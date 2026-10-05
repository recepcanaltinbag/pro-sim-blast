"""
ROAR yerel veritabanini kur (SQLite).

TASARIM ILKESI: ham komsu saklanir, siniflandirma SORGU ANINDA yapilir.

Mevcut pipeline is_matching_product() ile anahtar kelime filtresini CIKARIM
asamasinda uyguluyor -- "hypothetical protein" komsular veriye hic girmiyor.
Oysa hedef bilinmeyen genlere fonksiyon atfetmek. Ayrica regulator tanimini
degistirmek istedigin anda 17.000 gbk dosyasini bastan parse etmen gerekiyor.

Burada neighbor tablosu HAM tutulur, gene_category ayri bir tablodur.
Siniflandirmayi degistirmek = tek bir tabloyu yeniden doldurmak.

Tablolar:
    replicon        genom/plazmit kaydi, taksonomi, anotasyon durumu
    ro              dogrulanmis RO alpha subunit (koordinat + STRAND + kume)
    neighbor        +-10kb icindeki tum CDS'ler, filtresiz, STRAND'li
    gene_category   komsu -> kategori eslemesi (yeniden uretilebilir)
"""

import argparse
import csv
import os
import re
import sqlite3
import sys

csv.field_size_limit(10 ** 7)

SCHEMA = """
PRAGMA journal_mode = WAL;

CREATE TABLE IF NOT EXISTS replicon (
    nucleotide_id TEXT PRIMARY KEY,
    file          TEXT,
    organism      TEXT,
    taxonomy      TEXT,
    description   TEXT,
    length        INTEGER,
    mol_type      TEXT,
    is_plasmid    INTEGER,
    is_circular   INTEGER,
    cds_count     INTEGER,
    -- 'ok' | 'no_annotation'. Anotasyonsuz replikonlar SESSIZCE bos
    -- donmemeli; "komsusuz RO" ile "verisi eksik" ayri seylerdir.
    status        TEXT
);

CREATE TABLE IF NOT EXISTS ro (
    candidate_id  TEXT PRIMARY KEY,
    nucleotide_id TEXT REFERENCES replicon(nucleotide_id),
    start         INTEGER,
    end           INTEGER,
    strand        INTEGER,        -- 1 / -1 : orientation analizinin temeli
    product       TEXT,
    locus_tag     TEXT,
    protein_id    TEXT,
    gene          TEXT,
    -- HMM dogrulama sonuclari (annotate_ro.py dolduruyor)
    ro_cluster    TEXT,
    ro_group      TEXT,
    model_coverage REAL,
    hmm_score     REAL,
    hmm_evalue    REAL,
    rieske_intact    INTEGER,
    catalytic_intact INTEGER,
    is_confirmed     INTEGER
);

CREATE TABLE IF NOT EXISTS neighbor (
    neighbor_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    candidate_id  TEXT REFERENCES ro(candidate_id),
    nucleotide_id TEXT,
    start         INTEGER,
    end           INTEGER,
    strand        INTEGER,
    distance      INTEGER,    -- isaretli: negatif=sol, pozitif=sag, 0=ortusuyor
    same_strand   INTEGER,    -- RO ile ayni yonde mi
    gene_offset   INTEGER,    -- kac gen uzakta (-3 = uc gen solda)
    spans_origin  INTEGER,    -- dairesel replikonda orijini asiyor mu
    product       TEXT,
    locus_tag     TEXT,
    protein_id    TEXT,
    gene          TEXT
);

CREATE TABLE IF NOT EXISTS gene_category (
    neighbor_id INTEGER REFERENCES neighbor(neighbor_id),
    category    TEXT,
    method      TEXT
);

CREATE INDEX IF NOT EXISTS idx_nb_cand   ON neighbor(candidate_id);
CREATE INDEX IF NOT EXISTS idx_nb_dist   ON neighbor(distance);
CREATE INDEX IF NOT EXISTS idx_nb_offset ON neighbor(gene_offset);
CREATE INDEX IF NOT EXISTS idx_cat_nb    ON gene_category(neighbor_id);
CREATE INDEX IF NOT EXISTS idx_cat_cat   ON gene_category(category);
CREATE INDEX IF NOT EXISTS idx_ro_nuc    ON ro(nucleotide_id);
CREATE INDEX IF NOT EXISTS idx_ro_group  ON ro(ro_group);
"""

# Kategori tanimlari. Sirali degerlendirilir: ilk eslesen kazanir, boylece
# "transposase" iceren bir urun "hypothetical" kovasina dusmez.
# Tanimi degistirmek icin sadece burayi duzenleyip --recategorize calistir.
CATEGORIES = [
    ("transposon", r"transposase|transposon|insertion sequence|IS[0-9]{2,4}\b|"
                   r"integrase|recombinase|resolvase|invertase|mobile element"),
    ("regulator",  r"transcriptional regulator|regulatory protein|LysR|TetR|AraC|"
                   r"GntR|MarR|IclR|LuxR|XylR|NtrC|sigma factor|"
                   r"helix-turn-helix|DNA-binding response regulator|repressor|activator"),
    # ro_beta, ro_alpha'dan ONCE: "ring-hydroxylating dioxygenase subunit beta"
    # aksi halde ring.hydroxylating ile alpha kovasina duser (1.529 yanlis atama olcumlendi).
    ("ro_beta",    r"dioxygenase.*(beta|small subunit)|(beta|small subunit).*dioxygenase|"
                   r"oxygenase.*(beta|small) subunit"),
    ("ro_alpha",   r"ring.hydroxylating|dioxygenase.*(alpha|large subunit)|"
                   r"(alpha|large subunit).*dioxygenase|oxygenase.*large subunit"),
    ("ferredoxin", r"ferredoxin|2Fe-2S|\[2Fe-2S\]|rieske"),
    ("reductase",  r"reductase|oxidoreductase|NADH|FAD|flavoprotein"),
    ("transporter", r"transporter|permease|ABC.*binding|efflux|porin|MFS"),
    ("dehydrogenase", r"dehydrogenase"),
    ("hydrolase",  r"hydrolase|esterase|lipase|amidase"),
    ("ring_cleavage", r"catechol|muconate|muconolactone|protocatechuate|"
                      r"hydroxymuconate|semialdehyde|extradiol|intradiol"),
    ("hypothetical", r"hypothetical|DUF[0-9]+|uncharacterized|unknown"),
]
COMPILED = [(name, re.compile(pattern, re.I)) for name, pattern in CATEGORIES]

# Regulator AILESI -- 'regulator' kategorisinin alt kirilimi.
# Ayri bir boyut olarak saklanir (method='regulator_family_v1'), cunku
# "hangi RO grubu hangi regulator ailesiyle gidiyor" sorusu tek kovayla
# cevaplanamaz. LysR'lar klasik olarak divergent oturur; TetR'lar genelde
# ko-direksiyoneldir -- aileyi bilmek orientation analizini yorumlanabilir kilar.
REGULATOR_FAMILIES = [
    ("LysR",   r"LysR"),
    ("TetR",   r"TetR|AcrR"),
    ("AraC",   r"AraC|XylS"),
    ("MarR",   r"MarR"),
    ("IclR",   r"IclR"),
    ("GntR",   r"GntR|FadR"),
    ("LuxR",   r"LuxR"),
    ("ArsR",   r"ArsR"),
    ("Crp_Fnr", r"Crp|Fnr|CRP/FNR"),
    ("sigma54", r"sigma[- ]?54|NtrC|Fis family"),
    ("two_component", r"response regulator|two-component|sensor histidine"),
    ("sigma_factor",  r"sigma factor|sigma-70|RNA polymerase.*sigma"),
]
COMPILED_FAMILIES = [(name, re.compile(pattern, re.I))
                     for name, pattern in REGULATOR_FAMILIES]


def categorize(product):
    """Bir urun adini kategoriye ata. Eslesme yoksa 'other'."""
    for name, pattern in COMPILED:
        if pattern.search(product or ""):
            return name
    return "other"


def regulator_family(product):
    """Regulator ailesini cikar. Aile belirlenemezse 'regulator_unclassified'."""
    for name, pattern in COMPILED_FAMILIES:
        if pattern.search(product or ""):
            return name
    return "regulator_unclassified"


def to_int(value):
    """CSV'den gelen 'True'/'False'/'' degerlerini 0/1/None'a cevir."""
    if value in ("True", "true", "1"):
        return 1
    if value in ("False", "false", "0"):
        return 0
    if value in ("", None):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def load_csv(connection, path, table, columns, transform=None, batch=20000):
    """Bir CSV'yi tabloya yukle. Donen: yuklenen satir sayisi."""
    if not os.path.exists(path):
        print(f"[atlandi] {path} yok")
        return 0

    placeholders = ",".join("?" * len(columns))
    statement = f"INSERT OR REPLACE INTO {table} ({','.join(columns)}) VALUES ({placeholders})"
    rows, total = [], 0

    with open(path, newline="") as handle:
        for record in csv.DictReader(handle):
            if transform:
                record = transform(record)
                if record is None:
                    continue
            rows.append([record.get(column) for column in columns])
            if len(rows) >= batch:
                connection.executemany(statement, rows)
                total += len(rows)
                rows = []
                sys.stdout.write(f"\r  {table}: {total}")
                sys.stdout.flush()
    if rows:
        connection.executemany(statement, rows)
        total += len(rows)
    connection.commit()
    print(f"\r  {table}: {total} satir")
    return total


def build(context_dir, db_path):
    if os.path.exists(db_path):
        os.remove(db_path)
    connection = sqlite3.connect(db_path)
    connection.executescript(SCHEMA)

    print("[yukleniyor] replicon")
    load_csv(connection, os.path.join(context_dir, "replicons.csv"), "replicon",
             ["nucleotide_id", "file", "organism", "taxonomy", "description",
              "length", "mol_type", "is_plasmid", "is_circular", "cds_count", "status"],
             transform=lambda r: {**r,
                                  "is_plasmid": to_int(r.get("is_plasmid")),
                                  "is_circular": to_int(r.get("is_circular")),
                                  "length": to_int(r.get("length")),
                                  "cds_count": to_int(r.get("cds_count"))})

    print("[yukleniyor] ro")
    load_csv(connection, os.path.join(context_dir, "ro_candidates.csv"), "ro",
             ["candidate_id", "nucleotide_id", "start", "end", "strand",
              "product", "locus_tag", "protein_id", "gene"],
             transform=lambda r: {**r, "start": to_int(r.get("start")),
                                  "end": to_int(r.get("end")),
                                  "strand": to_int(r.get("strand"))})

    print("[yukleniyor] neighbor")
    load_csv(connection, os.path.join(context_dir, "neighbors.csv"), "neighbor",
             ["candidate_id", "nucleotide_id", "start", "end", "strand",
              "distance", "same_strand", "gene_offset", "spans_origin",
              "product", "locus_tag", "protein_id", "gene"],
             transform=lambda r: {**r,
                                  "start": to_int(r.get("start")),
                                  "end": to_int(r.get("end")),
                                  "strand": to_int(r.get("strand")),
                                  "distance": to_int(r.get("distance")),
                                  "same_strand": to_int(r.get("same_strand")),
                                  "gene_offset": to_int(r.get("gene_offset")),
                                  "spans_origin": to_int(r.get("spans_origin"))})

    categorize_all(connection)
    connection.close()
    print(f"[bitti] {db_path}")


def categorize_all(connection):
    """gene_category tablosunu ürün adlarindan yeniden uret."""
    print("[siniflandiriliyor] gene_category")
    connection.execute("DELETE FROM gene_category")
    cursor = connection.execute("SELECT neighbor_id, product FROM neighbor")
    rows, total = [], 0
    while True:
        chunk = cursor.fetchmany(50000)
        if not chunk:
            break
        for neighbor_id, product in chunk:
            category = categorize(product)
            rows.append((neighbor_id, category, "regex_v1"))
            # Regulator ise aileyi de AYRI bir satir olarak ekle. Boylece
            # kategori ve aile bagimsiz sorgulanabilir; biri digerini ezmez.
            if category == "regulator":
                rows.append((neighbor_id, regulator_family(product),
                             "regulator_family_v1"))
        connection.executemany(
            "INSERT INTO gene_category (neighbor_id, category, method) VALUES (?,?,?)", rows)
        total += len(rows)
        rows = []
        sys.stdout.write(f"\r  gene_category: {total}")
        sys.stdout.flush()
    connection.commit()
    print(f"\r  gene_category: {total} satir")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-dir", default="genomic_context")
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--recategorize", action="store_true",
                        help="Sadece gene_category'yi yeniden uret, veriyi tekrar yukleme")
    args = parser.parse_args()

    if args.recategorize:
        connection = sqlite3.connect(args.db)
        categorize_all(connection)
        connection.close()
        return

    build(args.context_dir, args.db)


if __name__ == "__main__":
    main()
