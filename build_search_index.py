"""
Arama indeksi (SQLite FTS5) + filtre alanlari.

MEVCUT DURUM: arama bes alanda LIKE '%...%' ile yapiliyordu. 11.422 kayitta
calisiyor ama (a) kelime sinirini bilmiyor, bu yuzden "nah" araminca
"Mycobacterium" gibi alakasiz eslesmeler donuyor, (b) siralama yok, (c) birden
fazla kelime yazinca ("Pseudomonas naphthalene") hicbir sey bulunamiyor cunku
tam dizgi araniyor.

BU SCRIPT: her dogrulanmis RO icin tek satirlik bir arama belgesi kurar ve
FTS5 ile indeksler. Ayrica arayuzde filtre olarak kullanilacak alanlari
(kanit duzeyi, yasam alani, plazmit, operon ortagi, kimyasal aile) tek tabloda
toplar ki sorgu anda join gerekmesin.

Tablolar:
    ro_search      filtre alanlari (her dogrulanmis RO icin bir satir)
    ro_fts         FTS5 sanal tablosu; icerik ro_search'ten okunur

Arama sozdizimi FTS5'in kendisidir: bosluk = VE, "tam ifade", onek*.
"""

import argparse
import csv
import os
import sqlite3

SCHEMA = """
DROP TABLE IF EXISTS ro_search;
CREATE TABLE ro_search (
    rowid         INTEGER PRIMARY KEY,
    candidate_id  TEXT UNIQUE,
    protein_id    TEXT,
    organism      TEXT,
    genus         TEXT,
    product       TEXT,
    cluster       TEXT,
    gene          TEXT,
    leaf_id       TEXT,
    substrate     TEXT,
    family        TEXT,     -- kimyasal aile (chemistry.csv)
    reaction      TEXT,     -- reaksiyon sinifi
    tier          TEXT,     -- kanit duzeyi
    domain        TEXT,     -- Bacteria / Eukaryota / ...
    is_plasmid    INTEGER,
    has_partner   INTEGER,  -- operonda beta/ferredoksin/reduktaz var mi
    has_regulator INTEGER,  -- operonun 5' ucunda divergent duzenleyici
    doc           TEXT      -- indekslenecek birlesik metin
);
CREATE INDEX idx_search_cluster ON ro_search(cluster);
CREATE INDEX idx_search_tier    ON ro_search(tier);
CREATE INDEX idx_search_family  ON ro_search(family);
CREATE INDEX idx_search_domain  ON ro_search(domain);

DROP TABLE IF EXISTS ro_fts;
CREATE VIRTUAL TABLE ro_fts USING fts5(
    doc,
    content='ro_search', content_rowid='rowid', tokenize='unicode61'
);
"""


def read_csv_map(path, key="cluster"):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r[key]: r for r in csv.DictReader(fh)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--ecology", default="cluster_ecology.csv")
    ap.add_argument("--chemistry", default="chemistry.csv")
    args = ap.parse_args()

    eco = read_csv_map(args.ecology)
    chem = read_csv_map(args.chemistry)

    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    con.executescript(SCHEMA)

    def has(table):
        return con.execute("SELECT name FROM sqlite_master WHERE name=?", (table,)).fetchone()

    rows = con.execute(f"""
        SELECT r.candidate_id, r.protein_id, r.locus_tag, r.gene, r.product, r.ro_cluster,
               p.organism, p.taxonomy, p.is_plasmid,
               {'rl.leaf_id' if has('ro_leaf') else 'NULL'} leaf_id,
               {'e.tier' if has('ro_evidence') else 'NULL'} tier,
               {'d.domain' if has('ro_domain') else 'NULL'} domain,
               {'o.has_beta + o.has_ferredoxin + o.has_reductase' if has('operon') else '0'} partners,
               {"g.architecture" if has('ro_regulation') else 'NULL'} architecture
        FROM ro r
        JOIN replicon p ON p.nucleotide_id = r.nucleotide_id
        {'LEFT JOIN ro_leaf rl ON rl.candidate_id = r.candidate_id' if has('ro_leaf') else ''}
        {'LEFT JOIN ro_evidence e ON e.candidate_id = r.candidate_id' if has('ro_evidence') else ''}
        {'LEFT JOIN ro_domain d ON d.candidate_id = r.candidate_id' if has('ro_domain') else ''}
        {'LEFT JOIN operon o ON o.candidate_id = r.candidate_id' if has('operon') else ''}
        {'LEFT JOIN ro_regulation g ON g.candidate_id = r.candidate_id' if has('ro_regulation') else ''}
        WHERE r.is_confirmed = 1""").fetchall()

    out = []
    for i, r in enumerate(rows, 1):
        cluster = r["ro_cluster"] or ""
        gene = cluster.split("_", 2)[-1] if cluster else ""
        ch = chem.get(cluster, {})
        substrate = ch.get("substrate_en") or eco.get(cluster, {}).get("substrate", "")
        organism = r["organism"] or ""
        genus = organism.split()[0] if organism else ""
        # Indekslenecek metin: arayan kisinin yazabilecegi her sey tek belgede.
        doc = " ".join(filter(None, [
            r["candidate_id"], r["protein_id"], r["locus_tag"], r["gene"], r["product"],
            cluster, gene, organism, r["taxonomy"], substrate,
            ch.get("product_en", ""), ch.get("reaction", ""),
            r["leaf_id"] or "", r["tier"] or "", r["domain"] or "",
            "plasmid" if r["is_plasmid"] else "chromosome",
        ]))
        out.append((i, r["candidate_id"], r["protein_id"], organism, genus, r["product"],
                    cluster, gene, r["leaf_id"], substrate, ch.get("family", ""),
                    ch.get("reaction_class", ""), r["tier"], r["domain"],
                    1 if r["is_plasmid"] else 0, 1 if (r["partners"] or 0) > 0 else 0,
                    1 if r["architecture"] == "divergent_regulator" else 0, doc))

    con.executemany("INSERT INTO ro_search VALUES (" + ",".join("?" * 18) + ")", out)
    con.execute("INSERT INTO ro_fts(ro_fts) VALUES('rebuild')")
    con.commit()

    n = con.execute("SELECT COUNT(*) FROM ro_search").fetchone()[0]
    print(f"[ro_search] {n} satir")
    for field in ("tier", "domain", "family"):
        dist = con.execute(f"SELECT {field}, COUNT(*) FROM ro_search GROUP BY 1 "
                           f"ORDER BY 2 DESC LIMIT 6").fetchall()
        print(f"   {field:8s}", ", ".join(f"{r[0] or '-'}:{r[1]}" for r in dist))
    for probe in ("naphthalene", "Pseudomonas putida", "plant OR Viridiplantae", "nah*"):
        try:
            hits = con.execute("SELECT COUNT(*) FROM ro_fts WHERE ro_fts MATCH ?",
                               (probe,)).fetchone()[0]
            print(f"   test '{probe}': {hits} sonuc")
        except sqlite3.OperationalError as exc:
            print(f"   test '{probe}': sorgu hatasi ({exc})")
    con.close()


if __name__ == "__main__":
    main()
