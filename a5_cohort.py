"""A5: `replicon_cohort` TURETILMIS tablosunu kurar -- sutun EKLEMEZ.

NEDEN BU TABLO VAR. `about.html` su satiri yayimliyor:
"17.074 annotated replicons scanned, N carrying at least one". Bu bir
YAYGINLIK iddiasidir ve paydasi "taranan genomlar"dir.

A5'te indirilen 13.079 kayit, tam da ICINDE BIR RO OLDUGU BILINDIGI ICIN
secildi. Ayni havuza katilirsa o oran ~%100'e cikar ve bir KESIF gibi
gorunur -- oysa degisen sey olcum degil ORNEKLEME. Bu, projenin tekrar eden
hata biciminin genom duzeyindeki hali: karar veren sey (ornekleme), karar
verdigini soyledigimiz sey (yaygınlık) degil.

Cozum: her replikon hangi kohorttan geldigini tasir, ve genom duzeyinde her
yaygınlık olcumu kohortunu BEYAN etmek zorundadir.

  survey            : ozgun tarama havuzu. Yaygınlık YALNIZCA burada anlamli.
  protein_recovery  : A5'te, bilinen bir dizi icin hedeflenerek cekilen kayit.
                      Yaygınlık paydasina GIRMEZ; operon ve komsuluk
                      olcumlerine girer, cunku onlar kayit ici olcumlerdir.

Sutun yerine ayri tablo olmasinin sebebi: `replicon` CSV sozlesmesinden
(`genomic_context/replicons.csv`) yukleniyor; oraya sutun eklemek cikarici,
build_db ve tazelik izleyicisini ayni anda degistirmek demekti. Turetilmis
tablo tazeligini kendi girdisine -- kohort bildirimine -- baglar.

Kullanim:
    python3 a5_cohort.py --db roar.sqlite          # tabloyu kur
    python3 a5_cohort.py --write-manifest          # hedef listesini yaz
"""

import argparse
import os
import sqlite3
import sys

MANIFEST = os.path.join("cohorts", "a5_protein_recovery.tsv")

SCHEMA = """
DROP TABLE IF EXISTS replicon_cohort;
CREATE TABLE replicon_cohort (
    nucleotide_id TEXT PRIMARY KEY,
    cohort        TEXT NOT NULL,   -- survey | protein_recovery
    fetched_for   INTEGER NOT NULL DEFAULT 0  -- kac dizi icin hedeflendi
);
CREATE INDEX idx_rc_cohort ON replicon_cohort(cohort);
"""


def write_manifest(con, path=MANIFEST):
    """Hedef listesini depoya yazar: hangi kayit hangi dizi icin cekildi.

    Liste DB'den TURETILIR, indirilen dosyalardan degil. Boylece indirme
    yarim kalsa bile kohortun TANIMI tam kalir; eksik olan sey indirme olur
    ve ledger'dan sayilabilir.
    """
    rows = con.execute("""
        SELECT DISTINCT u.seq_md5, p.embl_nuc
        FROM ro_union u
        JOIN ro_union_source s ON s.seq_md5 = u.seq_md5 AND s.source = 'uniprot'
        JOIN ro_protein_status p ON p.accession = s.accession
        WHERE u.has_genomic_context = 0
          AND p.embl_nuc IS NOT NULL AND TRIM(p.embl_nuc) <> ''
    """).fetchall()
    pairs = set()
    for md5, field in rows:
        for tok in str(field).replace(";", " ").replace(",", " ").split():
            acc = tok.strip()
            if acc:
                pairs.add((acc, md5))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("accession\tseq_md5\n")
        for acc, md5 in sorted(pairs):
            fh.write(f"{acc}\t{md5}\n")
    accs = {a for a, _ in pairs}
    print(f"[manifest] {path}: {len(accs):,} kayit, {len(pairs):,} kayit-dizi bagi")
    return path


def load_manifest(path=MANIFEST):
    targeted = {}
    if not os.path.exists(path):
        return targeted
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[0]:
                targeted[parts[0]] = targeted.get(parts[0], 0) + 1
    return targeted


SURVEY = os.path.join("cohorts", "survey.txt")


def load_survey(path=SURVEY):
    """Ozgun tarama havuzunun ANLIK GORUNTUSU.

    AYIRT EDICI SEY BU LISTE, A5 HEDEF LISTESI DEGIL -- ve bunu bir olcum
    soyletti. A5'in hedefledigi 13.079 kaydin 344'u `replicon`da ZATEN
    vardi: survey havuzunda taranmislar, ama tasiduklari dizi on filtrenin
    kor noktasina dustugu icin baglamsiz sayilmislar. Kohortu "A5 listesinde
    mi" diye atamak, mesru taranmis 344 genomu yaygınlık paydasindan
    CIKARIRDI -- ters yonde ayni hata.

    Liste A5 kayitlari girmeden once alindi ve BUYUMEZ.
    """
    out = set()
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                out.add(line)
    return out


def build(con, manifest=MANIFEST, survey_path=SURVEY):
    targeted = load_manifest(manifest)
    if not targeted:
        raise SystemExit(f"{manifest} yok ya da bos -- once --write-manifest")
    survey = load_survey(survey_path)
    if not survey:
        raise SystemExit(f"{survey_path} yok -- survey kohortu tanimsiz, "
                         f"yaygınlık paydasi hesaplanamaz")
    con.executescript(SCHEMA)
    # Surum numarasi tasiyan kayitlar (NZ_CP049045.1) ile tasimayanlar
    # (CP049045) ayni kaydi gosterir; eslesme surum ATILARAK yapilir.
    base = {}
    for acc, n in targeted.items():
        base[acc.split(".")[0]] = base.get(acc.split(".")[0], 0) + n
    rows = []
    blank = 0
    for (nid,) in con.execute("SELECT nucleotide_id FROM replicon"):
        # KIMLIGI BOS SATIR (2026-11-02 olcumu: 1 tane, NZ_QXDB01000002.1.gbk,
        # status=no_annotation, 0 CDS). `extract_genomic_context.py`
        # anotasyonsuz dosya icin nucleotide_id'yi bos birakiyor. Kohort
        # atanamaz, cunku kohort kimlige gore tanimli; ve atanirsa anlik
        # goruntu ile tablo bir satir ayrisir. Sayilir ve DISARIDA birakilir.
        if nid is None or not str(nid).strip():
            blank += 1
            continue
        key = str(nid).split(".")[0]
        alt = key[3:] if key.startswith("NZ_") else key
        n = base.get(key) or base.get(alt) or 0
        # Survey listesinde olan kayit her kosulda survey'dir.
        cohort = "survey" if nid in survey else (
            "protein_recovery" if n else "survey")
        rows.append((nid, cohort, n))
    con.executemany("INSERT INTO replicon_cohort VALUES (?,?,?)", rows)
    con.commit()
    counts = dict(con.execute("SELECT cohort, COUNT(*) FROM replicon_cohort "
                              "GROUP BY 1").fetchall())
    both = con.execute("SELECT COUNT(*) FROM replicon_cohort "
                       "WHERE cohort='survey' AND fetched_for>0").fetchone()[0]
    if blank:
        print(f"[not] {blank} replikon kaydinin nucleotide_id'si BOS "
              f"(anotasyonsuz dosya); kohort disinda birakildi")
    print(f"[replicon_cohort] {len(rows):,} replikon: "
          + ", ".join(f"{k}={v:,}" for k, v in sorted(counts.items())))
    if both:
        print(f"[not] {both} kayit hem survey havuzunda hem A5 hedefi: "
              f"taranmislar ama tasiduklari dizi on filtrenin kor noktasina "
              f"dusmus. Survey sayilirlar; yaygınlık paydasi bozulmaz.")
    # Mevcut kayitlarin ANAHTARLARI, atama ile AYNI kurala gore uretilir:
    # surum atilir ve NZ_ oneki de ayrica anahtarlanir. Ayni kural
    # kullanilmazsa `replicon`da zaten bulunan 344 kayit "indirilmemis"
    # sayilir ve kurtarma orani 344 kayit kadar kotu gorunur.
    present = set()
    for (n,) in con.execute("SELECT nucleotide_id FROM replicon"):
        key = str(n).split(".")[0]
        present.add(key)
        if key.startswith("NZ_"):
            present.add(key[3:])
    missing = len(set(base) - present)
    if missing:
        print(f"[not] hedeflenen {missing:,} kayit `replicon`da YOK -- "
              f"indirilmemis ya da anotasyonsuz; kurtarma orani bu sayiyi "
              f"SAYMAK zorunda")
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--manifest", default=MANIFEST)
    ap.add_argument("--write-manifest", action="store_true")
    args = ap.parse_args()
    con = sqlite3.connect(args.db)
    if args.write_manifest:
        write_manifest(con, args.manifest)
        return 0
    build(con, args.manifest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
