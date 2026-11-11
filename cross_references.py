"""Capraz referanslar ALAN olur: EC, UniProt, RefSeq, EMBL nukleotit.

NEDEN BU DOSYA VAR
  ROADMAP A6: "Bir toparlama veritabaninin isi ayni molekulun her yerdeki
  adini bir arada tutmaktir." Dort tanimlayici 2026-11-01'den beri
  `ro_protein_status`'ta DURUYOR ama sitede HICBIR YERDE gorunmuyordu:
  okuyucu bir diziyi UniProt'ta, ENA'da ya da RefSeq'te aramak istedigi
  anda elinde hicbir sey yok.

NEDEN DOGRUDAN SUTUN BASMIYORUZ
  Tanimlayicinin biri -- EC -- bir ISLEV IDDIASI, oteki ucu yalnizca
  adres. Ve bu iddianin kaynagi olculdu: EC tasiyan 3.880 accession'in
  3.788'i **UniProtKB unreviewed (TrEMBL)**, yalnizca 92'si Swiss-Prot.
  Yani EC'lerin %97,6'si elle dogrulanmis bir olcum degil, otomatik
  anotasyonun yaydigi bir TAHMIN. Bir tip sayfasina "EC 1.14.15.7" diye
  yazmak, bu veritabaninin kendi kuratorlu substrat etiketiyle ayni
  kefeye konmus 3.788 islev iddiasi uretirdi.

  O yuzden EC, KOSULUYLA yayimlanir: her satir reviewed/unreviewed
  etiketini tasir, ve dis anotasyonun BIZIM kuratorlu tipimizle ne kadar
  anlastigi tip basina olculur.

PAYDA (madde A3)
  Bu tanimlayicilarin tamami UniProt kaynakli dizilerden gelir. Payda
  dolayisiyla "UniProt kaynagi olan dizi" (68.747), TUM diziler (73.299)
  DEGIL. Iki payda JSON'da ayri ayri yazilir; birini otekinin yerine
  kullanan bir yuzde sessizce yanlis olur.

BIRIM
  Satir = DIZI (`seq_md5`), accession degil (madde A1). Ayni dizi birkac
  accession altinda gecebilir; tanimlayicilar o zaman birlesir.

Cikti:
    sqlite tablosu  ro_xref(seq_md5, db, identifier)
    analysis_out/cross_references.json
"""

import argparse
import json
import os
import re
import sqlite3
from collections import Counter, defaultdict

#  Rieske oksijenazin beklenen EC sinifi: 1.14.* (eslesmis donor, molekuler
#  oksijen katilimi). Bunun disindaki bir EC, dizinin kendisinin baska bir
#  sey oldugunu KANITLAMAZ -- UniProt kaydi cok alanli bir protein olabilir
#  ve EC kaydin tamamina aittir. O yuzden "yanlis" degil, "beklenen sinifin
#  disinda" diye raporlanir.
EXPECTED_EC_PREFIX = "1.14."

#  Tanimlayici bicimleri. Dogrulayici bunlari ayrica sinar; burada bozuk
#  bir deger TABLOYA GIRMEZ, cunku sayfada tiklanabilir bir baglanti olacak
#  ve kirik bir adres yayimlamak veriyi yayimlamamaktan kotudur.
PATTERNS = {
    "ec": re.compile(r"^\d+\.(\d+|-)\.(\d+|-)\.(\d+|-)$"),
    "refseq_protein": re.compile(r"^[A-Z]{2}_\d+(\.\d+)?$"),
    #  WGS accession'i 6 harf + 9 hane olabiliyor (JACBZF010000002); ilk
    #  desen \d{5,8} idi ve 5.358 GECERLI kaydi bozuk sayip atiyordu. Bicim
    #  kontrolunun atilan sayiyi RAPORLAMASI bunu yakaladi.
    "embl_nucleotide": re.compile(r"^[A-Z]{1,6}\d{5,10}(\.\d+)?$"),
    "uniprot": re.compile(r"^[A-Z0-9]{6,10}(-\d+)?$"),
}

URLS = {
    "ec": "https://enzyme.expasy.org/EC/{}",
    "refseq_protein": "https://www.ncbi.nlm.nih.gov/protein/{}",
    "embl_nucleotide": "https://www.ebi.ac.uk/ena/browser/view/{}",
    "uniprot": "https://www.uniprot.org/uniprotkb/{}",
}


def split_values(raw):
    """`ec` ve `refseq` alanlari virgulle cok degerli olabiliyor (7 satir)."""
    if not raw:
        return []
    return [v.strip() for v in re.split(r"[,;]", raw) if v.strip()]


def collect(con):
    """seq_md5 -> {db -> set(identifier)} ve EC'nin inceleme durumu."""
    xrefs = defaultdict(lambda: defaultdict(set))
    ec_review = defaultdict(set)      # "reviewed"/"unreviewed" -> seq_md5
    rejected = Counter()

    rows = con.execute("""
        SELECT s.seq_md5, s.accession, p.entry_type, p.refseq, p.embl_nuc, p.ec
          FROM ro_union_source s
          JOIN ro_protein_status p ON p.accession = s.accession
         WHERE s.source = 'uniprot'
    """)
    for seq_md5, accession, entry_type, refseq, embl, ec in rows:
        pairs = [("uniprot", [accession]),
                 ("refseq_protein", split_values(refseq)),
                 ("embl_nucleotide", split_values(embl)),
                 ("ec", split_values(ec))]
        for db, values in pairs:
            for value in values:
                if PATTERNS[db].match(value):
                    xrefs[seq_md5][db].add(value)
                else:
                    rejected[db] += 1
        if split_values(ec):
            #  Swiss-Prot ile TrEMBL arasindaki fark, EC'nin olcum mu tahmin
            #  mi oldugunu belirleyen TEK isaret; kaybedilemez.
            kind = "reviewed" if "reviewed" in (entry_type or "") \
                   and "unreviewed" not in (entry_type or "") else "unreviewed"
            ec_review[kind].add(seq_md5)
    return xrefs, ec_review, rejected


def write_table(con, xrefs):
    con.execute("DROP TABLE IF EXISTS ro_xref")
    con.execute("""
        CREATE TABLE ro_xref (
            seq_md5    TEXT NOT NULL,
            db         TEXT NOT NULL,   -- uniprot | refseq_protein |
                                        -- embl_nucleotide | ec
            identifier TEXT NOT NULL,
            PRIMARY KEY (seq_md5, db, identifier)
        )""")
    con.execute("CREATE INDEX idx_xref_md5 ON ro_xref(seq_md5)")
    con.execute("CREATE INDEX idx_xref_db  ON ro_xref(db, identifier)")
    con.executemany("INSERT OR IGNORE INTO ro_xref VALUES (?,?,?)",
                    ((md5, db, ident)
                     for md5, by_db in xrefs.items()
                     for db, idents in by_db.items()
                     for ident in idents))
    con.commit()


def ec_against_curation(con, xrefs):
    """Dis EC, BIZIM kuratorlu tipimizle anlasiyor mu -- tip basina."""
    cluster_of = dict(con.execute(
        "SELECT seq_md5, ro_cluster FROM ro_union "
        "WHERE ro_cluster IS NOT NULL AND ro_cluster <> ''"))
    per_cluster = defaultdict(Counter)
    for md5, by_db in xrefs.items():
        cluster = cluster_of.get(md5)
        if not cluster:
            continue
        for ec in by_db.get("ec", ()):
            per_cluster[cluster][ec] += 1

    out = []
    for cluster, counts in per_cluster.items():
        total = sum(counts.values())
        modal, modal_n = counts.most_common(1)[0]
        outside = sorted(ec for ec in counts
                         if not ec.startswith(EXPECTED_EC_PREFIX))
        out.append({
            "cluster": cluster,
            "ec_readings": total,
            "distinct_ec": len(counts),
            "modal_ec": modal,
            #  Tipin EN SIK EC'si, o tipteki EC okumalarinin kacini
            #  kapsiyor. 1,0 demek dis anotasyon tip icinde tutarli;
            #  kucuk deger demek ayni tipe birkac farkli islev yazilmis.
            "modal_share": round(modal_n / total, 3),
            "ec_outside_expected_class": outside,
        })
    out.sort(key=lambda r: -r["ec_readings"])
    return out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out", default=os.path.join("analysis_out",
                                                      "cross_references.json"))
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    xrefs, ec_review, rejected = collect(con)
    write_table(con, xrefs)

    all_sequences = con.execute("SELECT COUNT(*) FROM ro_union").fetchone()[0]
    uniprot_sourced = con.execute(
        "SELECT COUNT(DISTINCT seq_md5) FROM ro_union_source "
        "WHERE source = 'uniprot'").fetchone()[0]

    coverage = {}
    for db in PATTERNS:
        n = sum(1 for by_db in xrefs.values() if by_db.get(db))
        coverage[db] = {
            "sequences": n,
            "share_of_uniprot_sourced": round(n / uniprot_sourced, 4)
            if uniprot_sourced else None,
            "share_of_all_sequences": round(n / all_sequences, 4)
            if all_sequences else None,
            "rejected_malformed": rejected.get(db, 0),
            "url_template": URLS[db],
        }

    by_type = ec_against_curation(con, xrefs)
    ec_counts = Counter(ec for by_db in xrefs.values()
                        for ec in by_db.get("ec", ()))
    outside = sorted(((ec, n) for ec, n in ec_counts.items()
                      if not ec.startswith(EXPECTED_EC_PREFIX)),
                     key=lambda r: -r[1])
    n_reviewed = len(ec_review.get("reviewed", ()))
    n_unreviewed = len(ec_review.get("unreviewed", ()) - ec_review.get("reviewed", set()))

    result = {
        "method": {
            "question": ("which external identifiers exist for each sequence, "
                         "and whether the external functional annotation agrees "
                         "with this database's own curation"),
            "unit": "one amino acid sequence (seq_md5), not one accession",
            "denominator": ("every identifier here comes from a UniProt-sourced "
                            "sequence, so the denominator is the UniProt-sourced "
                            "set, not the whole union"),
            "ec_is_a_claim": ("EC is a statement about function, while the other "
                              "three are addresses. Almost every EC here is "
                              "carried by an unreviewed (TrEMBL) entry, which "
                              "means it was propagated by automatic annotation "
                              "rather than measured, so it is published with its "
                              "review status attached and is never treated as "
                              "equivalent to this database's curated substrate"),
            "expected_ec_class": (f"{EXPECTED_EC_PREFIX}* (paired donors, "
                                  "incorporation of molecular oxygen); an EC "
                                  "outside it does not prove the sequence is "
                                  "something else, because a UniProt entry can "
                                  "be multi-domain and the EC belongs to the "
                                  "whole record"),
            "malformed_identifiers": ("an identifier that fails its format "
                                      "pattern is NOT written to the table, "
                                      "because every one becomes a clickable "
                                      "link and a broken address is worse than "
                                      "a missing one"),
        },
        "denominators": {
            "all_sequences": all_sequences,
            "uniprot_sourced_sequences": uniprot_sourced,
        },
        "coverage": coverage,
        "ec_review_status": {
            "sequences_with_a_reviewed_ec": n_reviewed,
            "sequences_with_only_an_unreviewed_ec": n_unreviewed,
            "reviewed_share": round(n_reviewed / (n_reviewed + n_unreviewed), 4)
            if (n_reviewed + n_unreviewed) else None,
        },
        "ec_outside_expected_class": [{"ec": ec, "sequences": n}
                                      for ec, n in outside],
        "ec_by_type": by_type,
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=1)

    print("=" * 72)
    print("CAPRAZ REFERANSLAR")
    print("=" * 72)
    print(f"  tum dizi                      {all_sequences}")
    print(f"  UniProt kaynakli dizi         {uniprot_sourced}  <- PAYDA")
    for db, info in coverage.items():
        print(f"  {db:<28}{info['sequences']:>6}"
              f"  ({100 * (info['share_of_uniprot_sourced'] or 0):.1f} % paydadan)"
              + (f"  [{info['rejected_malformed']} bozuk atildi]"
                 if info["rejected_malformed"] else ""))
    print(f"  EC'si incelenmis dizi         {n_reviewed}")
    print(f"  EC'si yalnizca incelenmemis   {n_unreviewed}")
    print(f"  beklenen sinif disi EC        {len(outside)} farkli")
    for ec, n in outside[:6]:
        print(f"    {ec:<14}{n} dizi")
    print("  tip basina EC tutarliligi (en kalabalik 8):")
    for row in by_type[:8]:
        print(f"    {row['cluster']:<16}{row['ec_readings']:>5} okuma, "
              f"{row['distinct_ec']} farkli EC, "
              f"en sik {row['modal_ec']} (%{100 * row['modal_share']:.0f})")
    print(f"\n[yazildi] {args.out}  + ro_xref tablosu")
    con.close()


if __name__ == "__main__":
    main()
