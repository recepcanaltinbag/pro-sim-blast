"""INGEST ACIGI: kimligi belli ama indirilmemis nukleotid kayitlari.

NEDEN BU BETIK VAR. Kullanici `5_505_qxyA` icin "referansa daha yakin
homolog cikmasi lazim" dedi ve genel bir kusur olup olmadigini sordu.
Olculdu, uc katman cikti ve yalnizca ucuncusu gercek bir acik:

  1. KIMLIK HESABI SAGLAM. Dogrulanmis 15.960 dizide qxyA'ya en iyi
     isabet (kendisi haric) %45,5 -- `ro_evidence` ile birebir ayni.
  2. YAKIN HOMOLOGLAR VERIDE VAR ama uye SAYILMIYOR. Birlesik kumede
     %97,4 ve %97,1 kimlikli iki dizi kapiyi geciyor, Rieske saglam,
     ucgen 3/3, ve zaten `5_505_qxyA` atanmis -- ama
     `has_genomic_context=0` oldugu icin ne uye sayisina ne kademe
     tablosuna giriyorlar. Bu mimari bir karar, kusur degil; tip sayfasi
     birlesim sayisini zaten acikliyor.
  3. GENOMLARI INDIRILMEMIS -- GERCEK ACIK BU. %97,4'luk homolog
     `FNTJ01000003` contig'inde, `Pseudomonas saponiphila`'nin FNTJ01
     derlemesinde. Bizde o turden tek kayit var ve o BASKA bir
     derlemeden (`NZ_JAOXCB010000053.1`). Yani ingest CONTIG BAZLI:
     ilk aramada isabet veren contig alinmis, kardesleri ve alternatif
     derlemeler alinmamis.

NE OLCER. Bagiglami olmayan, kapiyi gecen ve bir tipe atanmis her dizinin
UniProt EMBL capraz referansi okunur; o nukleotid kaydi `replicon`
tablosunda yoksa INDIRME LISTESINE girer. Yani liste tahmin degil:
her satirin arkasinda, kimligi bilinen bir kayitta duran, kapimizi
gecmis ve tiplenmis en az bir dizi var.

NE YAPMAZ. Indirmez. 12.312 kayit indirmek yayimlanan her sayiyi
degistirir (yeni veri surumu) ve bu kuratorun karari.

Cikti: analysis_out/ingest_gap.csv + ingest_gap.json
"""

import argparse
import csv
import json
import os
import sqlite3
import sys
from collections import Counter, defaultdict


def base_accession(accession):
    """`NZ_FNTJ01000003.1` -> `FNTJ01000003`. Surum ve RefSeq oneki atilir."""
    return accession.split(".")[0].replace("NZ_", "")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-csv", default="analysis_out/ingest_gap.csv")
    parser.add_argument("--out-json", default="analysis_out/ingest_gap.json")
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    have = {base_accession(r[0])
            for r in con.execute("SELECT nucleotide_id FROM replicon")}

    rows = con.execute("""
        SELECT x.identifier, u.ro_cluster, u.seq_md5
          FROM ro_xref x
          JOIN ro_union u USING(seq_md5)
         WHERE x.db = 'embl_nucleotide'
           AND u.passes_gate = 1
           AND u.has_genomic_context = 0
           AND u.ro_cluster IS NOT NULL""").fetchall()

    per_record = defaultdict(set)
    record_types = defaultdict(set)
    per_type = Counter()
    already = 0
    for accession, cluster, md5 in rows:
        base = base_accession(accession)
        if base in have:
            already += 1
            continue
        per_record[base].add(md5)
        record_types[base].add(cluster)
        per_type[cluster] += 1

    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    with open(args.out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["nucleotide_accession", "sequences_gained", "types"])
        for record, seqs in sorted(per_record.items(),
                                   key=lambda kv: (-len(kv[1]), kv[0])):
            writer.writerow([record, len(seqs),
                             ";".join(sorted(record_types[record]))])

    report = {
        "question": ("which nucleotide records are named in our own "
                     "cross-references but were never downloaded?"),
        "why_it_matters": (
            "a sequence without genomic context has no operon, no "
            "neighbourhood and no replicon, so it is counted in no per-type "
            "member total and gets no evidence tier. When the record holding "
            "it is known by accession and simply absent, that is a gap in "
            "the ingest rather than in the data"),
        "how_the_list_is_built": (
            "every sequence that passes the gate, carries a type and has no "
            "genomic context is looked up through its UniProt EMBL "
            "cross-reference; if that nucleotide record is not in the "
            "replicon table it enters the list. Each row therefore stands on "
            "an identified record holding at least one sequence this "
            "database has already accepted and typed"),
        "does_not_download": (
            "nothing is fetched here. Taking 12,312 records changes every "
            "published number and is a data revision, which is a curator's "
            "decision"),
        "records_to_fetch": len(per_record),
        "sequences_that_would_gain_context": sum(
            len(v) for v in per_record.values()),
        "cross_references_already_held": already,
        "gain_by_type": dict(per_type.most_common()),
    }
    with open(args.out_json, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print(f"[kayit] indirilecek {len(per_record)}, "
          f"kazanilacak dizi {report['sequences_that_would_gain_context']}")
    print(f"[tip]   en cok: {list(per_type.most_common(5))}")
    print(f"[yazildi] {args.out_csv}  ve  {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
