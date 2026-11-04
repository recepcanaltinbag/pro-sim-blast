"""
A5 dizi eslesme kapisi -- bir kayit bir diziye hangi sartla baglam olur?

SORUN: A5'te her dizi icin NCBI capraz referansindan bir nukleotid kaydi
geldi. Capraz referans "bu dizi su kayitta gecer" der. BU, "su kayitta ayni
proteindir" DEMEZ. Kaydi acip oradaki bir adayi almak, dizinin komsulugunu
bir BASKA genin komsulugu olarak yayimlama riskidir.

KURAL (tek kural, gevsetilmemis): bir kayit bir diziye ancak AYNI KAYITTA
cevirisinin md5'i BIREBIR esit bir aday varsa baglam olur.

OLCUM (2026-11-04, 18.636 hedef cift / 18.532 tekil dizi):
  ayni kayitta birebir esleme : 18.139  %97,3   -> baglam kurulur
  baska kayitta esleme        :      9  %0,0    -> bu cift icin kurulmaz
  hic esleme yok              :    488  %2,6

Kalan 497 ciftin nedeni (bu betik yeniden uretir):
  A kayit indirilemedi            16   %3,2
  B kayitta hic CDS anotasyonu yok 50  %10,1   <- nukleotid VAR, ORF cagrimi gerek
  C anotasyonlu ama RO adayi yok   52  %10,5
  D adayi var, ceviri tutmuyor    379  %76,3

D'nin icinde neyin oldugu, kurali gevsetme kararinin dayanagidir:
  baska protein (ortak parca yok) 114  %30,1   <- KABUL EDILSE YANLIS KOMSULUK
  ic parca ortak, uclar farkli    105  %27,7
  aday hedefin icinde              89  %23,5
  hedef adayin icinde              71  %18,7

Yani "ayni kayit" kuralina guvenilseydi en az 114 dizi YANLIS komsulukla
yayimlanacakti. Icerme oranina dayali bir ikinci kural (ornegin %90 ortusme)
265 cifti kurtarabilir ama bu AYRI bir kuraldir: ilan edilmeden, olculmeden
ve kendi kapisi kurulmadan uygulanmaz. Bu surumde kural birebir md5'tir.

Cikti: analysis_out/a5_match_guard.json
"""

import argparse
import collections
import csv
import hashlib
import json
import os
import re
import sqlite3

# Kalan eslesmeyen cift butcesi: olculen degerin kendisi. Butce olcumden
# genis olursa diziler sessizce kaybolur; dar olursa cikarim degisince patlar.
RESIDUAL_BUDGET = 497
DIFFERENT_PROTEIN = 114


def strip_accession(acc):
    """Surum ve NZ_ onekini atar -- atama ile ayni kural (bkz. a5_cohort.py)."""
    acc = (acc or "").strip()
    acc = re.sub(r"\.\d+$", "", acc)
    return acc[3:] if acc.startswith("NZ_") else acc


def seq_md5(seq):
    """build_union.py ile ayni anahtar."""
    return hashlib.md5(seq.strip().upper().encode()).hexdigest()


def read_candidates(paths):
    """replikon -> aday ceviri listesi."""
    seqs = collections.defaultdict(list)
    for path in paths:
        if not os.path.exists(path):
            continue
        header, buf = None, []
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if header is not None:
                        seqs[strip_accession(header.split(":")[0])].append("".join(buf).upper())
                    header, buf = line[1:].strip(), []
                else:
                    buf.append(line.strip())
            if header is not None:
                seqs[strip_accession(header.split(":")[0])].append("".join(buf).upper())
    return seqs


def read_replicons(paths):
    """replikon -> (status, cds_count)."""
    info = {}
    for path in paths:
        if not os.path.exists(path):
            continue
        for row in csv.DictReader(open(path)):
            key = strip_accession(row["nucleotide_id"])
            if key:
                info[key] = (row["status"], int(row["cds_count"] or 0))
    return info


def classify_residual(target, rep, cand_seqs, info, db_seq):
    """Eslesmeyen bir cifti nedenine gore siniflar."""
    if rep not in info:
        return "A_record_missing", None
    if info[rep][1] == 0:
        return "B_no_cds_annotation", None
    if not cand_seqs:
        return "C_no_ro_candidate", None
    hedef = (db_seq.get(target) or "").strip().upper()
    if not hedef:
        return "D_translation_differs", "sequence_not_in_db"
    for cand in cand_seqs:
        if hedef in cand:
            return "D_translation_differs", "target_within_candidate"
        if cand in hedef:
            return "D_translation_differs", "candidate_within_target"
        probe = hedef[20:60]
        if probe and probe in cand:
            return "D_translation_differs", "shared_core_different_ends"
    return "D_translation_differs", "different_protein"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--targets", default="cohorts/a5_protein_recovery.tsv")
    ap.add_argument("--context-dir", action="append", default=None,
                    help="cikarim dizini; birden fazla verilebilir")
    ap.add_argument("--out", default="analysis_out/a5_match_guard.json")
    args = ap.parse_args()

    dirs = args.context_dir or ["/media/lin-bio/back2/roar_gc_a5",
                                "/media/lin-bio/back2/roar_gc_wide"]
    cand = read_candidates([os.path.join(d, "ro_candidates.fasta") for d in dirs])
    info = read_replicons([os.path.join(d, "replicons.csv") for d in dirs])
    hashed = {k: {seq_md5(s) for s in v} for k, v in cand.items()}
    global_md5 = set().union(*hashed.values()) if hashed else set()

    con = sqlite3.connect(args.db)
    db_seq = dict(con.execute("SELECT seq_md5, sequence FROM ro_union"))
    con.close()

    targets = list(csv.DictReader(open(args.targets), delimiter="\t"))
    same = other = 0
    linked_seqs = set()
    all_seqs = set()
    causes = collections.Counter()
    d_detail = collections.Counter()
    for row in targets:
        md5v, rep = row["seq_md5"], strip_accession(row["accession"])
        all_seqs.add(md5v)
        if md5v in hashed.get(rep, ()):
            same += 1
            linked_seqs.add(md5v)
            continue
        if md5v in global_md5:
            other += 1
        cause, detail = classify_residual(md5v, rep, cand.get(rep, []), info, db_seq)
        causes[cause] += 1
        if detail:
            d_detail[detail] += 1

    residual = sum(causes.values())
    out = {
        "rule": "context is linked only when a candidate translation in the SAME "
                "record has a byte-identical md5",
        "n_target_pairs": len(targets),
        "n_target_sequences": len(all_seqs),
        "matched_same_record": same,
        "matched_other_record_only": other,
        "residual_pairs": residual,
        "residual_budget": RESIDUAL_BUDGET,
        "sequences_linked": len(linked_seqs),
        "link_rate": round(len(linked_seqs) / max(1, len(all_seqs)), 4),
        "residual_causes": dict(causes),
        "translation_mismatch_detail": dict(d_detail),
        # Kuralin kazanci: "ayni kayit yeter" denseydi yanlis komsulukla
        # yayimlanacak dizi sayisi.
        "would_be_wrong_without_rule": d_detail.get("different_protein", 0),
    }

    print(f"[kural] {out['rule']}")
    print(f"[hedef] {len(targets):,} cift / {len(all_seqs):,} tekil dizi")
    print(f"[eslesme] ayni kayitta birebir {same:,} (%{100*same/max(1,len(targets)):.1f})")
    print(f"[baglam] {len(linked_seqs):,}/{len(all_seqs):,} dizi (%{100*out['link_rate']:.1f})")
    print(f"[kalan] {residual:,} cift (butce {RESIDUAL_BUDGET})")
    for cause, n in sorted(causes.items()):
        print(f"   {cause:<24} {n:>5} (%{100*n/max(1,residual):.1f})")
    for detail, n in d_detail.most_common():
        print(f"      - {detail:<30} {n:>5}")
    print(f"[kural kazanci] birebir sart olmasa {out['would_be_wrong_without_rule']} "
          f"dizi YANLIS komsulukla yayimlanirdi")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[yazildi] {args.out}")


if __name__ == "__main__":
    main()
