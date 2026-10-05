"""
Genomik baglamdan cikan RO adaylarini dogrula ve kume ata.

extract_genomic_context.py ucuz bir regex on filtresiyle aday cikarir
(Rieske ligand imzasi + uzunluk). Bu adim adaylari asil testlerden gecirir:

    1. hmmsearch RieskeDB71.hmm  -> model coverage + kume/grup atamasi
    2. hmmalign  ROmotif71.hmm   -> Rieske ligandlari + katalitik triad testi
    3. ro tablosunu doldur

Iki bagimsiz olcut (profil kaplamasi ve dogrudan kalinti testi) kalibrasyonda
NEG setinde sadece %0.9 ayristi -- birbirlerini dogruluyorlar. is_confirmed
ikisinin de gecmesini sart kosar.
"""

import argparse
import os
import sqlite3
import subprocess
import sys

from ro_filter import parse_domtblout, classify_target
from ro_motif import read_stockholm_matchcols, classify_motifs


def run(command):
    print("[calisiyor]", " ".join(command))
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-dir", default="genomic_context")
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--cluster-hmm", default="RieskeDB71.hmm")
    parser.add_argument("--motif-hmm", default="ROs_71_Clean/ROmotif71.hmm")
    parser.add_argument("--cpu", type=int, default=8)
    args = parser.parse_args()

    fasta = os.path.join(args.context_dir, "ro_candidates.fasta")
    domtbl = os.path.join(args.context_dir, "cand_dom.out")
    stockholm = os.path.join(args.context_dir, "cand_aln.sto")

    if not os.path.exists(domtbl) or os.path.getsize(domtbl) == 0:
        run(["hmmsearch", "--domtblout", domtbl, "--noali", "--cpu", str(args.cpu),
             "-E", "1e-5", args.cluster_hmm, fasta])
    else:
        print(f"[atlandi] {domtbl} zaten var")

    if not os.path.exists(stockholm) or os.path.getsize(stockholm) == 0:
        run(["hmmalign", "--trim", "--amino", "-o", stockholm, args.motif_hmm, fasta])
    else:
        print(f"[atlandi] {stockholm} zaten var")

    print("[parse] domtblout")
    hits = parse_domtblout(domtbl)
    print("[parse] hizalama")
    aligned = read_stockholm_matchcols(stockholm)

    print("[birlestiriliyor]")
    updates = []
    for candidate_id, target_hits in hits.items():
        result, status = classify_target(target_hits)
        if result is None:
            continue
        motifs = classify_motifs(aligned[candidate_id]) if candidate_id in aligned else {}
        confirmed = int(status == "RO_alpha" and motifs.get("is_RO_alpha_motif", False))
        updates.append((
            result["PredictedCluster"], result["Group"], result["ModelCoverage"],
            result["Score"], result["E-value"],
            int(motifs.get("rieske_intact", False)),
            int(motifs.get("catalytic_intact", False)),
            confirmed, candidate_id,
        ))

    connection = sqlite3.connect(args.db)
    # hmmsearch'te hic hit almayan adaylar NULL kalmasin: 'elendi' acikca yazilsin,
    # yoksa is_confirmed=0 sorgulari onlari gormez.
    connection.execute("UPDATE ro SET ro_cluster='N/A', is_confirmed=0, "
                       "rieske_intact=0, catalytic_intact=0 WHERE is_confirmed IS NULL")
    connection.executemany("""
        UPDATE ro SET ro_cluster=?, ro_group=?, model_coverage=?, hmm_score=?,
                      hmm_evalue=?, rieske_intact=?, catalytic_intact=?, is_confirmed=?
        WHERE candidate_id=?""", updates)
    connection.commit()

    total = connection.execute("SELECT COUNT(*) FROM ro").fetchone()[0]
    confirmed = connection.execute(
        "SELECT COUNT(*) FROM ro WHERE is_confirmed=1").fetchone()[0]
    print(f"\n[sonuc] aday          : {total}")
    print(f"[sonuc] DOGRULANMIS RO : {confirmed}  ({100*confirmed/max(1,total):.1f}%)")

    print("\nGrup dagilimi (dogrulanmis):")
    for group, count in connection.execute(
            "SELECT ro_group, COUNT(*) FROM ro WHERE is_confirmed=1 "
            "GROUP BY ro_group ORDER BY ro_group"):
        print(f"   grup {group}: {count}")

    print("\nEleme sebepleri:")
    for label, query in [
            ("coverage gecti, motif kaldi",
             "SELECT COUNT(*) FROM ro WHERE ro_cluster!='N/A' AND is_confirmed=0"),
            ("Rieske var, katalitik yok",
             "SELECT COUNT(*) FROM ro WHERE rieske_intact=1 AND catalytic_intact=0"),
            ("coverage kaldi",
             "SELECT COUNT(*) FROM ro WHERE ro_cluster='N/A'")]:
        print(f"   {label:32s}: {connection.execute(query).fetchone()[0]}")
    connection.close()


if __name__ == "__main__":
    sys.exit(main())
