"""
Tum PF00355 seti uzerinde RO alpha filtresini calistir.

Kullanim:
    python run_ro_filter.py --fasta combined_pfam.fasta --hmm RieskeDB.hmm

Mevcut mongo_access.py protein basina bir hmmsearch calistiriyor (189.657
subprocess, her biri 73 HMM profilini bastan yukluyor). Birlesik fasta uzerinde
tek cagri ayni sonucu verir ve saatler yerine dakikalar surer.

Cikti:
    ro_domtbl.out       ham hmmsearch --domtblout ciktisi (tekrar kullanilabilir)
    ro_alpha.csv        filtreyi gecen RO alpha subunitleri
    ro_rejected.csv     elenenler + eleme sebebi (coverage ile birlikte)
    ro_alpha.fasta      gecenlerin sekanslari -- downstream analiz icin
"""

import argparse
import csv
import os
import subprocess
import sys

from ro_filter import (parse_domtblout, classify_target, DEFAULT_MIN_COVERAGE,
                       DEFAULT_MAX_EVALUE, DEFAULT_MIN_LENGTH)


def run_hmmsearch(hmm_file, fasta_file, domtblout, cpu, evalue):
    """Tek hmmsearch cagrisi. Cikti dosyasi varsa yeniden calistirmaz."""
    if os.path.exists(domtblout) and os.path.getsize(domtblout) > 0:
        print(f"[atlandi] {domtblout} zaten var, hmmsearch tekrar calistirilmadi.")
        return

    command = [
        "hmmsearch",
        "--domtblout", domtblout,
        "--cpu", str(cpu),
        "-E", str(evalue),
        "--noali",          # hizalama bloklarini yazma, sadece tablo lazim
        hmm_file,
        fasta_file,
    ]
    print("[calisiyor]", " ".join(command))
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
    print(f"[bitti] {domtblout}")


def load_sequences(fasta_file):
    """accession -> sequence. Baslik formati: >tr|ACC|NAME description"""
    sequences = {}
    accession, chunks = None, []
    with open(fasta_file) as handle:
        for line in handle:
            if line.startswith(">"):
                if accession:
                    sequences[accession] = "".join(chunks)
                header = line[1:].strip()
                accession = header.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if accession:
        sequences[accession] = "".join(chunks)
    return sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default="combined_pfam.fasta")
    parser.add_argument("--hmm", default="RieskeDB.hmm")
    parser.add_argument("--domtblout", default="ro_domtbl.out")
    parser.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE)
    parser.add_argument("--max-evalue", type=float, default=DEFAULT_MAX_EVALUE)
    parser.add_argument("--min-length", type=int, default=DEFAULT_MIN_LENGTH)
    parser.add_argument("--cpu", type=int, default=22)
    parser.add_argument("--search-evalue", type=float, default=1e-5,
                        help="hmmsearch'un raporlama esigi. Filtre esiginden GEVSEK "
                             "olmali ki elenenler de kayda gecsin.")
    parser.add_argument("--out-prefix", default="ro")
    args = parser.parse_args()

    run_hmmsearch(args.hmm, args.fasta, args.domtblout, args.cpu, args.search_evalue)

    print("[calisiyor] domtblout parse ediliyor...")
    hits = parse_domtblout(args.domtblout)
    print(f"[bilgi] hmmsearch'te hit alan protein sayisi: {len(hits)}")

    from collections import defaultdict
    by_status = defaultdict(list)
    for target_hits in hits.values():
        result, status = classify_target(target_hits, args.min_coverage,
                                         args.max_evalue, args.min_length)
        by_status[status].append(result)

    accepted = by_status.get("RO_alpha", [])
    print(f"[sonuc] RO_alpha  (kabul)        : {len(accepted)}")
    print(f"[sonuc] fragment  (eksik kayit)  : {len(by_status.get('fragment', []))}")
    print(f"[sonuc] not_alpha (katalitik yok): {len(by_status.get('not_alpha', []))}")

    outputs = {
        f"{args.out_prefix}_alpha.csv": accepted,
        f"{args.out_prefix}_fragment.csv": by_status.get("fragment", []),
        f"{args.out_prefix}_not_alpha.csv": by_status.get("not_alpha", []),
    }
    for name, rows in outputs.items():
        if not rows:
            continue
        with open(name, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"[yazildi] {name}  ({len(rows)} satir)")

    # Gecenlerin sekanslarini ayri bir fastaya yaz
    print("[calisiyor] sekanslar yaziliyor...")
    sequences = load_sequences(args.fasta)
    fasta_out = f"{args.out_prefix}_alpha.fasta"
    written, missing = 0, 0
    with open(fasta_out, "w") as handle:
        for row in accepted:
            sequence = sequences.get(row["Accession Name"])
            if sequence is None:
                missing += 1
                continue
            handle.write(f">{row['Accession Name']} cluster={row['PredictedCluster']} "
                         f"cov={row['ModelCoverage']}\n{sequence}\n")
            written += 1
    print(f"[yazildi] {fasta_out}  ({written} sekans"
          + (f", {missing} sekans bulunamadi)" if missing else ")"))

    # Kume dagilimi -- ozet
    from collections import Counter
    groups = Counter(row["Group"] for row in accepted)
    print("\nGrup dagilimi (kabul edilenler):")
    for group, count in sorted(groups.items()):
        print(f"   grup {group}: {count}")


if __name__ == "__main__":
    sys.exit(main())
