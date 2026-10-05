"""
Kume ICI varyant analizi: "ayni enzim" diye siniflandirilanlar gercekten ayni mi?

Kume atamasi 71 model arasindan en yuksek bit skorunu alani secer. Ama bu
atama bir kumenin HOMOJEN oldugunu garanti etmez -- 3_313_KshA15 kumesinde
1.345 protein var ve 100'den fazla cinse yayilmis. Bunlarin tek bir enzim
olmasi beklenmez.

UC OLCUM

1. KUME ICI KIMLIK
   Tum uyeler ayni 426 kolonluk modele hizalanmis durumda, dolayisiyla ikili
   kimlik dogrudan match kolonlarindan okunabilir. Dusuk medyan kimlik =
   kume aslinda birden fazla enzimi barindiriyor.

2. ALT-AILE SAYISI (CD-HIT)
   Her kume %70 kimlikte kumelenir. Tek bir enzime karsilik gelen bir kume
   birkac alt kumeye ayrilir; heterojen bir kume onlarcaya.

3. SPESIFISITE BELIRLEYICI POZISYONLAR (SDP)
   Asil ilginc olan bu. Substrat spesifikligi katalitik cebi doseyen birkac
   kalinti tarafindan belirlenir (klasik ornek: BphA'da 335/336/376/377
   pozisyonlari PCB konjeneri secicigini belirler).

   SDP imzasi: kolon kume ICINDE korunmus ama kumeler ARASINDA farkli.
   Sadece korunmus kolonlar (her yerde ayni) yapisal; sadece degisken kolonlar
   gurultu. Ikisinin arasindaki kolonlar islevi tasir.

       SDP skoru = (kumeler arasi cesitlilik) - (kume ici ortalama cesitlilik)

   Yuksek skor = pozisyon her kumede kendi kalintisini korumus, ama kumeden
   kumeye degismis.
"""

import argparse
import csv
import math
import os
import random
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

from ro_motif import read_stockholm_matchcols, MODEL_COLUMNS, DEFAULT_MODEL

# Katalitik domain kabaca burada baslar (Rieske N-terminal, katalitik C-terminal)
CATALYTIC_DOMAIN_START = 180

# Ikili kimlik hesabinda kume basina en fazla kac uye orneklenecek
MAX_PAIRWISE_SAMPLE = 120

# CD-HIT alt-aile esigi
SUBFAMILY_IDENTITY = 0.70

RANDOM_SEED = 20260721


def identity(seq_a, seq_b):
    """Iki hizalanmis dizi arasinda kimlik (her ikisinde de bosluk olmayan kolonlar)."""
    match = total = 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        total += 1
        if a == b:
            match += 1
    return match / total if total else 0.0


def entropy(counter):
    total = sum(counter.values())
    if total <= 0:
        return 0.0
    return -sum((n / total) * math.log(n / total, 20) for n in counter.values() if n)


def run_cdhit(sequences, identity_threshold):
    """CD-HIT ile alt-aile sayisini bul. Donen: kume sayisi (veya None)."""
    if len(sequences) < 2:
        return 1
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = os.path.join(tmpdir, "in.fa")
        with open(fasta, "w") as handle:
            for index, sequence in enumerate(sequences):
                handle.write(">s%d\n%s\n" % (index, sequence.replace("-", "")))
        out = os.path.join(tmpdir, "out")
        word = 5 if identity_threshold >= 0.7 else 4
        try:
            subprocess.run(["cd-hit", "-i", fasta, "-o", out,
                            "-c", str(identity_threshold), "-n", str(word),
                            "-M", "2000", "-T", "4", "-d", "0"],
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            with open(out) as handle:
                return sum(1 for line in handle if line.startswith(">"))
        except (subprocess.CalledProcessError, FileNotFoundError, OSError):
            return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    parser.add_argument("--min-n", type=int, default=25)
    parser.add_argument("--out-dir", default="analysis_out")
    parser.add_argument("--skip-cdhit", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(RANDOM_SEED)

    print("[okunuyor] hizalama")
    aligned = read_stockholm_matchcols(args.alignment)
    print(f"  {len(aligned)} hizalanmis dizi")

    connection = sqlite3.connect(args.db)
    assignment = {row[0]: row[1] for row in connection.execute(
        "SELECT candidate_id, ro_cluster FROM ro WHERE is_confirmed=1 "
        "AND ro_cluster IS NOT NULL AND ro_cluster != 'N/A'")}
    connection.close()

    by_cluster = defaultdict(list)
    for candidate_id, cluster in assignment.items():
        sequence = aligned.get(candidate_id)
        if sequence:
            by_cluster[cluster].append(sequence)
    print(f"  {len(by_cluster)} kume, {sum(len(v) for v in by_cluster.values())} dizi")

    clusters = {k: v for k, v in by_cluster.items() if len(v) >= args.min_n}
    print(f"  n>={args.min_n} olan {len(clusters)} kume analiz edilecek")
    if not clusters:
        return 1

    length = len(next(iter(aligned.values())))

    # ---------------- 1 & 2: kume ici kimlik + alt-aile
    print("\n[hesaplaniyor] kume ici cesitlilik")
    rows = []
    for index, (cluster, sequences) in enumerate(
            sorted(clusters.items(), key=lambda x: -len(x[1])), 1):
        sample = (rng.sample(sequences, MAX_PAIRWISE_SAMPLE)
                  if len(sequences) > MAX_PAIRWISE_SAMPLE else sequences)
        identities = []
        for i in range(len(sample)):
            for j in range(i + 1, len(sample)):
                identities.append(identity(sample[i], sample[j]))
        identities.sort()
        median = identities[len(identities) // 2] if identities else 1.0
        low = identities[len(identities) // 20] if len(identities) >= 20 else median

        subfamilies = None if args.skip_cdhit else run_cdhit(sample, SUBFAMILY_IDENTITY)

        rows.append({
            "cluster": cluster, "n": len(sequences), "sampled": len(sample),
            "median_identity": round(median, 4), "p5_identity": round(low, 4),
            "subfamilies": subfamilies,
            "subfam_per_100": round(100 * subfamilies / len(sample), 1)
                              if subfamilies else None,
        })
        sys.stdout.write(f"\r  {index}/{len(clusters)} kume")
        sys.stdout.flush()
    print()

    rows.sort(key=lambda r: r["median_identity"])
    path = os.path.join(args.out_dir, "cluster_variance.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print("\n" + "=" * 72)
    print("KUME ICI CESITLILIK  (dusuk kimlik = kume aslinda heterojen)")
    print("=" * 72)
    print(f"{'kume':18s} {'n':>6s} {'medyan kimlik':>14s} {'%5 dilim':>10s} "
          f"{'alt-aile':>9s}")
    print("-" * 62)
    for row in rows[:16]:
        sub = str(row["subfamilies"]) if row["subfamilies"] else "-"
        print(f"{row['cluster']:18s} {row['n']:>6,} {row['median_identity']:>13.1%} "
              f"{row['p5_identity']:>9.1%} {sub:>9s}")
    print("  ...")
    for row in rows[-5:]:
        sub = str(row["subfamilies"]) if row["subfamilies"] else "-"
        print(f"{row['cluster']:18s} {row['n']:>6,} {row['median_identity']:>13.1%} "
              f"{row['p5_identity']:>9.1%} {sub:>9s}")
    print(f"\n[yazildi] {path}")

    # ---------------- 3: spesifisite belirleyici pozisyonlar
    print("\n[hesaplaniyor] spesifisite belirleyici pozisyonlar")
    within = [0.0] * length
    between = [0.0] * length
    coverage = [0] * length

    per_cluster_columns = {}
    for cluster, sequences in clusters.items():
        columns = []
        for position in range(length):
            counter = Counter(s[position] for s in sequences if s[position] != "-")
            columns.append(counter)
        per_cluster_columns[cluster] = columns

    total_clusters = len(clusters)
    for position in range(length):
        entropies, consensus = [], Counter()
        filled = 0
        for cluster, columns in per_cluster_columns.items():
            counter = columns[position]
            if sum(counter.values()) < 5:
                continue
            filled += 1
            entropies.append(entropy(counter))
            consensus[counter.most_common(1)[0][0]] += 1
        if filled < total_clusters * 0.5:
            continue
        within[position] = sum(entropies) / len(entropies)
        between[position] = entropy(consensus)
        coverage[position] = filled

    scored = [(between[i] - within[i], i + 1, within[i], between[i])
              for i in range(length) if coverage[i]]
    scored.sort(reverse=True)

    motif_positions = {c for c, _, _ in MODEL_COLUMNS[DEFAULT_MODEL]["rieske"]}
    motif_positions |= {c for c, _, _ in MODEL_COLUMNS[DEFAULT_MODEL]["catalytic"]}

    path = os.path.join(args.out_dir, "sdp_positions.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["column", "sdp_score", "within_cluster_entropy",
                         "between_cluster_entropy", "in_catalytic_domain",
                         "is_motif_position", "top_residues_by_cluster"])
        for score, column, w, b in scored:
            top = ";".join(
                f"{cluster}:{per_cluster_columns[cluster][column-1].most_common(1)[0][0]}"
                for cluster in list(clusters)[:8]
                if sum(per_cluster_columns[cluster][column-1].values()) >= 5)
            writer.writerow([column, round(score, 4), round(w, 4), round(b, 4),
                             int(column >= CATALYTIC_DOMAIN_START),
                             int(column in motif_positions), top])

    print("\n" + "=" * 72)
    print("SPESIFISITE BELIRLEYICI POZISYONLAR (SDP)")
    print("=" * 72)
    print("kume ICINDE korunmus, kumeler ARASINDA farkli kolonlar --")
    print("substrat secicigini tasimasi beklenen pozisyonlar.\n")
    print(f"{'kolon':>7s} {'SDP skoru':>10s} {'kume ici':>9s} {'kumeler arasi':>14s}  bolge")
    print("-" * 60)
    for score, column, w, b in scored[:20]:
        region = "KATALITIK" if column >= CATALYTIC_DOMAIN_START else "Rieske/N-term"
        flag = "  <- motif" if column in motif_positions else ""
        print(f"{column:>7d} {score:>10.3f} {w:>9.3f} {b:>14.3f}  {region}{flag}")

    catalytic_top = sum(1 for s, c, _, _ in scored[:30] if c >= CATALYTIC_DOMAIN_START)
    print(f"\nIlk 30 SDP'nin {catalytic_top}'i katalitik domainde "
          f"(kolon >= {CATALYTIC_DOMAIN_START}).")
    print(f"[yazildi] {path}")


if __name__ == "__main__":
    sys.exit(main() or 0)
