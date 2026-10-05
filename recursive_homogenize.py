"""
Ozyinelemeli homojenizasyon: her kumeyi homojen YAPRAKLARA bolunene kadar boler.

FIKIR (kullanicidan): "bir enzime baska atamalar olmussa ve homojen degilse,
onun alt dallarinda homojenize edene kadar gidebiliriz."

discover_subfamilies.py tek bir CD-HIT gecisi yapiyordu (%70). Ama bir alt-aile
hala heterojen olabilir. Burada islem OZYINELEMELI: bir dugum yeterince homojen
degilse (medyan ikili kimlik < hedef), daha yuksek bir esikte tekrar bolunur.
Sonuc bir AGAC -- her yaprak tek bir tutarli enzim tipi.

Neden onemli: substrat ataması ancak YAPRAK duzeyinde savunulabilir. "5_504_CntA
kumesine karnitin" demek %24,7 kimlikli 1210 protein icin gecersiz; ama o kumenin
bir homojen yapragina (ornegin %85 kimlikli 90 protein) atama anlamli.

ALGORITMA
  recurse(uyeler, derinlik):
    - cok kucukse veya medyan kimlik >= hedef  -> YAPRAK
    - degilse artan esikte CD-HIT ile bol, her cocuk icin recurse
    - esik bolemezse bir ust esige gec (derinlik+1)

Cikti:
    roar.sqlite   yeni tablolar: leaf, ro_leaf
    analysis_out/leaves.csv
"""

import argparse
import csv
import os
import random
import re
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

from ro_motif import read_stockholm_matchcols

# Bir dugumun "homojen yaprak" sayilmasi icin gereken medyan ikili kimlik
HOMOGENEITY_TARGET = 0.70
# Bundan kucuk dugumler bolunmez (asiri parcalanmayi onler)
MIN_SPLIT_SIZE = 8
# Artan CD-HIT esikleri -- her derinlikte bir sonraki kullanilir
THRESHOLDS = [0.50, 0.60, 0.70, 0.80, 0.90]
# Medyan kimlik hesabinda ornek buyuklugu
SAMPLE = 60
RANDOM_SEED = 20260722


def identity(a, b):
    match = total = 0
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue
        total += 1
        if x == y:
            match += 1
    return match / total if total else 0.0


def median_identity(sequences, rng):
    """Ornekten medyan ikili kimlik."""
    if len(sequences) < 2:
        return 1.0
    sample = (rng.sample(sequences, SAMPLE) if len(sequences) > SAMPLE else sequences)
    values = []
    for i in range(len(sample)):
        for j in range(i + 1, len(sample)):
            values.append(identity(sample[i], sample[j]))
    values.sort()
    return values[len(values) // 2] if values else 1.0


def cdhit_word(threshold):
    if threshold >= 0.70:
        return 5
    if threshold >= 0.60:
        return 4
    if threshold >= 0.50:
        return 3
    return 2


def cdhit(members, sequences, threshold):
    """CD-HIT ile bol. Donen: {candidate_id: alt_kume_index}."""
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = os.path.join(tmpdir, "in.fa")
        with open(fasta, "w") as handle:
            for candidate_id, sequence in zip(members, sequences):
                handle.write(">%s\n%s\n" % (candidate_id, sequence.replace("-", "")))
        out = os.path.join(tmpdir, "out")
        try:
            subprocess.run(["cd-hit", "-i", fasta, "-o", out, "-c", str(threshold),
                            "-n", str(cdhit_word(threshold)), "-M", "3000",
                            "-T", "4", "-d", "0"],
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        except (subprocess.CalledProcessError, FileNotFoundError, OSError):
            return {m: 0 for m in members}
        assignment, current = {}, None
        with open(out + ".clstr") as handle:
            for line in handle:
                if line.startswith(">Cluster"):
                    current = int(line.split()[1])
                else:
                    name = re.search(r">(.+?)\.\.\.", line)
                    if name:
                        assignment[name.group(1)] = current
    return assignment


def homogenize(members, aligned, rng):
    """Bir kumeyi homojen yapraklara boler. Donen: [(leaf_members, depth, median_id)]."""
    leaves = []
    # Yigin: (uyeler, derinlik)
    stack = [(members, 0)]
    while stack:
        node, depth = stack.pop()
        sequences = [aligned[m] for m in node]
        median = median_identity(sequences, rng)

        # Durma kosullari: kucuk, homojen, ya da esikler bitti
        if (len(node) < MIN_SPLIT_SIZE or median >= HOMOGENEITY_TARGET
                or depth >= len(THRESHOLDS)):
            leaves.append((node, depth, median))
            continue

        assignment = cdhit(node, sequences, THRESHOLDS[depth])
        children = defaultdict(list)
        for candidate_id in node:
            children[assignment.get(candidate_id, 0)].append(candidate_id)

        if len(children) <= 1:
            # Bu esik bolemedi -- bir ust esige gec, ayni uyelerle
            stack.append((node, depth + 1))
        else:
            for child in children.values():
                stack.append((child, depth + 1))
    return leaves


def ensure_schema(connection):
    connection.executescript("""
    DROP TABLE IF EXISTS leaf;
    DROP TABLE IF EXISTS ro_leaf;
    CREATE TABLE leaf (
        leaf_id        TEXT PRIMARY KEY,   -- "cluster#index"
        cluster        TEXT,
        size           INTEGER,
        depth          INTEGER,            -- kac bolme sonrasi homojenlesti
        median_identity REAL,
        is_homogeneous INTEGER,            -- medyan kimlik hedefi tutturdu mu
        representative TEXT,
        top_genera     TEXT
    );
    CREATE TABLE ro_leaf (
        candidate_id TEXT PRIMARY KEY,
        cluster      TEXT,
        leaf_id      TEXT,
        leaf_size    INTEGER,
        leaf_depth   INTEGER
    );
    CREATE INDEX idx_leaf_cluster ON leaf(cluster);
    CREATE INDEX idx_roleaf_leaf  ON ro_leaf(leaf_id);
    """)


def main():
    global HOMOGENEITY_TARGET
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    parser.add_argument("--out-dir", default="analysis_out")
    parser.add_argument("--target", type=float, default=HOMOGENEITY_TARGET)
    args = parser.parse_args()

    HOMOGENEITY_TARGET = args.target
    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(RANDOM_SEED)

    print("[okunuyor] hizalama")
    aligned = read_stockholm_matchcols(args.alignment)

    connection = sqlite3.connect(args.db)
    connection.row_factory = sqlite3.Row
    rows = connection.execute("""
        SELECT r.candidate_id, r.ro_cluster, rep.organism
        FROM ro r JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        WHERE r.is_confirmed = 1 AND r.ro_cluster IS NOT NULL AND r.ro_cluster != 'N/A'
    """).fetchall()
    by_cluster = defaultdict(list)
    organism_of = {}
    for row in rows:
        if row["candidate_id"] in aligned:
            by_cluster[row["ro_cluster"]].append(row["candidate_id"])
            organism_of[row["candidate_id"]] = row["organism"] or ""

    ensure_schema(connection)
    leaf_rows, member_rows = [], []
    clusters = sorted(by_cluster.items(), key=lambda x: -len(x[1]))

    for index, (cluster, members) in enumerate(clusters, 1):
        sys.stdout.write(f"\r  {index}/{len(clusters)} kume  ({cluster})      ")
        sys.stdout.flush()
        leaves = homogenize(members, aligned, rng)
        leaves.sort(key=lambda x: -len(x[0]))
        for leaf_index, (leaf_members, depth, median) in enumerate(leaves):
            leaf_id = f"{cluster}#{leaf_index}"
            genera = Counter(organism_of[m].split()[0] for m in leaf_members
                             if organism_of[m].split())
            leaf_rows.append((
                leaf_id, cluster, len(leaf_members), depth, round(median, 4),
                1 if median >= HOMOGENEITY_TARGET else 0, leaf_members[0],
                ";".join(f"{g}:{n}" for g, n in genera.most_common(5))))
            for candidate_id in leaf_members:
                member_rows.append((candidate_id, cluster, leaf_id,
                                    len(leaf_members), depth))
    print()

    connection.executemany("INSERT INTO leaf VALUES (?,?,?,?,?,?,?,?)", leaf_rows)
    connection.executemany("INSERT INTO ro_leaf VALUES (?,?,?,?,?)", member_rows)
    connection.commit()

    # --- Ozet
    total_leaves = len(leaf_rows)
    homogeneous = sum(1 for r in leaf_rows if r[5])
    big_leaves = [r for r in leaf_rows if r[2] >= 10]
    singletons = sum(1 for r in leaf_rows if r[2] == 1)
    covered = sum(r[2] for r in leaf_rows if r[2] >= 10)

    print("\n" + "=" * 70)
    print(f"OZYINELEMELI HOMOJENIZASYON (hedef medyan kimlik >= {HOMOGENEITY_TARGET:.0%})")
    print("=" * 70)
    print(f"  {len(clusters)} kume  ->  {total_leaves:,} yaprak")
    print(f"  homojen yaprak      : {homogeneous:,} ({100*homogeneous/total_leaves:.0f}%)")
    print(f"  >=10 uyeli yaprak   : {len(big_leaves):,}  ({covered:,} RO kapsiyor)")
    print(f"  tekil yaprak        : {singletons:,}")

    per_cluster = Counter(r[1] for r in leaf_rows)
    depth_of = defaultdict(int)
    for r in leaf_rows:
        depth_of[r[1]] = max(depth_of[r[1]], r[3])
    print("\n  En cok yaprak veren kumeler (heterojenlik = cok yaprak):")
    print(f"  {'kume':16s} {'uye':>6s} {'yaprak':>7s} {'max derinlik':>12s} "
          f"{'en buyuk yaprak':>15s}")
    print("  " + "-" * 62)
    size = {c: len(m) for c, m in by_cluster.items()}
    largest_leaf = defaultdict(int)
    for r in leaf_rows:
        largest_leaf[r[1]] = max(largest_leaf[r[1]], r[2])
    for cluster in sorted(per_cluster, key=lambda c: -per_cluster[c])[:14]:
        print(f"  {cluster:16s} {size[cluster]:>6,} {per_cluster[cluster]:>7d} "
              f"{depth_of[cluster]:>12d} {largest_leaf[cluster]:>15,}")

    path = os.path.join(args.out_dir, "leaves.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["leaf_id", "cluster", "size", "depth", "median_identity",
                         "is_homogeneous", "representative", "top_genera"])
        writer.writerows(leaf_rows)
    print(f"\n[yazildi] {path}  ({total_leaves} yaprak)")

    print("\n  En homojen buyuk yapraklar (temiz enzim tipleri):")
    for r in sorted(big_leaves, key=lambda x: -x[4])[:10]:
        genera = r[7].split(";")[0] if r[7] else "?"
        print(f"    {r[0]:20s} n={r[2]:>4d}  kimlik={r[4]:.2f}  {genera}")
    connection.close()


if __name__ == "__main__":
    sys.exit(main() or 0)
