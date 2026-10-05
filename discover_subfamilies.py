"""
Alt-aile kesfi + mutlak atama olcutu + novel aday siralamasi.

Uc geliştirmeyi birlestirir:

  (1) MUTLAK OLCUT
      Kume atamasi goreli en iyi modeli secer -- bir uyenin atandigi kumeye
      GERCEKTEN benzeyip benzemedigini soylemez. Burada model-uzunlugundan
      bagimsiz mutlak bir olcu kullanilir: uyenin, kendi kumesinin en kalabalik
      alt-ailesinin KONSENSUSUNA kimligi. Bu, "goreli kazanan" degil, dogrudan
      "ne kadar tipik bir uye" sorusunu cevaplar.

  (2) ALT-AILE BOLME
      Her kume CD-HIT ile %70 kimlikte alt-ailelere ayrilir. Heterojen bir kume
      (ornegin CntA, medyan kimlik %24,7) onlarca alt-aileye dagilir. En kalabalik
      alt-aile "cekirdek", geri kalanlar cevre.

  (3) NOVEL ADAY
      Novel RO tam olarak cekirdekten uzakta, kendi kucuk alt-ailesinde oturan
      uyedir. Bunlar su an "atanmis" gorundugu icin gorunmez; cekirdek kimligine
      gore siralaninca ortaya cikarlar. Gercek yeni aileler burada.

Cikti:
    roar.sqlite     yeni tablolar: subfamily, ro_subfamily
    analysis_out/subfamilies.csv
    analysis_out/novel_candidates.csv
"""

import argparse
import csv
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

from ro_motif import read_stockholm_matchcols

# CD-HIT alt-aile esigi ve cekirdek/cevre/novel siniri
SUBFAMILY_IDENTITY = 0.70
CORE_IDENTITY = 0.55        # bir REFERANS tipe bu ustu kimlik = "cekirdek"
NOVEL_IDENTITY = 0.40       # hicbir referans tipe bu kadar bile benzemiyor

MIN_CLUSTER_N = 10          # bundan kucuk kumelerde alt-aile analizi yapilmaz
MIN_SUBFAMILY_FOR_CONSENSUS = 3

# Bir alt-aile "referans tip" sayilmasi icin gereken minimum uye.
# NEDEN: heterojen bir kumede en kalabalik alt-aile bile %7-14 olabilir.
# Tek bir keyfi cekirdege kimlik olcmek novel'i sisirir -- 870 kisilik tutarli
# bir alt-aile "kumeye benzemiyor" diye novel sayilir, oysa o kumenin ikinci
# ana enzim tipidir (bolme sorunu, novellik degil). Bu yuzden uyeyi TUM buyuk
# alt-ailelere kiyaslar, en iyi uyumu aliriz. Novel = hicbir yerlesik tipe
# uymuyor VE kendisi de kucuk/tekil bir alt-ailede.
REFERENCE_TYPE_MIN_SIZE = 10


def consensus(sequences, length):
    """Hizalanmis dizilerden kolon-bazli konsensus (bosluk olmayan en sik kalinti)."""
    result = []
    for position in range(length):
        counter = Counter(s[position] for s in sequences if s[position] != "-")
        result.append(counter.most_common(1)[0][0] if counter else "-")
    return "".join(result)


def identity_to(sequence, reference):
    """Bir uyenin bir referans konsensusa kimligi (uyede bosluk olmayan kolonlar)."""
    match = total = 0
    for a, b in zip(sequence, reference):
        if a == "-":
            continue
        total += 1
        if a == b:
            match += 1
    return match / total if total else 0.0


def cdhit_subfamilies(members, sequences):
    """CD-HIT ile alt-aileler. Donen: {candidate_id: subfamily_index}, temsilciler.

    members: [candidate_id, ...] ; sequences: hizalanmis diziler (ayni sirada)
    """
    if len(members) < 2:
        return {members[0]: 0}, {0: members[0]}
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = os.path.join(tmpdir, "in.fa")
        with open(fasta, "w") as handle:
            for candidate_id, sequence in zip(members, sequences):
                handle.write(">%s\n%s\n" % (candidate_id, sequence.replace("-", "")))
        out = os.path.join(tmpdir, "out")
        try:
            subprocess.run(["cd-hit", "-i", fasta, "-o", out,
                            "-c", str(SUBFAMILY_IDENTITY), "-n", "5",
                            "-M", "3000", "-T", "4", "-d", "0"],
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        except (subprocess.CalledProcessError, FileNotFoundError, OSError):
            return {m: 0 for m in members}, {0: members[0]}

        assignment, representatives = {}, {}
        current = None
        with open(out + ".clstr") as handle:
            for line in handle:
                if line.startswith(">Cluster"):
                    current = int(line.split()[1])
                else:
                    name = re.search(r">(.+?)\.\.\.", line)
                    if not name:
                        continue
                    assignment[name.group(1)] = current
                    if line.rstrip().endswith("*"):
                        representatives[current] = name.group(1)
    return assignment, representatives


def ensure_schema(connection):
    connection.executescript("""
    DROP TABLE IF EXISTS subfamily;
    DROP TABLE IF EXISTS ro_subfamily;
    CREATE TABLE subfamily (
        cluster        TEXT,
        subfamily      INTEGER,          -- 0 = en kalabalik (cekirdek)
        subfamily_id   TEXT PRIMARY KEY, -- "cluster.subfamily"
        size           INTEGER,
        is_core        INTEGER,          -- kumenin cekirdek alt-ailesi mi
        representative TEXT,
        top_genera     TEXT,
        consensus_len  INTEGER
    );
    CREATE TABLE ro_subfamily (
        candidate_id       TEXT PRIMARY KEY,
        cluster            TEXT,
        subfamily          INTEGER,
        subfamily_id       TEXT,
        subfamily_size     INTEGER,
        core_identity      REAL,   -- en iyi uyan REFERANS tipe kimlik (MUTLAK olcut)
        assignment_class   TEXT    -- core | divergent | alt_type | novel_candidate
    );
    CREATE INDEX idx_rosub_cluster ON ro_subfamily(cluster);
    CREATE INDEX idx_rosub_class   ON ro_subfamily(assignment_class);
    """)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    parser.add_argument("--out-dir", default="analysis_out")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print("[okunuyor] hizalama")
    aligned = read_stockholm_matchcols(args.alignment)
    length = len(next(iter(aligned.values())))
    print(f"  {len(aligned)} dizi, {length} kolon")

    connection = sqlite3.connect(args.db)
    connection.row_factory = sqlite3.Row
    rows = connection.execute("""
        SELECT r.candidate_id, r.ro_cluster, rep.organism
        FROM ro r JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        WHERE r.is_confirmed = 1 AND r.ro_cluster IS NOT NULL
              AND r.ro_cluster != 'N/A'
    """).fetchall()

    by_cluster = defaultdict(list)
    organism_of = {}
    for row in rows:
        sequence = aligned.get(row["candidate_id"])
        if sequence:
            by_cluster[row["ro_cluster"]].append(row["candidate_id"])
            organism_of[row["candidate_id"]] = row["organism"] or ""

    ensure_schema(connection)

    subfamily_rows, member_rows = [], []
    class_totals = Counter()

    clusters = sorted(by_cluster.items(), key=lambda x: -len(x[1]))
    for index, (cluster, members) in enumerate(clusters, 1):
        sys.stdout.write(f"\r  {index}/{len(clusters)} kume isleniyor")
        sys.stdout.flush()

        sequences = [aligned[m] for m in members]

        def classify(best_fit, own_size):
            """MUTLAK olcut: en iyi uyan referans tipe kimlik + kendi alt-aile boyutu."""
            if best_fit >= CORE_IDENTITY:
                return "core"
            if best_fit >= NOVEL_IDENTITY:
                return "divergent"
            # Hicbir yerlesik tipe uymuyor. Kendi alt-ailesi buyukse bu bir
            # BOLME sorunu (tutarli ikinci tip); kucukse gercek novel aday.
            return "alt_type" if own_size >= REFERENCE_TYPE_MIN_SIZE else "novel_candidate"

        if len(members) < MIN_CLUSTER_N:
            core = consensus(sequences, length)
            for candidate_id, sequence in zip(members, sequences):
                ci = identity_to(sequence, core)
                cls = classify(ci, len(members))
                class_totals[cls] += 1
                member_rows.append((candidate_id, cluster, 0, f"{cluster}.0",
                                    len(members), round(ci, 4), cls))
            genera = Counter(organism_of[m].split()[0] for m in members
                             if organism_of[m].split())
            subfamily_rows.append((cluster, 0, f"{cluster}.0", len(members), 1,
                                   members[0], ";".join(f"{g}:{n}" for g, n
                                   in genera.most_common(5)), length))
            continue

        # CD-HIT alt-aileleri
        assignment, representatives = cdhit_subfamilies(members, sequences)
        subfam_members = defaultdict(list)
        for candidate_id in members:
            subfam_members[assignment.get(candidate_id, 0)].append(candidate_id)

        # Boyuta gore yeniden numaralandir: 0 = en kalabalik
        ordered = sorted(subfam_members.items(), key=lambda x: -len(x[1]))

        # REFERANS TIPLER: yeterince kalabalik alt-aileler. Hicbiri esigi
        # gecmiyorsa (kume tamamen dagilmis) en kalabalik alt-aileyi al.
        reference_consensuses = []
        for _, sub_members in ordered:
            if len(sub_members) >= REFERENCE_TYPE_MIN_SIZE:
                reference_consensuses.append(
                    consensus([aligned[m] for m in sub_members], length))
        if not reference_consensuses:
            reference_consensuses.append(
                consensus([aligned[m] for m in ordered[0][1]], length))

        for new_index, (old_index, sub_members) in enumerate(ordered):
            genera = Counter(organism_of[m].split()[0] for m in sub_members
                             if organism_of[m].split())
            subfamily_rows.append((
                cluster, new_index, f"{cluster}.{new_index}", len(sub_members),
                1 if new_index == 0 else 0,
                representatives.get(old_index, sub_members[0]),
                ";".join(f"{g}:{n}" for g, n in genera.most_common(5)), length))

            for candidate_id in sub_members:
                # En iyi uyan referans tipe kimlik -- tek keyfi cekirdek yerine
                best_fit = max(identity_to(aligned[candidate_id], ref)
                               for ref in reference_consensuses)
                cls = classify(best_fit, len(sub_members))
                class_totals[cls] += 1
                member_rows.append((candidate_id, cluster, new_index,
                                    f"{cluster}.{new_index}", len(sub_members),
                                    round(best_fit, 4), cls))
    print()

    connection.executemany(
        "INSERT INTO subfamily VALUES (?,?,?,?,?,?,?,?)", subfamily_rows)
    connection.executemany(
        "INSERT INTO ro_subfamily VALUES (?,?,?,?,?,?,?)", member_rows)
    connection.commit()

    # --- Ozet
    print("\n" + "=" * 68)
    print("MUTLAK ATAMA OLCUTU (cekirdek alt-aile konsensusuna kimlik)")
    print("=" * 68)
    total = sum(class_totals.values())
    labels = {
        "core": "cekirdek -- yerlesik bir tipe uyuyor",
        "divergent": "cevre -- uyuyor ama uzak",
        "alt_type": "alternatif tip -- tutarli, kume bolunmeli",
        "novel_candidate": "NOVEL -- hicbir tipe uymuyor, izole",
    }
    for cls in ("core", "divergent", "alt_type", "novel_candidate"):
        n = class_totals[cls]
        print(f"  {cls:16s}: {n:>7,}  ({100*n/total:4.1f}%)  {labels[cls]}")
    print(f"\n  novel_candidate = hicbir referans tipe kimlik < {NOVEL_IDENTITY:.0%}"
          f" VE kucuk alt-ailede.\n  alt_type = ayni durum ama BUYUK alt-ailede"
          f" -- gercek yeni enzim tipi, novellik degil bolme.")

    # --- Alt-aile sayilari
    per_cluster = Counter(row[0] for row in subfamily_rows)
    print("\n" + "=" * 68)
    print("EN COK ALT-AILEYE BOLUNEN KUMELER")
    print("=" * 68)
    print(f"  {'kume':18s} {'uye':>6s} {'alt-aile':>9s} {'cekirdek %':>11s} {'novel':>7s}")
    print("  " + "-" * 56)
    cluster_size = {c: len(m) for c, m in by_cluster.items()}
    novel_by_cluster = Counter()
    core_by_cluster = Counter()
    for row in member_rows:
        if row[6] == "novel_candidate":
            novel_by_cluster[row[1]] += 1
        if row[6] == "core":
            core_by_cluster[row[1]] += 1
    for cluster in sorted(per_cluster, key=lambda c: -per_cluster[c])[:15]:
        n = cluster_size[cluster]
        core_pct = core_by_cluster[cluster] / n if n else 0
        print(f"  {cluster:18s} {n:>6,} {per_cluster[cluster]:>9d} "
              f"{core_pct:>10.1%} {novel_by_cluster[cluster]:>7d}")

    # --- CSV ciktilari
    sub_path = os.path.join(args.out_dir, "subfamilies.csv")
    with open(sub_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cluster", "subfamily", "subfamily_id", "size", "is_core",
                         "representative", "top_genera", "consensus_len"])
        writer.writerows(subfamily_rows)

    # Novel adaylar: cekirdek kimligi en dusuk, kucuk alt-ailedekiler
    novel = connection.execute("""
        SELECT s.candidate_id, s.cluster, s.subfamily_id, s.subfamily_size,
               s.core_identity, rep.organism, r.hmm_score, r.model_coverage
        FROM ro_subfamily s
        JOIN ro r ON r.candidate_id = s.candidate_id
        JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        WHERE s.assignment_class = 'novel_candidate'
        ORDER BY s.core_identity ASC
    """).fetchall()
    novel_path = os.path.join(args.out_dir, "novel_candidates.csv")
    with open(novel_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["candidate_id", "assigned_cluster", "subfamily_id",
                         "subfamily_size", "core_identity", "organism",
                         "hmm_score", "model_coverage"])
        for row in novel:
            writer.writerow([row["candidate_id"], row["cluster"], row["subfamily_id"],
                             row["subfamily_size"], row["core_identity"],
                             row["organism"], row["hmm_score"], row["model_coverage"]])

    print(f"\n[yazildi] {sub_path}  ({len(subfamily_rows)} alt-aile)")
    print(f"[yazildi] {novel_path}  ({len(novel)} novel aday)")
    print("\nEn dusuk cekirdek kimlikli 12 novel aday:")
    for row in novel[:12]:
        organism = (row["organism"] or "?")[:32]
        print(f"  {row['core_identity']:.2f}  {row['cluster']:14s} "
              f"altaile={row['subfamily_size']:>3d}kisi  {organism}")
    connection.close()


if __name__ == "__main__":
    sys.exit(main() or 0)
