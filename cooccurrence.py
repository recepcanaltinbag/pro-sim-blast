"""
RO tiplerinin birlikte bulunmasi -- yol basamagi mi, gen duplikasyonu mu?

GOZLEM: operon dogrulamasinda (`operon_validation.py`) alfa alt biriminin
+-10 kb'inde 1.833 komsu genin KENDISI bir RO alfa alt birimi cikti. Bu iki
farkli sey olabilir:
    - ayni tipin ikinci kopyasi  -> gen duplikasyonu
    - farkli bir tip             -> ayni yolun iki basamagi ya da iki ayri yol

BU SCRIPT iki olcekte bakar:
    1. KOMSULUK: ayni +-10 kb icinde bulunan dogrulanmis RO ciftleri, ayni tip
       mi farkli tip mi; farkli tip ciftleri siklik sirasina gore.
    2. REPLIKON: ayni replikonda bulunan tip ciftleri, PERMUTASYON null'ina
       karsi. Null, her replikonun RO sayisini ve her tipin toplam uye sayisini
       KORUYARAK tipleri yeniden dagitir; boylece "iki tip de yaygin oldugu icin
       birlikte gorunuyor" aciklamasi elenir.

Null modeli neden boyle: tipler rastgele atanirsa kalabalik tipler her yerde
birlikte gorunur. Replikon basina sayi ve tip toplamlari sabit tutulunca test
"beklenenden FAZLA birlikte mi" sorusunu sorar.

Cikti: analysis_out/cooccurrence.json
"""

import argparse
import json
import os
import random
import sqlite3
from collections import Counter, defaultdict

NEIGHBOUR_WINDOW = 10000
PERMUTATIONS = 2000
RANDOM_SEED = 20261005
MIN_PAIR_COUNT = 5          # bu sayinin altindaki ciftler raporlanmaz


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/cooccurrence.json")
    ap.add_argument("--permutations", type=int, default=PERMUTATIONS)
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    rng = random.Random(RANDOM_SEED)

    ros = con.execute("""
        SELECT candidate_id, nucleotide_id, ro_cluster, start, end, strand
        FROM ro WHERE is_confirmed = 1""").fetchall()
    cluster_of = {r[0]: r[2] for r in ros}
    by_replicon = defaultdict(list)
    for cid, nuc, cluster, start, end, strand in ros:
        by_replicon[nuc].append((cid, cluster, start, end, strand))

    # --- 1. Komsuluk duzeyi: ayni pencerede iki dogrulanmis RO
    same_type, diff_type = 0, 0
    neighbour_pairs = Counter()
    same_strand_pairs = Counter()
    for nuc, items in by_replicon.items():
        items.sort(key=lambda x: x[2])
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                gap = items[j][2] - items[i][3]
                if gap > NEIGHBOUR_WINDOW:
                    break
                a, b = items[i], items[j]
                if a[1] == b[1]:
                    same_type += 1
                else:
                    diff_type += 1
                    neighbour_pairs[tuple(sorted((a[1], b[1])))] += 1
                    if a[4] == b[4]:
                        same_strand_pairs[tuple(sorted((a[1], b[1])))] += 1

    # --- 2. Replikon duzeyi: gozlenen tip ciftleri
    observed = Counter()
    replicon_types = {}
    for nuc, items in by_replicon.items():
        types = sorted({c for _, c, _, _, _ in items})
        replicon_types[nuc] = types
        for i in range(len(types)):
            for j in range(i + 1, len(types)):
                observed[(types[i], types[j])] += 1

    # --- permutasyon null'i: replikon basina RO sayisi ve tip toplamlari sabit
    counts_per_replicon = [len(v) for v in by_replicon.values()]
    pool = [c for _, c, _, _, _ in (x for v in by_replicon.values() for x in v)]
    null_counts = defaultdict(int)
    null_sq = defaultdict(int)
    for _ in range(args.permutations):
        shuffled = pool[:]
        rng.shuffle(shuffled)
        index = 0
        seen = Counter()
        for n in counts_per_replicon:
            types = sorted(set(shuffled[index:index + n]))
            index += n
            for i in range(len(types)):
                for j in range(i + 1, len(types)):
                    seen[(types[i], types[j])] += 1
        for pair, value in seen.items():
            null_counts[pair] += value
            null_sq[pair] += value * value

    pairs = []
    for pair, obs in observed.items():
        if obs < MIN_PAIR_COUNT:
            continue
        mean = null_counts[pair] / args.permutations
        var = max(0.0, null_sq[pair] / args.permutations - mean * mean)
        sd = var ** 0.5
        pairs.append({
            "type_a": pair[0], "type_b": pair[1], "observed": obs,
            "expected": round(mean, 2), "sd": round(sd, 2),
            "ratio": round(obs / mean, 2) if mean else None,
            "z": round((obs - mean) / sd, 2) if sd > 0 else None,
            "neighbours_within_10kb": neighbour_pairs.get(pair, 0),
        })
    pairs.sort(key=lambda p: -(p["z"] if p["z"] is not None else -99))

    out = {
        "neighbour_level": {
            "window_bp": NEIGHBOUR_WINDOW,
            "pairs_same_type": same_type,
            "pairs_different_type": diff_type,
            "top_different_type_pairs": [
                {"type_a": a, "type_b": b, "n": n,
                 "same_strand": same_strand_pairs.get((a, b), 0)}
                for (a, b), n in neighbour_pairs.most_common(15)],
        },
        "replicon_level": {
            "permutations": args.permutations,
            "replicons_with_an_ro": len(by_replicon),
            "replicons_with_more_than_one_type": sum(
                1 for t in replicon_types.values() if len(t) > 1),
            "pairs": pairs[:40],
        },
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"[komsuluk] ayni pencerede {same_type + diff_type:,} RO cifti: "
          f"{same_type:,} ayni tip (duplikasyon), {diff_type:,} farkli tip")
    print("[komsuluk] en sik farkli-tip ciftleri:")
    for p in out["neighbour_level"]["top_different_type_pairs"][:8]:
        print(f"   {p['type_a']:16s} + {p['type_b']:16s} {p['n']:4d}  "
              f"({p['same_strand']} ayni iplikte)")
    print(f"\n[replikon] {len(by_replicon):,} replikon, "
          f"{out['replicon_level']['replicons_with_more_than_one_type']:,} tanesinde birden fazla tip")
    print(f"[replikon] permutasyon null'i ({args.permutations} tekrar), en guclu birliktelikler:")
    for p in pairs[:10]:
        print(f"   {p['type_a']:16s} + {p['type_b']:16s} gozlenen {p['observed']:4d} "
              f"beklenen {p['expected']:7.2f}  {p['ratio']:>6}x  z={p['z']}")
    print(f"\n[yazildi] {args.out}")
    con.close()


if __name__ == "__main__":
    main()
