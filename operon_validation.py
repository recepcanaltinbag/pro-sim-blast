"""
Operon tanimi gercekten bir sey olcuyor mu? -- esigin ve kuralin sinanmasi.

SORUN: `build_operons.py` operonu "ayni iplik, genler arasi bosluk <=150 bp"
diye tanimlar. Bu bir KONVANSIYONDUR; transkripsiyon verisi olmadan dogrudan
dogrulanamaz. Dogrulanabilecek olan sudur: eger kural anlamliysa, dizi ile
dogrulanmis ortaklar (beta, ferredoksin, reduktaz) alfa alt biriminin yanina
RASTGELE BIR GENDEN DAHA SIK ve DAHA YAKIN yerlesmelidir.

BU SCRIPT ucunu olcer:
  1. Yon:    ortaklarin ayni iplikte olma orani, ayni penceredeki TUM genlerin
             orani ile karsilastirilir (Fisher).
  2. Konum:  ortaklarin |gen uzakligi| <= 2 olma orani, yine tum genlere karsi.
  3. Esik:   bosluk esigi 50-500 bp arasinda degistirilince operonda ortak
             bulunan RO sayisi nasil degisiyor -- 150 bp secimi sonucu ne kadar
             belirliyor.

Arka plan kumesi ayni pencerelerin kendisidir (her dogrulanmis RO'nun +-10 kb'i),
yani karsilastirma ayni genomlarda, ayni bolgelerde yapilir; genom bollugu veya
anotasyon yanliligi iki tarafi da ayni sekilde etkiler.

Cikti: analysis_out/operon_validation.json
"""

import argparse
import json
import os
import sqlite3
from collections import Counter, defaultdict

from scipy import stats

COMPONENTS = ["beta", "ferredoxin", "reductase", "alpha_other"]
GAP_THRESHOLDS = [50, 100, 150, 200, 300, 500]
NEAR = 2            # "yakin" tanimi: |gen uzakligi| <= 2
OFFSET_RANGE = 10


def fisher(a_yes, a_no, b_yes, b_no):
    odds, p = stats.fisher_exact([[a_yes, a_no], [b_yes, b_no]])
    ra = a_yes / max(1, a_yes + a_no)
    rb = b_yes / max(1, b_yes + b_no)
    return {"rate": round(ra, 4), "background_rate": round(rb, 4),
            "ratio": round(ra / rb, 3) if rb else None,
            "odds_ratio": round(float(odds), 3), "p": float(p),
            "n": a_yes + a_no, "n_background": b_yes + b_no}


def walk_counts(ro_rows, neighbors, components, max_gap):
    """build_operons.py ile ayni kurali uygulayip her esikte ortak sayisini verir."""
    found = Counter()
    for cid, strand, start, end in ro_rows:
        nbs = neighbors.get(cid, ())
        left = sorted((n for n in nbs if n[6] < 0), key=lambda n: -n[6])
        right = sorted((n for n in nbs if n[6] > 0), key=lambda n: n[6])
        seen = set()
        for side, direction in ((left, -1), (right, 1)):
            prev_start, prev_end = start, end
            for n in side:
                n_start, n_end, n_strand, key, spans = n[0], n[1], n[2], n[3], n[5]
                if n_strand != strand:
                    break
                gap = 0 if spans else (prev_start - n_end if direction < 0 else n_start - prev_end)
                if gap > max_gap:
                    break
                comp = components.get(key)
                if comp:
                    seen.add(comp)
                prev_start, prev_end = n_start, n_end
        for comp in seen:
            found[comp] += 1
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/operon_validation.json")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    components = dict(con.execute(
        "SELECT protein_key, component FROM neighbor_component WHERE component != 'none'"))

    rows = con.execute("""
        SELECT nb.candidate_id, nb.start, nb.end, nb.strand, nb.same_strand, nb.gene_offset,
               nb.distance, nb.spans_origin,
               nb.nucleotide_id||':'||nb.start||'-'||nb.end||':'||nb.strand AS pkey
        FROM neighbor nb JOIN ro r ON r.candidate_id = nb.candidate_id
        WHERE r.is_confirmed = 1""").fetchall()

    bg_same = Counter()
    bg_near = Counter()
    per_component = defaultdict(lambda: {"same": Counter(), "near": Counter(),
                                         "offsets": Counter(), "gaps": []})
    neighbors = defaultdict(list)
    for cid, start, end, strand, same, offset, distance, spans, pkey in rows:
        bg_same[bool(same)] += 1
        bg_near[abs(offset) <= NEAR] += 1
        neighbors[cid].append((start, end, strand, pkey, distance, spans, offset))
        comp = components.get(pkey)
        if comp in COMPONENTS:
            d = per_component[comp]
            d["same"][bool(same)] += 1
            d["near"][abs(offset) <= NEAR] += 1
            if abs(offset) <= OFFSET_RANGE:
                d["offsets"][offset] += 1
            if same and distance is not None:
                d["gaps"].append(abs(distance))

    out = {"n_neighbour_genes": len(rows),
           "background": {"same_strand_rate": round(bg_same[True] / max(1, len(rows)), 4),
                          "near_rate": round(bg_near[True] / max(1, len(rows)), 4),
                          "near_definition": f"|gene offset| <= {NEAR}"},
           "components": {}}

    print(f"[arka plan] {len(rows):,} komsu gen; ayni iplik "
          f"%{100 * out['background']['same_strand_rate']:.1f}, "
          f"yakin %{100 * out['background']['near_rate']:.1f}")
    for comp in COMPONENTS:
        d = per_component.get(comp)
        if not d:
            continue
        gaps = sorted(d["gaps"])
        entry = {
            "n": sum(d["same"].values()),
            "strand": fisher(d["same"][True], d["same"][False],
                             bg_same[True], bg_same[False]),
            "position": fisher(d["near"][True], d["near"][False],
                               bg_near[True], bg_near[False]),
            "offsets": {str(k): d["offsets"][k]
                        for k in range(-OFFSET_RANGE, OFFSET_RANGE + 1) if d["offsets"][k]},
            "gap_median": gaps[len(gaps) // 2] if gaps else None,
            "gap_quartiles": [gaps[len(gaps) // 4], gaps[3 * len(gaps) // 4]] if len(gaps) >= 4 else None,
            "gap_under_150": round(sum(1 for g in gaps if g <= 150) / len(gaps), 4) if gaps else None,
        }
        out["components"][comp] = entry
        print(f"[{comp:12s}] n={entry['n']:6,}  ayni iplik {100*entry['strand']['rate']:.1f}% "
              f"({entry['strand']['ratio']}x, p={entry['strand']['p']:.1e})  "
              f"yakin {100*entry['position']['rate']:.1f}% ({entry['position']['ratio']}x)  "
              f"medyan bosluk {entry['gap_median']} bp")

    ro_rows = con.execute(
        "SELECT candidate_id, strand, start, end FROM ro WHERE is_confirmed=1").fetchall()
    n_ro = len(ro_rows)
    sensitivity = {}
    for gap in GAP_THRESHOLDS:
        counts = walk_counts(ro_rows, neighbors, components, gap)
        sensitivity[str(gap)] = {c: counts.get(c, 0) for c in COMPONENTS}
        sensitivity[str(gap)]["any"] = 0
    # "en az bir ortak" ayrica hesaplanir (kume birlesimi gerekir)
    for gap in GAP_THRESHOLDS:
        total = 0
        for cid, strand, start, end in ro_rows:
            nbs = neighbors.get(cid, ())
            hit = False
            for side, direction in ((sorted((n for n in nbs if n[6] < 0), key=lambda n: -n[6]), -1),
                                    (sorted((n for n in nbs if n[6] > 0), key=lambda n: n[6]), 1)):
                prev_start, prev_end = start, end
                for n in side:
                    if n[2] != strand:
                        break
                    g = 0 if n[5] else (prev_start - n[1] if direction < 0 else n[0] - prev_end)
                    if g > gap:
                        break
                    if components.get(n[3]) in ("beta", "ferredoxin", "reductase"):
                        hit = True
                    prev_start, prev_end = n[0], n[1]
            total += hit
        sensitivity[str(gap)]["any"] = total
    out["threshold_sensitivity"] = {"n_ro": n_ro, "thresholds": GAP_THRESHOLDS,
                                    "counts": sensitivity}

    print(f"\n[esik duyarliligi] {n_ro:,} dogrulanmis RO")
    print(f"   {'esik':>6}  {'beta':>6} {'ferredoksin':>12} {'reduktaz':>9} {'en az biri':>11}")
    for gap in GAP_THRESHOLDS:
        s = sensitivity[str(gap)]
        print(f"   {gap:>4} bp  {s['beta']:>6} {s['ferredoxin']:>12} "
              f"{s['reductase']:>9} {s['any']:>11}")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\n[yazildi] {args.out}")
    con.close()


if __name__ == "__main__":
    main()
