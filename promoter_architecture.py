"""
Divergent duzenleyici gercekten daha kisa bir promotor bolgesi mi paylasiyor?

SORU. `analyze_regulation.py` her operonun 5' ucunun yukarisindaki intergenik
bolgeyi olcer. LysR tipi duzenleyiciler klasik olarak hedef operonla KAFA
KAFAYA oturur ve aradaki kisa bolgeyi paylasir. Eger bu mimari gercekse,
yukarisinda DUZENLEYICI olan divergent ciftlerin intergenik bolgesi,
yukarisinda baska bir gen olan divergent ciftlerinkinden KISA olmalidir.

NE OLCULMUYOR -- ve bu kasitli. Promotor DIZISI yok: veritabaninda komsu
genlerin protein cevirileri var, intergenik DNA yok. -35/-10 kutulari,
operator tekrarlari ve transkripsiyon baslangici icin nukleotid dizisi
gerekir. Dolayisiyla burada HICBIR promotor iddiasi yoktur; olculen sey
yalnizca "promotorun bulunmasi gereken bolgenin uzunlugu"dur. GC ve uzunluk
eslenmis bos modele karsi motif zenginlesmesi, dizinin cikarilmasini
bekliyor (ROADMAP 124.3'un kalan yarisi).

UC DUZEY, cunku 17.730 giriste p degeri neredeyse her seyi anlamli yapar:
  giris  -> her dogrulanmis RO bir gozlem
  dizi   -> ayni amino asit dizisi tek gozlem (tekrarlanan suslar cokertilir)
  cins   -> kume x cins basina tek gozlem; ASIL KARAR burasi

Ayrica CINS ICI esli sinama: her cinste iki grubun medyani ayri ayri alinir
ve yalnizca FARKIN ISARETI sinanir (Wilcoxon). Boylece "bazi cinslerin genomu
daha sikisik" turu bir kayma sonucu uretemez -- karsilastirma hep ayni cinsin
icinde kalir.

Cikti: analysis_out/promoter_architecture.json
"""

import argparse
import json
import os
import sqlite3
from collections import defaultdict

import numpy as np
from scipy import stats

MIN_PER_GENUS = 10          # cins ici esli sinamada her grup icin en az gozlem
ARCH_A = "divergent_regulator"
ARCH_B = "divergent_other"


def genus_of(organism):
    return (organism or "").split()[0] if organism else "?"


def cliffs_delta(a, b):
    """Parametrik olmayan etki buyuklugu: P(a>b) - P(a<b), -1..1."""
    a, b = np.asarray(a), np.asarray(b)
    # Siralama ile O(n log n): Mann-Whitney U'dan turetilir.
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2 * u / (len(a) * len(b)) - 1)


def compare(rows_a, rows_b):
    a = [r["gap"] for r in rows_a]
    b = [r["gap"] for r in rows_b]
    if len(a) < 3 or len(b) < 3:
        return {"n_regulator": len(a), "n_other": len(b), "usable": False}
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_regulator": len(a), "n_other": len(b), "usable": True,
        "median_regulator": float(np.median(a)),
        "median_other": float(np.median(b)),
        "median_difference": float(np.median(a) - np.median(b)),
        "cliffs_delta": round(cliffs_delta(a, b), 4),
        "p": float(u.pvalue),
    }


def collapse(rows, key):
    """Anahtar basina tek gozlem; bosluk degeri grubun MEDYANI."""
    groups = defaultdict(list)
    for r in rows:
        groups[key(r)].append(r)
    out = []
    for _k, items in groups.items():
        base = dict(items[0])
        base["gap"] = float(np.median([i["gap"] for i in items]))
        base["n_collapsed"] = len(items)
        out.append(base)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/promoter_architecture.json")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    rows = []
    for cid, gap, arch, cluster, seq, organism in con.execute("""
            SELECT g.candidate_id, g.intergenic_bp, g.architecture,
                   r.ro_cluster, r.sequence, p.organism
            FROM ro_regulation g
            JOIN ro r USING(candidate_id)
            JOIN replicon p ON p.nucleotide_id = r.nucleotide_id
            WHERE r.is_confirmed = 1
              AND g.intergenic_bp IS NOT NULL
              AND g.architecture IN (?, ?)""", (ARCH_A, ARCH_B)):
        rows.append({"candidate_id": cid, "gap": gap, "arch": arch,
                     "cluster": cluster, "sequence": seq or cid,
                     "genus": genus_of(organism)})

    levels = {}
    for name, data in (
            ("entry", rows),
            ("sequence", collapse(rows, lambda r: r["sequence"])),
            ("genus", collapse(rows, lambda r: (r["cluster"], r["genus"])))):
        levels[name] = compare([r for r in data if r["arch"] == ARCH_A],
                               [r for r in data if r["arch"] == ARCH_B])

    # CINS ICI ESLI SINAMA: her cins kendi icinde karsilastirilir.
    per_genus, paired = {}, []
    by_genus = defaultdict(lambda: {ARCH_A: [], ARCH_B: []})
    for r in rows:
        by_genus[r["genus"]][r["arch"]].append(r["gap"])
    for genus, d in sorted(by_genus.items()):
        if len(d[ARCH_A]) < MIN_PER_GENUS or len(d[ARCH_B]) < MIN_PER_GENUS:
            continue
        ma, mb = float(np.median(d[ARCH_A])), float(np.median(d[ARCH_B]))
        per_genus[genus] = {"n_regulator": len(d[ARCH_A]), "n_other": len(d[ARCH_B]),
                            "median_regulator": ma, "median_other": mb,
                            "difference": ma - mb}
        paired.append(ma - mb)

    if len(paired) >= 6:
        w = stats.wilcoxon(paired, alternative="two-sided")
        shorter = sum(1 for d in paired if d < 0)
        within = {"n_genera": len(paired), "genera_with_shorter_regulator_gap": shorter,
                  "median_of_differences": float(np.median(paired)),
                  "p": float(w.pvalue), "usable": True}
    else:
        within = {"n_genera": len(paired), "usable": False}

    report = {
        "question": "Does a divergent upstream regulator sit across a shorter "
                    "intergenic region than a divergent non-regulator?",
        "measures": "intergenic length in bp between the operon 5' end and the "
                    "first upstream gene. No promoter sequence is examined and "
                    "no promoter is predicted; the database holds protein "
                    "translations of neighbouring genes, not intergenic DNA.",
        "comparison": {"group_a": ARCH_A, "group_b": ARCH_B},
        "levels": levels,
        "within_genus": within,
        "per_genus": per_genus,
        "caveats": [
            "The regulator is identified by a product-name rule, so a gene whose "
            "product string does not say 'regulator' falls into the comparison "
            "group even if it regulates.",
            "A short intergenic region is consistent with a shared promoter but "
            "does not demonstrate one; nothing here measures transcription.",
            "Entries with no upstream gene on the replicon are absent from both "
            "groups, not counted as a long gap.",
        ],
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1)

    for name, lv in levels.items():
        if lv.get("usable"):
            print(f"[{name:8s}] duzenleyici medyan {lv['median_regulator']:6.0f} bp "
                  f"(n={lv['n_regulator']:5d})  digeri {lv['median_other']:6.0f} bp "
                  f"(n={lv['n_other']:5d})  fark {lv['median_difference']:+6.0f}  "
                  f"delta {lv['cliffs_delta']:+.3f}  p={lv['p']:.2e}")
    if within["usable"]:
        print(f"[cins ici] {within['n_genera']} cinsin "
              f"{within['genera_with_shorter_regulator_gap']}'inde duzenleyici "
              f"bosluğu daha kisa; farklarin medyani "
              f"{within['median_of_differences']:+.0f} bp, p={within['p']:.2e}")
    print(f"[yazildi] {args.out}")
    con.close()


if __name__ == "__main__":
    main()
