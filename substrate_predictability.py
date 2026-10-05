"""
Kimlik substrati ne kadar ongoruyor? -- veritabaninin merkezi uyarisinin olcumu.

Her yerde "kimlik substrat garantisi vermez" yaziyoruz ve ornek olarak EdoA1/CumA1
cifti (%99,8 kimlik, farkli substrat) veriliyor. Bu tek ornek. Bu script ayni
iddiayi SINIFLANDIRICI PERFORMANSI olarak nicelestirir:

    "iki protein %X'ten fazla benziyorsa ayni substrata etki eder" kuralini
    71 kuratorlu referansin TUM ciftlerinde sina.

Olculenler:
    - ROC egrisi altindaki alan (AUC): kimlik, "ayni substrat" ikilisini ne kadar
      ayirir. 0,5 = hic bilgi yok, 1,0 = kusursuz.
    - Her esik icin kesinlik (precision) ve duyarlilik (recall).
    - %90 kesinlige ulasan esik var mi (yoksa bu, hicbir esigin guvenli
      olmadigi anlamina gelir).
    - Ayni olcum SUBSTRAT SINIFI (ksenobiyotik / dogal aromatik / ...) icin de
      yapilir; sinif daha kaba oldugu icin ongorulebilirligi daha yuksek olmasi
      beklenir ve bu beklenti sinanir.

Girdi: analysis_out/reference_pairs.csv (evidence_tiers.py uretir),
       cluster_ecology.csv (substrat sinifi)
Cikti: analysis_out/substrate_predictability.json
"""

import argparse
import csv
import json
import os

UNKNOWN = {"bilinmiyor", "unknown", "", "?"}


def roc_auc(scores_labels):
    """Mann-Whitney U uzerinden AUC. scores_labels: [(skor, 0/1), ...]"""
    pos = [s for s, y in scores_labels if y == 1]
    neg = [s for s, y in scores_labels if y == 0]
    if not pos or not neg:
        return None
    # baglari 0,5 sayan siralamali toplam
    merged = sorted(((s, y) for s, y in scores_labels), key=lambda t: t[0])
    ranks = {}
    i = 0
    while i < len(merged):
        j = i
        while j + 1 < len(merged) and merged[j + 1][0] == merged[i][0]:
            j += 1
        average = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks.setdefault(k, average)
        i = j + 1
    rank_sum = sum(ranks[k] for k, (s, y) in enumerate(merged) if y == 1)
    return (rank_sum - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def curve(scores_labels, thresholds):
    """Her esik icin kesinlik/duyarlilik: "kimlik >= esik ise ayni substrat" kurali.

    Kesinlik TEK BASINA yanlis okunur: siniflar dengesizse (ornegin ciftlerin
    cogu zaten ayni substrat SINIFINA giriyorsa) hicbir bilgi tasimayan bir
    kural da yuksek kesinlik verir. Bu yuzden her satirda "lift" de yazilir:
    kesinligin taban orana (prevalence) bolumu. Lift ~1 ise esik bilgi
    katmiyordur.
    """
    out = []
    total_pos = sum(1 for _, y in scores_labels if y == 1)
    prevalence = total_pos / len(scores_labels) if scores_labels else 0
    for t in thresholds:
        tp = sum(1 for s, y in scores_labels if s >= t and y == 1)
        fp = sum(1 for s, y in scores_labels if s >= t and y == 0)
        out.append({
            "threshold": t,
            "pairs_at_or_above": tp + fp,
            "same": tp, "different": fp,
            "precision": round(tp / (tp + fp), 4) if tp + fp else None,
            "recall": round(tp / total_pos, 4) if total_pos else None,
            "lift": round((tp / (tp + fp)) / prevalence, 2)
                    if (tp + fp) and prevalence else None,
        })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default="analysis_out/reference_pairs.csv")
    ap.add_argument("--ecology", default="cluster_ecology.csv")
    ap.add_argument("--out", default="analysis_out/substrate_predictability.json")
    args = ap.parse_args()

    classes = {}
    with open(args.ecology) as fh:
        for row in csv.DictReader(fh):
            classes[row["cluster"]] = row.get("substrate_class", "unknown")

    rows = list(csv.DictReader(open(args.pairs)))
    exact, by_class = [], []
    pairs_kept = []
    for r in rows:
        a, b = r["substrate_a"].strip().lower(), r["substrate_b"].strip().lower()
        identity = float(r["identity"])
        if a in UNKNOWN or b in UNKNOWN:
            continue
        exact.append((identity, 1 if a == b else 0))
        ca, cb = classes.get(r["ref_a"], "unknown"), classes.get(r["ref_b"], "unknown")
        if ca not in UNKNOWN and cb not in UNKNOWN and "unknown" not in (ca, cb):
            by_class.append((identity, 1 if ca == cb else 0))
        pairs_kept.append(r)

    thresholds = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99]
    out = {
        "pairs_total": len(rows),
        "pairs_with_both_substrates_known": len(exact),
        "pairs_with_both_classes_known": len(by_class),
        "same_substrate_pairs": sum(y for _, y in exact),
        "exact_substrate": {
            "auc": roc_auc(exact), "curve": curve(exact, thresholds),
            "prevalence": round(sum(y for _, y in exact) / len(exact), 4) if exact else None},
        "substrate_class": {
            "auc": roc_auc(by_class), "curve": curve(by_class, thresholds),
            "prevalence": round(sum(y for _, y in by_class) / len(by_class), 4) if by_class else None},
    }

    # En yuksek kimlikli, FARKLI substratli cift: tek cumlelik uyari icin
    worst = max((r for r in pairs_kept if r["same_substrate_label"] == "0"),
                key=lambda r: float(r["identity"]), default=None)
    best_same = min((r for r in pairs_kept if r["same_substrate_label"] == "1"),
                    key=lambda r: float(r["identity"]), default=None)
    out["most_similar_different_substrate"] = worst
    out["least_similar_same_substrate"] = best_same

    # %90 ve %95 kesinlige ulasan en kucuk esik (varsa)
    for key in ("exact_substrate", "substrate_class"):
        for target in (0.9, 0.95):
            hit = next((p["threshold"] for p in out[key]["curve"]
                        if p["precision"] is not None and p["precision"] >= target
                        and p["pairs_at_or_above"] >= 10), None)
            out[key][f"threshold_for_precision_{int(target * 100)}"] = hit

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"[cift] {len(rows)} referans cifti, {len(exact)} tanesinde iki substrat da biliniyor "
          f"({out['same_substrate_pairs']} tanesi ayni substrat)")
    print(f"[AUC] tam substrat : {out['exact_substrate']['auc']:.3f}")
    print(f"[AUC] substrat sinifi: {out['substrate_class']['auc']:.3f}  (n={len(by_class)})")
    print("\n[egri] 'kimlik >= esik ise ayni substrat' kuralinin kesinligi:")
    print(f"   taban oran (prevalence): {out['exact_substrate']['prevalence']}")
    print(f"   {'esik':>5} {'cift':>6} {'ayni':>6} {'farkli':>7} {'kesinlik':>9} {'lift':>6} {'duyarlilik':>11}")
    for p in out["exact_substrate"]["curve"]:
        if p["pairs_at_or_above"]:
            print(f"   {p['threshold']:>4}% {p['pairs_at_or_above']:>6} {p['same']:>6} "
                  f"{p['different']:>7} {p['precision']:>9} {p['lift']:>6} {p['recall']:>11}")
    print(f"\n[sinif] taban oran {out['substrate_class']['prevalence']} -- sinif duzeyinde "
          f"yuksek kesinlik buyuk olcude bu taban orandan gelir, lift'e bakilmali:")
    for p in out["substrate_class"]["curve"]:
        if p["threshold"] in (30, 50, 70, 90) and p["pairs_at_or_above"]:
            print(f"   {p['threshold']:>4}% kesinlik {p['precision']} lift {p['lift']}")
    print(f"\n[esik] tam substrat icin %90 kesinlik: "
          f"{out['exact_substrate']['threshold_for_precision_90'] or 'HICBIR ESIKTE YOK'}")
    print(f"[esik] substrat sinifi icin %90 kesinlik: "
          f"{out['substrate_class']['threshold_for_precision_90'] or 'HICBIR ESIKTE YOK'}")
    if worst:
        print(f"[en kotu] {worst['ref_a']} / {worst['ref_b']} %{worst['identity']} kimlik, "
              f"{worst['substrate_a']} vs {worst['substrate_b']}")
    print(f"[yazildi] {args.out}")


if __name__ == "__main__":
    main()
