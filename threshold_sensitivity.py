"""
Dahil etme esiginin sonuclari NE KADAR degistirdigini olcer.

NEDEN BU SCRIPT VAR. Yontem sayfasi profil kapsama esigini 0,45 olarak veriyor
ve bu secimi KALIBRASYON setinde gerekcelendiriyor: 2.788 gercek alfa alt
birimi, 1.570 alfa olmayan protein. Ama kullanicinin sordugu soru farkliydi --
"o esik degisince burdaki bir suru sey degisebilir". Dogru cevap kalibrasyon
degil, VERITABANININ KENDI SONUCLARINI esik boyunca yeniden hesaplamaktir.

Hesap yeniden HMMER kosturmaz: kapsama, motif durumu ve karboksilat kimligi
giris basina saklandigi icin esik hareket ettirilip sayilar yeniden toplanabilir.

SINIR, ACIKCA. `ro_carboxylate` yalnizca ONAYLANMIS 11.422 giris icin olculdu.
Bu yuzden esigi 0,45'in ALTINA indirince yeni kabul edilen girislerin kalinti
kimligi BILINMIYOR ve o satirlarda kalinti sutunlari bos kalir. Esigi
yukseltmek bir alt kume secmek oldugundan orada her sey hesaplanabilir.

Cikti:
    analysis_out/threshold_sensitivity.json
"""

import argparse
import csv
import json
import os
import sqlite3
from collections import Counter, defaultdict

# Yayinlanan deger 0,45; cevresi hem asagi hem yukari taranir.
THRESHOLDS = (0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80)
PUBLISHED = 0.45


def cramers_v(table):
    """Chi-kare'den Cramer's V. scipy olmadan, kucuk tablolar icin."""
    rows = len(table)
    cols = len(table[0]) if rows else 0
    n = sum(sum(r) for r in table)
    if n == 0 or rows < 2 or cols < 2:
        return None
    row_sums = [sum(r) for r in table]
    col_sums = [sum(table[i][j] for i in range(rows)) for j in range(cols)]
    chi2 = 0.0
    for i in range(rows):
        for j in range(cols):
            expected = row_sums[i] * col_sums[j] / n
            if expected > 0:
                chi2 += (table[i][j] - expected) ** 2 / expected
    return (chi2 / (n * min(rows - 1, cols - 1))) ** 0.5


def load(connection, ecology_path):
    """Esikten BAGIMSIZ olan her seyi bir kez okur."""
    substrate_class = {}
    if os.path.exists(ecology_path):
        with open(ecology_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                substrate_class[row["cluster"]] = row.get("substrate_class", "")

    carboxylate = {cid: (cat, bridge) for cid, cat, bridge in connection.execute(
        "SELECT candidate_id, catalytic_residue, bridging_residue FROM ro_carboxylate")}

    rows = []
    for (cid, cluster, group, coverage, rieske, catalytic, plasmid,
         organism, domain) in connection.execute("""
            SELECT r.candidate_id, r.ro_cluster, r.ro_group, r.model_coverage,
                   r.rieske_intact, r.catalytic_intact, p.is_plasmid,
                   p.organism, d.domain
            FROM ro r JOIN replicon p USING(nucleotide_id)
            LEFT JOIN ro_domain d ON d.candidate_id = r.candidate_id
            WHERE r.model_coverage IS NOT NULL"""):
        rows.append({
            "id": cid, "cluster": cluster, "group": group,
            "coverage": coverage, "rieske": rieske, "catalytic": catalytic,
            "plasmid": bool(plasmid),
            "genus": (organism or "?").split()[0],
            "domain": domain,
            "sclass": substrate_class.get(cluster, ""),
            "carbox": carboxylate.get(cid),
        })
    return rows


def at_threshold(rows, threshold):
    """Esikteki kapi: motifler tam VE kapsama esigin uzerinde.

    Bu, yayinlanan kapidan biraz daha gevsek: `annotate_ro.py` ayrica bir durum
    etiketi uyguluyor ve 0,45'te 87 girisi daha disarida birakiyor. Fark
    raporda ACIKCA verilir, cunku burada amac mutlak sayiyi tekrar uretmek
    degil, esigin YONUNU ve BUYUKLUGUNU olcmek.
    """
    return [r for r in rows
            if r["rieske"] == 1 and r["catalytic"] == 1 and r["coverage"] >= threshold]


def metrics(kept, measured_only):
    """Bir esikteki ozet. measured_only: kalinti verisi olan girisler."""
    out = {
        "entries": len(kept),
        "types_with_members": len({r["cluster"] for r in kept
                                   if r["cluster"] and r["cluster"] != "N/A"}),
        "genera": len({r["genus"] for r in kept}),
        "plasmid_share": round(sum(r["plasmid"] for r in kept) / len(kept), 4) if kept else None,
        "eukaryote_share": round(sum(r["domain"] == "Eukaryota" for r in kept) / len(kept), 4)
                           if kept else None,
    }
    per_group = Counter(r["group"] for r in kept if r["group"])
    out["group_shares"] = {g: round(n / len(kept), 4) for g, n in sorted(per_group.items())} \
        if kept else {}

    # Ksenobiyotik girislerin plazmid orani -- sitenin yayinladigi bulgulardan
    # biri. Esik hareket edince yonu degisiyor mu?
    xen = [r for r in kept if r["sclass"] == "xenobiotic"]
    nat = [r for r in kept if r["sclass"] and r["sclass"] != "xenobiotic"]
    if xen and nat:
        a = sum(r["plasmid"] for r in xen) / len(xen)
        b = sum(r["plasmid"] for r in nat) / len(nat)
        out["plasmid_xenobiotic_share"] = round(a, 4)
        out["plasmid_natural_share"] = round(b, 4)
        out["plasmid_ratio"] = round(a / b, 3) if b else None

    # Koprüleyen karboksilat Asp/Glu -- en guclu yayinlanan iliski.
    carbox = [r for r in kept if r["carbox"] and r["carbox"][1] in ("D", "E") and r["group"]]
    if measured_only and carbox:
        groups = sorted({r["group"] for r in carbox})
        table = [[sum(1 for r in carbox if r["group"] == g and r["carbox"][1] == res)
                  for res in ("D", "E")] for g in groups]
        n_asp = sum(row[0] for row in table)
        n_glu = sum(row[1] for row in table)
        out["bridging_measured"] = len(carbox)
        out["bridging_asp"] = n_asp
        out["bridging_glu"] = n_glu
        out["bridging_asp_per_glu"] = round(n_asp / n_glu, 2) if n_glu else None
        out["bridging_cramers_v"] = round(cramers_v(table), 3) if cramers_v(table) else None
    else:
        out["bridging_measured"] = None
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--ecology", default="cluster_ecology.csv")
    parser.add_argument("--out", default="analysis_out/threshold_sensitivity.json")
    args = parser.parse_args()

    connection = sqlite3.connect(args.db)
    rows = load(connection, args.ecology)
    published_count = connection.execute(
        "SELECT COUNT(*) FROM ro WHERE is_confirmed=1").fetchone()[0]

    levels = []
    for threshold in THRESHOLDS:
        kept = at_threshold(rows, threshold)
        # Kalinti verisi yalnizca onaylanmis girisler icin olculdu; esik
        # yayinlanan degerin ALTINDA ise yeni girislerin kalinti kimligi yok.
        measured = threshold >= PUBLISHED
        entry = {"threshold": threshold, "residues_recomputable": measured}
        entry.update(metrics(kept, measured))
        if not measured:
            entry["unmeasured_entries"] = sum(1 for r in kept if r["carbox"] is None)
        levels.append(entry)

    baseline = next(l for l in levels if l["threshold"] == PUBLISHED)
    for level in levels:
        level["entries_vs_published"] = round(
            level["entries"] / baseline["entries"], 3) if baseline["entries"] else None

    result = {
        "method": {
            "question": ("how far do the database's own numbers move when the profile "
                         "coverage threshold moves"),
            "gate_recomputed": ("all four Rieske ligands present, catalytic triad present, "
                                "and profile coverage at or above the threshold"),
            "difference_from_published_gate": (
                f"the published pipeline applies one further status check and reports "
                f"{published_count} confirmed entries at coverage {PUBLISHED}, while this "
                f"recomputation keeps {baseline['entries']}; the gap is "
                f"{baseline['entries'] - published_count} borderline entries. The purpose "
                f"here is the direction and size of the change, not reproducing the "
                f"absolute count"),
            "residue_limit": ("ro_carboxylate was measured only for the confirmed set, so "
                              "below the published threshold the residue columns are left "
                              "empty rather than guessed"),
            "published_threshold": PUBLISHED,
            "published_confirmed_entries": published_count,
        },
        "levels": levels,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1, sort_keys=False)

    print("=" * 78)
    print("ESIK DUYARLILIGI")
    print("=" * 78)
    print(f"{'esik':>6} {'giris':>7} {'oran':>6} {'tip':>4} {'cins':>5} "
          f"{'plazmid':>8} {'okaryot':>8} {'Asp:Glu':>8} {'V':>6}")
    for level in levels:
        print(f"{level['threshold']:6.2f} {level['entries']:7,} "
              f"{level['entries_vs_published']:6.2f} {level['types_with_members']:4} "
              f"{level['genera']:5} "
              f"{(level['plasmid_share'] or 0) * 100:7.1f}% "
              f"{(level['eukaryote_share'] or 0) * 100:7.1f}% "
              f"{str(level.get('bridging_asp_per_glu') or '-'):>8} "
              f"{str(level.get('bridging_cramers_v') or '-'):>6}")
    print(f"\n[yazildi] {args.out}")


if __name__ == "__main__":
    main()
