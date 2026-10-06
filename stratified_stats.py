"""
Katalitik merkez kalintisinin iliskisini, UYENIN YAKINLIGINA gore yeniden olcer.

NEDEN BU SCRIPT VAR. Site "kopruleyen karboksilat ile enzim grubu arasinda
V=0,84 iliski var" diyordu, tek bir sayi olarak. Ama bu sayi 11.417 girisin
hepsini esit sayiyor, oysa girislerin cogu kuratorlu enzimden UZAK ve uzak bir
uyenin reaksiyonu bilinmiyor. Dolayisiyla sayinin ne kadari bilinen kimyadan,
ne kadari tanimsiz akrabalardan geliyor sorusu meşru.

Olcum cevabi veriyor ve rahatlaticidir: iliski yakin uyelerde ZAYIFLAMIYOR,
neredeyse kusursuzlasiyor (V=0,999) ve uzaklikla bozuluyor (V=0,805). Yani
uzak uyeler sinyali URETMIYOR, uzerine gurultu ekliyor.

Ayrica bir KURATORLUK IPUCU cikti: grup 2'nin uzak uyelerinden 62'si Glu
tasiyor, oysa ayni grubun yakin uyelerinin %0'i tasiyor. Bunlarin en yakin
referanslari agirlikla grup 3 ve grup 5'ten ve kimlik ortancasi %32. Yani
bunlar muhtemelen yanlis gruba atanmis.

Cikti:
    analysis_out/stratified_stats.json      ozet + tabakalar
    analysis_out/carboxylate_rows.json      tarayicida yeniden hesap icin kompakt tablo
"""

import argparse
import json
import os
import sqlite3
from collections import Counter, defaultdict

TIERS = ["characterized", "close_homolog", "family_member", "distant", "novel"]
CLASSES = ["core", "divergent", "novel_candidate"]
# "Yakin" tanimi: substrat etiketinin bilgilendirici sayildigi iki kademe.
# Sitenin baska yerlerinde de ayni esik (referansa >=%60 kimlik) kullaniliyor.
CLOSE = ("characterized", "close_homolog")


def chi2_and_v(table):
    """Chi-kare, serbestlik derecesi ve Cramer's V; scipy olmadan.

    Beklenen deger sifir olan satir/sutun ATILIR, yoksa bolme hatasi olur ve
    bir tabaka tek bir bos hucre yuzunden olculemez hale gelir.
    """
    table = [row for row in table if sum(row) > 0]
    if len(table) < 2:
        return None
    cols = len(table[0])
    col_sums = [sum(row[j] for row in table) for j in range(cols)]
    keep = [j for j in range(cols) if col_sums[j] > 0]
    if len(keep) < 2:
        return None
    table = [[row[j] for j in keep] for row in table]
    rows_n, cols_n = len(table), len(table[0])
    n = sum(sum(r) for r in table)
    row_sums = [sum(r) for r in table]
    col_sums = [sum(table[i][j] for i in range(rows_n)) for j in range(cols_n)]
    chi2 = 0.0
    for i in range(rows_n):
        for j in range(cols_n):
            expected = row_sums[i] * col_sums[j] / n
            if expected > 0:
                chi2 += (table[i][j] - expected) ** 2 / expected
    dof = (rows_n - 1) * (cols_n - 1)
    v = (chi2 / (n * min(rows_n - 1, cols_n - 1))) ** 0.5
    return {"chi2": round(chi2, 1), "dof": dof, "cramers_v": round(v, 3), "n": n}


def load(connection):
    rows = []
    for cluster, group, residue, catalytic, tier, identity, klass in connection.execute("""
            SELECT r.ro_cluster, r.ro_group, c.bridging_residue, c.catalytic_residue,
                   e.tier, e.ref_identity, s.assignment_class
            FROM ro_carboxylate c
            JOIN ro r USING(candidate_id)
            LEFT JOIN ro_evidence e USING(candidate_id)
            LEFT JOIN ro_subfamily s USING(candidate_id)
            WHERE r.ro_cluster IS NOT NULL AND r.ro_cluster <> 'N/A'"""):
        rows.append({"cluster": cluster, "group": group or "?",
                     "bridging": residue, "catalytic": catalytic,
                     "tier": tier or "unknown",
                     "identity": identity,
                     "class": klass or "unknown"})
    return rows


def association(rows, key="group", residue="bridging"):
    """key x (Asp/Glu) capraz tablosu ve V."""
    usable = [r for r in rows if r[residue] in ("D", "E")]
    keys = sorted({r[key] for r in usable})
    table = [[sum(1 for r in usable if r[key] == k and r[residue] == res)
              for res in ("D", "E")] for k in keys]
    stat = chi2_and_v(table)
    if stat is None:
        return None
    glu = sum(row[1] for row in table)
    stat.update({"glu": glu, "asp": stat["n"] - glu,
                 "glu_share": round(glu / stat["n"], 4) if stat["n"] else None,
                 "rows": keys,
                 "table": table})
    return stat


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default="analysis_out")
    args = parser.parse_args()

    connection = sqlite3.connect(args.db)
    rows = load(connection)

    strata = {}
    strata["all"] = association(rows)
    strata["close"] = association([r for r in rows if r["tier"] in CLOSE])
    strata["distant"] = association([r for r in rows if r["tier"] not in CLOSE
                                     and r["tier"] != "unknown"])
    for tier in TIERS:
        got = association([r for r in rows if r["tier"] == tier])
        if got:
            strata[f"tier:{tier}"] = got
    for klass in CLASSES:
        got = association([r for r in rows if r["class"] == klass])
        if got:
            strata[f"class:{klass}"] = got

    # Grup icinde tier'a gore Glu payi: tabloyu bozan yer burasi.
    within = defaultdict(dict)
    for group in sorted({r["group"] for r in rows}):
        for tier in TIERS:
            subset = [r for r in rows if r["group"] == group and r["tier"] == tier
                      and r["bridging"] in ("D", "E")]
            if len(subset) >= 10:
                glu = sum(1 for r in subset if r["bridging"] == "E")
                within[group][tier] = {"n": len(subset), "glu": glu,
                                       "glu_share": round(glu / len(subset), 4)}

    # Tip icinde tier'lar arasinda Glu payi ne kadar oynuyor?
    unstable = []
    for cluster in sorted({r["cluster"] for r in rows}):
        shares = {}
        for tier in TIERS:
            subset = [r for r in rows if r["cluster"] == cluster and r["tier"] == tier
                      and r["bridging"] in ("D", "E")]
            if len(subset) >= 10:
                shares[tier] = round(sum(1 for r in subset if r["bridging"] == "E")
                                     / len(subset), 4)
        if len(shares) >= 2:
            spread = max(shares.values()) - min(shares.values())
            if spread > 0.25:
                unstable.append({"cluster": cluster, "spread": round(spread, 3),
                                 "shares": shares})
    unstable.sort(key=lambda r: -r["spread"])

    # Kuratorluk ipucu: yakin uyeleri %0 Glu olan bir grupta Glu tasiyan uzak
    # uyeler. Bunlar muhtemelen yanlis gruba atanmis.
    suspects = []
    for group, tiers in within.items():
        close_share = max((tiers.get(t, {}).get("glu_share", 0) for t in CLOSE), default=0)
        far_share = max((tiers.get(t, {}).get("glu_share", 0)
                         for t in ("family_member", "distant", "novel")), default=0)
        if close_share <= 0.02 and far_share >= 0.05:
            found = connection.execute("""
                SELECT r.ro_cluster, e.nearest_ref, e.ref_identity
                FROM ro_carboxylate c JOIN ro r USING(candidate_id)
                JOIN ro_evidence e USING(candidate_id)
                WHERE c.bridging_residue='E' AND r.ro_group=?
                  AND e.tier IN ('family_member','distant','novel')""", (group,)).fetchall()
            if found:
                nearest_groups = Counter((f[1] or "?").split("_")[0] for f in found)
                identities = sorted(f[2] for f in found if f[2] is not None)
                suspects.append({
                    "group": group,
                    "close_glu_share": close_share,
                    "distant_glu_share": far_share,
                    "entries": len(found),
                    "assigned_types": Counter(f[0] for f in found).most_common(5),
                    "nearest_reference_group": nearest_groups.most_common(),
                    "median_identity_to_nearest": (identities[len(identities) // 2]
                                                   if identities else None),
                })
    suspects.sort(key=lambda r: -r["entries"])

    # Tarayicida yeniden hesap icin kompakt tablo. Dizgiler TEKRARLANMAZ:
    # 11.422 satirin her birinde tip, tier ve sinif adini duz metin yazmak
    # 508 KB ediyordu; sozluk + tamsayi indeksiyle ayni veri dortte bire
    # iniyor ve yayinlanan sayfa o kadar hizli aciliyor.
    def factorise(values):
        order = sorted(set(values))
        return order, {v: i for i, v in enumerate(order)}

    clusters, cluster_ix = factorise(r["cluster"] for r in rows)
    tiers_seen, tier_ix = factorise(r["tier"] for r in rows)
    classes_seen, class_ix = factorise(r["class"] for r in rows)
    groups_seen, group_ix = factorise(r["group"] for r in rows)
    residues = ["", "D", "E"]
    residue_ix = {v: i for i, v in enumerate(residues)}
    compact = [[group_ix[r["group"]], cluster_ix[r["cluster"]],
                residue_ix.get(r["bridging"] or "", 0),
                tier_ix[r["tier"]], class_ix[r["class"]]]
               for r in rows]

    result = {
        "method": {
            "question": ("how much of the bridging-carboxylate association comes from "
                         "members whose chemistry is actually known"),
            "close_definition": ("evidence levels " + " and ".join(CLOSE) +
                                 ", that is at least 60 % identity to a curated enzyme"),
            "why": ("a single association computed over every entry weights a distant "
                    "relative of unknown reaction the same as a characterised enzyme; "
                    "stratifying by evidence level says whether the signal depends on "
                    "the uncertain members"),
            "answer": ("it does not: the association is near-deterministic among close "
                       "members and degrades with distance, so distant members add noise "
                       "rather than creating the pattern"),
            "tiers": TIERS, "classes": CLASSES,
        },
        "strata": strata,
        "glu_share_within_group_by_tier": {g: v for g, v in within.items()},
        "types_unstable_across_tiers": unstable,
        "misassignment_suspects": suspects,
        "totals": {"entries": len(rows),
                   "with_bridging_residue": sum(1 for r in rows
                                                if r["bridging"] in ("D", "E"))},
    }
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "stratified_stats.json"), "w",
              encoding="utf-8") as fh:
        json.dump(result, fh, indent=1)
    with open(os.path.join(args.out_dir, "carboxylate_rows.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"columns": ["group", "cluster", "bridging", "tier", "class"],
                   "note": ("every column is an index into the matching list below; "
                            "the strings are stored once instead of 11,422 times"),
                   "groups": groups_seen, "clusters": clusters,
                   "residues": residues, "tiers": tiers_seen, "classes": classes_seen,
                   "rows": compact}, fh, separators=(",", ":"))

    print("=" * 74)
    print("KOPRULEYEN KARBOKSILAT -- YAKINLIGA GORE TABAKALANMIS")
    print("=" * 74)
    print(f"{'stratum':26} {'n':>7} {'%Glu':>7} {'V':>7} {'chi2':>9}")
    for name in ("all", "close", "distant"):
        s = strata.get(name)
        if s:
            print(f"{name:26} {s['n']:7} {100*s['glu_share']:6.1f}% "
                  f"{s['cramers_v']:7.3f} {s['chi2']:9.0f}")
    for name, s in strata.items():
        if name.startswith(("tier:", "class:")):
            print(f"{name:26} {s['n']:7} {100*s['glu_share']:6.1f}% "
                  f"{s['cramers_v']:7.3f} {s['chi2']:9.0f}")
    print("\nGlu payi, grup icinde tier'a gore:")
    for group, tiers in sorted(within.items()):
        parts = [f"{t[:5]}={100*v['glu_share']:.0f}%({v['n']})" for t, v in tiers.items()]
        print(f"  group {group}: " + "  ".join(parts))
    print(f"\ntier'lar arasinda oynak tip: {len(unstable)}")
    for item in unstable[:6]:
        print(f"  {item['cluster']:16} spread {100*item['spread']:.0f} pts  "
              + " ".join(f"{t[:5]}={100*v:.0f}%" for t, v in item['shares'].items()))
    print(f"\nyanlis atama supheli grup: {len(suspects)}")
    for item in suspects:
        print(f"  group {item['group']}: {item['entries']} uzak uye Glu tasiyor, "
              f"yakin uyelerde %{100*item['close_glu_share']:.0f}; "
              f"en yakin referanslarin grubu {item['nearest_reference_group']}, "
              f"kimlik ortancasi %{item['median_identity_to_nearest']}")
    print(f"\n[yazildi] {args.out_dir}/stratified_stats.json ve carboxylate_rows.json")


if __name__ == "__main__":
    main()
