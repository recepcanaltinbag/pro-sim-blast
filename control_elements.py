"""Duzenleyici aileleri ve IS aileleri: her birinin sayfasi icin veri.

Kullanici iki sey istedi ve ikisi de ayni bicimde cevaplaniyor:

  "TetR regulatorunun yonettikleri"   -- bir duzenleyici ailesinin yanindaki
                                        enzim tipleri, kimyasal aileler,
                                        reaksiyon siniflari
  "belli bi transpozonun iliskili     -- bir IS ailesinin yanindaki enzim
   oldugu enzimler"                     tipleri, ve asil soru: bu iliski
                                        ENZIME mi yoksa TURE mi bagli

Ucuncu soruyu kullanici kendisi sordu: "transposonlarin genel olarak tur
ozelinde mi enzim ozelinde mi oldugu". Olculebilir bir soru ve cevabi tek bir
sayida degil, iki sayinin KARSILASTIRILMASINDA: ayni IS ailesi tablosu bir kez
enzim tipine, bir kez cinse gore kurulur ve hangisinde daha guclu bir birliktelik
cikti karsilastirilir. Iki tablo farkli sekilde oldugu icin ham Cramér's V
karsilastirilamaz; her ikisi icin etiket karistirmasiyla bir null kurulur ve
NULL'IN UZERINDEKI fazla karsilastirilir.

Sayilan sey hakkinda iki uyari, ikisi de ciktiya yaziliyor:

  Duzenleyici KONUMLA belirlenir -- operonun 5' ucunun yukarisindaki ilk gen.
  Baglanma olculmedi. Uzaktan ya da kuresel bir duzenleyici bu yontemle
  gorunmez, ve yukarida duran her gen duzenleyici degildir.

  "Transpozon" etiketi regex ile konuyor ve icine rekombinasyon makinesi de
  giriyor: urun adlarinda 117 kez "Holliday junction resolvase RuvX" var, ki o
  bir mobil element degil. Bu yuzden IS AILELERI ayri ayri cikariliyor
  (urun adindaki "IS<n> family" kalibi) ve adlandirilmamis transpozazlardan
  ayri tutuluyor.

Cikti: analysis_out/control_elements.json
"""

import argparse
import csv
import json
import os
import re
import sqlite3
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

from stats_overview import bh_adjust, cramers_v, signal_map, trim_table

# "IS3 family transposase", "IS5/IS1182 family transposase" -> IS3, IS5/IS1182
IS_FAMILY_RE = re.compile(r"\b(IS\d+[A-Za-z]*(?:/IS\d+[A-Za-z]*)?)\s+family\b", re.I)

MIN_ENTRIES = 20        # bu sayinin altindaki aile icin sayfa kurulmaz
MIN_TYPE_ENTRIES = 10   # tabloda bir satir olmak icin gereken en az giris
PERMUTATIONS = 200
SEED = 7


def load(con, chemistry_path):
    chem = {}
    with open(chemistry_path) as fh:
        for row in csv.DictReader(fh):
            chem[row["cluster"]] = {
                "reaction_class": row.get("reaction_class") or "unknown",
                "chem_family": row.get("family") or "unknown",
            }
    rows = con.execute("""
        SELECT r.candidate_id, r.ro_cluster, r.ro_group, p.organism, p.is_plasmid,
               g.upstream_family, g.upstream_category, g.architecture, g.intergenic_bp
        FROM ro r
        JOIN replicon p USING(nucleotide_id)
        LEFT JOIN ro_regulation g USING(candidate_id)
        WHERE r.is_confirmed = 1
    """).fetchall()
    data = {}
    for (cid, cluster, group, organism, plasmid, fam, cat, arch, inter) in rows:
        c = chem.get(cluster, {})
        data[cid] = {
            "candidate_id": cid, "cluster": cluster, "group": group or "?",
            "genus": (organism or "?").split()[0],
            "species": " ".join((organism or "?").split()[:2]),
            "is_plasmid": int(bool(plasmid)),
            "regulator_family": fam if cat == "regulator" else None,
            "architecture": arch or "unknown",
            "intergenic_bp": inter,
            "reaction_class": c.get("reaction_class", "unknown"),
            "chem_family": c.get("chem_family", "unknown"),
            "is_families": set(),
            "unnamed_transposase": 0,
            "recombination_machinery": 0,
        }

    # mobil elemanlar: IS ailesi urun adindan cikarilir, adlandirilmamislar ve
    # rekombinasyon makinesi ayri sayilir.
    for cid, product in con.execute("""
        SELECT nb.candidate_id, nb.product
        FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id
                            AND c.category = 'transposon'
        JOIN ro r ON r.candidate_id = nb.candidate_id
        WHERE r.is_confirmed = 1 AND nb.product IS NOT NULL
    """):
        d = data.get(cid)
        if not d:
            continue
        m = IS_FAMILY_RE.search(product)
        if m:
            d["is_families"].add(m.group(1).upper())
        elif re.search(r"resolvase|holliday|ruv", product, re.I):
            d["recombination_machinery"] += 1
        else:
            d["unnamed_transposase"] += 1
    return list(data.values())


def enrichment(rows, key_field, keys, col_field, cols, min_expected=5.0):
    """Anahtar x sutun tablosu + hucre sinyalleri + cins duzeyinde tekrar."""
    table = [[sum(1 for d in rows if d[key_field] == k and d[col_field] == c)
              for c in cols] for k in keys]
    keys_t, cols_t, table_t = trim_table(keys, cols, table)
    if len(keys_t) < 2 or len(cols_t) < 2:
        return None
    seen, gsub = set(), []
    for d in rows:
        token = (d[key_field], d["genus"], d[col_field])
        if token in seen:
            continue
        seen.add(token)
        gsub.append(d)
    gtable = [[sum(1 for d in gsub if d[key_field] == k and d[col_field] == c)
               for c in cols_t] for k in keys_t]
    chi2, p, dof, _ = stats.chi2_contingency(table_t)
    return {"rows": keys_t, "cols": cols_t, "table": table_t,
            "genus_table": gtable,
            "chi2": float(chi2), "dof": int(dof), "p": float(p),
            "cramers_v": cramers_v(table_t),
            "signals": signal_map(keys_t, cols_t, table_t, gtable,
                                  min_expected=min_expected)}


def permuted_excess(rows, key_field, col_field, repeats=PERMUTATIONS, seed=SEED):
    """Olculen V ile etiket karistirmasinin urettigi V arasindaki FARK.

    Enzim tipi tablosu ~40 satir, cins tablosu ~200 satir olabilir. Cramér's V
    satir sayisiyla kendiliginden yukseldigi icin iki tabloyu dogrudan
    karsilastirmak yaniltir. Her tablo kendi null'una gore olculur ve yalnizca
    null'in UZERINDEKI fazla karsilastirilir.
    """
    import random as _random
    keys = sorted({d[key_field] for d in rows})
    cols = sorted({d[col_field] for d in rows})
    table = [[sum(1 for d in rows if d[key_field] == k and d[col_field] == c)
              for c in cols] for k in keys]
    keys_t, cols_t, table_t = trim_table(keys, cols, table)
    if len(keys_t) < 2 or len(cols_t) < 2:
        return None
    observed = cramers_v(table_t)
    rng = _random.Random(seed)
    labels = [d[col_field] for d in rows]
    kindex = {k: i for i, k in enumerate(keys_t)}
    cindex = {c: j for j, c in enumerate(cols_t)}
    nulls = []
    for _ in range(repeats):
        rng.shuffle(labels)
        t = [[0] * len(cols_t) for _ in keys_t]
        for d, lab in zip(rows, labels):
            i, j = kindex.get(d[key_field]), cindex.get(lab)
            if i is not None and j is not None:
                t[i][j] += 1
        v = cramers_v(trim_table(keys_t, cols_t, t)[2])
        if v is not None:
            nulls.append(v)
    if not nulls:
        return None
    mean = sum(nulls) / len(nulls)
    sd = (sum((x - mean) ** 2 for x in nulls) / len(nulls)) ** 0.5
    return {"observed": round(observed, 4), "null_mean": round(mean, 4),
            "null_sd": round(sd, 4), "excess": round(observed - mean, 4),
            "z": round((observed - mean) / sd, 2) if sd else None,
            "rows": len(keys_t), "cols": len(cols_t), "n": sum(sum(r) for r in table_t),
            "repeats": len(nulls)}


def family_profile(rows, all_rows, label, kind):
    """Tek bir aile icin sayfa icerigi."""
    types = Counter(d["cluster"] for d in rows)
    bg_types = Counter(d["cluster"] for d in all_rows)
    bg_total = len(all_rows)
    enriched = []
    praw, order = [], []
    for cluster, n in types.most_common():
        share = n / len(rows)
        base = bg_types[cluster] / bg_total if bg_total else 0.0
        # Fisher: bu tip, bu ailede mi yogunlasiyor
        a = n
        b = len(rows) - n
        c = bg_types[cluster] - n
        d_ = (bg_total - len(rows)) - c
        if min(a + b, c + d_) == 0:
            continue
        p = float(stats.fisher_exact([[a, b], [c, d_]])[1])
        praw.append(p)
        order.append({
            "cluster": cluster, "entries": n,
            "genera": len({x["genus"] for x in rows if x["cluster"] == cluster}),
            "share": round(share, 4), "baseline": round(base, 4),
            "lift": round(share - base, 4),
            "fold": round(share / base, 2) if base else None,
            "p": p,
        })
    for item, q in zip(order, bh_adjust(praw)):
        item["q"] = q
        item["flagged"] = bool(q < 0.05 and abs(item["lift"]) >= 0.05)
    enriched = sorted(order, key=lambda r: (0 if r["flagged"] else 1, -abs(r["lift"])))

    inter = [d["intergenic_bp"] for d in rows
             if d["intergenic_bp"] is not None and 0 <= d["intergenic_bp"] <= 2000]
    return {
        "family": label, "kind": kind,
        "entries": len(rows),
        "types": len(types),
        "genera": len({d["genus"] for d in rows}),
        "species": len({d["species"] for d in rows}),
        "plasmid_rate": round(sum(d["is_plasmid"] for d in rows) / len(rows), 4),
        "architecture": dict(Counter(d["architecture"] for d in rows).most_common()),
        "reaction_classes": dict(Counter(d["reaction_class"] for d in rows).most_common()),
        "chem_families": dict(Counter(d["chem_family"] for d in rows).most_common()),
        "groups": dict(Counter(d["group"] for d in rows).most_common()),
        "top_genera": Counter(d["genus"] for d in rows).most_common(8),
        "intergenic_median": round(float(np.median(inter)), 1) if inter else None,
        "intergenic_n": len(inter),
        "enriched_types": enriched[:20],
        "n_enriched": sum(1 for r in enriched if r["flagged"]),
    }


def sankey(rows, left_field, right_field, min_flow=10):
    """Sankey baglantilari: sol kategoriden saga akis.

    Kucuk akislar "other" altinda toplanir; bir Sankey'de otuz ince serit
    okunmuyor ve her birini ayri cizmek bilgi degil gurultu ekliyor.
    """
    flows = Counter((d[left_field], d[right_field]) for d in rows)
    kept = {k: v for k, v in flows.items() if v >= min_flow}
    dropped = sum(v for k, v in flows.items() if v < min_flow)
    lefts = sorted({k[0] for k in kept})
    rights = sorted({k[1] for k in kept})
    nodes = [f"{x}" for x in lefts] + [f"{x}" for x in rights]
    index = {}
    for i, name in enumerate(lefts):
        index[("L", name)] = i
    for j, name in enumerate(rights):
        index[("R", name)] = len(lefts) + j
    links = [{"source": index[("L", a)], "target": index[("R", b)], "value": v}
             for (a, b), v in sorted(kept.items(), key=lambda kv: -kv[1])]
    return {
        "nodes": nodes,
        "n_left": len(lefts),
        "links": links,
        "min_flow": min_flow,
        "entries_shown": sum(l["value"] for l in links),
        "entries_dropped_as_small": dropped,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--chemistry", default="chemistry.csv")
    ap.add_argument("--out", default="analysis_out/control_elements.json")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    data = load(con, args.chemistry)

    with_reg = [d for d in data if d["regulator_family"]]
    reg_counts = Counter(d["regulator_family"] for d in with_reg)
    regulators = []
    for fam, n in reg_counts.most_common():
        if n < MIN_ENTRIES:
            continue
        regulators.append(family_profile(
            [d for d in with_reg if d["regulator_family"] == fam], with_reg,
            fam, "regulator"))

    # IS aileleri: bir giris birden fazla IS ailesi tasiyabilir, bu yuzden
    # "aileye ait girisler" kumeleri ortusur ve bu ciktida yazilir.
    is_rows = defaultdict(list)
    for d in data:
        for fam in d["is_families"]:
            is_rows[fam].append(d)
    with_is = [d for d in data if d["is_families"]]
    is_families = []
    for fam, rows in sorted(is_rows.items(), key=lambda kv: -len(kv[1])):
        if len(rows) < MIN_ENTRIES:
            continue
        is_families.append(family_profile(rows, with_is, fam, "IS family"))

    # ------- asil soru: mobil eleman ENZIME mi TURE mi bagli
    flat = []
    for d in with_is:
        for fam in d["is_families"]:
            flat.append(dict(d, is_family=fam))
    big_is = {f for f, rows in is_rows.items() if len(rows) >= MIN_ENTRIES}
    flat = [d for d in flat if d["is_family"] in big_is]
    type_counts = Counter(d["cluster"] for d in flat)
    genus_counts = Counter(d["genus"] for d in flat)
    by_type = [d for d in flat if type_counts[d["cluster"]] >= MIN_TYPE_ENTRIES]
    by_genus = [d for d in flat if genus_counts[d["genus"]] >= MIN_TYPE_ENTRIES]

    specificity = {
        "question": ("Is a mobile element family tied to the ENZYME it sits next "
                     "to, or to the HOST it lives in?"),
        "how": ("the same insertion-sequence families are cross-tabulated twice, "
                "once against the enzyme type and once against the bacterial "
                "genus. The two tables have different shapes and Cramer's V rises "
                "on its own with the number of rows, so neither value is "
                "interpretable alone. Each is therefore measured against a null "
                "built by shuffling its own labels, and only the excess over that "
                "null is compared."),
        "by_enzyme_type": permuted_excess(by_type, "is_family", "cluster"),
        "by_host_genus": permuted_excess(by_genus, "is_family", "genus"),
    }
    a = specificity["by_enzyme_type"]
    b = specificity["by_host_genus"]
    if a and b:
        specificity["verdict"] = (
            "host genus" if b["excess"] > a["excess"] * 1.2 else
            "enzyme type" if a["excess"] > b["excess"] * 1.2 else
            "neither clearly")
        specificity["excess_ratio"] = (
            round(a["excess"] / b["excess"], 2) if b["excess"] else None)

    out = {
        "coverage": {
            "confirmed_entries": len(data),
            "with_an_upstream_regulator": len(with_reg),
            "with_a_named_IS_family": len(with_is),
            "entries_with_unnamed_transposase_only": sum(
                1 for d in data if d["unnamed_transposase"] and not d["is_families"]),
            "entries_whose_only_hit_is_recombination_machinery": sum(
                1 for d in data
                if d["recombination_machinery"] and not d["is_families"]
                and not d["unnamed_transposase"]),
        },
        "caveats": [
            "The regulator is identified by position: the first gene upstream of "
            "the operon 5' end whose product reads as a regulator. No binding was "
            "measured. A distal or global regulator is invisible to this, and a "
            "gene that merely sits upstream can be counted.",
            "The transposon label comes from a product-name regular expression, "
            "and it catches recombination machinery that is not a mobile element "
            "at all -- Holliday junction resolvase RuvX appears 117 times. Named "
            "insertion-sequence families are therefore extracted separately from "
            "the product string and kept apart from unnamed transposases.",
            "An entry can carry more than one insertion-sequence family, so the "
            "per-family sets overlap and their sizes do not sum to the total.",
        ],
        "regulators": regulators,
        "is_families": is_families,
        "mobile_specificity": specificity,
        "sankey_regulator_to_reaction": sankey(
            [d for d in with_reg if d["reaction_class"] != "unknown"],
            "regulator_family", "reaction_class"),
        "sankey_regulator_to_chem_family": sankey(
            [d for d in with_reg if d["chem_family"] != "unknown"],
            "regulator_family", "chem_family"),
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"[yazildi] {args.out}")
    print(f"  {len(with_reg)} entries have an upstream regulator, "
          f"{len(regulators)} families reach {MIN_ENTRIES} entries")
    for r in regulators:
        print(f"    {r['family']:24s} {r['entries']:5d} entries  {r['types']:3d} types  "
              f"{r['genera']:3d} genera  {r['n_enriched']} types enriched")
    print(f"  {len(with_is)} entries carry a named IS family, "
          f"{len(is_families)} families reach {MIN_ENTRIES}")
    for r in is_families:
        print(f"    {r['family']:24s} {r['entries']:5d} entries  {r['types']:3d} types  "
              f"{r['genera']:3d} genera  {r['n_enriched']} types enriched")
    if a and b:
        print(f"  mobile specificity: enzyme excess {a['excess']:+.3f} "
              f"(V={a['observed']:.3f} null={a['null_mean']:.3f}) vs "
              f"genus excess {b['excess']:+.3f} "
              f"(V={b['observed']:.3f} null={b['null_mean']:.3f}) "
              f"-> {specificity['verdict']}")


if __name__ == "__main__":
    main()
