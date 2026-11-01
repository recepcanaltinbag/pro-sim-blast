"""RO gruplari agacta GERCEKTEN birer dal mi -- olculur, varsayilmaz.

NEDEN BU DOSYA VAR
  Filogeni sayfasi su cumleyi yayimliyordu: "Deep nodes separating the colour
  blocks correspond to the historical RO groups, which this tree recovers from
  sequence alone." Bu bir SONUC beyaniydi ve hicbir yerde olculmemisti.
  2026-11-01'de kullanici "bizim olcumumuz bunu desteklemiyor" dedi; olcum
  kullaniciyi dogruladi.

  Olculen: her grup icin (a) grubun TUM yapraklarini iceren en kucuk dal
  (MRCA) ve o dalin safligi, (b) yalnizca o gruptan yaprak iceren en buyuk
  dal. Grup bir dal ise MRCA safligi %100'dur.

  Sonuc: iddia yalnizca grup 5 ve 6 icin dogru, ve yalnizca guvenle atanmis
  diziler arasinda. Grup 1-4 hicbir kesimde dal olusturmuyor; grup 2, 3, 4 ve
  5 icin tum yapraklar hesaba katildiginda "hepsini iceren en kucuk dal"
  AGACIN TAMAMI. Referanslar bile ayrik: grup 1'in 17 referansinin en buyuk
  tek dali 6 tanesini tutuyor, grup 3'un 19'unun 7'sini, grup 4'un 11'inin
  3'unu.

  Sayilar artik buradan okunuyor, duzyazidan degil. check_group_monophyly()
  sayfanin beyani ile bu olcumun ayrilmasini engelliyor.

NEDEN UC KESIM
  "Grup bir dal mi" sorusunun cevabi hangi dizileri saydigina bagli, ve bu
  bagimlilik SONUCU DEGISTIRIYOR -- o yuzden gizlenemez. `all_tips` cizilen
  agacin kendisidir. `references_only` kuratorlu 76 referanstir, yani
  gruplandirmanin tanimlandigi diziler. `confident` ise baskin kaniti
  novel/distant OLMAYAN yapraklar: uzak atamalar agacin her yerine dagildigi
  icin, grup yapisi ancak burada gorulebilir.

Kullanim:
    python3 group_monophyly.py --db roar.sqlite --out-dir analysis_out
"""

import argparse
import collections
import json
import os
import sqlite3
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "webapp"))
import atlas                                           # noqa: E402

FAR_TIERS = ("novel", "distant")


def group_of(name):
    """Yaprak adindan RO grubu: 'REF|5_505_qxyA' ve '5_505_qxyA#267' -> '5'."""
    base = name[4:] if name.startswith("REF|") else name
    return base.split("_")[0]


def leaf_sets(root):
    """Her dugum icin altindaki yaprak adlari. Donen: [frozenset, ...]."""
    out = []

    def walk(node):
        if not node.children:
            names = frozenset([node.name])
        else:
            acc = set()
            for child in node.children:
                acc |= walk(child)
            names = frozenset(acc)
        out.append(names)
        return names

    walk(root)
    return out


def measure(clades, pool):
    """Havuzdaki her grup icin MRCA safligi ve en buyuk saf dal.

    MRCA: grubun butun yapraklarini iceren en kucuk dal. Safligi
    |grup| / |MRCA| -- %100 ise grup bir daldir. En buyuk saf dal: yalnizca o
    gruptan yaprak iceren en genis dal; grup dagilmissa bu sayi uyeliginin
    kucuk bir kesridir ve dagilmanin OLCUSUDUR.
    """
    pool = set(pool)
    by_group = collections.defaultdict(set)
    for name in pool:
        by_group[group_of(name)].add(name)

    restricted = [(c & pool) for c in clades]
    restricted = [c for c in restricted if c]

    out = {}
    for group, mine in sorted(by_group.items()):
        mrca = pool
        for clade in restricted:
            if mine <= clade and len(clade) < len(mrca):
                mrca = clade
        largest_pure = 0
        for clade in restricted:
            if clade <= mine:
                largest_pure = max(largest_pure, len(clade))
        out[group] = {
            "n_leaves": len(mine),
            "mrca_size": len(mrca),
            "mrca_purity": round(len(mine) / len(mrca), 3),
            "largest_pure_clade": largest_pure,
            "largest_pure_fraction": round(largest_pure / len(mine), 3),
            "is_single_clade": len(mrca) == len(mine),
            "mrca_is_whole_pool": len(mrca) == len(pool),
        }
    return out


def dominant_tiers(con):
    """Yaprak -> uyelerinin BASKIN kanit katmani."""
    tier = dict(con.execute("SELECT candidate_id, tier FROM ro_evidence"))
    counts = collections.defaultdict(collections.Counter)
    for leaf_id, cid in con.execute("SELECT leaf_id, candidate_id FROM ro_leaf"):
        counts[leaf_id][tier.get(cid, "unknown")] += 1
    return {lid: c.most_common(1)[0][0] for lid, c in counts.items() if c}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out-dir", default="analysis_out")
    args = ap.parse_args()

    tree_path = os.path.join(args.out_dir, "tree_all.nwk")
    if not os.path.exists(tree_path):
        sys.exit(f"{tree_path} yok -- once build_phylogeny.py kosmali")

    root = atlas.parse_newick(open(tree_path, encoding="utf-8").read())
    clades = leaf_sets(root)
    all_tips = sorted(clades[-1])

    con = sqlite3.connect(args.db)
    dom = dominant_tiers(con)
    far = {lid for lid, t in dom.items() if t in FAR_TIERS}

    cuts = {
        "all_tips": all_tips,
        "references_only": [n for n in all_tips if n.startswith("REF|")],
        "confident": [n for n in all_tips
                      if n.startswith("REF|") or n not in far],
    }

    doc = {
        "tree": "tree_all.nwk",
        "n_tips": len(all_tips),
        "cut_definitions": {
            "all_tips": "every tip of the published tree",
            "references_only": "the curated references only, i.e. the "
                               "sequences on which the group scheme is defined",
            "confident": "references plus every variant leaf whose dominant "
                         "evidence tier is not novel or distant",
        },
        "far_tiers_excluded_from_confident": list(FAR_TIERS),
        "cuts": {},
    }
    for name, pool in cuts.items():
        doc["cuts"][name] = {"n_pool": len(pool), "groups": measure(clades, pool)}

    # Sayfanin tasiyacagi ozet: hangi gruplar dal, hangileri degil.
    headline = {}
    for name, block in doc["cuts"].items():
        single = sorted(g for g, v in block["groups"].items()
                        if v["is_single_clade"])
        scattered = sorted(g for g, v in block["groups"].items()
                           if not v["is_single_clade"])
        headline[name] = {
            "single_clade_groups": single,
            "scattered_groups": scattered,
            # Tek referansli grup icin soru anlamsiz: n=1 her zaman bir daldir.
            "single_clade_groups_nontrivial": sorted(
                g for g in single if block["groups"][g]["n_leaves"] > 1),
        }
    doc["headline"] = headline

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, "group_monophyly.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)

    print(f"[yazildi] {out_path}")
    for name, block in doc["cuts"].items():
        hl = headline[name]
        print(f"  {name:16s} havuz {block['n_pool']:5d} | tek dal: "
              f"{','.join(hl['single_clade_groups_nontrivial']) or '-'} | "
              f"dagilmis: {','.join(hl['scattered_groups']) or '-'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
