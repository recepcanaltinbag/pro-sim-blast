"""
Referans setinin kendi icindeki fazlaligi olcer ve ATAMAYA ne yaptigini gosterir.

NEDEN BU SCRIPT VAR. 71 kuratorlu referansin yalnizca 68'i farkli dizi. Birebir
ayni olan bir cift, HMM kutuphanesine IKI profil olarak giriyor ve bir uyenin
hangisine atandigi artik biyoloji degil, skorlardaki kucuk sayisal farklar
belirliyor. Sonuc gorunur: NarAa 30 uye toplarken birebir ikizi NDO(3_314) 2
uye topluyor, NahAc ile ikizi NDO(3_315) ise hicbir uye toplamiyor -- ucuncu
bir profil hepsini aliyor. Bu, okuyucuya bir tip sayfasinda "0 uye" olarak
gorunuyor ve hicbir aciklamasi yok.

Bu modul durumu DUZELTMEZ; duzeltme referans setini degistirip pipeline'i
bastan kosturmayi gerektiriyor ve hangi adin tutulacagi bir isimlendirme
karari. Modul durumu OLCER ve sayfalarda gorunur kilar, cunku bir sayinin
neden 0 oldugunu soyleyememek daha kotu.

Cikti:
    analysis_out/reference_redundancy.json
"""

import argparse
import csv
import json
import os
import sqlite3
from collections import defaultdict

NEAR_IDENTICAL = 99.0      # bu kimligin uzerindeki ciftler de bildirilir


def read_fasta(path, cluster_ids=()):
    """id -> dizi. Id'ler TIP kimligine indirilir.

    FASTA basliklari gen islevini de tasiyor
    ("3_309_NahAc_dioxygenase_pro") oysa veritabanindaki tip kimligi
    "3_309_NahAc". Eslesmeyi ad kesmekle yapmak kirilgan, cunku bazi gen
    adlari rakam ve harf karisik ("KshA15", "ROCH34"); bu yuzden kuratorlu
    tip listesiyle ONEK eslesemesi yapilir ve eslesmeyen baslik oldugu gibi
    birakilir, boylece sessizce kaybolmaz.
    """
    by_length = sorted(cluster_ids, key=len, reverse=True)

    def normalise(header):
        for cluster in by_length:
            if header == cluster or header.startswith(cluster + "_"):
                return cluster
        return header

    records, name, chunks = {}, None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    records[name] = "".join(chunks)
                name, chunks = normalise(line[1:].split()[0]), []
            elif line:
                chunks.append(line)
    if name:
        records[name] = "".join(chunks)
    return records


def member_counts(connection):
    return dict(connection.execute(
        "SELECT ro_cluster, COUNT(*) FROM ro WHERE is_confirmed=1 "
        "AND ro_cluster IS NOT NULL GROUP BY 1"))


def identical_groups(records):
    """Birebir ayni diziyi paylasan referanslar."""
    by_sequence = defaultdict(list)
    for name, sequence in records.items():
        by_sequence[sequence].append(name)
    return [sorted(names) for names in by_sequence.values() if len(names) > 1]


def contained_pairs(records):
    """Bir referansin dizisi digerinin ICINDE tam olarak geciyor mu?

    Bu bir hizalama sorusu degil, dizgi sorusu: NdmC'nin 355 kalintisi
    NdmB'nin 373 kalintisinin icinde oldugu gibi gecerse, kisa olan uzun
    olanin parcasidir ve ayri bir profil olarak durmasinin nedeni yok.
    """
    pairs = []
    items = sorted(records.items(), key=lambda kv: len(kv[1]))
    for i, (short_name, short_seq) in enumerate(items):
        for long_name, long_seq in items[i + 1:]:
            if len(short_seq) < len(long_seq) and short_seq in long_seq:
                pairs.append({"contained": short_name, "container": long_name,
                              "contained_length": len(short_seq),
                              "container_length": len(long_seq)})
    return pairs


def near_identical(pairs_path, threshold, exclude):
    """reference_pairs.csv'den yuksek kimlikli ama ayni OLMAYAN ciftler."""
    out = []
    if not os.path.exists(pairs_path):
        return out
    with open(pairs_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                identity = float(row["identity"])
            except (TypeError, ValueError):
                continue
            key = frozenset((row["ref_a"], row["ref_b"]))
            if identity >= threshold and row["ref_a"] != row["ref_b"] and key not in exclude:
                out.append({"ref_a": row["ref_a"], "ref_b": row["ref_b"],
                            "identity": identity,
                            "substrate_a": row.get("substrate_a", ""),
                            "substrate_b": row.get("substrate_b", ""),
                            "same_substrate": row.get("same_substrate_label") == "1"})
    out.sort(key=lambda r: -r["identity"])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--refs", default=os.path.join("ROs_71_Clean", "refs71.fasta"))
    parser.add_argument("--pairs", default=os.path.join("analysis_out", "reference_pairs.csv"))
    parser.add_argument("--chemistry", default="chemistry.csv")
    parser.add_argument("--out", default=os.path.join("analysis_out",
                                                      "reference_redundancy.json"))
    args = parser.parse_args()

    cluster_ids = []
    if os.path.exists(args.chemistry):
        with open(args.chemistry, newline="", encoding="utf-8") as handle:
            cluster_ids = [row["cluster"] for row in csv.DictReader(handle) if row.get("cluster")]
    records = read_fasta(args.refs, cluster_ids)
    unmatched = sorted(name for name in records if name not in set(cluster_ids))
    if unmatched:
        print(f"[uyari] {len(unmatched)} FASTA basligi tip listesinde yok: "
              f"{', '.join(unmatched[:4])}")
    connection = sqlite3.connect(args.db)
    counts = member_counts(connection)

    groups = []
    seen_pairs = set()
    for names in identical_groups(records):
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                seen_pairs.add(frozenset((a, b)))
        total = sum(counts.get(n, 0) for n in names)
        groups.append({
            "kind": "identical",
            "members": [{"id": n, "entries": counts.get(n, 0)} for n in names],
            "length": len(records[names[0]]),
            "entries_total": total,
            # Bolunmenin ne kadar dengesiz oldugu: 1,0 demek hepsi tek profile
            # gitti, 0,5 demek esit bolundu. Toplam 0 ise anlamsiz.
            "skew": round(max(counts.get(n, 0) for n in names) / total, 3) if total else None,
        })

    contained = contained_pairs(records)
    for pair in contained:
        pair["contained_entries"] = counts.get(pair["contained"], 0)
        pair["container_entries"] = counts.get(pair["container"], 0)
        seen_pairs.add(frozenset((pair["contained"], pair["container"])))

    near = near_identical(args.pairs, NEAR_IDENTICAL, seen_pairs)
    for pair in near:
        pair["entries_a"] = counts.get(pair["ref_a"], 0)
        pair["entries_b"] = counts.get(pair["ref_b"], 0)

    # Etkilenen her tip -> sayfada gosterilecek not.
    affected = {}
    for group in groups:
        names = [m["id"] for m in group["members"]]
        for name in names:
            others = [n for n in names if n != name]
            affected[name] = {"kind": "identical", "others": others,
                              "entries": counts.get(name, 0),
                              "entries_total": group["entries_total"]}
    for pair in contained:
        affected[pair["contained"]] = {
            "kind": "contained_in", "others": [pair["container"]],
            "entries": pair["contained_entries"],
            "entries_total": pair["contained_entries"] + pair["container_entries"]}
        affected.setdefault(pair["container"], {
            "kind": "contains", "others": [pair["contained"]],
            "entries": pair["container_entries"],
            "entries_total": pair["contained_entries"] + pair["container_entries"]})
    for pair in near:
        for name, other in ((pair["ref_a"], pair["ref_b"]), (pair["ref_b"], pair["ref_a"])):
            affected.setdefault(name, {
                "kind": "near_identical", "others": [other],
                "identity": pair["identity"],
                "entries": counts.get(name, 0),
                "entries_total": pair["entries_a"] + pair["entries_b"]})

    result = {
        "method": {
            "question": ("how much of the reference set is redundant, and what that does "
                         "to type assignment"),
            "identical": ("byte-identical amino acid sequences under different reference "
                          "names; each becomes its own profile, so which one recruits a "
                          "member is decided by small numerical differences in score "
                          "rather than by biology"),
            "contained": ("one reference sequence occurs exactly inside another, so the "
                          "shorter is a fragment of the longer"),
            "near_identical": f"distinct sequences at or above {NEAR_IDENTICAL} % identity",
            "not_fixed_here": ("this module measures the problem; repairing it means "
                               "changing the reference set and rebuilding the profile "
                               "library, and which name to keep is a nomenclature decision"),
            "near_identical_threshold": NEAR_IDENTICAL,
        },
        "totals": {
            "references": len(records),
            "distinct_sequences": len({s for s in records.values()}),
            "identical_groups": len(groups),
            "references_in_identical_groups": sum(len(g["members"]) for g in groups),
            "contained_pairs": len(contained),
            "near_identical_pairs": len(near),
            "affected_references": len(affected),
            "entries_under_an_ambiguous_reference": sum(
                counts.get(name, 0) for name in affected),
        },
        "identical_groups": groups,
        "contained_pairs": contained,
        "near_identical_pairs": near,
        "affected": affected,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=1)

    totals = result["totals"]
    print("=" * 72)
    print("REFERANS SETI FAZLALIGI")
    print("=" * 72)
    print(f"  referans                      {totals['references']}")
    print(f"  farkli dizi                   {totals['distinct_sequences']}")
    print(f"  birebir ayni grup             {totals['identical_groups']}")
    print(f"  bir digerinin parcasi         {totals['contained_pairs']}")
    print(f"  >= {NEAR_IDENTICAL}% kimlikli cift        {totals['near_identical_pairs']}")
    print(f"  etkilenen referans            {totals['affected_references']}")
    print(f"  bu referanslardaki giris      {totals['entries_under_an_ambiguous_reference']}")
    for group in groups:
        names = " == ".join(f"{m['id']} ({m['entries']})" for m in group["members"])
        print(f"    ayni: {names}  -> toplam {group['entries_total']}"
              f"  bolunme {group['skew']}")
    for pair in contained:
        print(f"    parca: {pair['contained']} ({pair['contained_entries']}) "
              f"icinde {pair['container']} ({pair['container_entries']})")
    for pair in near:
        print(f"    yakin: {pair['ref_a']} ({pair['entries_a']}) / "
              f"{pair['ref_b']} ({pair['entries_b']}) %{pair['identity']:.1f}"
              f"  ayni substrat: {pair['same_substrate']}")
    print(f"\n[yazildi] {args.out}")


if __name__ == "__main__":
    main()
