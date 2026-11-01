"""Dal destegi BOOTSTRAP ile -- hizalama kolonlari gercekten yeniden orneklenir.

NEDEN BU DOSYA VAR
  Yayimlanan agacin destek degerleri FastTree'nin VARSAYILAN destegidir:
  SH-benzeri YEREL destek. O test, her ic dali yalnizca kendi en yakin
  komsu yeniden duzenlemelerine karsi sinar ve hizalama kolonlarini
  YENIDEN ORNEKLEMEZ. Bootstrap ise tam bunu yapar, ve literaturun bir
  klad iddiasindan bekledigi sey budur.

  Fark somut bir karari engelliyordu. Olcum (group_monophyly.json) sunu
  soyluyor: bizim "ungrouped" dedigimiz bitki klorofil dali (PAO, CAO,
  PTC52) agacta TEK BIR DAL -- guvenle atanmis kesimde 26/26, kuratorlu
  referanslarda 3/3. Bu, "yeni bir grup onerelim" demek icin en cazip
  aday. Ama bir grup onerisi "proposed group VI (this study)" diye
  yazilacaksa, dalin ne kadar saglam oldugunu soyleyen sey DESTEKTIR ve
  SH-benzeri yerel destek o iddiayi tasimaz. Bu betik, tasiyip tasimadigini
  OLCER -- ve cevap "hayir" cikarsa o da bir sonuctur.

NASIL
  1. Secilen dizi kumesi hmmalign ile ayni profile hizalanir ve ayni kolon
     suzgeci uygulanir (build_phylogeny ile AYNI yol; farkli hizalama
     farkli agac demek olurdu).
  2. Referans agac FastTree ile cikarilir.
  3. N tekrar: kolonlar YERINE KOYARAK yeniden orneklenir, her tekrar icin
     agac cikarilir.
  4. Referans agacin her ic dali bir BOLUNMEDIR (bipartition). Destek, o
     bolunmeyi iceren tekrarlarin orani.

  Koksuz agacta bolunme, taksonlari iki kumeye ayirir ve hangisinin
  "ic" oldugu keyfidir; o yuzden her bolunme, sabit bir taksonu
  ICERMEYEN taraf olarak normallestirilir. Boylece ayni bolunme iki
  agacta ayni anahtari verir.

KUME SECIMI
  --set references : 76 kuratorlu referans. Grup semasinin tanimlandigi
                     diziler; bir grup onerisi once burada durmali.
  --set confident  : referanslar + baskin kaniti novel/distant OLMAYAN
                     yapraklar. Uzak atamalar agacin her yerine dagildigi
                     icin grup yapisi ancak burada gorulebilir.

Kullanim:
    python3 bootstrap_tree.py --set references --replicates 200
    python3 bootstrap_tree.py --set confident  --replicates 100
"""

import argparse
import collections
import json
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "webapp"))

import atlas                                           # noqa: E402
from build_phylogeny import MAX_GAP_FRAC, read_fasta    # noqa: E402
from ro_motif import read_stockholm_matchcols           # noqa: E402

FAR_TIERS = ("novel", "distant")


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def build_alignment(sequences, motif_hmm, work):
    """build_phylogeny ile AYNI hizalama yolu: hmmalign + kolon suzgeci."""
    with tempfile.TemporaryDirectory(dir=work) as tmp:
        fa = os.path.join(tmp, "in.fa")
        with open(fa, "w") as handle:
            for name, seq in sequences.items():
                handle.write(f">{name}\n{seq}\n")
        sto = os.path.join(tmp, "in.sto")
        run(["hmmalign", "--trim", "--amino", "-o", sto, motif_hmm, fa],
            stdout=subprocess.DEVNULL)
        aln = read_stockholm_matchcols(sto)
        if isinstance(aln, tuple):
            aln = aln[0]
    length = len(next(iter(aln.values())))
    keep = [i for i in range(length)
            if sum(1 for s in aln.values() if s[i] == "-") / len(aln) <= MAX_GAP_FRAC]
    clean = {k: "".join(v[i] if v[i] in "ACDEFGHIKLMNPQRSTVWY-" else "-"
                        for i in keep)
             for k, v in aln.items()}
    return clean, length, len(keep)


def fasttree(aln, path, fasttree_bin, extra=("-lg",)):
    """FastTree; cokerse build_phylogeny ile ayni yedekleme zinciri."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False,
                                     dir=os.path.dirname(path)) as handle:
        for name, seq in aln.items():
            handle.write(f">{name}\n{seq}\n")
        fa = handle.name
    try:
        for option in (list(extra), [], ["-noml"]):
            try:
                with open(path, "w") as out:
                    run([fasttree_bin, "-quiet", "-quote", "-nosupport"] + option
                        + [fa], stdout=out, stderr=subprocess.DEVNULL)
                return open(path, encoding="utf-8").read()
            except subprocess.CalledProcessError:
                continue
        raise RuntimeError("FastTree hicbir secenekte calismadi")
    finally:
        os.unlink(fa)


def splits_of(newick, pivot):
    """Agactaki ic bolunmeler. Her biri, pivot taksonu ICERMEYEN taraf.

    Yaprak sayisi 1 olan ya da butun taksonlari kapsayan bolunmeler
    atlanir: ikisi de her agacta vardir, yani hicbir sey ayirt etmezler.
    """
    root = atlas.parse_newick(newick)
    all_taxa = set()
    out = set()

    def walk(node):
        if not node.children:
            all_taxa.add(node.name)
            return frozenset([node.name])
        acc = set()
        for child in node.children:
            acc |= walk(child)
        out.add(frozenset(acc))
        return frozenset(acc)

    walk(root)
    normalised = set()
    for clade in out:
        if len(clade) < 2 or len(clade) >= len(all_taxa) - 1:
            continue
        side = clade if pivot not in clade else frozenset(all_taxa - clade)
        if 2 <= len(side) <= len(all_taxa) - 2:
            normalised.add(side)
    return normalised, all_taxa


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--reps", default="analysis_out/representatives_by_leaf.fasta")
    ap.add_argument("--motif-hmm", default="ROs_71_Clean/ROmotif71.hmm")
    ap.add_argument("--out-dir", default="analysis_out")
    ap.add_argument("--set", dest="which", default="references",
                    choices=("references", "confident"))
    ap.add_argument("--replicates", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261101)
    ap.add_argument("--fasttree", default="FastTree")
    ap.add_argument("--work", default=os.environ.get("ROAR_SCRATCH",
                                                     "/media/lin-bio/back2"))
    args = ap.parse_args()

    refs = {"REF|" + "_".join(k.split("_")[:3]): v
            for k, v in read_fasta(args.refs).items()}
    if args.which == "references":
        sequences = dict(refs)
    else:
        import sqlite3                                  # noqa: PLC0415
        con = sqlite3.connect(args.db)
        tier = dict(con.execute("SELECT candidate_id, tier FROM ro_evidence"))
        counts = collections.defaultdict(collections.Counter)
        for leaf_id, cid in con.execute("SELECT leaf_id, candidate_id FROM ro_leaf"):
            counts[leaf_id][tier.get(cid, "unknown")] += 1
        far = {lid for lid, c in counts.items()
               if c and c.most_common(1)[0][0] in FAR_TIERS}
        con.close()
        reps = read_fasta(args.reps)
        sequences = dict(refs)
        sequences.update({k: v for k, v in reps.items() if k not in far})

    os.makedirs(args.work, exist_ok=True)
    print(f"[kume] {args.which}: {len(sequences)} dizi")
    aln, n_cols_raw, n_cols = build_alignment(sequences, args.motif_hmm, args.work)
    print(f"[hizalama] {n_cols_raw} match kolonu -> {n_cols} tutuldu")

    names = list(aln)
    pivot = sorted(names)[0]
    columns = [[aln[n][i] for n in names] for i in range(n_cols)]

    with tempfile.TemporaryDirectory(dir=args.work) as tmp:
        ref_newick = fasttree(aln, os.path.join(tmp, "ref.nwk"), args.fasttree)
        ref_splits, all_taxa = splits_of(ref_newick, pivot)
        print(f"[referans agac] {len(ref_splits)} ic bolunme")

        rng = random.Random(args.seed)
        tally = collections.Counter()
        for replicate in range(1, args.replicates + 1):
            picked = [rng.randrange(n_cols) for _ in range(n_cols)]
            resampled = {}
            for index, name in enumerate(names):
                resampled[name] = "".join(columns[c][index] for c in picked)
            newick = fasttree(resampled, os.path.join(tmp, "rep.nwk"),
                              args.fasttree)
            got, _ = splits_of(newick, pivot)
            for split in ref_splits & got:
                tally[split] += 1
            if replicate % 25 == 0 or replicate == args.replicates:
                print(f"  tekrar {replicate}/{args.replicates}")

    support = {}
    for split in ref_splits:
        support[split] = tally[split] / args.replicates

    values = sorted(support.values())
    strong = sum(1 for v in values if v >= 0.70)
    very = sum(1 for v in values if v >= 0.95)

    # Agacin kendisi, dal desteklerini ic dugum adi olarak tasisin.
    def annotate(node):
        if not node.children:
            return frozenset([node.name])
        acc = set()
        for child in node.children:
            acc |= annotate(child)
        clade = frozenset(acc)
        side = clade if pivot not in clade else frozenset(all_taxa - clade)
        if side in support:
            node.name = f"{support[side]:.2f}"
        return clade

    tree = atlas.parse_newick(ref_newick)
    annotate(tree)

    def to_newick(node):
        if not node.children:
            return f"'{node.name}':{node.length:.6f}"
        inner = ",".join(to_newick(c) for c in node.children)
        return f"({inner}){node.name}:{node.length:.6f}"

    os.makedirs(args.out_dir, exist_ok=True)
    nwk_path = os.path.join(args.out_dir, f"bootstrap_{args.which}.nwk")
    with open(nwk_path, "w", encoding="utf-8") as handle:
        handle.write(to_newick(tree) + ";\n")

    # Ilgilendigimiz somut soru: hangi gruplar bootstrap destekli bir dal?
    def group_of(name):
        base = name[4:] if name.startswith("REF|") else name
        return base.split("_")[0]

    by_group = collections.defaultdict(set)
    for name in all_taxa:
        by_group[group_of(name)].add(name)
    group_support = {}
    for group, members in sorted(by_group.items()):
        if len(members) < 2:
            group_support[group] = {"n": len(members), "is_split": None,
                                    "bootstrap": None,
                                    "note": "single sequence: no branch to test"}
            continue
        side = members if pivot not in members else frozenset(all_taxa - members)
        side = frozenset(side)
        # Grup bir dal DEGILSE "ne kadar degil" sorusu kalir. En buyuk saf
        # dal, dagilmanin olcusudur ve group_monophyly.json ile ayni
        # anlamda okunur -- ama DIKKAT, o dosya YAYIMLANAN agaci (2.062 uc)
        # bir alt kumeye daraltarak olcer; burasi alt kumeden YENI bir agac
        # cikarir. Ayni soru degil, ve cevaplari farkli cikabilir.
        largest_pure = 0
        for split in ref_splits:
            for form in (split, frozenset(all_taxa - split)):
                if form and form <= members:
                    largest_pure = max(largest_pure, len(form))
        group_support[group] = {
            "n": len(members),
            "is_split": side in ref_splits,
            "bootstrap": round(support.get(side), 3) if side in support else None,
            "largest_pure_clade": largest_pure,
            "largest_pure_fraction": round(largest_pure / len(members), 3),
        }

    doc = {
        "set": args.which,
        "n_sequences": len(aln),
        "alignment_columns": n_cols,
        "replicates": args.replicates,
        "seed": args.seed,
        "support": "nonparametric bootstrap, alignment columns resampled with "
                   "replacement, tree re-inferred per replicate",
        "support_is_bootstrap": True,
        "n_internal_branches": len(ref_splits),
        "median": round(values[len(values) // 2], 3) if values else None,
        "at_least_0.70": strong,
        "at_least_0.95": very,
        "group_as_a_clade": group_support,
        "tree": os.path.basename(nwk_path),
        "scope_note": ("This tree is inferred from this sequence set ALONE. "
                       "It is not the published tree restricted to that set, "
                       "and the two can disagree: analysis_out/"
                       "group_monophyly.json measures clade structure inside "
                       "the published 2,062-tip tree, which is a different "
                       "question. Compare largest_pure_fraction between the "
                       "two rather than assuming they must agree."),
    }
    json_path = os.path.join(args.out_dir, f"bootstrap_{args.which}.json")
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(doc, handle, indent=1)

    print(f"\n[yazildi] {nwk_path}")
    print(f"[yazildi] {json_path}")
    print(f"  ic dal {len(ref_splits)} | medyan {doc['median']} | "
          f">=0,70 {strong} | >=0,95 {very}")
    for group, info in group_support.items():
        if info.get("is_split") is None:
            print(f"  grup {group}: {info['n']} dizi -- {info['note']}")
        else:
            print(f"  grup {group}: {info['n']:3d} dizi | tek dal: "
                  f"{'EVET' if info['is_split'] else 'hayir':5s} | bootstrap "
                  f"{info['bootstrap'] if info['bootstrap'] is not None else '-':>5} "
                  f"| en buyuk saf dal {info['largest_pure_clade']}/{info['n']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
