"""
Filogeni ve dizi benzerlik agi (SSN) -- referanslar + varyant temsilcileri.

Girdi:  ROs_71_Clean/refs71.fasta (71 deneysel RO)
        analysis_out/representatives_by_leaf.fasta (her yaprak/varyant icin 1 temsilci)
Adimlar:
  1. hmmalign (ROmotif71.hmm) -> match-state kolonlari -> %50'den fazla bosluklu
     kolonlar atilir -> FastTree -lg  => analysis_out/tree_all.nwk
     (ayni hizalamayi kume basina kesip kume agaclari: analysis_out/trees/<cluster>.nwk)
  2. diamond all-vs-all -> SSN kenarlari (kimlik >= --ssn-min, kaplama >= 70%)
     => analysis_out/ssn_edges.csv, ssn_nodes.csv
  3. kume x kume kimlik matrisi (temsilciler arasi en yuksek/medyan kimlik)
     => analysis_out/cluster_identity_matrix.csv

Agac yorumu: hizalama yalnizca Rieske+katalitik cekirdek match-state kolonlarina
dayandigi icin insert/fuzyon bolgeleri agaca girmez; dal uzunluklari cekirdek
domain evrimini olcer.
"""

import argparse
import csv
import os
import sqlite3
import subprocess
import sys
import tempfile
from collections import defaultdict

from ro_motif import read_stockholm_matchcols

MAX_GAP_FRAC = 0.5


def run(cmd, **kw):
    print("[calisiyor]", " ".join(cmd))
    subprocess.run(cmd, check=True, **kw)


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def write_phylip_free(aln, path):
    """FastTree fasta hizalama okur; isimlerde ':' sorun degil (newick'te tirnaklanir)."""
    with open(path, "w") as fh:
        for name, seq in aln.items():
            fh.write(f">{name}\n{seq}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--reps", default="analysis_out/representatives_by_leaf.fasta")
    ap.add_argument("--motif-hmm", default="ROs_71_Clean/ROmotif71.hmm")
    ap.add_argument("--out-dir", default="analysis_out")
    ap.add_argument("--ssn-min", type=float, default=30.0)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--fasttree", default="FastTree")
    args = ap.parse_args()

    refs = read_fasta(args.refs)
    reps = read_fasta(args.reps)
    # referans adi: model adi (ilk 3 parca) -> "REF|1_101_OxoO"
    refs = {"REF|" + "_".join(k.split("_")[:3]): v for k, v in refs.items()}
    allseq = {**refs, **reps}
    print(f"[bilgi] {len(refs)} referans + {len(reps)} temsilci")

    con = sqlite3.connect(args.db)
    leaf_of = dict(con.execute("SELECT candidate_id, leaf_id FROM ro_leaf"))
    cluster_of = dict(con.execute("SELECT candidate_id, ro_cluster FROM ro WHERE is_confirmed=1"))
    leaf_size = dict(con.execute("SELECT leaf_id, size FROM leaf"))
    organism = dict(con.execute(
        "SELECT r.candidate_id, p.organism FROM ro r JOIN replicon p USING(nucleotide_id)"))
    tier = dict(con.execute("SELECT candidate_id, tier FROM ro_evidence")) \
        if con.execute("SELECT name FROM sqlite_master WHERE name='ro_evidence'").fetchone() else {}
    con.close()

    os.makedirs(os.path.join(args.out_dir, "trees"), exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "all.fa")
        with open(fa, "w") as fh:
            for k, v in allseq.items():
                fh.write(f">{k}\n{v}\n")
        sto = os.path.join(tmp, "all.sto")
        run(["hmmalign", "--trim", "--amino", "-o", sto, args.motif_hmm, fa])
        aln = read_stockholm_matchcols(sto)
        # bosluklu kolonlari at
        L = len(next(iter(aln.values())))
        keep = [i for i in range(L)
                if sum(1 for s in aln.values() if s[i] == "-") / len(aln) <= MAX_GAP_FRAC]
        aln = {k: "".join(v[i] if v[i] in "ACDEFGHIKLMNPQRSTVWY-" else "-" for i in keep)
               for k, v in aln.items()}   # X/B/Z -> bosluk (FastTree sayisal hata)
        print(f"[hizalama] {len(aln)} dizi, {L} match kolonu -> {len(keep)} tutuldu")

        # --- 1. global agac
        aln_fa = os.path.join(tmp, "aln.fa")
        write_phylip_free(aln, aln_fa)
        tree_path = os.path.join(args.out_dir, "tree_all.nwk")
        # FastTree tek-duyarlikli ikili bazi hizalamalarda PairLogLk assertion'i ile
        # cokuyor; sirayla LG -> JTT -> ML'siz (NJ+ME) dene.
        for extra in (["-lg"], [], ["-noml"]):
            try:
                with open(tree_path, "w") as out:
                    run([args.fasttree, "-quiet", "-quote"] + extra + [aln_fa], stdout=out)
                print(f"[agac] model secenegi: {extra or ['JTT']}")
                break
            except subprocess.CalledProcessError:
                print(f"[uyari] FastTree {extra} basarisiz, sonraki secenek")
        else:
            raise RuntimeError("FastTree hicbir modelde calismadi")
        print(f"[yazildi] {tree_path}")

        # --- kume agaclari: kumenin yapraklari + kumenin referansi (+ en yakin 2 referans yok; sade)
        by_cluster = defaultdict(list)
        for k in reps:
            by_cluster[cluster_of.get(k, "?")].append(k)
        for cl, members in by_cluster.items():
            names = members + [n for n in refs if n == "REF|" + cl]
            if len(names) < 3:
                continue
            sub = {k: aln[k] for k in names if k in aln}
            sub_fa = os.path.join(tmp, "sub.fa")
            write_phylip_free(sub, sub_fa)
            for extra in (["-lg"], [], ["-noml"]):
                try:
                    with open(os.path.join(args.out_dir, "trees", f"{cl}.nwk"), "w") as out:
                        subprocess.run([args.fasttree, "-quiet", "-quote"] + extra + [sub_fa],
                                       check=True, stdout=out, stderr=subprocess.DEVNULL)
                    break
                except subprocess.CalledProcessError:
                    continue
        print(f"[yazildi] {len(by_cluster)} kume agaci -> {args.out_dir}/trees/")

        # --- 2. SSN (diamond all-vs-all)
        db = os.path.join(tmp, "all")
        run(["diamond", "makedb", "--in", fa, "-d", db, "--quiet"])
        hits = os.path.join(tmp, "hits.tsv")
        run(["diamond", "blastp", "-q", fa, "-d", db, "-o", hits, "--quiet", "-p", str(args.threads),
             "-k", "2000", "--max-hsps", "1", "--sensitive", "-e", "1e-5",
             "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "bitscore"])
        edges = {}
        best_pair = defaultdict(float)
        with open(hits) as fh:
            for line in fh:
                a, b, pid, qcov, bits = line.rstrip("\n").split("\t")
                if a >= b:
                    continue
                pid, qcov, bits = float(pid), float(qcov), float(bits)
                if qcov < 70:
                    continue
                key = (a, b)
                if pid > edges.get(key, (0,))[0]:
                    edges[key] = (pid, bits)
                ca = cluster_of.get(a, a[4:] if a.startswith("REF|") else "?")
                cb = cluster_of.get(b, b[4:] if b.startswith("REF|") else "?")
                if ca != cb:
                    k2 = tuple(sorted((ca, cb)))
                    best_pair[k2] = max(best_pair[k2], pid)

    nodes_path = os.path.join(args.out_dir, "ssn_nodes.csv")
    with open(nodes_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "kind", "cluster", "leaf_id", "leaf_size", "organism", "tier"])
        for k in allseq:
            if k.startswith("REF|"):
                w.writerow([k, "reference", k[4:], "", "", "", "characterized"])
            else:
                w.writerow([k, "representative", cluster_of.get(k, ""), leaf_of.get(k, ""),
                            leaf_size.get(leaf_of.get(k, ""), ""), organism.get(k, ""), tier.get(k, "")])
    edges_path = os.path.join(args.out_dir, "ssn_edges.csv")
    kept = 0
    with open(edges_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "target", "identity", "bitscore"])
        for (a, b), (pid, bits) in edges.items():
            if pid >= args.ssn_min:
                w.writerow([a, b, pid, bits]); kept += 1
    print(f"[yazildi] {nodes_path} ({len(allseq)} dugum), {edges_path} ({kept} kenar >= {args.ssn_min}%)")

    clusters = sorted({c for c in cluster_of.values()})
    mat_path = os.path.join(args.out_dir, "cluster_identity_matrix.csv")
    with open(mat_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cluster"] + clusters)
        for a in clusters:
            w.writerow([a] + [100.0 if a == b else round(best_pair.get(tuple(sorted((a, b))), 0.0), 1)
                              for b in clusters])
    print(f"[yazildi] {mat_path} (kumeler arasi en yuksek temsilci kimligi)")


if __name__ == "__main__":
    sys.exit(main())
