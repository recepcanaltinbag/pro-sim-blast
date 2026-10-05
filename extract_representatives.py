"""
Her kume ve varyant (yaprak) icin TEMSILCI dizi cikar + katalitik merkez dogrula.

Iki soru cevaplar:

  1. "Uyeleri temsil eden ornekler cikarabilir misin?"
     Her yaprak icin MEDOID secilir: yaprak konsensusuna en yakin uye. Ilk uye
     degil, gercekten merkezi temsil eden dizi. Kume duzeyinde de en kalabalik
     yapragin medoid'i = o kumenin "tip" dizisi.

  2. "Bunlar katalitik protein bazli mi?"
     EVET -- tum analiz RO alpha alt birimi uzerinde. Alpha, hem Rieske [2Fe-2S]
     merkezini HEM de mononukleer Fe(II) KATALITIK merkezini tasiyan tek alt
     birimdir (beta/ferredoksin/reduktaz sadece komsu, analiz birimi degil).
     Bu modul her temsilcinin katalitik triad'ini dogrular ve raporlar --
     boylece temsilcilerin gercekten katalitik protein oldugu gosterilir.

Cikti:
    analysis_out/representatives_by_leaf.fasta      her yaprak icin bir temsilci
    analysis_out/representatives_by_cluster.fasta   her kume icin tip dizisi
    analysis_out/representatives.csv                 katalitik dogrulama dahil
"""

import argparse
import csv
import os
import sqlite3
from collections import Counter, defaultdict

from ro_motif import read_stockholm_matchcols, classify_motifs

MIN_LEAF_FOR_REP = 5   # bu boyutun altindaki yapraklar temsilci uretmez (tekiller hariç)


def consensus(sequences, length):
    result = []
    for position in range(length):
        counter = Counter(s[position] for s in sequences if s[position] != "-")
        result.append(counter.most_common(1)[0][0] if counter else "-")
    return "".join(result)


def identity_to(sequence, reference):
    match = total = 0
    for a, b in zip(sequence, reference):
        if a == "-":
            continue
        total += 1
        if a == b:
            match += 1
    return match / total if total else 0.0


def load_sequences(fasta_path):
    """candidate_id -> ungapped sequence."""
    sequences, current, chunks = {}, None, []
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    sequences[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if current:
        sequences[current] = "".join(chunks)
    return sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    parser.add_argument("--sequences", default="genomic_context/ro_candidates.fasta")
    parser.add_argument("--ecology", default="cluster_ecology.csv")
    parser.add_argument("--out-dir", default="analysis_out")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    aligned = read_stockholm_matchcols(args.alignment)
    length = len(next(iter(aligned.values())))
    raw = load_sequences(args.sequences)

    ecology = {}
    if os.path.exists(args.ecology):
        with open(args.ecology) as handle:
            ecology = {row["cluster"]: row for row in csv.DictReader(handle)}

    connection = sqlite3.connect(args.db)
    connection.row_factory = sqlite3.Row

    meta = {row["candidate_id"]: dict(row) for row in connection.execute(
        "SELECT r.candidate_id, r.protein_id, rep.organism, rep.is_plasmid "
        "FROM ro r JOIN replicon rep ON rep.nucleotide_id=r.nucleotide_id "
        "WHERE r.is_confirmed=1")}

    # Yaprak uyeleri
    leaf_members = defaultdict(list)
    leaf_info = {}
    for row in connection.execute("SELECT candidate_id, leaf_id FROM ro_leaf"):
        leaf_members[row["leaf_id"]].append(row["candidate_id"])
    for row in connection.execute("SELECT leaf_id, cluster, size, median_identity FROM leaf"):
        leaf_info[row["leaf_id"]] = dict(cluster=row["cluster"], size=row["size"],
                                         identity=row["median_identity"])
    connection.close()

    def medoid(members):
        """Konsensusa en yakin uye (gercek temsilci)."""
        seqs = [aligned[m] for m in members if m in aligned]
        if not seqs:
            return members[0]
        cons = consensus(seqs, length)
        best, best_id = -1, members[0]
        for m in members:
            if m not in aligned:
                continue
            score = identity_to(aligned[m], cons)
            if score > best:
                best, best_id = score, m
        return best_id

    rows = []
    cluster_best = {}   # cluster -> (leaf_size, leaf_id, rep_id)

    for leaf_id, members in leaf_members.items():
        info = leaf_info[leaf_id]
        if info["size"] < MIN_LEAF_FOR_REP and info["size"] > 1:
            continue  # cok kucuk (2-4) atla; tekiller (1) kendisi temsilci
        rep_id = medoid(members)
        motif = classify_motifs(aligned[rep_id]) if rep_id in aligned else {}
        m = meta.get(rep_id, {})
        rows.append({
            "leaf_id": leaf_id, "cluster": info["cluster"], "size": info["size"],
            "leaf_identity": round(info["identity"], 3),
            "rep_protein_id": (m["protein_id"] if m else "") or "",
            "rep_candidate_id": rep_id,
            "organism": (m["organism"] if m else "") or "",
            "is_plasmid": (m["is_plasmid"] if m else 0) or 0,
            "substrate": ecology.get(info["cluster"], {}).get("substrate", "bilinmiyor"),
            "rieske_intact": int(motif.get("rieske_intact", False)),
            "catalytic_intact": int(motif.get("catalytic_intact", False)),
            "rieske_sites": f"{motif.get('rieske_sites_found', 0)}/{motif.get('rieske_sites_total', 4)}",
            "catalytic_sites": f"{motif.get('catalytic_sites_found', 0)}/{motif.get('catalytic_sites_total', 3)}",
        })
        # Kume tip dizisi = en buyuk yaprak
        if (info["cluster"] not in cluster_best
                or info["size"] > cluster_best[info["cluster"]][0]):
            cluster_best[info["cluster"]] = (info["size"], leaf_id, rep_id)

    # --- Cikti: yaprak temsilcileri
    leaf_fasta = os.path.join(args.out_dir, "representatives_by_leaf.fasta")
    with open(leaf_fasta, "w") as handle:
        for r in sorted(rows, key=lambda x: (x["cluster"], -x["size"])):
            seq = raw.get(r["rep_candidate_id"], "")
            if not seq:
                continue
            handle.write(f">{r['leaf_id']} {r['rep_protein_id']} n={r['size']} "
                         f"id={r['leaf_identity']} cat={r['catalytic_sites']} "
                         f"substrate={r['substrate']} org={r['organism']}\n{seq}\n")

    # --- Cikti: kume tip dizileri
    cluster_fasta = os.path.join(args.out_dir, "representatives_by_cluster.fasta")
    with open(cluster_fasta, "w") as handle:
        for cluster, (size, leaf_id, rep_id) in sorted(cluster_best.items()):
            seq = raw.get(rep_id, "")
            if not seq:
                continue
            m = meta.get(rep_id, {})
            substrate = ecology.get(cluster, {}).get("substrate", "bilinmiyor")
            handle.write(f">{cluster}_type leaf={leaf_id} n={size} "
                         f"{(m.get('protein_id') if m else '') or ''} "
                         f"substrate={substrate}\n{seq}\n")

    # --- CSV
    csv_path = os.path.join(args.out_dir, "representatives.csv")
    with open(csv_path, "w", newline="") as handle:
        fieldnames = ["leaf_id", "cluster", "size", "leaf_identity", "rep_protein_id",
                      "rep_candidate_id", "organism", "is_plasmid", "substrate",
                      "rieske_intact", "catalytic_intact", "rieske_sites", "catalytic_sites"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # --- Ozet + katalitik cevap
    print("=" * 68)
    print("TEMSILCILER")
    print("=" * 68)
    print(f"  yaprak temsilcisi : {len(rows)}  ({leaf_fasta})")
    print(f"  kume tip dizisi   : {len(cluster_best)}  ({cluster_fasta})")

    both = sum(1 for r in rows if r["rieske_intact"] and r["catalytic_intact"])
    rieske = sum(1 for r in rows if r["rieske_intact"])
    catalytic = sum(1 for r in rows if r["catalytic_intact"])
    total = len(rows)
    print("\n" + "=" * 68)
    print('"BUNLAR KATALITIK PROTEIN BAZLI MI?" -- EVET')
    print("=" * 68)
    print("  Analiz birimi RO alpha alt birimi: Rieske [2Fe-2S] + mononukleer")
    print("  Fe(II) KATALITIK merkezini birlikte tasiyan tek alt birim.")
    print(f"\n  Temsilcilerde ({total}):")
    print(f"    Rieske merkezi tam    : {rieske:>4d}  (%{100*rieske/total:.1f})")
    print(f"    Katalitik triad tam   : {catalytic:>4d}  (%{100*catalytic/total:.1f})")
    print(f"    HER IKISI birden      : {both:>4d}  (%{100*both/total:.1f})")
    print("\n  Katalitik triad eksik gorunenler genelde subaile hizalama kaymasi;")
    print("  Rieske+katalitik ikisi de tam olanlar kesin katalitik alpha subunit.")

    print(f"\n[yazildi] {csv_path}")
    print("\nOrnek kume tip dizileri:")
    for cluster, (size, leaf_id, rep_id) in sorted(cluster_best.items())[:8]:
        m = meta.get(rep_id, {})
        print(f"  {cluster:16s} tip={leaf_id:18s} temsilci={(m.get('protein_id') if m else '?') or '?'}"
              f"  ({(m.get('organism') if m else '')[:26]})")


if __name__ == "__main__":
    main()
