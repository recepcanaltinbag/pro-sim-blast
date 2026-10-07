"""
Uyesi olmayan tipler -- bos tip sayfasinin SEBEBINI yazar.

SORUN
  Kuratorlu referans seti 71 tip tanimliyor, ama ro tablosunda yalnizca 61
  tipin dogrulanmis uyesi var. Kalan 10 tipin sayfasi acildiginda neredeyse
  bos bir iskelet cikiyor: baslik, substrat kutusu, sonra hicbir sey. Kullanici
  bunu "bazi enzimlerin sayfalari eksik" diye bildirdi -- hakliydi, cunku sayfa
  SESSIZCE bostu. Eksik veri ile "bu tipin kendine ait uyesi yok" ayri seylerdir
  ve sayfanin hangisi oldugunu SOYLEMESI gerekiyor.

NEDEN BOSLAR
  Atama rekabetci: ro_filter.classify_target coverage barajini (>= 0.45) gecen
  modeller arasindan EN YUKSEK BIT SKORLU olani secer. Bir model hic uye
  almiyorsa bunun uc olasi sebebi var ve bu script ucunu ayirir:

    (a) reference_absorbed
        Referans proteinin KENDISI veritabaninda dogrulanmis olarak var, ama
        baska bir tipin kumesine atanmis. Iki kuratorlu ad tek modele cokmus
        demektir -- gercek bir BIRLESME.
    (b) gate_rejected
        Model coverage/E-value barajini hicbir adayda gecemiyor.
    (c) reference_not_in_genome_set
        Referansin tam dizisi taranan genom setinden hic cikmamis. Bu durumda
        model yine de yarisa giriyor ama her seferinde yakin bir kardes model
        tarafindan geriliyor; "absorbed_by" o kardes modeldir.

OLCULEN SEYLER
  - referansin tam amino asit dizisinin ro.sequence icinde birebir karsiligi
  - ro_evidence.nearest_ref: dizi kimligine gore bu referansa EN YAKIN olan
    girisler ve HMM'in onlari hangi kumeye koydugu
  - cand_dom.out uzerinde yeniden kurulan yarisma: modelin baraji kac adayda
    gectigi, ikinci geldigi adaylarda kimin kazandigi, en yuksek kendi skorunu
    aldigi adayi kimin aldigi ve skor farki
  - referans-referans diamond kimligi: en yakin kardes referans

Cikti: analysis_out/empty_types.json
"""

import argparse
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

# domtblout sutun indeksleri -- ro_filter.py ile ayni (0-tabanli).
_TARGET, _TLEN, _QUERY, _QLEN = 0, 2, 3, 5
_FULL_E, _FULL_SCORE = 6, 7
_HMM_FROM, _HMM_TO = 15, 16

# ro_filter.py'deki esikler. Burada TEKRAR yazilmiyor, oradan aliniyor:
# baraj degisirse bu teshis de degismeli.
try:
    from ro_filter import (DEFAULT_MAX_EVALUE, DEFAULT_MIN_COVERAGE,
                           DEFAULT_MIN_LENGTH, EXCLUDED_MODELS, _merge_intervals)
except ImportError:                                     # pragma: no cover
    DEFAULT_MIN_COVERAGE, DEFAULT_MAX_EVALUE, DEFAULT_MIN_LENGTH = 0.45, 1e-10, 300
    EXCLUDED_MODELS = frozenset({"2_205_IsoMO", "1_113_CdnD"})
    _merge_intervals = None

# gi|2317678|dbj|BAA21728.1| ... veya >WP_267253203.1 ... icindeki erisim no.
_ACC = re.compile(r"^[A-Z]{1,3}_?\d{5,}\.\d+$")
# PDB zincir basligi: >1NDO_1|Chains A, C, E|...
_PDB = re.compile(r"^(\d[A-Za-z0-9]{3})_\d+$")


def cluster_of(header):
    """refs71.fasta basligindan tip kimligi: ilk uc alt cizgi alani."""
    return "_".join(header.split("_")[:3])


def read_fasta(path):
    """{baslik: dizi} -- baslik ilk bosluga kadar."""
    out, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        out[name] = "".join(buf)
    return out


def accession_of(header_line):
    """Orijinal FASTA basligindan erisim numarasini cikar.

    Basliklar tek bicimde degil: bir kismi eski NCBI gi|...|db|acc| bicimi,
    bir kismi duz WP_ numarasi, biri de PDB zincir basligi. Tanimadigi bicimde
    tahmin URETMEZ, None doner -- yanlis bir erisim numarasi, hic olmamasindan
    daha kotudur.
    """
    tokens = [t.strip() for t in header_line.lstrip(">").split("|")]
    for token in tokens:
        if _ACC.match(token):
            return token
    pdb = _PDB.match(tokens[0]) if tokens else None
    if pdb:
        return "PDB " + pdb.group(1).upper()
    return None


def original_headers(directory):
    """{tip: ilk baslik satiri} -- kuratorlu tek-dizi FASTA dosyalarindan."""
    out = {}
    if not os.path.isdir(directory):
        return out
    for name in sorted(os.listdir(directory)):
        if not name.endswith((".fasta", ".fa", ".faa")):
            continue
        with open(os.path.join(directory, name)) as fh:
            for line in fh:
                if line.startswith(">"):
                    out.setdefault(cluster_of(os.path.splitext(name)[0]), line.strip())
                    break
    return out


def nearest_siblings(refs_path, threads):
    """Her referansin en yakin DIGER referanslari (diamond blastp).

    diamond yoksa bos doner: teshisin geri kalani bundan bagimsiz calisir.
    """
    if not shutil.which("diamond"):
        print("[uyari] diamond bulunamadi, referans-referans kimligi atlandi")
        return {}
    best = defaultdict(dict)
    with tempfile.TemporaryDirectory() as tmp:
        db, out = os.path.join(tmp, "refs"), os.path.join(tmp, "pairs.tsv")
        subprocess.run(["diamond", "makedb", "--in", refs_path, "-d", db, "--quiet"],
                       check=True)
        subprocess.run(["diamond", "blastp", "-q", refs_path, "-d", db, "-o", out,
                        "--quiet", "-p", str(threads), "-k", "200", "--max-hsps", "1",
                        "--ultra-sensitive", "-e", "1e-3", "--outfmt", "6",
                        "qseqid", "sseqid", "pident", "qcovhsp", "bitscore"], check=True)
        for line in open(out):
            query, subject, pid, qcov, bits = line.rstrip("\n").split("\t")
            a, b = cluster_of(query), cluster_of(subject)
            if a == b:
                continue
            row = {"type": b, "identity": round(float(pid), 1),
                   "qcov": round(float(qcov), 1), "bitscore": float(bits)}
            if row["bitscore"] > best[a].get(b, {"bitscore": -1})["bitscore"]:
                best[a][b] = row
    return {a: sorted(v.values(), key=lambda r: -r["bitscore"])
            for a, v in best.items()}


def read_domtbl(path, wanted):
    """cand_dom.out -> {aday: {model: {coverage, evalue, score, protein_length}}}.

    Yalnizca aranan modellerin bulundugu adaylar degil, HEPSI tutulmak zorunda:
    "bu modeli kim gerdi" sorusunun cevabi kazanan modelin skorunda.
    """
    envelopes = defaultdict(list)
    meta = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 19:
                continue
            model = fields[_QUERY]
            if model in EXCLUDED_MODELS:
                continue
            key = (fields[_TARGET], model)
            envelopes[key].append((int(fields[_HMM_FROM]), int(fields[_HMM_TO])))
            meta[key] = (int(fields[_QLEN]), int(fields[_TLEN]),
                         float(fields[_FULL_E]), float(fields[_FULL_SCORE]))
    by_target = defaultdict(dict)
    for (target, model), spans in envelopes.items():
        model_len, protein_len, evalue, score = meta[(target, model)]
        covered = _merge_intervals(spans) if _merge_intervals else 0
        by_target[target][model] = {
            "coverage": covered / model_len if model_len else 0.0,
            "evalue": evalue, "score": score, "protein_length": protein_len}
    print(f"[okundu] {path}: {len(by_target):,} aday, {len(envelopes):,} aday-model cifti")
    return by_target


def contest_stats(by_target, models):
    """Her bos model icin yarismayi yeniden kur.

    Donen (model basina):
        gate_passing       barajı gecen aday sayisi (tam boy)
        gate_rejected      hit alip barajı gecemeyen aday sayisi
        max_coverage       modelin herhangi bir adayda ulastigi en yuksek coverage
        runner_up          model IKINCI geldiginde kazanan modeller (sayim)
        best_contest       modelin en yuksek kendi skorunu aldigi aday ve kazanan
    """
    stats = {m: {"gate_passing": 0, "gate_rejected": 0, "max_coverage": 0.0,
                 "runner_up": Counter(), "best_contest": None, "hits": 0}
             for m in models}
    for target, hits in by_target.items():
        passing = [(model, h) for model, h in hits.items()
                   if h["coverage"] >= DEFAULT_MIN_COVERAGE
                   and h["evalue"] < DEFAULT_MAX_EVALUE
                   and h["protein_length"] >= DEFAULT_MIN_LENGTH]
        passing.sort(key=lambda p: -p[1]["score"])
        for model in models:
            hit = hits.get(model)
            if hit is None:
                continue
            row = stats[model]
            row["hits"] += 1
            row["max_coverage"] = max(row["max_coverage"], hit["coverage"])
            if not any(model == m for m, _ in passing):
                row["gate_rejected"] += 1
                continue
            row["gate_passing"] += 1
            winner, win_hit = passing[0]
            if winner == model:        # olmamasi gereken durum: kume bos degildi
                continue
            if len(passing) > 1 and passing[1][0] == model:
                row["runner_up"][winner] += 1
            best = row["best_contest"]
            if best is None or hit["score"] > best["own_score"]:
                row["best_contest"] = {
                    "candidate_id": target,
                    "own_score": round(hit["score"], 1),
                    "own_coverage": round(hit["coverage"], 4),
                    "winning_type": winner,
                    "winning_score": round(win_hit["score"], 1),
                    "score_gap": round(win_hit["score"] - hit["score"], 1)}
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/empty_types.json")
    ap.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--originals", default="Proteins_ROs_alphaSubUnits",
                    help="erisim numaralarini tasiyan orijinal tek-dizi FASTA dizini")
    ap.add_argument("--domtbl", default="genomic_context/cand_dom.out",
                    help="hmmsearch domtblout; yoksa yarisma analizi atlanir")
    ap.add_argument("--chemistry", default="chemistry.csv")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    refs = read_fasta(args.refs)
    curated = {cluster_of(h): (h, s) for h, s in refs.items()}
    print(f"[referans] {len(curated)} kuratorlu tip")

    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    counts = dict(con.execute("SELECT ro_cluster, COUNT(*) FROM ro WHERE is_confirmed=1 "
                              "GROUP BY ro_cluster"))
    empty = sorted(c for c in curated if not counts.get(c))
    print(f"[bos] {len(empty)} tipin dogrulanmis uyesi yok: {', '.join(empty)}")

    # Referans dizisinin veritabanindaki BIREBIR karsiligi. '*' ve '-' atilir,
    # cunku kuratorlu kayitlarda stop kodonu ve bosluk isareti bulunabiliyor.
    def canonical(seq):
        return (seq or "").upper().replace("*", "").replace("-", "")

    seq_index = defaultdict(list)
    for row in con.execute("SELECT candidate_id, protein_id, ro_cluster, is_confirmed, "
                           "model_coverage, hmm_score, sequence FROM ro "
                           "WHERE sequence IS NOT NULL"):
        seq_index[canonical(row["sequence"])].append(row)

    # ro_evidence.nearest_ref: dizi kimligine gore bu referansa en yakin olan
    # girisler. HMM onlari nereye koymus? Bagimsiz ikinci kanit.
    attributed = defaultdict(list)
    if con.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table' "
                   "AND name='ro_evidence'").fetchone()[0]:
        for row in con.execute("""
                SELECT e.nearest_ref, r.ro_cluster, COUNT(*) n,
                       MAX(e.ref_identity) max_identity, MAX(e.ref_qcov) max_qcov
                FROM ro_evidence e JOIN ro r USING(candidate_id)
                WHERE r.is_confirmed=1 AND e.nearest_ref IS NOT NULL
                GROUP BY 1, 2"""):
            attributed[row["nearest_ref"]].append(
                {"type": row["ro_cluster"], "entries": row["n"],
                 "max_identity": row["max_identity"], "max_qcov": row["max_qcov"]})
    for rows in attributed.values():
        rows.sort(key=lambda r: (-(r["max_identity"] or 0), -r["entries"]))

    chem = {}
    if os.path.exists(args.chemistry):
        import csv as _csv
        with open(args.chemistry) as fh:
            for row in _csv.DictReader(fh):
                chem[row["cluster"]] = row

    headers = original_headers(args.originals)
    siblings = nearest_siblings(args.refs, args.threads)
    contests = {}
    if os.path.exists(args.domtbl):
        contests = contest_stats(read_domtbl(args.domtbl, set(empty)), set(empty))
    else:
        print(f"[uyari] {args.domtbl} yok, yarisma analizi atlandi")

    records = []
    for cluster in empty:
        header, sequence = curated[cluster]
        canon = canonical(sequence)
        accession = accession_of(headers.get(cluster, ""))
        sibling = (siblings.get(cluster) or [None])[0]
        stat = contests.get(cluster)
        near = attributed.get(cluster) or []

        # 1. kanit: referansin tam dizisi veritabaninda var mi, nerede?
        exact = [r for r in seq_index.get(canon, []) if r["is_confirmed"]]
        exact_types = Counter(r["ro_cluster"] for r in exact)

        evidence = {
            "reference_length": len(canon),
            "reference_in_database": bool(exact),
            "reference_entries": len(exact),
            "reference_assigned_to": [{"type": t, "entries": n}
                                      for t, n in exact_types.most_common()],
            "reference_protein_id": exact[0]["protein_id"] if exact else None,
            "nearest_reference_sibling": sibling,
            "entries_nearest_to_this_reference": near[:4],
            "entries_nearest_to_this_reference_total": sum(r["entries"] for r in near),
        }
        if stat:
            evidence.update({
                "candidates_with_a_hit": stat["hits"],
                "candidates_passing_coverage_gate": stat["gate_passing"],
                "candidates_failing_coverage_gate": stat["gate_rejected"],
                "max_model_coverage_reached": round(stat["max_coverage"], 4),
                "coverage_gate": DEFAULT_MIN_COVERAGE,
                "runner_up_to": [{"type": t, "candidates": n}
                                 for t, n in stat["runner_up"].most_common(4)],
                "best_contest": stat["best_contest"],
            })

        # Siniflandirma. Barajı hic gecemediyse sebep barajdir; referansin
        # kendisi baska bir kumede dogrulanmissa gercek bir birlesme vardir;
        # geri kalanda referans genom setinden hic cikmamis ama model yine de
        # her yarismayi kaybediyor.
        if stat and stat["gate_passing"] == 0:
            case = "gate_rejected"
        elif exact:
            case = "reference_absorbed"
        else:
            case = "reference_not_in_genome_set"

        # absorbed_by sirasi: once referansin KENDI dustugu kume, sonra model
        # ikinci geldiginde en sik kazanan, sonra modelin en iyi adayini alan.
        absorbed_by = None
        if exact_types:
            absorbed_by = exact_types.most_common(1)[0][0]
        elif stat and stat["runner_up"]:
            absorbed_by = stat["runner_up"].most_common(1)[0][0]
        elif stat and stat["best_contest"]:
            absorbed_by = stat["best_contest"]["winning_type"]
        elif near:
            absorbed_by = near[0]["type"]

        gene = cluster.split("_", 2)[-1]
        absorbed_gene = absorbed_by.split("_", 2)[-1] if absorbed_by else None
        if case == "reference_absorbed":
            note = (f"The curated reference protein of {gene} is present in the database"
                    f" ({evidence['reference_protein_id'] or 'unnamed record'},"
                    f" {evidence['reference_entries']} genomic cop"
                    f"{'y' if evidence['reference_entries'] == 1 else 'ies'}), but every copy"
                    f" scores higher on the {absorbed_gene} profile, so the assignment step"
                    f" filed it under {absorbed_by}. Two curated names have collapsed onto one"
                    f" model: this is a merge, not missing data.")
        elif case == "gate_rejected":
            note = (f"No candidate protein aligned to the {gene} profile over at least"
                    f" {int(100 * DEFAULT_MIN_COVERAGE)} % of its length, so nothing reached"
                    f" the assignment step at all.")
        else:
            note = (f"The exact reference sequence of {gene} never came out of the scanned"
                    f" genome set. Its profile still competes"
                    + (f" on {evidence.get('candidates_passing_coverage_gate', 0):,}"
                       " candidates that clear the coverage gate" if stat else "")
                    + (f", but it is outscored every time, most often by {absorbed_gene}"
                       if absorbed_by else ", but it is outscored every time")
                    + ". The proteins this type would describe are therefore counted under"
                    + (f" {absorbed_by}." if absorbed_by else " a neighbouring type."))
        if sibling and sibling["identity"] >= 95:
            note += (f" The two references are {sibling['identity']:.1f} % identical, so the"
                     f" split between the names is a naming decision rather than a"
                     f" biological one.")

        records.append({
            "cluster": cluster,
            "gene": gene,
            "reference_accession": accession,
            "reference_header": headers.get(cluster),
            "substrate": (chem.get(cluster) or {}).get("substrate_en"),
            "case": case,
            "absorbed_by": absorbed_by,
            "evidence": evidence,
            "note": note,
        })

    out = {
        "generated": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
        "method": {
            "question": "why a curated enzyme type has no members of its own",
            "assignment": "a confirmed RO alpha subunit is filed under the reference profile "
                          "with the HIGHEST BIT SCORE among the profiles it covers at or above "
                          f"{DEFAULT_MIN_COVERAGE:.2f} of their length; assignment is therefore "
                          "competitive and a profile with a near-identical sibling can win "
                          "nothing at all",
            "cases": {
                "reference_absorbed": "the curated reference protein itself is in the database "
                                      "as a confirmed entry, but filed under another type -- a "
                                      "genuine merge of two curated names onto one model",
                "gate_rejected": "the profile never cleared the coverage / E-value gate, so no "
                                 "candidate reached the assignment step",
                "reference_not_in_genome_set": "the exact reference sequence was never retrieved "
                                               "from the scanned genomes; the profile still "
                                               "competes and loses every contest to a sibling",
            },
            "not_fixed_here": "this module explains the empty pages; removing the emptiness "
                              "means merging or renaming references and rebuilding the profile "
                              "library, which is a nomenclature decision",
        },
        "totals": {
            "curated_types": len(curated),
            "types_with_members": len([c for c in curated if counts.get(c)]),
            "types_without_members": len(empty),
            "by_case": dict(Counter(r["case"] for r in records)),
        },
        "types": records,
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)

    print()
    for record in records:
        print(f"[{record['case']:28s}] {record['cluster']:14s} "
              f"ref={record['reference_accession'] or '-':16s} "
              f"-> {record['absorbed_by'] or '(belirsiz)'}")
    print(f"\n[durum] {out['totals']['by_case']}")
    print(f"[yazildi] {args.out}")
    con.close()


if __name__ == "__main__":
    sys.exit(main())
