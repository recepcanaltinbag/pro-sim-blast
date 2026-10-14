"""Atama ne kadar guvenilir? Fikir yerine uc olcum.

NEDEN. Kullanici "hatalarin cogu HMM'den geliyor" dedi. Bu test edilebilir bir
iddia ve test edilmeden ne kabul edilmeli ne reddedilmeli. Bugune kadar
bulunan her sorunun sebebi baska cikti -- qxyA'nin uzak varyantlari genom hic
indirilmedigi icindi, dort bos tip kopya referans yuzundendi, sifir dogrulanmis
RO bayat bir ara dosya yuzundendi. Ucu de HMM degildi. Ama bu, atamanin saglam
oldugunu KANITLAMAZ; yalnizca o uc vakanin baska sebebi oldugunu gosterir.

Bir hipotez zaten olculdu ve CIKMADI: "buyuk profiller yarismayi kazanir".
Profil NSEQ'u 2 ile 900 arasinda degisiyor, ama NSEQ ile uye sayisi arasinda
Spearman rho yalnizca 0,23 (p=0,06) -- ve bu test sahte korelasyon LEHINE
yanli, cunku yaygin bir enzimin hem homologu bol hem uyesi cok olur.

Burada uc sey olculur:

  MARJ        Atamalarin yuzde kaci kil payiyla belirleniyor? Kazananla
              ikincinin bit farki kucukse atama kirilgandir: profil setinde
              kucuk bir degisiklik uyeleri savurur. Bilinen iki ornek var --
              PbaA1 tam beraberlikte (971,1'e 971,1) kaybetti, DntAc 1,5 bitle.
              Soru bunlarin istisna mi kural mi oldugu.

  KENDINI     Her kuratorlu referans kendi profiline mi atanyor? Dogrudan
  BULMA       dogruluk olcusu ve gecilememesi acik bir hatadir. Referansin
              kendisi profilin kuruldugu dizilerden biri oldugu icin bu KOLAY
              bir sinav; gecememek ciddi, gecmek ise yalnizca alt sinir.

  KAPSAMA     Kapsama barajina (0,45) takilan adaylar ne kadar? Ve takilanlar
              model uzunlugunun neresinde bitiyor? Profil katalitik loba
              ulasamiyorsa sorun esikte degil, modelin ERISIMINDE olur --
              karboksilat denetimi bunu ayri bir yoldan gostermisti.

Hicbiri "HMM iyi" ya da "kotu" demez; hangi hata turunun ne kadar oldugunu
soyler, ve duzeltmenin nereye yapilmasi gerektigini.

Cikti: analysis_out/assignment_quality.json
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

import numpy as np

# ro_filter ile AYNI esikler; buradan kopya tanim uretilmez.
try:
    from ro_filter import DEFAULT_MIN_COVERAGE
except ImportError:
    DEFAULT_MIN_COVERAGE = 0.45

# "Kil payi" esikleri. 2 bit, HMMER'da gurultu duzeyinde bir fark; 5 bit hala
# kucuk. Tek bir esik secmek yerine dagilim verilir.
MARGIN_BINS = [0.0, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0]


def parse_domtbl(path, min_coverage):
    """Aday basina (en iyi profil, skor) ve ikinci en iyi profil.

    Ayni profil bir adaya birden cok alan vurabilir; profil basina EN YUKSEK
    skor alinir, sonra profiller karsilastirilir -- ro_filter da boyle yapar.
    """
    best = defaultdict(dict)      # target -> {query: score}
    cov = defaultdict(dict)       # target -> {query: model coverage}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 23:
                continue
            target, query = p[0], p[3]
            qlen = int(p[5])
            score = float(p[13])
            hmm_from, hmm_to = int(p[15]), int(p[16])
            c = (hmm_to - hmm_from + 1) / qlen if qlen else 0.0
            if score > best[target].get(query, -1e9):
                best[target][query] = score
                cov[target][query] = c
    return best, cov


def margins(best, cov, min_coverage):
    """Kazanan ve ikinci arasindaki bit farki -- yalnizca barajı gecenler."""
    rows = []
    for target, scores in best.items():
        eligible = [(q, s) for q, s in scores.items()
                    if cov[target].get(q, 0) >= min_coverage]
        if not eligible:
            continue
        eligible.sort(key=lambda kv: -kv[1])
        winner, top = eligible[0]
        if len(eligible) == 1:
            rows.append((target, winner, top, None, None))
        else:
            runner, second = eligible[1]
            rows.append((target, winner, top, runner, top - second))
    return rows


def self_recovery(refs_fasta, library, workdir, cpu):
    """Her referans kendi profiline mi atanyor?"""
    out = os.path.join(workdir, "refs_vs_library.domtbl")
    try:
        subprocess.run(["hmmsearch", "--domtblout", out, "--noali",
                        "--cpu", str(cpu), "-E", "1e-5", library, refs_fasta],
                       check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        return {"status": "not run", "reason": str(exc)}
    best, cov = parse_domtbl(out, DEFAULT_MIN_COVERAGE)
    results, wrong = [], []
    for target, scores in best.items():
        # fasta basligi "1_101_OxoO_monooxygenase_pro" -> kume kimligi
        expected = "_".join(target.split("_")[:3])
        eligible = [(q, s) for q, s in scores.items()
                    if cov[target].get(q, 0) >= DEFAULT_MIN_COVERAGE]
        if not eligible:
            wrong.append({"reference": target, "expected": expected,
                          "got": None, "reason": "no profile cleared coverage"})
            continue
        eligible.sort(key=lambda kv: -kv[1])
        got, top = eligible[0]
        margin = (top - eligible[1][1]) if len(eligible) > 1 else None
        ok = (got == expected)
        results.append(ok)
        if not ok:
            wrong.append({"reference": target, "expected": expected, "got": got,
                          "winning_score": round(top, 1),
                          "own_score": round(scores.get(expected, float("nan")), 1)
                          if expected in scores else None,
                          "margin_bits": round(margin, 1) if margin else None})
    n = len(results)
    return {
        "status": "measured",
        "references_tested": n,
        "assigned_to_own_profile": sum(results),
        "rate": round(sum(results) / n, 4) if n else None,
        "failures": wrong,
        "what_this_is": ("an easy exam: the reference is one of the sequences "
                         "its own profile was built from, so passing is only a "
                         "lower bound. Failing it would be a clear fault."),
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domtbl", default="genomic_context/cand_dom.out",
                    help="aday x profil taramasi")
    ap.add_argument("--library", default="RieskeDB71.hmm")
    ap.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE)
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--out", default="analysis_out/assignment_quality.json")
    args = ap.parse_args()

    if not os.path.exists(args.domtbl):
        raise SystemExit(f"yok: {args.domtbl}")

    best, cov = parse_domtbl(args.domtbl, args.min_coverage)
    rows = margins(best, cov, args.min_coverage)

    # --- marj dagilimi
    contested = [(t, w, s, r, m) for t, w, s, r, m in rows if m is not None]
    only_one = len(rows) - len(contested)
    ms = sorted(m for *_, m in contested)
    bins = {}
    for lo, hi in zip(MARGIN_BINS, MARGIN_BINS[1:] + [float("inf")]):
        label = f"<{hi:g}" if lo == 0 else (f"{lo:g}-{hi:g}" if hi != float("inf")
                                            else f">={lo:g}")
        bins[label] = sum(1 for m in ms if lo <= m < hi)
    ties = sum(1 for m in ms if m == 0.0)

    # --- kapsama
    all_targets = len(best)
    passed = len(rows)
    failed_cov = all_targets - passed
    best_cov = [max(c.values()) for c in cov.values() if c]
    margin_block = {
        "question": ("how many assignments are decided by a hair? A small gap "
                     "between the winner and the runner-up means the call is "
                     "fragile: a small change to the profile set moves members."),
        "candidates_with_a_contest": len(contested),
        "candidates_with_only_one_eligible_profile": only_one,
        "exact_ties": ties,
        "margin_bits": {
            "median": round(float(np.median(ms)), 2) if ms else None,
            "p10": round(float(np.percentile(ms, 10)), 2) if ms else None,
            "p90": round(float(np.percentile(ms, 90)), 2) if ms else None,
        },
        "distribution": bins,
        "share_below_2_bits": round(sum(1 for m in ms if m < 2) / len(ms), 4)
        if ms else None,
        "share_below_5_bits": round(sum(1 for m in ms if m < 5) / len(ms), 4)
        if ms else None,
    }

    # hangi tip ciftleri birbirinden kil payiyla ayriliyor
    close_pairs = Counter()
    for t, w, s, r, m in contested:
        if m is not None and m < 5:
            close_pairs[tuple(sorted((w, r)))] += 1
    margin_block["closest_type_pairs"] = [
        {"types": list(k), "candidates_decided_by_under_5_bits": v}
        for k, v in close_pairs.most_common(12)]

    with tempfile.TemporaryDirectory(prefix="assignq_") as work:
        recovery = self_recovery(args.refs, args.library, work, args.cpu)

    out = {
        "what_was_measured": __doc__.strip().split("\n\n")[0],
        "parameters": {"min_coverage": args.min_coverage,
                       "library": args.library,
                       "domtbl": args.domtbl},
        "coverage": {
            "candidates_with_any_hit": all_targets,
            "cleared_the_coverage_gate": passed,
            "failed_the_coverage_gate": failed_cov,
            "share_failing": round(failed_cov / all_targets, 4) if all_targets else None,
            "best_model_coverage": {
                "median": round(float(np.median(best_cov)), 3) if best_cov else None,
                "p10": round(float(np.percentile(best_cov, 10)), 3) if best_cov else None,
            },
            "note": ("a candidate failing on coverage is not necessarily a "
                     "non-RO: it can also mean the profile does not reach far "
                     "enough into the protein. Those are different faults and "
                     "this number alone does not separate them."),
        },
        "decision_margins": margin_block,
        "self_recovery": recovery,
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    print(f"[yazildi] {args.out}")
    print()
    print(f"adaylar                     : {all_targets:,}")
    print(f"  kapsama barajini gecen    : {passed:,}")
    print(f"  takilan                   : {failed_cov:,} "
          f"({100 * failed_cov / all_targets:.1f} %)")
    print()
    print(f"yarisma olan aday           : {len(contested):,}")
    print(f"  tam beraberlik            : {ties:,}")
    print(f"  marj medyan               : {margin_block['margin_bits']['median']} bit")
    if margin_block["share_below_2_bits"] is not None:
        print(f"  2 bitten az farkla        : "
              f"{100 * margin_block['share_below_2_bits']:.1f} %")
        print(f"  5 bitten az farkla        : "
              f"{100 * margin_block['share_below_5_bits']:.1f} %")
    print()
    if recovery.get("status") == "measured":
        print(f"kendini bulma               : "
              f"{recovery['assigned_to_own_profile']}/{recovery['references_tested']} "
              f"({100 * (recovery['rate'] or 0):.1f} %)")
        for w in recovery["failures"][:8]:
            print(f"    {w['reference'][:34]:34s} beklenen {w['expected']:16s} "
                  f"-> {w.get('got')}")
    else:
        print(f"kendini bulma               : {recovery.get('reason')}")
    if margin_block["closest_type_pairs"]:
        print()
        print("5 bitten az farkla ayrilan tip ciftleri:")
        for p in margin_block["closest_type_pairs"][:6]:
            print(f"    {p['candidates_decided_by_under_5_bits']:5d}  "
                  f"{p['types'][0]} / {p['types'][1]}")


if __name__ == "__main__":
    sys.exit(main())
