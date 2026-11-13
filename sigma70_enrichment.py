"""Sigma70 kutusu ZENGINLESMESI -- eslenmis null'a karsi, uc duzeyde.

NEDEN BOYLE KURULDU. ROADMAP 124.3 bu olcumun sinirini BASTAN yaziyor:
sigma70 promotoru -10 (TATAAT) ve -35 (TTGACA) kutularindan olusur ve bu
kutular BILGI BAKIMINDAN ZAYIFTIR -- rastgele bir intergenik bolgede de
sik bulunur. Dolayisiyla operon basina "promotor burada" demek
yayimlanabilir bir iddia DEGILDIR ve burada oyle bir sey iddia edilmiyor.
Yayimlanabilir olan ZENGINLESME: gercek yukari akis bolgelerinde kutu
ciftinin siklig, ayni bolgelerden uretilmis bir null'dan yuksek mi?

IKI NULL, cunku tek null bir varsayimdir.
  A) MONONUKLEOTID KARISTIRMA: ayni bolgenin bazlari karistirilir.
     Uzunluk ve BAZ KOMPOZISYONU (dolayisiyla GC) tam korunur, kutularin
     konumu ve ciftlenmesi bozulur. Sinirli yani: dinukleotid
     kompozisyonu korunmaz.
  B) MARKOV-1 YENIDEN URETIM: ayni bolgeden gecis olasiliklari ogrenilir
     ve ayni uzunlukta yeni dizi uretilir. Dinukleotid kompozisyonu
     BEKLENTIDE korunur, yani (A)'nin sinirini kapatir.
Iddia ancak IKISINDE DE ayni yonu gosterirse yayimlanir; ayrisirlarsa
ayrisma BULGUNUN KENDISIDIR ve oyle yazilir.

UC DUZEY, KARARI CINS VERIR (bu projenin istatistik kurali). Ayni
genomdan gelen bolgeler bagimsiz degil; giris duzeyinde bir p degeri
tekrarli ornekleme yuzunden kendiliginden kucuk cikar.

IPLIK. Bolgeler `upstream_regions.py` tarafindan operonu SUREN iplikte
yaziliyor; burada yeniden cevirme YOK.

Cikti: analysis_out/sigma70_enrichment.json
"""

import argparse
import csv
import gzip
import json
import math
import os
import random
import sqlite3
import sys
from collections import Counter, defaultdict

import ro_taxonomy

#  Kanonik kutular ve araligi. Araliik literaturdeki 15-19 bp; kutular
#  birer mismatch'e kadar kabul edilir, cunku tam uzlasma dizisi
#  gercek promotorlarda da seyrek.
BOX35, BOX10 = "TTGACA", "TATAAT"
SPACER = (15, 19)
MAX_MISMATCH = 1
SHUFFLES = 20
SEED = 20261113


def mismatches(window, consensus):
    return sum(1 for a, b in zip(window, consensus) if a != b)


def count_pairs(sequence):
    """Bolgede kac (-35, -10) cifti var? Ayni -35 birden fazla -10 ile
    eslesebilir; SAYILAN sey ciftlerin sayisi degil, CIFT TASIYAN bolge
    olup olmadigi (asagida bool'a cevriliyor), cunku uzun bolge daha cok
    cift tasir ve uzunluk null'da da ayni."""
    hits = 0
    n = len(sequence)
    for i in range(n - 6):
        if mismatches(sequence[i:i + 6], BOX35) > MAX_MISMATCH:
            continue
        for gap in range(SPACER[0], SPACER[1] + 1):
            j = i + 6 + gap
            if j + 6 > n:
                break
            if mismatches(sequence[j:j + 6], BOX10) <= MAX_MISMATCH:
                hits += 1
    return hits


def mononucleotide_shuffle(sequence, rng):
    bases = list(sequence)
    rng.shuffle(bases)
    return "".join(bases)


def markov1_resample(sequence, rng):
    """Dinukleotid gecislerini ogrenip ayni uzunlukta dizi uretir."""
    if len(sequence) < 3:
        return sequence
    transitions = defaultdict(list)
    for a, b in zip(sequence, sequence[1:]):
        transitions[a].append(b)
    out = [rng.choice(sequence)]
    for _ in range(len(sequence) - 1):
        options = transitions.get(out[-1])
        out.append(rng.choice(options) if options else rng.choice(sequence))
    return "".join(out)


def load_regions(fasta, meta_csv):
    meta = {}
    with open(meta_csv, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            meta[row["candidate_id"]] = row
    regions = []
    with gzip.open(fasta, "rt", encoding="utf-8") as handle:
        name, chunks = None, []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name and name in meta:
                    regions.append((name, "".join(chunks), meta[name]))
                name, chunks = line[1:], []
            elif line:
                chunks.append(line)
        if name and name in meta:
            regions.append((name, "".join(chunks), meta[name]))
    return regions


def levels(con):
    """candidate_id -> (dizi anahtari, cins anahtari)."""
    genus = ro_taxonomy.genus_by_candidate(con)
    rows = con.execute(
        "SELECT candidate_id, sequence, ro_cluster FROM ro "
        "WHERE is_confirmed = 1").fetchall()
    out = {}
    for cid, seq, cluster in rows:
        got = genus.get(cid)
        out[cid] = (seq or cid,
                    (cluster, got.key if got else ro_taxonomy.UNKNOWN_KEY))
    return out


def binomial_tail(observed, n, probability):
    """Iki yonlu binom p degeri, LOG UZAYINDA.

    NEDEN LOG. Ilk yazim `math.comb(n, k)` kullaniyordu ve n=7.616'da
    tasti (`OverflowError: int too large to convert to float`): binom
    katsayisi float'a sigmiyor. Olasiliklar log uzayinda toplanir, yani
    buyuk n'de de tam deger hesaplanir ve bir yaklasima gecmek
    gerekmez.
    """
    if n <= 0 or not 0 < probability < 1:
        return None
    log_p, log_q = math.log(probability), math.log1p(-probability)
    log_fact_n = math.lgamma(n + 1)

    def log_pmf(k):
        return (log_fact_n - math.lgamma(k + 1) - math.lgamma(n - k + 1)
                + k * log_p + (n - k) * log_q)

    def tail(lo, hi):
        terms = [log_pmf(k) for k in range(lo, hi + 1)]
        if not terms:
            return 0.0
        top = max(terms)
        return math.exp(top + math.log(sum(math.exp(t - top)
                                           for t in terms)))

    upper = tail(observed, n)
    lower = tail(0, observed)
    return min(1.0, 2 * min(upper, lower))


def compare(rows, null_key):
    """GOZLENEN ile BEKLENEN bolge sayisini karsilastirir.

    NEDEN ISARET TESTI DEGIL -- ve ilk yazimda o kullanildi, YANLISTI.
    Bir bolgenin gercek degeri IKILI (cift var/yok), null degeri ise 20
    kopyanin ORTALAMASI (orn. 0,05). Ikisinin farkina isaret testi
    uygulamak, cift tasimayan bir bolgeyi "null lehte" sayiyor cunku null
    20'de 1 bulmus; yani -0,05 ile +0,95 ayni agirlikta sayiliyor. Sonuc
    apacik celiskiliydi: isaret testi null'u kazandiriyor, oran
    karsilastirmasi gercegi -- ayni veriden iki zit cevap. Birim soruya
    uymuyordu.

    Dogru karsilastirma SAYILAR uzerinde: kac bolgede cift GOZLENDI, ve
    null altinda kac bolgede BEKLENIR (null olasiliklarinin toplami).
    Binom testi bunu dogrudan sinar.
    """
    n = len(rows)
    observed = sum(1 for r in rows if r["real"])
    expected = sum(r[null_key] for r in rows)
    probability = expected / n if n else 0.0
    return {
        "regions": n,
        "observed_with_a_pair": observed,
        "expected_under_null": round(expected, 2),
        "observed_share": round(observed / n, 4) if n else None,
        "null_share": round(probability, 4),
        "ratio": round(observed / expected, 3) if expected else None,
        "binomial_two_sided_p": binomial_tail(observed, n, probability)
        if 0 < probability < 1 else None,
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--fasta",
                        default="analysis_out/upstream_regions.fasta.gz")
    parser.add_argument("--meta", default="analysis_out/upstream_regions.csv")
    parser.add_argument("--shuffles", type=int, default=SHUFFLES)
    parser.add_argument("--out",
                        default="analysis_out/sigma70_enrichment.json")
    args = parser.parse_args()

    regions = load_regions(args.fasta, args.meta)
    if not regions:
        sys.exit("bolge yok -- upstream_regions.py kosuldu mu?")
    print(f"[bolge] {len(regions)}", flush=True)

    con = sqlite3.connect(args.db)
    keys = levels(con)
    rng = random.Random(SEED)

    per_region = []
    for n, (cid, sequence, meta) in enumerate(regions, 1):
        real = 1 if count_pairs(sequence) else 0
        mono = sum(1 if count_pairs(mononucleotide_shuffle(sequence, rng))
                   else 0 for _ in range(args.shuffles)) / args.shuffles
        markov = sum(1 if count_pairs(markov1_resample(sequence, rng))
                     else 0 for _ in range(args.shuffles)) / args.shuffles
        per_region.append({"candidate_id": cid, "real": real,
                           "mono": mono, "markov": markov,
                           "architecture": meta["architecture"],
                           "length": int(meta["length"]),
                           "gc": float(meta["gc_fraction"])})
        if n % 1000 == 0:
            print(f"   {n}/{len(regions)}", flush=True)

    #  DUZEYLER. Birim degisince PAYDA degisir: giris duzeyinde her bolge
    #  bir birim; dizi ve cins duzeyinde ayni dizi/cinsten gelen bolgeler
    #  TEK birime indirilir (ortalama alinir), cunku ayni genomdan gelen
    #  bolgeler bagimsiz degil ve giris duzeyinde bir p degeri tekrarli
    #  ornekleme yuzunden kendiliginden kucuk cikar.
    def collapse(rows, key):
        grouped = defaultdict(list)
        for row in rows:
            grouped[key(row)].append(row)
        out = []
        for _k, group in grouped.items():
            out.append({
                "real": 1 if any(r["real"] for r in group) else 0,
                "mono": sum(r["mono"] for r in group) / len(group),
                "markov": sum(r["markov"] for r in group) / len(group)})
        return out

    result = {}
    for name, key in (("entry", None),
                      ("sequence", lambda r: keys.get(
                          r["candidate_id"], (r["candidate_id"],))[0]),
                      ("genus", lambda r: keys.get(
                          r["candidate_id"], (None, r["candidate_id"]))[1])):
        units = per_region if key is None else collapse(per_region, key)
        result[name] = {
            "vs_mononucleotide_shuffle": compare(units, "mono"),
            "vs_markov1_resample": compare(units, "markov"),
        }

    by_arch = defaultdict(list)
    for row in per_region:
        by_arch[row["architecture"]].append(row)

    report = {
        "question": ("are sigma70 -35/-10 box pairs more frequent in real "
                     "upstream intergenic regions than in a null built from "
                     "the same regions?"),
        "what_is_not_claimed": (
            "no promoter is located and none is predicted. The boxes are "
            "information-poor and occur often by chance, so a per-operon "
            "statement would not be publishable; only the enrichment is"),
        "method": {
            "boxes": {"minus_35": BOX35, "minus_10": BOX10,
                      "spacer_bp": list(SPACER),
                      "mismatches_allowed_per_box": MAX_MISMATCH},
            "nulls": {
                "mononucleotide_shuffle": ("bases of the same region "
                                           "shuffled: length and base "
                                           "composition, so GC, are exactly "
                                           "preserved; dinucleotide "
                                           "composition is not"),
                "markov1_resample": ("transition probabilities learned from "
                                     "the same region, a new sequence of the "
                                     "same length generated: dinucleotide "
                                     "composition preserved in expectation, "
                                     "which closes the first null's gap"),
                "why_two": ("one null is an assumption. The claim is "
                            "published only if both point the same way; if "
                            "they diverge, the divergence is the finding"),
            },
            "replicates_per_region": args.shuffles,
            "seed": SEED,
            "statistic": ("observed number of units carrying at least one pair "
                          "against the number expected under the null (the "
                          "sum of the null probabilities), tested with a "
                          "two-sided binomial. A sign test on per-unit "
                          "differences was tried first and was wrong: the "
                          "real value is binary while the null value is a "
                          "mean over replicates, so a region with no pair "
                          "counted against the real data whenever the null "
                          "found one in a single replicate"),
            "levels": ("entry, sequence and genus. Regions from one genome "
                       "are not independent, so an entry-level p value is "
                       "small by repeated sampling alone; the genus level "
                       "decides"),
            "strand": ("regions were written on the strand that drives the "
                       "operon by upstream_regions.py; nothing is "
                       "re-complemented here"),
        },
        "denominator": {
            "regions": len(regions),
            "note": ("only replicons whose nucleotide sequence is readable "
                     "could contribute; see replicon_sequence.csv and "
                     "ROADMAP 135.21. Whether the result generalises to the "
                     "unreadable majority is unmeasured"),
        },
        "levels": result,
        "by_architecture": {
            arch: {"regions": len(rows),
                   "share_with_a_pair": round(
                       sum(1 for r in rows if r["real"]) / len(rows), 4),
                   "null_mono": round(
                       sum(r["mono"] for r in rows) / len(rows), 4),
                   "null_markov": round(
                       sum(r["markov"] for r in rows) / len(rows), 4)}
            for arch, rows in sorted(by_arch.items())},
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    for name in ("entry", "sequence", "genus"):
        for null in ("vs_mononucleotide_shuffle", "vs_markov1_resample"):
            r = result[name][null]
            print(f"[{name:8s}] {null:28s} gozlenen {r['observed_with_a_pair']}"
                  f" / beklenen {r['expected_under_null']}  "
                  f"oran {r['ratio']}  p={r['binomial_two_sided_p']}")
    print(f"[yazildi] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
