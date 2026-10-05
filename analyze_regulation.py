"""
Regulasyon mimarisi: operonun 5' ucu, promotor bolgesi ve komsu duzenleyici.

NE OLCULUYOR
  Operon (build_operons.py) RO'yu iceren ayni iplikli gen dizisidir. Transkripsiyon
  bu dizinin 5' ucundan baslar; promotor o ucun yukarisindaki INTERGENIK BOLGEDE
  bulunur. Bu script her dogrulanmis RO icin:
    - operonun 5' ucunu (iplige gore) belirler,
    - yukari akistaki ilk geni ve aradaki intergenik uzunlugu olcer
      (= putatif promotor bolgesi),
    - o genin ayni iplikte mi (muhtemel okuma devami) yoksa ters iplikte mi
      (paylasilan divergent promotor bolgesi) oldugunu kaydeder,
    - duzenleyici ise ailesini (LysR, TetR, AraC, IclR ...) yazar.

NE OLCULMUYOR
  Promotor DIZISI tahmin edilmiyor. -35/-10 kutulari, operator tekrarlari ve
  transkripsiyon baslangici icin nukleotid dizisi gerekir; veritabaninda komsu
  genlerin protein cevirileri var, intergenik DNA yok. Burada raporlanan sey
  "promotorun bulunmasi gereken bolge ve onun genomik mimarisi"dir.

NEDEN DIVERGENT MIMARI ONEMLI
  LysR ailesi duzenleyiciler klasik olarak hedef operonla KAFA KAFAYA (divergent)
  oturur ve iki gen arasindaki kisa intergenik bolgeyi paylasir; duzenleyici hem
  kendi genini hem operonu ayni bolgeden kontrol eder. Bu mimari, duzenleyicinin
  operonla birlikte yatay transfer edilen islevsel bir modul olmasinin gostergesidir.
  Yon bilgisi olmadan bu analiz yapilamaz (eski pipeline'da strand kaydedilmiyordu).

Cikti: ro_regulation tablosu + analysis_out/regulation_by_cluster.csv
"""

import argparse
import csv
import os
import sqlite3
import statistics
from collections import Counter, defaultdict

# Bakterilerde divergent gen ciftlerinin paylastigi tipik promotor bolgesi
# genelde 60-400 bp'dir; bundan uzun bosluklar ayri transkripsiyon birimlerine
# isaret eder. Esik raporlamayi etkiler, veriyi degil (ham bosluk da saklanir).
MAX_SHARED_PROMOTER = 400


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out-dir", default="analysis_out")
    ap.add_argument("--max-shared", type=int, default=MAX_SHARED_PROMOTER)
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    con.executescript("""
        CREATE TABLE IF NOT EXISTS ro_regulation (
            candidate_id      TEXT PRIMARY KEY,
            operon_5p         INTEGER,  -- operonun 5' ucundaki koordinat
            upstream_gene_id  INTEGER,  -- neighbor_id (yukari akistaki ilk gen)
            intergenic_bp     INTEGER,  -- putatif promotor bolgesi uzunlugu
            upstream_strand   INTEGER,
            upstream_divergent INTEGER, -- ters iplikte mi (paylasilan promotor)
            upstream_category TEXT,
            upstream_family   TEXT,     -- duzenleyici ailesi (varsa)
            upstream_product  TEXT,
            architecture      TEXT      -- divergent_regulator | codirectional_regulator
                                        -- | divergent_other | codirectional_other | unknown
        );
        CREATE INDEX IF NOT EXISTS idx_reg_arch ON ro_regulation(architecture);
    """)

    families = defaultdict(dict)
    for nid, cat, method in con.execute(
            "SELECT neighbor_id, category, method FROM gene_category"):
        families[nid][method] = cat

    operons = {r["candidate_id"]: r for r in con.execute("SELECT * FROM operon")}
    neighbors = defaultdict(list)
    for r in con.execute("SELECT * FROM neighbor"):
        neighbors[r["candidate_id"]].append(r)

    rows = []
    for cid, op in operons.items():
        strand = op["strand"]
        nbs = [n for n in neighbors.get(cid, ()) if not n["spans_origin"]]
        if strand == 1:
            five_prime = op["start"]
            cands = [n for n in nbs if n["end"] <= five_prime]
            up = max(cands, key=lambda n: n["end"]) if cands else None
            gap = (five_prime - up["end"]) if up else None
        else:
            five_prime = op["end"]
            cands = [n for n in nbs if n["start"] >= five_prime]
            up = min(cands, key=lambda n: n["start"]) if cands else None
            gap = (up["start"] - five_prime) if up else None

        if up is None:
            rows.append((cid, five_prime, None, None, None, None, None, None, None, "unknown"))
            continue
        cat = families.get(up["neighbor_id"], {}).get("regex_v1")
        fam = families.get(up["neighbor_id"], {}).get("regulator_family_v1")
        divergent = int(up["strand"] != strand)
        if cat == "regulator":
            arch = "divergent_regulator" if divergent else "codirectional_regulator"
        else:
            arch = "divergent_other" if divergent else "codirectional_other"
        rows.append((cid, five_prime, up["neighbor_id"], gap, up["strand"], divergent,
                     cat, fam, up["product"], arch))

    con.execute("DELETE FROM ro_regulation")
    con.executemany("INSERT INTO ro_regulation VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit()

    arch = Counter(r[9] for r in rows)
    total = len(rows)
    print(f"[ro_regulation] {total} RO")
    for a, n in arch.most_common():
        print(f"   {a:26s} {n:6d}  ({100*n/total:.1f}%)")

    div_reg = [r[3] for r in rows if r[9] == "divergent_regulator" and r[3] is not None]
    other = [r[3] for r in rows if r[9] != "divergent_regulator" and r[3] is not None]
    if div_reg and other:
        print(f"   intergenik medyan: divergent duzenleyici {statistics.median(div_reg):.0f} bp, "
              f"digerleri {statistics.median(other):.0f} bp")
        shared = sum(1 for g in div_reg if g <= args.max_shared)
        print(f"   divergent duzenleyicilerin {100*shared/len(div_reg):.1f}%'i "
              f"<= {args.max_shared} bp (paylasilan promotor bolgesi olabilir)")

    fam_counts = Counter(r[7] for r in rows if r[9] == "divergent_regulator" and r[7])
    print("   divergent duzenleyici aileleri:",
          ", ".join(f"{k} {v}" for k, v in fam_counts.most_common(8)))

    cluster = dict(con.execute("SELECT candidate_id, ro_cluster FROM ro WHERE is_confirmed=1"))
    per = defaultdict(lambda: {"n": 0, "arch": Counter(), "fam": Counter(), "gaps": []})
    for r in rows:
        cl = cluster.get(r[0])
        if cl is None:
            continue
        d = per[cl]
        d["n"] += 1
        d["arch"][r[9]] += 1
        if r[7]:
            d["fam"][r[7]] += 1
        if r[3] is not None:
            d["gaps"].append(r[3])

    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, "regulation_by_cluster.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cluster", "n", "divergent_regulator", "codirectional_regulator",
                    "divergent_other", "codirectional_other", "unknown",
                    "divergent_regulator_rate", "median_intergenic_bp", "top_regulator_family"])
        for cl, d in sorted(per.items(), key=lambda kv: -kv[1]["n"]):
            w.writerow([cl, d["n"], d["arch"]["divergent_regulator"],
                        d["arch"]["codirectional_regulator"], d["arch"]["divergent_other"],
                        d["arch"]["codirectional_other"], d["arch"]["unknown"],
                        round(d["arch"]["divergent_regulator"] / d["n"], 3),
                        int(statistics.median(d["gaps"])) if d["gaps"] else "",
                        d["fam"].most_common(1)[0][0] if d["fam"] else ""])
    print(f"[yazildi] {path}")
    con.close()


if __name__ == "__main__":
    main()
