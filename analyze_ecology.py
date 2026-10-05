"""
Ekolojik/evrimsel analiz: substrat sinifi mobiliteyi ongoruyor mu?

HIPOTEZ
  Ksenobiyotik substratli RO'lar (PAH, PCB, BTEX, ftalat, nitroaromatik) son
  yuzyilda ortaya cikan bilesikleri okside eder. Bakteriler bu yetenegi hizla
  edinmis olmali -- yani beklenen desen: mobil element komsulugu, plazmit
  uzerinde tasinma, genis ve dagilmis taksonomik yayilim (yatay gen transferi).

  Dogal substratli RO'lar (vanillat/lignin, steroid, kafein, kolin, karnitin)
  ise eski metabolik yollarin parcasi. Beklenen: kromozomal, mobil elementsiz,
  dar ve tutarli taksonomik dagilim (dikey kalitim).

TEST EDILEBILIR TAHMINLER
  1. Ksenobiyotik kumelerde transpozon komsulugu > dogal kumelerde
  2. Ksenobiyotik kumelerde plazmit orani > dogal kumelerde
  3. Taksonomik yayginlik (kac farkli cins) ksenobiyotiklerde daha genis
  4. Ksenobiyotiklerde sinteni korunumu DAHA DUSUK (yeni, henuz oturmamis yollar)

UYARI: cluster_ecology.csv'deki substrat atamalari literaturden en iyi cabayla
yapildi ve 'confidence' sutunuyla isaretlendi. 72 referansi kuratorleyen kisi
bu atamalari gozden gecirmeli -- 'dusuk' guvenli satirlar analize dahil
edilmiyor ama yine de kontrol edilmeli.
"""

import argparse
import csv
import math
import os
import sqlite3
from collections import Counter, defaultdict


VALID_CLASSES = {"xenobiotic", "natural_aromatic", "natural_specialized", "unknown"}
VALID_CONFIDENCE = {"low", "medium", "high"}


def load_ecology(path):
    """cluster -> {substrate, substrate_class, note, confidence}"""
    mapping = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            if row.get("substrate_class") not in VALID_CLASSES or \
                    row.get("confidence") not in VALID_CONFIDENCE:
                raise ValueError(f"{path}: bozuk satir (virgul kacisi?) -> {row}")
            mapping[row["cluster"]] = row
    return mapping


def genus_of(organism):
    if not organism:
        return None
    parts = organism.split()
    return parts[0] if parts else None


def shannon(counts):
    """Kategori dagiliminin Shannon entropisi -- sinteni cesitliligi olcusu."""
    total = sum(counts.values())
    if total <= 0:
        return 0.0
    return -sum((n / total) * math.log(n / total) for n in counts.values() if n > 0)


def collect(connection, ecology, window=5000):
    """Her kume icin mobilite, taksonomi ve sinteni olculerini topla."""
    # --- RO basina: cins, plazmit mi, transpozon yakin mi
    rows = connection.execute("""
        SELECT r.candidate_id, r.ro_cluster, rep.organism, rep.is_plasmid,
               MAX(CASE WHEN c.category='transposon' AND ABS(nb.distance) <= ?
                        THEN 1 ELSE 0 END) has_tn
        FROM ro r
        JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        LEFT JOIN neighbor nb ON nb.candidate_id = r.candidate_id
        LEFT JOIN gene_category c ON c.neighbor_id = nb.neighbor_id
             AND c.method = 'regex_v1'
        WHERE r.is_confirmed = 1 AND rep.status = 'ok'
        GROUP BY r.candidate_id
    """, (window,)).fetchall()

    stats = defaultdict(lambda: {"n": 0, "tn": 0, "plasmid": 0,
                                 "genera": Counter(), "ids": []})
    for row in rows:
        cluster = row["ro_cluster"]
        if not cluster or cluster == "N/A":
            continue
        entry = stats[cluster]
        entry["n"] += 1
        entry["tn"] += row["has_tn"] or 0
        entry["plasmid"] += row["is_plasmid"] or 0
        genus = genus_of(row["organism"])
        if genus:
            entry["genera"][genus] += 1
        entry["ids"].append(row["candidate_id"])

    # --- Sinteni cesitliligi: kume icinde komsu kategori dagiliminin entropisi
    for row in connection.execute("""
        SELECT r.ro_cluster cluster, c.category, COUNT(*) n
        FROM ro r
        JOIN neighbor nb ON nb.candidate_id = r.candidate_id
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id
             AND c.method = 'regex_v1'
        WHERE r.is_confirmed = 1 AND ABS(nb.gene_offset) <= 4
        GROUP BY r.ro_cluster, c.category
    """):
        cluster = row["cluster"]
        if cluster in stats:
            stats[cluster].setdefault("syn", Counter())[row["category"]] = row["n"]

    # --- Operon tamligi: alpha'nin yaninda beta + ferredoksin + reduktaz var mi
    for row in connection.execute("""
        SELECT r.ro_cluster cluster,
               SUM(CASE WHEN parts >= 3 THEN 1 ELSE 0 END) complete,
               COUNT(*) total
        FROM (
          SELECT r.candidate_id, r.ro_cluster,
                 COUNT(DISTINCT CASE WHEN c.category IN ('ro_beta','ferredoxin','reductase')
                       THEN c.category END) parts
          FROM ro r
          LEFT JOIN neighbor nb ON nb.candidate_id = r.candidate_id
               AND nb.same_strand = 1 AND ABS(nb.distance) <= 2000
          LEFT JOIN gene_category c ON c.neighbor_id = nb.neighbor_id
               AND c.method = 'regex_v1'
          WHERE r.is_confirmed = 1
          GROUP BY r.candidate_id
        ) r GROUP BY r.ro_cluster
    """):
        if row["cluster"] in stats:
            stats[row["cluster"]]["operon_complete"] = row["complete"]
            stats[row["cluster"]]["operon_total"] = row["total"]

    # Ekoloji bilgisini ekle
    for cluster, entry in stats.items():
        info = ecology.get(cluster, {})
        entry["substrate"] = info.get("substrate", "bilinmiyor")
        entry["class"] = info.get("substrate_class", "unknown")
        entry["confidence"] = info.get("confidence", "low")
        entry["note"] = info.get("ecology_note", "")
    return stats


def load_euk_rates(path="analysis_out/domain_by_cluster.csv"):
    """Kume basina okaryot orani -- bakteriyel mobilite testini confound eder."""
    if not os.path.exists(path):
        return {}
    rates = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            try:
                rates[row["cluster"]] = float(row["eukaryota_rate"])
            except (KeyError, ValueError):
                pass
    return rates


def report(stats, min_n=30, out_csv=None, euk_threshold=0.25):
    print("=" * 78)
    print("SUBSTRAT SINIFI x MOBILITE")
    print("=" * 78)

    euk_rates = load_euk_rates()
    excluded_euk = []

    by_class = defaultdict(lambda: {"n": 0, "tn": 0, "plasmid": 0,
                                    "genera": set(), "clusters": 0,
                                    "syn_entropy": [], "op_c": 0, "op_t": 0})
    for cluster, entry in stats.items():
        # Dusuk guvenli substrat atamalari toplu istatistige girmez
        if entry["confidence"] == "low" or entry["class"] == "unknown":
            continue
        # OKARYOT-AGIRLIKLI kumeler bakteriyel plazmit/transpozon testini bozar
        # (okaryotta plazmit/mobil element farkli calisir) -- dislanir.
        if euk_rates.get(cluster, 0) >= euk_threshold:
            excluded_euk.append((cluster, euk_rates[cluster]))
            continue
        bucket = by_class[entry["class"]]
        bucket["n"] += entry["n"]
        bucket["tn"] += entry["tn"]
        bucket["plasmid"] += entry["plasmid"]
        bucket["genera"] |= set(entry["genera"])
        bucket["clusters"] += 1
        if entry.get("syn"):
            bucket["syn_entropy"].append(shannon(entry["syn"]))
        bucket["op_c"] += entry.get("operon_complete", 0)
        bucket["op_t"] += entry.get("operon_total", 0)

    if excluded_euk:
        print(f"\n[dislandi] okaryot-agirlikli {len(excluded_euk)} kume "
              f"(>=%{100*euk_threshold:.0f} okaryot, bakteriyel mobiliteyi confound eder):")
        for cluster, rate in sorted(excluded_euk, key=lambda x: -x[1]):
            print(f"    {cluster:16s} %{100*rate:.0f} okaryot")
        print()

    header = (f"{'substrat sinifi':22s} {'kume':>5s} {'RO':>7s} {'Tn yakin':>9s} "
              f"{'plazmit':>8s} {'cins':>6s} {'operon tam':>11s}")
    print(header)
    print("-" * len(header))
    for name in sorted(by_class, key=lambda k: -by_class[k]["n"]):
        b = by_class[name]
        tn_rate = b["tn"] / b["n"] if b["n"] else 0
        pl_rate = b["plasmid"] / b["n"] if b["n"] else 0
        op_rate = b["op_c"] / b["op_t"] if b["op_t"] else 0
        print(f"{name:22s} {b['clusters']:>5d} {b['n']:>7,} {tn_rate:>8.1%} "
              f"{pl_rate:>7.1%} {len(b['genera']):>6d} {op_rate:>10.1%}")

    # --- Ana karsilastirma
    xeno = by_class.get("xenobiotic")
    natural = {"n": 0, "tn": 0, "plasmid": 0, "genera": set()}
    for name in ("natural_aromatic", "natural_specialized"):
        if name in by_class:
            natural["n"] += by_class[name]["n"]
            natural["tn"] += by_class[name]["tn"]
            natural["plasmid"] += by_class[name]["plasmid"]
            natural["genera"] |= by_class[name]["genera"]

    if xeno and natural["n"]:
        print("\n" + "=" * 78)
        print("HIPOTEZ TESTI: ksenobiyotik vs dogal substrat")
        print("=" * 78)
        xt = xeno["tn"] / xeno["n"]
        nt = natural["tn"] / natural["n"]
        xp = xeno["plasmid"] / xeno["n"]
        npl = natural["plasmid"] / natural["n"]
        print(f"  {'olcut':28s} {'ksenobiyotik':>14s} {'dogal':>10s} {'oran':>9s}")
        print(f"  {'-'*28} {'-'*14} {'-'*10} {'-'*9}")
        print(f"  {'transpozon komsulugu':28s} {xt:>13.1%} {nt:>9.1%} "
              f"{xt/nt if nt else float('inf'):>8.2f}x")
        print(f"  {'plazmit uzerinde':28s} {xp:>13.1%} {npl:>9.1%} "
              f"{xp/npl if npl else float('inf'):>8.2f}x")
        print(f"  {'farkli cins sayisi':28s} {len(xeno['genera']):>13,} "
              f"{len(natural['genera']):>9,} "
              f"{len(xeno['genera'])/len(natural['genera']) if natural['genera'] else 0:>8.2f}x")
        print(f"  {'toplam RO':28s} {xeno['n']:>13,} {natural['n']:>9,}")

        # Iki oranli z-testi
        pooled = (xeno["tn"] + natural["tn"]) / (xeno["n"] + natural["n"])
        se = math.sqrt(pooled * (1 - pooled) * (1 / xeno["n"] + 1 / natural["n"]))
        z = (xt - nt) / se if se else 0
        print(f"\n  transpozon farki icin z = {z:.2f} "
              f"({'anlamli' if abs(z) > 2.58 else 'anlamli degil'}, p<0.01 esigi |z|>2.58)")

    # --- Kume duzeyi tablo
    print("\n" + "=" * 78)
    print(f"KUME DUZEYI (n >= {min_n})")
    print("=" * 78)
    header = (f"{'kume':16s} {'sinif':20s} {'RO':>6s} {'Tn':>6s} {'plaz':>6s} "
              f"{'cins':>5s} {'sinteni H':>10s}")
    print(header)
    print("-" * len(header))
    ranked = sorted((c for c in stats.items() if c[1]["n"] >= min_n),
                    key=lambda x: -(x[1]["tn"] / x[1]["n"]))
    rows_out = []
    for cluster, entry in ranked:
        tn_rate = entry["tn"] / entry["n"]
        pl_rate = entry["plasmid"] / entry["n"]
        entropy = shannon(entry.get("syn", Counter()))
        print(f"{cluster:16s} {entry['class'][:20]:20s} {entry['n']:>6,} "
              f"{tn_rate:>5.1%} {pl_rate:>5.1%} {len(entry['genera']):>5d} "
              f"{entropy:>10.2f}")
        rows_out.append([cluster, entry["substrate"], entry["class"],
                         entry["confidence"], entry["n"], entry["tn"],
                         round(tn_rate, 4), entry["plasmid"], round(pl_rate, 4),
                         len(entry["genera"]), round(entropy, 3),
                         entry.get("operon_complete", 0),
                         entry.get("operon_total", 0), entry["note"]])

    if out_csv:
        with open(out_csv, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["cluster", "substrate", "substrate_class", "confidence",
                             "ro_count", "transposon_adjacent", "transposon_rate",
                             "plasmid_count", "plasmid_rate", "genus_count",
                             "synteny_entropy", "operon_complete", "operon_total",
                             "ecology_note"])
            writer.writerows(rows_out)
        print(f"\n[yazildi] {out_csv}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--ecology", default="cluster_ecology.csv")
    parser.add_argument("--window", type=int, default=5000)
    parser.add_argument("--min-n", type=int, default=30)
    parser.add_argument("--out", default="analysis_out/cluster_ecology_stats.csv")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    connection = sqlite3.connect(args.db)
    connection.row_factory = sqlite3.Row
    ecology = load_ecology(args.ecology)
    stats = collect(connection, ecology, args.window)
    report(stats, args.min_n, args.out)
    connection.close()


if __name__ == "__main__":
    main()
