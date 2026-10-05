"""
ROAR veritabani analizleri: operon, orientation, transpozon, regulator.

Bu analizlerin hicbiri mevcut pipeline'da MUMKUN DEGILDI, cunku strand
kaydedilmiyordu (mongo_gbk_analysis.py sadece start/end aliyor). Strand
gelince su sorular sorulabilir hale geliyor:

  - Hangi genler RO ile ayni operonda tasiniyor?
  - Regulator RO'ya gore hangi yonde duruyor? (divergent = klasik LysR mimarisi)
  - Hangi taksonlarda RO kumeleri mobil element komsulugunda?
  - RO kumeleri plazmitte mi kromozomda mi zenginlesiyor?
"""

import argparse
import csv
import os
import sqlite3
from collections import Counter, defaultdict

# Ayni operonda kabul edilecek maksimum intergenik bosluk.
# Bakteriyel operonlarda gen arasi mesafe tipik olarak <50 bp; 100 bp
# muhafazakar bir ust sinir.
OPERON_MAX_GAP = 100

# Divergent regulator: RO'nun SOLUNDA, TERS yonde, promotor mesafesinde.
# LysR ailesi duzenleyiciler hedef operonun karsisina kafa kafaya oturur.
DIVERGENT_MAX_DISTANCE = 400


def connect(db_path):
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row
    return connection


def genus_of(organism):
    """Organizma adindan cins (genus) cikar."""
    if not organism:
        return "Unknown"
    parts = organism.split()
    return parts[0] if parts else "Unknown"


def section(title):
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")


def overview(connection):
    section("GENEL DURUM")
    queries = [
        ("replikon (toplam)", "SELECT COUNT(*) FROM replicon"),
        ("  anotasyonlu", "SELECT COUNT(*) FROM replicon WHERE status='ok'"),
        ("  ANOTASYONSUZ", "SELECT COUNT(*) FROM replicon WHERE status='no_annotation'"),
        ("  plazmit", "SELECT COUNT(*) FROM replicon WHERE is_plasmid=1"),
        ("RO adayi", "SELECT COUNT(*) FROM ro"),
        ("  DOGRULANMIS RO", "SELECT COUNT(*) FROM ro WHERE is_confirmed=1"),
        ("komsu kaydi", "SELECT COUNT(*) FROM neighbor"),
    ]
    for label, query in queries:
        print(f"  {label:22s}: {connection.execute(query).fetchone()[0]:>9,}")

    row = connection.execute(
        "SELECT COUNT(*) n FROM replicon WHERE status='no_annotation'").fetchone()
    total = connection.execute("SELECT COUNT(*) FROM replicon").fetchone()[0]
    if total:
        print(f"\n  UYARI: replikonlarin %{100*row['n']/total:.1f}'i anotasyonsuz. "
              f"Bunlar 'komsusuz RO' degil, 'verisi eksik' -- her "
              f"co-occurrence oranini asagi ceker, paydadan cikarilmali.")


def orientation(connection):
    section("ORIENTATION: komsular RO ile ayni yonde mi?")
    rows = connection.execute("""
        SELECT c.category,
               COUNT(*) total,
               SUM(n.same_strand) same,
               AVG(ABS(n.distance)) avg_dist
        FROM neighbor n
        JOIN gene_category c ON c.neighbor_id = n.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id = n.candidate_id AND r.is_confirmed = 1
        GROUP BY c.category
        HAVING total >= 50
        ORDER BY 1.0*same/total DESC
    """).fetchall()
    print(f"  {'kategori':16s} {'n':>8s} {'ayni yon':>10s} {'ort.mesafe':>11s}")
    for row in rows:
        fraction = row["same"] / row["total"] if row["total"] else 0
        print(f"  {row['category']:16s} {row['total']:>8,} {fraction:>9.1%} "
              f"{row['avg_dist']:>10.0f}bp")
    print("\n  Yorum: ayni-yon orani %50 civari = rastgele. Belirgin yuksek =")
    print("  operonik birliktelik. Belirgin dusuk = divergent (karsit) mimari.")


def operons(connection):
    """RO ile ayni yonde ve <OPERON_MAX_GAP boslukla duran komsular."""
    section(f"OPERON: RO ile ayni yonde, <{OPERON_MAX_GAP}bp bosluk")
    rows = connection.execute("""
        SELECT c.category, COUNT(*) n
        FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id = nb.candidate_id AND r.is_confirmed = 1
        WHERE nb.same_strand = 1 AND ABS(nb.distance) <= ?
        GROUP BY c.category ORDER BY n DESC LIMIT 12
    """, (OPERON_MAX_GAP,)).fetchall()
    for row in rows:
        print(f"  {row['category']:18s} {row['n']:>7,}")


def divergent_regulators(connection):
    section(f"DIVERGENT REGULATOR (RO'nun 5' tarafinda, ters yonde = kafa kafaya, <{DIVERGENT_MAX_DISTANCE}bp)")
    print("  LysR ailesi duzenleyicilerin klasik imzasi. STRAND olmadan")
    print("  tespit edilemez -- mevcut pipeline'da bu analiz mumkun degildi.\n")

    divergent = connection.execute("""
        SELECT COUNT(*) n FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id = nb.candidate_id AND r.is_confirmed = 1
        WHERE c.category='regulator' AND nb.same_strand = 0
          AND ((r.strand = 1 AND nb.distance < 0) OR (r.strand = -1 AND nb.distance > 0))
          AND ABS(nb.distance) <= ?
    """, (DIVERGENT_MAX_DISTANCE,)).fetchone()["n"]

    codirectional = connection.execute("""
        SELECT COUNT(*) n FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id = nb.candidate_id AND r.is_confirmed = 1
        WHERE c.category='regulator' AND nb.same_strand = 1
          AND ABS(nb.distance) <= ?
    """, (DIVERGENT_MAX_DISTANCE,)).fetchone()["n"]

    print("  NOT: distance genom koordinatinda isaretli; kafa kafaya geometri RO ipligine bagli\n"
          "       (+ iplik: solda, - iplik: sagda). Eski sorgu bunu yoksayiyordu.")
    print(f"  divergent (kafa kafaya) : {divergent:>7,}")
    print(f"  ko-direksiyonel         : {codirectional:>7,}")
    if divergent + codirectional:
        print(f"  divergent orani         : {divergent/(divergent+codirectional):>7.1%}")

    print("\n  En sik divergent regulator urunleri:")
    for row in connection.execute("""
        SELECT nb.product, COUNT(*) n FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id = nb.candidate_id AND r.is_confirmed = 1
        WHERE c.category='regulator' AND nb.same_strand = 0
          AND ((r.strand = 1 AND nb.distance < 0) OR (r.strand = -1 AND nb.distance > 0))
          AND ABS(nb.distance) <= ?
        GROUP BY nb.product ORDER BY n DESC LIMIT 10
    """, (DIVERGENT_MAX_DISTANCE,)):
        print(f"    {row['n']:>5,}  {row['product'][:60]}")


def transposon_by_taxon(connection, window=5000, min_ro=20, output_csv=None):
    """Hangi taksonlarda RO kumeleri mobil element komsulugunda?"""
    section(f"TRANSPOZON ILISKISI (RO'nun +-{window}bp icinde mobil element)")

    rows = connection.execute("""
        SELECT rep.organism, r.candidate_id,
               MAX(CASE WHEN c.category='transposon' AND ABS(nb.distance) <= ?
                        THEN 1 ELSE 0 END) has_tn
        FROM ro r
        JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        LEFT JOIN neighbor nb ON nb.candidate_id = r.candidate_id
        LEFT JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        WHERE r.is_confirmed = 1 AND rep.status = 'ok'
        GROUP BY r.candidate_id
    """, (window,)).fetchall()

    by_genus = defaultdict(lambda: [0, 0])   # genus -> [ro sayisi, transpozonlu]
    for row in rows:
        stats = by_genus[genus_of(row["organism"])]
        stats[0] += 1
        stats[1] += row["has_tn"]

    total_ro = sum(v[0] for v in by_genus.values())
    total_tn = sum(v[1] for v in by_genus.values())
    baseline = total_tn / total_ro if total_ro else 0
    print(f"  Genel taban oran: {total_tn:,}/{total_ro:,} = {baseline:.1%}\n")

    ranked = [(genus, count, tn, tn / count)
              for genus, (count, tn) in by_genus.items() if count >= min_ro]
    ranked.sort(key=lambda x: -x[3])

    print(f"  {'cins':28s} {'RO':>6s} {'Tn yakin':>9s} {'oran':>7s} {'zenginlesme':>12s}")
    print(f"  {'-'*28} {'-'*6} {'-'*9} {'-'*7} {'-'*12}")
    for genus, count, tn, fraction in ranked[:20]:
        enrichment = fraction / baseline if baseline else 0
        print(f"  {genus[:28]:28s} {count:>6,} {tn:>9,} {fraction:>6.1%} "
              f"{enrichment:>11.2f}x")

    if len(ranked) > 20:
        print(f"\n  ...en dusuk 5:")
        for genus, count, tn, fraction in ranked[-5:]:
            enrichment = fraction / baseline if baseline else 0
            print(f"  {genus[:28]:28s} {count:>6,} {tn:>9,} {fraction:>6.1%} "
                  f"{enrichment:>11.2f}x")

    if output_csv:
        with open(output_csv, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["genus", "ro_count", "transposon_adjacent",
                             "fraction", "enrichment_vs_baseline"])
            for genus, count, tn, fraction in ranked:
                writer.writerow([genus, count, tn, round(fraction, 4),
                                 round(fraction / baseline, 3) if baseline else 0])
        print(f"\n  [yazildi] {output_csv}")

    return ranked, baseline


def plasmid_enrichment(connection, window=5000):
    section("PLAZMIT vs KROMOZOM")
    row = connection.execute("""
        SELECT rep.is_plasmid, COUNT(DISTINCT r.candidate_id) n
        FROM ro r JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        WHERE r.is_confirmed = 1 GROUP BY rep.is_plasmid
    """).fetchall()
    for entry in row:
        label = "plazmit" if entry["is_plasmid"] else "kromozom"
        print(f"  {label:12s}: {entry['n']:>7,} RO")

    print("\n  Plazmit/kromozomda transpozon komsulugu:")
    for entry in connection.execute("""
        SELECT rep.is_plasmid,
               COUNT(DISTINCT r.candidate_id) total,
               COUNT(DISTINCT CASE WHEN c.category='transposon'
                     AND ABS(nb.distance) <= ? THEN r.candidate_id END) tn
        FROM ro r
        JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        LEFT JOIN neighbor nb ON nb.candidate_id = r.candidate_id
        LEFT JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        WHERE r.is_confirmed = 1 GROUP BY rep.is_plasmid
    """, (window,)):
        label = "plazmit" if entry["is_plasmid"] else "kromozom"
        fraction = entry["tn"] / entry["total"] if entry["total"] else 0
        print(f"  {label:12s}: {entry['tn']:>6,}/{entry['total']:<7,} = {fraction:.1%}")


def synteny_strings(connection, limit=15):
    """RO cevresindeki gen dizilimini yon oklariyla yazdir."""
    section("SINTENI ORNEKLERI (gen sirasi + yon)")
    candidates = connection.execute("""
        SELECT r.candidate_id, r.ro_cluster, rep.organism
        FROM ro r JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        WHERE r.is_confirmed = 1 AND rep.status='ok'
        ORDER BY r.candidate_id LIMIT ?
    """, (limit,)).fetchall()

    for candidate in candidates:
        neighbors = connection.execute("""
            SELECT nb.gene_offset, nb.same_strand, c.category
            FROM neighbor nb JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
            WHERE nb.candidate_id = ? AND ABS(nb.gene_offset) <= 4
            ORDER BY nb.gene_offset
        """, (candidate["candidate_id"],)).fetchall()
        parts = []
        for neighbor in neighbors:
            arrow = "->" if neighbor["same_strand"] else "<-"
            parts.append(f"{neighbor['category'][:4]}{arrow}")
            if neighbor["gene_offset"] == -1:
                parts.append(f"[{candidate['ro_cluster'] or 'RO'}]->")
        if not any("[" in p for p in parts):
            parts.append(f"[{candidate['ro_cluster'] or 'RO'}]->")
        print(f"  {(candidate['organism'] or '?')[:34]:34s} {' '.join(parts)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--window", type=int, default=5000)
    parser.add_argument("--out-dir", default="analysis_out")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    connection = connect(args.db)

    overview(connection)
    orientation(connection)
    operons(connection)
    divergent_regulators(connection)
    transposon_by_taxon(connection, args.window,
                        output_csv=os.path.join(args.out_dir, "transposon_by_genus.csv"))
    plasmid_enrichment(connection, args.window)
    synteny_strings(connection)

    connection.close()


if __name__ == "__main__":
    main()
