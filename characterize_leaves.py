"""
Yaprak (varyant) karakterizasyonu: bir kumeye atanmis farkli varyantlari ayirt et.

Ozyinelemeli homojenizasyon her kumeyi homojen yapraklara boldu. Ama "bu yaprak
farkli bir enzim mi" sorusu sadece dizi kimligiyle cevaplanamaz -- iki yaprak
%99 ici kimlikte olabilir ama BIRBIRINDEN cok farkli olabilir ve farkli
substratlara etki edebilir.

Bu modul her yaprak icin ayirt edici bir IMZA cikarir:

  1. TAKSONOMI      hangi cinslerde, ne kadar dagilmis
  2. MOBILITE       plazmit orani, transpozon komsulugu
  3. KOMSULUK       yaprak uyelerinin cevresinde en sik gecen genler
                    (yapisal ro_beta/ferredoksin/reduktaz haric -- bunlar her
                    yerde; ayirt edici olan tasiyici tipi, ozgul dioksijenazlar,
                    halka-acilim enzimleri gibi YOL-OZGU komsular)

Komsuluk imzasi en degerlisi: iki yaprak farkli genomik baglamda oturuyorsa
(biri ftalat tasiyicisi yaninda, digeri benzoat operonunda), ayni kumede
gorunseler bile buyuk olasilikla islevsel olarak farklidirlar.

Cikti:
    roar.sqlite   yeni tablo: leaf_profile
    analysis_out/leaf_profiles.csv
"""

import argparse
import csv
import os
import sqlite3
from collections import Counter, defaultdict

# Her yerde bulunan yapisal komsular -- varyanti AYIRT ETMEZLER, imzadan cikar
UBIQUITOUS = {"ro_alpha", "ro_beta", "ferredoxin", "reductase", "hypothetical",
              "other", "transporter"}

# Bir yaprak "profillenecek" minimum uye (daha kucukler sadece temel istatistik)
MIN_PROFILE_SIZE = 5

# Jenerik urun adlarini imzadan ele -- ayirt edici degiller
GENERIC_PRODUCT = ("hypothetical protein", "DUF", "uncharacterized",
                   "domain-containing protein")


def is_generic(product):
    if not product:
        return True
    low = product.lower()
    return any(g.lower() in low for g in GENERIC_PRODUCT)


def label_leaf(top_genus, genus_n, plasmid_rate, transposon_rate, signature):
    """Yaprak icin okunabilir otomatik etiket.

    METIN INGILIZCE: bu etiket varyant ve tip sayfalarinda dogrudan gosteriliyor,
    yani arayuz metnidir. Kod yorumlari Turkce kalir, kullanicinin okudugu her sey
    Ingilizce olmak zorunda.
    """
    parts = []
    if genus_n == 1:
        parts.append(f"confined to {top_genus}")
    elif genus_n <= 3:
        parts.append(f"mostly {top_genus}, narrow host range")
    else:
        parts.append(f"mostly {top_genus}")
    if plasmid_rate >= 0.30:
        parts.append("often plasmid-borne")
    if transposon_rate >= 0.30:
        parts.append("near a mobile element")
    if signature:
        parts.append(f"neighbour: {signature[0]}")
    return "; ".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default="analysis_out")
    parser.add_argument("--window", type=int, default=5000)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    connection = sqlite3.connect(args.db)
    connection.row_factory = sqlite3.Row

    # Uye -> yaprak
    member_leaf = {row["candidate_id"]: row["leaf_id"]
                   for row in connection.execute("SELECT candidate_id, leaf_id FROM ro_leaf")}

    # Yaprak temel bilgisi
    leaves = {row["leaf_id"]: dict(row) for row in connection.execute(
        "SELECT leaf_id, cluster, size, depth, median_identity FROM leaf")}

    # Uye basina: cins, plazmit, transpozon-yakin
    stats = defaultdict(lambda: {"genera": Counter(), "plasmid": 0, "tn": 0, "n": 0})
    for row in connection.execute("""
        SELECT r.candidate_id, rep.organism, rep.is_plasmid,
               MAX(CASE WHEN c.category='transposon' AND ABS(nb.distance) <= ?
                        THEN 1 ELSE 0 END) has_tn
        FROM ro r
        JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        LEFT JOIN neighbor nb ON nb.candidate_id = r.candidate_id
        LEFT JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        WHERE r.is_confirmed = 1
        GROUP BY r.candidate_id
    """, (args.window,)):
        leaf_id = member_leaf.get(row["candidate_id"])
        if not leaf_id:
            continue
        entry = stats[leaf_id]
        entry["n"] += 1
        entry["plasmid"] += row["is_plasmid"] or 0
        entry["tn"] += row["has_tn"] or 0
        organism = row["organism"] or ""
        if organism.split():
            entry["genera"][organism.split()[0]] += 1

    # Komsuluk imzasi: yaprak uyelerinin +-4 gen komsu URUNLERI
    leaf_products = defaultdict(Counter)
    leaf_categories = defaultdict(Counter)
    for row in connection.execute("""
        SELECT nb.candidate_id, nb.product, c.category
        FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        WHERE ABS(nb.gene_offset) <= 4
    """):
        leaf_id = member_leaf.get(row["candidate_id"])
        if not leaf_id:
            continue
        if row["category"] not in UBIQUITOUS and not is_generic(row["product"]):
            leaf_products[leaf_id][row["product"]] += 1
        if row["category"] not in ("ro_alpha", "other", "hypothetical"):
            leaf_categories[leaf_id][row["category"]] += 1

    # Kume ortalamasi -- bir yapragin komsulugu kume genelinden farkli mi
    cluster_cat = defaultdict(Counter)
    cluster_n = Counter()
    for leaf_id, entry in stats.items():
        cluster = leaves[leaf_id]["cluster"]
        cluster_n[cluster] += entry["n"]
        for category, count in leaf_categories[leaf_id].items():
            cluster_cat[cluster][category] += count

    connection.executescript("""
        DROP TABLE IF EXISTS leaf_profile;
        CREATE TABLE leaf_profile (
            leaf_id        TEXT PRIMARY KEY,
            cluster        TEXT,
            size           INTEGER,
            median_identity REAL,
            depth          INTEGER,
            genus_n        INTEGER,
            top_genera     TEXT,
            plasmid_rate   REAL,
            transposon_rate REAL,
            neighbor_signature TEXT,  -- ayirt edici komsu urunleri
            enriched_categories TEXT, -- kume ortalamasina gore zengin kategoriler
            label          TEXT
        );
        CREATE INDEX idx_lp_cluster ON leaf_profile(cluster);
    """)

    rows_out = []
    for leaf_id, info in leaves.items():
        entry = stats.get(leaf_id, {"genera": Counter(), "plasmid": 0, "tn": 0, "n": 0})
        n = entry["n"] or info["size"]
        top_genera = entry["genera"].most_common(5)
        top_genus = top_genera[0][0] if top_genera else "?"
        genus_n = len(entry["genera"])
        plasmid_rate = entry["plasmid"] / n if n else 0
        tn_rate = entry["tn"] / n if n else 0

        # Ayirt edici komsu urunleri (en sik 4, jenerik olmayan)
        signature = [p for p, _ in leaf_products[leaf_id].most_common(4)]

        # Kume ortalamasina gore zengin kategoriler
        cluster = info["cluster"]
        enriched = []
        total_cluster = cluster_n[cluster] or 1
        for category, count in leaf_categories[leaf_id].most_common():
            leaf_rate = count / max(1, n)
            cluster_rate = cluster_cat[cluster][category] / total_cluster
            if cluster_rate > 0 and leaf_rate / cluster_rate >= 1.5 and count >= 3:
                enriched.append(f"{category}({leaf_rate/cluster_rate:.1f}x)")

        label = label_leaf(top_genus, genus_n, plasmid_rate, tn_rate, signature)
        rows_out.append((
            leaf_id, cluster, n, round(info["median_identity"], 4), info["depth"],
            genus_n, ";".join(f"{g}:{c}" for g, c in top_genera),
            round(plasmid_rate, 4), round(tn_rate, 4),
            " | ".join(signature[:4]), "; ".join(enriched[:4]), label))

    connection.executemany(
        "INSERT INTO leaf_profile VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rows_out)
    connection.commit()

    path = os.path.join(args.out_dir, "leaf_profiles.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["leaf_id", "cluster", "size", "median_identity", "depth",
                         "genus_n", "top_genera", "plasmid_rate", "transposon_rate",
                         "neighbor_signature", "enriched_categories", "label"])
        writer.writerows(rows_out)

    # --- Ozet: en cok varyanta bolunen kumeler ve imzalarindaki farklar
    big = [r for r in rows_out if r[2] >= MIN_PROFILE_SIZE]
    print(f"[yazildi] {path}  ({len(rows_out)} yaprak, {len(big)} profillendi)")

    print("\n" + "=" * 72)
    print("VARYANT ORNEKLERI: ayni kume, farkli komsuluk imzasi")
    print("=" * 72)
    # Komsuluk imzasi olan buyuk yapraklardan ornek kumeler
    by_cluster = defaultdict(list)
    for r in rows_out:
        if r[2] >= 10 and r[9]:   # size>=10 ve imza var
            by_cluster[r[1]].append(r)
    shown = 0
    for cluster, leaf_rows in sorted(by_cluster.items(),
                                     key=lambda x: -len(x[1])):
        if len(leaf_rows) < 2:
            continue
        print(f"\n{cluster}  ({len(leaf_rows)} buyuk varyant):")
        for r in sorted(leaf_rows, key=lambda x: -x[2])[:4]:
            print(f"  {r[0]:20s} n={r[2]:>4d} id={r[3]:.2f}  {r[11]}")
            if r[9]:
                print(f"       neighbour: {r[9][:64]}")
        shown += 1
        if shown >= 6:
            break
    connection.close()


if __name__ == "__main__":
    main()
