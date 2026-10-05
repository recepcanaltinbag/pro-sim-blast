"""
Komsuluk zenginlesmesi icin NULL MODEL: ayni replikonda rastgele pencereler.

NEDEN GEREKLI: "RO'larin %40'inin yaninda bir regulator var" tek basina hicbir
sey soylemez. Regulatorler bakteriyel genomlarda zaten yaygin -- rastgele bir
20 kb pencerede de buyuk olasilikla bir tane vardir. Iddianin savunulabilir
olmasi icin karsilastirma noktasi gerekir.

Mevcut pipeline'in (ve ALL.csv'deki Score hesabinin) en buyuk metodolojik
acigi bu: co-occurrence sayiliyor ama arka plana gore normalize edilmiyor.

BURADAKI KONTROL: dogrulanmis RO tasiyan her replikondan, RO penceresiyle
AYNI GENISLIKTE rastgele pencereler ornekleniyor ve ayni kategori siniflandirmasi
uygulaniyor. Boylece "RO yaninda" frekansi "ayni genomda rastgele yerde"
frekansiyla kiyaslanabilir. Genom kompozisyonu, GC, gen yogunlugu ve anotasyon
kalitesi otomatik olarak kontrol edilmis olur -- cunku karsilastirma AYNI
replikon icinde yapiliyor.
"""

import argparse
import csv
import os
import random
import sqlite3
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

from Bio import SeqIO

from build_db import categorize
from extract_genomic_context import feature_spans, NEIGHBOR_WINDOW

# Her replikondan kac rastgele pencere ornekleneceği. Yuksek deger daha duzgun
# arka plan tahmini verir; 5 pratikte yeterli (replikon sayisi zaten binlerce).
WINDOWS_PER_REPLICON = 5

# Tekrarlanabilirlik: sabit tohum. Rastgelelik pencere KONUMUNDA, sonucta degil.
RANDOM_SEED = 20260721


def sample_windows(args):
    """Bir gbk dosyasindan rastgele pencereler ornekle, kategori say.

    Donen: (nucleotide_id, [Counter, ...]) -- her pencere icin kategori sayimlari
    """
    path, ro_positions, seed = args
    try:
        records = list(SeqIO.parse(path, "genbank"))
    except Exception:
        return None, []

    results = []
    rng = random.Random(seed)

    for record in records:
        cds_list = []
        for feature in record.features:
            if feature.type != "CDS":
                continue
            spans = feature_spans(feature)
            cds_list.append((min(s for s, _ in spans), max(e for _, e in spans),
                             feature.qualifiers.get("product", ["hypothetical protein"])[0]))
        if len(cds_list) < 5:
            continue

        length = len(record.seq)
        half = NEIGHBOR_WINDOW
        # RO'lardan uzak duran rastgele merkezler sec -- ayni bolgeyi tekrar
        # ornekleyip "kontrol"u kirletmemek icin.
        forbidden = ro_positions.get(record.id, [])
        for _ in range(WINDOWS_PER_REPLICON):
            for _attempt in range(20):
                center = rng.randint(half, max(half + 1, length - half))
                if all(abs(center - p) > 2 * half for p in forbidden):
                    break
            window_start, window_end = center - half, center + half
            counts = Counter()
            for start, end, product in cds_list:
                if end < window_start or start > window_end:
                    continue
                counts[categorize(product)] += 1
            counts["__genes__"] = sum(v for k, v in counts.items() if k != "__genes__")
            results.append(counts)

    return (records[0].id if records else None), results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--gbk-dir", default="gbk_files")
    parser.add_argument("--out", default="analysis_out/null_model.csv")
    parser.add_argument("--processes", type=int, default=6)
    parser.add_argument("--max-replicons", type=int, default=4000,
                        help="Ornekleme icin kac replikon kullanilacak")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    connection = sqlite3.connect(args.db)
    connection.row_factory = sqlite3.Row

    # Dogrulanmis RO tasiyan replikonlar -- kontrol AYNI genomlarda yapilmali
    rows = connection.execute("""
        SELECT rep.nucleotide_id, rep.file, r.start, r.end
        FROM ro r JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        WHERE r.is_confirmed = 1 AND rep.status = 'ok'
    """).fetchall()
    if not rows:
        print("[hata] dogrulanmis RO yok -- once annotate_ro.py calistirilmali")
        return 1

    by_file = defaultdict(lambda: defaultdict(list))
    for row in rows:
        by_file[row["file"]][row["nucleotide_id"]].append(
            (row["start"] + row["end"]) // 2)

    files = sorted(by_file)
    rng = random.Random(RANDOM_SEED)
    if len(files) > args.max_replicons:
        files = rng.sample(files, args.max_replicons)
    print(f"[bilgi] {len(files)} replikon ornekleniyor, "
          f"her birinden {WINDOWS_PER_REPLICON} pencere")

    tasks = [(os.path.join(args.gbk_dir, f), dict(by_file[f]), RANDOM_SEED + i)
             for i, f in enumerate(files)]

    # --- Kontrol pencerelerini topla
    background = Counter()
    window_count = 0
    with Pool(args.processes) as pool:
        for index, (_, windows) in enumerate(
                pool.imap_unordered(sample_windows, tasks, chunksize=8), 1):
            for counts in windows:
                background.update(counts)
                window_count += 1
            if index % 500 == 0:
                sys.stdout.write(f"\r  {index}/{len(tasks)} replikon | "
                                 f"{window_count} kontrol penceresi")
                sys.stdout.flush()
    print(f"\r  {len(tasks)} replikon | {window_count} kontrol penceresi")

    # --- Gercek RO pencereleri
    observed = Counter()
    # Kontrol pencereleri yalnizca >=5 CDS'li replikonlardan ornekleniyor;
    # payda da ayni popülasyon olmali, yoksa oranlar ~%9 asagi cekilir.
    ro_windows = connection.execute(
        "SELECT COUNT(*) FROM ro r JOIN replicon p USING(nucleotide_id) "
        "WHERE r.is_confirmed=1 AND p.cds_count >= 5").fetchone()[0]
    for row in connection.execute("""
        SELECT c.category, COUNT(*) n
        FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        JOIN ro r ON r.candidate_id = nb.candidate_id AND r.is_confirmed = 1
        JOIN replicon p ON p.nucleotide_id = r.nucleotide_id AND p.cds_count >= 5
        GROUP BY c.category
    """):
        observed[row[0]] = row[1]
    observed["__genes__"] = sum(v for k, v in observed.items() if k != "__genes__")

    # --- Zenginlesme
    print(f"\n{'kategori':16s} {'RO/pencere':>11s} {'kontrol/pencere':>16s} "
          f"{'zenginlesme':>12s}")
    print("-" * 60)
    results = []
    for category in sorted(set(observed) | set(background)):
        if category == "__genes__":
            continue
        ro_rate = observed.get(category, 0) / max(1, ro_windows)
        null_rate = background.get(category, 0) / max(1, window_count)
        enrichment = ro_rate / null_rate if null_rate else float("inf")
        results.append((category, observed.get(category, 0), ro_rate,
                        background.get(category, 0), null_rate, enrichment))

    results.sort(key=lambda x: -x[5])
    for category, obs_n, ro_rate, null_n, null_rate, enrichment in results:
        marker = "  <<<" if enrichment >= 2 else ("  (dusuk)" if enrichment < 0.5 else "")
        display = f"{enrichment:>11.2f}x" if enrichment != float("inf") else "        inf"
        print(f"{category:16s} {ro_rate:>11.3f} {null_rate:>16.3f} {display}{marker}")

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["category", "observed_count", "per_ro_window",
                         "null_count", "per_null_window", "enrichment"])
        for row in results:
            writer.writerow([row[0], row[1], round(row[2], 4), row[3],
                             round(row[4], 4), round(row[5], 3)
                             if row[5] != float("inf") else ""])
    print(f"\n[yazildi] {args.out}")
    print(f"[bilgi] {ro_windows} RO penceresi vs {window_count} kontrol penceresi")
    print("\nYorum: zenginlesme >1 = RO yaninda beklenenden SIK. ~1 = sadece")
    print("genomda yaygin oldugu icin gorunuyor, RO ile ozel bir iliskisi yok.")
    connection.close()


if __name__ == "__main__":
    sys.exit(main() or 0)
