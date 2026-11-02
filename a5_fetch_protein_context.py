"""A5: protein duzeyinde duran RO alpha'lar icin GENOM KAYDINI ceker.

NEDEN. Bugun `ro_union`'da 73.299 tekil dizi var, bunlarin yalnizca 15.960'i
genomik baglam tasiyor. Geri kalan 57.339'un operonu, komsulugu, replikonu,
izolasyon kaynagi YOK -- yani sitenin konuma bagli her olcumu onlari
gormuyor. Ama 18.636 tanesi icin UniProt capraz referansi bir EMBL/GenBank
NUKLEOTID kaydi veriyor; o kayit indirilirse baglam KURULABILIR.

OLCUM (2026-11-02): 18.636 dizi -> 13.079 TEKIL nukleotid kaydi. Yani kayit
basina ortalama 1,4 dizi; en yogun kayit CP000699 tek basina 52 dizi tasiyor.
Indirme biriminin dizi degil KAYIT olmasi isi 18.636'dan 13.079'a indiriyor.

BU SCRIPT NE YAPAR, NE YAPMAZ. Yalnizca indirir ve her kaydin sonucunu bir
siciIe yazar. Veritabanina DOKUNMAZ, pipeline kosturmaz, sitedeki hicbir
sayiyi degistirmez. Baglamin gercekten kurulmasi ayri bir adim
(`a5_link_protein_context.py`) ve orada bir DIZI ESLESME SARTI var: indirilen
kayitta cevirisinin md5'i bizim dizimizle BIREBIR ayni olan bir CDS yoksa o
dizi baglamsiz kalir. Capraz referans "ayni kayitta gecer" demektir, "ayni
proteindir" demez; bu ayrim yapilmazsa yanlis genoma bagli operon yayimlanir.

Yeniden calistirilabilir: gecerli dosya varsa atlanir, yarim indirme .part
uzantisinda kalir ve tamamlanmadan yerine konmaz.

Kullanim:
    python3 a5_fetch_protein_context.py --limit 30        # pilot
    python3 a5_fetch_protein_context.py                   # tamami
"""

import argparse
import concurrent.futures
import gzip
import json
import os
import sqlite3
import sys
import threading
import time

from refetch_genomes import (RATE_NO_KEY, RATE_WITH_KEY, TRANSIENT, RateLimiter,
                             fetch, inspect)

# Projenin NCBI iletisim adresi (verify_reference_sources.py ile ayni).
DEFAULT_EMAIL = "recocanextra@gmail.com"


def targets(con, limit=None):
    """Indirilecek tekil nukleotid kayitlari, tasidiklari dizi sayisina gore.

    Bir UniProt girisi birden fazla EMBL kaydina isaret edebiliyor (alanda
    noktali virgulle ayrilmis); hepsi ayri hedef sayilir, cunku hangisinin
    dogru CDS'i tasidigini ancak indirip dizi esleyerek anlariz.
    """
    rows = con.execute("""
        SELECT DISTINCT u.seq_md5, p.embl_nuc
        FROM ro_union u
        JOIN ro_union_source s ON s.seq_md5 = u.seq_md5 AND s.source = 'uniprot'
        JOIN ro_protein_status p ON p.accession = s.accession
        WHERE u.has_genomic_context = 0
          AND p.embl_nuc IS NOT NULL AND TRIM(p.embl_nuc) <> ''
    """).fetchall()
    per_acc = {}
    for md5, field in rows:
        for tok in str(field).replace(";", " ").replace(",", " ").split():
            acc = tok.strip()
            if acc:
                per_acc.setdefault(acc, set()).add(md5)
    order = sorted(per_acc.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    out = [(acc, len(md5s)) for acc, md5s in order]
    return out[:limit] if limit else out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out-dir", default="/media/lin-bio/back2/roar_gbk_a5",
                    help="indirilen kayitlar buraya; kok bolumde yer yok")
    ap.add_argument("--ledger", default=None,
                    help="sonuc sicili (varsayilan: <out-dir>/ledger.jsonl)")
    ap.add_argument("--email", default=DEFAULT_EMAIL)
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    ap.add_argument("--limit", type=int, help="yalnizca ilk N kayit (pilot)")
    ap.add_argument("--retries", type=int, default=4)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    ledger_path = args.ledger or os.path.join(args.out_dir, "ledger.jsonl")
    con = sqlite3.connect(args.db)
    todo = targets(con, args.limit)

    print(f"[hedef] {len(todo)} tekil nukleotid kaydi, "
          f"{sum(n for _, n in todo):,} dizi-kayit bagi")
    print(f"[dizin] {args.out_dir}")
    print(f"[sicil] {ledger_path}")
    print(f"[hiz]   saniyede {'~9 (api anahtari var)' if args.api_key else '~3 (anahtar yok)'}")
    print(flush=True)

    limiter = RateLimiter(RATE_WITH_KEY if args.api_key else RATE_NO_KEY)
    counts = {"done": 0, "failed": 0, "skipped": 0, "stubs": 0, "notfound": 0,
              "bytes": 0, "cds": 0, "seen": 0}
    lock = threading.Lock()
    ledger = open(ledger_path, "a", buffering=1)
    t0 = time.time()

    def handle(target):
        acc, n_seq = target
        path = os.path.join(args.out_dir, acc + ".gbk.gz")
        if os.path.exists(path) and os.path.getsize(path) > 512:
            return "skipped", acc, n_seq, None

        blob, last = None, None
        for attempt in range(1, args.retries + 1):
            limiter.take()
            try:
                blob = fetch(acc, args.email, args.api_key,
                             tool="ROAR-DB-a5")
                break
            except TRANSIENT as exc:
                last = f"{type(exc).__name__}: {exc}"
                # 400/404 kalici: accession gecersiz ya da cekilemez.
                code = getattr(exc, "code", None)
                if code in (400, 404, 410):
                    return "notfound", acc, n_seq, {"error": last}
                time.sleep(min(60, 2 ** attempt))
        if blob is None:
            return "failed", acc, n_seq, {"error": last}
        if not blob.lstrip().startswith(b"LOCUS"):
            head = blob[:120].decode("utf-8", "replace").replace("\n", " ")
            return "notfound", acc, n_seq, {"error": f"GenBank degil: {head}"}

        info = inspect(blob)
        tmp = path + f".part{threading.get_ident()}"
        with gzip.open(tmp, "wb") as fh:
            fh.write(blob)
        os.replace(tmp, path)
        return ("done" if info["usable"] else "stub"), acc, n_seq, info

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for kind, acc, n_seq, info in pool.map(handle, todo):
            with lock:
                counts["seen"] += 1
                counts[{"done": "done", "stub": "stubs", "skipped": "skipped",
                        "notfound": "notfound", "failed": "failed"}[kind]] += 1
                if info and info.get("usable"):
                    counts["bytes"] += info["bytes"]
                    counts["cds"] += info["cds"]
                i = counts["seen"]
                ledger.write(json.dumps({"accession": acc, "outcome": kind,
                                         "sequences": n_seq,
                                         "bytes": (info or {}).get("bytes"),
                                         "cds": (info or {}).get("cds"),
                                         "error": (info or {}).get("error")}) + "\n")
            if i % 50 == 0 or i == len(todo):
                el = time.time() - t0
                rate = i / el if el else 0
                left = (len(todo) - i) / rate if rate else 0
                print(f"[{i}/{len(todo)}] tam={counts['done']} "
                      f"iskelet={counts['stubs']} yok={counts['notfound']} "
                      f"atlanan={counts['skipped']} hata={counts['failed']}  "
                      f"{counts['bytes'] / 1e9:.2f} GB  "
                      f"kalan ~{left / 60:.0f} dk", flush=True)
    ledger.close()

    print()
    print(f"[bitti] {counts['done']} anotasyonlu, {counts['stubs']} iskelet, "
          f"{counts['notfound']} cekilemedi (kalici), "
          f"{counts['skipped']} zaten vardi, {counts['failed']} gecici hata")
    if counts["done"]:
        print(f"        {counts['bytes'] / 1e9:.2f} GB, {counts['cds']:,} CDS, "
              f"kayit basina ~{counts['bytes'] / counts['done'] / 1e3:.0f} KB "
              f"ve {counts['cds'] // counts['done']:,} CDS")
    if counts["failed"]:
        print(f"        UYARI: {counts['failed']} kayit GECICI hatayla dustu. "
              f"Bunlar eksik, 'yok' DEGIL -- script yeniden kosturulmali.")
    print()
    print("Bu script veritabanina DOKUNMADI. Baglamin kurulmasi icin "
          "a5_link_protein_context.py kosturulmali.")


if __name__ == "__main__":
    sys.exit(main())
