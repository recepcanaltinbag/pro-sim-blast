"""Anotasyonsuz indirilmis genomlari YENIDEN ceker.

SORUN. `replicon` tablosundaki 2.975 kayit (%17,4) `no_annotation` durumunda:
dosya var ama icinde ozellik yok. Hepsi ayni sekilde ~5 KB'lik bir CON
iskeleti -- NCBI'nin "contig" kaydi, yani dizi baska kayitlara isaret eden bir
yonlendirme listesi. `rettype=gb` bu kayitlar icin iskeleti dondurur;
`rettype=gbwithparts` parcalari birlestirip tam kaydi verir.

NEDEN ONEMLI. Eksiklik rastgele degil. En cok etkilenen cinsler Streptomyces
(248), Pseudomonas (191) ve Burkholderia (114) -- yani Rieske oksijenaz
bakimindan EN ZENGIN olanlar. Yani uye sayilari sistematik olarak eksik ve
eksiklik tam da en cok uye beklenen yerlerde.

DOGRULANDI. `NZ_CP049045.1` (Pseudomonas sp. BIOMIG1BAC, kullanicinin kendi
makalesi doi:10.1128/MRA.00309-20) `replicon`'da var ama o genomdan tek bir
aday yok. `gbwithparts` ile yeniden cekildiginde 19,5 MB ve 7.119 CDS geliyor,
ve qxyA referans dizisi makalenin verdigi 6.844.288..6.845.439
koordinatlarinda birebir bulunuyor (WP_068587002.1).

BU SCRIPT NE YAPAR, NE YAPMAZ. Yalnizca INDIRIR ve dogrular. Veritabanina
dokunmaz, pipeline'i kosturmaz, sitedeki hicbir sayiyi degistirmez. Indirilen
dosyalar --out-dir altina yazilir (varsayilan harici diskte, cunku kok
bolumde yer yok). Yeniden calistirilabilir: gecerli bir dosya zaten varsa
atlanir, yarim kalan indirme .part uzantisinda kalir ve tamamlanmadan
yerine konmaz.

Kullanim:
    python3 refetch_genomes.py --email ADRES [--api-key ANAHTAR] --limit 20
    python3 refetch_genomes.py --email ADRES --api-key ANAHTAR
"""

import argparse
import concurrent.futures
import gzip
import http.client
import threading
import os
import re
import socket
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# NCBI'nin ilan ettigi sinir: anahtarsiz saniyede 3, anahtarla 10. Biraz
# altinda kaliniyor, cunku sinir asilinca IP engellenir ve bu is uzun surer.
RATE_NO_KEY = 3.0
RATE_WITH_KEY = 9.0

# Bir kaydin "gercek" sayilmasi icin en az bu kadar ozellik tasimasi gerekir.
# CON iskeleti tipik olarak tek bir source ozelligi tasir ve hic CDS tasimaz.
MIN_CDS = 1

FEATURE_RE = re.compile(rb"^     CDS             ", re.M)
CONTIG_ONLY_RE = re.compile(rb"^CONTIG\s", re.M)
LOCUS_RE = re.compile(rb"^LOCUS\s+(\S+)", re.M)


def targets(con, limit=None):
    """Anotasyonsuz replikonlar, en cok etkilenen cinsler once.

    Sira onemli: is yarida kesilirse en degerli kayitlar cekilmis olur.
    """
    rows = con.execute("""
        SELECT nucleotide_id, organism,
               (SELECT COUNT(*) FROM ro r WHERE r.nucleotide_id = p.nucleotide_id) ro_count
        FROM replicon p
        WHERE (cds_count IS NULL OR cds_count = 0)
          AND nucleotide_id IS NOT NULL AND trim(nucleotide_id) <> ''
        ORDER BY ro_count DESC, organism, nucleotide_id
    """).fetchall()
    return rows[:limit] if limit else rows


# OLCULDU: `gbwithparts` ile buyuk bir kayit istendiginde NCBI akisi zaman
# zaman yarida kapatiyor ve http.client.IncompleteRead firlatiyor. Bu sinif
# OSError'dan TUREMEZ (HTTPException + ValueError), yani alelade bir
# `except OSError` onu yakalamaz ve script ilk kayitta cokuyordu.
TRANSIENT = (urllib.error.URLError, urllib.error.HTTPError, OSError,
             socket.timeout, http.client.HTTPException)


def fetch(accession, email, api_key, tool="ROAR-DB-refetch", timeout=900):
    params = {"db": "nuccore", "id": accession, "rettype": "gbwithparts",
              "retmode": "text", "tool": tool, "email": email}
    if api_key:
        params["api_key"] = api_key
    url = EFETCH + "?" + urllib.parse.urlencode(params)
    # Parca parca okunur: tek seferlik read() koptugunda okunan her sey
    # kayboluyor, oysa IncompleteRead'in partial'i cogu zaman tam kayittir.
    chunks = []
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            while True:
                part = resp.read(1 << 20)
                if not part:
                    break
                chunks.append(part)
    except http.client.IncompleteRead as exc:
        if exc.partial:
            chunks.append(exc.partial)
        blob = b"".join(chunks)
        # Kayit bitis isaretini tasiyorsa kopma zararsizdir.
        if blob.rstrip().endswith(b"//"):
            return blob
        raise
    return b"".join(chunks)


def inspect(blob):
    """Kaydin gercekten anotasyonlu olup olmadigini soyler."""
    n_cds = len(FEATURE_RE.findall(blob))
    locus = LOCUS_RE.search(blob)
    return {
        "bytes": len(blob),
        "cds": n_cds,
        "locus": locus.group(1).decode() if locus else None,
        "still_a_contig_stub": bool(CONTIG_ONLY_RE.search(blob)) and n_cds < MIN_CDS,
        "usable": n_cds >= MIN_CDS,
    }


class RateLimiter:
    """Butun is parcaciklarinin PAYLASTIGI hiz sinirlayici.

    NCBI'nin siniri saniyedeki ISTEK sayisinda, es zamanli baglanti
    sayisinda degil. Bir kayit 20 saniyede indigi icin tek parcacikla
    saniyede 0,05 istek yapiliyor -- sinirin cok altinda ve is 20 saate
    yayiliyor. Birkac parcacik ayni kovadan jeton alirsa hem sinir
    asilmaz hem is birkac saate iner.
    """

    def __init__(self, per_second):
        self.interval = 1.0 / per_second
        self.lock = threading.Lock()
        self.next_at = 0.0

    def take(self):
        with self.lock:
            now = time.monotonic()
            wait = max(0.0, self.next_at - now)
            self.next_at = max(now, self.next_at) + self.interval
        if wait:
            time.sleep(wait)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out-dir", default="/media/lin-bio/back2/roar_gbk_refetch",
                    help="indirilen kayitlar buraya; kok bolumde yer yok")
    ap.add_argument("--email", required=True,
                    help="NCBI E-utilities iletisim adresi (zorunlu tutuluyor: "
                         "anonim toplu indirme IP engeline yol acabilir)")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                    help="varsa saniyede 10 istege cikar")
    ap.add_argument("--limit", type=int, help="yalnizca ilk N kayit (pilot icin)")
    ap.add_argument("--gzip", action="store_true", default=True,
                    help="diskte sikistirilmis sakla (varsayilan)")
    ap.add_argument("--no-gzip", dest="gzip", action="store_false")
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--workers", type=int, default=4,
                    help="es zamanli indirme. Hiz siniri ORTAK bir kovadan "
                         "uygulanir, yani parcacik sayisi istek hizini "
                         "artirmaz; yalnizca aktarim suresini ortuyor")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    con = sqlite3.connect(args.db)
    todo = targets(con, args.limit)
    delay = 1.0 / (RATE_WITH_KEY if args.api_key else RATE_NO_KEY)

    print(f"[hedef] {len(todo)} anotasyonsuz replikon")
    print(f"[dizin] {args.out_dir}")
    print(f"[hiz]   saniyede {'~9 (api anahtari var)' if args.api_key else '~3 (anahtar yok)'}")
    print()

    ext = ".gbk.gz" if args.gzip else ".gbk"
    limiter = RateLimiter(RATE_WITH_KEY if args.api_key else RATE_NO_KEY)
    opener = gzip.open if args.gzip else open
    counts = {"done": 0, "failed": 0, "skipped": 0, "stubs": 0,
              "bytes": 0, "cds": 0, "seen": 0}
    lock = threading.Lock()
    t0 = time.time()

    def handle(target):
        acc, organism, _ro_count = target
        path = os.path.join(args.out_dir, acc + ext)
        if os.path.exists(path) and os.path.getsize(path) > 2048:
            return "skipped", acc, None

        blob = None
        for attempt in range(1, args.retries + 1):
            limiter.take()
            try:
                blob = fetch(acc, args.email, args.api_key)
                break
            except TRANSIENT as exc:
                wait = min(60, 2 ** attempt)
                print(f"  [yeniden] {acc} deneme {attempt}/{args.retries}: "
                      f"{type(exc).__name__} -- {wait} s", flush=True)
                time.sleep(wait)
        if blob is None:
            return "failed", acc, None
        if not blob.lstrip().startswith(b"LOCUS"):
            head = blob[:120].decode("utf-8", "replace").replace("\n", " ")
            print(f"  [GENBANK DEGIL] {acc}: {head}", flush=True)
            return "failed", acc, None

        info = inspect(blob)
        tmp = path + f".part{threading.get_ident()}"
        with opener(tmp, "wb") as fh:
            fh.write(blob)
        os.replace(tmp, path)   # yarim dosya asla yerine konmaz
        if not info["usable"]:
            print(f"  [HALA ISKELET] {acc} {organism} "
                  f"{info['bytes']} bayt, {info['cds']} CDS", flush=True)
            return "stub", acc, info
        return "done", acc, info

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for kind, acc, info in pool.map(handle, todo):
            with lock:
                counts["seen"] += 1
                if kind == "done":
                    counts["done"] += 1
                    counts["bytes"] += info["bytes"]
                    counts["cds"] += info["cds"]
                elif kind == "stub":
                    counts["stubs"] += 1
                elif kind == "skipped":
                    counts["skipped"] += 1
                else:
                    counts["failed"] += 1
                i = counts["seen"]
            if i % 25 == 0 or i == len(todo):
                el = time.time() - t0
                rate = i / el if el else 0
                left = (len(todo) - i) / rate if rate else 0
                print(f"[{i}/{len(todo)}] tam={counts['done']} "
                      f"iskelet={counts['stubs']} atlanan={counts['skipped']} "
                      f"hata={counts['failed']}  "
                      f"{counts['bytes'] / 1e9:.1f} GB  "
                      f"kalan ~{left / 60:.0f} dk", flush=True)

    done, stubs, skipped, failed = (counts["done"], counts["stubs"],
                                    counts["skipped"], counts["failed"])
    total_bytes, total_cds = counts["bytes"], counts["cds"]

    print()
    print(f"[bitti] {done} kayit anotasyonlu olarak indirildi, "
          f"{stubs} tanesi yeniden cekildiginde de bos, "
          f"{skipped} zaten vardi, {failed} basarisiz")
    if done:
        print(f"        toplam {total_bytes / 1e9:.1f} GB (sikistirilmamis), "
              f"{total_cds:,} CDS, kayit basina ortalama "
              f"{total_bytes / done / 1e6:.1f} MB ve {total_cds // done:,} CDS")
    print(f"        dizin: {args.out_dir}")
    print()
    print("Bu script veritabanina DOKUNMADI. Sayilarin guncellenmesi icin "
          "pipeline'in yeniden kosturulmasi gerekir ve bu sitedeki her sayiyi "
          "degistirir.")


if __name__ == "__main__":
    sys.exit(main())
