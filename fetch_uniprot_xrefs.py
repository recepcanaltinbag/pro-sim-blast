"""UniProt'tan InterPro / KEGG / MetaCyc(BioCyc) tanimlayicilarini ceker.

NEDEN BU BETIK VAR. ROADMAP A6: "bir toparlama veritabaninin isi ayni
molekulun her yerdeki adini bir arada tutmaktir." `cross_references.py`
dort tanimlayiciyi yayimliyor (UniProt, RefSeq, EMBL nukleotit, EC) ama
dordu de `ro_protein_status`'ta ZATEN duran sutunlardan geliyor. InterPro,
KEGG ve MetaCyc hic cekilmemisti.

TASARIM OLCUMLE SECILDI, TAHMINLE DEGIL -- ve ilk iki tasarim YANLISTI.
Dokumanlar duz metin olarak alinamadigi icin servisin davranisi kucuk
ollcekli problarla olculdu, ve su uc olgu tasarimi belirledi:

  1. GONDERIM BOYUTU SINIRLI, ve sinir 100.000 kimlik DEGIL. Olculdu:
     1.500 kimlik (gövde 19,5 KB) TAM calisiyor; 5.000 (65 KB) ve 20.000
     (260 KB) SIFIR sonuc donduruyor ve is BIR SANIYEDE "bitiyor" --
     yani gonderimde reddediliyor ama yine de bir jobId veriyor. 73.699
     kimlikli tek is de (810 KB) ayni sekilde sifir dondu. Dolayisiyla
     parti boyu OLCULEN degerdir: 1.500.
  2. `/results/stream/` uc noktasi KESIYOR: 5.000 kimlikli bir isten tam
     1000 satir dondu, 3.699 kimlikli partiden de 1000. Yuvarlak sayi bir
     sayfa sinirinin imzasi. Kullanilmiyor.
  3. Sayfali `/results/` + `Link: rel="next"` TAM donuyor: 1.500 kimlikle
     3 sayfa, 1.500 satir, 1.500 tekil kimlik.
  4. ISTEK HIZI ONEMLI. 15 partiyi aralik vermeden gondermek durum
     yoklamalarina HTTP 400 getirdi; 400 akista degil DURUM uc
     noktasinda geliyordu, yani sorun kodda degil hizdaydi. Partiler TEK
     TEK, aralarinda bekleyerek gonderilir.

BOS VE KESIK SONUC BASARI DEGILDIR. Ilk yazim bir alanin YOKLUGUNU
("jobStatus" yok) "is bitti" saymisti; hata yanitlari da o alani
tasimadigi icin is hazir olmadan akis cekildi ve yalnizca BASLIK satiri
olan dosyalar "fetched" diye bildirildi. Simdi "bitti"nin isareti ya
`FINISHED`, ya sonuclara YONLENME, ya da yanitta `results` olmasi; ve
cekim sonunda GONDERILEN ile DONEN kimlik sayisi karsilastirilip
onbellegin kendi meta dosyasina yaziliyor -- cunku eksik bir cekim "bu
dizinin InterPro kaydi yok" gibi gorunur, yani eksikligini GIZLEYEN
yonde hata uretir.

DEVAM EDEBILIR. Her sayfa ayri dosyaya yaziliyor ve siradaki cursor bir
durum dosyasinda duruyor; betik yeniden kosturuldugunda kaldigi sayfadan
devam eder.

TLS DOGRULAMASI ACIK. Bu depodaki eski `download_pfam.py`
`ssl._create_unverified_context()` kullaniyor; buraya kopyalanmadi.

NE YAPMAZ. Tanimlayiciyi YORUMLAMAZ. InterPro bir domain IDDIASI, KEGG ve
MetaCyc birer yol IDDIASI; hangisinin ne kadar guvenilir oldugu
`cross_references.py`nin isi. Burada yalnizca cekim ve onbellek var.

Cikti: analysis_out/uniprot_xrefs/page_NNNN.tsv.gz + _fetch.json
"""

import argparse
import glob
import gzip
import hashlib
import json
import os
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

#  DUZ ARAMA UC NOKTASI -- ID-mapping yerine. Olculdu (2026-11-13):
#  100 accession'lik bir `accession:(A OR B ...)` sorgusu 100/100 kimlik
#  donduruyor ve SAYFALAMA bile gerekmiyor; 200 accession'da URL uzunlugu
#  HTTP 400 veriyor. Bu yol ID-mapping'in BUTUN basarisizlik yuzeyini
#  ortadan kaldiriyor: asenkron is yok (dolayisiyla sessiz-sifir yok),
#  cursor yok (dolayisiyla kesilme yok), kuyruk kotasi yok. ID-mapping ile
#  73.699 accession'in ancak 27.000'i cekilebildi; gerisinde servis isi
#  kabul edip BOS donduruyordu.
SEARCH = "https://rest.uniprot.org/uniprotkb/search"
CHUNK = 100

RUN = "https://rest.uniprot.org/idmapping/run"
STATUS = "https://rest.uniprot.org/idmapping/status/{job}"
RESULTS = ("https://rest.uniprot.org/idmapping/uniprotkb/results/{job}"
           "?format=tsv&fields=accession,xref_interpro,xref_kegg,xref_biocyc"
           "&size={size}")
UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
OUT_DIR = "analysis_out/uniprot_xrefs"
PAGE = 500
#  OLCULEN gonderim siniri: 1.500 tam calisiyor, 5.000 sifir donduruyor.
BATCH = 1500


def open_url(url, data=None, tries=6):
    """Yeniden denemeli istek. `Retry-After` varsa sunucunun suresi kullanilir."""
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=data,
                                         headers={"User-Agent": UA})
            return urllib.request.urlopen(req, timeout=300)
        except urllib.error.HTTPError as err:
            if attempt == tries - 1:
                raise
            wait = err.headers.get("Retry-After") if err.headers else None
            time.sleep(int(wait) if wait and wait.isdigit()
                       else 20 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError):
            if attempt == tries - 1:
                raise
            time.sleep(20 * (attempt + 1))


def submit(ids):
    body = urllib.parse.urlencode({"from": "UniProtKB_AC-ID",
                                   "to": "UniProtKB",
                                   "ids": ",".join(ids)}).encode()
    with open_url(RUN, body) as res:
        return json.loads(res.read())["jobId"]


def wait_for_job(job, poll=5.0, limit=3600):
    """Is BITENE KADAR bekler. Bir alanin YOKLUGUNU basari SAYMAZ."""
    waited = 0.0
    while True:
        with open_url(STATUS.format(job=job)) as res:
            url, payload = res.geturl(), json.loads(res.read())
        if payload.get("jobStatus") == "ERROR" or payload.get("errors"):
            raise RuntimeError(f"ID-mapping isi hata verdi: {payload}")
        if (payload.get("jobStatus") == "FINISHED" or payload.get("results")
                or "/results/" in url):
            return
        if waited > limit:
            raise TimeoutError(f"is {limit} s icinde bitmedi: {job}")
        time.sleep(poll)
        waited += poll


def next_link(header):
    if not header or 'rel="next"' not in header:
        return None
    return header.split(";")[0].strip().lstrip("<").rstrip(">")


def accessions(db):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    rows = con.execute("SELECT DISTINCT accession FROM ro_union_source "
                       "WHERE source = 'uniprot' ORDER BY accession").fetchall()
    con.close()
    return [r[0] for r in rows]


def accessions_in_file(path):
    """Bir onbellek dosyasindaki ACCESSION'lar. Sutun BASLIKTAN bulunur.

    NEDEN BASLIKTAN. Dizinde IKI duzen bir arada: ID-mapping doneminden
    kalan dosyalar `From, Entry, InterPro, ...`, arama dosyalari ise
    `Entry, InterPro, ...`. Ilk yazim ikisini `parts[:2]` ile idare etmeye
    calisiyordu ve arama dosyalarinda IKINCI SUTUNU (InterPro) accession
    sanip `have` kumesine IPR kimlikleri katiyordu -- sonucta "10.000
    accession eksik" diye YUVARLAK bir sayi uretti, ki yuvarlak sayi bu
    projede kesilme/sayim hatasi imzasi. Sutunu addan bulmak varsayimi
    tamamen kaldiriyor.
    """
    out = set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        cols = [c.strip().lower() for c in header.rstrip("\n").split("\t")]
        if "entry" not in cols:
            return out
        index = cols.index("entry")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > index and parts[index].strip():
                out.add(parts[index].strip())
    return out


def fetch_chunk(ids, path, pause=0.4):
    """Bir parcayi DUZ ARAMA ile ceker. Eksik/bos sonuc ONBELLEGE GIRMEZ."""
    query = "accession:(" + " OR ".join(ids) + ")"
    url = SEARCH + "?" + urllib.parse.urlencode({
        "query": query, "format": "tsv", "size": 500,
        "fields": "accession,xref_interpro,xref_kegg,xref_biocyc"})
    with open_url(url) as res:
        text = res.read().decode("utf-8", "replace")
    rows = text.splitlines()
    header = rows[0] if rows else ""
    lines = [r for r in rows[1:] if r.strip()]
    seen = {ln.split("\t")[0].strip() for ln in lines}
    seen.discard("")
    if not seen:
        raise RuntimeError(f"{len(ids)} accession icin SIFIR satir dondu; "
                           "onbellege yazilmadi")
    payload = ("\n".join([header] + lines) + "\n").encode()
    tmp = path + ".partial"
    with gzip.open(tmp, "wb") as handle:
        handle.write(payload)
    os.replace(tmp, path)
    time.sleep(pause)
    return len(seen)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default=OUT_DIR)
    parser.add_argument("--chunk", type=int, default=CHUNK,
                        help="bir sorgudaki accession sayisi (OLCULEN sinir "
                             "100-200 arasi)")
    parser.add_argument("--pause", type=float, default=0.4)
    args = parser.parse_args()

    ids = accessions(args.db)
    os.makedirs(args.out_dir, exist_ok=True)

    #  ID-mapping doneminden kalan `page_*` dosyalari GECERLI VERI; hangi
    #  accession'larin zaten elimizde oldugu onlardan okunur ve o
    #  accession'lar yeniden istenmez.
    have = set()
    for path in sorted(glob.glob(os.path.join(args.out_dir, "*.tsv.gz"))):
        have |= accessions_in_file(path)
    todo = [a for a in ids if a not in have]
    chunks = [todo[i:i + args.chunk] for i in range(0, len(todo), args.chunk)]
    print(f"[accession] {len(ids)} toplam, {len(have & set(ids))} elimizde, "
          f"{len(todo)} eksik -> {len(chunks)} sorgu", flush=True)

    #  DOSYA ADI ICERIKTEN TURETILIR, KONUMDAN DEGIL.
    #
    #  Ilk yazim `chunk_{n:05d}` kullaniyordu, yani adi `todo` listesindeki
    #  KONUM belirliyordu -- ama `todo` her kosuda degisiyor (elimizde olan
    #  accession'lar cikarildigi icin). Sonuc olculdu: ikinci kosu, birinci
    #  kosunun 1-100 numarali dosyalarini "zaten var" diye ATLADI, oysa o
    #  dosyalar BASKA accession'lari tasiyordu. Tam 100 parca x 100 = 10.000
    #  accession sessizce cekilmemis kaldi ve sayim "10.001 eksik" dedi.
    #
    #  Ad artik parcanin KENDI iceriginin ozeti: ayni accession kumesi ayni
    #  dosyaya gider, degisen bir liste eski dosyalari gecersiz kilmaz.
    done = failed = 0
    for n, chunk in enumerate(chunks, 1):
        tag = hashlib.sha1(",".join(chunk).encode()).hexdigest()[:12]
        path = os.path.join(args.out_dir, f"chunk_{tag}.tsv.gz")
        if os.path.exists(path) and os.path.getsize(path) > 100:
            done += 1
            continue
        try:
            seen = fetch_chunk(chunk, path, args.pause)
        except Exception as err:                      # noqa: BLE001
            print(f"[sorgu {n}/{len(chunks)}] DUSTU: {err}",
                  file=sys.stderr, flush=True)
            failed += 1
            continue
        done += 1
        if n % 50 == 0 or n == len(chunks):
            print(f"[sorgu {n}/{len(chunks)}] {seen}/{len(chunk)}", flush=True)

    seen = set()
    for path in sorted(glob.glob(os.path.join(args.out_dir, "*.tsv.gz"))):
        seen |= accessions_in_file(path)
    covered = len(seen & set(ids))
    state = {
        "submitted": len(ids),
        "returned_ids": covered,
        "missing_ids": len(ids) - covered,
        "chunks_ready": done,
        "chunks_failed": failed,
        "chunk_size": args.chunk,
        "finished": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "endpoint": ("plain /uniprotkb/search with accession:(A OR B ...), "
                     "100 accessions per query -- measured: 100 returns "
                     "100/100 with no paging, 200 exceeds the URL limit and "
                     "returns HTTP 400. ID-mapping was abandoned: it accepts "
                     "a job and returns zero results once many jobs have been "
                     "submitted, which is a silent refusal"),
    }
    with open(os.path.join(args.out_dir, "_fetch.json"), "w",
              encoding="utf-8") as handle:
        json.dump(state, handle, indent=1)
    print(f"[bitti] {done} sorgu, {failed} dustu; {covered}/{len(ids)} "
          f"accession elimizde, {len(ids) - covered} eksik")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default=OUT_DIR)
    parser.add_argument("--size", type=int, default=PAGE,
                        help="sonuc sayfasi boyu")
    parser.add_argument("--batch", type=int, default=BATCH,
                        help="bir isteki kimlik sayisi (OLCULEN sinir 1500)")
    parser.add_argument("--gap", type=float, default=20.0,
                        help="isler arasi bekleme; 0 vermek 400 firtinasi "
                             "uretir (olculdu)")
    args = parser.parse_args()

    ids = accessions(args.db)
    os.makedirs(args.out_dir, exist_ok=True)
    batches = [ids[i:i + args.batch] for i in range(0, len(ids), args.batch)]
    print(f"[accession] {len(ids)} -> {len(batches)} parti "
          f"(parti {args.batch}, sayfa {args.size}, aralik {args.gap}s)",
          flush=True)

    done = failed = total_seen = 0
    for n, chunk in enumerate(batches, 1):
        path = os.path.join(args.out_dir, f"page_{n:04d}.tsv.gz")
        if os.path.exists(path) and os.path.getsize(path) > 100:
            done += 1
            continue
        try:
            seen, pages = fetch_batch(chunk, path, args.size)
        except Exception as err:                      # noqa: BLE001
            print(f"[parti {n}/{len(batches)}] DUSTU: {err}",
                  file=sys.stderr, flush=True)
            failed += 1
            time.sleep(args.gap)
            continue
        done += 1
        total_seen += seen
        if n % 5 == 0 or n == len(batches):
            print(f"[parti {n}/{len(batches)}] {seen}/{len(chunk)} kimlik, "
                  f"{pages} sayfa", flush=True)
        time.sleep(args.gap)

    #  TAMLIK, diskteki dosyalardan yeniden sayilir -- betigin kendi
    #  sayacindan degil (bir onceki kosunun partileri de sayilmali).
    seen = set()
    for path in sorted(glob.glob(os.path.join(args.out_dir, "page_*.tsv.gz"))):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            fh.readline()
            for line in fh:
                first = line.split("\t")[0].strip()
                if first:
                    seen.add(first)
    state = {
        "submitted": len(ids),
        "returned_ids": len(seen),
        "missing_ids": len(ids) - len(seen),
        "batches_ready": done,
        "batches_failed": failed,
        "batch_size": args.batch,
        "finished": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "endpoint": ("paged /idmapping/uniprotkb/results with Link rel=next. "
                     "The /results/stream variant was measured to cap at 1000 "
                     "rows, and a submission above about 1500 ids returns zero "
                     "results while still handing back a job id, so the batch "
                     "size is a measured value rather than the documented "
                     "maximum"),
    }
    with open(os.path.join(args.out_dir, "_fetch.json"), "w",
              encoding="utf-8") as handle:
        json.dump(state, handle, indent=1)
    print(f"[bitti] {done} parti hazir, {failed} dustu; gonderilen {len(ids)}, "
          f"donen {len(seen)}, eksik {len(ids) - len(seen)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
