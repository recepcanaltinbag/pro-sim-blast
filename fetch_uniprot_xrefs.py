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
import json
import os
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

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


def page_through(job, size):
    """Bir isin sonuclarini sayfalayarak toplar. Dondurur: (satirlar, sayfa)."""
    url = RESULTS.format(job=job, size=size)
    lines, pages, header = [], 0, None
    while url:
        with open_url(url) as res:
            text, link = res.read().decode("utf-8", "replace"), res.headers.get("Link")
        rows = text.splitlines()
        if rows:
            if header is None:
                header = rows[0]
            lines.extend(r for r in rows[1:] if r.strip())
        pages += 1
        url = next_link(link)
    return header, lines, pages


def fetch_batch(ids, path, size):
    """Bir partiyi ceker. SIFIR sonuc ONBELLEGE YAZILMAZ.

    Eksik kimlik (silinmis accession) olabilir ve bu bir hata degil, ama
    BEYAN edilir. Sifir sonuc ise hatadir: parti onbellege girmez ve bir
    sonraki kosuda yeniden denenir -- cunku onbellege giren bos bir parti,
    1.500 accession'i sessizce "InterPro kaydi yok" diye yayimlar.
    """
    job = submit(ids)
    wait_for_job(job)
    header, lines, pages = page_through(job, size)
    seen = {ln.split("\t")[0].strip() for ln in lines if ln.strip()}
    seen.discard("")
    if not seen:
        raise RuntimeError(f"{len(ids)} kimlik icin SIFIR satir dondu "
                           f"({pages} sayfa); onbellege yazilmadi")

    #  KESIK SONUC DA ONBELLEGE GIRMEZ -- ve bu, sifir kontrolunden AYRI
    #  bir kontrol olmak zorundaydi. Olculdu: 19 partinin 14'u 1500/1500
    #  dondu, BESI 1000/1500 -- yani sayfali uc nokta `Link: rel="next"`i
    #  iki sayfadan sonra ARALIKLI olarak dusuruyor. Sifiri reddeden
    #  koruma bunlari KABUL ediyordu, ve kesik bir parti onbellege girince
    #  500 accession sessizce "InterPro kaydi yok" diye yayimlanirdi.
    #
    #  AYIRT EDICI: gercek bir eksiklik (silinmis accession) rastgele bir
    #  sayidir; KESILME sayfa boyunun TAM KATIDIR. O yuzden yuvarlak bir
    #  eksiklik kesilme sayilip yeniden denenir, yuvarlak olmayan eksiklik
    #  beyan edilip kabul edilir.
    if len(seen) < len(ids) and len(seen) % size == 0:
        raise RuntimeError(
            f"KESIK: {len(seen)}/{len(ids)} kimlik, {pages} sayfa "
            f"(sayfa boyu {size}'in tam kati); onbellege yazilmadi")
    payload = ("\n".join([header or ""] + lines) + "\n").encode()
    tmp = path + ".partial"
    with gzip.open(tmp, "wb") as handle:
        handle.write(payload)
    os.replace(tmp, path)
    return len(seen), pages


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
