"""`ingest_gap.csv`deki nukleotid kayitlarini indirir.

NEDEN. 135.28'de olculdu: kimligi BELLI ama hic indirilmemis 12.312
nukleotid kaydi var ve her biri kapimizi gecmis, tiplenmis en az bir
Rieske oksijenazi tasiyor. Bu betik o listeden secilen kayitlari ceker.

DOGRULAMA INDIRMEDEN ONCE YAPILDI (ROADMAP 135.28): ornek kayit
`FNTJ01000003` cekildi ve icinde qxyA referansiyla **%97,4** ayni 383 aa
bir ceviri bulundu -- ustelik yaninda toluat dioksijenazin alfa, beta ve
elektron transfer bilesenleriyle birlikte, ve ORIGIN dizisiyle. Yani
kayitlarin aradigimiz seyi tasidigi VARSAYILMADI, olculdu.

HER KAYIT YAZILMADAN ONCE SINANIR: LOCUS satiri ve en az bir CDS olmali.
Yarim ya da hatali bir indirme `gbk_files/`e girerse, pipeline onu
"anotasyonsuz replikon" diye okur ve sessizce paydaya katar.

DEVAM EDEBILIR: zaten inmis kayitlar atlanir.

Cikti: gbk_files/<accession>.gbk.gz
"""

import argparse
import csv
import gzip
import os
import re
import sys
import time
import urllib.error
import urllib.request

EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text")
UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
VERSION = re.compile(r"^VERSION\s+(\S+)", re.MULTILINE)


def fetch(accession, tries=4):
    for attempt in range(tries):
        try:
            request = urllib.request.Request(EFETCH.format(acc=accession),
                                             headers={"User-Agent": UA})
            with urllib.request.urlopen(request, timeout=300) as response:
                return response.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as err:
            if attempt == tries - 1:
                raise
            wait = err.headers.get("Retry-After") if err.headers else None
            time.sleep(int(wait) if wait and wait.isdigit()
                       else 10 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError):
            if attempt == tries - 1:
                raise
            time.sleep(10 * (attempt + 1))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--list", default="analysis_out/ingest_gap.csv")
    parser.add_argument("--gbk", default="gbk_files")
    parser.add_argument("--type", default="",
                        help="yalnizca bu tipi kazandiran kayitlar")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--pause", type=float, default=0.4)
    args = parser.parse_args()

    with open(args.list, newline="", encoding="utf-8") as handle:
        rows = [r for r in csv.DictReader(handle)
                if not args.type or args.type in r["types"]]
    rows.sort(key=lambda r: -int(r["sequences_gained"]))
    if args.limit:
        rows = rows[:args.limit]
    print(f"[liste] {len(rows)} kayit"
          + (f" ({args.type} kazandiranlar)" if args.type else ""), flush=True)

    os.makedirs(args.gbk, exist_ok=True)
    done = skipped = failed = 0
    for n, row in enumerate(rows, 1):
        accession = row["nucleotide_accession"]
        existing = [f for f in os.listdir(args.gbk)
                    if f.startswith(accession + ".")]
        if existing:
            skipped += 1
            continue
        try:
            text = fetch(accession)
        except Exception as err:                      # noqa: BLE001
            print(f"[{accession}] DUSTU: {err}", file=sys.stderr, flush=True)
            failed += 1
            continue
        #  SINAMA: yarim bir kayit pipeline'a "anotasyonsuz replikon" diye
        #  girer ve paydayi sessizce bozar.
        if "LOCUS" not in text or "     CDS  " not in text:
            print(f"[{accession}] GECERSIZ: LOCUS ya da CDS yok, yazilmadi",
                  file=sys.stderr, flush=True)
            failed += 1
            continue
        match = VERSION.search(text)
        name = match.group(1) if match else accession
        path = os.path.join(args.gbk, f"{name}.gbk.gz")
        tmp = path + ".partial"
        with gzip.open(tmp, "wt", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(tmp, path)
        done += 1
        if n % 25 == 0 or n == len(rows):
            print(f"   {n}/{len(rows)}  inen {done}, atlanan {skipped}, "
                  f"dusen {failed}", flush=True)
        time.sleep(args.pause)

    print(f"[bitti] inen {done}, atlanan {skipped}, dusen {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
