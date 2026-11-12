"""Hangi replikonun NUKLEOTID DIZISI elimizde?

NEDEN BU BETIK VAR. ROADMAP 124.3'un kalan yarisi -- sigma70 kutusu
zenginlesmesini uzunluk ve GC eslenmis bir null'a karsi olcmek --
intergenik DNA gerektiriyor, ve madde "gbk_files/ten cikarilmasi
gerekiyor" diyor. O varsayim SINANMAMISTI.

Olculdu: `gbk_files/` kayitlarinin COGU CON tipi, yani anotasyonu ve
protein cevirilerini tasiyor ama DIZIYI TASIMIYOR; dosya
`CONTIG join(...)` ile bitiyor ve 184.457 bp beyan eden bir kayit 231 KB
olabiliyor. Dolayisiyla intergenik DNA her replikon icin cikarilamaz.

Bu betik kimde dizi oldugunu SAYAR ve LISTELER, cunku zenginlesme
olcumunun paydasi tam olarak bu: dizisi olmayan bir replikonun yukari
akisi okunamaz, ve okunamayan bolgeyi paydada tutmak "motif yok" gibi
gorunur -- bu projenin tekrar eden hata bicimi.

Cikti: analysis_out/replicon_sequence.csv  (nucleotide_id, has_sequence,
       declared_bp, file_kind)
"""

import argparse
import csv
import gzip
import os
import re
import sqlite3
import sys
from collections import Counter

LOCUS_BP = re.compile(r"^LOCUS\s+\S+\s+(\d+)\s+bp", re.MULTILINE)


def open_record(gbk_dir, nucleotide_id):
    """Kaydin dosyasini acar; dizin duz ve gz bicimleri bir arada tutuyor."""
    plain = os.path.join(gbk_dir, nucleotide_id + ".gbk")
    if os.path.exists(plain):
        return open(plain, "rt", encoding="utf-8", errors="replace"), "plain"
    packed = plain + ".gz"
    if os.path.exists(packed):
        return gzip.open(packed, "rt", encoding="utf-8",
                         errors="replace"), "gzip"
    return None, "missing"


def scan(handle):
    """(dizi var mi, beyan edilen bp). ORIGIN satirini ARAR, dizi OKUMAZ."""
    declared, has_origin = None, False
    for line in handle:
        if declared is None and line.startswith("LOCUS"):
            match = LOCUS_BP.match(line)
            if match:
                declared = int(match.group(1))
        if line.startswith("ORIGIN"):
            has_origin = True
            break
    return has_origin, declared


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--gbk", default="gbk_files")
    parser.add_argument("--out", default="analysis_out/replicon_sequence.csv")
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    nucs = [r[0] for r in con.execute(
        "SELECT nucleotide_id FROM replicon ORDER BY nucleotide_id")]

    rows, tally = [], Counter()
    for nuc in nucs:
        handle, kind = open_record(args.gbk, nuc)
        if handle is None:
            tally["missing_file"] += 1
            rows.append((nuc, 0, "", kind))
            continue
        with handle:
            has_origin, declared = scan(handle)
        tally["with_sequence" if has_origin else "annotation_only"] += 1
        tally[f"{kind}_{'seq' if has_origin else 'noseq'}"] += 1
        rows.append((nuc, int(has_origin), declared or "", kind))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["nucleotide_id", "has_sequence", "declared_bp",
                         "file_kind"])
        writer.writerows(rows)

    #  PROMOTOR ANALIZININ PAYDASI. Zenginlesme olcumu yalnizca dizisi
    #  olan replikonlarda yapilabilir; kac DOGRULANMIS giris oraya
    #  dusuyor, simdiden yazili olmali.
    have = {nuc for nuc, flag, _bp, _k in rows if flag}
    placeholders = ",".join("?" * len(have)) if have else "''"
    confirmed_on_seq = con.execute(
        f"SELECT COUNT(*) FROM ro WHERE is_confirmed = 1 "
        f"AND nucleotide_id IN ({placeholders})",
        tuple(have)).fetchone()[0] if have else 0
    confirmed_total = con.execute(
        "SELECT COUNT(*) FROM ro WHERE is_confirmed = 1").fetchone()[0]

    print(f"[replikon] {len(nucs)} kayit")
    print(f"[dizi]     {tally['with_sequence']} tanesinde ORIGIN, "
          f"{tally['annotation_only']} yalnizca anotasyon, "
          f"{tally['missing_file']} dosya yok")
    print(f"[bicim]    {dict(tally)}")
    print(f"[giris]    dogrulanmis {confirmed_on_seq}/{confirmed_total} giris "
          f"dizisi okunabilir bir replikonda")
    print(f"[yazildi]  {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
