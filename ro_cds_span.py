"""
RO adayinin GENOMIK AYAK IZI -- min/max ile gercek parcalar ayri seylerdir.

SORUN (2026-11-04). `extract_genomic_context.py` icindeki `feature_spans`
dairesel replikonda orijini asan genler icin DOGRU davranir: parcalari
ayri ayri dondurur, ve komsu mesafeleri o parcalar uzerinden hesaplanir.
Ama adayin KENDI `start`/`end` degerleri parcalarin MIN ve MAX'i olarak
yazilir. Kaynaktaki yorum bunu "raporlama icin" diye isaretler; asagi
akista ise `build_operons.py`, `operon_validation.py` ve `analyze.py` bosluk
aritmetigini TAM O iki degerle yapar. Yani ayni tuzak komsular icin
kapatilmis, adayin kendisi icin kapatilmamis -- ve `ro` tablosu komsu
tablosunun tasidigi `spans_origin` alanini tasimiyor, yani bilgi kayboluyor.

OLCULDU (72.315 aday):
  span, cevirinin uc katindan buyuk   4.631  %6,40
  span > 20 kb                        1.319  %1,82
    bunlarin 1.312'si Eukaryota       -> intron; operon mantigi zaten gecersiz
    7'si Bacteria                     -> gercek vaka
  span == replikon uzunlugu, dairesel      2  CP063373.1 (9,96 Mb / 318 aa),
                                            CP060714.1 (5,17 Mb / 905 aa)

Boyle bir RO'da replikondaki HER komsu "ortusuyor" sayilir, bosluk hep 0
cikar, operon yurumesi hic kirilmaz ve sonuc TUM KROMOZOMU kapsayan bir
operon olur. Canli veritabaninda onaylanmis hicbir giriste yok (en buyuk
span 3.513 bp), yani bugun yayimlanmis bir sayiyi bozmuyor; A5 kayitlari
girdiginde girecekti.

NE YAPIYOR. Aciklanamaz ayak izi olan her aday icin kaynak GenBank yeniden
okunur ve GERCEK parcalar `ro_cds_span` turetilmis tablosuna yazilir:
parca sayisi, parcalarin toplam uzunlugu ve kenarlar. Bosluk aritmetigi
bundan sonra bu tabloyu okur.

NEDEN SUTUN DEGIL TABLO. Yeni tip-basi alan `chemistry.csv`e sutun olarak
eklenmez; turetilmis tablo + tazelik kontrolu kurali gecerli (sema
genisletme kurali). Ayrica yeniden cikarim gerektirmez: aciklanamaz ayak
izi olan 470 replikon okunur, 29.799'un tamami degil.

Cikti: `ro_cds_span` tablosu + analysis_out/ro_cds_span.json
"""

import argparse
import collections
import gzip
import json
import os
import re
import sqlite3
import sys

# Ayak izi, cevirinin bu katindan buyukse ACIKLANAMAZ sayilir ve kaynaktan
# yeniden okunur. 1.5: bir CDS'in genomik izi cevirisinin 3*(aa+1) bp'sini
# ancak intron/parca ile asar; tek parcali bir CDS'te oran 1,0'a cok yakindir.
UNEXPLAINED_RATIO = 1.5

SCHEMA = """
CREATE TABLE IF NOT EXISTS ro_cds_span (
  candidate_id TEXT PRIMARY KEY,
  n_parts      INTEGER NOT NULL,
  cds_length   INTEGER NOT NULL,   -- parcalarin TOPLAM uzunlugu (gercek ayak izi)
  reported_span INTEGER NOT NULL,  -- ro.end - ro.start (min/max)
  left_edge    INTEGER NOT NULL,   -- sola bakan en yakin parca kenari
  right_edge   INTEGER NOT NULL,   -- saga bakan en yakin parca kenari
  wraps_origin INTEGER NOT NULL,   -- parcalar replikonun iki ucuna dayaniyor mu
  parts        TEXT NOT NULL       -- "start-end;start-end" (denetlenebilir olsun)
)
"""


def translations(fasta_paths):
    """candidate_id -> ceviri uzunlugu (aa)."""
    out = {}
    for path in fasta_paths:
        if not os.path.exists(path):
            continue
        header, n = None, 0
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if header:
                        out[header] = n
                    header, n = line[1:].strip(), 0
                else:
                    n += len(line.strip())
            if header:
                out[header] = n
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--fasta", action="append", required=True,
                    help="ro_candidates.fasta; birden fazla verilebilir")
    ap.add_argument("--gbk-dir", action="append", required=True)
    ap.add_argument("--out", default="analysis_out/ro_cds_span.json")
    args = ap.parse_args()

    from Bio import SeqIO
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from extract_genomic_context import feature_spans

    aa = translations(args.fasta)
    con = sqlite3.connect(args.db)
    con.executescript(SCHEMA)

    rows = con.execute("SELECT candidate_id, nucleotide_id, start, end FROM ro").fetchall()
    # Aciklanamaz ayak izi olanlari sec: ceviri uzunlugundan beklenen bp'nin
    # UNEXPLAINED_RATIO katindan buyuk olanlar.
    todo = []
    for cid, nid, start, end in rows:
        n_aa = aa.get(cid)
        if n_aa is None:
            continue
        expected = 3 * (n_aa + 1)
        if (end - start) > UNEXPLAINED_RATIO * expected:
            todo.append((cid, nid, start, end, n_aa))
    print(f"[tarama] {len(rows):,} aday, aciklanamaz ayak izi {len(todo):,} "
          f"(%{100*len(todo)/max(1,len(rows)):.2f})")

    # Replikon basina gruplandir: her dosya BIR kez okunur.
    by_rep = collections.defaultdict(list)
    for item in todo:
        by_rep[item[1]].append(item)
    files = dict(con.execute("SELECT nucleotide_id, file FROM replicon"))
    lengths = dict(con.execute("SELECT nucleotide_id, length FROM replicon"))
    print(f"[kaynak] {len(by_rep):,} replikon okunacak")

    written = missing = unmatched = 0
    stats = collections.Counter()
    for nid, items in by_rep.items():
        fname = files.get(nid)
        path = None
        for d in args.gbk_dir:
            if fname and os.path.exists(os.path.join(d, fname)):
                path = os.path.join(d, fname)
                break
        if not path:
            missing += len(items)
            continue
        try:
            handle = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
            records = list(SeqIO.parse(handle, "genbank"))
        except Exception as exc:
            print(f"[hata] {nid}: {exc}", file=sys.stderr)
            missing += len(items)
            continue

        # (start, end, strand) -> gercek parcalar
        index = {}
        for rec in records:
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                spans = feature_spans(feat)
                key = (min(s for s, _ in spans), max(e for _, e in spans),
                       feat.location.strand or 0)
                index[key] = spans

        rlen = lengths.get(nid) or 0
        for cid, _nid, start, end, n_aa in items:
            m = re.match(r".*:(-?\d+)$", cid)
            strand = int(m.group(1)) if m else 0
            spans = index.get((start, end, strand))
            if spans is None:
                unmatched += 1
                continue
            parts = sorted(spans)
            cds_len = sum(e - s for s, e in parts)
            # Orijini asma: parcalar replikonun iki ucuna dayaniyor mu.
            wraps = bool(rlen and parts[0][0] <= 1 and parts[-1][1] >= rlen - 1
                         and len(parts) > 1)
            # Sola/saga bakan en yakin kenar: orijini asan bir gende min/max
            # DEGIL, gene komsu olan parcanin kenari dogrudur.
            left_edge = parts[-1][0] if wraps else parts[0][0]
            right_edge = parts[0][1] if wraps else parts[-1][1]
            con.execute(
                "INSERT OR REPLACE INTO ro_cds_span VALUES (?,?,?,?,?,?,?,?)",
                (cid, len(parts), cds_len, end - start, left_edge, right_edge,
                 int(wraps), ";".join(f"{s}-{e}" for s, e in parts)))
            written += 1
            stats["wraps_origin" if wraps else "multi_part"] += 1
    con.commit()

    out = {
        "n_candidates": len(rows),
        "unexplained_ratio": UNEXPLAINED_RATIO,
        "n_unexplained": len(todo),
        "n_resolved": written,
        "n_source_missing": missing,
        "n_feature_unmatched": unmatched,
        "by_kind": dict(stats),
    }
    print(f"[yazildi] ro_cds_span: {written:,} satir "
          f"({stats.get('wraps_origin', 0)} orijin asan, "
          f"{stats.get('multi_part', 0)} cok parcali)")
    print(f"          kaynak yok {missing:,}, feature eslesmedi {unmatched:,}")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[yazildi] {args.out}")
    con.close()


if __name__ == "__main__":
    main()
