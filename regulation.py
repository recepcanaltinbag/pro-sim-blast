"""
Regulasyon baglami: her dogrulanmis RO icin en yakin transkripsiyon regulatoru
ve (divergent ciftlerde) aradaki promoter bolgesi.

Yon siniflari (RO ipligine gore, genom koordinatina gore DEGIL):
    divergent      regulator RO'nun 5' tarafinda, ters iplikte  (kafa kafaya; LysR klasigi)
    convergent     regulator RO'nun 3' tarafinda, ters iplikte  (kuyruk kuyruga)
    upstream_co    ayni iplik, 5' tarafta (ayni operonda olabilir)
    downstream_co  ayni iplik, 3' tarafta

Divergent ciftlerde regulator ile RO arasindaki DNA (<= --max-intergenic bp)
gbk'dan okunur, RO yonune cevrilir ve basit bir sigma70 taramasi yapilir:
-35 (TTGACA) ve -10 (TATAAT) kutulari, her birinde <= --mismatch uyumsuzluk,
15-19 bp aralik. Bu TAHMINDIR; deneysel promoter degil, inceleme icin ipucu.

Cikti: ro_regulation tablosu + analysis_out/regulation_by_cluster.csv
"""

import argparse
import csv
import gzip
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq

MAX_DIST = 1000          # regulator-RO arasi en fazla bp (herhangi yon)
BOX35, BOX10 = "TTGACA", "TATAAT"


def orientation(ro_strand, nb_strand, distance):
    """distance: isaretli genom koordinati (negatif = solda)."""
    upstream = (distance < 0) if ro_strand == 1 else (distance > 0)
    if nb_strand != ro_strand:
        return "divergent" if upstream else "convergent"
    return "upstream_co" if upstream else "downstream_co"


def mismatches(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def scan_sigma70(seq, max_mm):
    """Donen: [(pos35, spacer, pos10, mm35, mm10)] -- en iyi 3 aday."""
    seq = seq.upper()
    hits = []
    for i in range(len(seq) - 6):
        m35 = mismatches(seq[i:i + 6], BOX35)
        if m35 > max_mm:
            continue
        for spacer in range(15, 20):
            j = i + 6 + spacer
            if j + 6 > len(seq):
                break
            m10 = mismatches(seq[j:j + 6], BOX10)
            if m10 <= max_mm:
                hits.append((m35 + m10, i, spacer, j, m35, m10))
    hits.sort()
    return [h[1:] for h in hits[:3]]


def _extract(job):
    """Bir gbk dosyasindan istenen aralik dizilerini cek. job = (path, [(cid, start, end, strand)])"""
    path, wants = job
    out = []
    try:
        handle = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
        with handle:
            rec = next(SeqIO.parse(handle, "genbank"))
            for cid, s, e, strand in wants:
                s, e = max(0, s), min(len(rec.seq), e)
                if e <= s:
                    out.append((cid, ""))
                    continue
                frag = rec.seq[s:e]
                if strand == -1:
                    frag = frag.reverse_complement()
                out.append((cid, str(frag)))
    except Exception as exc:
        sys.stderr.write(f"[hata] {path}: {exc}\n")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db",