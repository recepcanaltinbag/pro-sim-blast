# -*- coding: utf-8 -*-
"""RefSeq taramasinin adaylarini birlesime al.

scan_ncbi_proteins.py ham adaylari toplar; bu betik onlari KARARA baglar ve
`ro_union`a yazar.

AYNI KURAL, YENI KURAL DEGIL. Tip atamasi `ro_filter.classify_target` ile,
yani UniProt tarafinda kullanilan kuralin AYNISI ile yapilir (RieskeDB71'in
76 referans modeli, kapsama barajini gecenler arasindan en yuksek BIT
SKORU -- en kucuk E-degeri degil). Triad kapisi da ayni:
ingest_protein_entries.py'nin `site_residue`'u ve `ro_motif`in kolonlari
dogrudan import edilir. Bu dosyada kopyalanmis tek bir esik yoktur; olsaydi
iki kaynak zamanla ayrisir ve "NCBI'dan gelenler baska bir kurala gore
sayilmis" olurdu.

NE YAZILIR, NE YAZILMAZ
  `ro`ya DOKUNULMAZ. O tablonun her satiri bir GENOM KONUMUDUR ve sitedeki
  konuma bagli her istatistik bu varsayima dayanir (madde: ingest_protein_
  entries.py'nin basindaki gerekce). RefSeq protein kayitlarinin konumu
  yoktur.

  BU BETIK `ro_union`A DA YAZMAZ. Kendi KAYNAK tablosunu yazar
  (`ro_ncbi_hit`), birlesimi ise tek bir yer kurar: build_union.py. Ilk
  tasarimda buradan dogrudan `ro_union`a yaziliyordu ve bu SESSIZ BIR VERI
  KAYBI olurdu -- build_union.py tabloyu DROP edip yeniden kuruyor, yani
  bir sonraki kurulum bu satirlari silerdi ve kimse fark etmezdi. Kaynak
  tablosu + tek birlestirici, `ro_protein` / `ro_protein_status`'ta
  kullanilan ayni dizen.
"""
import argparse
import glob
import hashlib
import os
import subprocess
import sqlite3
import sys

import ro_motif
from ro_motif import (BRIDGING_SITE, CATALYTIC_SITES, RIESKE_SITES,
                      read_stockholm_matchcols)
from ingest_protein_entries import site_residue
from ro_filter import (parse_domtblout, classify_target,
                       DEFAULT_MIN_COVERAGE, DEFAULT_MAX_EVALUE,
                       DEFAULT_MIN_LENGTH)

SCAN_DIR = "/media/lin-bio/back2/roar_ncbi"


def md5(seq):
    return hashlib.md5(seq.strip().upper().encode()).hexdigest()


def gather(scan_dir, out_fasta):
    """slot dosyalarini tek fasta'da birlestir, accession'a gore tekille."""
    seen, n_dup = set(), 0
    desc = {}
    with open(out_fasta, "w") as dst:
        for path in sorted(glob.glob(os.path.join(scan_dir, "hits_slot*.faa"))):
            keep, chunks = None, []
            with open(path, errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if keep:
                            dst.write(">%s\n%s\n" % (keep, "".join(chunks)))
                        head = line[1:].rstrip("\n")
                        acc = head.split(None, 1)[0]
                        if acc in seen:
                            keep, n_dup = None, n_dup + 1
                        else:
                            seen.add(acc)
                            keep = acc
                            desc[acc] = (head.split(None, 1)[1]
                                         if " " in head else "")
                        chunks = []
                    elif keep:
                        chunks.append(line.strip())
                if keep:
                    dst.write(">%s\n%s\n" % (keep, "".join(chunks)))
    print("[toplandi] %s tekil aday (%s tekrar atildi)"
          % (format(len(seen), ","), format(n_dup, ",")), flush=True)
    return desc


def read_fasta(path):
    out, name, chunks = {}, None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name:
                chunks.append(line.strip())
    if name:
        out[name] = "".join(chunks)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--scan-dir", default=SCAN_DIR)
    ap.add_argument("--work", default=os.path.join(SCAN_DIR, "ingest"))
    ap.add_argument("--library-hmm", default="RieskeDB71.hmm",
                    help="76 referans modeli; tip atamasi bununla yapilir")
    ap.add_argument("--motif-hmm", default="ROs_71_Clean/ROmotif71.hmm",
                    help="603 kolonluk profil; triad kapisi bununla okunur")
    ap.add_argument("--cpu", type=int, default=14)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    hits_fa = os.path.join(args.work, "ncbi_hits.faa")
    desc = gather(args.scan_dir, hits_fa)
    if not desc:
        print("[hata] hic aday yok -- tarama bitmis mi?", file=sys.stderr)
        return 2
    seqs = read_fasta(hits_fa)

    # ---- 1. TIP ATAMASI: UniProt tarafiyla AYNI kural ----
    domtbl = os.path.join(args.work, "ncbi_domtbl.out")
    if not (os.path.exists(domtbl) and os.path.getsize(domtbl) > 0
            and os.path.getmtime(domtbl) >= os.path.getmtime(hits_fa)):
        print("[calisiyor] hmmsearch: 76 referans modeli", flush=True)
        subprocess.run(
            ["hmmsearch", "--cpu", str(args.cpu), "-E", "1e-5",
             "--domtblout", domtbl, "--noali",
             args.library_hmm, hits_fa],
            check=True, stdout=subprocess.DEVNULL)
    else:
        print("[atlandi] %s guncel" % domtbl)

    by_status = {}
    assigned = {}
    for target, target_hits in parse_domtblout(domtbl).items():
        result, status = classify_target(
            target_hits, DEFAULT_MIN_COVERAGE, DEFAULT_MAX_EVALUE,
            DEFAULT_MIN_LENGTH)
        by_status[status] = by_status.get(status, 0) + 1
        if status == "RO_alpha":
            assigned[target] = result
    print("[atama] " + " | ".join("%s %s" % (k, format(v, ","))
                                  for k, v in sorted(by_status.items())))
    print("[atama] RO_alpha: %s" % format(len(assigned), ","))
    if not assigned:
        print("[hata] hicbiri RO_alpha degil -- durduruluyor", file=sys.stderr)
        return 2

    # ---- 2. TRIAD KAPISI: ingest_protein_entries ile AYNI okuma ----
    alpha_fa = os.path.join(args.work, "ncbi_alpha.faa")
    with open(alpha_fa, "w") as fh:
        for acc in sorted(assigned):
            fh.write(">%s\n%s\n" % (acc, seqs[acc]))
    sto = os.path.join(args.work, "ncbi_aln.sto")
    if not (os.path.exists(sto) and os.path.getsize(sto) > 0
            and os.path.getmtime(sto) >= os.path.getmtime(alpha_fa)):
        print("[calisiyor] hmmalign --trim", flush=True)
        subprocess.run(["hmmalign", "--trim", "--amino", "-o", sto,
                        args.motif_hmm, alpha_fa], check=True)
    else:
        print("[atlandi] %s guncel" % sto)
    aligned = read_stockholm_matchcols(sto)
    print("[hizalandi] %s" % format(len(aligned), ","))

    measured = {}
    gate_pass = 0
    for acc, al in aligned.items():
        rieske_ok = all(site_residue(al, c, r)[1] is not None
                        for c, r, _ in RIESKE_SITES)
        triad = [site_residue(al, c, r) for c, r, _ in CATALYTIC_SITES]
        n_triad = sum(1 for _, off in triad if off is not None)
        passes = rieske_ok and n_triad >= ro_motif.MIN_CATALYTIC_SITES
        gate_pass += int(passes)
        measured[acc] = (int(rieske_ok), n_triad, int(passes))
    print("[kapi] Rieske 4/4 + triad >=%d gecen: %s / %s"
          % (ro_motif.MIN_CATALYTIC_SITES, format(gate_pass, ","),
             format(len(measured), ",")))

    # ---- 3. KAYNAK TABLOSUNU YAZ (birlesimi build_union.py kurar) ----
    con = sqlite3.connect(args.db)
    rows = []
    for acc, res in assigned.items():
        seq = seqs.get(acc, "")
        if not seq or acc not in measured:
            continue
        rieske_ok, n_triad, passes = measured[acc]
        rows.append((
            acc, md5(seq), seq.strip().upper(), len(seq.strip()),
            res.get("PredictedCluster"), res.get("Group"),
            res.get("ModelCoverage"), res.get("Score"), res.get("E-value"),
            rieske_ok, n_triad, passes,
            (desc.get(acc) or "").strip() or None))

    print("\n[kaynak] RO_alpha ve olculmus: %s" % format(len(rows), ","))
    uniq = len({r[1] for r in rows})
    print("   tekil dizi: %s" % format(uniq, ","))
    if args.dry_run:
        print("[dry-run] tablo yazilmadi")
        return 0

    con.executescript("""
        DROP TABLE IF EXISTS ro_ncbi_hit;
        CREATE TABLE ro_ncbi_hit (
            accession         TEXT PRIMARY KEY,
            seq_md5           TEXT,
            sequence          TEXT,
            protein_length    INTEGER,
            ro_cluster        TEXT,
            ro_group          TEXT,
            model_coverage    REAL,
            hmm_score         REAL,
            hmm_evalue        REAL,
            rieske_intact     INTEGER,
            triad_sites_found INTEGER,
            passes_current_gate INTEGER,
            product           TEXT
        );
        CREATE INDEX idx_rnh_md5     ON ro_ncbi_hit(seq_md5);
        CREATE INDEX idx_rnh_cluster ON ro_ncbi_hit(ro_cluster);
    """)
    con.executemany(
        "INSERT INTO ro_ncbi_hit (accession, seq_md5, sequence, "
        "protein_length, ro_cluster, ro_group, model_coverage, hmm_score, "
        "hmm_evalue, rieske_intact, triad_sites_found, passes_current_gate, "
        "product) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit()
    n = con.execute("SELECT COUNT(*) FROM ro_ncbi_hit").fetchone()[0]
    print("[yazildi] ro_ncbi_hit %s satir" % format(n, ","))
    print("[sonraki] python3 build_union.py  -- birlesimi yeniden kurar")
    return 0


if __name__ == "__main__":
    sys.exit(main())
