# -*- coding: utf-8 -*-
"""RefSeq'in TUM protein uzayini ROmotif71 ile tara.

NEDEN: madde 95'in yazdigi gibi, bu DB'nin protein evreni UniProt'tan geldi,
NCBI'DAN DEGIL. Kullanicinin istedigi "NCBI'daki her RO olabilecek sey"
bugune kadar HIC aranmamisti. Bu betik o boslugu kapatir.

ESIK KEYFI DEGIL, OLCULDU: onaylanmis 19.125 girisin TAMAMI ROmotif71'e
E <= 5,6e-15 ile eslesiyor, en kotu bit skoru 55,5. Bu yuzden on suzgec
E <= 1e-3'te tutuluyor -- en kotu gercek girisin on kat mertebe altinda,
yani hicbir gercek RO on suzgecte dusmez. Suzgec KARAR VERMEZ; karari
ingest asamasindaki kapi (Rieske 4/4 + triad >=2/3) verir. Burada
toplanan sey ADAY'dir, uye degil.

AKITARAK CALISIR: 871 dosya, ~36 GB sikistirilmis, ~143 milyon protein.
Hicbiri diskte birikmez -- indir, tara, eslesenlerin dizisini al, sil.
Durum dosyasi tutuldugu icin kesilirse kaldigi yerden devam eder.
"""
import argparse
import gzip
import json
import os
import subprocess
import sys
import urllib.request

BASE = "https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/"
HMM = "ROs_71_Clean/ROmotif71.hmm"
WORK = "/media/lin-bio/back2/roar_ncbi"
EVALUE = "1e-3"


def listing():
    """Surumdeki protein parcalarinin adlari."""
    with urllib.request.urlopen(BASE, timeout=120) as fh:
        html = fh.read().decode("utf-8", "replace")
    import re
    names = sorted(set(re.findall(
        r"complete\.\d+\.protein\.faa\.gz", html)),
        key=lambda s: int(s.split(".")[1]))
    return names


def release():
    with urllib.request.urlopen(
            "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER",
            timeout=60) as fh:
        return fh.read().decode().strip()


def fetch(name, dest):
    urllib.request.urlretrieve(BASE + name, dest)


def scan_one(name, work, hmm, slot):
    """Bir parcayi indir, tara, eslesen dizileri dondur, dosyayi sil."""
    gz = os.path.join(work, "tmp", "slot%d.faa.gz" % slot)
    faa = os.path.join(work, "tmp", "slot%d.faa" % slot)
    tbl = os.path.join(work, "tmp", "slot%d.tbl" % slot)
    for path in (gz, faa, tbl):
        if os.path.exists(path):
            os.remove(path)
    fetch(name, gz)

    n_seq = 0
    with gzip.open(gz, "rt", errors="replace") as src, \
            open(faa, "w") as dst:
        for line in src:
            if line.startswith(">"):
                n_seq += 1
            dst.write(line)
    os.remove(gz)

    subprocess.run(
        ["hmmsearch", "--cpu", "2", "-E", EVALUE, "--noali",
         "--tblout", tbl, hmm, faa],
        check=True, stdout=subprocess.DEVNULL)

    hits = {}
    with open(tbl) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            # ad, tam-dizi E, tam-dizi bit skoru
            hits[f[0]] = (float(f[4]), float(f[5]))

    # Eslesenlerin dizisini ve tanimini ayni dosyadan cek.
    records = []
    if hits:
        keep, desc, chunks = None, None, []
        with open(faa, errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if keep:
                        records.append((keep, desc, "".join(chunks)))
                    head = line[1:].rstrip("\n")
                    acc = head.split(None, 1)[0]
                    keep = acc if acc in hits else None
                    desc = head.split(None, 1)[1] if " " in head else ""
                    chunks = []
                elif keep:
                    chunks.append(line.strip())
            if keep:
                records.append((keep, desc, "".join(chunks)))
    os.remove(faa)
    os.remove(tbl)
    return n_seq, hits, records


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=WORK)
    ap.add_argument("--hmm", default=HMM)
    ap.add_argument("--limit", type=int, default=0,
                    help="yalnizca ilk N parca (deneme icin)")
    ap.add_argument("--slot", type=int, default=0,
                    help="paralel kosuda bu islemin numarasi")
    ap.add_argument("--of", type=int, default=1,
                    help="toplam paralel islem sayisi")
    ap.add_argument("--parallel", type=int, default=0,
                    help="N alt islem baslat ve hepsini bekle; yeniden "
                         "uretilebilirlik icin zincirin cagirdigi bicim")
    args = ap.parse_args()

    # Darbogaz INDIRME, tarama degil: tek parca 5 saniyede iniyor ve 0,8
    # saniyede taraniyor. O yuzden paralellik slot sayisiyla saglaniyor,
    # hmmsearch'un --cpu'suyla degil.
    if args.parallel > 1:
        os.makedirs(os.path.join(args.work, "tmp"), exist_ok=True)
        kids = []
        for slot in range(args.parallel):
            cmd = [sys.executable, os.path.abspath(__file__),
                   "--work", args.work, "--hmm", args.hmm,
                   "--slot", str(slot), "--of", str(args.parallel)]
            if args.limit:
                cmd += ["--limit", str(args.limit)]
            log = open(os.path.join(args.work, "log_slot%d.txt" % slot), "a")
            kids.append(subprocess.Popen(cmd, stdout=log, stderr=log))
        bad = [k.wait() for k in kids]
        if any(bad):
            print("[hata] slot cikis kodlari: %s" % bad, file=sys.stderr)
            return 1
        tot = hit = 0
        for slot in range(args.parallel):
            st = os.path.join(args.work, "state_slot%d.json" % slot)
            if os.path.exists(st):
                d = json.load(open(st))
                tot += d["n_seq"]
                hit += d["n_hit"]
        print("[bitti] %s protein tarandi, %s aday"
              % (format(tot, ","), format(hit, ",")))
        return 0

    os.makedirs(os.path.join(args.work, "tmp"), exist_ok=True)
    names = listing()
    if args.limit:
        names = names[:args.limit]
    mine = [n for i, n in enumerate(names) if i % args.of == args.slot]
    rel = release()
    print("[surum] RefSeq release %s | %d parca, bu islem %d tane"
          % (rel, len(names), len(mine)), flush=True)

    state_path = os.path.join(args.work, "state_slot%d.json" % args.slot)
    state = {"release": rel, "done": [], "n_seq": 0, "n_hit": 0}
    if os.path.exists(state_path):
        state = json.load(open(state_path))
        if state.get("release") != rel:
            print("[uyari] durum dosyasi surum %s icin, simdiki %s -- sifirla"
                  % (state.get("release"), rel), flush=True)
            state = {"release": rel, "done": [], "n_seq": 0, "n_hit": 0}
    done = set(state["done"])

    fa_path = os.path.join(args.work, "hits_slot%d.faa" % args.slot)
    tb_path = os.path.join(args.work, "hits_slot%d.tsv" % args.slot)
    fa = open(fa_path, "a")
    tb = open(tb_path, "a")

    for i, name in enumerate(mine, 1):
        if name in done:
            continue
        try:
            n_seq, hits, records = scan_one(name, args.work, args.hmm,
                                            args.slot)
        except Exception as exc:                      # agdan kopma olagan
            print("[hata] %s: %s -- atlanmadi, tekrar denenecek"
                  % (name, exc), flush=True)
            continue
        for acc, desc, seq in records:
            ev, sc = hits[acc]
            fa.write(">%s %s\n%s\n" % (acc, desc, seq))
            tb.write("%s\t%s\t%.3g\t%.1f\t%d\t%s\n"
                     % (acc, name, ev, sc, len(seq), desc.replace("\t", " ")))
        fa.flush()
        tb.flush()
        state["n_seq"] += n_seq
        state["n_hit"] += len(hits)
        state["done"].append(name)
        json.dump(state, open(state_path, "w"))
        if i % 10 == 0 or i == len(mine):
            print("  %d/%d | %s | taranan %s dizi | aday %s"
                  % (i, len(mine), name, format(state["n_seq"], ","),
                     format(state["n_hit"], ",")), flush=True)

    fa.close()
    tb.close()
    print("[bitti] slot %d: %d dizi tarandi, %d aday"
          % (args.slot, state["n_seq"], state["n_hit"]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
