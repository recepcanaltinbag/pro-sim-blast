"""
PROTEIN DUZEYINDE GIRISLER -- genomik baglami olmayan RO alpha'lar.

NEDEN BU MODUL VAR
  Pipeline'in ilk adimi (run_ro_filter.py) UniProt'un PF00355 setinde
  73.699 RO alpha buluyor (68.747 tekil dizi). Yayinlanan veritabaninda ise
  19.114 giris (15.953 tekil dizi) var ve bunlarin 14.449'u o protein-duzeyi
  setle birebir ayni. Yani kendi ilk adimimizin RO alpha dedigi tekil
  dizilerin ~%79'u (54.298) veritabanina HIC GIRMIYOR.

  Sebep yontemseldi, biyolojik degil: `ro` tablosu GENOM kayitlarindan
  kuruluyor ve bir protein ancak annotasyonlu bir GenBank kaydina
  baglanabildiyse giriyor. Genom baglami operon, komsuluk, plazmit ve
  ekoloji analizleri icin gerekli -- ama bir proteinin RO alpha OLMASI icin
  gerekli degil. Kullanici karari (2026-10-30): bu proteinler de
  veritabaninda olmali.

NEDEN AYRI TABLO, `ro`'YA EK SATIR DEGIL
  `ro` tablosunun her satiri bir genom konumudur ve sitedeki her giris-bazli
  istatistik (plazmit orani, operon kurallari, komsuluk zenginlesmesi,
  izolasyon kaynagi, regulasyon mimarisi) bu varsayima dayanir. 54.298
  baglamsiz satir `ro`'ya eklenirse bu istatistiklerin hepsi SESSIZCE
  anlam degistirir: payda 3,4 katina cikar ve pay degismez. Ayri tablo
  bunu imkansiz kilar -- baglam isteyen bir sorgu bu satirlari bulamaz,
  yani yanlis sayi yerine EKSIK sayi uretir. Eksik sayi fark edilir,
  yanlis sayi fark edilmez.

NE OLCULUYOR
  Her protein ayni kapidan gecirilir: `ro_motif.classify_motifs` ile Rieske
  ligandlari 4/4 ve katalitik triad kolonlari. Triad'in UC uyesi de ayri
  ayri kaydedilir (yalnizca "kac tanesi" degil), cunku kapinin 2/3'ten 3/3'e
  cikarilmasi tartisiliyor ve o karar her iki tabloda ayni olcume dayanmak
  zorunda.

  Dizisi bir genom-baglamli girisle BIREBIR ayni olan proteinler ayrica
  isaretlenir (`linked_candidate_id`): bunlar yeni enzim degil, ayni enzimin
  baglamsiz kopyasidir ve hicbir sayimda iki kez sayilmamalidir.

Kullanim:
    python3 ingest_protein_entries.py --work /media/lin-bio/back2/roar_protein
"""

import argparse
import os
import sqlite3
import subprocess
import sys

import ro_motif
from ro_motif import (BRIDGING_SITE, CATALYTIC_SITES, RIESKE_SITES,
                      read_stockholm_matchcols)

WINDOW = ro_motif.COLUMN_WINDOW


def read_fasta_headers(path):
    """Donen: {dizi_adi: (cluster, coverage)} -- fasta basligindan."""
    out = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split()
            name = parts[0]
            meta = dict(p.split("=", 1) for p in parts[1:] if "=" in p)
            out[name] = (meta.get("cluster"), meta.get("cov"))
    return out


def site_residue(aligned, column, accepted):
    """Kolonda (+-WINDOW) kabul edilen kalinti var mi. Donen: (kalinti, kayma)."""
    length = len(aligned)
    for delta in sorted(range(-WINDOW, WINDOW + 1), key=lambda x: (abs(x), x)):
        idx = column - 1 + delta
        if 0 <= idx < length and aligned[idx].upper() in accepted:
            return aligned[idx].upper(), delta
    idx = column - 1
    return (aligned[idx].upper() if 0 <= idx < length else "?"), None


def run_hmmalign(hmm, fasta, out_sto):
    if os.path.exists(out_sto) and os.path.getsize(out_sto) > 0:
        if os.path.getmtime(out_sto) >= max(os.path.getmtime(hmm),
                                            os.path.getmtime(fasta)):
            print("[atlandi] %s guncel" % out_sto)
            return
    print("[calisiyor] hmmalign --trim -> %s" % out_sto)
    subprocess.run(["hmmalign", "--trim", "--amino", "-o", out_sto, hmm, fasta],
                   check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--fasta", default="ro_alpha.fasta")
    ap.add_argument("--hmm", default="ROs_71_Clean/ROmotif71.hmm")
    ap.add_argument("--work", default="/media/lin-bio/back2/roar_protein")
    ap.add_argument("--dry-run", action="store_true",
                    help="tabloyu yazma, yalnizca sayilari bas")
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    sto = os.path.join(args.work, "protein_aln.sto")
    run_hmmalign(args.hmm, args.fasta, sto)

    print("[okunuyor] %s" % sto)
    aligned = read_stockholm_matchcols(sto)
    print("[bilgi] hizalanan protein: %d" % len(aligned))
    meta = read_fasta_headers(args.fasta)

    con = sqlite3.connect(args.db)
    by_sequence = {}
    for cid, seq in con.execute(
            "SELECT candidate_id, sequence FROM ro "
            "WHERE is_confirmed=1 AND sequence IS NOT NULL"):
        by_sequence.setdefault(seq, cid)
    print("[bilgi] genom-baglamli onaylanmis tekil dizi: %d" % len(by_sequence))

    raw = {}
    name = None
    with open(args.fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                raw[name] = []
            elif name:
                raw[name].append(line.strip())
    raw = {k: "".join(v) for k, v in raw.items()}

    rows = []
    gate_pass = linked = triad3 = 0
    for pname, al in aligned.items():
        cluster, cov = meta.get(pname, (None, None))
        rieske = [site_residue(al, c, r)[1] is not None
                  for c, r, _ in RIESKE_SITES]
        triad = [site_residue(al, c, r) for c, r, _ in CATALYTIC_SITES]
        carbox_res, carbox_off = triad[2]
        bridge_res, bridge_off = site_residue(al, BRIDGING_SITE[0],
                                              BRIDGING_SITE[1])
        n_triad = sum(1 for _, off in triad if off is not None)
        rieske_ok = all(rieske)
        passes = rieske_ok and n_triad >= ro_motif.MIN_CATALYTIC_SITES
        gate_pass += int(passes)
        triad3 += int(rieske_ok and n_triad == 3)
        seq = raw.get(pname, "")
        link = by_sequence.get(seq)
        linked += int(link is not None)
        rows.append((
            pname,
            pname.split("|")[1] if "|" in pname else pname,
            cluster, float(cov) if cov else None, len(seq),
            int(rieske_ok), n_triad,
            int(triad[0][1] is not None), int(triad[1][1] is not None),
            carbox_res, carbox_off, bridge_res, bridge_off,
            int(passes), link,
        ))

    print("\n[sonuc] protein-duzeyi giris: %d" % len(rows))
    print("   Rieske 4/4 + triad >=2/3 (mevcut kapi) : %d" % gate_pass)
    # Bu sayi SABIT KOLON okumasidir ve dolayisiyla SISMIS bir yokluk
    # bildirir. `ro` tablosunda ayni okuma 2388 giriste "karboksilat yok"
    # dedi; tip-ici capali okuma o sayiyi 24'e dusurdu. Buradaki 3/3 sayisi
    # bir KAPI DEGERI olarak kullanilmamalidir -- protein-duzeyi girislerde
    # henuz capali okuma YOK (tip basina pozitif capa kurulmasi gerekiyor).
    print("   Rieske 4/4 + triad 3/3 SABIT KOLONDA     : %d" % triad3)
    print("      (SISMIS: sabit kolon okumasi; `ro`da ayni okuma 2388 yokluk "
          "dedi, capali okuma 24'e dusurdu -- kapi degeri DEGIL)")
    print("   genom-baglamli bir girisle AYNI dizi   : %d" % linked)
    print("   yeni, baglamsiz dizi                   : %d"
          % (len(rows) - linked))

    if args.dry_run:
        print("[dry-run] tablo yazilmadi")
        return 0

    con.executescript("""
        DROP TABLE IF EXISTS ro_protein;
        CREATE TABLE ro_protein (
            protein_name        TEXT PRIMARY KEY,
            accession           TEXT,
            ro_cluster          TEXT,
            model_coverage      REAL,
            protein_length      INTEGER,
            rieske_intact       INTEGER,
            triad_sites_found   INTEGER,
            catalytic_his1      INTEGER,
            catalytic_his2      INTEGER,
            carboxylate_residue TEXT,
            carboxylate_offset  INTEGER,
            bridging_residue    TEXT,
            bridging_offset     INTEGER,
            passes_current_gate INTEGER,
            linked_candidate_id TEXT
        );
        CREATE INDEX idx_rop_cluster ON ro_protein(ro_cluster);
        CREATE INDEX idx_rop_linked  ON ro_protein(linked_candidate_id);
    """)
    con.executemany("INSERT INTO ro_protein VALUES (%s)" % ",".join(["?"] * 15),
                    rows)
    con.commit()
    con.close()
    print("[yazildi] ro_protein: %d satir" % len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
