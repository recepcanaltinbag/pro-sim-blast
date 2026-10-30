#!/usr/bin/env python3
"""Triad kapisini SABIT NOKTAYA getir: ro_triad artik degismeyene kadar don.

SORUN. `triad_gate.py` ile `triad_anchor_his.py` birbirini besliyor:

    triad_anchor_his.py  ro_triad'i OKUR  -> hangi girislerin His'i sabit
                         kolonda okunamadigini oradan ogrenir
    triad_gate.py        ro_triad'i YAZAR -> capalamanin JSON ciktisini okuyup
                         3/3 kararini tabloya basar

Yani tek gecis YETMEZ. 2026-10-30'da olculdu: yeniden kurulum onaylanmis kumeyi
19.114'ten 19.125'e cikardi; capalama BAYAT ro_triad'dan is listesi okudugu
icin yeni 12 girisi hic gormedi, ve birinci gecisten sonra ro_triad'da 12 giris
`his_needs_anchoring` kaldi. Ikinci gecis 0'a indirdi.

Bu betik o turu kodluyor, elle dondurmek yerine:

    gate --apply  ->  his_needs_anchoring kac?  ->  0 ise BITTI
                                                ->  degilse anchor, basa don

Yakinsama, secimin `his_source` uzerinden yapilmasina DAYANIR (triad_status
uzerinden DEGIL): `his_source != 'fixed_column'` kosudan kosuya sabittir, oysa
triad_status'u gate'in kendisi degistirir ve duruma bakan bir secim ikinci
turda BOS kume dondurur. Ayrintisi triad_anchor_his.py'de.

NEDEN BU AYRI BIR ADIM. rebuild.py bir adimin yazdigi tablolari, kendisinden
ONCEKI adimlarin tazelik karsilastirmasinda yok sayar (`ignore`). ro_triad ve
ro_triad_meta'yi yazan betikler zincirde adim OLMADIGI icin hicbir adimin
`db_writes` kumesinde degillerdi, dolayisiyla hicbirinin `ignore`una girmiyor,
ve `roar.sqlite`i girdi sayan HER adim onlari "degismis ust girdi" saniyordu.
2026-10-30'da olculdu: 51 artefaktin 45'i bu yuzden BAYAT gorunuyordu, oysa
52 adim betiginden yalnizca validate_curation.py bu tablolari okuyor. Sonuc:
`rebuild.py --check` yayin oncesi kapi olarak ise yaramaz hale geliyordu --
45 adim her zaman "bayat" diyorsa gercek bayatligi sahtesinden ayirt edemezsin.
Cozum tazelik kontrolunu zayiflatmak (ro_triad'i izlemekten cikarmak) DEGIL,
bu yazmalarin zincirdeki gercek yerini dogru yazmaktir. Bu betik o yer.
"""

import argparse
import os
import sqlite3
import subprocess
import sys

MAX_PASSES = 5


def log(msg):
    print(msg, flush=True)


def needs_anchoring(db):
    """ro_triad'da His okumasi hala capalanmamis giris sayisi."""
    con = sqlite3.connect("file:%s?mode=ro" % db, uri=True)
    try:
        have = con.execute(
            "SELECT name FROM sqlite_master WHERE name='ro_triad'").fetchone()
        if not have:
            return None
        return con.execute("SELECT COUNT(*) FROM ro_triad "
                           "WHERE triad_status = 'his_needs_anchoring'").fetchone()[0]
    finally:
        con.close()


def run(cmd):
    log("    $ " + " ".join(cmd))
    res = subprocess.run(cmd)
    if res.returncode != 0:
        log("[HATA] cikis %d: %s" % (res.returncode, " ".join(cmd)))
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--work", default=None,
                    help="triad_anchor_his.py'nin hizalama dizini")
    ap.add_argument("--mafft-threads", type=int, default=1)
    ap.add_argument("--max-passes", type=int, default=MAX_PASSES)
    args = ap.parse_args()

    gate = [sys.executable, "triad_gate.py", "--db", args.db, "--apply"]
    anchor = [sys.executable, "triad_anchor_his.py", "--db", args.db,
              "--mafft-threads", str(args.mafft_threads)]
    if args.work:
        anchor += ["--work", args.work]

    for i in range(1, args.max_passes + 1):
        log("\n=== gecis %d/%d: kapiyi uygula ===" % (i, args.max_passes))
        if not run(gate):
            return 2
        left = needs_anchoring(args.db)
        if left is None:
            log("[HATA] ro_triad yazilmadi.")
            return 2
        log("[durum] his_needs_anchoring: %d giris" % left)
        if left == 0:
            log("\n[yakinsadi] %d gecisde sabit noktaya ulasildi." % i)
            return 0
        log("=== gecis %d: His kolonlarini capala ===" % i)
        # --reuse-alignments KULLANILMAZ: is listesi her gecisde degisiyor,
        # yani onbellekli hizalamalar yeni girisleri icermez.
        if not run(anchor):
            return 2

    log("\n[HATA] %d gecisde yakinsamadi; his_needs_anchoring hala %d. "
        "Bu bir dongu isaretidir, kapi SIKILASTIRILMAMALI."
        % (args.max_passes, needs_anchoring(args.db)))
    return 1


if __name__ == "__main__":
    sys.exit(main())
