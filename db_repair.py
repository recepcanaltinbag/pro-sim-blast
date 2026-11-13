"""Veri duzeltmelerini BILDIRIMLE uygular, elle SQL ile degil.

NEDEN. `curate.py` ayni isi `chemistry.csv` icin yapiyor: her duzeltme eski
degeri de yazar, uygulamadan once DOGRULAR, ve iki kez kosturulabilir. Bu
dosya ayni sozlesmeyi VERITABANI satirlari icin kurar. Elle atilan bir
`DELETE` uc seyi kaybeder: neyin silindigi, NEDEN silindigi, ve o satirin
geri gelip gelmedigini anlamanin yolu.

KURAL. Bir duzeltme ancak SEBEBI KODDA kapatildiktan sonra buraya girer.
Yoksa bu dosya, pipeline'in her kosusunda ayni satiri yeniden silen bir
temizlik adimina doner -- yani hatayi gizler. Her duzeltme, hatayi ureten
yerin de kapandigini soyleyen bir `fixed_in` alani tasir.

GRAFIKTE OLMAK ZORUNDA -- ve 2026-11-12'ye kadar DEGILDI. Sebep kodda
kapandi (`extract_genomic_context.py` artik satiri uretmiyor), ama ARA
ARTEFAKT eski kodla uretilmis: `genomic_context/replicons.csv` o hayalet
satiri hala tasiyor (13087. satir, nucleotide_id bos). `build_db.py` o
CSV'yi okuyor, yani bir sonraki tam yeniden kurulum hayalet satiri GERI
GETIRIR ve onu silecek hicbir sey yok. Bir duzeltmeyi elle uygulamak,
duzeltmeyi bir kerelik yapar; grafige koymak, duzeltmeyi bir KURAL yapar.

Kullanim:
    python3 db_repair.py              # yalnizca rapor
    python3 db_repair.py --apply      # uygula
"""

import argparse
import json
import os
import sqlite3
import sys
import time


#  Her duzeltme: (ad, neden, fixed_in, on_kosul SQL, beklenen sayi, uygula SQL)
REPAIRS = [
    (
        "phantom_replicon",
        "A replicon row with no identifier. `gbk_files/NZ_QXDB01000002.1.gbk` "
        "is 0 bytes (a download that stopped half way in December 2024), so "
        "no record could be parsed and the row that was written carried "
        "nucleotide_id=''. Something with no identifier cannot be the UNIT "
        "of the replicon table, and this row made every per-replicon "
        "denominator one too large (17,074 instead of 17,073). No "
        "ro/neighbour/operon/taxon row depends on it, so deleting it loses "
        "no observation. The FACT that a requested record never arrived is "
        "not deleted: it is now declared in "
        "`genomic_context/unreadable_files.csv`.",
        "extract_genomic_context.py: a file whose identifier cannot be "
        "parsed is now given `status='unreadable'` and NO replicon row is "
        "written",
        "SELECT COUNT(*) FROM replicon "
        "WHERE nucleotide_id IS NULL OR TRIM(nucleotide_id) = ''",
        1,
        "DELETE FROM replicon "
        "WHERE nucleotide_id IS NULL OR TRIM(nucleotide_id) = ''",
    ),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--out", default="analysis_out/db_repair.json",
                    help="rapor dosyasi; adimin grafikteki ciktisi")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    pending = done = 0
    log = []
    for name, why, fixed_in, precondition, expected, statement in REPAIRS:
        found = con.execute(precondition).fetchone()[0]
        if found == 0:
            #  Idempotent: zaten uygulanmis. Bu bir hata degil.
            print(f"[uygulanmis] {name}: on kosul 0 satir buluyor")
            log.append({"repair": name, "state": "already_applied",
                        "rows_found": 0, "why": why, "fixed_in": fixed_in})
            done += 1
            continue
        if found != expected:
            #  SESSIZCE UYGULAMAZ. Beklenen sayi tutmuyorsa durum
            #  tanimlandigi gibi degil, ve bir duzeltmenin bilmedigi bir
            #  durumu "duzeltmesi" en kotu sonuctur.
            print(f"[DURDU] {name}: {expected} satir bekleniyordu, "
                  f"{found} bulundu. Hicbir sey yazilmadi.")
            sys.exit(1)
        print(f"[bekliyor] {name}: {found} satir")
        print(f"    neden    : {why}")
        print(f"    kapatildi: {fixed_in}")
        if args.apply:
            con.execute(statement)
            con.commit()
            kalan = con.execute(precondition).fetchone()[0]
            if kalan:
                print(f"[DURDU] {name}: uygulamadan sonra {kalan} satir kaldi")
                sys.exit(1)
            print(f"[uygulandi] {name}")
        log.append({"repair": name,
                    "state": "applied" if args.apply else "pending",
                    "rows_found": found, "why": why, "fixed_in": fixed_in})
        pending += 1

    if pending and not args.apply:
        print(f"\n{pending} duzeltme bekliyor. Uygulamak icin: "
              f"python3 db_repair.py --apply")
    elif not pending:
        print(f"\n{done}/{len(REPAIRS)} duzeltme zaten uygulanmis.")

    report = {
        "question": "Which data repairs does this database carry, and why?",
        "rule": ("a repair enters this file only after its cause is closed in "
                 "code; every entry names the place that was fixed. Applying a "
                 "repair by hand fixes it once; putting it in the rebuild graph "
                 "makes it a rule, which matters because the intermediate CSV "
                 "still carries the row the new code no longer produces."),
        "checked_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "repairs": log,
        "pending": pending if not args.apply else 0,
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    con.close()


if __name__ == "__main__":
    main()
