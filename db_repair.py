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

Kullanim:
    python3 db_repair.py              # yalnizca rapor
    python3 db_repair.py --apply      # uygula
"""

import argparse
import sqlite3
import sys


#  Her duzeltme: (ad, neden, fixed_in, on_kosul SQL, beklenen sayi, uygula SQL)
REPAIRS = [
    (
        "phantom_replicon",
        "Kimligi olmayan replikon satiri. `gbk_files/NZ_QXDB01000002.1.gbk` "
        "0 bayt (Ara 2024'te yarim kalmis indirme), yani hicbir kayit "
        "ayristirilamadi; uretilen satir nucleotide_id='' tasiyordu. "
        "Kimligi olmayan bir sey replikon tablosunun BIRIMI olamaz, ve bu "
        "satir her replikon-basi paydayi bir fazla gosteriyordu "
        "(17.074 yerine 17.073). Hicbir ro/neighbor/operon/taxon satiri "
        "ona bagli degil -- silmek hicbir gozlem kaybetmiyor. Istenen bir "
        "kaydin hic gelmedigi OLGUSU silinmiyor: artik "
        "`genomic_context/unreadable_files.csv` icinde beyan ediliyor.",
        "extract_genomic_context.py: kimlik ayristirilamayan dosya artik "
        "`status='unreadable'` alir ve replicon satiri URETILMEZ",
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
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    pending = done = 0
    for name, why, fixed_in, precondition, expected, statement in REPAIRS:
        found = con.execute(precondition).fetchone()[0]
        if found == 0:
            #  Idempotent: zaten uygulanmis. Bu bir hata degil.
            print(f"[uygulanmis] {name}: on kosul 0 satir buluyor")
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
        pending += 1

    if pending and not args.apply:
        print(f"\n{pending} duzeltme bekliyor. Uygulamak icin: "
              f"python3 db_repair.py --apply")
    elif not pending:
        print(f"\n{done}/{len(REPAIRS)} duzeltme zaten uygulanmis.")
    con.close()


if __name__ == "__main__":
    main()
