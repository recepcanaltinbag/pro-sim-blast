"""
Ara dosya atlama kurali -- TEK kaynak.

ATLAMA KURALI: bir ara dosyanin VAR OLMASI onu guncel yapmaz. Girdilerinden
herhangi biri daha yeniyse yeniden uretilmek ZORUNDADIR.

NEDEN BU DOSYA VAR. Ayni tuzak bu depoda IKI KEZ kurban aldi:

1. `annotate_ro.py` (2026-10): iki ara dosya 21 Temmuz'dan kalmisti, aday
   fastasi 9 Ekim'de yeniden uretilmisti. Script "zaten var" gorup
   hmmsearch'u atladi, 26.693 yeni adayi 14.229 adaylik eski taramanin
   sonucuyla eslestirmeye calisti ve DOGRULANMIS RO: 0 yazdi. Hicbir sey
   hata vermedi. Orada `stale()` yazildi ve sorun cozuldu.

2. `build_operons.py` (2026-11-04'te bulundu): ayni satir duzeltilmemisti --
   `if force or not os.path.exists(domtbl) or os.path.getsize(domtbl) == 0`.
   `operon_work/neighbor_components.domtbl` 2026-10-05'ten,
   `operon_work/neighbor_proteins.fasta` ise 2026-10-14'ten kalmisti. Yani
   komsu proteinleri yenilenmis, bilesen taramasi yenilenmemisti. Sonuc:
   2026-10-05'ten SONRA eklenen her komsu protein bilesen taramasi GORMEDI
   ve hepsi `none` olarak siniflandi -- beta, ferredoksin ve reduktaz
   dahil. Kanit: anotasyonun "dioxygenase small subunit" dedigi 628 protein
   domtbl'de HIC gecmiyor, ama Pfam'in kendi Ring_hydroxyl_B modeli o
   628'in 628'ini E <= 1e-5 ile tutuyor.

Yani duzeltme BIR dosyaya uygulanmis, AYNI kurali uygulayan otekine
uygulanmamisti. Bu, bu projenin tekrar eden hata bicimidir: karar veren
kural, iddia edilen kural degil. O yuzden kural artik kopyalanmaz, BURADAN
alinir.

NOT. `rebuild.py:staleness` ve `provenance.py:stale_inputs` AYRI isler
yapar (yeniden kurma grafigi ve koken raporu) ve bilincli olarak
birlestirilmedi; burasi yalnizca "ara dosya girdisinden eski mi" sorusudur.
"""

import os


def stale(output, *inputs):
    """Ara dosya yeniden uretilmeli mi?

    Dondurur: None (guncel) ya da NEDENI aciklayan bir dizgi. Neden
    dondurulur cunku atlama karari LOGA yazilabilir olmali -- sessiz bir
    atlama tam olarak yukaridaki iki olayin olusma bicimiydi.
    """
    if not os.path.exists(output) or os.path.getsize(output) == 0:
        return "yok"
    out_age = os.path.getmtime(output)
    newer = [p for p in inputs
             if os.path.exists(p) and os.path.getmtime(p) > out_age]
    if not newer:
        return None
    return "girdi daha yeni: " + ", ".join(os.path.basename(p) for p in newer)
