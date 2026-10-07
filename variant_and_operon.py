"""
Varyant karsilastirmasini CEBE tasiyan, ve operon/transpozon iliskisini
birbirinden ayirmaya calisan uc olcum.

NEDEN BU MODUL VAR -- UC AYRI SORU, HEPSI KURATORUN.

1. SORU: "varyantlarin amino asit karsilastirmasinda hatali bir kiyas var;
   TM-align mi daha iyi olur, yoksa katalitik demire yakin kalintilardaki
   farklar gibi aktif-bolge temelli bir sey mi?"

   Kusku hakli ve iki ayri hata birbirine karismis durumda:
     (a) Sitede varyantin yanina yazilan `leaf.median_identity` VARYANT ICI
         medyan ikili kimliktir (recursive_homogenize.py, HOMOGENEITY_TARGET
         0,70): "bu varyantin uyeleri birbirine ne kadar benziyor" demek.
         VARYANTLAR ARASI bir sayi DEGIL. Iki varyanti bu iki sayiyi yan yana
         koyarak karsilastirmak dogrudan yanlistir, cunku olculen sey ayni
         nicelik bile degil.
     (b) Varyantlar arasi global kimlik dogru hesaplansa bile, bu ailede
         substrat secimini global benzerlik belirlemiyor: reference_pairs.csv
         farkli substratli kuratorlu ciftlerin %99,8 kimlige kadar ciktigini
         gosterdi. Yani global kimlik soruya zayif bir cevap.

   YAPISAL SUPERPOZISYON (TM-align) BILINCLI OLARAK YAPILMADI. Elimizdeki
   yapilarin 46'si AlphaFold ongorusu ve bu ailenin omurgasi tipler arasinda
   bile neredeyse ayni; bir superpozisyon agirlikla bu ortak omurgayi olcer,
   islevi olcmez. Ustelik varyantlarin kendi modelleri yok -- superpozisyon
   ancak TIP temsilcileri arasinda yapilabilirdi, soru ise varyantlar arasi.

   Bunun yerine KOLON temelli olculdu: active_site.py'nin her tip icin
   kristalden (17 tip) ya da dogrulanmis ongorulen modelden (46 tip)
   cikardigi, katalitik demirin 8 A komsulugundaki hizalama kolonlari alinir;
   her varyantin bu kolonlardaki konsensus kalintisi okunur; kimlik YALNIZCA
   bu kolonlar uzerinden hesaplanir ve ayni yontemle hesaplanmis global
   kimligin YANINA konur. Boylece iki sayi gercekten karsilastirilabilir
   olur (ayni hizalama, ayni bosluk kurali, ayni konsensus).

   OLCULEN SONUC. 35 tipte 646 varyant ve 13.636 varyant cifti olculdu.
   Ham haliyle ceb protein ortalamasindan daha korunmus gorunuyor: ortalama
   ceb uyusmazligi global uyusmazligin 0,73 katinda ve 35 tipin 35'inde oran
   1'in altinda. AMA BU SAYI YANILTICI VE SEBEBI OLCULDU. 8 A cebi iki
   katalitik histidini, demir karboksilatini ve kopru aspartatini iceriyor;
   bu kolonlar TANIM GEREGI degismez, cunku onlari tasimak bir girisi
   "dogrulanmis RO alpha" yapan testin kendisi (ro_motif). Bu dort kolon
   cikarilinca oran 0,73'ten 0,91'e ciktigi gibi, 35 tipin yalnizca 23'unde
   1'in altinda kaliyor (isaret testi p=0,09, yani anlamsiz). Ayni buyuklukte
   RASTGELE kolon kumelerine karsi da ayni sonuc: ligandlarla ceb null'un
   %1,5'lik diliminde, ligandsiz ceb %37'lik diliminde -- yani siradan bir
   kolon kumesi.

   AMA HAVUZLANMIS SAYI DA YANILTICI, VE SEBEBI OLCULDU
   (by_global_identity_band). Korunma kontrasti, karsilastirilan iki dizi
   birbirinden uzaklastikca SONUYOR; tek bir oran olcumun calismadigi
   rejimleri de ortalamaya katiyor. Bu veride olculen ciftlerin %96'si
   %60'in ALTINDA global kimlikte, %56'si ise ampirik olarak olculmus
   kimlik tabaninin (%27,4; farkli tipten ciftlerin 95. yuzdeligi) altinda
   -- yani havuzlanmis 0,91'in cogu, akraba olmayan enzimlerin bile ayni
   kadar benzestigi bir bolgede hesaplanmis. Farkli tipten ciftlerde oran
   dogrudan 1,03 olcularak doygunluk seviyesi GOSTERILDI.
   Ayni oran global kimlik bantlarinda yeniden hesaplandiginda sonuc
   tersine doniyor: 0,80-0,90 bandinda 0,15 (10/10 tip), 0,70-0,80'de 0,29
   (10/10), 0,60-0,70'te 0,35 (18/18); ucu de BH duzeltmesinden sonra
   q<0,05. Bu, ayni bantlarda cekilen ayni buyuklukteki RASTGELE kolon
   kumelerinde GORULMUYOR (onlar 0,87-0,95 veriyor), ve bantlarin farkli
   tipler icermesinden de gelmiyor: iki bolgede de cifti olan 19 tipin
   19'unda kendi olculebilir ciftlerindeki oran kendi doygun ciftlerindeki
   orandan dusuk (p=3,8e-6). SONUC: ligand olmayan ceb kalintilari,
   karsilastirma gorulebilecek kadar yakin oldugunda proteinin geri
   kalanindan iki ila yedi kat daha yavas degisiyor. Havuzlanmis null sonuc
   bir GUC problemi degil, bir KOMPOZISYON problemiydi.
   Cekince: ayakta kalan bantlar alti bandin arasindan secildi; BH yanlis
   kesif oranini kontrol eder ama bunu on-kayitli bir teste cevirmez.
   Iddiayi tasiyan sey uc KOMSU bandin ayni yonde olmasi, siralamanin
   argumanin ongordugu sonumleme ile ayni gitmesi, rastgele kolon
   kumelerinin bunu uretmemesi ve ayni tiplerin kendilerine karsi ayni seyi
   gostermesi.
   Global ve ceb kimligi arasindaki iliski gercek ve her tipte ayni yonde
   (havuzlanmis Spearman 0,69; tip basina medyan 0,74; 24/24 tipte pozitif).
   Kuratorun tarif ettigi "global benzer ama cepte farkli" cifti VAR ama
   NADIR: %85 ve uzeri global kimlikteki 32 ciftin yalnizca 1'i ceb
   kolonlarinin >=3'unde ayriliyor (KshA15#24 x KshA15#36, global %92,7,
   16 kolonda 3 fark). Sebebi de olculdu ve onemli: varyantlar zaten ham
   dizide CD-HIT ile ayrildigi icin, ayni tipin iki FARKLI varyantinin
   global olarak cok benzer olmasi tanim geregi nadirdir. Yani bu quadrant
   bos degil, ORNEKLEM tarafindan kurutulmus.

2. SORU: "transpozon tur ozelinde mi enzim ozelinde mi?" (ROADMAP 53)

   Ikili sonuc degiskeni (10 kb icinde regex transpozonu var/yok) icin iki
   faktorun acikladigi VARYANS ayristirildi. kxN capraz tabloda Cramer's
   V'nin KARESI tam olarak ikili degiskenin aciklanan varyans oranidir
   (eta-kare), yani iki faktor ayni birimde karsilastirilabilir. Her olcum
   PERMUTASYON null'una karsi okunur, cunku 508 cinsli bir faktor hicbir
   iliski olmasa da yuksek eta-kare uretir.

   OLCULEN SONUC: ENZIM TIPI, CINSTEN daha cok sey anlatiyor, ve bu uc
   duzeyde de ayni: giris duzeyinde sans ustu 0,117'ye 0,083; cins sabit
   tutuldugunda tipin etkisi 0,244, tip sabit tutuldugunda cinsin etkisi
   0,186; tip x cins hucresine cokertildiginde 0,080'e 0,039. Iki faktor
   KISMEN birbirine karisik (V=0,36) ve bazi tipler tek bir cinsten
   geliyor (GbcA'nin %97'si Pseudomonas), bu yuzden tam ayristirma
   mumkun degil.

   ZORUNLU CEKINCE, operon_relations.py'nin olcumu: bu 'transposon'
   cagrisi GenBank /product METNINE uygulanan bir regex'tir (gene_category,
   method='regex_v1'), dizi kaniti degil. Olctugu sey "metin regex'e uyuyor
   mu", "orada gercekten hareketli bir element var mi" DEGIL. Ve apacik
   fazlaligin yaklasik dortte biri (ortak OR 1,41 -> 1,31) anotasyon
   yogunlugundan geliyor. Burada anotasyon yogunlugunun kendi etkisi de
   ayrica olculdu (sans ustu eta-kare 0,015) ve iki faktorun ikisinden de
   kucuk cikti; ayrica iki etki yogunluk ceyrekligi icinde yeniden
   hesaplandi. Yine de olculen sey "nerede transpozon ANOTE EDILMIS"tir.

3. SORU: "birbirine cok benzeyen enzimlerin operonlari da benzer olmali;
   operon benzemiyorsa eslesme rastgele olabilir, ayni gen adini tasiyip
   anlamsiz olabilir." (ROADMAP 54)

   Tip ICINDEKI giris ciftleri icin enzim kimligi ile operon gen SIRASI
   benzerligi iliskilendirildi. Gen sirasi olcutu olarak kisa operonun
   uzunluguna normalize edilmis EN UZUN ORTAK ALT DIZI (LCS) secildi:
   sira duyarli, araya giren fazladan geni cezalandirmaz (operon siniri
   150 bp bosluk kuralindan ve anotasyon derinliginden geliyor, biyolojiden
   degil), simetrik ve 0-1 arasi. Duzenleme mesafesi (Levenshtein)
   REDDEDILDI cunku eklemeyi cezalandirir ve burada ekleme agirlikla
   anotasyon derinligidir. Kompozisyon ayri olcuTLDU (operon_gene uzerinden
   Jaccard) ki sira ile kompozisyon birbirinden ayirt edilebilsin.

   OLCULEN SONUC: KURATOR HAKLI. Havuzlanmis Spearman 0,54; 54 tipin
   54'unde pozitif (isaret testi p=1,1e-16); cins carprazi ciftlerde 0,45,
   ayni cinste 0,70. Ve bulgu anotasyon metnine DAYANMIYOR: tek tamamen
   dizi-kanitli olcut olan bilesen ortusmesi (profil HMM ile cagrilan
   beta/ferredoksin/reduktaz) kimlikle 0,70 korelasyon veriyor, yani gen
   sirasindan DAHA GUCLU. Kuratorun korktugu durum GERCEK AMA COK NADIR:
   >=%95 kimlikli yaklasik 2.200 ciftin yalnizca 1'i tipler arasi SANS
   duzeyinde ya da altinda operon benzerligi gosteriyor. Bu nadirlik iki
   anlama geliyor: operon uyusmazligi goruldugunde bakmaya deger, ama genel
   bir SINAMA olarak agirlik tasimaz cunku neredeyse hicbir cift bu testi
   gecemiyor.

NE OLCULMUYOR -- bastan yazmak gerekiyor:
  * Ceb bir MESAFE kumesidir, bagli substratla temas eden kalintilarin
    deneysel listesi degil. 46 tipte yapi ongorudur ve ceb tam olarak yan
    zincir kumesidir, yani en kirilgan kisim.
  * Ceb kimligi varyant KONSENSUSLARI uzerinden hesaplanir; >=3 uyesi
    olmayan varyantin konsensusu yok ve olcume girmez.
  * Operon layout'unun yarisi dizi kaniti (beta/ferredoksin/reduktaz,
    profil HMM), yarisi anotasyon metni (other/dehydrogenase/...). Bu yuzden
    ayni olcum yalnizca DIZI kanitli jetonlarla da tekrarlandi.
  * Nedensellik yok. Hicbir blokta yon iddiasi yoktur.

Cikti: analysis_out/variant_and_operon.json

Kullanim:
    python3 variant_and_operon.py --db roar.sqlite --out-dir analysis_out
"""

import argparse
import datetime
import json
import os
import random
import sqlite3
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

from ro_motif import DEFAULT_MODEL, MODEL_COLUMNS, read_stockholm_matchcols
# Bant duzeltmesi icin AYNI BH uygulamasi: istatistik sayfasi ne kullaniyorsa
# burada da o kullanilir, iki yerde iki ayri duzeltme olmasin.
from stats_overview import bh_adjust

# ---------------------------------------------------------------- SABITLER
# Her esik adlandirilmis ve gerekcesi yaninda. Sessiz kesim yok.

# --- Ortak

# Permutasyon sayisi ve tohum: operon_relations.py ve cooccurrence.py ile
# AYNI, boylece uc modulun null'lari karsilastirilabilir ve kosu yeniden
# uretilebilir.
PERMUTATIONS = 2000
RANDOM_SEED = 20261005

# Sacilim grafigi icin JSON'a yazilan nokta ustu (operon_relations.py ile
# ayni). Tamami yazilsa dosya gereksiz buyur; secim sabit tohumla yapilir.
MAX_SCATTER_POINTS = 3000

# --- 1. soru: ceb kimligi

# Birincil yaricap. active_site.py'nin PRIMARY_RADIUS_A degeri ile AYNI
# olmak zorunda: 5,0 A "kesin olan" (demir ligandlari), 8,0 A "ilgili olan"
# (substrat cebi). Kuratorun verdigi ornek de 18 ceb pozisyonundan
# bahsediyor, ki bu 8 A cebinin buyuklugu (16-24 kolon); 5 A cebi 6-7
# kalintidir ve uzerinden kimlik hesaplamak cok kaba olur.
POCKET_RADIUS_A = 8.0
# Ikinci yaricap yalnizca DUYARLILIK kontrolu olarak hesaplanir.
POCKET_RADIUS_SENSITIVITY_A = 5.0

# Bir varyantin konsensusundan soz edebilmek icin gereken en az uye sayisi.
# variant_signature.py ve active_site.py ile AYNI esik.
MIN_LEAF_MEMBERS = 3

# Bir kolonun varyant konsensusunda kalinti tasimasi icin, varyant uyelerinin
# en az bu kadarinda dolu olmasi gerekir. Altinda kolon '-' sayilir: azinlik
# bir uyenin kalintisini butun varyantin konsensusu gibi sunmak, tam olarak
# bu modulun duzeltmeye calistigi tur bir hatadir.
MIN_CONSENSUS_OCCUPANCY = 0.5

# Bir varyant cifti icin ceb kimligi hesaplanabilsin diye, iki tarafta da
# kalinti tasiyan en az bu kadar ceb kolonu gerekir. 16-24 kolonluk bir
# cepte 10'un altina dusmek "cebin kimligi" degil "cebin gorulen parcasinin
# kimligi" olurdu.
MIN_POCKET_COLUMNS_COMPARED = 10

# Kuratorun tarif ettigi iki ilginc durumun kapilari.
# "Global benzer": kuratorun kendi ornegi %85. Projenin baska yerinde
# kullanilan kimlik bandi kenarlari 95/80/60 (operon_relations.IDENTITY_BANDS)
# oldugu icin 0,85 bu bandin icinde kaliyor ve ayri bir kesim noktasi
# yaratmiyor. Daha gevsek bir kapi da ayrica raporlanir, cunku 0,85 ustunde
# cift sayisi azdir.
HIGH_GLOBAL_IDENTITY = 0.85
RELAXED_GLOBAL_IDENTITY = 0.70   # recursive_homogenize.HOMOGENEITY_TARGET:
                                 # varyant icinde hedeflenen medyan kimlik.
                                 # Iki varyant bunun altindaysa, kumeleme
                                 # onlari zaten "ayni varyant degil" saymistir.
# "Cepte ayrisik" kapisi SAYI uzerinden, oran uzerinden DEGIL. Cep
# buyuklugu tipe gore 16-24 kolon arasinda degisiyor; sabit bir oran
# kapisi kucuk cepte buyuk cepten daha az fark istemek anlamina gelir ve
# kuratorun gercek ornegini (KshA15, 16 kolonda 3 fark = %81 kimlik) tam
# olarak disarida birakiyordu. Uc fark esigi secildi cunku bir ya da iki
# farkin tamami tek bir hizalama kaymasindan cikabilir (motif testi zaten
# +-2 kolon penceresi kullaniyor).
MIN_POCKET_MISMATCHES = 3
# Oran yalnizca ALT SAYIM olarak raporlanir, kapi olarak degil.
POCKET_DIVERGENT_MAX_IDENTITY = 0.80

# "Ceb protein ortalamasindan daha mi korunmus" sorusunun DOGRU null'u.
# Protein ortalamasiyla karsilastirmak yetmez: 20 kolonluk herhangi bir alt
# kume, 426 kolonun ortalamasindan farkli dagilim gosterir. Bu yuzden her tip
# icin ayni buyuklukte RASTGELE kolon kumeleri cekilir ve cebin bu null
# icindeki yuzdeligi raporlanir.
RANDOM_COLUMN_SETS = 200
# Rastgele kolon kumesi yalnizca tipin varyant konsensuslarinda DOLU olan
# kolonlardan cekilir; bos kolonlari havuza koymak null'u yapay olarak
# kolaylastirirdi.
MIN_NULL_COLUMN_OCCUPANCY = 0.9
# "Global ayrisik ama ceb ayni" tarafi icin global kapi: evidence_tiers.py'nin
# "ayni reaksiyon cok olasi" kapisi %60.
LOW_GLOBAL_IDENTITY = 0.60

# Tip basina Spearman hesaplanabilmesi icin gereken en az cift sayisi.
MIN_PAIRS_FOR_CORRELATION = 10

# Profil-kolon kimliginin GERCEK ikili hizalama kimligine denk oldugunu
# sinamak icin ornek buyuklugu. Sabit tohum. Bu kontrol zorunlu: profil
# hizalamasi insertion kolonlarini atar, ve atilan kisim cifte gore
# degisirse kimlik sistematik kayar.
IDENTITY_VALIDATION_PAIRS = 300

# --- 2. soru: transpozon

# gene_category tablosundaki metin tabanli yontemin adi ve kategori.
# operon_relations.py ile AYNI dizgiler; bu modulun butun 2. sorusu bu iki
# dizginin anotasyon metninden geliyor olmasina bagli.
REGEX_METHOD = "regex_v1"
TRANSPOSON_CATEGORY = "transposon"

# Komsuluk penceresi, extract_genomic_context.py ile ayni; yalnizca
# raporlamak icin tekrarlandi.
NEIGHBOUR_WINDOW_BP = 10000

# Hic anote komsusu olmayan girisler analizden CIKARILIR. operon_relations.py
# olctu: 945 giris (%8,3) hic komsusu olmayan kayitlarda ve hepsi okaryot
# mRNA kaydi. Bu girislerde transpozon yakinligi sifir OLMAK ZORUNDA, yani
# paydada durmalari her iki faktorun etkisini de seyreltir.
MIN_ANNOTATED_NEIGHBOURS = 1

# Bir cins ya da tipin "tipin etkisi cins icinde" / "cinsin etkisi tip
# icinde" tabakasina girmesi icin gereken en az giris sayisi. Daha azi
# havuzlanmis eta-kareye gurultuden baska sey katmaz.
MIN_STRATUM_ENTRIES = 30

# Anotasyon yogunlugu tabakalari VERIDEN cikarilir (ceyreklikler), elle
# yazilmaz: pencere doluluk dagilimi veri setine ozgudur ve sabit bir kesim
# noktasi baska bir surumde anlamini yitirir. operon_relations.py ile ayni.
DENSITY_QUANTILES = (25, 50, 75)

# Bir tipin ya da cinsin betimleyici tabloda gosterilmesi icin en az giris.
MIN_CELL_ENTRIES = 20
# Ve en az kac FARKLI cins: 40 girisi olan ama tek cinsten gelen bir hucre
# tek bir dizileme projesidir, bulgu degil (operon_relations.py ile ayni).
MIN_CELL_GENERA = 3

# --- 3. soru: operon benzerligi

# Gen SIRASI karsilastirilabilsin diye operonda en az bu kadar gen olmali.
# 1 genli operon yalnizca RO'nun kendisidir (11.422'nin 2.479'u), 2 genli
# operonda capaya gore tek bir sira vardir. Uc genden once "sira" diye bir
# sey yok.
MIN_OPERON_GENES = 3

# Tip icinde cift ornekleme ustu. operon_relations.py ile AYNI deger ve ayni
# gerekce: amac tip ICINDEKI cift dagilimini kestirmek, duzgun ornekleme
# bunu yanlilik katmadan yapar.
MAX_PAIRS_PER_TYPE = 2000

# Tipler ARASI cift sayisi: operon benzerliginin SANS duzeyi icin. Bu taban
# zorunlu, cunku her layout '[alpha]' iceriyor ve 'other' jetonu operonlarin
# %90'inda var -- yani LCS'in sifirdan cok yukarida bir tabani var.
CHANCE_LEVEL_PAIRS = 20000

# "Protein olarak neredeyse ayni" kapisi: %95, evidence_tiers.py'deki
# 'characterized' kademesinin kapisi.
HIGH_ENZYME_IDENTITY = 0.95
RELAXED_ENZYME_IDENTITY = 0.90

# Dizi kanitli jetonlar. Bunlar komsu CDS'in PROTEIN DIZISI profil HMM'e
# taranarak bulundu (build_operons.py); geri kalan jetonlar GenBank urun
# metninin regex'i. Ayni olcum bir kez de yalnizca bu jetonlarla yapilir ki
# sonucun anotasyon metnine ne kadar dayandigi gorulebilsin.
SEQUENCE_EVIDENCE_TOKENS = frozenset({
    "[alpha]", "beta", "ferredoxin", "reductase", "alpha_other",
    "rieske_other"})
# Metin kanitli jetonlarin yerine konan joker.
TEXT_TOKEN_PLACEHOLDER = "x"

# Capa cevresinde dogrudan pozisyon uyusmasi icin pencere: RO'dan iki yone
# kac gen. +-3 secildi cunku operonlarin medyan uzunlugu 3 gendir ve daha
# genis bir pencere agirlikla bos pozisyon sayar.
ANCHOR_WINDOW = 3

# --- 1. soru: GLOBAL KIMLIK TABAKALARI
#
# NEDEN TABAKA ZORUNLU. Havuzlanmis oran, olcumun CALISAMADIGI iki rejimi de
# ortalamaya katiyor:
#   * cok benzer ciftlerde hicbir yerde neredeyse hicbir sey degismiyor; ceb
#     ve zemin uyusmazligi ikisi de sifira yakin ve oran iki kucuk sayinin
#     kararsiz bolumu oluyor;
#   * cok uzak ciftlerde ikisi de zemin uyusmazlik oranina DOYUYOR ve oran
#     tanim geregi 1'e suruluyor;
#   * ikisinin arasindaki bantta hem ayirt edilecek kadar fark var hem de
#     doygunluk yok -- etki gercekse YALNIZCA orada gorunur.
# Bu modulun kendi docstring'i doyma problemini yarim gormustu ("kimlikler
# tavanina yakin duruyor, oranlari farki sikistirir" dedigi icin uyusmazlik
# oranlarini kullaniyor) ama tabakalamayi hic yapmamisti. Ayni oran, ayni
# ceb, ayni global kimlik her bantta yeniden hesaplanir; boylece bant
# sonucu havuzlanmis sayiyla DOGRUDAN karsilastirilabilir olur.
#
# BANT SINIRLARI keyfi degil: 0,95 ve 0,90 bu modulde zaten tanimli olan
# HIGH_ENZYME_IDENTITY / RELAXED_ENZYME_IDENTITY kapilari (evidence_tiers.py
# 'characterized' kademesi), 0,60 ise LOW_GLOBAL_IDENTITY, yani
# evidence_tiers.py'nin "ayni reaksiyon cok olasi" kapisi. Aradaki 0,80 ve
# 0,70 esit genislikte ara basamaklar. Alti bant TEK BIR AILE olarak
# Benjamini-Hochberg ile duzeltilir (stats_overview.bh_adjust): alti banttan
# en iyisini duzeltmesiz bildirmek, bu veritabaninin istatistik sayfasinin
# uyardigi hatanin tam kendisi olurdu.
GLOBAL_IDENTITY_BANDS = (
    (HIGH_ENZYME_IDENTITY, 1.01, ">= 0.95"),
    (RELAXED_ENZYME_IDENTITY, HIGH_ENZYME_IDENTITY, "0.90 - 0.95"),
    (0.80, RELAXED_ENZYME_IDENTITY, "0.80 - 0.90"),
    (0.70, 0.80, "0.70 - 0.80"),
    (LOW_GLOBAL_IDENTITY, 0.70, "0.60 - 0.70"),
    (0.00, LOW_GLOBAL_IDENTITY, "< 0.60"))

# EN ALT BANDIN ICI, YALNIZCA BETIMLEYICI. Ciftlerin buyuk cogunlugu 0,60'in
# altinda oldugu icin tek bir '< 0,60' satiri doygunlugun NEREDE bastigini
# gizler. Bu satirlar BH ailesine GIRMEZ: bant sayisini artirmak duzeltmeyi
# zayiflatir, ve bu alt bolme hipotez sinamak icin degil tabani GOSTERMEK
# icin var. Hicbirine p ya da q yazilmaz.
GLOBAL_IDENTITY_SUB_BANDS = (
    (0.50, 0.60, "0.50 - 0.60"),
    (0.40, 0.50, "0.40 - 0.50"),
    (0.30, 0.40, "0.30 - 0.40"),
    (0.00, 0.30, "< 0.30"))

# Bir bandin TAVANA carptiginin olcutu: cift basina, zemin uyusmazlik
# oraniyla BEKLENEN ceb uyusmazligi sayisi (karsilastirilan ceb kolonu x
# global uyusmazlik orani). 16-20 kolonluk bir cepte bu sayi 1'in altindaysa
# gozlenen sifir uyusmazlik korunmanin kaniti DEGIL -- sans eseri zaten en
# olasi sonuc odur ve oran hicbir sey olcmuyordur.
MIN_EXPECTED_POCKET_MISMATCHES = 1.0

# Bir bandin TABANA oturdugunun olcutu: ciftlerinin yarisindan fazlasi
# ampirik olculmus kimlik tabaninin altindaysa, bandin ortalamasi doygun
# olcumlerin ortalamasidir.
MAX_SATURATED_SHARE_FOR_A_USABLE_BAND = 0.5

# Doygunluk tavani ASSERT EDILMEZ, OLCULUR -- operon_relations.py'nin kimlik
# tabanini olctugu yontemin aynisi. Orada 'akraba olmadigi bilinen' karsiligi
# farkli tip VE farkli aileden duzenleyici ciftiydi; burada FARKLI ENZIM
# TIPINDEN iki varyant konsensusu, yani bu veride mumkun olan en uzak
# karsilastirma. Gercek ciftlerle AYNI bicimde puanlanir. 95. yuzdelik
# "akraba olmayan iki RO alpha bu kadar benzesebiliyor" demektir; altindaki
# her olcum doygundur.
SATURATION_SAMPLE_PAIRS = 4000
SATURATION_PERCENTILE = 95

# Tabakalanmis blok KENDI uretecini kullanir, ortak py_rng'den cekmez.
# Sebebi bir kez olculdu: question_one'in SONUNDA cagrilan profil kimligi
# dogrulamasi ayni uretecten 300 cift cekiyor, ve bu blok araya girince o
# ornek kayiyordu (rho 0,860 -> 0,849). Yeni bir blok eski bir sayiyi
# degistirmemeli; ayni gerekce main()'deki soru basina tohumlarin
# gerekcesiyle birebir ayni.
BAND_SEED_OFFSET = 11


# ------------------------------------------------------------- YARDIMCILAR

def column_identity(seq_a, seq_b, indices):
    """Verilen kolonlarda yuzde kimlik; BOSLUK KOLONLARI PAYDAYA GIRMEZ.

    Bosluk kurali recursive_homogenize.identity ile birebir ayni, boylece
    burada hesaplanan global kimlik `leaf.median_identity`'nin olculdugu
    nicelikle ayni birimdedir (biri varyantlar ARASI, oteki varyant ICI).
    Donen: (kimlik, karsilastirilan kolon sayisi, uyusmayan kolon sayisi).
    """
    matches = compared = 0
    for i in indices:
        a, b = seq_a[i], seq_b[i]
        if a in "-." or b in "-.":
            continue
        compared += 1
        if a == b:
            matches += 1
    if not compared:
        return None, 0, 0
    return matches / compared, compared, compared - matches


def spearman(x, y):
    """Spearman rho + p; sabit vektorde (tanimsiz) None doner."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or len(set(x.tolist())) < 2 or len(set(y.tolist())) < 2:
        return None, None
    rho, p = stats.spearmanr(x, y)
    if not np.isfinite(rho):
        return None, None
    return float(rho), float(p)


def sign_test(values, reference=0.0, quantity=None):
    """Tip basina olculen bir degerin isaret testi.

    Tip basina bir gozlem = ciftlerin bagimli olmasina karsi en basit ve en
    savunulabilir denetim. Etki buyuklugu olarak medyan ve pozitif oran
    birlikte doner. `quantity` ciktiya yazilir cunku bazi yerlerde test
    edilen sey dogrudan olcum degil, ondan turetilmis bir sayidir ve
    okuyanin neye bakildigini bilmesi gerekir.
    """
    values = [v for v in values if v is not None]
    if not values:
        return None
    above = sum(1 for v in values if v > reference)
    result = stats.binomtest(above, len(values), 0.5)
    out = {"n_types": len(values),
           "n_above_reference": above,
           "reference": reference,
           "median": round(float(np.median(values)), 4),
           "min": round(float(np.min(values)), 4),
           "max": round(float(np.max(values)), 4),
           "p_sign_test": float(result.pvalue)}
    if quantity:
        out["tested_quantity"] = quantity
    return out


def eta_squared(codes, outcome, n_levels):
    """Ikili sonucun bir faktor tarafindan aciklanan varyans orani.

    kx2 capraz tabloda Cramer's V'nin karesi ile ayni sayidir, ama buradaki
    yazim neyin olculdugunu dogrudan soyluyor: gruplar arasi kareler toplami
    / toplam kareler toplami. Iki faktor boylece AYNI birimde okunur.
    """
    n = len(outcome)
    mean = outcome.mean()
    if n == 0 or mean in (0.0, 1.0):
        return None
    sums = np.bincount(codes, weights=outcome, minlength=n_levels)
    counts = np.bincount(codes, minlength=n_levels)
    group_means = np.divide(sums, counts, out=np.zeros(n_levels),
                            where=counts > 0)
    between = float((((group_means - mean) ** 2) * counts).sum())
    return between / (n * mean * (1.0 - mean))


def factor_effect(values, outcome, rng, permutations=PERMUTATIONS):
    """Bir faktorun etkisi, PERMUTASYON null'u ile birlikte.

    Null zorunlu: 508 seviyeli bir faktor hicbir iliski olmasa da yuksek
    eta-kare uretir (beklenen taban kabaca (k-1)/n). Okunabilir olan sayi
    eta-karenin kendisi degil, null'u ASAN kismi.
    """
    levels = sorted(set(values))
    index = {v: i for i, v in enumerate(levels)}
    codes = np.array([index[v] for v in values])
    observed = eta_squared(codes, outcome, len(levels))
    if observed is None:
        return None
    null = np.empty(permutations)
    for i in range(permutations):
        null[i] = eta_squared(codes, rng.permutation(outcome), len(levels))
    mean, sd = float(null.mean()), float(null.std())
    return {"eta_squared": round(observed, 4),
            "permuted_eta_squared_mean": round(mean, 4),
            "permuted_eta_squared_sd": round(sd, 5),
            "eta_squared_above_chance": round(observed - mean, 4),
            "z_against_permutation": (round((observed - mean) / sd, 1)
                                      if sd else None),
            "n_levels": len(levels),
            "n": int(len(outcome))}


def pooled_within_stratum_effect(records, stratum_key, factor_key, rng,
                                 permutations=PERMUTATIONS // 2):
    """Tabaka ICINDE havuzlanmis eta-kare: oteki faktor sabit tutuldugunda.

    "Cinsi biliyorsam, ustune tipi bilmek ne katiyor" sorusunun dogrudan
    karsiligi. Null, sonucu HER TABAKA ICINDE karistirarak kurulur; boylece
    tabaka etkisi null'da da aynen durur ve karsilastirma yalnizca faktorun
    kendi katkisini olcer.
    """
    strata = defaultdict(list)
    for i, record in enumerate(records):
        strata[record[stratum_key]].append(i)
    # Tabaka basina faktor kodlari BIR KEZ cikarilir. Permutasyonun icinde
    # sozluk kurmak 500 seviyeli bir faktorde olcumu dakikalara cikariyordu;
    # matematik ayni, yalnizca kodlama disariya alindi.
    usable = []
    for indices in strata.values():
        outcomes = [records[i]["outcome"] for i in indices]
        labels = [records[i][factor_key] for i in indices]
        if (len(indices) >= MIN_STRATUM_ENTRIES and len(set(labels)) >= 2
                and 0 < sum(outcomes) < len(outcomes)):
            index = {v: i for i, v in enumerate(sorted(set(labels)))}
            usable.append((np.array(indices),
                           np.array([index[v] for v in labels]),
                           len(index)))
    if not usable:
        return None
    outcome = np.array([r["outcome"] for r in records], dtype=float)

    def statistic(values):
        between = total = 0.0
        for indices, codes, n_levels in usable:
            subset = values[indices]
            mean, n = subset.mean(), len(subset)
            if mean <= 0.0 or mean >= 1.0:
                continue
            sums = np.bincount(codes, weights=subset, minlength=n_levels)
            counts = np.bincount(codes, minlength=n_levels)
            group_means = np.divide(sums, counts, out=np.zeros(n_levels),
                                    where=counts > 0)
            between += float((((group_means - mean) ** 2) * counts).sum())
            total += n * mean * (1.0 - mean)
        return between / total if total else None

    observed = statistic(outcome)
    if observed is None:
        return None
    null = np.empty(permutations)
    for b in range(permutations):
        shuffled = outcome.copy()
        for indices, _, _ in usable:
            shuffled[indices] = rng.permutation(shuffled[indices])
        null[b] = statistic(shuffled)
    mean, sd = float(null.mean()), float(null.std())
    return {"held_fixed": stratum_key,
            "varying_factor": factor_key,
            "pooled_eta_squared": round(observed, 4),
            "permuted_mean": round(mean, 4),
            "permuted_sd": round(sd, 5),
            "eta_squared_above_chance": round(observed - mean, 4),
            "z_against_permutation": (round((observed - mean) / sd, 1)
                                      if sd else None),
            "n_strata_used": len(usable),
            "n_entries_used": int(sum(len(indices)
                                      for indices, _, _ in usable))}


def cramers_v(table):
    """Iki kategorik alan arasindaki iliski (karisiklik denetimi icin)."""
    table = np.asarray(table, dtype=float)
    if min(table.shape) < 2 or table.sum() == 0:
        return 0.0
    chi2 = stats.chi2_contingency(table, correction=False)[0]
    return float(np.sqrt(chi2 / (table.sum() * (min(table.shape) - 1))))


def longest_common_subsequence(a, b):
    """Iki jeton dizisinin en uzun ortak alt dizisinin uzunlugu."""
    previous = [0] * (len(b) + 1)
    for token in a:
        current = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if token == b[j - 1]:
                current[j] = previous[j - 1] + 1
            else:
                current[j] = max(current[j - 1], previous[j])
        previous = current
    return previous[len(b)]


def order_similarity(a, b):
    """Kisa operona normalize edilmis LCS: 0-1, sira duyarli, simetrik.

    NEDEN LCS. Operonun siniri 150 bp bosluk kuralindan ve kaydi gonderenin
    kac geni anote ettiginden geliyor, biyolojiden degil. Duzenleme mesafesi
    araya giren geni CEZALANDIRIR, yani agirlikla anotasyon derinligini
    olcer. LCS ise ortak sirayi arar ve fazladan geni atlayabilir.
    Normalizasyon KISA operona gore yapilir, cunku uzun operonun fazlasi
    kisa olanda zaten gozlenemez.
    """
    shorter = min(len(a), len(b))
    if not shorter:
        return None
    return longest_common_subsequence(a, b) / shorter


def jaccard(a, b):
    """Kompozisyon benzerligi. IKI TARAF DA BOSSA None doner, 1 DONMEZ.

    "Iki operonun da dizi kanitli bileseni yok" ifadesi onlarin benzer
    oldugunu soylemez; 1 dondurmek bos veriyi mukemmel uyum gibi gosterirdi.
    """
    union = a | b
    return len(a & b) / len(union) if union else None


def anchor_agreement(a_map, b_map, window=ANCHOR_WINDOW):
    """Capa (RO = pozisyon 0) cevresinde dogrudan pozisyon uyusmasi.

    LCS'e tamamlayici: LCS kaymayi bagislar, bu olcum bagislamaz. Ikisi
    birlikte "ayni genler ayni sirada mi" ile "ayni genler ayni YERDE mi"
    sorularini ayirir.
    """
    compared = agreed = 0
    for position in range(-window, window + 1):
        if position == 0:
            continue
        left, right = a_map.get(position), b_map.get(position)
        if left is None or right is None:
            continue
        compared += 1
        if left == right:
            agreed += 1
    if not compared:
        return None, 0
    return agreed / compared, compared


def sample_pairs(items, max_pairs, rng):
    """Tip icinde cift ornekleme; tumu az ise tumu, degilse duzgun ornek."""
    items = sorted(items)
    n = len(items)
    if n < 2:
        return []
    total = n * (n - 1) // 2
    if total <= max_pairs:
        return [(items[i], items[j])
                for i in range(n) for j in range(i + 1, n)]
    chosen = set()
    while len(chosen) < max_pairs:
        i, j = rng.randrange(n), rng.randrange(n)
        if i != j:
            chosen.add((min(i, j), max(i, j)))
    return [(items[i], items[j]) for i, j in sorted(chosen)]


def scatter_sample(points, rng, limit=MAX_SCATTER_POINTS):
    """JSON'a yazilacak nokta kumesi; sabit tohumla secilir."""
    if len(points) <= limit:
        return points
    return rng.sample(points, limit)


def describe(values):
    """Bir sayi dizisinin ozeti; JSON'da tekrar eden blok."""
    array = np.asarray([v for v in values if v is not None], dtype=float)
    if not array.size:
        return None
    return {"n": int(array.size),
            "mean": round(float(array.mean()), 4),
            "median": round(float(np.median(array)), 4),
            "p25": round(float(np.percentile(array, 25)), 4),
            "p75": round(float(np.percentile(array, 75)), 4),
            "min": round(float(array.min()), 4),
            "max": round(float(array.max()), 4)}


# -------------------------------------------------- 1. SORU: CEB KIMLIGI

def load_pocket_columns(active_site_path, radius):
    """active_site.json'dan tip basina ceb kolonlari.

    Kristal varsa KRISTAL kazanir; yoksa ongorulen model kullanilir. Kaynak
    her tip icin ciktida yazar, cunku ikisi ayni kanit degildir.
    """
    with open(active_site_path, encoding="utf-8") as handle:
        payload = json.load(handle)
    key = "%.1f" % radius
    columns, source = {}, {}
    for record in payload.get("structures", []):
        if record.get("status") == "ok":
            columns[record["type"]] = list(
                record["pocket"][key]["alignment_columns"])
            source[record["type"]] = {
                "evidence": "crystal structure, iron-centred pocket",
                "pdb_id": record.get("pdb_id"),
                "reference_name": record.get("reference_name")}
    for record in payload.get("predicted_structures", []):
        if record.get("status") == "ok" and record["type"] not in columns:
            columns[record["type"]] = list(
                record["pocket"][key]["alignment_columns"])
            source[record["type"]] = {
                "evidence": ("predicted model, pocket around the catalytic "
                             "triad centroid because the model carries no "
                             "metal"),
                "uniprot_accession": record.get("uniprot_accession"),
                "global_mean_plddt": record.get("global_mean_plddt"),
                "mean_pocket_plddt": record["pocket"][key].get(
                    "mean_pocket_plddt"),
                "reference_name": record.get("reference_name")}
    validation = payload.get("predicted_pocket_transfer_validation", {})
    return columns, source, validation.get("summary", {})


def ligand_columns():
    """ro_motif'in metal ligandi kolonlari -- neredeyse degismezler."""
    model = MODEL_COLUMNS[DEFAULT_MODEL]
    cols = {c for c, _, _ in model["rieske"]}
    cols |= {c for c, _, _ in model["catalytic"]}
    cols.add(model["bridging"][0])
    return cols


def variant_consensus(members, aligned, length):
    """Varyantin hizalama konsensusu; seyrek kolonlar '-' kalir."""
    out = []
    for i in range(length):
        counter = Counter(aligned[m][i] for m in members
                          if aligned[m][i] not in "-.")
        if sum(counter.values()) / len(members) < MIN_CONSENSUS_OCCUPANCY:
            out.append("-")
        else:
            out.append(counter.most_common(1)[0][0])
    return "".join(out)


def validate_profile_identity(con, aligned, leaves_by_cluster, rng):
    """Profil-kolon kimligi GERCEK ikili hizalama kimligine denk mi.

    Profil hizalamasi insertion kolonlarini atar. Atilan kisim cifte gore
    degisirse "kimlik" sistematik kayar ve butun 1. soru bu sayinin uzerine
    kuruldugu icin sinanmasi zorunlu. Ornek, ayni tipin FARKLI
    varyantlarindan alinir, yani olcumun gercekten kullanildigi yerden.
    """
    try:
        from Bio import Align
    except ImportError:
        return {"status": "skipped",
                "reason": "biopython is not installed, so the pairwise "
                          "alignment cross-check could not be run"}
    aligner = Align.PairwiseAligner(scoring="blastp")
    aligner.mode = "global"
    sequences = dict(con.execute(
        "SELECT candidate_id, sequence FROM ro "
        "WHERE is_confirmed = 1 AND sequence IS NOT NULL AND sequence <> ''"))
    candidates = []
    for cluster, leaves in leaves_by_cluster.items():
        leaf_ids = sorted(leaves)
        if len(leaf_ids) < 2:
            continue
        for i in range(len(leaf_ids)):
            for j in range(i + 1, len(leaf_ids)):
                a = leaves[leaf_ids[i]][0]
                b = leaves[leaf_ids[j]][0]
                if a in sequences and b in sequences:
                    candidates.append((a, b))
    if not candidates:
        return {"status": "skipped", "reason": "no comparable pair found"}
    sample = (rng.sample(candidates, IDENTITY_VALIDATION_PAIRS)
              if len(candidates) > IDENTITY_VALIDATION_PAIRS else candidates)
    indices = list(range(len(next(iter(aligned.values())))))
    profile, pairwise = [], []
    for a, b in sample:
        column_value = column_identity(aligned[a], aligned[b], indices)[0]
        if column_value is None:
            continue
        alignment = aligner.align(sequences[a], sequences[b])[0]
        matches = compared = 0
        for x, y in zip(alignment[0], alignment[1]):
            if x != "-" and y != "-":
                compared += 1
                if x == y:
                    matches += 1
        if not compared:
            continue
        profile.append(column_value)
        pairwise.append(matches / compared)
    if len(profile) < 10:
        return {"status": "skipped", "reason": "too few scoreable pairs"}
    differences = np.array(profile) - np.array(pairwise)
    rho, p = spearman(profile, pairwise)
    # Karar METINDEN degil SAYIDAN cikar: olcum degisirse cumle de degisir.
    median_gap = float(np.median(np.abs(differences)))
    if rho is not None and rho >= 0.9 and median_gap <= 0.03:
        verdict = ("the profile-column identity and the pairwise alignment "
                   "identity agree closely, so the column measure can stand "
                   "in for a real alignment")
    elif rho is not None and rho >= 0.75 and median_gap <= 0.06:
        verdict = ("the profile-column identity tracks the pairwise alignment "
                   "identity well enough to rank pairs, with a typical "
                   "absolute gap of %.1f points and a bias of %+.1f points. "
                   "It is usable for comparisons and should not be quoted as "
                   "an exact identity." % (100 * median_gap,
                                           100 * float(differences.mean())))
    else:
        verdict = ("the profile-column identity does NOT reproduce the "
                   "pairwise alignment identity well (rho %s, typical "
                   "absolute gap %.1f points); every number built on it "
                   "should be treated as provisional"
                   % ("%.2f" % rho if rho is not None else "undefined",
                      100 * median_gap))
    return {
        "status": "measured",
        "what_this_checks": ("whether identity computed over the profile "
                             "match columns equals identity from an "
                             "independent pairwise global alignment of the "
                             "raw sequences; the profile drops insertion "
                             "columns, so a pair-dependent loss would shift "
                             "the measure systematically"),
        "n_pairs": len(profile),
        "sampling": ("one member from each variant, for variant pairs inside "
                     "the same type, drawn with seed %d" % RANDOM_SEED),
        "spearman_rho": round(rho, 4) if rho is not None else None,
        "p": p,
        "mean_difference_profile_minus_pairwise": round(
            float(differences.mean()), 4),
        "median_absolute_difference": round(
            float(np.median(np.abs(differences))), 4),
        "max_absolute_difference": round(float(np.abs(differences).max()), 4),
        "verdict": verdict}


def measured_divergence_ceiling(consensus, owner, band_inputs, all_indices,
                                py_rng):
    """Doygunluk tavaninin AMPIRIK olcumu -- iddia degil, olcum.

    operon_relations.measured_identity_floor ile AYNI mantik ve ayni
    gerekce: taban literaturden alinmis bir sayi olarak varsayilmaz, bu
    veride olculur. 'Akraba olmadigi bilinen cift'in buradaki karsiligi
    FARKLI enzim tipinden iki varyant konsensusudur; bu veri setinde
    yapilabilecek en uzak karsilastirma budur.

    Puanlama gercek ciftlerle birebir ayni: global kimlik ayni butun match
    kolonlari uzerinden, ceb kimligi iki tipin KENDI ligandsiz ceb
    kolonlari uzerinden. Her cift iki kez puanlanir (bir kez her tipin
    cebiyle), cunku tek bir tipin cebini secmek asimetrik ve keyfi olurdu.
    """
    types = sorted(band_inputs)
    leaves = sorted(leaf for leaf in consensus if leaf in owner)
    if len(types) < 2 or len(leaves) < 2:
        return None
    globals_, pockets, seen = [], [], set()
    attempts = 0
    while (len(globals_) < SATURATION_SAMPLE_PAIRS
           and attempts < 100 * SATURATION_SAMPLE_PAIRS):
        attempts += 1
        a, b = py_rng.choice(leaves), py_rng.choice(leaves)
        if owner[a] == owner[b]:
            continue
        key = (a, b) if a < b else (b, a)
        if key in seen:
            continue
        seen.add(key)
        glob = column_identity(consensus[a], consensus[b], all_indices)[0]
        if glob is None:
            continue
        scored = False
        for cluster in (owner[a], owner[b]):
            value, compared, _ = column_identity(
                consensus[a], consensus[b],
                band_inputs[cluster]["non_ligand_idx"])
            if value is not None and compared >= MIN_POCKET_COLUMNS_COMPARED:
                pockets.append(value)
                scored = True
        if scored:
            globals_.append(glob)
    if not globals_ or not pockets:
        return None
    global_divergence = 1.0 - float(np.mean(globals_))
    pocket_divergence = 1.0 - float(np.mean(pockets))
    floor = float(np.percentile(np.asarray(globals_, dtype=float),
                                SATURATION_PERCENTILE))
    return {
        "what_was_sampled": (
            "variant consensus pairs drawn from DIFFERENT enzyme types, "
            "scored exactly like the real within-type pairs: global identity "
            "over the same alignment columns, pocket identity over each "
            "type's own pocket columns with the metal ligands removed. These "
            "are the most distant comparisons this data allows, so they "
            "measure where both numbers stop carrying information."),
        "n_cross_type_pairs": len(globals_),
        "n_pocket_scores": len(pockets),
        "global_identity": describe(globals_),
        "pocket_identity_excluding_metal_ligands": describe(pockets),
        "mean_global_divergence": round(global_divergence, 4),
        "mean_pocket_divergence_excluding_metal_ligands": round(
            pocket_divergence, 4),
        "divergence_ratio_at_saturation": (
            round(pocket_divergence / global_divergence, 4)
            if global_divergence > 0 else None),
        "floor_percentile": SATURATION_PERCENTILE,
        "measured_global_identity_floor": round(floor, 4),
        "how_to_read": (
            "The divergence ratio between unrelated enzymes is the value the "
            "statistic is driven to when the measure saturates. A real band "
            "whose ratio equals it has measured nothing. The floor is the "
            "%dth percentile of cross-type global identity: a within-type "
            "pair below it is no more similar than two enzymes of different "
            "types, so neither its global nor its pocket divergence carries "
            "information." % SATURATION_PERCENTILE)}


def _per_type_band_ratios(records):
    """Bant icinde TIP BASINA uyusmazlik orani.

    Toplama sirasi havuzlanmis satirla birebir ayni tutuldu: bant icindeki
    ciftler de bagimsiz degil, bu yuzden once tip basina tek bir oran
    hesaplanir, isaret testi ancak ondan sonra tipler uzerinde yapilir.
    Donen: {tip: oran} ve {tip: o bandaki cift sayisi}.
    """
    grouped = defaultdict(list)
    for record in records:
        grouped[record["type"]].append(record)
    ratios, counts = {}, {}
    for cluster, rows in grouped.items():
        global_divergence = 1.0 - float(np.mean(
            [r["global_identity"] for r in rows]))
        free = [r["pocket_identity_excluding_metal_ligands"] for r in rows
                if r["pocket_identity_excluding_metal_ligands"] is not None]
        counts[cluster] = len(rows)
        if not free or global_divergence <= 0:
            continue
        ratios[cluster] = ((1.0 - float(np.mean(free))) / global_divergence)
    return ratios, counts


def _band_row(records, label, lower, upper, band_inputs, consensus, floor,
              py_rng, with_null=True):
    """Tek bir bandin satiri. Hicbir yeni olcum YOK, ayni oran yeniden.

    Her satir yalnizca oran degil, orani olusturan IKI sayiyi da ayri ayri
    tasiyor (ortalama ceb uyusmazligi ve ortalama global uyusmazlik), cunku
    tavan ve taban iddia edilmemeli, okuyanin dogrudan gorebilmesi gerekir.
    `with_null` False ise ayni buyuklukte rastgele kolon kumesi null'u
    atlanir (betimleyici alt bantlar icin).
    """
    if not records:
        return None
    ratios, counts = _per_type_band_ratios(records)
    values = sorted(ratios.values())
    free = [r["pocket_identity_excluding_metal_ligands"] for r in records
            if r["pocket_identity_excluding_metal_ligands"] is not None]
    mismatches = [r["pocket_mismatches_excluding_metal_ligands"]
                  for r in records
                  if r["pocket_mismatches_excluding_metal_ligands"]
                  is not None]
    # Beklenen ceb uyusmazligi: cebin de zemin oraninda degistigi varsayimi.
    # Bu sayi bandin TAVANA carpip carpmadigini dogrudan gosterir.
    scored_columns = [r["pocket_columns_compared_excluding_metal_ligands"]
                      * (1.0 - r["global_identity"]) for r in records
                      if r["pocket_columns_compared_excluding_metal_ligands"]]
    expected = float(np.mean(scored_columns)) if scored_columns else 0.0
    saturated_share = float(np.mean(
        [1.0 if r["global_identity"] <= floor else 0.0 for r in records])
    ) if floor is not None else None
    row = {
        "band": label,
        "global_identity_from": round(lower, 4),
        "global_identity_to": round(min(upper, 1.0), 4),
        "n_pairs": len(records),
        "n_types_in_band": len(counts),
        "n_types_with_a_ratio": len(values),
        "n_types_with_at_least_5_pairs_in_the_band": sum(
            1 for n in counts.values() if n >= 5),
        "mean_global_divergence": round(
            1.0 - float(np.mean([r["global_identity"] for r in records])), 4),
        "mean_pocket_divergence_excluding_metal_ligands": (
            round(1.0 - float(np.mean(free)), 4) if free else None),
        "per_type_divergence_ratio": describe(values),
        "expected_pocket_mismatches_per_pair_at_the_background_rate": round(
            expected, 2),
        "observed_pocket_mismatches_per_pair": (
            round(float(np.mean(mismatches)), 2) if mismatches else None),
        "share_of_pairs_with_an_identical_non_ligand_pocket": (
            round(float(np.mean([1.0 if m == 0 else 0.0
                                 for m in mismatches])), 4)
            if mismatches else None),
        "share_of_pairs_at_or_below_the_measured_identity_floor": (
            round(saturated_share, 4) if saturated_share is not None
            else None)}

    # Isaret testi, havuzlanmis satirla AYNI sekilde: test edilen nicelik
    # 1 eksi oran, yani pozitif deger "ceb daha yavas degisiyor" demek.
    summary = describe(values)
    test = sign_test([1.0 - v for v in values]) if values else None
    if test and summary:
        row["sign_test_against_one"] = dict(
            test,
            tested_quantity="1 minus the divergence ratio computed inside "
                            "this band without the metal ligand columns")
        # Isaret testinin kendi etki buyuklugu: 1'in altindaki tiplerin
        # orani ve TAM binom guven araligi. Tek basina medyan oran yeterli
        # degil, cunku ince bir bantta medyan tek bir tipten gelebilir.
        attempt = stats.binomtest(test["n_above_reference"], test["n_types"],
                                  0.5)
        interval = attempt.proportion_ci(confidence_level=0.95,
                                         method="exact")
        row["effect_size"] = {
            "median_divergence_ratio": summary["median"],
            "median_slowdown_percent": round(
                100.0 * (1.0 - summary["median"]), 1),
            "share_of_types_where_the_pocket_changes_more_slowly": round(
                test["n_above_reference"] / test["n_types"], 4),
            "share_ci95_exact": [round(float(interval.low), 4),
                                 round(float(interval.high), 4)],
            "note": ("the share and its exact binomial interval are the "
                     "effect size the sign test actually tests; the median "
                     "ratio is the size of the difference")}

    # Bu bant hangi rejimde? Tavan ve taban testi SAYIYLA, metinle degil.
    reasons = []
    if expected < MIN_EXPECTED_POCKET_MISMATCHES:
        reasons.append(
            "ceiling: at this identity the background rate predicts only "
            "%.2f mismatched pocket columns per pair, so an unchanged pocket "
            "is the most likely outcome by chance alone and the ratio is a "
            "quotient of two near-zero numbers" % expected)
    if (saturated_share is not None
            and saturated_share > MAX_SATURATED_SHARE_FOR_A_USABLE_BAND):
        reasons.append(
            "floor: %.0f %% of the pairs in this band are at or below the "
            "measured identity floor, where two enzymes of different types "
            "already score as high, so both divergences are saturated"
            % (100 * saturated_share))
    row["usable"] = not reasons
    row["why_not_usable"] = reasons or None

    if with_null:
        row["same_size_random_column_set_null"] = _band_null(
            records, band_inputs, consensus, ratios, py_rng)
    return row


def _band_null(records, band_inputs, consensus, ratios, py_rng):
    """Bandin KENDI null'u: ayni bantta ayni buyuklukte rastgele kolonlar.

    Bu kontrol zorunlu, cunku bant SECIMI global kimlige gore yapiliyor ve
    ceb kolonlari global kimligin de icinde. Rastgele ayni buyuklukte kolon
    kumeleri AYNI bantta ayni orani veriyorsa bant etkisi secimden gelir,
    cebin kendisinden gelmez. Null havuzu ligand kolonlarini icermez, cunku
    karsilastirilan ceb de icermiyor.
    """
    grouped = defaultdict(list)
    for record in records:
        grouped[record["type"]].append(record)
    observed, null_medians, percentiles = [], [], []
    for cluster, rows in grouped.items():
        if cluster not in ratios:
            continue
        inputs = band_inputs[cluster]
        size = len(inputs["non_ligand_idx"])
        pool = inputs["null_pool"]
        global_divergence = 1.0 - float(np.mean(
            [r["global_identity"] for r in rows]))
        if size < 1 or len(pool) <= size or global_divergence <= 0:
            continue
        drawn_ratios = []
        for _ in range(RANDOM_COLUMN_SETS):
            drawn = [c - 1 for c in py_rng.sample(pool, size)]
            values = [column_identity(consensus[r["variant_a"]],
                                      consensus[r["variant_b"]], drawn)[0]
                      for r in rows]
            values = [v for v in values if v is not None]
            if values:
                drawn_ratios.append(
                    (1.0 - float(np.mean(values))) / global_divergence)
        if not drawn_ratios:
            continue
        array = np.asarray(drawn_ratios)
        observed.append(ratios[cluster])
        null_medians.append(float(np.median(array)))
        percentiles.append(float((array <= ratios[cluster]).mean()))
    if not observed:
        return None
    return {
        "n_types_with_a_null": len(observed),
        "n_random_column_sets_per_type": RANDOM_COLUMN_SETS,
        "median_random_set_ratio": round(float(np.median(null_medians)), 4),
        "median_pocket_ratio": round(float(np.median(observed)), 4),
        "median_pocket_percentile_in_null": round(
            float(np.median(percentiles)), 4),
        "n_types_with_pocket_below_the_5th_percentile_of_its_null": sum(
            1 for value in percentiles if value <= 0.05),
        "how_to_read": ("if the random sets reproduce the pocket's ratio "
                        "inside the band, the band effect is an artefact of "
                        "selecting pairs on global identity; if they stay "
                        "near 1 while the pocket does not, the effect "
                        "belongs to the pocket")}


def conservation_by_identity_band(pairs, band_inputs, consensus, ceiling,
                                  py_rng):
    """Ayni oran, GLOBAL KIMLIGE gore tabakalanmis.

    Havuzlanmis sonucun (0,91, 35 tipin 23'u, p=0,09) bir GUC problemi mi
    yoksa bir KOMPOZISYON problemi mi oldugunu bu blok ayirir. Yeni bir
    olcum tanimlanmaz: ayni ceb kimligi, ayni global kimlik, ayni uyusmazlik
    orani, ayni tip-basina toplama, ayni isaret testi. Degisen tek sey,
    testin hangi cift kumesinde yapildigi.
    """
    floor = (ceiling or {}).get("measured_global_identity_floor")
    rows = []
    for lower, upper, label in GLOBAL_IDENTITY_BANDS:
        records = [r for r in pairs if lower <= r["global_identity"] < upper]
        row = _band_row(records, label, lower, upper, band_inputs, consensus,
                        floor, py_rng)
        if row:
            rows.append(row)
    # Alti bant TEK AILE olarak duzeltilir. Bir banti duzeltmesiz bildirmek,
    # alti deneme yapip en iyisini secmek olurdu.
    tested = [r for r in rows if r.get("sign_test_against_one")]
    for row, q in zip(tested, bh_adjust(
            [r["sign_test_against_one"]["p_sign_test"] for r in tested])):
        row["benjamini_hochberg_q"] = float(q)
        row["survives_bh_at_0.05"] = bool(q < 0.05)
    for row in rows:
        row.setdefault("benjamini_hochberg_q", None)
        row.setdefault("survives_bh_at_0.05", None)

    # Betimleyici alt bantlar: tabanin NEREDE bastigini gostermek icin.
    sub_rows = []
    for lower, upper, label in GLOBAL_IDENTITY_SUB_BANDS:
        records = [r for r in pairs if lower <= r["global_identity"] < upper]
        row = _band_row(records, label, lower, upper, band_inputs, consensus,
                        floor, py_rng, with_null=False)
        if row:
            row.pop("benjamini_hochberg_q", None)
            sub_rows.append(row)

    # Referans satiri: AYNI cift kumesinde tabakalanmamis test. Havuzlanmis
    # ile tabakalanmis karsilastirmasi boylece birebir ayni kumede olur.
    pooled = _band_row(pairs, "all pairs, unstratified", 0.0, 1.0,
                       band_inputs, consensus, floor, py_rng, with_null=False)
    if pooled:
        pooled.pop("benjamini_hochberg_q", None)
        pooled.pop("survives_bh_at_0.05", None)

    # TIP ICINDE eslesmis kontrol: hangi tiplerin hangi banta dustugu
    # sonucu suruklememis olsun. Ayni tipin olculebilir penceredeki orani ile
    # doygun bolgedeki orani karsilastirilir, yani tip kimligi tamamen sabit.
    window_low, window_high = LOW_GLOBAL_IDENTITY, HIGH_ENZYME_IDENTITY
    window_ratios, _ = _per_type_band_ratios(
        [r for r in pairs
         if window_low <= r["global_identity"] < window_high])
    saturated_ratios, _ = _per_type_band_ratios(
        [r for r in pairs if r["global_identity"] < window_low])
    shared = sorted(set(window_ratios) & set(saturated_ratios))
    paired = None
    if len(shared) >= 5:
        differences = [saturated_ratios[c] - window_ratios[c]
                       for c in shared]
        test = sign_test(differences)
        try:
            wilcoxon_p = float(stats.wilcoxon(differences).pvalue)
        except ValueError:
            wilcoxon_p = None
        paired = {
            "question": ("inside one type, is the ratio lower among its "
                         "measurable pairs than among its saturated pairs?"),
            "why": ("the bands hold different types, so a band difference "
                    "could be a difference between types. This comparison "
                    "holds the type fixed and only changes which of its own "
                    "pairs are counted."),
            "measurable_window": "%.2f <= global identity < %.2f" % (
                window_low, window_high),
            "saturated_region": "global identity < %.2f" % window_low,
            "n_types_with_pairs_in_both": len(shared),
            "median_ratio_in_the_window": round(
                float(np.median([window_ratios[c] for c in shared])), 4),
            "median_ratio_in_the_saturated_region": round(
                float(np.median([saturated_ratios[c] for c in shared])), 4),
            "sign_test": dict(test or {},
                              tested_quantity="saturated ratio minus window "
                                              "ratio, so a positive value "
                                              "means the pocket signal is "
                                              "stronger in the window"),
            "wilcoxon_p": wilcoxon_p,
            "per_type": [
                {"type": c,
                 "n_pairs_in_the_window": sum(
                     1 for r in pairs
                     if r["type"] == c
                     and window_low <= r["global_identity"] < window_high),
                 "ratio_in_the_window": round(window_ratios[c], 4),
                 "n_pairs_saturated": sum(
                     1 for r in pairs
                     if r["type"] == c
                     and r["global_identity"] < window_low),
                 "ratio_saturated": round(saturated_ratios[c], 4)}
                for c in sorted(shared, key=lambda c: window_ratios[c])]}

    usable_rows = [r for r in rows if r["usable"]]
    survivors = [r for r in rows if r.get("survives_bh_at_0.05")]
    reading = _band_reading(rows, pooled, ceiling, survivors, usable_rows,
                            paired)
    return {
        "question": ("the pooled test averages over regimes where the "
                     "measurement cannot work. Does the same statistic "
                     "behave differently when the variant pairs are "
                     "stratified by how similar the two sequences are "
                     "overall?"),
        "why_stratify": (
            "Conservation contrast decays as the two sequences being "
            "compared get less similar. One number over all variant pairs "
            "therefore mixes three regimes: pairs so similar that nothing "
            "varies anywhere, where the ratio is an unstable quotient of two "
            "near-zero numbers; pairs so distant that both measures have "
            "saturated at the background mismatch rate, where the ratio is "
            "driven to 1 by construction; and a middle band with enough "
            "variation to detect a difference and not enough to saturate. If "
            "the effect is real it is visible only in the middle."),
        "statistic": ("exactly the statistic of the row above: divergence "
                      "ratio = (1 - mean pocket identity) / (1 - mean global "
                      "identity), with the metal ligand columns removed, one "
                      "ratio per enzyme type, then a sign test against 1 "
                      "across types. Only the set of pairs entering it "
                      "changes, so every band is directly comparable to the "
                      "pooled number."),
        "band_edges_why": (
            "0.95 and 0.90 are this module's own HIGH_ENZYME_IDENTITY and "
            "RELAXED_ENZYME_IDENTITY gates, taken from the 'characterized' "
            "tier in evidence_tiers.py; 0.60 is LOW_GLOBAL_IDENTITY, the "
            "'same reaction very likely' gate in the same file. 0.80 and "
            "0.70 are equal steps between them. The edges were not chosen "
            "after looking at the result."),
        "multiple_testing": (
            "The six bands are treated as one family and corrected with "
            "Benjamini-Hochberg (stats_overview.bh_adjust). The descriptive "
            "sub-bands below 0.60 are deliberately outside that family and "
            "carry no p or q: they exist to show where the floor begins, not "
            "to test anything. A band selected from six is a weaker claim "
            "than a pre-registered one even after correction."),
        "measured_saturation": ceiling,
        "reference_row_unstratified": pooled,
        "bands": rows,
        "descriptive_detail_inside_the_lowest_band": sub_rows,
        "within_type_paired_check": paired,
        "n_bands_tested": len(rows),
        "n_bands_surviving_bh": len(survivors),
        "reading": reading}


def _band_reading(rows, pooled, ceiling, survivors, usable_rows, paired):
    """Bandin KARARI, metinden degil sayidan.

    Bu cumle elle yazilmaz: veri degisirse cumle de degismek zorunda.
    """
    scored = [r for r in rows if r["per_type_divergence_ratio"]]
    if not scored or not pooled or not pooled["per_type_divergence_ratio"]:
        return "no band held enough pairs to be measured"
    lowest = min(scored, key=lambda r: r["n_pairs"])
    biggest = max(scored, key=lambda r: r["n_pairs"])
    parts = []
    parts.append(
        "The pooled test is not underpowered; it is dominated by pairs the "
        "measure cannot read. %d of the %d variant pairs, %.0f %% of them, "
        "sit in the %s band, and the pooled ratio (%.2f) is that band's "
        "ratio (%.2f) almost exactly."
        % (biggest["n_pairs"], pooled["n_pairs"],
           100 * biggest["n_pairs"] / pooled["n_pairs"], biggest["band"],
           pooled["per_type_divergence_ratio"]["median"],
           biggest["per_type_divergence_ratio"]["median"]))
    if ceiling:
        parts.append(
            "The saturation level was measured, not assumed: variant "
            "consensuses from different enzyme types align at %.0f %% global "
            "identity and %.0f %% pocket identity, a divergence ratio of "
            "%.2f, so unrelated enzymes already produce the ratio of 1 that "
            "the statistic is supposed to be tested against. The %dth "
            "percentile of that distribution is %.0f %% identity, and %.0f %% "
            "of all the variant pairs measured here lie at or below it."
            % (100 * ceiling["global_identity"]["mean"],
               100 * ceiling["pocket_identity_excluding_metal_ligands"][
                   "mean"],
               ceiling["divergence_ratio_at_saturation"],
               ceiling["floor_percentile"],
               100 * ceiling["measured_global_identity_floor"],
               100 * (pooled.get(
                   "share_of_pairs_at_or_below_the_measured_identity_floor")
                   or 0.0)))
    parts.append(
        "At the other end the measure is blind rather than saturated: in the "
        "%s band the background rate predicts %.2f mismatched pocket columns "
        "per pair, so the %.2f ratio there is a quotient of two near-zero "
        "numbers and the sign test cannot reach significance on %d types "
        "whatever the biology."
        % (lowest["band"],
           lowest["expected_pocket_mismatches_per_pair_at_the_background_"
                  "rate"],
           lowest["per_type_divergence_ratio"]["median"],
           lowest["n_types_with_a_ratio"]))
    if survivors:
        strongest = min(survivors,
                        key=lambda r: r["per_type_divergence_ratio"]["median"])
        powered = max(survivors, key=lambda r: r["n_pairs"])
        ratios = [r["per_type_divergence_ratio"]["median"] for r in survivors]
        parts.append(
            "In between, the effect is present and large. %d of the %d bands "
            "survive Benjamini-Hochberg at q < 0.05: %s. The best powered of "
            "them is %s, with %d pairs over %d types, where the non-ligand "
            "pocket diverges at %.2f of the background rate, %d of %d types "
            "agree and q = %.1e; the lowest ratio of the three is %.2f in "
            "%s. Across the surviving bands the pocket changes between %.1f "
            "and %.1f times more slowly than the rest of the protein."
            % (len(survivors), len(rows),
               ", ".join(r["band"] for r in survivors), powered["band"],
               powered["n_pairs"],
               powered["sign_test_against_one"]["n_types"],
               powered["per_type_divergence_ratio"]["median"],
               powered["sign_test_against_one"]["n_above_reference"],
               powered["sign_test_against_one"]["n_types"],
               powered["benjamini_hochberg_q"],
               strongest["per_type_divergence_ratio"]["median"],
               strongest["band"],
               1.0 / max(ratios) if max(ratios) else float("inf"),
               1.0 / min(ratios) if min(ratios) else float("inf")))
        # SIRALAMA IDDIA EDILMEZ, SINANIR. Tavana carpan bantlar disarida
        # birakilir (orada oran iki sifira yakin sayinin bolumu, siralamaya
        # girmesi anlamsiz). Sonra artan siralamanin NEREDEN itibaren
        # bozulmadan devam ettigi bulunur; monotonluk yalnizca gercekten
        # varsa yazilir.
        ordered = [r for r in scored
                   if r["expected_pocket_mismatches_per_pair_at_the_"
                        "background_rate"] >= MIN_EXPECTED_POCKET_MISMATCHES]
        ordered.sort(key=lambda r: -r["global_identity_from"])
        series = [r["per_type_divergence_ratio"]["median"] for r in ordered]
        start = 0
        for index in range(len(series)):
            rest = series[index:]
            if all(a <= b for a, b in zip(rest, rest[1:])):
                start = index
                break
        tail = ordered[start:]
        listing = "; ".join("%s %.2f" % (
            r["band"], r["per_type_divergence_ratio"]["median"])
            for r in tail)
        if start == 0 and len(tail) > 2:
            parts.append(
                "The ordering is the decay itself, not a single lucky band: "
                "the ratio rises without exception as identity falls across "
                "every band the measure can read (%s)." % listing)
        elif len(tail) > 2:
            out_of_order = ordered[:start]
            parts.append(
                "The ordering mostly follows the decay the argument "
                "predicts: the ratio rises without exception from %s down to "
                "%s (%s). The exception is %s, which holds %d pairs over %d "
                "types, and a band that thin is expected to sit out of order."
                % (tail[0]["band"], tail[-1]["band"], listing,
                   ", ".join(r["band"] for r in out_of_order),
                   sum(r["n_pairs"] for r in out_of_order),
                   max(r["n_types_with_a_ratio"] for r in out_of_order)))
        nulls = [r["same_size_random_column_set_null"] for r in survivors
                 if r.get("same_size_random_column_set_null")]
        if nulls:
            parts.append(
                "The band effect is not an artefact of selecting pairs on "
                "global identity. Random column sets of the same size, drawn "
                "inside the same bands, give ratios of %s, while the pocket "
                "gives %s and sits at percentile %s of its own null."
                % (", ".join("%.2f" % n["median_random_set_ratio"]
                             for n in nulls),
                   ", ".join("%.2f" % n["median_pocket_ratio"]
                             for n in nulls),
                   ", ".join("%.3f" % n["median_pocket_percentile_in_null"]
                             for n in nulls)))
        if paired and paired["sign_test"]:
            parts.append(
                "And it is not a difference between types either. For the %d "
                "types that have pairs in both regions, the ratio is lower "
                "among their own measurable pairs than among their own "
                "saturated pairs in %d of %d cases (median %.2f against "
                "%.2f, sign test p = %.1e), with the type held fixed."
                % (paired["n_types_with_pairs_in_both"],
                   paired["sign_test"]["n_above_reference"],
                   paired["sign_test"]["n_types"],
                   paired["median_ratio_in_the_window"],
                   paired["median_ratio_in_the_saturated_region"],
                   paired["sign_test"]["p_sign_test"]))
        parts.append(
            "So stratifying rescues the effect: the non-ligand pocket of "
            "these enzymes IS measurably more conserved than the rest of the "
            "protein wherever the comparison is close enough to see it. The "
            "honest caveat is that the "
            "surviving bands were selected from six, so the correction "
            "controls the false discovery rate but does not make this a "
            "pre-registered test; what carries the claim is that three "
            "adjacent bands agree, that the ordering follows the decay the "
            "argument predicts, that random column sets do not reproduce it, "
            "and that the same types show it against themselves.")
    else:
        parts.append(
            "No band survives Benjamini-Hochberg at q < 0.05, so the pooled "
            "null result was not a power problem: there is no identity range "
            "in which the non-ligand pocket is measurably more conserved "
            "than the rest of the protein.")
    return " ".join(parts)


def question_one(con, aligned, active_site_path, py_rng):
    """Varyantlari CEB kolonlarinda karsilastir."""
    length = len(next(iter(aligned.values())))
    columns, column_source, transfer_summary = load_pocket_columns(
        active_site_path, POCKET_RADIUS_A)
    columns_inner, _, _ = load_pocket_columns(
        active_site_path, POCKET_RADIUS_SENSITIVITY_A)
    ligands = ligand_columns()

    stored = {leaf_id: {"size": size, "median_identity": median,
                        "top_genera": genera}
              for leaf_id, size, median, genera in con.execute(
                  "SELECT leaf_id, size, median_identity, top_genera "
                  "FROM leaf")}
    leaves_by_cluster = defaultdict(dict)
    for candidate_id, cluster, leaf_id in con.execute(
            "SELECT candidate_id, cluster, leaf_id FROM ro_leaf"):
        if candidate_id in aligned:
            leaves_by_cluster[cluster].setdefault(leaf_id, []).append(
                candidate_id)

    # Konsensusu olan varyantlar ve karsilastirilabilir tipler.
    usable = {}
    for cluster, leaves in leaves_by_cluster.items():
        big = {leaf_id: members for leaf_id, members in leaves.items()
               if len(members) >= MIN_LEAF_MEMBERS}
        if cluster in columns and len(big) >= 2:
            usable[cluster] = big

    consensus = {}
    for cluster, leaves in usable.items():
        for leaf_id, members in leaves.items():
            consensus[leaf_id] = variant_consensus(members, aligned, length)

    all_indices = list(range(length))
    per_type, per_variant, pairs = {}, [], []
    # Tabakalanmis test ayni ceb kolonlarini ve AYNI null havuzunu yeniden
    # kullanir; burada tip basina saklanir ki paralel bir uygulama yazilmasin.
    band_inputs = {}
    for cluster in sorted(usable):
        pocket = sorted(columns[cluster])
        inner = sorted(columns_inner.get(cluster, []))
        pocket_idx = [c - 1 for c in pocket]
        inner_idx = [c - 1 for c in inner]
        non_ligand = [c - 1 for c in pocket if c not in ligands]
        leaf_ids = sorted(usable[cluster])

        type_pairs = []
        for i in range(len(leaf_ids)):
            for j in range(i + 1, len(leaf_ids)):
                a, b = leaf_ids[i], leaf_ids[j]
                glob, glob_n, _ = column_identity(
                    consensus[a], consensus[b], all_indices)
                pock, pock_n, pock_mismatch = column_identity(
                    consensus[a], consensus[b], pocket_idx)
                if glob is None or pock is None:
                    continue
                if pock_n < MIN_POCKET_COLUMNS_COMPARED:
                    continue
                free, free_n, free_mismatch = column_identity(
                    consensus[a], consensus[b], non_ligand)
                tight = column_identity(consensus[a], consensus[b],
                                        inner_idx)[0] if inner_idx else None
                genus_a = (stored.get(a, {}).get("top_genera") or "?"
                           ).split(":")[0]
                genus_b = (stored.get(b, {}).get("top_genera") or "?"
                           ).split(":")[0]
                record = {
                    "type": cluster, "variant_a": a, "variant_b": b,
                    "global_identity": round(glob, 4),
                    "global_columns_compared": glob_n,
                    "pocket_identity": round(pock, 4),
                    "pocket_columns_compared": pock_n,
                    "pocket_columns_total": len(pocket),
                    "pocket_mismatches": pock_mismatch,
                    "pocket_identity_excluding_metal_ligands": (
                        round(free, 4) if free is not None else None),
                    # Bant tabakalamasi bu iki sayiyi KULLANIR: bir bandin
                    # tavana carpip carpmadigi, cift basina beklenen ceb
                    # uyusmazligi sayisindan okunur.
                    "pocket_columns_compared_excluding_metal_ligands": (
                        free_n if free is not None else None),
                    "pocket_mismatches_excluding_metal_ligands": (
                        free_mismatch if free is not None else None),
                    "pocket_identity_at_5A": (
                        round(tight, 4) if tight is not None else None),
                    "dominant_genus_a": genus_a,
                    "dominant_genus_b": genus_b,
                    "same_dominant_genus": int(genus_a == genus_b)}
                type_pairs.append(record)
        if not type_pairs:
            continue
        pairs.extend(type_pairs)

        # Kolon basina varyant konsensus durumlari: cebin ne kadar oynak
        # oldugunun dogrudan olcumu.
        per_column = {}
        for column in pocket:
            states = Counter()
            for leaf_id in leaf_ids:
                residue = consensus[leaf_id][column - 1]
                if residue not in "-.":
                    states[residue] += 1
            per_column[str(column)] = {
                "n_variants_with_residue": sum(states.values()),
                "variant_consensus_states": dict(
                    sorted(states.items(), key=lambda kv: (-kv[1], kv[0]))),
                "n_distinct_states": len(states),
                "invariant_across_variants": (len(states) == 1
                                              if states else None),
                "is_metal_ligand_column": column in ligands}
        invariant = [int(c) for c, v in per_column.items()
                     if v["invariant_across_variants"] is True]
        variable = [int(c) for c, v in per_column.items()
                    if v["invariant_across_variants"] is False]

        globals_ = [r["global_identity"] for r in type_pairs]
        pockets = [r["pocket_identity"] for r in type_pairs]
        free_values = [r["pocket_identity_excluding_metal_ligands"]
                       for r in type_pairs
                       if r["pocket_identity_excluding_metal_ligands"]
                       is not None]
        global_divergence = 1.0 - float(np.mean(globals_))
        pocket_divergence = 1.0 - float(np.mean(pockets))
        ratio = (pocket_divergence / global_divergence
                 if global_divergence > 0 else None)
        free_ratio = None
        if free_values and global_divergence > 0:
            free_ratio = (1.0 - float(np.mean(free_values))) / global_divergence
        rho, p = spearman(globals_, pockets)

        # Ayni buyuklukte rastgele kolon kumeleri: cebin gercek null'u.
        # Havuz, tipin varyant konsensuslarinda yeterince dolu kolonlar.
        pool = [c for c in range(1, length + 1)
                if (sum(1 for leaf_id in leaf_ids
                        if consensus[leaf_id][c - 1] not in "-.")
                    / len(leaf_ids)) >= MIN_NULL_COLUMN_OCCUPANCY]
        # Ligandsiz karsilastirma icin AYRI bir null: boyut da havuz da
        # eslesmeli. Ligandsiz ceb daha az kolon tasiyor ve ligand kolonlari
        # neredeyse degismez oldugu icin havuzda kalmalari null'u yapay
        # olarak korunmus gosterir.
        pool_free = [c for c in pool if c not in ligands]
        band_inputs[cluster] = {"non_ligand_idx": non_ligand,
                                "null_pool": pool_free}

        def null_ratios_for(size, draw_pool):
            out = []
            if size < 1 or len(draw_pool) <= size or global_divergence <= 0:
                return out
            for _ in range(RANDOM_COLUMN_SETS):
                drawn = [c - 1 for c in py_rng.sample(draw_pool, size)]
                values = [column_identity(consensus[r["variant_a"]],
                                          consensus[r["variant_b"]], drawn)[0]
                          for r in type_pairs]
                values = [v for v in values if v is not None]
                if values:
                    out.append(
                        (1.0 - float(np.mean(values))) / global_divergence)
            return out

        null_ratios = null_ratios_for(len(pocket), pool)
        null_free = null_ratios_for(len(non_ligand), pool_free)
        null_block = None
        if null_ratios and ratio is not None:
            array = np.asarray(null_ratios)
            free_block = None
            if null_free and free_ratio is not None:
                free_array = np.asarray(null_free)
                free_block = {
                    "n_random_column_sets": len(null_free),
                    "columns_drawn_per_set": len(non_ligand),
                    "pool_size": len(pool_free),
                    "random_set_divergence_ratio_mean": round(
                        float(free_array.mean()), 4),
                    "random_set_divergence_ratio_sd": round(
                        float(free_array.std()), 4),
                    "pocket_percentile_in_null": round(
                        float((free_array <= free_ratio).mean()), 4)}
            null_block = {
                "n_random_column_sets": len(null_ratios),
                "columns_drawn_per_set": len(pocket),
                "pool_size": len(pool),
                "random_set_divergence_ratio_mean": round(
                    float(array.mean()), 4),
                "random_set_divergence_ratio_sd": round(
                    float(array.std()), 4),
                "pocket_percentile_in_null": round(
                    float((array <= ratio).mean()), 4),
                "excluding_metal_ligands": free_block,
                "note": ("the percentile is the share of same-size random "
                         "column sets whose divergence ratio is at or below "
                         "the pocket's, so a small number means the pocket "
                         "is unusually conserved for a set of this size. The "
                         "second block repeats it for the pocket without its "
                         "metal ligand columns, drawing the random sets from "
                         "a pool that also excludes them, so both size and "
                         "pool are matched.")}

        # Degismezligin ne kadari BILGI, ne kadari ORNEK AZLIGI. Iki varyanti
        # olan bir tipte bir kolonun degismez cikmasi neredeyse kacinilmazdir;
        # bu yuzden karar varyant sayisina bakmadan verilemez.
        if not variable:
            reading = ("the pocket is invariant across every variant of this "
                       "type, so pocket identity is 1 for every pair and the "
                       "measure carries no information here")
        elif len(leaf_ids) < 5:
            reading = ("%d of %d pocket columns vary, but the type has only "
                       "%d variants with a consensus, so an invariant column "
                       "here means mostly that there was little opportunity "
                       "to vary; this type should not be read as evidence of "
                       "pocket conservation"
                       % (len(variable), len(pocket), len(leaf_ids)))
        elif len(variable) <= 2:
            reading = ("only %d of %d pocket columns vary across %d variants, "
                       "so pocket identity moves in very few steps and small "
                       "differences should not be over-read"
                       % (len(variable), len(pocket), len(leaf_ids)))
        else:
            reading = ("%d of %d pocket columns vary across %d variants, so "
                       "the measure is informative for this type"
                       % (len(variable), len(pocket), len(leaf_ids)))

        per_type[cluster] = {
            "pocket_column_source": column_source[cluster],
            "n_pocket_columns": len(pocket),
            "pocket_columns": pocket,
            "metal_ligand_pocket_columns": sorted(
                c for c in pocket if c in ligands),
            "n_variants_with_consensus": len(leaf_ids),
            "n_variants_in_type": len(leaves_by_cluster[cluster]),
            "n_pairs": len(type_pairs),
            "global_identity": describe(globals_),
            "pocket_identity": describe(pockets),
            "pocket_identity_excluding_metal_ligands": describe(free_values),
            "invariant_pocket_columns": sorted(invariant),
            "variable_pocket_columns": sorted(variable),
            "n_invariant_pocket_columns": len(invariant),
            "n_variable_pocket_columns": len(variable),
            "pocket_divergence_ratio": (round(ratio, 4)
                                        if ratio is not None else None),
            "pocket_divergence_ratio_excluding_metal_ligands": (
                round(free_ratio, 4) if free_ratio is not None else None),
            "spearman_global_vs_pocket": (round(rho, 4)
                                          if rho is not None else None),
            "spearman_p": p,
            "pairs_scored_for_correlation": len(type_pairs),
            "random_column_set_null": null_block,
            "information_content": reading,
            "per_column": per_column}

        for leaf_id in leaf_ids:
            signature = []
            present = 0
            for column in pocket:
                residue = consensus[leaf_id][column - 1]
                if residue not in "-.":
                    present += 1
                signature.append("%d%s" % (column, residue))
            mates = [r for r in type_pairs
                     if r["variant_a"] == leaf_id or r["variant_b"] == leaf_id]
            per_variant.append({
                "leaf_id": leaf_id,
                "type": cluster,
                "n_members_in_alignment": len(usable[cluster][leaf_id]),
                "size_in_database": stored.get(leaf_id, {}).get("size"),
                "within_variant_median_identity_stored": stored.get(
                    leaf_id, {}).get("median_identity"),
                "dominant_genus": (stored.get(leaf_id, {}).get("top_genera")
                                   or "?").split(":")[0],
                "pocket_signature": " ".join(signature),
                "n_pocket_columns_with_consensus": present,
                "mean_global_identity_to_other_variants": (
                    round(float(np.mean([r["global_identity"]
                                         for r in mates])), 4)
                    if mates else None),
                "mean_pocket_identity_to_other_variants": (
                    round(float(np.mean([r["pocket_identity"]
                                         for r in mates])), 4)
                    if mates else None)})

    # --- Doygunluk tavani ve global kimlik tabakalari. SIRASI onemli:
    # ikisi de py_rng'den cekiyor ve tip basina null'lardan SONRA cagriliyor,
    # boylece onceki sayilar birebir yeniden uretilebilir kalir.
    owner = {leaf_id: cluster for cluster in band_inputs
             for leaf_id in usable[cluster]}
    band_rng = random.Random(RANDOM_SEED + BAND_SEED_OFFSET)
    ceiling = measured_divergence_ceiling(consensus, owner, band_inputs,
                                          all_indices, band_rng)
    bands_block = conservation_by_identity_band(pairs, band_inputs, consensus,
                                                ceiling, band_rng)

    # --- Global ve ceb kimligi arasindaki iliski, uc duzeyde
    globals_ = [r["global_identity"] for r in pairs]
    pockets = [r["pocket_identity"] for r in pairs]
    rho_all, p_all = spearman(globals_, pockets)
    cross = [r for r in pairs if not r["same_dominant_genus"]]
    rho_cross, p_cross = spearman([r["global_identity"] for r in cross],
                                  [r["pocket_identity"] for r in cross])
    same = [r for r in pairs if r["same_dominant_genus"]]
    rho_same, p_same = spearman([r["global_identity"] for r in same],
                                [r["pocket_identity"] for r in same])
    per_type_rho = [v["spearman_global_vs_pocket"] for v in per_type.values()
                    if v["spearman_global_vs_pocket"] is not None
                    and v["n_pairs"] >= MIN_PAIRS_FOR_CORRELATION]
    ratios = [v["pocket_divergence_ratio"] for v in per_type.values()
              if v["pocket_divergence_ratio"] is not None]
    free_ratios = [v["pocket_divergence_ratio_excluding_metal_ligands"]
                   for v in per_type.values()
                   if v["pocket_divergence_ratio_excluding_metal_ligands"]
                   is not None]

    # --- "Ceb daha mi korunmus" sorusunun KARARI, metinden degil sayidan.
    # Bu blogun sonucu analizi degistirdigi icin cumleyi elle yazmak
    # tehlikeli olurdu: veri degisirse cumle de degismeli.
    ratio_test = sign_test([1.0 - r for r in ratios])
    free_test = sign_test([1.0 - r for r in free_ratios])
    percentiles = [v["random_column_set_null"]["pocket_percentile_in_null"]
                   for v in per_type.values() if v["random_column_set_null"]]
    free_percentiles = [
        v["random_column_set_null"]["excluding_metal_ligands"][
            "pocket_percentile_in_null"]
        for v in per_type.values()
        if v["random_column_set_null"]
        and v["random_column_set_null"]["excluding_metal_ligands"]]
    conservation_reading = (
        "Read the two versions together and the answer changes. With the "
        "metal ligands in, the pocket looks clearly more conserved than the "
        "protein average: median divergence ratio %.2f, below 1 in %d of %d "
        "types, and around percentile %.2f of same-size random column sets. "
        "With the ligands taken out, and with the random sets matched for "
        "both size and pool, the effect mostly disappears: the median ratio "
        "rises to %.2f, it is below 1 in only %d of %d types (sign test "
        "p = %.2g, so not significant), and the pocket now sits around "
        "percentile %.2f of the random sets instead of %.2f. So most "
        "of the apparent extra conservation of the active site comes from the "
        "three or four positions that the pipeline already required to be "
        "present. "
        "This block used to stop there, and stopping there was wrong. A "
        "single ratio over every variant pair averages over identity ranges "
        "in which the measurement cannot work at all, and almost all of the "
        "pairs here lie in one of them. by_global_identity_band repeats the "
        "same statistic inside bands of global identity. %s"
        % (describe(ratios)["median"], ratio_test["n_above_reference"],
           ratio_test["n_types"],
           float(np.median(percentiles)) if percentiles else float("nan"),
           describe(free_ratios)["median"], free_test["n_above_reference"],
           free_test["n_types"], free_test["p_sign_test"],
           float(np.median(free_percentiles))
           if free_percentiles else float("nan"),
           float(np.median(percentiles)) if percentiles else float("nan"),
           (bands_block or {}).get("reading", "")))

    # --- Kuratorun iki ilginc durumu
    def gate(records, label, note):
        return {"label": label, "definition": note, "n_pairs": len(records),
                "n_types": len({r["type"] for r in records}),
                "n_genera": len({r["dominant_genus_a"] for r in records}
                                | {r["dominant_genus_b"] for r in records}),
                "n_also_below_pocket_identity_%.2f"
                % POCKET_DIVERGENT_MAX_IDENTITY: sum(
                    1 for r in records
                    if r["pocket_identity"] <= POCKET_DIVERGENT_MAX_IDENTITY),
                "pairs": sorted(records,
                                key=lambda r: (-r["global_identity"],
                                               r["pocket_identity"]))[:60]}

    strict_similar = [r for r in pairs
                      if r["global_identity"] >= HIGH_GLOBAL_IDENTITY]
    relaxed_similar = [r for r in pairs
                       if r["global_identity"] >= RELAXED_GLOBAL_IDENTITY]
    strict_divergent = [r for r in strict_similar
                        if r["pocket_mismatches"] >= MIN_POCKET_MISMATCHES]
    relaxed_divergent = [r for r in relaxed_similar
                         if r["pocket_mismatches"] >= MIN_POCKET_MISMATCHES]
    identical_pocket_strict = [r for r in pairs
                               if r["global_identity"] <= LOW_GLOBAL_IDENTITY
                               and r["pocket_mismatches"] == 0]
    identical_pocket_relaxed = [r for r in pairs
                                if r["global_identity"]
                                <= RELAXED_GLOBAL_IDENTITY
                                and r["pocket_mismatches"] == 0]

    scatter = scatter_sample(
        [{"type": r["type"], "a": r["variant_a"], "b": r["variant_b"],
          "global": r["global_identity"], "pocket": r["pocket_identity"],
          "mismatches": r["pocket_mismatches"],
          "cross_genus": 1 - r["same_dominant_genus"]} for r in pairs],
        py_rng)

    skipped = sorted(set(leaves_by_cluster) - set(usable))
    big_variants = {c: sum(1 for m in leaves.values()
                           if len(m) >= MIN_LEAF_MEMBERS)
                    for c, leaves in leaves_by_cluster.items()}
    no_columns_but_comparable = [c for c in skipped
                                 if c not in columns
                                 and big_variants.get(c, 0) >= 2]
    no_columns_and_single_variant = [c for c in skipped
                                     if c not in columns
                                     and big_variants.get(c, 0) < 2]
    too_few = [c for c in skipped if c in columns]

    return {
        "question": ("The curator wrote that comparing variants by whole-"
                     "protein amino-acid identity is the wrong comparison, "
                     "and asked whether TM-align or an active-site measure "
                     "would be better. This block compares variants only at "
                     "the alignment columns that line the catalytic iron."),
        "why_not_structural_superposition": (
            "No superposition was attempted, on purpose. 46 of the 63 "
            "pocket definitions come from predicted models, and the backbone "
            "of this family is nearly identical even between types, so a "
            "superposition would mostly measure that shared fold rather than "
            "function. There are also no models for the variants themselves, "
            "only for the reference of each type, so a superposition could "
            "not answer a question about variants at all."),
        "what_was_wrong_with_the_old_comparison": (
            "Two separate faults were conflated. First, the number shown next "
            "to a variant, leaf.median_identity, is the median pairwise "
            "identity WITHIN that variant, that is how similar its own "
            "members are. It is not a between-variant quantity, so placing "
            "two such numbers side by side does not compare the two variants "
            "at all. Second, even a correct between-variant global identity "
            "is a weak measure in this family: curated pairs with different "
            "substrates reach 99.8 % identity, so substrate choice is not "
            "set by global similarity. Both numbers are reported here, "
            "computed on the same alignment with the same gap rule, so they "
            "can finally be read against each other."),
        "method": {
            "pocket_definition": (
                "alignment columns within %.1f A of the catalytic "
                "mononuclear iron, taken from active_site.json: the crystal "
                "structure where one exists, otherwise the AlphaFold model "
                "with the pocket centred on the catalytic triad centroid"
                % POCKET_RADIUS_A),
            "identity_rule": (
                "percent identity over the chosen columns only; columns "
                "where either side has a gap are left out of the "
                "denominator, the same rule recursive_homogenize.py uses, so "
                "the pocket number and the global number are in the same "
                "units"),
            "variant_representation": (
                "each variant is represented by the consensus residue of its "
                "aligned members at every column, requiring at least %d "
                "members and at least %.0f %% occupancy for a column to "
                "carry a residue"
                % (MIN_LEAF_MEMBERS, 100 * MIN_CONSENSUS_OCCUPANCY)),
            "pocket_transfer_validation_from_active_site": transfer_summary},
        "coverage": {
            "types_with_pocket_columns": len(columns),
            "types_with_pocket_from_crystal": sum(
                1 for v in column_source.values()
                if v["evidence"].startswith("crystal")),
            "types_with_pocket_from_predicted_model": sum(
                1 for v in column_source.values()
                if v["evidence"].startswith("predicted")),
            "types_measured_here": len(per_type),
            "variants_measured": len(per_variant),
            "variant_pairs_measured": len(pairs),
            "pocket_columns_per_type": describe(
                [v["n_pocket_columns"] for v in per_type.values()]),
            "types_lost_because_no_pocket_columns_though_comparable": {
                "note": ("these types have at least two variants with a "
                         "consensus but neither a crystal structure nor a "
                         "usable predicted model, so they could have been "
                         "measured and could not be; BmoA is the largest "
                         "loss, with 140 variants over 610 members"),
                "types": no_columns_but_comparable},
            "types_without_pocket_columns_and_without_two_variants":
                no_columns_and_single_variant,
            "types_skipped_fewer_than_two_variants_with_consensus": too_few,
            "why_large_types_are_missing": (
                "The biggest types contribute no variant pair at all, because "
                "their members fall into a single variant: XylX (422 "
                "members), GbcA (336), AntA (184), HcaE (178), PhtAa (126) "
                "each have one variant, so there is nothing to compare "
                "within them. This is a property of the variant clustering, "
                "not of the pocket.")},
        "profile_identity_validation": validate_profile_identity(
            con, aligned, usable, py_rng),
        "relation_between_global_and_pocket_identity": {
            "question": ("does identity at the pocket track identity over the "
                         "whole protein, or is it a different measurement"),
            "levels": [
                {"level": "every variant pair",
                 "n": len(pairs),
                 "n_types": len(per_type),
                 "spearman_rho": round(rho_all, 4) if rho_all else None,
                 "p": p_all,
                 "note": ("pairs inside a type are not independent, so this "
                          "number is the weakest of the three")},
                {"level": "one correlation per type, then a sign test",
                 "detail": sign_test(per_type_rho),
                 "note": ("this is the level that controls for pairs not "
                          "being independent: every type contributes one "
                          "observation")},
                {"level": "only pairs whose two variants have different "
                          "dominant genera",
                 "n": len(cross),
                 "n_types": len({r["type"] for r in cross}),
                 "n_genera": len({r["dominant_genus_a"] for r in cross}
                                 | {r["dominant_genus_b"] for r in cross}),
                 "spearman_rho": (round(rho_cross, 4)
                                  if rho_cross is not None else None),
                 "p": p_cross,
                 "note": ("a pair of variants that share their dominant "
                          "genus can be two halves of one sequencing effort; "
                          "this level removes them")},
                {"level": "only pairs whose two variants share their "
                          "dominant genus",
                 "n": len(same),
                 "spearman_rho": (round(rho_same, 4)
                                  if rho_same is not None else None),
                 "p": p_same}],
            "reading": ("The two measures move together, in every type and at "
                        "every level, so the pocket measure is not "
                        "independent of global identity. It is however "
                        "consistently higher, which is the next block.")},
        "is_the_pocket_more_conserved_than_the_protein": {
            "statistic": ("divergence ratio = (1 - mean pocket identity) / "
                          "(1 - mean global identity), computed per type. "
                          "Below 1 means pocket positions change more slowly "
                          "than the protein average. The ratio is built on "
                          "mismatch rates rather than identities because the "
                          "identities themselves sit close to their ceiling "
                          "and their ratio would compress the difference."),
            "per_type_ratio": describe(ratios),
            "sign_test_against_one": dict(
                ratio_test,
                tested_quantity="1 minus the divergence ratio, so a positive "
                                "value means the pocket changes more slowly "
                                "than the protein average"),
            "per_type_ratio_excluding_metal_ligand_columns":
                describe(free_ratios),
            "sign_test_against_one_excluding_metal_ligand_columns": dict(
                free_test,
                tested_quantity="1 minus the divergence ratio computed "
                                "without the metal ligand columns"),
            "why_the_exclusion_matters": (
                "The 8 A pocket contains the two catalytic histidines, the "
                "iron carboxylate and the bridging aspartate. Those columns "
                "are invariant by definition: carrying them is the test that "
                "made an entry a confirmed RO alpha subunit in the first "
                "place. Leaving them in the pocket makes the pocket look more "
                "conserved for a reason that is circular, so the ratio is "
                "reported both ways and the second number is the one to "
                "believe."),
            "same_size_random_column_set_null": {
                "why": ("Comparing a 20-column subset against the mean of all "
                        "426 columns is not a fair test: any small subset has "
                        "a different sampling distribution. For every type, "
                        "%d random column sets of the same size were drawn "
                        "from the columns that are filled in at least %.0f %% "
                        "of that type's variant consensuses, and the pocket's "
                        "position inside that null is reported."
                        % (RANDOM_COLUMN_SETS,
                           100 * MIN_NULL_COLUMN_OCCUPANCY)),
                "pocket_percentile": describe(
                    [v["random_column_set_null"]["pocket_percentile_in_null"]
                     for v in per_type.values()
                     if v["random_column_set_null"]]),
                "pocket_percentile_excluding_metal_ligands": describe(
                    [v["random_column_set_null"]["excluding_metal_ligands"][
                        "pocket_percentile_in_null"]
                     for v in per_type.values()
                     if v["random_column_set_null"]
                     and v["random_column_set_null"][
                         "excluding_metal_ligands"]]),
                "n_types_with_a_null": sum(
                    1 for v in per_type.values()
                    if v["random_column_set_null"]),
                "how_to_read": ("a percentile near 0 means the pocket is "
                                "more conserved than almost any random set "
                                "of the same size; a percentile near 0.5 "
                                "means the pocket is an ordinary set of "
                                "columns")},
            "pooled_pocket_identity": describe(pockets),
            "pooled_global_identity": describe(globals_),
            "share_of_pairs_with_an_identical_pocket": round(
                float(np.mean([1.0 if r["pocket_mismatches"] == 0 else 0.0
                               for r in pairs])), 4),
            "invariant_column_budget": {
                "total_pocket_columns_across_types": sum(
                    v["n_pocket_columns"] for v in per_type.values()),
                "invariant_across_the_type_variants": sum(
                    v["n_invariant_pocket_columns"]
                    for v in per_type.values()),
                "metal_ligand_columns_across_types": sum(
                    len(v["metal_ligand_pocket_columns"])
                    for v in per_type.values()),
                "types_with_at_least_five_variants": {
                    "n_types": sum(1 for v in per_type.values()
                                   if v["n_variants_with_consensus"] >= 5),
                    "invariant_share": round(float(np.mean(
                        [v["n_invariant_pocket_columns"]
                         / v["n_pocket_columns"]
                         for v in per_type.values()
                         if v["n_variants_with_consensus"] >= 5])), 4)},
                "note": ("An invariant column in a type that has only two or "
                         "three variants mostly records that there was no "
                         "opportunity to vary, so the share is also given for "
                         "the types with at least five variants. Every type's "
                         "own reading is in per_type.information_content.")},
            "by_global_identity_band": bands_block,
            "reading": conservation_reading},
        "discordant_pairs": {
            "question": ("which variants are globally similar but differ at "
                          "the pocket, and which differ globally but have the "
                          "same pocket"),
            "why_these_are_the_interesting_ones": (
                "A pair that is globally similar but differs at several "
                "pocket positions is a candidate for differing specificity. A "
                "pair that differs globally but has an identical pocket "
                "probably performs the same chemistry. Neither statement is "
                "tested here against any measured activity; they are "
                "candidate readings only."),
            "globally_similar_pocket_divergent_strict": gate(
                strict_divergent,
                "global identity at least %.2f and at least %d mismatched "
                "pocket columns" % (HIGH_GLOBAL_IDENTITY,
                                    MIN_POCKET_MISMATCHES),
                "the curator's own example, taken literally"),
            "globally_similar_pocket_divergent_relaxed": gate(
                relaxed_divergent,
                "global identity at least %.2f and at least %d mismatched "
                "pocket columns" % (RELAXED_GLOBAL_IDENTITY,
                                    MIN_POCKET_MISMATCHES),
                "the same idea with the global gate lowered to the variant "
                "homogeneity target, because very few variant pairs reach "
                "0.85"),
            "globally_divergent_pocket_identical_strict": gate(
                identical_pocket_strict,
                "global identity at most %.2f and no mismatched pocket column"
                % LOW_GLOBAL_IDENTITY,
                "the reverse case, with the global gate at the threshold "
                "evidence_tiers.py uses for 'same reaction very likely'"),
            "globally_divergent_pocket_identical_relaxed": gate(
                identical_pocket_relaxed,
                "global identity at most %.2f and no mismatched pocket column"
                % RELAXED_GLOBAL_IDENTITY,
                "the same with the gate at the variant homogeneity target"),
            "n_pairs_at_each_global_gate": {
                "at_least_0.95": sum(1 for r in pairs
                                     if r["global_identity"] >= 0.95),
                "at_least_0.85": len(strict_similar),
                "at_least_0.70": len(relaxed_similar),
                "at_least_0.60": sum(1 for r in pairs
                                     if r["global_identity"] >= 0.60)},
            "why_the_first_quadrant_is_nearly_empty": (
                "Variants are produced by splitting a type with CD-HIT until "
                "the median identity inside a variant reaches 0.70. Two "
                "different variants of the same type are therefore, by "
                "construction, rarely very similar: the clustering split them "
                "precisely because they were not. The quadrant the curator "
                "described is not empty because the biology is absent, it is "
                "nearly empty because the variant partition removes most of "
                "its candidates before the comparison starts.")},
        "sensitivity_to_the_radius": {
            "why": ("active_site.py defines the pocket at two radii and this "
                    "module uses the larger one. The smaller one is reported "
                    "here so the reader can see how much the conclusion "
                    "depends on that choice."),
            "primary_radius_A": POCKET_RADIUS_A,
            "second_radius_A": POCKET_RADIUS_SENSITIVITY_A,
            "pocket_identity_at_the_primary_radius": describe(pockets),
            "pocket_identity_at_the_second_radius": describe(
                [r["pocket_identity_at_5A"] for r in pairs
                 if r["pocket_identity_at_5A"] is not None]),
            "n_pairs_scoreable_at_the_second_radius": sum(
                1 for r in pairs if r["pocket_identity_at_5A"] is not None),
            "share_of_pairs_identical_at_the_second_radius": (round(
                float(np.mean([1.0 if r["pocket_identity_at_5A"] == 1.0
                               else 0.0 for r in pairs
                               if r["pocket_identity_at_5A"] is not None])), 4)
                if any(r["pocket_identity_at_5A"] is not None
                       for r in pairs) else None),
            "columns_per_type_at_the_second_radius": describe(
                [len(v) for v in columns_inner.values()]),
            "reading": ("The 5 A shell is six or seven residues and most of "
                        "them are the metal ligands themselves, so pocket "
                        "identity there is close to 1 for almost every pair "
                        "and the measure has almost no resolution. That is "
                        "why the 8 A shell is used throughout: it is the "
                        "radius at which a bound aromatic substrate could "
                        "still touch a side chain, and it is the only one of "
                        "the two with enough varying positions to separate "
                        "variants.")},
        "per_type": per_type,
        "per_variant": sorted(per_variant,
                              key=lambda r: (r["type"], r["leaf_id"])),
        "scatter": {"note": ("a seeded sample of at most %d variant pairs for "
                             "plotting; the full set stays in the per-type "
                             "summaries" % MAX_SCATTER_POINTS),
                    "seed": RANDOM_SEED,
                    "points": scatter},
        "what_this_cannot_settle": [
            "The pocket is a set of residues within a distance of the metal, "
            "not an experimentally mapped substrate contact list. For 46 of "
            "the 63 types the structure is an AlphaFold model and the pocket "
            "is exactly a set of side chains, which is the part of a "
            "prediction that is least trustworthy.",
            "A pocket difference is not a specificity difference. The "
            "database holds no measured activity for any member, so every "
            "discordant pair listed here is a candidate for an experiment, "
            "not a result.",
            "An identical pocket is not proof of identical chemistry either: "
            "the curated set already contains pairs that differ by a single "
            "residue 15.4 A from the iron and are labelled with different "
            "substrates, which the pocket measure cannot see.",
            "Variants with fewer than %d aligned members have no consensus "
            "and are absent from every number here, and the types with the "
            "most members contribute no pair at all because their members "
            "form a single variant." % MIN_LEAF_MEMBERS,
            "Pocket columns were defined on the reference of each type and "
            "read on its members. For a distant member the column mapping is "
            "only as good as the profile alignment, and the motif test itself "
            "allows a plus-or-minus two column window, so a residue can sit "
            "one column away from where the pocket expects it."]}


# ------------------------------------------- 2. SORU: TRANSPOZON KIMIN OZELLIGI

def question_two(con, rng):
    """Transpozon yakinligi cinsin mi tipin mi ozelligi."""
    rows = con.execute("""
        SELECT r.candidate_id, r.ro_cluster, r.ro_group, p.organism, p.status,
               (SELECT COUNT(*) FROM neighbor nb
                 WHERE nb.candidate_id = r.candidate_id) AS n_neighbours,
               (SELECT COUNT(*) FROM neighbor nb
                  JOIN gene_category c ON c.neighbor_id = nb.neighbor_id
                   AND c.method = ?
                 WHERE nb.candidate_id = r.candidate_id
                   AND c.category = ?) AS transposons
        FROM ro r JOIN replicon p USING(nucleotide_id)
        WHERE r.is_confirmed = 1""",
        (REGEX_METHOD, TRANSPOSON_CATEGORY)).fetchall()

    everything = []
    for (candidate_id, cluster, group, organism, status, n_neighbours,
         transposons) in rows:
        words = (organism or "?").split()
        everything.append({
            "candidate_id": candidate_id,
            "type": cluster or "?",
            "group": group or "?",
            "genus": words[0] if words else "?",
            "species": " ".join(words[:2]) if words else "?",
            "replicon_status": status,
            "n_neighbours": n_neighbours or 0,
            "outcome": int((transposons or 0) > 0)})
    scored = [r for r in everything
              if r["n_neighbours"] >= MIN_ANNOTATED_NEIGHBOURS]
    outcome = np.array([r["outcome"] for r in scored], dtype=float)

    # Anotasyon yogunlugu tabakalari veriden.
    counts = np.array([r["n_neighbours"] for r in scored])
    cuts = np.percentile(counts, DENSITY_QUANTILES)
    for record in scored:
        record["density_quartile"] = "q%d" % int(
            np.searchsorted(cuts, record["n_neighbours"], side="right"))

    effects = {}
    for key, label in (("genus", "genus"), ("species", "species"),
                       ("type", "enzyme type"), ("group", "RO group"),
                       ("density_quartile", "annotation density quartile")):
        effect = factor_effect([r[key] for r in scored], outcome, rng)
        if effect:
            effect["factor"] = label
            effects[key] = effect
    joint = factor_effect([r["type"] + "|" + r["genus"] for r in scored],
                          outcome, rng, permutations=PERMUTATIONS // 4)
    if joint:
        joint["factor"] = "type and genus together, as one cell per pair"
        joint["caution"] = ("there are about five entries per cell, so this "
                            "model overfits badly; it is reported only to "
                            "show that the two factors together explain more "
                            "than the sum of their separate excesses, which "
                            "is what an interaction looks like")
        effects["type_x_genus_cell"] = joint

    held = {
        "type_within_genus": pooled_within_stratum_effect(
            scored, "genus", "type", rng),
        "genus_within_type": pooled_within_stratum_effect(
            scored, "type", "genus", rng),
        "type_within_density_quartile": pooled_within_stratum_effect(
            scored, "density_quartile", "type", rng),
        "genus_within_density_quartile": pooled_within_stratum_effect(
            scored, "density_quartile", "genus", rng)}

    # Karisiklik denetimi: tip ile cins birbirini ne kadar belirliyor.
    types = sorted({r["type"] for r in scored})
    genera = sorted({r["genus"] for r in scored})
    type_index = {t: i for i, t in enumerate(types)}
    genus_index = {g: i for i, g in enumerate(genera)}
    table = np.zeros((len(types), len(genera)))
    for record in scored:
        table[type_index[record["type"]], genus_index[record["genus"]]] += 1
    confounding_v = cramers_v(table)

    by_type = defaultdict(list)
    for record in scored:
        by_type[record["type"]].append(record)
    type_rows = []
    for cluster, items in by_type.items():
        genus_counts = Counter(r["genus"] for r in items)
        dominant, dominant_n = genus_counts.most_common(1)[0]
        if len(items) < MIN_CELL_ENTRIES:
            continue
        type_rows.append({
            "type": cluster,
            "n_entries": len(items),
            "n_genera": len(genus_counts),
            "transposon_rate": round(float(np.mean(
                [r["outcome"] for r in items])), 4),
            "dominant_genus": dominant,
            "dominant_genus_share": round(dominant_n / len(items), 4),
            "reportable": bool(len(genus_counts) >= MIN_CELL_GENERA),
            "rate_outside_the_dominant_genus": (
                round(float(np.mean([r["outcome"] for r in items
                                     if r["genus"] != dominant])), 4)
                if len(items) > dominant_n else None)})
    type_rows.sort(key=lambda r: -r["transposon_rate"])

    by_genus = defaultdict(list)
    for record in scored:
        by_genus[record["genus"]].append(record)
    genus_rows = []
    for genus, items in by_genus.items():
        if len(items) < MIN_CELL_ENTRIES:
            continue
        type_counts = Counter(r["type"] for r in items)
        genus_rows.append({
            "genus": genus,
            "n_entries": len(items),
            "n_types": len(type_counts),
            "transposon_rate": round(float(np.mean(
                [r["outcome"] for r in items])), 4),
            "dominant_type": type_counts.most_common(1)[0][0],
            "dominant_type_share": round(
                type_counts.most_common(1)[0][1] / len(items), 4),
            "reportable": bool(len(type_counts) >= 2)})
    genus_rows.sort(key=lambda r: -r["transposon_rate"])

    # Tip x cins hucresine cokertilmis duzey: projenin her biyolojik iddia
    # icin kosulu. Ikili ozellik hucre icinde cogunluk oyuyla cokertilir
    # (stats_overview.py / operon_relations.py deseni).
    cells = defaultdict(list)
    for record in scored:
        cells[(record["type"], record["genus"])].append(record["outcome"])
    collapsed = [{"type": key[0], "genus": key[1], "n_entries": len(values),
                  "outcome": int(float(np.mean(values)) >= 0.5)}
                 for key, values in cells.items()]
    collapsed_outcome = np.array([r["outcome"] for r in collapsed],
                                 dtype=float)
    collapsed_effects = {}
    for key, label in (("genus", "genus"), ("type", "enzyme type")):
        effect = factor_effect([r[key] for r in collapsed],
                              collapsed_outcome, rng)
        if effect:
            effect["factor"] = label
            collapsed_effects[key] = effect

    density_rows = []
    for quartile in sorted({r["density_quartile"] for r in scored}):
        items = [r for r in scored if r["density_quartile"] == quartile]
        density_rows.append({
            "quartile": quartile,
            "n_entries": len(items),
            "mean_annotated_neighbours": round(float(np.mean(
                [r["n_neighbours"] for r in items])), 1),
            "transposon_rate": round(float(np.mean(
                [r["outcome"] for r in items])), 4)})

    genus_excess = effects["genus"]["eta_squared_above_chance"]
    type_excess = effects["type"]["eta_squared_above_chance"]
    verdict = (
        "At every level the enzyme type explains more of the transposon call "
        "than the genus does. Over single entries the type sits %.3f above "
        "its permuted chance level and the genus %.3f; holding the genus "
        "fixed, the type still adds %.3f, while holding the type fixed the "
        "genus adds %.3f; and after collapsing to one observation per type "
        "and genus the type is %.3f against %.3f. So the pattern is more "
        "enzyme-specific than species-specific. Two reservations are real: "
        "the two factors are partly confounded (Cramer's V %.2f between type "
        "and genus, and some types sit almost entirely in one genus), and "
        "what is being predicted is where a transposon has been ANNOTATED, "
        "not where one is." % (
            type_excess, genus_excess,
            held["type_within_genus"]["eta_squared_above_chance"]
            if held["type_within_genus"] else float("nan"),
            held["genus_within_type"]["eta_squared_above_chance"]
            if held["genus_within_type"] else float("nan"),
            collapsed_effects["type"]["eta_squared_above_chance"],
            collapsed_effects["genus"]["eta_squared_above_chance"],
            confounding_v))

    return {
        "question": ("Is transposon association a property of the species or "
                     "of the enzyme? (ROADMAP item 53)"),
        "critical_caveat_carried_from_operon_relations": {
            "what_the_call_actually_is": (
                "The transposon label is a regular expression over the "
                "GenBank /product text of the neighbouring genes "
                "(gene_category, method='%s'), not sequence evidence. Keys "
                "such as integrase, recombinase and resolvase also match "
                "proteins that are not transposases. What is measured is "
                "whether the annotation text matches the expression, not "
                "whether a mobile element is there." % REGEX_METHOD),
            "measured_size_of_the_annotation_bias": (
                "operon_relations.py stratified the substrate-class "
                "association by how many genes are annotated in the window "
                "and combined the strata with Mantel-Haenszel: the common "
                "odds ratio fell from 1.41 to 1.31, so roughly a quarter of "
                "the apparent excess is annotation density and three "
                "quarters are not. The window-occupancy rate itself runs "
                "from 11.3 % to 25.2 % across quartiles, a factor of 2.2."),
            "how_it_is_handled_here": (
                "Annotation density is entered as a third factor measured on "
                "the same scale as genus and type, and both effects are "
                "recomputed inside density quartiles. That does not repair "
                "the regular expression; it only bounds how much of the "
                "ranking could be an annotation artefact."),
            "not_presented_as_verified_biology": (
                "Every number in this block describes annotated transposon "
                "proximity. Sequence verification would need a transposase "
                "profile HMM run over the neighbour proteins and has not "
                "been done.")},
        "method": {
            "outcome": ("1 if at least one neighbour within %d bp carries the "
                        "regex category '%s', else 0"
                        % (NEIGHBOUR_WINDOW_BP, TRANSPOSON_CATEGORY)),
            "effect_size": (
                "eta squared, the share of variance of the binary outcome "
                "explained by a factor. For a k-by-2 table this is exactly "
                "the square of Cramer's V, so genus, enzyme type and "
                "annotation density are all read on one scale."),
            "why_a_permutation_null": (
                "A factor with 508 levels produces a large eta squared even "
                "with no association at all, roughly (k-1)/n. The number to "
                "read is therefore not eta squared but how far it stands "
                "above its permuted level, with %d permutations and seed %d."
                % (PERMUTATIONS, RANDOM_SEED)),
            "excluded_entries": (
                "%d of %d confirmed entries have no annotated neighbour at "
                "all and were dropped. In those records the transposon call "
                "is forced to zero, so keeping them in the denominator would "
                "dilute both factors equally and flatter neither."
                % (len(everything) - len(scored), len(everything))),
            "density_quartile_cuts": [int(c) for c in cuts]},
        "coverage": {
            "confirmed_entries": len(everything),
            "entries_scored": len(scored),
            "entries_dropped_without_annotated_neighbours":
                len(everything) - len(scored),
            "transposon_rate_all_confirmed_entries": round(
                float(np.mean([r["outcome"] for r in everything])), 4),
            "transposon_rate_scored_entries": round(
                float(outcome.mean()), 4),
            "n_types": len(types),
            "n_genera": len(genera),
            "n_species": len({r["species"] for r in scored})},
        "factor_effects_every_entry": effects,
        "factor_effects_one_factor_held_fixed": held,
        "factor_effects_collapsed_to_one_observation_per_type_and_genus": {
            "n_cells": len(collapsed),
            "rule": ("one observation per (enzyme type, genus) cell, the "
                     "binary outcome collapsed by majority vote, the same "
                     "rule stats_overview.py and operon_relations.py use"),
            "transposon_rate": round(float(collapsed_outcome.mean()), 4),
            "effects": collapsed_effects},
        "are_the_two_factors_confounded": {
            "cramers_v_type_vs_genus": round(confounding_v, 4),
            "reading": ("The two factors are not independent, so no split "
                        "between them is exact. The most extreme case is "
                        "GbcA, whose entries are almost all one genus; for "
                        "such a type the type effect and the genus effect "
                        "cannot be told apart at all."),
            "single_genus_types": sorted(
                [{"type": r["type"], "dominant_genus": r["dominant_genus"],
                  "dominant_genus_share": r["dominant_genus_share"],
                  "n_entries": r["n_entries"]}
                 for r in type_rows if r["dominant_genus_share"] >= 0.5],
                key=lambda r: -r["dominant_genus_share"])},
        "annotation_density_strata": density_rows,
        "per_type": type_rows,
        "per_genus": genus_rows,
        "verdict": verdict,
        "what_this_cannot_settle": [
            "Whether a mobile element is actually present. The call is a "
            "regular expression over annotation text; its own accuracy was "
            "not measured here and cannot be measured from this table.",
            "Whether the enzyme-over-genus ranking would survive a sequence-"
            "verified transposase call. The annotation-density control bounds "
            "one bias, not the regular expression's own error rate.",
            "Direction. Nothing here says whether transposons moved these "
            "enzymes or accumulated beside them.",
            "The genus effect is measured on the genera that were sequenced. "
            "ROADMAP item 61 records that 2,975 replicons (17.4 %) were "
            "downloaded as empty skeletons and that the loss concentrates in "
            "Streptomyces, Pseudomonas and Burkholderia, which is exactly "
            "where a genus effect would show, so the genus side of this "
            "comparison is measured on an incomplete sample."]}


# -------------------------------------- 3. SORU: BENZER ENZIM, BENZER OPERON?

def question_three(con, aligned, py_rng):
    """Enzim kimligi ile operon gen sirasi benzerligi arasindaki iliski."""
    entries = {}
    for candidate_id, cluster, organism, n_genes, layout in con.execute("""
            SELECT r.candidate_id, r.ro_cluster, p.organism,
                   o.n_genes, o.layout
            FROM ro r JOIN replicon p USING(nucleotide_id)
            JOIN operon o USING(candidate_id)
            WHERE r.is_confirmed = 1"""):
        if candidate_id not in aligned or not layout:
            continue
        if not n_genes or n_genes < MIN_OPERON_GENES:
            continue
        words = (organism or "?").split()
        entries[candidate_id] = {
            "type": cluster or "?",
            "genus": words[0] if words else "?",
            "n_genes": n_genes,
            "layout": layout.split(" > ")}
    components, categories, positions = (defaultdict(set), defaultdict(set),
                                         defaultdict(dict))
    for candidate_id, position, component, category in con.execute(
            "SELECT candidate_id, position, component, category "
            "FROM operon_gene"):
        if candidate_id not in entries:
            continue
        # Pozisyon 0 RO'nun kendisidir ve her operonda 'alpha'/'ro_alpha'
        # etiketiyle durur. Kompozisyon kumesine katilsa iki operon KOMSUSU
        # olmasa bile Jaccard 1 verirdi; olculen sey komsularin kompozisyonu
        # oldugu icin capa disarida tutulur.
        if position == 0:
            positions[candidate_id][position] = "[alpha]"
            continue
        if component and component != "none":
            components[candidate_id].add(component)
        if category:
            categories[candidate_id].add(category)
        label = component if component and component != "none" else (
            category or "other")
        positions[candidate_id][position] = label

    def sequence_only(layout):
        return [t if t in SEQUENCE_EVIDENCE_TOKENS else TEXT_TOKEN_PLACEHOLDER
                for t in layout]

    def measures(a, b):
        left, right = entries[a]["layout"], entries[b]["layout"]
        return {
            "order_lcs": order_similarity(left, right),
            "order_lcs_sequence_evidence_only": order_similarity(
                sequence_only(left), sequence_only(right)),
            "component_jaccard": jaccard(components[a], components[b]),
            "category_jaccard": jaccard(categories[a], categories[b]),
            "anchor_agreement": anchor_agreement(positions[a],
                                                 positions[b])[0]}

    all_indices = list(range(len(next(iter(aligned.values())))))
    by_type = defaultdict(list)
    for candidate_id, record in entries.items():
        by_type[record["type"]].append(candidate_id)

    pairs = []
    sampling_log = []
    for cluster in sorted(by_type):
        members = by_type[cluster]
        total = len(members) * (len(members) - 1) // 2
        selected = sample_pairs(members, MAX_PAIRS_PER_TYPE, py_rng)
        if not selected:
            continue
        sampling_log.append({"type": cluster, "n_entries": len(members),
                             "possible_pairs": total,
                             "pairs_used": len(selected),
                             "exhaustive": total <= MAX_PAIRS_PER_TYPE})
        for a, b in selected:
            identity = column_identity(aligned[a], aligned[b],
                                       all_indices)[0]
            if identity is None:
                continue
            record = {"type": cluster, "a": a, "b": b,
                      "enzyme_identity": round(identity, 4),
                      "genus_a": entries[a]["genus"],
                      "genus_b": entries[b]["genus"],
                      "same_genus": int(entries[a]["genus"]
                                        == entries[b]["genus"]),
                      "min_operon_genes": min(entries[a]["n_genes"],
                                              entries[b]["n_genes"])}
            record.update(measures(a, b))
            pairs.append(record)

    # SANS duzeyi: tipler ARASI ciftler. Zorunlu, cunku her layout
    # '[alpha]' iceriyor ve 'other' jetonu neredeyse her operonda var.
    candidates = sorted(entries)
    chance = []
    attempts = 0
    while len(chance) < CHANCE_LEVEL_PAIRS and attempts < 20 * CHANCE_LEVEL_PAIRS:
        attempts += 1
        a, b = py_rng.choice(candidates), py_rng.choice(candidates)
        if a == b or entries[a]["type"] == entries[b]["type"]:
            continue
        identity = column_identity(aligned[a], aligned[b], all_indices)[0]
        if identity is None:
            continue
        record = {"enzyme_identity": round(identity, 4)}
        record.update(measures(a, b))
        chance.append(record)

    fields = ("order_lcs", "order_lcs_sequence_evidence_only",
              "component_jaccard", "category_jaccard", "anchor_agreement")
    identities = np.array([r["enzyme_identity"] for r in pairs])

    def correlate(records, field):
        x = [r["enzyme_identity"] for r in records if r[field] is not None]
        y = [r[field] for r in records if r[field] is not None]
        rho, p = spearman(x, y)
        return {"n": len(x),
                "spearman_rho": round(rho, 4) if rho is not None else None,
                "p": p}

    measure_block = {}
    for field in fields:
        within = describe([r[field] for r in pairs])
        between = describe([r[field] for r in chance])
        per_type_rho = []
        for cluster in sorted(by_type):
            subset = [r for r in pairs if r["type"] == cluster]
            if len(subset) < MIN_PAIRS_FOR_CORRELATION:
                continue
            got = correlate(subset, field)
            if got["spearman_rho"] is not None:
                per_type_rho.append(got["spearman_rho"])
        measure_block[field] = {
            "within_type_pairs": within,
            "cross_type_chance_level": between,
            "chance_level_median": (between or {}).get("median"),
            "correlation_with_enzyme_identity": {
                "every_pair": correlate(pairs, field),
                "one_correlation_per_type_then_sign_test":
                    sign_test(per_type_rho),
                "cross_genus_pairs_only": correlate(
                    [r for r in pairs if not r["same_genus"]], field),
                "same_genus_pairs_only": correlate(
                    [r for r in pairs if r["same_genus"]], field)}}

    # Sira kompozisyonun UZERINE bir sey katiyor mu: kismi korelasyon.
    # Spearman'lar uzerinden birinci dereceden kismi korelasyon; amaci
    # "ayni genler" ile "ayni sirada" arasindaki farki ayirmak.
    def partial_spearman(records, field, control):
        usable = [r for r in records
                  if r[field] is not None and r[control] is not None]
        if len(usable) < 50:
            return None
        xy = spearman([r["enzyme_identity"] for r in usable],
                      [r[field] for r in usable])[0]
        xz = spearman([r["enzyme_identity"] for r in usable],
                      [r[control] for r in usable])[0]
        yz = spearman([r[field] for r in usable],
                      [r[control] for r in usable])[0]
        if None in (xy, xz, yz):
            return None
        denominator = ((1 - xz ** 2) * (1 - yz ** 2)) ** 0.5
        if denominator <= 0:
            return None
        return {"n": len(usable),
                "spearman_rho": round(xy, 4),
                "spearman_rho_controlling_for_" + control: round(
                    (xy - xz * yz) / denominator, 4),
                "correlation_between_the_two_measures": round(yz, 4)}

    # Kuratorun supheli durumu: protein neredeyse ayni, operon benzemiyor.
    chance_median = measure_block["order_lcs"]["chance_level_median"]
    def meaningless(records, identity_gate):
        out = [r for r in records
               if r["enzyme_identity"] >= identity_gate
               and r["order_lcs"] is not None
               and r["order_lcs"] <= chance_median]
        return out

    suspect_strict = meaningless(pairs, HIGH_ENZYME_IDENTITY)
    suspect_relaxed = meaningless(pairs, RELAXED_ENZYME_IDENTITY)
    high_pairs = [r for r in pairs
                  if r["enzyme_identity"] >= HIGH_ENZYME_IDENTITY]
    relaxed_pairs = [r for r in pairs
                     if r["enzyme_identity"] >= RELAXED_ENZYME_IDENTITY]

    def detail(records):
        return [{"type": r["type"], "a": r["a"], "b": r["b"],
                 "enzyme_identity": r["enzyme_identity"],
                 "order_lcs": round(r["order_lcs"], 4),
                 "component_jaccard": (round(r["component_jaccard"], 4)
                                       if r["component_jaccard"] is not None
                                       else None),
                 "genus_a": r["genus_a"], "genus_b": r["genus_b"],
                 "layout_a": " > ".join(entries[r["a"]]["layout"]),
                 "layout_b": " > ".join(entries[r["b"]]["layout"])}
                for r in sorted(records, key=lambda r: r["order_lcs"])[:40]]

    scatter = scatter_sample(
        [{"type": r["type"], "identity": r["enzyme_identity"],
          "order_lcs": r["order_lcs"],
          "category_jaccard": r["category_jaccard"],
          "cross_genus": 1 - r["same_genus"]} for r in pairs], py_rng)

    primary = measure_block["order_lcs"]["correlation_with_enzyme_identity"]
    component = measure_block["component_jaccard"][
        "correlation_with_enzyme_identity"]
    verdict = (
        "The curator's expectation holds. Enzyme identity and operon gene "
        "order move together: Spearman %.2f over all %d within-type pairs, "
        "and positive in %d of %d types taken one at a time (sign test "
        "p = %.1e, median rho %.2f). It survives the genus control (%.2f on "
        "cross-genus pairs). It is not an artefact of the annotation text "
        "either: the only purely sequence-based measure, the overlap of the "
        "HMM-called operon components, correlates with enzyme identity at "
        "%.2f over %d pairs, more strongly than gene order does, while the "
        "order measure built only from sequence-evidence tokens gives %.2f. "
        "The case the curator feared is real but rare: of %d pairs at %.0f %% "
        "identity or above, %d sit at or below the cross-type chance level of "
        "operon similarity. That rarity is what makes the measure useful as a "
        "flag, and it also means an operon mismatch cannot carry much weight "
        "as a general test, because almost nothing fails it." %
        (primary["every_pair"]["spearman_rho"], primary["every_pair"]["n"],
         primary["one_correlation_per_type_then_sign_test"][
             "n_above_reference"],
         primary["one_correlation_per_type_then_sign_test"]["n_types"],
         primary["one_correlation_per_type_then_sign_test"]["p_sign_test"],
         primary["one_correlation_per_type_then_sign_test"]["median"],
         primary["cross_genus_pairs_only"]["spearman_rho"],
         component["every_pair"]["spearman_rho"],
         component["every_pair"]["n"],
         measure_block["order_lcs_sequence_evidence_only"][
             "correlation_with_enzyme_identity"]["every_pair"][
                 "spearman_rho"],
         len(high_pairs), 100 * HIGH_ENZYME_IDENTITY, len(suspect_strict)))

    return {
        "question": ("The curator wrote that enzymes very similar to each "
                     "other should also have similar operons, and that a pair "
                     "whose operons are not similar may be a random match, "
                     "the same gene name carrying no meaning. (ROADMAP item "
                     "54)"),
        "method": {
            "unit": ("pairs of confirmed entries inside the same enzyme type"),
            "enzyme_identity": (
                "percent identity over the match-state columns of the motif "
                "profile, gaps excluded from the denominator. This is the "
                "same quantity the variant clustering is built on, and it was "
                "cross-checked against an independent pairwise global "
                "alignment in the first block of this file."),
            "gene_order_measure": (
                "longest common subsequence of the two layout token "
                "sequences, divided by the shorter layout's length. Chosen "
                "over an edit distance because an operon boundary here comes "
                "from a 150 bp gap rule and from how many genes the submitter "
                "annotated, so an inserted gene is mostly an annotation "
                "difference; an edit distance would charge for it while the "
                "longest common subsequence can step over it. Normalising by "
                "the shorter layout is deliberate: the extra genes of the "
                "longer operon could not have been observed in the shorter "
                "one. The measure is symmetric, order-aware and bounded by 0 "
                "and 1."),
            "why_a_chance_level_is_required": (
                "Every layout contains the anchor token [alpha] and about "
                "nine operons in ten contain the token 'other', so the "
                "longest common subsequence has a high floor. Cross-type "
                "pairs give that floor empirically; a within-type similarity "
                "is only informative against it."),
            "composition_measures": (
                "Jaccard over the HMM-called components and over the regex "
                "categories from operon_gene, so that sharing the same genes "
                "can be told apart from sharing the same order."),
            "anchor_measure": (
                "fraction of positions from -%d to +%d around the RO where "
                "both operons have a gene and the labels agree. The longest "
                "common subsequence forgives a shift, this does not."
                % (ANCHOR_WINDOW, ANCHOR_WINDOW)),
            "minimum_operon_size": (
                "both operons must have at least %d genes. A one-gene operon "
                "is the RO alone, which is 2,479 of the 11,422 entries, and a "
                "two-gene operon has only one possible order relative to the "
                "anchor, so there is no order to compare below three."
                % MIN_OPERON_GENES),
            "sampling": (
                "at most %d pairs per type, drawn uniformly without "
                "replacement with seed %d, the same cap and the same seed "
                "operon_relations.py uses; %d cross-type pairs were drawn the "
                "same way for the chance level"
                % (MAX_PAIRS_PER_TYPE, RANDOM_SEED, CHANCE_LEVEL_PAIRS)),
            "half_of_the_layout_is_annotation_text": (
                "Layout tokens are a mixture: beta, ferredoxin, reductase, "
                "alpha_other and rieske_other come from scanning the "
                "neighbour protein sequence with profile HMMs, while other, "
                "dehydrogenase, transporter, hydrolase, regulator, "
                "ring_cleavage, hypothetical and transposon come from a "
                "regular expression over the GenBank product text. The whole "
                "measurement is therefore repeated with the text tokens "
                "replaced by a single wildcard.")},
        "coverage": {
            "minimum_operon_genes": MIN_OPERON_GENES,
            "entries_with_a_long_enough_operon": len(entries),
            "confirmed_entries_in_database": con.execute(
                "SELECT COUNT(*) FROM ro WHERE is_confirmed = 1").fetchone()[0],
            "within_type_pairs_scored": len(pairs),
            "n_types": len({r["type"] for r in pairs}),
            "n_genera": len({r["genus_a"] for r in pairs}
                            | {r["genus_b"] for r in pairs}),
            "cross_genus_pairs": sum(1 for r in pairs if not r["same_genus"]),
            "same_genus_pairs": sum(1 for r in pairs if r["same_genus"]),
            "chance_level_pairs": len(chance),
            "sampling_per_type": sorted(sampling_log,
                                        key=lambda r: -r["n_entries"])},
        "measures": measure_block,
        "does_order_add_anything_beyond_composition": {
            "question": ("is the relation between enzyme identity and operon "
                         "order just a restatement of the two operons sharing "
                         "the same genes"),
            "order_controlling_for_component_composition": partial_spearman(
                pairs, "order_lcs", "component_jaccard"),
            "order_controlling_for_category_composition": partial_spearman(
                pairs, "order_lcs", "category_jaccard"),
            "reading": ("A first-order partial correlation, so it assumes a "
                        "monotone relation and nothing more; it says whether "
                        "order carries any signal once shared gene content is "
                        "held fixed, not how much.")},
        "which_measure_to_use": {
            "recommended": "order_lcs",
            "why": ("It is defined for every pair, it is order-aware, it is "
                    "symmetric, it is bounded, and it does not charge for an "
                    "extra annotated gene. Its weakness is a high floor, "
                    "which is why the cross-type chance level is reported "
                    "beside every value and why the flag for a suspicious "
                    "pair is set at that floor rather than at a round "
                    "number."),
            "what_the_other_measures_add": (
                "component_jaccard is the only measure built purely on "
                "sequence evidence, the beta, ferredoxin, reductase and "
                "related subunits called by profile HMMs on the neighbour "
                "proteins, and it correlates with enzyme identity more "
                "strongly than gene order does. That is the most important "
                "line in this block: the relation does not depend on the "
                "annotation text. Its limits are that it is undefined for "
                "the pairs where neither operon has any called component, "
                "and that it ignores order entirely. anchor_agreement also "
                "correlates more strongly than order_lcs, because it refuses "
                "to forgive a shift, but it too is undefined when the two "
                "operons have no overlapping positions. category_jaccard is "
                "literally 'do the two operons carry the same gene names', "
                "the thing the curator distrusted, and it correlates with "
                "enzyme identity about as strongly as gene order does, which "
                "is a reminder that part of the pooled relation lives in the "
                "annotation text. order_lcs is kept as the headline measure "
                "because it is the only one that is defined for every pair "
                "and still sensitive to order, not because it is the "
                "strongest."),
            "why_not_the_strongest_measure": (
                "A measure that is undefined for part of the data cannot be "
                "the number printed next to a pair, because the pairs where "
                "it is undefined are exactly the thinly annotated ones and "
                "dropping them silently would flatter every summary. The "
                "stronger measures are reported in full beside it."),
            "floor_of_the_order_measure": (
                "Both layouts always contain the anchor token, so the "
                "smallest possible value is 1 divided by the shorter "
                "layout's length, that is 0.33 for the shortest operons "
                "accepted here.")},
        "suspected_meaningless_matches": {
            "question": ("which pairs are nearly identical as proteins yet "
                         "share no operon structure"),
            "definition": (
                "enzyme identity at least %.2f and operon order similarity at "
                "or below %.3f, the median of the cross-type chance level, "
                "which means the two operons resemble each other no more than "
                "two operons from different enzyme types do"
                % (HIGH_ENZYME_IDENTITY, chance_median)),
            "n_pairs_at_or_above_the_identity_gate": len(high_pairs),
            "n_pairs_flagged": len(suspect_strict),
            "share_flagged": (round(len(suspect_strict) / len(high_pairs), 5)
                              if high_pairs else None),
            "types_involved": sorted({r["type"] for r in suspect_strict}),
            "pairs": detail(suspect_strict),
            "relaxed_gate": {
                "definition": ("the same with the identity gate at %.2f"
                               % RELAXED_ENZYME_IDENTITY),
                "n_pairs_at_or_above_the_identity_gate": len(relaxed_pairs),
                "n_pairs_flagged": len(suspect_relaxed),
                "types_involved": sorted({r["type"]
                                          for r in suspect_relaxed}),
                "pairs": detail(suspect_relaxed)},
            "reading": (
                "The curator's suspicion is correct in kind and small in "
                "size. Almost every near-identical protein pair also has a "
                "recognisable operon in common, so an operon that does not "
                "match is unusual enough to be worth looking at one by one. "
                "What the flag does not say is which of the two explanations "
                "applies: a genuinely different genomic context, or a "
                "neighbour that nobody annotated.")},
        "scatter": {"note": ("a seeded sample of at most %d within-type pairs "
                             "for plotting" % MAX_SCATTER_POINTS),
                    "seed": RANDOM_SEED,
                    "points": scatter},
        "verdict": verdict,
        "what_this_cannot_settle": [
            "Whether a low operon similarity means a wrong enzyme assignment. "
            "It can equally mean an incompletely annotated window, a contig "
            "that ends inside the operon, or a real difference in genomic "
            "context.",
            "Half of the layout vocabulary is annotation text, so the "
            "headline order measure inherits part of the weakness described "
            "for the transposon call. Replacing the text tokens with a "
            "wildcard lowers the correlation, so those tokens do carry some "
            "of the signal; the purely sequence-based component overlap "
            "carries more of it, which is why the conclusion survives, but "
            "the exact size of the pooled correlation does depend on "
            "annotation quality.",
            "Pairs inside a type are not independent and neither are entries "
            "from the same genus. The per-type sign test and the cross-genus "
            "subset are the controls, and they are weaker than the pooled "
            "number, which is the expected direction.",
            "The entries without an operon of at least %d genes are absent, "
            "and they are not a random subset: a short operon usually means a "
            "short contig or a thin annotation, so the measurement is made on "
            "the better-annotated part of the database." % MIN_OPERON_GENES,
            "Operon direction and gene strand are already folded into the "
            "layout by build_operons.py, so this measure cannot see a "
            "rearrangement that preserves the component order."]}


# ---------------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default="analysis_out")
    parser.add_argument("--alignment", default="genomic_context/cand_aln.sto",
                        help="uye hizalamasi (match-state kolonlari)")
    parser.add_argument("--active-site",
                        default="analysis_out/active_site.json",
                        help="ceb kolonlarinin kaynagi")
    parser.add_argument("--only", choices=["1", "2", "3"], action="append",
                        help="yalnizca bu sorulari kosur (tekrarlanabilir)")
    args = parser.parse_args()

    if not os.path.exists(args.active_site):
        print("[hata] %s bulunamadi; once active_site.py kosturulmali."
              % args.active_site, file=sys.stderr)
        return 2

    wanted = set(args.only or ["1", "2", "3"])
    connection = sqlite3.connect(args.db)
    # Her sorunun KENDI tohumu var. Ortak bir uretec kullanilsa, --only ile
    # tek soru kosturuldugunda sayilar degisirdi (ureticinin durumu bir
    # onceki sorunun kac cekilis yaptigina bagli olurdu) ve ciktinin yeniden
    # uretilebilirligi hangi sorularin kosturuldugune baglanirdi.
    seeds = {"1": RANDOM_SEED + 1, "2": RANDOM_SEED + 2, "3": RANDOM_SEED + 3}

    aligned = {}
    if {"1", "3"} & wanted:
        print("[okunuyor] uye hizalamasi")
        aligned = read_stockholm_matchcols(args.alignment)
        print("[bilgi] hizalamada %d dizi, %d match kolonu"
              % (len(aligned), len(next(iter(aligned.values())))))

    payload = {
        "module": "variant_and_operon.py",
        "generated_utc": datetime.datetime.utcnow().strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "inputs": {"database": os.path.basename(args.db),
                   "member_alignment": args.alignment,
                   "pocket_columns_from": args.active_site},
        "what_was_measured": (
            "Three questions from the project owner. First, variants are "
            "compared at the alignment columns that line the catalytic iron "
            "rather than over the whole protein, because a whole-protein "
            "identity is a weak measure in a family where curated pairs with "
            "different substrates reach 99.8 % identity, and because the "
            "number previously shown next to a variant was a within-variant "
            "quantity rather than a between-variant one. Second, the variance "
            "of the annotated-transposon call is decomposed between the genus "
            "and the enzyme type, each against a permutation null, with "
            "annotation density entered as a third factor on the same scale. "
            "Third, enzyme identity is related to operon gene-order "
            "similarity for pairs inside a type, to test whether similar "
            "enzymes really do sit in similar operons."),
        "parameters": {
            "pocket_radius_A": POCKET_RADIUS_A,
            "pocket_radius_sensitivity_A": POCKET_RADIUS_SENSITIVITY_A,
            "minimum_members_per_variant": MIN_LEAF_MEMBERS,
            "minimum_consensus_occupancy": MIN_CONSENSUS_OCCUPANCY,
            "minimum_pocket_columns_compared": MIN_POCKET_COLUMNS_COMPARED,
            "high_global_identity_gate": HIGH_GLOBAL_IDENTITY,
            "relaxed_global_identity_gate": RELAXED_GLOBAL_IDENTITY,
            "pocket_divergent_max_identity": POCKET_DIVERGENT_MAX_IDENTITY,
            "minimum_pocket_mismatches": MIN_POCKET_MISMATCHES,
            "low_global_identity_gate": LOW_GLOBAL_IDENTITY,
            "minimum_annotated_neighbours": MIN_ANNOTATED_NEIGHBOURS,
            "minimum_stratum_entries": MIN_STRATUM_ENTRIES,
            "minimum_cell_entries": MIN_CELL_ENTRIES,
            "minimum_cell_genera": MIN_CELL_GENERA,
            "minimum_operon_genes": MIN_OPERON_GENES,
            "max_pairs_per_type": MAX_PAIRS_PER_TYPE,
            "chance_level_pairs": CHANCE_LEVEL_PAIRS,
            "high_enzyme_identity_gate": HIGH_ENZYME_IDENTITY,
            "global_identity_bands": [b[2] for b in GLOBAL_IDENTITY_BANDS],
            "descriptive_sub_bands": [b[2]
                                      for b in GLOBAL_IDENTITY_SUB_BANDS],
            "minimum_expected_pocket_mismatches_for_a_usable_band":
                MIN_EXPECTED_POCKET_MISMATCHES,
            "maximum_saturated_share_for_a_usable_band":
                MAX_SATURATED_SHARE_FOR_A_USABLE_BAND,
            "saturation_sample_pairs": SATURATION_SAMPLE_PAIRS,
            "saturation_floor_percentile": SATURATION_PERCENTILE,
            "global_identity_band_seed": RANDOM_SEED + BAND_SEED_OFFSET,
            "anchor_window": ANCHOR_WINDOW,
            "permutations": PERMUTATIONS,
            "random_column_sets": RANDOM_COLUMN_SETS,
            "minimum_null_column_occupancy": MIN_NULL_COLUMN_OCCUPANCY,
            "relaxed_enzyme_identity_gate": RELAXED_ENZYME_IDENTITY,
            "random_seed": RANDOM_SEED,
            "per_question_seeds": {"question_1": RANDOM_SEED + 1,
                                   "question_2": RANDOM_SEED + 2,
                                   "question_3": RANDOM_SEED + 3},
            "neighbour_window_bp": NEIGHBOUR_WINDOW_BP,
            "regex_method": REGEX_METHOD,
            "transposon_category": TRANSPOSON_CATEGORY},
        "non_independence_rule": (
            "Entries in this database are not independent: many are "
            "near-identical strains. Every claim here is therefore reported "
            "at more than one level. For the transposon question the second "
            "level is one observation per enzyme type and genus, the rule the "
            "rest of the pipeline uses. For the two pair-based questions a "
            "collapse to one observation per pair is meaningless, so the "
            "second level is one correlation per type followed by a sign "
            "test, and the third is the subset of pairs whose two sides come "
            "from different genera. Where the levels disagree, the pooled "
            "number is the one to distrust."),
    }

    if "1" in wanted:
        print("[1. soru] ceb kimligi hesaplaniyor")
        payload["question_1_pocket_identity_between_variants"] = question_one(
            connection, aligned, args.active_site,
            random.Random(seeds["1"]))
    if "2" in wanted:
        print("[2. soru] transpozon varyans ayristirmasi")
        payload["question_2_transposon_species_or_enzyme"] = question_two(
            connection, np.random.default_rng(seeds["2"]))
    if "3" in wanted:
        print("[3. soru] operon benzerligi")
        payload["question_3_operon_similarity_vs_enzyme_identity"] = \
            question_three(connection, aligned,
                           random.Random(seeds["3"]))

    payload["limitations"] = [
        "None of the three blocks measures activity. A pocket difference, a "
        "transposon call and an operon mismatch are all descriptions of "
        "annotation and sequence, and the database holds no measured reaction "
        "for any of the 11,422 entries.",
        "The transposon call is a regular expression over GenBank product "
        "text. operon_relations.py measured that roughly a quarter of its "
        "apparent excess is annotation density, and the regular expression's "
        "own error rate was never measured. Nothing in this file should be "
        "read as verified presence of a mobile element.",
        "Pocket columns come from one structure per type and are read on "
        "every member of that type through a profile alignment. For 46 of 63 "
        "types that structure is a prediction, and side chains are the least "
        "reliable part of a prediction.",
        "The curated reference set has known errors that are recorded but not "
        "yet applied (ROADMAP items 43 and 60), including one duplicated "
        "enzyme and two wrong PDB identifiers. Type-level numbers inherit "
        "them.",
        "17.4 % of replicons were downloaded as annotation-free skeletons "
        "(ROADMAP item 61) and the loss concentrates in the genera richest in "
        "these enzymes, so every count and every genus effect here is "
        "measured on an incomplete sample.",
        "With 11,422 entries almost any comparison reaches a small p value, "
        "so no p value is reported on its own anywhere in this file; each one "
        "stands next to an effect size and, where a factor has many levels, "
        "next to a permutation null."]

    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, "variant_and_operon.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)

    print_summary(payload)
    print("\n[yazildi] %s" % path)
    connection.close()
    return 0


def print_summary(payload):
    """Konsol ozeti. Her blok icin etki buyuklugu + duzeyler."""
    line = "=" * 74
    one = payload.get("question_1_pocket_identity_between_variants")
    if one:
        print("\n" + line)
        print("1. SORU -- VARYANTLARI CEPTE KARSILASTIRMA")
        print(line)
        cov = one["coverage"]
        print("ceb kolonu olan tip: %d (kristal %d, ongoru %d)"
              % (cov["types_with_pocket_columns"],
                 cov["types_with_pocket_from_crystal"],
                 cov["types_with_pocket_from_predicted_model"]))
        print("olculen: %d tip, %d varyant, %d varyant cifti"
              % (cov["types_measured_here"], cov["variants_measured"],
                 cov["variant_pairs_measured"]))
        validation = one["profile_identity_validation"]
        if validation.get("status") == "measured":
            print("kimlik olcutu dogrulamasi: %d ciftte ikili hizalamayla "
                  "rho=%.3f, ortalama fark %+.3f"
                  % (validation["n_pairs"], validation["spearman_rho"],
                     validation["mean_difference_profile_minus_pairwise"]))
        rel = one["relation_between_global_and_pocket_identity"]
        pooled = rel["levels"][0]
        per_type = rel["levels"][1]["detail"]
        cross = rel["levels"][2]
        print("\nglobal x ceb kimligi:")
        print("  tum ciftler          rho=%.3f  (n=%d, %d tip)"
              % (pooled["spearman_rho"], pooled["n"], pooled["n_types"]))
        print("  tip basina           medyan rho=%.3f, %d/%d tipte pozitif, "
              "isaret testi p=%.2g"
              % (per_type["median"], per_type["n_above_reference"],
                 per_type["n_types"], per_type["p_sign_test"]))
        print("  cins carprazi        rho=%.3f  (n=%d, %d cins)"
              % (cross["spearman_rho"], cross["n"], cross["n_genera"]))
        cons = one["is_the_pocket_more_conserved_than_the_protein"]
        ratio = cons["per_type_ratio"]
        free = cons["per_type_ratio_excluding_metal_ligand_columns"]
        test = cons["sign_test_against_one"]
        free_test = cons["sign_test_against_one_excluding_metal_ligand_columns"]
        print("\nceb protein ortalamasindan daha mi korunmus:")
        print("  uyusmazlik orani     medyan %.3f (%.3f-%.3f), "
              "%d/%d tipte 1'in altinda, p=%.2g"
              % (ratio["median"], ratio["min"], ratio["max"],
                 test["n_above_reference"], test["n_types"],
                 test["p_sign_test"]))
        if free and free_test:
            print("  metal ligandlari HARIC  medyan %.3f (%.3f-%.3f), "
                  "%d/%d tipte 1'in altinda, p=%.2g   <-- asil sayi"
                  % (free["median"], free["min"], free["max"],
                     free_test["n_above_reference"], free_test["n_types"],
                     free_test["p_sign_test"]))
        null = cons["same_size_random_column_set_null"]
        if null["pocket_percentile"]:
            print("  ayni buyuklukte rastgele kolon kumesine karsi yuzdelik: "
                  "medyan %.3f (ligandlar haric %.3f)"
                  % (null["pocket_percentile"]["median"],
                     (null["pocket_percentile_excluding_metal_ligands"]
                      or {}).get("median", float("nan"))))
        budget = cons["invariant_column_budget"]
        print("  ceb kolonu toplami %d, varyantlar arasinda degismez %d "
              "(bunun %d'si metal ligandi)"
              % (budget["total_pocket_columns_across_types"],
                 budget["invariant_across_the_type_variants"],
                 budget["metal_ligand_columns_across_types"]))
        print("  >=5 varyantli %d tipte degismez kolon payi %.1f%%"
              % (budget["types_with_at_least_five_variants"]["n_types"],
                 100 * budget["types_with_at_least_five_variants"][
                     "invariant_share"]))
        print("  ceb birebir ayni olan cift orani: %.1f%%"
              % (100 * cons["share_of_pairs_with_an_identical_pocket"]))
        bands = cons.get("by_global_identity_band")
        if bands:
            sat = bands.get("measured_saturation")
            if sat:
                print("\n  OLCULEN DOYGUNLUK (farkli tipten %d cift): global "
                      "kimlik %.3f, ceb kimligi %.3f, doygunlukta oran %.2f; "
                      "%d. yuzdelik taban %.3f"
                      % (sat["n_cross_type_pairs"],
                         sat["global_identity"]["mean"],
                         sat["pocket_identity_excluding_metal_ligands"][
                             "mean"],
                         sat["divergence_ratio_at_saturation"],
                         sat["floor_percentile"],
                         sat["measured_global_identity_floor"]))
            print("\n  GLOBAL KIMLIK BANTLARI (metal ligandlari haric):")
            print("    %-14s %6s %5s %6s %7s %7s %7s %9s %9s %s"
                  % ("bant", "cift", "tip", "oran", "cebUyus", "globUyus",
                     "beklUyus", "p", "q", "durum"))
            rows = list(bands.get("bands") or [])
            reference = bands.get("reference_row_unstratified")
            if reference:
                rows = rows + [reference]
            for row in rows:
                ratio = row.get("per_type_divergence_ratio") or {}
                test = row.get("sign_test_against_one") or {}
                flag = ("kullanilabilir" if row.get("usable")
                        else "OLCMUYOR: " + "; ".join(
                            r.split(":")[0] for r in
                            (row.get("why_not_usable") or [])))
                print("    %-14s %6d %2d/%-2d %6.3f %7.3f %7.3f %7.2f "
                      "%9.2g %9.2g %s"
                      % (row["band"], row["n_pairs"],
                         test.get("n_above_reference", 0),
                         row["n_types_with_a_ratio"],
                         ratio.get("median", float("nan")),
                         (row["mean_pocket_divergence_excluding_metal_"
                              "ligands"] if row["mean_pocket_divergence_"
                              "excluding_metal_ligands"] is not None
                          else float("nan")),
                         row["mean_global_divergence"],
                         row["expected_pocket_mismatches_per_pair_at_the_"
                             "background_rate"],
                         test.get("p_sign_test", float("nan")),
                         (row.get("benjamini_hochberg_q")
                          if row.get("benjamini_hochberg_q") is not None
                          else float("nan")),
                         flag))
            for row in bands.get("bands") or []:
                null = row.get("same_size_random_column_set_null")
                if null:
                    print("      %-14s rastgele kolon kumesi orani %.2f, "
                          "cebin null icindeki yuzdeligi %.3f"
                          % (row["band"], null["median_random_set_ratio"],
                             null["median_pocket_percentile_in_null"]))
            paired = bands.get("within_type_paired_check")
            if paired and paired.get("sign_test"):
                print("    TIP ICINDE eslesmis: %d tipte pencere orani %.3f, "
                      "doygun bolge orani %.3f, %d/%d tipte pencere daha "
                      "dusuk, p=%.2g"
                      % (paired["n_types_with_pairs_in_both"],
                         paired["median_ratio_in_the_window"],
                         paired["median_ratio_in_the_saturated_region"],
                         paired["sign_test"]["n_above_reference"],
                         paired["sign_test"]["n_types"],
                         paired["sign_test"]["p_sign_test"]))
            print("    BH sonrasi ayakta kalan bant: %d/%d"
                  % (bands["n_bands_surviving_bh"], bands["n_bands_tested"]))
        disc = one["discordant_pairs"]
        gates = disc["n_pairs_at_each_global_gate"]
        print("\nkuratorun iki ilginc durumu:")
        print("  global kimlik >=0,95 / >=0,85 / >=0,70 cift sayisi: "
              "%d / %d / %d" % (gates["at_least_0.95"], gates["at_least_0.85"],
                                gates["at_least_0.70"]))
        for key in ("globally_similar_pocket_divergent_strict",
                    "globally_similar_pocket_divergent_relaxed",
                    "globally_divergent_pocket_identical_strict",
                    "globally_divergent_pocket_identical_relaxed"):
            block = disc[key]
            print("  %-48s %4d cift, %2d tip, %2d cins"
                  % (key, block["n_pairs"], block["n_types"],
                     block["n_genera"]))
        for row in disc["globally_similar_pocket_divergent_relaxed"][
                "pairs"][:6]:
            print("    %-16s %s / %s  global %.3f, ceb %.3f (%d/%d fark)"
                  % (row["type"], row["variant_a"], row["variant_b"],
                     row["global_identity"], row["pocket_identity"],
                     row["pocket_mismatches"], row["pocket_columns_compared"]))
        uninformative = [t for t, v in one["per_type"].items()
                         if v["n_variable_pocket_columns"] <= 2
                         or v["n_variants_with_consensus"] < 5]
        print("\nolcumun bilgi tasimadigi ya da zayif tasidigi tip: %d/%d"
              % (len(uninformative), len(one["per_type"])))
        for name in sorted(uninformative):
            v = one["per_type"][name]
            print("  %-16s %d varyant, %d/%d kolon oynak"
                  % (name, v["n_variants_with_consensus"],
                     v["n_variable_pocket_columns"], v["n_pocket_columns"]))

    two = payload.get("question_2_transposon_species_or_enzyme")
    if two:
        print("\n" + line)
        print("2. SORU -- TRANSPOZON: TUR OZELINDE MI ENZIM OZELINDE MI")
        print(line)
        cov = two["coverage"]
        print("olculen: %d/%d giris (%d girisin hic anote komsusu yok), "
              "%d tip, %d cins"
              % (cov["entries_scored"], cov["confirmed_entries"],
                 cov["entries_dropped_without_annotated_neighbours"],
                 cov["n_types"], cov["n_genera"]))
        print("regex transpozon orani: tum girisler %.1f%%, olculenler %.1f%%"
              % (100 * cov["transposon_rate_all_confirmed_entries"],
                 100 * cov["transposon_rate_scored_entries"]))
        print("\naciklanan varyans (eta-kare), permutasyon null'una karsi:")
        print("  %-34s %8s %8s %8s %6s"
              % ("faktor", "eta2", "null", "sans ustu", "k"))
        for key in ("genus", "species", "type", "group", "density_quartile",
                    "type_x_genus_cell"):
            e = two["factor_effects_every_entry"].get(key)
            if e:
                print("  %-34s %8.4f %8.4f %8.4f %6d"
                      % (e["factor"], e["eta_squared"],
                         e["permuted_eta_squared_mean"],
                         e["eta_squared_above_chance"], e["n_levels"]))
        print("\nbir faktor sabit tutuldugunda:")
        for key, block in two["factor_effects_one_factor_held_fixed"].items():
            if block:
                print("  %-34s %8.4f %8.4f %8.4f  (tabaka %d, n=%d)"
                      % (key, block["pooled_eta_squared"],
                         block["permuted_mean"],
                         block["eta_squared_above_chance"],
                         block["n_strata_used"], block["n_entries_used"]))
        collapsed = two[
            "factor_effects_collapsed_to_one_observation_per_type_and_genus"]
        print("\ntip x cins hucresine cokertilmis (%d hucre):"
              % collapsed["n_cells"])
        for key, e in collapsed["effects"].items():
            print("  %-34s %8.4f %8.4f %8.4f %6d"
                  % (e["factor"], e["eta_squared"],
                     e["permuted_eta_squared_mean"],
                     e["eta_squared_above_chance"], e["n_levels"]))
        print("\nkarisiklik: tip x cins Cramer's V = %.3f"
              % two["are_the_two_factors_confounded"]["cramers_v_type_vs_genus"])
        single = two["are_the_two_factors_confounded"]["single_genus_types"]
        for row in single[:5]:
            print("  %-16s %%%.0f %s (n=%d)"
                  % (row["type"], 100 * row["dominant_genus_share"],
                     row["dominant_genus"], row["n_entries"]))
        print("\nanotasyon yogunlugu ceyreklikleri:")
        for row in two["annotation_density_strata"]:
            print("  %s n=%5d ortalama komsu %.1f oran %.1f%%"
                  % (row["quartile"], row["n_entries"],
                     row["mean_annotated_neighbours"],
                     100 * row["transposon_rate"]))
        print("\nen yuksek ve en dusuk oranli tipler (>=%d giris, >=%d cins):"
              % (MIN_CELL_ENTRIES, MIN_CELL_GENERA))
        reportable = [r for r in two["per_type"] if r["reportable"]]
        for row in reportable[:5] + reportable[-5:]:
            print("  %-16s n=%5d oran %.1f%%  baskin cins %s (%%%.0f), "
                  "o cins disinda %.1f%%"
                  % (row["type"], row["n_entries"],
                     100 * row["transposon_rate"], row["dominant_genus"],
                     100 * row["dominant_genus_share"],
                     100 * (row["rate_outside_the_dominant_genus"] or 0.0)))
        print("\nKARAR: " + two["verdict"])

    three = payload.get("question_3_operon_similarity_vs_enzyme_identity")
    if three:
        print("\n" + line)
        print("3. SORU -- BENZER ENZIMLERIN OPERONLARI DA BENZER MI")
        print(line)
        cov = three["coverage"]
        print("olculen: %d giris (>=%d genli operon, %d dogrulanmis giristen), "
              "%d tip ici cift, %d tip, %d cins"
              % (cov["entries_with_a_long_enough_operon"], MIN_OPERON_GENES,
                 cov["confirmed_entries_in_database"],
                 cov["within_type_pairs_scored"], cov["n_types"],
                 cov["n_genera"]))
        print("  cins carprazi cift %d, ayni cins %d, sans duzeyi cifti %d"
              % (cov["cross_genus_pairs"], cov["same_genus_pairs"],
                 cov["chance_level_pairs"]))
        print("\n%-38s %7s %7s %7s %7s %7s"
              % ("olcut", "tipici", "sans", "rho", "medrho", "+/tip"))
        for field, block in three["measures"].items():
            corr = block["correlation_with_enzyme_identity"]
            test = corr["one_correlation_per_type_then_sign_test"]
            print("%-38s %7.3f %7.3f %7.3f %7.3f %3d/%-3d"
                  % (field, block["within_type_pairs"]["median"],
                     block["cross_type_chance_level"]["median"],
                     corr["every_pair"]["spearman_rho"],
                     test["median"], test["n_above_reference"],
                     test["n_types"]))
        primary = three["measures"]["order_lcs"][
            "correlation_with_enzyme_identity"]
        print("\ngen sirasi (order_lcs) duzeyler:")
        print("  tum ciftler   rho=%.3f (n=%d)"
              % (primary["every_pair"]["spearman_rho"],
                 primary["every_pair"]["n"]))
        test = primary["one_correlation_per_type_then_sign_test"]
        print("  tip basina    medyan rho=%.3f, %d/%d pozitif, p=%.2g"
              % (test["median"], test["n_above_reference"], test["n_types"],
                 test["p_sign_test"]))
        print("  cins carprazi rho=%.3f (n=%d) | ayni cins rho=%.3f (n=%d)"
              % (primary["cross_genus_pairs_only"]["spearman_rho"],
                 primary["cross_genus_pairs_only"]["n"],
                 primary["same_genus_pairs_only"]["spearman_rho"],
                 primary["same_genus_pairs_only"]["n"]))
        partial = three["does_order_add_anything_beyond_composition"][
            "order_controlling_for_component_composition"]
        if partial:
            print("  kompozisyon sabit tutulunca rho=%.3f -> %.3f"
                  % (partial["spearman_rho"],
                     partial["spearman_rho_controlling_for_"
                             "component_jaccard"]))
        suspect = three["suspected_meaningless_matches"]
        print("\nkuratorun suphelendigi durum (protein ~ayni, operon benzemiyor):")
        print("  >=%.0f%% kimlikli cift: %d, bunlardan sans duzeyinde ya da "
              "altinda: %d (%.2f%%)"
              % (100 * HIGH_ENZYME_IDENTITY,
                 suspect["n_pairs_at_or_above_the_identity_gate"],
                 suspect["n_pairs_flagged"],
                 100 * (suspect["share_flagged"] or 0.0)))
        print("  ilgili tipler: %s" % (", ".join(suspect["types_involved"])
                                       or "-"))
        relaxed = suspect["relaxed_gate"]
        print("  gevsek kapi (>=%.0f%%): %d ciftten %d isaretlendi, tipler: %s"
              % (100 * RELAXED_ENZYME_IDENTITY,
                 relaxed["n_pairs_at_or_above_the_identity_gate"],
                 relaxed["n_pairs_flagged"],
                 ", ".join(relaxed["types_involved"]) or "-"))
        for row in suspect["pairs"][:6]:
            print("    %-16s %s / %s  kimlik %.3f, sira %.3f"
                  % (row["type"], row["genus_a"], row["genus_b"],
                     row["enzyme_identity"], row["order_lcs"]))
        print("\nKARAR: " + three["verdict"])


if __name__ == "__main__":
    sys.exit(main())
