# predictions/ -- TAHMIN katmani

Bu dizin, **dizi analiziyle bulunmus ama deneysel olarak karakterize
EDILMEMIS** Rieske oksijenaz adaylarini tutar.

## Neden ayri bir katman

`types/` karakterize edilmis enzimlerin seti. Bir referansin oradaki anlami
sudur: biri bu proteini saflastirdi, substratini olctu, urununu belirledi.
Deneysel kaniti olmayan bir diziyi oraya koymak o anlami bozar -- ve sette
19.000 girisin siniflandirmasi o referanslara dayandigi icin, bozulma tek bir
satirda kalmaz.

O yuzden tahminler AYRI durur ve asagidaki kurallar KATIDIR.

## Kurallar

Bu katmandaki hicbir sey:

- `types/` dizinine girmez
- `chemistry.csv`'ye girmez
- HMM referans kutuphanesine (`RieskeDB71.hmm`, `ROmotif71.hmm`) girmez
- hicbir istatistige, orana, q-degerine veya sekile girmez
- bir `ro_group`a atanmaz

Sitede her zaman "tahmin, karakterize edilmemis" olarak etiketlenir.

Bir tahmin `types/`'a ancak **deneysel kanit yayimlandiginda** tasinir: o
noktada normal `add_type.py` kapilarindan gecer ve burada degil orada yasar.

## Kanitin okunmasi

Her tahmin icin uc sayi birlikte raporlanir:

| alan | ne demek |
|---|---|
| Rieske 4/4 | dort Rieske ligandi modelin TURETILMIS kolonunda mi |
| katalitik 3/3 | mononukleer Fe 2-His-1-karboksilat kolonunda mi |
| kapsama | dizinin modelin kac kolonuna hizalandigi |

**Kapsama olmadan motif skoru okunmaz.** Dusuk kapsama NEGATIF DEGILDIR:
hizalama katalitik kolonlara ulasmamissa "triad yok" degil "olculemedi"
denir. Bu kural AerC'de olculdu -- 717 aa'lik protein modele %18.9 kapsamayla
hizalandi, son hizali kolon 231 oldu ve katalitik kolonlara (295/300/479) hic
ulasilmadi; sonuc 0/3 cikti ve YANLISTI. Parca parca hizalandiginda ligandlar
yerindeydi. Ayrintisi: `analysis_out/aerc_resolution.json`.

Ayni ders karboksilat denetiminde de olculmustu: "eksik" denen 1.414 girisin
1.388'i hizalama artefaktiydi.

## Testin kalibrasyonu

Bu katmani uretirken kullanilan test, ayni kosuda pozitif ve negatif
kontrollerle sinanmistir:

| kontrol | n | 4/4 Rieske | 3/3 katalitik | gecen |
|---|---|---|---|---|
| 75 kuratorlu referans (pozitif) | 75 | 75 | 73 | **75** |
| sitokrom b6-f kompleks ISP (negatif) | 55 | **55** | **0** | **0** |
| Rieske ferredoksin (negatif) | 4 | 0 | 0 | 0 |
| NirD (negatif) | 7 | 1 | 0 | 0 |

Negatif kontroller onemli: 66 gercek non-RO'nun 56'sinda dort Rieske ligandi
VAR, hicbirinde katalitik triad yok, hicbiri gecmiyor. Ayrimin biyokimyasal
temeli tam olarak bu, ve test onu yakaliyor. Yanlis pozitif 0/66.

## Regex TEK BASINA kanit degil

Ham dizide desen aramak zayif kanittir ve bu olculdu. Diziyi 10.000 kez
karistirinca mononukleer deseni `[DE]-x(2)-H-x(3,5)-H` hala uyuyor:

- ScAerC (717 aa): **%16.6**
- NsSuoC (720 aa): %11.0
- MsmA (414 aa): **%27.8**

Rieske deseni daha iyi (p~0.004) ama o da tek basina kesin degil. Agirligi
KOLON-TABANLI test tasir. Bir tahmin yalnizca regex'e dayaniyorsa, bu
kabul edilmez.

## Icerik

| dosya | ne |
|---|---|
| `archaeal.tsv` | 109 arkeal aday, tam katalitik donanimli |
| `fungal.tsv` | 654 fungal aday, tam katalitik donanimli |
| `../analysis_out/archaeal_fungal_ro_prediction.json` | yontem, kontroller, cekinceler |
| `../analysis_out/archaeal_clade/` | filogeni (clade.nwk), HGT testi, kosu betigi |

## Arkeal bulgunun ozeti

109 aday, **331-505 aa** (medyan 381), kopru Asp 109/109. Uc filum:
Halobacteria 92, Thermoproteota 15, Thermoplasmata 2. Hepsi "Predicted" ya da
"Inferred from homology" -- deneysel kanit SIFIR.

En yakin kuratorlu referanslar COGUNLUKLA GRUP 5: 61/109, basta CntA (59).
Grup 3 ikinci (44). Filum kirilimi: Halobacteria grup 5'e (61) ve grup 3'e
(31); Thermoproteota agirlikli grup 3'e (11/15).

22.531 bakteriyel Rieske proteinine karsi en yakin akraba: medyan %33.4,
maksimum %51.7. %70 ustu yok. Yani YAKIN ZAMANLI HGT DEGIL.

Ama TEK BIR SOY da degiller: agacta 7 saf kume var, her biri kendi icinde
guclu destekli (0,985-1,0) ve her birinin kardesi bakteriyel. Yalniz
Halobacteria'da en az iki bagimsiz edinim gorunuyor. Derin dugum destekleri
0,52-0,97 arasi oldugu icin VERICI GRUPLAR IDDIA EDILEMEZ.

## Fungal bulgunun ozeti

654 aday, kimlik medyan %26.9, maksimum %38.5. Arkeal kumeden FARKLI bir yere
dusuyorlar: en yakin referanslar CntA 493, CmoB 88, CmoS 50 -- yani grup 5,
kuaterner aminler. 16 tanesi Nvd/DAF-36 sterol 7-desaturazlarina dusuyor.

Not: fungal "Choline monooxygenase, chloroplastic" anotasyonlari bitkilerden
transfer artefaktidir. Bu test anotasyona DEGIL kolona baktigi icin bulgu o
artefaktten bagimsizdir.

**DIKKAT.** Arkeal kume de cogunlukla grup 5'e dusuyor (61/109). Yani "arkea
grup 3'e, fungi grup 5'e dusuyor, bu ayrisma ortak bir artefakta karsi kanit"
bicimindeki arguman GECERSIZDIR. Ortak artefakt olasiligi bu olcumle
dislanmamistir; ayri bir kontrol gerekir.
