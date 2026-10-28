# Varyant kimlikleri: icerikten turetilmis, atamadan bagimsiz

## Sorun

Varyant (yaprak) kimligi bugun soyle kuruluyor (`recursive_homogenize.py`):

```python
by_cluster[row["ro_cluster"]].append(...)   # kumeleme TIP ICINDE
leaves.sort(key=lambda x: -len(x[0]))       # buyukluge gore
leaf_id = f"{cluster}#{leaf_index}"         # kimlik = tip adi + konum
```

Iki kirilganlik var:

1. **Kimlik tip adini tasiyor.** `5_505_qxyA#23`. Atama degisirse kimlik de
   degisir, oysa kume ayni dizilerden olusmaya devam ediyor olabilir.
2. **Indeks konumsal ve buyukluge gore siralanmis.** Bir tipin uyeligi
   degisince O TIPIN butun yapraklari yeniden numaralanir.

OLCULEN ORNEK. AerC eklendiginde (rev-0004 -> rev-0005) dort tipin uyeligi
degisti: 7_701_AerC +7, 1_106_TsaM -4, 1_113_CdnA +1, 1_105_DdmC +1. Bu dort
tipin yapraklari yeniden numaralandi; yaprak sayisi 2.769 -> 2.783 oldu.

KAPSAM, DURUSTCE. `leaf_index` TIP ICINDE sayildigi icin yalnizca uyeligi
degisen tiplerin kimlikleri kayar; dokunulmayan 64 tip etkilenmedi. Yani sorun
gercek ama sinirli. Yine de yayimlanmis bir URL'in anlaminin sessizce degismesi
kabul edilemez: biri `/leaf/1_106_TsaM-7.html` adresini alintiladiysa, sonraki
kosuda o adres baska bir kume gosterebilir.

## Neden simdi onemli -- kullanicinin gerekcesi

Tip adlari ve gruplar LITERATURDEN alindi, ve literaturde hata var. Bu oturumda
olculenler: bir onbaskinin ferredoksin iddiasi curutuldu, bu projenin kendi
"regulators diverge faster" iddiasi kismen geri cekildi, 3_313_KshA15'in ADI
yanlis (dizi M. tuberculosis'in, ad R. rhodochrous paralogunu dusunduruyor),
grup 6'nin aile etiketi "chlorophyll_pathway" iken kimyasal siniflandirmaya
uymadigi icin "tetrapyrroles" yapildi.

Dolayisiyla kume kimligini literaturden gelen bir ADA baglamak IKI KAT
kirilgan: hem atama degisince, hem makale duzeltilince kayiyor. Icerikten
turetilmis kimlik veriyi adlandirmadan AYIRIR -- bir makalenin iddiasi
duzeltildiginde etiket degisir, kimlik ayakta kalir.

## Neden simdi onemli (teknik)

Veritabaninin amaci, yeni Rieske oksijenazlar kesfedildikce buyumek. Her yeni
referans atamayi oynatiyor, her oynama kimlikleri kaydiriyor. Modulerlik tam
burada kiriliyor: **veri guncellenmeli, kimlikler durmali.**

## Tasarim

### 1. Kimlik icerikten turetilir

```
RO-<12 hex> = sha256(uye protein dizilerinin siralanmis ozetleri)[:12]
```

- Ayni dizi kumesi -> ayni kimlik, kacinci kosu oldugu onemsiz
- Uye girip cikarsa kimlik degisir -- OLMASI GEREKEN BU, artik baska bir kume
- Dokunulmayan kumeler kimligini KORUR (bugun ayni tipteki tum kardesleri
  kaydiriyor)
- Tip adi kimligin ICINDE degil, YANINDA bir alan

### 2. Tip adi bir etiket, kimlik degil

`leaf` tablosuna:

| alan | anlam |
|---|---|
| `variant_id` | `RO-7f3a9c21`, kalici |
| `nearest_type` | o anki en yakin kuratorlu tip (ETIKET, degisebilir) |
| `label_kind` | `named` (referansa benziyor) / `placeholder` (benzemiyor) |

### 3. Yer tutucu ad uzayi (GTDB mantigi)

Olculdu: 2.081 yaprak (10.093 giris) TAMAMI uzak/novel uyelerden olusuyor, yani
hicbir uyesi referansa benzemiyor -- ama hepsi `5_505_qxyA#23` gibi karakterize
bir enzimin adi altinda duruyor. O kume qxyA degil; yalnizca qxyA profiline en
yuksek skoru vermis.

GTDB adlandirilmis taksonla yer tutucu taksonu ayri ad uzayinda tutar. Burada
karsiligi: bu kumeler `label_kind=placeholder` olur ve sayfada
"adlandirilmamis kume, en yakin profil qxyA" diye sunulur -- "qxyA varyanti"
diye DEGIL.

Esik, kanit duzeyinden gelir (yeni esik icat edilmez): bir yaprakta
characterized / close_homolog / family_member duzeyinde en az bir uye varsa
`named`, yoksa `placeholder`.

### 4. Eski URL'ler kirilmaz

`freeze.py` zaten yeniden adlandirilan tipler icin 63 yonlendirme uretiyor;
mekanizma var. Eski `leaf_id -> variant_id` eslemesi bir dosyada tutulur ve
eski adresler yeni adrese yonlendirilir.

## Dokunulan yerler

| dosya | ne |
|---|---|
| `recursive_homogenize.py` | kimlik uretimi, yeni alanlar |
| `leaf`, `ro_leaf` semasi | `variant_id`, `nearest_type`, `label_kind` |
| `characterize_leaves.py`, `variant_signature.py`, `build_phylogeny.py` | leaf_id okuyan yerler |
| `webapp/app.py`, `templates/leaf.html` | sayfa, baslik, cerceve |
| `freeze.py` | yonlendirmeler |
| `validate_curation.py` | KARARLILIK KONTROLU (asagi) |

## Kabul olcutu

Tasarimin isledigini soyleyecek tek sinama: **ayni girdiyle iki kez kosuldugunda
butun kimlikler ayni kalmali**, ve **bir referans eklendiginde yalnizca uyeligi
gercekten degisen kumelerin kimligi degismeli.**

`validate_curation.py`'ye eklenecek kontrol: `analysis_out/variant_ids.json`
onceki kosunun kimlik -> uye-ozeti eslemesini tutar; yeni kosuda ayni uye
kumesinin ayni kimligi aldigini dogrular. Kimlik uretimi bozulursa kontrol
kirilir. (Bu depoda ayni sinif icin zaten uc sinama var: md5 parca anahtari,
SVG portu, iki render yolunun metni.)
