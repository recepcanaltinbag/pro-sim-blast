# Katalitik karboksilatın okunması: metot, doğrulama ve bulunan hata

**30 Ekim 2026** · ROAR-DB · hazırlayan: Claude (kurator denetimi için)

Bu rapor iki soruya cevap verir: **katalitik karboksilat nasıl okunuyor** ve
**doğru okunuyor mu**. İkinci soruyu denetlerken okumanın yanlış kolondan
yapıldığı ortaya çıktı; rapor bunu da, düzeltmesini de içerir.

---

## 1. Neden bu soru var

Rieske oksijenazının kanonik tanımı mononükleer demirde tam bir
**2-His-1-karboksilat** yüzeysel üçlüsü ister. ROAR-DB'nin kapısı ise Rieske
ligandlarının 4/4'ünü zorunlu tutarken üçlünün yalnızca **2/3'ünü** istiyor —
yani bir giriş karboksilat okunmadan da kabul edilebiliyor. Soru şuydu: o
üçüncü ligand gerçekten yok mu, yoksa **okuması** mı bozuk?

## 2. Metot: her adımda neye bakılıyor

Tam çıktı: `python3 carboxylate_method_card.py`
(kayıtlı kopya `analysis_out/carboxylate_method_card.txt`, İngilizce — yayına
gidebilecek biçimde).

| adım | neye bakılıyor | neden bu adım var |
|---|---|---|
| 1 | tek sabit hizalama kolonunda oturan kalıntı; saklı tablo hizalamadan yeniden hesaplanıp karşılaştırılıyor | "yok" denen şeyin ne olduğunu ikiye ayırmak: **boş kolon** (hizalama olayı) ile **ikame** (biyolojik iddia) aynı şey değil |
| 2 | her girişin ilk/son **dolu** match kolonu | kolon girişin hizalanmış bölgesinin dışında kalıyorsa orada okunacak bir şey yoktur |
| 3 | aynı kolon ±2 / ±5 / ±10 / ±20 / ±40 | ucuz çözümü (pencereyi gevşetmek) elemek; boş kolonun etrafını taramak boş kolon taramaktır |
| 4 | His-2'ye kalibre edilmiş ham uzaklıkta D/E | bu okumanın **kanıt olmadığını** göstermek: pencere ~90 kalıntı geniş olmak zorunda ve şans beklentisi isabetle neredeyse aynı |
| 5 | referanstan kesilen 3 alt-model (Rieske / katalitik lob / karboksilat çevresi), 3 kümeye hmmsearch | kusurun **dizide mi profilde mi** olduğunu ayırmak; Rieske alt-modeli pozitif kontrol |
| 6 | **asıl ölçüm:** tip içi mafft; karboksilat kolonu o tipin pozitiflerinin **oyu** | sabit sayı yerine tipin kendi koordinatı; ±3 kolon pencere |
| 7 | çapa kolonunun D/E'si, grubun **kendi** D/E sıklığına karşı binom (BH düzeltmeli) + aynı hizalamadaki bütün kolonların D/E dağılımı | bulunan korunmuş kolonun şans olup olmadığını sınamak |
| 8 | bulunan karboksilatın His-2'ye ham uzaklığının sapması | sabit kolonun **neden** kaydığını göstermek (His çifti ile karboksilat arasındaki ~140 kalıntılık parçada biriken indel) |

Adım 6'nın kendini denetleyen üç yeri var ve çıktıda görünür: kolonun
**oybirliği oranı**, tipin kuratörlü **referansının** çapadaki kalıntısı, ve
yöntemin aynı hizalamadaki **geri kazanımı**.

## 3. Metodun kendi doğrulaması (yeni ölçüm)

Adım 6'da bir döngüsellik vardı: kolonu pozitiflerin oyu seçiyor, sonra
doğruluk **aynı** pozitiflerde ölçülüyordu. Bunu kırmak için
`carboxylate_anchor_holdout.py` yazıldı: pozitifler ikiye bölünüyor, kolonu
yalnızca **yarısı** oyluyor, doğruluk oyu olmayan diğer yarıda ölçülüyor.

| ölçüm | sonuç |
|---|---|
| dışarıda bırakılan pozitiflerde geri kazanım | **490/490 = 1,000** |
| yarı veriyle seçilen kolon, tam oyun kolonuyla aynı mı | **10/10 tip** |
| en küçük oy marjı | 8 oy (NdmA), diğerleri 41–60 |

Yani **çapalama makinesi doğru çalışıyor**; döngüsellik şüphesi ölçümle
kalktı. Bu, raporun geri kalanındaki hatayı metodun kendisinden ayırmayı
mümkün kılıyor: metot sağlam, **hedefi** yanlıştı.

## 4. Literatür denetimi

| iddia | literatür | bizdeki durum |
|---|---|---|
| üçlü = 2 His + **Asp veya Glu** | Barry & Challis 2013 (ACS Catal., doi:10.1021/cs400087p): "2-His-1-carboxylate triad … an aspartate/glutamate residue" | `CARBOX_OK = "DE"` — **doğru** |
| köprü Asp bir demir ligandı **değil**, elektron transferi aracısı | aynı derleme; ayrıca 2025 derlemesi (PMC12445004): "a highly conserved aspartate residue … forms a hydrogen-bonding bridge that links the mononuclear iron site of one subunit to the Rieske cluster of a neighboring subunit" | ligand/röle ayrımımız ve 8/9 yapısal ölçümümüz — **doğru**, hatta derlemenin ifadesiyle birebir |
| RedG'de köprü Asp'ı **Glu** ile yer değiştirmiş görünüyor | aynı derleme: "the 'bridging' aspartate residue implicated in electron transfer … appears to be mutated to glutamate in RedG" — çekingen ifade, **işlevsel yorum yok** | sitede "An open question, not a result." dediğimiz yer **literatürle aynı duruşta**; dahası bizde emsal var: RedG tipi 44 D / 20 E, CntA 1693 E, Stc2 1361 E — derlemede tek cümlelik gözlem, bizde ölçülmüş bir dağılım |
| üçlüyü **tamamen** kaybetmiş RO kataloğu | yok; aranan derlemelerde böyle bir liste bulunmadı | o halde bizim "gerçekten karboksilatsız" dediğimiz girişler **şüpheli sayılmalı**, yeni biyoloji iddiası olarak sunulmamalı |
| karboksilat genellikle **Asp** | aynı derleme + NDO/KF707 numaralandırması His208, His213, **Asp362** | 4. bölümdeki hatanın ortaya çıktığı yer |

## 5. BULUNAN HATA: karboksilat yanlış kolondan okunuyordu

### Ne oldu

Model 603 kolonluk `ROmotif71`. Üç farklı kolon numarası dolaşımda:

| kolon | nereden geliyor | 76 referansta ne var | durum |
|---|---|---|---|
| 355 | `ro_motif.py` içindeki gömülü eski sabit (71 referanslı, 426 kolonluk ESKİ hizalamaya ait) | – 25, G 19, K 8, D 6 | **anlamsız** |
| 471 | `motif_columns.json` — türetme betiğinin seçtiği | E 42, D 29 (D/E %93,4) | **yanlış site** |
| **431** | NDO'nun yapıdan bilinen demir ligandı | **D 68** (D/E %89,5) | **doğru** |

Kanıt zinciri:

1. 1NDO'da demir ligandları His208 / His213 / **Asp362** (literatür; ayrıca
   kendi yapısal ölçümümüz Asp362'yi demirden 2,27 Å'da buldu).
2. `refs71` içindeki NahAc dizisi PDB numaralandırmasından **3 kalıntı**
   geride başlıyor: ham 202 / 205 / 210 → PDB 205 / 208 / 213. Köprü Asp205 ve
   iki His bu haritalamayla **tam** oturuyor.
3. Aynı haritalamada ham **359 = PDB Asp362**, ve ham 359 match kolonu
   **431**'e düşüyor. Ham 407 ise kolon 471'e — yani 471 başka bir korunmuş
   karboksilat.
4. 19.114 onaylanmış giriş üzerinde: kolon 431'de D 16.287 / E 439 (±2 pencere
   ile **%87,5** D/E), kolon 471'de E 7.583 / D 6.008 (**%71,1**). Literatürün
   "genellikle Asp" ifadesine uyan kolon 431; 471'in E/D karışımı o tarife
   uymuyor.

### Kök neden

`build_reference_models.py` kolonu **korunma oranını maksimize ederek**
seçiyordu: `his[h1] + his[h2] + acid[d]`. Kolon 471'in D/E korunması (%93,4)
kolon 431'den (%89,5) yüksek, yani doğru kolonu **3 referans farkla**
kaybediyor. "En korunmuş D/E kolonu" ile "demiri bağlayan D/E kolonu" aynı şey
değil — ama betik bunları aynı sayıyordu.

### Hata neden sessiz kaldı

Yanlış kolon hiçbir yerde patlamadı: tablo doldu, istatistikler üretildi, site
yayınlandı. Dahası, yanlış kolon **kendi açıklamasını da üretti**: boş kolonda
okunan "karboksilat yok" çağrıları "üyelerin %12'sinde karboksilat eksik"
bulgusuna, o da 8 adımlık bir denetime ve "kapı 2/3'e toleranslı olmak zorunda"
gerekçesine dönüştü. **Bulgunun kendisi artefakttı.**

### Etkisi

| ne | önce (kolon 471) | sonra (kolon 431) |
|---|---|---|
| üçlüsü 3/3 tam olan giriş | 13.543 | **16.685** |
| onaylanmışlarda karboksilat okunma oranı | %71,1 | **%87,5** |
| **üyelik** değişimi (Rieske'si tam olanlar) | — | **+12 / −1** (19.114 içinde) |

Yani üyelik pratikte değişmiyor; değişen şey **aktif bölge hakkında
söylediğimiz her şey**.

### Düzeltme (yapıldı)

1. `build_reference_models.py`: karboksilat kolonu artık **seçilmiyor,
   yapıdan okunuyor** (`STRUCTURAL_ANCHOR` = 1NDO Asp362 → NahAc ham 359).
   Korunma yalnızca saplama olarak raporlanıyor; çapa doğrulanamazsa betik
   sessizce devam etmiyor, **hata verip duruyor**.
2. `ROs_71_Clean/motif_columns.json` yeniden üretildi; karboksilat **431**,
   ve dosya artık `carboxylate_source` ile kolonun **neyle** seçildiğini
   kaydediyor. Yeniden üretim, farkı ekrana basıyor:
   `[UYARI] korunma en cok kolon 471'i isaret ediyor (D/E 0.9342) ama yapisal
   capa kolon 431 (D/E 0.8947). CAPA KAZANIR: kolon 431.`
3. `ro_motif.py`'deki gömülü eski sabitler (355 dahil) ölçülmüş değerlerle
   değiştirildi ve neden yanlış oldukları yazıldı — JSON silinse bile sessiz
   yanlışlık dönmesin.
4. Kapıya **iki yeni kontrol** eklendi (`validate_curation.py`, artık 107
   kontrol): kullanılan kolon yapısal çapayla aynı mı, ve saklı kolon
   korunmayla mı çapayla mı seçilmiş. İkisi de **PASS**.

## 6. Bu düzeltmeyle geçersizleşen sayılar

Bunları hiçbir yerde kullanmayalım, çünkü hepsi yanlış kolondan:

- `analysis_out/carboxylate_audit.json` (7 Ekim): **1.414** karboksilatsız,
  **1.388** çapalamayla kurtarılan, **26** gerçek yokluk, %12,4, 8,9× vb.
- Dün ürettiğim **metot kartındaki bütün sayılar** (kolon 355 satırı dahil).
  Kartın **yapısı** geçerli — adımlar ve neye bakıldığı doğru — **sayıları**
  yeniden ölçülmeli.
- Kapının "triad 3/3 tam olan üye" istatistiği (13.543) ve ona dayanan site
  metni.
- Karboksilat üzerinden kurulmuş klad kayması / pencere genişletme
  karşılaştırmaları.

**Geçerli kalanlar:** köprü karboksilat ölçümleri (kolon 286 doğrulandı, 8/9
yapısal ölçüm), ligand/röle ayrımı, oksijen kaderi ve metabolik rol tabloları,
reaksiyon doğrulamaları — hiçbiri bu kolona bağlı değil.

## 7. Sırada ne var

1. **Yeniden kurma** (zorunlu): `ro_carboxylate` ve `catalytic_intact` kolon
   431 ile yeniden üretilmeli. Bu, pahalı `rebuild.py --run --include-heavy`
   adımıdır ve zaten bekleyen 3 bayat `ro_search` kaydını da düzeltir.
2. Denetimi **kolon 431'de yeniden koşmak** (`carboxylate_audit.py`) ve metot
   kartını yeniden üretmek. Beklenti: "karboksilatsız" küme büyük ölçüde
   erir, çünkü artık boş kolon okunmuyor.
3. Kalan gerçek yoklukları yeni sayılarla incelemek. Şimdiden bilinen üç
   sorulu nokta: referans **PobA** kolon 431'de K taşıyor, referans
   **KshA15** orada boşluk veriyor, ve **3_309_NahAc** hiç onaylanmış üye
   toplamıyor (kapı WARN'ı) — NahAc neredeyse tıpatıp komşularına
   (NBDO/3NTDO) kaydığı için.
4. Site metninin güncellenmesi ve ardından yayın.

## 8. Dürüst değerlendirme

Sen "hata yapmayalım, hepsini güncel literatürden bakalım" dediğinde haklıydın
ve bu rapor bunun neden haklı olduğunun kanıtı. Metot iyiydi — sekiz adım,
null modeli, kendi doğruluğunu ölçmesi, hatta dışarıda bırakma testinde
490/490. **Ama sekiz adımın hepsi yanlış kolona bakıyordu.** Bu, bir ölçümün
ne kadar titiz olduğunun, neyi ölçtüğünü garanti etmediğinin örneğidir; ve
hatayı yakalayan şey istatistik değil, **yapıya ve literatüre geri dönmek**
oldu: NDO'nun Asp362'si.

Düzeltmenin kalıcı kısmı da bu: kolon artık korunmadan değil **yapıdan**
geliyor ve her kapı koşusunda yeniden doğrulanıyor.
