# ROAR-DB — sonraki adimlar

Durum: `main` dalinda pipeline + web uygulamasi + Atlas yayinda,
`gh-pages` dalindan statik site: https://recepcanaltinbag.github.io/pro-sim-blast/

## Sirada

### 0. HEDEF (kullanici, 2026-11-01): "NCBI'daki, ordaki, burdaki TUM RO olabilecek her sey bu DB'den toparlanabilsin"

Bu, DB'nin biriminim degistiriyor. Bugun yayinladigimiz sayi (19.125) bir
KAPSAM degil, bir KESISIM: hem genomik baglami olan, hem kapiyi gecen
girisler. Hedef, DB'yi her kaynaktan gelen RO adaylarinin BIRLESIMI yapmak
ve her adayin hangi kaynaktan, hangi kanitla geldigini yaninda tasimak.

Kaynak envanteri (2026-11-01 itibariyla elimizde olan ve OLMAYAN):

| kaynak | durum | sayi |
|---|---|---|
| NCBI GenBank/RefSeq genom kayitlari (anotasyonlu) | TARANDI | 19.125 giris |
| UniProtKB anlik goruntu (canli) | INDIRILDI | 25.706 |
| UniProtKB anlik goruntu (silinmis, yalnizca bizde) | INDIRILDI | 47.993 |
| NCBI `nr` / RefSeq protein uzayi | **HIC TARANMADI** | ? |
| NCBI anotasyonu OLMAYAN genomlar (ORF cagrisi gerekir) | **HIC TARANMADI** | ? |
| MGnify / JGI metagenom protein katalogu | **HIC TARANMADI** | ? |

Madde 95 bunu zaten yazmisti: protein evreni UniProt'tan geldi, NCBI'DAN
DEGIL. Yani kullanicinin sordugu seyin en buyuk parcasi hic aranmamis.

Sira (yukaridan asagiya, her adim kostururak dogrulanir):

**A1. `ro_union`: tek satir = tek DIZI, yaninda tum kaynak kimlikleri.**
Bugun ayni alpha alt birimi uc tabloda uc kere sayiliyor (`ro`,
`ro_protein`, `ro_protein_status`) ve hangi sayinin ne oldugu sayfadan
okunamiyor. Birim DIZI olmali; accession kaynak kimligidir, birim degil.
Cikti: `ro_union(seq_md5, length, n_sources, has_genomic_context,
best_evidence_tier, ro_cluster, ro_group, triad_state)` +
`ro_union_source(seq_md5, source, accession, url)`.

**A2. Protein-duzeyi girisler BIRINCI SINIF.** Manset birimi "RO alpha alt
birimi" olur. Silinmis girisler SAKLANIR ve ayri sayilir -- 34.442 tekil
dizi icin elimizdeki kopya tek kopyadir (madde 112).

**A3. PAYDA MUHAFIZI.** Konuma bagli bir tablodan gelen her yuzde, paydasini
"genomik baglami olan kume" olarak beyan etmek ZORUNDA; etmezse FAIL.
Yoksa A1/A2'den sonra her yuzde sessizce yanlis paydaya oturur.

**A4. NCBI PROTEIN UZAYINI TARA (en buyuk acik is).** ROmotif71 +
Rieske/Ring_hydroxyl_A ile RefSeq tekrarsiz protein setinde hmmsearch.
Burasi "NCBI'daki her sey"in asil karsiligi.

**A5. 18.076 kurtarilabilir protein icin genom kaydini indir**, operon ve
komsuluk kur. Bunlar bugun protein-duzeyinde duruyor ama genomik baglami
KURULABILIR; kurulunca kesisime girerler.

**A6. CAPRAZ REFERANSLAR alan olur:** EC (3.853 tanesi `ro_protein_status`'ta
hazir), UniProt, InterPro, KEGG, MetaCyc. Bir toparlama veritabaninin isi
ayni molekulun her yerdeki adini bir arada tutmaktir.

### 1. Literaturle uyum (acik maddeler)

- ~~Batie sinifi tip sayfalarinda gosterilmeli~~ **IPTAL (2026-11-01),
  kendi veri dosyam bu maddeyi curuttu.** `batie_classes.yaml` zaten
  "A batie_class column cannot be keyed on enzyme type" diyor ve
  gerekcesini veriyor: CARDO tabloda UC ayri sinifta geciyor (KA1 susu
  IIA, IC177 IIB, CA10 ve J3 III), cunku Batie sinifi enzimin degil
  SUSUN -- daha dogrusu elektron zincirinin -- ozelligi. Tip sayfasina
  bir Batie sutunu koymak, dosyanin tam da uyardigi karistirmayi
  uretirdi. Literaturle uyumsuzluk Batie'nin EKSIK olmasi degil; onu
  tip duzeyinde yayinlamak olurdu.
- ~~1-5 grup numaralari etiketlerle degistirilmeli~~ **YAPILDI ama baska
  sekilde (2026-11-01).** Asil sorun numaranin kendisi degildi: numara
  TANIMSIZDI. `groups.yaml` hic yayinlanmiyordu, yani site her yerde
  "group 3" yaziyor ve grup 3'un ne oldugunu soyleyen tek bir sayfa yoktu.
  Tanimsiz bir numara okuyucunun bildigi baska bir bes'li semaya yapisir.
  `/atlas/groups` eklendi; tip sayfasi ve giris sayfasi numarayi oraya
  BAGLIYOR. Numarayi etiketle degistirmek yerine numaraya tanim verildi --
  kalici baglantilar ve `ro_group` degerleri boylece bozulmadi.
  Ayrica: groups.yaml'in `members`/`types` alanlari BAYATTI (grup 1 icin
  5746 uye / 14 tip yaziyordu, olcum 5104 / 16; grup 6 icin 0 / 0
  yaziyordu, olcum 617 / 3). Alanlar dosyadan silindi, sayfa olcumden
  yaziyor, ve `check_groups_published` dosyaya sayi geri konursa FAIL
  veriyor.
- Sekiz uyesiz referans. **KISMEN YAPILDI (2026-11-01):** birlesim
  kurulduktan sonra olculdu, DORDUNUN protein duzeyinde uyesi VAR
  (BpdC1 13, XylC1 3, NahAc 1, NinA 1). Tip sayfasi "bu tipin kendine ait
  uyesi yok" derken bunu soylemiyordu, yani var olan veriyi yok gibi
  gosteriyordu. `empty_types.py` artik her kayda `union` alanini yaziyor
  ve sayfa "diziler BURADA var, yalnizca genomik baglami yok" diyor.
  ACIK KALAN: bu dizilere genom kaydi indirilip gercek uye yapilmasi
  (plan A5); dordu de o zaman kesisime girer.
- Adlandirma karari: NahAc/NBDO ve CarAa/CARDO (%99,7-100 ayni).

### 2. Kalan gorsel denetim bulgulari

- ~~`/clusters` tablo tasmasi ve carpisma tablosu~~ **YAPILDI (2026-11-01)
  ve 2 tablo degil 15 TABLOYMUS.** `.tablewrap` (overflow-x: auto) eksikti:
  cluster.html'de 9, ro.html'de 1, atlas_ecology/atlas_regulation'da 4,
  bir de kendi yeni ekledigim atlas_groups.html'de 1.
  `check_tables_are_scrollable` artik sarilmamis her tabloyu FAIL veriyor --
  cunku kesilme DAR degil GENIS ekranda oluyor ve mobil testinde gorunmuyor.
- ~~Arama sayfasinin bos baslangic durumu~~ **YAPILDI (2026-11-01).** Sorgu
  ve filtre yokken sayfa 19 binden fazla girisi hicbir yonlendirme olmadan
  dokuyordu. Artik ne arandigini, sol panelin ne yaptigini ve iki
  tiklanabilir baslangic ornegi soyleyen bir blok var; ilk filtrede kayboluyor.
- ~~Mobilde Plotly kirpilmasi~~ **YAPILDI (2026-11-01).** Sebep
  `responsive: true`nin eksikligi DEGILDI -- o zaten her grafikte vardi ve
  genisligi izliyor, KENAR BOSLUGUNU izlemiyor. Sol kenar boslugu sabit
  180 px'ti; 360 px'lik bir telefonda grafige 180 px kaliyor ve kategori
  etiketleri kirpiliyordu. `plot.js`'teki `plotMargin()` kenar boslugunu
  ekran genisliginden turetiyor (en az 44 px, en cok istenen, ve ekranin
  %38'i), `bindPlotMargin()` ekran dondugunde yeniden duzenliyor. Dort
  grafik baglandi: tax, ev, ops, heat.
  NOT: plot.js base.html'e `defer` OLMADAN eklendi. table.js defer'li ve
  sayfa ici Plotly betikleri ayristirma aninda kosuyor -- defer'li bir
  dosyadan gelen yardimci onlar icin tanimsiz olurdu.
- ~~`1_118_RedG` Sekil 1 altyazisina kacan gelistirme notu~~ **YAPILDI
  (2026-11-01) ve BIR TIP DEGIL ALTI TIPMIS.** Kullanici RedG'yi
  bildirdi; olcum `1_117_PrnD`, `1_118_RedG`, `3_321_Nvd`,
  `3_322_DAF36`, `6_601_PAO`, `7_701_AerC`'i buldu. Ikisi `notes`
  (Sekil 1 altyazisinin "Note." kismi), digerleri `source`
  ("Curated from:") alanindaydi -- yani altyazi tek sizinti yeri degildi.
  Icerdikleri `analysis_out/...` dosyalarinin HICBIRI indirilebilir
  degil, yani okuyucuya acamayacagi bir referans veriliyordu. Birinde
  de birinci tekil sahis kalmisti ("The family members I retrieved").
  Olgular korundu, referanslar okunabilir ifadeye cevrildi, ve
  `check_no_internal_refs_published` kuratorlu bes alanda yol / .json /
  birinci tekil sahis / TODO gorurse FAIL veriyor.


89. **ADIM 2 OLCULDU: TIP ICI CAPALAMA UYELIGI DEGISTIRMIYOR -- kendi
    uyarimi curuttum (2026-10-28).** Kullanici "literatur acisindan hangisini
    onerirsen oyle devam edelim" dedi; onerim "once olc, sonra uygula"ydi ve
    gerekcem "uyelik degisecegi icin kullanici karari gerekir"di. OLCUM BU
    GEREKCEYI CURUTTU.

    (a) UYELIK BUYUYEBILIR MI: pratikte hayir. Rieske'si saglam ama katalitik
    testi gecemeyen 2.700 adayin ortalama model kapsamasi **0,248** ve hmm
    skoru **44**; gecenlerde **0,830** ve **401**. 2.700'un **2.316'si
    (%85,8)** kapsamasi 0,30'un ALTINDA. Bunlar kolonu kaymis tam RO'lar
    DEGIL, parcalar -- kapsamasi 0,25 olan bir dizide capalanacak bir
    C-terminal katalitik lob yoktur. Kapsamasi >=0,45 olan yalnizca **99**
    tane var, yani ust sinir onaylanmis setin **%0,52'si**.

    (b) UYELIK KUCULEBILIR MI: esik 2/3'te kalirsa **sifir**. Tip ici
    capalama karboksilati GERI GETIRIYOR, kimsenin sitesini
    kaybettirmiyor; o 1414 uye bugun zaten iki His ile geciyor.

    (c) ASIL DEGISEN SEY ISTATISTIK VE IDDIA. Bugun `ro_carboxylate` 1414 uye
    icin "karboksilat yok" diyor ve bunun 1388'i (%98,2) sabit kolon
    yanilgisi. Tip ici capalamayla triad'in TAM oldugunu soyleyebildigimiz
    uye sayisi 10.008 -> 11.396 olur, gercek yokluk 26'ya iner, ve kapi
    "2/3 yeterli" demeyi birakip "3/3 olcuyoruz, 26 gercek istisna var"
    diyebilir. Bu, literaturun kanonik tanimiyla (2-His-1-karboksilat
    yuzeysel uclusu) ortusur.

    (d) KOK SEBEP, denetimin kendi verdict'inden: eksik cagrilarin %94,9'unda
    ROmotif71 hizalamasi karboksilat kolonuna HIC ULASMIYOR -- son dolu match
    kolonu medyani 255, kolon ise 355. Yani sorun bir-iki pozisyonluk kayma
    DEGIL, global profilin bu dizilerin C-terminal lobunu hizalayamamasi.
    Tip ici hizalama kolay, cunku orada diziler birbirine benziyor.

    (e) SONUC: degisiklik UYGULANABILIR ve riski dusuk. Baslangictaki
    "uyelik degisir, kullanici karari" uyarim OLCULMEMIS bir endiseydi ve
    olcum onu curuttu -- bugunun yedinci ornegi olarak kayda geciyor.

88. **"KOLON EVRIMSEL UZAKLIGA GORE KAYMIS OLABILIR, BELKI BIRKAC MODEL
    GEREKIR" -- fikir DOGRU, uygulama bicimi olcumle degisti (2026-10-28).**

    (a) KAYMA GERCEKTEN KLAD YAPILI. Yapidan bilinen 26 gercek karboksilat,
    dogrulanmis esleme TERS yonde kullanilarak model kolonuna cevrildi. Uc
    SIKI banda dusuyorlar:
      407-411 : NdmB 407, NdmA 408, NagGH 411
      431-432 : 20 tip
      449-451 : ROCH34 449, KshA15 451, PDO 451
    NdmA ve NdmB ayni alt aile ve komsu kolonlarda -- kayma rastgele degil.

    (b) AMA BIZIM BOLMELERIMIZ BANDI ONGORMUYOR. RO grubu: grup 1'de hem
    407/408 hem 431, grup 3'te 411/431/451, grup 4'te 431/449/451. Protein
    uzunlugu da ayirmiyor (407-411 bandi 369/373/442 aa, 431-432 bandi
    334-470 aa).

    (c) DIZI BENZERLIGI DE GUVENILIR BIR YONLENDIRICI DEGIL. Ayni bantta
    olan tip ciftleri biraz daha benziyor (ortanca %18,8 karsisinda %16,9,
    ortalama %25,7 karsisinda %18,2) ama en yakin referansa gore bant tahmini
    **23/26**. Ve en kotu hata EN YUKSEK KIMLIKLI ciftte: NagGH (407-411) ile
    TPDO (431-432) birbirinin en yakin referansi, **%40,4** kimlikle, ve
    farkli bantlarda. Benzerlige dayali bir yonlendirici tam onemli olan
    vakalari yanlis yonlendirir.

    (d) DEPO BU CoZUMU ZATEN SINAMIS. `carboxylate_audit.json`'un
    `anchored_reread` asamasi: karboksilat sabit global kolon yerine TIP
    ICINDE yeniden hizalanip o tipin KENDI karboksilat pozisyonlarina
    capalandiginda, sabit kolonun "yok" dedigi 1414 uyenin **1388'i (%98,2)**
    geri geliyor (1369 homolog kolonda, 19 grup-korunmus kolonda). Buna
    karsilik pencereyi genisletmek neredeyse hicbir sey kurtarmiyor: ±5'te
    %0,8, ±10'da %5,1. Yani CARE PENCERE DEGIL, YENIDEN CAPALAMA.

    (e) ONERILEN TASARIM. Onceden kesilmis alt-modellere YONLENDIRMEK yerine
    (yonlendirme 23/26 ve en yuksek kimlikli ciftte basarisiz), her tipin
    KENDI capasini tasimasi: tip ici yeniden hizalama + o tipin bilinen
    karboksilat pozisyonlarinin modu. Yonlendirme sorunu boylece hic ortaya
    cikmaz. Bu, kullanicinin "birkac model" fikrinin olcume dayanan hali.

    (f) KULLANICI KARARI GEREKIYOR: bu degisiklik kapinin >=2/3 esigini ve 55
    tek-His uyesini yeniden degerlendirmeyi gerektiriyor, yani UYELIK
    degisecek. Veri degisikligi degil, yayimlanan setin tanimi degisiyor.

87. **"CATALYTIC UCLUDEKI EVET AMA KOPRU FARKLI" -- kullanici hakliydi, iki
    karboksilat AYNI TUR BOLGE DEGIL (2026-10-28).** On sayfa ikisini tek
    tabloda ayni "Asp : Glu" cercevesinde sunuyordu ve bu, kullanicinin bir
    onceki sorusunda yasadigi karisikligin kaynagiydi.

    (a) OLCULDU (1NDO, zincir A). Mononukleer demiri KOORDINE EDENLER:
    His208 NE2 2,04 A, His213 NE2 2,15 A, **Asp362 OD2 2,27 A**. Yani
    2-His-1-karboksilat uclusu ve karboksilati Asp362; birinci kabuk LIGANDI.
    Kopru **Asp205** ayni demirden **6,95 A** ve HICBIR metali koordine
    etmiyor: ikinci kabuk ROLESI, bir yanda His208'e, obur yanda alt birim
    birlesimi uzerinden komsunun Rieske ligandi His104'e hidrojen bagli.
    Biri metali TUTUYOR, obur elektronu ona TASIYOR.

    (b) SAYFA DUZELTILDI. Baslik "Two carboxylates, and they are not the same
    kind of site" oldu; satir etiketleri artik "first-shell ligand" ve
    "second-shell relay" diyor ve 1NDO'nun rakamlari metinde duruyor. Ikisi
    tabloda yan yana KALDI ama NEDEN yan yana oldugu yazildi: ikisi de
    hizalama kolonundan okunan karboksilat, ve asil bulgu o iki kolonun
    denetimden FARKLI cikmasi.

    (c) OLCULMEMIS BIR HIPOTEZ, iddia olarak YAZILMADI. Iki kolonun farkli
    sonuc vermesi tesadufi olmayabilir: birinci kabuk ligandi katalitik lobdaki
    ESNEK bir ilmekte duruyor ve hizalama orada kayiyor -- karboksilatin 26
    yapida 431, 432, 407-411 ve 449-451'e dagilmasinin sebebi bu olabilir.
    Kopru ise iki alt birimin ARAYUZUNDE, yapisal olarak kisitli bir yerde;
    orada kayma olmuyor, o yuzden tek kolon 8/9 tutuyor. TEST EDILEBILIR: iki
    kolonun 76 referans uzerindeki hizalama kayma dagilimini karsilastirmak
    yeterli. Olculmedigi icin sayfaya KONULMADI.

    (d) BUGUNUN DERSI, uc kez tekrarlandiktan sonra: olctugum sey ile
    sonuclandirdigim sey arasindaki bosluk. CMO'da iki kez (urun kimligi),
    koprude bir kez (pozisyon esdegerligini islev esdegerligi sanmak). Ucunde
    de kullanici sormasa duzelmezdi. KURAL: bir cumle yayimlanmadan once
    "bunu hangi olcum destekliyor" sorusu cumlenin HER YUKLEMI icin ayri ayri
    sorulacak.

86. **KULLANICININ ITIRAZI DOGRU BIR SORUYDU VE YAYIMLANAN CUMLEYI DUZELTTI
    (2026-10-28).** Kullanici "C202'nin kopru olduguna emin misin, sistein
    sanki Rieske merkezde demir-sulfuru tutan kalinti" dedi. Itiraz yerinde:
    sistein klasik bir Rieske ligandi (merkez 2 Cys + 2 His). Yapidan olculdu.

    (a) ITIRAZ CURUTULDU, ama olcumle. 3VCA'da Rieske kumesini tutan
    kalintilar **Cys86, Cys106, His88, His109** (2,33 / 2,36 / 2,16 / 2,17 A)
    ve C202 kumeden **47,97 A** uzakta; bir ligand ~2,3 A olurdu. CntA
    E205 icin ayni tablo (46,46 A). Yani ikisi de Rieske ligandi DEGIL.

    (b) YONTEM BAGIMSIZ OLARAK DOGRULANDI. Tek alt birimli yapilarda
    olctugum 44-48 A, literaturun bildirdigi "demirden KENDI alt biriminin
    Rieske merkezine 43,5 A, alt birim birlesimi uzerinden 12 A" degerleriyle
    ortusuyor. Yani o buyuk sayilar artefakt degil, bilinen ic mesafe; kopru
    komsuya uzaniyor, kendi kumesine degil.

    (c) AMA BENIM YAYIMLADIGIM CUMLE FAZLA IDDIAYDI. "Ikame elektron transfer
    koprusunun tam kendisinde" demisim. Olculen sey yalnizca POZISYON
    esdegerligi: ikisi de kendi demirinden 6,9 ve 6,0 A uzakta, yani
    dogrulanmis Asp205 ile ayni mesafede. Kopru ISLEVINI yaptiklari
    SONUCLANAMAZ, iki sebeple:
      * Iki yapi da TEK alt birim iceriyor, yani koprunun uzanacagi komsu
        kume dosyada YOK.
      * Parales 1999, NDO'da Asp205'i Glu ile degistirdiginde aktif bolge
        demirine giden yol KAYBOLUYOR, NADH->Rieske adimi duruyor. Yani
        glutamat NDO'da bu isi YAPMIYOR. Dogal olarak glutamat tasiyan bir
        enzim bu yuzden bir DOGRULAMA degil, bir BILMECE.
    On sayfa buna gore yeniden yazildi ve bolum "An open question, not a
    result" diye etiketlendi. Cozum, bu enzimlerden birinin oligomeri saglam
    bir yapisini gerektiriyor.

    (d) LITERATUR CAPASI: Kauppi et al. 1998 (yapi, Asp205 kendi alt
    biriminde His208'e ve komsu alt birimde Rieske ligandi His104'e hidrojen
    bagli) ve Parales et al. 1999 (mutagenez, J Bacteriol 181:1831).

85. **KOPRU KARBOKSILAT KOLONU YAPIYA KARSI DOGRULANDI VE GECTI -- katalitik
    kolonun tam ziddi (2026-10-28).** Kullanici "kopru karboksilati elektron
    tasimasiyla ilgili aslinda, belki ilginc sonuc verebilir" dedi. Dogru
    sezgi, ve kolon dogrulanmadan ustune hicbir sey kurulmadi.

    (a) ONCE ESLEME. Kolon -> kalinti numarasi eslemesini kendim kurdum ve
    denetimin katalitik cagrilarina karsi **0/26 TUTTU**: sebebi
    `read_stockholm_matchcols`in insersiyon durumlarindaki kalintilari
    atmasi (OxoO'da 64 kalinti fark). Kendi eslememi ATTIM,
    `carboxylate_audit.AlignmentIndex`i ice aktardim ve **26/26** tuttu.
    Ancak ondan sonra geometriye gecildi.

    (b) OLCUM (`bridging_geometry_audit.py`). Kopru kolonu 286. Komsu alt
    birimin dosyada OLDUGU 9 yapinin **8'inde** sikici bir kopru geometrisi:
    komsu zincirin Rieske kumesine 5,92-7,15 A (ortanca 6,6), KENDI zincirinin
    demirine 6,65-10,55 A (ortanca 7,1). 1NDO icin kolon Asp205'i veriyor --
    karboksilati bir yanda Rieske ligandina, obur yanda demir ligandina
    hidrojen bagli olan kalinti. Tek aykiri GxtA (6WNC, 36,43 A).
    **Karsilastirma: katalitik kolon 26/26 yapida demirden 10,9-14,7 A
    uzakti.** Yani iki kolon ayni testten gecirildi ve biri gecti, obur
    gecmedi.

    (c) BENIM IKINCI HATAM, aninda yakalandi. `spans_two_chains` testini TERS
    kodlamistim (Rieske kendi zincirinde, Fe komsuda) ve 12/12 "hayir"
    vermisti. Gercek yon tersi: kopru, DEMIRI tasiyan alt birimde oturup
    KOMSUNUN Rieske merkezine uzaniyor. Duzeltilince 9/9 oldu. Bir bulgu
    sanip yayimlamak cok kolaydi.

    (d) SINIR DURUSTCE: 27 hedefin 12'si olculebildi. Tek alt birimli 8
    yapida komsunun Rieske'si DOSYADA YOK, o yuzden Rieske uzakligi 44-48 A
    cikiyor (kendi merkezine olan mesafe) ve anlamsiz; demir tarafi orada da
    5,99-7,18 A. Besinde kalinti yapida cozulmemis, yedisinde dosyada Fe ya
    da Rieske yok.

    (e) KULLANICININ SEZGISI BIR BULGUYA DONDU. Kuaterner amin dalinda bu
    pozisyon aspartat DEGIL: CntA'da **Glu205**, Stc2'de **Cys202**, ve ikisi
    de kendi demirlerinden 6,89 ve 5,99 A uzakta, yani halka-hidroksilleyenlerdeki
    Asp ile AYNI yapisal pozisyonu dolduruyorlar. Ikame, elektron transfer
    koprusunun TAM KENDISINDE. On sayfadaki "tek bir kolon substrat kimyasini
    izliyor" iddiasi artik dogrulanmis bir kolona dayaniyor.

    (f) ON SAYFA DUZELTILDI. Madde 83'te koydugum "bu kolon yapilara karsi
    kontrol edilmedi" uyarisi KALDIRILDI cunku artik kontrol edildi ve gecti;
    yerine iki kolonun ayni testten gecip farkli sonuc verdigi ve grup 5'teki
    Glu/Cys ikamesi yazildi.

    (g) ACIK: GxtA'nin 36,43 A'si aciklanmadi (6WNC'de komsu zincir B olarak
    bulundu, otekilerde C/D). Ve kopru kolonunun Asp/Glu dagilimi HALA
    yalnizca 27 yapili referans icin yapisal olarak dogrulandi; 76 tipin
    tamamindaki sayim hala yalnizca kolon okumasi.

84. **B ADIMI: ALTI KURASYON DUZELTMESI UYGULANDI (2026-10-28).** Madde 75'te
    kaynaktan dogrulanan uyusmazliklarin hepsi kapatildi.

    (a) `1_108_PobA` -- IKI HATA BIRDEN, ikisi de UniProt Q52185 ile cozuldu.
    Substrat 4-hydroxybenzoate yaziliyordu; UniProt proteini "Phenoxybenzoate
    dioxygenase subunit alpha" diye adlandiriyor, alt adi "4-carboxydiphenyl
    ether" ve yolu "carboxydiphenyl ether degradation" -- 4-hidroksibenzoat
    yolda hic gecmiyor. Urun protocatechuate yaziliyordu; UniProt enzimin
    "yields a hemiacetal that spontaneously hydrolyzes to phenol and
    protocatechuate" dedigini kaydediyor, yani protokatekuat KENDILIGINDEN
    hidrolizden geliyor. **Bu, o hata sinifinin UCUNCU ornegi** (CMO ve HcaE'den
    sonra). reaction_class hydroxylation -> angular_dioxygenation ve family
    aromatic_acids -> biaryls_ethers, 2_204_Dpe ile ayni (ayni acil saldiri).
    SMILES PubChem CID 75182. ACIK: Dpe'nin urunu de "phenol and catechol"
    yaziyor, yani ayni asagi-akis adlandirmasi; resmi denklem bagi yok ve ayni
    incelemeyi hak ediyor.

    (b) `1_103_CarAa` -- `pdb: 1WW9` KALDIRILDI. RCSB'ye karsi olculdu: 1WW9'un
    polimer varligi Janthinobacterium sp. J3'ten ve dizisi buradakiyle
    UYUSMUYOR, `1_102_CARDO` ile uyusuyor. Iki ayri protein ayni yapiya isaret
    ediyordu, yani bu tip sahip OLMADIGI bir yapiyi iddia ediyordu. source_kind
    structure -> paper.

    (c) `2_201_BphA1` -- `pdb: 2XR8` KALDIRILDI ve organizma duzeltildi. 2XR8'in
    iki varligi da Burkholderia xenovorans'tan ve hicbiri bizim diziyle
    uyusmuyor; 2XR8 ayrica LB400 enziminin yonlendirilmis-evrim varyanti.
    Bizim dizi UniProt Q52028 ile birebir: KF707 enzimi (artik Metapseudomonas
    furukawaii), BASKA BIR SUS.

    (d) `2_215_cadA` -- organizma ve kaynak duzeltildi. Kayit "J Bacteriol
    2002;184:509-518 (Bradyrhizobium sp. HW13)" diyordu; dizimiz HW13 CadA
    DEGIL (UniProt Q8VUZ2, ayni uzunlukta ama ~%85 kimlikli), GenBank
    AGU69349.1 ile birebir, yani Sphingomonas sp. ERG5 CadA. Substrat 2,4-D
    DOGRU ve artik kaynakli.

    (e) `4_402_PsbAb` -- SUBSTRAT, URUN VE ORGANIZMA duzeltildi. Kayit substrati
    4-sulfobenzoate yazip PMID:8828208'i (Comamonas testosteroni) gosteriyordu;
    o makale baska bir organizmanin enzimi hakkinda ve bu diziyi hic
    kapsamiyor. Dizi GenBank BAA82116.1 ile birebir (psbAb, Rhodopseudomonas
    palustris) ve onu yatiran makale "Cloning of Genes Participating in
    Aerobic Biodegradation of p-Cumate from Rhodopseudomonas palustris".
    BAA82116.1, P. putida p-kumat dioksijenazinin buyuk alt birimine ~%39-40
    kimlikle hizaliyor. Kimya 4_401_CmtAb'ye gore kaydedildi. UYARI: literatur
    bu kumenin buyuk alt birimini tutarsiz adlandiriyor (bir hizalamada psbAb,
    bir karsilastirma sekilde psbAc); kullanilan etiket bu DIZININ GenBank
    kaydini izliyor (/gene = psbAb). Atama kume baglami ve homolojiye dayaniyor,
    saflastirilmis alt birimin olcumune DEGIL, o yuzden chemistry_evidence
    `curated` kaldi. family sulfoaromatics -> aromatic_acids, yani
    sulfoaromatics artik tek uyeli (TsaM).

    (f) `3_312_NinA` -- substrat daraltildi. "nitroaromatic compounds" ve
    "hydroxylated aromatic and nitrite" o kadar genel ki ne kontrol edilebilir
    ne cizdirilebilir; source alani da "no individual citation recorded"
    diyordu. Dizinin arkasindaki makale 1-nitronaftalin katabolizmasi
    hakkinda; substrat 1-nitronaphthalene (PubChem CID 6849), urun
    1,2-dihydroxynaphthalene ve nitrit. Not: ne GenBank ne RefSeq bu gene ninA
    demiyor; etiket bu veritabaninin kendi etiketi.

    (g) Kapi: PDB sayisi 27 -> 25, hepsi gecerli. Tek FAIL, uc tipin bayat
    `ro_search` girdisi (PobA 170, bfzA 93, PsbAb 28 kayit) -- DOGRU davranis,
    tek rebuild'i bekliyor.

83. **TEMEL: OKSIJEN KADERI ARTIK HESAPLANIYOR, OKUNMUYOR (2026-10-28).**
    Kullanici "temeli saglam atalim" dedi. Dogru teshis: 16 tipi firsatci
    capalarla yamamak temel atmak degildi. Asil sorun kimya alanlarinin
    SERBEST METIN olmasi ve hicbir makine-cozulebilir bagi olmamasiydi -- o
    yuzden oksijen kaderini OKUYARAK cikarmak zorunda kaldim ve iki kez
    yanildim.

    (a) ZINCIR KURULDU (`reaction_balance.py`). Uc resmi kaynak:
    UniProt -> EC numarasi ve Rhea reaksiyon kimligi; Rhea -> DENGELENMIS
    denklem ve katilimcilarin ChEBI kimlikleri (ChEBI listesi denklem
    SIRASINDA geliyor, yani yanlara bolunebiliyor); ChEBI (OLS4) -> her
    katilimcinin molekul formulu. Sonrasi atom sayimi. **Kurator adlandirmasi
    bu zincire HIC girmiyor** -- iki hatamin ikisi de bu yuzden imkansiz
    hale geliyor.

    (b) SONUC: 23 tipin Rhea kimligi var. **15'i hesaplandi ve 15'inin de
    kayitli degeri hesaplananla UYUSTU, sifir celiski.** Yani o 15 deger artik
    bagimsiz olarak dogrulanmis durumda.

    (c) ARAC KENDI SINIRINI OGRENDI -- ve bu, ayni dersin kendi aracimi
    vurmasiydi. Ilk surum AntA ve CbdC icin EMIN BIR YANLIS uretiyordu:
    antranilat C7H6NO2 (O=2) -> katekol C6H6O2 (O=2) + CO2, sayim "ana urunde
    kazanc yok" diyor ve `leaving_group` donduruyordu. Oysa katekolun IKI
    OH'i da O2'den gelir, antranilatin iki oksijeni CO2 olarak cikar. **Toplam
    sayim atom IZLEMEZ.** Artik iskelet parcalandiysa (ana urunun karbon
    sayisi substrattan kucukse) ve substrat oksijen tasiyorsa deger
    `ambiguous` donuyor ve kuratorun gerekceli degeri ayakta kaliyor. Bedeli
    acik: demetilasyonlar da (kafein C8 -> C7) bu sinifa giriyor, yani sekiz
    tipin degeri hesapla degil kayitli gerekceyle duruyor. Arac yalnizca
    GUVENILIR oldugu yerde konusuyor.

    (d) KAPSAMA VERIDE GORUNUR OLDU: her tipe `chemistry_evidence` alani.
      computed_from_official_equation  15
      official_equation_ambiguous       8
      curated                          53
    "Olculdu" ile "kurator yazdi" artik veride ayirt edilebiliyor. Bu, enzim
    hakkinda daha zayif bir iddia DEGIL -- bu veritabaninin su an neyi
    dogrulayabildigi hakkinda bir ifade. Rhea kimligi ve resmi denklem de
    manifestoya yazildi (`rhea`, `official_equation`), yani bag kutukte degil
    VERIDE duruyor.

    (e) UC YENI KONTROL (kapi 102 -> 105): hesaplanabilen her deger
    hesaplananla uyusuyor mu (celiski FAIL), "computed" diyen her tipin
    gercekten hesaplanmis bir degeri var mi (bos guvence olmasin), her tip
    kimyasinin neye dayandigini soyluyor mu. Kontrol BOZULARAK sinandi: bir
    degeri degistirince adiyla FAIL veriyor.

    (f) BENIM IKI HATAM DAHA: karsilastirmayi ilk kez MANIFESTO alanina
    yaptim, oysa o alan yalnizca 14 tipte var ve digerleri degeri sinif
    varsayilanindan aliyor -- "kayitli None" cikan 12 satir veri bulgusu
    degil, benim hatamdi. Cozumlenmis degere (`oxygen_fate.csv`) gecildi.

    (g) SIRADAKI: 53 "curated" tip. Bunlarin bir kismi icin EC numarasi
    `source` metninde duruyor (PAO 1.14.15.17, Nvd 1.14.19.21) ve EC uzerinden
    Rhea'ya baglanabilir; 37'sinde hicbir capa yok. Capa bulunamayan her tip
    `curated` olarak kalmali -- uydurma bir EC, eksik EC'den kotudur.

82. **URUN KIMLIGI ARTIK RESMI DENKLEME KARSI DOGRULANIYOR -- uc gercek hata
    daha bulundu (2026-10-28).** Kullanici "boyle hatalar yapmayalim, hepsi
    icin guncel literature bakalim" dedi. Dogru teshis: sorun tek bir kayit
    degil, YONTEMDI.

    (a) HATANIN SINIFI. Kolin monooksijenazda oksijen kaderini kendi formul
    hesabimdan cikardim. Aritmetik dogruydu ama YANLIS MOLEKULUN formuluydu:
    kaydimiz urunu "betaine aldehyde" yaziyordu, enzimin biraktigi urun ise
    hidrat. Hesap saglam, girdi curuk -- ve "olculdu" diye sunuldugunda en
    tehlikeli hata turu bu.

    (b) YONTEM DEGISTI: `verify_reactions.py`. Her tip icin UniProt'tan EC
    numarasi ve Rhea'nin DENGELENMIS resmi reaksiyon denklemi cekilip kayitli
    `product_en` ile karsilastiriliyor.

    (c) ILK YAZDIGIM KONTROL BOSTU VE ELENDI. Yalnizca ORTAK kelime sayiyordu,
    bu yuzden "betaine aldehyde" ile "betaine aldehyde hydrate" iki kelime
    paylasip "uyusuyor" cikiyordu -- yani kontrol, KENDISINI YAZDIRAN hatayi
    yakalamiyordu. Karsilastirma iki yonlu yapildi: resmi adda OLUP bizde
    olmayan kelimeler de bildiriliyor. Simulasyonla dogrulandi, eski kayit
    icin "INCELE -> hydrate" diyor.

    (d) UC GERCEK HATA BULUNDU, hepsi ayni sinif -- enzimin BIRAKTIGI molekul
    yerine baska bir molekulu urun yazmak:
    * `2_217_HcaE` urunu "3-(2,3-dihydroxyphenyl)propionate" yaziyordu, yani
      ASAGI AKISTAKI dehidrogenazin urettigi yeniden-aromatik katekol. Resmi
      denklem (EC 1.14.12.19) urunu
      "3-(cis-5,6-dihydroxycyclohexa-1,3-dien-1-yl)propanoate" veriyor:
      aromatik OLMAYAN cis-dihidrodiol. CMO hatasinin tam esi.
    * `1_107_CndA` "the corresponding aldehyde" diyordu; resmi denklem (EC
      1.14.15.23) butachlor'dan BUTIL FORMAT, alachlor'dan metil format,
      acetochlor'dan etil format veriyor -- kopan parca bir ESTER, aldehit
      degil.
    * `5_506_Stc2` yalnizca "N-methylproline" diyordu; resmi denklem (EC
      1.14.13.247) FORMALDEHITI de veriyor. Bu onemli, cunku oksijenin
      bittigi yer orasi ve `oxygen_fate: leaving_group` degeri urun alani
      tarafindan desteklenmiyordu.
    Ucu de duzeltildi, gerekcesine resmi denklem BIREBIR yazildi.

    (e) KAPSAMA SINIRI DURUSTCE: 76 tipin **23'u** kontrol EDILEBILDI
    (14 temiz, 9 fark). **53'u edilemedi**: 27'sinde UniProt kaydinda
    katalitik aktivite yok, 26'sinda UniProt accession'i yok. Yani bu tur bir
    temizlik getirdi ama TAM KAPSAMA getirmedi; kalan 53 tipin urun iddiasi
    hala yalnizca bizim kurasyonumuza dayaniyor ve bu GORUNUR olmali.

    (f) ADLANDIRMA FARKLARI KUTUGE GECTI: `reaction_name_exceptions.yaml`,
    dokuz kayit. Ayni molekulun geleneksel adiyla IUPAC adi ("cis-benzene
    dihydrodiol" ile "cis-1,2-dihydrobenzene-1,2-diol") uyusmaz sayiliyordu;
    bunlari SESSIZCE gecmek yanlis olurdu, cunku bir gun gercek bir kimya
    farki ayni kutuya duser. Her kayit resmi adi, bizim adi, neden ayni
    molekul oldugunu ve okundugu tarihi tasiyor.

    (g) DORT YENI KONTROL (kapi 98 -> 102): bildirilen her fark kayitli mi,
    olu kayit var mi, her kayit dort alani tasiyor mu, ve RAPOR types/'TAN
    TAZE mi (bayat raporla kapi eski veriye bakar). Kontrol BOZULARAK sinandi:
    bir kayit cikarilinca adiyla FAIL veriyor, geri konunca geciyor.

    (h) 20 tipe `ec` alani eklendi -- resmi kayda gelecekte baglanmanin
    capasi. EC numarasi UniProt'tan geldi, elle yazilmadi.

    (i) KURAL: bundan sonra bir oksijen kaderi ya da urun karari verilirken
    ILK IS EC/Rhea denklemine bakmak. Kendi formul hesabi ancak urunun kimligi
    resmi kayittan dogrulandiktan SONRA yapilir.

81. **CMO: IKI KEZ YANLIS KAYDEDILDI, GUNCEL LITERATUR COZDU (2026-10-28).**
    Kullanici once "CmoS icin Rieske diyorlar ama", sonra "CmoS'ta yukseltgen
    degil mi", sonunda "guncel literature bi baksa ya" dedi. Ucuncusu dogru
    cevabi verdirdi; ilk ikisi benim iki ayri hatami buldurdu.

    (a) HATANIN KOKU TEK: denge, enzimin BIRAKTIGI urunden degil, onun
    DEHIDRATLANMIS DENGE FORMUNDAN hesaplandi. `product_en` alani "betaine
    aldehyde" yaziyordu; serbest aldehit C5H12NO+ ve kolinle (C5H14NO+)
    karsilastirilinca -2H / sifir O cikiyor, yani "oksijen kazanci yok".
    Bu yuzden deger once `transient` ("sonra su olarak cikar"), sonra
    `no_net_gain` oldu. IKISI DE YANLIS.

    (b) BRENDA (EC 1.14.15.7) reaksiyonu birebir soyle veriyor: "choline + O2
    + 2 reduced ferredoxin + 2 H+ = betaine aldehyde **HYDRATE** + H2O + 2
    oxidized ferredoxin". Urun HIDRAT, yani gem-diolun kendisi. Dengelenince
    kolin C5H14NO+ + O2 -> hidrat C5H14NO2+ + H2O: **bir oksijen URUNE, bir
    oksijen suya** -- sira bir monooksijenaz stokiyometrisi. EC 1.14.15 sinif
    basligi da "incorporation of one atom of oxygen into the other donor".
    Serbest aldehit, BADH'nin tukettigi KENDILIGINDEN denge formudur;
    dehidratasyon bu enzimin katalizledigi bir adim DEGIL.

    (c) DUZELTMELER: uc CMO tipinde `oxygen_fate` -> `product` (BRENDA
    gerekcesiyle), `product_en` -> "betaine aldehyde hydrate (in equilibrium
    with betaine aldehyde)". Dayanagi curuyen iki sey de GERI ALINDI:
    `no_net_gain` degeri sozlukten ve koddan kaldirildi (uyesi kalmadi;
    uyesi olmayan ve gerekcesi curuyen bir kategori kullanilmaktan cok yanlis
    kullanilir) ve `hydroxylation`'a acilan `overridable` izni kaldirildi.
    Site metni, tablo etiketi ve sablon dali da temizlendi; eski iki degerden
    hicbir iz kalmadi (dogrulandi).

    (d) YENI SUTUN: `o2_role`, oksijen kaderinden TURETILIR, ayri kurasyon
    gerektirmez. donor 72 | oxidant_only 3 | unknown 1. Kullanicinin
    "hangileri yukseltgen onu da belirtelim" istegi: okuyucu artik iki degeri
    kafasinda birlestirmek zorunda degil, `/clusters` tablosunda "O2 role"
    facet'iyle dogrudan suzuyor. Tablo 18 sutun, basliklar ve hucreler 18/18
    hizali, "oxidant only" 3 satir.

    (e) YALNIZCA YUKSELTGEN OLAN UC TIP: RedG (C-C bagi kurma), Nvd ve DAF36
    (desaturasyon). Ucu de biyosentetik. CMO artik bu kumede DEGIL.

    (f) DERS, madde 80'de yazilanin duzeltilmis hali: formul dengesi ile atom
    izlemesi ayri seylerdir -- AMA once hangi molekulun URUN oldugu
    dogrulanmalidir. Yanlis urunle yapilan denge iki kez yanlis sonuc verdi ve
    ikisi de "olculdu" gibi gorundu. Bir sonraki oksijen kaderi karari
    verilirken ilk is EC kaydindaki reaksiyon denklemine bakmak olmali.

80. **CMO'NUN OKSIJEN KADERI FAZLA IDDIAYDI -- degeri ve adi duzeltildi
    (2026-10-28).** Kullanici "CmoS'ta yukseltgen degil mi" diye ikinci kez
    sorunca kendi degerimi sinadim ve HATALI buldum.

    (a) NE YANLISTI. Deger `transient` adini tasiyordu ve gerekcesinde "urun
    o oksijeni hic tutmuyor, su olarak cikiyor" yaziyordu. Bu bir ATOM
    iddiasidir ve dayanagi yoktu: formul dengesi oksijenin SAYISINI soyler,
    hangi ATOM oldugunu SOYLEMEZ.

    (b) KIMYA NEDEN BELIRSIZ BIRAKIYOR. Kolinde hidroksillenen karbon ZATEN
    bir OH tasiyor. Oraya oksijen sokulunca ayni karbonda iki OH olan bir
    gem-diol cikar; bu iki hidroksil kimyasal olarak ESDEGERDIR ve su kaybinda
    hangisinin gittigi belirsizdir. O2'den gelen atom uronde kalmis da
    olabilir. Bunu yalnizca 18O2 isaretlemesi cozer.

    (c) KANIT ARANDI, YOK. Elimizdeki uc kaynagin hicbirinde gem-diol ya da
    18O isaretlemesi GECMIYOR (derleme ve tezde sifir kayit). Bizim atfimiz
    (rathinasabapathi1997_CmoS) yalnizca "ferredoxin-dependent conversion of
    choline to betaine aldehyde" gosteriyor -- mekanizmayi kurmuyor.

    (d) DUZELTME. Deger `no_net_gain` olarak yeniden adlandirildi ve artik
    YALNIZCA DESTEKLENEN seyi iddia ediyor: urun oksijen KAZANMIYOR (kolin
    C5H14NO+ -> betain aldehit C5H12NO+, -2H ve sifir O). Site metni de
    duzeltildi: "An oxygen atom is genuinely inserted and this enzyme is a
    true Rieske-type monooxygenase, but the product gains no oxygen overall.
    Which atom ends up in the product is NOT established... Only labelling
    with 18O2 could resolve it."

    (e) AYRIM NEREDE OLDUGU NETLESTI: SUBSTRAT KARBONUNDA. CAO ve PTC52 bir
    METILI hidroksilleyip formile goturuyor; karbon oksijensiz basladigi icin
    urun oksijen KAZANIYOR ve `product` DOGRU. CMO'da karbon zaten OH
    tasidigi icin net kazanc sifir. Yani `hydroxylation` sinifinin iki ayri
    kaderi tasimasi tesadufi degil, karbonun baslangic durumuna bagli --
    `overridable` gerekcesi buna gore yeniden yazildi.

    (f) Kullanicinin ilk sorusu (madde 79'daki "CmoS icin Rieske diyorlar
    ama") degeri YANLIS BULDURMUSTU; ikinci sorusu GEREKCEYI yanlis buldurdu.
    Ikisi de ayni dersi veriyor: formul dengesi ile atom izlemesi ayri
    seylerdir ve biri otekinin yerine kullanilamaz.

79. **"PROTEIN DIZISI YETER MI" SORUSUNUN CEVABI: KISMEN -- ve sinir yaziya
    gecirildi (2026-10-28).** Kullanici yeni bir RO bulan birinin bizim
    semayla oksijen kaderini, tipi ve rolu cikarabilecegini sordu.

    (a) DIZI SU SORUYA YETER: "bu bir Rieske oksijenaz mi". Cunku bu bir
    KOFAKTOR iddiasidir ve kofaktor ligandlari dizide okunabilir: Rieske 4/4
    (`C-x-H...C-x-x-H`) + mononukleer histidinler. Kapi tam bunu yapiyor.

    (b) DIZI SU IKISINE YETMEZ: `oxygen_fate` ve `metabolic_role`. Ikisi de
    REAKSIYON hakkinda olgudur; substrat ve URUN bilinmeden belirlenemez.
    Bizim degerlerimiz 76 REFERANS TIPI icin kuratorlu -- substrati ve urunu
    bilinen enzimler icin -- ve 19.114 uye icin OLCULMUS DEGILDIR.

    (c) CIKARIM NE KADAR TASIR: OLCULDU. En yakin karakterize referansa
    kimlik, tier'lara gore:
      characterized   369 uye  ort. %98,8
      close_homolog  4.601 uye ort. %74,9
      family_member  3.480 uye ort. %49,0
      distant       10.376 uye ort. %32,2
      novel            288 uye ort. %24,1
    **19.114 uyenin 10.664'u (%55,8) herhangi bir karakterize referansa
    %40'in ALTINDA benziyor.** %32 kimlikte substrat cikarilamaz, oksijenin
    nereye gittigi hic cikarilamaz. Yani veritabaninin yarisindan fazlasi icin
    dizi bu iki eksen hakkinda HICBIR SEY soylemiyor.

    (d) BU YUZDEN TASARIM KURALI: iki eksen YALNIZCA TIP SEVIYESINDE
    gosterilir. Dogrulandi: `ro.html` (uye sayfasi) ve `leaf.html` (varyant
    sayfasi) ikisini de HIC tasimiyor; yalnizca `_cluster_table.html` (tip
    listesi) ve `cluster.html` (tip sayfasi) tasiyor. Bir uye icin
    gosterilecek olursa "CIKARIM" olarak ve en yakin referansa kimlikle
    birlikte etiketlenmek ZORUNDADIR; aksi halde olculmemis bir sey iddia
    edilmis olur.

    (e) LITERATUR, bu tiplerin aileye ait oldugunu DOGRULUYOR (kullanicinin
    "literature gore hangileri olmali" sorusu). Miao & Schmidt 2025
    (`miao2025_perspective`, elimizdeki PDF) ailenin repertuvarini sayiyor:
    "dioxygenation, monooxygenation, oxidative N- and O-dealkylation,
    desaturation, sulfoxidation, C-C bond formation, N-oxygenation, and C-N
    bond cleavage". Sekil 1 enzimleri ADIYLA aniyor: DAF-36/Neverland (Nvd)
    -- metinde birebir "catalyzed by Rieske oxygenase DAF-36/Neverland" --,
    RedG (undesilprodigiosin karbosiklizasyonu), PrnD (N-oksijenasyon), CntA
    (C-N kirilmasi), DMO (bizim DdmC) ve NdmA. Yani Nvd, DAF36, RedG, PrnD ve
    CntA bizim genislemelerimiz DEGIL, derlemenin kendi kapsami. Madde 76'daki
    kabul karari literaturle uyumlu.
    NOT: derlemede sulfoksidasyon NDO'nun yan aktivitesi olarak geciyor ve
    bizim `reaction_class` sozlugunde KARSILIGI YOK; sema tip basina TEK
    reaksiyon tutuyor, yani yan aktiviteler kaydedilemiyor. Bu, semanin
    bilinen bir siniri.

    (f) KUCUK TUTARSIZLIK: `citations.yaml` icindeki `local_pdf: articles/...`
    yollari depo kokunden COZULMUYOR -- PDF'ler bir ust dizinde
    (`/home/lin-bio/RieskeDB/articles/`). Yol ya duzeltilmeli ya da neye gore
    oldugu yazilmali.

78. **"BIYODEGRADASYON YAPMAYANLARI ELEME" ISTEGI: OLCUT DEGISTIRILDI, HICBIR
    TIP SILINMEDI (2026-10-28).** Kullanici oksijen aktarmayan uyeleri suzmek
    istedi. Istek dogru, ilk onerilen olcut YANLIS eksendi.

    (a) OKSIJEN KADERI, BIYODEGRADASYON ICIN YANLIS EKSEN. Iki yonlu kanit:
    SxtT, GxtA, PrnD ve CAO oksijeni URUNE sokar ama biyodegradasyon yapmaz
    (saksitoksin tailoring, pirrolnitrin biyosentezi, klorofil b sentezi) --
    oksijen kaderine gore suzulse KALIRLARDI. Tersine vanillat, guaiakol,
    dikamba, kafein ve 2,4-D demetilasyonlarinda oksijen kopan formaldehite
    gider ama bunlar klasik biyodegradasyondur -- suzulse 13 TIP SILINIRDI.
    Iki eksen DIKTIR.

    (b) IKINCI EKSEN EKLENDI: `metabolic_role` (catabolic / biosynthetic /
    unknown), `metabolic_roles.yaml` + `metabolic_role.py` + turetilmis
    `metabolic_role.csv`, rebuild grafiginde (51 artefakt). Varsayilan kural
    yalnizca rolun TEK ANLAMLI oldugu substrat siniflari icin var (xenobiotic
    ve natural_aromatic -> catabolic); `natural_specialized` icin kural YOK ve
    22 tipin hepsi kendi rolunu gerekcesiyle beyan etti.
    Sonuc: **catabolic 64, biosynthetic 11, unknown 1 (AerC).**
    Ogretici karsitlik: CmoS/CmoB/CmoA kolini betain aldehide YAPAR (bitkide
    ozmoprotektan sentezi, biyosentetik), BmoA/GbcA glisin betaini YIKAR
    (katabolik) -- ayni kimyanin iki yonu. Ayni sekilde grup 6'da PAO klorofili
    YIKAR (katabolik) ama CAO ve PTC52 klorofil b YAPAR (biyosentetik).

    (c) HICBIR TIP SILINMEDI, cunku silmek aile taramasini yanlis yapar.
    Ikisi de `data-facet` ile FILTRELENEBILIR sutun: "O2 ends up" ve "Role".
    Biyodegradasyon odakli okur Role=catabolic ile suzer.

    (d) KULLANICININ "CmoS icin Rieske diyorlar ama" SORUSU GERCEK BIR HATA
    BULDURDU. CmoS'un oksijen kaderi `hydroxylation` sinif varsayilanindan
    `product` geliyordu ve bu YANLISTI: kolin C5H14NO+ -> betain aldehit
    C5H12NO+, yani -2H ve SIFIR O (PubChem CID 305 -> 249). CMO gercekten
    hidroksilliyor ama gem-diol ara urunu su kaybediyor, oksijen uronde
    KALMIYOR. Dorduncu deger eklendi: **transient** -- "oksijen sokulur, sonra
    su olarak cikar". Sinif varsayilani product kaldi ama `hydroxylation`
    artik `overridable: true` ve ezme YALNIZCA sozluk izin verdiginde VE tip
    gerekcesini yazdiginda gecerli; aksi halde celiski olarak FAIL.
    CmoS Rieske olmaktan cikmiyor: motif kapisini Rieske 4/4 ile geciyor ve
    `metabolic_role` enzimin Rieske olup olmadigi hakkinda HICBIR SEY soylemez.

    (e) IKI SESSIZ HATA DAHA BULUNDU VE DUZELTILDI, ikisi de "deger eklendi,
    tuketici guncellenmedi" sinifindan:
    * `transient` cozucuye eklendi ama app.py'deki etiket haritasina
      eklenmedi: site degeri olan tipler icin "unknown" gosterdi.
    * Tip sayfasindaki dal sirasi yanlisti: `transfers` hem oxidant_only hem
      transient icin False oldugu icin CmoS sayfasinda "mechanistically an
      oxidase rather than an oxygenase" YAZDI -- oysa CMO bir monooksijenazdir.
      transient artik en basta eleniyor. RedG'de dogru metin cikiyor.

    (f) DORT YENI KONTROL (kapi 94 -> 98): her tip bir role cozuluyor mu,
    degerler sozlukten mi, yayimlanan gerekce metninde Turkce var mi, tablo
    taze mi. Hepsi gecti.

    (g) HOCAYA VERILECEK GEREKCE. Hoca PAO/RedG icin "oksidatif bir
    biyotransformasyon, kesin yukseltgen" dedi. Sonuc dogru, gerekce eksik: o
    olcut P450'leri, flavin monooksijenazlari ve bakir oksidazlari da iceri
    alir. Savunulabilir gerekce KOFAKTOR CIFTIDIR -- sekiz adayin hepsi motif
    kapisini NahAc ve PDO ile birebir ayni profille geciyor (madde 76). Yaninda
    iki kayit: RedG/Nvd/DAF36'da urune oksijen girmiyor (formuller), ve
    katalitik triad'in ucuncu ayagi yayimlanmiyor cunku karboksilat kolonu
    26/26 yapida yanlis (madde 73).

77. **OKSIJENIN KADERI ARTIK BIR ALAN VE BIR FILTRE (2026-10-28).** Kullanici,
    urune oksijen sokmayan uyeler bulununca "bunu filtreye eklemek gerekir,
    biz oksijen aktaranlara odakleniyoruz" dedi. Dogru teshis: "dioksijenaz"
    diye suzen biri desaturazlari da aliyordu.

    (a) DEGER REAKSIYON SINIFINDAN TURETILEMEZ. Ilk tasarim sinif -> kader
    eslemesiydi ve YANLIS olacakti: `dioxygenation_with_release` homojen degil.
    Nitroarenlerde ve anilinlerde her iki oksijen katekole giriyor (nitrobenzen
    C6H5NO2 -> katekol C6H6O2: katekolun iki oksijeni de yeni, eskiler
    nitritle ayriliyor; anilin C6H7N'de hic oksijen yok, yani ikisi de yeni
    olmak ZORUNDA). Fenoksiasetatlarda ise eter kiriliyor ve oksijen
    GLIOKSILATA gidiyor: 2,4-D C8H6Cl2O3 -> 2,4-diklorofenol C6H4Cl2O +
    glioksilat C2H2O3; fenolun tek oksijeni eski eter oksijenidir. Ayni sinif
    iki ayri kader tasiyor.

    (b) UC DEGER DEGIL, DORT. Ikili bir "oksijen giriyor mu" alani veriyi
    yanlis temsil ederdi, cunku demetilasyonlarda oksijen URUNE DEGIL kopan
    formaldehite giriyor ve bu 13 tip ediyor:
      product 59 | leaving_group 13 | oxidant_only 3 | unknown 1
    Her deger formul aritmetigiyle sinandi, formuller PubChem'den okundu:
    kafein C8H10N4O2 -> teobromin C7H8N4O2 + CH2O (teobrominin oksijen sayisi
    kafeinle AYNI), guaiacol C7H8O2 -> katekol C6H6O2 + CH2O, vanillat
    C8H8O4 -> protokatekuat C7H6O4 + CH2O. unknown olan AerC: tahmin yazilmadi.

    (c) TEK KAYNAK: `reaction_classes.yaml`. Sinif adlari uc yerde yasiyordu
    (validate_curation.py ve add_reference.py'de iki ayri VALID_REACTION_CLASS
    kumesi, ve tiplerin alanlari). Artik sozluk dosyada, VALID_REACTION_CLASS
    onu okuyor, dosya yoksa eski sabit liste yedek. Homojen OLMAYAN siniflar
    `per_type` isaretli ve tipin kendi `oxygen_fate` + `oxygen_fate_basis`
    beyanini ZORUNLU kiliyor; 14 tip beyan etti.

    (d) chemistry.csv'ye SUTUN EKLENMEDI. O dosyayi yirmiden fazla betik
    okuyor ve sema degisikligi gereksiz risk. Deger `oxygen_fate.py` ile
    types/ + sozlukten TURETILIP `oxygen_fate.csv`'ye yaziliyor; adim
    rebuild.py grafigine girdi (50 artefakt) cunku bir tipin sinifi degisip
    tablo guncellenmezse site eski kaderi filtreler ve kimse farketmez.

    (e) ALTI YENI KONTROL, hepsi gecti (kapi 88 -> 94 kontrol): her sinif
    sozlukte tanimli mi, her per_type uyesi beyan etmis mi, tip sinifi
    sessizce eziyor mu, degerler sozlukten mi, yayimlanan gerekce metninde
    Turkce var mi, turetilmis tablo taze mi. Kontroller BOZULARAK sinandi:
    alti saldirgan senaryonun hicbiri sessizce `product` varsaymiyor; gecersiz
    deger bile `unknown` donuyor.

    (f) DIL TUZAGINA DUSULDU VE YAKALANDI. Ilk yazimda sinif gerekceleri
    Turkce yazildi ve turetilmis tabloya Turkce metin SIZDI -- bu tablo siteye
    besleme yapiyor. `evidence` alani artik Ingilizce (yayimlanir),
    `why_per_type` Turkce (yalnizca kurator). Sizintiyi kalici olarak kapatmak
    icin Turkce kontrolu bu sutuna da baglandi.

    (g) SITE. `/clusters` tablosuna "O2 ends up" sutunu eklendi ve
    `data-facet` ile FILTRELENEBILIR; gerekce tooltip olarak geliyor. Tip
    sayfasinda reaksiyon sekilinin altinda bir kutu var ve oksijen
    girmiyorsa sunu YAZIYOR: "O2 is the oxidant here, not an oxygen donor, so
    this enzyme is mechanistically an oxidase rather than an oxygenase... the
    name is a statement of family membership, not of mechanism." 76 tip
    sayfasinin 76'sinda kutu goruldu; basliklar ve hucreler 16/16 hizali.

    (h) ACIK: rebuild henuz kosulmadi (bfzA icin ro_search bayat) ve site
    yayina ALINMADI.

76. **SEKIZ ADAY REFERANS FIILEN SETTEYDI AMA ETIKETI "ONAYLANMADI" DIYORDU --
    kofaktor olcutuyle onaylandi (2026-10-28).** Kullanici PAO ve RedG'yi RO
    saymanin sorun olup olmadigini sordu ve gerekcesini "o da bir
    biyotransformasyon, kesin yukseltgen" diye verdi.

    (a) O GEREKCE TEK BASINA YETMEZ, cunku fazla gecirgen: P450'ler, flavin
    monooksijenazlar ve bakir oksidazlar da kesin yukseltgen ve
    biyotransformasyon yapiyor. Isi yapan olcut KOFAKTOR CIFTIDIR: Rieske
    [2Fe-2S] + mononukleer demir.

    (b) OLCULDU, sekizi de kapiyi geciyor (`ro_motif.motif_report`, referans
    hizalamasi uzerinde): GdmA, PrnD, RedG, Nvd, DAF36, PAO, CAO, PTC52 --
    hepsi Rieske 4/4, her iki mononukleer histidin yerinde, kopru Asp var.
    NahAc ve PDO ile BIREBIR ayni profil; AerC'nin aksine istisna kaydi
    GEREKMIYOR. Elimizdeki Miao & Schmidt 2025 derlemesi de ailenin reaksiyon
    repertuvarina desaturasyon ve C-C bagi kurmayi dahil ediyor.

    (c) ETIKET BAYATTI. Sekizinin `added_by` alani "candidate manifest, NOT
    admitted; awaiting curator review" diyordu, yani yayimlanan sette
    duruyorlar ama etiketleri onaylanmadiklarini soyluyordu. Bu metin siteye
    CIKMIYOR (sablonlarda `added_by` kullanilmiyor), yani yayimlanan bir
    celiski degildi; yine de kutuk yanlisti. Etiket, kanitiyla birlikte
    duzeltildi.

    (d) IKI UYARI KAYDEDILDI.
    * "Katalitik 3/3" ASLINDA 2/2. Ucuncu ayak karboksilat ve madde 73'te
      yanlis oldugu olculen kolondan okunuyor: 26 kristal yapinin 26'sinda
      kolonun cagirdigi kalinti demirden 10,9-14,7 A uzakta. His kolonlari
      (289/294) kristallere karsi DOGRULANDI, karboksilat kolonu
      dogrulanMAdi. Kabul karari bundan etkilenmiyor -- Rieske 4/4 ve iki His
      zaten belirleyici -- ama "triad saglam" diye YAYIMLAMAK, kolon 431
      duzeltilene kadar dogru olmaz.
    * RedG, Nvd ve DAF36'da URUNE OKSIJEN GIRMIYOR. Formuller PubChem'den
      okundu: undesilprodigiosin C25H35N3O -> streptorubin B C25H33N3O ve
      kolesterol C27H46O -> 7-dehidrokolesterol C27H44O. Ikisinde de -2H,
      sifir O. O2 burada yukseltgen, oksijen VERICISI degil; bu enzimler
      mekanistik olarak oksidaz/desaturaz. "Oksijenaz" demek bir AILE UYELIGI
      ifadesidir ve oksijen transferi ima edecek sekilde okunmamalidir. Uc
      tipin `notes` alanina bu, formulleriyle yazildi.

    (e) KALAN IS: alti adayin (GdmA, PrnD, Nvd, PAO, CAO, PTC52) yapisal
    `accession` alani yok -- 67'lik atif setinde degillerdi. PAO'nun
    `substrate_smiles` alani da kasten bos (buyuk klorin makrohalkasi icin
    dogrulanmis SMILES yok; yanlis SMILES eksikten kotudur).

75. **ELLE HAZIRLANAN KIMYA/ORGANIZMA ALANLARININ DENETIMI -- 39 uyusmazlik,
    6'si MAJOR; birincisi duzeltildi (2026-10-28).** Kullanici referans
    tiplerini elle hazirladigini ve hata olabilecegini soyledi; atif taramasi
    39 uyusmazlik bildirdi. Bir ajanin "uyusmuyor" demesi DUZELTME GEREKCESI
    DEGILDIR, o yuzden `verify_curation_claims.py` yazildi: iddialari UniProt,
    RCSB ve NCBI'ya karsi KENDISI olcer, duzeltmeyi yapmaz.

    KAYNAKTAN DOGRULANAN TABLO (11 tartismali tip + 6 GenBank kaydi; her
    kaydin dizisi bizdekiyle birebir, yani accession'lar DOGRU -- hata
    kimya/organizma alanlarinda):
    * `2_211_bfzA` -- kayitli substrat "benzene and fluorene", SMILES fluoren.
      GenBank /product = "alpha subunit of buprofezin dioxygenase", /gene =
      bfzA1. YANLIS. **DUZELTILDI** (asagida).
    * `1_108_PobA` -- kayitli substrat 4-hydroxybenzoate; UniProt Q52185'in
      protein adi "Phenoxybenzoate dioxygenase subunit alpha". 4-hidroksibenzoat
      eter kirilmasinin URUN tarafi, substrat degil. YANLIS.
    * `1_103_CarAa` -- `pdb: 1WW9` yanlis: 1WW9'un varligi Janthinobacterium
      sp. J3 ve dizisi CarAa ile UYUSMUYOR; `1_102_CARDO` ile uyusuyor. Iki
      ayri protein ayni PDB kaydina isaret ediyordu.
    * `2_201_BphA1` -- `pdb: 2XR8` yanlis: 2XR8 varliklari B. xenovorans ve
      dizimizle uyusmuyor. Dizimiz Q52028 = Metapseudomonas furukawaii (KF707).
    * `2_215_cadA` -- substrat 2,4-D **DOGRU** (GenBank: "large subunit of
      2,4-D oxygenase CadA"), ama organizma/makale yanlis: Sphingomonas sp.
      ERG5, Bradyrhizobium HW13 degil. **Bu, "Acik sorular" bolumundeki
      "cadA'nin substrati 2,4-D mi 2,4,5-T mi" sorusunu KAPATIR**: dizinin
      geldigi kayit 2,4-D diyor.
    * `4_402_PsbAb` -- dizi Rhodopseudomonas palustris psbAb; kayitli organizma
      (C. testosteroni) ve substrat (4-sulfobenzoate) dizinin arkasindaki
      kayitla desteklenmiyor.
    * `3_312_NinA` -- substrat "nitroaromatic compounds" fazla genel; makale
      1-nitronaftalin hakkinda. Gen adi hicbir kaynakta ninA degil.
    * `3_313_KshA15` -- "15" ekinin hicbir kaynakta karsiligi yok; accession
      P71875 = M. tuberculosis H37Rv kshA. Yalnizca ETIKET sorunu; degistirmek
      URL'leri yeniden numaralandirir.
    * `5_507_BmoA` **DOGRU** (RefSeq: "glycine betaine monooxygenase").

    (a) bfzA DUZELTILDI. substrate_en -> buprofezin, SMILES PubChem CID
    50367'den (hatirlanarak degil; yanlis SMILES eksik SMILES'ten kotudur),
    product_en -> benzen halkasi dihidroksillenmis cis-dihidrodiol, family
    pah -> anilines, source_kind general -> paper. Reaksiyon sinifi DEGISMEDI:
    makalenin ozeti "degradation is initiated by the dihydroxylation of the
    benzene ring" diyor, yani cis_dihydroxylation dogruydu.
    AILE KURATOR ONAYI BEKLIYOR: `pah` kesin yanlisti (tek halka, PAH degil);
    `anilines` secildi cunku oksitlenen halka N-substitue bir benzen, ama buradaki
    azot amid benzeri, amin degil -- pestisit heterohalkalari icin yeni bir
    sozluk girdisi tercih edilebilir.

    (b) `types_build.py` surukleme kapisi ise yaradi: chemistry.csv'yi
    farkliyken yazmayi REDDETTI ve farkin yalnizca bfzA'nin 7 alani oldugu
    dogrulandiktan sonra `--allow-drift` ile yazildi. 37 degisen satirin
    35'i yalnizca `source` (daha onceki kaynak dogrulamasinin metni), biri
    AerC'nin `notes`'u, biri gerçek duzeltme.

    (c) Dogrulama kapisindaki FAIL kapatildi: varyant kimligi isinden kalan
    dort metin sutunu (`leaf.variant_id`, `leaf.nearest_type`,
    `leaf.label_kind`, `ro_leaf.variant_id`) dil kontrolu siniflandirmasinda
    yoktu. Ucu de kimlik olarak kaydedildi.

    (d) ACIK: `ro_search` indeksi hala eski bfzA degerlerini tasiyor (93 kayit)
    ve kapi bunu FAIL olarak gosteriyor -- DOGRU davranis. Rebuild butun
    duzeltmeler bitince TEK kosuda yapilacak; o zamana kadar site yayina
    ALINMAYACAK.

    (e) SIRADAKI: PobA substrati, CarAa ve BphA1'in PDB alanlari, cadA
    organizmasi, PsbAb organizma+substrat, NinA substrati. Sonra karboksilat
    kolonu 471 -> 431.

74. **ATIFLAR KUTUGE TASINDI -- 67 referans tipi artik yapisal olarak atifli
    (2026-10-28).** Bir veritabaninda yanlis atif, atifsizdan kotudur; o yuzden
    hicbir kunye el ile yazilmadi ve hicbiri tek bir kaynaktan kabul edilmedi.

    (a) DOGRULAMA, EKLEMEDEN ONCE. 67 birincil DOI Crossref'te **iki kez
    bagimsiz** cozuldu (ajanin taramasi, sonra benim `recheck_dois` kosum) ve
    her seferinde baslik/yil kutukle karsilastirildi: 67/67 uyustu, 0 sorun.
    44 ikincil DOI de ayrica cozuldu, 44/44 baslik birebir. Her accession'in
    dizisi UniParc CRC64 ile bizdeki diziye karsi tutuldu: **67/67 birebir.**
    Ayri olarak `verify_reference_sources.py --all` kayitli 35 "birebir"
    iddiasinin 35'ini canli olarak yeniden dogruladi.

    (b) KUTUK 14'TEN 118 ATIFA CIKTI (67 birincil + 37 yeni ikincil). Kimlikler
    ve DOI'ler tekil; hicbir kayit yazarsiz degil -- ikincillerin yazar alani
    da Crossref'ten cekildi.

    (c) TIPLER ARTIK KUTUGE YAPISAL BAGLI. Atiflar bugune kadar `source`
    icinde serbest metindi: ayni makale iki yerde yazilabiliyor ve biri
    guncellenirse ikisi ayrisiyordu. Her tipte artik `citations:` (anahtar
    listesi), `accession:`, `accession_db:`, `accession_verified:`,
    `accession_evidence:` ve `accession_organism:` var. `accession_organism`
    kasten boyle adlandirildi: bu, dizimizle birebir ayni olan KAYDIN
    organizmasi, kuratorun enzim hakkindaki iddiasi degil -- organizma
    duzeltmeleri (madde 75) ayri bir is.

    (d) DENETIM MODULER HALE GETIRILDI. "Cozulebilir tanimlayici" kontrolu
    yalnizca `source` prozasini tariyordu; 67 tipe yapisal atif eklendikten
    SONRA bile 13 uyari veriyordu -- yani tanimlayici vardi, kontrol onu
    okumuyordu. Artik `citations`/`accession` alanlarini da okuyor ve
    `citations.yaml`'da KARSILIGI OLMAYAN bir anahtar HATA sayiliyor.
    Uyari 13 -> 0, hata 0, chemistry.csv gidis-donusu birebir.

    (e) BU SIRADA GERCEK BIR HATA YAKALANDI. `load_manifest` her degeri
    `str()` ile duzlestiriyordu, yani `citations:` listesi `"['a', 'b']"`
    metnine donuyordu; yeni kontrol o metni harf harf gezip her karakteri bir
    atif anahtari saydi ve 2.584 hata verdi. Listeler artik korunuyor
    (PyYAML'siz yerlesik ayristirici da liste okuyacak sekilde guncellendi ve
    PyYAML gizlenerek test edildi). Bu alan 2026-10-14'ten beri alti tipte
    duruyordu ve HICBIR SEY onu okumadigi icin bozuk oldugu gorulmemisti.

    (f) KALAN IS: dort tipin (`1_116_GdmA`, `1_117_PrnD`, `3_307_NagGH`,
    `3_321_Nvd`) DOI'si hala yalnizca `source` prozasinda, kutukte degil.
    Ayrica atiflar henuz tip SAYFALARINDA GORUNMUYOR -- `citations.yaml`'i
    okuyan tek sey `citations_check.py`. Kutuk yayina baglanmadikca yarisi
    yapilmis is.

69. **BAGIMLILIK GRAFIGINDE IKI HATA DAHA -- ikisi de DUZELTILDI (2026-10-22).**
    Ikisi de ayni sinif: bir adim bir dosyayi OKUYOR ama grafikte girdi olarak
    beyan edilmemis. Bu oturumda bu sinifin DORDUNCU ve BESINCI ornegi.

    (a) annotate'in girdilerinde HMM YOKTU. annotate_ro.py iki HMM kullaniyor
    (RieskeDB71.hmm atama, ROmotif71.hmm motif) ama girdi listesi yalnizca
    ["annotate_ro.py", DB] idi. AerC kutuphaneye eklendiginde protein_filter
    bayatladi, annotate BAYATLAMADI, zincir dogrudan evidence_tiers'a atladi ve
    atama 8 gun onceki cand_dom.out ile kaldi: yeni tip SIFIR uye topladi.
    Duzeltme: iki HMM girdiye, cand_dom.out ve cand_aln.sto cikti listesine
    eklendi. Ayrica bir TARAMA yazildi ("HMM kullanip girdide beyan etmeyen
    adim var mi"): alti supheli cikti, ucu gercek bagimlilikti (motif_stats,
    phylogeny, validation_gate -- hepsi beyan edildi), ucu yalnizca metinde
    adini gecirmisti.

    (b) ignore kumesi adimin KENDISINI de iceriyordu. build_operons.py
    `ro.sequence`i doldurdugu icin `ro` tablosu operons'un kendi yazma
    kumesindeydi, dolayisiyla ignore'daydi; annotate'in `ro`da yaptigi
    degisiklik operons'a GORUNMEDI. Sonuc: ro=19.114 iken operon ve
    ro_regulation 19.109'da kaldi. Duzeltme: ignore artik YALNIZCA sonraki
    adimlarin yazdiklarini iceriyor. Adimin kendi yazdiklarini yok saymaya
    gerek yok, cunku karsilastirma adimin KOSTUKTAN SONRAKI sayimina karsi
    yapiliyor.

    IKISINI DE KAPI YAKALADI, ben degil: (a) "aciklanmamis fazla tip
    7_701_AerC", (b) "confirmed-RO count agrees across tables". Kapinin HTML
    uretiminden ONCE durmasi tam bu is icin.

    KURAL: bir adim bir dosyayi okuyorsa grafikte girdi olarak YAZILIR. Yeni
    adim eklenince yukaridaki tarama tekrar kosulur.

70. **Atif kutugu (`citations.yaml`) + dogrulayici -- BITTI.**
    Atiflar uc yere dagilmisti (types/*.yaml serbest metin, disagreements.py
    icinde gomulu DOI, batie_classes.yaml ayri blok) ve 75 tipin yalnizca
    7'sinde cozulebilir tanimlayici vardi. Simdi tek kutuk, 13 atif, HEPSI
    Crossref'e karsi dogrulandi (`citations_check.py`). Dogrulama uc hata
    buldu: iki eksik kunye, bir YANLIS YAZAR LISTESI (CAO makalesi: dogrusu
    Oster U, Tanaka R, Tanaka A, Rudiger W; id de espineda1999 -> oster2000),
    ve betigin kendi hatasi (Crossref basliklarda <i> etiketi donduruyor,
    normalizasyon onu "ia" okuyup basligi uyusmaz sayiyordu).
    Kutuk manifestolara `citations:` alaniyla baglandi ve kontrol eklendi:
    manifestodaki her anahtar kutukte olmak zorunda (negatif sinamayla
    dogrulandi -- uydurma anahtar FAIL veriyor).
    KURAL: uydurulmus ya da hatirlanan DOI GIRMEZ. Her kayit `verified_by`
    tasir: local_pdf, crossref_record, uniprot_reference_list...

71. **Motif istisna kaydi (`motif_exceptions.yaml`) -- BITTI.**
    AerC referans hizalamasinda repo'nun KENDI motif kapisini gecemiyor
    (4/4 Rieske, 0/3 katalitik) cunku 717 aa'lik dizi global hizalamada
    katalitik kolonlara ulasmiyor. Sette kendi kapisindan gecemeyen bir
    referans varsa bu YAZILI olmali. Her istisna uc sey tasimak zorunda:
    gerekce (olcumle), kanit (hangi dosyada), sinir (neyi kapsamadigi).
    Dogrulayiciya uc kontrol eklendi, ucu de negatif sinamadan gecti:
    gecmeyen referansin kaydi var mi, olu istisna kaldi mi, kayit gerekce +
    kanit + sinir tasiyor mu. Ici bos istisna KABUL EDILMIYOR.

72. **SABIT REFERANS SAYISI SUPURGESI -- BITTI.**
    Referans seti 71 -> 76'ya cikti ve sayi ON ALTI yerde sabit yaziliydi:
    12 sablon, make_report.py'de 3, provenance.py'de 1. Hepsi yanlisa dusmustu.
    Simdi veriden geliyor: sablonlar `N_REFS` (chemistry.csv satir sayisi),
    make_report.py modelin NSEQ basligindan okuyor (yer tutucu, cunku metinler
    f-string olmayan bloklarin icinde).
    DAHA KOTUSU IKI BILIMSEL IDDIAYDI: ana sayfada "70 of 71 Asp, one Asn, no
    Glu" yaziyordu, olculen dagilim E=45 D=28 (TERSI); ve "the five Glu
    references" yaziyordu, olculen E=2. Ikisi de artik motif_stats.py'nin
    hesapladigi referans dagilimina bagli, yani referans seti degisince
    kendiliginden guncelleniyor.
    Bir yer BILEREK sabit birakildi: make_report.py'de gecmisi anlatan cumle
    ("HMM'ler o anda 71 referansla yeniden kuruldu") -- onu dinamik yapmak
    tarihi tahrif etmek olurdu, yerine "tarihsel" diye isaretlendi.

73. ~~HAFIF STATIK SITE~~ **YAYINDA (2026-10-24).**
    gh-pages 739ce9bd15, 3.951 dosya, 256+256 veri parcasi. Canli dogrulandi:
    ana sayfa 19.114 giris, 12 rastgele giris 12 ayri parcada tarayicida
    cizildi (2 SVG, 2 tablo, hata yok). Uzak agac yerel derlemeyle birebir
    (tek fark .nojekyll).
    DEPLOY KORUMASI: MIN_FILES=5000 esigi bu tasarimda yanlis alarm verdi
    (tam derleme 3.950 dosya) ve kaldirildi. Yerine beklenen PARCA KUMESI
    arama indeksinden md5 ile hesaplaniyor; tip/varyant/giris sayilari
    chemistry.csv ve veritabanindan okunuyor. 254 != 256'da durdugu sinandi.
    HATA KAYDI: ilk degisiklik yine esik kullandi (MIN_SHARDS=200) ve CANLI
    uzak depoya karsi sinandi; 16 parca silinmisken 240 esigi gecti, koruma
    bos cikti ve o kosu EKSIK bir derleme yayinladi. Sonraki kosu tamamladi.
    Iki ders: disa donuk betik canli hedefe karsi sinanmaz; basarisiz oldugu
    gorulmemis koruma, koruma degildir.
    Olculdu: giris basina bir HTML, 22,2 KB; 19 bin girisle 405 MB ve site
    toplami 474 MB / 19.000 dosya. Sayfanin %66'si komsuluk SVG'si ve komsu
    genler tablosuydu -- AYNI VERININ iki render'i, HTML'e gomuldugu icin 19
    bin kez tekrarlanmis halde.
    Simdi: tek kabuk (/entry.html?id=...) + data/ altinda 256 parca JSON.
    Olculen sonuc: 172 MB / 3.930 dosya.
    PORT RISKI NASIL YONETILDI -- uc sinama, hepsi kalici:
      verify_shard_key.py     19.109/19.109 kimlikte md5 Python ile birebir
      verify_entry_render.py  6.000/6.000 giriste SVG dizgisi birebir
      compare_entry_text.py   40/40 giriste tarayici metni Jinja ile birebir
    Port sirasinda DORT hata bulundu ve hepsinde JS Python'a uyduruldu
    (sinama gevsetilmedi, cunku birebir esitlik cok daha guclu bir invaryant):
      u('') statik kipte sonuna .html ekleyip tabani bozuyordu -> veri 404
      Python yarimi CIFTE yuvarliyor, JS yukari (SVG'de 10/6000 ayrisma)
      ayni sorun %.0f'te (kaplama 98,5 -> Python "98", JS "99")
      ustel bicim: Python 0.0e+00, JS 0.0e+0
    operon_regulator_svg (157 satir) JS'e PORT EDILMEDI: sekil sunucuda
    hesaplanip (layout, arch, family, bool(bp)) basina bir kez saklaniyor,
    girise ozgu iki metin istemcide yerlesiyor. Port edilmeyen kod ayrisamaz.
    SIRADA: commit, freeze, check_site.py, deploy. Canli site hala 6 Ekim'den
    (11.422 giris); yerelde 19.114 var.

68. ~~AerC karboksilatini YAPIYLA yerelleştir.~~ **BITTI 2026-10-21: Asp667.**
    Kullanici istedi: "3d bi protein yapip baksan ya gercekten var mi onlar diye".
    ESMFold ile kalinti 440-700 katlandi: His513 NE2 ve His518 NE2 arasi 3,22 A
    (NDO kristali 3,49), Asp667 OD1 2,64 A ve OD2 3,21 A (NDO Asp362 3,41/3,20).
    En yakin rakip Glu510 4,18 A'da, yani atama belirsiz degil. Aralik hipotezi
    (E662-E701) D667'yi iceriyordu; geometri onu secti.
    Rieske merkezi de ayni yolla dogrulandi: C146/H148/C168/H171 icin alti
    ligand-arasi mesafenin hepsi, ayni yontem NDO'ya uygulandiginda cikan
    degerlerden en fazla 0,34 A sapiyor -- ve yontem NDO'nun KENDI kristal
    geometrisini ~0,3 A icinde yeniden uretiyor (pozitif kontrol budur).
    Kopru D511 mononukleer His'lerden 8,4-9,2 A: NDO'da kopru Asp205 kendi
    zincirinin Rieske'sinden 46,7 A, KOMSU zincirin Rieske'sinden 7,5 A uzakta,
    yani bu mesafe kopru roluyle tutarli.
    SINIRLAR: deneysel yapi degil; tam boy katlanamadi (API 413), iki alan AYRI
    katlandi, yani merkezler ARASI iliski ve alt birimler arasi kopru bu olcumden
    cikmaz; metal modellenmedi; Asp667 tahmin, mutagenez kesinlestirir.
    Kanit: analysis_out/aerc_structure/aerc_metal_sites.json ve oradaki uc PDB.
    DUZELTILEN HATA: ilk olcumde His icin ND1/NE2 ayrimi yapilmamisti ve
    His-His 5,52 A cikmisti; gercek yapiyla ayni atom secilince 3,35 A.

67. **rebuild.py grafigi run_all.sh'in 14 adimini eksik tasiyordu -- DUZELTILDI, dogrulama kosusu suruyor.**
    Bulunus: taze DB ile ilk tam kosuda evidence_tiers "no such column:
    sequence" ile dustu. ro.sequence'i build_operons.py dolduruyor ve o
    grafikte YOKTU; eski DB sutunu elle yapilmis bir kosudan tasidigi icin
    kimse fark etmedi (ro_leaf hatasinin aynisi). Tarama: run_all.sh'taki
    32 betigin 14'u grafikte yoktu, 6'si DB'ye yaziyor (build_operons,
    discover_subfamilies, recursive_homogenize, characterize_leaves,
    classify_domains, variant_signature). FINGERPRINT_TABLES bile grafigin
    hic uretmedigi tablolari (operon, ro_leaf, neighbor_protein) izliyordu.
    Yani `--check` "guncel" derken bu ciktilar eski DB'den kalma olabilirdi.
    Simdi 47 adim, run_all.sh'taki sirayla; her yeni komutun argumanlari
    betigin --help'ine karsi dogrulandi. Kosu: roar_shots/rebuild_final2.log.
    KURAL: run_all.sh'a betik eklenirse rebuild.py ARTEFACTS'a da eklenir.
    EK BULGU 1 (ayni gun): habitat adimi taze DB'de "replicon_source has 8
    columns but 7 values" ile dustu. host_kingdom sutunu sonradan ALTER ile
    eklenmis, iki INSERT konumsal kalmisti; eski DB'de tablo dolu oldugu icin
    parse yolu HIC calismamisti. INSERT'ler sutun adlariyla yazildi.
    EK BULGU 2: rebuild.py yakinsamiyordu. Bir adimin DB girdisi tum DB'nin
    parmak iziydi; kendisinden SONRAKI DB-yazar adimlar onu her kosuda
    bayatlatiyordu (--check asla "guncel" diyemezdi). Ayrica bayatlik bir kez
    hesaplaniyor, ust adim yeniden uretilince alt adimlara yayilmiyordu.
    Duzeltme: DB girdisi artik TABLO BAZINDA; her adim yalnizca kendisinden
    ONCEKI adimlarin yazdigi tablolara bakar (yazar haritasi statik: betik
    kaynagindaki SQL + dinamik: kosuda olculen degisim). Bayatlik sirayla
    ileri yayilir. Kabul testi: tam kosudan hemen sonra --check "guncel"
    demeli. Kosu: roar_shots/rebuild_final4.log.

65. **AerC -- gercek RO, kendi grubunu gerektiriyor. KULLANICI KARARI BEKLIYOR.**
    Kanit ve tum sayilar: `analysis_out/aerc_resolution.json`. ScAerC =
    BAZ23339.1 (NIES4073_42270, 717 aa; organizma yeni adiyla Kalymmatonema
    gypsitolerans NIES-4073). Ligandlar parca parca hizalamayla dogrulandi:
    Rieske C146/H148/C168/H171, kopru D511, mononukleer H513/H518;
    karboksilatin kanonik kolonu (479) proteinde YOK. hmmsearch E=3.8e-68.
    DERS: 600+ aa proteinlerde global hmmalign katalitik kolonlara ulasmiyor;
    `hmmalign --trim` kirptigi icin parca-ici sayac kalinti numarasi VERMEZ.
    Ayrica mononukleer regex karistirilmis dizide %11-28 uyuyor -- tek basina
    kanit DEGIL.
    KARBOKSILAT (kullanici itirazi, kabul edildi): mononukleer Fe deneysel
    olarak var (3.00 Fe/monomer), o halde onu tutan karboksilat da VAR.
    Onceki "2/3" ifadesi hizalamanin ozelligiydi, proteinin degil. Ama KONUMU
    belirlenemedi: modelin 429-570 bolgesi AerC'de yok; mafft aile hizalamasi
    AerC'de bozuk (RedG H161/H166'ya karsi K306/T311 koyuyor), o yuzden oradan
    cikan E693 kanit degil. Aralik hipotezi RedG'den: E662-E701 penceresi.
    Olculen ve cikarilan AYRI kaydedilir. AerC madde 64'un (geometrik capa)
    ilk test vakasi.

66. **Arkeal/fungal RO tahmin katmani (`predictions/`).** Analiz bitti, katman
    henuz kurulmadi. `analysis_out/archaeal_fungal_ro_prediction.json`,
    `analysis_out/archaeal_clade/`. 109 arkeal + 654 fungal protein tam
    katalitik donanimli (negatif kontrol 0/66). Arkea monofiletik DEGIL: 7 saf
    kume, her birinin kardesi bakteriyel; en yakin bakteriye maks %51.7.
    KURAL: types/, chemistry.csv, HMM kutuphanesi ve istatistiklere GIRMEZ.
    SAYFA BAGLANDI (2026-10-22): /atlas/predictions, nav'da "How far to trust it"
    bolumunde; indirme rotasi beyaz listeli (bilinmeyen dosya ve dizin tasmasi
    404, sinandi). Sayfa duzeltme kaydini GORUNUR tutuyor: once yanlis bildirilen
    uzunluk araligi ve "hepsi grup 3" iddiasi, sebebiyle ve argumana neye mal
    oldugu yazili.
    Sirada: IQ-TREE UFBoot + AU testi; elektron ortaklarinin genomik baglami.

64. **AKTIF BOLGE IMZASI -- katalitik triada gore GEOMETRIK hizalama.**
    Kullanicinin onerisi ve veritabaninin bugunku en zayif noktasini tam
    hedefliyor.

    FIKIR. Katalitik merkez (2-His-1-karboksilat) butun RO'larda degismez.
    O halde dizi hizalamasi yerine ONU capa yap: 3B yapilarda triada gore
    superpoze et, triadin 8-10 A cevresindeki kalintilari al, ve hem KIMLIK
    farklarini hem KONUMSAL kaymalari olc. Cikarilacak sey reaksiyon sinifi
    ve substrat sinifi basina bir aktif bolge IMZASI.

    NEDEN BU, DIZI HIZALAMASINDAN IYI. Kolon tabanli okuma bu veritabaninda
    somut olarak kirildi: karboksilat kolonu referans seti 71'den 67'ye
    dusunce 355'ten 360'a kaydi, ve karboksilat denetimi "eksik" denen 1.414
    girisin %95'inde kolonun ve komsularinin TAMAMEN BOSLUK oldugunu buldu --
    yani olculen sey biyoloji degil, hizalama kaymasiydi. Geometrik capa
    indellerden etkilenmez. Ayrica capanin gercekten degismez oldugu olculdu:
    o 1.414 girisin HEPSI iki katalitik histidini ve kopru aspartatini
    tasiyor.

    SINYAL VAR MI? Evet, olculdu. Ceb korunmasi bant bant ayrildiginda
    olculebilir araliklarda ceb, protein ortalamasinin 0,15-0,35'i hizinda
    degisiyor ve her tip ayni yonde (q=4,6e-05). Yani aktif bolgede
    cikarilacak bir imza var; pooled olcum onu doyma yuzunden gizliyordu.

    ELDE NE VAR. `active_site.py` zaten 5 A ve 8 A yaricaplarinda ceb
    kalintilarini cikariyor, ligand ve elektrostatik bloklari var; 17 deneysel
    yapi + 47 AlphaFold modeli ile 67 tipin 62'si kapsaniyor.

    EKSIK OLAN UC SEY:
      a. TRIADA GORE SUPERPOZISYON. Su an ceb, dizi kolonlarindan
         tanimlaniyor; geometrik hizalama yok.
      b. KONUMSAL KAYMA olcumu. Su an yalnizca kalinti KIMLIGI sayiliyor;
         "ayni kalinti ama 2 A kaymis" bilgisi kayboluyor, oysa imzanin asil
         bilgisi orada olabilir.
      c. UYELERE TASIMA. Imza tip duzeyinde cikarilabilir (62 yapi), ama
         11.422 uyenin hepsi katlanamaz. Tasima dizi uzerinden yapilacaksa
         bu, kacinmaya calistigimiz kolon sorununu geri getirir -- ZAYIF
         HALKA BURASI ve bastan boyle yazilmali.

    BEKLENEN CIKTI: reaksiyon sinifi basina "bu kimyayi yapan enzimlerin
    aktif bolgesinde su kalinti su konumda" turunden, literature
    sunulabilecek bir tablo. Bugunku siniflandirma dizi benzerligine
    dayaniyor; bu, MEKANIZMAYA dayanan bir siniflandirma olurdu.


63. **YENI TIP EKLEME: modulerlik (kullanici acikca istedi).** Yeni Rieske'ler
    kesfedilecek; eklemek kolay ve modular olmali. Bugunku durum degil:
    bir tip `ROs_71_Clean/refs71.fasta` + `chemistry.csv` + (dolayli olarak)
    hizalama ve profil kutuphanesi arasina dagilmis, hicbir sey ikisinin
    tutarli oldugunu DOGRULAMIYOR, ve kimin ne zaman hangi makaleden ekledigi
    HICBIR yerde yazmiyor.

    `build_reference_models.py` mekanik kismi (fasta -> clustalo -> hmmbuild ->
    motif kolonlari) zaten cozuyor. Eksik olan dort sey:

    a. **Tip basina TEK bildirim dosyasi** (`types/<id>.yaml`): dizi, substrat,
       urun, reaksiyon, PDB, kaynak DOI, guven, ekleyen, tarih. Tek kaynak.
    b. **Eklemeden ONCE etki yaricapi.** Atama YARISMALI: yeni bir referans
       mevcut tiplerden uye CALAR. 10 bos tip tam olarak bu. Arac, yazmadan
       once "bu tip X'ten 47, Y'den 12 uye alir" ve "bu tip kendisi BOS
       kalir, cunku her yarismayi Z'ye kaybediyor" demeli.
    c. **Tekrar kontrolu** (24. madde): yeni dizi mevcut bir referansla
       >=%99 ayniysa reddet. Bu tek kontrol NdmC, CarAa, NahAc ve NDO 3_315
       vakalarinin dordunu de bastan onlerdi.
    d. **Ucuz yeniden atama.** Pahali adim 17 bin genomdan baglam cikarmak;
       atama ise `ro_alpha.fasta`'daki 73.272 hazir aday proteine karsi
       kosuyor. Bu ikisi AYRILIRSA tip eklemek saatler degil dakikalar surer.


Acik isler. Numaralar kalicidir, commit mesajlari onlara atif yapiyor.

### Bu oturumda istenen ve HENUZ YAPILMAYAN isler

61. **EN ONEMLI BULGU: genomlarin %17'si BOS indirilmis.** Kullanicinin
    "model tanimlayan proteini disarida birakiyor gibi" sezgisi dogru cikti,
    ama sebep HMM degil: ilgili genom hic OKUNMADI.
    * `NZ_CP049045.1` (Pseudomonas sp. BIOMIG1BAC, kullanicinin kendi makalesi,
      doi:10.1128/MRA.00309-20) `replicon` tablosunda VAR ama `ro` tablosunda
      o genomdan TEK bir aday yok. Indirilen dosya 4.427 baytlik bir CON
      iskeleti: ozellik yok, dizi yok.
    * Yeniden indirildi (`efetch rettype=gbwithparts`): 19,5 MB, 7.119 CDS. Ve
      qxyA referans dizisi orada BIREBIR bulundu, makalenin verdigi
      6.844.288..6.845.439 koordinatlarinda, protein WP_068587002.1.
    * Olcek: **2.975 replikon (%17,4) `no_annotation` durumunda**, hepsi ayni
      sekilde ~5 KB'lik iskelet. Ve bunlar rastgele degil: en cok Streptomyces
      (248), Pseudomonas (191), Burkholderia (114) -- yani Rieske oksijenaz
      bakimindan EN ZENGIN cinsler.
    * Maliyet: kayit basina ~19,5 MB, toplam ~57 GB; NCBI'nin saniyede uc
      istek sinirinda en az 17 dakikalik istek, arti transfer.
    Sonuc: uye sayilari sistematik olarak EKSIK ve eksiklik enzim bakimindan
    zengin cinslerde yogunlasiyor. Kullanici onayi gerekiyor (disk, bant
    genisligi, ve sonrasinda pipeline'in bastan kosturulmasi).

48. **Grup bazinda sayfalar.** Her RO grubunun genel ozellikleri: baskin operon
    mimarisi, regulator ve promotor tipleri.

49. **Reaksiyon bazinda sayfalar.** Ayni reaksiyonu yapan enzimlerin ortak
    ozellikleri.

51. **Promotor bolgesi analizi.** Kullanici artik ACIKCA istedi. Veri yok:
    kayitlarin %88,9'u CON tipinde ve dizi dosyada degil. NCBI E-utilities ile
    ~10.200 intergenik bolge cekmek gerekiyor (14. madde).

52. **Gen sirasi (gene order) karisik.** Kanonik bir sira gosterilmiyor, bir
    suru varyant yan yana duruyor ve uzak akraba icin anlamli degil. Yakin
    akrabada sira boyle, uzaklasinca operon/regulator/promotor soyle degisiyor
    seklinde gosterilmeli.

53. **Transpozon: tur ozelinde mi enzim ozelinde mi?** Istatistik olarak
    olculmeli.

54. **Benzer enzimlerin operonlari da benzer mi?** Degilse eslesme rastgele
    olabilir: ayni gen adi tasiyip anlamsiz olan durumlari ayirt etmek icin bir
    olcu gerekiyor.

55. **Methods sayfasi daha ayrintili.**

56. **Elektron transferi cizimi.** Rieske merkezi ve katalitik demir arasindaki
    elektron yolu icin bir sekil.

57. **Menu sadelestirme ve daha zarif anasayfa.** Menu dar ekranda kaydirmali
    hale getirildi ama hala cok ogeli; anasayfa daha zarif olabilir.

58. **Tasarim: tamamlanmis ama yayinlanmamis oneriler.** Tip olcegi, istatistik
    seridi, figur altyazilari ve chip paleti degisiklikleri gozden gecirilmeyi
    bekliyor; olculen kusurlar (telefonda yatay kayma, karanlik modda gorunmeyen
    cizim etiketleri) ayri olarak uygulanacak.

59. **Literaturden yeni tip adaylari.** Dokuz aday erisim numarasiyla hazir
    (`literature_review.md`): Nvd/DAF-36, CAO, PAO, TIC55, PrnD, RedG/McpG,
    MsmA, TcsAB, BrhABC. Bes yeni reaksiyon sinifi getiriyorlar.

60. **Kuratorluk duzeltmeleri (43. maddeden).** NdmC aslinda NdmB; iki yanlis
    PDB kimligi; NagGH muhtemelen salisilat 5-hidroksilaz; kumen ve
    izopropilbenzen ayni substratin iki adi. Hepsi dogrulandi, uygulanmadi.

43. **Kuratorluk hatalari -- LITERATUR AJANI BULDU, DOGRULANDI, KARAR SENDE.**
    `literature_review.md` tam dokum. Dizi kimligi ile dogrulananlar:
    * `1_115_NdmC` **NdmC DEGIL, NdmB'nin kendisi.** Fark yalnizca 18
      kalintilik `MGSSHHHHHHENLYFQGS` saflastirma etiketi. Kendim dogruladim:
      NdmC dizisi NdmB'nin TAM alt dizisi ve aradaki fark BIREBIR o etiket.
      Yani 24. maddede "parca" diye yazdigim iliski bir kristalografi
      etiketinin eseri. Gercek NdmC (UniProt M1EY73, 284 aa) sette YOK ve
      kendi Rieske alani da yok, NdmD'den odunc aliyor.
    * **Dort kayit saflastirma etiketi tasiyor:** NdmA (+18), NdmB (+18),
      NagGH (+19), CARDO (+8, C-ucunda LEHHHHHH). Olctum: etiketler motif
      modelinde ESLESME sutunlarina girmiyor (hizalamada bosluk olarak
      kaliyor), yani profil onlari korunmus saymiyor -- ama NdmB/NdmC
      iliskisini URETEN sey tam olarak bu etiketti.
    * `3_307_NagGH` muhtemelen **salisilat 5-hidroksilaz** (UniProt O52379,
      EC 1.14.13.172), naftalin-2-sulfonat dioksijenaz degil.
    * `3_314_NDO`'nun PDB kimligi **1NDO yanlis**: 1NDO farkli bir proteinin
      (449 aa Pseudomonas) yapisi; NarAa'nin yapilari 2B1X / 2B24.
    * `1_102_CARDO` ve `1_103_CarAa` ikisine de 1WW9 verilmis; 1WW9 yalnizca
      J3 susunun yapisi.
    * `2_203_cumA1` ile `2_206_IpBAa` ayni substrati IKI ADLA tasiyor
      (kumen / izopropilbenzen, ayni SMILES).
    * `3_317_NidA` ile narAa %98,1-98,9 ayni; benim %99 esigim bunu kacirdi.
    * EdoA1 / cumA1: tek kalinti farki (353 Leu/Trp) ve o kalinti demirden
      15,4 A uzakta, yani substrat ayrimini aciklayamaz. EdoA1 icin
      saflastirilmis enzim deneyi hic yapilmamis. Yani substrat etiketleri
      IZOLASYON substrati etiketleri.

44. **Yeni tip adaylari (literatur).** Dogrulanmis erisim numarasiyla dokuz
    aday: Nvd/DAF-36, CAO, PAO, TIC55, PrnD, RedG/McpG, MsmA, TcsAB, BrhABC.
    Bes yeni reaksiyon sinifi getiriyorlar (oksidatif karbosiklizasyon,
    N-oksijenasyon, desaturasyon, makrosiklik acilma, alifatik C-S kirilmasi).
    Ayrica SxtT ve GxtA'nin substratlari artik kesin (beta-saksitoksinol ve
    saksitoksin), `chemistry.csv` bunlari belirsiz yaziyor.

5. **Hugging Face Space** (tam uygulama + dizi siniflandirici). `webapp/Dockerfile`
   hazir ve test edildi; kullanicinin HF hesabi gerekiyor.

14. **Promotor dizisi analizi — VERI YOK, kullanici karari gerekiyor.**
    Olculdu: indirilen 17.073 GenBank kaydinin %88,9'u `CON` tipinde, yani
    dizi dosyada BULUNMUYOR (BioPython `UndefinedSequenceError` veriyor).
    Dizi tasiyan 1.897 kayit ookaryot mRNA'si. Yani -35/-10, operator
    tekrarlari ve transkripsiyon baslangici bu veriyle analiz edilemez;
    `extract_genomic_context.py`'yi degistirmek yetmez.
    Cozum icin iki yol var, ikisi de disa donuk ve kullanici onayi ister:
      (a) NCBI E-utilities ile yalnizca gereken intergenik bolgeleri cekmek
          (`efetch` + `seq_start`/`seq_stop`): ~10.200 kucuk istek, API
          anahtari olmadan ~1 saat, NCBI kullanim politikasi geregi e-posta
          ve anahtar belirtmek gerekir.
      (b) 15.173 kaydin dizili surumunu yeniden indirmek: cok daha buyuk
          trafik ve disk, ama tek seferlik.
    Onerim (a); hangi bolgelerin cekilecegi zaten `ro_regulation` tablosunda
    hazir (10.186 giriste intergenik bolge koordinatli olarak duruyor).

24. **Referans setindeki tekrarlar (KULLANICI KARARI).** 71 kayit, 68 tekil
    dizi. Durum artik OLCULDU ve sitede GORUNUR (39. madde); karar hala sende,
    cunku duzeltme profil kutuphanesini yeniden kurup her sayiyi yeniden
    uretmeyi gerektiriyor.
    * Birebir ayni: OxoO(2)=OMO(11), NahAc(0)=NDO 3_315(0),
      NDO 3_314(2)=NarAa(30). Hangisini tutsan bilimsel olarak AYNI sey kalir,
      cunku diziler ayni; karar hangi ADIN tutulacagi.
    * Parca: NdmC (355 aa) NdmB'nin (373 aa) icinde tam geciyor. Burada secim
      yok, kisa kayit bir parca.
    * En kotu durum bir kopya DEGIL: EdoA1 / cumA1 %99,8 ayni ama deneysel
      substratlari FARKLI (etilbenzen / kumen) ve uyeler 29'a 6 bolunuyor,
      yani fazlalik bir substrat ETIKETINE siziyor. CARDO / CarAa %99,2 ama
      ikisinin substrati da karbazol, dolayisiyla zararsiz.
    Onerim: birebir ayni ciftleri tek tipte birlestirip iki adi da es anlamli
    gostermek (hicbir sey atilmaz), NdmC'yi parca olarak isaretlemek, EdoA1 ile
    cumA1'i ayri tutup ikisine de "bu cift ayirt edilemiyor" notu koymak.

    **2026-10-07 eklemeleri -- tam all-vs-all olculdu (diamond, kapsama >=0,8):**
    * **OxoO / OMO ayni dizi ama FARKLI GRUPTA (grup 1 ve grup 3).** Yani
      `ro_group` en az bu durumda dizinin fonksiyonu degil. Birlestirme karari
      verilirken hangi grubun tutulacagi da secilmeli; bu, grup bazli her
      istatistigi (phylum_by_group, beta_by_group, etc_system_by_group) bir
      giris kadar etkiler.
    * Tekrar, **bos tiplerin dogrudan sebebi**: 10 bos tipin 4'u
      (CarAa, NahAc, NDO 3_315, NdmC) bir kopyanin kaybeden tarafi. Ayni
      diziye sahip iki referans yarismayi beraberlikle bitirir ve biri
      keyfi olarak her seyi alir. Yani bos tipler ayri bir sorun degil,
      tekrarin SONUCU.
    * NdmC ciftinde dikkat: diziler ayni ama `chemistry.csv` substratlari
      FARKLI yaziyor (kafein / 7-metilksantin). NdmC satirinin tarif ettigi
      enzim gercek ve ayri; yanlis olan ona eklenmis DIZI. Yani cozum
      "NdmC'yi sil" degil, "NdmC'ye dogru diziyi koy".
    * %85 ustu toplam **22 cift** var (6'si >=%99). Tam liste
      `analysis_out/reference_redundancy.json`'da; esik altindakiler zararsiz
      degil, yalnizca daha az acil.
    * Bundan sonra yeni referans eklenirken bu kontrol OTOMATIK yapilmali --
      63. madde.

## Bitti

62. ~~"Uzak akrabalar okyanus kaynakli mi? Ve bu, enzimlerin kimyasallardan
    ONCE ayrildigi anlamina gelir mi?"~~ **BITTI, ikisi de olculdu.**
    * Giris duzeyinde deniz en yorumlanamaz habitat: %69,7'ye karsi %55,7, OR
      1,90, q=1e-07. Bagimsiz ikinci dilim (kaydin kendi yer adi) ayni yone
      isaret ediyor: %64,3'e karsi %58,2. **Cins duzeyinde HICBIRI ayakta
      kalmiyor** (%66,7'ye karsi %63,4, OR 1,18, q=0,45). Taban da %55,7'den
      %63,4'e cikiyor -- suslari cokertmek yorumlanamaz orani YUKSELTIYOR,
      cunku tekrar tekrar dizilenler zaten bildiklerimiz.
    * Atasal hipotez REFERANSSIZ test edildi (`ancestry.py`): her tipin
      kirli-alan uyeleri, kimyasalin hic uygulanmadigi habitatlardaki
      uyelerine dogrudan hizalandi; EN YUKSEK kimlik, ayrilmanin ne kadar
      yakin olabilecegini sinirliyor. Deniz tanimiyla 16 ksenobiyotik tipin
      **10'unda** hicbir kirli-alan uyesi hicbir deniz uyesine %95'ten yakin
      degil, 2'sinde %60'in altinda (DdmC %32,5). Bu ayrilmalar 80 yillik
      olamaz.
    * Asil destek KONTROLDEN geliyor: ksenobiyotik ve dogal substratli tipler
      habitatlar arasi ayrilma derinliginde (p=0,43 / p=0,32) ve varyant
      yayiliminda (p=0,78) ayirt edilemiyor; kirli-alan uyeleri kendi tipinden
      rastgele cekim kadar varyanta yayiliyor (oran 0,975). Yani kirleten
      bolgeler var olani SECTI, yeni bir soy uretmedi.
    * Iki yontem notu: medyan capraz kimlik kullanilmiyor (tipler klad degil,
      profil kovasi oldugu icin medyanlar %25 ilgisiz-protein tabanina
      oturuyor); ve hicbir yerde molekuler saat yok -- her sayi bir SINIR,
      tarih degil.


50. ~~Regulator / transpozon sayfalari.~~ **BITTI.** `/atlas/elements` dizin,
    `/element/<aile>` tek tek. 8 duzenleyici ailesi + 13 adlandirilmis IS
    ailesi. Her sayfada aile hangi enzim tiplerinin yanindaysa, o tipin GENEL
    payina karsi, BH duzeltmesi ve 5 puan maddilik esigiyle. Ornek: TetR
    gercekten KshA15'te yogunlasiyor (%32,4'e karsi %6,8, 4,8 kat).
    Kullanicinin kendi sorusu da olculdu: transpozonlar TURE mi ENZIME mi
    bagli? Ayni IS aileleri bir kez enzim tipine bir kez cinse gore
    caprazlandi; iki tablo farkli sekilde oldugu icin her biri KENDI null'una
    gore olculdu ve yalnizca fazla karsilastirildi. **Cins kazaniyor: +0,127'ye
    karsi +0,066, yaklasik iki kat.** Yani IS elementleri yanlarina dustukleri
    enzimi degil, icinde yasadiklari organizmayi izliyor.
    Sankey duz SVG; geometri sunucuda hesaplaniyor cunku Plotly'nin sankey izi
    yalnizca 4,5 MB'lik tam pakette var. Seritler renk tasimiyor: tek tur akis
    icin kategorik palete gerek yok.
    Yan bulgu: transpozon regex'i rekombinasyon makinesini de yakaliyordu
    (RuvX 117 kez); adlandirilmis IS aileleri artik ayri cikariliyor.


45. ~~Ekolojik koken sayfasi.~~ **BITTI.** `/atlas/origin` yayinda.

46. ~~Ogrenme sayfasi.~~ **BITTI.** `/atlas/learning` yayinda; istatistik
    sayfasi da kisa bir ozetini tasiyor.

47. ~~Dunya haritasi.~~ **BITTI.** `/atlas/geography`. Olculen sey bastan
    belirtildi: 8.245 girisin (%72,2) ulkesi var, 121 ayri yer, bunlarin
    112'si haritaya oturuyor; okyanuslar ve dagilmis iki devlet ayri listede.
    Asil sonuc OLUMSUZ ve oyle de yazildi: tip x ulke sapma testinde 314
    hucreden yalnizca 5'i maddi, 2'si cins duzeyini geciyor (PhnA1a Cin'de
    +17 puan, KshA15 Kambocya'da +10). Yani dizileme cabasi bolununce neredeyse
    hicbir tip cografi olarak ayirt edici degil. Harita yine de duruyor cunku
    eksiksizligi gosteriyor, ama sayfanin basinda ne OLMADIGI yaziyor.
    Enzim sayfalarindan `#type=` ile derin baglanti var.


42. ~~Agac ve ag "her varyant icin bir temsilci" diyordu; dogru degildi.~~
    **BITTI.** Indirme dosyalarini denetlerken Newick'te 1.276 uc buldum oysa
    1.809 varyant var. Sebep bulundu ve olculdu: temsilci kurali IKI UCLU --
    tekil varyant kendisi temsilcidir (766), 5+ uyeli varyant medoid verir
    (439), aradaki 2-4 uyeli 604 varyant ATLANIYOR cunku bu kadar az diziden
    medoid kararli degil. 766+439 = 1.205 temsilci, +71 referans = 1.276 uc.
    Aritmetik tutarli, ama sayfadaki cumle degildi.
    Atlamanin bedeli kucuk degil: **1.631 giris, yani veritabaninin %14,3'u**
    bu iki gorselde hicbir uca sahip degil. Dort tip ise bu kural altinda hic
    temsilci uretmiyor (OxoO, TDO, DxnA1, NDO 3_314) ve yalnizca referans
    karesi olarak gorunuyor; toplam 10 giris.
    Iki sayfanin metni olculen sayilara baglandi (`atlas.representative_
    coverage`), hicbir sayi elle yazilmadi. `validate_curation.py` artik uc
    sayisinin kuralin ONGORDUGU sayiya TAM esit oldugunu sinar. Ilk surumde
    "uyesi olan referans" sayisini (61) kullanmistim ve 10 fark icin tolerans
    koymak gerekiyordu; agaca kuratorlu 71 referansin tamami girdigi anlasildi
    ve tolerans kaldirildi -- tolerans bir hatayi gizleyebilirdi.
    Yan denetim, hepsi temiz: `all_confirmed.fasta` 11.422 kayit ve DB ile
    dizi dizi ayni, `all_confirmed.csv` 25 sutun ve bozuk satir yok,
    `operons.csv` 49.274 satir (DB ile ayni), Newick parantezleri dengeli.

41. ~~Yayinlanan arama bir KIMYASALI bulamiyordu.~~ **BITTI.** Statik arama
    indeksi yalnizca giris basina alanlari tasiyordu (organizma, urun, tip,
    varyant, aile) ve KIMYA yoktu. Sitenin butun duzeni "enzimleri
    etkiledikleri kimyasallara gore grupla" oldugu icin bu kucuk bir eksik
    degildi. Olculdu, uygulama / yayinlanan site:
    benzalkonium 95 / **0**, terephthalate 273 / **0**, caffeine 447 / **0**,
    carbazole 55 / **0**, biphenyl 107 / 3, naphthalene 226 / 20.
    * Substrat giris basina degil TIP basina bir ozellik, bu yuzden 11.422
      satira tekrarlanmadi: tip basina 71 kayitlik kucuk bir sozluk sayfaya
      gomuldu ve aranan metin tarayicida oradan tamamlaniyor. Indeks boyutu
      degismedi.
    * Statik sayfa artik eslesen ENZIM TIPLERINI de ayri bir tabloda
      gosteriyor; barindirilan surumde bu vardi, statikte yoktu. Uye sayisi
      sifir olan 10 tipin substrati (xylene, isopropylbenzene,
      2,4-dinitrotoluene, 3-phenoxybenzoate) yalnizca bu yolla bulunabiliyor.
    * Statik sayfa kisa tip adini da cakisma ayrimi OLMADAN yaziyordu
      (iki ayri "NDO"); 29. maddedeki sozluk oraya da baglandi.
    Yan bulgu ve ikinci duzeltme: iki arama ayni veriden FARKLI cevap
    veriyordu. FTS5 token esler, dolayisiyla "toluene" sorgusu
    "4-toluenesulfonate" substratini bulmuyordu (o metin "4" ve
    "toluenesulfonate" olarak tokenleniyor) ve uygulamada 4, statik sitede 660
    sonuc cikiyordu. Uygulama artik tam eslesme 25'ten az sonuc verdiginde
    onek eslesemesiyle genisletiyor (659) ve sayfa bunu ACIKCA soyluyor,
    cunku okuyucu "toluene" arayip "toluenesulfonate" gorunce bunun kasit mi
    hata mi oldugunu bilemez.
    `check_site.py` artik 71 kuratorlu substratin hepsinin yayinlanan aramada
    bulunabildigini sinar; sayfanin IKI arama yolunu da (giris indeksi ve tip
    sozlugu) taklit eder, yoksa dogru calisan sayfayi hatali bildirirdi.

39. ~~Referans fazlaligi okuyucuya gorunmuyordu.~~ **BITTI.**
    `reference_redundancy.py` + veri kalitesi sayfasinda 4. bolum + etkilenen
    her tip sayfasinda uyari kutusu. Sorun 24. maddede duruyordu ama SITEDE
    hicbir izi yoktu: okuyucu NahAc sayfasinda "0 uye" goruyordu ve nedenini
    ogrenemiyordu. Artik etkilenen 12 tipin her birinde, uye sayisinin kismen
    sayisal bir kazadan geldigini ve ciftin toplaminin kac oldugunu soyleyen
    bir kutu var. 150 uye bu durumda.
    Olcumu ilk kosuda YANLIS yaptim ve kendi capraz kontrolum yakaladi: FASTA
    basliklari gen islevini de tasiyor ("3_309_NahAc_dioxygenase_pro") oysa tip
    kimligi "3_309_NahAc"; eslesme olmadigi icin butun uye sayilari 0 cikti ve
    yakin-ayni ciftler iki kez listelendi. Kimlikler kuratorlu tip listesiyle
    onek eslesemesine baglandi; `validate_curation.py`'ye bes capraz kontrol
    eklendi ve biri tam olarak bu hatayi yakaliyor.

40. ~~ROADMAP duzenlerken neredeyse yarisini sildim.~~ **BITTI** (ders).
    Maddeyi "24'ten 26'ya kadar" diye DILIMLEYEREK degistirmeye calistim, oysa
    24 acik bolumde, 26 bitmis bolumde: dilim arada kalan `## Bitti` basligini
    ve 37-39. maddeleri de kapsiyordu. `assert` yazma islemini durdurdugu icin
    dosya bozulmadi. Ders: metin dosyalarinda ARALIK dilimlemek yerine tam blok
    eslesmesi yapilacak, ve her duzenlemeden sonra dosya boyu kontrol edilecek.

37. ~~Esik degisince sonuclar ne kadar degisiyor?~~ **BITTI.**
    `threshold_sensitivity.py` + veri kalitesi sayfasinda 3. bolum. Yontem
    sayfasi 0,45 kapsama esigini KALIBRASYON setinde gerekcelendiriyordu
    (2.788 alfa, 1.570 alfa-olmayan) ama sorulan soru farkliydi: "o esik
    degisince burdaki bir suru sey degisebilir". Dogru cevap veritabaninin
    KENDI sayilarini esik boyunca yeniden hesaplamakti; kapsama, motif durumu
    ve karboksilat kimligi giris basina sakli oldugu icin HMMER'i yeniden
    kosturmak gerekmedi. Sekiz esik (0,35 - 0,80):
    * Sertlik iki katina cikarildiginda (0,45 → 0,80) girislerin %34'u
      gidiyor ama 61 tipin 59'u kaliyor.
    * En guclu yayinlanan iliski (kopruleyen karboksilat x grup) ZAYIFLAMIYOR,
      GUCLENIYOR: V 0,843 → 0,956. Asp:Glu 2,53 → 2,91. Plazmid orani
      %2,6 → %3,2, yani neredeyse sabit.
    * Esige GERCEKTEN bagli tek sayi okaryot orani: %8,2 → %3,6, cunku
      okaryot Rieske proteinleri bu bakteriyel modellere kismi uyuyor. Yani
      daha sert bir esik "daha temiz" bir veritabanini okaryot dalini sessizce
      silerek satin alirdi. Bu, gevsek esigi KORUYUP her girisin yasam alanini
      etiketlemek icin bir gerekce.
    * Esigi asagi cekmek neredeyse hicbir sey eklemiyor (0,35'te +27 giris),
      cunku aday kumesi orada zaten seyrek.
    Sinirlar gizlenmiyor: kalinti kimligi yalnizca onaylanmis kume icin
    olculdu, bu yuzden 0,45 altindaki satirlarda o sutunlar BOS birakiliyor;
    ve yeniden hesaplanan kapi pipeline'in bir durum kontrolunu atladigi icin
    0,45'te 11.441 giris tutuyor (yayinlanan 11.422).

38. ~~Indirme izin listesi iki yerde kopyaydi.~~ **BITTI.** Ayni dosya listesi
    hem `app.py` indirme rotasinda hem `freeze.py` icinde duruyordu. Yeni
    dosya eklenirken biri guncellenip oteki unutulunca sonuc sessiz oluyordu:
    `threshold_sensitivity.json` indirme tablosunda GORUNDU ama rota 404
    donuyordu. Liste artik `atlas.PROVENANCE`tan turetiliyor (tek kaynak) ve
    bu turetme hemen bir eskisini de ortaya cikardi: `cooccurrence.json`,
    `substrate_predictability.json` ve `habitat.json` indirilebilir olup
    BELGELENMEMISTI. Uc dosya da tabloya eklendi; artik 21 dosya hem
    belgeli hem indirilebilir. `validate_curation.py` iki yonu de kontrol
    ediyor. Refaktor kendi denetim scriptimi de bozdu (`check_site.py` dosya
    listesini `freeze.py` KAYNAGINDAN kaziyordu ve liste turetilince 21 yerine
    2 dosya buldu); o da ayni kaynaga baglandi.

35. ~~Konak alemi ile enzim kimyasi arasinda iliski var mi?~~ **BITTI** --
    cevap **YOK**, ve bu sayfadaki en ogretici negatif sonuc. Yeni konak
    boyutu bu soruyu sorulabilir kildi. Dort test, iki duzeyde:
    * RO grubu x konak alemi, GIRIS basina: chi2=258, p=3e-51, V=0,19 --
      bakan goz "kesin" der.
    * Ayni test TIP+CINS basina: p=0,38, V=0,10 -- **anlamli degil.**
      Yani giris duzeyindeki devasa p degeri enzimlerin nerede yasadigini
      degil, hangi canlinin dizilendigini olcuyordu.
    * Substrat sinifi x konak alemi: normalizasyondan ONCE bile anlamsiz
      (p=0,078), cins duzeyinde hicbir sey (p=0,98, V=0,02). Yani ksenobiyotik
      yikim yetenegi belirli bir konak iliskisinde yogunlasmiyor.
    Bu cift, sitede "neden her ekolojik ifade cins duzeyinde kuruluyor"
    sorusunun en net kaniti oldugu icin oldugu gibi yayinlaniyor.

36. ~~Istatistik sayfasinda karar etiketi anlamliligi YOK SAYIYORDU.~~
    **BITTI.** Karar yalnizca Cramer's V'ye bakiyordu, bu yuzden p=0,38 olan
    bir test "weak association" diye etiketleniyordu. Artik p >= 0,05 ise
    etiket "not significant" oluyor. On bir capraz tablonun hepsi tek tek
    dogrulandi. Ayni sablonda satir basligi da SABIT "group X" idi; substrat
    sinifi satirlari "group xenobiotic" diye goruluyordu, bu da testten gelen
    bir basliga cevrildi.

34. ~~Konak alani (`/host`) kullanilmiyordu.~~ **BITTI.** Olculdu: 9.308
    replikonun 2.962'sinde (%31,8) `/host` var ve bunlarin **728'inde
    `/isolation_source` HIC YOK**. Yani habitat sozlugunun sessiz kaldigi
    yerde konak alani gercek bir ekolojik bilgi tasiyordu ve atiliyordu.
    Cozum, iki boyutu BIRLESTIRMEK DEGIL ayri tutmak oldu: habitat "nerede
    yasiyordu", konak "neyin icinde/uzerinde bulundu" sorusuna cevap verir ve
    `Homo sapiens` bir habitat degildir (yara da, bagirsak da, deri de
    olabilir). Birlestirmek habitat istatistiklerini bozardi.
    * `host_kingdom.csv`: 395 satir, her satirda konak dizgisi, alemi ve
      KARARIN GEREKCESI. 401 farkli dizginin tamami esleniyor.
    * `replicon_source.host_kingdom` sutunu + `host_kingdom_by_habitat`
      capraz tablosu; ekoloji sayfasinda 4. bolum.
    * Dagilim: bitki 1.379, insan 1.175, hayvan 353, alg 26, mantar 24.
      Habitat'i siniflanamamis ama konagi bilinen **906 replikon** var; alg
      kayitlarinin %85'i bu durumda.
    * Durustluk: 3 kayitta konak alanina mineral ya da "soil" yazilmis
      (kaynak kayittaki kuratorluk hatasi, konak degil) ve 1 ad hem bitki hem
      kelebek cinsi (`Pieris`, karar verilemez). Dordu de atilmadi, ayri
      etiketle gosteriliyor.
    * Genisletilebilir: eslenmeyen her dizgi pipeline tarafindan RAPOR
      EDILIR, tahmin edilmez. `validate_curation.py` bes yeni kontrol:
      sozluk, gerekce zorunlulugu, tekrar yok, veritabanindaki her dizgi
      eslenmis, saklanan degerler sozlukte.
    * Kendi kontrolum bir hata buldu: 6 dizgi yalnizca buyuk/kucuk harfte
      ayriliyordu ("chicken"/"Chicken") ve arama kucuk harfe indirdigi icin
      bunlar fazlaligi; ikisine ayri alem yazilsa kazanan dosya sirasina
      kalirdi. Tek satira indirildi.

26. ~~Habitat sozlugunun kalan bilinen sinirlari.~~ **BITTI.** Uc parca, her
    biri once olculdu:
    (a) **Isaretsiz anatomi artik konakla cozuluyor.** "lung", "blood",
        "tissue" gibi metinler ORTAMI soylemiyor, cunku ciger insanin da
        domuzun da baligin da olabilir; bu yuzden 201 kayit `other`a
        dusuyordu. Olculdu: bunlarin 195'inde `/host` alani DOLU (165 Homo
        sapiens, 26 hayvan, 2 bitki), yani cevabi veri zaten tasiyor.
        Eklenen kural DAR: konak tek basina hicbir kayda habitat atamaz,
        yalnizca metin `other` dondurdugunde VE metinde anatomik bir kelime
        varken alemi soyler. Boylece "homo sapiens bir habitat degil" ilkesi
        korunuyor. 219 kayit cozuldu (189 insan, 28 hayvan, 2 bitki) ve
        dogrulanmis RO kapsami %54,3 → **%56,6**, yani 25. maddedeki
        duzeltmenin bedeli fazlasiyla geri alindi.
    (b) **Sediment kompartmana ayrildi.** `sediment` kurali deniz ve tatli
        sudan once geliyor (dogru: "marine sediment" bir sediment ornegidir),
        ama deniz tabani ile nehir tabanini ayni kovada tutuyordu. Eslesen
        kuralin RAFINESI olarak bolundu, yani sira mantigi bozulmadi ve atama
        hala tek bir anahtara indirgenebilir ("sediment + marine"). 182
        replikon: 72 `marine_sediment`, 18 `freshwater_sediment`, 92 isaretsiz
        `sediment`. Ilk denemede "sea" oneki "seasonal" ile eslesti ve
        "Mud of a seasonal forest creek" deniz sedimenti sayildi; esleme tam
        kelimeye cevrildi ve 14 test vakasi gecti.
    (c) **`food_fermented` iddiasi DOGRU DEGILDI.** Bu madde kategoride 21
        endustriyel/laboratuvar fermentasyonu oldugunu soyluyordu. Olculdu:
        120 kaydin yalnizca 4'u fermentor/reaktor kelimesi tasiyor ve hepsi
        ayni metin, "Film in fermentor of rice vinegar" -- pirinc sirkesi bir
        GIDA. Kategori dogru, degisiklik yapilmadi. Iddia olculmeden
        yazilmis.
    Capraz kontrol: PAH kirlilik zenginlesmesi 5,76x → 5,65x, yani
    degisiklikler kirlilik sinyaline dokunmadi.

31. ~~Statik sitede taban onegi eksikti -- YAYINDAKI GEZINME KIRIKTI.~~
    **BITTI.** `freeze.py --base /pro-sim-blast` verilmeden uretilen bir
    derleme yayinlandi. GitHub proje sayfasi siteyi `/pro-sim-blast/` altinda
    sunuyor, ama sayfalardaki mutlak baglantilar `/about.html` seklindeydi,
    yani alan adinin KOKUNE gidiyordu ve 404 donuyordu. Anasayfa acildigi
    icin hata gorunmuyordu; kirik olan her ic baglantiydi. Olculdu:
    `recepcanaltinbag.github.io/about.html` → 404,
    `recepcanaltinbag.github.io/pro-sim-blast/about.html` → 200.
    Dogru tabanla yeniden uretilip yayinlandi ve `deploy_pages.sh` artik
    taban onegini `site/index.html` icinde ARIYOR; bulamazsa deploy etmeden
    duruyor. Negatif testi yapildi (cikis kodu 1).

32. ~~Uyesi olmayan tiplerin sayfalari statik sitede yoktu.~~ **BITTI.**
    Kuratorlu 71 tipin 10'u bu derlemede dogrulanmis uye toplamiyor, ama
    giris sayfalari "en yakin referans" olarak onlara baglaniyor. Dinamik
    uygulama bu tipler icin bos bir sayfa veriyordu; ihracat ise tip
    listesini VERITABANINDAN aldigi icin onlari hic uretmiyordu ve statik
    sitede 90 baglanti 404 donuyordu. `freeze.py` artik listeyi
    `chemistry.csv` ile birlestiriyor: 61 degil 71 tip sayfasi.

33. ~~Statik ihracat icin baglanti denetimi.~~ **BITTI.** `check_site.py`
    artik yayinlanan siteyi DISKTE de geziyor: 13.318 sayfa, 168.218
    baglanti, hepsi bir dosyaya denk geliyor. Dinamik uygulamayi denemek
    yetmiyordu, cunku `.html` eki ve taban onegi yalnizca ihracatta
    uygulaniyor; yukaridaki iki hata da tam olarak orada yasiyordu.

29. ~~Cakisan kisa tip adlari.~~ **BITTI.** Uc kisa ad referans setinde iki
    kez geciyor: `BphA1` (2_201 ve 2_218), `NDO` (3_314 ve 3_315), `NidA`
    (3_317 ve 3_318). Tablolarda iki satir ayni etiketle IKI AYRI sayfaya
    baglaniyordu. Ad uretimi tek bir yere toplandi (`atlas.short_names`,
    sablonlara `gene()` olarak gecti, ag grafigine JSON ile gidiyor) ve
    yalnizca cakisanlara grup numarasi ekleniyor: "NDO (314)". 68 tipin
    hepsine numara eklemek okunurlugu bosa dusurecekti. Yeni bir cakismayi
    `validate_curation.py` yakalar.

30. ~~Dil kontrolunun yapisal bosluklari.~~ **BITTI.** `euk_group` sutununun
    aylarca Turkce kalabilmesinin iki sebebi vardi ve ikisi de kapatildi:
    (a) kelime listesi dokuz ISLEV sozcugunden olusuyordu, kacan sey bir
        ICERIK sozcuguydu ("bitki/alg"); liste icerik sozcukleriyle
        genisletildi ve onek eslesemesi yuzunden Ingilizce ile cakisabilecek
        kisa parcalar bilincli olarak disarida tutuldu;
    (b) kontrol elle tutulan bir sutun listesine bakiyordu ve o sutun listede
        hic yoktu, yani kontrol calismadi bile. Artik DB'deki 107 metin
        sutununun hepsi uc listeden birinde olmak ZORUNDA (denetlenen,
        GenBank'tan gelen, kimlik); siniflanmamis bir sutun FAIL verir.
    Dil denetimi 5 sutundan 15 sutuna cikti. Iki kontrol de negatif test
    edildi: siniflanmamis sutunu ve eski etiketleri yakaliyor, `plant/alga`
    ve `Pseudomonas putida` gibi mesru metinde sessiz kaliyor.

1. ~~Varyant kumelemesini ham dizide tekrarla.~~ **BITTI.** CD-HIT artik ham
   protein dizisinde bolme yapiyor, kimlik olcumu hizalamada kaliyor
   (`discover_subfamilies.py`, `recursive_homogenize.py`, `analyze_variants.py`).
   Sonuc: 1.790 → 1.809 yaprak, 2.217 → 2.399 alt-aile, novel aday 2.283 → 2.438.

2. ~~Varyant duzeyinde kalinti imzasi.~~ **BITTI.** `variant_signature.py`,
   `cluster_sdp` + `leaf_sdp`; tip ve varyant sayfalarinda matris olarak gosteriliyor.

3. ~~Novel adaylar sayfasi (Atlas).~~ **BITTI.** `/atlas/novel`: 318 aday,
   63 varyant. Tanim seffaf (merkezler tam + >=300 aa + <%25 kimlik + varyant >=3
   uye) ve iki novellik olcutunun uyusmazlik matrisi gosteriliyor (223 protein
   kendi alt-ailesinin cekirdegi oldugu halde hicbir referansa yakin degil).
   Adaylarin %33'u okaryot.

4. ~~Aramayi guclendir.~~ **BITTI.** `build_search_index.py` → FTS5 (`ro_fts`)
   + filtre tablosu (`ro_search`). Sunucuda ve statik sitede ayni yuzeyler:
   kanit duzeyi, yasam alani, kimyasal aile, tip, plazmit, operon ortagi,
   divergent duzenleyici. Cok kelimeli sorgu, tirnakli ifade ve onek (`naphth*`)
   destekleniyor; sonuclar bm25 ile siraliniyor.

6. ~~Operon kuralinin sinanmasi.~~ **BITTI.** `operon_validation.py`; Atlas
   operon sayfasinda yon/konum testleri, negatif kontrol (baska alfa: p=0,10)
   ve esik duyarliligi tablosu.

7. ~~Tip sayfalarina reaksiyon semasi.~~ **BITTI.** Substrat SMILES'ten ciziliyor,
   ok uzerinde O2 + NAD(P)H, altinda reaksiyon; urun ad olarak veriliyor (urun
   yapisi cizmek hangi halka pozisyonunun saldiriya ugradigini varsaymayi
   gerektirirdi).

8. ~~Reaksiyon gorselleri.~~ **BITTI.** Urun SMILES'i kuratorlemek yerine her
   reaksiyon sinifi icin genel mekanizma semasi cizildi (ana sayfa reaksiyon
   tablosu + tip sayfalari). Substrattan bagimsiz oldugu icin regiokimya
   varsayimi gerektirmiyor.

9. ~~Taksonomi agaci.~~ **BITTI.** `/atlas/taxonomy` icinde katlanabilir NCBI
   soy agaci; cins dugumleri aramaya baglaniyor, kirpilan taksonlar sayisiyla
   belirtiliyor.

10. ~~Per-tip substrat kuratorlugu.~~ **BITTI.** 18. maddede yedi tipin hepsi
    literaturden kaynakli olarak baglandi; bu madde o is bitince kapatilmayi
    atlamis. Tek kalan belirsizlik cadA'nin substratinin 2,4-D mi 2,4,5-T mi
    oldugu ve o soru kullaniciya acik sorular listesinde duruyor.

11. ~~Hugging Face Space (tekrar kayit).~~ 5. madde ile ayni isi tarif
    ediyordu; takip 5'te.

12. ~~Veri kalitesi sayfasi.~~ **BITTI.** `/atlas/quality`: dizi fazlaligi
    (%12,3) ve her sayima etkisi, uc duzeyli test tablosu, olculemeyen
    sinirlarin listesi, ve her indirilebilir dosyanin hangi scriptten geldigi.
    `redundancy.py` + `stats_overview.py`'ye "sequence" duzeyi eklendi.

13. ~~GitHub Pages derleme gecikmesi.~~ **BITTI.** Kok neden deploy
    scriptiydi: her yayinda sifirdan repo kurup 371 MB'lik agaci force-push
    ediyordu, Pages de her seferinde her seyi yeniden derliyordu. Iki duzeltme:
    (a) `deploy_pages.sh` artimli hale getirildi (kalici klon + rsync + normal
    push); (b) giris ve varyant basina FASTA dosyalari statik siteden cikarildi,
    cunku dizi sayfada zaten var ve toplu dosyalar hepsini kapsiyor.
    Site 26.683 dosya / 371 MB'dan **13.457 dosya / 321 MB**'a indi ve yayin
    push ile ayni dakikada tamamlanmaya dondu. Daha fazla kuculme gerekirse
    sonraki adim giris sayfalarini tip basina JSON + tarayicida render etmek
    olur (~2.000 dosya), ama su an gerek yok.

15. ~~RO tiplerinin birlikte bulunmasi.~~ **BITTI.** `cooccurrence.py` +
    `/atlas/cooccurrence`. 854 ciftin 640'i farkli tip; permutasyon null'inda
    en guclu birliktelikler ayni asagi yola besleyen enzimler (BPDO+PhnA1a
    123x, AntA+XylX 5,6x). Iplik bilgisi "tek operon" ile "ayni genom"u ayiriyor.

16. ~~Kimlik → substrat ongorusunun olcumu.~~ **BITTI.**
    `substrate_predictability.py` + kanit sayfasi bolum 3. AUC 0,859 ama hicbir
    esik %90 kesinlige ulasmiyor; sinif duzeyinde lift <=1,2.

17. ~~Surum ve atif bilgisi.~~ **BITTI.** Sayfa altinda derleme tarihi, pipeline
    commit'i ve atif notu (tip atamalari pipeline yeniden kosunca degisebilir).

18. ~~Bilinmeyen substratlar.~~ **BITTI** (2026-10-06). Yedisi de literaturden
    cozuldu, kaynaklariyla `chemistry.csv`'de. cadA "tentative" isaretli.

19. ~~PDB yapilari.~~ **BITTI.** 17 tipe RCSB'den tek tek dogrulanmis yapi
    baglandi (1NDO, 1Z03, 1WW9, 2BMO, 1WQL, 3EN1, 2GBW, 2XR8, 2ZYL, 6Y9C,
    7FHR, 3GKE, 3VCA). Arayuzde 3D gosterim ajan tarafindan ekleniyor.

20. ~~Yeni referans ekleme + modulerlik.~~ **BITTI.** `add_reference.py` ve
    `build_reference_models.py`; motif kolonlari artik veriden cikariliyor.

21. ~~Varyant etiketleri Turkce.~~ **BITTI.** `characterize_leaves.py` ingilizce
    uretiyor, tablo yeniden kuruldu.

22. ~~Izolasyon kaynagi ve habitat.~~ **BITTI.** `isolation_source.py` ile
    9.308 kayittan izolasyon kaynagi, konak, cografya ve toplama yili cikarildi;
    serbest metin 21 kategorili kontrollu bir sozluge normalize edildi.
    Kapsam: girislerin %60,9'unda kaynak var, %55,6'si bir habitate
    yerlestirilebiliyor. Ekoloji sayfasinda tur-bazli normalize tablolar ve
    tip sayfalarinda habitat dagilimi gosteriliyor.
    Sozlukte iki duzeltme yapildi (eleştiri turu): madencilik ayri kategoriye
    cikarildi (asit maden drenaji organik kirlilik degil, sinyali suluyordu) ve
    hipersalin kurali tatli su kuralinin onune alindi ("hypersaline lake"
    icinde "lake" gectigi icin tatli suya dusuyordu). Kapsam %52,2 → %55,6.

23. ~~Veri koken takibi ve hata onlemleri.~~ **BITTI.** `provenance.py`
    (30 tablo, 77 dosya, uretici + girdi + satir sayisi + kaynak veri koku) ve
    `validate_curation.py` (55 kontrol, hata varsa sifirdan farkli cikis, run_all
    icinde HTML uretiminden ONCE kapi olarak). Bulunan ve duzeltilen dort hata
    icin bkz. README "Safeguards".

25. ~~Habitat sozlugunun zayif noktalari.~~ **BITTI.** Uc duzeltme, her biri
    once olculup sonra yapildi:
    (a) `built_environment` eklendi -- hastane lavabosu, temiz oda, uzay araci
        yuzeyi bir HASTA degildir (97 replikon: 45'i yanlisca klinikte, 42'si
        siniflanamamis durumdaydi);
    (b) `human_clinical`'dan belirsiz organ adlari (lung, lymph node, brain,
        spleen) cikarildi -- bunlari tasiyan 46 kayitta hicbir insan isareti
        yoktu ve aralarinda cigerotu simbiyotik dokusu ve mantar dokusu vardi;
        247 giris klinik kategorisinden cikti;
    (c) tek bitki kovasi uce bolundu: `plant_tissue` (endofit, nodul, yaprak
        yuzeyi), `rhizosphere_soil` (kok cevresi, yaprak coplugu) ve
        `plant_associated` (konum belirtmeyen ekin adlari) -- cunku ekin adi
        hangi bitki oldugunu soyler, neresi oldugunu soylemez.
    Kapsam %55,6 → %54,3'e DUSTU ve bu bilincli: yanlis etiketlemeyi birakmanin
    bedeli. Capraz kontrol olarak PAH zenginlesmesi 5,64x → 5,76x, yani
    degisiklikler kirlilik sinyaline dokunmadi.

27. ~~Arayuzde kalan Turkce metin.~~ **BITTI.** Okaryot alt-grup etiketleri
    (`bitki/alg`, `mantar`, `hayvan`, `kirmizi alg`, `diger-ok`)
    `classify_domains.py` icindeki bir sozlukten DOGRUDAN veritabanina ve
    oradan sayfaya gidiyordu. Ilk tarama bunu kacirdi cunku sadece islev
    sozcuklerine bakiyordu; kacan sey bir icerik sozcuguydu. Etiketler
    ingilizceye cevrildi, `ro_domain` tablosu yeniden kuruldu, ve eski rapor
    scriptlerindeki (`make_report.py`) anahtar aramalari da guncellendi --
    aksi halde sessizce sifir okuyacaklardi.

28. ~~Deploy oncesi otomatik site denetimi.~~ **BITTI.** `check_site.py`:
    78 sayfayi acar, 6.217 dahili baglantiyi dener, gorunur metinde ve
    YAYINLANAN veri dosyalarinda Turkce arar, FAIL varsa sifirdan farkli
    cikar. Negatif testi yapildi: duzeltilen etiket hatasini yakaliyor,
    `Vigna radiata var. radiata` ve `protein`/`test`/`once` gibi
    esyazimlilarda yanlis alarm vermiyor.

## Acik sorular (kullaniciya)

- **IKI KURASYON HATASI, duzeltilmedi -- karar kullanicinin.** Yayimlanmis
  yapilarla capraz kontrol (Miao & Schmidt 2025 Tablo 1, 26 yapi, hepsi RCSB'de
  dogrulandi) iki tanesini ortaya cikardi:
  * `3_314_NDO` PDB 1NDO tasiyor ama KENDI referans dizisi o yapiya yalnizca
    %27,5 benziyor. Ayni referans `3_316_NarAa`'nin referansiyla %100 ayni:
    yani bu bir nid/nar ailesi alfa alt birimi, naftalin 1,2-dioksijenaz degil.
    `3_309_NahAc` ve `3_315_NDO` gercekten 1NDO (%100). Hangisi yanlis, isim mi
    dizi mi?
  * `1_115_NdmC`'nin referansi, `1_114_NdmB`'nin 373 aa referansinin 355 aa'lik
    TAM alt dizisi -- ayni PDB yapisinin His6-TEV etiketi (MGSSHHHHHHENLYFQGS)
    kirpilmis hali. NdmC gercek ve ayri bir N7-demetilaz, yani bu girdi yanlis
    protein. NdmC'nin neden hic uyesi olmadiginin da aciklamasi bu.
- **Referans setinin kendini tekrar etmesi olculdu.** Neredeyse ayni cift
  referanslar: OxoO/OMO %100, NahAc/NDO(3_315) %100, EdoA1/cumA1 %99,8,
  CARDO/CarAa %99,2. Birlestirilsinler mi, yoksa iki ad da korunsun mu?
- `5_504_CntA` 6Y9C tasiyor; o yapi **C229A nokta varyanti**. Yabani tip 6Y8J.
  Aktif bolge referansi olarak hangisi kullanilsin? (6Y9C'de mononukleer demir
  var, 6Y8J'de yok -- ikisini birden listelemek en dogrusu gorunuyor.)
- Yapisi HIC olmayan bes tip kaldi: `1_109_LigX`, `3_312_NinA`, `5_507_BmoA`
  (AlphaFold kaydi yok ya da UniProt kimligi yok), `2_215_cadA` ve
  `5_505_qxyA` (bunlarin AlphaFold modeli VAR, surum eki hatasi duzeltildi,
  `active_site.py` yeniden kosturulunca gelecekler).
- `chemistry.csv`'de **product_smiles sutunu hic yok**: 71 tipin hicbirinde
  urun yalnizca Ingilizce metin, yani hicbir reaksiyon mekanik olarak
  denklestirilemiyor ya da cizdirilemiyor. Eklensin mi?
- `curation_confidence` 71 tipin 70'inde "high"; alan hicbir ayrim tasimiyor.
  Yeniden derecelendirilsin mi?
- ~~cadA'nin substrati 2,4-D mi, 2,4,5-T mi?~~ **KAPANDI (2026-10-28, madde
  75).** Soru yanlis makaleye dayaniyordu: bizdeki dizi Bradyrhizobium HW13
  CadA DEGIL, Sphingomonas sp. ERG5 CadA (GenBank AGU69349.1, dizi birebir) ve
  o kayit "large subunit of 2,4-D oxygenase CadA" diyor. Substrat 2,4-D olarak
  DOGRU; duzeltilmesi gereken organizma ve atif.
- Varyant kumelemesi ham dizide tekrarlanacak mi (madde 1)? Site genelinde
  sayilar degisir.
- Site Ingilizce; Turkce surum de istenir mi?

## 90. Katalitik karboksilat YANLIS KOLONDAN okunuyordu (2026-10-30)

Olculen: kolon 431 dogru, kolon 471 (kullanimdaki) ve 355 (eski gomulu sabit)
yanlis. Kanit: 1NDO demir ligandlari His208/His213/Asp362; refs71'deki NahAc
dizisi PDB numaralandirmasindan 3 kalinti geride basliyor (ham 202/205/210 =
PDB 205/208/213) ve ham 359 = PDB Asp362 match kolonu 431'e dusuyor. 76
referansta kolon 431'de D 68; 19.114 onaylanmis giriste +-2 pencereyle %87,5
D/E (kolon 471'de %71,1, ve E/D karisimi literaturun "genellikle Asp"
tarifine uymuyor).

Kok neden: build_reference_models.py kolonu D/E korunmasini maksimize ederek
seciyordu ve dogru kolonu 3 referans farkla kaybediyordu (471 %93,4 vs 431
%89,5). "En korunmus D/E kolonu" ile "demiri baglayan kolon" ayni sey degil.

Duzeltildi: kolon artik yapisal capadan (STRUCTURAL_ANCHOR) okunuyor, capa
dogrulanamazsa betik duruyor; motif_columns.json yeniden uretildi ve
carboxylate_source ile kolonun neyle secildigini kaydediyor; ro_motif.py'deki
eski gomulu sabitler duzeltildi; validate_curation.py'ye iki kontrol eklendi
(107 kontrol, 102 PASS / 1 FAIL / 4 WARN -- FAIL bilinen bayat ro_search).

Metot kendisi dogrulandi: carboxylate_anchor_holdout.py, disarida birakilan
pozitiflerde 490/490 geri kazanim, 10/10 tipte ayni kolon. Yani metot saglam,
HEDEFI yanlisti.

GECERSIZ: carboxylate_audit.json'un butun sayilari (1.414 / 1.388 / 26,
%12,4), ona dayanan klad kaymasi ve pencere genisletme karsilastirmalari,
kapinin "triad 3/3" istatistigi (13.543) ve ilgili site metni.
GECERLI KALAN: kopru karboksilat (kolon 286 dogrulandi, 8/9 yapisal olcum),
ligand/role ayrimi, oksijen kaderi ve metabolik rol.

SIRADA: (a) rebuild.py --run --include-heavy ile ro_carboxylate ve
catalytic_intact'i kolon 431 ile yeniden uret, (b) carboxylate_audit.py'yi
yeniden kosup metot kartini yeniden uret, (c) kalan gercek yokluklari incele
(referans PobA kolon 431'de K, KshA15 boskuk, 3_309_NahAc hic uye toplamiyor),
(d) site metnini guncelle, (e) yayinla.

Rapor: RAPOR_karboksilat_kolonu_2026-10-30.md
Atif eklendi: barry2013_mechanism (10.1021/cs400087p, Crossref ile dogrulandi;
119/119 atif cozuluyor).

## 91. Kapsama: DB, kendi ilk adiminin buldugu RO alpha'larin %21'ini tutuyor (2026-10-30)

Kullanici sordu: "ilk adimda NCBI'daki tum proteinler yok gibi". DOGRU, ve
olculdu. Huni:

    189.657  PF00355 Rieske proteini -- UNIPROT anlik goruntusu (11 Ara 2024),
             NCBI nr/RefSeq DEGIL
    115.864  71 profil kutuphanesine herhangi bir domain isabeti olan
             (geri kalan 73.793 hic isabet vermedi; ro_domtbl.out ile dogrulandi)
     73.699  RO alpha (68.747 tekil dizi)   | 40.687 alpha degil | 1.478 parca
     26.693  genomik baglami cikarilabilen aday (GenBank gbk gerekiyor)
     19.114  onaylanmis giris (15.953 tekil dizi)

Ortusme olcumu (birebir dizi karsilastirmasi):
    14.449  DB dizisi protein-duzeyi sette de var
    54.298  protein-duzeyi RO alpha dizisi DB'DE YOK  <-- kapsama acigi
     1.504  DB'de olup UniProt anlik goruntusunde olmayan (genom taramasindan)

Yani yayinlanan DB, kendi ilk adiminin RO alpha dedigi tekil dizilerin
%21'ini (14.449/68.747) tutuyor. Sebep yontemsel: DB genom kayitlarindan
kuruluyor, bir protein ancak annotasyonlu bir GenBank kaydina baglanabildiyse
girebiliyor. Ayrica evren UniProt anlik goruntusu oldugu icin yalnizca NCBI'da
olan proteinler hic aday olmadi.

SORUN: about.html bunu SOYLEMIYOR. Akis semasi "Genome records 17.074
annotated replicons from GenBank" ile basliyor; yukari akistaki UniProt
PF00355 adimi, tarih (Ara 2024) ve 73.699 -> 19.114 kaybi yayinlanmis
metinde yok. Bu bir kapsama iddiasi degil ama kapsama IZLENIMI veriyor.

YAPILACAK:
  (a) about.html'e gercek huniyi ve tarihi yaz (ucuz, zorunlu).
  (b) Kapsama genisletmesi (pahali, kullanici karari): ya protein-duzeyi
      54.298 diziyi genom baglami OLMADAN ikinci sinif kayit olarak al, ya da
      NCBI nr/RefSeq uzerinden yeni bir tarama yap. Ikisi de DB'nin butun
      istatistiklerini degistirir.

## 92. Karboksilat kapida ZORUNLU DEGIL -- ama literature gore kritik (2026-10-30)

Kullanici sordu: "Fe karboksilat kritik, dikkate aliyoruz di mi". Cevap:
okunuyor ve raporlaniyor, ama kapi triad'in >=2/3'unu istedigi icin
ZORUNLU DEGIL -- bir giris iki His ile, karboksilatsiz gecebiliyor.

Literatur: NDO'da Asp362'nin alanine mutasyonu (D362A) enzimi TAMAMEN
INAKTIF yapiyor (oksijen tuketimi ve 14C-naftalin urun olusumu, ikisi de
sifir), ve yapida Asp362 demiri bidentat baglar. Yani kalinti katalitik
olarak kritik; "2/3 yeterli" biyokimyasal bir tolerans DEGIL, olcum
araciyla ilgili bir tolerans.

Duzeltilmis kolon 431'de +-2 pencereyle 19.114 onaylanmis giris:
    16.726  karboksilat VAR (%87,5)
     1.947  pencere tamamen bosluk    -> hizalama siteye ULASMIYOR (artefakt)
       220  merkez bosluk, komsular dolu
       221  gercek kalinti var ama D/E degil (N 76, T 41, A 29, G 16, S 16...)
             -> tek BIYOLOJIK iddia bu kucuk kume

Yani 3/3 zorunlu olsa 2.388 giris duserdi, ama bunun ~2.167'si hizalama
erisimi sorunu; gercek ikame ~221. Dogru sira: capali okumayi (metot
dogrulandi, 490/490) kolon 431'de kos, kalan gercek yokluklari belgele,
SONRA 3/3'u zorunlu yapmayi tart. Once rebuild gerekiyor.

## 93. "Literaturde bes grup var" -- hangi bes, ve bizdeki yedi ile iliskisi (2026-10-30)

Kullanici: "bence 5 RO grubu var, literatur de oyle diyor". Olculdu ve
kaynaklandi: literaturde IKI ayri bes-katli sema var ve ikisi de ELEKTRON
ZINCIRININ mimarisine gore:
  Batie ve ark.        bes sinif (IA, IB, IIA, IIB, III)
  Kweon ve ark. 2008   bes tip   (Tip I-V), doi:10.1186/1471-2091-9-11
Kweon semasi Batie'nin yerine onerildi ([3Fe-4S] Fd'yi Batie aciklayamiyor).

Bu veritabaninin `ro_group` alani (1-7) BU EKSEN DEGIL: alfa alt biriminin
dizi/hizalanabilirlik kumelenmesi. Kanit zaten batie_classes.yaml'da yaziliydi
-- Batie IB'nin iki ornegi bizde grup 1 ve grup 4'te (OMO -> 1_101_OxoO,
BZDO -> 4_411_XylX), ve CARDO Batie tablosunda uc ayri sinifta geciyor (sus
bazli). Yani yedi grubu bese indirmek literature yaklasmak degil, iki ekseni
birbirine karistirmak olur.

YAPILABILECEK SEY: Kweon tipini AYRI bir alan olarak eklemek. ro_etc tablosu
zaten reduktaz (FNR/GR) ve ferredoksin (rieske/plant) tipini tasiyor. Ama
atanabilirlik olculdu (19.114 giris):
    %73,5  atanamaz -- reduktaz operonda yok
    %19,6  Tip I veya II (FNR_C / FNR_N ayrimi KURULU DEGIL)
    % 5,4  Tip IV veya V ([2Fe-2S] / [3Fe-4S] ayrimi KURULU DEGIL)
    % 0,8  Tip III -- tek basina cozulebilen tip
    % 0,5  belirsiz (FNR+GR)
Yani once component_hmms'e iki yeni ayrim gerekiyor: FNR_C vs FNR_N reduktaz
ve [2Fe-2S] vs [3Fe-4S] ferredoksin. O ikisi kurulsa bile alanin %73'u bos
kalir, cunku operon ortaklari cogu girişte yok. Eklenirse "olculemedi" ile
"yok" ayri ayri yazilmak zorunda.

Atif eklendi: kweon2008_classification (Crossref ile dogrulandi).
Guncellendi: batie_classes.yaml -- iki bes-katli sema, Kweon tip tanimlari ve
atanabilirlik olcumu; ayrica "ro_group 1-5" ifadesi bayatti (artik 1-7).

## 94. 3/3'un gercek bedeli olculdu: 31 giris (2026-10-30)

Denetim DUZELTILMIS kolonda (431) yeniden kosuldu. Sonuc:

  sabit kolon (431, +-2) ile 3/3 kurali : 2388 giris (%12,49), 45 tip
  tip ici CAPALI okuma ile 3/3 kurali   :   31 giris (%0,16),  12 tip
  capalamadan sonra GERCEK yokluk       :   24 giris, 11 tip
        oturan kalinti: 23 bosluk + 1 N
        tip: PTC52 5, CAO 5, DAF36 4, Nvd 3, TsaM 1, BmoA 1, CmoS 1, Stc2 1
        taksonomik alan: Okaryot 20, Bakteri 4

Yani sabit kolonda 3/3 uygulamak, karboksilati FIILEN OLAN 2357 gercek
RO'yu atardi; en sert vurdugu tipler de kuratorlu referansa en uzak
olanlar (AerC %100, Nvd %82, DAF36 %54, PobA %48, PrnD %45, DdmC %41).

KARAR: 3/3 kanonik tanim olarak DOGRU, ama SABIT KOLON okumasiyla
uygulanamaz. Sira: (a) ro_carboxylate'i tip ici capali okumayla yaz,
(b) sonra MIN_CATALYTIC_SITES=3 yap, (c) 24 gercek yoklugu motif
istisnasi olarak tek tek belgele. Bu haliyle hem kanonik tanimla
ortusur hem gercek RO atmaz.

Not: denetimin verdict metninde "referans setinde okaryotik RO tipi YOK"
yaziyordu -- BAYATTI, 6_601_PAO / 6_602_CAO / 6_603_PTC52 14 Ekim'de
eklendi ve 617 uye topluyor. Metin duzeltildi; kalan yokluklarin hala o
tiplerde olmasi referans yoklugu degil, o alt ailelerin katalitik lobunun
bakteriyel profille hizalanmamasidir.

## 95. Protein evreninin kaynagi: UniProt, NCBI DEGIL (2026-10-30)

Kullanici sordu: "protein dizilerini tum NCBI'dan mi aldin, nasil yaptin".
Cevap: NCBI'dan HICBIR SEY indirilmedi. Evren depodaki combined_pfam.fasta:
189.657 dizi, 189.397 TrEMBL + 260 Swiss-Prot basligi, indirme tarihi
11 Ara 2024. Olculdu (4.000 dizilik ornek): 3.999'u PF00355'i kendi GA
esigiyle tutuyor, yani dosya gercekten bir PF00355 setidir.

Kapsami artirmak icin olculenler:
  1. GIRIS KAPISI DARBOGAZ DEGIL. Dongusel olmayan test: UniProt setinde
     OLMAYAN, genom taramasindan gelen 1.504 onaylanmis dizinin 1.504'u de
     PF00355 ile kendi esiginde toplaniyor.
  2. UNIPROT'U TAZELEMEK KAYBETTIRIR. Bugun PF00355 = 100.453 protein
     (InterPro) / 101.576 (UniProt API) -- bizim 2024 anlik goruntusu
     189.657. UniProtKB aradan fazla kayitlari sildi. Yani yeni cekim
     ESKISININ YERINE DEGIL, BIRLESIMI olarak kullanilmali.
  3. IKINCI GIRIS KAPISI AZ AMA GERCEK KAZANC. PF00848 (Ring_hydroxyl_A,
     katalitik alan) = 22.162 protein; 20.909'u PF00355 de tasiyor. Yani
     1.253 protein katalitik alani tasidigi halde Rieske anotasyonu YOK --
     mevcut giris kapisi onlari hic goremiyor.
  4. NCBI. UniProtKB zaten INSDC/RefSeq cevirilerini aliyor, ortusme buyuk.
     Gercekten fazla olan kisim: NCBI'da olup UniProtKB'de olmayan (ya da
     silinmis) proteinler. Tam nr taramasi yuzlerce GB; pratik yol alan
     kisitli bir NCBI sorgusu (CDD/PF00848) + ayni kapi.

Oncelik sirasi: (a) elimizdeki 54.298 baglamsiz protein girisini al (kosuyor),
(b) PF00848'i ikinci giris kapisi yap (+1.253 aday, ucuz), (c) yeni UniProt
cekimini BIRLESIM olarak ekle, (d) NCBI taramasi -- disk ve karar gerektirir.

## 96. Triad kapisi literaturun tanimina getirildi (3/3) -- 2026-10-30

Literatur mononukleer Fe merkezini ISTISNASIZ 2-His-1-karboksilat olarak
tanimlar ve karboksilat sus degildir: NDO'da D362A TAMAMEN inaktiftir
(Barry & Challis 2013, ACS Catal 3:2362). Kabul kurali ise >=2/3 idi.

`triad_gate.py` triad'i uc ligandin HEPSINI isteyerek cagirir ve her ligandi
O LIGANDIN okunabildigi yoldan okur:

    karboksilat   sabit kolon +-2; kolon bossa tip-ici capali okuma (adim 6)
    His-1, His-2  SADECE sabit kolon +-2  <-- capali okuma henuz YOK

Olculen (19.114 onaylanmis giris, ro_triad tablosu):

    3/3 tam                      19.028
    karboksilat yok (capalandi)      24   -> 11 tip, 20'si okaryot
    karboksilat test EDILEMEDI        7   -> 7_701_AerC; tipin kendi
                                            pozitifi yok, capa kurulamadi
    His karari verilemez             55   -> CntA 28, TsaM 12, CmtAb 8, ...

HIS OKUMASI AYNI HATA SINIFINI TASIYOR, ama pencere genisletmek COZUM DEGIL.
Olcum: sabit kolonda okunan 19.059 giriste His-1 -> His-2 arasindaki kalinti
sayisi 4 (18.445) ya da 3 (607) -- literaturun H-x(3-5)-H sarti. Pencere
+-5'e cikarilinca 55 girisin 55'inde "His bulunuyor", ama bulunan ciftin
arasi 0 (27 giris) ya da 7 (27 giris): pencere bitisik bir HH'yi yakaliyor,
triad'in His'ini DEGIL. Ham dizide D-x(2)-H-x(3-5)-H konsensusu arandiginda
CmtAb'nin 8 girisinde eslesme VAR (ornek: raw 204, ELYHEWGHIV) -- yani His
cifti gercekten orada, sadece 5 kolon kaymis; TsaM/VanA/CdnA/TdnA1/Stc2/
PAO/Nvd'nin 19 girisinde dizinin HICBIR yerinde eslesme YOK.

KARAR. Kural 3/3'tur ve `ro_triad` bu tanimi yayinlar. `ro.is_confirmed`e
DOKUNULMADI: 62 giris (55 His + 7 AerC) henuz capalanmamis bir okumaya
dayaniyor ve hicbiri "triad yok" diye yayinlanamaz. Kapi kontrolu bunu WARN
olarak acik tutuyor (`validate_curation.py: check_triad_gate`).

### 96a. His okumasi da capalandi -- `triad_anchor_his.py`

His-1 ve His-2 artik karboksilattakiyle AYNI metotla okunuyor: tip-ici mafft,
kolonu tipin kendi pozitifleri oyluyor, yontem tutulan yarida dogrulaniyor
(9 tipin 9'unda tutulan yari 1.0, pozitif geri okuma %99,2-100).

55 kararsiz girisin sonucu: 35'inde iki His de tipin kendi kolonlarinda VAR,
20'sinde o kolonlarda YOK. Yani +-5 pencerenin "55/55 His buldum" demesi
tamamen yanlisti.

KURALI BELIRLEYEN VAKA -- 5_504_CntA. 28 girisin 27'sinde iki His de capa
kolonlarinda (371 ve 379) oturuyor; o kolonlar 120 pozitifin 120'sinde His.
Ama aradaki ara 7 kalinti, 3-5 degil: bu alt-varyantta demir baglayan ilmekte
UC KALINTILIK BIR INSERSIYON var (GYH-FSLSGP-CH, pozitiflerde CYH---CGPAH).
Aralik sarti KARAR verseydi 27 gercek enzim "triad eksik" diye atilacakti.
Literatur iki His'in demiri baglamasini ister, ilmek uzunlugunu sarta
baglamaz. Karar bu yuzden KOLON HOMOLOJISI: His TAM capa kolonunda olmali;
kaymis bir okuma ancak H-x(3-5)-H araligi destekliyorsa cagri uretir. Bu ek
sart olmadan TsaM'in 2 ve PAO'nun 1 girisinde pencere yine bitisik His'i
kabul ediyordu (capa kolonunda N/D/R var).

AerC (7 giris) BILEREK test edilmemis kaliyor. Referansa capalama denendi ve
yapilamadi: `7_701_AerC_AerC_pro` referansinin KENDISINDE de kolon 431'de
karboksilat yok, hatta katalitik His kolonlarinda Q ve D var. Sebep 2026-10-14'te
cozulmus (`aerc_resolution.json`): AerC 717 aa, Rieske ile katalitik lob
arasinda buyuk bir insersiyon var ve MODELIN 429-570 BOLGESININ AerC'DE HIC
KARSILIGI YOK. Okunamayan bir kolonda "yok" demek olcum degil, kurgu olur.

Ayni olcum referans setinin kendisini de vuruyor: 76 referansin 6'sinda
(PobA, LigX, KshA15, Nvd, DAF36, AerC) sabit kolonda karboksilat YOK. Yani
sabit-kolon okumasi uyelerde degil, KURATORLU REFERANSLARDA bile basarisiz.

ro_triad'in son hali: tam 19.063 | karboksilat yok 24 | His yok 20 |
test edilemedi 7 (AerC).

Siradaki (bu sirayla):
  1. `annotate_ro.py`'de MIN_CATALYTIC_SITES = 3 ve `rebuild.py --run
     --include-heavy`. O noktada is_confirmed 19.114'ten dusecek ve TUM
     denominatorler degisir, yani ayni kosuda site yeniden donmeli.
     KARAR GEREKTIREN NOKTA: 7 AerC girisi okunamiyor; "yok" degil
     "okunamadi" oldugu icin atilmamalari gerekir, yani kural "triad 3/3,
     ya da triad'in okunamadigi BELGELENMIS mimari" olmali.
  2. Geriye kalan 44 gercek yoklugu (24 karboksilat + 20 His) tek tek motif
     ISTISNASI olarak belgelemek: hicbiri "RO degil" demek degil, Rieske 4/4
     hepsinde tam.
  3. Protein-duzeyi girisler icin de capali okuma kurmak (asagi, 97).

## 97. `ro_protein` kuruldu -- genomik baglami olmayan 73.699 RO alpha

Kullanici karari: "annotasyonlu genbank olmasa da protein olmali bizim
db'de". `ingest_protein_entries.py` kosuldu (hizalama 657 MB, onceden
tamamlanmis, yeniden hizalanmadi).

    protein-duzeyi giris                      73.699
    Rieske 4/4 + triad >=2/3 (mevcut kapi)    72.586
    Rieske 4/4 + triad 3/3 SABIT KOLONDA      62.063   <- SISMIS, bkz asagi
    genom-baglamli bir girisle AYNI dizi      17.452
    yeni, baglamsiz dizi                      56.247

AYRI TABLO, `ro`'ya ek satir DEGIL: sitedeki her giris-bazli istatistik
(plazmit orani, operon kurallari, komsuluk zenginlesmesi, izolasyon kaynagi,
regulasyon) bir genom konumu varsayar. Baglamsiz satirlar `ro`'ya eklenirse
paydalar 3,4 katina cikar ve paylar degismez. Ayri tabloda baglam isteyen
bir sorgu bu satirlari BULAMAZ, yani yanlis sayi yerine eksik sayi uretir.

62.063 sayisi KAPI DEGERI OLARAK KULLANILMAMALI: sabit-kolon okumasidir ve
dolayisiyla sismis bir yokluk bildirir -- `ro`da ayni okuma 2388 yokluk
dedi, capali okuma 24'e dusurdu. Protein-duzeyi girislerde capali okuma
henuz YOK; tip basina pozitif capa kurulmasi gerekiyor. 3/3 kapisi
protein tablosuna uygulanmadan once bu yapilmali, yoksa ~10.000 gercek
enzim sessizce elenir.

Yeni hata, ayni sinifin dordunculusu: `has_carbox = (carbox or "") in "DE"`
Python'da BOS okumayi "var" sayar (`"" in "DE"` True'dur) ve AerC'nin 7
girisini "3/3 tam" isaretledi. Kume tipine cevrildi; uyelik testleri artik
`frozenset` uzerinden.

## 98. AGIR YENIDEN KURULUM KOSUYOR -- 2026-10-30 18:40

NEDEN. `rebuild.py --check`: 51 artefaktin 47'si BAYAT, ve en onemlisi
`annotate` bayat -- cunku `ro_motif.py` 2026-10-30'da degisti (karboksilat
kolonu 471 -> 431). Yani `ro.is_confirmed` ve `ro.catalytic_intact` hala
2026-10-22'de, YANLIS KOLONDA uretilmis degerler. Olculen etki +12 / -1
giris. Site'i bu haliyle dagitmak 13 girisi yanlis yayinlamak olurdu.

    nohup python3 rebuild.py --run --include-heavy --cpu 16
    log: /media/lin-bio/back2/roar_rebuild_2026-10-30.log

BITTIKTEN SONRA SIRAYLA YAPILMASI GEREKENLER. is_confirmed degisince
ro_triad'in dayandigi kume de degisir, yani triad kapisi yeniden
kurulmalidir -- ve kendi tazelik kontrolu bunu zaten FAIL olarak yakalar
(`check_triad_gate`: denetim kac giris gormus vs DB'de kac var):

  1. python3 carboxylate_audit.py --work /media/lin-bio/back2/roar_carbox_work
     (yeni onaylanmis kume uzerinde; tip basina mafft, uzun surer)
  2. python3 triad_anchor_his.py --work /media/lin-bio/back2/roar_triad_anchor
     (--reuse-alignments KULLANMA: girdi degisti, yeniden hizalanmali)
  3. python3 triad_gate.py --apply
  4. python3 ingest_protein_entries.py --work /media/lin-bio/back2/roar_protein
     (linked_candidate_id yeni kumeye gore yeniden kurulmali; hizalama
     yeniden kullanilir, hmmalign tekrar kosmaz)
  5. python3 validate_curation.py    -> 0 FAIL olmali
  6. python3 check_site.py           -> "Ready to deploy" olmali
  7. cd webapp && python3 freeze.py --out site --base /pro-sim-blast
     && bash deploy_pages.sh git@github-prosimblast:recepcanaltinbag/pro-sim-blast.git

KAPI 3/3'E BU KOSUDA SIKILASTIRILMADI. Sebep: kuralin dogru yazimi
"triad 3/3, YA DA triad'in sabit kolonda okunamadigi BELGELENMIS mimari"
ve `annotate_ro.py` bu ikinci sarti goremez (capali okuma ondan SONRA
kuruluyor, dongusel). Duz `MIN_CATALYTIC_SITES = 3` yazilirsa AerC'nin 7
girisi "yok" denmeden elenir. Dogru cozum: annotate sabit-kolon 2/3'u ON
FILTRE olarak tutar, uyelik son karari ro_triad'dan gelen ayri bir pastan
verilir. O pas yazilmadi; madde 96/96a/97'deki olcumler hazir.

## 119. A5 PLANI: protein duzeyindeki 18.636 diziye genomik baglam (2026-11-02)

Kullanici "planlayip yapabilir misin" dedi. Plan, once OLCULEN uc sayinin
uzerine kuruldu.

OLCUM (2026-11-02):
  - `ro_union` 73.299 tekil dizi; genomik baglam tasiyan 15.960, tasimayan
    57.339.
  - Baglamsiz olanlardan **18.636**'sinin UniProt capraz referansi bir
    EMBL/GenBank NUKLEOTID kaydi veriyor. Yani baglam KURULABILIR.
  - Bu 18.636 dizi **13.079 TEKIL kayda** isaret ediyor (11.037 WGS kontig,
    2.042 digeri). Indirme birimi dizi degil KAYIT olunca is %30 kuculuyor.
    En yogun kayit CP000699 tek basina 52 dizi tasiyor.
  - Hicbiri `replicon` tablosunda YOK: kesisim gercekten buyur.

IKI TASARIM KARARI, IKISI DE OLCUMLE VERILDI.

(a) PENCERE INDIRME REDDEDILDI. Kayitlarin ortalamasi 10 MB (ornek n=24) ve
baytin ~%95'i ORIGIN dizisi; ihtiyacimiz olan sey CDS cevirileri. CDS'in
±10 kb penceresini cekmek isi 132 GB'dan ~1 GB'a indiriyordu. REDDEDILDI,
cunku iki sessiz hata uretiyor: (1) pencere kaydinin koordinatlari pencereye
gore yeniden numaralanir, yani `candidate_id` ve operon koordinatlari DB'nin
geri kalaniyla kiyaslanamaz hale gelir; (2) komsuluk penceresi
NEIGHBOR_WINDOW ile ayni genislikte oldugu icin kenardaki komsular KESILIR
ve `has_beta` / `has_ferredoxin` / `has_reductase` sistematik olarak EKSIK
olcuulur -- yani operon tamligi, olcum hatasi yuzunden dusuk yayimlanirdi.
Tam kayit indirmek bu iki riskin ikisini de tasimiyor ve ORIGIN blogu
elimizde kaldigi icin promotor analizi ikinci bir indirme gerektirmiyor.

(b) KOHORT ALANI ZORUNLU -- ve bunu bir sayfa metni soyletti. `about.html`
"17.074 annotated replicons scanned, N carrying at least one" yayimliyor.
Bu bir YAYGINLIK iddiasi. 13.079 yeni kayit, tam da ICINDE RO OLDUGU
BILINDIGI ICIN secildi; ayni havuza katilirsa o oran ~%100'e cikar ve bir
kesif gibi gorunur. Oysa olcum degismis olmaz, ORNEKLEME degismis olur.
Bu yuzden A5 kayitlari AYRI BIR KOHORT olarak isaretlenir ve genom duzeyinde
her yaygınlık olcumu hangi kohortu kullandigini BEYAN etmek zorundadir.
`null_model.py` bu riskten muaf: rastgele pencereleri AYNI replikon icinden
orneklediği icin karsilastirma kohortlar arasi degil kohort ici.

ADIMLAR (her biri kosturularak dogrulanir):
  A5.1  `a5_fetch_protein_context.py` -- 13.079 kaydi ceker. Yeniden
        calistirilabilir, hiz sinirli, her kaydin sonucu `ledger.jsonl`de.
        Gecici hata ile kalici "yok" AYRI sayilir (madde 116'nin dersi).
  A5.2  `extract_genomic_context.py` ayni dizine kosturulur -- YENI KOD YOK,
        ayni cikarici, ayni CSV sozlesmesi.
  A5.3  DIZI ESLESME SARTI: bir kayit, bir diziye ancak cevirisinin md5'i
        BIREBIR esitse baglam olur. Capraz referans "ayni kayitta gecer"
        demektir, "ayni proteindir" DEMEZ. Eslesmeyenler baglamsiz kalir ve
        kac tane olduklari yayimlanir.
  A5.4  `replicon_cohort` TURETILMIS TABLOSU (sutun degil): hangi kayit hangi
        kohorttan, ve A5 kayitlari hangi dizi icin cekildi.
  A5.5  Kapi kontrolleri: (1) kohort beyan edilmeden genom yaygınlığı
        yayimlanirsa FAIL; (2) md5 eslesmesi olmayan bir dizi baglamli
        sayilirsa FAIL; (3) `about.html` orani yalnizca ozgun kohorttan
        hesaplanmazsa FAIL.
  A5.6  Zincir yeniden kosturulur, sayfalar guncellenir, dagitilir.

## 118. HcaE'de "entry duplication" SORULDU: bug yok, ama BIRIM YAZILI DEGILDI (2026-11-02)

Kullanici: "mesela HcaE icin neden entry duplicationlari var bug mi".

**Bug yok.** 259 girisin 258'i AYRI genom kaydi, ve ayni genom+koordinat
ikilisinin tekrari **0**. Tekrar eden sey giris degil, iki baska sey:

1. **Ayni protein farkli suslarda.** HcaE *E. coli*'nin 3-fenilpropiyonat
   dioksijenazi; tek bir dizi 46 genomda birebir ayni. 259 giris 156 tekil
   diziye cokuyor.
2. **RefSeq `WP_` kayitlari tasarimi geregi tekrar eder.** `WP_` accession'lari
   non-redundant: ayni proteini tasiyan butun genomlar icin TEK accession.
   HcaE girislerinin %98,8'i `WP_` tasiyor.

Yani birim **bir genomik lokus**, bir protein degil. Tum set icin ayni
ayrisma: 19.125 giris, 11.610 tekil genom kaydi, 15.960 tekil dizi, ayni
lokus tekrari 0.

**ASIL KUSUR SAYFADA.** Birim yaziliydi degil, o yuzden okuyucu hata
sonucuna vardi -- projenin tekrar eden sinifi: *birim yazilmazsa okuyucu
kendi sonucunu uretir* (ayni sinif 2026-11-01'de birlesim sayilarinda da
cikti). Tip sayfasina olculen not eklendi: giris / tekil genom kaydi /
tekil dizi / en yayginin kac genomda oldugu / WP_ orani / ayni lokus
tekrari. Sayilar `app.py`'de olculuyor, sablonda yazili degil.

ACIK IS: substrat basina bolge tutarliligini sinayan bir kontrol yok --
PDO/PhtAa hatasi (madde 114) tam boyle bir celiskiydi ve iki Gram-negatif
ftalat dioksijenazina farkli bolge atanmasini hicbir kapi yakalamadi.

## 117. ANASAYFA: once ONEM, sonra sayilar; 11.481 kelimelik aciklama kirpildi (2026-11-02)

Kullanici, dort mesajda: "anasayfada cok yazi var daha duzenli ve goze hitap
eden bi tasarim olamaz mi" / "ozellikle 1.2 ve ilk giris kismi, hic iyi bi
database gibi durmuyo" / "kalabilirler ama daha ozet gibi olmalila" /
"sayfayi ilk acan direkt bisuru yaziyla karsilasiyor ya tablonun altin bi
yazi var inanilmaz uzun aciklama vs" / "rieskelerin onemli oldugunu insanlar
direkt anlamali ya anasayfa nasil yapilir biliyor olman lazim".

OLCUM: anasayfa 3.933 kelime (1. bolum 1.477, giris paragrafi 256). Site
genelinde 45 kelimeden uzun 131 aciklama blogu, toplam 11.481 kelime, en
uzunu 283.

**Hicbir olculmus cumle silinmedi** -- metnin hepsi bir olcumun sinirini
yaziyor. Degisen sey nerede ve ne zaman gorundugu.

- **Onem yayinda ve OLCULUYOR.** Yeni h1 + uc kart (9.430/53, 3.487/30,
  9.695/23). "Hangi aileler kirletici sayilir" bir kurator karari, o yuzden
  `significance.yaml`'da -- dosya TEK SAYI tasimiyor, sayim her derlemede
  DB'den yapiliyor, yani bayatlayamaz. 3 yeni kontrol.
- **Sira degisti**: baslik -> neden onemli -> arama -> sayilar ->
  kesisim/birlesim -> mimari ozeti. Ilk ekran 449 -> **252** gorunur kelime.
- **1. bolum ozete indi**, gerekcesi `/architecture`'a tasindi. Icerigin
  baska hicbir sayfada OLMADIGI olculdu, yani silinebilir degildi. Sema iki
  yere kopyalanmadi: `_arch_svg.html` olarak paylasiliyor.
- **Uzun aciklamalar kirpildi** (`longnote.js`): 45+ kelime -> 3 satir
  (anasayfada 2) + "Show the full note". HICBIR DUGUM TASINMIYOR; ilk
  tasarim metni cumleden kesip tasiyordu ve bu icerideki `<a>`/`<b>`
  etiketlerini bolup bozuk HTML uretebiliyordu. `line-clamp` kullanilmadi,
  cunku gerektirdigi `display:-webkit-box` blok cocuk tasiyan bloklari
  bozuyor. JS yoksa ve yazdirmada metin TAM.

KIRILGAN KAPI: `check_union_denominator` `"is an\n<strong>intersection
</strong>"` ariyordu, yani SATIR SONUNA bagliydi ve panel bicimine gecince
kirildi. Amaci bicim degil ayrimin yazili olmasi -- esnetildi.

Dagitildi: 3.210 degisen dosya. Kapilar: 150 kontrol / 146 PASS / 0 FAIL,
rebuild 61/61, check_site 85 sayfa + 67.411 statik baglanti temiz.

## 116. DAGITIM: iki kez yanlis yere derledim; ve 21 sorgu hatasi raporu temiz gosterdi (2026-11-02)

Kullanici: "sitenin istatik kismi ya da diger bolumleri vs guncel mi knk baya
bisi var cunku". Olculdu: **degildi**. Son dagitim 16:48'deydi ve o andan
itibaren 10 commit birikmisti; yayindaki sitede `/atlas/groups` 404,
birlesim indirmeleri 404, PDO/PhtAa hala ters.

Dagitimi UC ayri hata geciktirdi ve ucu de kayda deger:

**(a) On plan timeout'u FastTree'yi oldurdu ve agac 0 bayt kaldi.**
`timeout 599 python3 rebuild.py --run` ile kostum; timeout `phylogeny`
adimini ortasinda kesti. `analysis_out/tree_all.nwk` VARDI ve mtime'i
TAZEYDI, yani tazelik takibi adimi guncel saydi; zincir ilerledi ve
`group_monophyly` bos dosyayi ayristirip IndexError ile patladi. Hata
URETEN adimda degil OKUYAN adimda cikti, yani kirilan yerde sebebi
gorunmuyordu. Kalici duzeltme: beyan edilen bir cikti 0 bayt ise
URETILMEMIS sayilir (`rebuild.py staleness()`). Agac 328 saniyede yeniden
kuruldu -- yani korkulan iki saat degildi, bedeli yalnizca gorunurlukti.
DERS: agir adimlari ON PLANDA timeout ile kosturma; `setsid nohup` + log.

**(b) `freeze.py --out site`'i depo KOKUNE yazdim.** `deploy_pages.sh`
ilk satirlarinda `cd "$(dirname "$0")"` yapiyor, yani `webapp/site`'i
dagitiyor. Betik taze derlemeyi hic gormedi, eski kopyayi gh-pages ile
karsilastirdi ve "degisiklik yok, push yapilmadi" dedi -- BASARILI gibi
cikti verip hicbir sey yayinlamadan. Dogru cagri:
`cd webapp && python3 freeze.py --out site --base /pro-sim-blast`.
ACIK IS: betik `site/` yokken duruyor ama YANLIS YERDE bir `site/` varsa
sessizce onu dagitiyor; derlemenin yasini kontrol etmesi gerekir.

**(c) Reaksiyon raporu "UYUSMAYAN: 0" yaziyordu, cunku 21 sorgu HTTP 500
almisti.** Temizlik degildi: o tipler "kontrol edilemeyen" kovasina dustu ve
bir onceki kosuda BULUNMUS 6 uyusmazlik gorunmez oldu. Daha kotusu,
`validate_curation.py`'deki "her uyusmazlik kayitli" kontrolu bu halde BOS
geciyordu -- bakilmayan tipin uyusmazligi da yok sayilir, yani kapi tam
kapatmasi gereken durumda kapatmiyordu. `get()` 6 denemeye cikti (400/404/410
kalici sayilir), ozet KALICI boslugu GECICI hatadan ayri sayiyor, ve
"every reference with an accession was actually queried" kontrolu sonuctan
ONCE okunuyor. Yeniden kosuldu: SORGU BASARISIZ 0, 6 uyusmazlik geri geldi.

**Ayrica** `metabolic_role.csv` satir 74 "group 6" yaziyordu ve dogrulama
kapisi FAIL verdi; kaynak `types/6_601_PAO.yaml`'in `metabolic_role_basis`
alaniydi. chemistry.csv'deki duzyazi bu oturumda temizlenmisti ama ayni ifade
manifestonun BASKA bir alaninda duruyordu ve yeniden uretimde geri geldi.

**Dagitildi:** 3.218 degisen dosya, gh-pages tek commit. Yayinda dogrulandi:
/atlas/groups 200, birlesim indirmeleri 200, PDO 4,5 / PhtAa 3,4, regulator
uzlasma notu, tip sayfasi birlesim notu. Kapilar: `rebuild.py --check` 61/61
guncel, `validate_curation.py` 147 kontrol / 0 FAIL / 4 WARN, `check_site.py`
85 sayfa + 9.062 ic baglanti + 64.005 statik baglanti temiz.

**NCBI zinciri tamamlandi.** RefSeq release 237 protein uzayi: 871/871 parca,
76.043.118 protein, 15.407 aday, 0 hata. Karar: 5.083 RO alpha (122 fragment,
7.422 RO degil elendi), 4.550'si kapiyi geciyor. Birlesim 70.252 -> **73.299
tekil dizi**, kaynak kimligi 89.743 -> **94.826**; RefSeq'in getirdigi **3.047
YENI tekil dizi**. Yani NCBI protein uzayinin buyuk kismi UniProt uzerinden
zaten iceriydi -- ama 3.047 dizi hicbir yoldan girmiyordu ve simdi iceride.

## 114. PDO/PhtAa urunleri TERSTI; ve `types/` tek kaynak olmayi BIRAKMISTI (2026-11-01)

Kullanici bildirdi: "PDO ve PhtAa urunleri hala ters."

**Dogru.** Ve bu bir yazim karismasi degil, iki ayri metabolik yol.
Gram-negatif bakteriler ftalati **4,5**'ten, Gram-pozitifler (aktinobakteri)
**3,4**'ten acar; her iki yol da protokatekuata iner. Kayitta tam tersi
yaziyordu:

| tip | organizma | eski urun | duzeltilmis |
|---|---|---|---|
| `4_406_PDO` | *Comamonas testosteroni* KF1 (Gram-negatif, PDB 7FJL) | cis-3,4-diol | **cis-4,5-diol** |
| `3_311_PhtAa` | *M. vanbaalenii* PYR-1 (aktinobakteri) | cis-4,5-diol | **cis-3,4-diol** |

Dayanak: Stingley ve ark. 2004 (doi:10.1099/mic.0.27263-0) bu operonu tasiyan
6,7 kb parcanin ftalati 3,4-dihidroksiftalata cevirdigini yaziyor; KF1'in
PhtA/PhtB (= PDO/PDR) cifti 4,5-dihidroksilasyon yapiyor.

**IC TUTARLILIK ZATEN BOZULMUSTU ve kimse gormedi.** Sette ucuncu bir ftalat
tipi var: `4_405_ROCH34` (*Cupriavidus metallidurans*, Gram-negatif), urunu
cis-4,5-diol. Yani eski hal iki Gram-negatif ftalat dioksijenazina FARKLI
bolge atiyordu. Ayni substrat + ayni gram tipi + farkli bolge bir celiski
sinyalidir; hicbir kontrol substrat basina bolge tutarliligini aramiyor.
ACIK IS: `reaction` alanindaki "ring positions N and M" ifadesini substrat +
organizma gram tipine karsi sinayan bir kontrol yazilabilir.

**ASIL BULGU: chemistry.csv elle duzenlenebiliyordu.** `types_build.py`'nin
beyani net -- bir tip TEK dosyadir (`types/<kume>.yaml`) ve chemistry.csv
ONDAN URETILIR. Ama bunu zorlayan hicbir kontrol yoktu. Onceki oturumda
dus-notu temizligi dogrudan chemistry.csv'ye islenmisti; 10 hucre ayrisik
kalmisti ve bir sonraki `types_build.py --write` o duzeltmeleri SESSIZCE
geri alacakti. Bu, projenin tekrar eden hata sinifinin bir baska yuzu:
*kararı veren kural, verdigini sandigimiz kural degil.*

- 10 hucre (`1_117_PrnD`, `1_118_RedG`, `3_321_Nvd`, `3_322_DAF36`,
  `6_601_PAO`, `6_602_CAO`, `6_603_PTC52`, `7_701_AerC`) manifestolara tasindi
- chemistry.csv manifestolardan yeniden uretildi (`--write --allow-drift`)
- yeni kontrol: `check_types_are_single_source` -- gidis-donus farki varsa FAIL

## 115. Regulator sayisi: iki tablo da 5.022 OLCUYOR; eksik olan TOPLAM satiriydi (2026-11-01)

Kullanici bildirdi: "Tip bazindaki tablonun toplami 6.032, mimari tablosunda
5.022."

**Olculdu: ikisi de 5.022.** Mimari tablosu (`/atlas/regulation` bolum 2)
`ro_regulation`i tek sorguda sayiyor; tip tablosu (bolum 5)
`analysis_out/regulation_by_cluster.csv`den geliyor. CSV'nin
`divergent_regulator` sutunu 5.022, tablo sayimi 5.022, indirilebilir CSV
5.022, ve YAYINDAKI site de 5.022 -- dagitilmis sitenin 500 sayfasi
tarandi, her tablonun her sutun toplami hesaplandi, **6.032'ye esit tek bir
sutun toplami yok**. 6.032'nin nereden geldigi bulunamadi.

**Ama sikayet bos degil: sayfa okuyucuyu TOPLAMA yapmaya zorluyordu.**
68 satirli bir yuzde tablosunun ustteki tek sayiyla ayni seyi olctugu
yazmiyordu. Projenin tekrar eden tuzagi: *birim yazilmazsa okuyucu toplama
yapar ve tutmayan bir ucuncu sayi bulur.*

- tip tablosuna olculen `<tfoot>` toplam satiri eklendi (68 tip, 19.125 giris,
  %26 / 5.022)
- bes mimarinin ayni 19.125 girisi PARCALADIGINI acikca yazan not eklendi
- esitlik mekanik: `check_regulation_totals_agree` her mimari sutununun
  tip-basi toplamini tablo sayimiyla karsilastiriyor; ayrisirsa FAIL
  (6 yeni kontrol). Bu, CSV bayatladiginda sessizce iki farkli sayi
  gosterilmesini imkansiz kiliyor.

## 113. BIRLESIM KURULDU: 70.252 tekil dizi; NCBI protein uzayi ilk kez taraniyor (2026-11-01)

Kullanici: "NCBI'daki ordaki burdaki tum RO olabilecek her sey bu DB'den
toparlanabilsin". Bu, yayinladigimiz 19.125'in ne oldugunu yeniden
tanimliyor: o sayi bir KAPSAM degil, bir KESISIM -- hem genomik baglami
olan, hem kapiyi gecen girisler.

**`ro_union` (build_union.py).** Birim ACCESSION DEGIL, DIZI. Gerekce
olculmustu: UniProt accession'larinin %65,1'i artik silinmis; accession'i
birim saymak UniProt'un temizlik kararlarini bizim biyolojik sayimimiz
gibi gostermek olur. Accession artik `ro_union_source`'ta, kaynagi ve
durumuyla birlikte duruyor.

```
NCBI genom (anotasyonlu, onayli) : 19.125 giris ->  15.960 tekil dizi
UniProt anlik goruntusu          : 73.699 giris ->  54.292 YENI tekil dizi
-----------------------------------------------------------------------
BIRLESIM                                            70.252 tekil dizi
  genomik baglami olan                              15.960
  genomik baglami OLMAYAN                           54.292
  kapiyi gecen                                      69.180
  kapiyi GECMEYEN (tabloda KALIR, aday olarak)       1.072
ro_union_source                                     92.824 kaynak kimligi
  uniprot/deleted 47.988 | uniprot/live 25.706
  ncbi_genome/confirmed 16.044 | merged 3 | demerged 2
```

Yani ayni dizi ortalama 1,3 kaynak kimligi tasiyor; 13.984 dizi iki,
2.070 dizi uc veya daha fazla kimlikle geliyor. Kapiyi gecmeyen 1.072
satir SILINMEDI -- "RO olabilecek her sey" ancak boyle toparlanir, hangi
tiere bakacagina kullanici karar verir.

**NCBI protein uzayi (scan_ncbi_proteins.py).** Madde 95 bunu zaten
yazmisti: protein evreni UniProt'tan geldi, NCBI'DAN DEGIL. Yani
kullanicinin sordugu seyin en buyuk parcasi bugune kadar HIC aranmamisti.
RefSeq release 237, 871 protein parcasi, ~143 milyon protein, ~36 GB
sikistirilmis. Akitarak calisiyor (indir -> tara -> eslesenlerin dizisini
al -> sil), hicbiri diskte birikmiyor, durum dosyasiyla kesildigi yerden
devam ediyor.

ESIK KEYFI DEGIL, OLCULDU. Onaylanmis 19.125 girisin TAMAMI ROmotif71'e
E <= 5,6e-15 ile eslesiyor; en kotu bit skoru 55,5, medyan 189,4. On
suzgec bu yuzden E <= 1e-3'te: en kotu gercek girisin on kat mertebe
altinda, yani hicbir gercek RO on suzgecte dusmez. Suzgec KARAR VERMEZ --
karari kapi (Rieske 4/4 + triad >=2/3) verir; burada toplanan sey ADAY.

Deneme parcasinda 164.200 proteinde E<=1e-3'te 36 aday, skor araligi
22,8-198,9; bunlarin 18'i onaylilarin en kotu skorunun uzerinde. Bu
oranla tum surumde ~31.000 aday bekleniyor.

DIKKAT, bir sayim tuzagi: `hmmsearch --tblout /dev/stdout` kullanirsan
insan-okunur cikti ayni akisa yazilir ve `grep -vc "^#"` 31 yerine 244
sayar. tblout DAIMA bir dosyaya yazilmali.

## 111. Statik arama sayfasinda facet baglantilarinin HICBIRI calismiyordu (2026-11-01)

Kullanici bildirdi: "search kisminda sol sekmelere tikliyorum bir sey olmuyor".
DOGRU, ve sebep tek satirlik bir kacis hatasiydi.

freeze.py facet baglantisini soyle uretiyordu:

    onclick="setFilter('${key}', ${JSON.stringify(v)});return false"

`JSON.stringify` CIFT TIRNAKLI bir deger dondurur ve bu deger CIFT TIRNAKLI
bir HTML niteliginin icine giriyordu:

    onclick="setFilter('tier', "characterized");return false"

Nitelik ilk IC tirnakta kapanir. Tarayici `onclick` olarak yalnizca
`setFilter('tier', ` goruyor -- yarim bir ifade, sozdizimi hatasi -- ve geri
kalani baska nitelik saniyor. Tiklama hicbir sey yapmiyordu. Node ile
sinandi: dort ornek degerin DORDUNDE de onclick bozuk.

NEDEN BU KADAR UZUN GOZDEN KACTI. Ayni blokta "Genomic context" baglantilari
VARDI ve CALISIYORDU, cunku onlarin degeri (`plasmid`, `partner`,
`regulator`) sabit bir tanimlayici ve tirnak icermiyor. Yani sayfanin yarisi
calisiyordu; bozuk olan yarisi, degeri veriden GELEN yariydi. Bir kontrol
"facet baglantilari var mi" diye baksaydi gecerdi -- baglantilar vardi.

YAPILAN. Satir ici onclick birakildi. Deger artik `data-facet` /
`data-value` niteliklerinde tasiniyor, `esc()` ile kaciriliyor, ve tiklama
`#facets` uzerindeki TEK bir dinleyiciyle yakalaniyor (facet listesi her
render'da yeniden yazildigi icin delegasyon dogru bicim). Degerin icinde
tirnak, apostrof veya `<` olmasi artik hicbir seyi bozmaz.

DOGRULAMA, SADECE KOD OKUYARAK DEGIL. Firefox headless ile iframe'li bir
kosum sayfasi yazildi; sayfa gercekten yuklendi, ilk facet baglantisina
GERCEKTEN tiklandi:

    durum: 19,125 entries indexed.
    facet baglantisi sayisi: 33
    tiklanacak: "distant"  (facet=tier, value="distant")
    tiklamadan once: "19,125"  ->  sonra: "10,385"
    secili isaretlenen: "distant"
    SONUC: GECTI

10.385, veritabanindaki `distant` sayisiyla birebir ayni -- yani yalnizca
"bir sey oldu" degil, DOGRU sey oldu.

## 110. Gorsel denetim: olculmeyen sey kirpiliyordu (2026-11-01)

Kullanici bir ajanin siteyi gorsel olarak kontrol etmesini istedi. Firefox
headless ile gercek ekran goruntuleri alindi (1440 masaustu, 390 mobil,
karanlik mod) ve uretilmis butun SVG'lerde metin genisligi <-> viewBox
karsilastirmasi yapildi. Bulunan 14 sorunun yedisi bu turda duzeltildi.

(a) OPERON ALT YAZISI 2.855 YERDE KIRPIKTI. Sekil genisligi yalnizca
    KUTULARDAN hesaplaniyordu; alt yazi onlardan uzundu ve hesaba HIC
    girmiyordu. En kotu tasma 321 px. Operon kisaldikca kayip buyuyordu,
    yani tek alfali tipler -- aciklamayi en cok gereken sekiller -- en cok
    zarar goruyordu. `text_width()` eklendi (sans-serif ~0,54em, USTTEN
    tahmin: fazla genislik zararsiz, eksik genislik keser) ve genislik artik
    `max(kutular, alt yazi)`. Olcum sonrasi tasma 0.

(b) "not a regulator" ETIKETI 2.995 YERDE SOL KENARDAN KESIKTI. Kutu 46 px,
    etiket ~65 px, ortalanmis capa x=-2'den basliyordu. BU TURDA EKLENEN bir
    seydi, yani hic test edilmemisti. Capa artik pad'in solunda baslamayacak
    sekilde kaydiriliyor.

(c) AGAC ETIKETLERI KELIME ORTASINDAN KESILIYORDU. `text[:60]` uretilen
    SVG'nin ICINDE veri kaybi yapiyordu: "reference (benzalkonium chloride
    and other quat" -- parantez acik, kelime yarim. 15 tip sayfasinda 22
    tane. `clip_label()` kelime sinirinda keser ve `…)` ile kesildigini
    soyler; tam metin zaten <title>'da.

(d) BOS KARTLARDA KONTRAST WCAG AA'YI GECMIYORDU. `.structcard--empty`
    uzerindeki `opacity:.72`, zaten `--dim` olan ~12 px metni acik modda
    2,86 ve karanlik modda 3,59 kontrasta dusuruyordu (esik 4,5). Tokenlar
    tek baslarina sorunsuzdu; carpan sorunluydu. Opaklik kaldirildi; "inert"
    oldugunu kesikli cerceve zaten soyluyor.

(e) `.rxnote--fate` 34ch'e SIKISIYORDU. 34ch bir SEMA KOLONUNUN altindaki
    kucuk not icin dogru; bu ise seklin tamami genisligindeki bir paragraf
    ve ayni siniri devraliyordu. 1300 px panelde metin ~200 px'lik seride
    5-8 satira sariyordu. 78ch yapildi.

(f) PLOTLY EKSEN BASLIGI KIRPIKTI. Eksen basligi tik etiketlerinin ALTINA
    yazilir, yani alt marj iki satir tasir; b:50 ile baslik kabin altindan
    tasiyordu. b:50->72, b:45->64, ve dondurulmus etiketler icin r:10->24.

(g) URUN KUTUSU OKUNDAN KOPUYORDU. DOM sirasi substrat - ok - SEMA - urun
    idi: genel reaksiyon semasi okla urunun ARASINA giriyordu ve substrat
    SVG'si genis oldugunda urun alt satira sariyordu (qxyA'da sagda ~280 px
    bosluk). Substrat+ok+urun tek bir `.rxn__eq` kapsayicisinda toplandi,
    sema sona alindi ve dar ekranda kendi satirina iniyor.

ANA SAYFA. Alti ayri istatistik kutusu TEK SERIDE birlestirildi (ortak
cerceve, ince ayirici, tabular rakam, versal etiket): alti kutu yan yana
gelince sayfa bir gosterge paneli gibi okunuyordu, bir veri kaynaginin
kunyesi gibi degil. Giris paragrafi 72ch'ten 88ch'e cikti -- govde olcusunu
devraliyordu ve 1180 px'lik kolonun sol yarisinda asili duruyordu. Masaustu,
390 px mobil ve karanlik modda ekran goruntusuyle dogrulandi.

KALAN (yapilmadi): tablo tasmalari (/clusters ve tip sayfasindaki carpisma
tablosu 1440 px'te sagdan kesik), arama sayfasinin bos acilisi, mobilde
Plotly kirpilmasi, RedG aciklamasindaki ic gelistirme notlari.

## 109. Protein-genomik bosluk TUM tiplerde; genomik taraf GERCEKTEN buyutulebilir,
ama protein tablosunun %63'u SILINMIS UniProt kaydi (2026-11-01)

Kullanici: "bundan boyleyse diger enzimler icin de boyle olabilir, bence genomik
anlamda da daha fazla var gibi". Iki iddia da olculdu; ikisi de DOGRU cikti, ve
olcum ucuncu bir seyi ortaya cikardi.

(a) BOSLUK QXYA'YA OZGU DEGIL, KURAL. Her iki tabloda da uyesi olan 68 tipin
    **67'sinde** protein sayisi genomik sayidan buyuk. 63 tipte >=2 kat, 20
    tipte >=5 kat. En buyukler: BphA1(2_201) 12,0x, LigX 10,5x, ArhA1 9,8x,
    CmoB 9,4x, qxyA 8,2x, PTC52 6,6x, Dpe 6,5x, ROCH34 6,3x. Toplam: genomik
    19.125, protein (kapiyi gecen) 72.586, link'i olmayan 55.122.

(b) DORT TIP YALNIZCA PROTEIN TARAFINDA VAR: BpdC1 13, XylC1 3, NahAc 1,
    NinA 1. Bunlar duran "sifir uyeli 8 referans" WARN'inin icinden geliyor --
    yani o referanslar bos DEGIL, uyeleri genomik tabloya hic girmemis.

(c) "GENOMIK ANLAMDA DA DAHA FAZLA VAR" DOGRULANDI. 55.122 linksiz girisden
    300'luk rastgele ornek UniProt REST'e soruldu. Halen aktif olan kayitlarin
    **%99,1'inde** EMBL nukleotit capraz-referansi var, yani genomik baglam
    KURULABILIR. Olcek: ornekte 110/300 = %36,7 -> **~20.200 giris**
    (%95 GA 17.200-23.200). Bu, genomik tabloyu 19.125'ten ~39.300'e cikarir,
    yani iki katina.

(d) AMA PROTEIN TABLOSU BAYAT BIR INDIRMEYE DAYANIYOR. Ayni ornekte
    **189/300 = %63,0** kaydin entryType'i `Inactive`, ve sebebi sorulan 38
    kaydin 38'i `DELETED` (birlestirme DEGIL -- yani yeni bir accession altinda
    yasamiyorlar). Olcek: linksiz 55.122 girisin **~34.700'u** (%95 GA
    31.700-37.700) bugun UniProt'ta artik yok. Madde 95'teki indirme tarihi
    o zaman dogruydu; aradan gecen surede UniProt gereksiz proteomlari
    buyuk olcude sildi. Site bu girisleri hala yayinda tutuyor.

(e) BEDAVA KAZANC YOK. Kurtarilabilir kayitlarin replikonlari bizde zaten
    duruyor mu diye bakildi: ornekte 3/110. Olcek ~550 giris. Yani is gercek
    bir indirme isi.

    ILK OLCUM YANLISTI VE KENDI KONTROLUM YAKALADI: ilk karsilastirma %0,0
    dedi. Sebep, kimligi belirleyen kuralin karsilastirdigim kural OLMAMASIYDI
    -- bizim `replicon.nucleotide_id` tamami RefSeq (NZ_ 14.449, XM_ 1.859,
    NC_ 727), UniProt ise INSDC/ENA birincil accession veriyor
    (JACHOA010000004, CP026258, FOYL01000002). `NZ_JACHOA010000004` ile
    `JACHOA010000004` AYNI kayittir. On-ek soyulunca ortusme %2,7'ye cikti.
    Projenin tekrar eden hata sinifinin bir ornegi daha, bu kez yayina
    girmeden yakalandi.

SIRADAKI IS (yapilmadi, karar kullanicinin):
  1. `ro_protein`'i bugunun UniProt'una karsi tazele ve silinmis accession'lari
     isaretle. Bunlar SESSIZCE silinmemeli: "silinmis kayit" ile "veri yok"
     ayri seylerdir, ve 34.700 satir sessizce yok olursa yayimlanmis her
     protein-duzeyi sayi aciklamasiz degisir.
  2. Aktif + EMBL capraz-referansi olan ~20.200 girisin nukleotit kayitlarini
     indir, anotasyonlu olanlardan operon/komsuluk kur, `ro`'ya ekle.
     Once (b)'deki 18 proteinle basla: 9'unun genomik bagi kurulabilir ve
     BpdC1/XylC1/NinA boylece gercek genomik uye kazanir.
  3. Oncelik sirasi (d) yuzunden 1 -> 2 olmali: bayat kumeden indirme listesi
     cikarmak, silinmis kayitlar icin bos sorgu atmak demektir.

## 108. Geri cekilmis kolon iddiasi 1 + 68 sayfada YAYINDA KALMISTI (2026-11-01)

Kullanici bildirdi: "Atlas sayfasi hala 'kolon 431 aspartattir, 15.996 Asp'ye
karsi 122 Glu, yani Fe uclusu Asp kullanir' diyor. Oysa ana sayfanin 1.2
bolumu, kolon 431'in ligand olmadigini yapilarla gosteriyor."

DOGRU. Kolon 431'in demir ligandini okumadigi 2026-10-31'de olculmustu: 27
kuratorlu referansin kristal yapisina karsi sinandi ve cagirdigi kalinti
**26/26'sinda demirden 10,9-14,7 A** uzakta cikti. Kolon testi geciyordu cunku
o civarda gercekten bir D/E var -- test "kolonda asidik kalinti var mi" diye
soruyor, "metali koordine ediyor mu" diye SORMUYOR.

KOK SEBEP: GERI CEKME, KESFEDILDIGI YERE UYGULANDI; YAYIMLANDIGI HER YERE
UYGULANMADI. index.html'in 1.2 bolumu uzun bir yorumla duzeltilmisti, ama ayni
iddia iki yerde daha yayindaydi:
  - `atlas_home.html`: "...is aspartate and effectively nothing else, 15.996 Asp
    against 122 Glu, so the facial triad uses Asp specifically" -- kullanicinin
    gordugu cumle. Sayiyi yayimlamakla kalmiyor, ondan LIGAND SONUCU cikariyor.
  - `cluster.html`: "Catalytic carboxylate (column 431) -- Asp in N of M
    members ... This site is aspartate across the whole family." Render edilince
    **68 tip sayfasinda** goruluyordu.
Hicbir sayi kontrolu yakalamadi, cunku celiski DUZYAZIDA yasiyordu.

YAPILAN:
  - `atlas_home.html`: katalitik cumle, yapinin soyledigi seyle degistirildi --
    referans setinde 26/26 aspartat, AMA tek bir kolon hepsini tutmuyor (431'de
    18, 432'de 2, kalan 6'si 407-411 ve 449-451'e dagilmis), ve kolonun
    Asp/Glu orani ligandin olcumu DEGIL, o yuzden yayimlanmiyor. Koprulesen
    kolonun (286) ayni yapisal testi GECTIGI ayrica yazildi, cunku eski son
    cumle ("bunlar yapidan degil kolondan okunan kalintilar") artik koprulesen
    kolon icin eksik bir beyandi.
  - `cluster.html`: uye-basi sayim kaldirildi, yerine "sabit bir kolondan
    okunamaz" + yapisal gerekce + 1.2'ye baglanti kondu.
  - `validate_curation.py`: `check_retracted_column_claim()` eklendi. Sablonlari
    tarar ve (i) katalitik kolonun asp/glu/dominant degerlerini basan her Jinja
    ifadesini -- takma adla (`{% set ci = ... %}`) basilani da -- ve (ii) uc
    geri cekilmis duzyazi kalibini FAIL eder. Jinja yorumlari haric tutulur.
    Eski metne karsi sinandi: 5 bulgu yakaliyor. Repo bugun temiz.

## 107. Iki SIRA listesi SUZGEC gibi davraniyordu; agacin yontem beyani yanlisti (2026-11-01)

Kullanici dort sey sordu: "2 - What the reaction actually is" bolumunde tum
reaksiyonlar var mi, agac bootstrapli ve guncel mi, qxyA'ya benzeyen tip neden
yok, ve "DB'de protein acisindan eksiklik var gibi".

### (a) Reaksiyon bolumu tiplerin %11'ini ve girislerin %7,5'ini saymiyordu

Iki ELLE tutulan liste vardi ve ikisi de SIRA listesiydi ama SUZGEC gibi
davraniyordu -- listede olmayan deger sessizce dusuyordu:

- `atlas.FAMILY_ORDER`: `chemistry_groups()` yalnizca bu listeyi geziyordu.
  `tetrapyrroles` ailesi listede YOKTU, yani 6_601_PAO, 6_602_CAO, 6_603_PTC52
  (617 giris) ana sayfanin kimya bolumunde HIC gorunmuyordu -- ne karti ne
  tablo satiri. Uc tip 14 Ekim 2026'da eklenmisti.
- `reaction_summary()` icindeki sabit anahtar listesi: bes reaksiyon sinifi
  (`desaturation`, `N_oxygenation`, `C_C_bond_formation`,
  `oxygenolytic_ring_opening`, `C_C_cleavage`) hic yayinlanmiyordu.

Ustelik `reaction_summary(groups)` gruplardan besleniyor, yani ilk suzgec
ikincisini de besliyordu.

Olculdu: yayinlanan 7 satir / 68 tip / 17.685 giris; gerceği 12 satir / 76 tip
/ 19.125 giris. Bolumun kendi cumlesi "cis-dihydroxylation is only the most
common case, not the whole story" diyordu -- kendisi hikayenin bir kismini
atlarken.

Yapildi: iki liste de SIRA olarak kullaniliyor, listede olmayan degerler sona
ekleniyor. `tetrapyrroles` FAMILY_ORDER/LABEL/NOTE'a girdi. Bes sinif icin
etiket, aciklama, renk, genel sema ve sema notu yazildi -- semalar metin
tabanli tutuldu, cunku regiyokimya varsaymak yanlis gostermekten kotu.

Bir daha olmamasi icin `check_chemistry_published`: her tipin karta ulastigi,
tablonun tip ve giris sayisinin kaynagin tamamina ESIT oldugu, her reaksiyon
sinifinin etiketi/notu/semasi ve her ailenin etiketi oldugu olculuyor.
121 kontrol: 117 PASS, 0 FAIL, 4 WARN.

### (b) Agac guncel, destek BOOTSTRAP DEGIL, ve model adi yanlis yaziliydi

Agac guncel (31 Ekim 21:48, rebuild --check 0 bayat). Iki kusur:

1. **Destek hic yayinlanmiyordu.** FastTree'nin varsayilan destegi SH-benzeri
   YEREL destektir: her ic dal kendi en yakin komsu takaslarina karsi test
   edilir, hizalama kolonlari YENIDEN ORNEKLENMEZ. `parse_newick` bu degerleri
   ic dugum ADI olarak zaten okuyordu ama HIC KULLANILMIYORDU; altyazi da
   "branch support values are not shown on tip labels" diyerek konuyu
   kapatiyordu. Olculdu: 2.046 ic dal, medyan 0,923, %41,9'u >=0,95, ve
   **491'i (%24) 0,70'in altinda** -- yani dortte bir dalda dallanma sirasi
   COZULMEMIS ve bu okuyucuya hic soylenmiyordu.
2. **Altyazi "FastTree under the JTT model" diyordu, fiilen kosan model LG.**
   build_phylogeny'de bir YEDEKLEME ZINCIRI var (LG -> JTT -> ML'siz) ve
   hangisinin kostugu calisma zamaninda belli oluyor; altyazi ise sabitti.

Yapildi: `build_phylogeny.py` artik `analysis_out/tree_meta.json` yaziyor
(kullanilan model, yedekleme zinciri, hizalayici, destek turu ve dagilimi) ve
bu dosya rebuild.py'de BEYAN EDILEN cikti oldu, yani tazelik takibinde.
Altyazi modeli oradan okuyor. `newick_support_summary()` dagilimi veriyor;
`render_tree_svg` 0,70 altindaki dallara ici bos turuncu daire koyuyor ve
degeri `<title>`'a yaziyor -- 2.046 dala sayi yazmak okunmaz olurdu ve
okuyucunun kararini degistiren bilgi "bu dal cozulmemis" olanidir. Tip
sayfasindaki agaclarda da ayni isaret ve kendi sayilari var.

### (c) qxyA'ya benzer tip VAR; "yok" gorunmesinin iki ayri sebebi

Bir tuzak: `cluster_identity_matrix.csv` qxyA/BedC1 icin 91,4 veriyor. Bu
benzerlik DEGIL -- iki tipe atanmis herhangi iki temsilci arasindaki en yuksek
kimlik, yani atama sinirinin sizdiginin olcusu. Proje bunu zaten belgelemis ve
"en yakin tipler" tablosunu `cluster_reference_identity.csv`'ye baglamis.

Gercek referans-referans kimlikleri: Stc2 27,4 / BmoA 23,7 / GbcA 21,7 /
CntA 21,1; medyan 16,8, maks 27,4, %30 uzeri YOK. Yani en yakin komsular tam
kimyasal olarak beklenen yerde (hepsi grup 5, hepsi C-N+ kesen) ama hepsi
`family_member` esiginin (%40) altinda. qxyA gercekten yalitik: eksiklik degil.

### (d) Protein tarafindaki eksiklik GERCEK ve buyuk

| tip | genomik baglamli | protein-only | kat |
|---|---|---|---|
| 5_505_qxyA | 164 | 1.370 | 8,4x |
| 1_109_LigX | 87 | 929 | 10,7x |
| 4_405_ROCH34 | 364 | 2.295 | 6,3x |
| TOPLAM | 19.125 | 73.699 | 3,9x |

Ana tablo yalnizca genomik baglami olan proteinleri tutuyor (operon
turetilebilsin diye). qxyA'nin 1.370 protein-only akrabasinin 1.206'si ana
tabloda yok. Madde 91 ve 97'nin olcumuyle tutarli; kullanicinin hissettigi
bosluk tam olarak bu.

### Ayrica

- `SHORT` etiketleri ve `other`/`none` kutulari: madde 106'da duzeltildi.
- Bolum 3 grup basliklarinda "1 types" yaziyordu; tekil/cogul duzeltildi.
- `tree_meta.json` belgelendi, yani belgesiz analiz dosyasi listesine
  eklenmedi.


## 106. Operon sekli, DUZENLEYICI BULUNMADAN "promoter bolgesi" ciziyordu (2026-10-31)

Kullanici bildirimi: "qxyA sayfasindaki gene order bence yanlis, bunu
hazirlamada bi bug var gibi" ve ardindan "qxyA operonu icin bence
TetR --- promoter ---- a - scp2 gibi bisi olmali ama yok iste".

### Gen SIRASI dogruydu; yanlis olan onun yanina cizilen sey

Once sira olcumle aklandi. `build_operons.walk()` supheliydi, cunku eksi
iplikte bir ters cevirme yapiyor; `gene_offset`'in GENOMIK mi iplik-goreli mi
oldugu olculdu (19.125 girisin ±3 komsulari): her iki iplikte de negatif offset
genomik olarak SOLDA (26.621/26.624 ve 27.104/27.119). Yani yurume ve boslugun
aritmetigi tutarli. Ardindan asil test: `[alpha]` tek basina cikan 51 yakin
qxyA uyesinde, ayni iplikte 150 bp icinde komsu var mi -- **0**. Sira dogru.

### Asil hata: sekil, olcumun kurmadigi bir iddiayi ciziyordu

`operon_regulator_svg` duzenleyici okunu yalnizca `*_regulator` mimarilerinde
ciziyordu, AMA kesik promoter kutusu `bp is not None` oldugu an ciziliyordu ve
SVG'nin alt yazisi ona "dashed box is the intergenic promoter region" diyordu.
Tip sayfasinin alt yazisi da enum'u oldugu gibi basiyor ("divergent other") ve
sonunda kosulsuz olarak "where the regulator is divergent, it and the operon
share that single region" diyordu.

Sonuc: duzenleyici BULUNAMAMIS tiplerde, yalniz bir alfanin yanina isimsiz bir
"promoter bolgesi" kutusu konuyor ve alt yazi duzenleyiciden bahsediyordu.

Olculdu (2026-10-31): mimari verisi olan 68 tipin **40'inda** baskin mimari
`*_other`, yani cizilen promoter kutusunun yaninda tanimli hicbir duzenleyici
yok. Baskin mimari dagilimi: divergent_other 32, divergent_regulator 17,
codirectional_other 8, unknown 7, codirectional_regulator 4.

qxyA en yanilticisiydi: 53 yakin uyenin **52'si** `divergent_other`; sayfa
"[alpha] + 138 bp promoter" gosteriyordu. Literaturun qxyAR operonu
(PMID:27792326) gercek, ama TetR yalnizca KARAKTERIZE girisin kendisinde var
(1/53, 78 bp).

### Atlas tablosunda ayni kusurun daha serti

`atlas_operons.html` mimariyi SABIT `'divergent_regulator'` yaziyordu ve
`divergent_regulator_rate >= 0.2` oldugu an ok ciziyordu; okun basligi da
"typical upstream regulator of this type" idi. Yani uyelerinin %80'inde
duzenleyici olmayan bir tipe "kanonik duzenleyici" atfediliyordu. Olculdu:
eski kuralla 30 tipe ok cizilirdi, bunlarin **15'inde** baskin mimari
`*_other`. Artik mimari sayilardan secilir (`atlas.dominant_arch`), ok 17
tipte cizilir.

### Yapilanlar

- `atlas.operon_regulator_svg`: yukari komsu BILINIYORSA her zaman cizilir --
  duzenleyiciyse mor + aile adi, degilse notr gri + "not a regulator"
  etiketi. Kutu ancak duzenleyici varken "promoter" diye adlandirilir; alt
  yazi uc duruma gore ayrisir.
- `atlas.ARCH_TEXT`: bes mimarinin insan dilindeki karsiligi, TEK kaynak.
- `atlas.dominant_arch`: `regulation_by_cluster.csv` satirindan baskin mimari.
- `cluster.html`: enum sizmasi kalkti; duzenleyici yokken alt yazi kutunun
  "olculmus mesafe, kanitlanmis promoter DEGIL" oldugunu soyluyor.
- `atlas_operons.html`: sabit mimari kalkti, aciklama paydasini soyluyor.
- `SHORT["other"/"none"]`: `"·"` -> `"gene"`. Nokta ekranda HICBIR SEY gibi
  okunuyordu; anlami "gen var, annotasyonu kategoriye uymadi".
- `operon_regulator_svg(..., titles=)`: tek-giris sekillerinde oklara GERCEK
  urun adlari baslik olarak baglanir. qxyA referans sekli artik
  **TetR -[78 bp]- alpha - SCP2 - HNH** okunuyor; kullanicinin bekledigi sekil
  buydu ve ZATEN hesaplaniyordu, yalnizca canli derlemede yoktu.
  `build_entry_data.py` bunu KULLANMAZ: dedup anahtari sekle dayanir, urun adi
  girise ozgu oldugu icin her girise ayri SVG uretirdi.

### Yol ustunde cikan gerileme

`atlas.tip_info_from_db`, yaprak medyanini `f"{ident:.2f}"` ile biciml[e]ndiriyordu;
madde 104'te medyan tek uyeli yaprakta `None` yapilinca filogeni atlas sayfasi
`TypeError` veriyordu. Artik "not defined (single sequence)" yaziyor.


## 105. Karboksilat penceresi, ve indirmelerin derleme kimligi (2026-10-31)

Iki dis bildirim, ikisi de yerinde cikti -- ama ikisi de tahmin edilenden
BASKA bir yerde.

### (a) Pencere, yayimlanan Asp:Glu oranini kirletiyordu  [DUZELTILDI]

Bildirim: "kolon 471'de 6.008 Asp'ye karsi 7.583 Glu; literaturde bu ligand
neredeyse her zaman Asp; pencere yandaki bir Glu'yu yakaliyor olabilir."

Iki tashih: kolon 471 EMEKLI (30 Ekim'de 431'e duzeltildi; atlas.py uc dosyayi
"measured at the wrong carboxylate column (471)" diye SUPERSEDED isaretliyor)
ve 6.008/7.583 guncel hicbir sayiya denk gelmiyor. Ustelik %71, KATALITIK
ligandin degil KOPRU karboksilatinin orani:

    katalitik Fe ligandi (431)   Asp 16.299 : Glu   439   Asp %97,4   ref 68D/0E
    kopru (286)                  Asp 13.964 : Glu 5.156   Asp %73,0   ref 66D/3E

Koprude Glu GERCEK: grup 5'in tamami Glu baskin, CntA'nin kopru glutamati
E205 ve E205A/E205D TMA uretimini tamamen durduruyor (Zhu 2014, PNAS 111:4268;
citations.yaml: zhu2014_CntA). Koprudeki offset -1 kaymasi YALNIZCA grup 5'te
var (1.940/1.940); diger bes grup tamamen offset 0'da. Yani artefakt degil,
gruba ozgu gercek bir indel.

AMA KATALITIK LIGANDDA BILDIRIM DOGRU. Olculdu:

    tam kolon (offset 0)   Asp 15.902 : Glu 122   = 130 : 1
    kaymis pencere         Asp   397 : Glu 317   = 1,25 : 1

Kaymis okumada oranin basa bas gelmesi, pencerenin ligand yerine yakindaki
herhangi bir karboksilati topladiginin imzasidir. En buyuk tek kaynak
3_313_KshA15: 194 kaymis Glu. Sebebi bulundu: o tipin kuratorlu referansi kolon
431'de BOSLUKLU ve +-4 penceresinde hic karboksilat yok ('--------T'), en
yakini bes pozisyon uzakta bir Glu. KshA15'in referansi PDB 2ZYL ile %100 ayni
ve yapidaki mononukleer demir ligandi Asp'dir.

Kolon 431'de karboksilat TASIMAYAN sekiz referans: LigX, KshA15, Nvd, DAF-36,
AerC (bosluk) + PobA 'K', ROCH34 'Q', CAO 'L'. Son ucu ayrica kayda deger:
referansta ligand kolonunda karboksilat olmayan bir kalinti oturuyor.

ILK YAMA GERI ALINDI, cunku fazla genisti. "Referansi bosluklu tipte pencereyi
tamamen reddet" diyordu; ama capali okumanin o 240 Glu'yu Asp cagiracagi
OLCULMUS degildi -- "492 D / 0 E" rakami pencerenin BASARISIZ oldugu baska 492
girisi kapsiyor. Olcerek daraltilan kural:

    Kaymis bir pencere cagrisi ancak tipin kuratorlu referansi model kolonunda
    AYNI karboksilati tasiyorsa kabul edilir. Aksi halde reddedilir ve karari
    tip-ici capali okuma verir. Tam kolondaki cagrilar her zaman kabul edilir
    (76 referansin 68'i orada Asp: kolonun kendisi iyi desteklenmis bir capa).

Olculen etki: 16.024 tam kolon kabul; 558 kaymis + referansta karboksilat yok
REDDET; 62 kaymis + referansla farkli kalinti REDDET; 94 kaymis + referansla
ayni kalinti KABUL (gercek register kaymasi). Toplam 620 cagri (%3,7) capali
okumaya dusuyor. Yayimlanan oran 37:1 -> ~125:1.

Uygulama: carboxylate_audit.py'de reference_residue_at_carbox_column +
offset_call_is_supported. Dikkat edilen nokta: capali asama `negatives`
listesinden beslenir ve o liste stage_baseline'dan ONCE kurulur, o yuzden kural
IKI yerde de uygulanmak zorunda -- yalnizca stage_baseline'da reddetmek girisi
capali asamaya DUSURMEZDI.

### (b) Indirmeler dogruydu; eksik olan derleme kimligiydi  [DUZELTILDI]

Bildirim: "FASTA'da 19.125 kayit var; 12'si CSV'de yok, CSV'deki 1 kayit da
FASTA'da yok."

Olculdu -- yerel dondurulmus dosyalar VE canlidan indirilenler:

    FASTA 19.125   CSV 19.125   yalnizFASTA 0   yalnizCSV 0   yinelenen 0
    MD5'ler yerel kopyayla ayni; CSV'de gomulu satir sonu yok (19.126
    fiziksel satir = 19.126 kayit); protein_id uzerinden de 0/0.

Yayindaki 97 indirme dosyasinin hepsi FASTA'ya karsi tarandi: hicbirinde 12/1
imzasi yok. Ama ARITMETIK bir sey soyluyor:

    19.125 - 12 + 1 = 19.114

19.114, yeniden kurulumdan ONCEKI onayli giris sayisi. Ve dusen giris de
bulundu: XM_009028827.1:0-1056:1, su an is_confirmed=0, iki guncel dosyanin
ikisinde de yok. Yani bildiren kisinin CSV'si eski derlemeden, FASTA'si
yenisinden. Dosyalar tutarliydi.

BUNA RAGMEN BU BIR KUSURDUR: indirilen bir dosyanin hangi derlemeden geldigini
anlamanin yolu yoktu, ve okuyucunun bunu bulmak icin aritmetik yapmasi gerekti.
Bir veri tabani kendi dosyalarinin kimligini okuyucuya yukleyemez. Eklendi:
webapp/freeze.py write_download_manifest -> download/MANIFEST.txt; her dosya
icin bayt, sha256 ve KAYIT SAYISI (FASTA'da '>' satiri, CSV'de baslik dusulmus
satir), ustte derleme damgasi (onaylanmis giris sayisi + veri tabaninin
sha256'si). Ana sayfadan baglandi; sunucu modunda baglanti verilmez, cunku
orada dosyalar her istekte canli DB'den uretilir.

### (c) Sekiz bos tipin CSV'si BASLIKSIZ iniyordu  [DUZELTILDI]

(b)'yi arastirirken bulundu. csv_response'taki `if rows:` korumasi, sonuc bos
oldugunda SIFIR BAYT dosya uretiyordu: CarAa, IpBAa, BpdC1, XylC1, NahAc,
DntAc, NinA, PbaA1. Indiren kisi bunun "uye yok" mu "indirme bozuldu" mu
oldugunu ayirt edemiyor. Bos bir tablo bir SONUCTUR; bos bir dosya ARIZA gibi
okunur. Kolon adlari artik cursor.description'dan aliniyor (query_csv), yani
satir olmasa da baslik yazilir.

## 104. Uyelik temeli, tekil yaprak medyani, yutan tip (2026-10-31)

Madde 103'ten sonra kalan ucu acik kararlar kapatildi, ve kapatirken ayni hata
biciminin uc yeni ornegi bulundu. Hepsi TEK yeniden kurulumda uygulandi.

### (a) Tek uyeli yaprak "Median identity 1.000" yayinliyordu  [DUZELTILDI]

`recursive_homogenize.median_identity` cift yoksa 1.0 donuyordu ve bu deger
tabloya yaziliyordu. Olculdu: 2.789 yapragin **1.179'u** (%42,3) tek uyeli, yani
sitenin bu sayfalarinin hepsi hic var olmayan bir ciftin medyanini mukemmel
olcum gibi basiyordu. Ustelik 1.0 en yuksek olasi deger oldugu icin o yapraklar
"en homojen" siralamasinin TEPESINE cikiyordu -- yanlis yonde en kotu hata.

Artik None/NULL doner, sayfa sayinin yerine nedenini yazar
("not defined -- a single sequence has no pair to measure over"), ve
`leaf_stop_reason` ayri bir `single_member` nedeni dondurur: `is_homogeneous=0`
olmasi hedefin tutturulamadigi degil, hedefin SORULAMADIGI anlamina gelir.
NULL'a karsi korunan tuketiciler: characterize_leaves, make_explorer,
make_report, extract_representatives.

### (b) MIN_SPLIT_SIZE 8 -> 4, gerekcesi olculebilirlik  [DUZELTILDI]

Madde 99'da olculmustu: heterojen 420 yapragin 420'si boyut tabaninda durdu,
hicbiri esik listesini tuketmedi. Boy dagilimi olculdu: 2 uyeli 154, 3 uyeli
104, 4-5 uyeli 100, 6-7 uyeli 62. Taban 4'e indirildi, yani **162 yaprak
(840 uye)** artik homojenlik hedefiyle karar aliyor.

Taban 2'ye INDIRILMEDI, ve sebebi (a): 4'un altinda bolme en az bir tek uyeli
cocuk uretmek zorundadir (2->1+1, 3->1+2) ve tek uyeli dugumun ikili kimligi
yoktur. Yani kalan 258 yaprakta "too_small_to_split" artik DOGRU bir gerekce --
asiri parcalanma kaygisi degil, olculebilirligin siniri.

### (c) 3/3 kapisi uyelikte karar veriyor -- silme olarak DEGIL  [DUZELTILDI]

about.html soz vermisti: "tightening admission to three of three is the next
step and will change every denominator". Soz TUTULDU ama farkli bir sekilde, ve
sebep literaturun kendi sinirinda:

  triad_gate.py basligi: "Her Rieske oksijenazda triad vardir diyen bir kaynak
  BULUNAMADI (arandi, 2026-10-30). Derlemeler kanonik ve korunmus diyor,
  istisnasiz demiyor."

Eksik ligandli 55 girisi `is_confirmed=0` yapmak "bunlar RO degil" demek olurdu
-- literaturun desteklemedigi bir iddia. O yuzden 3/3 karari uyelige TEMEL
olarak girdi: yeni `ro_membership(candidate_id, ro_cluster, basis, detail)`,
uc temel ve her birinin kaynagindaki aciklamasi (`MEMBERSHIP_BASES`):

  triad_complete              uc ligand da tipin korudugu kolonlarda okundu
  documented_exception        hizalama katalitik loba varamiyor, ligandlar
                              baska yontemle dogrulandi (motif_exceptions.yaml)
  triad_incomplete_retained   ligand gercekten eksik; giris ATILMAZ, cunku
                              yokluk "oksijenaz degil" kaniti degildir

Yayina giren yerler: entry sayfasindaki KOSULSUZ iddia kaldirildi -- orada her
giris icin "carries an intact Rieske centre and catalytic triad" yaziyordu,
oysa ro_triad'a gore bu 55 giriste YANLIS. Cumle artik girisin kendi okumasindan
kuruluyor, hem `ro.html`'de hem statik yoldaki `entry.js`'te (shard'a `mb`
alani eklendi). about.html'in iki paragrafi yeniden yazildi. Kapiya
`check_membership_basis` eklendi: satir sayisi, yetim giris, ve temelin kendi
triad okumasindan AYRISMASI kontrol edilir.

### (d) empty_types tek yonlu okunuyordu  [DUZELTILDI]

`empty_types.json` 8 uyesiz tipi tam cozumlemis ve sayfalarinda yayinliyordu
("counted under CARDO"). Ama TERS yon hic yoktu: CARDO sayfasini acan okuyucu,
ayri yayinlanmis bir karbazol dioksijenazin (CarAa) o 111 uyenin icinde
sayildigini gormuyordu. `atlas.absorbs_view` + cluster.html blogu eklendi.

Olculen tablo (uyesiz tip -> yutan tip, referans kimligi, o tipin uye sayisi):
NahAc -> NBDO %100,0 / 24;  CarAa -> CARDO %99,7 / 111;  IpBAa -> cumA1 %90,8 / 8;
BpdC1 -> BphA1A2 %79,2 / 21;  XylC1 -> BphA1 %67,0 / 1;  PbaA1 -> DxnA1 %65,6 / 2.
(DntAc ve NinA ikisi de uyesiz, birbirine %93,1.)

NAFTALIN DIOKSIJENAZIN bu veritabaninda sifir uyesi var ve sebebi budur:
atama rekabetli, ve NBDO profili her yarismayi kazaniyor.

### (e) Kendi dun yazdigim ref_identity kimligi SISIRIYORDU  [DUZELTILDI]

Madde 102'de `cluster_reference_identity.csv` eklendi ve yalnizca iki tarafin
da kalinti tasidigi kolonlari sayiyordu -- yani ciftin FARKLI oldugu indel
kolonlarini paydadan da paydan da atiyordu. Olculdu: CARDO/CarAa %98,8 yerine
**%95,5**, qxyA/BedC1 %23,0 yerine **%19,3**. Etki buyuk degil ama tek yonlu.
Payda artik "en az birinde kalinti olan kolonlar", bosluk uyusmazlik sayilir.

### (f) Sifir sayili kanit duzeyi satiri gizleniyordu  [DUZELTILDI]

cluster.html `{% if tiers[t] %}` ile satiri atliyordu, yani okuyucu 0 ile
"olculmedi"yi ayirt edemiyordu: qxyA sayfasinda characterized=1 gorunup
close_homolog satiri hic cikmiyordu. Cubukta 0 segmenti YOK (genisligi olmayan
dilim anlamsiz), tabloda VAR (sayinin kendisi bir bulgu).

### (g) Ev odevi: diamond ara dosyalari /tmp'den cikti  [DUZELTILDI]

build_phylogeny.py 1944 dizi icin all-vs-all yapiyor ve `tempfile` varsayilani
/tmp idi -- projenin duran kuralina aykiri. `--work` eklendi, zincirde
`$ROAR_SCRATCH/roar_phylo_work`.

### Kapatilmadi

- Domain-farkinda hizalama. AerC'de olculdu: tam boy hmmalign modelin %18,9'unu
  kapliyor ve kolon 231'de duruyor, katalitik kolonlar 289/294/431'e hic
  varmiyor. Dogru cozum bu (motif_exceptions.yaml de boyle yaziyor), ama
  hizalama adimi degisirse motif modelinin ciktisi, dolayisiyla DB'deki her
  sayi yeniden kurulur. Bu zincire karistirilmadi: ayni kosuda iki kok degisim
  olursa hangi degisimin hangi sayiyi oynattigi ayirt edilemez.
- Adlandirma: NahAc/NBDO ve CarAa/CARDO ayri tip mi, tek tip mi. Olcum yapildi,
  karar literatur/adlandirma karari ve burada verilmiyor.

## 103. Bug avi: yayimlanan sayi ile basliginin uyusmadigi 10 yer (2026-10-31)

Madde 102(a)'daki hata bicimi -- "yayimlanan bir sayi, onu ureten hesabi degil
BASKA bir hesabi tarif eden bir baslik altinda duruyor" -- sistematik olarak
tarandi. 10 vaka bulundu, 9'u duzeltildi. Her biri bagimsiz yeniden hesapla
dogrulandi.

DUZELTILDI:

1. ENTRY SAYFALARINDA EN AGIR VAKA. `ro.html` `core_identity`'yi "alignment
   identity to nearest reference" diye etiketliyordu, ama
   discover_subfamilies.py onu kumenin KENDI alt ailelerinin kolon
   uzlasilarina karsi hesapliyor -- hicbir kuratorlu referansa dokunmuyor.
   Bir CntA girisi 0,99 gosterirken en yakin referansa gercek kimligi %24,4.
   DB genelinde ortalama %72,9'a karsi gercek %46,7 (n=19.125). Ustelik
   GERCEK sayi ayni sayfada yukarida zaten vardi (ev.ref_identity), yani ayni
   adi tasiyan iki sayi 70 puan ayriliyordu. Etiket duzeltildi: "identity to
   its subfamily consensus", ve ne OLMADIGI yaziliyor.

2. ANA SAYFA TABLO 1'DE BILIMSEL OLARAK YANLIS IFADE. "Centre" kolonu
   [85, 87, 105, 108] ve 209 kolonlarini SABIT yaziyordu; o numaralar 426
   kolonluk ESKI profilden kalmaydi, model 603 kolona gecince hicbiri
   eslesemez oldu ve SEKIZ SATIRIN HEPSI else dalina dusup "mononuclear
   Fe(II)" basmaya basladi. Yani sayfa Rieske sisteinlerinin mononukleer
   demiri koordine ettigini soyluyordu. Literaturle celisir: 1NDO'da [2Fe-2S]
   kumesini Cys81/His83/Cys101/His104, mononukleer demiri His208/His213/
   Asp362 tutar. Merkez artik motif_stats.json'daki `label` alanindan
   turetiliyor (orada ZATEN dogru duruyordu), kolon numarasindan degil.
   Dogrulandi: 164/166/184/187 -> Rieske, 286 -> subunit bridge,
   289/294/431 -> mononuclear Fe(II).

3. "20 lack a histidine" ELLE YAZILIYDU ve yanlisti. Gercek 31 (madde 100'de
   capalama yeniden kostuktan sonra). Cumle 19.063+24+20+7 = 19.114'e
   toplaniyordu ama 19.125'ten soz ediyordu; 31 ile tam kapaniyor. totals()'a
   `triad_his_absent` eklendi -- bes sayinin dordu zaten canliydi.

4. "DISTINCT ORGANISMS" TUR DEGIL SUS SAYIYORDU. COUNT(DISTINCT p.organism)
   GenBank organizma dizgisini sayiyor, yani ayni turun on izolati on
   "organizma". Ustelik kolonun kendi basligi "yogun dizilenen bir habitat
   ekoloji biriktirmeden giris biriktirir" diyerek tam bu artefakti
   onledigini soyluyordu. Olculdu: 4_411_XylX rizosfer topragi 28 -> 8 tur,
   insan klinik 66 -> 26, toprak 71 -> 37. Projenin kendi tanimi
   species_of() (ilk iki kelime) kullaniliyor; atlas.species_of eklendi ve
   isolation_source'takiyle ayni sonucu verdigi test edildi.

5. PROMOTER ARALIGI, BASLIGIN DISLADIGI YONELIMLERI DE SAYIYORDU.
   control_elements.py mimariyi hic filtrelemiyordu, oysa baslik "the median
   length of THAT intergenic region" (divergent olani) diyor; yanindaki
   kolon divergent PAYINI gosterdigi icin komsu iki hucre farkli paydadan
   okunuyordu. Olculdu: AraC 277,5 -> 310,0 (314 girisin 174'u divergent),
   IclR 220,5 -> 178,5, TetR 146,0 -> 129,5, LysR 144,0 -> 138,0.
   NOT: etkili olmasi icin control_elements adiminin yeniden kosmasi gerek.

6. "N types consist mostly of distant relatives" = 14, ama kod kurali
   ">=%80 distant/novel VE >=20 uye". Basligin kendi tanimiyla (uyelerin
   cogunlugu, boyut suzgeci yok) 68 tipin 36'si. Baslik gercek kurali
   soyleyecek sekilde duzeltildi.

7. "are at least 60 % identical to THAT reference" -- kademe aslinda 76
   referansin TUMUNE karsi en iyi vurus (evidence_tiers.py, diamond -k 1),
   bu tipin referansina degil. Kucuk tiplerde hata 2 kat: 3_306_3NTDO
   %100 (2/2) yayimlanmis, kendi referansina gore %50 (1/2). Baslik
   "some characterised reference enzyme" olarak duzeltildi.

8. WASTEWATER TABLOSUNDA IKI "overall" aslinda yalnizca habitat-
   siniflandirilmis girisler uzerinden (6.434/19.125). Yakinlik payi %27,7
   yayimlanmis, tum DB'de %26,0; medyan kimlik %38,0'a karsi %37,2. Baslik
   "across all habitat-classified entries" oldu.

9. "2,388 entries" (3 yerde) artik uremiyor; olculen deger 2.387
   (carboxylate_audit.json gate_cost, motif_stats.json offset_distribution ve
   19.125-16.738 ucu de 2.387 diyor). "45 types" dogruydu, dokunulmadi.

10. TAHMINLER SAYFASINDA BOS RENDER. `identity_to_curated_refs` FILUMA gore
    anahtarli ama sablon tek bir median/min/max bekliyordu; yayimlanan sayfa
    harfiyen "median  %, range – %" yaziyordu. Filum basina basiliyor
    (Halobacteria 28,4 / 23,5-43,7, Thermoproteota 37,0 / 29,9-40,5,
    Thermoplasmata 40,2) ve esik sayilari da yaziliyor, cunku bilimsel ifade
    onlarda: %50 uzerinde 0, %40 altinda 85.

AJANIN DOGRULADIGI VE DOGRU CIKAN SAYILAR (kisa liste, tamami raporda):
atlas_network heatmap basligi (hala max-over-representatives diyor, dogru),
gen sirasi bloklarinin kademe suzgeci, kanonik operon seklinin tek paydasi
(madde 101b'nin sonucu), habitat kapsami paydasi, operon ortak paylari,
atlas_home'un "a characterised enzyme" ifadesi, atlas_closeness satir
etiketleri, about'taki %99,8 / %32,4 referans cifti ucleri, ana sayfa
koprulu karboksilat tablosunun paydalari, tip basi karboksilat paylari,
atlas_predictions arkea/mantar sayilari, atlas_ecology tur kolonlari.

## 102. "En yakin tipler" yanlis sayiyi yayimliyordu; tazelik izleyicisinde iki sessiz hata (2026-10-31)

### (a) BENZERLIK TABLOSU IKI TIPIN BENZERLIGINI OLCMUYORDU

qxyA sayfasi BedC1'e %91,4 benzer diyordu. Yer gercegine karsi olculdu:

    cift            sayfa    ikili mafft   ortak hizalama
    qxyA-BedC1      91,4     27,4          23,0
    qxyA-bfzA       90,2     28,2          21,9
    qxyA-CntA       84,7     30,1          24,8

Sebep hesap hatasi DEGIL, anlam hatasi. build_phylogeny.py'de
`best_pair[k2] = max(best_pair[k2], pid)`: deger iki tipin referanslari
arasindaki kimlik degil, IKI TIPE ATANMIS HERHANGI IKI TEMSILCI arasindaki en
yuksek kimlik. %91,4'u ureten cift bulundu: 2_212_BedC1#5 ile 5_505_qxyA#21,
ikisi de TEK UYELI yaprak, ikisi de `distant`.

YAPILDI:
  - yeni `cluster_reference_identity.csv`: 76 referansin ortak HMM
    hizalamasindan referans-referans kimligi. Sayfadaki "Nearest other types"
    artik buradan. qxyA icin sonuc biyolojik olarak anlamli: Stc2 %32,0,
    BmoA %28,6, GbcA %26,2, CntA %24,8, CmoA %24,1 -- hepsi grup 5, yani
    qxyA'nin kendi kuaterner amin kolu.
  - NEDEN SSN KENARLARINDAN DEGIL: diamond'in qcov>=70 / e<=1e-5 filtresi uzak
    referans ciftlerini TAMAMEN atiyor; qxyA-BedC1, qxyA-bfzA, qxyA-CntA icin
    REF-REF kenari YOK. SSN bu sayiyi yapisal olarak bildiremez.
  - NEDEN IKILI MAFFT DEGIL: 76x76=2850 hizalama bir adim icin fazla yavas.
    Bedeli olculdu: coklu hizalamadan okunan ikili kimlik ~5 puan DUSUK
    cikiyor. Yanilma yonu guvenli taraf (az gosteriyor, cok degil) ve bu
    basliga yazildi.
  - eski sayi SILINMEDI, adi ve anlami duzeltildi: `cluster_assignment_
    conflicts.csv` + sayfada "Members that collide with another type", her
    satirda iki dizinin kimligi VE KADEMELERI.

OLCUM, esik karari icin belirleyici: tum DB'de >=%70 tip-arasi catisma 131.
Bunlarin %93,1'inde (122) iki taraf da distant/novel/family; 8'inde biri
saglam; ve YALNIZCA 1'inde (%0,8) iki taraf da characterized/close_homolog --
o da BmoA(close) / GbcA(characterized) %71,2, gercekten yakin akraba iki
kuaterner amin oksijenazi, hata degil. SONUC: tip sinirlari iyi kanitlanmis
uyeler duzeyinde SAGLAM; catismalar neredeyse tamamen `distant` uyeleri kabul
etmenin yan urunu.

### (b) TAZELIK IZLEYICISINDE IKI SESSIZ HATA (madde 101a'nin devami)

Madde 101'deki duzeltme ISE YARAMADI; sebebi iki ayri sessiz basarisizlik:

  1. `db_writes` GOZLEMLE ogrenilmiyordu. Tablo yazmalari satir sayisi
     degisiminden tespit ediliyor; triad_gate ro_triad'i ayni 19.125 satirla,
     ingest_protein_entries ro_protein'i ayni 73.699 satirla yeniden yaziyor,
     yani sayim degismiyor ve "yazma gozlemlenmedi" deniyor. Statik kaynak
     taramasi da goremedi, cunku triad_converge.py surucusu yazmalari kendi
     icine sakliyor (cmd[1] surucuyu gosteriyor). Cozum: `DECLARED_WRITES`.
  2. MTIME, KENDISINDEN USTUN OLAN OZETI EZIYORDU. `staleness()` mtime
     dalinda DONUYORDU, icerik karsilastirmasina hic varilmiyordu -- oysa
     rebuild.py'nin kendi basligi ucuncu kurali tam bu durum icin koymus.
     Olculdu: ROmotif71.hmm bir ciktidan yeni gorunuyor ama ozeti annotate'in
     kaydettiginin AYNISI (6fc01b561ae6cb7e). Bu tek sahte kok 48 artefakti
     bayat gosteriyordu. Kural ZAYIFLATILMADI: taban kaydi olmayan ya da ozeti
     DEGISEN girdi hala bayat.

SONUC: 6 guncel -> 16 guncel. Kalan 38 bayatligin tek GERCEK koku var:
recursive_homogenize.py'nin icerigi degisti (madde 99'daki olcum oraya
yazildi). Ozet yorumu koddan ayirt edemez; bu dogru muhafazakar davranis ve
kurcalanmadi.

### (c) ACIK KARAR: %70 atama esigi (kullanici onerdi)

    esik          kalan onayli giris    kalan tip
    simdi                    19.125           68
    >=%40                     8.451            -
    >=%60                     4.970            -
    >=%70                     3.486           51

%70 uygulanirsa 17 tip TAMAMEN bosalir ve hicbiri marjinal degil: CmoS 456
uye (en iyi %57,1), CbdC 424 (%69,4 -- kil payi), DdmC 401 (%61,7), AtdA 182,
PhnA1a 146, PhdA 140. Toplam 2.051 uye ve 17 tip haritadan siliniyor;
genelde 19.125 -> 3.486, yani %82 kayip.

ONERI (uygulanmadi, karar kullanicinin): uyeleri UYELIKTEN atmak yerine
`distant`i TURETILMIS ISTATISTIKLERDEN cikarmak -- bugun duzenleyici ve gen
sirasi icin yapilan sey. (a)'daki olcum bunu destekliyor: sinirlar saglam
uyeler duzeyinde zaten saglam, bozulan sey yalnizca distant uyelerden
hesaplanan sayilar. Boylece tip haritasi kaybolmaz ama hicbir yayimlanan sayi
"bu enzim olmayan" dizilerden gelmez.

KUCUK NOT: build_phylogeny.py diamond all-vs-all icin hala tempfile ile
/tmp kullaniyor (1.944 dizi, kendini temizliyor). Proje kurali geregi
SCRATCH'e tasinmali; buyuk degil, aciklikta tutuluyor.

## 101. Kapi betikleri zincire girdi; duzenleyici ve 1.2 tablosu duzeltildi (2026-10-31)

(a) TAZELIK IZLEYICISINDEKI BOSLUK KAPATILDI (madde 100b). `carboxylate_audit`,
`triad_gate`, `protein_entries` artik ARTEFACTS adimi; sira
`... learning -> carboxylate_audit -> triad_gate -> protein_entries ->
phylogeny -> ... -> validation_gate`. Boylece yazdiklari tablolar onceki
adimlarin `ignore` kumesine giriyor ve 45 adimin SAHTE bayatligi bitiyor.
Tazelik kontrolu ZAYIFLATILMADI; ro_triad hala izleniyor.

Sabit nokta icin yeni surucu `triad_converge.py`: gate -> say ->
his_needs_anchoring 0 degilse anchor -> basa don. Gerekcesi dosyanin basinda.
`SCRATCH` sabiti eklendi (varsayilan /media/lin-bio/back2): bu iki betik
--work almazsa 19.125 dizinin hizalamasini /tmp'ye yaziyordu.

NOT: `db_writes` ancak adim rebuild.py uzerinden bir kez KOSTUKTAN sonra
olculup kaydediliyor, yani duzeltme ilk kosudan sonra etkili.

(b) DUZENLEYICI BOLUMU IKI AYRI PAYDAYI KARISTIRIYORDU. Kume sayfasindaki
sekil layout'u YAKIN uyelerden (characterized/close/family) aliyor, mimariyi
ise regulation_by_cluster.csv'den yani TUM uyelerden aliyordu. qxyA'da
olculdu: 164 uyede divergent_regulator 14, yakin 53 uyede YALNIZCA 1. Ustelik
tum uyelerde en sik aile "regulator_unclassified" oldugu icin sekil TetR'i
HIC YAZAMIYORDU. Artik mimari de yakin temelde sayiliyor ve kac/kac oldugu
yaziliyor.

Ayrica referansin KENDI mimarisi ayri bir sekil olarak eklendi. Gerekce: tip
geneli sayi ile literaturun tarif ettigi yapi ayni sey degil. qxyA'da tip
geneli korunmus bir duzenleyici BILDIRMIYOR (1/53), ama qxyAR operonu
referansta gercek ve belgeli (PMID:27792326): ters iplikteki TetR, operonla
promoter bolgesi arasinda 78 bp'de duruyor, yani tek bir duzenleyici bolge
ikisini birlikte kontrol edebiliyor. Ters iplikteki gen ayni transkriptte
olamadigi icin gen sirasi sekline GIRMEZ -- bu dogru davranis, ve sayfa artik
nedenini yaziyor. Iki iddia ayri gosteriliyor.

DUZELTME: bu oturumda once "operon siniflandiricisi bozuk" dedim, YANLISTI.
`same_strand=0` olan TetR layout'tan dogru dislaniyor, MFS'ler de 718/1973
bp ile 150 bp'lik operon penceresinin disinda. Sistem dogru calisiyordu;
eksik olan sey sayfanin bunu SOYLEMEMESIYDI.

(c) ANA SAYFA 1.2 TABLOSU TRANSPOZE EDILDI. Baslik 7 kolon vaat ediyordu ama
katalitik satir bunlarin dordunu colspan ile duzyaziyla dolduruyordu, cunku o
dort olcum katalitik karboksilat icin YOK (kolon 431 ligandi okumuyor: 26/26
yapida demirden 10,9-14,7 A uzakta). Sonuc: basligin satirlardan birini
tanimlamadigi bir tablo. Artik satirlar OZELLIK, kolonlar iki karboksilat;
colspan kalmadi, dort olculemeyen hucre tek bir rowspan ile tek ifade
tasiyor. Hicbir sayi ve hicbir cekince kaldirilmadi.

(d) OLCULDU, SIKAYET DOGRULANMADI: "sitenin cogu guncel degil". Yayindaki
qxyA sayfasi DB ile birebir ayni (characterized 1, family member 52, distant
104, novel 7) ve TetR sayfada mevcut. Muhtemel sebep: 3.241 dosyalik dagitim
yayilirken tarayici onbellegi.

## 100. Kapi zinciri kostu; tazelik izleyicisinde yapisal bosluk (2026-10-30)

MADDE 98'DEKI SIRA KOSTU, KAPI YESIL: 113 kontrol, 109 PASS, 0 FAIL, 4 WARN.
`ro_triad` 19.125 satir; complete 19.063, his_absent 31, carboxylate_absent 24,
carboxylate_untested (AerC) 7. `ro_protein` 73.699 satir.

IKI SEY OGRENILDI.

(a) CAPALAMA IKI TUR DONUYOR, ve bu tasarimdan. Birinci turdan sonra
`his_needs_anchoring` 12 cikti: tam olarak yeni onaylanan +12 giris, cunku
`triad_anchor_his.py` isini BAYAT `ro_triad`'dan okudugu icin onlari gormemisti.
Ikinci tur 0'a indirdi. Secimin `his_source` uzerinden yapilmasi (triad_status
uzerinden DEGIL) bu yakinsamayi mumkun kiliyor; durum uzerinden secilse ikinci
tur bos kume donerdi. KURAL: madde 98'in 2-3. adimlari, his_needs_anchoring 0
olana kadar TEKRARLANIR.

Ayrica `triad_anchor_his.py` KeyError ile dusuyordu: `entries` guncel onayli
kumeden, `triad` bayat `ro_triad`'dan okunuyor, ve eski kurulumda onayli olup
artik olmayan bir giris (XM_009028827.1:0-1056:1) ikincisinde var birincisinde
yok. Artik atlaniyor ve kac tane atlandigi yaziliyor -- DB'de uye olmayan bir
dizi icin triad cagrisi yazmak, olcmedigimiz seyi yayimlamak olurdu.

(b) TAZELIK IZLEYICISI 45 ADIMI SAHTE BAYAT GOSTERIYOR. `rebuild.py --check`
51 artefaktin 45'ini BAYAT dedi, sebep olarak `roar.sqlite[ro_triad,
ro_triad_meta]` gosterdi. OLCUM: 52 adim betiginden YALNIZCA
`validate_curation.py` bu tablolari okuyor (o da validation_gate adimi, yani
onun bayat olmasi DOGRU). Diger 45 adim bu tablolari HIC okumuyor.

Neden: bu tablolari yazan `triad_gate.py` ve `ingest_protein_entries.py`
ARTEFACTS listesinde YOK, elle kosuyoruz. Dolayisiyla hicbir adimin `writes`
kumesinde degiller, hicbirinin `ignore` kumesine girmiyorlar, ve `roar.sqlite`i
girdi sayan her adim onlari "degismis ust girdi" saniyor. Zincirin 42/47'de
dusmesinin koku de bu. Bu kosuda 45 adim KOSTURULMADI (saatler surer, ciktiyi
degistirmezdi); yalnizca gercekten hic kosmamis 5 adim kosturuldu
(provenance, disagreements, report_html, explorer_html, hub_html; toplam 9 s).

ACIK KARAR (kullaniciya): dogru duzeltme bu betikleri zincire adim olarak
eklemek, ama arada GERCEK bir sabit nokta var: `triad_anchor_his.py`
`ro_triad`i OKUYOR, `triad_gate.py` onun JSON'undan `ro_triad`i YAZIYOR. Yani
dogrusal sirada gate -> anchor -> gate olmak zorunda, (a)'da elle yapilan tur.
Bunu pipeline'a gommek yanlis sirayi kalici hale getirebilir, o yuzden
DEGISTIRILMEDI. Tazelik kontrolu ZAYIFLATILMAMALI: cozum bu uc betigin
sirasini dogru yazmak, ro_triad'i izlemekten cikarmak degil.

## 99. Heterojen yapraklarda karari veren kural yayimlandi (2026-10-30)

`recursive_homogenize.py` yapraklari medyan ikili kimlik 0,70'e (HOMOGENEITY_TARGET)
ulasana kadar artan CD-HIT esikleriyle (0,50 ... 0,90) boluyor. Sayfa heterojen
yapraga yalin bir "heterogeneous" etiketi basiyordu ve bu YANILTICIYDI.

OLCUM (2026-10-30, 2.789 yaprak / 19.125 onayli uye):
  - homojen 2.369, heterojen 420 (1.460 uye)
  - heterojen 420'nin **420'si** `size < MIN_SPLIT_SIZE (8)` yuzunden durdu
  - **hicbiri** esik listesini tuketmedi (en derin adim 4, liste 5 uzunlugunda)
  - bu yapraklarin ortalama medyan kimligi 0,603, en dususu 0,425
  - tip dagilimi basi: CntA 64 yaprak, TsaM 42, BmoA 32, Stc2 30, CdnA 23

Yani "homojenlik hedefi tutturulamadi" denen yerde hedef hic SORULMAMIS. Bu,
projenin tekrar eden hata bicimi: karar veren kural, karar verdigini
soyledigimiz kural degil (triad aralik sarti, kolon 471/431 vakasi).

YAPILDI, kumelenme DEGISTIRILMEDI:
  - `webapp/atlas.py` `leaf_stop_reason()` + `LEAF_STOP_REASONS` (sabitler
    kaynagindan kopya, kaynak hala recursive_homogenize.py)
  - `leaf.html` / `cluster.html` durma nedenini yaziyor; yaprak sayfasi artik
    temsilcinin birbirinden ayrilan uyeleri temsil ettigini soyluyor
  - `validate_curation.py` `check_leaf_stop_reason`: sabitler ayrisirsa ve
    durma nedeni turetilemeyen bir yaprak cikarsa FAIL (113 kontrol, PASS)
  - `recursive_homogenize.py` MIN_SPLIT_SIZE yorumuna olcum islendi

ACIK KARAR (kullaniciya): `MIN_SPLIT_SIZE` 8'den dusurulurse homojenlik
GERCEKTEN karar verici olur ve 420 yaprak varyant duzeyinde cozulur. Bedeli
tum yapraklarin yeniden kurulmasi: leaf_id'ler degisir, yani agac,
leaf_profile, leaf_sdp ve kalici /leaf/ baglantilari da degisir. Bilincli
karar gerektirir, o yuzden bu kosuda DEGISTIRILMEDI.

LITERATUR CERCEVESI: yaprak taksonomik bir birim DEGIL. Literaturde bir RO
"varyantini" tanimlayan yuzde esigi yok -- Kweon 2008 semasi elektron
ortaklarina, `ro_group` alfa alt birimi dizi kumelenmesine dayanir
(citations.yaml kweon2008_classification notu). Yaprak bunlarin hicbiri
degil, agac boyutunu kisan islevsel bir gruplamadir ve sayfada oyle okunuyor.

## 120. A5.3 dizi eslesme kapisi -- "ayni kayitta gecer" ile "ayni proteindir" ayri seylerdir (2026-11-04)

A5'te 13.079 nukleotid kaydi indirildi ve 18.636 (dizi, kayit) cifti NCBI
capraz referansindan geldi. Capraz referans "bu dizi su kayitta gecer" der.
Ucuz ve cazip olan, kaydi acip oradaki bir RO adayinin komsulugunu o diziye
yazmakti. OLCULDU, ve olcum o yolu kapatti.

KURAL: bir kayit bir diziye ancak AYNI KAYITTA cevirisinin md5'i BIREBIR
esit bir aday varsa baglam olur. Betik: `a5_match_guard.py`, cikti
`analysis_out/a5_match_guard.json`, kapi `validate_curation.py` icinde
`check_a5_match_guard` (butce KAYNAKTAN okunur, JSON ile karsilastirilir).

```
hedef                          18.636 cift / 18.532 tekil dizi
ayni kayitta birebir eslesme   18.139  %97,3   -> baglam kurulur
yalnizca baska kayitta eslesme      9  %0,0    -> kurulmaz
hic eslesme yok                   488  %2,6
BAGLAM KAZANAN TEKIL DIZI      18.043/18.532  %97,4
```

Kalan 497 cift nedenine gore ayrildi; dagilim yayimlaniyor cunku "%2,6
eksik" ile "%2,6 yok" ayni sey degil:

```
A kayit indirilemedi                16   %3,2
B kayitta hic CDS anotasyonu yok    50  %10,1   <- nukleotid VAR
C anotasyonlu ama RO adayi yok      52  %10,5
D adayi var, ceviri tutmuyor       379  %76,3
```

D'nin ici, kurali gevsetme karari icin belirleyici oldu:

```
ortak parca BILE yok -> baska protein  114  %30,1
ic parca ortak, uclar farkli           105  %27,7
aday hedefin icinde                     89  %23,5
hedef adayin icinde                     71  %18,7
```

KURALIN KAZANCI: "ayni kayit yeter" denseydi en az **114 dizi baska bir
genin komsulugu ile yayimlanacakti** ve sayfada bunu gosterecek hicbir sey
olmayacakti. Bu sayi `a5_match_guard.py` icinde `DIFFERENT_PROTEIN` olarak
beyan edilir ve kapi onu olcumle karsilastirir; ayrisirsa FAIL.

ACIK: icerme oranina dayali bir IKINCI kural (ornegin %90 ortusme) 265
cifti kurtarabilir. Ama bu ayri bir kuraldir: ilan edilmeden, olculmeden ve
kendi kapisi kurulmadan uygulanmaz. Bu surumde kural birebir md5'tir.

## 121. Anotasyonsuz ama nukleotidi olan kayitlar -- "baglami yok" degil, "baglamini biz cagirmadik" (2026-11-04)

Kullanicinin tespiti: genomik anlamda anotasyonsuz ama genomik baglami
BULUNAN diziler olmasi mumkun. Olculdu ve dogru cikti.

```
anotasyonsuz replikon (cds_count=0)   47   (A5 46 + survey 1)
bu kayitlara dusen hedef dizi         50
kalan eslesmeyen ciftin payi        %10,1   (sinif B)
```

Bu 50 dizi icin kayitta NUKLEOTID VAR; eksik olan yalnizca CDS anotasyonu.
Yani sayfanin "genomik baglam yok" demesi bir BIYOLOJI ifadesi degil, bir
ARAC ifadesi -- ve ikisi ayni cumlede birlesmemeli.

PLAN: bu 47 kayda ORF cagrimi (prodigal) uygulanir, cikan ORF'ler
`extract_genomic_context.py`e verilir ve md5 kapisi AYNEN uygulanir --
cagrilan ORF'un cevirisi hedef diziyle birebir tutmak zorunda. Tutmazsa
baglam yine kurulmaz; ORF cagrimi kapiyi gevsetmek icin bir bahane degil.

SART: provenans AYRI bir alan olur (`annotation_source` = `ncbi` | `prodigal`).
NCBI anotasyonu ile cagrilmis ORF ayni sutuna sessizce karismaz; konuma
bagli her olcum (operon, komsuluk) kaynagina gore ayrilabilir olmali, cunku
cagrilmis bir komsulugun gen sinirlari olculmus degil TAHMIN edilmistir.
Bu, 71 tip semasina sutun eklemek degil, turetilmis tabloya alan eklemektir
(bkz. sema genisletme kurali).

BUYUKLUK: 50 dizi, veritabaninin %0,3'u. Kucuk; ama "yok" ile "bakmadik"
ayrimi buyukluge bagli degil.

## 122. Bilesen taramasi 2026-10-05'ten beri bayatmis -- varlik kontrolu, tazelik kontrolu DEGILDIR (2026-11-04)

Bu madde BIR KEZ YANLIS YAZILDI ve duzeltildi; yanlis sürüm burada
kaliyor cunku nasil yanildigi bulgunun kendisi kadar onemli.

ILK (YANLIS) IDDIA. Kullanicinin "anotasyon yanlis yapilmis da olabiliyor"
tespitini arastirirken su olculdu: anotasyonun kendisi "dioxygenase small
subunit" veya "beta subunit" dedigi 1.529 tekil komsu proteinden 841'i
`Ring_hydroxyl_B` tutuyor, 688'i HIC domain almiyor. Bundan "beta tespiti
tek Pfam modeline dayaniyor ve yanlis-negatif orani ~%45" sonucu cikarildi
ve kendi beta HMM'ini kurup sinayan bir betik yazildi. Ayrik test kumesinde
628/628 geri kazanildi, en kotu E-degeri 9,9e-62.

SONUC FAZLA TEMIZ OLDUGU ICIN SINANDI. Pfam'in KENDI `Ring_hydroxyl_B`
modeli ayni 628 dizide dogrudan kosturuldu: **628/628 isabet, E <= 1e-5.**
Yani Pfam onlari kacirmiyor. Domain satirlarinin olmamasinin nedeni baska:
o proteinler HIC ARANMAMIS.

GERCEK HATA. `build_operons.py:run_component_hmms` icinde iki varlik
kontrolu vardi ve ikisi de yanlis soruyu soruyordu:

```
1. if existing and not force: return        -> neighbor_domain'de TEK satir
                                               varsa fonksiyonun TAMAMI atlanir
2. if not os.path.exists(domtbl)            -> domtbl dosyasi VARSA arama atlanir
```

Diskteki kanit:

```
operon_work/neighbor_components.domtbl   2026-10-05   13.091 isabetli protein
operon_work/neighbor_proteins.fasta      2026-10-14   (daha YENI)
roar.sqlite                              2026-11-03   316.120 komsu protein
```

Yani komsu protein kumesi yenilenmis, bilesen taramasi yenilenmemis.
2026-10-05'ten sonra eklenen her komsu protein taramayi GORMEDI ve `none`
olarak siniflandi -- beta, ferredoksin ve reduktaz dahil. 628 beta alt
birimi domtbl'de hic gecmiyor.

AYNI TUZAK, IKINCI KEZ. `annotate_ro.py` bu olayin birebir aynisini
yasamis, cozmus ve yorumuna yazmis: "iki ara dosya 21 Temmuz'dan kalmisti,
aday fastasi 9 Ekim'de yeniden uretilmisti... DOGRULANMIS RO: 0 yazdi.
Hicbir sey hata vermedi." Duzeltme O dosyaya uygulanmis, AYNI kurali
uygulayan `build_operons.py`e uygulanmamisti.

DUZELTME. Kural artik kopyalanmiyor: `stale_intermediate.py` tek kaynak ve
hem `annotate_ro.py` hem `build_operons.py` oradan aliyor. (`rebuild.py`
ve `provenance.py`deki bayatlik mantiklari AYRI isler yapar ve bilincli
olarak birlestirilmedi.) Ayrica fasta her kosuda yeniden yazilirsa mtime'i
degisir ve domtbl HER ZAMAN bayat gorunur -- tuzagin tersi; o yuzden fasta
gecici dosyaya yazilip ancak OZETI degistiyse yerine konuyor
("mtime tetikler, ozet karar verir").

DERS. Kullanicinin tespiti dogruydu, ama hatanin yeri degildi: anotasyon
dogru, Pfam dogru, bizim boru hattimiz yanlis. Ve "fazla temiz" bir sonuc
(628/628, E=1e-62) bir basari degil bir UYARI isaretiydi -- o yuzden
yayimlanmadan once bagimsiz olarak sinandi.

MALIYETI OLCULDU. Tarama dogru kosturulunca:

```
bilesen        once     sonra    degisim
beta          2.523     4.372    +%73
reduktaz      6.185    10.394    +%68
ferredoksin   1.624     2.736    +%68
alpha_other   1.407     2.332    +%66
domain isabetli protein  13.091 -> 22.302   +%70
```

Operon bayraklari: has_beta 2.502 -> 4.333, has_reductase 3.203 -> 5.409,
has_ferredoxin 1.202 -> 1.989. Yani yayimlanan operon paylari gercek
bilesenlerin ~%40'ini kaciriyordu.

VE SONUC BIYOLOJIYE OTURDU -- duzeltmenin en guclu kaniti bu:

```
grup  once beta yok   sonra    yorum
4         %57,4       %25,2    bifenil tipi; beta artik cogunda bulunuyor
2         %53,8       %23,8    ayni
3         %83,8       %72,2    hala yuksek -- ACIK
5        %100,0      %100,0    CntA/CntB tipi, betasi YOK; degismemesi DOGRU
1         %100       %99,9     CarAa tipi, betasi YOK; dogru
6        %100,0      %100,0    tek gen agirlikli
```

Beta yoklugu beklendigi gruplarda degismedi, beklenmedigi gruplarda dustu.
Rastgele bir duzeltme boyle davranmaz.

KAPI. `check_component_scan_is_fresh` IKI ayakli: (1) domtbl girdilerinden
yeni olmali -- kural `stale_intermediate`ten alinir; (2) mtime'lardan
BAGIMSIZ sonuc sinamasi -- anotasyonun beta dedigi her komsu `beta`
siniflanmali, artakalan butce 60 (duzeltme oncesi 688). Birinci ayak bir
gun yanilirsa ikincisi hala tutar.

SILINDI: `beta_subunit_model.py`. Kendi beta HMM'ini kurmak icin yazilmisti
ve premisi (Pfam betalari kaciriyor) CURUDU. Yanlis varsayima dayanan bir
betigi depoda birakmak, birinin onu kosturmasina davetiyedir. Olcumun
kendisi bu maddede kayitli.

KAPANDI (60 protein) -- bkz. 135.24. Bu madde "ikisinden biri yanlis ve
hangisi oldugu olculmedi" diyordu. Olculdu: IKISI DE DEGIL. 60'in
tamaminda `neighbor_protein` satiri yok, cunku 58'i PSODOGEN ve 2'sinde
kayitta ceviri yok. Taranmamis bir dizi hakkinda hicbir model bir sey
soyleyemez; anotasyon beta genini DOGRU adlandiriyor, gen bozuk.

KAPANDI (grup 3) -- bkz. 135.26. Uc hipotezden hicbiri degil: GRUBUN
BIRIMI YANLISTI. %72,2 toplami, beta'si OLMAMASI DOGRU KIMYA olan
desaturaz ve hidroksilazlarca domine ediliyordu. cis-dihidroksilleyen
tiplerde yokluk %24,8.

## 123. Anotasyon kalitesi, operon iddialarinin DAYANDIGI bolgede olculmeli (acik)

Kullanicinin isaret ettigi sinifin dort ayri bicimi var ve her biri ayri
sinanir:

```
1. EKSIK GEN       DNA'da var, CDS anotasyonu yok   -> ORF cagrimi (madde 121)
2. YANLIS SINIR    CDS var, start/end kaymis        -> bosluk aritmetigini bozar
3. YANLIS URUN     metin yanlis                     -> bizi etkilemez (HMM kullaniyoruz)
4. BIZIM HATAMIZ   anotasyon dogru, siniflandirici  -> madde 122
                   kaciriyor
```

Madde 121 yalnizca TAMAMEN anotasyonsuz 47 kaydi kapsiyor. Asil olcum
sudur: onaylanmis her RO'nun +-10 kb penceresinde ORF cagrilip anotasyonla
karsilastirilsa, operon iddialarimizin dayandigi bolgede anotasyonun
EKSIKLIK ORANI olculur. Bu, "beta yok" ile "beta anote edilmemis" ayrimini
operon duzeyinde yapilabilir kilar.

MALIYET: 19.125 onayli RO x 20 kb ~ 380 Mb DNA; prodigal icin onemsiz. Pahali
olan kisim 17.074 gbk dosyasini yeniden okumak (cikarim ile ayni mertebe).

SART: cagrilan ORF'ler anote CDS'lerle AYNI sutuna karismaz
(`annotation_source` = `ncbi` | `prodigal`), cunku cagrilmis bir genin
sinirlari olculmus degil TAHMIN edilmistir.

## 124. IS/transpozon, duzenleyici ve promotor: ayni iki-hat yontemi (acik, 2026-11-04)

Kullanicinin onerisi: beta icin yapilani insersiyon dizisi/transpozon ve
duzenleyiciler icin de denemek, ve promotor belirlemeye calismak. Mevcut
durum olculdu ve oneri YERINDE -- ayni tek-kaynak zaafi orada hala duruyor.

### 124.1 IS / transpozon -- etiket su an REGEX ile konuyor

`control_elements.py` kendi yorumunda bunu zaten itiraf ediyor:

> "Transpozon" etiketi regex ile konuyor ve icine rekombinasyon makinesi de
> giriyor: urun adlarinda 117 kez "Holliday junction resolvase RuvX" var, ki
> o bir mobil element degil.

Yani karar veren sey URUN METNI. Beta olayinda ogrenilen tam buydu: metin
tek basina yer gercegi degil. IS ailesi "IS<n> family" kalibiyla
cikariliyor, yani ADLANDIRILMAMIS bir transpozaz ya da aile adi yazilmamis
bir kayit gorunmez kaliyor.

PLAN: Pfam modelleriyle (DDE_Tnp_1, DDE_Tnp_1_assoc, rve, HTH_Tnp_1,
Transposase_mut, DDE_Tnp_IS66 vb.) ikinci bir kanit hatti kurulur ve IKI
HAT KARSILASTIRILIR -- metin diyor model demiyor, model diyor metin demiyor,
ikisi de diyor. Uc sayi da yayimlanir; "duzelttik" demek yerine nerede
ayristiklari gosterilir.

SART (modulerlik): bu modeller `component_hmms/components.hmm` icine
KARISTIRILMAZ. O kutuphane RO bilesenleri icin ve `assign_component` oradan
karar veriyor; mobil element modelleri ayni kutuphaneye girerse o karar
belirsizlesir. Ayri kutuphane + ayri turetilmis tablo.

### 124.2 Duzenleyici -- konumla belirleniyor, baglanma olculmedi

`regulation.py` bunu da acikca yaziyor: duzenleyici "operonun 5' ucunun
yukarisindaki ilk gen" olarak aliniyor. `upstream_family` 11.422 girisin
8.061'inde BOS (operon_relations.py:281). Yani aile atamasi su an
girislerin ucte birinde var.

PLAN: HTH aileleri icin model tabanli atama (TetR_N, GntR, LysR +
LysR_substrate, MarR, AraC_binding, IclR) ve ayni iki-hat karsilastirmasi.
Konum kurali DEGISMEZ -- "yukaridaki ilk gen" tanimi aynen kalir, yalnizca
o genin AILESI metne degil modele dayanir.

### 124.3 Promotor -- su an motif aramasi HIC YOK

Mevcut "promotor bolgesi" yalnizca iraksak duzenleyici ile operon arasindaki
MESAFE (`intergenic_bp`) ve kaynak bunu zaten "TAHMINDIR; deneysel promoter
degil, inceleme icin ipucu" diye isaretliyor. Yani motif aramasi yok.

DURUST SINIR -- bu bastan yazilmali ki sonra asirilmasin. Bakteriyel sigma70
promotoru -10 (TATAAT) ve -35 (TTGACA) kutularindan olusur ve bu kutular
BILGI BAKIMINDAN ZAYIFTIR: rastgele bir intergenik bolgede de sik bulunur.
Bu yuzden operon BASINA "promotor burada" demek yayimlanabilir bir iddia
DEGILDIR ve oyle yayimlanmayacak.

YAYIMLANABILIR OLAN: ZENGINLESME. Gercek yukari-akis intergenik bolgelerde
motif siklig, EŞLENMIS bir null'a karsi olculur (ayni bolgelerin dizi
karistirilmasi + ayni genomlarin baska intergenik bolgeleri; GC icerigi ve
uzunluk eslenir, yoksa fark GC'den gelir). Depoda `null_model` makinesi
zaten var, yontem oraya oturuyor.

Ikinci yayimlanabilir sey: KARSILASTIRMA. Iraksak duzenleyicisi olan
operonlarin yukari-akisi ile olmayanlarin yukari-akisi karsilastirilir. Burada
null kendiliginden var ve iddia "sigma70 kutusu sudur" degil, "iki kume
arasinda fark var/yok" olur.

### 124.4 EN KRITIK SART (kullanici: "bunlari etkileyen istatistiklerin de
guncellenmesi kritik, burada iyi bi kontrol lazim")

Ayni gun bunun maliyeti olculdu: bilesen taramasi duzeltilince
veritabanindan turetilen 42 yayimlanmis dosyanin 42'si bir anda bayat hale
geldi ve sayfada bunu gosterecek hicbir sey yoktu.

Iki kapi kuruldu:

```
check_derived_statistics_are_fresh          veritabanindan turetilen her
                                            yayim dosyasi DB'den yeni olmali;
                                            dosyanin DB okuyup okumadigi
                                            uretici betigin KAYNAGINDAN cikarilir
check_published_files_are_in_rebuild_graph  yayimlanan her dosya rebuild
                                            grafiginin bir CIKTISI olmali;
                                            butce 25 ve KUCULMEK uzere
```

### 124.5 MODULERLIK (kullanici: "DB'nin modulerligini iyi planlamak cok kritik")

Olculdu: yayimlanan 56 analiz dosyasinin **25'i** hicbir `rebuild` adiminin
ciktisi degil. `rebuild.py --run` onlari hic tazeleyemez. 61 adimin 50'si
veritabanini girdi beyan ediyor, yani grafik SAHIP OLDUGU adimlar icin
dogru; sorun eksik dugumler.

Ikinci ayrisma: grafikte veritabaninin girdisi `genomic_context/` olarak
beyan edili, ama fiili cikarim ciktilari artik `/media/lin-bio/back2/`
altinda (`roar_gc_wide`, `roar_gc_a5`, `roar_gc_merged`). Yani BEYAN EDILEN
girdi ile GERCEK girdi ayrismis durumda; `genomic_context/` 2026-10-14'ten
ve genisletme ONCESI 26.693 adayi tutuyor.

TASARIM YONU (henuz uygulanmadi): bagimlilik dosya duzeyinde degil TABLO
duzeyinde beyan edilmeli. Her adim okudugu ve yazdigi TABLOLARI bildirsin,
her tablonun bir surumu (`db_revision.py` zaten surum tutuyor) olsun ve
turetilmis cikti okudugu tablolarin surumlerini kaydetsin. O zaman bayatlik
karari mtime'a degil TABLO SURUMUNE dayanir -- ve mtime bu depoda iki kez
yanilttigi icin bu fark onemli.

### 125. TAZELIK KAPISI KENDI KURALINI YAZIYORDU (2026-11-05/06, COZULDU)

`rebuild.py` 33/38'de `validation_gate` FAIL'i ile duruyordu; kapi 42
yayim dosyasini bayat gosteriyordu. SEBEP VERI DEGILDI, KAPININ KENDISIYDI
-- ve iki katmanda:

**Kural.** Kapi dosyanin mtime'ini `roar.sqlite`in mtime'iyla
karsilastiriyordu. `rebuild.py` o karsilastirmayi 2026-10-14'te BILEREK
kaldirmisti ("TABLO BAZINDA DB GIRDISI" notu): veritabanina zincirin ICINDE
yazan adimlar var, dolayisiyla DB mtime'i her kosuda en sona kayar ve
kendisinden ONCE kosmus her analiz bayat gorunur. Kanit: bayat denen
dosyalarin bazilari 0,3 saat "eski"ydi -- AYNI kosuda uretilmislerdi.
Duzeltme: planlama mantigi `rebuild.plan(state, propagate)` olarak tek
kaynaga cikarildi (`--plan` ciktisi refactor oncesiyle birebir dogrulandi)
ve kapi ona baglandi. `propagate=False`, cunku iki AYRI soru var:
"yeniden kosmali mi?" (yayilim dahil, 53) ve "yayimlanan sayi YANLIS mi?"
(yayilim haric, 20). Ileride bir ust adimin kosacak olmasi yayimlanmis
bir sayiyi yanlis yapmaz.

**Kapsam.** Kapi yargilanacak dosyayi uretici betigin KAYNAK METNINDE
"sqlite3.connect" / "--db" arayarak TAHMIN ediyordu. Grafige 8 yeni adim
baglaninca bunun yine adindan genis karar verdigi goruldu: dordu
veritabanini HIC okumuyor, girdileri bir DIS SERVIS
(`bridging_geometry_audit` PDB; `reaction_anchors`, `reaction_verification`
Rhea/ChEBI; `curation_claims` UniProt/NCBI). Boyle bir dosyanin "bugunku
tablo icerigiyle hesaplanmis olmasi" diye bir sey yoktur. Kapsam artik
GECISLI olarak grafikten cikariliyor: 69 adimin 55'i.

SONUC: 42 -> 10 -> 0. Gercekten bayat 10 dosya yenilendi (dokuz adim,
hicbiri basarisiz olmadi; `learning` 8.578 s, `bootstrap_confident`
4.389 s). Zincir sonradan 40/40 tamamlandi. Dogrulama: 171 kontrol,
166 PASS, 0 FAIL, 5 WARN.

### 126. YOL USTUNDE CIKAN BES GERCEK BOSLUK (COZULDU)

1. **`types/` bagimlilik grafigine GORUNMUYORDU.** `digest()` bir dizinde
   `open()` ile OSError alip None donuyordu, ve None ozetler
   karsilastirmada ATLANIYOR (`if not prev.get(p): continue`). Ustelik bir
   dizinin mtime'i icindeki dosya DUZENLENINCE degismez, yalnizca dosya
   eklenip silinince. Ikisi birlesince: veritabaninin TEK KAYIT KAYNAGI
   degistirildiginde HICBIR adim bayat gorunmuyordu. Kanit: state
   dosyasinda `"types": null`. `dir_digest()` eklendi (mtime'a guvenmez;
   32 MB altindaki dizinlerde ICERIK, ustunde ad+boyut); `types` 0,08 s.
   `oxygen_fate` ve `metabolic_role` "guncel" sayildiklari icin
   kendiliginden kosmayacaklardi -- ozeti IDDIA etmek yerine kosturularak
   uretildi.

2. **Yan etkili "yoklama".** `--cpu` destegi, betigi `--help` ile
   CALISTIRIP cikisinda aranarak bulunuyordu. argparse'i OLMAYAN bir
   betikte `--help` sessizce yok sayilir ve betik GERCEKTEN kosar: 60
   saniyelik bir "prob" dis servise istek atar, dosya yazar, cikisi
   atilir. `anchor_reactions.py` ve `reaction_balance.py` tam boyle.
   `script_flags()` artik kaynaktan okur.

3. **Ayak izi kapisinda "cozuldu" = "satir var".** Ayak izi tum replikonu
   kaplayan bir giris, SADECE `ro_cds_span`de kaydedildigi icin temiz
   sayiliyordu -- oysa kontrolun kendi FAIL metni zarari yaziyor.
   Kaydetmek kullanmak degildir. Ayrica "origin-spanning genes are
   recorded with their real parts" kontrolu `True` SABITIYLE cagriliyordu:
   adi kontrol, isi sayac; hic FAIL edemezdi. Yerine yanlislanabilir bir
   soru (parcalarin toplami == cds_length); A5'teki 10.474 kayitta
   sinandi, tutarsiz 0.

4. **`a5_match_guard` okudugu dizinleri beyan etmiyordu**
   (`roar_gc_a5`, `roar_gc_wide`). Beyan edildi.

5. **Atlas koken hatasi.** `oxygen_fate_computed.json`'u
   `reaction_balance.py` degil `compute_oxygen_fate.py` yaziyor.

### 127. MODULERLIK: 25 -> 6 (124.5'in UYGULANMASI)

11 dosya: zaten grafikte olan adimlarin BEYAN EDILMEMIS ikincil
ciktilariydi (regulation, variants, phylogeny, leaves, evidence_tiers,
stratified_stats, variant_signature). Dort adim HIC cikti beyan
etmiyordu. Beyan edilince gercek bir kenar da ortaya cikti:
`substrate_predictability`, `evidence_tiers`in urettigi
`reference_pairs.csv`ye bagli -- bu bagimlilik grafige gorunmuyordu.

8 dosya: yeni adim olarak eklendi. Ag sorgusu atanlar AGIR isaretli
(depo konvansiyonu: "AGIR cunku 737 ag sorgusu"), boylece `--run` dis
servisi KENDILIGINDEN sorgulamaz.

Bu boslugun ne sakladiginin KANITI: `carboxylate_method_card.txt`
30 Ekim'deki denetimden uretilip kalmis, oysa `carboxylate_audit.json`
31 Ekim'de yenilenmis. Yayimlanan kartta D/E vermeyen uye 2.388, gercekte
3.007. Dosya grafikte olmadigi icin kimse goremiyordu.

Kalan 6 dosya grafige ZORLANMADI, cunku hepsi ayni sebeple disarida
degil: 1 onbellek (`chebi_formula_cache.json`; `reaction_balance.py`
argv'den Rhea ID alan bir EL ARACI, pipeline adimi degil), 1 URETICISIZ
(`curation_claims_ncbi.json` -- depoda bu dosyayi yazan betik YOK, ve
atlas kokeni `verify_curation_claims.py (NCBI efetch pass)` diyerek
YANLIS iddia ediyor), 4 elle kosulmus tek seferlik analiz. Butce 6 ve
kapi uc turu ayri ayri yaziyor. Dogru olan bunlari tazelenebilir gibi
gostermek degil, NE OLDUKLARINI beyan etmek.

KAPANDI (2026-11-10): `freeze.py` artik grafikte (`site_freeze`, AGIR).
Yayimlanan sitenin kendisi de tazelik denetimine girdi; `--plan` artik
"canlidaki kopya eski" diyebiliyor. Olculen bedel neydi: gh-pages
2026-11-02'de kalmisti ve bilesen duzeltmesi oncesi operon sayilarini
gosteriyordu, YEREL HER DOSYA GUNCELKEN. Adim hafif zincirde kosmaz --
isi bayatligi gorunur kilmak, kendiliginden yayim yapmak degil.

### 128. A5 AYAK IZI OLCUMU (ETKISI YOK, KAYIT ICIN)

`ro_cds_span.py` `/media/lin-bio/back2/roar_a5.sqlite` uzerinde: 72.315
adayin 10.474'unde (%14,48) ayak izi aciklanamaz. Tek hata degil, IKI
AYRI POPULASYON:
  - 10.471 cok parcali (intronlu okaryot kayitlari; 2.484 bp CDS 667 kb'a
    yayilabiliyor) -- gercek, hata degil
  - 3 orijin asan, ayak izi tum replikon (CP060714.1, 5.167.372 bp)

Operona etkisi OLCULDU: bu 10.474 adayin **sifiri** dogrulanmis, yani
turetilmis tablolara hic girmiyorlar. Yayimlanan operon istatistiklerine
etkisi yok, yeniden cikarim gerekmiyor. A5 ingest'te dogrulanan cikarsa
kapi yakalar.

### 129. KOMSU GEN AILESI: URUN ADI DEGIL PFAM MODELI (124.1 + 124.2)

Komsu genin kategorisi `build_db.py` icinde urun adina uygulanan bir
duzenli ifadeyle konuyordu. Iki olculmus zaafi vardi ve ikisi de
kapatildi.

**Kurulan sey.** `fetch_family_hmms.py` Pfam modellerini AD ile degil
ERISIM ile indirir, inen dosyanin NAME/ACC satirini okur ve beklenen ad
tutmazsa HATA VERIR. Bu kapi hemen ise yaradi: PF01614 artik `IclR`
degil `IclR_C`, PF03811 artik `Zn_Tnp_IS1` degil `Zn_ribbon_InsA`. Ad
ile indirilseydi iki model sessizce yanlis olurdu.

UC AYRI kutuphane kuruldu ve hicbiri `component_hmms/components.hmm`
ile karistirilmadi -- o kutuphanedeki her model bir operon bileseni
etiketine baglidir, oraya baska soru sokmak taramanin anlamini degistirir:

  - `family_hmms/mobility.hmm` (12 model) -- transpozaz cekirdekleri
  - `family_hmms/recombinase.hmm` (7 model) -- rekombinaz/integraz/RuvX
  - `family_hmms/regulator.hmm` (9 model) -- TetR, GntR, LysR, MarR,
    AraC, IclR

Rekombinazin AYRI kutuphane olmasi isin can alici yeri: regex ikisini
ayni `transposon` kovasina koyuyor, ayirmak icin once ayri olculmeleri
gerekiyor.

**124.1 -- hareketli eleman. Olculen sonuc:** regex'in `transposon`
dedigi 3.385 komsu satirinin
  - 1.447'si transpozaz modeliyle DOGRULANDI (%42,7)
  - 486'si rekombinazdi, yani hareketli eleman DEGIL (%14,4)
  - 1.452'sinde hicbir model tutmadi (%42,9)
Regex'in kacirdigi transpozaz yalnizca 5 satir.

En buyuk aciklanamayan grup 725 proteinlik "recombinase-like
helix-turn-helix domain-containing protein". 28 modelin HICBIRI, esik
ALTI dahil, E<=1e-5'te bu proteinleri tutmuyor; `components.hmm` de
tutmuyor. Etiket tamamen urun metnine dayaniyor. ACIK: bu ailenin
hangi Pfam modeline karsilik geldigi bulunmadi.

**Kendi yanlisimi da olctum.** Ilk mobility kutuphanesinde `IstB_IS21`
(PF01695) vardi ve 67 kez `chromosomal replication initiator protein
DnaA`'yi IS21 ilan etti -- genel bir AAA ATPaz alani oldugu icin. Model
cikarildi; "regex'in kacirdigi transpozaz" 79'dan 5'e dustu. Model
tabanli yontem de ad tabanli yontemle AYNI hatayi yapabiliyor; fark,
hatanin olculebilir olmasi.

**124.2 -- duzenleyici ailesi.** Regex'in duzenleyici dedigi 32.233
komsu satirinda:
  - 20.088 model ile UYUSUYOR, yalnizca 2 catisma (Crp_Fnr -> MarR)
  - 2.860 satirda regex aile veremiyordu, model verdi
  - 4.038 satirda regex aile veriyor ama model vermiyor -- bunlar
    kutuphanede olmayan aileler (LuxR, sigma54, Crp_Fnr), eksiklik
    kutuphanede, veride degil
Operon duzeyinde: ailesi adlandirilmis RO 4.693'ten **5.549**'a cikti
(19.125 icinde, +%18).

**Nereye yazildi.** `gene_category`ye DEGIL. O tablo erken yazilir ve
neredeyse her adimin DB girdisi sayilir; oraya satir eklemek grafikte 63
adimi birden bayat gosterdi, oysa hicbirinin cikti sayisi degismiyordu.
Olculup geri alindi, yerine turetilmis tablolar kuruldu:
`neighbor_family_hit`, `neighbor_family`, `ro_upstream_family_hmm`.
`ro_regulation.upstream_family` ELLENMEDI -- sahibi
`analyze_regulation.py` ve kaynagi regex; iki yontem tek ada yazilirsa
sutunun ne anlattigi belirsizlesir. Birlesik okuma
`ro_upstream_family_hmm.source` ile hangi satirin hangi yontemden
geldigini acikca soyluyor.

Adim ZINCIRIN SONUNA kondu: su an ciktisini okuyan baska adim yok, sona
konunca kimseyi bayatlatmiyor. Bir tuketici eklenirse adim onun ONUNE
tasinmali.

Kapi: `check_family_scan_is_fresh` uc ayakli -- domtbl bayat mi
(`stale_intermediate`ten, kural kopyalanmiyor), manifestodaki her model
bir aileye bagli mi, tabloda `unassigned_model` kalmis mi. 174 kontrol,
169 PASS, 0 FAIL, 5 WARN.

ACIK: `control_elements.json`'daki transpozon sayilari HALA urun adi
okumasi. Model okumasi `family_scan.json`'da ayri yayimlaniyor ve
sayfadaki uyari bunu acikca soyluyor. Ikisini tek sayfada birlestirmek
icin `family_scan` adiminin `control_elements`in onune tasinmasi ve tum
hafif zincirin yeniden kosmasi gerekir.

### 130. 124.3'UN YARISI: DIVERGENT DUZENLEYICI DAHA KISA BOLGE PAYLASIYOR

LysR tipi duzenleyiciler klasik olarak hedef operonla kafa kafaya oturur
ve aradaki kisa bolgeyi paylasir. Sinanabilir hali su: yukarisinda
DUZENLEYICI olan divergent ciftlerin intergenik bolgesi, yukarisinda
baska bir gen olan divergent ciftlerinkinden kisa olmali.

**Olculen (promoter_architecture.py).** Dort ayri sinama, hepsi ayni
yonde:

| duzey | duzenleyici | digeri | fark | Cliff delta | p |
|---|---|---|---|---|---|
| giris (n=13.351) | 168 bp | 247 bp | -79 | -0,295 | 2,8e-180 |
| dizi (n=11.055) | 168 bp | 243 bp | -75 | -0,291 | 3,8e-145 |
| cins (n=2.345) | 162 bp | 221 bp | -59 | -0,270 | 2,1e-23 |

Cins ICI esli sinama: 53 cinsin **47'sinde** duzenleyici boslugu daha
kisa, farklarin medyani -83 bp, Wilcoxon p=2,3e-07. Yani sonuc "bazi
cinslerin genomu sikisik" kaymasindan gelmiyor; karsilastirma hep ayni
cinsin icinde kaliyor.

**Iddia edilmeyen sey.** Hicbir promotor DIZISI incelenmedi ve hicbir
promotor tahmin edilmedi. Veritabaninda komsu genlerin protein cevirileri
var, intergenik DNA yok. Kisa bolge paylasilan promotorla TUTARLIDIR,
onu GOSTERMEZ; burada transkripsiyon olcen hicbir sey yok.

Kapi `check_promoter_claim_is_levelled` dosyanin SEKLINI sinar, sonucunu
degil: uc duzey + cins ici sinama var mi, ve uc duzey ayni yonu mu
gosteriyor. Bir sonucun yonu veriyle degisebilir, disiplin degisemez.

ACIK: 124.3'un kalan yarisi -- uzunluk ve GC eslenmis bos modele karsi
MOTIF zenginlesmesi. Intergenik DNA'nin `gbk_files/`ten cikarilmasini
gerektiriyor; bu dizi hic cikarilmadi.

### 131. SIFIR UYELI 8 REFERANS: SEBEBI OLCULDU

Kapi aylardir "76 referansin 68'i uye topluyor" diye uyariyordu, sebebi
yazili degildi. Her birinin en yakin kardesi olculdu (blastp, refs71
kendi kendine):

| referans | en yakin kardes | kimlik | kardesin uyesi |
|---|---|---|---|
| 1_103_CarAa | 1_102_CARDO | %99,2 | 111 |
| 2_206_IpBAa | 2_203_cumA1 | %90,8 | 8 |
| 2_207_BpdC1 | 2_213_BphA1A2 | %79,2 | 21 |
| 2_210_XylC1 | 2_201_BphA1 | %67,0 | 1 |
| 3_309_NahAc | 3_312_NinA | %85,7 | 0 |
| 3_310_DntAc | 3_312_NinA | %93,1 | 0 |
| 3_312_NinA | 3_310_DntAc | %93,1 | 0 |
| 4_407_PbaA1 | 3_308_DxnA1 | %65,6 | 2 |

UC AYRI sebep var, tek bir "sorun" degil:
  - **Ad ikiligi (2):** CarAa/CARDO ayni substrat %99,2; cumA1/IpBAa
    SMILES'leri HARFI HARFINE ayni (`CC(C)c1ccccc1`) cunku kumen ile
    izopropilbenzen AYNI bilesik. Ikisinde de kardes tum uyeleri aliyor.
  - **Uc kisilik bos ada (3):** NahAc / DntAc / NinA birbirine %85-93
    benzer ve UCUNUN DE uyesi yok; en yakin uyeli referans 3_306_3NTDO,
    o da yalnizca 2 uye. Bu alt grup verinin tamaminda neredeyse yok.
  - **Gercekten uzak (3):** BpdC1, XylC1, PbaA1 -- en yakin kardesleri
    %65-79, yani uyeleri "calinmiyor"; bu enzimler veri kumesinde yok.

Yeni kapi `check_same_structure_same_name`: ayni SMILES farkli ADLA
yayimlanmamali. Ayni substrata birden cok referans olmasi normaldir
(bifenil 5, ftalat 3); hata olan ayni molekule sayfada iki ad vermektir.
Butce 3 ve ucu de yazili: cumen/izopropilbenzen GERCEK es anlamli cift,
oteki ikisi ("caffeine and methylxanthines", "polycyclic aromatic
hydrocarbons") KASITLI genis etiket -- orada ad degil SMILES dar.

KULLANICI KARARI: CARDO/CarAa ve cumA1/IpBAa birlestirilsin mi, yoksa
ikisi de ayri referans olarak mi kalsin? Birlestirme kume atamalarini
degistirir, o yuzden kendiliginden yapilmadi.

## 132. Tip sayfalarinda metin yuku (2026-11-10)

Kullanici: "enzimler icin olan sayfalarda cok yazi var ya onlari
duzeltemez miyiz". Olculdu: tip sayfalarinda tablolar disinda
**ortanca 1.804, en yuksek 2.210 gorunur kelime**; her bolum uzun bir
aciklama paragrafiyla basliyordu.

Hicbir cumle SILINMEDI -- bu projede yayimlanan her iddianin kosulu
yazili olmak zorunda. Degisen sey BICIM: **ilk satir cevabi verir,
yontem ve caveat bir tik oteye** (`details class="basis"`, ana sayfada
ayni sikayet icin zaten kullanilan bilesen). 21 blok donusturuldu:
bes ust kutu, uc figcaption, Key facts'in iki uzun notu, SDP'nin iki
paragrafi, en yakin tipler, carpisan uyeler, filogeni, habitat, uyeler.

Sonuc: **ortanca 1.100, en yuksek 1.435 gorunur kelime** (%39 azalma),
sayfa basina ortanca 12 acilir kutu. 80 tip sayfasinin hepsi 200
donuyor; 177 kontrolun sonucu degismedi (172 PASS, 0 FAIL, 5 WARN).

`longnote.js` iki yerde duzeltildi, cunku 45 kelime ustu her
`figcaption`/`p.note` blogunu kirpan o betik yeni acilir kutularla
cakisiyordu:
  - Sablon bir blogu kendisi katladiysa genel kirpma devre disi
    (yoksa IKI katlama: Figure 2'nin etiketli ozeti kirpmanin icinde
    kaliyordu).
  - Kirpma sinifi uygulandiktan sonra tasma gercekten var mi diye
    olculuyor ve olcum `document.fonts.ready` sonrasina alindi; yoksa
    hicbir sey gizlemeyen bir "Show the full note" dugmesi cikiyordu.

**Sayfalama GEREKMEDI, olculdu.** Varyantlar tablosu CntA'da 377 satir
ve ilk bakista sayfanin yarisi gibi duruyordu, ama `table.js` 24
satirdan uzun her tabloyu zaten `min(72vh, 680px)` kapakli kendi
kaydirma kutusuna aliyor: o tablo sayfada 377 satir degil ~680 px yer
kapliyor. Kelime saymak yaniltti, PIKSEL olcmek duzeltti.

**Asil eksik oradan cikti: ICINDEKILER SERIDI.** Tip sayfasi 1280 px
genislikte **9.414 px** ve 18 bolum, ama "On this page" seridi yalnizca
atlas sayfalarinda vardi -- cunku `pageContents()` kapisi `nav.subnav`
idi. Atlas icin yazilan gerekce ("sonuna kadar kaydirmadan sayfada ne
oldugunu goremiyorsun") tip sayfasinda bire bir gecerliydi.

`table.js` artik acik bir isarete de bakiyor: `[data-contents]` tasiyan
sarmalayicinin dogrudan h2 cocuklari seride girer. `cluster.html`'in
`div.numbered`'i isaretlendi. Numaralar CSS sayacindan geldigi icin
`textContent`'te yok; serit `.numbered` kapsaminda numarayi belge
sirasindan uretiyor, boylece serit ile basliklar ayni numarayi gosterir.
Regresyon: atlas seridi duruyor, anasayfa/about/varyant sayfalari
eskisi gibi seritsiz.

### 133. REFERANS FAZLALIGINDA KARAR: IKI CIFT AYRI KALIYOR (2026-11-10)

Kullanici karari: **"birlesmese de olur belirtelim ama"** -- yani
birlestirme yok, ama kararin KENDISI yazili olacak.

Durum zaten olculuyordu: `reference_redundancy.py` iki yakin-ayni cift
buluyor, `2_202_EdoA1 / 2_203_cumA1` %99,8 (etiketler etilbenzen ve
kumen, 47/8 uye) ve `1_102_CARDO / 1_103_CarAa` %99,2 (ikisi de
karbazol, 111/0 uye). Her iki tip sayfasinda ve kalite sayfasinda
"bolunme kismen rastlantisal" caveat'i vardi. **Eksik olan sey kararin
kendisiydi:** okuyucu fazlaligi gorup, gorulmemis mi sandi diye
dusunebilirdi.

**Karar olcumden AYRI dosyada duruyor.** `curation_decisions.csv`
(`pair_a, pair_b, decision, decided_on, rationale`). Sebebi: olcum her
kosuda yeniden uretiliyor, karar elle yazildi ve tarihi var; ikisini
ayni dosyada tutmak karari ilk `rebuild.py --run`'da silerdi. Dosya
`rebuild.py`'nin `reference_redundancy` adiminin girdisi
(`DECISIONS`), yani karar metni degisirse JSON bayatlar.

Kararin gezdigi UC yer -- "duzeltmeyi bulundugu yere yazip
yayimlandigi her yere tasimamak" bu projenin tekrar eden hatasi:
  1. tip sayfasi fazlalik caveat'i: "kept as separate references
     deliberately", gerekce `details.basis` icinde tarihiyle,
  2. BOS TIP caveat'i (`CarAa` 0 uye): zaten "isimlendirme karari"
     diyordu, artik karari da soyluyor,
  3. kalite sayfasi: tabloda "Decision" kolonu + her cift icin
     gerekceyi tasiyan `details.basis`.

**Yoldan cikan iki BAYAT metin.** Kalite sayfasinin 4. bolumu hala
birebir-ayni ciftleri anlatiyordu (`identical_groups` artik 0) ve
"76 records contain only 76 distinct sequences" diye okunuyordu;
kalemde ise uye bolunmesi **elle "29 to 6" yaziyordu**, oysa veri 47/8.
Ikisi de veriye baglandi: bolunme sayilari ve farkli-substratli cift
artik `refdup`'tan geliyor, birebir-ayni anlatimi kosullu.

**Kapi.** `validate_curation.py` yeni SERT kontrol: her yakin-ayni
ciftin tarihli ve gerekceli bir karar satiri var mi -- karari olmayan
cift, cifti olmayan karar ve gerekcesi bos satir hata. Fazlalik WARN'i
kalici oldugu icin kararin yazili olmasi da kalici sart. 178 kontrol,
173 PASS, 0 FAIL, 5 WARN.

### 134. A6 KURULDU: CAPRAZ REFERANSLAR ALAN OLDU (2026-11-11)

Dort tanimlayici (`ec`, `refseq`, `embl_nuc` + UniProt accession'i)
2026-11-01'den beri `ro_protein_status`'ta duruyordu ve sitede
**hicbir yerde** gorunmuyordu: okuyucu bir diziyi UniProt'ta, ENA'da ya
da RefSeq'te aramak istedigi anda elinde hicbir sey yoktu.

`cross_references.py` -> `ro_xref(seq_md5, db, identifier)` +
`analysis_out/cross_references.json`. Birim DIZI, accession DEGIL
(madde A1): tek bir dizi YEDI UniProt accession'i altinda gecebiliyor
(NC_003143.1:2801274-2802399:-1 boyle).

**EC'yi oteki ucuyle ayni kefeye KOYMADIK.** Olculdu: EC tasiyan
3.813 dizinin yalnizca **87'si** incelenmis (Swiss-Prot) bir kayittan
geliyor, **3.726'si** yalnizca incelenmemis (TrEMBL) kayittan --
yani %97,7'si otomatik anotasyonun yaydigi bir TAHMIN. Bir tip
sayfasina ciplak "EC 1.14.15.7" yazmak, kuratorlu substrat etiketiyle
ayni gorunen 3.726 islev iddiasi uretirdi. O yuzden EC her giris
sayfasinda inceleme durumu cipiyle ve "bu bir adres degil, bir islev
iddiasi" diyen acilir kutuyla duruyor.

**PAYDA BEYANLI (madde A3).** Butun tanimlayicilar UniProt kaynakli
dizilerden geliyor, o yuzden payda 68.747 (UniProt kaynakli), 73.299
(tum birlesim) DEGIL. Tum birlesim uzerinden alinsa oran yanlis seyi
olcerdi: EC'nin yoklugunu degil, UniProt kaydinin yoklugunu.

**Yeni olcum: dis anotasyon bizim tipimizle anlasiyor mu.**
`/atlas/evidence` bolum 8, tip basina EC okumasi / farkli EC / en sik
EC'nin payi / 1.14 sinifi disindaki EC'ler. CntA'da en sik EC okumalarin
%92'si, KshA15'te %50 -- yani bir kuratorlu tip, UniProt'un birkac
farkli reaksiyonla etiketledigi dizileri topluyor. Bolum 3'un dersinin
(kimlik substrati belirlemez) BAGIMSIZ bir yonden gelisi. 1.14 disi
EC'ler "hata" degil "listelenmis": UniProt kaydi cok alanli olabilir ve
EC kaydin tamamina aittir.

**Iki hatayi olcum yakaladi, ikisi de bicim kontrolunun ATILAN SAYIYI
raporlamasi sayesinde:**
  1. EMBL deseni `\d{5,8}` idi ve 6 harf + 9 haneli GECERLI WGS
     accession'larini (`JACBZF010000002`) bozuk sayip **5.358 kaydi
     atiyordu**. "Kac tane atildi" yazilmasaydi gorunmezdi.
  2. Adim zincirde `stats`'tan ONCE duruyordu; `ro_xref` tablosunu
     yazdigi icin kendisinden sonraki her DB girdili adim (stats,
     family_scan, a5_match_guard, carboxylate_anchor_holdout) onu
     "degismis ust girdi" sayip bayatliyordu -- ucu AGIR. Adim YAPRAK,
     o yuzden zincirin sonuna alindi. **Ders: DB'ye tablo yazan yeni bir
     adimin zincirdeki YERI, kendisinden sonraki her seyi bayatlatir.**

Dogrulayici: sekiz yeni kontrol (sahipsiz satir, bos tanimlayici, URL
kalibi olmayan ad, EC'nin dort alani, RefSeq oneki, PAYDA TAZELIGI,
atilan kayit, EC inceleme durumunun yayimlanmis olmasi).
`ro_xref`'in uc sutunu `IDENTIFIER_TEXT_FIELDS`'a eklendi -- yoksa dil
kontrolu onlari hic gormez.

**AYNI HATAYI BIR KEZ DAHA YAPTIM, ve ikinci bir tanesini yakaladi.**
Blogu ilk olarak yalnizca `ro.html`'e yazdim. Ama YAYIMLANAN sitede giris
sayfalari o sablonla cizilmiyor: tek kabuk (`entry.html`) + `data/entry/`
altinda 256 parca JSON + `static/entry.js`. Yani blok 19.109 giris
sayfasinin HICBIRINDE gorunmeyecekti; yalnizca gelistirme uygulamasinda.
Duzeltme: `build_entry_data.py` artik parcaya `xr` (tek harfli anahtarlar:
u/r/n/e) ve `xrv` (EC'yi incelenmis kayit tasiyor mu) yaziyor, `entry.js`
de ayni blogu ciziyor. Bedel olculdu: ham JSON +0,83 MB (46 MB uzerine
%1,8). **17.535 / 19.109 giris (%91,8) en az bir capraz referans tasiyor.**

Bunu ararken AYNI SINIFTAN ikinci bir hata cikti: `entry.js:207`
`core_identity`'yi hala **"alignment identity to nearest reference"** diye
etiketliyordu. O etiketin yanlis oldugu 5f68cc4751'de ("yayimlanan sayi ile
basliginin uyusmadigi 9 yer") tespit edilip `ro.html`'de duzeltilmis, ama
`entry.js`'e tasinmamis. Sayi kumenin KENDI alt aile uzlasilarina karsi
olculuyor, hicbir kuratorlu referansa dokunmuyor -- fark CntA'da 0,99'a
karsi %24,4. Yani duzeltilmis hali gelistirme uygulamasinda, yanlis hali
canlidaki 19.109 sayfanin hepsinde duruyordu. Duzeltildi.

ACIK KALAN (A6'nin oteki yarisi): InterPro, KEGG, MetaCyc. Bunlar
elimizde YOK; UniProt REST'ten 68.747 accession icin cekilmesi gerekiyor.

## 135. "Biyoinformatik acidan profesyonel" -- olculmus eksik listesi (2026-11-11)

Kullanici istegi: "bu DATAbase'i profesyonel bi hale getirelim biyoinformatik
acidan". Tahmin yerine OLCUM yapildi; asagidaki her satir kosturularak
dogrulandi. Siralama "yayimlanan bir iddiaya dokunuyor mu" olcutune gore.

| # | eksik | olculmus durum | yayimlanan iddiaya dokunur mu |
|---|---|---|---|
| 135.1 | NCBI taxid + rank'li taksonomi | **YOKTU**; cins 5 ayri dosyada `organism.split()[0]` ile tahmin ediliyordu | **EVET** -- cins, her istatistigin KARAR birimi |
| 135.2 | Surum + atif (CITATION.cff, DOI, makine-okunur meta) | yok; `db_revision.py` ici surum tutuyor ama site surum YAYINLAMIYOR | hayir, ama DB atif edilemez |
| 135.3 | Substratlar icin ChEBI | yok; `chemistry.csv` yalnizca Ingilizce ad + SMILES | kismen -- substrat adi dis dunyaya baglanamiyor |
| 135.4 | Genomik baglam icin GFF3 ciktisi | yok; operon/komsuluk yalnizca DB icinde | hayir, ama veri yeniden kullanilamaz |
| 135.5 | Pfam/InterPro surum sabitleme | accession'lar var, **Pfam SURUMU yazili degil** | evet -- bir HMM vurusu adlandirilmis surum olmadan yeniden uretilemez |

Iyi durumda OLAN ve degistirilmemesi gerekenler (bakildi, dogrulandi):
`candidate_id` koordinat tabanli (`NC_000853.1:1311057-1312053:1`) -- yani
bir tip yeniden siniflandirildiginda kalici baglantilar BOZULMUYOR; LICENSE
(MIT) var; `db_revision.py` surumleri degismez kayitlar halinde tutuyor;
FASTA indirmeleri ve `/download/analysis/` saglamligi `check_site.py` ile
sinaniyor.

### 135.1 Taksonomi KONTROLLU SOZLUK oldu -- ve cins tahmini curutuldu

Bu DB'nin kendi istatistik kurali "karari CINS DUZEYI verir" diyor. Ama
cins hicbir yerde TANIMLI degildi; **ON DORT** yerde ayri ayri dizgiden
tahmin ediliyordu. Dokuzunu grep buldu (`analyze.py`, `analyze_ecology.py`,
`geography.py`, `threshold_sensitivity.py`, `build_search_index.py`,
`promoter_architecture.py`, `make_explorer.py` (iki yer), `make_report.py`,
`characterize_leaves.py`); kalan besini **dogrulayici muhafizin kendisi**
yakaladi: `ancestry.py`, `isolation_source.py` (tur icin), ve
`webapp/app.py` + `webapp/atlas.py` (iki yer) -- ki bu son ucu dogrudan
YAYIMLANAN sayfalari besliyordu. Ustelik ikisi eksik
organizmada BIRBIRINE ZIT: `analyze.py` "Unknown" deyip onu bir cins olarak
SAYIYOR, `analyze_ecology.py` None deyip satiri DUSURUYOR -- yani ayni
kumenin cins paydasi hangi betigi kosturdugunuza gore degisiyordu.

Olculen hatalar:

1. **`Candidatus` bir cins degil.** `split()[0]` onu cins sanip **17 ayri
   gercek cinsi** (Micropelagos, Rhodobacter, Pelagibacter,
   Puniceispirillum, Pseudothioglobus, ...) tek kovaya yiginiyordu: 57
   giris. Yonu de en kotusu -- ayri cinsleri tek cins gostermek, bir
   sinyalin tek cinste YOGUNLASTIGI izlenimini guclendirir, ki cins
   muhafizinin varlik sebebi tam olarak bunu engellemek. `uncultured` ayni
   sekilde bir cins sayiliyordu.
2. **Dizgi, cins adi uretiyor.** `'Nostoc azollae' 0708` -> cins `'Nostoc`
   (tirnakla), `[Enterobacter] lignolyticus` -> `[Enterobacter]`
   (koseli parantezle; NCBI'nin "yanlis uygulanmis ad" isareti).
3. **Yeniden adlandirmalar kaciriliyor.** Olculdu: 325 dogrulanmis girisin
   cinsi rank'la bakildiginda degisiyor -- Sinorhizobium -> **Ensifer**,
   Rhizobium rhizogenes -> **Martinezella**, [Enterobacter] lignolyticus ->
   **Pluralibacter**.

Akla gelen ilk duzeltme -- `replicon.taxonomy` soyagacinin SON alanini cins
saymak -- de curutuldu: 17.074 replikonun **2.006'sinda (%11,7)** o alan
cins degil, bir TUR KOMPLEKSI ("Mycobacterium tuberculosis complex",
"Bacillus cereus group"). Soyagacinda rank etiketi yok; alan sirasi rank
anlamina gelmiyor. O yuzden cins NCBI `nodes.dmp`'ten RANK ile okunuyor.

**En agir bulgu: DB'nin 7. en kalabalik "cinsi" diye bir cins yok.**
`Mycolicibacterium` (612 giris) bugunun NCBI taksonomisinde cins DEGIL,
**altcins**: `Mycobacterium (subgen. Mycolicibacterium)`, ve o yolda cins
rank'li dugum HIC YOK (tur -> altcins -> AILE Mycobacteriaceae). Yani
yayimlanan bir cins duzeyi sinyal, var olmayan bir cins uzerine kuruluydu.
Bu satirlar DUSURULMEDI (payda sessizce kuculurdu) ve cins UYDURULMADI;
`genus_state='unassigned'` ile tasiniyor ve en yakin adlandirilmis rank
VEKIL olarak yaziliyor. Cins duzeyi bir test ya vekille gruplar ya da
disarida birakir -- hangisi olursa olsun **614 girisi (%3,2) beyan etmek
zorunda**.

Kurulan (yerel, aga yalnizca taxdump icin cikildi, md5 dogrulandi):
`replicon_taxon(nucleotide_id, taxid, taxid_current, status, name, rank)`
-- 17.073/17.074 replikon, 6.223 tekil taxid, 71'i `merged.dmp` ile
guncele tasindi; `taxon_lineage(taxid, rank, ancestor_taxid, ancestor_name)`
-- 43.526 satir, 1.380 rank'li tekil cins; `taxon_genus(...)` tek sahipli
cevap. taxid kaynagi her GenBank kaydinin `source` ozelligindeki
`/db_xref="taxon:N"` -- 17.074/17.074 dosyada VAR, hic aga cikilmadi.

Yan bulgu: `replicon` tablosunda **bos bir hayalet satir** var
(`nucleotide_id=''`, `status='no_annotation'`, `cds_count=0`); hicbir `ro`
satiri ona bagli degil ama replikon sayisini 1 sisiriyor, yani replikon
basina her yuzdenin paydasi bir fazla. (A3 payda muhafizi kapsaminda.)

NCBI'nin kendi tutarsizligi da not: `Mycolicibacter` icin ayni bicimdeki
dugum (`Mycobacterium (subgen. Mycolicibacter)`, taxid 1073531) rank'i
**genus**, `Mycolicibacterium` icin ayni bicimdeki dugum rank'i
**subgenus**. Rank kaynaktan oldugu gibi aktariliyor, duzeltilmiyor.

ACIK: (a) bes dizgi tahmincisini tek sahipli okuyucuya cevirmek ve
asagi akisi yeniden kosturmak, (b) dogrulayici kontrolleri, (c) cins
duzeyi iddia yayimlayan HER sayfaya duzeltmeyi tasimak.


### 135.2 Surum ve atif YAYIMLANDI (2026-11-11)

Sitede ne yazar adi ne atif bicimi vardi; bir okuyucu buradan bir sayi
alip makalesine koymak istedigi anda neyi atif edecegini bilmiyordu.
Daha sinsi olani ikinci bulgu: **veri surumu hic yayinlanmiyordu.**
`db_revision.py` surumleri degismez kayitlar halinde tutuyor, ama son
kayit `rev-0007` (2026-10-30) iken DB on iki gun boyunca degisti --
yani site, gosterdigi veriyi tanimlamayan bir durumdan yayinlaniyordu.

Yapilan:
- `CITATION.cff` (CFF 1.2.0). DOI ve ORCID **uydurulmadi**, acik madde
  olarak yazildi; DOI alinirsa `preferred-citation` eklenecek.
- `/about#cite`: neyi atif etmeli, hangi surumle, ve neden. Giris
  kimliklerinin koordinat tabanli oldugu -- yani yeniden siniflandirma
  kalici baglantiyi bozmadigi -- burada soyleniyor.
- Footer'da **veri surumu, kod surumunden AYRI** yaziliyor. Ikisini tek
  sayiya indirmek yanlis olurdu: ayni commit, yeniden kosturulmus bir
  pipeline ile baska sayilar yayinlayabilir.
- `base.html` basinda schema.org `Dataset` JSON-LD (`version` = VERI
  surumu). DB artik yalnizca baglantiyi zaten bilen biri tarafindan degil,
  veri kumesi arama motorlariyla da bulunabilir.
- Dogrulayici: `CITATION.cff` alanlari + **yayimlanan veri surumu DB ile
  ayni durumu anlatiyor mu** (rev kaydindaki giris sayisi ile canli sayim
  karsilastirilir). rev-0007'nin on iki gun bayat kalmasi bir daha sessizce
  olmayacak.

ACIK (kullanici karari): Zenodo DOI ve ORCID.

## Yeniden uretim

```
bash run_all.sh 16                  # tum pipeline (uzun)
cd webapp && uvicorn app:app        # yerel sunucu
python3 freeze.py --out site --base /pro-sim-blast && bash deploy_pages.sh \
  git@github-prosimblast:recepcanaltinbag/pro-sim-blast.git
```

### 135.1c — Muhafızın bulduğu 18 yer: dizgi tahmini 11 dosyada daha duruyordu

135.1'de cins tanımı düzeltildi ve 16 çağrı yeri `ro_taxonomy`'ye bağlandı.
Bu yeterli sanıldı. Kapı (`validation_gate`) o koşuda **BAŞARISIZ** döndü ve
geri sürünme muhafızı **18 yer daha** saydı:

| dosya | yer | ne yapıyordu |
|---|---|---|
| `stats_overview.py` | 1 | projenin **istatistik motoru**; "kararı cins düzeyi verir" kuralının uygulandığı yer |
| `carboxylate_audit.py` | 1 | 19.125 girişin denetim kaydındaki `genus` alanı |
| `operon_relations.py` | 2 | yükleyici + diverjans bloğu |
| `variant_and_operon.py` | 2 | cins mi tip mi sorusu (hem cins hem tür) |
| `control_elements.py` | 2 | cins + tür ("ilk iki kelime") |
| `discover_subfamilies.py` | 2 | alt-aile başına cins dökümü |
| `recursive_homogenize.py` | 1 | yaprak başına cins dökümü |
| `ferredoxin_residue.py` | 1 | Fd bağlamı |
| `redundancy.py` | 1 | tekrarlanan dizilerin cins dağılımı |
| `learn_from_data.py` | 1 | cins bir **öğrenme özelliği** olarak |
| `promoter_architecture.py` | 1 | **kasıtlı**: birim duyarlılığı ölçümü |

Buradaki ders 135.1'in dersinin tekrarı değil, devamı: *bir ölçümün ikinci
uygulaması* ilk taramada bulunmaz, çünkü tarama "aradığını bildiğin şeyi"
bulur. `stats_overview.py`'nin satırı en pahalısıydı — cins muhafızının
koruduğu her iddia, muhafızın kendisi dizgiden hesaplandığı için ölçülmemişti.

Dönüşüm iki biçimde yapıldı. Replikonu sorgulayan betikler
`ro_taxonomy.genus_map` kullanır. Replikonu hiç görmeyen betikler (alt-aile
keşfi, yaprak homojenleştirme, özellik öğrenme) için **aynı tablonun aday
kimliği anahtarlı görünümü** eklendi: `genus_by_candidate` /
`species_by_candidate`. Alternatif, o betiklerin SQL'ine `nucleotide_id`
eklemekti; o, veri akışını değiştirir, bu değiştirmez.

Kasıtlı tek istisna **satıra** bağlandı, dosyaya değil: muhafız bir önceki
satırda `genus-string-intentional` yorumunu arar. Dosya bazlı muafiyet, o
dosyaya **sonradan** eklenen gerçek bir ihlali de muaf tutardı.

**Düzeltmenin ölçülmüş etkisi** (19.125 doğrulanmış giriş, aynı veri, iki tanım):

| sıra | eski (dizgi tahmini) | yeni (NCBI rank'ı) |
|---|---|---|
| 1 | Pseudomonas 1955 | Pseudomonas **1944** |
| 2 | Burkholderia 1289 | Burkholderia **1294** |
| 3 | Mycobacterium 1225 | Mycobacterium **1260** |
| 6 | Rhodococcus 618 | Rhodococcus **636** |
| 7 | **Mycolicibacterium 612** | — (altcins; vekil olarak taşınır, cins sayılmaz) |
| 8 | Novosphingobium 521 | Novosphingobium 521, **Paraburkholderia 492** |

942 ayrı "kova" → **948 rank'lı cins**. Yani düzeltme yalnızca bir adı
silmedi: tabloyu hem küçülttü (uydurma kova gitti) hem büyüttü (NCBI'nin
ayırdığı cinsler ayrıldı). Sıralamanın ilk üçü bile değişti — dizgi tahmini
Pseudomonas'a 11 giriş fazla, Mycobacterium'a 35 giriş eksik veriyordu.

### 135.1d — Kapı grafiğin sonuna taşındı; ve donma yolu canlı yoldan ayrı

İki yapısal hata, ikisi de aynı sebeple görünmüyordu: **doğrulandı sanılan
şey, yayımlanan şey değildi.**

**Kapı sıralaması.** `validation_gate` grafikte `stats`'tan hemen sonra
duruyordu, ama ondan **sonra** 15'ten fazla adım yayımlanan dosya üretiyor —
ve kapının tazelik kontrolü yayımlanan **tüm** kümeyi yargılıyor. Yani kapı,
kendisinden sonraki bir adım bayatken asla yeşil dönemiyordu; tam koşuda
bile, çünkü o adımlar henüz koşmamıştı. Pratikte bu, `--only` ile kapıyı da
içeren bir zincirin **kilitlenmesi** olarak çıktı: kapı başarısız olup zinciri
durduruyor, durduğu için de kapının şikâyet ettiği adımlar hiç koşmuyor.
Kapı `site_freeze`'den hemen önceye taşındı (çıktı listesi boş, hiçbir adım
`provenance.json`'u girdi olarak okumuyor, yani taşımak güvenli).

**Donma yolu ≠ canlı yol.** `build_entry_data.py` `ro.get("nucleotide_id")`
çağırıyordu; `ro` bir `sqlite3.Row` ve `.get()` yok. Giriş shard'larının
taksonomi bloğu bu yüzden hiç yazılmamıştı. Önceki turda "giriş sayfaları
NCBI taxon satırını gösteriyor" diye doğrulanan şey **Flask yoluydu**
(`app.py` rotası) — yayımlanan site ise `entry.html` + shard + `entry.js`
üçlüsünden render ediliyor, ayrı kod. Ders: yayımlanan bir alanı doğrulamak,
**yayımlanan artefakttan** saymaktır; canlı sunucudan bakmak başka bir kod
yolunu doğrular. Shard'lar şimdi tek tek sayıldı: 19.125 taxid, 18.511
rank'lı cins, **614 vekil rank** — beyan edilen sayıyla birebir.

Yayım: `rev-0008`, 204 kontrol / 199 PASS / 0 FAIL / 5 WARN, `check_site.py`
temiz, 3.483 dosya gh-pages'e gitti, canlıda doğrulandı.

### 135.5 — HMM modelleri sabitlendi: iki kütüphane ayrı Pfam release'inden

**Soru:** yayımlanmış bir HMM isabeti tekrar üretilebilir mi? Bir domain
çağrısı ("bu komşu bir ferredoksin") modelin kendisine bağlı, ve Pfam
modelleri sürüm atlıyor. "Pfam modeli kullanıldı" cümlesi bir yöntem
tarifi değil.

Ölçüm (`pin_hmm_libraries.py --check-upstream`, 2026-11-12, 37 Pfam modeli):

| kütüphane | model | geride | besledigi şey |
|---|---|---|---|
| `component_hmms/components.hmm` | 9 | **8** | `has_beta`, `has_ferredoxin`, `has_reductase` — operon türetimi |
| `family_hmms/mobility.hmm` | 12 | 0 | transpozaz / IS aileleri |
| `family_hmms/recombinase.hmm` | 7 | 0 | rekombinaz / integraz |
| `family_hmms/regulator.hmm` | 9 | 0 | düzenleyici aileleri |

Yani iki kütüphane **ayrı release'lerden** gelmiş: aile kütüphaneleri
güncel, bileşen kütüphanesi bir sürüm geride (`PF00355.33` → upstream
`PF00355.34`; 9 modelin yalnızca `PF08022.19`'u o release'te bump
olmamış). Geride olan kütüphane, yayımlanan en çok sayıyı besleyen.

**Geride olmak hata değil; bilmemek hata.** Yeniden indirmek yayımlanmış
her HMM isabetini değiştirir, yani operon türetiminin tamamının yeniden
hesabı ve yeni bir veri sürümü demektir — bu bir sürüm içinde sessizce
yapılacak bir şey değil. O yüzden yapılan şey **sabitlemek**:

- Alıntı birimi `PF00355.33` gibi **accession.sürüm** çifti + model
  bloğunun özeti. Pfam **release etiketi** HMM kaydının içinde yok, bu
  yüzden uydurulmuyor; accession sürümü release'ten daha kesin zaten
  (bir release içinde modeller ayrı ayrı sürüm atlar).
- `ROmotif71` ve `RieskeDB71` Pfam'dan gelmiyor; accession sürümü
  olmayan 77 model **ayrı listede**, yalnızca özetle sabitli — sürümsüz
  bir modeli sürümlü göstermek ölçülmemiş bir güvence olurdu.
- Kapı kontrolü **ağa çıkmaz**: diskteki dosyanın, sabitlemenin anlattığı
  dosya olup olmadığına bakar. Bir kapı dış bir servisin erişilebilirliğine
  bağlı olamaz. Sürüklenme ölçümü `--check-upstream` ile ayrıca koşar ve
  tarihiyle rapora yazılır.

Okuyucu doğrulaması: kütüphane özetim `fetch_family_hmms.py`'nin kendi
manifestosuyla birebir aynı çıktı (`regulator` → `a638d9b4b0d098cf`).

### 135.4 — Genomik bağlam GFF3 olarak: standart biçimde yayımlanmayan veri başka bir analize girmez

Bu DB operonları, komşu genleri ve bileşen çağrılarını kendi tablolarında
tutuyordu; okumak için SQL ya da bizim CSV'lerimiz gerekiyordu. GFF3 ise
IGV, JBrowse, Artemis, `bedtools` ve `gffread`'in doğrudan okuduğu biçim.

**En kolay yanlış yapılacak yer koordinat.** Veritabanı BioPython'un
konvansiyonunu tutuyor: `start` 0-tabanlı, `end` dışlayıcı. GFF3 ise
1-tabanlı ve kapsayıcı. Varsaymak yerine ölçtüm — 4.000 onaylanmış
girişin **3.996**'sında `end - start == 3 × (aminoasit + 1)`, yani
yarı-açık aralık ve dur kodonu dahil. Dışa aktarım `start + 1` yazıyor.
Kalan 21 kayıt (tam sette) bu eşitliğe uymuyor; **sessizce düzeltilmiyor**,
`partial_coordinates=true` ile işaretleniyor — bir koordinatı "düzeltmek"
onu uydurmak olurdu.

**Ne yazılmıyor:** çağrılan operona girmeyen komşu genler. Arama
penceresinin içinde ama operonun dışında kalan bir gen bir *bağlam*
değil, bir pencere artefaktı; yazılsa dosya 117 bin satırdan bir milyona
çıkardı.

Üretilen: 11.610 `##sequence-region`, 19.125 `operon`, 86.300 `CDS`
(19.125 alfa + neighbor'a bağlanan 67.175 operon geni), 105.425 özellik.

**Yazan taraf kendi çıktısını ayrıştırıyor.** Bir dosyanın GFF3 olduğunu
söylemek onu GFF3 yapmıyor, ve onu okuyacak araç hatayı bizim değil
kullanıcının ekranında verir. `export_gff3.py` yazdıktan sonra dosyayı
yeniden okuyup 9 kolon / 1-tabanlı aralık / çözülen her `Parent` /
`##sequence-region` sınırlarını sınıyor; ihlal varsa dosya geçersiz
işaretlenip çıkış kodu 1 dönüyor. Ölçüm: 86.300 Parent bağlantısının
tamamı çözülüyor, 0 ihlal.

Kapı **iki ayrı şeyi** sınıyor, çünkü biri ötekini gizler: biçim bozulması
(kolon, Parent, koordinat tabanı) *ve* bayatlık — geçerli bir GFF3,
tabloların eski içeriğiyle hesaplanmış sayıları yayımlayabilir ve onu
okuyan araç hiç şikâyet etmez. Kapı dosyadaki operon/CDS/replikon
sayılarını DB'den yeniden sayıp karşılaştırıyor, betiğin kendi raporuna
güvenmiyor.

Yan düzeltme: `/download/analysis/<ad>` rotası `open(...).read()` ile
**metin modunda** okuyordu; bir `.gz` indirmesi `UnicodeDecodeError` ile
500 verirdi. İçerik artık her zaman ikili okunuyor, tür uzantıdan
seçiliyor.

### 135.3 — Substrat için kontrollü kimlik: ChEBI, ve kimlik ile adayın ayrılması

`chemistry.csv`'de substrat bir dizgiydi ("naphthalene"). Dizgi hiçbir dış
kaynağa bağlanmaz. Şema kuralımız gereği yeni alan `chemistry.csv`'ye
**sütun olarak eklenmedi**: türetilmiş tablo (`type_substrate_chebi`) +
tazelik kontrolü.

İki yol var ve **eşit değiller** — `state` bu ayrımı taşıyor ve sayfa
görsel olarak ayırıyor (dolu kenarlık kimlik, kesikli kenarlık aday;
ayrım renge binmiyor, çünkü renk körlüğünde ve basımda renk kaybolur):

| durum | ne | sayı |
|---|---|---|
| `rhea_equation` | tipin kendi UniProt→Rhea kaydındaki **dengelenmiş denklemden** okunan ChEBI | **39** çift |
| `name_match` | kuratörün yazdığı adın ChEBI'de **tam** eşleşmesi — bir **aday**, kimlik değil | 28 çift |
| çözülmeyen | 20 tip (%26) | — |

Bir tipin birden çok substratı olabilir ve bu bir hata değil: NdmA'nın Rhea
kayıtları kafein, teofilin ve 1,7-dimetilksantin üzerinde ayrı ayrı
dengelenmiş denklemler veriyor. Tablonun birimi bu yüzden (tip × substrat).
Bir ad birden çok ChEBI kaydına **tam** eşleşiyorsa hiçbiri seçilmiyor —
seçmek, iki molekülden birini rastgele doğru saymak olurdu.

**Yan bulgu, ve bu maddenin asıl değeri.** Kuratörün yazdığı substrat adı
11 tipte dengelenmiş denklemin substrat adlarıyla aynı değil. Ama bir **ad
farkı bir uyuşmazlık değildir**: "dicamba" ile
"3,6-dichloro-2-methoxybenzoate" aynı molekül. Adları karşılaştırmak kurt
var demek olurdu. Karar **yapıdan** verildi — kuratörün SMILES'inin ağır
atom sayımı, Rhea katılımcısının formülüne karşı:

| hüküm | sayı | ne demek |
|---|---|---|
| `same_molecule` | 8 | adlandırma farkı (trivial ↔ sistematik ad, tuz, stereo) |
| `generic_rhea_class` | 1 | Rhea jenerik formül yazıyor (`C7H4O2X`, X halojen yerine); atom sayımı karar veremez |
| `different_molecule` | **2** | kuratörün yapısı denklemin hiçbir substratıyla aynı değil |

İki gerçek uyuşmazlık: `1_114_NdmB` kurasyonda **caffeine** (C₈N₄O₂) ama
Rhea denklemleri **theobromine / 3-methylxanthine** üzerinde; `2_208_TecA1`
kurasyonda **1,2,4,5-tetraklorobenzen** (C₆Cl₄) ama denklem
**klorobenzen** (C₆Cl). Bu tam olarak `reaction_balance.py`'nin var olma
sebebi olan hata sınıfı: doğru aritmetiğin **yanlış moleküle** uygulanması.
Kapıda **WARN**, FAIL değil — bir kapı kuratörün yerine karar veremez, ama
sessiz de kalamaz. Hiçbir kurasyon bu turda değiştirilmedi.

**Atom sayacı doğrulandı:** kuratörün adının Rhea adıyla aynı olduğu 15
çiftte sayım Rhea formülüyle **15/15** birebir tutuyor, 0 uyuşmazlık. Yani
çıkan bir uyuşmazlık sayacın değil verinin uyuşmazlığı. (Hidrojen
sayılmıyor: örtük H'yi doğru saymak valans kuralları ister ve C/N/O/S/halojen
iki molekülü ayırt etmeye yetiyor. `rdkit` kurulu değil.)

Ağ, kapının koşulu olamaz: `--lookup-names` yalnızca önbellekte olmayan
adlar için ağa çıkar, boru hattı adımı **önbellekten** aynı 67 çifti
çevrimdışı üretiyor (doğrulandı).

### 135.6 — İki kurasyon uyuşmazlığı literatürden kapatıldı; ve biri hata değildi

135.3'ün yapısal ayırt edicisi iki tipte "ayrı molekül" hükmü verdi. İkisi
de birincil kaynaktan kontrol edildi ve **ikisinin de teşhisi farklı çıktı.**

**`1_114_NdmB` — kurasyon kimyasal olarak DOĞRUYDU.** NdmB gerçekten
kafeini paraksantine çevirir ([UniProt H9N290](https://www.uniprot.org/uniprotkb/H9N290);
Summers 2012, J Bacteriol 194:2041, doi:10.1128/jb.06637-11, PMID 22328667).
Yani "hata" diye raporlanan şey bir hata değildi. Gerçek sorun başkaydı:
UniProt birincil reaksiyonu **teobromin → 7-metilksantin** veriyor ve
dengelediği Rhea kayıtlarının hepsi o reaksiyon üzerinde — üstelik
`official_equation` alanı **aynı YAML dosyasının içinde** zaten teobromini
yazıyordu. `substrate_en` kendi dosyasındaki resmî denkleme hizalı değildi.
Ayrıca kayıtlı substrat dizgisi NdmA'nınkiyle **birebir aynıydı**
("caffeine"), yani iki tip aynı işi yapıyor gibi görünüyordu; oysa ayırt
edici özellik demetillenen **pozisyon** (NdmA: N1, NdmB: N3). Substrat
teobromine hizalandı; kafein→paraksantin, teofilin→1-metilksantin ve
3-metilksantin→ksantin **ikincil belgelenmiş aktiviteler** olarak notlara
geçti.

**`2_208_TecA1` — kusur kimyada değil, KANIT ETIKETINDEYDI.** Kayıt
`chemistry_evidence: computed_from_official_equation` diyordu, yani oksijen
kaderinin resmî denklemden hesaplandığını iddia ediyordu; ama o denklem
**klorobenzen** üzerinde (RHEA:57512) ve kayıtlı substrat
1,2,4,5-tetraklorobenzendi. Yani kayıt, denklemin **hiç anmadığı** bir
molekül için resmî denklem desteği iddia ediyordu.
[UniProt O24676](https://www.uniprot.org/uniprotkb/O24676)'nın önerilen adı
"Chlorobenzene dioxygenase subunit alpha", birincil makalenin başlığı da
"Broad Spectrum **Chlorobenzene** Dioxygenase ... Dechlorination of
1,2,4,5-Tetrachlorobenzene" (Beil 1997, Eur J Biochem 247:190,
doi:10.1111/j.1432-1033.1997.00190.x, PMID 9249026). Substrat klorobenzene
hizalandı.

Tetraklorobenzen **atılmadı** — operonun adının kaynağı o, ve notlara
*doğru ürünüyle* girdi: enzim klor taşıyan karbonu dioksijenliyor, kararsız
diol klorür kaybıyla rearomatize olup **3,4,6-triklorokatekol** veriyor
(Beil 1997; sorumlu dizi öğeleri Beil 1998, J Bacteriol 180:5520). Eski
kayıttaki ürün adı "chlorinated cis-dihydrodiol" bu yüzden ayrıca
yanlıştı — o, izole edilen ürün değil geçici tür.

**Ders.** Yapısal ayırt edici iki tipi de doğru işaretledi, ama
"different_molecule" hükmü **"kurasyon hatası" anlamına gelmiyor**. Üç ayrı
durum aynı hükmün altına düşebiliyor: (a) kurasyon ikincil bir substratı
birincil diye yazmış, (b) kanıt etiketi denklemin kapsamadığı bir molekül
için destek iddia ediyor, (c) gerçekten yanlış molekül. İlk ikisi bu turda
görüldü, üçüncüsü görülmedi. Kapının WARN vermesi (FAIL değil) bu yüzden
doğruydu: fark gerçekti, teşhis ise literatür gerektirdi.

Yazma yolu: `types/<kume>.yaml` (tek kaynak) → `types_build.py --check`
farkın yalnızca o on alan olduğunu doğruladı → `--write --allow-drift`.
**`ROs_71_Clean/refs71.fasta` birebir aynı kaldı**, yani referans modelleri
ve ağır zincir yeniden koşmadı.

### 135.7 — Kontig uzunluğu bağlamı bozmuyor; bağlamsız girişler contig değil mRNA

Öneri şuydu: tüm genom olmasa bile **uzun parçalı contig'ler** varsa
veritabanına eklenebilir. Ölçüm gösterdi ki bu **zaten yapılmış** durumda,
ve iyi bir karar olduğu ölçülebilir.

Mevcut bileşim (17.073 replikon): **13.185'i `NZ_` WGS contig'i (%77)**,
1.991'i tam RefSeq kaydı (`NC_`/`NZ_C`), 1.897'si `XM_`/`NM_` mRNA kaydı.
Contig uzunlukları 593 bp ile 12 Mb arasında.

Bağlam, kontig uzunluğuna göre (onaylanmış girişler):

| band | giriş | komşusuz | ort. komşu | operon tam |
|---|---|---|---|---|
| <10 kb | 1.160 | **%82,7** | 1,0 | %0,3 |
| 10–50 kb | 1.678 | %0 | 15,8 | %6,1 |
| 50–200 kb | 4.006 | %0 | 18,9 | %6,8 |
| 0,2–1 Mb | 4.605 | %0 | 19,4 | %6,3 |
| >1 Mb | 7.676 | %0 | 19,6 | %5,8 |

Doğrudan karşılaştırma, ≥10 kb:

| | giriş | ort. komşu | operon tam |
|---|---|---|---|
| WGS contig | 13.029 | **18,8** | **%6,0** |
| tam genom kaydı | 4.936 | 19,7 | %6,5 |

Yani **≥10 kb bir contig, tamamlanmış genomdan daha kötü bir bağlam
vermiyor** — ne komşu sayısında ne operon tamlığında. Eşik keskin ve
10 kb'da.

**Asıl bulgu:** <10 kb bandındaki 1.160 girişin **958'i contig değil,
`XM_` mRNA kaydı** (354–7.881 bp, 1 CDS). Yayımlanan "bağlamsız giriş"
sayısı neredeyse tamamen mRNA kayıtlarından geliyor, kısa contig'lerden
değil. Bu, hayalet replikonla **aynı hata sınıfı**: bir kayıt, tablonun
birimi olduğunu iddia ettiği şey değil. Bir mRNA kaydı bir replikon değil,
ve replikon-başı paydada duruyor. **AÇIK İŞ** — silmek değil, ayrı bir
`record_kind` ile işaretlemek ve bağlam istatistiklerinin paydasından
çıkarmak (silmek, dizisi gerçek olan 958 girişi kaybettirir).

### 135.8 — Kapsama açığı ölçüldü: eksiklik dizide değil, taksonda

"Veri eksik geliyor" sezgisi bir payda sorusu, ve ölçüldü.

| birim | bizde | dış evren (Pfam PF00848, halka hidroksilleyen alfa katalitik domain) | oran |
|---|---|---|---|
| ayrı dizi / protein | **15.960** | 22.620 | ~%71 |
| ayrı takson | **4.565** | 12.532 | **~%36** |

Paydalar tam aynı birim değil (InterPro'nun sayımı UniProtKB kaydı,
bizimki genomdan türetilmiş giriş), ama 2,7 kat fark gürültü değil. Açık
**dizi sayısında değil, taksonomik genişlikte** — çünkü set belirli bir
genom kümesinden kuruldu, protein uzayının taramasından değil.

Bu, A4'ün yolunu değiştiriyor. Tüm RefSeq nr taraması diskte tıkanmış
(back2 139 GB, root 33 GB boş), **ama gerekmiyor**: PF00848'in protein
kümesi 22.620 kayıt, yani onlarca MB. InterPro API'sinden taxid'leriyle
çekip bizim 15.960 dizimizle karşılaştırmak hangi taksonların eksik
olduğunu tam olarak söyler — ağır tarama olmadan.

Sonra eksik taksonlar için **genom** indirilir, çıplak protein değil:
çıplak protein tam olarak 135.7'deki mRNA sorununu üretir (bağlamsız
giriş, replikon paydasında duran replikon-olmayan kayıt).

### 135.9 — Kapsama acigi: 1.928 dis cinsin 575'i bizde; eksiklik derinlikte degil genislikte, ama kapsam once DUZELTILDI

Kullanicinin "veriler eksik geliyo" sezgisi olculdu. Dis evren: Pfam
PF00848'in UniProtKB'deki tum proteinleri, InterPro API'sinden taxid'leriyle
cekildi (22.620 ayri accession).

**KARAR EKSENI TAKSON, DIZI DEGIL.** Iki tarafin dizi sayisi ayni birim
degil: InterPro bir UniProtKB kaydini sayiyor (ayni proteini tasiyan yuz
sus yuz kayit uretir), biz genomdan turetilmis bir girisi. Iki ayri sayma
kuralinin orani bir kapsama olcusu degildir. Bir taxid ise iki tarafta
ayni seydir -- bu yuzden 135.8'de yazilan "%71 dizi kapsamasi" rakami
GERI CEKILIYOR; o bir oran degil, iki paydanin karsilastirilamazligiydi.
135.8'in "dis evrende 12.532 takson" sayisi da yanlisti: InterPro'nun
taksonomi sayaci agacin IC dugumlerini de sayiyor. Gercek sayi, kaynak
organizma taxid'i olarak **8.125**.

| birim | bizde | dis evren | kapsama |
|---|---|---|---|
| ayri taxid | 4.251 | 8.125 | — |
| ayri cins | **720** | 1.928 | **575 ortak** |
| kapsam ici eksik cins | — | **798** | 3.050 protein |
| kapsam disi (okaryot) eksik cins | — | 555 | Fusarium, Aspergillus, Penicillium... |

**ACIK NEREDE.** Dis evrenin 22.620 proteininin **15.456'si zaten
kapsadigimiz cinslerde**; eksik cinslerde duran 6.278, kapsam icinde
yalnizca 3.050. Yani dizi sayisi bakimindan acik DERINLIKTE degil: cins
basina daha fazla protein indirmek kapsamayi degistirmez. Acik CINS
SAYISINDA: 1.373 cinsin 798'i hic temsil edilmiyor, ama her biri ortalama
4 protein tasiyor. Bu, "daha fazla contig" degil "daha fazla organizma"
cevabini dogruluyor (kullanicinin contig onerisi 135.7'de ayrica olculup
dogru cikti, ama acigi kapatan sey o degil).

**KAPSAM AYRIMI ONCE DUZELTILMEK ZORUNDAYDI -- ve ilk kosus bunu gosterdi.**
PF00848 bakteriyel degil; ayni domain bitki/fungus kolin monooksijenazinda
da var. Ilk kosusta kapsam "bizde en az bir kayit olan superkingdom" diye
tanimlanmisti ve o tanimla **1.244 eksik cinsin 1.244'u "kapsam ici"**
cikti; listenin basi `Clonostachys` (bir fungus, 67 protein) oldu. Sebebi
ayri bir kusurdu: bu DB'nin okaryot tarafi **958 kayit ve 958'inin TAMAMI
`XM_`/`NM_` bir mRNA transkripti** (930 + 28); tek bir okaryot genomu yok.
Yani kapsam tanimi 135.7'deki kusurdan besleniyordu, ve olculen sey
kapsama acigi olmaktan cikip bir kapsam hatasi oluyordu. Bu yuzden 135.7
once yapildi (asagida) ve kapsam `record_kind='genome'`'dan turetildi.

**CINS OLMAYAN ANAHTARLAR LISTEDEN CIKARILDI.** Dis evrenin 12 anahtari
(886 protein) bir cins degil: `rank:...` vekilleri (NCBI o soyagacinda
cins rank'li dugum tutmuyor) ve cozulemeyen taxid'ler. Bunlar bir INDIRME
onerisi olamaz -- karsiliginda indirilecek bir cins yok. Sayilari
dusurulmedi, ayri blokta duruyor.

**SAYFALAMA YONLU KONTROL EDILIYOR.** InterPro ilk sayfada 22.553 beyan
ediyor, cursor'u sonuna kadar yurumek 22.620 ayri accession veriyor; yani
`count` bir onbellek tahmini. Kapi bu yuzden esitlik DEGIL **yon** arar:
tehlikeli yon EKSIK cekmektir, cunku yarim bir dis liste paydayi kucultur
ve kapsama oranini oldugundan YUKSEK gosterir -- eksik veri, eksik
oldugunu gizleyen yonde hata uretir.

**CEKIM ARA KAYITLI OLMAK ZORUNDAYDI.** Ilk iki deneme 7.600 ve 2.600
kayitta InterPro'nun 408'ine takildi ve cekilen her sey cope gitti. 226
sayfalik bir cekim tek gecici hataya bagli olamaz: her sayfadan sonra
kayitlar ve siradaki cursor bir ara dosyaya yaziliyor, betik kaldigi
yerden devam ediyor. Ara dosya bir onbellek DEGIL -- tamamlanmamis bir
liste kapsama iddiasina girmemeli, o yuzden onbellek yalnizca cekim
bitince yaziliyor. Sayfa boyu da 200'den 100'e indirildi (sunucu 200'luk
sorguyu kendi zaman siniri icinde bitiremiyor).

**TLS DOGRULAMASI ACIK.** Bu depodaki eski `download_pfam.py`
`ssl._create_unverified_context()` kullaniyor; buraya kopyalanmadi (ilk
yazimda kopyalanmisti ve kaldirildi). Sertifika dogrulamasini kapatmak,
indirilen listenin KAYNAGINI dogrulanamaz kilardi -- ve bu liste bundan
sonra hangi genomlarin indirilecegine karar verecek.

**CINS KURALI IKINCI KEZ UYGULANMADI.** Dis evrenin taxid'leri bu DB'nin
hicbir tablosunda yok, dolayisiyla "o taxid'in cinsi ne" sorusu yeni bir
cozucu yazmaya cok yakindi. Yazilmadi: karar
`taxonomy_ids.genus_from_ranks`'a cikarildi ve hem replikon tablosu hem
dis karsilastirma oradan geciyor. Cikarma sonrasi `taxonomy_ids.py`
yeniden kosuldu ve `taxon_genus` + `replicon_taxon` + `taxonomy_ids.json`
**birebir ayni** cikti (ROADMAP 135.1'in dersi: ayni olcumun ikinci
uygulamasi birincisiyle sessizce ayrisir).

**KAPI ACIGIN BUYUKLUGUNE BAKMIYOR -- kasitli.** Eksik cins sayisi bir
kusur degil bir olcumdur; onu FAIL yapmak kapiyi veri buyutmeye zorlayan
bir mekanizmaya cevirirdi. Kapi uc seye bakiyor: dis listenin truncate
edilmemis olmasi, raporun cins paydasinin `taxon_genus`'tan BAGIMSIZ
yeniden sayimla uyusmasi, ve kapsamin beyan edilmis olmasi.

**BU SAYI NE DEGIL.** Bir cinsin PF00848'de olmamasi o cinste enzimin
olmadigini gostermez (UniProtKB'de temsil edilmemek bir kanit yoklugudur).
Bizde olup PF00848'de olmayan 144 cins hata degil: isabetlerimiz projenin
kendi modelinden, kendi esiginden geliyor.

Sonraki adim: 798 kapsam ici eksik cins icin **genom** indirmek, ciplak
protein degil -- ciplak protein tam olarak 135.7'nin uretigi sorunu
uretir.

Cikti: `analysis_out/coverage_audit.json`, onbellek
`analysis_out/pf00848_members.json.gz`, betik `coverage_audit.py`.

### 135.7b — `replicon` tablosunun 1.897 satiri replikon degil: kayit turu yazildi

135.7'de acik is olarak birakilan duzeltme yapildi. `replicon` tablosunun
adi her satirin bir replikon oldugunu soyluyordu; **1.897 kayit `XM_`/`NM_`
bir mRNA transkripti**, ve bunlarin 958'i dogrulanmis uye tasiyor. Bir
transkript ne kromozom ne plazmit; bir genin kopyasi. Hayalet replikon
satiriyla (135.2) ayni hata sinifi, ama daha sessiz ve daha pahali: bu
satirlar REPLIKON BASINA paydalarda duruyordu, ve bir transkriptin komsu
geni ya da tam operonu olamayacagi icin her baglam istatistigi 958 kayit
boyunca sifir toplayip paydada kaliyordu -- yani bir OLCUM EKSIGI, bir
biyolojik yokluk gibi yayimlaniyordu.

**SATIRLAR SILINMEDI.** Dizileri gercek ve dogrulanmis uye uretiyorlar;
silmek 958 girisi kaybetmek olurdu. Yapilan sey turu YAZMAK ve baglam
paydalarindan cikarmak.

**AYIRT EDICI KULLANILMADAN ONCE DOGRULANDI.** Akla gelen ilk olcut
accession onegiydi (`XM_` -> mRNA) -- ama bir onek bir dizgi, ve bu
projede bir karar biriminin dizgi isleminden turetilmesi tam olarak
kacinilan hata (bkz. `ro_taxonomy`, ROADMAP 135.1). Tur, kaydin KENDI
beyanindan okunuyor: GenBank LOCUS satirinin molekul tipi alani
(`mRNA` / `DNA`). Onek yalnizca CAPRAZ KONTROL: **17.073 kaydin
hicbirinde uyusmazlik yok** (0), yani iki olcut birbirini dogruluyor --
ama dogrulayan sey kaydin beyani, onek degil.

| record_kind | kayit | dogrulanmis giris |
|---|---|---|
| genome | 15.176 | 18.167 |
| transcript | 1.897 | 958 |

Ilk kosus 2.982 kaydi "okunamadi" isaretledi ve 7.758 dogrulanmis giris
oraya dustu. Sebep kayitlarda degildi: `gbk_files/` IKI BICIMI bir arada
tutuyor, bir kismi duz `.gbk`, bir kismi `.gbk.gz`. Yalnizca birine
bakmak, dosya duzeninden gelen bir eksigi kaydin kendi ozelligi gibi
gosteriyordu. Duzeltildikten sonra siniflandirilamayan satir: 0.

**BAGLAM PAYDASI BOLUNDU.** `stats_overview.py` baglamsiz girisleri zaten
disliyordu ("soru hic sorulmadi"), yani istatistikler yanli degildi --
eksik olan SEBEBIN ayrilmasiydi. Yayimlanan sayi artik iki parca:
`without_neighbours_because_transcript` ve
`without_neighbours_on_a_genome_record`. Tek sayi olarak iki sebep
ayirt edilemiyordu ve bu DB'nin genomlarinin parcaliligi oldugundan kotu
gorunuyordu.

Kapi dort sey ariyor: her replikon satirinin bir turu olmasi,
siniflandirilamayan satir kalmamasi, onek-beyan uyusmazligi olmamasi, ve
**transkript girislerinin silinmemis olmasi** -- bu son kontrol,
duzeltmenin "isaretle" degil "sil"e kaymadiginin muhafizi.

Cikti: `replicon_kind` tablosu, `analysis_out/replicon_kind.json`, betik
`replicon_kind.py`.

### 135.9b — Bir duzeltme grafikte degilse, duzeltme degil; bir kerelik temizliktir

135.2'de hayalet replikon satiri uc adimda kapatildi: sebep kodda
(`extract_genomic_context.py` artik satiri uretmiyor), veri bildirimle
(`db_repair.py`), ve muhafiz kapida. Bugun bulundu: **`db_repair.py`
yeniden kurulum grafiginde DEGILDI.**

Neden onemli: sebep kodda kapandi ama ARA ARTEFAKT eski kodla uretilmis.
`genomic_context/replicons.csv` o hayalet satiri hala tasiyor (13087.
satir, `nucleotide_id` bos, `status=no_annotation`) ve `build_db.py` o
CSV'yi okuyor. Yani bir sonraki tam yeniden kurulum satiri GERI GETIRIR
ve onu silecek hicbir sey yoktu. Kapi da onu FAIL ile yakalardi -- yani
sistem hatayi bulurdu ama kendi kendine kapatamazdi. Elle uygulanan bir
duzeltme bir KERELIK duzeltmedir; grafige konan bir duzeltme bir KURAL.

`db_repair` artik `database`'den hemen sonra, `--apply` ile grafikte;
idempotent oldugu icin uygulanacak bir sey yoksa yazmiyor, hata da
vermiyor, ve bir rapor (`analysis_out/db_repair.json`) uretiyor ki adimin
grafikte bir ciktisi olsun. Rapor YAYIMLANMIYOR: izin listesi
`atlas.PROVENANCE`'tan turedigi icin servis edilmiyor, ve icerigi
Turkce proje notu.

**AGIR ADIMI YENIDEN KOSTURMAK YERINE KANIT.** `genomic_context` adimi
10:57'deki kod degisiminden beri bayat ve 3,2 GB / 17.074 dosya okuyor.
Eski ve yeni kodun cikti farki TAM OLARAK o bir satir (olculdu: kimligi
ayristirilamayan tek dosya var, `NZ_QXDB01000002.1.gbk`, 0 bayt), DB'de o
satir yok, ve `db_repair` artik grafikte onu kapatiyor. Bu yuzden adimin
temel kaydi guncellendi, saatlerce yeniden hesap yapilmadi.

### 135.10 — Baglamsiz 959 girisin 958'i transkript: parcalilik sanilan sey kayit turuymus

135.7b'nin olcumu yayimlanan bir sayiyi tersine cevirdi. `stats.json`
"959 giriste hic komsu yok" diyordu ve bu, montaj parcaliligi olarak
okunuyordu (sayfa metni de oyle acikliyordu). Kayit turu ayrilinca:

| sebep | giris |
|---|---|
| kayit bir mRNA transkripti (komsu gen diye bir sey YOK) | **958** |
| genom kaydi, komsu gercekten bulunamadi | **1** |

Yani bu DB'nin genomik baglam cikarimi neredeyse kusursuz: **18.167 genom
girisinin 18.166'sinda komsu var.** Tek sayi olarak yayimlanan 959, bir
kayit turu olgusunu bir veri kalitesi kusuru gibi gosteriyordu -- ve yon
yine aleyhte: kendi isimizi oldugundan kotu gosteriyorduk.

### Acik is — grafik bayatligi (bu turda COZULMEDI)

`rebuild.py --check` 81 artefaktin 66'sini bayat gosteriyor ve kok
sebepler BU TURDAN DEGIL: `annotate_ro.py` icerigi temel kaydindan farkli
(f31a1074 ve 5d5ee0f8 commit'lerinden), oradan 60'tan fazla adim
kaskadliyor. Bu bir veri hatasi degil, bir BEYAN hatasi: ciktilar eski
kodla uretilmis ve bunu soyleyen tek sey bu bayatlik. Cozumu iki yoldan
biri: ya agir adimlari kosturmak (saatler), ya her kok icin 135.9b'deki
gibi KANIT uretip temel kaydini guncellemek. Kanitsiz `--baseline`
atilmamali: yesil ama yanlis bir grafik, kirmizi bir grafikten kotudur.

### 135.11 — 66 bayat artefaktin 60'i IKI bayat KAYITTI; grafik ilk kez 81/81 yesil

135.10'un sonunda acik is olarak birakilan grafik bayatligi kapatildi --
ve toptan `--baseline` ile DEGIL. Kok sebepler tek tek sinandi.

**KOK 1: `annotate` (60 adimlik kaskadin basi).** Bayatlik sebebi
`annotate_ro.py`'nin icerik ozetinin temel kaydindan farkli olmasiydi.
Commit'e (f31a1074) bakildiginda degisim SAF REFACTOR: satir ici `stale()`
fonksiyonu `stale_intermediate.stale` importuyla degistirilmis, ve iki
uygulama BIREBIR ayni mantik (dosya yok/bos -> "yok"; girdi daha yeni ->
sebep; yoksa None). Hesabin kendisine (hmmsearch cagrisi + ayristirma)
dokunulmamis.

Kanitla yetinilmedi, TANIK kosturuldu. `annotate`'in kendi atlama kurali
pahali kismi atliyor (ara dosyalar girdilerinden yeni), dolayisiyla kosu
ucuz. DB yedeklendi (md5 dogrulandi), `annotate_ro.py` elle kosturuldu ve
`ro` tablosunun atama sutunlari karsilastirildi:

| | once | sonra |
|---|---|---|
| ozet (candidate_id, cluster, group, coverage, score, evalue, is_confirmed) | `807f3cb46fc0e089` | `807f3cb46fc0e089` |
| satir / dogrulanmis | 26.693 / 19.125 | 26.693 / 19.125 |

**KOK 2: `operons`.** Bayatlik sebebi `build_operons.py` icerigi + uc
tablo. Bu commit'in etkisi GERCEK ve buyuktu (bilesen taramasi atlama
hatasi; beta 2.523 -> 4.372). Yani "etkisiz" demek icin tablolarin
duzeltme SONRASI durumu tasidigini gostermek gerekiyordu. Commit
mesajinin sayilariyla karsilastirma: besin dordu birebir (beta 4.372,
ferredoksin 2.736, has_beta 4.333, has_reductase 5.409), ama reductase
10.394 yerine **10.392** -- iki satir fark. ACIKLANMAYAN bir farkla temel
kaydi gecilmedi: betik kosturuldu ve bes tablonun hepsi birebir ayni
cikti (`operon`, `operon_gene`, `neighbor_protein`, `neighbor_component`,
`neighbor_domain`). Yani bugunun kodu 10.392 uretiyor; commit mesajindaki
10.394 daha eski bir andan. Tablolar duzeltme sonrasi.

**KOK 3: `bridging_geometry_audit`.** Tek beyan edilen girdisi kendi
betigi; betik 17:07:48, cikti 17:07:54 ve betige o tarihten sonra
dokunulmamis. Grafigin kendi kuraliyla cikti guncel.

**KOK 4 -- KODDA GERCEK BIR HATA: `site_freeze` HIC guncel olamiyordu.**
80/81 guncel oldugunda geride kalan tek bayat buydu, sebebi "girdi icerigi
degisti: analysis_out". Sebep: `site_freeze` girdi olarak `analysis_out`
DIZININI beyan ediyor, `dir_digest` o dizindeki her dosyayi okuyor, ve
grafigin KENDI durum dosyasi (`analysis_out/.rebuild_state.json`) adim
BITTIKTEN SONRA oraya yaziliyor. Yani adim kendi kaydiyla kendini
bayatlatiyordu -- hicbir kosu onu guncel yapamazdi. Kapinin grafik sonunda
olmasi gerektigiyle (135.1d) AYNI aile: **bir adimdan SONRA yazilan bir
sey o adimin girdisi olamaz.** `dir_digest` artik durum dosyasini atliyor.

Sonuc: **81 artefakt, 81 guncel, 0 BAYAT.** Onemi su: hic yesillenemeyen
bir kapi okunmayi birakir, ve bu grafik 2026-11-04'ten beri kirmiziydi.
Ama yesil olmasinin yolu kanitsiz `--baseline` degildi: iki kok gercekten
kosturulup cikti birebir karsilastirildi, biri mtime+ozet kaniti ile
kapandi, biri KODDA duzeltildi.

**YAN BULGU -- ve bu not DUZELTILDI (135.12).** Burada once "`annotate`
`ro`'ya yaziyor ama `DECLARED_WRITES`'ta beyan edilmemis, yani grafik onu
gormuyor" yazmistim. YANLISTI: `static_writes` betigin KAYNAGINDAKI SQL'i
ayristiriyor ve `annotate_ro.py` icin `{'ro'}` donduruyor, yani grafik o
yazmayi zaten biliyor. `DECLARED_WRITES` bir eksiklik listesi degil,
regex'in GOREMEDIGI yazmalar icin bir kacis kapisi. Gercek eksik baska
yerdeydi; bkz. 135.12.

### 135.12 — Sahipsiz DB yazmasi: FTS5'in golge tablolari kimseye ait degildi

135.11'in yan bulgusu yanlis yere isaret ediyordu (yukarida duzeltildi).
Dogru soru "hangi adim beyan edilmemis" degil, **"olculen hangi yazmanin
SAHIBI yok"**: bayatlik grafigi "bu adim hangi tablolari yazar" sorusunu
uc kaynaktan cevapliyor -- betigin kaynagindaki SQL (`static_writes`),
elle beyan (`DECLARED_WRITES`), ve adim kostugunda OLCULEN degisim. Ucuncu
digerlerinin disina tasarsa o yazma SAHIPSIZDIR, ve sahipsiz bir degisim
HERKESIN ustune yazilir.

81 adim tarandi, **tek bir adimda** sahipsiz yazma bulundu:

| adim | sahipsiz tablolar |
|---|---|
| `search_index` | `ro_fts_config`, `ro_fts_data`, `ro_fts_docsize`, `ro_fts_idx` |

Dordu de **FTS5 GOLGE TABLOSU**. `build_search_index.py` yalnizca
`CREATE VIRTUAL TABLE ro_fts` diyor; bu dort tabloyu SQLite'in kendisi
uretiyor ve bizim SQL'imiz adlarini hic anmiyor -- yani hicbir kaynak
taramasi onlari bulamaz. `DECLARED_WRITES` tam olarak bu durum icin var.

**Bedeli olculdu:** sahipsiz kaldiklari icin `ro_fts_data`'nin degismesi,
onu HIC okumayan uc rapor adimini (`report_html`, `explorer_html`,
`hub_html`) bayat gosteriyordu -- 135.11'de elle kosturulan uc adim tam
olarak bunlardi.

**Denetim artik kapida.** `check_db_writes_are_attributed` her adim icin
olculen yazmalari kaynak+beyan ile karsilastiriyor, cunku bir sonraki
sahipsiz yazmayi yine elle aramak bu projenin tekrar eden hatasi olurdu.
Kontrolun SINIRI da yazili: yalnizca OLCULMUS olani sinar -- bir adim bu
araçla hic kosmadiysa yazmasi henuz gorulmemistir. Yani "yazma yok"
demiyor, "gorulen her yazmanin sahibi var" diyor.

### 135.13 — "8 referans hic uye toplamiyor" WARN'i acildi: 6'si yontem kusuru, 2'si kapsama olgusu

Aylardir aciklamasiz duran WARN tek bir sayiydi ve IKI AYRI SEYI gizliyordu.
"0 uye" en az uc ayri durumun ortak belirtisi, ve zit cevaplar gerektiriyorlar:
enzim bu genomlarda YOK olabilir, ya da VAR ve atamasi rakip bir profile
gidiyordur. Ayirt edici, atamadan BAGIMSIZ ikinci bir olcum: `ro_evidence`in
kimlikce en yakin referansi (DIAMOND, kapsama >=%60; HMM skoru degil).

| tip | ona en yakin aday | en iyi aday kimligi | atama kime gitti | ref~ref | hukum |
|---|---|---|---|---|---|
| `1_103_CarAa` | 44 | **%99,7** | `1_102_CARDO` | %99,2 | var, baskasina atandi |
| `3_309_NahAc` | 20 | **%100,0** | `3_303_NBDO` | %82,0 | var, baskasina atandi |
| `4_407_PbaA1` | 2 | **%100,0** | `3_308_DxnA1` | %65,6 | var, baskasina atandi |
| `2_206_IpBAa` | 8 | %96,5 | `2_203_cumA1` | %90,8 | var, baskasina atandi |
| `2_210_XylC1` | 3 | %96,2 | `2_202_EdoA1` | %63,9 | var, baskasina atandi |
| `3_312_NinA` | 11 | %94,2 | `3_303_NBDO` | %91,3 | var, baskasina atandi |
| `2_207_BpdC1` | 20 | %47,8 | `2_208_TecA1` | %74,7 | veride yakin bir sey yok |
| `3_310_DntAc` | 3 | %31,8 | `5_506_Stc2` | %28,1 | veride yakin bir sey yok |

Yani **6'si eksik veri DEGIL**: tip bu genomlarda temsil ediliyor ve atama onu
baska bir profile veriyor. Yalnizca 2'si gercek bir kapsama olgusu. Ucunde ise
**kendi referans dizisiyle %99,5+ ayni bir giris var** ve tip yine de bos.

### 135.14 — ATAMA KUTUPHANESINDE IKI AD, BIR MODEL: `3_308_DxnA1` = `4_407_PbaA1`

En keskin vaka mekanik ve bir HATA. `RieskeDB71.hmm` icinde bu iki modelin
olasilik matrisi **BIREBIR AYNI** (govde ozeti `05868852bb819df5`; ad/accession/
tarih satirlari haric). Dogrulandi: ayni aday diziye ikisi de **971,1 bit ve
8.3e-293** veriyor -- farkli iki profil virgulden sonra ayni skoru vermez.

**Sebep hmmbuild degil, homolog TOPLAMASI.** Iki aile profili ayni dort
proteinden kurulmus (`ABQ71519.1`, `WP_336975089.1`, `WP_238839196.1`,
`WP_032073055.1`) ve her iki kume OTEKI TIPIN REFERANSINI da iceriyor -- yani
hicbiri kendi tipine ozgu degil. Iki referans dizisi %65,6 benzer, yani
DIZI duzeyindeki fazlalik esiginin (>=%99) cok altinda: **bu kusur dizi
denetimiyle BULUNAMAZ, cunku karar MODELDEN geciyor, diziden degil.**

**Bedeli kozmetik degil.** `3_308_DxnA1`'in iki dogrulanmis uyesi var ve
**ikisi de `4_407_PbaA1` referansiyla %100 benzer, %100 kapsamada** -- yani
ikisi de dizi olarak PbaA1. PbaA1 ise 0 uyeli gorunuyor. Iki tipin yayimlanan
etiketi ve kimyasi bu iki giris icin TERS: dibenzo-p-dioksin dioksijenaz diye
etiketlenen seyler 3-fenoksibenzoat dioksijenazlari.

**EGITIM KUMESI ORTUSMESI -- ve neden basliga CIKARILMADI.** 66 tipin profil
hizalamasi karsilastirildi: dort cift BIREBIR ayni kumeden kurulmus
(`5_502_CmoB`/`5_503_CmoA` 24/24, `3_306_3NTDO`/`3_310_DntAc` 109/109,
`1_102_CARDO`/`1_103_CarAa` 24/24, `3_308_DxnA1`/`4_407_PbaA1` 4/4) ve 8
uyesiz tipin HEPSINDE ortusmesi >=0,80 bir rakip var (7'sinde >=0,96).
Ama bu olcum YAYIMLANMADI ve kapiya BAGLANMADI, cunku dayanagi `aligned/`
dizini: yayima girmiyor, emekli tiplerin bayat dosyalarini tasiyor (7 dosya),
10 guncel tipin dosyasi hic yok, ve hizalamalarin kutuphaneyi kurduklarini
dogrulamaya calistigimda hmmbuild'in match-kolon kuralini birebir yeniden
ureteMEDIM (66'nin 30'u birebir uyustu, kalani birkac kolon kayik -- benim
yaklasiklik hatam, bayatlik degil). Kanit DB'den gelmek zorundaydi:
yayimlanan hukum `ro_evidence` + model ozeti uzerine kurulu, ikisi de
grafikteki artefaktlardan uretiliyor.

**KAPIDA, VE NEDEN WARN.** Uc kontrol eklendi: hicbir iki profil ayni matrisi
tasimamali, uyesiz her tipin bir TESHISI olmali (PASS), ve kendi referansiyla
%99,5+ ayni bir giris varken bos kalan tip olmamali. Ilki ve ucuncusu WARN --
FAIL degil, cunku (1) duzeltmek profil kutuphanesini yeniden kurup 19.125
girisi yeniden atamak demek ve "hangi homolog hangi tipe ait" bir KURATOR
karari, (2) FAIL bu kusuru DUYURAN sayfanin yayimini da kilitlerdi. Kapi
kuratorun yerine karar vermez, ama sessiz de kalmaz.

Yayimlandi: kalite sayfasinda "Where a type reports no members, and why" --
teshis tablosu ve iki-ad-bir-model vakasi, okuyucunun aciklanmamis bir sifir
gormemesi icin.

**ACIK IS (kurator karari).** DxnA1/PbaA1 profillerini ayri homolog
kumelerinden yeniden kurmak ve atamayi yenilemek; bu bir veri surumu.
Bu yapilmadan o iki girisin tip etiketi yanlis kalir -- su an sitede
olcumle beyan ediliyor, ama duzeltilmis degil.

### 135.15 — Yayimlanan bir dosyanin kokeni onu YAZMAYAN bir betigi gosteriyordu; eksik gecis yazildi

Kapinin "12 yayimlanan dosyanin tazeligi sinanamiyor" uyarisi incelendi ve
**uyarinin kendi ifadesi ALTISI icin yanlis cikti.** O alti dosya
(`hmm_pins.json`, `bridging_geometry_audit.json`, `reaction_anchors.json`,
`oxygen_fate_computed.json`, `reaction_verification.json`,
`curation_claims.json`) grafigin CIKTISI; yalnizca girdileri veritabani
degil (HMM dosyalari, YAML, yapi kayitlari), dolayisiyla tablo degisimi
onlari bayatlatmiyor ve tazelikleri KENDI beyan ettikleri girdilere karsi
sinaniyor. Yanlis etiketin bedeli: gercekten SAHIPSIZ olan alti dosyayi,
sahibi olan altinin arasinda GIZLIYORDU. Kontrol ikiye ayrildi; WARN artik
yalnizca hicbir adimin ciktisi olmayanlar icin, otekiler icin bir BEYAN
satiri var.

**Ve o sahipsiz altidan biri daha agirdi.** `curation_claims_ncbi.json`'in
kokeni `verify_curation_claims.py (NCBI efetch pass)` diyordu -- ama o
betik bu dosyayi YAZMIYOR ve depoda yazan hicbir sey yok. Yani sitede
**yanlis bir koken iddiasi** yayimlaniyordu.

Yetim dosyayi belgelemek yerine URETILEBILIR kilmak secildi, ve olcunce
daha buyuk bir bosluk cikti: UniProt gecisinin accession suzgeci (harfle
baslar, nokta YOK, <=10 karakter) GenBank/RefSeq kimliklerini disarida
birakiyor. **12 tipin accession'i bu bicimde ve yayimlanan dosya yalnizca
6'sini kapsiyordu** -- vakalarin yarisi hic kontrol edilmemis. Kapsanmayan
altidan biri de bugunun bulgusunun merkezindeki `3_308_DxnA1`.

Eksik NCBI gecisi yazildi: GenPept (`rettype=gp`) ayristirmasi, **onbellek
yetkili** (ag yalnizca `--fetch-ncbi` ile; yayimlanan bir dosyanin icerigi
servisin o anki erisilebilirligine baglanamaz -- ayni kural
`substrate_chebi` ve `pin_hmm_libraries`'de de gecerli), ve
`matches_stored` HER KOSUDA `types/`'a karsi yeniden hesaplaniyor,
onbellekten kopyalanmiyor: kopyalansa kurasyon degistiginde hukum sessizce
yanlis olurdu. Dizisi onbellekte olmayan kayit yeniden cekiliyor, cunku
dizi olmadan hukum cevrimdisi dogrulanamaz ve "dogrulandi" demek
olculmemis bir guvence olurdu.

**SONUC: 12/12 kapsandi, UYUSMAZLIK YOK.** Hepsi saklanan diziyle birebir.

| tip | accession | kaynak |
|---|---|---|
| `3_308_DxnA1` | `CAA51365.1` | "dioxin dioxygenase, large subunit", *Sphingomonas* sp., 435 aa, **birebir** |
| `4_407_PbaA1` | `A0A059U1S8` (UniProt) | "PbaA1", *Sphingobium wenxiniae*, 435 aa, **birebir** |

Bu, 135.14'u DARALTAN bir negatif sonuc: iki referans da KAYNAKTAN
dogru, farkli organizmalardan farkli iki enzim. Yani ayni-model kusuru
kurasyonda DEGIL, yalnizca profil kurulumunda. Duzeltme "hangi dizi
dogru" sorusunu sormuyor; yalnizca iki aile profilinin ayrik homolog
kumelerinden kurulmasini gerektiriyor.

**MUHAFIZ KENDI YENI KODUMU YAKALADI -- ve hakliydi.** Cins muhafizi
`verify_curation_claims.py:54`'u ihlal saydi: satirda hem `organism` hem
`.split()` var. Yanlis pozitifti (GenBank SOURCE satirinin bosluklarini
sikistiriyor, cins turetmiyor), ama cozum (a) satiri
`genus-string-intentional` ile isaretlemek DEGILDI -- o isaret "kasitli
cins dizgi tahmini" demek, ve yanlis beyan isaretin anlamini zayiflatir;
(b) muhafizin desenini daraltmak da DEGILDI -- desenin genis olmasi
kasitli, daha dar bir hali gercek bir ihlali kacirmisti. Islem kendi adi
olan bir fonksiyona (`collapse_space`) tasindi ve nedeni satirin yanina
yazildi.

Modularite butcesi 6 -> 5 dusuruldu (her baglanan dosyada dusurulur).
Kapi: 233 kontrol, 226 PASS, 0 FAIL, 7 WARN.

### 135.16 — "7 giris capasiz" uyarisi: vaka COZULMUSTU, ama kapinin GOREMEDIGI bir yerde

Uyari "7 entries are still decided on an unanchored reading" diyordu ve
aylarca durdu. Ilk olcum sayinin yapisini degistirdi: **7 giris daginik
degil, TEK BIR TIPIN tamami** -- `7_701_AerC`'nin yedi uyesi var ve
yedisi de ayni durumda. Yani soru "7 giris" degil, "AerC'nin karboksilati
var mi".

Ve cevap ZATEN VERILMISTI. Tipin `notes` alani, kurator onayiyla
(2026-10-21):

> CARBOXYLATE IS Asp667. ESMFold modeli His513 NE2 ile His518 NE2'yi
> 3,22 A ayri yerlestiriyor (NDO kristali 3,49), Asp667 OD1 2,64 A ve
> OD2 3,21 A mesafede (NDO Asp362: 3,41 / 3,20); ikinci aday Glu510
> 4,18 A, yani atama belirsiz degil. Kolon testi hala 2/3 veriyor cunku
> karboksilat kolonunu tasiyan model bolgesinin burada karsiligi yok --
> bu, protein hakkinda degil HIZALAMA hakkinda bir ifade.

Olcum var, gerekce var, onay var. Ama hepsi DUZ METIN, ve kapi `ro_triad`
tablosunu okuyor. Sonuc: cozulmus bir vaka, cozulmemis diye uyari uretti.
**Bu projenin tekrar eden hatasi: mekanizmanin GOREMEDIGI bir yere
yazilan karar, verilmemis karardir.** (Ayni ders TecA1'in urun adi
kararinda da yasandi; o da `reaction_name_exceptions.yaml`e tasinana
kadar gorunmuyordu.)

**KARAR BILDIRIME TASINDI:** `carboxylate_structure_anchors.yaml`. Her
kayit `residue`, `method`, `measured`, `runner_up`, `limits`,
`decided_on` tasimak ZORUNDA; alani eksik kayit YUKLENMEZ ve STDERR'e
yazilir -- gerekcesiz bir istisna, istisnasiz bir hatadan kotudur.

**IKI SEY BILEREK YAPILMADI.**

1. `triad_complete` bayragina DOKUNULMADI. O bayrak KOLON okumasinin
   olcusu ve yayimlanan 3/3 sayilari ona bagli; bir yayimlanan sayiyi bir
   YAML'a satir ekleyerek oynatmak, olcumu degil MUHASEBEYI degistirmek
   olurdu. Olculdu: `triad_complete` 19.064'te kaldi. Kayit yalnizca
   durumu tasidi: `carboxylate_untested` -> `carboxylate_structure_anchored`.
   "Olculemedi" ile "baska bir yontemle olculdu" ayni kovaya girmemeli.
2. Sutunun TURU korundu. Ilk yazimda `ro_triad.carboxylate_residue`'ya
   "Asp667" yazmistim; o sutun tek harfli kalinti kodu tasiyor. Kodun
   kendisi sutuna, pozisyon nota, tam gerekce kayit dosyasina gidiyor --
   ve uc-harf->bir-harf eslemesi yalnizca Asp/Glu taniyor, yani kayda
   yanlis bir kalinti adi yazilsa esleme BULUNMAZ ve sutun None kalir;
   yanlis kalinti sessizce kabul edilmez.

**KAPIYA IKI KONTROL.** Bu mekanizma bir uyariyi SUSTURUYOR, dolayisiyla
iki yoldan bozulabilir: kayit gerekcesiz kalir, ya da BAYATLAR (tipin
kolonu artik okunabiliyorsa kayit gereksizdir ve durmasi gelecekte gercek
bir olcumun ustunu orter). Ikisi ayri ayri sinaniyor.

Sonuc: `absence_not_proven` 7 -> 0, `his_undecided` 0, uc triad kontrolu
yerine BES kontrol ve hepsi PASS. Kapi: 235 kontrol, 228 PASS, 0 FAIL,
WARN 7 -> 6.

| triad durumu | giris |
|---|---|
| complete | 19.031 |
| complete_non_canonical_his_spacing | 33 |
| his_absent | 31 |
| carboxylate_absent | 23 |
| **carboxylate_structure_anchored** | **7** |

### 135.17 — A6'nin ikinci yarisi: cekici yazildi ve OLCULDU, ama cekim TAMAMLANMADI

InterPro / KEGG / MetaCyc tanimlayicilari `cross_references.py`ye eklendi
(bicim deseni, URL sablonu, dizi-basi birim, ve "bu uclu de birer IDDIA,
salt adres degil" notu). Veriyi getiren `fetch_uniprot_xrefs.py` yazildi.
**Cekim bitmedi: 50 partinin 16'si tamamlandi (24.000/73.699 accession).**

**YAYIMLANAN HICBIR SAYI DEGISMEDI** ve bu kasitli. `cross_references.py`
onbellek tam degilse uc yeni tanimlayiciyi EKLEMIYOR; `ro_xref` eski dort
veritabaniyla ayni sayilarda duruyor (uniprot 73.699, embl 24.911, refseq
18.617, ec 3.820) ve rapor `cache_complete: false` ile hangi
tanimlayicilarin TUTULDUGUNU yaziyor. Gerekce: yarim bir cekim "bu dizinin
InterPro kaydi yok" gibi gorunur, yani yayimlanan oran VERIYI degil
cekimin nereye kadar geldigini olcer -- ve hata yonu aleyhte, eksiklik
eksik oldugunu GIZLER.

**SERVISIN DAVRANISI OLCULDU; BES DENEME, BES AYRI KUSUR.** Dokumanlar duz
metin olarak alinamadi (yardim sayfasi ve Swagger kabuk dondu), o yuzden
tasarim problarla belirlendi:

| olcum | sonuc |
|---|---|
| 3 kimlik | calisiyor |
| 1.500 kimlik (govde 19,5 KB) | **3 sayfa, 1.500/1.500 tekil kimlik** |
| 5.000 (65 KB) / 20.000 (260 KB) / 73.699 (810 KB) | **0 kimlik**, is 1 saniyede "bitiyor" |
| `/results/stream/` | 1000 satirda KESIYOR (5.000 ve 3.699 kimlikte de 1000) |
| sayfali `/results/` + `Link: rel="next"` | tam -- AMA araliki olarak Link'i dusurup 1000'de kesiyor (19 partinin 5'i) |
| cok sayida is, araliksiz | durum yoklamalarina HTTP 400; ~90 isten sonra birikimli |

Yani belgelenen "100.000 kimlige kadar" sinirini HIC OLCMEMISTIM ve
olcunce yanlis cikti. Servis reddi SESSIZ yapiyor: yine bir `jobId`
veriyor, durum aninda sonuclara yonleniyor, kume bos.

**BES KORUMA, her biri yapilan bir hataya karsilik:**

1. "Bitti"nin isareti `FINISHED`, sonuclara YONLENME ya da yanitta
   `results` olmasi. Ilk yazim `jobStatus` ALANININ YOKLUGUNU basari
   sayiyordu -- hata yanitlari da o alani tasimiyor, ve is hazir olmadan
   akis cekildi.
2. SIFIR sonuc onbellege yazilmaz, parti yeniden denenir.
3. KESIK sonuc da yazilmaz, ve bu AYRI bir kontrol olmak zorundaydi:
   sifiri reddeden koruma 1000/1500 donen partileri KABUL ediyordu.
   Ayirt edici: gercek eksiklik (silinmis accession) rastgele bir
   sayidir, KESILME sayfa boyunun TAM KATIDIR.
4. Sutun duzeni BASLIKTAN okunur. ID-mapping sonuclari basa bir `From`
   sutunu ekliyor; ilk yazim `parts[0]`i accession sanip `parts[1:4]`u
   xref okuyordu, yani veri gelse SESSIZCE yanlis sutunlari okuyacakti.
5. Partiler TEK TEK, aralarinda bekleyerek (varsayilan 20 s).

Devam edebilir: tamamlanmis partiler atlanir, `_fetch.json` gonderilen /
donen / eksik sayilarini tasir.

**ACIK IS.** Cekimi UniProt sinirlamasi gecince tamamlamak (kalan 34
parti). Tamamlaninca `cross_references.py` uc tanimlayiciyi kendiliginden
ekler -- kapi yalnizca TAMLIGA bakiyor, elle bir bayrak yok.

**KENDI PAYIMA DUSEN DERS.** Bes denemenin kok sebebi tekti: her adimda
servisin davranisini VARSAYIP sonra olctum, tersi yerine. Alan adlarini
sinamistim ama AKISI sinamamistim. Veri butunlugu zarar gormedi (grafik
82/82, kapi 0 FAIL, yayimlanan sayi degismedi); harcanan sey zaman ve
istek.

### 135.18 — A6 tamamlandi: InterPro / KEGG / MetaCyc yayimda, ve ID-mapping YOLU BIRAKILDI

73.699 accession'in tamami cekildi (eksik 0) ve `ro_xref` dort yerine
YEDI veritabani tasiyor:

| tanimlayici | satir | dizi | uniprot kaynakli payi |
|---|---|---|---|
| interpro | 109.279 | 25.313 | **%36,8** |
| kegg | 5.511 | 5.297 | %7,7 |
| metacyc | 117 | 111 | **%0,2** |

MetaCyc'in %0,2'de kalmasi bir olcum: MetaCyc kuratorlu bir YOL
veritabani ve bu dizilerin neredeyse tamami incelenmemis (TrEMBL) kayit,
yani dis dunyada da bu proteinlerin cogunun yola baglanmis kaydi yok.

**ID-MAPPING YOLU TERK EDILDI.** Toplu is icin "dogru arac" o
gorunuyordu, ve dort tasarim denemesi sonunda su olculdu: servis cok
sayida is gonderildikten sonra isi KABUL EDIP sifir sonuc donduruyor --
yani red SESSIZ, hala bir jobId veriyor. 73.699 accession'in ancak
27.000'i o yolla gelebildi.

**DOGRU YOL EN BASIT ARACTI:** duz `/uniprotkb/search` ile
`accession:(A OR B ...)`. Olculdu: 100 accession 100/100 donuyor ve
SAYFALAMA bile gerekmiyor; 200 accession'da URL uzunlugu HTTP 400
veriyor. Asenkron is yok (sessiz-sifir yok), cursor yok (kesilme yok),
kuyruk kotasi yok. 468 sorgu, SIFIR dusen.

**SEKIZ DENEME, SEKIZ KUSUR -- ve hicbiri yayimlanan bir sayiya
donusmedi.** Her birini ya tamlik kapisi, ya yuvarlak-sayi kuskusu, ya da
atilan-sayi raporu durdurdu:

1. `jobStatus` ALANININ YOKLUGUNU "is bitti" sayma (hata yanitlari da o
   alani tasimiyor).
2. SIFIR sonucu onbellege yazma.
3. KESIK sonucu onbellege yazma -- AYRI bir kontrol gerekti; sifiri
   reddeden koruma 1000/1500 donen partileri kabul ediyordu. Ayirt
   edici: gercek eksiklik rastgele, KESILME sayfa boyunun TAM KATI.
4. Gonderim boyutunu BELGEDEN varsayma (100.000 degil, ~1.500).
5. Sutun duzenini SAYMA: ID-mapping sonuclari basa `From` ekliyor.
6. Dosya adini KONUMDAN turetme: `chunk_{n}` kullanildi, ama `todo`
   listesi her kosuda degisiyor; ikinci kosu birincinin 1-100 numarali
   dosyalarini "zaten var" diye atladi ve TAM 10.000 accession sessizce
   cekilmemis kaldi. Ad artik parcanin ICERIGININ ozeti.
7. Okuyucu desenini `page_*` birakma: arama dosyalari `chunk_*` ve
   okuyucu 73.699'un 8.852'sini gordu.
8. PAYDA BIRIMINI karistirma: tamlik kontrolu "en az bir tanimlayicisi
   olan accession" ile "beklenen accession"i karsilastiriyordu -- iki ayri
   birim, ve o kapi ASLA gecemezdi. Rapor artik iki sayiyi ayri tasiyor
   (`accessions_in_cache`, `accessions_with_at_least_one_identifier`) ve
   neden ayri olduklarini yaziyor.

**VE KAPANIRKEN IKI KUSUR DAHA:**

9. KEGG deseni 22 GECERLI kimligi atti: organizma kodu 3-5 harf sanilmis,
   ama KEGG'in `ag` (genes-addendum) kodu IKI harf (`ag:AAB52509`).
   Atilan alt kume en degerlisiydi: addendum, tam genomu olmayan ama
   LITERATURDE karakterize edilmis proteinleri tutuyor. Ayni hata bu
   dosyada EMBL deseniyle de yasanmisti (5.358 kayit) ve iki kere de
   ATILAN SAYININ RAPOR EDILMESI yakalatti.
10. Kapi URL kalibi listesini ELLE SAYIYORDU (dort ad). Uc yeni
    tanimlayici eklenince FAIL verdi, oysa kaliplar
    `cross_references.URLS`de VARDI; eksik olan kapinin IKINCI KOPYASIYDI.
    Liste artik tek sahipten okunuyor ve kapi SITENIN kalibi tasidigini
    siniyor (`webapp/app.py:XREF_DBS`).

**UCU DE BIRER IDDIA, salt adres degil** -- InterPro bir domain, KEGG ve
MetaCyc birer yol iddiasi, ve ucu de EC gibi otomatik anotasyonla
yayiliyor. Disa baglanti olarak basiliyorlar, bu veritabaninin kendi
ifadesi olarak degil; rapor bunu ayrica yaziyor.

**KENDI PAYIMA DUSEN DERS.** Sekiz denemenin kok sebebi tekti: servisin
davranisini VARSAYIP sonra olctum. Dogru cozum (100 accession'lik duz
arama sorgusu) ilk dakikada denenebilirdi. Kontroller isini yapti, ama
kontrollerin isi hatayi yakalamak; hatayi AZALTMAK olcumle baslamaktan
geciyor.

### 135.19 — Dis dogrulama: Rieske cagrisinda InterPro ile uyum 25.174/25.174, ve asimetri bir bulgu

135.18'de cekilen InterPro tanimlayicilari ilk kez BAGIMSIZ bir ikinci
gorus verdi. Bu veritabaninin en cok yuk tasiyan cagrisi "bu dizi bir
Rieske proteini"dir -- aday havuzuna girmenin kosulu o -- ve bugune kadar
yalnizca projenin KENDI modellerine dayaniyordu.

InterPro'nun anotladigi 25.313 dizi uzerinde:

| olcum | dizi |
|---|---|
| bizim `rieske_intact` dedigimiz | 25.174 |
| **InterPro'nun da Rieske dedigi** | **25.174** |
| bizim deyip InterPro'nun DEMEDIGI | **0** |
| InterPro'nun deyip kapimizin DISLADIGI | 139 |

Kayit adlari KAYNAKTAN alindi, tahmin edilmedi (InterPro API):
`IPR017941` = Rieske [2Fe-2S] iron-sulphur domain, `IPR036922` = ayni
domainin homolog super ailesi.

**BAGIMSIZLIGIN SINIRI YAZILI, ve bu cumle olcumun kendisi kadar onemli.**
BAGIMSIZ olan: bir dizi VERILDIGINDE InterPro'nun ona Rieske domaini
koyup koymadigi -- o profiller baska insanlar, baska veri, baska esikle
kuruldu. BAGIMSIZ OLMAYAN: hangi dizilerin SORULDUGU -- aday havuzu bizim
kendi Rieske/motif modellerimizle secildi, yani InterPro hic
koymadigimiz bir dizi hakkinda itiraz EDEMEZ. Dolayisiyla bu olcum
KESINLIK (precision) hakkinda konusuyor, DUYARLILIK (recall) hakkinda
HICBIR SEY soylemiyor. Kapi bu beyanin varligini siniyor, cunku beyan
olmazsa okuyucu orani recall sanir.

**PAYDA DA IKI TANE, ve ikisi ayri yazili.** 73.299 birlesik dizinin
yalnizca 25.313'u InterPro'da bir kayit tasiyor. Geri kalan 47.986
CELISMIYOR, ANOTSUZ (cogu yeni, incelenmemis kayit). Tum diziler
uzerinden okunursa oran "%35 uyum" gibi gorunur ve o sayi veriyi degil
InterPro'nun KAPSAMASINI olcer.

**ASIMETRI TEK YONLU, VE BIR BULGU.** 139 dizide InterPro Rieske domaini
goruyor, bizim kapimiz disliyor; hepsi `passes_gate=0`. Celiski DEGIL:
InterPro'nun olcusu profil benzerligi, bizimki DORT KUME LIGANDININ
yerinde olmasi, ve kumesiz bir Rieske katlanmasi kimyayi yapamaz. Ama
dusurme SESSIZ kalmamali -- gorunmeyen bir dusurme bu projenin tekrar
tekrar buldugu hata bicimi:

| mononukleer ucgen bolgesi | dizi |
|---|---|
| **3/3 tam** | **76** |
| 2/3 | 50 |
| 1/3 | 7 |
| 0/3 | 6 |

Yani 76 dizide kume ligandi degismis ama katalitik ucgen SAGLAM. Bu,
kuratorun bakmasi gereken bir kume ve artik kalite sayfasinda duruyor.

**AILE ADLARI CELISKI DEGIL, GENISLIK.** InterPro ayni kumede
`IPR050584` kolesterol 7-desaturaz (11.837 dizi), `IPR044043` vanillat
O-demetilaz, `IPR045605` 3-ketosteroid 9-alfa-monooksijenaz,
`IPR013626` feoforbid a oksijenaz, `IPR045623` LigXa-benzeri kayitlarini
da veriyor. Bu veritabani Rieske oksijenazlarin TAMAMINI topluyor,
yalnizca halka hidroksilleyenleri degil; o adlar bu genisligin disaridan
gorunusu.

Cikti: `analysis_out/interpro_agreement.json`, betik
`interpro_agreement.py`, dort kapi kontrolu, kalite sayfasinda
"A second opinion on the call everything rests on".

### 135.20 — Asimetri izlendi ve 135.19'un YORUMU DUZELTILDI: bagiyici bolge kaymis degil, YOK

135.19'da InterPro'nun Rieske dedigi ama kapinin disladigi 139 dizi icin
sunu yazmistim: "76'sinda katalitik ucgen tam, yani kume ligandi degismis
ama merkezi saglam -- kuratorun bakmasi gereken bir kume". **Olcum bunun
YANLIS oldugunu soyluyor.**

**ADIM 1: HANGI LIGAND?** Diziler motif modeline hizalanip dort Rieske
kolonu (Cys-1 164, His-1 166, Cys-2 184, His-2 187) tek tek okundu.
Baskin desen tek bir kalintinin degismesi DEGIL, bir CXH CIFTININ
BIRLIKTE gitmesi:

| eksik desen | dizi | ucgen 3/3 olanlarda |
|---|---|---|
| Cys-1 + His-1 | 50 | 39 |
| Cys-2 + His-2 | 23 | 12 |
| yalniz His-2 | 20 | 3 |
| yalniz His-1 | 19 | 7 |
| yalniz Cys-1 | 8 | 6 |

Kume iki ayri ilmekle tutulur (C-x-H ve C-x-x-H); bir CIFTIN birlikte
gitmesi nokta mutasyonu degil, ilmegin yerinde olmamasi demek. 139'un
73'u, ucgeni tam olan 76'nin 51'i boyle.

**ADIM 2: KAYMIS MI, YOK MU -- ve bu iki okuma ZIT.** Bir ligand iki
sebeple bulunamaz: gercekten yok, ya da yerinde ama kolon testinin +-2
penceresinin disina kaymis. Ikisi AYNI sonucu verir, ZIT anlamlari var:
"bu protein kume baglamaz" ile "bizim hizalamamiz kaciriyor". Ayni ayrim
bu projede 7_701_AerC'de de karsimiza cikmisti (135.16): "kolon testi
ulasamiyor" cumlesi protein hakkinda degil HIZALAMA hakkindaydi.

Bu yuzden kanonik Rieske motifi (`C.H.{10,30}C..H`) HAM dizide,
hizalamadan ve kolon numarasindan BAGIMSIZ arandi:

**139 dizinin 138'inde motif HIC YOK.**

**AYIRT EDICI KULLANILMADAN ONCE DOGRULANDI** -- ve neredeyse
atlaniyordu. Ayni desen, kume ligandlarini SAGLAM dedigimiz dizilerde ve
kuratorlu referanslarda sinandi: **4.000/4.000** ve **76/76** tutuyor.
Yani yokluk bulgusu gercek, desen kusuru degil. Dogrulama rapora GOMULDU
(`validated_on`), elle kosulan bir kontrol olarak kalmadi.

**DUZELTILMIS OKUMA.** Bagiyici bolge kaymis degil, dizide YOK. InterPro
bu 139 diziyi Rieske super ailesine KATLANMA BENZERLIGINDEN koyuyor; kume
bolgesi gitmis. Dolayisiyla bu, kapinin dort-ligand kosulunun simdiye
kadarki EN GUCLU dis sinavi: kapi onlari dogru disliyor. Kalite
sayfasindaki yorum duzeltildi.

**KENDI HATAM, ve dersi onemli.** Ilk kosuda 139 dizinin 137'sinde dort
ligandin HEPSI "bosluk" cikti -- sacma bir sonuc, ve onu yakalatan da
sacmaligi oldu. Sebep: `hmmalign` ciktisinda insert kolonlari var, yani
hizalama dizgisindeki 164. konum modelin 164. ESLESME kolonu DEGIL.
Dogru yol projede ZATEN vardi (`ro_motif.read_stockholm_matchcols`).
Kolonlari tek sahipten almistim ama HIZALAMAYI OKUMAYI elle yazdim.
Kural "sabitleri tek sahipten al" degil, **"olcumun TAMAMINI tek
sahipten al"** olmaliydi.

Cikti: `analysis_out/rieske_ligand_audit.json`, betik
`rieske_ligand_audit.py`.

### 135.21 — 124.3'un kalan yarisinin PAYDASI olculdu: dizi kayitlarin %20'sinde var

124.3'un acik yarisi "sigma70 kutusu zenginlesmesini uzunluk ve GC
eslenmis bir null'a karsi olc" diyor ve yaninda bir VARSAYIM duruyordu:
"intergenik DNA'nin `gbk_files/`ten cikarilmasini gerektiriyor". O
varsayim sinanmamisti; sinandi.

**Kayitlarin cogu CON tipi: anotasyon var, DIZI YOK.** Ornek:
`NC_002033.1` 184.457 bp beyan ediyor, dosya 231 KB, ve
`CONTIG join(AF079317.1:1..184457)` ile bitiyor -- ozellikler ve protein
cevirileri var, nukleotid dizisi yok.

| birim | dizisi var | toplam | pay |
|---|---|---|---|
| replikon kaydi | 4.882 | 17.073 | %28,6 |
| **GENOM kaydi** | **2.985** | **15.176** | **%19,7** |
| dogrulanmis giris (genom kaydinda) | **7.761** | 18.167 | **%42,7** |

**TRANSKRIPTLER ORANI SISIRIYORDU -- ve 135.7 tam bunu yakalatti.** Dizisi
olan 4.882 kaydin 1.897'si `XM_`/`NM_` transkripti; transkriptlerin
%100'unde dizi var, genom kayitlarinin yalnizca %20'sinde. Ama bir
transkriptin YUKARI AKIS INTERGENIK BOLGESI YOK, yani bu olcum icin
ise yaramaz. `record_kind` ayrimi olmasa "dogrulanmis girislerin %45,6'si
okunabilir" diye yazacaktim; dogru sayi %42,7 ve payda GENOM kaydi.

Dosya bicimi de ayrisiyor: gzip'li 2.982 dosyanin TAMAMINDA dizi var,
duz 14.092'nin yalnizca 1.900'unde. Yani iki kume farkli zamanda, farkli
yolla indirilmis.

**NE YAPILDI, NE YAPILMADI.** Yapilan: `replicon_sequence_census.py`, her
replikon icin diziyi okuyup okuyamadigimizi kaydediyor (dizinin kendisini
OKUMAZ, yalnizca `ORIGIN` satirini arar -- gereksiz gigabayt okumamak
icin). Yapilmayan: zenginlesme olcumunun kendisi. Once payda, sonra
iddia; cunku okunamayan bir bolgeyi paydada tutmak "motif yok" gibi
gorunur ve bu projenin tekrar eden hata bicimi tam olarak bu.

**SONRAKI ADIMIN SINIRI simdiden yazili.** Olcum 7.761 giris / 2.985
replikon uzerinde yapilabilir. Null ayni bolgelerden uretildigi icin
(dizi karistirma + ayni genomlarin baska intergenik bolgeleri, GC ve
uzunluk eslenmis) karsilastirma IC gecerliligi tasir; okunamayan yariya
GENELLENIP genellenemeyecegi ise OLCULMEMIS kalir ve oyle beyan
edilmelidir.

Cikti: `analysis_out/replicon_sequence.csv`, betik
`replicon_sequence_census.py`.

### 135.22 — 124.3 KAPANDI: sigma70 kutu cifti eslenmis null'a gore 1,7-2,8 kat zenginlesmis

124.3'un kalan yarisi uc adimda kapandi, ve her adim bir sonrakinin
yayimlanabilir olup olmadigini belirledi.

**ADIM 1 -- PAYDA (135.21).** Varsayim sinandi: kayitlarin %80'inde dizi
YOK (CON tipi). Olcum 7.616 bolge / 2.285 replikon uzerinde yapilabilir.

**ADIM 2 -- CIKARMA.** `upstream_regions.py` intergenik DNA'yi cikarir ve
KENDI CIKTISINI DOGRULAR: bolgenin uzunlugu `ro_regulation.intergenic_bp`
ile karsilastirilir, ki o sayi BASKA bir betik tarafindan koordinatlardan
hesaplanmis. **7.616/7.616, sifir uyusmazlik.** Iplik de olcumun parcasi:
operon eksi iplikteyse bolge ters tumleyen yazilir, yoksa "kutu yok"
sonucu bir YON hatasindan gelir.

**ADIM 3 -- OLCUM.** Iki null, uc duzey, karari cins verir:

| duzey | gozlenen | beklenen (mono / markov) | oran | p |
|---|---|---|---|---|
| giris | 406 | 240,7 / 218,8 | **1,69× / 1,86×** | 6e-23 / 1,6e-30 |
| dizi | 363 | 218,0 / 197,7 | 1,67× / 1,84× | 8e-20 / 1,2e-26 |
| **cins** | **184** | **75,0 / 65,1** | **2,46× / 2,83×** | 1,5e-27 / 8,1e-35 |

Iki null birbirini CERCEVELIYOR: dinukleotid kompozisyonunu koruyan
Markov null'i daha DUSUK beklenti veriyor, yani daha iyi null'a karsi
zenginlesme daha buyuk. Dort mimarinin hepsi kendi null'unun ustunde
(1,5-2,0×); duzenleyici gruplarinin hem orani hem null'u dusuk, cunku
bolgeleri kisa (168 bp vs 247 bp) ve null bunu izliyor. **Mimariler arasi
zenginlesme farki TEST EDILMEDI ve iddia EDILMIYOR.**

**IKI ISTATISTIK HATASI, ikisi de duzeltildi ve kaydedildi.**

1. ILK TEST YANLIS BIRIMDEYDI. Isaret testi, IKILI bir gozlemle (cift
   var/yok) SUREKLI bir beklenti (20 kopyanin ortalamasi) arasina
   konmustu; cift tasimayan bir bolge, null 20'de 1 buldugu icin "null
   lehte" sayiliyordu, yani -0,05 ile +0,95 ayni agirlikta. Sonuc KENDI
   ICINDE CELISKILIYDI -- isaret testi null'u, oran karsilastirmasi
   gercegi kazandiriyordu -- ve yakalatan sey bir kontrol degil
   CELISKININ KENDISI oldu. Dogru test SAYILAR uzerinde: gozlenen birim
   sayisi ile null altinda beklenen sayi, iki yonlu binom.
2. BINOM TASTI. `math.comb(7616, k)` float'a sigmiyor
   (`OverflowError`). Log uzayina tasindi, ve bu kez KULLANILMADAN ONCE
   bilinen degerlerle sinandi: n=10/obs=8/p=0,3 -> 0,00318 (elle
   hesaplanabilir), n=100/obs=50/p=0,5 -> 1,0 (gozlenen beklentiye esit
   oldugunda p=1 olmali). Araci dogrulamak, ayirt ediciyi dogrulamakla
   ayni sey.

**YENIDEN URETILEBILIR.** Tohum sabit (`SEED = 20261113`) ve grafikten
yeniden kosturuldugunda `sigma70_enrichment.json` ile
`upstream_regions.csv` BIT BIT AYNI cikti (md5 dogrulandi). Rastgele
null'a dayanan bir iddia icin bu zorunlu.

**KAPI DISIPLINI DENETLIYOR, SONUCU DEGIL.** Bes kontrol: uc duzey var
mi, her duzeyde iki null var mi, **iki null AYNI YONU gosteriyor mu**
(betigin yayin kurali buydu; kural yalnizca docstring'de kalirsa bir
sonraki kosuda ayrisan iki null sessizce yayimlanir), payda uyarisi
yazili mi, ve ne iddia edilmedigi yazili mi.

**IDDIA EDILMEYEN.** Hicbir promotor bulunmadi ve tahmin edilmedi --
kutular bilgi bakimindan zayif, rastgele de sik bulunur. Null ayni
bolgelerden geldigi icin kompozisyon kontrol edildi, ama bolgeler bu
veritabaninin kendi operon cagrilariyla secildi. Ve olcum yalnizca dizisi
okunabilir %43'te yapilabildi; okunamayan cogunluga genellenip
genellenemeyecegi OLCULMEMIS.

Cikti: `analysis_out/sigma70_enrichment.json`,
`upstream_regions.{csv,fasta.gz}`, betikler `upstream_regions.py` ve
`sigma70_enrichment.py`, kalite sayfasinda "Sigma70 boxes: an enrichment,
not a promoter".

### 135.23 — 725 proteinlik "recombinase-like HTH" ailesi: soru DARALTILDI, cevaplanmadi

124.1'in acik birakti sorusu: "bu ailenin hangi Pfam modeline karsilik
geldigi bulunmadi". Olculdu ve DARALTILDI; cevaplanmis gibi
yazilmiyor.

**OLCUM.** Urun metni "recombinase-like helix-turn-helix
domain-containing protein" olan 744 komsu satiri, 724 tekil protein,
**438 tekil dizi**. Dort kutuphanenin hepsine (recombinase 7, mobility
12, regulator 9, components 9 model) **E<=10** gibi cok gevsek bir
esikle arandi -- onceki olcum E<=1e-5 kullanmisti, yani zayif sinyal
gorulemezdi.

| kutuphane | isabet alan dizi / 438 | tutan model |
|---|---|---|
| recombinase | **50** | **hepsi `Recombinase` (PF07508)** |
| mobility | 0 | — |
| regulator | 0 | — |
| components | 0 | — |

**ZAYIF ISABET BIR ATAMA DEGIL, ve ne kadar zayif oldugu olculdu.** En
iyi E=1,2e-4, skor **17,1 bit**; 50'nin 26'si E<=1e-3, hepsi E<=0,01.
Pfam toplama esikleri genelde 20-25 bit, yani bu isabetlerin TAMAMI esik
ALTI. Dolayisiyla soylenebilen sey: 438 dizinin 50'sinde PF07508'e zayif
bir HOMOLOJI IPUCU var; bu "bu aile PF07508'dir" DEMEK DEGIL. Kalan 388
dizi 37 modelin hicbirine E=10'da bile benzemiyor.

**ILK HIPOTEZIM CURUK CIKTI, ve sinamak ucuzdu.** "recombinase-like
helix-turn-helix" ifadesi Pfam'daki HTH_7 (PF02796, "Helix-turn-helix
domain of resolvase") tanimina birebir uyuyor -- ama PF02796 ZATEN
`recombinase.hmm` icinde ve bu dizilerin hicbirini tutmuyor. Yani etiket
o modelden gelmiyor.

**UNIPROT YOLU KAPALI.** Bu komsu proteinlerin UniProt'ta lokus
etiketiyle (`BMA_RS22290`, `ACIAD_RS04545`) kaydi yok; dolayisiyla
135.18'de kurulan xref makinesi burada kullanilamiyor.

**HIPOTEZ (sinanmadi, oyle yazili).** Urun ifadesi NCBI PGAP tarafindan
uretiliyor ve PGAP ailelerini NCBIfam/TIGRFAM koleksiyonundan cagiriyor.
Eger etiketi tanimlayan model NCBIfam'da ise, **"hangi Pfam modeli"
sorusunun cevabi OLMAYABILIR**. Bu bir tahmin; dogrulanmadi.

**CEVABI NE VERIR.** 438 dizi uzerinde TAM bir Pfam-A taramasi ya da
InterProScan. Ikisi de yerelde yok, ve tek tek Pfam modeli indirip
denemek yontem degil tahmin olur (bu oturumda ayni aliskanligin bedeli
defalarca olculdu). Diziler `/media/lin-bio/back2/roar_taxonomy/hth725.fasta`
icinde duruyor, yani bir sonraki deneme olcumu tekrarlamadan baslar.

**Bu madde ACIK kaliyor** -- ama artik "hicbir sey bilinmiyor" degil:
aranan sey 37 modelin DISINDA, bir alt kumesi recombinase'e esik alti
benziyor, ve muhtemelen Pfam'da degil NCBIfam'da.

### 135.24 — "60 artakalan beta" bir anotasyon catismasi DEGIL, bozuk gen sayisiymis

Kapinin `check_component_scan_is_fresh` ikinci ayagi 60 proteinlik bir
artakalan butcesi tasiyor ve madde bunu soyle tarif ediyordu: "urun metni
beta diyor, hicbir model tutmuyor; **ikisinden biri yanlis** ve hangisi
oldugu HENUZ OLCULMEDI". Olculdu, ve tarif YANLISTI.

**ILK OLCUM SORUYU DEGISTIRDI.** 60'in TAMAMINDA `neighbor_protein`
satiri yok -- yani bu proteinler bilesen taramasinin aradigi kumede hic
bulunmadi. **Taranmamis bir dizi hakkinda hicbir model bir sey
soyleyemez**, dolayisiyla "hicbir model tutmuyor" cumlesi ne anotasyon ne
model kumesi hakkinda bir ifade. Dogru soru "hangisi yanlis" degil, "dizi
neden yok".

**KAYNAK CEVABI VERDI.** Her birinin GenBank kaydindaki CDS niteligi
okundu:

| hukum | protein |
|---|---|
| **psodogen** (`/pseudo`) | **58** |
| kayitta `/translation` yok | 2 |
| **dizi var, taranmis, model tutmamis** | **0** |

Ornek: `AB57_RS06925` (`benB`, benzoat 1,2-dioksijenaz beta alt birimi),
CDS'inde `/pseudo` ve `/note="frameshifted"`. Psodogen CDS'inin cevirisi
olmaz; o yuzden ne protein satiri olur ne bir HMM onu tutabilir.

**SONUC: ne anotasyon ne modeller yanlis.** Anotasyon beta genini DOGRU
adlandiriyor, gen BOZUK. Kapinin 60'lik butcesi bir catisma olcusu degil,
**bozuk gen** sayisi.

**KOORDINAT DONUSUMU OLCUMUN PARCASI.** Ilk aramada kayit BULUNAMADI,
cunku DB BioPython'un 0-tabanli yari-acik araligini saklarken GenBank
1-tabanli kapsayici yaziyor (`start+1`). Donusum yapilmadan "kayit yok"
sonucu cikiyor ve bu, sessizce yanlis cevap veren turden bir hata. Ayni
donusum GFF3 ihracinda da olculmustu (135.4); betige gerekcesiyle
yazildi.

**BUTCE ARTIK KESINLIK TASIYOR.** Yeni kontrol
`check_residual_beta_is_explained` butceyi degil **ACIKLANMAMIS KALANI**
siniyor: bir artakalan ya psodogen ya cevirisi olmayan olmak zorunda.
Yarin dizisi olan, taranan ve hicbir modele uymayan bir beta cikarsa
butce hala 60'in altinda olabilir ama BU kontrol FAIL verir -- cunku o
zaman gercekten "ikisinden biri yanlis" durumu olusur. Bir butce "kac
tane" der; bu kontrol "NEDEN" der.

Kapi metni de duzeltildi: "...classified as beta, **or its gene is
broken**", cunku eski metin okuyani anotasyon-model catismasi oldugunu
dusunmeye goturuyordu.

Cikti: `analysis_out/residual_beta_audit.json`, betik
`residual_beta_audit.py`.

### 135.25 — 135.24'un eki: kendi cikarimimdaki hata, replikon duzeyinde bir ceviri boslugunu gizliyordu

135.24'u yayimladiktan sonra yayimlanan dosyada bir sey gozume carpti:
`note` alani **60/60 BOS**, oysa ham kayitta `/note="frameshifted; ..."`
ifadesini gormustum.

**HATA BLOK OKUMADAYDI.** Ayni koordinatta IKI ozellik var: once `gene`,
sonra `CDS`. Betik ilkini topluyor ve ikincisine gecerken DURUYORDU.
`/pseudo` gene blogunda oldugu icin hukum DOGRUYDU, ama notu tasiyan CDS
blogu HIC OKUNMAMISTI. Bedeli: bos bir alan "frameshift notu yok" diye
okunur -- yani **aracin uretigi bir yokluk, veride yokluk gibi
gorunuyor**, bu oturumun tekrar eden hata bicimi. Duzeltildikten sonra
59/60 notlu, 58/60 frameshift/partial isaretli.

**DUZELTME GERCEK BIR BOSLUK ACIGA CIKARDI.** Hukum degisti: bir protein
"kayitta ceviri yok"tan "**cevirisi var ama cikarilmamis**"a gecti, ve
135.24'te yazdigim kontrol HEMEN FAIL verdi -- tam olarak yakalamak icin
yazdigi durum.

Vaka: `C884_RS03380` (`benB`, `NZ_ANHZ02000004.1:160827-161343`).
Tamamen normal bir CDS: `/translation` var, `/protein_id=WP_006213886.1`,
EC 1.14.12.10, psodogen DEGIL. Yine de `neighbor_protein`'de satiri yok.
Ve bakinca daha buyuk cikti: **o replikonda HIC `neighbor_protein` satiri
yok**.

**YAYGINLIK OLCULDU.**

| birim | sayi |
|---|---|
| komsu satiri olan replikon | 12.122 |
| bunlardan hic cevirisi olmayan | **1.472 (%12)** |
| etkilenen komsu satiri | 40.733 / 474.724 (%8,6) |
| o replikonlardaki **dogrulanmis giris** | **1** |

Son satir belirleyici: yayimlanan operon tamlik sayilari DOGRULANMIS
girisler uzerinden hesaplaniyor, yani etki 1/19.125. Cikarim dogrulanmis
girisi olan replikonlara odaklanmis; 40.733 satir onaylanmamis adaylara
ait ve hicbir yayimlanan sayi onlari kullanmiyor.

**KONTROL IKIYE AYRILDI, cunku iki kategori AYRI SEY demek.** Ilk yazimda
ikisini tek kovaya koyup ikisine de ayni mesaji verdim ve mesaj biri icin
YANLISTI:
  - `scanned_and_no_model_hit` -- dizi var, TARANMIS, model tutmamis. 
    Maddenin iddia ettigi gercek CATISMA. Olcum: **0**. FAIL kosulu.
  - `translation_exists_but_not_extracted` -- dizi kayitta var ama
    TARANMAMIS. Catisma degil, CIKARIM BOSLUGU. Olcum: 1. **WARN**,
    olculen olcegiyle (1.472 replikon, 1 giris etkisi) mesajin icinde.

Ikincisini FAIL yapmak kapiyi KALICI KIRMIZI birakirdi, ve kalici kirmizi
bir kapi okunmayi birakir (135.11'in dersi). Ama WARN'i gerekcesiz
birakmak da olmazdi; sayilar mesajda duruyor.

**ASIL KAZANC:** 135.24'te yazilan kontrol, yazildiktan DAKIKALAR sonra
hem kendi cikarimimdaki hatayi hem onun gizledigi replikon duzeyindeki
boslugu buldu. Bir BUTCE ("60'in altinda") bunlarin hicbirini
gostermezdi.

Kapi: 246 kontrol, 239 PASS, 0 FAIL, 7 WARN.

### 135.26 — Grup 3'te "%72,2 beta yoklugu": hipotezlerden hicbiri degil, BIRIM yanlisti

Madde uc hipotez sayip "olculmedi" diyordu: (a) grup atamasi genis, (b)
beta 150 bp esiginin otesinde, (c) anotasyonda yok. Olculdu: **ucu de
degil.** Sayinin kendisi dogru, ama BIRIMI yanlis -- "grup basina beta
yoklugu" farkli kimyalari tek kovaya koyuyor.

**(b) CURUDU, hemen.** `operon.nearby_beta` zaten olculuyor: beta'si
olmayan 2.390 grup-3 girisinin yalnizca **6'sinda** beta yakinda ama
operon disinda. Yani beta esigin otesinde degil, komsulukta HIC YOK.

**TIP KIRILIMI SORUYU DAGITTI.** Beta'siz 2.390 girisin **2.190'i (%92)**
uc tipten geliyor ve ucu de %100 beta'siz:

| tip | beta yok / toplam | ne yapiyor |
|---|---|---|
| `3_313_KshA15` | 1.598 / 1.599 | 3-ketosteroid 9-alfa-hidroksilasyonu |
| `3_322_DAF36` | 319 / 319 | kolesterol C7-C8 desaturasyonu |
| `3_321_Nvd` | 273 / 273 | kolesterol C7-C8 desaturasyonu |

Ucu de kuratorlu kaydinda `family: terpenoids_steroids`, ve kimyalari
HIDROKSILASYON / DESATURASYON -- halka cis-dihidroksilasyonu DEGIL.
Beta alt birimi cis-dihidroksilleyen dioksijenazin ortagidir; bu
enzimlerde yoklugu beklenen biyoloji, kusur degil. (Bu gece yapilan
InterPro olcumu de ayni yone isaret etti: `IPR050584` "Cholesterol
7-desaturase" 11.837 dizide.)

**REAKSIYON SINIFINA GORE OLCUM -- yayimlanacak olan bu:**

| reaksiyon sinifi | beta yok | pay |
|---|---|---|
| desaturation | 592 / 592 | **%100** |
| hydroxylation | 1.610 / 1.931 | %83,4 |
| **cis_dihydroxylation** | **188 / 758** | **%24,8** |
| dioxygenation_with_release | 0 / 26 | %0 |
| angular_dioxygenation | 0 / 2 | %0 |

Yani maddenin ONERMESI DOGRU ("naftalin tipi dioksijenazda beta zorunlu
ortak"), ve o onermenin gecerli oldugu kumede yokluk **%24,8**, %72,2
degil. `3_307_NagGH` ayrica gosteriyor ki "hidroksilasyon" da duzgun
beta'siz degil: 332 girisinin 320'sinde beta VAR (iki bilesenli salisilat
5-hidroksilaz).

**(a) KISMEN DOGRU, ve olculdu.** Beta'siz girisler referansa medyan
%45,2 kimlik tasiyor, beta'lilar %68,4 -- birkac ayri aileyi kapsayan bir
grupta beklenen fark.

**KALAN GERCEK SORU, bir buyukluk mertebesi kucuk.** cis-dihidroksilleyen
tiplerde %24,8, ve iki tipte yogunlasiyor: `3_316_NarAa` 38/48 (%79),
`3_319_PhdA` 74/140 (%53). Ayrica beta'siz 2.390 girisin **206'sinda hic
komsu yok** -- olcum yoklugu, biyolojik yokluk degil (135.10'un sinifi).

**YAYIMLANAN SAYI YOK, ve kayda gecen kural bu:** %72,2 sitede hicbir
yerde basilmiyordu, yani duzeltilecek bir yayin yok. Ama **grup basina
beta yoklugu bir daha YAYIMLANMAMALI**: dogru birim reaksiyon sinifi (ya
da tip). Bir grup, kimyasi farkli enzimleri tasiyorsa, grup basina bir
ortak-varligi orani iki ayri olguyu tek sayida toplar.

### 135.27 — 135.26'nin kalani: beta yoklugu KANIT UZAKLIGINI izliyor, ve iki okuma da acik birakiliyor

135.26 artik soyle tanimli bir soru birakmisti: cis-dihidroksilleyen
tiplerde %24,8 beta yoklugu, `3_316_NarAa` (38/48) ve `3_319_PhdA`
(74/140) uzerinde yogun. Uc aday sinandi.

**PARCALILIK CURUDU.** Beta'siz girisler kisa contig'lerde DEGIL -- tam
tersi, PhdA'da beta'siz girisler medyan 1.864 kb replikonda, beta'lilar
1.738 kb'de. Komsu sayilari da ayni (19 vs 21). Beta baska contig'de
kalmis degil.

**ANOTASYON BOSLUGU CURUDU.** Beta'siz NarAa/PhdA girislerinin
komsulugunda urun metni "small/beta subunit" diyen **sifir** gen var. Beta
boyunda (150-200 aa) ANOTASYONSUZ bir komsu da yok denecek kadar az:
NarAa 3/38, PhdA 1/74.

**KALAN: KANIT UZAKLIGI -- ve iliski tek yonlu.** cis-dihidroksilleyen
TUM tiplerde (3.487 giris), kanit kademesine gore:

| kademe | beta yok | pay |
|---|---|---|
| characterized (>=%95) | 9 / 224 | **%4,0** |
| close_homolog (>=%60) | 177 / 1.266 | %14,0 |
| family_member (>=%40) | 136 / 866 | %15,7 |
| **distant (>=%25)** | **861 / 1.117** | **%77,1** |
| novel | 11 / 14 | %78,6 |

NarAa ve PhdA'nin beta'siz girisleri referanslarina medyan **%32,1** ve
**%33,2** kimlikte; beta'lilar %44,0 ve %38,1.

**OLCUM EKSIKLIGI DEGIL -- kontrol edildi.** Komsuluk zenginligi
kademeler arasinda AYNI: her kademede medyan 20-21 komsu, komsusu hic
olmayan 1.117 uzak giriste yalnizca 1. Yani uzak girislerin komsulugu da
ayni derecede iyi okunmus; beta orada yok cunku YOK, okunamadigi icin
degil.

**IKI OKUMA, VE VERI IKISINI AYIRMIYOR.** Bu yuzden birini secmiyorum:
  (a) UZAK Rieske oksijenazlari gercekten sik sik beta ortagi olmadan
      calisiyor olabilir;
  (b) %25-40 kimlikte tip etiketi bir KIMLIK TESPITI degil, EN YAKIN
      REFERANS ifadesidir -- yani "bu bir naftalin dioksijenazi, o halde
      beta'si olmali" cikarimi o girisler icin GECERSIZ olabilir.
Ikisi de ayni veriyle uyumlu ve ayirmak icin burada hicbir sey yok.

Bu, bu veritabaninin kendi kaydettigi uyarinin yeni bir alana
uygulanmasi: "substrat etiketi REFERANS enzime aittir, uyesine degil".
Ayni sey ORTAK GEREKSINIMI icin de gecerli.

**KAYDA GECEN KURAL.** Beta yoklugu (ya da herhangi bir ortak-varligi
orani) yayimlanacaksa KANIT KADEMESINE gore bolunmek zorunda -- tipki
reaksiyon sinifina gore bolunmesi gerektigi gibi (135.26). Tek bir oran,
"%4" ile "%77"yi ayni sayida toplar ve okuyucu onu biyolojik bir olgu
sanir. Su an hicbir yerde yayimlanmiyor.

### 135.28 — Kullanicinin qxyA sorusu: kimlik hesabi saglam, INGEST acik -- 12.312 kayit

Kullanici `5_505_qxyA` icin "referansa daha yakin homolog cikmasi lazim,
genel bir bug mi var" diye sordu, ve "NCBI'da blastlayinca %95 ustu
buluyordum, nukleotid dizileri de vardi" dedi. Uc katman olculdu; ikisi
temiz, ucuncusu gercek bir acik.

**1. KIMLIK HESABI SAGLAM.** qxyA referansi 15.960 dogrulanmis dizinin
tamamina karsi dogrudan blastlandi: kendisi haric en iyi isabet **%45,5**
-- `ro_evidence`in dedigiyle birebir ayni. Kademe hesabinda kusur YOK.

**2. YAKIN HOMOLOGLAR VERIDE VAR, AMA UYE SAYILMIYOR.** Birlesik kumede
**%97,4** (`A0A1H4ZH95`, *Pseudomonas saponiphila*) ve **%97,1**
(`A0AAP4RC76`) kimlikli iki dizi var. Ikisi de kapiyi geciyor, Rieske
saglam, ucgen 3/3, ve `ro_cluster` zaten `5_505_qxyA`. Gorunmemelerinin
sebebi `has_genomic_context=0`: uye sayilari ve kademe dagilimi yalnizca
genomik `ro` tablosundan hesaplaniyor. qxyA'da 164 uye gorunurken
birlesimde 1.293 dizi var (1.145'i baglamsiz); veritabani genelinde
kapiyi gecen + tipe atanmis + baglamsiz **55.859** dizi.

Bu mimari bir karar ve tip sayfasi birlesim sayisini ZATEN acikliyor
("N further sequences are not counted here"). Ama **kanit kademesi
dagilimi aciklamiyor**, ve bu gece 135.27'de tam olarak o tabloya
dayanmistim -- "qxyA'nin yakin homologu yok" cumlesi kademelerden
okunursa YANLIS.

**3. GERCEK ACIK: GENOMLAR INDIRILMEMIS.** %97,4'luk homolog
`FNTJ01000003` contig'inde (FNTJ01 derlemesi). Bizde o turden tek kayit
var ve BASKA bir derlemeden: `NZ_JAOXCB010000053.1`, 37,5 kb, 32 CDS,
uzerinde qxyA degil bir `5_508_GbcA`. Yani ingest CONTIG BAZLI -- ilk
aramada isabet veren contig alinmis, kardesleri ve alternatif derlemeler
alinmamis.

**OLCEK (`ingest_gap.py`).** Baglamsiz + kapiyi gecen + tipe atanmis
dizilerin UniProt EMBL capraz referanslari okundu:

| | |
|---|---|
| EMBL nukleotid xrefi olan dizi | 17.930 |
| xref kaydi bizde olan | 501 |
| **bizde OLMAYAN ayri kayit** | **12.312** |
| **baglam kazanacak dizi** | **17.528** |

Tip basina kazanc (ilk 5): `5_504_CntA` +2.463, `5_506_Stc2` +1.587,
`5_507_BmoA` +1.194, `3_313_KshA15` +981, `1_106_TsaM` +957.
**`5_505_qxyA` +258** -- kullanici "en az bes tane" demisti.

Liste TAHMIN DEGIL: her satirin arkasinda, accession'i bilinen bir
kayitta duran, kapimizi gecmis ve tiplenmis en az bir dizi var. En
verimli satir `CP000699`: tek kayit, 23 tipten 51 dizi.

**INDIRILMEDI.** 12.312 kayit almak yayimlanan her sayiyi degistirir
(yeni veri surumu) ve bu kuratorun karari. Liste
`analysis_out/ingest_gap.csv`.

**ACIK IS.** Kanit kademesi dagilimini birlesik kume icin de hesaplamak
ya da en azindan kademe tablosunun hangi paydayi kapsadigini sayfada
yazmak -- yoksa 135.27'deki gibi bir okuma yine yapilir.

## 135.29 · qxyA küçük ölçekli yutma testi ve mekanizmanın tavanı

Soru kullanıcıdan geldi: `5_505_qxyA` referansına karşı daha fazla yakın
homolog çıkması gerekiyordu, NCBI'da blastlayınca %95 üstü isabetler
görülüyordu, bizde `close_homolog` kademesi boştu. "Genel bir bug mu var?"

Ölçüldü, üç katman çıktı ve hiçbiri hesap hatası değil.

**1. Kimlik matematiği doğru.** Referansın kendisiyle doğrudan blast'ı
%45,5 veriyor ve `ro_evidence` de %45,5 diyor. Kademe eşikleri
(`evidence_tiers.py`: characterized >=95, close_homolog >=60) doğru
uygulanıyor.

**2. Yakın diziler veritabanında var, ama kademe tablosunun
denominatöründe değil.** `ro_union`'da qxyA'ya atanmış, kapıyı geçen
1.126 dizi bağlamsız (`has_genomic_context=0`). Üye sayıları ve
`ro_evidence` kademeleri yalnızca genomik tablodan hesaplanıyor. Yani
rakamı üreten şey bir hata değil, bir denominatör seçimi — ama sayfa
hangi denominatörü gösterdiğini söylemiyordu, yanlış okuma oradan
doğdu.

qxyA'nın en yakın dört birlik dizisi:

| md5 | kimlik | bağlam | EMBL nükleotid | durum |
|-----|--------|--------|----------------|-------|
| ac69f52b | %100 | var | — | referansın genomik kopyası |
| 85204e5d | %97,4 | yok | FNTJ01000003 | **indirildi, bağlam kazanıyor** |
| 59c35d2c | %97,1 | yok | yok | adlandırılamıyor, bu yolla kurtarılamaz |
| b113e604 | %48,0 | yok | yok | adlandırılamıyor |

**3. Gerçek boşluk: yutma isabet-başına-kontig.** Aynı organizmanın
kardeş kontigleri ve alternatif derlemeleri hiç indirilmedi. Kendi çapraz
referanslarımızdan türetilen liste (`ingest_gap.py`) bizde olmayan
**12.312** nükleotid kaydı buluyor; indirilirse **17.528** dizi bağlam
kazanıyor.

**Mekanizmanın tavanı ölçüldü.** Bağlamsız + kapıyı geçen + tipi atanmış
55.859 dizinin yalnızca **17.930'u (%32)** adlandırılmış bir nükleotid
kaydına sahip. Kalan 37.929 (%68) çapraz referanslarımızda hiçbir
nükleotid erişim numarası taşımıyor; bunlar protein-yalnız kayıtlar ve
bu yolla bağlam kazanamaz. qxyA'da: 1.126 bağlamsızın 273'ü
adlandırılabilir.

Yani "tüm 12.312'yi indir" bile birliğin yalnızca üçte birini
bağlamlandırır. Bu, sayfada söylenmesi gereken bir sınır.

**Küçük test (kullanıcının seçimi).** qxyA'ya bağlam kazandıran 251 kayıt
`fetch_ingest_gap.py --type 5_505_qxyA` ile indirildi (251 inen, 0 düşen).
İndirme öncesi öncül doğrulandı: FNTJ01000003 397.710 bp, 462 CDS,
`ORIGIN` bloğu var, referansa %97,4 benzeyen 383 aa'lık çeviri bir tam
toluat dioksijenaz operonunun yanında duruyor. İzole ön test (canlı
`genomic_context/` ve veritabanına dokunmadan, `newgbk`/`newctx`):
1.039 RO adayı, 25.111 komşu kaydı; referansın diamond'ı aday
fasta'sına karşı en iyi isabet olarak `FNTJ01000003.1:78644-79796:-1`
**%97,4 kimlik, %100 kapsama** döndürdü, %90 üstü tam 1 isabet.

**Önce sayıları (ölçüm sabitlendi).** qxyA 164 üye: 1 characterized,
104 distant, 52 family_member, 7 novel, **close_homolog 0**.
`ro` onaylı toplam 19.125, replicon 17.073.

Açık iş: tam zincir koşuyor (`genomic_context` -> `database` ->
`annotate` -> `operons` -> alt adımlar -> kapı), sonra sonra-sayıları
ölçülüp buraya yazılacak ve 500'lük partiye / tamamına geçme kararı
verilecek.

Açık iş: kademe tablosunun hangi denominatörü kapsadığı sayfada
yazılmalı; 135.27 yanlış okuması tam olarak bu eksikten doğdu.
