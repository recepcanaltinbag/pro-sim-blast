# ROAR-DB — sonraki adimlar

Durum: `main` dalinda pipeline + web uygulamasi + Atlas yayinda,
`gh-pages` dalindan statik site: https://recepcanaltinbag.github.io/pro-sim-blast/

## Sirada

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

## Yeniden uretim

```
bash run_all.sh 16                  # tum pipeline (uzun)
cd webapp && uvicorn app:app        # yerel sunucu
python3 freeze.py --out site --base /pro-sim-blast && bash deploy_pages.sh \
  git@github-prosimblast:recepcanaltinbag/pro-sim-blast.git
```
