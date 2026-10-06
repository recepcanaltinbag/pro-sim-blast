# ROAR-DB — sonraki adimlar

Durum: `main` dalinda pipeline + web uygulamasi + Atlas yayinda,
`gh-pages` dalindan statik site: https://recepcanaltinbag.github.io/pro-sim-blast/

## Sirada

Acik isler. Numaralar kalicidir, commit mesajlari onlara atif yapiyor.

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

## Bitti

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

- cadA'nin substrati 2,4-D olarak isaretlendi ama KESIN DEGIL; Bradyrhizobium
  HW13 cadABC makalesi hem 2,4-D hem 2,4,5-T aktivitesinden soz ediyor.
  Elindeki orijinal referans hangisiyse soyler misin?
- Varyant kumelemesi ham dizide tekrarlanacak mi (madde 1)? Site genelinde
  sayilar degisir.
- Site Ingilizce; Turkce surum de istenir mi?

## Yeniden uretim

```
bash run_all.sh 16                  # tum pipeline (uzun)
cd webapp && uvicorn app:app        # yerel sunucu
python3 freeze.py --out site --base /pro-sim-blast && bash deploy_pages.sh \
  git@github-prosimblast:recepcanaltinbag/pro-sim-blast.git
```
