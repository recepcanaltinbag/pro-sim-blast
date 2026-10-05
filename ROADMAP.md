# ROAR-DB — sonraki adimlar

Durum: `main` dalinda pipeline + web uygulamasi + Atlas yayinda,
`gh-pages` dalindan statik site: https://recepcanaltinbag.github.io/pro-sim-blast/

## Sirada (oncelik sirasina gore)

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

5. **Hugging Face Space** (tam uygulama + dizi siniflandirici). `webapp/Dockerfile`
   hazir ve test edildi; kullanicinin HF hesabi gerekiyor.

## Sirada

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

## Sirada

10. **Per-tip substrat kuratorlugu.** 7 tip hala `bilinmiyor` (OxoO, CndA, cadA,
    OMO, PsbAb, ROCH34, OxyA) — orijinal makaleler gerekiyor, kullanicidan
    bekleniyor.
11. **Hugging Face Space.** Docker imaji hazir ve test edildi; kullanicinin
    hesabi gerekiyor.
12. ~~Veri kalitesi sayfasi.~~ **BITTI.** `/atlas/quality`: dizi fazlaligi
    (%12,3) ve her sayima etkisi, uc duzeyli test tablosu, olculemeyen
    sinirlarin listesi, ve her indirilebilir dosyanin hangi scriptten geldigi.
    `redundancy.py` + `stats_overview.py`'ye "sequence" duzeyi eklendi.

13. **GitHub Pages derleme gecikmesi (COZULUYOR).** Site 26.683 dosya / 371 MB'a
    cikinca Pages derlemesi push'un bir saat gerisinde kaldi. Iki duzeltme
    yapildi: (a) `deploy_pages.sh` artik artimli (kalici klon + rsync + normal
    push), eskiden her yayinda sifirdan repo kurup tum agaci force-push
    ediyordu; (b) giris ve varyant basina FASTA dosyalari statik siteden
    cikarildi (-13.200 dosya, ~-52 MB) cunku dizi sayfada zaten var ve toplu
    dosyalar hepsini kapsiyor. Hala geride kalirsa sonraki adim: giris
    sayfalarini 11.422 HTML yerine tip basina JSON + tarayicida render etmek
    (dosya sayisi ~2.000'e duser).

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

## Sirada

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

## Sirada

24. **Referans setindeki tekrarlar (KULLANICI KARARI).** 71 kayit, 68 tekil dizi.
    OxoO=OMO, NahAc=NDO(3_315), NDO(3_314)=NarAa birebir ayni; NdmC, NdmB'nin
    alt dizisi. HMM ayni profili iki kez icerdigi icin atama ikizler arasinda
    keyfi bolunuyor. Hangi ikiz tutulacak?
25. **Habitat sozlugunun kalan zayif noktalari.** Ajan sekiz tanesini isaretledi;
    en onemlileri: `hospital` anahtari lavabo/yuzey orneklerini de klinige
    sokuyor, `lymph node` fare deneylerini insan kliniğine sokuyor, ve
    `rhizosphere_plant` rizosfer/endofit/yaprak/colemen hepsini tek kovada
    tutuyor. Insan ve hayvan dokusu hic ayrilmiyor.

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
