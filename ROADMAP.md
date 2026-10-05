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

## Acik sorular (kullaniciya)

- 7 kume substrati bilinmiyor: OxoO, CndA, cadA, OMO, PsbAb, ROCH34, OxyA.
  Orijinal makaleler ya da enzim/organizma adlari gerekiyor.
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
