# ROAR-DB — yeniden kurulan pipeline

Rieske oksijenaz (RO) alpha alt birimlerini PF00355 Rieske proteinlerinden ayıran,
genomik bağlamlarını GenBank'tan çıkaran, yerel bir veritabanı kuran ve varyant/novel
düzeyine kadar analiz eden pipeline.

Bu, eski MongoDB tabanlı akışın (`mongo_*.py`, `filter_hmm_out.py`) yerine geçer.
Eski kodun neden değiştiği ve bulunan hatalar için interaktif rapora bakın.

## Neden yeniden kuruldu

PF00355 ile çekilen protein seti yalnızca RO alpha değil; ferredoksin, sitokrom
*bc*₁ Rieske ISP, NirD gibi **Rieske merkezi taşıyan ama katalitik merkezi olmayan**
proteinleri de içerir. Ölçüldü: **E-value bu ayrımı yapamıyor** (her iki grup da
%100 geçiyor). Ayırt edici olan Rieske değil, ona eşlik eden **mononükleer Fe(II)
katalitik merkezi**dir.

## Kurulum

```
# Python 3.8+ ve: biopython, pandas, numpy
# Harici: hmmer (hmmsearch/hmmalign/hmmbuild/hmmpress), cd-hit
pip install biopython pandas numpy
# Web uygulaması için ayrıca: pip install -r webapp/requirements.txt
```

Veri: `combined_pfam.fasta` (189.657 PF00355 proteini), `gbk_files/` (17.074 GenBank),
`ROs_71_Clean/` (71 küratörlü referans), `RieskeDB71.hmm`.

## Pipeline (sırayla)

| # | Script | Girdi → Çıktı |
|---|--------|---------------|
| 1 | `run_ro_filter.py` | `combined_pfam.fasta` → `ro_alpha.csv/.fasta` — coverage ≥0,45 + katalitik motif |
| 2 | `extract_genomic_context.py` | `gbk_files/` → `genomic_context/` — RO + komşular (strand'li) |
| 3 | `build_db.py` | `genomic_context/` → `roar.sqlite` — replicon/ro/neighbor/gene_category |
| 4 | `annotate_ro.py` | hmmsearch + hmmalign → `ro` tablosu doğrulanır (11.422 RO) |
| 4b | `build_operons.py` | komşu protein dizileri (gbk) + Pfam HMM (`component_hmms/`) → `neighbor_protein`, `neighbor_component`, `operon`, `operon_gene`; `ro.sequence` |
| 4c | `evidence_tiers.py` | diamond blastp → `ro_evidence` (kanit kademesi: characterized / close_homolog / family_member / distant / novel) |
| 4d | `etc_types.py` | komsu Pfam domainleri → `ro_etc` (FNR/GR reduktaz, Rieske/bitki ferredoksin, α3 vs α3β3) |
| 4e | `analyze_regulation.py` | operon 5' ucu + yukari intergenik bolge + regulator ailesi → `ro_regulation` |
| 5 | `analyze.py` | orientation, divergent regülatör, transpozon×takson |
| 6 | `null_model.py` | komşuluk zenginleşmesi — rastgele pencere kontrolü |
| 7 | `analyze_variants.py` | küme içi kimlik, alt-aile, SDP pozisyonları |
| 8 | `discover_subfamilies.py` | CD-HIT alt-aileler + mutlak ölçüt + novel adaylar |
| 9 | `recursive_homogenize.py` | özyinelemeli bölme → homojen yapraklar (`leaf` tablosu) |
| 10 | `characterize_leaves.py` | yaprak = varyant, komşuluk imzasıyla (`leaf_profile`) |
| 11 | `extract_representatives.py` | her yaprak/küme için temsilci dizi + katalitik doğrulama |
| 12 | `classify_domains.py` | yaşam alanı (bakteri/ökaryot/arke) — `ro_domain` |
| 13 | `analyze_ecology.py` | substrat sınıfı × mobilite hipotez testi |
| 12b | `variant_signature.py` | varyantlari ayiran hizalama kolonlari → `cluster_sdp`, `leaf_sdp` |
| 12c | `motif_stats.py` | 8 tanimlayici kolonun korunma olcumu → `motif_stats.json` |
| 12c2 | `operon_validation.py` | operon kuralinin sinanmasi + esik duyarliligi → `operon_validation.json` |
| 12c3 | `redundancy.py` | dizi fazlaligi ve her sayima etkisi → `redundancy.json` |
| 12c4 | `cooccurrence.py` | tip birliktelikleri + permutasyon null'i → `cooccurrence.json` |
| 12c5 | `substrate_predictability.py` | kimlik → substrat ongorusunun ROC/kesinlik olcumu → `substrate_predictability.json` |
| 12c6 | `isolation_source.py` | `gbk_files/` `source` nitelikleri + `cluster_ecology.csv`/`chemistry.csv` → `replicon_source` tablosu (izolasyon kaynagi, konak, cografya, yil, habitat) + `analysis_out/habitat.json`; habitat × substrat sinifi tur-normalize capraz tablolari |
| 12d | `build_search_index.py` | FTS5 tam metin indeksi + filtre alanlari → `ro_search`, `ro_fts` |
| 13b | `build_phylogeny.py` | hmmalign + FastTree → `tree_all.nwk`, `trees/<kume>.nwk`; diamond all-vs-all → `ssn_*.csv`, `cluster_identity_matrix.csv` |
| 13c | `stats_overview.py` | hipotez testleri (entry + genus duzeyi, etki buyuklugu) → `stats.json` |
| 13d | `validate_curation.py` | **KAPI**: DB + iki kuratorlu CSV denetimi → PASS/FAIL satirlari; bir FAIL varsa sifir-disi cikar ve pipeline durur (bkz. Safeguards) |
| 13e | `provenance.py` | her tablo/dosya icin ureten script, girdiler, kayit sayisi, zaman, kaynak zinciri → `analysis_out/provenance.json` |
| 14 | `make_report.py` / `make_explorer.py` / `make_hub.py` | üç HTML sayfa |

Tümünü koşmak için: `bash run_all.sh` (bkz. aşağı).

## Filtre mantığı (aşama 1)

Üç katman, her biri bağımsız:

1. **Kapı** — `ROmotif71.hmm` (71 küratörlü referanstan), HMM model kaplaması ≥0,45.
   Kaplama 2.788 tam-boy alpha / 1.570 kesin non-alpha üzerinde kalibre edildi.
2. **Doğrulama** — `ro_motif.py`, hizalama kolonlarında doğrudan kalıntı testi:
   Rieske ligandları 4/4 + katalitik triad ≥2/3. Kolonlar referans hizalamasından çıkarıldı.
3. **Atama** — `RieskeDB71.hmm` 71 model, en yüksek **bit skoru** (E-value değil).

Sonuç: 189.657 → **73.272 RO alpha** (eski yöntem 117.537 kabul ediyordu; ~44.000
kontaminant temizlendi).

## Temel bulgular

- **Referans setinde 2 kontaminant**: `2_205_IsoMO` (Rieske merkezi yok), `1_113_CdnD`
  (elektron transfer bileşeni) — çıkarıldı, modeller 71 referansla yeniden kuruldu.
- **Null model**: regülatör komşuluğu RO'ya özel değil (1,08x, arka plan); mevcut
  co-occurrence skorlaması genom bolluğunu ölçüyor. Gerçek zenginleşenler: ro_beta
  (4,12x), ring_cleavage (2,81x), ferredoksin (2,05x).
- **Kümeler heterojen**: doğrulanan RO'ların ~%67'si medyan kimliği <%35 olan
  kümelerde. Özyinelemeli homojenizasyon 61 kümeyi 1.809 yapraga böldü.
- **Varyantlar genomik bağlamla ayrışıyor** (KshA'nın bir varyantı steroid dehidrogenaz
  yanında; VanA'nın mobil *Burkholderia* soyu IclR+transpozon yanında).
- **Novel**: 318 yuksek-guven aday (güvenle RO, bilinen tipe uymuyor) + bitki/alg dalı.
- **Ekoloji**: ksenobiyotik RO'lar plazmit uzerinde 4,0x daha sik (tur bazinda 3,3x, cins bazinda 4,4x) (yatay transfer);
  transpozon farkı sınırda (z=2,57); taksonomik yayılım ters (doğal substratlar daha
  çok soyda — kadim dikey kalıtım).
- **Veri kompozisyonu**: %91,2 bakteri, %8,3 ökaryot (483 bitki/alg, 219 mantar,
  201 hayvan). Katalitik filtre yaşam alanını ayırmaz.

## Bulunan hatalar (eski kodda)

1. `mongo_gbk_analysis.py:174` — gen/CDS çift sayımı (prokaryotta ikisi aynı gen)
2. `mongo_gbk_analysis.py:175` — strand hiç okunmuyor (orientation analizi imkânsızdı)
3. `mongo_gbk_analysis.py:175` — dairesel orijini aşan genler tüm replikonu kaplıyor
4. `mongo_gbk_analysis.py` — anotasyonsuz dosya (%17,4) sessizce boş dönüyor
5. `mongo_gbk_analysis.py:186` — keyword filtresi çıkarımda (`hypothetical` komşular kaybolur)
6. `filter_hmm_out.py:41` — küme ataması E-value ile (bit skoru daha kararlı)

## Bulunan ve düzeltilen hatalar (yeni pipeline, 2026-10-05 incelemesi)

1. **Divergent regülatör sayımı RO ipliğini yoksayıyordu** (`analyze.py`,
   `make_report.py`): `distance<0` genom koordinatında "sol" demek; eksi iplikli
   RO'da kafa-kafaya geometri sağda. 960 → **1.782** (ve 92 konvergent yanlış sayılmıştı).
2. **`ro_beta` regex sırası**: "ring-hydroxylating dioxygenase subunit beta" önce
   `ro_alpha`'ya düşüyordu (1.529 yanlış kategori). Sıra değişti; null model
   ro_beta zenginleşmesi 4,1× → 6,1×.
3. **Null model paydası**: kontrol pencereleri ≥5 CDS'li replikonlardan örneklenirken
   payda tüm RO'ları sayıyordu (~%9 aşağı çekme). Payda aynı popülasyona çekildi.
4. **`analyze_ecology.py` operon paydası** INNER JOIN yüzünden komşusuz RO'ları düşürüyordu → LEFT JOIN.
5. **`cluster_ecology.csv`** `2\,4` kaçışı CSV'de geçersizdi; TftA/DntAc satırları kaymıştı. Alıntılandı + yükleyici doğrulama yapıyor.
6. **`make_report.py --db` yoksayılıyordu** (novel/domain bölümleri sabit `roar.sqlite` açıyordu); boş CSV'de çöküyordu; sabit yazılmış 960 / "%95+ kimlik" / "4,3×" metinleri DB ile çelişiyordu.
7. **`annotate_ro.py`**: hmmsearch hiti olmayan 945 aday `is_confirmed=NULL` kalıyordu → 0.
8. **`extract_genomic_context.py`**: çok kayıtlı gbk'da başka kontigin genleri komşu sayılabilirdi (kontig kontrolü eklendi); dairesel replikonda `gene_offset` sarılmıyordu (90 satır); artık komşu dizileri de yazıyor.
9. Explorer: null kimlikte kararsız sıralama, CSV'de tırnak kaçışı yok, DB metinleri HTML'e kaçışsız basılıyordu. Hub: boş DB'de "None".
10. `run_all.sh`: gereksiz `clustalo` ön koşulu; eksik modelde uyarıp devam ediyordu (artık çıkar).

## Operonlar (adım 4b, `build_operons.py`)

Komşu CDS'lerin protein dizileri gbk'dan çekilir (`neighbor_protein`), Pfam
modelleriyle taranır (`component_hmms/components.hmm`: Ring_hydroxyl_B, Fer2,
Rieske, Rieske_2, NAD_binding_1, FAD_binding_6, FAD_binding_8, Pyr_redox_2,
Pyr_redox_dim) ve beta / ferredoksin (≤250 aa) / redüktaz / alpha_other etiketi
verilir. Operon = RO'dan iki yöne aynı iplikte, genler arası boşluk ≤150 bp
(`--max-gap`). `operon` tablosu operon içi ve ±10 kb içi bileşen varlığını ayrı
tutar. Sonuç (11.422 RO): operonda beta %21,9, ferredoksin %10,5, redüktaz %28,0,
üçü birden %4,4. Beta yalnızca beklenen tiplerde (TdnA1, XylX, NBDO…) çıkıyor;
KshA/CntA/VanA'da yok — biyolojiyle tutarlı.

## Yeni referans ekleme ve modulerlik

Referans seti kapali bir liste degil; literaturde karakterize edilmis ve burada
olmayan RO'lar eklenebilir. Eskiden bu, bes dosyayi elle tutarli tutmayi ve
**motif kolon numaralarini elle yeniden cikarmayi** gerektiriyordu: kolonlar
`ro_motif.py` icinde SABIT yaziliydi, referans seti degisince hizalama kayar ve
kod sessizce yanlis pozisyonlara bakardi.

Simdi iki script var:

- **`add_reference.py`** — tek komutla ekler, once dogrular: ID bicimi ve
  cakismasi, dizi karakterleri ve uzunlugu, ayni dizinin zaten var olup olmadigi,
  SMILES'in ayraç dengesi, PDB kimliginin bicimi, ve **kaynak alaninin bos
  olmamasi** (bu veri nereden geliyor sorusu cevapsiz birakilamaz). Yazmadan
  once `.bak-<tarih>` yedegi alir, sonra uc dosyayi birlikte guncelleyip
  sonraki adimlari ekrana yazar. `--dry-run` ve `--list-vocabularies` var.
- **`build_reference_models.py`** — hizalama + hmmbuild + hmmalign zincirini
  kosar ve **motif kolonlarini VERIDEN cikarir**: kolonlari mutlak pozisyonla
  degil bagil yapiyla arar (Rieske C-x-H...C-x-x-H, katalitik H..H..D/E,
  aralik kisitlari modul basinda sabit olarak tanimli) ve korunmayi en yuksek
  yapan kombinasyonu secer. Sonuc `ROs_71_Clean/motif_columns.json`'a yazilir;
  `ro_motif.py` bu dosya varsa onu okur, yoksa gomulu varsayilanlari kullanir.
  **Dogrulama:** mevcut 71-referans seti uzerinde cikarim, koda elle yazilmis
  sekiz kolonun tamamini birebir yeniden uretti (85/87/105/108, 212/217/355,
  kopru 209), yani hem kod hem algoritma karsilikli dogrulanmis oldu.

Istatistik adimlarinin hepsi bagimsiz kosulabilir ve hepsi `--db` alir; referans
seti degismediyse yalnizca istatistikleri tazelemek icin pipeline'in basina
donmek gerekmez.

## Kanit duzeyi (adim 4c) — en onemli metodolojik nokta

HMM atamasi bir proteini **en yakin kuratorlu referansa** koyar; bu fonksiyon
atamasi DEGILDIR. `evidence_tiers.py` her uyenin 71 referans proteine diamond
kimligini olcer ve kademelendirir (varsayilan esikler, `--cut-*` ile degisir):

| kademe | kimlik | ne soylenebilir |
|---|---|---|
| characterized | ≥95% (+%90 kaplama) | referans enzimin kendisi / sus varyanti |
| close_homolog | ≥60% | ayni reaksiyon cok olasi |
| family_member | 40–60% | ayni tip, substrat belirsiz |
| distant | 25–40% | RO alpha, tip yalnizca homolojiyle atandi |
| novel | <25% | yakin referans yok, yeni tip adayi |

Olculen dagilim: characterized %1,9 · close_homolog %23,2 · family_member %16,7 ·
distant %54,9 · novel %3,3. Yani substrat etiketi uyelerin yalnizca **%25'ine**
aktarilabilir. 16 kume ≥%80 distant/novel uyeden olusuyor. Uyelerin %29,3'unde
en yuksek skorlu profil ile en yakin referans protein FARKLI tipe ait.

**Esik kalibrasyonu** (`analysis_out/reference_pairs.csv`): 71 referansin kendi
aralarinda farkli substratli ciftler %99,8 kimlige kadar cikiyor (EdoA1 etilbenzen /
CumA1 kumen; TDO toluen / BedC1 benzen %92; NarAa naftalen / NidA piren %91), ayni
substratli ciftlerin en dusugu %32,4. Sonuc: **hicbir global kimlik esigi "ayni
substrat" garantisi veremez**; substrat secimi birkac aktif-bolge kalintisiyla
belirlenir. Kademeler "ayni enzim TIPI" guvenini olcer, substrat garantisi degil.

## Regulasyon mimarisi (adim 4e)

Operonun 5' ucu belirlenir, yukari akistaki ilk gene kadarki intergenik bolge
(putatif promotor bolgesi) olculur, o genin yonu ve duzenleyici ailesi kaydedilir.
Sonuc: divergent duzenleyici %25,2 · ayni yonde duzenleyici %4,2 · divergent
baska gen %42,9 · ayni yonde baska gen %16,9 · pencerede gen yok %10,8.
Divergent duzenleyicilerde intergenik medyan 165 bp, digerlerinde 248 bp; divergent
duzenleyicilerin %94,6'si ≤400 bp (paylasilan promotor bolgesi araligi).
Aileler: LysR 878, TetR 487, MarR 281, IclR 275, GntR 167 — LysR baskinligi
aromatik katabolizma literaturuyle ortusuyor.
**Sinir — ve asil nedeni:** promotor DIZISI tahmin edilmiyor. Bunun sebebi
sadece veritabani tasarimi degil, KAYNAK VERININ KENDISI: indirilen 17.073
GenBank kaydinin **%88,9'u `CON` tipinde**, yani dizi dosyada yok, referansla
kuruluyor (olculdu: `awk '/^LOCUS/ && / CON /'`). BioPython bu kayitlarda
`len(record.seq)` icin bildirilen uzunlugu verir ama icerige erisilemez
(`UndefinedSequenceError`). Dizi tasiyan 1.897 kayit ise ookaryot mRNA'lari
(`XM_`/`NM_`), yani operon analizi icin uygun degil.
Sonuc: -35/-10 kutulari, operator tekrarlari ve transkripsiyon baslangici bu
veriyle ANALIZ EDILEMEZ; bunun icin 15.000 kaydin dizili surumunun yeniden
indirilmesi gerekir (bkz. ROADMAP).

## Filogeni ve dizi uzayi (adim 13b)

71 referans + her yaprak icin 1 temsilci (1.195) hmmalign ile katalitik cekirdege
hizalanir, >%50 bosluklu kolonlar atilir (426 → 350 kolon), FastTree ile agac
kurulur. Ayrica diamond all-vs-all ile SSN (30.935 kenar ≥%30) ve kume×kume
en yuksek temsilci kimligi matrisi uretilir. 44 kume icin ayri agac.

## Operon kuralinin sinanmasi (adim 12c2)

Operon tanimi bir konvansiyondur (ayni iplik, bosluk ≤150 bp) ve transkripsiyon
verisi olmadan dogrudan dogrulanamaz. Dogrulanabilen sey, dizi ile dogrulanmis
ortaklarin rastgele bir genden daha yakin ve daha sik ayni iplikte olup
olmadigidir. Arka plan ayni pencerelerdeki 192.541 genin tamamidir
(ayni iplik %56,6, ±2 gen icinde %21,5):

| bilesen | n | ayni iplik | ±2 gen | medyan bosluk |
|---|---|---|---|---|
| beta | 3.050 | %90,7 (1,60x) | %83,5 (3,89x) | 0 bp |
| ferredoksin | 1.849 | %80,9 (1,43x) | %64,9 (3,02x) | 514 bp |
| reduktaz | 6.492 | %70,2 (1,24x) | %47,8 (2,23x) | 770 bp |
| **baska alfa (kontrol)** | 1.833 | %54,7 (0,97x, p=0,10) | %25,1 (1,17x) | 2.380 bp |

Son satir negatif kontroldur: pencerede bulunan baska bir RO alfa alt birimi
arka plandan FARKSIZDIR. Yontem "alfanin yanindaki her seyi" buluyor olsaydi o
satir da zenginlesmis gorunurdu. Beta'nin medyan boslugunun 0 bp olmasi
(stop ve start kodonlarinin bitismesi) translasyonel eslesmenin klasik imzasidir.

**Esik duyarliligi** (50→500 bp): beta 2.371→2.518 (%6 degisim, esikten bagimsiz),
ferredoksin 921→1.317 (%43), reduktaz 2.214→3.411 (%54). Yani beta sonucu saglam,
ferredoksin/reduktaz yuzdeleri 150 bp konvansiyonuna bagli okunmali.

## Reaksiyon semalari ve taksonomi agaci (web)

Tip sayfalarinda substrat SMILES'ten cizilir, urun ADLA verilir. Urun yapisi
cizilmedi: bu her tipte hangi halka pozisyonunun saldiriya ugradigini varsaymayi
gerektirir ve yanlis regiokimya gostermek hic gostermemekten kotudur. Onun yerine
her REAKSIYON SINIFI icin genel mekanizma semasi cizilir (`atlas.reaction_scheme_svg`);
bunlar substrattan bagimsiz ve kesindir. Taksonomi sayfasinda yigili cubugun
yaninda NCBI soyu uzerinde katlanabilir agac var (`atlas.taxonomy_tree`), kirpilan
dugum sayisi acikca yaziliyor.

## Veri fazlaligi (adim 12c3)

11.422 giris, 10.019 tekil dizi → **%12,3 fazlalik**. Bir dizinin en fazla 29
kopyasi var; 2.322 giris kopyali bir dizi tasiyor. Kume boyutlari tekillestirmede
%5-16 dusuyor (VanA -%15, CntA -%12, KshA -%10). Fazlalik en cok dizilenmis
cinslerde birikiyor (Pseudomonas 415, Acinetobacter 304).

Girisler SILINMIYOR: ayni protein farkli genomda farkli komsulukta bulunur ve
genomik baglam bu veritabaninin asil konusudur. Onun yerine `stats_overview.py`
her ikili testi **uc duzeyde** kosar: giris basina, tekil dizi basina, tip×cins
basina. Plazmit bulgusu ucunde de ayakta (4,20x / 3,48x / 4,46x), yani tekrarlanan
suslarin eseri degil.

## Kimlik substrati ne kadar ongoruyor (adim 12c5)

"Kimlik substrat garantisi vermez" iddiasi tek ornekle (EdoA1/CumA1) degil,
871 referans cifti uzerinde SINIFLANDIRICI olarak olculdu. Ciftlerin yalnizca
%3,1'i ayni substrati paylasiyor; kesinlik bu taban orana karsi okunmali.

| esik | cift | ayni | kesinlik | lift | duyarlilik |
|---|---|---|---|---|---|
| %40 | 171 | 18 | %10,5 | 3,4x | %67 |
| %70 | 42 | 8 | %19,1 | 6,1x | %30 |
| %90 | 15 | 4 | %26,7 | 8,6x | %15 |
| %95 | 5 | 4 | %80,0 | 25,8x | %15 |

**AUC = 0,859**: kimlik siralama sinyali olarak gercekten bilgi tasiyor
(0,5 = hicbir bilgi). Ama KARAR KURALI olarak basarisiz: %90 kimlikte 15 ciftin
yalnizca 4'u ayni substrati paylasiyor ve **hicbir esik (n>=10 iken) %90 kesinlige
ulasmiyor**. En benzer farkli-substratli cift EdoA1/CumA1 %99,8; ayni substrati
paylasan en uzak cift %32,4.

Substrat SINIFI duzeyinde (ksenobiyotik/dogal) kesinlik %99'a cikiyor ama bu
yanıltici: ciftlerin **%83'u zaten ayni sinifta** cunku karakterize enzimlerin
cogu ksenobiyotik uzerine. Lift hicbir esikte 1,2'yi gecmiyor, yani esik bilgi
katmiyor; dogru ozet AUC = 0,706.

## Tip birliktelikleri (adim 12c4)

Operon dogrulamasinda ortaya cikan 1.833 "komsu alfa alt birimi" takip edildi.
Ayni +-10 kb icinde 854 dogrulanmis RO cifti var: **214 ayni tip** (duplikasyon),
**640 farkli tip**. Yani bu agirlikla gen duplikasyonu degil.

Iplik bilgisi iki deseni ayiriyor:
- `CntA + BmoA` 54 cift, 44'u ayni iplikte → tek transkripsiyon birimi adayi;
  ikisi de metilamin birakir, metabolik olarak tutarli.
- `KshA15 + CntA` 61 cift, 0'i ayni iplikte → ayni genomda ama ayri transkribe
  ediliyor; buyuk bir genomda iki alakasiz yetenek.

Replikon duzeyinde permutasyon null'i (2.000 tekrar; replikon basina RO sayisi ve
tip toplamlari SABIT, boylece "ikisi de yaygin" aciklamasi elenir):

| cift | gozlenen | beklenen | oran | z |
|---|---|---|---|---|
| BPDO + PhnA1a | 16 | 0,13 | 123,6x | 43,2 |
| BPDO + TPDO | 19 | 0,39 | 48,2x | 30,2 |
| PhnA1a + TPDO | 31 | 1,39 | 22,3x | 25,6 |
| TPDO + CmoS | 42 | 4,14 | 10,1x | 19,0 |
| AntA + XylX | 23 | ~4,1 | 5,6x | 10,6 |
| TdnA1 + AntA | 25 | ~4,3 | 5,8x | 10,3 |

En guclu birliktelikler **urunleri ayni asagi yola giren** enzimler: anilin,
antranilat ve benzoat oksijenazlarinin ucu de katekol veriyor (beta-ketoadipat
yolu). Yani genomlar tek enzim degil, tum huni ediniyor gibi gorunuyor.
**Sinir:** ayni replikonu paylasmak tek organizmada ortak yol KANITI degildir;
bunlar dizilenmis genomlar uzerinde sayimlar ve PAH yikan izolatlar iyi
calisilmis bir grup oldugu icin orneklem yanliligi da ayni yone iter.

## Izolasyon kaynagi ve habitat (adim 12c6)

GenBank kayitlarinin `source` ozelliklerinde kullanilmayan ekolojik veri vardi:
9.591 dosyada `/isolation_source`, 4.911'inde `/host`, 11.485'inde cografi konum.
`isolation_source.py` bunlari cikarir, serbest metni 21 kategorili siralı ve
denetlenebilir bir kelime haritasiyla normalize eder (`replicon_source` tablosu +
`analysis_out/habitat.json`).

**Kapsam durust verilir:** girislerin %60,9'unda kaynak var, **%55,6'si** bir
habitate yerlestirilebiliyor. 4.493 girişte kaynak hic yok, 578'inde metin
siniflandirilamiyor (214 farkli dizgi: etiketsiz ontoloji numaralari, bitki ve
hayvan orneklerinde de gecen ciplak anatomik kelimeler, "culture"/"tissue").
Yanlis siniflamaktansa siniflamamak yegdir.

**Tur bazinda normalizasyon neden zorunlu:** giris sayisi dizilenmis suslari
sayar. Insan klinigi ve bitki iliskili kaynaklarda tur basina 6,3 ve 6,0 giris
var, denizde 3,0. Bu yuzden her ekolojik ifade TUR uzerinden kurulur.

**Sonuc.** Kirli/sanayi sahasi arka plani tum tur-habitat gozlemlerinin %2,4'u.
Tur bazinda zenginlesme:

| | tur payi | kat |
|---|---|---|
| alkilbenzenler | %17,9 | 7,6x |
| nitroaromatikler | %16,7 | 7,1x |
| PAH | %13,3 | 5,6x |
| biaril/eterler | %11,5 | 4,9x |
| ksenobiyotik (sinif) | %6,1 | 2,6x |

Yani ilişki **sinifa degil belirli kimyasal ailelere** ait. Alkilbenzen ve
nitroaromatik satirlari birkac ture dayaniyor, oran olarak okunmamali; PAH
16 turle en guvenilir olani.

**Bir sozluk karari sonucu degistirdi.** Madencilik ve asit maden drenaji
baslangicta kirli/sanayi kategorisinin icindeydi. Kirliliktir ama metaliktir,
organik degil, ve kategoriyi suluyordu: ayirinca PAH zenginlesmesi ~3x'ten
5,6x'e cikti. Madencilik artik ayri bir habitat ve FARKLI, daha zayif bir aile
kumesini zenginlestiriyor — ayirmanin dogru oldugunun capraz kontrolu bu.

## Arama (adim 12d)

Eski arama bes alanda `LIKE '%...%'` yapiyordu: kelime siniri yok, siralama yok,
iki kelime yazinca hic sonuc yok. Simdi her dogrulanmis RO icin tek bir arama
belgesi (`ro_search.doc`) kurulur ve FTS5 ile indekslenir; sonuclar bm25 ile
siralanir. Filtre alanlari ayni tabloda: kanit duzeyi, yasam alani, kimyasal
aile, tip, plazmit, operon ortagi, divergent duzenleyici.
Statik sitede ayni filtreler tarayicida calisir (`search_index.json`, 2,3 MB).

## Web uygulaması

`webapp/` — FastAPI + Jinja2. Sayfalar: ana sayfa, tip listesi, tip sayfasi
(kanit kademesi, operon ortaklari, regulasyon, kume agaci, varyantlar, uyeler),
varyant sayfasi, giris sayfasi (gen komsulugu SVG, operon tablosu, promotor
bolgesi, kanit kademesi), arama, dizi siniflandirici ve **Atlas** bolumu:
filogeni, dizi uzayi (SSN + kimlik matrisi + referans kalibrasyonu), taksonomi,
operonlar/ETC, regulasyon, ekoloji, kanit duzeyleri, istatistikler.
Her sekilde yontem + uyari notu var. Yayınlama: `webapp/README_DEPLOY.md`
(GitHub Pages statik dışa aktarım `freeze.py`; tam uygulama için Docker /
Hugging Face Spaces).

## Bilinen sınırlar

- **Operon tamlığı** artık dizi tabanlı ölçülüyor (adım 4b). `analyze_ecology.py`
  içindeki anotasyon-metni tabanlı `operon_complete` eski ölçüttür; karşılaştırma
  için tutuldu.
## Varyant kalinti imzasi (adim 12b)

Her kume icin varyantlari AYIRAN hizalama kolonlari bulunur: kolon uyelerin
≥%90'inda dolu, varyantlar arasi en fazla 4 farkli kalinti, en sik iki kalinti
varyantlarin ≥%70'ini kapsiyor ve kalinti varyant ICINDE korunmus (entropi ≤0,5).
Bu kapilar olmadan metrik yalnizca en kotu hizalanmis ilmekleri buluyordu
(olculdu: CntA'da her varyant farkli kalinti, ~11 durum).
Sonuc: 35 kumede imza, 1.749 varyant; kolonlarin %48'i katalitik bolgede,
varyantlar kume cogunlugundan ortalama 4,8/12 kolonda ayriliyor.
**Sinir:** bunlar istatistiksel olarak ayirt edici kolonlardir, yapisal olarak
dogrulanmis baglanma cebi kalintilari degildir; veritabaninda yapi yok.

## Korunmus merkez olcumu (adim 12c)

11.422 dogrulanmis RO'da: Rieske ligandlarinin 4/4'u %100, kopru Asp/Glu %100,
katalitik triad 3/3 %87,3 (kalan %12,7'de 2/3). En oynak pozisyon Fe(II)
karboksilati (kolon 355, %87,6) — hem Asp hem Glu kabul ediyor.
- ~~7 küme substratı bilinmiyor~~ **COZULDU** (2026-10-06): hepsi literaturden
  belirlendi ve `chemistry.csv`'ye kaynakla birlikte islendi. OxoO ve OMO ayni
  enzim (2-oksokinolin 8-monooksijenaz, PDB 1Z03) — referans seti bu enzimi iki
  kez iceriyor, ki ikisinin %100 kimlikli olmasi bunu zaten gosteriyordu.
  CndA kloroasetanilid herbisit N-dealkilazi, PsbAb 4-sulfobenzoat
  3,4-dioksijenazi, ROCH34 ftalat 4,5-dioksijenazi (PDB 7FHR), OxyA/qxyA
  benzalkonyum klorur (QAC) oksijenazi. cadA 2,4-D oksijenazi olarak
  **kesin degil** diye isaretlendi.
- **Substrat ataması küme düzeyinde**: heterojen kümelerde tek substrat tüm üyeler
  için geçerli değil; yaprak (varyant) düzeyinde yapılmalı.

## Çıktılar

- `roar.sqlite` — ana veritabanı (replicon, ro, neighbor, gene_category, subfamily,
  ro_subfamily, leaf, ro_leaf, leaf_profile, ro_domain)
- `analysis_out/` — tüm analiz CSV'leri + novel/temsilci FASTA'ları
- `analysis_out/provenance.json` — köken manifestosu (`provenance.py`): her tablo
  ve her analiz dosyası için üreten script, tükettiği girdiler, satır/kayıt
  sayısı, değişiklik zamanı, kısa İngilizce açıklama ve zincirin dayandığı
  kaynak veri
- `report.html`, `explorer.html`, `hub.html` — üç HTML sayfa (DB'den canlı üretilir)

Sınıflandırmayı değiştirmek için gbk'ları yeniden parse etmek gerekmez:
`python3 build_db.py --recategorize` yalnızca `gene_category` tablosunu yeniden üretir.

## Referans setinde bulunan tekrarlar (2026-10-06)

`validate_curation.py` iki SERT hata veriyor ve ikisi de kuratorluk karari
gerektiriyor, bu yuzden kod tarafindan sessizce duzeltilmedi:

1. **71 referans, 68 tekil dizi.** Uc cift birebir ayni: `1_101_OxoO` = `3_304_OMO`,
   `3_309_NahAc` = `3_315_NDO`, `3_314_NDO` = `3_316_NarAa`. Sonuc kozmetik degil:
   `RieskeDB71.hmm` ayni profili iki kez iceriyor ve atama ikizler arasinda
   keyfi bolunuyor (OxoO 2 / OMO 11, NDO 2 / NarAa 30, NahAc 0 / NDO 0).
   "En iyi profil ile en yakin referans protein uyusmuyor" oraninin %29,3
   olmasinin bir kismi bundan.
2. **`1_115_NdmC`, `1_114_NdmB`'nin alt dizisi** (355 aa, 373 aa icinde).
   Gercek NdmB ve NdmC ~%65 benzer ayri demetilazlardir, yani ayni protein iki
   kez girilmis, biri N-ucundan kisaltilmis. NdmC sifir uye topluyor, NdmB 15.
3. Uyari duzeyinde: `1_102_CARDO` / `1_103_CarAa` %99,2 ve CarAa sifir uye
   topluyor. Ayni desen ama sert kontrolu bir kalinti farkla gecmiyor; otomatik
   kural bunu EdoA1/cumA1 (%99,8, farkli substrat, gercek bir bulgu) vakasindan
   ayirt edemez, bu yuzden karar insana birakildi.

Hangi ikizin tutulacagi ve NdmC'nin kendi dizisiyle yeniden alinip alinmayacagi
kuratore aittir; degistirildiginde `RieskeDB71.hmm` yeniden kurulmali ve
pipeline 1. adimdan itibaren kosulmalidir.

## Safeguards (adim 13d/13e) — ne yakalanir, ne yakalanmaz

Bu pipeline birkac SESSIZ hata yayinladi: yutulan cd-hit cokmesi, kolonlari
kaydiran `2\,4` CSV kacisi, en kotu hizalanmis ilmekleri "ayirt edici kolon"
sanan metrik, `is_confirmed=NULL` kalan 945 aday. Hepsinin ortak yani ayni:
hata ciktiyi BOZDU ama hicbir sey durmadi. Iki script bunu adres aliyor ve
`run_all.sh` icinde HTML uretiminden ONCE kosuyor.

### `validate_curation.py` — kapi (sifir-disi cikar)

Her kontrol icin tek satir `PASS` / `FAIL` / `WARN` ve bir sayi yazar. Bir
`FAIL` varsa sifir-disi cikar, yani `set -e` altinda pipeline durur ve bozuk
veri web sitesine gitmez. Kapsam:

| grup | ne sinaniyor |
|---|---|
| **referans seti** (`ROs_71_Clean/refs71.fasta`) | iki referans dizisi birbirinin AYNISI olmasin; bir dizi bir baskasinin ICINDE gecmesin (parca referans); her kimlik FASTA + `cluster_ecology.csv` + `chemistry.csv` ucunde de tam bir kez olsun. Ayrica WARN olarak: %99 uzeri kimlikli ciftler (substratlariyla) ve hic uye toplamayan referanslar (ikiz/parca etiketiyle) |
| sozluk | `add_reference.py`'nin YAZARKEN kullandigi izinli deger listeleri ile dogrulayicinin KONTROL ETTIGI listeler ayni mi (ayrisirsa yeni referans sessizce gecer) |
| CSV bicimi | her satirin alan sayisi baslikla ayni mi — `2\,4` kacisi hata sinifi tam olarak bu |
| `ro` | dogrulanan her RO'nun gercek bir `ro_cluster`'i (`NULL`/`'N/A'` degil), dizisi ve `rieske_intact=1`'i var mi; `is_confirmed` hic `NULL` kalmis mi |
| kapsama | `ro`'daki her kume iki kuratorlu CSV'de var mi; CSV'lerde olup `ro`'da uyesi olmayan tip var mi (bilinen 10 uyesiz referans `REFERENCE_ONLY_CLUSTERS` ile muaf, muafiyet bayatlarsa `WARN`); iki CSV ayni tip kumesini kapsiyor mu; tekrarlanan satir var mi |
| `cluster_ecology.csv` | `substrate_class` ∈ {xenobiotic, natural_aromatic, natural_specialized, unknown}, `confidence` ∈ {low, medium, high} |
| `chemistry.csv` | `reaction_class` ve `family` bilinen kumede, `source_kind` bilinen kumede, `curation_confidence` ∈ {low, medium, high}; her bos olmayan `substrate_smiles` makul (parantez/koseli parantez dengesi, yalnizca mesru atom-bag karakterleri, en az bir atom, eslesen halka kapanis rakamlari); her bos olmayan `pdb` tam dort alfanumerik; adlandirilmis substratin kaynagi yazili mi |
| tutarlilik | iki CSV bir tipin substratinin bilinip bilinmedigi konusunda celisiyor mu (chemistry substrati adlandirmis ama ecology hala `unknown` diyorsa ekoloji testi o tipi disarida birakir ve web sayfasi iki farkli sey soyler) |
| dil | `leaf_profile.label`, `leaf.top_genera`, `ro_etc.etc_profile`, `operon.layout` ve iki CSV'nin HER alani; ayrica `analysis_out/` icinde web'in yayinladigi metin dosyalari. Iki tuzak ayri ayri araniyor: Turkce'ye ozgu harfler **ve** ASCII'ye sadelestirilmis Turkce kelimeler (`agirlikli`, `komsu`, `kume`, `yaprak`, `bilinmiyor`, `dusuk`, `orta`, `yuksek`, `degil`) — ikincisi daha sinsi, cunku kodlama kontrolleri onlari gormez |
| referans butunlugu | `neighbor.candidate_id` → `ro`; `ro_leaf.leaf_id` → `leaf`; `leaf_sdp.cluster` → `cluster_sdp`; `leaf_sdp.leaf_id` → `leaf`; `gene_category.neighbor_id` → `neighbor`; `ro.nucleotide_id` → `replicon`; `ro_subfamily.subfamily_id` → `subfamily`; `leaf.size` o yapragin `ro_leaf` uye sayisina esit mi ve `sum(leaf.size)` = `ro_leaf` satir sayisi mi; `subfamily.size` ayni sekilde |
| aritmetik | dogrulanan RO sayisi `ro`, `ro_search`, `ro_evidence`, `ro_domain`, `operon`, `ro_leaf`, `ro_etc`, `ro_regulation`, `ro_subfamily` arasinda ayni mi (ayrisirsa farkli sayilar yazdirilir); her dogrulanan RO her tureyen tabloda var mi; tureyen tablolar yalnizca dogrulanmis RO tutuyor mu |
| tazelik | `ro_search`'un substrat/aile/reaksiyon alanlari kuratorlu CSV'lerin SU ANKI halini mi yansitiyor (SQLite tablo basina zaman damgasi tutmadigi icin bu ancak icerik karsilastirmasiyla gorulebilir) |

**Referans seti kontrolu neden sert bir kapi.** Pipeline'in TAMAMI
`ROs_71_Clean/refs71.fasta` uzerine kuruluyor: tip atamasi, katalitik motif
kolonlari, her kanit kademesi. Ayni enzim iki adla girildiginde
`RieskeDB71.hmm` birbirinin aynisi iki profil tasir ve atama ikizler arasinda
KEYFI boluur — olculdu: OxoO 2 / OMO 11, NDO 2 / NarAa 30, NahAc 0 / NDO 0.
"En yuksek skorlu profil ile en yakin referans protein farkli tipte" orani
(%29,3) buyuk olcude bundan sisiyor. Parca referans ayni sorunun diger yuzu:
`1_115_NdmC` (355 aa) `1_114_NdmB`'nin (373 aa) harfi harfine bir parcasi, oysa
gercek NdmB ve NdmC ~%65 kimlikli ayri demetilazlar; parcadan kurulan profil
sistematik olarak zayif ve NdmC hic uye toplamiyor.

**%99 uzeri kimlik neden FAIL degil WARN.** Bu veritabaninin merkezi bulgusu
tam olarak bu: `2_202_EdoA1` ve `2_203_cumA1` %99,8 kimlikte FARKLI
substratlara (etilbenzen / kumen) etki ediyor. Bunu kusur saymak bulguyu
silmek olur. Duplikasyon kusuru ise AYRI ve sert kontrol; tam kopya ve parca
ciftleri WARN listesinden cikarilir, yoksa gercek bulgu gurultuye karisir.

### `provenance.py` — koken manifestosu (yuksek sesle uyarir)

`analysis_out/provenance.json` yazar. Script adlari README'nin pipeline
tablosundan capraz kontrol edilir ama YETKILI olan `provenance.py` icindeki
modul duzeyi `TABLE_PROVENANCE` / `FILE_PROVENANCE` dict'leridir. Uyari
siniflari:

- **haritada olmayan tablo/dosya** — haritalanmamis bir artefakt, sessizce
  bayatlayan seyin ta kendisidir (bu kontrol `ro_carboxylate` tablosunu
  bulmustur: `motif_stats.py` kuruyor, `stats_overview.py` okuyor, hicbir
  belgede yoktu);
- **haritada olup diskte/veritabaninda olmayan artefakt** — o adim kosulmamis;
- **ureticisi dogrulanamayan artefakt** — haritadaki script'in KOD GOVDESINDE
  (docstring ve yorumlar atilarak) cikti adi hic gecmiyorsa, ya harita
  yanlistir ya da artefakt yetimdir;
- **girdisi kendisinden yeni olan artefakt** — bayat zincir.

Her artefakt icin zincir kaynak veriye kadar cozulur: `combined_pfam.fasta`,
`gbk_files/`, `ROs_71_Clean/`, `cluster_ecology.csv`, `chemistry.csv`
(+ `RieskeDB71.hmm`, `component_hmms/`).

`provenance.py` varsayilan olarak sifir doner (yalnizca uyarir); CI'da sert
kapi istenirse `--strict`.

### Yakalanamayanlar — durust liste

Bunlar kontrol edilmiyor ve bu dogrulayici gectigi icin veri "dogru" olmus
sayilmaz:

1. **Kuratorlu bilginin DOGRULUGU.** Bir substrat yanlis makaleden alinmissa,
   SMILES yanlis molekulu cizse, PDB kodu baska bir yapiya isaret etse
   kontroller gecer. Sinanan sey bicim ve ic tutarliliktir, literatur degil.
2. **Referans setinde YAKLASIK duplikasyon.** Sert kontroller yalnizca TAM
   kopyayi ve harfi harfine PARCAYI gorur. Tek kalintisi degisen bir kopya
   ikisinden de kacar; o ancak %99 uzeri kimlik WARN'inda gorunur ve karari
   insan verir. Ornegin `1_102_CARDO` / `1_103_CarAa` %99,2 kimlikte ve AYNI
   substrata (karbazol) atanmis, CarAa hic uye toplamiyor — bu desen
   duplikasyona benziyor ama otomatik olarak boyle ilan edilmiyor. Ayrica o
   WARN'in kimlik olcumu `analysis_out/reference_pairs.csv`'den okunur: o dosya
   bayatsa veya yoksa kontrol sessizce zayiflar (yoksa WARN yazar).
3. **SMILES kimyasi.** Kontrol parantez dengesi ve mesru karakter duzeyinde;
   valans, aromatiklik ve stereokimya sinanmiyor (RDKit kurulu degil). Bicimsel
   olarak gecerli ama kimyasal olarak sacma bir SMILES gecer.
4. **SQLite tablo basina tazelik.** SQLite tablolar icin zaman damgasi tutmaz.
   Dosya duzeyinde bayatlik yakalanir; tablo duzeyinde yalnizca ozel olarak
   yazilmis icerik karsilastirmalari yakalar (su an sadece `ro_search` ↔ iki
   CSV). Baska bir tablo bir CSV'den kopyalama yapmaya baslarsa o kontrol elle
   eklenmeli. Dosya duzeyindeki olcut de KABA: mtime karsilastirmasi bir
   girdinin NE kadarinin degistigini bilmez, bu yuzden bir CSV'de tek alan
   duzeltilince o CSV'ye bagli TUM ciktilar "bayat" gorunur. Uyari her zaman
   "yeniden kos" demek degil, "yeniden kosulup kosulmayacagina KARAR VER"
   demektir.
5. **Istatistiksel iddialar.** `stats.json`, `null_model.csv`,
   `cooccurrence.json`, `substrate_predictability.json` icindeki sayilarin
   dogrulugu sinanmiyor; yalnizca dosyalarin var oldugu ve girdilerinden yeni
   oldugu. Yanlis bir null model sessizce gecer.
6. **Biyolojik anlam.** Operon tanimi (ayni iplik, bosluk ≤150 bp) bir
   konvansiyondur; yaprak = varyant esitligi bir karardir; kume duzeyi substrat
   atamasi heterojen kumelerde uyelerin cogu icin gecerli degildir. Hicbiri
   dogrulanabilir bir onerme degil.
7. **Kaynak verinin kendisi.** `gbk_files/` kayitlarinin %88,9'u `CON` tipinde,
   yani dizi dosyada yok; bu eksiklik kontrollerle giderilemez (bkz. yukarida
   regulasyon mimarisi bolumu).
8. **Ingilizce metnin kalitesi.** Dil kontrolu Turkce harf ve sabit bir Turkce
   kelime listesi arar. Listede olmayan bir Turkce kelime, bozuk Ingilizce ya
   da yarim cumle gecer. Ters yonde de kusurlu: arama ONEK eslesmesi yapar
   (Turkce eklemeli oldugu icin gerekli) ve bu mesru metni yakalayabilir —
   GenBank'tan gelen `New Zealand: Kumeu` yer adi 'kume' onekine takiliyor.
   Bilinen yanlis pozitifler `LANGUAGE_FALSE_POSITIVES` ile muaf tutuluyor; o
   liste uzadikca kontrolun degeri duser.
9. **`webapp/` sablonlari.** Dogrulayici veritabani ve iki CSV'ye bakar; HTML
   sablonlarindaki sabit yazilmis metinler ve sayilar kapsamda degildir —
   daha once tam bu sinifta hata cikmisti (`make_report.py` icinde sabit
   yazilmis 960 / "4,3x" degerleri DB ile celisiyordu).
