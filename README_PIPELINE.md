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
| 12d | `build_search_index.py` | FTS5 tam metin indeksi + filtre alanlari → `ro_search`, `ro_fts` |
| 13b | `build_phylogeny.py` | hmmalign + FastTree → `tree_all.nwk`, `trees/<kume>.nwk`; diamond all-vs-all → `ssn_*.csv`, `cluster_identity_matrix.csv` |
| 13c | `stats_overview.py` | hipotez testleri (entry + genus duzeyi, etki buyuklugu) → `stats.json` |
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
  kümelerde. Özyinelemeli homojenizasyon 61 kümeyi 1.790 yaprağa böldü.
- **Varyantlar genomik bağlamla ayrışıyor** (KshA'nın bir varyantı steroid dehidrogenaz
  yanında; VanA'nın mobil *Burkholderia* soyu IclR+transpozon yanında).
- **Novel**: 155 yüksek-güven aday (güvenle RO, bilinen tipe uymuyor) + bitki/alg dalı.
- **Ekoloji**: ksenobiyotik RO'lar plazmit üzerinde 4,3x daha sık (yatay transfer);
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
**Sinir:** promotor DIZISI tahmin edilmiyor (-35/-10, operator). DB'de intergenik
DNA yok, yalnizca komsu protein cevirileri var.

## Filogeni ve dizi uzayi (adim 13b)

71 referans + her yaprak icin 1 temsilci (1.195) hmmalign ile katalitik cekirdege
hizalanir, >%50 bosluklu kolonlar atilir (426 → 350 kolon), FastTree ile agac
kurulur. Ayrica diamond all-vs-all ile SSN (30.935 kenar ≥%30) ve kume×kume
en yuksek temsilci kimligi matrisi uretilir. 44 kume icin ayri agac.

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
- **7 küme substratı bilinmiyor** (`OxoO/CndA/cadA/OMO/PsbAb/ROCH34/OxyA`): standart
  veritabanlarında bulunamadı, orijinal referans makaleleri gerekir.
- **Substrat ataması küme düzeyinde**: heterojen kümelerde tek substrat tüm üyeler
  için geçerli değil; yaprak (varyant) düzeyinde yapılmalı.

## Çıktılar

- `roar.sqlite` — ana veritabanı (replicon, ro, neighbor, gene_category, subfamily,
  ro_subfamily, leaf, ro_leaf, leaf_profile, ro_domain)
- `analysis_out/` — tüm analiz CSV'leri + novel/temsilci FASTA'ları
- `report.html`, `explorer.html`, `hub.html` — üç HTML sayfa (DB'den canlı üretilir)

Sınıflandırmayı değiştirmek için gbk'ları yeniden parse etmek gerekmez:
`python3 build_db.py --recategorize` yalnızca `gene_category` tablosunu yeniden üretir.
