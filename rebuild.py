"""Hangi analiz BAYAT, ve hangi sirayla yeniden uretilmeli.

SORUN. `run_all.sh` dogrusal bir betik ve atlama kurali dosyanin VAR OLUP
olmamasina bakiyor, GUNCEL olup olmamasina degil. Bir kuratorluk duzeltmesi
`chemistry.csv`'ye islendiginde `stats.json` yerinde durmaya devam eder, site
onu okur, ve sayfa artik veriyle uyusmayan bir sayi gosterir. Hicbir sey
bagirmaz. Ayrica bu oturumda yazilan on alti analiz `run_all.sh` icinde HIC
gecmiyor, yani onlar hicbir zaman kendiliginden yenilenmiyor.

BU DOSYA NE YAPAR. Her uretilmis artefakt icin uc seyi bilir: hangi komut
uretir, neye BAGLIDIR, ve son uretildiginde girdilerin icerigi neydi. Bayatlik
uc sekilde olusur ve ucu de yakalanir:

    cikti yok
    bir girdi ciktidan YENI
    bir girdinin ICERIGI son uretimden beri degisti (mtime degismemis olsa bile)

Ucuncusu sart: `git checkout` ve kopyalama mtime'i ileri atar ama icerik ayni
olabilir; tersine bir dosya ayni saniyede degistirilebilir. Icerik ozeti
ikisini de dogru cozer.

UC KIP:
    --check    hicbir sey calistirmaz, bayat olanlari listeler, bayat varsa
               1 doner. Yayin oncesi kapi icin budur.
    --plan     ne calisacagini sirayla gosterir, yine calistirmaz
    --run      bayat olanlari BAGIMLILIK SIRASINA gore uretir

Agir adimlar (17 bin genomdan baglam cikarimi, 189 bin protein filtresi) ayri
isaretlidir ve --include-heavy verilmeden calistirilmaz: bunlar saatler surer
ve bir kuratorluk duzeltmesinden sonra gerekmezler.
"""

import argparse
import functools
import hashlib
import json
import glob
import os
import subprocess
import sys
import time

STATE = "analysis_out/.rebuild_state.json"
OUT = "analysis_out"
# Yer isteyen ara dosyalar sistem diskine DEGIL buraya yazilir (proje kurali).
# Yoksa carboxylate_audit ve His capalamasi tempfile.mkdtemp ile /tmp'ye
# gider; ikisi de 19.125 dizinin tip-ici hizalamalarini uretiyor.
SCRATCH = os.environ.get("ROAR_SCRATCH", "/media/lin-bio/back2")

# Cekirdek girdiler. Bunlardan biri degisirse asagidaki her sey bayatlar.
DB = "roar.sqlite"
CHEM = "chemistry.csv"
ECO = "cluster_ecology.csv"
# Kuratorun fazlalik ciftleri icin verdigi karar: olcumun degil, insanin ciktisi.
DECISIONS = "curation_decisions.csv"
REFS = "ROs_71_Clean/refs71.fasta"
MOTIF_HMM = "ROs_71_Clean/ROmotif71.hmm"
# RefSeq taramasinin is dizini. SCRATCH altinda, cunku ara dosyalar 40 MB'lik
# parcalar halinde inip siliniyor ve sistem diskine yazilmamalari gerekiyor.
NCBI = os.path.join(SCRATCH, "roar_ncbi")

# Her kayit: (ad, ciktilar, komut, girdiler, agir mi)
#
# Girdilere BETIGIN KENDISI de dahildir: analiz kodu degistiginde ciktisi
# bayattir, ve bu oturumda tam olarak bu oldu -- stats_overview.py degisti,
# stats.json elle yeniden uretilmeseydi sayfa eski sayilari gosterecekti.
ARTEFACTS = [
    # --- cekirdek pipeline (agir olanlar isaretli)
    # reference_models: types/ -> chemistry.csv + refs71.fasta zincirinin SON
    # halkasi. GRAFIKTE YOKTU (2026-10-28'de bulundu): biri types/'a referans
    # ekleyip types_build.py --write kostuktan sonra rebuild.py --run
    # calistirsa, motif modeli YENIDEN KURULMAZDI ve butun zincir eski modelle
    # koscaktI. Bugune kadar yalnizca elle kosuldugu icin sorun cikmadi --
    # yani dogru sonuc bir insanin hatirlamasina bagliydi. Bu, bu oturumda
    # bulunan ALTINCI eksik kenar.
    #
    # RieskeDB71.hmm (atama kutuphanesi) BILEREK burada DEGIL: her tip icin
    # profili build_type_profile.py kuruyor ve bu, aday havuzuna ve kurator
    # kararina (kayma sinamasi, yarisma onizlemesi) bagli bir adim.
    ("reference_models", [MOTIF_HMM, "ROs_71_Clean/ROmotif71_hmmaln.sto",
                          "ROs_71_Clean/motif_columns.json"],
     ["python3", "build_reference_models.py", "--refs", REFS,
      "--out-dir", "ROs_71_Clean", "--name", "ROmotif71"],
     ["build_reference_models.py", REFS], True),
    # oxygen_fate: types/ + reaction_classes.yaml -> oxygen_fate.csv. Hafif
    # ama GRAFIKTE OLMAK ZORUNDA: bir tipin reaksiyon sinifi degisip tablo
    # guncellenmezse site eski kaderi filtreler ve kimse farketmez.
    # chemistry.csv BILEREK girdi degil -- deger tiplerden ve sozlukten
    # turetilir, chemistry.csv'ye sutun eklenmedi (o dosyayi yirmiden fazla
    # betik okuyor).
    ("oxygen_fate", ["oxygen_fate.csv"],
     ["python3", "oxygen_fate.py", "--write"],
     ["oxygen_fate.py", "reaction_classes.yaml", "types"], False),
    # metabolic_role: types/ + metabolic_roles.yaml + cluster_ecology.csv ->
    # metabolic_role.csv. Ekoloji dosyasi GIRDIDIR cunku varsayilan kural
    # substrate_class'a bakiyor; o dosya degisirse rol yeniden cozulmelidir.
    ("metabolic_role", ["metabolic_role.csv"],
     ["python3", "metabolic_role.py", "--write"],
     ["metabolic_role.py", "metabolic_roles.yaml", "types",
      "cluster_ecology.csv"], False),
    ("protein_filter", ["ro_alpha.csv", "ro_alpha.fasta", "ro_domtbl.out"],
     ["python3", "run_ro_filter.py", "--fasta", "combined_pfam.fasta",
      "--hmm", "RieskeDB71.hmm"],
     ["run_ro_filter.py", "ro_filter.py", "RieskeDB71.hmm"], True),
    # gbk_files GIRDIDIR ve 2026-11-10'a kadar BEYAN EDILMEMISTI. Dort adim
    # (genomic_context, operons, null_model, habitat) o dizinden okuyor;
    # hicbiri onu girdi saymiyordu, yani genom eklense ya da bir kayit
    # yenilense bu adimlarin HICBIRI bayat gorunmuyordu. Dizin ozeti
    # 2026-11-05'te eklendiginde bu beyan mumkun hale geldi (bkz. dir_digest);
    # 3,2 GB / 17.074 dosya icin ozet ~2,8 s suruyor ve surec basina bir kez
    # hesaplaniyor (lru_cache).
    ("genomic_context", ["genomic_context/neighbors.csv",
                         "genomic_context/ro_candidates.fasta"],
     ["python3", "extract_genomic_context.py", "--gbk-dir", "gbk_files",
      "--out-dir", "genomic_context"],
     ["extract_genomic_context.py", "gbk_files"], True),
    ("database", [DB],
     ["python3", "build_db.py", "--context-dir", "genomic_context", "--db", DB],
     ["build_db.py", "genomic_context/neighbors.csv"], True),
    # db_repair: BILDIRIMLE uygulanan veri duzeltmeleri. `database`'den
    # HEMEN SONRA olmak zorunda. Sebep kodda kapandi ama ARA ARTEFAKT eski
    # kodla uretilmis: `genomic_context/replicons.csv` hayalet satiri hala
    # tasiyor ve `build_db.py` o CSV'yi okuyor -- yani grafikte olmasa bir
    # sonraki tam yeniden kurulum satiri geri getirir ve onu silecek
    # hicbir sey olmazdi (ROADMAP 135.9b). Idempotent: uygulanacak bir sey
    # yoksa yazmaz, hata da vermez.
    ("db_repair", [f"{OUT}/db_repair.json"],
     ["python3", "db_repair.py", "--db", DB, "--apply",
      "--out", f"{OUT}/db_repair.json"],
     ["db_repair.py", DB], False),
    # annotate'in ciktilari genomic_context/cand_dom.out ve cand_aln.sto'dur ve
    # bunlari IKI HMM uretir. Ikisi de girdi listesinde YOKTU: 2026-10-22'de
    # AerC kutuphaneye eklendiginde protein_filter bayatladi ama annotate
    # bayatlamadi, zincir dogrudan evidence_tiers'a atladi ve atama 10-14
    # tarihli cand_dom.out ile kaldi -- yeni tip SIFIR uye topladi, dogrulama
    # kapisi da "aciklanmamis fazla tip" diye bunu yakaladi. Ciktilari da
    # yaziyorum, yoksa "cikti yok" diye her kosuda bayat gorunurdu.
    ("annotate", ["genomic_context/cand_dom.out", "genomic_context/cand_aln.sto"],
     ["python3", "annotate_ro.py", "--context-dir", "genomic_context", "--db", DB],
     ["annotate_ro.py", DB, "RieskeDB71.hmm", MOTIF_HMM], True),

    # build_operons: ro.sequence, operon, operon_gene, neighbor_protein,
    # neighbor_component, neighbor_domain. Bu adim grafikte YOKTU; eski DB
    # sutunu tasidigi icin fark edilmedi, taze DB'de evidence_tiers
    # "no such column: sequence" ile dustu (2026-10-14).
    ("operons", [],
     ["python3", "build_operons.py", "--db", DB, "--gbk-dir", "gbk_files"],
     ["build_operons.py", DB, "component_hmms/components.hmm", "gbk_files"], True),

    # --- veritabanina yazan ara adimlar
    ("evidence_tiers", [f"{OUT}/evidence_by_cluster.csv", f"{OUT}/reference_pairs.csv"],
     ["python3", "evidence_tiers.py", "--db", DB, "--out-dir", OUT],
     ["evidence_tiers.py", DB, CHEM], False),
    ("etc_types", [f"{OUT}/etc_by_cluster.csv"],
     ["python3", "etc_types.py", "--db", DB, "--out-dir", OUT],
     ["etc_types.py", DB], False),
    ("regulation", [f"{OUT}/regulation_by_cluster.csv"],
     ["python3", "analyze_regulation.py", "--db", DB, "--out-dir", OUT],
     ["analyze_regulation.py", DB], False),
    # taxonomy_ids: taksonomi KONTROLLU SOZLUK olur (ROADMAP 135.1).
    # Uc tablo YAZAR (replicon_taxon, taxon_lineage, taxon_genus), yani
    # kendisinden sonraki her DB girdili adimi bayatlatir. Buraya konuldu
    # cunku CINS SORAN DOKUZ ADIMIN HEPSI bundan sonra geliyor
    # (promoter_architecture, characterize_leaves, analyze_ecology,
    # threshold_sensitivity, geography, build_search_index, make_report,
    # make_explorer) ve hepsi cinsi artik `taxon_genus`'tan okuyor.
    # Daha erkeye konsa yayilim buyur, daha geriye konsa tuketiciler
    # tabloyu bulamaz. Girdisi: replicon (build_db) + gbk_files + taxdump.
    ("taxonomy_ids", [f"{OUT}/taxonomy_ids.json"],
     ["python3", "taxonomy_ids.py", "--db", DB, "--gbk", "gbk_files",
      "--out", f"{OUT}/taxonomy_ids.json"],
     ["taxonomy_ids.py", DB, "gbk_files"], False),

    # replicon_kind: her kaydin NE OLDUGU (genom / transkript). Kaydin
    # kendi LOCUS beyanindan okunur, accession onekinden degil; onek
    # yalnizca capraz kontrol. `taxonomy_ids`'dan SONRA ve bu tabloyu
    # okuyan adimlardan (stats_overview, coverage_audit) ONCE olmali
    # (ROADMAP 135.7).
    ("replicon_kind", [f"{OUT}/replicon_kind.json"],
     ["python3", "replicon_kind.py", "--db", DB, "--gbk", "gbk_files",
      "--out", f"{OUT}/replicon_kind.json"],
     ["replicon_kind.py", DB, "gbk_files"], False),
    # Tablo YAZMAZ, yalnizca okur: nereye konursa konsun kimseyi
    # bayatlatmaz (bkz. family_scan'daki not).
    ("promoter_architecture", [f"{OUT}/promoter_architecture.json"],
     ["python3", "promoter_architecture.py", "--db", DB,
      "--out", f"{OUT}/promoter_architecture.json"],
     ["promoter_architecture.py", DB], False),

    # --- run_all.sh 5-12: grafikte eksik olan analiz ve DB-yazar adimlari
    ("analyze", [],
     ["python3", "analyze.py", "--db", DB, "--out-dir", OUT],
     ["analyze.py", DB], False),
    ("null_model", [f"{OUT}/null_model.csv"],
     ["python3", "null_model.py", "--db", DB, "--gbk-dir", "gbk_files",
      "--processes", "16", "--out", f"{OUT}/null_model.csv"],
     ["null_model.py", DB, "gbk_files"], True),
    ("variants", [f"{OUT}/sdp_positions.csv"],
     ["python3", "analyze_variants.py", "--db", DB, "--out-dir", OUT],
     ["analyze_variants.py", DB], False),
    ("subfamilies", [],
     ["python3", "discover_subfamilies.py", "--db", DB, "--out-dir", OUT],
     ["discover_subfamilies.py", DB], False),
    ("homogenize", [],
     ["python3", "recursive_homogenize.py", "--db", DB, "--out-dir", OUT],
     ["recursive_homogenize.py", DB], False),
    ("leaves", [f"{OUT}/leaf_profiles.csv"],
     ["python3", "characterize_leaves.py", "--db", DB, "--out-dir", OUT],
     ["characterize_leaves.py", DB], False),
    ("representatives", [],
     ["python3", "extract_representatives.py", "--db", DB, "--out-dir", OUT],
     ["extract_representatives.py", DB], False),
    ("domains", [],
     ["python3", "classify_domains.py", "--db", DB, "--out-dir", OUT],
     ["classify_domains.py", DB], False),

    # --- kuratorlu veriye bagli analizler
    ("ecology_stats", [f"{OUT}/cluster_ecology_stats.csv"],
     ["python3", "analyze_ecology.py", "--db", DB,
      "--out", f"{OUT}/cluster_ecology_stats.csv"],
     ["analyze_ecology.py", DB, ECO], False),
    ("variant_signature", [f"{OUT}/variant_signatures.csv"],
     ["python3", "variant_signature.py", "--db", DB, "--out-dir", OUT],
     ["variant_signature.py", DB], False),
    ("motif_stats", [f"{OUT}/motif_stats.json"],
     ["python3", "motif_stats.py", "--db", DB, "--out", f"{OUT}/motif_stats.json"],
     ["motif_stats.py", DB, MOTIF_HMM, "ROs_71_Clean/ROmotif71_hmmaln.sto"], False),
    ("operon_validation", [f"{OUT}/operon_validation.json"],
     ["python3", "operon_validation.py", "--db", DB,
      "--out", f"{OUT}/operon_validation.json"],
     ["operon_validation.py", DB], False),
    # reference_redundancy: dogrulama kapisi bu dosyayi OKUYOR ama betigi ne
    # run_all.sh'ta ne burada vardi. Sonuc: reference_redundancy.json 6 Ekim'den
    # kalmisti (71 referans) ve kapi onu 75'lik FASTA ile karsilastirip uc
    # kontrolde FAIL verdi. Artik grafikte (2026-10-16).
    ("reference_redundancy", [f"{OUT}/reference_redundancy.json"],
     ["python3", "reference_redundancy.py", "--db", DB, "--refs", REFS,
      "--chemistry", CHEM, "--decisions", DECISIONS,
      "--out", f"{OUT}/reference_redundancy.json"],
     ["reference_redundancy.py", DB, CHEM, REFS, DECISIONS], False),
    ("redundancy", [f"{OUT}/redundancy.json"],
     ["python3", "redundancy.py", "--db", DB, "--out", f"{OUT}/redundancy.json"],
     ["redundancy.py", DB], False),
    ("cooccurrence", [f"{OUT}/cooccurrence.json"],
     ["python3", "cooccurrence.py", "--db", DB, "--out", f"{OUT}/cooccurrence.json"],
     ["cooccurrence.py", DB], False),
    ("substrate_predictability", [f"{OUT}/substrate_predictability.json"],
     ["python3", "substrate_predictability.py", "--pairs", f"{OUT}/reference_pairs.csv",
      "--ecology", ECO, "--out", f"{OUT}/substrate_predictability.json"],
     ["substrate_predictability.py", f"{OUT}/reference_pairs.csv", ECO], False),
    ("habitat", [f"{OUT}/habitat.json"],
     ["python3", "isolation_source.py", "--db", DB, "--gbk-dir", "gbk_files",
      "--ecology", ECO, "--chemistry", CHEM, "--out", f"{OUT}/habitat.json"],
     ["isolation_source.py", DB, ECO, CHEM, "gbk_files"], False),
    ("threshold_sensitivity", [f"{OUT}/threshold_sensitivity.json"],
     ["python3", "threshold_sensitivity.py", "--db", DB,
      "--out", f"{OUT}/threshold_sensitivity.json"],
     ["threshold_sensitivity.py", DB, ECO], False),
    ("stratified_stats", [f"{OUT}/stratified_stats.json", f"{OUT}/carboxylate_rows.json"],
     ["python3", "stratified_stats.py", "--db", DB, "--out-dir", OUT],
     ["stratified_stats.py", DB, ECO], False),

    # --- bu oturumda yazilanlar; run_all.sh'de HIC gecmiyorlardi
    ("geography", [f"{OUT}/geography.json"],
     ["python3", "geography.py", "--db", DB, "--out", f"{OUT}/geography.json"],
     ["geography.py", "stats_overview.py", DB], False),
    ("control_elements", [f"{OUT}/control_elements.json"],
     ["python3", "control_elements.py", "--db", DB, "--chemistry", CHEM,
      "--out", f"{OUT}/control_elements.json"],
     ["control_elements.py", "stats_overview.py", DB, CHEM], False),
    ("ecological_origin", [f"{OUT}/ecological_origin.json"],
     ["python3", "ecological_origin.py", "--db", DB, "--out-dir", OUT],
     ["ecological_origin.py", DB, ECO], False),
    ("ancestry", [f"{OUT}/ancestry.json"],
     ["python3", "ancestry.py", "--db", DB, "--ecology", ECO,
      "--out", f"{OUT}/ancestry.json"],
     ["ancestry.py", DB, ECO], False),
    ("ferredoxin_residue", [f"{OUT}/ferredoxin_residue.json"],
     ["python3", "ferredoxin_residue.py", "--db", DB,
      "--out", f"{OUT}/ferredoxin_residue.json"],
     ["ferredoxin_residue.py", "stats_overview.py", DB], False),
    ("empty_types", [f"{OUT}/empty_types.json"],
     ["python3", "empty_types.py", "--db", DB, "--out", f"{OUT}/empty_types.json"],
     ["empty_types.py", DB, CHEM], False),
    ("completeness_audit", [f"{OUT}/completeness_audit.json"],
     ["python3", "completeness_audit.py", "--db", DB, "--chemistry", CHEM,
      "--out", f"{OUT}/completeness_audit.json", "--alphafold-probe", "none"],
     ["completeness_audit.py", DB, CHEM], False),
    ("reference_structures", [f"{OUT}/reference_structures.json"],
     ["python3", "verify_reference_structures.py", "--chemistry", CHEM,
      "--out", f"{OUT}/reference_structures.json", "--offline"],
     ["verify_reference_structures.py", CHEM, "reference_structures.csv"], False),
    ("operon_relations", [f"{OUT}/operon_relations.json"],
     ["python3", "operon_relations.py", "--db", DB, "--out-dir", OUT],
     ["operon_relations.py", DB, CHEM, ECO], True),
    ("variant_and_operon", [f"{OUT}/variant_and_operon.json"],
     ["python3", "variant_and_operon.py", "--db", DB, "--out-dir", OUT],
     ["variant_and_operon.py", DB], True),
    ("active_site", [f"{OUT}/active_site.json"],
     ["python3", "active_site.py", "--db", DB, "--out-dir", OUT, "--offline"],
     ["active_site.py", DB, CHEM], True),
    ("learning", [f"{OUT}/learning.json"],
     ["python3", "learn_from_data.py", "--db", DB, "--out-dir", OUT],
     ["learn_from_data.py", DB, CHEM, ECO], True),

    # --- triad kapisi: uyelik kesinlestikten SONRA, ozetlerden ONCE ---------
    # Bu uc adim eskiden ELLE kosturuluyordu ve zincirde yoktu. Sonucu:
    # yazdiklari tablolar (ro_triad, ro_triad_meta, ro_protein) hicbir adimin
    # `db_writes` kumesinde olmadigi icin hicbirinin `ignore`una girmiyor, ve
    # `roar.sqlite`i girdi sayan HER adim onlari "degismis ust girdi"
    # saniyordu. 2026-10-30'da olculdu: 51 artefaktin 45'i bu yuzden sahte
    # BAYAT gorunuyordu, oysa 52 adim betiginden yalnizca validate_curation.py
    # bu tablolari OKUYOR. Boylece `rebuild.py --check` yayin oncesi kapi
    # olarak ise yaramaz hale geliyordu: 45 adim her zaman bayat diyorsa
    # gercek bayatligi sahtesinden ayirt edemezsin, ve zincir 42/47'de
    # validation_gate'te duruyordu.
    #
    # Sira onemli: kapi onaylanmis kumeyi okur, yani uyeligi belirleyen
    # adimlardan (annotate ve sonrasi) SONRA gelmeli; validation_gate ise bu
    # tablolari okudugu icin BUNLARDAN sonra gelmeli -- zaten oyle.
    ("carboxylate_audit", [f"{OUT}/carboxylate_audit.json"],
     ["python3", "carboxylate_audit.py", "--db", DB, "--out",
      f"{OUT}/carboxylate_audit.json", "--work", f"{SCRATCH}/roar_carbox_work"],
     ["carboxylate_audit.py", DB, MOTIF_HMM], True),
    # Tek gecis yetmez: anchor ro_triad'i OKUR, gate onu YAZAR. Sabit noktaya
    # kadar donen surucu triad_converge.py; gerekcesi o dosyanin basinda.
    ("triad_gate", [f"{OUT}/triad_gate.json", f"{OUT}/triad_anchor_his.json"],
     ["python3", "triad_converge.py", "--db", DB,
      "--work", f"{SCRATCH}/roar_triad_anchor"],
     ["triad_converge.py", "triad_gate.py", "triad_anchor_his.py", DB,
      f"{OUT}/carboxylate_audit.json"], True),
    ("protein_entries", [],
     ["python3", "ingest_protein_entries.py", "--db", DB,
      "--work", f"{SCRATCH}/roar_protein"],
     ["ingest_protein_entries.py", DB, MOTIF_HMM], True),
    # --- ozet katmani: HER SEYDEN sonra
    # tree_meta.json BEYAN EDILIR: sayfanin altyazisi kullanilan modeli ve
    # destek turunu buradan okuyor, yani dosya tazelik takibine girmek
    # zorunda -- yoksa altyazi eski bir modelin adini yayinlayabilir.
    ("phylogeny", [f"{OUT}/tree_all.nwk", f"{OUT}/tree_meta.json", f"{OUT}/cluster_assignment_conflicts.csv", f"{OUT}/cluster_identity_matrix.csv", f"{OUT}/cluster_reference_identity.csv", f"{OUT}/ssn_edges.csv", f"{OUT}/ssn_nodes.csv"],
     ["python3", "build_phylogeny.py", "--db", DB, "--out-dir", OUT,
      "--fasttree", "FastTree", "--work", f"{SCRATCH}/roar_phylo_work"],
     ["build_phylogeny.py", DB, MOTIF_HMM], False),
    # Gruplarin agacta dal olusturup olusturmadigi OLCULUR. Filogeni sayfasi
    # "this tree recovers the historical RO groups from sequence alone" diye
    # olculmemis bir sonuc yayimliyordu; 2026-11-01'de olculdu ve iddia
    # yalnizca grup 5 ve 6 icin dogru cikti. Ciktisi tazelik takibinde, cunku
    # sayfanin cumlesi artik bu dosyadan okunuyor.
    ("group_monophyly", [f"{OUT}/group_monophyly.json"],
     ["python3", "group_monophyly.py", "--db", DB, "--out-dir", OUT],
     ["group_monophyly.py", DB, f"{OUT}/tree_all.nwk"], False),
    # UniProt anlik goruntusunun BUGUNKU durumu. AGIR cunku 737 ag sorgusu.
    # Silinmis kayitlar atilmaz, ayri tabloya durum yazilir -- "silinmis
    # kayit" ile "veri yok" ayri seylerdir.
    ("protein_freshness", [f"{OUT}/protein_freshness.json"],
     ["python3", "protein_freshness.py", "--db", DB, "--out-dir", OUT],
     ["protein_freshness.py", DB], True),
    # Bootstrap destegi: SH-benzeri yerel destegin tasiyamadigi klad
    # iddialarini sinar. AGIR cunku her tekrar icin agac yeniden cikarilir.
    ("bootstrap_references", [f"{OUT}/bootstrap_references.json",
                              f"{OUT}/bootstrap_references.nwk"],
     ["python3", "bootstrap_tree.py", "--set", "references",
      "--replicates", "200", "--db", DB, "--out-dir", OUT],
     ["bootstrap_tree.py", DB, MOTIF_HMM], True),
    ("bootstrap_confident", [f"{OUT}/bootstrap_confident.json",
                             f"{OUT}/bootstrap_confident.nwk"],
     ["python3", "bootstrap_tree.py", "--set", "confident",
      "--replicates", "100", "--db", DB, "--out-dir", OUT],
     ["bootstrap_tree.py", DB, MOTIF_HMM], True),
    # SIRALAMA ONEMLI: bu uc adim phylogeny'den SONRA gelir. Once oncesine
    # konmuslardi ve tazelik izleyicisi haklı olarak agaci BAYAT saydi --
    # "adim i'nin DB girdisi, kendisinden oncekilerin yazdiklaridir" kurali
    # geregi. Oysa phylogeny ro_union'i hic okumuyor; sonuc iki saatlik bir
    # agac isinin bosuna tetiklenmesi olacakti.
    # RefSeq'in TUM protein uzayi. AGIR: 871 dosya indirilir (~36 GB akitarak,
    # diskte birikmez) ve ~143 milyon protein taranir. Madde 95'in yazdigi
    # boslugu kapatir -- protein evrenimiz UniProt'tan gelmisti, NCBI'dan
    # DEGIL. Ciktisi slot durum dosyalari; kesilirse kaldigi yerden devam
    # eder, surum degisirse kendini sifirlar.
    ("ncbi_scan", [f"{NCBI}/hits_slot0.faa"],
     ["python3", "scan_ncbi_proteins.py", "--work", NCBI, "--parallel", "6"],
     ["scan_ncbi_proteins.py", MOTIF_HMM], True),
    # Taramanin adaylarini KARARA baglar. Tip atamasi ro_filter ile (UniProt
    # tarafiyla ayni kural), triad kapisi ingest_protein_entries'in okumasiyla.
    # Kendi kaynak tablosunu yazar; birlesimi build_union.py kurar.
    ("ncbi_hits", [],
     ["python3", "ingest_ncbi_hits.py", "--db", DB, "--scan-dir", NCBI],
     ["ingest_ncbi_hits.py", "ro_filter.py", "ingest_protein_entries.py",
      MOTIF_HMM, f"{NCBI}/hits_slot0.faa"], True),
    # Birlesim: her kaynaktan gelen adaylarin tek tabloda toplanmasi. ro ve
    # ro_protein'den SONRA gelmek zorunda, cunku ikisini de okur. Anasayfa
    # kesisim/birlesim ayrimini bu tablodan yaziyor; bayat kalirsa iki sayi
    # sessizce ayrisir ve validate_curation check_union_denominator FAIL
    # verir -- ama once burada tazelenmesi gerekir.
    ("union", [],
     ["python3", "build_union.py", "--db", DB, "--fasta", "ro_alpha.fasta"],
     ["build_union.py", DB, "ro_alpha.fasta"], False),
    ("stats", [f"{OUT}/stats.json"],
     ["python3", "stats_overview.py", "--db", DB, "--out", f"{OUT}/stats.json"],
     ["stats_overview.py", DB, CHEM, ECO], False),
    ("search_index", [],
     ["python3", "build_search_index.py", "--db", DB, "--ecology", ECO,
      "--chemistry", CHEM],
     ["build_search_index.py", DB, ECO, CHEM], False),
    # disagreements.py yalnizca baska JSON'lari OKUR, bu yuzden en sonda ve
    # girdileri o dosyalardir.
    ("disagreements", [f"{OUT}/disagreements.json"],
     ["python3", "disagreements.py", "--analysis-dir", OUT,
      "--out", f"{OUT}/disagreements.json"],
     ["disagreements.py", f"{OUT}/stats.json", f"{OUT}/operon_relations.json",
      f"{OUT}/ferredoxin_residue.json"], False),
    # --- run_all.sh 14: yerel HTML sayfalar
    ("report_html", ["report.html"],
     ["python3", "make_report.py", "--db", DB, "-o", "report.html"],
     ["make_report.py", DB], False),
    ("explorer_html", ["explorer.html"],
     ["python3", "make_explorer.py", "--db", DB, "-o", "explorer.html"],
     ["make_explorer.py", DB], False),
    ("hub_html", ["hub.html"],
     ["python3", "make_hub.py", "--db", DB, "--out-dir", OUT, "-o", "hub.html"],
     ["make_hub.py", DB], False),
    # --- 2026-11-05: YAYIMLANAN AMA GRAFIKTE OLMAYAN artefaktlar.
    # Olculdu: atlas'in indirmeye sundugu 56 dosyanin 25'i hicbir adimin
    # ciktisi degildi, yani `rebuild.py --run` onlari ASLA tazeleyemezdi ve
    # tazelik kapisi onlari YARGILAYAMIYORDU. Kanit: carboxylate_method_card
    # 30 Ekim'deki denetimden uretilmis ve kalmisti, oysa denetim 31 Ekim'de
    # yenilenmisti -- yayimlanan kart bayatti ve bunu kimse goremiyordu.
    #
    # AGIR OLANLAR DIS SERVISE ISTEK ATIYOR (ChEBI, Rhea, UniProt, NCBI,
    # PDB). Depodaki konvansiyon bunu zaten `heavy` ile isaretliyor
    # ("AGIR cunku 737 ag sorgusu", yukarida). Yani `--run` onlari
    # KENDILIGINDEN kosturmaz; dis servis ancak bilerek sorgulanir.
    ("a5_match_guard", [f"{OUT}/a5_match_guard.json"],
     ["python3", "a5_match_guard.py", "--db", DB],
     # Betigin GERCEKTEN okudugu cikarim dizinleri de girdi: a5_match_guard.py
     # icindeki varsayilanlarla ayni. Beyan edilmezse bu dizinler degistiginde
     # adim bayat gorunmez (dizin ozeti 2026-11-05'te eklendi, bak dir_digest).
     ["a5_match_guard.py", DB, "cohorts/a5_protein_recovery.tsv",
      "/media/lin-bio/back2/roar_gc_a5", "/media/lin-bio/back2/roar_gc_wide"], True),
    ("carboxylate_anchor_holdout", [f"{OUT}/carboxylate_anchor_holdout.json"],
     ["python3", "carboxylate_anchor_holdout.py", "--db", DB],
     ["carboxylate_anchor_holdout.py", DB, "genomic_context/cand_aln.sto",
      "ROs_71_Clean/refs71.fasta"], True),
    ("carboxylate_method_card", [f"{OUT}/carboxylate_method_card.txt"],
     ["python3", "carboxylate_method_card.py",
      "--out", f"{OUT}/carboxylate_method_card.txt"],
     ["carboxylate_method_card.py", f"{OUT}/carboxylate_audit.json"], False),
    ("oxygen_fate_computed", [f"{OUT}/oxygen_fate_computed.json"],
     ["python3", "compute_oxygen_fate.py"],
     ["compute_oxygen_fate.py", "types", "oxygen_fate.csv"], False),
    ("reaction_anchors", [f"{OUT}/reaction_anchors.json"],
     ["python3", "anchor_reactions.py"],
     ["anchor_reactions.py", "types"], True),
    ("reaction_verification", [f"{OUT}/reaction_verification.json"],
     ["python3", "verify_reactions.py"],
     ["verify_reactions.py", "types"], True),
    ("curation_claims", [f"{OUT}/curation_claims.json"],
     ["python3", "verify_curation_claims.py"],
     ["verify_curation_claims.py", "types"], True),
    ("bridging_geometry_audit", [f"{OUT}/bridging_geometry_audit.json"],
     ["python3", "bridging_geometry_audit.py"],
     ["bridging_geometry_audit.py"], True),

    # Komsu geni URUN ADINDAN degil Pfam modelinden siniflandirir: hareketli
    # eleman, rekombinaz ve duzenleyici ailesi. Uc AYRI kutuphane kullanir ve
    # hicbiri component_hmms/components.hmm ile karistirilmaz.
    #
    # NEDEN ZINCIRIN SONUNDA. Yeni tablo yazan bir adim nereye konursa
    # KENDISINDEN SONRAKI her adimin DB girdisi degisir ve hepsi bir kez
    # bayat gorunur. Bu adimin ciktisini su an okuyan baska adim yok;
    # sona konunca kimseyi bayatlatmaz. Bir tuketici eklenirse adim onun
    # ONUNE tasinmali -- yoksa tuketici eski tabloyu okur.
    ("family_hmms", ["family_hmms/mobility.hmm", "family_hmms/recombinase.hmm",
                     "family_hmms/regulator.hmm",
                     "family_hmms/mobility.manifest.json",
                     "family_hmms/recombinase.manifest.json",
                     "family_hmms/regulator.manifest.json"],
     ["python3", "fetch_family_hmms.py", "--library", "all"],
     ["fetch_family_hmms.py"], True),
    #  substrate_chebi: substrat icin KONTROLLU kimlik. `--lookup-names`
    #  BURADA YOK -- ag erisimi bir kapi adiminin kosulu olamaz; onbellek
    #  yetkili ve ayni tabloyu cevrimdisi uretiyor (ROADMAP 135.3).
    ("substrate_chebi", [f"{OUT}/substrate_chebi.json"],
     ["python3", "substrate_chebi.py", "--db", DB, "--types", "types",
      "--out", f"{OUT}/substrate_chebi.json"],
     ["substrate_chebi.py", DB, "types",
      f"{OUT}/oxygen_fate_computed.json"], False),
    #  coverage_audit: dis evrene gore KAPSAMA acigi. `--fetch` BURADA
    #  YOK -- ag erisimi bir kapi adiminin kosulu olamaz (bkz.
    #  substrate_chebi); onbellek yetkili. Girdisi onbellek + DB, cikti
    #  bir rapor; tablo yazmaz, kimseyi bayatlatmaz (ROADMAP 135.9).
    ("coverage_audit", [f"{OUT}/coverage_audit.json"],
     ["python3", "coverage_audit.py", "--db", DB,
      "--out", f"{OUT}/coverage_audit.json"],
     ["coverage_audit.py", DB, f"{OUT}/pf00848_members.json.gz"], False),
    #  gff3: genomik baglam STANDART bicimde. Girdisi DB oldugu icin
    #  tablolar degisince bayat gorunur; kendi ciktisini ayristirip
    #  dogruladigi icin gecersiz bir degisim dosyasi yayimlanamaz.
    ("gff3", [f"{OUT}/roar_genomic_context.gff3.gz",
              f"{OUT}/roar_genomic_context.json"],
     ["python3", "export_gff3.py", "--db", DB,
      "--out", f"{OUT}/roar_genomic_context.gff3.gz"],
     ["export_gff3.py", DB], False),
    #  hmm_pins: hangi modelin hangi SURUMU yayimlanan isabetleri uretti.
    #  Kutuphane dosyalarini girdi olarak alir, yani bir model degisirse
    #  sabitleme bayat gorunur -- model kutuphanesi sessizce degisemez.
    #  `--check-upstream` KAPIDA KOSMAZ: ag erisimi gerektirir ve bir kapi
    #  dis bir servisin erisilebilirligine bagli olamaz.
    ("hmm_pins", [f"{OUT}/hmm_pins.json"],
     ["python3", "pin_hmm_libraries.py", "--out", f"{OUT}/hmm_pins.json"],
     ["pin_hmm_libraries.py", "component_hmms/components.hmm",
      "family_hmms/mobility.hmm", "family_hmms/recombinase.hmm",
      "family_hmms/regulator.hmm", "ROs_71_Clean/ROmotif71.hmm",
      "RieskeDB71.hmm"], False),
    ("family_scan", [f"{OUT}/family_scan.json"],
     ["python3", "family_scan.py", "--db", DB, "--out", f"{OUT}/family_scan.json"],
     ["family_scan.py", DB, "family_hmms/mobility.hmm",
      "family_hmms/recombinase.hmm", "family_hmms/regulator.hmm"], False),

    # cross_references: A6. Girdisi ro_union_source + ro_protein_status, yani
    # union'dan SONRA. YAPRAK olarak en sona konuldu: `ro_xref` tablosunu
    # yazdigi icin zincirin ortasinda dururken kendisinden SONRAKI her DB
    # girdili adim (stats, family_scan, a5_match_guard,
    # carboxylate_anchor_holdout) onu "degismis ust girdi" sayip bayatliyordu
    # -- ucu AGIR. Bu tabloyu okuyan tek sey webapp, yani adimin burada
    # durmasinin hicbir bedeli yok.
    ("cross_references", [f"{OUT}/cross_references.json"],
     ["python3", "cross_references.py", "--db", DB,
      "--out", f"{OUT}/cross_references.json"],
     ["cross_references.py", DB], False),

    # YAYIMLANAN SITENIN KENDISI. 2026-11-09'a kadar grafigin disindaydi:
    # her sey taze olsa bile canlidaki kopyanin eski olup olmadigini hicbir
    # sey sormuyordu. Olculen bedel: gh-pages 2026-11-02'de kalmisti ve
    # bilesen duzeltmesi oncesi operon sayilarini gosteriyordu (beta 2.523
    # yerine 4.372 olmaliydi) -- yerel her dosya guncelken.
    #
    #  KAPI EN SONDA, site_freeze'den hemen once. Onceki yerinde (stats'tan
    #  hemen sonra) grafigin 15'ten fazla adimi KAPIDAN SONRA yayimlanan
    #  dosya uretiyordu, ve kapinin tazelik kontrolu yayimlanan TUM kumeyi
    #  yargiliyor. Yani kapi, kendisinden sonraki bir adim bayatken asla
    #  yesil donemiyordu -- tam kosuda bile, cunku o adimlar henuz
    #  kosmamisti. Tasimak guvenli: `validation_gate`'in cikti listesi BOS
    #  ve hicbir adim `provenance.json`'u girdi olarak okumuyor.
    ("validation_gate", [],
     ["python3", "validate_curation.py", "--db", DB, "--ecology", ECO,
      "--chemistry", CHEM, "--out-dir", OUT],
     ["validate_curation.py", DB, ECO, CHEM, "ROs_71_Clean/ROmotif71_hmmaln.sto"], False),
    ("provenance", [f"{OUT}/provenance.json"],
     ["python3", "provenance.py", "--db", DB, "--out-dir", OUT],
     ["provenance.py", DB], False),

    # AGIR: ~16 dakika, 3.282 statik sayfa. Hafif zincirde kosmaz; isi
    # bayatligi GORUNUR kilmak, kendiliginden yayim yapmak degil. Yayim
    # sirasi degismedi: check_site.py -> freeze.py -> deploy_pages.sh.
    ("site_freeze", ["webapp/site/index.html", "webapp/site/search_index.json",
                     "webapp/site/download/MANIFEST.txt"],
     ["python3", "webapp/freeze.py", "--out", "webapp/site",
      "--base", "/pro-sim-blast"],
     ["webapp/freeze.py", "webapp/app.py", "webapp/atlas.py",
      "webapp/build_entry_data.py", "webapp/templates", "webapp/static",
      DB, OUT], True),
]


# Veritabani OZEL. mtime ile karsilastirmak ise yaramaz: pipeline'in uc adimi
# (evidence_tiers, etc_types, analyze_regulation) roar.sqlite'a YAZAR, yani her
# kosuda mtime ileri gider ve ondan sonraki her analiz sonsuza dek bayat
# gorunur -- olculdu, --run iki kez ust uste kosuldugunda ikincisi de her seyi
# yeniden uretiyordu. Dosya baytlari da ise yaramaz, cunku SQLite ayni mantiksal
# icerigi farkli sayfa duzeniyle yazabilir.
#
# Bunun yerine MANTIKSAL parmak izi: onemli tablolarin satir sayilari. Pipeline
# ayni veriyi yeniden turetirse parmak izi degismez ve hicbir sey bayatlamaz;
# veri gercekten degisirse degisir.
FINGERPRINT_TABLES = ["ro", "replicon", "neighbor", "operon", "operon_gene",
                      "ro_evidence", "ro_etc", "ro_regulation", "ro_leaf",
                      "neighbor_protein", "gene_category", "replicon_source"]


def db_fingerprint(path):
    import sqlite3
    try:
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    except sqlite3.Error:
        return None
    parts = []
    try:
        have = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        for t in FINGERPRINT_TABLES:
            if t in have:
                n = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
                parts.append(f"{t}={n}")
        # dogrulanmis RO sayisi ve tip sayisi: kuratorluk degisikligi bunlari
        # dogrudan oynatir
        if "ro" in have:
            parts.append("confirmed=%d" % con.execute(
                "SELECT COUNT(*) FROM ro WHERE is_confirmed=1").fetchone()[0])
            parts.append("types=%d" % con.execute(
                "SELECT COUNT(DISTINCT ro_cluster) FROM ro "
                "WHERE is_confirmed=1").fetchone()[0])
    except sqlite3.Error:
        return None
    finally:
        con.close()
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# TABLO BAZINDA DB GIRDISI (2026-10-14)
#
# Eski davranis: bir adimin DB girdisi tum DB'nin parmak iziydi. Ama DB'ye
# yazan adimlar ZINCIRIN ICINDE: annotate'ten sonra operons, evidence_tiers,
# etc_types... hepsi tablo yaziyor. Sonuc: annotate (40 dk) ve operons her
# --run'da yeniden "bayat" goruluyordu, cunku kendisinden SONRA gelenlerin
# yazdigi tablolar parmak izini oynatiyordu. Olculdu: tam bir kosudan hemen
# sonra --plan yine annotate'ten basliyordu.
#
# Simdi: her adimin DB girdisi, yalnizca KENDISINDEN ONCEKI adimlarin (ve
# temel kurulumun) yazdigi tablolardir. Bir adimin yazdigi tablolar iki
# kaynaktan birlestirilir:
#   statik : betigin kaynagindaki CREATE/INSERT/UPDATE/ALTER/DELETE/DROP
#   dinamik: kosu sirasinda olculen, oncesi/sonrasi sayisi degisen tablolar
# ---------------------------------------------------------------------------
import re as _re

_WRITE_RX = _re.compile(
    r"(?:CREATE\s+(?:VIRTUAL\s+)?TABLE(?:\s+IF\s+NOT\s+EXISTS)?|"
    r"INSERT(?:\s+OR\s+\w+)?\s+INTO|UPDATE|ALTER\s+TABLE|DELETE\s+FROM|"
    r"DROP\s+TABLE(?:\s+IF\s+EXISTS)?)\s+[\"`]?([A-Za-z_][A-Za-z0-9_]*)",
    _re.IGNORECASE)
_SQL_WORDS = {"set", "or", "if", "into", "from", "table", "where", "select"}


# Bir adimin yazdigi tablolar normalde IKI yoldan ogreniliyor: betigin
# kaynagini tarayarak (static_writes) ve adim kostuktan sonra tablo
# sayimlarini karsilastirarak (db_writes). 2026-10-31'de olculdu: bu uc adimda
# IKISI DE SESSIZCE BASARISIZ OLUYOR.
#
#   - Kaynak taramasi: `triad_gate` adiminin komutu triad_converge.py, ve o
#     surucu tabloyu KENDISI yazmiyor, triad_gate.py'yi cagiriyor. Tarayici
#     surucunun icinden gecemiyor.
#   - Sayim karsilastirmasi: iki betik de idempotent. triad_gate ro_triad'i
#     ayni 19.125 satirla, ingest_protein_entries ro_protein'i ayni 73.699
#     satirla yeniden yaziyor. Satir sayisi degismedigi icin "yazma
#     gozlemlenmedi" deniyor ve db_writes BOS kaliyor.
#
# Sonucu agir: bu tablolar hicbir adimin ignore kumesine girmiyor, ve
# roar.sqlite'i girdi sayan HER onceki adim onlari "degismis ust girdi"
# saniyor -- 48 artefaktin bayat gorunmesinin sebebi bu. Tazelik kapisi
# boylece ise yaramaz hale geliyor. O yuzden burada BEYAN ediliyorlar.
# Beyan, gozlemin yerine gecmez; gozlemle BIRLESTIRILIR (writes[name]).
DECLARED_WRITES = {
    "triad_gate": {"ro_triad", "ro_triad_meta", "ro_membership"},
    "protein_entries": {"ro_protein"},
    "ncbi_hits": {"ro_ncbi_hit"},
    "union": {"ro_union", "ro_union_source"},
    "taxonomy_ids": {"replicon_taxon", "taxon_lineage", "taxon_genus"},
    "substrate_chebi": {"type_substrate_chebi"},
    "replicon_kind": {"replicon_kind"},
    "db_repair": {"replicon"},
}


#  OKUDUGU TABLOLAR. Bir adimin DB girdisi, "kendisinden once yazilan TUM
#  tablolar" olarak hesaplaniyordu. Bu asiri genis, ve uc kez ayni yanlis
#  alarmi uretti: zincire tablo yazan yeni bir adim eklendiginde, o
#  tablolari HIC OKUMAYAN onlarca adim bayat gorundu (2026-11-01
#  cross_references'ta 4 adim; 2026-11-11 taxonomy_ids'te 32 adim, 8'i
#  AGIR). Saatlerce yeniden hesap, birebir ayni cikti -- ve dogrulama
#  kapisi yanlis yere FAIL verdi.
#
#  Bir adim, adini HIC ANMADIGI bir tablodan etkilenemez. O yuzden
#  karsilastirma, adimin kaynaginda (ve yerel olarak import ettigi
#  modullerin kaynaginda) adi GECEN tablolarla sinirlanir. Kapanis
#  gereklidir: cins tuketicileri `taxon_genus`'u kendi SQL'lerinde
#  anmiyor, `ro_taxonomy`yi import edip ona soruyorlar.
#
#  KURAL ZAYIFLATILMIYOR, cunku belirsizlikte TUM tablolar karsilastirilir:
#    - kaynak okunamiyorsa,
#    - kaynak tablo adini DINAMIK kuruyorsa (f-string ya da dizgi toplama
#      icinde FROM/JOIN/UPDATE/INTO).
#  Yani kacirma yonunde degil, fazla sayma yonunde hata yapar.
#
#  `sqlite_master` TEK BASINA belirsizlik SAYILMAZ, ve bunun sebebi var:
#  depodaki kullanimlarin tamami VARLIK KONTROLU -- ya `WHERE name=?` ile
#  tek tabloyu soruyor, ya da `WHERE type='table'` ile adlari toplayip
#  kaynakta HARFI HARFINE yazili adlara karsi uyelik bakiyor. Iki durumda
#  da tablo adi kaynakta geciyor, yani okuma gorunur. SQL'de tanimlayici
#  parametre olamaz; dolayisiyla adi gezip o tabloyu gercekten OKUYAN bir
#  betik, adi dizgiyle kurmak ZORUNDA -- ve o durumu asagidaki desen
#  yakaliyor. Olcum: bu ayrim 18 adimi yanlis alarmdan cikardi.
_IDENT_RX = _re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_DYNAMIC_SQL_RX = _re.compile(
    r"(?:FROM|JOIN|INTO|UPDATE|TABLE)\s*(?:[\"']?\{|[\"']\s*\+|%s)",
    _re.IGNORECASE)


def _local_modules():
    """Depo kokundeki ve webapp/ altindaki modul adi -> dosya yolu.

    Bu iki dizin, depodaki TUM import edilebilir modulleri kapsiyor
    (`find . -name '*.py'`: 129 kok + 7 webapp; otekiler test ve cikti
    dizinleri, hicbir adim onlari import etmiyor). Yeni bir paket dizini
    eklenirse buraya da eklenmek ZORUNDA, yoksa o modulun okudugu
    tablolar gorunmez olur -- ve gorunmez bir okuma, bayatligi KACIRIR.
    """
    out = {}
    for pattern in ("*.py", os.path.join("webapp", "*.py")):
        for path in glob.glob(pattern):
            out.setdefault(os.path.basename(path)[:-3], path)
    return out


@functools.lru_cache(maxsize=None)
def source_closure(script):
    """Betigin kaynagi + yerel olarak import ettigi modullerin kaynagi.

    Doner: (metin, kesin_mi). `kesin_mi` False ise cagiran TUM tablolari
    karsilastirmak zorundadir.
    """
    modules = _local_modules()
    seen, stack, chunks, certain = set(), [script], [], True
    while stack:
        path = stack.pop()
        if path in seen:
            continue
        seen.add(path)
        try:
            text = open(path, encoding="utf-8", errors="replace").read()
        except OSError:
            return "", False
        chunks.append(text)
        if _DYNAMIC_SQL_RX.search(text):
            certain = False
        for match in _re.finditer(r"^\s*(?:import|from)\s+([A-Za-z_][\w.]*)",
                                  text, _re.M):
            head = match.group(1).split(".")[0]
            if head in modules:
                stack.append(modules[head])
    return "\n".join(chunks), certain


def unread_tables(script, tables):
    """`tables` icinde betigin HIC ANMADIGI olanlar (karsilastirma disi)."""
    if not script or not os.path.exists(script):
        return set()
    text, certain = source_closure(script)
    if not certain or not text:
        return set()
    words = set(_IDENT_RX.findall(text))
    return {t for t in tables if t not in words}


def static_writes(script):
    """Betigin kaynaginda yazdigi gorulen tablo adlari."""
    try:
        src = open(script).read()
    except OSError:
        return set()
    return {m.lower() for m in _WRITE_RX.findall(src)} - _SQL_WORDS


def db_counts(path):
    """Her tablonun satir sayisi + dogrulanmis RO ve tip sayisi."""
    import sqlite3
    try:
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    except sqlite3.Error:
        return {}
    out = {}
    try:
        names = [r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name NOT LIKE 'sqlite_%'")]
        for t in names:
            try:
                out[t.lower()] = con.execute(f'SELECT COUNT(*) FROM "{t}"').fetchone()[0]
            except sqlite3.Error:
                pass
        if "ro" in out:
            out["ro#confirmed"] = con.execute(
                "SELECT COUNT(*) FROM ro WHERE is_confirmed=1").fetchone()[0]
            out["ro#types"] = con.execute(
                "SELECT COUNT(DISTINCT ro_cluster) FROM ro WHERE is_confirmed=1"
            ).fetchone()[0]
    finally:
        con.close()
    return out


def _table_of(key):
    return key.split("#", 1)[0]


def db_changed(prev, cur, ignore):
    """prev ve cur sayim sozluklerini, ignore'daki tablolar HARIC karsilastir."""
    if not isinstance(prev, dict):
        return True            # eski bicim (tek dizgi): bir kez yeniden uret
    keys = (set(prev) | set(cur))
    return sorted(k for k in keys
                  if _table_of(k) not in ignore and prev.get(k) != cur.get(k))


# Bir dizinin icerik ozeti icin: bu boyutun ALTINDAKI dizinlerde dosya
# ICERIKLERI okunur, ustunde yalnizca ad+boyut. `types/` kucuk yaml'lar
# oldugu icin tam icerik, `genomic_context/` (GB'lar) ad+boyut.
DIR_CONTENT_LIMIT = 32 << 20


@functools.lru_cache(maxsize=None)
def dir_digest(path):
    """DIZIN OZETI (2026-11-05).

    NEDEN. `digest()` bir dizinde `open()` ile OSError aliyor ve None
    donuyordu. None ise karsilastirma onu ATLIYOR
    (`if not prev.get(p): continue`). Ustelik bir dizinin mtime'i icindeki
    bir dosya DUZENLENINCE degismez -- yalnizca dosya eklenip silinince
    degisir. Ikisi birlesince su cikti: `types/` iki adimda girdi olarak
    BEYAN EDILMIS ama bir `types/*.yaml` duzenlendiginde HICBIR adim bayat
    gorunmuyordu. `types/` bu veritabaninin TEK KAYIT KAYNAGI; yani kayit
    kaynagi bagimlilik grafigine gorunmuyordu.

    Olculdu: state dosyasinda `"types": null`.

    mtime'a GUVENMEZ (bu dosyanin ucuncu kurali): ad + boyut + icerik ozeti.
    """
    entries, total = [], 0
    #  GRAFIGIN KENDI KAYDI dizin ozetine GIRMEZ. Olculdu (2026-11-12):
    #  `site_freeze` girdi olarak `analysis_out` dizinini beyan ediyor ve
    #  grafigin durum dosyasi (`analysis_out/.rebuild_state.json`) adim
    #  BITTIKTEN SONRA oraya yaziliyor. Yani adim kendi kaydiyla kendini
    #  bayatlatiyordu ve HICBIR kosu onu guncel yapamiyordu -- 80/81 guncel
    #  oldugunda geride kalan tek bayat buydu. Kapinin grafik sonunda olmasi
    #  gerektigiyle ayni aile: bir adimdan SONRA yazilan bir sey o adimin
    #  girdisi olamaz.
    state_file = os.path.abspath(STATE)
    for root, dirs, files in os.walk(path):
        dirs.sort()
        for fn in sorted(files):
            p = os.path.join(root, fn)
            if os.path.abspath(p) == state_file:
                continue
            try:
                size = os.path.getsize(p)
            except OSError:
                continue
            entries.append((os.path.relpath(p, path), size, p))
            total += size
    h = hashlib.sha256()
    read_content = total <= DIR_CONTENT_LIMIT
    for rel, size, p in entries:
        h.update(f"{rel}:{size}|".encode())
        if read_content:
            try:
                with open(p, "rb") as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(chunk)
            except OSError:
                h.update(b"<okunamadi>")
    return h.hexdigest()[:16]


def digest(path):
    """Icerik ozeti. Buyuk dosyalarda bas + son + boyut yeter ve hizlidir."""
    if path.endswith(".sqlite"):
        return db_counts(path)
    if os.path.isdir(path):
        return dir_digest(path)
    try:
        size = os.path.getsize(path)
    except OSError:
        return None
    h = hashlib.sha256()
    h.update(str(size).encode())
    try:
        with open(path, "rb") as fh:
            h.update(fh.read(1 << 20))
            if size > (2 << 20):
                fh.seek(-(1 << 20), os.SEEK_END)
                h.update(fh.read())
    except OSError:
        return None
    return h.hexdigest()[:16]


def load_state():
    try:
        with open(STATE) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def save_state(state):
    os.makedirs(os.path.dirname(STATE), exist_ok=True)
    with open(STATE, "w") as fh:
        json.dump(state, fh, indent=1, sort_keys=True)


def staleness(name, outputs, inputs, state, ignore=None, counts=None):
    """Bayatlik sebebi, ya da guncelse None."""
    missing_out = [p for p in outputs if not os.path.exists(p)]
    if missing_out:
        return f"cikti yok: {', '.join(missing_out)}"
    # BOS CIKTI "GUNCEL" SAYILAMAZ. 2026-11-01: `phylogeny` adimi FastTree'yi
    # kosarken dis bir timeout tarafindan oldurüldu ve `tree_all.nwk` 0 bayt
    # kaldi. Dosya VARDI ve mtime'i TAZEYDI, yani bu islev adimi guncel
    # sayiyor; zincir ilerliyor ve bir sonraki adim (group_monophyly) bos
    # dosyayi ayristirmaya calisip IndexError ile patliyor. Kirilan yerde
    # sebebi gorunmuyor, cunku hata agaci URETEN adimda degil OKUYAN adimda
    # cikiyor. Beyan edilen bir cikti 0 bayt ise URETILMEMIS sayilir.
    empty_out = [p for p in outputs
                 if not os.path.isdir(p) and os.path.getsize(p) == 0]
    if empty_out:
        return f"cikti BOS (0 bayt, yarim kalmis uretim): {', '.join(empty_out)}"
    missing_in = [p for p in inputs if not os.path.exists(p)]
    if missing_in:
        return f"girdi yok: {', '.join(missing_in)}"
    prev = (state.get(name) or {}).get("inputs") or {}
    if outputs:
        oldest_out = min(os.path.getmtime(p) for p in outputs)
        # Veritabani mtime karsilastirmasinin DISINDA; yukaridaki nota bak.
        #
        # MTIME TETIKLER, OZET KARAR VERIR. Eskiden bu dal mtime'a bakip
        # DONUYORDU, yani asagidaki icerik karsilastirmasina hic varilmiyordu.
        # Oysa bu dosyanin kendi basligi ucuncu kurali tam bu durum icin
        # koymus: "git checkout ve kopyalama mtime'i ileri atar ama icerik
        # ayni olabilir". 2026-10-31'de olculdu: ROs_71_Clean/ROmotif71.hmm
        # bir cikti dosyasindan yeni gorunuyor ama ozeti annotate'in
        # kaydettiginin AYNISI (6fc01b561ae6cb7e), ve bu tek sahte kok
        # annotate'ten baslayip 48 artefakti bayat gosteriyordu -- saatlerce
        # yeniden hesap, birebir ayni cikti.
        #
        # Kural ZAYIFLATILMADI: taban kaydi olmayan ya da ozeti DEGISEN bir
        # girdi hala bayatlik sayilir. Yalnizca "yeni ama icerigi kanitlanmis
        # sekilde ayni" durumu disari alindi.
        newer = []
        for p in inputs:
            if p.endswith(".sqlite") or os.path.getmtime(p) <= oldest_out:
                continue
            baseline = prev.get(p)
            if baseline and baseline == digest(p):
                continue        # mtime ileri atmis, icerik ayni
            newer.append(p)
        if newer:
            return f"girdi daha yeni: {', '.join(sorted(newer)[:3])}"
    changed = []
    for p in inputs:
        if not prev.get(p):
            continue
        if p.endswith(".sqlite"):
            diff = db_changed(prev[p], counts if counts is not None else db_counts(p),
                              ignore or set())
            if diff:
                changed.append(p if diff is True else
                               f"{p}[{','.join(diff[:3])}]")
        elif prev[p] != digest(p):
            changed.append(p)
    if changed:
        return f"girdi icerigi degisti: {', '.join(sorted(changed)[:3])}"
    if not prev:
        return "daha once bu araçla uretilmedi (temel kaydi yok)"
    return None


def script_flags(script):
    """Betik --cpu ya da --threads kabul ediyor mu -- KAYNAKTAN okunur.

    Calistirarak sormak yan etki uretiyordu; yukaridaki nota bak.
    """
    if not script or not os.path.exists(script):
        return None
    try:
        text = open(script, errors="replace").read()
    except OSError:
        return None
    if '"--cpu"' in text or "'--cpu'" in text:
        return "--cpu"
    if '"--threads"' in text or "'--threads'" in text:
        return "--threads"
    return None


def plan(state, propagate=True):
    """Her artefakt icin bayatlik sebebi (ya da None) -- TEK KAYNAK.

    NEDEN AYRI ISLEV (2026-11-05). Bu mantik `main()` icinde gomuluydu ve
    `validate_curation.py` "yayimlanan istatistik guncel mi?" sorusunu KENDI
    kuraliyla cevapliyordu: analiz dosyasinin mtime'ini `roar.sqlite`
    mtime'i ile karsilastiriyordu. O karsilastirma bu dosyanin 2026-10-14'te
    BILEREK kaldirdigi karsilastirmanin aynisi (yukaridaki "TABLO BAZINDA DB
    GIRDISI" notuna bak): DB'ye zincirin ICINDE yazan adimlar var, dolayisiyla
    DB mtime'i her kosuda en sona kayar ve kendisinden ONCE kosmus her analiz
    "bayat" gorunur. Olculdu: kapi 42 dosyanin 42'sini bayat gosterdi ve
    bazilari 0,3 saat "eski"ydi -- yani AYNI kosuda uretilmislerdi.

    Soru ayniysa kural da ayni olmali. Artik iki yerde iki kural yok: kapi
    da bu islevi cagiriyor.

    propagate=False ise bayatlik ILERI YAYILMAZ: her adim yalnizca kendi
    girdilerine karsi sinanir. Iki ayri soru vardir ve karistirilmamalidir:

      "yeniden kosmali mi?"      -> yayilim DAHIL (plan/run bunu sorar; ust
                                    adim kosacaksa alt adim da kosmalidir)
      "yayimlanan sayi yanlis mi?" -> yayilim HARIC (dogrulama kapisi bunu
                                    sorar; bir sayi, OKUDUGU tablo icerigi
                                    o gunden beri degistiyse yanlistir --
                                    ileride bir ust adimin kosacak olmasi
                                    yayimlanmis sayiyi yanlis yapmaz)

    Doner: [(name, outputs, cmd, inputs, heavy, why), ...], ARTEFACTS sirasinda.
    """
    # Her adimin yazdigi tablolar: statik (kaynak) + dinamik (olculen).
    writes = {}
    for name, outputs, cmd, inputs, heavy in ARTEFACTS:
        w = static_writes(cmd[1]) if len(cmd) > 1 else set()
        w |= DECLARED_WRITES.get(name, set())
        w |= set((state.get(name) or {}).get("db_writes") or [])
        writes[name] = w
    # Adim i'nin DB girdisi = kendisinden ONCEKILERIN yazdiklari. Yani i ve
    # sonrasinin yazdigi tablolar karsilastirmada yok sayilir.
    order = [a[0] for a in ARTEFACTS]
    ignore = {}
    # ignore = YALNIZCA bu adimdan SONRAKILERIN yazdigi tablolar. Adimin
    # KENDI yazdiklari ignore'a GIRMEZ.
    #
    # NEDEN. Onceki surumde kendisi de dahildi ve su hataya yol acti: 
    # build_operons.py `ro.sequence`i doldurdugu icin `ro` tablosu operons'un
    # kendi yazma kumesindeydi, dolayisiyla ignore'daydi; annotate'in `ro`da
    # yaptigi degisiklik operons'a GORUNMEDI ve adim atlandi. Sonuc: ro=19.114
    # iken operon=19.109 kaldi ve dogrulama kapisi tutarsizligi yakaladi
    # (2026-10-22).
    #
    # Kendi yazdiklarini ignore etmeye GEREK YOK, cunku karsilastirma adimin
    # KOSTUKTAN SONRAKI sayimina karsi yapiliyor; adimin kendi etkisi o kayda
    # zaten dahil.
    acc = set()
    for name in reversed(order):
        ignore[name] = set(acc)
        acc = acc | writes[name]

    counts_now = db_counts(DB) if os.path.exists(DB) else {}

    #  ... ve adimin adini hic anmadigi tablolar da karsilastirma disi
    #  (gerekcesi `unread_tables`'in basinda). Tablo listesi CANLI DB'den
    #  gelir, yani yeni bir tablo eklenince kendiliginden kapsanir.
    live_tables = {t for t in counts_now if not t.startswith("__")}
    for name, _o, cmd, _i, _h in ARTEFACTS:
        script = cmd[1] if len(cmd) > 1 else None
        ignore[name] |= unread_tables(script, live_tables)

    # Bayatligi SIRAYLA ILERI YAY: bir adim bayatsa ciktilari ve yazdigi
    # tablolar "kirli" sayilir; onlari okuyan sonraki adimlar da bayattir.
    # Eskiden bayatlik bir kez hesaplaniyordu ve bir ust adim yeniden
    # uretilse bile, basta guncel gorunen alt adimlar calismiyordu.
    rows = []
    dirty_files, dirty_tables = set(), set()
    for name, outputs, cmd, inputs, heavy in ARTEFACTS:
        why = staleness(name, outputs, inputs, state,
                        ignore=ignore[name], counts=counts_now)
        if not why and propagate:
            up_files = [p for p in inputs if p in dirty_files]
            up_tabs = (sorted(dirty_tables - ignore[name])
                       if DB in inputs else [])
            if up_files:
                why = f"ust adim yeniden uretilecek: {', '.join(up_files[:3])}"
            elif up_tabs:
                why = f"ust adim tablo yazacak: {', '.join(up_tabs[:3])}"
        if why:
            dirty_files.update(outputs)
            dirty_tables.update(writes[name])
        rows.append((name, outputs, cmd, inputs, heavy, why))
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true",
                   help="hicbir sey calistirma; bayat varsa 1 don")
    g.add_argument("--plan", action="store_true", help="ne calisacagini goster")
    g.add_argument("--run", action="store_true", help="bayat olanlari uret")
    g.add_argument("--baseline", action="store_true",
                   help="hicbir sey calistirmadan mevcut durumu guncel say "
                        "(ilk kurulum icin)")
    ap.add_argument("--include-heavy", action="store_true",
                    help="saatler suren adimlari da calistir")
    ap.add_argument("--only", action="append",
                    help="yalnizca bu artefakt(lar); tekrarlanabilir")
    ap.add_argument("--cpu", type=int, default=8)
    args = ap.parse_args()

    state = load_state()
    rows = plan(state)
    if args.only:
        rows = [r for r in rows if r[0] in args.only]

    if args.baseline:
        for name, outputs, cmd, inputs, heavy, _ in rows:
            prevw = (state.get(name) or {}).get("db_writes") or []
            state[name] = {"inputs": {p: digest(p) for p in inputs},
                           "db_writes": prevw,
                           "at": time.strftime("%Y-%m-%dT%H:%M:%S")}
        save_state(state)
        print(f"[temel] {len(rows)} artefakt guncel sayildi -> {STATE}")
        return 0

    stale = [r for r in rows if r[5]]
    fresh = len(rows) - len(stale)
    print(f"{len(rows)} artefakt: {fresh} guncel, {len(stale)} BAYAT")
    if stale:
        print()
        for name, outputs, cmd, inputs, heavy, why in stale:
            tag = " [AGIR]" if heavy else ""
            print(f"  BAYAT  {name}{tag}")
            print(f"         sebep: {why}")
    if args.check:
        print()
        if stale:
            print("Yayinlamadan once bunlari uret: python3 rebuild.py --run")
            return 1
        print("Her sey guncel.")
        return 0

    runnable = [r for r in stale if args.include_heavy or not r[4]]
    skipped = [r for r in stale if r not in runnable]
    if args.plan:
        print()
        print("Sira:")
        for name, outputs, cmd, inputs, heavy, why in runnable:
            print(f"  {name:24s} {' '.join(cmd)}")
        for name, *_ in skipped:
            print(f"  {name:24s} [AGIR -- atlanir; --include-heavy ile calisir]")
        return 0

    db_before = db_fingerprint(DB) if os.path.exists(DB) else None
    failed = []
    for i, (name, outputs, cmd, inputs, heavy, why) in enumerate(runnable, 1):
        print(f"\n==> [{i}/{len(runnable)}] {name}  ({why})")
        full = list(cmd)
        # --cpu kabul eden betiklere gecir.
        #
        # ESKIDEN betigi `--help` ile CALISTIRIP cikisinda "--cpu" arardi.
        # Bu, argparse'i OLMAYAN bir betikte `--help`in sessizce yok sayilip
        # betigin GERCEKTEN kosmasi demekti: 60 saniyelik bir "prob" ic
        # ag sorgulari atar, dosya yazar, ve cikisi atilir. Yan etkisi olan
        # bir yoklama yoklama degildir. (2026-11-05; `anchor_reactions.py` ve
        # `reaction_balance.py` grafige eklenirken bulundu -- ikisi de
        # argparse'siz ve ikisi de dis servise istek atiyor.)
        #
        # Artik KAYNAKTAN okunur: calistirmak gerekmez, yan etkisi yoktur ve
        # hizlidir.
        flag = script_flags(full[1] if len(full) > 1 else "")
        if flag:
            full += [flag, str(args.cpu)]
        t0 = time.time()
        before = db_counts(DB) if os.path.exists(DB) else {}
        res = subprocess.run(full)
        if res.returncode != 0:
            failed.append(name)
            print(f"    BASARISIZ ({res.returncode}) -- zincir burada duruyor, "
                  f"cunku sonraki adimlar bunun ciktisini okuyor")
            break
        after = db_counts(DB) if os.path.exists(DB) else {}
        observed = {_table_of(k) for k in set(before) | set(after)
                    if before.get(k) != after.get(k)}
        prevw = set((state.get(name) or {}).get("db_writes") or [])
        # DB girdisinin kaydi: adim CALISTIKTAN SONRAKI sayim. Karsilastirmada
        # kendisinin ve sonrakilerin tablolari zaten yok sayiliyor.
        state[name] = {"inputs": {p: (after if p == DB else digest(p))
                                  for p in inputs},
                       "db_writes": sorted(prevw | observed),
                       "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                       "seconds": round(time.time() - t0, 1)}
        save_state(state)
        print(f"    tamam, {time.time() - t0:.1f} s")

    # Veritabani gercekten degistiyse bir SURUM kaydet. Her --run'da degil:
    # yalnizca mantiksal parmak izi oynadiginda, yani sayilar degistiginde.
    # Boylece "bu tip neden 94'ten 310'a cikti" sorusunun cevabi sonradan
    # aranabilir hale gelir.
    if db_before and db_fingerprint(DB) != db_before:
        ran = ", ".join(n for n, *_ in runnable if n not in failed) or "-"
        subprocess.run(["python3", "db_revision.py", "--db", DB, "--note",
                        f"rebuild.py --run: {ran}"])

    print()
    if failed:
        print(f"DURDU: {', '.join(failed)} basarisiz. Duzeltip tekrar calistir.")
        return 1
    if skipped:
        print(f"Atlanan agir adimlar: {', '.join(n for n, *_ in skipped)}")
        print("Bunlar gerekiyorsa: python3 rebuild.py --run --include-heavy")
    print("Bitti. Yayin oncesi: python3 rebuild.py --check && python3 check_site.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
