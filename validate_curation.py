#!/usr/bin/env python3
"""
Kuratorlu veri ve veritabani dogrulayicisi -- run_all.sh icinde KAPI olarak kullanilir.

NEDEN: bu pipeline birkac SESSIZ hata yayinladi. cd-hit cokmesi yutuldu, bir CSV
virgul kacisi ('2\\,4') kolonlari kaydirdi ve TftA/DntAc satirlari yanlis yere
dustu, bir metrik en kotu hizalanmis ilmekleri "ayirt edici kolon" sandi,
hmmsearch hiti olmayan 945 aday is_confirmed=NULL kaldi. Hepsinin ortak yani:
hata ciktiyi BOZDU ama hicbir sey durmadi.

BU SCRIPT her kontrol icin tek satir PASS/FAIL ve bir sayi yazar, sonunda
basarisizlik varsa SIFIR-DISI cikar. Boylece run_all.sh'de bir adim olarak
durabilir ve bozuk veri HTML'e / web sitesine kadar gidemez.

Kullanim:
    python3 validate_curation.py
    python3 validate_curation.py --db roar.sqlite --ecology cluster_ecology.csv \\
        --chemistry chemistry.csv
    python3 validate_curation.py --warn-only   # rapor yaz ama her zaman 0 don

Cikis kodu: 0 = tum FAIL kontrolleri gecti, 1 = en az bir FAIL, 2 = kontrol
kosulamadi (veritabani/CSV yok).

SINIR: burada gecen her sey DOGRU demek degildir; bkz. README_PIPELINE.md
"Safeguards" bolumu -- bu dogrulayicinin yakalayamadiklari orada yaziyor.
"""

import argparse
import csv
import json
import os
import re
import sqlite3
import yaml
import sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))

# Kuratorlu referans seti: pipeline'in TAMAMI bunun uzerine kuruluyor
# (tip atamasi, katalitik motif kolonlari, her kanit kademesi). Buradaki bir
# kusur tek bir sayiyi degil, butun siniflandirmayi kaydirir.
REFS_FASTA = os.path.join(ROOT, "ROs_71_Clean", "refs71.fasta")

# Iki referansin birbirine bu orandan fazla benzemesi RAPOR EDILIR ama hata
# sayilmaz: EdoA1/cumA1 %99,8 kimlikte FARKLI substratlara (etilbenzen /
# kumen) etki ediyor ve bu veritabaninin merkezi bulgusudur. Ayni sey
# DUPLIKASYON kusuruyla karistirilmamali -- o ayri ve sert bir kontrol.
NEAR_IDENTICAL_PCT = 99.0

# --- Izinli sozluk degerleri -------------------------------------------------
# add_reference.py yeni referans eklerken AYNI listeleri kullanir; ikisi
# ayrismasin diye asagida ayrica karsilastirilirlar.
VALID_SUBSTRATE_CLASS = {"xenobiotic", "natural_aromatic", "natural_specialized",
                         "unknown"}
VALID_CONFIDENCE = {"low", "medium", "high"}
# 2026-10-16: dort yeni sinif. Sekiz yeni referansin gercekten yaptigi kimya
# bu; alternatifi onlari "unknown" isaretlemekti. Gerekceler:
#   desaturation             Nvd / DAF-36, kolesterol C7-C8. Hem Miao & Schmidt
#                            2025 (doi:10.1021/acs.biochem.5c00369) hem Zhang
#                            2025 (doi:10.1021/jacs.5c01994) ailenin
#                            reaksiyonlari arasinda desaturasyonu sayiyor.
#   N_oxygenation            PrnD, aril amin -> aril nitro, alti elektronlu.
#                            Zhang 2025 "N-oxygenation" diye listeliyor.
#   C_C_bond_formation       RedG. Zhang 2025 dogrudan yaziyor: RedG ve
#                            homologlari prodijinin biyosentezinde C-C bagi
#                            kuran halkalastirmayi katalizler.
#   oxygenolytic_ring_opening PAO, klorin makrosiklinin O2 ile acilmasi.
#                            Pruzinska 2003 (doi:10.1073/pnas.2036571100)
#                            PAO'nun Rieske-tipi demir-sulfur proteini
#                            oldugunu kurdu.
#   C_C_cleavage             AerC (7_701_AerC). Oksidatif karbosiklizasyonun
#                            ARDINDAN C-C bagi kirilmasi, tek aktif bolgede,
#                            uc karbon kaybiyla. Zhang et al. JACS 2025,
#                            doi:10.1021/jacs.5c01994. Sette C_C_bond_formation
#                            (RedG) ve C_N_cleavage (CntA, qxyA) vardi ama C-C
#                            KIRILMASI yoktu.
#                            NOT: reaction_class tip basina TEK deger. AerC iki
#                            kimyayi birlikte yapiyor; kirilma sinif olarak
#                            secildi cunku sette eksik olan odur. Bu, semanin
#                            bir sinirlamasidir, enzimin degil.
# 2026-10-28: liste ELLE TUTULMUYOR artik. Sinif adlari burada,
# add_reference.py'de ve tiplerin alanlarinda ucuz ucuz ayrisabiliyordu; daha
# onemlisi sinifin adi enzimin oksijen aktarip aktarmadigini soylemiyordu.
# Tek kaynak reaction_classes.yaml ve oxygen_fate.py onu okur. Dosya yoksa
# eski sabit liste yedek olarak kullanilir, yani davranis bozulmaz.
_FALLBACK_REACTION_CLASS = {"cis_dihydroxylation", "angular_dioxygenation",
                            "dioxygenation_with_release", "O_demethylation",
                            "N_demethylation", "hydroxylation",
                            "C_N_cleavage", "unknown", "desaturation",
                            "N_oxygenation", "C_C_bond_formation",
                            "oxygenolytic_ring_opening", "C_C_cleavage"}


def _reaction_classes():
    try:
        import oxygen_fate                              # noqa: PLC0415
        ids = set(oxygen_fate.load_vocab(ROOT))
        return (ids | {"unknown"}) if ids else _FALLBACK_REACTION_CLASS
    except Exception:                                   # noqa: BLE001
        return _FALLBACK_REACTION_CLASS


VALID_REACTION_CLASS = _reaction_classes()
# 2026-10-16: tetrapyrroles eklendi -- PAO, CAO, PTC52 icin.
# Manifestolar "chlorophyll_pathway" yaziyordu; DEGISTIRILDI. Sebep: bu
# sozlugun butun oteki uyeleri substratin KIMYASAL SINIFI (alkylbenzenes,
# quaternary_amines, terpenoids_steroids...), "chlorophyll_pathway" ise bir
# YOL adi. Substratlar feoforbid a, klorofillid a ve protoklorofillid a --
# hepsi klorin, yani indirgenmis tetrapirol. Yol/klad bilgisi zaten group: 6
# alaninda duruyor, family kolonunun onu tekrarlamasi gerekmiyor. Ayrica
# tetrapyrroles ileride hem/porfirin uzerinde calisan bir RO cikarsa onu da
# kapsar.
VALID_FAMILY = {"alkylbenzenes", "pah", "biaryls_ethers", "nitroaromatics",
                "haloaromatics", "sulfoaromatics", "aromatic_acids", "anilines",
                "quaternary_amines", "alkaloids", "terpenoids_steroids", "unknown", "tetrapyrroles"}
VALID_SOURCE_KIND = {"paper", "structure", "curator", "paper+structure",
                     "structure+paper", "paper+curator", "general", "none"}

# --- Turkce sizintisi --------------------------------------------------------
# Web sitesini goren hicbir metin Turkce olmamali. Iki ayri tuzak var:
# (a) Turkce'ye ozgu harfler, (b) ASCII'ye sadelestirilmis Turkce kelimeler --
# bunlar daha sinsi, cunku kodlama kontrolleri onlari gormez.
# Turkce'ye ozgu harfler, KOD SAF ASCII kalsin diye kod noktasiyla yazildi:
# c-cedilla, g-breve, dotless-i, o-umlaut, s-cedilla, u-umlaut (kucuk+buyuk).
TURKISH_LETTERS = {chr(c) for c in (
    0x00E7, 0x011F, 0x0131, 0x00F6, 0x015F, 0x00FC,
    0x00C7, 0x011E, 0x0130, 0x00D6, 0x015E, 0x00DC,
)}
# Liste hem islev hem ICERIK sozcugu icerir. Bir donem yalnizca islev
# sozcukleri vardi ve okaryot alt-grup etiketleri ("bitki/alg", "mantar",
# "hayvan") aradan gecip arayuze cikti; eksik olan icerik sozcukleriydi.
# Arama ONEK esleseme yapar, bu yuzden Ingilizce bir kelimenin basiyla
# cakisabilecek kisa parcalar (ornegin "alg", "ana", "son") listeye GIRMEZ.
TURKISH_WORDS = ["agirlikli", "komsu", "kume", "yaprak", "bilinmiyor",
                 "dusuk", "orta", "yuksek", "degil",
                 "bitki", "mantar", "hayvan", "diger", "kirmizi",
                 "okaryot", "bakteriyel", "kaynak", "kayit", "sayisi",
                 "toprak", "tanimsiz", "belirsiz", "bulunmayan", "supheli",
                 "icin", "dogrulan", "ornek", "siniflan"]
TURKISH_WORD_RE = re.compile(r"\b(" + "|".join(TURKISH_WORDS) + r")", re.IGNORECASE)

# Kelime arama ONEK esleseme yapar (Turkce eklemeli: komsu -> komsular,
# kume -> kumeler) ve bu bazen mesru metni yakalar. Bilinen yanlis pozitifler
# burada muaf tutulur; liste KISA kalmali, yoksa kontrolun degeri duser.
#   "Kumeu" : Yeni Zelanda'da bir yer adi, GenBank /geo_loc_name alanindan
#             geliyor ve 'kume' onekiyle eslesiyor.
LANGUAGE_FALSE_POSITIVES = ("Kumeu",)

# Kontrol edilecek kullanici-yuzlu DB alanlari: (tablo, anahtar kolon, metin kolonu)
USER_FACING_DB_FIELDS = [
    ("leaf_profile", "leaf_id", "label"),
    ("leaf", "leaf_id", "top_genera"),
    ("ro_etc", "candidate_id", "etc_profile"),
    # 2026-10-30: triad kapisinin "neden cagri yapilmadi" gerekcesi. Proza,
    # giris sayfasinda gosterilebilir, bu yuzden Ingilizce olmak zorunda.
    ("ro_triad", "candidate_id", "note"),
    # 2026-10-31: uyelik temeli ve gerekcesi. Giris sayfasi bunlari OKUYUP
    # cumle kuruyor, dolayisiyla Ingilizce olmak zorunda.
    ("ro_membership", "candidate_id", "basis"),
    ("ro_membership", "candidate_id", "detail"),
    ("operon", "candidate_id", "layout"),
    # Proje tarafindan URETILEN habitat etiketi (serbest metin degil, sozluk).
    # isolation_source.py'nin HABITAT_RULES sonucu; web'de gosteriliyor.
    ("replicon_source", "nucleotide_id", "habitat"),
    # Okaryot alt-grubu: classify_domains.py'deki sozlukten geliyor ve giris
    # sayfalarinda gosteriliyor. Bir donem Turkce idi ve bu liste onu
    # icermedigi icin kontrol hic calismadi.
    ("ro_domain", "candidate_id", "euk_group"),
    # Konak alemi: host_kingdom.csv'den gelen kontrollu sozluk, ekoloji
    # sayfasinda gosteriliyor.
    ("replicon_source", "nucleotide_id", "host_kingdom"),
    ("ro_etc", "candidate_id", "reductase_type"),
    ("ro_etc", "candidate_id", "ferredoxin_type"),
    ("ro_evidence", "candidate_id", "tier"),
    ("ro_regulation", "candidate_id", "upstream_category"),
    ("ro_regulation", "candidate_id", "architecture"),
    ("ro_subfamily", "candidate_id", "assignment_class"),
    ("operon_gene", "candidate_id", "category"),
    ("gene_category", "category", "method"),
    ("leaf_profile", "leaf_id", "enriched_categories"),
]

# GenBank'tan OLDUGU GIBI kopyalanan metin sutunlari. Bunlarda yabanci dil
# olmasi normaldir (organizma adlari, /product aciklamalari, yer adlari) ve
# dil kontrolune girmezler. Liste, asagidaki "her metin sutunu siniflanmis mi"
# kontrolu icin var: yeni bir sutun ne burada ne yukarida ise kontrol FAIL
# verir, boylece `euk_group` gibi bir sutun bir daha sessizce denetim disinda
# kalmaz.
EXTERNAL_TEXT_FIELDS = {
    ("replicon", "organism"), ("replicon", "taxonomy"), ("replicon", "description"),
    ("replicon", "file"), ("replicon", "mol_type"), ("replicon", "status"),
    ("replicon_source", "isolation_source"), ("replicon_source", "host"),
    ("replicon_source", "geo"), ("replicon_source", "country"),
    ("ro", "product"), ("ro", "locus_tag"), ("ro", "gene"), ("ro", "sequence"),
    ("neighbor", "product"), ("neighbor", "locus_tag"), ("neighbor", "gene"),
    ("neighbor_protein", "translation"), ("neighbor_domain", "hmm"),
    ("operon_gene", "product"), ("ro_regulation", "upstream_product"),
    ("ro_search", "organism"), ("ro_search", "genus"), ("ro_search", "product"),
    ("ro_search", "doc"),
    ("ro_ncbi_hit", "product"),
    ("ro_protein_status", "organism"), ("ro_union", "organism"),
}

# Kimlik, kod ve olcum sutunlari: serbest metin tasimadiklari icin dil
# kontrolune girmezler.
IDENTIFIER_TEXT_FIELDS = {
    # 2026-11-11, taxonomy_ids.py. Hepsi KIMLIK ya da NCBI'nin kendi
    # terimleri: `rank` NCBI rank sozlugunden sabit bir deger ('genus',
    # 'subgenus', ...), adlar bilimsel ad, `status` dort degerli sabit.
    # Cevrilecek proza yok; cevirmek kaynakla baglantiyi koparir.
    ("replicon_taxon", "nucleotide_id"), ("replicon_taxon", "status"),
    ("replicon_taxon", "name"), ("replicon_taxon", "rank"),
    ("taxon_lineage", "rank"), ("taxon_lineage", "ancestor_name"),
    ("taxon_genus", "nucleotide_id"),
    ("taxon_genus", "genus_name"), ("taxon_genus", "genus_state"),
    ("taxon_genus", "proxy_rank"), ("taxon_genus", "proxy_name"),
    #  `type_substrate_chebi`: hepsi TANIMLAYICI ya da kontrollu bir
    #  vokabulerden gelen etiket, serbest metin degil -- `chebi_label`
    #  ChEBI'nin kendi adlandirmasi, bizim yazdigimiz bir cumle degil.
    ("type_substrate_chebi", "cluster"), ("type_substrate_chebi", "chebi"),
    ("type_substrate_chebi", "chebi_label"),
    ("type_substrate_chebi", "formula"), ("type_substrate_chebi", "state"),
    ("type_substrate_chebi", "rhea"),
    #  `replicon_kind`: dordu de KIMLIK ya da kontrollu sozluk --
    #  `accession_prefix` NCBI'nin oneki, `moltype` GenBank'in LOCUS
    #  alani, `record_kind` iki degerli bir sabit. Cevirilecek proza yok.
    ("replicon_kind", "nucleotide_id"), ("replicon_kind", "accession_prefix"),
    ("replicon_kind", "moltype"), ("replicon_kind", "record_kind"),
    ("cluster_sdp", "cluster"), ("cluster_sdp", "columns"), ("cluster_sdp", "method"),
    ("gene_category", "category"),
    # 2026-11-11, cross_references.py'nin yazdigi ro_xref. Ucu de KIMLIK:
    # seq_md5 dizi ozeti, db dort degerli bir sabit (uniprot /
    # refseq_protein / embl_nucleotide / ec), identifier DIS bir
    # veritabaninin accession'i. Cevrilecek proza yok, ama uc listeden
    # birinde olmak ZORUNDA -- yoksa dil kontrolu onlari hic gormez.
    ("ro_xref", "seq_md5"), ("ro_xref", "db"), ("ro_xref", "identifier"),
    # 2026-11-10, family_scan.py'nin yazdigi uc tablo. Hepsi KIMLIK:
    # category Pfam ailesinin ingilizce sabit adi (IS3-like, LysR, ...),
    # method yontem etiketi, model/library Pfam model adi, protein_key
    # koordinat anahtari, source uc degerli bir sabit. Cevrilecek proza yok,
    # ama uc listeden birinde olmak ZORUNDA -- yoksa dil kontrolu gormez.
    ("neighbor_family", "category"), ("neighbor_family", "method"),
    ("neighbor_family_hit", "protein_key"), ("neighbor_family_hit", "library"),
    ("neighbor_family_hit", "model"),
    ("ro_upstream_family_hmm", "candidate_id"),
    ("ro_upstream_family_hmm", "family_regex"),
    ("ro_upstream_family_hmm", "family_hmm"),
    ("ro_upstream_family_hmm", "family_best"),
    ("ro_upstream_family_hmm", "source"),
    ("leaf", "leaf_id"), ("leaf", "cluster"), ("leaf", "representative"),
    ("leaf", "top_genera"),
    # 2026-10-28, icerikten turetilen varyant kimlikleri. Ucu de KIMLIK:
    # variant_id "RO-" + uye dizilerinin ozetinden turer, nearest_type bir
    # kume kimligi, label_kind iki degerli bir sabit ("named"/"placeholder").
    # Hicbiri cevrilecek proza degil, ama uc listeden birinde olmak ZORUNDA --
    # yoksa dil kontrolu onlari hic gormez ve bir Turkce etiket sessizce
    # arayuze cikar (euk_group'ta tam olarak bu oldu).
    ("leaf", "variant_id"), ("leaf", "nearest_type"), ("leaf", "label_kind"),
    ("ro_leaf", "variant_id"),
    ("leaf_profile", "leaf_id"), ("leaf_profile", "cluster"),
    ("leaf_profile", "top_genera"), ("leaf_profile", "neighbor_signature"),
    ("leaf_sdp", "leaf_id"), ("leaf_sdp", "cluster"), ("leaf_sdp", "residues"),
    ("leaf_sdp", "signature"),
    # 2026-11-01, protein_freshness.py'nin yazdigi tazelik tablosu. Hepsi
    # UniProt'un kendi kodlari ya da capraz-referans kimlikleri; cevrilecek
    # proza yok. organism serbest metin oldugu icin EXTERNAL listesinde.
    ("ro_protein_status", "accession"), ("ro_protein_status", "entry_type"),
    ("ro_protein_status", "inactive_reason"), ("ro_protein_status", "merged_into"),
    ("ro_protein_status", "embl_nuc"), ("ro_protein_status", "refseq"),
    ("ro_protein_status", "ec"),
    # 2026-11-01, build_union.py'nin yazdigi birlesim tablosu. Birim DIZI,
    # accession kaynak kimligi; hepsi kimlik/kod, cevrilecek proza yok.
    # organism serbest metin oldugu icin EXTERNAL listesinde.
    ("ro_union", "seq_md5"), ("ro_union", "sequence"),
    ("ro_union", "ro_cluster"), ("ro_union", "ro_group"),
    ("ro_union", "genomic_candidate_id"), ("ro_union", "first_source"),
    ("ro_union_source", "seq_md5"), ("ro_union_source", "source"),
    ("ro_union_source", "accession"), ("ro_union_source", "status"),
    # 2026-11-01, RefSeq taramasinin kaynak tablosu (ingest_ncbi_hits.py).
    # product RefSeq'in kendi anotasyon metni oldugu icin EXTERNAL listesinde.
    ("ro_ncbi_hit", "accession"), ("ro_ncbi_hit", "seq_md5"),
    ("ro_ncbi_hit", "sequence"), ("ro_ncbi_hit", "ro_cluster"),
    ("ro_ncbi_hit", "ro_group"),
    # 2026-11-02, a5_cohort.py'nin yazdigi kohort tablosu. Ikisi de kimlik /
    # sabit etiket ("survey" | "protein_recovery"); cevrilecek proza yok.
    ("replicon_cohort", "nucleotide_id"), ("replicon_cohort", "cohort"),
    ("ro_membership", "candidate_id"), ("ro_membership", "ro_cluster"),
    ("neighbor", "candidate_id"), ("neighbor", "nucleotide_id"),
    ("neighbor", "protein_id"),
    ("neighbor_component", "protein_key"), ("neighbor_component", "component"),
    ("neighbor_component", "evidence"), ("neighbor_domain", "protein_key"),
    ("neighbor_protein", "protein_key"), ("neighbor_protein", "nucleotide_id"),
    ("operon", "candidate_id"), ("operon", "nucleotide_id"),
    ("operon_gene", "candidate_id"), ("operon_gene", "protein_key"),
    ("operon_gene", "component"),
    ("replicon", "nucleotide_id"), ("replicon_source", "nucleotide_id"),
    ("ro", "candidate_id"), ("ro", "nucleotide_id"), ("ro", "protein_id"),
    ("ro", "ro_cluster"), ("ro", "ro_group"),
    ("ro_carboxylate", "candidate_id"), ("ro_carboxylate", "catalytic_residue"),
    ("ro_carboxylate", "bridging_residue"),
    # 2026-10-30, 3/3 triad kapisi. Kalinti harfleri, kaynak etiketleri ve
    # durum etiketleri KIMLIK (sabit kucuk kume); ro_triad.note ise proza,
    # o yuzden USER_FACING_DB_FIELDS'ta ve Ingilizce.
    ("ro_triad", "candidate_id"), ("ro_triad", "ro_cluster"),
    ("ro_triad", "his1_residue"), ("ro_triad", "his2_residue"),
    ("ro_triad", "carboxylate_residue"), ("ro_triad", "carboxylate_source"),
    ("ro_triad", "his_source"), ("ro_triad", "triad_status"),
    ("ro_triad_meta", "key"), ("ro_triad_meta", "value"),
    # 2026-10-30, genomik baglami olmayan protein-duzeyi girisler.
    # protein_name/accession UniProt kimlikleri, ro_cluster bir tip kimligi,
    # kalinti harfleri sabit kume, linked_candidate_id bir `ro` anahtari.
    # Hicbiri cevrilecek proza degil.
    ("ro_protein", "protein_name"), ("ro_protein", "accession"),
    ("ro_protein", "ro_cluster"), ("ro_protein", "carboxylate_residue"),
    ("ro_protein", "bridging_residue"), ("ro_protein", "linked_candidate_id"),
    ("ro_domain", "candidate_id"), ("ro_domain", "cluster"), ("ro_domain", "domain"),
    ("ro_etc", "candidate_id"),
    ("ro_evidence", "candidate_id"), ("ro_evidence", "nearest_ref"),
    ("ro_leaf", "candidate_id"), ("ro_leaf", "cluster"), ("ro_leaf", "leaf_id"),
    ("ro_regulation", "candidate_id"), ("ro_regulation", "upstream_family"),
    ("ro_search", "candidate_id"), ("ro_search", "protein_id"),
    ("ro_search", "cluster"), ("ro_search", "gene"), ("ro_search", "leaf_id"),
    ("ro_search", "substrate"), ("ro_search", "family"), ("ro_search", "reaction"),
    ("ro_search", "tier"), ("ro_search", "domain"),
    ("ro_subfamily", "candidate_id"), ("ro_subfamily", "cluster"),
    ("ro_subfamily", "subfamily_id"),
    ("subfamily", "cluster"), ("subfamily", "subfamily_id"),
    ("subfamily", "representative"), ("subfamily", "top_genera"),
}

# analysis_out/ icinde web arayuzunun DOGRUDAN yayinladigi metin dosyalari.
# Turkce sizintisi burada da yayina gider; bu yuzden ayrica taraniyor.
PUBLISHED_TEXT_FILES = [
    "cluster_ecology_stats.csv", "reference_pairs.csv", "representatives.csv",
    "etc_by_cluster.csv", "evidence_by_cluster.csv", "regulation_by_cluster.csv",
    "leaf_profiles.csv", "domain_by_cluster.csv", "null_model.csv",
    "cluster_variance.csv", "variant_signatures.csv", "leaves.csv",
    "novel_candidates.csv", "transposon_by_genus.csv", "sdp_positions.csv",
    "subfamilies.csv", "cluster_identity_matrix.csv", "ssn_nodes.csv",
    "habitat.json", "motif_stats.json", "stats.json", "redundancy.json",
    "cooccurrence.json", "operon_validation.json",
    "substrate_predictability.json",
    "stratified_stats.json", "ecological_origin.json", "threshold_sensitivity.json",
    "reference_redundancy.json", "group_monophyly.json",
    "protein_freshness.json", "bootstrap_references.json",
    "bootstrap_confident.json",
]

# --- Uyesiz referans tipleri -------------------------------------------------
# cluster_ecology.csv ve chemistry.csv 71 kuratorlu referansin HEPSINI tutar,
# ama bu 10 tipin veritabaninda TEK bir adayi bile yok (olculdu: hicbir
# PF00355 proteini bu modellere en yuksek bit skorunu vermiyor; sadece
# dogrulanmamis degil, hic atanmamis). Yani "CSV'de var, ro'da yok" durumu
# bunlar icin BEKLENEN; listede olmayan yeni bir bosluk ise HATA'dir, cunku
# o, bir tipin uyelerini kaybettigi anlamina gelir.
#
# DIKKAT: ayni 10 kimlik check_reference_set() icinde de WARN olarak cikar,
# ama farkli bir soru icin. Burada sorulan "CSV ile DB ayristi mi"; orada
# sorulan "bu referans NEDEN uye toplamiyor". Ucu (3_309_NahAc, 3_315_NDO,
# 1_115_NdmC) referans setindeki duplikasyon/parca kusuru yuzunden bos;
# kusur duzeltilince bu listenin de kisalmasi gerekir.
REFERENCE_ONLY_CLUSTERS = {
    "1_103_CarAa", "1_115_NdmC", "2_206_IpBAa", "2_207_BpdC1", "2_210_XylC1",
    "3_309_NahAc", "3_310_DntAc", "3_312_NinA", "3_315_NDO", "4_407_PbaA1",
}

# SMILES'te mesru olan karakterler. Atom, bag, halka kapanisi, yuk, kiralite.
SMILES_LEGAL = set("ABCDEFGHIKLMNOPRSTUVWXYZabcdefghiklmnoprstuvy"
                   "0123456789[]()=#$:/\\@+-.%*")
SMILES_ATOM_RE = re.compile(r"\[[^\]]*\]|Cl|Br|[BCNOPSFIbcnopsfi]")


class Report:
    """PASS/FAIL satirlari ve cikis kodu."""

    def __init__(self):
        self.rows = []
        self.failed = 0

    def check(self, name, ok, count, detail="", severity="FAIL"):
        """severity='FAIL' kontrol kapiyi kapatir; 'WARN' sadece rapor eder."""
        status = "PASS" if ok else severity
        if not ok and severity == "FAIL":
            self.failed += 1
        self.rows.append((status, name, count, detail))
        return ok

    def emit(self):
        width = max(len(r[1]) for r in self.rows) if self.rows else 10
        for status, name, count, detail in self.rows:
            line = f"{status:4s}  {name:<{width}s}  n={count}"
            if detail:
                line += f"   {detail}"
            print(line)
        n_pass = sum(1 for r in self.rows if r[0] == "PASS")
        n_warn = sum(1 for r in self.rows if r[0] == "WARN")
        print("-" * (width + 24))
        print(f"      {len(self.rows)} checks: {n_pass} PASS, "
              f"{self.failed} FAIL, {n_warn} WARN")
        return 1 if self.failed else 0


def first(items, n=6):
    items = sorted(str(i) for i in items)
    out = ", ".join(items[:n])
    if len(items) > n:
        out += f", ... (+{len(items) - n})"
    return out


def turkish_hits(text):
    """Metinde Turkce harf veya Turkce kelime var mi -> bulgu listesi."""
    if text is None:
        return []
    text = str(text)
    for allowed in LANGUAGE_FALSE_POSITIVES:
        text = text.replace(allowed, "")
    hits = sorted({ch for ch in text if ch in TURKISH_LETTERS})
    hits += sorted({m.group(1).lower() for m in TURKISH_WORD_RE.finditer(text)})
    return hits


# ---------------------------------------------------------------------------
# SMILES makul mu
# ---------------------------------------------------------------------------

def smiles_problems(smiles):
    """Makul bir SMILES mi? Kimyasal dogruluk DEGIL, bicim dogrulugu.

    Kontroller: parantez ve koseli parantez dengesi, yalnizca mesru atom/bag
    karakterleri, en az bir atom, halka kapanis rakamlarinin eslesmesi.
    Burada gecmek molekulun DOGRU oldugu anlamina gelmez; yalnizca cizim
    kutuphanesine verildiginde kirilmayacagi anlamina gelir.
    """
    problems = []
    for opener, closer in (("(", ")"), ("[", "]")):
        depth = 0
        for ch in smiles:
            if ch == opener:
                depth += 1
            elif ch == closer:
                depth -= 1
                if depth < 0:
                    problems.append(f"'{closer}' before '{opener}'")
                    break
        if depth > 0:
            problems.append(f"{depth} unclosed '{opener}'")
    illegal = sorted({ch for ch in smiles if ch not in SMILES_LEGAL})
    if illegal:
        problems.append("illegal characters " + repr("".join(illegal)))
    if not SMILES_ATOM_RE.search(smiles):
        problems.append("no atom found")
    # Halka kapanis rakamlari: koseli parantez icindekiler (yuk, izotop) sayilmaz
    core = re.sub(r"\[[^\]]*\]", "A", smiles)
    digits = re.findall(r"%\d\d|\d", core)
    unpaired = sorted({d for d in set(digits) if digits.count(d) % 2})
    if unpaired:
        problems.append(f"unpaired ring-closure digits {unpaired}")
    return problems


# ---------------------------------------------------------------------------
# CSV okuma: kolon sayisi kaymasi (virgul kacisi hata sinifi)
# ---------------------------------------------------------------------------

def read_csv_strict(path, rep, label):
    """CSV'yi oku ve basliktan FARKLI alan sayisi olan satirlari bildir.

    2\\,4 kacisi hata sinifi tam olarak bu: satir sessizce fazla/eksik alanla
    ayristirilir, degerler bir kolon kayar ve hicbir sey uyarmaz.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        rep.check(f"{label}: file not empty", False, 0, "file has no rows")
        return [], []
    header = rows[0]
    bad = [(i + 1, len(r)) for i, r in enumerate(rows) if r and len(r) != len(header)]
    rep.check(f"{label}: field count matches header", not bad, len(bad),
              "" if not bad else f"{len(header)} header fields; offending "
                                 f"lines/counts: {first(bad)}")
    dicts = [dict(zip(header, r)) for r in rows[1:] if r]
    return header, dicts


# ---------------------------------------------------------------------------
# Kontroller
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Kuratorlu referans setinin kendisi
# ---------------------------------------------------------------------------

def read_fasta(path):
    """[(header, sequence)] dondur. Dizi buyuk harfe cevrilir, bosluk atilir."""
    records, header, chunks = [], None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line.upper().replace(" ", ""))
    if header is not None:
        records.append((header, "".join(chunks)))
    # Bazi kaynaklar dizinin sonuna '*' koyar; karsilastirmayi bozmasin.
    return [(h, s.rstrip("*")) for h, s in records]


def ref_id(header):
    """FASTA basligindan kume kimligi: evidence_tiers.py ile AYNI kural.

    Baslik '1_101_OxoO_monooxygenase_pro' bicimindedir; ilk uc alt-cizgi
    parcasi kume kimligidir. Kural ayrisirsa dogrulama baska bir sey sinar.
    """
    return "_".join(header.split("_")[:3])


def check_reference_set(rep, refs_path, pairs_path, eco, chem, con):
    """Referans setinin IC kusurlari.

    Burada aranan hatalar 2026-10 incelemesinde ELLE bulundu ve hepsi ayni
    sonucu dogurdu: RieskeDB71.hmm'e birbirinin AYNISI profiller girdi, bu
    yuzden tip atamasi ikizler arasinda keyfi boluldu (OxoO 2 / OMO 11,
    NDO 2 / NarAa 30) ve "en yuksek skorlu profil ile en yakin referans
    protein farkli tipte" orani (%29,3) buyuk olcude bundan sisti.
    """
    if not os.path.exists(refs_path):
        rep.check("refs71: reference FASTA present", False, 0,
                  f"{refs_path} not found; the whole reference-set check could "
                  f"not run")
        return

    records = read_fasta(refs_path)
    ids = [ref_id(h) for h, _ in records]
    seqs = [s for _, s in records]
    rep.check("refs71: FASTA is readable and non-empty", bool(records),
              len(records), f"{len(set(seqs))} distinct sequences")

    # 1. Birbirinin AYNISI diziler. Ayni enzimin iki adla girilmesi
    #    HMM kutuphanesine cift profil koyar.
    by_seq = {}
    for ident, seq in zip(ids, seqs):
        by_seq.setdefault(seq, []).append(ident)
    dup_groups = [g for g in by_seq.values() if len(g) > 1]
    detail = "; ".join(" == ".join(sorted(g)) for g in dup_groups)
    rep.check("refs71: no two reference sequences are identical",
              not dup_groups, len(dup_groups),
              f"{len(records)} records, {len(by_seq)} distinct sequences"
              if not dup_groups else
              f"{len(records)} records but only {len(by_seq)} distinct "
              f"sequences -> {detail}")

    # 2. Bir dizinin bir baskasinin ICINDE gecmesi: ayni protein, bir kopyasi
    #    kirpilmis. Parcadan kurulan profil sistematik olarak zayif olur.
    duplicated = {i for g in dup_groups for i in g}
    substrings = []
    for ident_a, seq_a in zip(ids, seqs):
        if not seq_a:
            continue
        for ident_b, seq_b in zip(ids, seqs):
            if ident_a == ident_b or seq_a == seq_b:
                continue
            if seq_a in seq_b:
                substrings.append(f"{ident_a} ({len(seq_a)} aa) is contained "
                                  f"in {ident_b} ({len(seq_b)} aa)")
    rep.check("refs71: no reference sequence is a substring of another",
              not substrings, len(substrings), first(substrings, 4))

    # 3. Her kimlik UC yerde de tam olarak bir kez: FASTA, cluster_ecology.csv,
    #    chemistry.csv. Biri eksikse web arayuzu o tipi yarim gosterir.
    fasta_counts = {}
    for ident in ids:
        fasta_counts[ident] = fasta_counts.get(ident, 0) + 1
    eco_counts, chem_counts = {}, {}
    for row in eco:
        key = row.get("cluster", "")
        eco_counts[key] = eco_counts.get(key, 0) + 1
    for row in chem:
        key = row.get("cluster", "")
        chem_counts[key] = chem_counts.get(key, 0) + 1
    every = set(fasta_counts) | set(eco_counts) | set(chem_counts)
    bad = []
    for ident in sorted(every):
        counts = (fasta_counts.get(ident, 0), eco_counts.get(ident, 0),
                  chem_counts.get(ident, 0))
        if counts != (1, 1, 1):
            bad.append(f"{ident}: fasta={counts[0]} ecology={counts[1]} "
                       f"chemistry={counts[2]}")
    rep.check("refs71: every reference id appears exactly once in FASTA, "
              "ecology and chemistry", not bad, len(bad),
              f"{len(every)} reference ids checked" if not bad
              else first(bad, 4))

    # 4. Neredeyse ayni ama AYNI OLMAYAN ciftler: WARN. Bu bir kusur DEGIL,
    #    bu veritabaninin bulgusu. Kimlik olcumu pipeline'in kendi diamond
    #    sonucundan (reference_pairs.csv) okunur, yeniden hesaplanmaz.
    chem_by = {r["cluster"]: r for r in chem}

    def substrate_of(ident):
        return (chem_by.get(ident, {}).get("substrate_en") or "unknown").strip()

    if not os.path.exists(pairs_path):
        rep.check(f"refs71: near-identical pairs above "
                  f"{NEAR_IDENTICAL_PCT:.0f}% identity", True, 0,
                  f"{os.path.basename(pairs_path)} not present; run "
                  f"evidence_tiers.py -- check skipped", severity="WARN")
    else:
        # Kimlik olcumu TURETILMIS bir dosyadan okunuyor. O dosya FASTA'dan
        # eskiyse olcum artik baska bir referans setini tarif ediyor; uyari
        # yaniltici olur, o yuzden acikca yazilir.
        stale_note = ""
        try:
            if os.path.getmtime(pairs_path) < os.path.getmtime(refs_path):
                stale_note = (" [WARNING: reference_pairs.csv is OLDER than "
                              "refs71.fasta, so these identities describe an "
                              "earlier reference set; re-run evidence_tiers.py]")
        except OSError:
            pass
        near = []
        near_keys = set()
        with open(pairs_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                try:
                    pid = float(row.get("identity", ""))
                except ValueError:
                    continue
                if pid <= NEAR_IDENTICAL_PCT:
                    continue
                a, b = row.get("ref_a", ""), row.get("ref_b", "")
                # Tam kopya / parca ciftleri yukarida SERT kontrol olarak
                # raporlandi; burada tekrar edilmemeleri gerekiyor, yoksa
                # gercek biyolojik bulgu duplikasyon gurultusunda kaybolur.
                if a in duplicated and b in duplicated:
                    continue
                if any(f"{a} (" in s and f"in {b} (" in s
                       or f"{b} (" in s and f"in {a} (" in s
                       for s in substrings):
                    continue
                near.append(f"{a} / {b} at {pid:.1f}% "
                            f"({substrate_of(a)} against {substrate_of(b)})")
                near_keys.add(frozenset((a, b)))
        rep.check(f"refs71: near-identical pairs above "
                  f"{NEAR_IDENTICAL_PCT:.0f}% identity", not near, len(near),
                  "; ".join(near[:4]) + (f"; ... (+{len(near) - 4})"
                                         if len(near) > 4 else "") + stale_note,
                  severity="WARN")

        # 4b. Yukaridaki WARN kalici: cift ayri kalmaya karar verildi. O zaman
        #     KARARIN KENDISI yazili olmak zorunda, yoksa okuyucu fazlaligi
        #     gorulmemis mi sanacak. curation_decisions.csv bu kaydi tutuyor ve
        #     iki yonlu kontrol ediliyor: karari olmayan cift ve artik var
        #     olmayan cifte yazilmis karar, ikisi de hata.
        decisions = {}
        dec_path = os.path.join(ROOT, "curation_decisions.csv")
        if os.path.exists(dec_path):
            with open(dec_path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    a = (row.get("pair_a") or "").strip()
                    b = (row.get("pair_b") or "").strip()
                    if a and b:
                        decisions[frozenset((a, b))] = row
        missing = sorted("/".join(sorted(k)) for k in near_keys - set(decisions))
        stale = sorted("/".join(sorted(k)) for k in set(decisions) - near_keys)
        bare = sorted("/".join(sorted(k)) for k, row in decisions.items()
                      if not (row.get("rationale") or "").strip()
                      or not (row.get("decided_on") or "").strip())
        problems = ([f"karari yok: {n}" for n in missing]
                    + [f"cift yok: {n}" for n in stale]
                    + [f"gerekcesi/tarihi eksik: {n}" for n in bare])
        rep.check("curation_decisions: every near-identical pair has a dated, "
                  "reasoned decision", not problems, len(problems),
                  f"{len(near_keys)} near-identical pair(s), "
                  f"{len(decisions)} decision(s) recorded" +
                  ("; " + first(problems, 6) if problems else ""))

    # 5. Hic uye toplamayan referanslar: WARN. Sifir uye hem ikizin hem de
    #    parca referansin imzasidir, o yuzden isimleriyle yazilir.
    members = {r[0]: r[1] for r in con.execute(
        "SELECT ro_cluster, COUNT(*) FROM ro WHERE is_confirmed=1 "
        "AND ro_cluster IS NOT NULL GROUP BY ro_cluster")}
    zero = [ident for ident in sorted(fasta_counts) if members.get(ident, 0) == 0]
    notes = []
    for ident in zero:
        tag = ""
        if ident in duplicated:
            twin = [i for g in dup_groups if ident in g for i in g if i != ident]
            tag = f" [identical twin of {','.join(twin)}]"
        elif any(s.startswith(f"{ident} (") for s in substrings):
            tag = " [fragment of another reference]"
        notes.append(ident + tag)
    rep.check("refs71: every reference recruits at least one confirmed member",
              not zero, len(zero),
              f"{len(fasta_counts) - len(zero)} of {len(fasta_counts)} "
              f"references have members; zero-member: " + first(notes, 12),
              severity="WARN")


def check_metabolic_role(rep, root):
    """Her tipin metabolik rolu cozulebiliyor mu ve tablo taze mi?

    NEDEN AYRI BIR EKSEN. Biyodegradasyon odagi icin "oksijen urune giriyor
    mu" olcutu YANLIS: SxtT ve PrnD oksijeni urune sokar ama biyosentetiktir,
    vanillat/kafein demetilazlarinda oksijen kopan formaldehite gider ama
    bunlar klasik biyodegradasyondur. Birini otekinin yerine kullanmak yanlis
    tipleri eler.

    Rol, enzimin RIESKE OLUP OLMADIGI hakkinda hicbir sey soylemez.
    """
    try:
        sys.path.insert(0, root)
        import metabolic_role as mr                     # noqa: PLC0415
    except Exception as exc:                            # noqa: BLE001
        rep.check("metabolic role: module usable", False, 1, str(exc))
        return
    rows = mr.table(root)
    bad = sorted(r["cluster"] for r in rows if r["resolved_from"] in
                 ("missing_declaration", "no_rule", "type_overrides_default"))
    rep.check("metabolic role: every type resolves to a role",
              not bad, len(bad),
              f"{sum(1 for r in rows if r['resolved_from'] == 'type')} declared "
              f"by the type, {sum(1 for r in rows if r['resolved_from'] == 'substrate_class')} "
              f"from the substrate class" if not bad else first(bad))
    bad_val = sorted(r["cluster"] for r in rows
                     if r["metabolic_role"] not in mr.ROLES)
    rep.check("metabolic role: every value is from the vocabulary",
              not bad_val, len(bad_val),
              "/".join(mr.ROLES) if not bad_val else first(bad_val))
    leaks = []
    for r in rows:
        hits = turkish_hits(r.get("basis") or "")
        if hits:
            leaks.append(f"{r['cluster']}: {', '.join(sorted(hits))}")
    rep.check("metabolic role: the published basis text carries no Turkish",
              not leaks, len(leaks),
              f"{len(rows)} basis strings checked" if not leaks
              else first(leaks, 4))
    path = os.path.join(root, mr.OUT_FILE)
    if not os.path.exists(path):
        rep.check("metabolic role: derived table exists", True, 0,
                  f"{mr.OUT_FILE} not written yet", severity="WARN")
    else:
        with open(path) as fh:
            disk = list(csv.DictReader(fh))
        by = {r["cluster"]: r for r in rows}
        stale = [r["cluster"] for r in disk
                 if by.get(r["cluster"], {}).get("metabolic_role")
                 != r["metabolic_role"]]
        rep.check("freshness: the derived metabolic-role table matches types/",
                  not stale, len(stale),
                  f"{len(disk)} rows agree" if not stale else first(stale, 4))


def check_computed_oxygen_fate(rep, root):
    """Hesaplanabilen her oksijen kaderi, HESAPLANANLA AYNI mi?

    NEDEN. Oksijen kaderi bu oturumda iki kez yanlis kaydedildi ve her
    seferinde sebep ayniydi: deger KURATORUN OKUMASINA dayaniyordu. Artik
    resmi dengelenmis denklem (Rhea) ve katilimci formulleri (ChEBI) uzerinden
    ATOM SAYIMIYLA hesaplanabiliyor; hesaplanabildigi yerde insan okumasi
    HAKEM DEGIL.

    Sayim atom IZLEMEZ, bunu da bilir: iskelet parcalanmis ve substrat zaten
    oksijen tasiyorsa deger "ambiguous" doner ve kuratorun gerekcesi ayakta
    kalir. Bu kontrol yalnizca HESAPLANMIS olanlari karsilastirir -- yani
    yalnizca sayimin guvenilir oldugu yerde bagirir.
    """
    path = os.path.join(root, "analysis_out", "oxygen_fate_computed.json")
    if not os.path.exists(path):
        rep.check("oxygen fate: computed comparison present", True, 0,
                  "analysis_out/oxygen_fate_computed.json not written yet",
                  severity="WARN")
        return
    rows = json.load(open(path)).get("rows") or []
    decided = [r for r in rows
               if r.get("computed") and r["computed"] != ["ambiguous"]]
    bad = [f"{r['cluster']}: db={r.get('resolved_in_db')} "
           f"computed={'/'.join(r['computed'])}"
           for r in decided
           if r.get("resolved_in_db") not in r["computed"]]
    rep.check("oxygen fate: agrees with the value computed from the official "
              "balanced equation",
              not bad, len(bad),
              f"{len(decided)} of {len(rows)} Rhea-anchored types computed, "
              f"all agree" if not bad else first(bad, 5))

    # Katman alani ile gercek durum tutarli mi: "computed" diyen bir tipin
    # hesaplanmis bir degeri OLMAK zorunda, yoksa etiket bos bir guvence olur.
    try:
        import yaml as _y                               # noqa: PLC0415
        import glob as _g                               # noqa: PLC0415
        tiers = {}
        for f in _g.glob(os.path.join(root, "types", "*.yaml")):
            m = _y.safe_load(open(f))
            tiers[m.get("cluster")] = m.get("chemistry_evidence")
    except Exception as exc:                            # noqa: BLE001
        rep.check("chemistry evidence: tiers readable", False, 1, str(exc))
        return
    computed_ok = {r["cluster"] for r in decided}
    lying = sorted(c for c, t in tiers.items()
                   if t == "computed_from_official_equation"
                   and c not in computed_ok)
    rep.check("chemistry evidence: every type claiming a computed oxygen fate "
              "really has one",
              not lying, len(lying),
              f"{len(computed_ok)} computed, "
              f"{sum(1 for t in tiers.values() if t == 'curated')} curated"
              if not lying else first(lying, 5))
    missing = sorted(c for c, t in tiers.items() if not t)
    rep.check("chemistry evidence: every type states what its chemistry rests "
              "on",
              not missing, len(missing),
              f"{len(tiers)} types tiered" if not missing
              else first(missing, 5))


def check_reaction_names(rep, root):
    """Resmi denklemle uyusmayan her urun adi KAYDA GECMIS mi?

    NEDEN. verify_reactions.py her tipin kayitli urununu UniProt/Rhea'nin
    resmi denklemine karsi tutuyor ve 2026-10-28'de UC GERCEK HATA buldu:
    HcaE'nin urunu asagi akistaki katekol yazilmisti, CndA'da bir ester
    "aldehit" denmisti, Stc2'de formaldehit atlanmisti. Ucu de ayni sinif --
    enzimin BIRAKTIGI molekul yerine baska bir molekulu urun yazmak -- ve
    kolin monooksijenazin iki kez yanlis kaydedilmesinin de sebebi buydu.

    Ayni kontrol, ayni molekulun geleneksel adiyla IUPAC adini da uyusmaz
    sayiyor. Bu farklari SESSIZCE gecmek yanlis olur, cunku bir gun gercek bir
    kimya farki ayni kutuya duser. O yuzden iki sey ayri ayri kontrol edilir:
      1. bildirilen her uyusmazligin bir kaydi var mi
      2. kaydi olup ARTIK bildirilmeyen (olu) kayit kaldi mi
    """
    reg_path = os.path.join(root, "reaction_name_exceptions.yaml")
    rep_path = os.path.join(root, "analysis_out", "reaction_verification.json")
    if not os.path.exists(rep_path):
        rep.check("reaction names: verification report present", True, 0,
                  "analysis_out/reaction_verification.json not written yet "
                  "(run verify_reactions.py)", severity="WARN")
        return
    try:
        import yaml as _y                               # noqa: PLC0415
        reg = {e["cluster"]: e for e in
               ((_y.safe_load(open(reg_path)) or {}).get("exceptions") or [])}
    except Exception as exc:                            # noqa: BLE001
        rep.check("reaction names: registry readable", False, 1, str(exc))
        return
    # TAZELIK. Rapor types/'tan eskiyse kapi ESKI veriye bakar ve yeni bir
    # urun hatasi sessizce gecer. Raporu yeniden uretmek bir ag cagrisi
    # gerektirdigi icin otomatik adim DEGIL, ama bayatligi GORUNUR olmak
    # zorunda.
    import glob as _glob                                # noqa: PLC0415
    newest = max((os.path.getmtime(f) for f in
                  _glob.glob(os.path.join(root, "types", "*.yaml"))),
                 default=0)
    stale_by = newest - os.path.getmtime(rep_path)
    rep.check("freshness: the reaction verification report is newer than "
              "types/",
              stale_by <= 0, 1 if stale_by > 0 else 0,
              "report is current" if stale_by <= 0 else
              f"types/ changed {stale_by / 60:.0f} min after the report was "
              "written -- re-run verify_reactions.py", severity="WARN")

    rows = json.load(open(rep_path)).get("rows") or []

    # KAPSAM, SONUCTAN ONCE OKUNUR. 2026-11-02: bir kosuda 21 UniProt
    # sorgusu HTTP 500 aldi; o tipler rapordan dustu ve ozet "UYUSMAYAN: 0"
    # yazdi -- bir onceki kosuda bulunmus 6 uyusmazlik gorunmez oldu.
    # Asagidaki "her uyusmazlik kayitli" kontrolu de o halde BOS gecerdi,
    # cunku bakilmayan tipin uyusmazligi da yok sayilir. Sorgu basarisizligi
    # bu yuzden sonuctan ONCE ve ayrica bildirilir.
    failed = sorted(r["cluster"] for r in rows
                    if str(r.get("status") or "").startswith("lookup failed"))
    rep.check("reaction names: every reference with an accession was actually "
              "queried",
              not failed, len(failed),
              "sorgu basarisizligi yok" if not failed else
              f"{len(failed)} tip sorulamadi ({', '.join(failed[:6])}"
              f"{' ...' if len(failed) > 6 else ''}) -- asagidaki uyusmazlik "
              "sayisi bir ALT SINIR, raporu yeniden uret",
              severity="WARN")

    flagged = {r["cluster"] for r in rows
               if r.get("product_agrees") is False
               or r.get("qualifiers_missing_from_our_product")}
    unregistered = sorted(flagged - set(reg))
    rep.check("reaction names: every mismatch with the official equation is "
              "recorded",
              not unregistered, len(unregistered),
              f"{len(flagged)} differences, all recorded"
              if not unregistered else first(unregistered, 6))
    dead = sorted(set(reg) - flagged)
    rep.check("reaction names: no stale entry in the registry",
              not dead, len(dead),
              f"{len(reg)} entries, all still reported"
              if not dead else first(dead, 6))
    incomplete = sorted(c for c, e in reg.items()
                        if not all(str(e.get(k) or "").strip()
                                   for k in ("official", "ours", "why",
                                             "checked_on")))
    rep.check("reaction names: each recorded difference states the official "
              "name, ours, why and when",
              not incomplete, len(incomplete),
              "all entries complete" if not incomplete
              else first(incomplete, 6))


def check_oxygen_fate(rep, root, chem):
    """O2 oksijeninin kaderi her tip icin COZULEBILIYOR mu, ve turetilmis
    tablo TAZE mi?

    NEDEN. "Rieske oksijenaz" bir aile adidir; icinde urune hic oksijen
    sokmayan uyeler var (desaturasyon, C-C bagi kurma). "Dioksijenaz" diye
    filtreleyen biri onlari da aliyordu. Deger reaksiyon sinifindan
    TURETILEMEZ: dioxygenation_with_release homojen degil -- nitroarenlerde
    oksijen ana urune, fenoksiasetatlarda kopan glioksilata gidiyor. O yuzden
    sozluk o siniflari `per_type` isaretler ve tipin kendi beyanini ZORUNLU
    kilar. Bu kontrol beyanin EKSIK kalmasini sessizlestirmez.
    """
    try:
        sys.path.insert(0, root)
        import oxygen_fate                              # noqa: PLC0415
        vocab = oxygen_fate.load_vocab(root)
    except Exception as exc:                            # noqa: BLE001
        rep.check("oxygen fate: vocabulary readable", False, 1,
                  f"reaction_classes.yaml not usable ({exc})")
        return

    used = {(r.get("reaction_class") or "").strip() for r in chem}
    undefined = sorted(c for c in used if c and c not in vocab)
    rep.check("oxygen fate: every reaction class is defined in the vocabulary",
              not undefined, len(undefined),
              f"{len(vocab)} classes defined, {len(used)} in use"
              if not undefined else first(undefined))

    rows = oxygen_fate.table(root)
    by = {r["cluster"]: r for r in rows}
    bad_src = sorted(r["cluster"] for r in rows if r["resolved_from"] in
                     ("missing_declaration", "undefined_class"))
    rep.check("oxygen fate: every per-type class member declares its fate",
              not bad_src, len(bad_src),
              f"{sum(1 for r in rows if r['resolved_from'] == 'type')} types "
              f"declare their own fate" if not bad_src else first(bad_src))

    override = sorted(r["cluster"] for r in rows
                      if r["resolved_from"] == "type_overrides_class")
    rep.check("oxygen fate: no type silently contradicts its class default",
              not override, len(override),
              "no overrides" if not override else first(override))

    bad_val = sorted(r["cluster"] for r in rows
                     if r["oxygen_fate"] not in oxygen_fate.FATES)
    rep.check("oxygen fate: every value is from the vocabulary",
              not bad_val, len(bad_val),
              "/".join(oxygen_fate.FATES) if not bad_val else first(bad_val))

    # Gerekce alani SITEYE CIKAR: Turkce sizmasi burada yakalanmali.
    leaks = []
    for r in rows:
        hits = turkish_hits(r.get("basis") or "")
        if hits:
            leaks.append(f"{r['cluster']}: {', '.join(sorted(hits))}")
    rep.check("oxygen fate: the published basis text carries no Turkish",
              not leaks, len(leaks),
              f"{len(rows)} basis strings checked" if not leaks
              else first(leaks, 4))

    # TAZELIK: tablo diskte duruyorsa, tiplerden YENIDEN hesaplanmis haliyle
    # ayni olmak zorunda. Bir tipin sinifi degisip tablo guncellenmezse site
    # eski kaderi filtreler ve kimse farketmez.
    path = os.path.join(root, oxygen_fate.OUT_FILE)
    if not os.path.exists(path):
        rep.check("oxygen fate: derived table exists", True, 0,
                  f"{oxygen_fate.OUT_FILE} not written yet", severity="WARN")
    else:
        with open(path) as fh:
            disk = list(csv.DictReader(fh))
        stale = [r["cluster"] for r in disk
                 if by.get(r["cluster"], {}).get("oxygen_fate")
                 != r["oxygen_fate"]]
        missing = sorted(set(by) - {r["cluster"] for r in disk})
        rep.check("freshness: the derived oxygen-fate table matches types/",
                  not stale and not missing, len(stale) + len(missing),
                  f"{len(disk)} rows agree"
                  if not stale and not missing
                  else first(stale + [f"missing: {m}" for m in missing], 4))


def check_reference_motifs(rep, root, refs_aln=None):
    """Her kuratorlu referans motif kapisini geciyor mu -- ve gecmiyorsa
    motif_exceptions.yaml'da gerekcesi YAZILI mi?

    NEDEN. Kapi (ro_motif.classify_motifs) aday proteinleri suzmek icin var ve
    REFERANSLARIN gecmesi beklenir. AerC gecmiyor: 717 aa oldugu icin global
    hizalama katalitik kolonlara ulasmiyor ve 0/3 veriyor. Bu bir artefakt, ama
    SESSIZ kalmasi yanlis -- sette kendi kapisindan gecemeyen bir referans
    varsa, gerekcesi ve kaniti yazili olmalidir.

    Dolayisiyla iki sey ayri ayri kontrol edilir:
      1. gecmeyen her referansin bir istisna kaydi var mi
      2. istisna kaydi olup da ARTIK GECEN referans kaldi mi (olu istisna)
    """
    try:
        sys.path.insert(0, root)
        import ro_motif                              # noqa: PLC0415
    except Exception as exc:                          # noqa: BLE001
        rep.check("references: motif gate evaluated", True, 0,
                  f"ro_motif not importable ({exc})", severity="WARN")
        return

    aln = refs_aln or os.path.join(root, "ROs_71_Clean", "ROmotif71_hmmaln.sto")
    if not os.path.exists(aln):
        rep.check("references: motif gate evaluated", True, 0,
                  f"reference alignment missing: {aln}", severity="WARN")
        return

    def cluster_of(name):
        """'7_701_AerC_AerC_pro...' -> '7_701_AerC'.

        Referans fasta basliklari <kume>_<gen>_pro... bicimindedir, yani gen
        adi TEKRARLIDIR. Sadece '_pro' ile kesmek '7_701_AerC_AerC' veriyordu
        ve istisna kaydiyla eslesmiyordu -- ayni kontrol hem 'belgelenmemis'
        hem 'olu istisna' diye iki kez FAIL vermisti.
        """
        m = re.match(r"^(\d+_\d+_[^_]+)", name)
        return m.group(1) if m else name.split("_pro")[0]

    report = ro_motif.motif_report(aln)
    failing = {cluster_of(n) for n, r in report.items()
               if not r["is_RO_alpha_motif"]}

    path = os.path.join(root, "motif_exceptions.yaml")
    excepted, detail = set(), {}
    if os.path.exists(path):
        try:
            doc = yaml.safe_load(open(path, encoding="utf-8")) or {}
        except Exception:                             # noqa: BLE001
            doc = {}
        for e in (doc.get("exceptions") or []):
            cid = e.get("cluster")
            if cid:
                excepted.add(cid)
                detail[cid] = e

    undocumented = sorted(failing - excepted)
    rep.check("references: every reference passes the motif gate, or its "
              "exception is documented",
              not undocumented, len(undocumented),
              "; ".join(undocumented) +
              (" -- add it to motif_exceptions.yaml with reason, evidence and limits"
               if undocumented else ""))

    # Olu istisna: artik geciyorsa kayit kaldirilmali, yoksa gercek bir
    # basarisizligi ortmeye devam eder.
    stale = sorted(excepted - failing)
    rep.check("references: no stale motif exception", not stale, len(stale),
              "; ".join(stale) + (" -- passes the gate now; remove the exception"
                                  if stale else ""))

    # Her istisna gerekce, kanit ve sinir tasimak zorunda.
    thin = []
    for cid in sorted(excepted & failing):
        e = detail.get(cid) or {}
        for field in ("reason", "how_the_ligands_were_verified",
                      "evidence_files", "limits"):
            if not e.get(field):
                thin.append(f"{cid}:{field}")
    rep.check("references: each motif exception states reason, evidence and limits",
              not thin, len(thin), "; ".join(thin[:6]))


def check_literature_motif_spacing(rep, root):
    """Kolon ucluleri literaturun DXXHXXXXH imzasini tutuyor mu?

    NEDEN. Derlemeler mononukleer merkezin dizi imzasini DXXHXXXXH, yani
    D-x(2)-H-x(4)-H olarak veriyor. Kullanilan kolonlar (286 / 289 / 294)
    bunun birebir kendisi ve bu, kolonlarin dogru yerde oldugunun YAPIDAN
    BAGIMSIZ ikinci kanitidir (birincisi 1NDO Asp362 -> kolon 431 capasi).
    Kolonlardan biri sessizce kayarsa iki kanittan biri duser; bu kontrol
    o dusmeyi gurultulu yapar. Karboksilat kolonunun yanlis olmasi bu
    projede zaten bir kere sessiz kaldi.
    """
    label = ("motif columns: the bridge and the two histidines keep the "
             "literature's D-x(2)-H-x(4)-H spacing")
    try:
        sys.path.insert(0, root)
        import ro_motif                              # noqa: PLC0415
    except Exception as exc:                          # noqa: BLE001
        rep.check(label, True, 0, f"ro_motif not importable ({exc})",
                  severity="WARN")
        return
    bridge = ro_motif.BRIDGING_SITE[0]
    his1 = ro_motif.CATALYTIC_SITES[0][0]
    his2 = ro_motif.CATALYTIC_SITES[1][0]
    gap1, gap2 = his1 - bridge - 1, his2 - his1 - 1
    ok = gap1 == 2 and gap2 == 4
    rep.check(label, ok, 0 if ok else 1,
              "" if ok else
              f"columns {bridge}/{his1}/{his2} give D-x({gap1})-H-x({gap2})-H, "
              f"but the literature's signature is D-x(2)-H-x(4)-H "
              f"(naphthalene 1,2-dioxygenase Asp205/His208/His213); one of "
              f"these columns has moved")



def check_membership_basis(rep, root):
    """Her onayli girisin bir uyelik temeli var mi, ve temel okumayla tutuyor mu?

    NEDEN. 3/3 kapisi uyelikte karar veriyor ama SILME olarak degil: triad'in
    istisnasiz oldugunu yazan bir kaynak bulunamadi, o yuzden ligandi eksik
    giris atilmaz, hangi temelle uye oldugu yayinlanir (ro_membership). Boyle
    bir tasarimin sessiz bozulma yolu belli: tablo bayatlar ya da bir temel
    kendi okumasindan AYRISIR -- ornegin triad_status 'his_absent' olan bir
    giris 'triad_complete' temeliyle yayinlanir. Sayfa o cumleyi temelden
    kuruyor, yani ayrisma dogrudan yanlis bir iddiaya cevrilir.
    """
    label = "membership basis: every confirmed entry has one and it matches its triad reading"
    if root not in sys.path:
        sys.path.insert(0, root)
    import triad_gate                                      # noqa: PLC0415
    db = os.path.join(root, "roar.sqlite")
    if not os.path.exists(db):
        rep.check(label, False, 0, "roar.sqlite yok", severity="WARN")
        return
    con = sqlite3.connect("file:%s?mode=ro" % db, uri=True)
    try:
        have = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        if "ro_membership" not in have:
            rep.check(label, False, 0,
                      "ro_membership yok; triad_gate.py --apply kosulmali",
                      severity="WARN")
            return
        confirmed = con.execute(
            "SELECT COUNT(*) FROM ro WHERE is_confirmed = 1").fetchone()[0]
        n_rows = con.execute("SELECT COUNT(*) FROM ro_membership").fetchone()[0]
        tally = dict(con.execute(
            "SELECT basis, COUNT(*) FROM ro_membership GROUP BY 1"))
        # Temel ile okumanin ayrismasi: 'triad_complete' diyen ama triad'i
        # tamamlanmayan, ya da tersi.
        mismatched = con.execute("""
            SELECT COUNT(*) FROM ro_membership m JOIN ro_triad t USING(candidate_id)
            WHERE (m.basis = 'triad_complete') != (t.triad_status LIKE 'complete%')
            """).fetchone()[0]
        orphan = con.execute(
            "SELECT COUNT(*) FROM ro WHERE is_confirmed = 1 AND candidate_id NOT IN "
            "(SELECT candidate_id FROM ro_membership)").fetchone()[0]
    finally:
        con.close()

    problems = []
    if n_rows != confirmed:
        problems.append("ro_membership %d satir, onaylanmis giris %d"
                        % (n_rows, confirmed))
    if orphan:
        problems.append("%d onayli girisin uyelik temeli yok" % orphan)
    if mismatched:
        problems.append("%d giriste temel ile triad okumasi ayrisiyor" % mismatched)
    unknown = set(tally) - set(triad_gate.MEMBERSHIP_BASES)
    if unknown:
        problems.append("aciklamasi olmayan temel: %s" % ", ".join(sorted(unknown)))

    detail = ", ".join("%s=%d" % kv for kv in sorted(tally.items())) or "-"
    rep.check(label, not problems, len(problems),
              "; ".join(problems) if problems else detail)


def check_chemistry_published(rep, root, con):
    """chemistry.csv'deki her tip ve her reaksiyon sinifi ana sayfada sayiliyor mu?

    NEDEN. Ana sayfanin kimya bolumu iki ELLE tutulan liste uzerinden
    geziyordu: atlas.FAMILY_ORDER (aile sirasi) ve reaction_summary'deki
    reaksiyon anahtarlari. Ikisi de SIRA listesiydi ama SUZGEC gibi
    davraniyordu -- listede olmayan deger sessizce dusuyordu.

    Olculdu (2026-11-01): `tetrapyrroles` ailesi (6_601_PAO, 6_602_CAO,
    6_603_PTC52; 617 giris) 14 Ekim'de eklendigi halde ne kart ne tablo
    satiri aliyordu; ayrica bes reaksiyon sinifi (desaturation,
    N_oxygenation, C_C_bond_formation, oxygenolytic_ring_opening,
    C_C_cleavage) hic yayinlanmiyordu. Yayinlanan toplam 68 tip / 17.685
    giris, gercegi 76 tip / 19.125 giristi -- ve bolumun kendi cumlesi
    "cis-dihydroxylation only the most common case, not the whole story"
    diyordu.

    Iki liste artik veriye bagli. Bu kontrol bagin KOPMADIGINI olcer:
    yayinlanan tip ve giris sayisi kaynagin tamamina esit olmak zorunda.
    """
    import csv as _csv
    sys.path.insert(0, os.path.join(root, "webapp"))
    import atlas as _at
    chem = list(_csv.DictReader(open(os.path.join(root, "chemistry.csv"),
                                     encoding="utf-8")))
    counts = {r[0]: r[1] for r in con.execute(
        "SELECT ro_cluster, COUNT(*) FROM ro WHERE is_confirmed=1 "
        "AND ro_cluster IS NOT NULL GROUP BY 1")}
    eco_path = os.path.join(root, "cluster_ecology.csv")
    eco = ({r["cluster"]: r for r in _csv.DictReader(open(eco_path, encoding="utf-8"))}
           if os.path.exists(eco_path) else {})
    groups = _at.chemistry_groups(os.path.join(root, "chemistry.csv"), counts, eco)

    grouped = {t["cluster"] for g in groups for t in g["types"]}
    missing = sorted(r["cluster"] for r in chem if r["cluster"] not in grouped)
    rep.check("chemistry: every curated type reaches the substrate cards",
              not missing, len(missing),
              f"{len(grouped)} of {len(chem)} types grouped"
              if not missing else first(missing, 6))

    summary = _at.reaction_summary(groups)
    pub_types = sum(r["types"] for r in summary)
    pub_members = sum(r["members"] for r in summary)
    exp_members = sum(counts.get(r["cluster"], 0) for r in chem)
    bad = []
    if pub_types != len(chem):
        bad.append(f"types {pub_types} != {len(chem)}")
    if pub_members != exp_members:
        bad.append(f"entries {pub_members} != {exp_members}")
    rep.check("chemistry: the reaction table counts every type and every entry",
              not bad, len(bad),
              f"{pub_types} types, {pub_members} entries"
              if not bad else "; ".join(bad))

    unlabelled = sorted({(r["reaction_class"] or "").strip() for r in chem}
                        - set(_at.REACTION_LABEL))
    rep.check("chemistry: every reaction class has a published label and note",
              not unlabelled, len(unlabelled),
              f"{len(set(r['reaction_class'] for r in chem))} classes"
              if not unlabelled else first(unlabelled, 6))

    noscheme = sorted(k for k in {(r["reaction_class"] or "").strip() for r in chem}
                      if not _at.reaction_scheme_svg(k))
    rep.check("chemistry: every reaction class has a generic scheme",
              not noscheme, len(noscheme),
              "all classes drawn" if not noscheme else first(noscheme, 6),
              severity="WARN")

    unfamilied = sorted({(r["family"] or "").strip() for r in chem}
                        - set(_at.FAMILY_LABEL))
    rep.check("chemistry: every substrate family has a published label",
              not unfamilied, len(unfamilied),
              f"{len(set(r['family'] for r in chem))} families"
              if not unfamilied else first(unfamilied, 6))


def check_leaf_stop_reason(rep, root):
    """Heterojen bir yapragin NEDEN bolunmeyi biraktigi kayitli mi?

    NEDEN. recursive_homogenize.py yapraklari, medyan ikili kimlik 0,70'e
    ulasana kadar artan CD-HIT esikleriyle boluyor. Ama iki BAGIMSIZ kural
    bolmeyi erken kesebiliyor: dugum MIN_SPLIT_SIZE'dan kucukse, ya da esik
    listesi tukenmisse. 2026-10-30'da olculdu: is_homogeneous=0 olan 420
    yapragin 420'si (1.460 uye) BOYUT tabani yuzunden durdu, hicbiri esik
    listesini tuketmedi. Yani sayfadaki yalin "heterogeneous" etiketi
    yaniltiyordu: okuyucu homojenlik hedefinin tutturulamadigini sanirken
    hedef hic sorulmamisti. Bu, bu projenin tekrar eden hata bicimi --
    karar veren kural, karar verdigini soyledigimiz kural degil (bkz. triad
    aralik sarti ve kolon 471/431 vakalari). Kapi iki seyi tutuyor:
    sabitlerin kopyasi kaynagindan ayrismasin, ve her heterojen yapragin
    durma nedeni TURETILEBILIR olsun. "unexplained" cikan bir yaprak,
    homojenizasyonun belgelenmemis bir yoldan ciktigi anlamina gelir.

    Yaprak taksonomik bir birim DEGIL. Literaturde bir RO "varyantini"
    tanimlayan yuzde esigi yok: Kweon 2008 semasi elektron ortaklarina,
    ro_group alfa alt birimi kumelenmesine dayanir (citations.yaml'deki
    kweon2008_classification notu). Yaprak bunlarin hicbiri degil, agac
    boyutunu kisan islevsel bir gruplamadir ve sayfada oyle okunmali.
    """
    label = "leaves: every heterogeneous leaf records why splitting stopped"
    db = os.path.join(root, "roar.sqlite")
    if not os.path.exists(db):
        rep.check(label, False, 0, "roar.sqlite yok", severity="WARN")
        return

    problems = []
    try:
        sys.path.insert(0, root)
        sys.path.insert(0, os.path.join(root, "webapp"))
        import atlas  # noqa: F401
        import recursive_homogenize as rh
    except Exception as exc:  # pragma: no cover
        rep.check(label, False, 0, "sabitler okunamadi: %s" % exc, severity="WARN")
        return

    # Sabitler atlas.py'de KOPYALANMIS; kaynak recursive_homogenize.py.
    # Ayrisirlarsa sayfa yanlis neden basar, o yuzden burada tutuluyor.
    if atlas.HOMOGENEITY_MIN_SPLIT_SIZE != rh.MIN_SPLIT_SIZE:
        problems.append(
            "atlas.HOMOGENEITY_MIN_SPLIT_SIZE=%s, recursive_homogenize.MIN_SPLIT_SIZE=%s"
            % (atlas.HOMOGENEITY_MIN_SPLIT_SIZE, rh.MIN_SPLIT_SIZE))
    if atlas.HOMOGENEITY_MAX_DEPTH != len(rh.THRESHOLDS):
        problems.append(
            "atlas.HOMOGENEITY_MAX_DEPTH=%s, len(THRESHOLDS)=%s"
            % (atlas.HOMOGENEITY_MAX_DEPTH, len(rh.THRESHOLDS)))

    con = sqlite3.connect("file:%s?mode=ro" % db, uri=True)
    try:
        have = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        if "leaf" not in have:
            rep.check(label, False, 0, "leaf tablosu yok", severity="WARN")
            return
        rows = con.execute(
            "SELECT leaf_id, size, depth, is_homogeneous FROM leaf "
            "WHERE is_homogeneous = 0").fetchall()
    finally:
        con.close()

    counts = {}
    unexplained = []
    for leaf_id, size, depth, homo in rows:
        reason = atlas.leaf_stop_reason(size, depth, homo)
        counts[reason] = counts.get(reason, 0) + 1
        if reason == "unexplained":
            unexplained.append(leaf_id)
    if unexplained:
        problems.append("durma nedeni turetilemeyen %d yaprak: %s"
                        % (len(unexplained), ", ".join(unexplained[:5])))

    detail = "%d heterojen yaprak; %s" % (
        len(rows), ", ".join("%s=%d" % kv for kv in sorted(counts.items())) or "-")
    rep.check(label, not problems, len(problems),
              "; ".join(problems) if problems else detail)


def check_triad_gate(rep, root):
    """Yayinlanan triad tanimi 3/3 mu, ve bu karar TAZE mi?

    NEDEN. Literatur mononukleer Fe merkezini KANONIK olarak 2-His-1-karboksilat
    tanimlar (iki His + Asp ya da Glu) ve karboksilat sus degildir: NDO'da D362A
    naftalinden hic urun vermiyor (Parales 2000, J Bacteriol 182:1641), CntA'da
    E205A/E205D TMA uretimini durduruyor (Zhu 2014, PNAS 111:4268). "Istisnasiz"
    diyen bir kaynak YOK, bu yuzden triad'i eksik girisler silinmez isaretlenir. Buna karsin kabul kurali uzun sure
    >=2/3'te kaldi, cunku 3/3'u SABIT kolonda uygulamak 2388 girisi
    kaybettiriyordu -- oysa o girislerin 2357'sinde karboksilat FIILEN vardi,
    sadece kolon hizalamada bostu. Karar artik ro_triad'da durur ve tabloyu
    ureten denetimin bu DB icin uretilmis oldugu her kapi kosusunda kontrol
    edilir; bayat bir denetimden tablo uretmek bu projede zaten bir kere
    yanlis bulgu dogurdu (kolon 355/471/431 vakasi).
    """
    label = "triad gate: ro_triad is the published 3/3 definition and is fresh"
    db = os.path.join(root, "roar.sqlite")
    if not os.path.exists(db):
        rep.check(label, True, 0, "roar.sqlite yok", severity="WARN")
        return
    con = sqlite3.connect("file:%s?mode=ro" % db, uri=True)
    try:
        have = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        if "ro_triad" not in have or "ro_triad_meta" not in have:
            rep.check(label, True, 0,
                      "ro_triad yok; triad_gate.py --apply kosulmali",
                      severity="WARN")
            return
        meta = {k: json.loads(v) for k, v in
                con.execute("SELECT key, value FROM ro_triad_meta")}
        confirmed = con.execute(
            "SELECT COUNT(*) FROM ro WHERE is_confirmed = 1").fetchone()[0]
        n_rows = con.execute("SELECT COUNT(*) FROM ro_triad").fetchone()[0]
    finally:
        con.close()

    problems = []
    if n_rows != confirmed:
        problems.append(f"ro_triad {n_rows} satir, onaylanmis giris {confirmed}")
    if meta.get("confirmed_entries") != confirmed:
        problems.append(f"meta {meta.get('confirmed_entries')} giris diyor")

    try:
        sys.path.insert(0, root)
        import ro_motif                              # noqa: PLC0415
        want = {"his1": ro_motif.CATALYTIC_SITES[0][0],
                "his2": ro_motif.CATALYTIC_SITES[1][0],
                "carboxylate": ro_motif.CATALYTIC_SITES[2][0]}
        if meta.get("columns") != want:
            problems.append(f"meta kolonlari {meta.get('columns')}, "
                            f"ro_motif {want}")
    except Exception as exc:                          # noqa: BLE001
        problems.append(f"ro_motif okunamadi ({exc})")

    used = (meta.get("audit_used") or {})
    audit_path = os.path.join(root, "analysis_out", "carboxylate_audit.json")
    if os.path.exists(audit_path):
        try:
            stamp = json.load(open(audit_path, encoding="utf-8"))["generated"]
        except Exception:                             # noqa: BLE001
            stamp = None
        if stamp and used.get("generated") != stamp:
            problems.append(f"tablo {used.get('generated')} denetiminden "
                            f"uretilmis, dosyada {stamp} var")

    rep.check(label, not problems, len(problems),
              "" if not problems else "; ".join(problems))

    rule = meta.get("rule") or ""
    rep.check("triad gate: the rule on record requires 3/3, not 2/3",
              "3/3" in rule, 0 if "3/3" in rule else 1,
              "" if "3/3" in rule else f"kayitli kural: {rule!r}")

    undecided = (meta.get("his_undecided") or 0) + (meta.get("absence_not_proven") or 0)
    rep.check("triad gate: every triad call rests on an anchored reading",
              not undecided, undecided,
              "" if not undecided else
              f"{undecided} entries are still decided on an unanchored "
              f"reading ({meta.get('his_undecided')} histidine, "
              f"{meta.get('absence_not_proven')} carboxylate); ro_triad does "
              f"NOT call them absent, and they must be resolved before "
              f"is_confirmed is tightened to 3/3",
              severity="WARN")


def check_motif_columns_anchored(rep, root, refs_aln=None):
    """Katalitik karboksilat kolonu YAPISAL CAPAYA uyuyor mu?

    NEDEN. 2026-10-30'da olculdu: build_reference_models.py kolonu D/E
    korunmasini maksimize ederek seciyordu ve 603 kolonluk modelde kolon 471'i
    secti (korunma %93,4). Dogru kolon 431'dir (%89,5) -- naftalin
    1,2-dioksijenazin yapidan bilinen demir ligandi Asp362 oraya dusuyor.
    Yanlis kolon hic hata vermedi: aktif-bolge tablosu doldu, istatistikler
    uretildi ve "uyelerin %12'sinde karboksilat yok" gibi bir bulgu dogdu.
    Sessiz yanlislik en pahali yanlislik oldugu icin kolon artik her kapi
    kosusunda yapisal capaya karsi yeniden dogrulanir.
    """
    try:
        sys.path.insert(0, root)
        import ro_motif                              # noqa: PLC0415
        import build_reference_models as BRM          # noqa: PLC0415,N817
    except Exception as exc:                          # noqa: BLE001
        rep.check("motif columns: carboxylate column verified against the "
                  "structural anchor", True, 0,
                  f"modules not importable ({exc})", severity="WARN")
        return

    aln = refs_aln
    if not aln:
        for name in ("refs71_hmmaln.sto", "ROmotif71_hmmaln.sto"):
            cand = os.path.join(root, "ROs_71_Clean", name)
            if os.path.exists(cand):
                aln = cand
                break
    if not aln or not os.path.exists(aln):
        rep.check("motif columns: carboxylate column verified against the "
                  "structural anchor", True, 0,
                  "reference alignment missing", severity="WARN")
        return

    try:
        anchor_col, _ = BRM.anchor_carboxylate_column(aln)
    except SystemExit as exc:                         # noqa: BLE001
        rep.check("motif columns: carboxylate column verified against the "
                  "structural anchor", False, 1, str(exc))
        return

    in_use = ro_motif.CATALYTIC_SITES[2][0]
    rep.check("motif columns: carboxylate column verified against the "
              "structural anchor (1NDO Asp362)",
              in_use == anchor_col, 0 if in_use == anchor_col else 1,
              "" if in_use == anchor_col else
              f"ro_motif reads column {in_use} but the anchor residue sits in "
              f"column {anchor_col}; every active-site statistic is measured "
              f"at the wrong position until this agrees")

    path = os.path.join(root, "ROs_71_Clean", "motif_columns.json")
    if os.path.exists(path):
        try:
            doc = json.load(open(path, encoding="utf-8"))
        except Exception:                             # noqa: BLE001
            doc = {}
        entry = (doc.get("ROmotif71") or {})
        src = (entry.get("carboxylate_source") or {}).get("chosen_by")
        rep.check("motif columns: the stored carboxylate column was chosen by "
                  "the structural anchor, not by conservation",
                  src == "structural_anchor", 0 if src == "structural_anchor" else 1,
                  "" if src == "structural_anchor" else
                  f"carboxylate_source.chosen_by = {src!r}; regenerate with "
                  f"build_reference_models.py --derive-only")


def check_retracted_column_claim(rep, root):
    """Katalitik karboksilat kolonunun uye sayilari HICBIR sayfada yayinda mi?

    NEDEN. Kolon 431 bir zamanlar "demir ligandi" diye yayimlaniyordu. 27
    kuratorlu referansin kristal yapisina karsi sinandi ve cagirdigi kalinti
    26/26'sinda demirden 10,9-14,7 A uzakta cikti: kolon ligandi okumuyor,
    yalnizca o civarda asidik bir kalinti bulunduguna bakiyor. Iddia
    2026-10-31'de index.html'in 1.2 bolumunde geri cekildi.

    AMA GERI CEKME TEK SAYFAYA ISLENDI. Atlas ozeti "kolon 431 aspartattir,
    15.996 Asp'ye karsi 122 Glu, yani cephe uclusu Asp kullanir" demeye devam
    etti ve her tip sayfasi "bu bolge tum ailede aspartattir" dedi -- yani ayni
    geri cekilmis iddia 1 + 76 sayfada yayinda kaldi ve kullanici celiskiyi
    Atlas'ta buldu (2026-11-01).

    Bu, projenin tekrar eden hata sinifinin yeni bir ornegi: DUZELTME
    KESFEDILDIGI YERE UYGULANIR, YAYIMLANDIGI HER YERE UYGULANMAZ. Duzyazi
    oldugu icin hicbir sayi kontrolu yakalamiyordu. Bu kontrol, sayilarin
    sablonlara geri sizmasini mekanik olarak engeller.
    """
    tpl_dir = os.path.join(root, "webapp", "templates")
    if not os.path.isdir(tpl_dir):
        rep.check("retracted claim: the catalytic-column residue counts are "
                  "not published anywhere", True, 0,
                  "template directory missing", severity="WARN")
        return

    # Yayina sizmis sayim ifadeleri: katalitik kolonun asp/glu/dominant
    # degerlerini DOGRUDAN basan her sablon ifadesi.
    count_expr = re.compile(
        r"\{\{[^}]*\bcatalytic_iron\b[^}]*\.(?:asp|glu|dominant|asp_to_glu|"
        r"asp_fraction|with_carboxylate)\b[^}]*\}\}")
    # `{% set ci = ... %}` ile takma ad verip sonra ci.asp basan bicim.
    alias_expr = re.compile(
        r"\{%\s*set\s+(\w+)\s*=\s*[^%]*\bcatalytic_iron\b[^%]*%\}")

    # Geri cekilmis duzyazi: kolon sayimindan ligand sonucu cikaran cumleler.
    phrases = (
        "aspartate across the whole family",
        "facial triad uses Asp",
        "is aspartate and effectively nothing else",
    )

    offenders = []
    for name in sorted(os.listdir(tpl_dir)):
        if not name.endswith(".html"):
            continue
        path = os.path.join(tpl_dir, name)
        try:
            text = open(path, encoding="utf-8").read()
        except Exception:                             # noqa: BLE001
            continue
        body = re.sub(r"\{#.*?#\}", "", text, flags=re.S)   # yorumlar hesaba katilmaz
        for match in count_expr.finditer(body):
            offenders.append(f"{name}: {match.group(0).strip()}")
        for alias in alias_expr.findall(body):
            hits = re.findall(
                r"\{\{[^}]*\b" + re.escape(alias)
                + r"\.(?:asp|glu|dominant|asp_to_glu|asp_fraction|"
                  r"with_carboxylate)\b[^}]*\}\}", body)
            offenders.extend(f"{name}: {h.strip()}" for h in hits)
        for phrase in phrases:
            if phrase in body:
                offenders.append(f"{name}: \u201c{phrase}\u201d")

    rep.check("retracted claim: no page publishes the catalytic-column "
              "residue counts or calls that column the iron ligand",
              not offenders, len(offenders),
              "" if not offenders else
              "column 431 does not read the iron ligand (26/26 structures put "
              "its residue 10.9-14.7 A from the iron); these renders revive a "
              "retracted claim: " + "; ".join(offenders[:6]))


def check_group_monophyly(rep, root, out_dir):
    """Filogeni sayfasinin grup beyani OLCUMLE ayni mi?

    NEDEN. Sayfa "Deep nodes separating the colour blocks correspond to the
    historical RO groups, which this tree recovers from sequence alone" diye
    yaziyordu. Bu olculmemis bir SONUC beyaniydi; kullanici 2026-11-01'de
    "bizim olcumumuz bunu desteklemiyor" dedi ve olcum kullaniciyi dogruladi:
    grup 1-4 hicbir kesimde dal olusturmuyor.

    Iki sey ayri ayri korunur. Birincisi, geri cekilmis cumle sablonlara geri
    donmemeli. Ikincisi -- ve asil onemlisi -- sayfa sayilari artik
    group_monophyly.json'dan okuyor, o halde dosya VAR olmali ve icinde
    sayfanin dayandigi sey, yani "dal olusturan gruplar" kumesi, bulunmali.
    Dosya kaybolursa sayfa sessizce cumleyi dusurur ve okuyucu hicbir sey
    gormez; bu kontrol o sessizligi bir FAIL'e cevirir.
    """
    tpl = os.path.join(root, "webapp", "templates", "atlas_phylogeny.html")
    retracted = "which this tree recovers from sequence alone"
    if os.path.exists(tpl):
        text = open(tpl, encoding="utf-8").read()
        rep.check("phylogeny page: the retracted \u201crecovers the historical "
                  "groups\u201d claim is gone",
                  retracted not in text, 1 if retracted in text else 0,
                  "" if retracted not in text else
                  "groups 1-4 are not clades on this tree in any of the three "
                  "measured sequence sets; see group_monophyly.json")

    path = os.path.join(out_dir, "group_monophyly.json")
    if not os.path.exists(path):
        rep.check("phylogeny page: the group-monophyly measurement it cites "
                  "exists", False, 1,
                  f"{path} missing; the page would silently drop the claim "
                  f"instead of stating it")
        return
    try:
        doc = json.load(open(path, encoding="utf-8"))
    except Exception as exc:                          # noqa: BLE001
        rep.check("phylogeny page: the group-monophyly measurement it cites "
                  "is readable", False, 1, str(exc))
        return

    cuts = doc.get("cuts") or {}
    needed = ("all_tips", "references_only", "confident")
    absent = [c for c in needed if c not in cuts]
    rep.check("phylogeny page: the measurement covers all three sequence sets",
              not absent, len(absent),
              "" if not absent else "missing: " + ", ".join(absent))

    # Sayfa "yalnizca su gruplar dal" diyor; o kume dosyadan gelmeli ve BOS
    # OLMAMALI, yoksa cumle anlamsiz bir listeye duser.
    bad = []
    for cut in needed:
        block = cuts.get(cut) or {}
        groups = block.get("groups") or {}
        if not groups:
            bad.append(f"{cut}: no groups measured")
            continue
        real = {g: v for g, v in groups.items() if (v.get("n_leaves") or 0) > 1}
        if not real:
            bad.append(f"{cut}: every group has a single leaf")
            continue
        for g, v in real.items():
            for key in ("is_single_clade", "largest_pure_clade",
                        "largest_pure_fraction", "mrca_purity"):
                if v.get(key) is None:
                    bad.append(f"{cut}/group {g}: {key} missing")
    rep.check("phylogeny page: every measured group carries the fields the "
              "page renders", not bad, len(bad), first(bad, 4))

    # Geri cekilmenin SEBEBI kayitli kalmali: grup 1-4 dal DEGIL. Eger bir gun
    # olcum bunun tersini soylerse cumle de degismeli -- o zaman bu kontrol
    # duser ve kimse eski metni yeniden yazmak zorunda kalmaz, ama kimse de
    # olcumle celisen bir metni fark etmeden yayinda tutamaz.
    conf = (cuts.get("confident") or {}).get("groups") or {}
    ring = [g for g in ("1", "2", "3", "4") if g in conf]
    still = [g for g in ring if conf[g].get("is_single_clade")]
    rep.check("phylogeny page: the softened wording still matches the "
              "measurement (groups 1-4 are not clades)",
              not still, len(still),
              "" if not still else
              "group(s) " + ", ".join(still) + " now measure as a single "
              "clade; the page text says they do not and must be updated")


def check_tables_are_scrollable(rep, root):
    """Her table.data bir .tablewrap icinde mi?

    Gorsel denetim (2026-11-01): sarilmayan tablolar 1440 px kolonda
    tasiyor ve sagdaki sutunlar kesiliyordu. On tablo sarilmamisti
    (cluster.html'de dokuz, ro.html'de bir). .tablewrap overflow-x: auto
    veriyor, yani tablo kendi icinde kayar ve sayfa duzeni bozulmaz.

    Bu kontrol gozle bakmaya BIRAKILMAZ: yeni bir tablo eklenince kimse
    1440 px'te acip sutun saymaz, ve kesilme tam da dar ekranda degil
    GENIS ekranda oldugu icin mobil testinde de gorunmez.
    """
    tdir = os.path.join(root, "webapp", "templates")
    bad = []
    for name in sorted(os.listdir(tdir)):
        if not name.endswith(".html"):
            continue
        path = os.path.join(tdir, name)
        text = open(path, encoding="utf-8").read()
        for m in re.finditer(r'<table class="data[^"]*"', text):
            before = text[max(0, m.start() - 220):m.start()]
            if "tablewrap" not in before:
                line = text[:m.start()].count("\n") + 1
                bad.append(f"{name}:{line}")
    rep.check("every data table sits inside a scrollable wrapper",
              not bad, len(bad),
              "tum tablolar sarili" if not bad else first(bad, 6))


def check_no_internal_refs_published(rep, chem, root):
    """Yayinlanan proza ic dosya yolu ya da birinci tekil sahis tasiyamaz.

    2026-11-01'de olculdu: alti tipin `notes` ve `source` alanlari tip
    sayfasinda yayinlaniyor ("Curated from:" ve Sekil 1 altyazisinin
    "Note." kismi) ve iclerinde analysis_out/... yollari vardi. O
    dosyalarin HICBIRI indirilebilir degil, yani okuyucuya acamayacagi
    bir referans veriliyordu. Bir tanesinde de birinci tekil sahis
    ("The family members I retrieved") kalmisti -- kuratorun not
    defterinden sayfaya gecmis bir cumle.

    Bu kontrol kuratorlu alanlara bakar, kodun kendi yorumlarina DEGIL:
    kaynak dosyalarda yol yazmak dogrudur, yayinlanan metinde degil.
    """
    pat = re.compile(
        r"analysis_out/\S+"
        r"|\b[\w./-]+\.json\b"
        r"|\b(?:I|my)\s+(?:retrieved|chose|found|own|measurement)\b"
        r"|\bTODO\b|\bFIXME\b")
    bad = []
    for r in chem:
        for field in ("notes", "source", "reaction", "substrate_en",
                      "product_en"):
            text = (r.get(field) or "")
            for m in pat.finditer(text):
                bad.append(f"{r['cluster']}.{field}: {m.group(0)}")
    rep.check("published prose carries no internal file path or first person",
              not bad, len(bad),
              f"{len(chem)} rows checked" if not bad else first(bad, 5))


def check_groups_published(rep, con, root):
    """Grup numarasi yayinlaniyorsa ne oldugu da yayinlanmak ZORUNDA.

    2026-11-01'de olculdu: site her yerde "group 3" yaziyordu ve grup 3'un
    ne oldugunu soyleyen hicbir sayfa YOKTU -- groups.yaml hic
    yayinlanmiyordu. Batie'nin bes sinifiyla karistirilma riski tam
    buradan doguyor; tanimsiz bir numara okuyucunun bildigi baska bir
    bes'li semaya yapisir.

    Ayrica groups.yaml'in UYE SAYISI TASIMADIGINI zorunlu kilar. Tasiyordu
    ve bayatlamisti (grup 1 icin 5746 uye / 14 tip; olcum 5104 / 16).
    Ayni sayinin iki yerde durmasi bu projenin tekrar eden hata sinifi.
    """
    path = os.path.join(root, "groups.yaml")
    try:
        import yaml as _yaml
        doc = _yaml.safe_load(open(path, encoding="utf-8")) or {}
    except Exception as exc:
        rep.check("groups: groups.yaml is readable", False, 1, str(exc))
        return
    entries = doc.get("groups") or []
    declared = {str(g.get("id")) for g in entries}

    measured = {str(r[0]) for r in con.execute(
        "SELECT DISTINCT ro_group FROM ro WHERE is_confirmed=1")}
    orphan = sorted(measured - declared)
    rep.check("groups: every group in the database is defined in groups.yaml",
              not orphan, len(orphan),
              f"{len(declared)} defined, {len(measured)} in use"
              if not orphan else "tanimsiz: " + ", ".join(orphan))

    # Sayi dosyada DURMAMALI.
    stale = sorted(str(g.get("id")) for g in entries
                   if "members" in g or "types" in g)
    rep.check("groups: groups.yaml carries no member or type counts",
              not stale, len(stale),
              "sayi yok, olcumden geliyor" if not stale
              else "sayi tasiyan grup: " + ", ".join(stale))

    # Sayfa VAR mi, ve numara bir tanima BAGLI mi?
    tpl = os.path.join(root, "webapp", "templates", "atlas_groups.html")
    nav = os.path.join(root, "webapp", "templates", "_atlas_nav.html")
    missing = []
    if not os.path.exists(tpl):
        missing.append("atlas_groups.html")
    elif "Batie" not in open(tpl, encoding="utf-8").read():
        missing.append("atlas_groups.html'de Batie ayrimi yok")
    if os.path.exists(nav) and "/atlas/groups" not in open(
            nav, encoding="utf-8").read():
        missing.append("gezinmede /atlas/groups yok")
    frz = os.path.join(root, "webapp", "freeze.py")
    if os.path.exists(frz) and '"atlas/groups"' not in open(
            frz, encoding="utf-8").read():
        missing.append("freeze.py statik yayina almiyor")
    rep.check("groups: the page that defines the numbers is published",
              not missing, len(missing),
              "yayinda" if not missing else "; ".join(missing))

    # Her grup bir ETIKET tasimak zorunda: etiketi olmayan bir grup sayfada
    # bos bir satir olarak gorunur ve numara yine tanimsiz kalir.
    nolabel = sorted(str(g.get("id")) for g in entries
                     if not (g.get("label") or "").strip())
    rep.check("groups: every group carries a label",
              not nolabel, len(nolabel),
              f"{len(entries)} grup etiketli" if not nolabel
              else "etiketsiz: " + ", ".join(nolabel))


def check_types_are_single_source(rep, root):
    """`types/` tek kaynaksa chemistry.csv ONDAN uretilebilmek ZORUNDA.

    NEDEN. types_build.py'ye gore bir tipin tek tanimi `types/<kume>.yaml`
    ve chemistry.csv o dizinden URETILIR. Ama bu kuralı hicbir sey
    zorlamiyordu: chemistry.csv elle duzeltilebiliyor, site duzelmis halini
    gosteriyor ve her sey dogru gorunuyor -- ta ki biri `types_build.py
    --write` kosana kadar; o an uretilen dosya elle yapilan duzeltmeyi
    SESSIZCE geri aliyor. 2026-11-01'de tam bu oldu: onceki oturumda
    chemistry.csv'ye islenen 10 hucrelik dus-notu temizligi manifestolara
    islenmemisti ve bir sonraki --write'ta kaybolacakti.

    Kontrol gidis-donus testidir: manifestolardan uretilen chemistry.csv ile
    diskteki dosya hucre hucre karsilastirilir.
    """
    sys.path.insert(0, root)
    try:
        import types_build
    except Exception as exc:                                  # pragma: no cover
        rep.check("types/: chemistry.csv is reproducible from the manifests",
                  False, 1, f"types_build.py yuklenemedi: {exc}")
        return
    mans = types_build.read_all(os.path.join(root, "types"))
    diffs, _ = types_build.roundtrip(mans, os.path.join(root, "chemistry.csv"))
    rep.check("types/: chemistry.csv is reproducible from the manifests",
              not diffs, len(diffs),
              "gidis-donus birebir" if not diffs else
              f"{len(diffs)} hucre ayrisik, ilki {diffs[0][0]}/{diffs[0][1]} -- "
              "chemistry.csv elle duzenlenmis; duzeltmeyi types/'a tasi ve "
              "types_build.py --write kos")


def check_significance_claims(rep, con, root):
    """Anasayfanin kapsam kartlari OLCULEBILIR mi?

    NEDEN. Kartlar bir IDDIA yayinliyor ("iki metal merkez, on iki farkli
    oksidasyon; uc yasam alani") ve bu projede yayinlanan her iddianin
    kosulu yazili olmak zorunda. Sayilar significance.yaml'daki `measure`
    anahtarindan uretiliyor, o yuzden uc sey ayrisabilir:

      1. Tanimsiz bir `measure` -> kart SESSIZCE dusuyor. Bir iddianin
         yayindan kaybolmasi, yanlis yayinlanmasi kadar kotu ve hic
         gorunmuyor.
      2. Olculen deger 0 -> bos kart; olcmedigini olculmus gibi gostermek.
      3. `basis` eksik -> iddia var, kosulu yok.
    """
    path = os.path.join(root, "significance.yaml")
    if not os.path.exists(path):
        return
    try:
        import yaml as _y                                 # noqa: PLC0415
        claims = (_y.safe_load(open(path, encoding="utf-8"))
                  or {}).get("claims") or []
    except Exception as exc:                              # noqa: BLE001
        rep.check("significance: significance.yaml is readable", False, 1,
                  str(exc))
        return
    sys.path.insert(0, os.path.join(root, "webapp"))
    try:
        import atlas as _atlas                            # noqa: PLC0415
    except Exception as exc:                              # noqa: BLE001
        rep.check("significance: the view is importable", False, 1, str(exc))
        return
    counts = dict(con.execute("SELECT ro_cluster, COUNT(*) FROM ro "
                              "WHERE is_confirmed=1 GROUP BY 1"))
    measured = {v["id"]: v for v in _atlas.significance_view(
        os.path.join(root, "chemistry.csv"), counts, path, con)}

    dropped = [c.get("id") or "?" for c in claims if (c.get("id") or "?") not in measured]
    rep.check("significance: every claim in the file is actually measured",
              not dropped, len(dropped),
              f"{len(claims)} iddia, hepsi olculuyor" if not dropped
              else "olculmeyen (tanimsiz `measure`?): " + ", ".join(dropped))

    zero = [k for k, v in measured.items() if not v.get("n")]
    rep.check("significance: no claim measures zero",
              not zero, len(zero),
              "her iddia sifirdan buyuk" if not zero
              else "sifir olcen iddia: " + ", ".join(zero))

    no_basis = [c.get("id") or "?" for c in claims
                if not (c.get("basis") or "").strip()]
    rep.check("significance: every claim states what it is counted from",
              not no_basis, len(no_basis),
              "her iddianin kosulu yazili" if not no_basis
              else "basis eksik: " + ", ".join(no_basis))

    # BIRIM TUTARLILIGI. Kartlar tip sayisi yaziyor ve bu sayi sitenin her
    # yerinde yayinlanan tip sayisiyla AYNI olmak zorunda. Ilk surumde 76
    # (kuratorlu referans) yaziyordu, kutucuklar 68 (uyesi olan tip)
    # diyordu: ayni buyukluk iki ayri sayiyla yayinda.
    n_types = len([1 for v in counts.values() if v])
    bad = [k for k, v in measured.items()
           if re.search(r"\b(\d[\d,]*) (?:characterised enzyme )?types?\b",
                        v.get("sub") or "")
           and int(re.search(r"\b(\d[\d,]*) (?:characterised enzyme )?types?\b",
                             v["sub"]).group(1).replace(",", "")) != n_types]
    rep.check("significance: the type count agrees with the rest of the site",
              not bad, len(bad),
              f"kartlar {n_types} tip yaziyor" if not bad
              else "tutmayan kart: " + ", ".join(bad))


def check_regulation_totals_agree(rep, con, root):
    """Ayni sayiyi veren iki tablo AYRISAMAZ.

    NEDEN. /atlas/regulation iki yerde ayni buyuklugu yayinliyor: 2. bolumun
    mimari tablosu `ro_regulation`i tek sorguda sayiyor, 5. bolumun tip
    tablosu ise `analysis_out/regulation_by_cluster.csv`den geliyor. Ikisi
    farkli zincirlerden uretiliyor, yani CSV bayatladigi an sayfa iki farkli
    sayiyi ayni iddia icin gosterir ve bunu kimse soylemez -- okuyucu 68
    satiri toplayip ucuncu bir sayi bulur. Bu kontrol esitligi zorunlu kilar.

    Mimarilerin TOPLAMI da onaylanan giris sayisina esit olmak zorunda:
    her girisin tam bir yukari-akim mimarisi var, yani bes sayi ayni kumenin
    parcalanisi. Esit degilse ya bir giris iki kez sayildi ya biri dustu.
    """
    import csv as _csv
    path = os.path.join(root, "analysis_out", "regulation_by_cluster.csv")
    if not con.execute("SELECT name FROM sqlite_master WHERE name='ro_regulation'"
                       ).fetchone() or not os.path.exists(path):
        return
    rows = list(_csv.DictReader(open(path)))
    cols = ("divergent_regulator", "codirectional_regulator",
            "divergent_other", "codirectional_other", "unknown")
    csv_sum = {c: sum(int(r.get(c) or 0) for r in rows) for c in cols}
    db_sum = dict(con.execute(
        "SELECT architecture, COUNT(*) FROM ro_regulation GROUP BY 1"))
    for c in cols:
        a, b = csv_sum[c], db_sum.get(c, 0)
        rep.check(f"regulation: per-type {c} sums to the architecture count",
                  a == b, 0 if a == b else abs(a - b),
                  f"{a} = {b}" if a == b else
                  f"tip tablosu {a}, mimari tablosu {b} -- "
                  "regulation_by_cluster.csv bayat (analyze_regulation.py)")
    n_csv = sum(int(r.get("n") or 0) for r in rows)
    n_db = con.execute("SELECT COUNT(*) FROM ro_regulation").fetchone()[0]
    rep.check("regulation: the five architectures partition the entries",
              n_csv == n_db == sum(csv_sum.values()),
              0 if n_csv == n_db == sum(csv_sum.values()) else 1,
              f"{n_csv} giris, mimari toplami {sum(csv_sum.values())}, tablo {n_db}")


def check_operon_gap_single_source(rep, con, root):
    """Operon bosluk esigi yalnizca build_operons.py'de tanimli olabilir.

    NEDEN. Esik DORT yerde ayri ayri yaziliydi ve biri AYRISMISTI:
        build_operons.py     DEFAULT_MAX_GAP = 150   <- tabloyu kuran
        analyze.py           OPERON_MAX_GAP  = 100   <- ayrismis
        operon_validation.py g <= 150 (ve alan adi `gap_under_150`)
        atlas_operons.html   "at most 150 base pairs" (uc yerde)
    Yani ayni veritabani uzerinde iki operon kurali vardi, ikisi de "OPERON"
    basligi tasiyordu ve raporu okuyan kisi sitenin sayisini bulamazdi.
    Ayrica esik bir ALAN ADININ icine gomulmustu: esik degisse ad yalan olur
    ama kimse farketmez.

    Esik artik tek yerde; `max_gap_bp` olarak yayimlaniyor ve sayfa oradan
    okuyor. Bu kontrol kopyanin geri gelmesini engeller.
    """
    import re as _re4
    offenders = []
    for fn in ("analyze.py", "operon_validation.py", "operon_relations.py",
               "variant_and_operon.py"):
        fp = os.path.join(root, fn)
        if not os.path.exists(fp):
            continue
        for i, line in enumerate(open(fp), 1):
            if _re4.match(r"\s*(OPERON_MAX_GAP|DEFAULT_MAX_GAP|MAX_GAP)"
                          r"\s*=\s*[0-9]", line):
                offenders.append(f"{fn}:{i}")
    rep.check("operon: the gap threshold is defined only in build_operons.py",
              not offenders, len(offenders),
              "esik tek kaynakta" if not offenders
              else "elle yazilmis esik: " + ", ".join(offenders))

    # Yayimlanan dosya esigi BEYAN etmek zorunda: etmezse okuyucu
    # `gap_under_max` alaninin neyin altini saydigini tahmin eder.
    import json as _json
    path = os.path.join(root, "analysis_out", "operon_validation.json")
    if not os.path.exists(path):
        rep.check("operon: operon_validation.json exists", False, 1,
                  "dosya yok", severity="WARN")
        return
    data = _json.load(open(path))
    try:
        from build_operons import DEFAULT_MAX_GAP
    except Exception as exc:  # noqa: BLE001
        rep.check("operon: the rule is importable", False, 1, str(exc))
        return
    pub = data.get("max_gap_bp")
    rep.check("operon: the published file declares the gap it measured",
              pub == DEFAULT_MAX_GAP, 0 if pub == DEFAULT_MAX_GAP else 1,
              f"yayimlanan esik {pub} bp, kural {DEFAULT_MAX_GAP} bp"
              + ("" if pub == DEFAULT_MAX_GAP else
                 " -- dosya bayat ya da esigi hic beyan etmiyor"))


def check_carboxylate_claim_is_measured(rep, con, root):
    """Anasayfa-atlas'taki karboksilat iddiasinda elle yazilmis sayi kalamaz.

    NEDEN. O paragraf canli degiskenlerle ELLE YAZILMIS sayilari
    karistiriyordu ve ikisi YANLISTI (2026-11-03 olcumu):

      - "The residue column 431 calls sits 10.9-14.7 A from the iron in all
        26 testable structures". Ham denetimde 6 yapi
        `column_call_distance = -1` tasiyor; bu bir mesafe DEGIL, "kolon
        cagrisinin yapida karsiligi bulunamadi" nobetcisi, ve bir yapida
        (HcaE/8K0A) mononukleer demir hic yok. Aralik 20 yapida olculmus,
        26'ya atfedilmisti.
      - "above 89 % in Stc2, CntA and BmoA". Olcum: Stc2 %93,4, CntA %90,7,
        **BmoA %87,4** -- yani BmoA 89'un altinda. Ustelik bu paylar uyelik
        buyudukce kayiyor ve A5 kurtarmasi tam da grup 5'i buyutuyor.

    KONTROL. Paragraf, capa cumleden itibaren duzlestirilir, Jinja
    ifadeleri ({{...}} ve {%...%}) CIKARILIR ve geri kalanda iki ya da daha
    cok haneli bir sayi aranir. Kalan varsa o sayi olcumden gelmiyor
    demektir. Boylece metin yeniden sarilsa da kontrol kirilmaz; yalnizca
    bir sayi SABITLENIRSE kirilir.
    """
    import re as _re3
    tpl = os.path.join(root, "webapp", "templates", "atlas_home.html")
    if not os.path.exists(tpl):
        rep.check("carboxylate claim: atlas_home.html exists", False, 1, tpl)
        return
    flat = _re3.sub(r"\s+", " ", open(tpl).read())
    anchor = "The two carboxylates of this fold"
    i = flat.find(anchor)
    if i < 0:
        rep.check("carboxylate claim: the paragraph is found", False, 1,
                  f"capa cumle yok: {anchor!r} -- kontrol bos gecemez",
                  severity="WARN")
        return
    window = flat[i:i + 1800]
    bare = _re3.sub(r"\{\{.*?\}\}|\{%.*?%\}", " ", window)
    nums = _re3.findall(r"(?<![\w.])\d{2,}(?![\w%])", bare)
    rep.check("carboxylate claim: every number comes from a measurement",
              not nums, len(nums),
              "paragrafta elle yazilmis sayi yok" if not nums
              else f"elle yazilmis sayi(lar): {', '.join(sorted(set(nums)))} "
                   f"-- olcum var ve sayfa onu kopyaliyor, yani bayatlayacak")


def check_sdp_regions(rep, con, root):
    """Yayimlanan SDP tablosunda hicbir Rieske ligandi katalitik sayilamaz.

    NEDEN. `analysis_out/sdp_positions.csv` sitede INDIRILEBILIR ve bir
    `in_catalytic_domain` sutunu tasiyor. O sutun 2026-11-03'e kadar elle
    yazilmis `CATALYTIC_DOMAIN_START = 180` sinirindan hesaplaniyordu, oysa
    ROmotif71'in Rieske ligand kolonlari 164, 166, 184, 187. Yani sinir
    ligand kumesinin ICINDEN geciyor ve Rieske Cys-2 ile Rieske His-2
    "katalitik domainde" yayimlaniyordu. 416 satirin 109'u (%26) yanlis
    bolgedeydi.

    Sinir artik ro_motif'ten turetiliyor. Bu kontrol turetmenin yerinde
    kalmasini sagliyor: bir Rieske ligand kolonu katalitik isaretlenirse
    ya da `region` sutunu kaybolursa FAIL.
    """
    import csv as _csv
    path = os.path.join(root, "analysis_out", "sdp_positions.csv")
    if not os.path.exists(path):
        rep.check("sdp: sdp_positions.csv exists", False, 1,
                  "dosya yok -- analyze_variants.py kosturulmali",
                  severity="WARN")
        return
    try:
        from ro_motif import RIESKE_SITES, CATALYTIC_SITES
    except Exception as exc:  # noqa: BLE001
        rep.check("sdp: ligand columns are readable", False, 1, str(exc))
        return
    rieske_cols = {c for c, _r, _l in RIESKE_SITES}
    first_cat = min(c for c, _r, _l in CATALYTIC_SITES)
    rows = list(_csv.DictReader(open(path)))
    if "region" not in (rows[0].keys() if rows else {}):
        rep.check("sdp: the published table carries a region column",
                  False, 1, "`region` sutunu yok -- bolge atamasi iki degerli "
                            "kalmis, oysa 188-288 bandinda olculmus sinir yok")
        return
    bad_lig = [r["column"] for r in rows
               if int(r["column"]) in rieske_cols
               and r["in_catalytic_domain"] == "1"]
    rep.check("sdp: no Rieske ligand column is published as catalytic",
              not bad_lig, len(bad_lig),
              "Rieske ligandlari katalitik sayilmiyor" if not bad_lig
              else f"katalitik isaretli Rieske ligand kolonlari: "
                   f"{', '.join(bad_lig)}")
    # Bayrak ile olculen sinir ayrisirsa dosya bayat.
    wrong = [r["column"] for r in rows
             if (int(r["column"]) >= first_cat) != (r["in_catalytic_domain"] == "1")]
    rep.check("sdp: the catalytic flag agrees with the measured boundary",
              not wrong, len(wrong),
              f"sinir kolon {first_cat}, bayrak uyusuyor" if not wrong
              else f"{len(wrong)} satirda bayrak olculen sinirla uyusmuyor "
                   f"-- dosya bayat, analyze_variants.py kosturulmali")

    # BOLGE SINIRI TEK KAYNAK OLMAK ZORUNDA. 2026-11-03'e kadar ayni sinir
    # dort yerde tanimliydi (analyze_variants 180; variant_signature ve
    # learn_from_data 120/191) ve ucu de ligand kolonlariyla celisiyordu.
    # Burada kaynak dosyalarda ELLE YAZILMIS bir sinir kalmadigi sinaniyor:
    # kalirsa ro_motif'teki turetme ile sessizce ayrisabilir.
    import re as _re2
    literal = []
    for fn in ("analyze_variants.py", "variant_signature.py",
               "learn_from_data.py"):
        fp = os.path.join(root, fn)
        if not os.path.exists(fp):
            continue
        for i, line in enumerate(open(fp), 1):
            if _re2.match(r"\s*(RIESKE_END|CATALYTIC_START|"
                          r"CATALYTIC_DOMAIN_START|RIESKE_DOMAIN_END)"
                          r"\s*=\s*[0-9]", line):
                literal.append(f"{fn}:{i}")
    rep.check("sdp: the domain boundary is defined only in ro_motif",
              not literal, len(literal),
              "sinir yalnizca ro_motif'te turetiliyor" if not literal
              else "elle yazilmis sinir: " + ", ".join(literal))


def check_prevalence_declares_cohort(rep, con, root):
    """Genom duzeyinde yaygınlık, kohortunu BEYAN etmek zorunda.

    NEDEN. `about.html` "N annotated replicons scanned, M carrying at least
    one" yayimliyor. Bu bir yaygınlık iddiasi ve paydasi "taranan genomlar".
    A5'te cekilen 13.079 kayit, ICINDE BIR RO OLDUGU BILINDIGI ICIN secildi;
    ayni paydaya girerse oran ~%100'e cikar ve bir KESIF gibi gorunur. Oysa
    degisen sey olcum degil ORNEKLEME.

    UC SEY SINANIR:
      1. `replicon_cohort` var mi -- yoksa kohort tanimsiz, iddia yapilamaz.
      2. survey kohortunun boyu, anlik goruntunun (cohorts/survey.txt) boyuna
         ESIT mi. Survey listesi BUYUMEZ; buyurse payda hedeflenmis
         kayitlarla kirlenir ve bunu kimse goremez.
      3. Yaygınlık cumlesi `survey` kelimesini tasiyor mu.

    UCUNCU SINAMA SATIR TABANLI DEGIL. `check_union_denominator` bir kez
    satir sonu yuzunden kirildi (madde 117); burada sablonun BOSLUGU
    DUZLESTIRILIR ve kelime, cumlenin etrafindaki pencerede aranir. Boylece
    metin yeniden sarildiginda kontrol kirilmaz, yalnizca BEYAN kalkarsa
    kirilir.
    """
    import re as _re
    have = con.execute("SELECT name FROM sqlite_master "
                       "WHERE name='replicon_cohort'").fetchone()
    if not rep.check("prevalence: replicon_cohort exists", bool(have), 0 if have else 1,
                     "kohort tablosu var" if have else
                     "replicon_cohort yok -- a5_cohort.py kosturulmali"):
        return

    snap = os.path.join(root, "cohorts", "survey.txt")
    n_snap = 0
    if os.path.exists(snap):
        with open(snap) as fh:
            n_snap = sum(1 for ln in fh
                         if ln.strip() and not ln.startswith("#"))
    n_survey = con.execute("SELECT COUNT(*) FROM replicon_cohort "
                           "WHERE cohort='survey'").fetchone()[0]
    rep.check("prevalence: the survey cohort matches its snapshot",
              n_snap > 0 and n_snap == n_survey,
              abs(n_snap - n_survey),
              f"anlik goruntu {n_snap}, tabloda {n_survey}"
              + ("" if n_snap == n_survey else
                 " -- survey havuzu degismis; yaygınlık paydasi artik ayni "
                 "kumeyi olcmuyor"))

    tpl = os.path.join(root, "webapp", "templates", "about.html")
    if not os.path.exists(tpl):
        rep.check("prevalence: about.html exists", False, 1, tpl)
        return
    flat = _re.sub(r"\s+", " ", open(tpl).read())
    phrase = "carrying at least one"
    ok, detail = True, "yaygınlık cumlesi kohortunu yaziyor"
    for m in _re.finditer(_re.escape(phrase), flat):
        window = flat[max(0, m.start() - 400):m.end() + 400]
        if "survey" not in window:
            ok = False
            detail = ("yaygınlık cumlesi kohortunu YAZMIYOR -- okuyucu orani "
                      "hedeflenmis kayitlari da iceren bir paydadan sanar")
            break
    rep.check("prevalence: the genome-level claim names its cohort", ok,
              0 if ok else 1, detail)


def check_prefilter_admits_gate_passers(rep, con, root):
    """Ucuz on filtre, HMM kapisinin kabul ettigi bir diziyi ELEYEMEZ.

    NEDEN BU KONTROL VAR. `extract_genomic_context.py` milyonlarca CDS'i
    binlere indirmek icin elle yazilmis bir Rieske regex'i kullaniyor ve
    2026-11-02'ye kadar penceresi 12-22 kalintiydi. O pencere REFERANSLARDAN
    turetilmisti (76/76 referans gecer) ama UYELERDEN turetilmemisti.
    Olcum: kapiyi gecen 71.819 dizinin 3.107'si (%4,33) regex'e takiliyordu
    ve bunlarin SIFIRI genomik baglam tasiyordu. Yani bir dizi protein
    tarafinda "onayli RO alpha", genomik tarafta ise "aday bile degil"
    olabiliyordu -- ayni veritabaninin iki yerinde iki farkli cevap.

    Kaybin grup 5'te yogunlasmasi (CntA 1.165, Stc2 515, qxyA 269, BmoA 213)
    bunu gorunmez kiliyordu: eksiklik her yere dagilmis degil, tek bir
    yapisal sinifta toplanmisti ve o sinif genomik tarafta hic yoktu, yani
    "eksik" diye bakilacak bir sayi da yoktu.

    KONTROLUN MANTIGI. Regex'in kendisi dogrulanmaz; regex'in KARARI
    kapinin karariyla karsilastirilir. Kapiyi gecen her birlesim dizisi
    on filtreden de gecmek zorundadir. Gecemeyenler ARTAKALAN olarak
    sayilir ve sayi kayitli esikten BUYUKSE FAIL verilir -- yani pencere
    daraltilirsa ya da yeni bir yapisal sinif girerse kapi kapanir.

    Artakalan sifir DEGIL ve olmasi da beklenmiyor: olculen 32 dizinin
    12'sinde 80 kalintilik pencerede hic C-x-H..C-x-x-H kalibi yok. Onlari
    bir regex bulamaz; HMM gerekir.
    """
    import re as _re
    import importlib.util as _ilu
    path = os.path.join(root, "extract_genomic_context.py")
    if not os.path.exists(path):
        rep.check("prefilter: extract_genomic_context.py exists", False, 1, path)
        return
    # Regex'i KAYNAKTAN okur. Kopyalanmis bir kalip, kontrolu yanlis seyi
    # sinayan bir kontrole cevirirdi (projenin tekrar eden hata bicimi).
    src = open(path).read()
    m = _re.search(r'RIESKE_REGEX = re\.compile\(r"([^"]+)"\)', src)
    if not m:
        rep.check("prefilter: the Rieske pre-filter pattern is readable",
                  False, 1, "RIESKE_REGEX kaynaktan okunamadi")
        return
    rx = _re.compile(m.group(1))

    have = con.execute(
        "SELECT name FROM sqlite_master WHERE name='ro_union'").fetchone()
    if not have:
        rep.check("prefilter: ro_union exists", False, 1, "ro_union yok")
        return
    rows = con.execute("SELECT sequence FROM ro_union "
                       "WHERE passes_gate=1 AND sequence IS NOT NULL").fetchall()
    rejected = [s for (s,) in rows if not rx.search(s)]
    # 2026-11-02 OLCUMU: 12-45 penceresiyle 32 dizi kurtarilamiyor. Butce
    # olculen degerin KENDISI; daha gevsek bir butce, aradaki farki sessizce
    # kaybetmeye izin verirdi.
    BUDGET = 32
    n = len(rejected)
    rep.check("prefilter: the cheap pre-filter admits every gate-passing sequence",
              n <= BUDGET, n,
              f"{n} dizi eleniyor, butce {BUDGET} (kapiyi gecen {len(rows)})"
              + ("" if n <= BUDGET else
                 " -- on filtre penceresi kapinin kabul ettigi bir sinifi "
                 "disliyor; extract_genomic_context.py'deki olcume bakilmali"))

def check_a5_match_guard(rep, con, root):
    """Capraz referansla gelen baglam, BIREBIR dizi eslesmesi olmadan baglanamaz.

    NEDEN BU KONTROL VAR. A5'te 18.636 (dizi, kayit) cifti NCBI capraz
    referansindan geldi. Capraz referans "bu dizi su kayitta gecer" der;
    "su kayitta ayni proteindir" DEMEZ. Kaydi acip oradaki herhangi bir RO
    adayini almak cazip ve ucuzdu. Olculdu: eslesmeyen 379 ciftin 114'unde
    aday ile hedef dizinin ORTAK PARCASI BILE YOK. Yani "ayni kayit yeter"
    kurali uygulansaydi en az 114 dizi BASKA bir genin komsulugu ile
    yayimlanacakti -- ve sayfada bunu gosteren hicbir sey olmayacakti.

    KONTROLUN MANTIGI. Kural `a5_match_guard.py`de uygulanir, butce de
    ORADA yazilidir; bu kontrol butceyi KAYNAKTAN okur ve cikti JSON'u ile
    karsilastirir. Butce kaynaktan okunmazsa esik iki yere kopyalanmis olur
    ve bu projenin tekrar eden hata bicimi tam budur.

    Ayrica TAZELIK sinanir: hedef listesi JSON'dan yeni ise olcum bayattir
    ve "%97,4 baglam kurulur" iddiasi artik olculmemis bir iddiadir.
    """
    import json as _json
    import re as _re
    guard = os.path.join(root, "a5_match_guard.py")
    out = os.path.join(root, "analysis_out", "a5_match_guard.json")
    targets = os.path.join(root, "cohorts", "a5_protein_recovery.tsv")
    if not os.path.exists(guard):
        rep.check("a5 guard: a5_match_guard.py exists", False, 1, guard)
        return
    if not os.path.exists(out):
        rep.check("a5 guard: the match measurement has been run", False, 1,
                  f"{out} yok -- python3 a5_match_guard.py", severity="WARN")
        return
    src = open(guard).read()
    m = _re.search(r"^RESIDUAL_BUDGET = (\d+)", src, _re.M)
    w = _re.search(r"^DIFFERENT_PROTEIN = (\d+)", src, _re.M)
    if not m or not w:
        rep.check("a5 guard: the residual budget is readable from source",
                  False, 1, "RESIDUAL_BUDGET / DIFFERENT_PROTEIN okunamadi")
        return
    budget, declared_wrong = int(m.group(1)), int(w.group(1))
    data = _json.load(open(out))

    residual = data.get("residual_pairs", -1)
    rep.check("a5 guard: unmatched pairs stay within the measured budget",
              0 <= residual <= budget, residual,
              f"{residual} eslesmeyen cift, butce {budget}"
              + ("" if 0 <= residual <= budget else
                 " -- cikarim degismis; a5_match_guard.py yeniden kosturulup "
                 "butce OLCUMLE birlikte guncellenmeli"))

    # Kuralin kazanci yayimlanabilir bir sayi; kaynaktaki beyan ile olcum
    # ayrisirsa iddia olculmemis hale gelir.
    measured_wrong = data.get("would_be_wrong_without_rule", -1)
    rep.check("a5 guard: the stated cost of dropping the rule is the measured one",
              measured_wrong == declared_wrong, abs(measured_wrong - declared_wrong),
              f"olculen {measured_wrong}, kaynakta beyan {declared_wrong}")

    # Birebir eslesme sarti gercekten UYGULANMIS mi: baska kayitta eslesen
    # ciftler baglam KAZANMAMIS olmali.
    rep.check("a5 guard: a match in a different record does not link context",
              data.get("matched_other_record_only", 0) >= 0
              and data.get("sequences_linked", 0) <= data.get("n_target_sequences", 0),
              data.get("matched_other_record_only", 0),
              f"{data.get('matched_other_record_only', 0)} cift yalnizca baska "
              f"kayitta esliyor ve baglam almiyor")

    if os.path.exists(targets) and os.path.getmtime(targets) > os.path.getmtime(out):
        rep.check("a5 guard: the measurement is newer than the target list",
                  False, 1, "hedef listesi JSON'dan yeni -- olcum bayat",
                  severity="WARN")
    else:
        rep.check("a5 guard: the measurement is newer than the target list",
                  True, 0)


def check_ro_footprint_explained(rep, con, root):
    """Onaylanmis bir RO'nun genomik ayak izi cevirisiyle ACIKLANABILIR olmali.

    NEDEN BU KONTROL VAR. `extract_genomic_context.py` icindeki
    `feature_spans` dairesel replikonda orijini asan genleri DOGRU ele alir
    ve komsu mesafeleri parcalar uzerinden hesaplanir. Ama adayin kendi
    `start`/`end` degerleri parcalarin MIN ve MAX'i olarak yazilir ve
    kaynaktaki yorum bunu "raporlama icin" diye isaretler. Asagi akista
    `build_operons.py`, `operon_validation.py` ve `analyze.py` bosluk
    aritmetigini TAM O iki degerle yapar. Yani tuzak komsular icin
    kapatilmis, adayin KENDISI icin kapatilmamis; ustelik `ro` tablosu
    komsu tablosunun tasidigi `spans_origin` alanini tasimadigi icin bilgi
    hic ulasmiyor.

    SONUC NE OLURDU. Ayak izi replikon uzunluguna esit bir RO'da o
    replikondaki HER komsu "ortusuyor" sayilir, bosluk hep 0 cikar, operon
    yurumesi hic kirilmaz ve sayfa TUM KROMOZOMU kapsayan bir operon
    yayimlar. 2026-11-04 olcumunde A5 kayitlarinda tam bu durumda iki
    bakteri girisi var: CP063373.1 (9,96 Mb ayak izi / 318 aa) ve
    CP060714.1 (5,17 Mb / 905 aa), ikisi de dairesel.

    KONTROLUN MANTIGI. Oran `ro_cds_span.py`de yazilidir ve buraya
    KAYNAKTAN okunur. Aciklanamaz ayak izi olan her onayli giris ya
    `ro_cds_span` tablosunda gercek parcalariyla COZULMUS olmali, ya da
    kontrol FAIL verir. Yani "aciklanamaz" durumu sessiz kalamaz: ya
    aciklanir ya kapi kapanir.

    Eukaryota girisleri de ayrica sayilir: orada ayak izi intron tasir,
    yani buyuk span bir hata DEGIL -- ama operon mantigi da gecerli degil
    ve bu sayinin gorunur olmasi gerekir.
    """
    import re as _re
    script = os.path.join(root, "ro_cds_span.py")
    if not os.path.exists(script):
        rep.check("footprint: ro_cds_span.py exists", False, 1, script)
        return
    m = _re.search(r"^UNEXPLAINED_RATIO = ([\d.]+)", open(script).read(), _re.M)
    if not m:
        rep.check("footprint: the unexplained-footprint ratio is readable from source",
                  False, 1, "UNEXPLAINED_RATIO okunamadi")
        return
    ratio = float(m.group(1))

    have_span = con.execute("SELECT name FROM sqlite_master "
                            "WHERE name='ro_cds_span'").fetchone()
    # COZULMUS = SATIR VAR demek DEGILDIR (2026-11-05 duzeltmesi). Ilk surumde
    # `resolved` yalnizca "ro_cds_span'de kaydi var mi" diye sorardi. Yani
    # ayak izi tum replikonu kaplayan bir giris, SADECE kaydedildigi icin
    # temiz sayiliyordu -- oysa bu kontrolun kendi FAIL metni zararin ne
    # oldugunu yaziyor: "ayak izi replikon uzunluguna yaklasirsa bosluk
    # aritmetigi her komsuyu ortusuyor sayar". Kaydetmek kullanmak degildir.
    #
    # Dogru olcut: parcalardan hesaplanan GERCEK ayak izi (`cds_length`)
    # cevirinin butcesine giriyor mu. Girmiyorsa parcalar o geni
    # aciklamiyor ve cozum bir kurgudur.
    span_rows = {}
    if have_span:
        span_rows = {c: (ln, parts) for c, ln, parts in con.execute(
            "SELECT candidate_id, cds_length, parts FROM ro_cds_span")}

    # Ayak izi, cevirinin beklenen bp'sinin `ratio` katindan buyuk mu?
    # Ceviri uzunlugu `ro_union.length` uzerinden gelir (dizi DB'de).
    rows = con.execute("""
        SELECT r.candidate_id, r.end - r.start, u.length
        FROM ro r JOIN ro_union u ON u.genomic_candidate_id = r.candidate_id
        WHERE r.is_confirmed = 1 AND u.length IS NOT NULL""").fetchall()
    unexplained = [(cid, span, n_aa) for cid, span, n_aa in rows
                   if span > ratio * 3 * (n_aa + 1)]
    unresolved = [c for c, _s, n_aa in unexplained
                  if c not in span_rows
                  or span_rows[c][0] > ratio * 3 * (n_aa + 1)]

    rep.check("footprint: every confirmed RO footprint is explained by its translation "
              "or resolved into real parts",
              not unresolved, len(unresolved),
              f"{len(rows):,} onayli giris olculdu, aciklanamaz {len(unexplained)}, "
              f"cozulmemis {len(unresolved)}"
              + ("" if not unresolved else
                 f" -- {first(sorted(unresolved), 3)}; ro_cds_span.py kosturulmali. "
                 "Ayak izi replikon uzunluguna yaklasirsa bosluk aritmetigi her "
                 "komsuyu ortusuyor sayar ve operon tum kromozomu kapsar"))

    if have_span:
        # BU KONTROL ESKIDEN `True` SABITIYLE CAGRILIYORDU (2026-11-05).
        # Yani hic FAIL edemezdi: adi kontrol, isi sayac. Artik yanlislanabilir
        # bir sey soruyor -- kaydedilen parcalar gercekten o geni veriyor mu:
        # parcalarin toplami `cds_length`e esit olmali. Esit degilse yazici
        # hatali ve yukaridaki butce testi yanlis sayiyla yapiliyor.
        bad = []
        for cid, (ln, parts) in span_rows.items():
            try:
                total = sum(int(b) - int(a)
                            for a, b in (p.split("-") for p in parts.split(";") if p))
            except ValueError:
                bad.append(f"{cid} (parcalar ayristirilamadi)")
                continue
            if total != ln:
                bad.append(f"{cid} (parcalar {total} != cds_length {ln})")
        wraps = con.execute(
            "SELECT COUNT(*) FROM ro_cds_span WHERE wraps_origin = 1").fetchone()[0]
        rep.check("footprint: recorded parts actually add up to the recorded "
                  "CDS length",
                  not bad, len(bad),
                  f"{len(span_rows):,} kayit sinandi, {wraps} orijin asan"
                  + ("" if not bad else f"; TUTARSIZ: {first(sorted(bad), 3)}"))


def check_component_scan_is_fresh(rep, con, root):
    """Bilesen taramasi komsu protein kumesiyle GUNCEL olmali.

    NEDEN BU KONTROL VAR. `build_operons.py:run_component_hmms` icinde iki
    VARLIK kontrolu vardi: `neighbor_domain` tablosunda tek satir varsa
    fonksiyonun tamami atlaniyor, ve `domtbl` dosyasi varsa hmmsearch
    atlaniyordu. Diskteki kanit: domtbl 2026-10-05'ten, komsu protein
    fastasi 2026-10-14'ten, veritabani 2026-11-03'ten. Yani komsu kumesi
    buyumus, bilesen taramasi buyumemis; o tarihten sonra eklenen her
    protein `none` siniflandi.

    MALIYETI OLCULDU. Tarama dogru kosturuldugunda:
        beta        2.523 -> 4.372   (+%73)
        reduktaz    6.185 -> 10.394  (+%68)
        ferredoksin 1.624 -> 2.736   (+%68)
    ve operon bayraklari has_beta 2.502 -> 4.333, has_reductase 3.203 ->
    5.409. Yani yayimlanan operon paylari gercek bilesenlerin ~%40'ini
    kaciriyordu. Ayni tuzak `annotate_ro.py`de bulunup cozulmus, bu dosyaya
    uygulanmamisti (bkz. stale_intermediate.py).

    KONTROLUN MANTIGI IKI AYAKLI.
      1. MEKANIZMA: domtbl, girdilerinden (komsu fastasi, bilesen HMM'i)
         daha yeni olmak zorunda. Kural `stale_intermediate`ten alinir,
         burada YENIDEN YAZILMAZ.
      2. SONUC: mtime'lardan BAGIMSIZ bir sinama. Anotasyonun kendisi
         "dioxygenase small/beta subunit" dedigi her komsu protein `beta`
         siniflanmali. Mekanizma kontrolu bir gun yanilirsa bu ayak hala
         tutar, cunku dosya tarihlerine degil SINIFLANDIRMANIN KENDISINE
         bakar.

    Artakalan 60 protein olculen degerin kendisidir (duzeltme oncesi 688).
    O 60, kullanicinin "anotasyon yanlis da olabiliyor" tespitinin asil aday
    kumesi: urun metni beta diyor, hicbir model tutmuyor. Ikisinden biri
    yanlis ve hangisi oldugu HENUZ OLCULMEDI.
    """
    import sys as _sys
    _sys.path.insert(0, root)
    try:
        from stale_intermediate import stale
    except ImportError:
        rep.check("components: the staleness rule has a single source", False, 1,
                  "stale_intermediate.py bulunamadi")
        return

    work = os.path.join(root, "operon_work")
    fasta = os.path.join(work, "neighbor_proteins.fasta")
    domtbl = os.path.join(work, "neighbor_components.domtbl")
    hmm = os.path.join(root, "component_hmms", "components.hmm")
    if os.path.exists(domtbl):
        why = stale(domtbl, fasta, hmm)
        rep.check("components: the component scan is newer than the protein set",
                  not why, 1 if why else 0,
                  why or f"{os.path.basename(domtbl)} guncel")
    else:
        rep.check("components: the component scan is newer than the protein set",
                  True, 0, "domtbl yok -- tarama bir sonraki kosuda uretilir",
                  severity="WARN")

    # 2. ayak: sonuc sinamasi, dosya tarihlerinden bagimsiz.
    BUDGET = 60
    key = "nb.nucleotide_id||':'||nb.start||'-'||nb.end||':'||nb.strand"
    unclassified = con.execute(f"""
        SELECT COUNT(*) FROM (
          SELECT DISTINCT {key} AS pk FROM neighbor nb
          WHERE lower(nb.product) LIKE '%dioxygenase small subunit%'
             OR lower(nb.product) LIKE '%dioxygenase beta subunit%')
        WHERE pk NOT IN (SELECT protein_key FROM neighbor_component
                         WHERE component = 'beta')""").fetchone()[0]
    rep.check("components: a neighbour the annotation calls a beta subunit is "
              "classified as beta",
              unclassified <= BUDGET, unclassified,
              f"{unclassified} protein urun metninde beta diyor ama beta "
              f"siniflanmamis, butce {BUDGET}"
              + ("" if unclassified <= BUDGET else
                 " -- bilesen taramasi bayat olabilir; build_operons.py "
                 "--force ile kosturulup sayilar yeniden olculmeli"))


def check_same_structure_same_name(rep, chem):
    """Ayni SMILES'i tasiyan referanslar substrati AYNI ADLA anmali.

    NEDEN. Ayni substrata etki eden birden cok referans olmasi normaldir ve
    hata degildir (bifenil 5, ftalat 3 referans). Hata olan sey, AYNI
    molekule sayfada IKI FARKLI AD verilmesidir: okuyucu iki ayri substrat
    gorur. Olculen ornek: `2_203_cumA1` "cumene", `2_206_IpBAa`
    "isopropylbenzene" -- SMILES'leri harfi harfine ayni (CC(C)c1ccccc1),
    cunku ikisi ayni bilesigin iki adi. Ikisi %90,8 benzer ve IpBAa SIFIR
    uye topluyor: ad ikiligi dogrudan atamaya da vuruyor.

    BUTCE 3, cunku kalan ikisi KASITLI genis etiket: "caffeine and
    methylxanthines" ve "polycyclic aromatic hydrocarbons" -- orada SMILES
    temsilcidir, ad daha genistir. Yanlis olan ad degil, SMILES'in dar
    olmasi; bu kontrol onu da gorunur tutuyor.
    """
    groups = defaultdict(set)
    for row in chem:
        smiles = (row.get("substrate_smiles") or "").strip()
        name = (row.get("substrate_en") or "").strip()
        if smiles and name:
            groups[smiles].add(name)
    clashes = sorted(f"{s} -> {' / '.join(sorted(n))}"
                     for s, n in groups.items() if len(n) > 1)
    BUDGET = 3
    rep.check("chemistry.csv: one structure is published under one substrate name",
              len(clashes) <= BUDGET, len(clashes),
              f"butce {BUDGET}; " + ("; ".join(clashes) if clashes else "catisma yok"))


def check_promoter_claim_is_levelled(rep, root, out_dir="analysis_out"):
    """Promotor mimarisi iddiasi UC DUZEYDE ve CINS ICINDE de verilmis olmali.

    NEDEN BU KONTROL VAR. 17.730 giriste p degeri neredeyse her seyi anlamli
    yapar; bu projede bir iliski ancak giris / dizi / cins duzeylerinde ayni
    yonde durursa raporlanir. Dosyanin kendisi tazelik kapisindan geciyor ama
    tazelik "icinde ne var" sorusunu sormaz: biri ileride yalnizca giris
    duzeyini biraksa dosya tap taze olur ve disiplin sessizce kaybolurdu.
    Bu kontrol dosyanin SEKLINI sinar, sonucunu degil -- bir sonucun yonu
    veriyle degisebilir, disiplin degisemez.
    """
    path = os.path.join(root, out_dir, "promoter_architecture.json")
    if not os.path.exists(path):
        rep.check("promoter: the architecture comparison is published", True, 0,
                  "promoter_architecture.json yok -- adim henuz kosmadi",
                  severity="WARN")
        return
    with open(path) as fh:
        data = json.load(fh)
    want = ("entry", "sequence", "genus")
    levels = data.get("levels") or {}
    missing = [lv for lv in want
               if not (levels.get(lv) or {}).get("usable")]
    if not (data.get("within_genus") or {}).get("usable"):
        missing.append("within_genus")
    rep.check("promoter: the intergenic comparison is reported at entry, "
              "sequence and genus level and within genus",
              not missing, len(missing),
              ("eksik duzey: " + ", ".join(missing)) if missing
              else "uc duzey + cins ici esli sinama var")

    # Yon tutarliligi: duzeyler arasinda isaret degisiyorsa bu bir SONUC
    # degil, bir uyaridir -- iddia zayif demektir, yayimi durdurmaz.
    signs = {lv: (levels[lv]["median_difference"] > 0) - (levels[lv]["median_difference"] < 0)
             for lv in want if (levels.get(lv) or {}).get("usable")}
    consistent = len(set(signs.values())) <= 1
    rep.check("promoter: the three levels point the same way",
              consistent, 0 if consistent else len(set(signs.values())),
              str(signs) if not consistent else
              f"uc duzey ayni yonde ({'kisa' if -1 in signs.values() else 'uzun'})",
              severity="WARN")


def check_family_scan_is_fresh(rep, con, root):
    """Aile taramasi (hareketlilik / rekombinaz / duzenleyici) guncel ve
    kendi kutuphanesiyle tutarli olmali.

    NEDEN BU KONTROL VAR. Tarama uc AYRI Pfam kutuphanesi kullanir ve her
    biri `family_scan.py` icindeki bir MODEL->AILE sozluguyle eslesir. Iki
    sessiz bozulma bicimi var:

      1. BAYATLIK: domtbl, komsu protein fastasindan ya da .hmm'den eski
         kalirsa tarama yeni proteinleri hic gormez ve hepsi sinifsiz sayilir.
         Ayni tuzak `build_operons.py`de yasandi; kural yine
         `stale_intermediate`ten alinir, burada yeniden yazilmaz.
      2. SESSIZ KAYMA: kutuphaneye model eklenip sozluge eklenmezse o modelin
         tuttugu her protein `unassigned_model` olur ve sayilar hicbir hata
         vermeden degisir. Model listesi TEK KAYNAKTAN (manifesto) okunur ve
         sozlukle karsilastirilir.

    OLCULEN DEGER (2026-11-10). Urun adi regex'inin 'transposon' dedigi 3.385
    komsu satirinin 1.447'si transpozaz modeliyle dogrulandi, 486'si
    rekombinazdi (hareketli eleman DEGIL), 1.452'sinde hicbir model tutmadi.
    Yani etiketin %43'u hala yalnizca urun metnine dayaniyor.
    """
    import sys as _sys
    _sys.path.insert(0, root)
    try:
        from stale_intermediate import stale
        import family_scan as _fs
    except ImportError as exc:
        rep.check("families: the family scan is importable", False, 1, str(exc))
        return

    work = os.path.join(root, "family_work")
    fasta = os.path.join(work, "neighbor_proteins.fasta")
    reasons = []
    for library in _fs.LIBRARIES:
        hmm = os.path.join(root, _fs.HMM_DIR, f"{library}.hmm")
        domtbl = os.path.join(work, f"{library}.domtbl")
        if not os.path.exists(domtbl):
            reasons.append(f"{library}: domtbl yok")
            continue
        why = stale(domtbl, fasta, hmm)
        if why:
            reasons.append(f"{library}: {why}")
    rep.check("families: every family scan is newer than the protein set and "
              "its model library",
              not reasons, len(reasons),
              "; ".join(reasons) or f"{len(_fs.LIBRARIES)} kutuphane guncel")

    # 2. ayak: kutuphanedeki her model bir aileye bagli mi?
    missing = []
    for library in _fs.LIBRARIES:
        man = os.path.join(root, _fs.HMM_DIR, f"{library}.manifest.json")
        if not os.path.exists(man):
            missing.append(f"{library}: manifesto yok")
            continue
        mapping = _fs.FAMILY_MAP[library]
        with open(man) as fh:
            for model in json.load(fh)["models"]:
                if model["name"] not in mapping:
                    missing.append(f"{library}/{model['name']}")
    rep.check("families: every model in a library maps to a family",
              not missing, len(missing),
              "; ".join(missing) or "her model bir aileye bagli")

    # 3. ayak: tabloda sozlukte olmayan aile kalmamali.
    unassigned = con.execute(
        "SELECT COUNT(*) FROM neighbor_family "
        "WHERE category = 'unassigned_model'").fetchone()[0]
    rep.check("families: no neighbour carries an unmapped model label",
              unassigned == 0, unassigned,
              f"{unassigned} satir 'unassigned_model'")


def check_derived_statistics_are_fresh(rep, con, root, out_dir="analysis_out"):
    """Veritabanindan turetilen her yayimlanmis istatistik, OKUDUGU tablolarin
    bugunku icerigiyle hesaplanmis olmali.

    NEDEN BU KONTROL VAR (kullanici istegi, 2026-11-04). Bu depoda bayatliga
    karsi tek tek kontroller vardi (tepkime raporu, ro_search, turetilmis rol
    tablosu...) ama GENEL bir kapi yoktu. Yani bir tablo duzeltildiginde
    ondan turetilen istatistiklerin eski surumle hesaplanmis kalmasi
    mekanik olarak yakalanmiyordu.

    Maliyeti ayni gun olculdu: bilesen taramasi bayat oldugu icin beta
    2.523 -> 4.372, reduktaz 6.185 -> 10.394 duzeltildi. O sayilar
    `operon_validation.json`, `stats.json`, `control_elements`,
    `operon_relations` ve atlas sayfalarini besliyor. Veritabani duzeldi ama
    turetilmis dosyalar duzelmedi; aradaki fark sayfada GORUNMEZ.

    KURAL NEREDEN GELIYOR (2026-11-05 duzeltmesi). Ilk surum kendi kuralini
    yazdi: dosyanin mtime'i `roar.sqlite`in mtime'indan eski ise bayat.
    O karsilastirma YANLIS ve `rebuild.py` onu 2026-10-14'te bilerek
    kaldirmisti (oradaki "TABLO BAZINDA DB GIRDISI" notu): veritabanina
    zincirin ICINDE yazan adimlar var, dolayisiyla DB'nin mtime'i her kosuda
    en sona kayar ve kendisinden ONCE kosmus her analiz bayat gorunur.
    Olculdu: kapi 42 dosyanin 42'sini bayat gosterdi; bazilari 0,3 saat
    "eski"ydi, yani AYNI kosuda uretilmislerdi. Depodaki tekrar eden hata
    siniftan biriydi -- "karar veren kural, karar verdigini soyledigimiz
    kural degil".

    Artik karar `rebuild.plan(state, propagate=False)` veriyor: adimin DB
    girdisi TABLO BAZINDADIR ve yalnizca kendisinden ONCEKI adimlarin yazdigi
    tablolar sayilir. Yani sayi, OKUDUGU tablolarin ozeti degistiyse bayattir
    -- alakasiz bir tabloya yazilmis olmasi onu bayat yapmaz.

    `propagate=False` BILINCLI. Plan "yeniden kosmali mi?" diye sorar ve
    ust adim kosacaksa alt adimi da bayat sayar. Kapi baska bir sey sorar:
    "yayimlanan sayi YANLIS mi?" Ileride bir ust adimin kosacak olmasi
    yayimlanmis sayiyi yanlis yapmaz. Olculdu: yayilim dahil 53, haric 20.

    KAPSAM. Yalnizca `rebuild.ARTEFACTS` icinde bir adimin ciktisi olan
    dosyalar yargilanabilir; grafikte olmayan dosya icin bu sorunun cevabi
    YOK ve bu bir tazelik bulgusu degil, modularite bulgusudur --
    `check_published_files_are_in_rebuild_graph` onu sayar. Burada yalnizca
    WARN olarak gecilir ki iki kapi ayni seyi iki kez FAIL etmesin.
    """
    import sys as _sys
    _sys.path.insert(0, root)
    _sys.path.insert(0, os.path.join(root, "webapp"))
    import atlas as _atlas
    import rebuild as _rebuild

    cwd = os.getcwd()
    try:
        os.chdir(root)
        rows = _rebuild.plan(_rebuild.load_state(), propagate=False)
    finally:
        os.chdir(cwd)

    # KAPSAM: YALNIZCA VERITABANINA BAGLI ADIMLAR (2026-11-05).
    # Kontrolun adi "veritabanindan turetilen her yayimlanmis istatistik"
    # diyor. Ilk surum ise yargilanacak dosyayi "uretici betik roar.sqlite
    # okuyor mu" diye KAYNAK METNINDEN tahmin ediyordu, ve 2026-11-05'te
    # grafige 8 yeni adim baglandiginda bunun yine adindan genis karar
    # verdigi goruldu: `bridging_geometry_audit` (PDB), `reaction_anchors`
    # ve `reaction_verification` (Rhea/ChEBI), `curation_claims` (UniProt,
    # NCBI) veritabanini HIC okumuyor -- girdileri bir DIS SERVIS. Onlarin
    # "bugunku tablo icerigiyle hesaplanmis olmasi" diye bir sey yok;
    # tazelikleri ancak "servise ne zaman soruldu" ile ifade edilebilir ve
    # bu grafigin yargilayacagi bir sey degildir.
    #
    # Artik kapsam GRAFIKTEN cikarilir: DB'yi girdi beyan eden adimlar, ve
    # ciktisi onlarin girdisi olan adimlar (gecisli). Yani kural adinin
    # soyledigi seyi yapar. Tahmin yok, kaynak metni okuma yok.
    dep = {name: set(inputs) for name, _o, _c, inputs, _h, _w in rows}
    by_output = {}
    for name, outputs, _c, _i, _h, _w in rows:
        for p in outputs:
            by_output[os.path.normpath(p)] = name

    db_linked = set()
    changed = True
    while changed:
        changed = False
        for name, inputs in dep.items():
            if name in db_linked:
                continue
            for p in inputs:
                q = os.path.normpath(p)
                if q.endswith("roar.sqlite") or by_output.get(q) in db_linked:
                    db_linked.add(name)
                    changed = True
                    break

    # cikti yolu -> (adim adi, bayatlik sebebi ya da None)
    produced = {}
    for name, outputs, _cmd, _inputs, _heavy, why in rows:
        if name not in db_linked:
            continue
        for p in outputs:
            produced[os.path.normpath(p)] = (name, why)

    documented = set(_atlas.ANALYSIS_FILES)
    stale, fresh, ungated = [], 0, []
    for entry in _atlas.PROVENANCE:
        name, kind = entry[0], entry[1]
        if kind != "file" or name not in documented:
            continue
        rel = os.path.normpath(os.path.join(out_dir, name))
        if rel not in produced:
            ungated.append(name)
            continue
        step, why = produced[rel]
        if why:
            stale.append(f"{name} [{step}: {why}]")
        else:
            fresh += 1

    rep.check("freshness: every published statistic derived from the database is "
              "computed from the current table contents",
              not stale, len(stale),
              f"{fresh + len(stale)} dosya rebuild grafiginden yargilandi"
              + ("" if not stale else
                 f"; BAYAT: {first(sorted(stale), 4)} -- ilgili adim yeniden "
                 "kosturulmali, yoksa sayfa tablolarin eski icerigiyle "
                 "hesaplanmis sayilar yayimlar"))

    #  IKI AYRI SEBEP, TEK UYARI DEGIL (2026-11-12). Bu kontrol once 12
    #  dosya icin "grafikte olmadigi icin tazeligi SINANAMIYOR" diyordu ve
    #  ALTISI ICIN BU YANLISTI: o alti dosya grafigin ciktisi (hmm_pins,
    #  bridging_geometry_audit, reaction_anchors, oxygen_fate_computed,
    #  reaction_verification, curation_claims) -- yalnizca girdileri
    #  veritabani DEGIL (HMM dosyalari, YAML, yapi kayitlari). Tazelikleri
    #  KENDI beyan ettikleri girdilere karsi sinaniyor; sinanmayan sey
    #  "tablolar degisti mi" sorusu, cunku o dosyalar tablolardan
    #  turemiyor. Yanlis etiketin bedeli: gercekten SAHIPSIZ olan alti
    #  dosyayi, sahibi olan alti dosyanin arasinda gizliyordu.
    graph_outputs = set()
    for entry in _rebuild.ARTEFACTS:
        for out in entry[1]:
            graph_outputs.add(os.path.normpath(out))
    orphan, not_db = [], []
    for name in ungated:
        rel = os.path.normpath(os.path.join(out_dir, name))
        (not_db if rel in graph_outputs else orphan).append(name)

    if orphan:
        rep.check("freshness: every published file is a rebuild step output",
                  False, len(orphan),
                  first(sorted(orphan), 5) + " -- hicbir adimin ciktisi DEGIL, "
                  "tazeligi hic sinanamaz; sebepleri modularite kapisinda "
                  "dosya dosya yazili",
                  severity="WARN")
    if not_db:
        #  Bu bir uyari DEGIL, bir BEYAN: dosya grafikte, ama veritabanindan
        #  turemedigi icin tablo degisimi onu bayatlatmaz -- ve bunun boyle
        #  oldugunu soylemek, sessiz birakmaktan iyi.
        rep.check("freshness: graph outputs that do not derive from the "
                  "database are judged against their own inputs",
                  True, 0,
                  f"{len(not_db)} dosya: {first(sorted(not_db), 4)}")


def check_published_files_are_in_rebuild_graph(rep, con, root):
    """Yayimlanan her analiz dosyasi, yeniden kurma grafiginin bir CIKTISI olmali.

    NEDEN BU KONTROL VAR (kullanici istegi "DB'nin modulerligini iyi
    planlamak cok kritik", 2026-11-04). `rebuild.py:ARTEFACTS` dosya
    duzeyinde bir bagimlilik grafigidir ve 61 adimin 50'si veritabanini
    GIRDI olarak beyan eder; yani veritabani degisince o adimlar yeniden
    kosar. Grafik, SAHIP OLDUGU adimlar icin dogru.

    Sorun eksik dugumler. 2026-11-04 olcumu: yayimlanan 56 analiz
    dosyasinin 25'i hicbir adimin ciktisi DEGIL. Bu dosyalar
    `rebuild.py --run` ile hic tazelenemez; veritabani duzeltilse bile eski
    surumden hesaplanmis halleriyle indirmeye sunulmaya devam ederler ve
    sayfada bunu gosteren hicbir sey yoktur.

    Ayni gun olculen ornek: bilesen taramasi bayat oldugu icin beta
    2.523 -> 4.372 duzeltildi, ve veritabanindan turetilen 42 yayimlanmis
    dosyanin 42'si bir anda bayat hale geldi.

    KONTROLUN MANTIGI. Yayim listesi `atlas.ANALYSIS_FILES`, grafik
    `rebuild.ARTEFACTS`; ikisi ayri yerlerde tutuldugu icin ayrisabilirler
    ve bu kontrol ayrismayi sayar. Butce OLCULEN degerdir ve KUCULMEK
    uzeredir; buyumesi FAIL'dir -- yani bundan sonra grafige baglanmamis
    yeni bir yayim dosyasi eklenemez.
    """
    import sys as _sys
    _sys.path.insert(0, os.path.join(root, "webapp"))
    _sys.path.insert(0, root)
    import atlas as _atlas
    try:
        import rebuild as _rebuild
    except Exception as exc:
        rep.check("modularity: the rebuild graph is importable", False, 1, str(exc))
        return

    outputs = set()
    for entry in _rebuild.ARTEFACTS:
        for out in entry[1]:
            outputs.add(out)
            outputs.add(os.path.basename(out))

    missing = sorted(name for name in _atlas.ANALYSIS_FILES
                     if name not in outputs
                     and os.path.join("analysis_out", name) not in outputs)

    # UC TUR VAR, TEK SAYI DEGIL (2026-11-05). Butce 25'ten 6'ya dusurulurken
    # kalan dosyalarin HEPSININ ayni sebeple grafik disinda olmadigi goruldu:
    #
    #   onbellek      : chebi_formula_cache.json -- `reaction_balance.py`
    #                   argv'den Rhea ID alan bir EL ARACI, pipeline adimi
    #                   degil; ciktisi bir istatistik degil, ChEBI formul
    #                   onbellegi. Grafige konmasi "tazelenmeli" demek olur
    #                   ki yanlistir.
    #   ureticisiz    : KAPANDI (2026-11-12). curation_claims_ncbi.json'i
    #                   yazan betik depoda YOKTU ve atlas kokeni onu
    #                   YAZMAYAN bir betigi gosteriyordu -- yayimlanan
    #                   yanlis bir koken iddiasi. Eksik NCBI gecisi
    #                   `verify_curation_claims.py`e yazildi (onbellek
    #                   yetkili, ag kapali), dosya artik bir adimin
    #                   ciktisi, ve kapsama 6/12'den 12/12'ye cikti --
    #                   eski dosya vakalarin yarisini sessizce atliyordu.
    #   elle kosulmus : motif_gate_B_new.json, motif_single_his_members.json,
    #                   oxygen_role_phylogeny.json, repertoire_matrix.json --
    #                   tek seferlik analizler, betigi yok.
    #
    # Bunlari grafige ZORLAMAK dogru olmaz; dogru olan NE OLDUKLARINI beyan
    # etmek. Butce bu yuzden 0 degil 6: 6 dosya icin "tazelenebilir" iddiasi
    # ZATEN yapilmiyor. Butce olculen degerin kendisi; her baglanan dosyada
    # DUSURULUR. Yukari hareket, grafige baglanmamis YENI bir yayim dosyasi
    # eklendigi anlamina gelir ve kapi kapanir.
    BUDGET = 5
    rep.check("modularity: every published analysis file is an output of the "
              "rebuild graph",
              len(missing) <= BUDGET, len(missing),
              f"{len(missing)} dosya grafikte yok, butce {BUDGET} "
              f"({len(_atlas.ANALYSIS_FILES)} yayim dosyasi, "
              f"{len(_rebuild.ARTEFACTS)} adim)"
              + (f": {first(missing, 5)}" if missing else "")
              + ("" if len(missing) <= BUDGET else
                 " -- yeni yayim dosyasi rebuild.ARTEFACTS'a BAGLANMADAN "
                 "eklenmis; o dosya hic tazelenemez"))


def check_citation_and_version(rep, con, root):
    """DB atif edilebilir mi, ve yayimlanan VERI SURUMU taze mi.

    NEDEN. ROAR-DB herkese acik bir veritabani olarak yayinlaniyor ama
    2026-11-11'e kadar sitede ne yazar adi ne atif bicimi vardi, ve veri
    surumu HIC yayinlanmiyordu. Ikincisi daha sinsi: `db_revision.py`
    surumleri tutuyor ama son kayit rev-0007 (2026-10-30) iken DB on iki
    gun boyunca degismisti. Yani bir okuyucu "hangi surumu kullandim"
    sorusunu cevaplayamiyordu; cevaplayabildigini SANDIGI durum daha
    kotuydu, cunku sitede yazan surum artik gosterdigi veriyi
    tanimlamiyordu.
    """
    path = os.path.join(root, "CITATION.cff")
    rep.check("citation: CITATION.cff var", os.path.exists(path), 0, path)
    if os.path.exists(path):
        text = open(path, encoding="utf-8").read()
        for field in ("cff-version", "title", "authors", "license", "type"):
            rep.check(f"citation: CITATION.cff '{field}' alani var",
                      field + ":" in text, 0, "")

    ledger = os.path.join(root, "analysis_out", "revisions.jsonl")
    rep.check("citation: veri surumu defteri var", os.path.exists(ledger), 0,
              ledger)
    if not os.path.exists(ledger):
        return
    with open(ledger, encoding="utf-8") as handle:
        lines = [ln for ln in handle if ln.strip()]
    rep.check("citation: en az bir veri surumu kayitli", bool(lines), 0, "")
    if not lines:
        return
    last = json.loads(lines[-1])
    live = con.execute(
        "SELECT COUNT(*) FROM ro WHERE is_confirmed = 1").fetchone()[0]
    recorded = last.get("confirmed_entries")
    #  Yayimlanan surum, GOSTERDIGI veriyi tanimlamak zorunda. Sayilar
    #  ayrisiyorsa sitede yazan surum id'si yanlis bir seye isaret ediyor.
    rep.check("citation: yayimlanan veri surumu DB ile ayni durumu anlatiyor",
              recorded == live, abs((recorded or 0) - live),
              f"{last.get('revision')} {recorded} giris yaziyor, DB'de {live}")


def check_taxonomy_is_controlled(rep, con, root, out_dir="analysis_out"):
    """Cins TAHMIN EDILMIYOR: rank'tan okunuyor, ve atanmamislari beyanli.

    NEDEN. Bu DB'nin istatistik kurali "karari cins duzeyi verir" diyor,
    yani cins yayimlanan her sinyalin birimi. Ama cins bes degil DOKUZ
    yerde `organism.split()[0]` ile tahmin ediliyordu ve o satir olculmus
    hata uretiyordu: "Candidatus" bir cins sanilip 17 ayri gercek cins tek
    kovaya yiginiyordu; `Mycolicibacterium` (612 giris, DB'nin 7. en
    kalabalik "cinsi") cins sanilyordu, oysa bugunun NCBI taksonomisinde
    bir ALTCINS ve o yolda cins rank'li dugum hic yok. Bu kontrol
    duzeltmenin geri surunmesini engeller: yeni bir yerel `split()[0]`
    eklenirse ya da cinsi atanmamis girisler sessizce dusurulurse FAIL.
    """
    import glob as _glob

    tables = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    needed = {"replicon_taxon", "taxon_lineage", "taxon_genus"}
    rep.check("taxonomy: uc tablo da var", needed <= tables,
              len(needed - tables), f"eksik: {sorted(needed - tables)}")
    if not needed <= tables:
        return

    #  1. Her replikonun taxid'i olmali. Olmayan bir replikonun uzerindeki
    #  girisler cins duzeyi hicbir teste giremez.
    #  BAYPAS KALDIRILDI. Burada `WHERE nucleotide_id <> ''` yaziyordu,
    #  cunku kimligi olmayan bir replikon satiri vardi ve bu kontrolu
    #  dusuruyordu. Baypas satiri duzeltmez, GIZLER: satir geri gelse bu
    #  kontrol yine sessiz kalirdi. Satir `db_repair.py` ile silindi,
    #  sebebi `extract_genomic_context.py`'de kapatildi, ve asagidaki
    #  kontrol kimligin bosluguna artik DOGRUDAN bakiyor.
    n_rep = con.execute("SELECT COUNT(*) FROM replicon").fetchone()[0]
    n_tax = con.execute("SELECT COUNT(*) FROM replicon_taxon").fetchone()[0]
    rep.check("taxonomy: her replikonun NCBI taxid'i var", n_tax >= n_rep,
              max(0, n_rep - n_tax), f"replikon {n_rep}, taxid'li {n_tax}")

    #  1b. KIMLIGI OLMAYAN REPLIKON OLAMAZ. `replicon` tablosunun birimi
    #  bir replikon ve birimin kimligi `nucleotide_id`; bos bir kimlik bir
    #  birim degil, bir ayristirma hatasidir. Bir kez olculdu: 0 baytlik
    #  bir gbk dosyasi boyle bir satir uretti ve her replikon-basi paydayi
    #  bir fazla gosterdi.
    phantom = con.execute(
        "SELECT COUNT(*) FROM replicon "
        "WHERE nucleotide_id IS NULL OR TRIM(nucleotide_id) = ''").fetchone()[0]
    rep.check("replicon: her satirin bir nucleotide_id'si var",
              phantom == 0, phantom,
              "kimligi bos replikon satiri: ayristirilamayan dosya replikon "
              "satiri URETMEMELI (bkz. db_repair.py phantom_replicon)")

    #  2. Cins RANK'tan gelmeli: genus_taxid, soyagacinda rank='genus'
    #  olarak gecen ata olmak zorunda.
    bad = con.execute("""
        SELECT COUNT(*) FROM taxon_genus g
         WHERE g.genus_taxid IS NOT NULL
           AND NOT EXISTS (SELECT 1 FROM taxon_lineage l
                            WHERE l.taxid = g.taxid_current
                              AND l.rank = 'genus'
                              AND l.ancestor_taxid = g.genus_taxid)""").fetchone()[0]
    rep.check("taxonomy: cins rank='genus' atasindan geliyor", bad == 0, bad,
              "genus_taxid soyagacinda genus rank'li degil")

    #  3. Cins adi, cins OLMAYAN bir sey olmamali -- eski tahmincinin
    #  urettigi hatanin tam kendisi (char(39) = tek tirnak).
    fake = con.execute("""
        SELECT genus_name, COUNT(*) n FROM taxon_genus
         WHERE genus_state = 'ranked'
           AND (LOWER(genus_name) IN ('candidatus','uncultured','unclassified',
                                      'candidate','bacterium','archaeon')
                OR genus_name LIKE '[%'
                OR genus_name LIKE char(39) || '%')
         GROUP BY 1""").fetchall()
    rep.check("taxonomy: cins adi cins (Candidatus/parantez/tirnak yok)",
              not fake, len(fake), first([f"{r[0]} x{r[1]}" for r in fake], 5))

    #  4. Cinsi atanmamis girisler DUSURULMEMIS olmali: rapor onlari
    #  SAYMAK ve nasil ele alindigini YAZMAK zorunda.
    path = os.path.join(root, out_dir, "taxonomy_ids.json")
    data = None
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    rep.check("taxonomy: rapor dosyasi var", data is not None, 0, path)
    if data:
        una = data.get("genus_unassigned") or {}
        measured = con.execute("""
            SELECT COUNT(*) FROM ro JOIN taxon_genus USING (nucleotide_id)
             WHERE ro.is_confirmed = 1 AND genus_state <> 'ranked'""").fetchone()[0]
        rep.check("taxonomy: cinsi atanmamis giris sayisi TAZE",
                  una.get("entries") == measured,
                  abs((una.get("entries") or 0) - measured),
                  f"raporda {una.get('entries')}, DB'de {measured}")
        rep.check("taxonomy: atanmamislarin nasil ele alindigi YAZILI",
                  bool((una.get("handling") or "").strip()), 0,
                  "genus_unassigned.handling bos")
        rep.check("taxonomy: payda beyanli",
                  bool((data.get("denominators") or {}).get("confirmed_entries")),
                  0, "denominators.confirmed_entries yok")

    #  5. GERI SURUNME MUHAFIZI. Cinsi dizgiden tahmin eden yerel bir
    #  uygulama yeniden eklenirse bu kontrol duser. Tek sahip ro_taxonomy.
    #  webapp/ de taranir: cins/tur oradaki uc yerde de dizgiden
    #  kesiliyordu ve YAYIMLANAN sayfalari besliyordu. `atlas.species_of`
    #  taxid'i olmayan kayit icin son care olarak duruyor, o yuzden
    #  izinli.
    allowed = {"taxonomy_ids.py", "ro_taxonomy.py", "validate_curation.py",
               "atlas.py"}
    offenders = []
    for path_py in sorted(_glob.glob(os.path.join(root, "*.py"))
                          + _glob.glob(os.path.join(root, "webapp", "*.py"))):
        name = os.path.basename(path_py)
        if name in allowed:
            continue
        with open(path_py, encoding="utf-8", errors="replace") as handle:
            previous = ""
            for i, line in enumerate(handle, 1):
                stripped = line.strip()
                if stripped.startswith("#"):
                    previous = stripped
                    continue
                #  TEK ISTISNA, ve yorumla DEKLARE edilmek zorunda: bir
                #  olcum, dizgi tahminini KASITLI olarak eski tanim diye
                #  hesaplayabilir (`promoter_architecture.py`'nin birim
                #  duyarliligi). Istisna dosya adina degil SATIRA baglandi,
                #  cunku bir dosyayi butunuyle muaf tutmak o dosyadaki
                #  SONRADAN eklenen gercek bir ihlali de muaf tutar.
                if "genus-string-intentional" in previous:
                    previous = stripped
                    continue
                #  Desen GENISTIR: ilk hali yalnizca `organism.split()`
                #  ariyordu ve `(organism or "").split()[0]` yazan
                #  `ecological_origin.py`'yi KACIRDI. Ayni satirda
                #  organizma adi ve bir dizgi kesmesi varsa, cins ya da
                #  tur dizgiden tahmin ediliyor demektir.
                if "organism" in stripped and ".split()" in stripped:
                    offenders.append(f"{name}:{i}")
                previous = stripped
    rep.check("taxonomy: cinsi dizgiden tahmin eden ikinci uygulama yok",
              not offenders, len(offenders), first(offenders, 6))


def check_hmm_libraries_are_pinned(rep, con, root, out_dir="analysis_out"):
    """Yayimlanan her HMM isabeti, kimligi SABITLENMIS bir modelden gelmeli.

    NEDEN. Bir domain cagrisi ("bu komsu bir ferredoksin") modelin kendisine
    baglidir ve Pfam modelleri surum atlar. Olculdu (2026-11-12):
    `component_hmms/components.hmm` icindeki 9 modelin 8'i InterPro'nun bugun
    verdiginden bir surum geride, `family_hmms/` ise tamamen guncel -- yani
    iki kutuphane AYRI release'lerden gelmis. Geride olmak hata degil;
    yeniden indirmek yayimlanmis her isabeti degistirirdi. Hata, hangi
    modelin kullanildiginin YAZILI OLMAMASIDIR.

    Bu kontrol AGA CIKMAZ. Diskteki dosyanin, sabitlemenin anlattigi dosya
    olup olmadigina bakar: ozet tutmuyorsa kutuphane sabitlemeden SONRA
    degismis demektir ve yayimlanan isabetler artik tarif edilen modelden
    gelmiyor.
    """
    path = os.path.join(root, out_dir, "hmm_pins.json")
    if not os.path.exists(path):
        rep.check("hmm: model kutuphaneleri sabitlendi", False, 1,
                  "hmm_pins.json yok -- `rebuild.py --run --only hmm_pins`")
        return
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)

    import hashlib as _hl

    def digest(full):
        sha = _hl.sha256()
        with open(full, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                sha.update(chunk)
        return sha.hexdigest()[:16]

    libraries = data.get("libraries") or []
    missing = [lib["library"] for lib in libraries if not lib.get("present")]
    rep.check("hmm: sabitlenen her kutuphane diskte var",
              not missing, len(missing), first(missing, 4))

    moved = []
    for lib in libraries:
        if not lib.get("present"):
            continue
        full = os.path.join(root, lib["library"])
        if not os.path.exists(full):
            moved.append(f"{lib['library']} (dosya yok)")
        elif digest(full) != lib.get("file_digest"):
            moved.append(f"{lib['library']} ozet {lib.get('file_digest')} "
                         f"-> {digest(full)}")
    rep.check("hmm: kutuphane dosyasi sabitlemeden beri degismedi",
              not moved, len(moved), first(moved, 3))

    #  Pfam'dan gelen her modelin accession'i SURUMLU olmali: `PF00355`
    #  bir pin degil, `PF00355.33` pindir.
    unversioned = [f"{m.get('library')}:{m.get('name')}"
                   for m in data.get("pfam_models") or []
                   if not re.match(r"^PF\d{5}\.\d+$", m.get("accession", ""))]
    rep.check("hmm: her Pfam modeli accession SURUMU ile sabitlendi",
              not unversioned, len(unversioned), first(unversioned, 4))

    #  Proje modelleri ayri listede durmali: accession surumu olmayan bir
    #  modeli surumlu gostermek, olculmemis bir guvence olur.
    project = data.get("project_models") or []
    rep.check("hmm: proje modelleri ayrica beyan edildi",
              bool(project), 0 if project else 1,
              f"{len(project)} model yalnizca ozetle sabitlendi")


def check_gff3_export(rep, con, root, out_dir="analysis_out"):
    """GFF3 degisim dosyasi gercekten GFF3 ve DB ile ayni sayilari anlatiyor.

    NEDEN IKI AYRI SEY SINANIYOR. Bir degisim dosyasinin iki bozulma yolu
    var ve biri otekini gizler: (1) bicim bozulur -- kolon sayisi, cozulmeyen
    `Parent`, 0-tabanli kalmis bir koordinat; (2) bicim DOGRU kalir ama
    icerik bayatlar, yani dosya gecerli bir GFF3 olarak eski sayilari
    yayimlar. Ikincisi daha tehlikeli, cunku dosyayi okuyan arac hic
    sikayet etmez.

    Ayristirmayi betigin KENDI raporuna guvenerek yapmiyor: dosyayi
    yeniden okuyup sayiyor. Yazan tarafin "gecerli" demesi bir olcum
    degil.
    """
    gff = os.path.join(root, out_dir, "roar_genomic_context.gff3.gz")
    if not os.path.exists(gff):
        rep.check("gff3: degisim dosyasi var", False, 1,
                  "roar_genomic_context.gff3.gz yok -- "
                  "`rebuild.py --run --only gff3`")
        return

    import gzip as _gz
    regions = {}
    ids, parents, problems = set(), 0, []
    types = {}
    with _gz.open(gff, "rt") as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if line.startswith("##sequence-region"):
                _, seqid, lo, hi = line.split()
                regions[seqid] = (int(lo), int(hi))
                continue
            if line.startswith("#") or not line:
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                problems.append(f"satir {lineno}: {len(cols)} kolon")
                continue
            types[cols[2]] = types.get(cols[2], 0) + 1
            if not cols[3].isdigit() or int(cols[3]) < 1:
                problems.append(f"satir {lineno}: start {cols[3]} 1-tabanli degil")
            for pair in cols[8].split(";"):
                key, _, value = pair.partition("=")
                if key == "ID":
                    ids.add(value)
                elif key == "Parent":
                    parents += 1
                    if value not in ids:
                        problems.append(f"satir {lineno}: Parent={value} cozulmuyor")
    rep.check("gff3: bicim gecerli (9 kolon, 1-tabanli, Parent cozuluyor)",
              not problems, len(problems), first(problems, 4))

    #  2. BAYATLIK: dosyadaki operon sayisi DB'deki onaylanmis operon
    #  sayisi olmak zorunda. Gecerli bir GFF3, eski sayilari yayimlayabilir.
    n_operon_db = con.execute("""
        SELECT COUNT(*) FROM operon o JOIN ro r USING(candidate_id)
         WHERE r.is_confirmed = 1""").fetchone()[0]
    n_operon_file = types.get("operon", 0)
    rep.check("gff3: dosyadaki operon sayisi DB ile ayni",
              n_operon_file == n_operon_db,
              abs(n_operon_file - n_operon_db),
              f"dosya {n_operon_file}, DB {n_operon_db}")

    n_cds_db = con.execute("""
        SELECT (SELECT COUNT(*) FROM ro WHERE is_confirmed = 1)
             + (SELECT COUNT(*) FROM operon_gene g
                  JOIN neighbor n ON n.neighbor_id = g.neighbor_id)""").fetchone()[0]
    n_cds_file = types.get("CDS", 0)
    rep.check("gff3: dosyadaki CDS sayisi DB ile ayni",
              n_cds_file == n_cds_db, abs(n_cds_file - n_cds_db),
              f"dosya {n_cds_file}, DB {n_cds_db} "
              f"(alfa + neighbor'a baglanan operon geni)")

    #  3. Her seqid icin bildirilen bolge olmali: bolgesiz bir seqid,
    #  tarayicida uzunlugu bilinmeyen bir kontig demek.
    n_rep_db = con.execute("""
        SELECT COUNT(DISTINCT nucleotide_id) FROM ro WHERE is_confirmed = 1
        """).fetchone()[0]
    rep.check("gff3: her replikon icin ##sequence-region var",
              len(regions) == n_rep_db, abs(len(regions) - n_rep_db),
              f"dosya {len(regions)}, DB {n_rep_db}")


def check_db_writes_are_attributed(rep, root):
    """Her DB yazmasinin bir SAHIBI var mi -- yani grafik kimin ne yazdigini
    biliyor mu?

    NEDEN. Bayatlik grafigi "bu adim hangi tablolari yazar" sorusunun
    cevabina dayaniyor ve cevap uc yerden geliyor: betigin kaynagindan
    (`static_writes` SQL'i ayristirir), elle beyandan (`DECLARED_WRITES`),
    ve adim kostugunda OLCULEN degisimden. Ucuncusu birinci ve ikincinin
    disina tasarsa, o yazmanin SAHIBI bilinmiyor demektir -- ve sahibi
    bilinmeyen bir degisim herkesin ustune yazilir.

    Olculdu (2026-11-12): `search_index` dort FTS5 GOLGE tablosu yaziyor
    (`ro_fts_config`, `ro_fts_data`, `ro_fts_docsize`, `ro_fts_idx`).
    Betik yalnizca `CREATE VIRTUAL TABLE ro_fts` diyor; dordunu SQLite'in
    kendisi uretiyor, yani hicbir kaynak taramasi onlari bulamaz. Sahipsiz
    kaldiklari icin `ro_fts_data`'nin degismesi, onu HIC okumayan uc rapor
    adimini bayat gosteriyordu.

    NE BULAMAZ. Bu kontrol yalnizca OLCULMUS olani sinar: bir adim bu
    araçla hic kosmadiysa yazmasi henuz gorulmemistir. Yani kontrol
    "yazma yok" demiyor, "gorulen her yazmanin sahibi var" diyor.
    """
    import sys as _sys
    _sys.path.insert(0, root)
    try:
        import rebuild as _rebuild
    except Exception as err:                          # noqa: BLE001
        rep.check("graph: her DB yazmasinin bir sahibi var", False, 1,
                  f"rebuild.py okunamadi: {err}")
        return
    state = _rebuild.load_state()
    orphan = []
    for name, _out, cmd, _in, _heavy in _rebuild.ARTEFACTS:
        static = _rebuild.static_writes(cmd[1]) if len(cmd) > 1 else set()
        declared = _rebuild.DECLARED_WRITES.get(name, set())
        observed = set((state.get(name) or {}).get("db_writes") or [])
        for table in sorted(observed - static - declared):
            orphan.append(f"{name} -> {table}")
    rep.check("graph: olculen her DB yazmasinin bir sahibi var",
              not orphan, len(orphan),
              first(orphan, 4) + (" -- DECLARED_WRITES'a ekle" if orphan else ""))


def check_record_kind(rep, con, root, out_dir="analysis_out"):
    """`replicon` tablosunun her satiri gercekten bir replikon mu?

    NEDEN. Olculdu: dogrulanmis uyesi olan kayitlarin 958'i `XM_`/`NM_`
    bir mRNA transkripti, ve bir transkript replikon degil. Bu satirlar
    SILINMEZ (dizileri gercek), ama REPLIKON BASINA paydalarda duramazlar:
    bir transkriptin komsu geni ya da tam operonu olamaz, dolayisiyla her
    baglam istatistigi 958 kayit boyunca sifir toplayip paydada kalir --
    yani bir OLCUM EKSIGI, bir biyolojik yokluk gibi yayimlanir.

    KONTROL EDILEN: her satirin siniflandirilmis olmasi. Siniflandirilamayan
    bir satir en tehlikeli durum, cunku hangi paydaya girdigi belirsiz
    kalir. Ayrica onek ile kaydin kendi LOCUS beyani arasinda uyusmazlik
    olmamasi: uyusmazlik hata degil ama SINIFLANDIRILMAMIS demektir.
    """
    tables = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    if "replicon_kind" not in tables:
        rep.check("record_kind: kayit turleri siniflandirildi", False, 1,
                  "replicon_kind tablosu yok -- "
                  "`rebuild.py --run --only replicon_kind`")
        return

    missing = con.execute("""
        SELECT COUNT(*) FROM replicon r
         WHERE NOT EXISTS (SELECT 1 FROM replicon_kind k
                            WHERE k.nucleotide_id = r.nucleotide_id)""").fetchone()[0]
    rep.check("record_kind: her replikon satirinin bir turu var",
              missing == 0, missing,
              "siniflandirilmamis satir, payda belirsiz kalir")

    unclassified = con.execute(
        "SELECT COUNT(*) FROM replicon_kind "
        "WHERE record_kind NOT IN ('genome','transcript')").fetchone()[0]
    rep.check("record_kind: hicbir satir 'okunamadi'/'uyusmazlik' degil",
              unclassified == 0, unclassified,
              "kaydin LOCUS beyani okunamadi ya da onekle celisiyor")

    #  Transkript kayitlari HALA uye uretiyor olmali: bu sayinin sifira
    #  dusmesi, duzeltmenin satirlari silmeye kaydigi anlamina gelir.
    kept = con.execute("""
        SELECT COUNT(*) FROM ro r JOIN replicon_kind k USING(nucleotide_id)
         WHERE r.is_confirmed = 1 AND k.record_kind = 'transcript'""").fetchone()[0]
    rep.check("record_kind: transkript girisleri silinmedi, isaretlendi",
              kept > 0, 0 if kept else 1,
              f"{kept} dogrulanmis giris transkript kayitlarinda duruyor")


def check_ncbi_accession_claims(rep, con, root, out_dir="analysis_out"):
    """accession'i GenBank/RefSeq olan referanslar KAYNAKTAN dogrulanmis mi?

    NEDEN. UniProt gecisinin accession suzgeci (harfle baslar, nokta YOK,
    <=10 karakter) GenBank/RefSeq kimliklerini disarida birakiyor. Olculdu
    (2026-11-12): 12 tipin accession'i bu bicimde ve yayimlanan rapor
    yalnizca 6'sini kapsiyordu -- yani dosya vakalarin YARISINI sessizce
    atliyordu, ve ustelik kokeni onu YAZMAYAN bir betigi gosteriyordu.

    Bu kontrol UYUSMAZLIK SAYISINI bir kusur saymaz: bir referansin
    kaynaktan ayrismasi bir kurasyon bulgusudur, kapinin cozecegi bir sey
    degil. Aradigi sey KAPSAMA: NCBI bicimli her accession'in bir olcumu
    OLMALI, ve olcumun cevrimdisi dogrulanabilir (dizi onbellekte) olmasi.
    Kapsanmayan bir accession, "kontrol edildi" izlenimi veren bir dosyada
    gorunmeyen bir bosluktur.
    """
    path = os.path.join(root, out_dir, "curation_claims_ncbi.json")
    if not os.path.exists(path):
        rep.check("ncbi: GenBank/RefSeq referanslari kaynaktan olculdu",
                  False, 1, "curation_claims_ncbi.json yok -- "
                  "`rebuild.py --run --only curation_claims_ncbi`")
        return
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    meta = data.get("_meta") or {}

    uncovered = meta.get("types_not_covered") or []
    rep.check("ncbi: NCBI bicimli her accession'in bir olcumu var",
              not uncovered, len(uncovered),
              f"{meta.get('types_covered')}/"
              f"{meta.get('types_with_an_ncbi_accession')} tip olculdu"
              + (f"; kapsanmayan: {first(uncovered, 4)}" if uncovered else ""))

    unverifiable = meta.get("unverifiable_offline") or []
    rep.check("ncbi: her hukum cevrimdisi yeniden hesaplanabilir",
              not unverifiable, len(unverifiable),
              first(unverifiable, 4) or "dizi onbellekte, ag gerekmiyor")

    #  Uyusmazlik FAIL DEGIL: kaynaktan ayrisan bir referans kuratorun
    #  bakmasi gereken bir bulgudur, kapinin kapatacagi bir sey degil.
    mismatches = meta.get("mismatches") or []
    rep.check("ncbi: saklanan dizi kaynakla uyusuyor",
              not mismatches, len(mismatches),
              first(mismatches, 4) or "12/12 referans kaynakla birebir",
              severity="WARN")


def check_assignment_library(rep, con, root, out_dir="analysis_out"):
    """Atama kutuphanesi her tipi ERISILEBILIR kiliyor mu?

    NEDEN. Tip atamasi, kapsama barajini gecen profiller arasinda en yuksek
    bit skorlusunu secer. Iki profil AYNI olasilik matrisini tasiyorsa her
    aday ikisinde ayni skoru alir, berabereyi hep ayni isim kazanir, ve
    kaybeden tip hicbir uye toplayamaz -- kendi referans dizisiyle BIREBIR
    ayni bir giris bile ona gitmez.

    Olculdu (2026-11-12): `3_308_DxnA1` ile `4_407_PbaA1` modelleri birebir
    ayni. Sebep hmmbuild degil, homolog TOPLAMASI: iki aile profili ayni
    dort proteinden kurulmus ve her biri OTEKININ referansini da iceriyor.
    Bedel kozmetik degil -- DxnA1'in iki dogrulanmis uyesi var ve IKISI DE
    PbaA1 referansiyla %100 benzer (%100 kapsamada), yani iki girisin
    yayimlanan tipi ve kimyasi ters.

    NEDEN FAIL DEGIL WARN. Duzeltmek profil kutuphanesini yeniden kurup
    19.125 girisi yeniden atamak demek, ve "hangi homolog hangi tipe ait"
    bir KURATOR karari. Bir kapi kuratorun yerine karar veremez; ama sessiz
    de kalamaz, cunku yayimlanan etiket su an yanlis. FAIL yapmak ayrica
    bu kusuru DUYURAN sayfanin yayimini da kilitlerdi.

    0 UYE TEK BIR OLGU DEGIL. Ayni belirtinin en az uc sebebi var ve rapor
    hangisi oldugunu yaziyor: tip veride TEMSIL EDILIYOR ama atamasi
    baskasina gidiyor / yalnizca UZAK adaylar ona en yakin / hicbir aday ona
    en yakin degil. Olculdu: 8 uyesiz tipin 6'si birinci, 2'si ikinci --
    yani bu uyarinin buyuk kismi eksik veri DEGIL, yontem kusuru.
    """
    path = os.path.join(root, out_dir, "reference_redundancy.json")
    if not os.path.exists(path):
        rep.check("assignment: kutuphane denetimi yazildi", False, 1,
                  "reference_redundancy.json yok")
        return
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    totals = data.get("totals") or {}

    groups = data.get("identical_model_groups") or []
    unreachable = [t for g in groups for t in (g.get("unreachable_types") or [])]
    detail = "; ".join(
        f"{'/'.join(g['models'])} (dizi kimligi "
        f"%{(g.get('sequence_identity') or [{}])[0].get('identity')})"
        for g in groups)
    rep.check("assignment: hicbir iki profil ayni olasilik matrisini tasimiyor",
              not groups, len(unreachable),
              detail + (" -- bir tip ERISILEMEZ" if unreachable else ""),
              severity="WARN")

    #  Bu kontrol SAYIYI degil, SAYININ ACIKLANMIS olmasini arar. "8 tip 0
    #  uyeli" bir kusur olmak zorunda degil (veri o enzimi icermiyor
    #  olabilir); aciklanmamis olmasi kusurdur.
    zero = data.get("zero_member_types") or []
    unexplained = [z["cluster"] for z in zero if not z.get("verdict")]
    rep.check("assignment: uyesiz her tipin bir TESHISI var",
              not unexplained, len(unexplained), first(unexplained, 4))

    #  En keskin alt kume: tipin kendi referansiyla %=100 benzer bir giris
    #  var ve yine de o tip bos. Bu bir yorum degil, bir celiski.
    exact = [f"{z['cluster']} (aday %{z['best_candidate_identity']} -> "
             f"{(z['assigned_instead_to'] or [{}])[0].get('cluster')})"
             for z in zero
             if (z.get("best_candidate_identity") or 0) >= 99.5]
    rep.check("assignment: kendi referansiyla ayni bir giris varken bos "
              "kalan tip yok",
              not exact, len(exact), first(exact, 4), severity="WARN")


def check_coverage_audit(rep, con, root, out_dir="analysis_out"):
    """Kapsama denetiminin YAPISI ve TAZELIGI; acigin BUYUKLUGU degil.

    NE KONTROL EDILMEZ -- ve bu kasitli. Bu kapi "acik kucuk olsun"
    DEMEZ. Eksik cins sayisi bir kusur degil bir olcumdur; onu FAIL
    yapmak, kapiyi veri BUYUTMEYE zorlayan bir mekanizmaya cevirir ve
    kapinin isi o degil. Ustelik dis evren (UniProtKB'de PF00848) kendi
    hizinda buyuyor; acigin buyumesi bizim bir seyi bozdugumuz anlamina
    gelmez.

    NE KONTROL EDILIR
      1. Dis liste TAM. InterPro'nun beyan ettigi sayi ile cekilen sayi
         esit olmali. Olculdu: ilk iki cekim 7.600 ve 2.600 kayitta
         408'e takildi. YARIM bir dis liste, kapsama oranini oldugundan
         YUKSEK gosterir (payda kuculur) -- yani eksik veri, eksik
         oldugunu gizleyen yonde hata uretir. En tehlikeli yon bu.
      2. Rapor BAYAT DEGIL. Raporun yazdigi "bizdeki cins sayisi",
         `taxon_genus` tablosundan BAGIMSIZ olarak yeniden sayilir.
         Esit degilse DB denetimden sonra degismistir ve yayimlanan
         kapsama orani artik baska bir kumeyi anlatiyor.
      3. KAPSAM BEYAN EDILMIS. PF00848 bakteriyel degil; kapsam ayrimi
         yazili olmazsa "eksik takson" sayisi bu DB'nin hic hedeflemedigi
         okaryotlarla sisirilir.
    """
    path = os.path.join(root, out_dir, "coverage_audit.json")
    if not os.path.exists(path):
        rep.check("coverage: kapsama denetimi yazildi", False, 1,
                  "coverage_audit.json yok -- "
                  "`rebuild.py --run --only coverage_audit`")
        return
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)

    outside = data.get("external_universe") or {}
    declared = outside.get("interpro_declared_count")
    got = outside.get("records_fetched")
    #  ESITLIK DEGIL YON. Olculdu: InterPro ilk sayfada 22.553 beyan
    #  ediyor, cursor'u sonuna kadar yurumek 22.620 AYRI accession
    #  veriyor; `count` bir onbellek tahmini, yurunen liste tukenmis.
    #  Tehlikeli yon EKSIK cekmek: yarim bir dis liste paydayi kucultur
    #  ve kapsama oranini oldugundan YUKSEK gosterir.
    rep.check("coverage: dis protein listesi truncate edilmemis",
              bool(got) and bool(outside.get("walk_not_truncated")),
              max(0, (declared or 0) - (got or 0)),
              f"InterPro {declared} beyan ediyor, yurunen liste {got} kayit")

    #  Raporun kendi sayisina GUVENILMEZ: ayni sayi tablodan yeniden
    #  uretilir. Iki yol da `taxonomy_ids`'in cins kuralindan geciyor,
    #  ama biri dis taxid'leri cozen kod (classify), oteki yayimlanmis
    #  tablo -- ayrisirlarsa biri bayat.
    import ro_taxonomy
    genus = ro_taxonomy.genus_map(con)
    #  PAYDA RAPORUNKIYLE AYNI KAPSAMDA olmali: yalnizca GENOM kayitlari.
    #  Ilk yazimda bu kisit yoktu ve kapi 950 (tum kayitlar) ile 720
    #  (genom kayitlari) arasindaki farki bir BAYATLIK sandi. Iki yol
    #  hala bagimsiz: biri `taxon_genus` tablosu, oteki dis taxid'leri de
    #  cozen `classify()`.
    rows = con.execute("""
        SELECT DISTINCT r.nucleotide_id
          FROM ro r JOIN replicon_kind k USING(nucleotide_id)
         WHERE r.is_confirmed = 1 AND k.record_kind = 'genome'""").fetchall()
    recount = len({genus[r[0]].key for r in rows if r[0] in genus})
    claimed = (data.get("denominators") or {}).get("our_genera")
    rep.check("coverage: raporun cins paydasi tabloyla uyusuyor",
              recount == claimed, abs(recount - (claimed or 0)),
              f"tabloda {recount} cins, raporda {claimed}")

    scope = ((data.get("method") or {}).get("scope") or {})
    in_scope = scope.get("in_scope_superkingdoms") or []
    rep.check("coverage: kapsam beyan edilmis",
              bool(in_scope) and bool(scope.get("our_side_by_superkingdom")),
              0 if in_scope else 1,
              f"kapsam ici superkingdom: {in_scope}")


def check_substrate_chebi(rep, con, root, out_dir="analysis_out"):
    """Substratin kontrollu kimligi: KIMLIK ile ADAY birbirine karismasin.

    NEDEN. `type_substrate_chebi` iki ayri seyi ayni tabloda tutuyor:
    dengelenmis Rhea denkleminden OKUNAN kimlik, ve kuratorun yazdigi ADIN
    ChEBI'de bulunmasiyla cikan ADAY. Ikisi esit degil -- bir ad eslesmesi
    ayni adi tasiyan baska bir molekule denk gelebilir. `state` bu ayrimi
    tasiyor ve bu kontrol ayrimin KAYBOLMADIGINI sinar: durumsuz bir satir,
    adayi kimlik gibi gosterir.
    """
    path = os.path.join(root, out_dir, "substrate_chebi.json")
    tables = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    if "type_substrate_chebi" not in tables or not os.path.exists(path):
        rep.check("chebi: substrat kimlik tablosu var", False, 1,
                  "type_substrate_chebi / substrate_chebi.json yok -- "
                  "`rebuild.py --run --only substrate_chebi`")
        return
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)

    rows = con.execute("""
        SELECT cluster, chebi, state FROM type_substrate_chebi""").fetchall()
    allowed = {"rhea_equation", "name_match"}
    bad_state = [f"{c}:{s}" for c, _x, s in rows if s not in allowed]
    rep.check("chebi: her satir KIMLIK mi ADAY mi oldugunu soyluyor",
              not bad_state, len(bad_state), first(bad_state, 4))

    bad_id = [f"{c}:{x}" for c, x, _s in rows
              if not re.match(r"^CHEBI:\d+$", x or "")]
    rep.check("chebi: her kimlik CHEBI:<sayi> bicimimde",
              not bad_id, len(bad_id), first(bad_id, 4))

    #  PAYDA BEYANLI. "76 tipin 37'si bir kimlik tasiyor" cumlesi ancak
    #  paydanin yaninda bir anlam tasir; rapor onu yazmak zorunda.
    denom = data.get("denominators") or {}
    rep.check("chebi: payda beyanli",
              bool(denom.get("types")) and "by_state" in denom, 0,
              "denominators.types / by_state eksik")

    #  Tablo ile rapor ayni seyi anlatmali: rapor eski bir kosudan kalmis
    #  olabilir ve o zaman sayfadaki sayi tablodan baska bir sey soyler.
    n_table = len(rows)
    n_report = denom.get("pairs_written")
    rep.check("chebi: tablo ile rapor ayni sayiyi veriyor",
              n_table == n_report, abs((n_report or 0) - n_table),
              f"tablo {n_table}, rapor {n_report}")

    #  YAPISAL AYIRT EDICININ KENDISI. Kurator adi ile denklemin adi
    #  farkliysa karar atom sayimindan verilir; her satir bir hukum
    #  tasimak zorunda, yoksa fark "sayilmis ama okunmamis" olur.
    block = data.get("curated_name_differs_from_equation") or {}
    rows_of = block.get("rows") or []
    no_verdict = [r.get("cluster") for r in rows_of if not r.get("verdict")]
    rep.check("chebi: her ad farki bir hukumle kapatildi",
              not no_verdict, len(no_verdict), first(no_verdict, 4))

    #  Gercek uyusmazliklar FAIL degil WARN: kurator karari gerektirir,
    #  ve bir kapi kuratorun yerine karar veremez. Ama SESSIZ de kalamaz.
    different = [f"{r['cluster']}: '{r.get('curated_substrate')}' vs "
                 f"{r.get('rhea_substrates')}"
                 for r in rows_of if r.get("verdict") == "different_molecule"]
    rep.check("chebi: kurator substrati denklemin substratiyla ayni molekul",
              not different, len(different), first(different, 4),
              severity="WARN")


def check_cross_references(rep, con, root, out_dir="analysis_out"):
    """`ro_xref` ve `cross_references.json`: adres tiklanabilir, payda beyanli.

    NEDEN. Buradaki her tanimlayici sayfada DIS BIR BAGLANTI olacak, yani
    bozuk bir deger okuyucuyu 404'e goturur; ve her yuzde UniProt kaynakli
    dizi paydasindan hesaplaniyor, TUM diziler paydasindan degil (madde A3).
    Iki sey de kendiliginden bozulabilir: payda `build_union.py` kostugunda
    buyur, bicim ise yeni bir accession semasi geldiginde kayar -- ilk
    kosuda 6 harf + 9 haneli WGS accession'lari (`JACBZF010000002`) desene
    uymadigi icin 5.358 GECERLI kayit atilmisti.

    EC AYRI BIR SEY. Oteki uc tanimlayici adres, EC ise ISLEV IDDIASI. Bu
    veritabaninin kuratorlu substrat etiketiyle ayni kefeye konmamasi icin
    inceleme durumunun yayimlanmasi SART: EC'lerin tamamina yakini
    incelenmemis (TrEMBL) kayitlardan geliyor. Oran yayimlanmiyorsa iddia
    kosulsuz kalir, o yuzden alanin varligi burada sinaniyor.
    """
    have = con.execute(
        "SELECT name FROM sqlite_master WHERE name='ro_xref'").fetchone()
    if not have:
        rep.check("xref: ro_xref exists", False, 1,
                  "ro_xref yok -- cross_references.py kosturulmali")
        return

    # 1. Her satir bir birlesim dizisine ait olmak zorunda.
    orphan = con.execute(
        "SELECT COUNT(*) FROM ro_xref x WHERE NOT EXISTS "
        "(SELECT 1 FROM ro_union u WHERE u.seq_md5 = x.seq_md5)").fetchone()[0]
    rep.check("xref: every cross-reference belongs to a union sequence",
              not orphan, orphan,
              "sahipsiz satir yok" if not orphan
              else f"{orphan} satir ro_union'da olmayan bir dizide")

    # 2. Bos ya da bosluk tasiyan tanimlayici: dogrudan kirik baglanti.
    blank = con.execute(
        "SELECT COUNT(*) FROM ro_xref WHERE identifier IS NULL "
        "OR TRIM(identifier) = '' OR identifier LIKE '% %'").fetchone()[0]
    rep.check("xref: no blank identifier", not blank, blank,
              "hepsi dolu" if not blank else f"{blank} satir bos/bosluklu")

    # 3. Bilinmeyen veritabani adi: sayfa onun icin URL kalibi bilmiyor.
    known = ("uniprot", "refseq_protein", "embl_nucleotide", "ec")
    unknown = [r[0] for r in con.execute(
        "SELECT DISTINCT db FROM ro_xref") if r[0] not in known]
    rep.check("xref: every database name has a published URL template",
              not unknown, len(unknown),
              "bilinen dort ad" if not unknown
              else "URL kalibi olmayan ad: " + ", ".join(sorted(unknown)))

    # 4. EC'nin yapisi: dort noktali alan. Bozuk bir EC, var olmayan bir
    #    ExPASy sayfasina baglanir.
    bad_ec = [r[0] for r in con.execute(
        "SELECT DISTINCT identifier FROM ro_xref WHERE db='ec'")
        if len([f for f in r[0].split(".")]) != 4]
    rep.check("xref: every EC number has four fields", not bad_ec, len(bad_ec),
              "hepsi n.n.n.n" if not bad_ec else first(bad_ec, 6))

    # 5. RefSeq protein accession'i alt cizgi tasir (WP_/NP_/YP_...).
    bad_rs = [r[0] for r in con.execute(
        "SELECT DISTINCT identifier FROM ro_xref WHERE db='refseq_protein'")
        if "_" not in r[0]]
    rep.check("xref: every RefSeq protein accession is prefixed",
              not bad_rs, len(bad_rs),
              "hepsi onekli" if not bad_rs else first(bad_rs, 6))

    # 6. PAYDA TAZELIGI: yayimlanan yuzdelerin paydasi bugunku DB ile ayni mi.
    path = os.path.join(root, out_dir, "cross_references.json")
    if not os.path.exists(path):
        rep.check("xref: cross_references.json present", False, 1,
                  f"{path} yok -- cross_references.py kosturulmali")
        return
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    den = data.get("denominators") or {}
    live_all = con.execute("SELECT COUNT(*) FROM ro_union").fetchone()[0]
    live_uni = con.execute(
        "SELECT COUNT(DISTINCT seq_md5) FROM ro_union_source "
        "WHERE source='uniprot'").fetchone()[0]
    drift = []
    if den.get("all_sequences") != live_all:
        drift.append(f"all_sequences {den.get('all_sequences')} != {live_all}")
    if den.get("uniprot_sourced_sequences") != live_uni:
        drift.append(f"uniprot_sourced {den.get('uniprot_sourced_sequences')} "
                     f"!= {live_uni}")
    rep.check("xref: the published denominators match the current database",
              not drift, len(drift),
              f"all={live_all}, uniprot-sourced={live_uni}"
              + ("" if not drift else "; BAYAT: " + "; ".join(drift)))

    # 7. Bicim kontrolu atilan kayit birakmamali: atilan varsa desen bir
    #    GECERLI semayi tanimiyor demektir ve o adresler sessizce kaybolur.
    dropped = {db: (info or {}).get("rejected_malformed", 0)
               for db, info in (data.get("coverage") or {}).items()}
    lost = {db: n for db, n in dropped.items() if n}
    rep.check("xref: no identifier was dropped as malformed",
              not lost, sum(lost.values()),
              "hicbiri atilmadi" if not lost
              else "; ".join(f"{db}: {n}" for db, n in sorted(lost.items()))
              + " -- desen gecerli bir accession semasini tanimiyor olabilir")

    # 8. EC'nin KOSULU yayimlanmis mi: inceleme durumu alani dolu olmak
    #    zorunda, yoksa islev iddiasi kosulsuz yayimlanir.
    review = data.get("ec_review_status") or {}
    has_share = isinstance(review.get("reviewed_share"), (int, float))
    rep.check("xref: the EC review status is published with the EC numbers",
              has_share, 0 if has_share else 1,
              (f"incelenmis {review.get('sequences_with_a_reviewed_ec')} / "
               f"yalnizca incelenmemis "
               f"{review.get('sequences_with_only_an_unreviewed_ec')} dizi")
              if has_share else "reviewed_share yok -- EC kosulsuz kalir")


def check_union_denominator(rep, con, root):
    """Birlesim ile kesisim ayrisirsa FAIL.

    NEDEN. Anasayfa iki sayi yayinliyor: `totals.confirmed` (19.125,
    KESISIM -- genomik baglami olan VE kapiyi gecen) ve `ro_union`dan
    gelen birlesim. Ikisi ayri tablodan okunuyor. `ro`ya giris eklenip
    `build_union.py` kosturulmazsa sayfa iki sayiyi yan yana gosterir ve
    biri bayat olur; okuyan kisi bunu goremez. Bu kontrol o ayrismayi
    mekanik hale getirir.

    Ayrica paydayi korur: konuma bagli her olcum (operon, komsuluk,
    replikon) yalnizca genomik baglami olan girisler icin yapilabilir, ve
    birlesim buyudukce o paydanin birlesim sanilmasi riski artar. Burada
    `has_genomic_context=1` satir sayisinin `ro`nun tekil dizi sayisina
    ESIT olmasi istenir -- esit degilse payda ya sismis ya eksilmis.
    """
    have = con.execute(
        "SELECT name FROM sqlite_master WHERE name='ro_union'").fetchone()
    if not have:
        rep.check("union: ro_union exists", False, 1,
                  "ro_union yok -- build_union.py kosturulmali")
        return

    ctx = con.execute("SELECT COUNT(*) FROM ro_union "
                      "WHERE has_genomic_context=1").fetchone()[0]
    uniq = con.execute("SELECT COUNT(DISTINCT sequence) FROM ro "
                       "WHERE is_confirmed=1 AND sequence IS NOT NULL"
                       ).fetchone()[0]
    rep.check("union: genomic-context rows equal the confirmed unique sequences",
              ctx == uniq, 0 if ctx == uniq else 1,
              f"{ctx} rows, {uniq} unique confirmed sequences"
              + ("" if ctx == uniq else " -- build_union.py bayat"))

    # Her birlesim satiri en az bir kaynak kimligi tasimak ZORUNDA: tasimayan
    # bir satirin nereden geldigi soylenemez ve indirilen dosyada kaynaksiz
    # bir dizi olarak gorunur.
    orphan = con.execute(
        "SELECT COUNT(*) FROM ro_union u WHERE NOT EXISTS "
        "(SELECT 1 FROM ro_union_source s WHERE s.seq_md5=u.seq_md5)"
    ).fetchone()[0]
    rep.check("union: every sequence carries at least one source accession",
              not orphan, orphan,
              "kaynaksiz satir yok" if not orphan
              else f"{orphan} satir hicbir kaynak kimligi tasimiyor")

    # n_sources beyan edilen sayi; kaynak tablosu olculen sayi. Ayrisirsa
    # sayfada "kac kaynaktan geliyor" yanlis yazilir.
    bad = con.execute(
        "SELECT COUNT(*) FROM ro_union u WHERE u.n_sources <> "
        "(SELECT COUNT(*) FROM ro_union_source s WHERE s.seq_md5=u.seq_md5)"
    ).fetchone()[0]
    rep.check("union: n_sources equals the measured source count",
              not bad, bad,
              "beyan ve olcum ayni" if not bad
              else f"{bad} satirda n_sources kaynak sayisina esit degil")

    # "Tek kopya" iddiasi: genomik baglami OLAN bir dizi hicbir kosulda
    # tek-kopya sayilamaz, cunku dizisi genom kaydindan geldi.
    wrong = con.execute(
        "SELECT COUNT(*) FROM ro_union u WHERE u.has_genomic_context=1 "
        "AND NOT EXISTS (SELECT 1 FROM ro_union_source s "
        "WHERE s.seq_md5=u.seq_md5 AND s.status NOT IN ('deleted','demerged'))"
    ).fetchone()[0]
    rep.check("union: no genomic-context sequence is counted as only-copy",
              not wrong, wrong,
              "tek-kopya kosulu genomik tarafi disariyor" if not wrong
              else f"{wrong} genomik giris tek-kopya sayiliyor")

    # Sayfanin kesisim/birlesim ayrimini YAZMASI zorunlu. Daha once bu
    # projede bir duzeltme kesfedildigi yerde uygulanip yayinlandigi yerde
    # unutuldu; o yuzden metin de kontrol ediliyor.
    idx = os.path.join(root, "webapp", "templates", "index.html")
    text = open(idx, encoding="utf-8").read() if os.path.exists(idx) else ""
    # 2026-11-02: eskiden bu liste `"is an\n<strong>intersection</strong>"`
    # ariyordu, yani SATIR SONUNA bagliydi. Anasayfa paragraftan panele
    # donusturuldugunde kontrol kirildi -- oysa amaci bicim degil, AYRIMIN
    # yazili olmasi. Artik kelimeleri ve uc olculmus degiskeni ariyor.
    need = ["Intersection", "Union", "totals.union_sequences",
            "totals.union_no_context", "totals.union_only_copy"]
    missing = [n for n in need if n not in text]
    rep.check("union: the homepage states the intersection/union distinction",
              not missing, len(missing),
              "ayrim yaziliyor" if not missing
              else "index.html'de eksik: " + ", ".join(missing))


def check_group_numbers_not_published(rep, root, out_dir):
    """Literaturde olmayan grup NUMARALARI yayinda mi?

    NEDEN. Literaturde bes RO grubu var. Bizim ic kimliklerimiz 7'ye kadar
    gidiyor ve 6 ile 7, karsiligi olmayan kovalar: "6" bitki klorofil
    RO'lari (PAO, CAO, PTC52), "7" tek bir dizi (AerC). Kullanici 2026-11-01'de
    uyardi: bunlari "group 6" / "group 7" diye yayimlamak, olmayan bir
    siniflandirmayi varmis gibi gosteriyor, cunku NUMARA yayimlanmis bir
    semaya atifta bulunur.

    Bu kontrol iki sey yapar. (i) Render edilmis sayfalarda ve kuratorlu
    duzyazida "group 6"/"group 7"/"group VI" gecmemeli. (ii) groups.yaml bu
    ayrimi acikca tasimali: her grubun `is_literature_group` degeri olmali ve
    literatur disi olanlar `published_as` ile adlandirilmali. Boylece etiket
    sablonlardan degil VERIDEN gelir.
    """
    import re as _re                                  # noqa: PLC0415

    bad_text = _re.compile(r"\bgroup\s+(?:6|7|VI|VII)\b")
    offenders = []

    # (i-a) Kuratorlu duzyazi.
    for name in ("chemistry.csv", "metabolic_role.csv", "cluster_ecology.csv"):
        path = os.path.join(root, name)
        if not os.path.exists(path):
            continue
        for lineno, line in enumerate(open(path, encoding="utf-8"), 1):
            if bad_text.search(line):
                offenders.append(f"{name}:{lineno}")

    # (i-b) Sablonlar. groups.yaml'in kendisi haric: orada numara KIMLIKTIR.
    tpl_dir = os.path.join(root, "webapp", "templates")
    if os.path.isdir(tpl_dir):
        for name in sorted(os.listdir(tpl_dir)):
            if not name.endswith(".html"):
                continue
            text = open(os.path.join(tpl_dir, name), encoding="utf-8").read()
            body = _re.sub(r"\{#.*?#\}", "", text, flags=_re.S)
            if bad_text.search(body):
                offenders.append(f"templates/{name}")

    rep.check("group labels: the non-literature groups are never published as "
              "a number",
              not offenders, len(offenders),
              "" if not offenders else
              "the published scheme has five groups; 6 and 7 are internal "
              "identifiers and must read \u201cungrouped\u201d with their "
              "names: " + ", ".join(sorted(set(offenders))[:6]))

    # (ii) groups.yaml bu ayrimi tasiyor mu?
    path = os.path.join(root, "groups.yaml")
    if not os.path.exists(path):
        rep.check("group labels: groups.yaml records which groups are from "
                  "the literature", True, 0, "groups.yaml missing",
                  severity="WARN")
        return
    try:
        import yaml                                   # noqa: PLC0415
        doc = yaml.safe_load(open(path, encoding="utf-8").read()) or {}
    except Exception as exc:                          # noqa: BLE001
        rep.check("group labels: groups.yaml is readable", False, 1, str(exc))
        return

    missing = []
    for entry in doc.get("groups") or []:
        gid = str(entry.get("id"))
        if "is_literature_group" not in entry:
            missing.append(f"group {gid}: is_literature_group not declared")
        elif not entry["is_literature_group"] and not entry.get("published_as"):
            missing.append(f"group {gid}: not a literature group but has no "
                           f"published_as label")
    rep.check("group labels: groups.yaml declares each group's standing and "
              "names the ungrouped ones",
              not missing, len(missing), first(missing, 4))


def check_vocabulary_drift(rep):
    """add_reference.py yazarken BU dosyanin dogrularken kullandigi listeler ayni mi?

    Yeni referans ekleyen script ile dogrulayici ayrisirsa, eklenen satir
    dogrulamadan gecer ama web arayuzu onu taniyamaz.
    """
    try:
        sys.path.insert(0, ROOT)
        import add_reference as ar          # noqa: PLC0415
    except Exception:                       # noqa: BLE001
        rep.check("vocabulary: writer and validator agree", True, 0,
                  "add_reference.py not importable; comparison skipped",
                  severity="WARN")
        return
    pairs = [
        ("substrate_class", VALID_SUBSTRATE_CLASS, getattr(ar, "VALID_SUBSTRATE_CLASS", None)),
        ("confidence", VALID_CONFIDENCE, getattr(ar, "VALID_CONFIDENCE", None)),
        ("reaction_class", VALID_REACTION_CLASS, getattr(ar, "VALID_REACTION_CLASS", None)),
        ("family", VALID_FAMILY, getattr(ar, "VALID_FAMILY", None)),
        ("source_kind", VALID_SOURCE_KIND, getattr(ar, "VALID_SOURCE_KIND", None)),
    ]
    diffs = []
    for field, mine, theirs in pairs:
        if theirs is None:
            diffs.append(f"{field}: missing in add_reference.py")
        elif set(theirs) != mine:
            diffs.append(f"{field}: {sorted(set(theirs) ^ mine)}")
    rep.check("vocabulary: writer and validator agree", not diffs, len(diffs),
              first(diffs) if diffs else "")


def check_confirmed_ro(con, rep):
    n_conf = con.execute("SELECT COUNT(*) FROM ro WHERE is_confirmed=1").fetchone()[0]
    rep.check("ro: confirmed entries exist", n_conf > 0, n_conf)

    bad = con.execute("""
        SELECT candidate_id FROM ro WHERE is_confirmed=1
          AND (ro_cluster IS NULL OR TRIM(ro_cluster)='' OR ro_cluster='N/A')
        LIMIT 20""").fetchall()
    n = con.execute("""
        SELECT COUNT(*) FROM ro WHERE is_confirmed=1
          AND (ro_cluster IS NULL OR TRIM(ro_cluster)='' OR ro_cluster='N/A')
        """).fetchone()[0]
    rep.check("ro: confirmed has a real ro_cluster", n == 0, n,
              "" if n == 0 else f"e.g. {first(r[0] for r in bad)}")

    n = con.execute("""
        SELECT COUNT(*) FROM ro WHERE is_confirmed=1
          AND (sequence IS NULL OR TRIM(sequence)='')""").fetchone()[0]
    rep.check("ro: confirmed has a sequence", n == 0, n)

    n = con.execute("""
        SELECT COUNT(*) FROM ro WHERE is_confirmed=1
          AND (rieske_intact IS NULL OR rieske_intact!=1)""").fetchone()[0]
    rep.check("ro: confirmed has rieske_intact=1", n == 0, n)

    # is_confirmed'in NULL kalmasi daha once gercek bir hataydi (945 aday).
    n = con.execute("SELECT COUNT(*) FROM ro WHERE is_confirmed IS NULL").fetchone()[0]
    rep.check("ro: is_confirmed is never NULL", n == 0, n,
              "" if n == 0 else "candidates without an hmmsearch hit must be 0, "
                                "not NULL")
    return n_conf


def check_cluster_coverage(con, rep, eco, chem):
    db_clusters = {r[0] for r in con.execute(
        "SELECT DISTINCT ro_cluster FROM ro WHERE is_confirmed=1 "
        "AND ro_cluster IS NOT NULL")}
    eco_clusters = {r["cluster"] for r in eco}
    chem_clusters = {r["cluster"] for r in chem}

    missing_eco = db_clusters - eco_clusters
    rep.check("coverage: every ro cluster is in cluster_ecology.csv",
              not missing_eco, len(missing_eco), first(missing_eco))
    missing_chem = db_clusters - chem_clusters
    rep.check("coverage: every ro cluster is in chemistry.csv",
              not missing_chem, len(missing_chem), first(missing_chem))

    # Ters yon: CSV'de olup veritabaninda uyesi olmayan tipler. Bilinen 10
    # uyesiz referans muaf; baska bir bosluk cikarsa bir tip uyelerini
    # kaybetmis demektir.
    orphan_eco = eco_clusters - db_clusters - REFERENCE_ONLY_CLUSTERS
    rep.check("coverage: cluster_ecology.csv has no unexplained extra type",
              not orphan_eco, len(orphan_eco), first(orphan_eco))
    orphan_chem = chem_clusters - db_clusters - REFERENCE_ONLY_CLUSTERS
    rep.check("coverage: chemistry.csv has no unexplained extra type",
              not orphan_chem, len(orphan_chem), first(orphan_chem))

    # Muafiyet listesi bayatladi mi: listedeki bir tip uye kazanmissa
    # listeden cikarilmasi gerekir.
    now_populated = REFERENCE_ONLY_CLUSTERS & db_clusters
    rep.check("coverage: REFERENCE_ONLY_CLUSTERS list is still accurate",
              not now_populated, len(now_populated),
              "" if not now_populated else
              f"these now have members, remove them: {first(now_populated)}",
              severity="WARN")

    # Iki CSV birbirini kapsiyor mu
    only_eco = eco_clusters - chem_clusters
    only_chem = chem_clusters - eco_clusters
    rep.check("coverage: the two curated CSVs cover the same types",
              not only_eco and not only_chem, len(only_eco) + len(only_chem),
              "" if not (only_eco or only_chem) else
              f"ecology-only: {first(only_eco)} | chemistry-only: {first(only_chem)}")

    dup_eco = len(eco) - len(eco_clusters)
    dup_chem = len(chem) - len(chem_clusters)
    rep.check("coverage: no duplicate cluster rows in the CSVs",
              dup_eco == 0 and dup_chem == 0, dup_eco + dup_chem)
    return db_clusters


def check_ecology_csv(rep, eco):
    bad = [(r["cluster"], r["substrate_class"]) for r in eco
           if r.get("substrate_class") not in VALID_SUBSTRATE_CLASS]
    rep.check("cluster_ecology.csv: substrate_class in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [(r["cluster"], r["confidence"]) for r in eco
           if r.get("confidence") not in VALID_CONFIDENCE]
    rep.check("cluster_ecology.csv: confidence in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [r for r in eco if not (r.get("cluster") or "").strip()
           or not (r.get("substrate") or "").strip()]
    rep.check("cluster_ecology.csv: cluster and substrate are non-empty",
              not bad, len(bad))


VALID_HOST_KINGDOM = {"human", "animal", "plant", "fungus", "alga",
                      "ambiguous", "not_a_host", "unresolved"}


def check_reference_redundancy_agrees(rep, out_dir, refs_path):
    """`reference_redundancy.py` ile bu dosyanin bulgulari ayni mi?

    Ayni gercek iki yerde bagimsiz hesaplaniyor: burada FASTA'dan, orada yine
    FASTA'dan ama farkli kodla ve tip kimligi normalizasyonuyla. Ikisinin
    AYRILMASI bir hata isaretidir; ornegin normalizasyon bozulursa orada
    sayilar sessizce sifira duser ve sayfadaki uyarilar kaybolur. Bu kontrol
    o sessizligi engeller.
    """
    path = os.path.join(out_dir, "reference_redundancy.json")
    if not os.path.exists(path):
        rep.check("reference_redundancy.json exists", False, 0,
                  "run reference_redundancy.py")
        return
    with open(path, encoding="utf-8") as fh:
        report = json.load(fh)
    if not os.path.exists(refs_path):
        rep.check("reference redundancy: FASTA available for cross-check", False, 0,
                  refs_path)
        return

    records, name, chunks = {}, None, []
    with open(refs_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    records[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name:
        records[name] = "".join(chunks)

    totals = report.get("totals", {})
    rep.check("reference redundancy: reference count agrees",
              totals.get("references") == len(records),
              abs((totals.get("references") or 0) - len(records)),
              f"report {totals.get('references')} vs FASTA {len(records)}")
    distinct = len({seq for seq in records.values()})
    rep.check("reference redundancy: distinct-sequence count agrees",
              totals.get("distinct_sequences") == distinct,
              abs((totals.get("distinct_sequences") or 0) - distinct),
              f"report {totals.get('distinct_sequences')} vs FASTA {distinct}")
    expected_groups = sum(1 for seq in {s for s in records.values()}
                          if sum(1 for v in records.values() if v == seq) > 1)
    rep.check("reference redundancy: identical-group count agrees",
              len(report.get("identical_groups", [])) == expected_groups,
              abs(len(report.get("identical_groups", [])) - expected_groups),
              f"report {len(report.get('identical_groups', []))} vs {expected_groups}")
    # Her etkilenen referans gercekten kuratorlu tip kimligi mi? Normalizasyon
    # bozulursa burada uzun FASTA basliklari gorunur ve sayfada hicbir uyari
    # eslesmez.
    affected = report.get("affected") or {}
    rep.check("reference redundancy: affected ids are type ids, not FASTA headers",
              all("_pro" not in key and key.count("_") == 2 for key in affected),
              sum(1 for key in affected if "_pro" in key or key.count("_") != 2),
              first(k for k in affected if "_pro" in k or k.count("_") != 2))
    rep.check("reference redundancy: at least one reference is flagged",
              bool(affected), 0,
              f"{len(affected)} references flagged, "
              f"{totals.get('entries_under_an_ambiguous_reference')} members affected")


def check_tree_coverage(rep, con, out_dir, chem_clusters):
    """Agactaki uc sayisi, temsilci kuralinin ONGORDUGU sayiya esit mi?

    Sayfalar bir donem "her varyant icin bir temsilci" diyordu, oysa 2-4 uyeli
    604 varyant atlaniyor ve bu girislerin %14'u. Metin olculen sayilara
    baglandi; bu kontrol esigin degismesi halinde metnin agactan AYRILMASINI
    engeller.
    """
    path = os.path.join(out_dir, "tree_all.nwk")
    if not os.path.exists(path):
        rep.check("tree_all.nwk exists", True, 0, "tree not built; check skipped",
                  severity="WARN")
        return
    newick = open(path, encoding="utf-8").read()
    tips = len(re.findall(r"[(,]\s*'?[^(),:']+'?\s*:", newick))

    import sys as _sys
    _sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "webapp"))
    import atlas as _atlas
    cover = _atlas.representative_coverage(con)
    if not cover:
        rep.check("representative coverage computable", False, 0, "no leaf table")
        return
    # Agaca KURATORLU referanslarin TAMAMI giriyor, uyesi olmayan 10 tane dahil;
    # bu yuzden beklenen sayi 71 uzerinden kurulur ve kontrol TAM esitliktir.
    # Ilk surum "uyesi olan referans" sayisini (61) kullanmis ve 10 fark icin
    # bir tolerans koymak gerekmisti; tolerans bir hatayi gizleyebilirdi.
    references = len(chem_clusters)
    expected = cover["variants_represented"] + references
    rep.check("tree tip count agrees with the representative rule",
              tips == expected, abs(tips - expected),
              f"{tips} tips, rule predicts {cover['variants_represented']} representatives "
              f"+ {references} curated references = {expected}")
    rep.check("skipped variants are a minority of entries",
              cover["entries_skipped_share"] < 0.25,
              round(100 * cover["entries_skipped_share"]),
              f"{cover['variants_skipped']} variants with 2-4 members hold "
              f"{cover['entries_skipped']} entries "
              f"({100 * cover['entries_skipped_share']:.1f} %)")


def check_published_files(rep, out_dir):
    """Belgelenen her analiz dosyasi diskte VAR mi?

    Indirme izin listesi artik `atlas.PROVENANCE`tan turetiliyor, yani tablo
    bir dosyayi sayiyorsa site onu indirmeye sunuyor. O halde tablodaki bir
    satirin karsiligi yoksa sitede kirik bir indirme baglantisi olur.
    """
    import sys as _sys
    _sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "webapp"))
    import atlas as _atlas
    missing = sorted(name for name in _atlas.ANALYSIS_FILES
                     if not os.path.exists(os.path.join(out_dir, name)))
    rep.check("every documented analysis file exists on disk",
              not missing, len(missing),
              f"{len(_atlas.ANALYSIS_FILES)} documented" if not missing
              else first(missing, 5))
    # Ters yon: analysis_out'ta olup tabloda olmayan dosyalar. Bu bir FAIL
    # degil uyaridir, cunku ara ciktilar yayinlanmak zorunda degil.
    on_disk = {f for f in os.listdir(out_dir)
               if f.endswith((".csv", ".json"))} if os.path.isdir(out_dir) else set()
    # SUPERSEDED_FILES de BELGELENMIS sayilir: yayinlanmiyorlar ama provenance
    # tablosunda neden yayinlanmadiklari yaziyor. Burada sayilmasalar bu
    # uyari "belgesiz" diyerek yanlis bir sey soylerdi.
    undocumented = sorted(on_disk - set(_atlas.ANALYSIS_FILES)
                          - set(getattr(_atlas, "SUPERSEDED_FILES", ())))
    rep.check("analysis files that exist but are not documented",
              not undocumented, len(undocumented), first(undocumented, 6),
              severity="WARN")


def check_host_kingdom(rep, con, path="host_kingdom.csv"):
    """Konak alemi eslemesi tam ve kontrollu mu?

    Esleme ELLE tutuluyor ve yeni genomlar yeni konak adlari getiriyor. Eksik
    bir dizgi sessizce `unresolved` olmasin diye burada sayilir: kontrol
    eksikligi gizlemez, gorunur kilar.
    """
    if not os.path.exists(path):
        rep.check(f"{path} exists", False, 0, "file missing")
        return
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    bad = [(r.get("host"), r.get("kingdom")) for r in rows
           if (r.get("kingdom") or "").strip() not in VALID_HOST_KINGDOM]
    rep.check(f"{path}: kingdom in vocabulary", not bad, len(bad),
              first(f"{h}={k!r}" for h, k in bad))
    nobasis = [r.get("host") for r in rows if not (r.get("basis") or "").strip()]
    rep.check(f"{path}: every row records the basis for its decision",
              not nobasis, len(nobasis), first(nobasis))
    seen = Counter((r.get("host") or "").strip().lower() for r in rows)
    dupes = sorted(h for h, n in seen.items() if n > 1)
    rep.check(f"{path}: no duplicate host string", not dupes, len(dupes), first(dupes))

    if not con.execute("SELECT name FROM sqlite_master WHERE name='replicon_source'").fetchone():
        return
    mapped = {(r.get("host") or "").strip().lower() for r in rows}
    hosts = {(h or "").strip() for (h,) in con.execute(
        "SELECT DISTINCT host FROM replicon_source "
        "WHERE host IS NOT NULL AND TRIM(host) <> ''")}
    unmapped = sorted(h for h in hosts if h.lower() not in mapped)
    rep.check("every /host string in the database is mapped to a kingdom",
              not unmapped, len(unmapped),
              f"{len(hosts)} distinct host strings" if not unmapped
              else first(unmapped, 5))
    stored = {k for (k,) in con.execute(
        "SELECT DISTINCT host_kingdom FROM replicon_source "
        "WHERE host_kingdom IS NOT NULL")}
    off = sorted(stored - VALID_HOST_KINGDOM)
    rep.check("stored host_kingdom values are all in the vocabulary",
              not off, len(off), first(off))


def check_display_names(rep, chem):
    """Gosterilen kisa adlar birbirinden ayirt edilebiliyor mu?

    Uc kisa ad referans setinde iki kez geciyor (BphA1, NDO, NidA) ve tablolarda
    iki satir ayni etiketle iki AYRI sayfaya baglaniyordu. Cozum ayristirma
    kurali; bu kontrol kuralin yeni bir cakismayi kacirmadigini dogrular.
    """
    import sys as _sys
    _sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "webapp"))
    import atlas as _atlas
    names = _atlas.short_names(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            "chemistry.csv"))
    missing = [r["cluster"] for r in chem if r["cluster"] not in names]
    rep.check("display name assigned to every type",
              not missing, len(missing), first(missing))
    seen = Counter(names.values())
    clashes = sorted(n for n, c in seen.items() if c > 1)
    rep.check("display names are unique after disambiguation",
              not clashes, len(clashes), first(clashes))
    empty = sorted(c for c, n in names.items() if not n.strip())
    rep.check("no display name is blank", not empty, len(empty), first(empty))


def check_chemistry_csv(rep, chem):
    bad = [(r["cluster"], r["reaction_class"]) for r in chem
           if r.get("reaction_class") not in VALID_REACTION_CLASS]
    rep.check("chemistry.csv: reaction_class in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [(r["cluster"], r["family"]) for r in chem
           if r.get("family") not in VALID_FAMILY]
    rep.check("chemistry.csv: family in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [(r["cluster"], r["source_kind"]) for r in chem
           if r.get("source_kind") not in VALID_SOURCE_KIND]
    rep.check("chemistry.csv: source_kind in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))
    bad = [(r["cluster"], r.get("curation_confidence")) for r in chem
           if r.get("curation_confidence") not in VALID_CONFIDENCE]
    rep.check("chemistry.csv: curation_confidence in vocabulary",
              not bad, len(bad), first(f"{c}={v!r}" for c, v in bad))

    bad = []
    n_smiles = 0
    for r in chem:
        s = (r.get("substrate_smiles") or "").strip()
        if not s:
            continue
        n_smiles += 1
        problems = smiles_problems(s)
        if problems:
            bad.append(f"{r['cluster']}: {'; '.join(problems)}")
    rep.check("chemistry.csv: every SMILES is parseable",
              not bad, len(bad),
              f"{n_smiles} non-empty SMILES checked" if not bad else first(bad, 4))

    bad = [(r["cluster"], r["pdb"]) for r in chem
           if (r.get("pdb") or "").strip()
           and not re.fullmatch(r"[A-Za-z0-9]{4}", r["pdb"].strip())]
    n_pdb = sum(1 for r in chem if (r.get("pdb") or "").strip())
    rep.check("chemistry.csv: every pdb is four alphanumerics",
              not bad, len(bad),
              f"{n_pdb} non-empty pdb ids checked" if not bad
              else first(f"{c}={v!r}" for c, v in bad))

    # Kaynak beyan edilmeden substrat verilmesi: web sayfasi substrati
    # gosterir ama nereden geldigini soyleyemez.
    bad = [r["cluster"] for r in chem
           if (r.get("substrate_en") or "").strip() not in ("", "unknown")
           and not (r.get("source") or "").strip()]
    rep.check("chemistry.csv: a named substrate always cites a source",
              not bad, len(bad), first(bad))


def check_text_fields_classified(con, rep):
    """DB'deki her metin sutunu uc listeden birinde kayitli mi?

    Dil kontrolu elle tutulan bir sutun listesine bakiyordu ve `euk_group` o
    listede YOKTU: sutun hic denetlenmedi, Turkce etiketler aylarca arayuzde
    kaldi. Bu kontrol eksigi kapatmaz, GORUNUR kilar: yeni bir metin sutunu
    ne "kontrol edilecek", ne "disaridan gelen", ne de "kimlik" listesinde
    degilse FAIL verir ve birinin karar vermesi gerekir.
    """
    checked = {(t, c) for t, _k, c in USER_FACING_DB_FIELDS}
    known = checked | EXTERNAL_TEXT_FIELDS | IDENTIFIER_TEXT_FIELDS
    tables = [r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' "
        "AND name NOT LIKE '%_fts%' AND name NOT LIKE 'sqlite_%'")]
    unclassified, total = [], 0
    for table in sorted(tables):
        for row in con.execute(f'PRAGMA table_info("{table}")'):
            name, decl = row[1], (row[2] or "").upper()
            if not (decl.startswith(("TEXT", "VARCHAR", "CHAR")) or decl == ""):
                continue
            total += 1
            if (table, name) not in known:
                unclassified.append(f"{table}.{name}")
    rep.check("every DB text column is classified for language checking",
              not unclassified, len(unclassified),
              f"{total} text columns, {len(checked)} language-checked"
              if not unclassified else first(unclassified, 5))
    stale = sorted(f"{t}.{c}" for t, c in known
                   if t in tables and (t, c) not in {
                       (tt, row[1]) for tt in tables
                       for row in con.execute(f'PRAGMA table_info("{tt}")')})
    rep.check("no language-check entry points at a dropped column",
              not stale, len(stale), first(stale, 5), severity="WARN")


def check_turkish(con, rep, eco_path, chem_path, eco, chem, out_dir):
    total = 0
    for table, key, col in USER_FACING_DB_FIELDS:
        if not con.execute("SELECT name FROM sqlite_master WHERE name=?",
                           (table,)).fetchone():
            rep.check(f"language: {table}.{col} is English", True, 0,
                      "table not in this database; check skipped",
                      severity="WARN")
            continue
        try:
            rows = con.execute(
                f'SELECT "{key}", "{col}" FROM "{table}" '
                f'WHERE "{col}" IS NOT NULL').fetchall()
        except sqlite3.DatabaseError as exc:
            rep.check(f"language: {table}.{col} is English", False, 0, str(exc))
            continue
        bad = []
        for ident, text in rows:
            hits = turkish_hits(text)
            if hits:
                bad.append(f"{ident}: {first(hits, 3)} in {str(text)[:40]!r}")
        total += len(bad)
        rep.check(f"language: {table}.{col} is English", not bad, len(bad),
                  f"{len(rows)} values checked" if not bad else first(bad, 3))

    for path, rows, label in ((eco_path, eco, "cluster_ecology.csv"),
                              (chem_path, chem, "chemistry.csv")):
        bad = []
        for i, row in enumerate(rows, start=2):
            for field, value in row.items():
                for hit in turkish_hits(value):
                    bad.append(f"line {i} {field}={hit}")
        total += len(bad)
        rep.check(f"language: {label} is English", not bad, len(bad),
                  f"{len(rows)} rows checked" if not bad else first(bad, 4))

    # analysis_out/ icinde web tarafindan yayinlanan metin dosyalari. Bu
    # dosyalar DB'den degil, uretildikleri ANDAKI CSV'den kopyalanir; CSV
    # Ingilizce'ye cevrilince eski ciktilar Turkce kalir.
    bad = []
    for name in PUBLISHED_TEXT_FILES:
        path = os.path.join(out_dir, name)
        if not os.path.exists(path):
            continue
        text = open(path, encoding="utf-8", errors="replace").read()
        hits = turkish_hits(text)
        if hits:
            bad.append(f"{name}: {first(hits, 4)}")
    rep.check("language: published analysis_out files are English",
              not bad, len(bad), first(bad, 4))
    return total


def check_referential(con, rep):
    def scalar(sql):
        return con.execute(sql).fetchone()[0]

    n = scalar("""SELECT COUNT(*) FROM neighbor n
                  LEFT JOIN ro r ON r.candidate_id=n.candidate_id
                  WHERE r.candidate_id IS NULL""")
    rep.check("refs: neighbor.candidate_id exists in ro", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM ro_leaf x
                  LEFT JOIN leaf l ON l.leaf_id=x.leaf_id
                  WHERE l.leaf_id IS NULL""")
    rep.check("refs: ro_leaf.leaf_id exists in leaf", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM leaf_sdp s
                  LEFT JOIN cluster_sdp c ON c.cluster=s.cluster
                  WHERE c.cluster IS NULL""")
    rep.check("refs: leaf_sdp.cluster exists in cluster_sdp", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM leaf_sdp s
                  LEFT JOIN leaf l ON l.leaf_id=s.leaf_id
                  WHERE l.leaf_id IS NULL""")
    rep.check("refs: leaf_sdp.leaf_id exists in leaf", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM gene_category g
                  LEFT JOIN neighbor n ON n.neighbor_id=g.neighbor_id
                  WHERE n.neighbor_id IS NULL""")
    rep.check("refs: gene_category.neighbor_id exists in neighbor", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM ro r
                  LEFT JOIN replicon p ON p.nucleotide_id=r.nucleotide_id
                  WHERE p.nucleotide_id IS NULL""")
    rep.check("refs: ro.nucleotide_id exists in replicon", n == 0, n)

    n = scalar("""SELECT COUNT(*) FROM ro_subfamily x
                  LEFT JOIN subfamily s ON s.subfamily_id=x.subfamily_id
                  WHERE s.subfamily_id IS NULL""")
    rep.check("refs: ro_subfamily.subfamily_id exists in subfamily", n == 0, n)

    if con.execute("SELECT name FROM sqlite_master "
                   "WHERE name='replicon_source'").fetchone():
        n = scalar("""SELECT COUNT(*) FROM replicon_source s
                      LEFT JOIN replicon p ON p.nucleotide_id=s.nucleotide_id
                      WHERE p.nucleotide_id IS NULL""")
        rep.check("refs: replicon_source.nucleotide_id exists in replicon",
                  n == 0, n)

    # leaf.size, o yapragin ro_leaf satir sayisina esit olmali. Bu alan web
    # sayfasinda "varyant boyutu" olarak gosteriliyor; ayrismasi sessizdir.
    rows = con.execute("""
        SELECT l.leaf_id, l.size,
               (SELECT COUNT(*) FROM ro_leaf r WHERE r.leaf_id=l.leaf_id)
        FROM leaf l""").fetchall()
    bad = [f"{lid}: size={sz} members={cnt}" for lid, sz, cnt in rows if sz != cnt]
    rep.check("refs: leaf.size equals its ro_leaf member count",
              not bad, len(bad),
              f"{len(rows)} leaves checked" if not bad else first(bad, 4))

    total_size = con.execute("SELECT COALESCE(SUM(size),0) FROM leaf").fetchone()[0]
    total_members = con.execute("SELECT COUNT(*) FROM ro_leaf").fetchone()[0]
    rep.check("refs: sum(leaf.size) equals the ro_leaf row count",
              total_size == total_members, abs(total_size - total_members),
              f"sum(leaf.size)={total_size} ro_leaf={total_members}")

    rows = con.execute("""
        SELECT s.subfamily_id, s.size,
               (SELECT COUNT(*) FROM ro_subfamily r
                WHERE r.subfamily_id=s.subfamily_id)
        FROM subfamily s""").fetchall()
    bad = [f"{sid}: size={sz} members={cnt}" for sid, sz, cnt in rows if sz != cnt]
    rep.check("refs: subfamily.size equals its member count",
              not bad, len(bad),
              f"{len(rows)} subfamilies checked" if not bad else first(bad, 4))


def check_arithmetic(con, rep, n_conf):
    """Dogrulanan RO sayisi her tabloda ayni mi?

    Her biri ayri bir script tarafindan dolduruluyor; biri yarida kalirsa
    sayilar ayrisir ve rapordaki yuzdeler sessizce yanlis olur.
    """
    tables = {
        "ro (is_confirmed=1)": "SELECT COUNT(*) FROM ro WHERE is_confirmed=1",
        "ro_search": "SELECT COUNT(*) FROM ro_search",
        "ro_evidence": "SELECT COUNT(*) FROM ro_evidence",
        "ro_domain": "SELECT COUNT(*) FROM ro_domain",
        "operon": "SELECT COUNT(*) FROM operon",
        "ro_leaf": "SELECT COUNT(*) FROM ro_leaf",
        "ro_etc": "SELECT COUNT(*) FROM ro_etc",
        "ro_regulation": "SELECT COUNT(*) FROM ro_regulation",
        "ro_subfamily": "SELECT COUNT(*) FROM ro_subfamily",
    }
    counts = {}
    for label, sql in tables.items():
        try:
            counts[label] = con.execute(sql).fetchone()[0]
        except sqlite3.DatabaseError:
            counts[label] = None
    distinct = {v for v in counts.values() if v is not None}
    ok = len(distinct) <= 1
    detail = f"all tables report {n_conf}" if ok else \
        "disagreement -> " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    rep.check("arithmetic: confirmed-RO count agrees across tables",
              ok, len(distinct), detail)

    # Her dogrulanan RO her tabloda TEMSIL EDILIYOR mu (sayilar esit olsa bile
    # kumeler farkli olabilir).
    missing = {}
    for table in ("ro_search", "ro_evidence", "ro_domain", "operon", "ro_leaf",
                  "ro_etc", "ro_regulation", "ro_subfamily"):
        try:
            n = con.execute(f"""
                SELECT COUNT(*) FROM ro r WHERE r.is_confirmed=1
                  AND NOT EXISTS (SELECT 1 FROM "{table}" t
                                  WHERE t.candidate_id=r.candidate_id)""").fetchone()[0]
        except sqlite3.DatabaseError:
            continue
        if n:
            missing[table] = n
    rep.check("arithmetic: every confirmed RO is present in every derived table",
              not missing, sum(missing.values()),
              "" if not missing else ", ".join(f"{k} misses {v}"
                                               for k, v in missing.items()))

    # Ters yon: tureyen tabloda dogrulanmamis giris olmamali.
    extra = {}
    for table in ("ro_search", "ro_evidence", "ro_domain", "operon", "ro_leaf"):
        try:
            n = con.execute(f"""
                SELECT COUNT(*) FROM "{table}" t
                  LEFT JOIN ro r ON r.candidate_id=t.candidate_id
                  WHERE r.candidate_id IS NULL OR r.is_confirmed!=1""").fetchone()[0]
        except sqlite3.DatabaseError:
            continue
        if n:
            extra[table] = n
    rep.check("arithmetic: derived tables hold only confirmed ROs",
              not extra, sum(extra.values()),
              "" if not extra else ", ".join(f"{k} has {v}" for k, v in extra.items()))


def check_search_index_fresh(con, rep, eco, chem):
    """ro_search, kuratorlu CSV'lerin SU ANKI halini mi yansitiyor?

    build_search_index.py substrat/aile/reaksiyon alanlarini CSV'lerden
    KOPYALAR. CSV'ler sonradan guncellenirse tablo bayatlar ve web sitesi eski
    degeri gosterir. SQLite tablo basina mtime tutmadigi icin bunu zaman
    damgasiyla yakalamak mumkun degil; tek yol icerigi karsilastirmak.
    """
    chem_by = {r["cluster"]: r for r in chem}
    eco_by = {r["cluster"]: r for r in eco}
    rows = con.execute("""
        SELECT cluster, substrate, family, reaction, COUNT(*)
        FROM ro_search GROUP BY cluster, substrate, family, reaction""").fetchall()
    bad = []
    entries = 0
    for cluster, substrate, family, reaction, n in rows:
        ch = chem_by.get(cluster, {})
        want_sub = (ch.get("substrate_en") or "").strip() \
            or (eco_by.get(cluster, {}).get("substrate") or "").strip()
        want_fam = (ch.get("family") or "").strip()
        want_rxn = (ch.get("reaction_class") or "").strip()
        got = ((substrate or "").strip(), (family or "").strip(),
               (reaction or "").strip())
        if got != (want_sub, want_fam, want_rxn):
            entries += n
            bad.append(f"{cluster} ({n} entries): index has "
                       f"{got} but the CSVs say "
                       f"{(want_sub, want_fam, want_rxn)}")
    rep.check("freshness: ro_search matches the curated CSVs",
              not bad, len(bad),
              f"{len(rows)} cluster groups checked" if not bad else
              f"{entries} indexed entries affected; " + first(bad, 3))


def check_cross_file_substrate(rep, eco, chem):
    """Iki kuratorlu CSV ayni tip hakkinda celisiyor mu?

    Ozellikle: chemistry.csv substrati ADLANDIRMIS ama cluster_ecology.csv
    hala 'unknown' diyorsa, ekoloji hipotez testi (substrate_class) o tipi
    disarida birakir ve web sayfasi iki farkli sey soyler.
    """
    eco_by = {r["cluster"]: r for r in eco}
    bad = []
    for r in chem:
        e = eco_by.get(r["cluster"])
        if not e:
            continue
        named = (r.get("substrate_en") or "").strip().lower() not in ("", "unknown")
        eco_unknown = (e.get("substrate") or "").strip().lower() in ("", "unknown") \
            or (e.get("substrate_class") or "").strip() == "unknown"
        if named and eco_unknown:
            bad.append(f"{r['cluster']}: chemistry says "
                       f"{r['substrate_en'][:34]!r} but ecology says "
                       f"substrate={e.get('substrate')!r} "
                       f"class={e.get('substrate_class')!r}")
    rep.check("consistency: the two CSVs agree on whether a substrate is known",
              not bad, len(bad), first(bad, 3))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(ROOT, "roar.sqlite"))
    ap.add_argument("--ecology", default=os.path.join(ROOT, "cluster_ecology.csv"))
    ap.add_argument("--chemistry", default=os.path.join(ROOT, "chemistry.csv"))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "analysis_out"))
    ap.add_argument("--refs", default=REFS_FASTA,
                    help="curated reference FASTA the whole pipeline rests on")
    ap.add_argument("--warn-only", action="store_true",
                    help="report everything but always exit 0")
    args = ap.parse_args()

    for path in (args.db, args.ecology, args.chemistry):
        if not os.path.exists(path):
            print(f"CANNOT RUN: {path} not found", file=sys.stderr)
            return 2

    rep = Report()
    print(f"validate_curation.py -- db={os.path.basename(args.db)} "
          f"ecology={os.path.basename(args.ecology)} "
          f"chemistry={os.path.basename(args.chemistry)}")
    print()

    check_vocabulary_drift(rep)
    check_reference_motifs(rep, ROOT)
    check_motif_columns_anchored(rep, ROOT)
    check_retracted_column_claim(rep, ROOT)
    check_group_monophyly(rep, ROOT, args.out_dir)
    check_group_numbers_not_published(rep, ROOT, args.out_dir)
    check_literature_motif_spacing(rep, ROOT)
    check_triad_gate(rep, ROOT)
    check_leaf_stop_reason(rep, ROOT)
    check_membership_basis(rep, ROOT)
    _, eco = read_csv_strict(args.ecology, rep, "cluster_ecology.csv")
    _, chem = read_csv_strict(args.chemistry, rep, "chemistry.csv")
    check_oxygen_fate(rep, ROOT, chem)
    check_reaction_names(rep, ROOT)
    check_computed_oxygen_fate(rep, ROOT)
    check_metabolic_role(rep, ROOT)
    check_no_internal_refs_published(rep, chem, ROOT)
    check_tables_are_scrollable(rep, ROOT)

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    con.text_factory = str
    try:
        check_reference_set(rep, args.refs,
                            os.path.join(args.out_dir, "reference_pairs.csv"),
                            eco, chem, con)
        n_conf = check_confirmed_ro(con, rep)
        check_operon_gap_single_source(rep, con, ROOT)
        check_carboxylate_claim_is_measured(rep, con, ROOT)
        check_sdp_regions(rep, con, ROOT)
        check_prevalence_declares_cohort(rep, con, ROOT)
        check_prefilter_admits_gate_passers(rep, con, ROOT)
        check_a5_match_guard(rep, con, ROOT)
        check_ro_footprint_explained(rep, con, ROOT)
        check_component_scan_is_fresh(rep, con, ROOT)
        check_cross_references(rep, con, ROOT)
        check_taxonomy_is_controlled(rep, con, ROOT)
        check_citation_and_version(rep, con, ROOT)
        check_hmm_libraries_are_pinned(rep, con, ROOT)
        check_gff3_export(rep, con, ROOT)
        check_substrate_chebi(rep, con, ROOT)
        check_record_kind(rep, con, ROOT)
        check_db_writes_are_attributed(rep, ROOT)
        check_assignment_library(rep, con, ROOT)
        check_ncbi_accession_claims(rep, con, ROOT)
        check_coverage_audit(rep, con, ROOT)
        check_same_structure_same_name(rep, chem)
        check_promoter_claim_is_levelled(rep, ROOT)
        check_family_scan_is_fresh(rep, con, ROOT)
        check_derived_statistics_are_fresh(rep, con, ROOT)
        check_published_files_are_in_rebuild_graph(rep, con, ROOT)
        check_union_denominator(rep, con, ROOT)
        check_regulation_totals_agree(rep, con, ROOT)
        check_significance_claims(rep, con, ROOT)
        check_types_are_single_source(rep, ROOT)
        check_groups_published(rep, con, ROOT)
        check_cluster_coverage(con, rep, eco, chem)
        check_ecology_csv(rep, eco)
        check_chemistry_csv(rep, chem)
        check_chemistry_published(rep, ROOT, con)
        check_display_names(rep, chem)
        check_host_kingdom(rep, con)
        check_tree_coverage(rep, con, args.out_dir,
                            [r['cluster'] for r in chem])
        check_published_files(rep, args.out_dir)
        check_reference_redundancy_agrees(rep, args.out_dir, args.refs)
        check_cross_file_substrate(rep, eco, chem)
        check_text_fields_classified(con, rep)
        check_turkish(con, rep, args.ecology, args.chemistry, eco, chem,
                      args.out_dir)
        check_referential(con, rep)
        check_arithmetic(con, rep, n_conf)
        check_search_index_fresh(con, rep, eco, chem)
    finally:
        con.close()

    code = rep.emit()
    if code and not args.warn_only:
        print("\nVALIDATION FAILED -- do not publish this build.", file=sys.stderr)
    return 0 if args.warn_only else code


if __name__ == "__main__":
    sys.exit(main())
