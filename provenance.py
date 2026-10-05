#!/usr/bin/env python3
"""
Koken (provenance) manifestosu -- her tablo ve her analiz dosyasi nereden geliyor.

NEDEN: bu pipeline'da en pahali hatalar sessiz olanlar oldu (cd-hit hatasinin
yutulmasi, CSV'de gecersiz virgul kacisi, hizalama gurultusunu secen metrik).
Bunlarin ortak ozelligi, bir ciktinin URETILDIGI ANDAN sonra girdisinin
degismesi ve kimsenin fark etmemesi. Bir artefakt haritada YOKSA ya da onu
uretmesi gereken script artik o dosyayi hic yazmiyorsa, o artefakt sessizce
bayatlar ve yine de web sitesinde yayinlanir.

BU SCRIPT uc sey yapar:
  1. Her SQLite tablosu ve analysis_out/ icindeki her dosya icin ureten script,
     tukettigi girdiler, satir/kayit sayisi, degisiklik zamani ve kisa bir
     Ingilizce aciklama yazar.
  2. Her zinciri kaynak veriye kadar coz: combined_pfam.fasta, gbk_files/,
     ROs_71_Clean/ ve iki el-kuratorlu CSV.
  3. YUKSEK SESLE uyarir:
       - haritada olmayan tablo/dosya (bayatlamanin klasik girisi),
       - haritada olup diskte olmayan artefakt,
       - ureten script'in icinde cikti adinin HIC GECMEDIGI artefakt
         (yani kodda artik uretici yok -- dosya eski bir kosudan kalmistir),
       - girdisi kendisinden YENI olan artefakt (bayat zincir).

Kullanim:
    python3 provenance.py                       # analysis_out/provenance.json yaz
    python3 provenance.py --strict               # uyari varsa sifir-disi cik
    python3 provenance.py --quiet                # sadece uyarilari yaz

SINIR (bilinmesi gereken): SQLite tablo basina degisiklik zamani TUTMAZ. Tablo
satirlarinin mtime'i yoktur; burada raporlanan zaman roar.sqlite DOSYASININ
zamanidir. Yani "su tablo girdisinden eski" tespiti dosyalar icin calisir,
tablolar icin CALISMAZ. Tablolarin icerik tutarliligini validate_curation.py
denetler (ornegin ro_search ile chemistry.csv'nin ayni seyi soyleyip
soylemedigi).
"""

import argparse
import csv
import datetime as dt
import fnmatch
import json
import os
import re
import sqlite3
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))

# --- Kaynak veri: zincirlerin en altinda duran, pipeline'in URETMEDIGI girdiler.
SOURCE_DATA = {
    "combined_pfam.fasta": {
        "kind": "file",
        "description": "All PF00355 (Rieske) proteins downloaded from Pfam; "
                       "the raw protein input of the whole pipeline.",
        "hand_curated": False,
    },
    "gbk_files/": {
        "kind": "dir",
        "description": "Downloaded GenBank records; the source of every genomic "
                       "coordinate, neighbour gene and neighbour protein sequence.",
        "hand_curated": False,
    },
    "ROs_71_Clean/": {
        "kind": "dir",
        "description": "The 71 hand-curated reference RO alpha subunits with their "
                       "HMM models and motif columns; defines type assignment, the "
                       "catalytic motif test and every evidence tier.",
        "hand_curated": True,
    },
    "cluster_ecology.csv": {
        "kind": "file",
        "description": "Hand-curated substrate, substrate class, ecology note and "
                       "confidence per reference type.",
        "hand_curated": True,
    },
    "chemistry.csv": {
        "kind": "file",
        "description": "Hand-curated substrate name, SMILES, product, reaction, "
                       "reaction class, chemical family, PDB entry and literature "
                       "source per reference type.",
        "hand_curated": True,
    },
    "RieskeDB71.hmm": {
        "kind": "file",
        "description": "Pressed HMM library of the 71 clean references, used for "
                       "type assignment by bit score.",
        "hand_curated": True,
    },
    "component_hmms/": {
        "kind": "dir",
        "description": "Pfam models for operon components (beta, ferredoxin, "
                       "reductase); input of operon derivation.",
        "hand_curated": False,
    },
}

# Pipeline'in kendi urettigi ara dizinler: kaynak degil, ama girdi olarak anilir.
INTERMEDIATE = {
    "ro_alpha.csv": ("run_ro_filter.py", ["combined_pfam.fasta", "RieskeDB71.hmm",
                                          "ROs_71_Clean/"]),
    "ro_alpha.fasta": ("run_ro_filter.py", ["combined_pfam.fasta", "RieskeDB71.hmm",
                                            "ROs_71_Clean/"]),
    "ro_not_alpha.csv": ("run_ro_filter.py", ["combined_pfam.fasta", "RieskeDB71.hmm"]),
    "ro_fragment.csv": ("run_ro_filter.py", ["combined_pfam.fasta", "RieskeDB71.hmm"]),
    "genomic_context/": ("extract_genomic_context.py", ["gbk_files/", "ro_alpha.csv"]),
}


def P(script, inputs, description, producer_pattern=None):
    """Tek bir artefaktin koken kaydi."""
    return {
        "script": script,
        "inputs": list(inputs),
        "description": description,
        # Ureticinin dogrulanmasi icin script icinde aranacak dizgi
        # (varsayilan: artefaktin kendi adi).
        "producer_pattern": producer_pattern,
    }


# --- YETKILI HARITA: tablolar -------------------------------------------------
# Script adlari README_PIPELINE.md'deki pipeline tablosundan turetilebilir, ama
# yetkili olan BU dict'tir; README ile celisirse uyari verilir.
TABLE_PROVENANCE = {
    "replicon": P("build_db.py", ["genomic_context/"],
                  "One row per GenBank replicon: organism, taxonomy, length, "
                  "plasmid and circularity flags, CDS count and whether the "
                  "record carried annotation at all."),
    "ro": P("build_db.py", ["genomic_context/"],
            "One row per RO alpha candidate with its genomic coordinates and "
            "strand. annotate_ro.py later fills the HMM verification columns "
            "(ro_cluster, model_coverage, rieske_intact, is_confirmed) and "
            "build_operons.py fills ro.sequence."),
    "neighbor": P("build_db.py", ["genomic_context/"],
                  "Neighbouring CDS of each RO candidate with signed distance, "
                  "same-strand flag, gene offset and origin-spanning flag."),
    "gene_category": P("build_db.py", ["genomic_context/", "table:neighbor"],
                       "Functional category assigned to each neighbour by "
                       "annotation-text rules (regulator, transposon, ro_beta, "
                       "ring_cleavage, ferredoxin and so on)."),
    "ro_domain": P("classify_domains.py", ["table:ro", "table:replicon"],
                   "Domain of life per confirmed RO (Bacteria, Eukaryota, "
                   "Archaea) with a eukaryotic subgroup where applicable."),
    "neighbor_protein": P("build_operons.py", ["gbk_files/", "table:neighbor"],
                          "Protein sequence of each neighbouring CDS, pulled back "
                          "out of the GenBank records and keyed by coordinates."),
    "neighbor_domain": P("build_operons.py", ["component_hmms/",
                                              "table:neighbor_protein"],
                         "Pfam component-model hits per neighbour protein with "
                         "score, E-value and coverage."),
    "neighbor_component": P("build_operons.py", ["table:neighbor_domain"],
                            "Component label per neighbour protein: beta, "
                            "ferredoxin, reductase, alpha_other, rieske_other "
                            "or none."),
    "operon": P("build_operons.py", ["table:neighbor_component", "table:ro"],
                "Derived operon per confirmed RO: extent, gene count, which "
                "components sit inside the operon versus merely within 10 kb, "
                "and a readable 5'->3' layout string."),
    "operon_gene": P("build_operons.py", ["table:operon", "table:neighbor"],
                     "One row per gene of each derived operon with its position "
                     "relative to the RO, component label and intergenic gap."),
    "ro_evidence": P("evidence_tiers.py", ["ROs_71_Clean/", "table:ro"],
                     "Evidence tier per confirmed RO from diamond identity to the "
                     "nearest curated reference: characterized, close_homolog, "
                     "family_member, distant or novel."),
    "ro_etc": P("etc_types.py", ["table:neighbor_component", "table:operon"],
                "Electron transport chain composition per RO: reductase type "
                "(FNR / GR), ferredoxin type (Rieske / plant) and whether a beta "
                "subunit is present."),
    "ro_regulation": P("analyze_regulation.py", ["table:operon", "table:neighbor",
                                                 "table:gene_category"],
                       "Regulatory architecture per RO: operon 5' end, the "
                       "intergenic region up to the first upstream gene, that "
                       "gene's orientation and regulator family."),
    "subfamily": P("discover_subfamilies.py", ["table:ro"],
                   "CD-HIT subfamilies within each type, with size, "
                   "representative and top genera."),
    "ro_subfamily": P("discover_subfamilies.py", ["table:ro", "ROs_71_Clean/"],
                      "Subfamily membership per RO plus identity to the best "
                      "matching curated reference and an assignment class "
                      "(core, divergent, alt_type, novel_candidate)."),
    "leaf": P("recursive_homogenize.py", ["table:ro", "table:ro_subfamily"],
              "Homogeneous leaves produced by recursively splitting "
              "heterogeneous types; a leaf is the variant level of this "
              "database."),
    "ro_leaf": P("recursive_homogenize.py", ["table:leaf"],
                 "Leaf membership per confirmed RO."),
    "cluster_sdp": P("variant_signature.py", ["table:ro", "table:leaf"],
                     "Alignment columns that separate the variants of a type, "
                     "with the gating statistics that qualified each column."),
    "leaf_sdp": P("variant_signature.py", ["table:cluster_sdp", "table:leaf"],
                  "Residue signature of each leaf at its type's discriminating "
                  "columns, and how far it departs from the type consensus."),
    "leaf_profile": P("characterize_leaves.py",
                      ["table:leaf", "table:neighbor", "table:gene_category"],
                      "Variant profile per leaf: median identity, genus spread, "
                      "plasmid and transposon rates, the distinctive neighbour "
                      "products and a human-readable label."),
    "ro_search": P("build_search_index.py",
                   ["table:ro", "table:ro_leaf", "table:ro_evidence",
                    "table:ro_domain", "table:operon", "table:ro_regulation",
                    "cluster_ecology.csv", "chemistry.csv"],
                   "One search document per confirmed RO plus the filter fields "
                   "used by the web interface (evidence tier, domain, chemical "
                   "family, reaction, plasmid, operon partner, regulator)."),
    "ro_fts": P("build_search_index.py", ["table:ro_search"],
                "FTS5 full-text index over ro_search.doc, ranked with bm25."),
    # Bu tablo README_PIPELINE.md'nin "Ciktilar" listesinde YOKTU; haritaya
    # eklenmesi gerekti (olculdu: motif_stats.py kuruyor, stats_overview.py
    # okuyor).
    "ro_carboxylate": P("motif_stats.py", ["table:ro", "ROs_71_Clean/"],
                        "Per-entry identity of the Fe(II) carboxylate and the "
                        "bridging acid residue at the catalytic columns, with "
                        "the alignment offset at which each was found; the "
                        "evidence behind the statement that the carboxylate "
                        "position accepts both Asp and Glu."),
    "replicon_source": P("isolation_source.py", ["gbk_files/"],
                         "Isolation metadata per replicon read from the GenBank "
                         "source feature: free-text isolation source, host, "
                         "geography, country, collection year, and the habitat "
                         "label assigned by the ordered keyword rules."),
}

# FTS5'in kendi golge tablolari ve SQLite'in ic tablolari: ayri ele alinir,
# "haritada yok" uyarisi uretmemeleri gerekir.
SHADOW_TABLES = {
    "ro_fts_data": "FTS5 internal storage for ro_fts.",
    "ro_fts_idx": "FTS5 internal term index for ro_fts.",
    "ro_fts_docsize": "FTS5 internal document size table for ro_fts.",
    "ro_fts_config": "FTS5 internal configuration table for ro_fts.",
    "sqlite_sequence": "SQLite internal AUTOINCREMENT bookkeeping.",
    "sqlite_stat1": "SQLite internal query planner statistics (ANALYZE).",
}

# --- YETKILI HARITA: analysis_out/ dosyalari ---------------------------------
# Anahtar analysis_out/ icine goreli yoldur; glob kalibi da olabilir.
FILE_PROVENANCE = {
    "null_model.csv": P("null_model.py", ["table:ro", "table:neighbor", "gbk_files/"],
                        "Neighbourhood enrichment of each gene category against "
                        "random same-size windows drawn from the same replicon "
                        "population."),
    "transposon_by_genus.csv": P("analyze.py", ["table:ro", "table:gene_category",
                                                "table:replicon"],
                                 "Transposon adjacency per genus crossed with RO "
                                 "type, used for the mobility analysis."),
    "cluster_variance.csv": P("analyze_variants.py", ["table:ro"],
                              "Within-type identity spread per RO type: median, "
                              "quartiles and the fraction of members below the "
                              "homogeneity target."),
    "sdp_positions.csv": P("analyze_variants.py", ["table:ro"],
                           "Candidate specificity-determining alignment positions "
                           "per type from the first-pass variant analysis."),
    "subfamilies.csv": P("discover_subfamilies.py", ["table:subfamily"],
                         "Flat export of the CD-HIT subfamily table."),
    "novel_candidates.csv": P("discover_subfamilies.py", ["table:ro_subfamily"],
                              "RO members whose identity to every curated "
                              "reference falls below the novel-type threshold."),
    "leaves.csv": P("recursive_homogenize.py", ["table:leaf"],
                    "Flat export of the leaf table (one row per variant)."),
    "leaf_profiles.csv": P("characterize_leaves.py", ["table:leaf_profile"],
                           "Flat export of the per-leaf variant profiles."),
    "representatives.csv": P("extract_representatives.py",
                             ["table:leaf", "table:ro", "cluster_ecology.csv"],
                             "One representative sequence per leaf and per type "
                             "with its catalytic-site verification result."),
    "representatives_by_leaf.fasta": P("extract_representatives.py", ["table:leaf",
                                                                      "table:ro"],
                                       "Representative protein sequences, one per "
                                       "leaf; the input of the phylogeny."),
    "representatives_by_cluster.fasta": P("extract_representatives.py", ["table:ro"],
                                          "Representative protein sequences, one "
                                          "per RO type."),
    "domain_by_cluster.csv": P("classify_domains.py", ["table:ro_domain"],
                               "Domain-of-life composition per RO type."),
    "cluster_ecology_stats.csv": P("analyze_ecology.py",
                                   ["table:ro", "table:gene_category",
                                    "table:replicon", "table:operon",
                                    "cluster_ecology.csv"],
                                   "Per-type ecology and mobility table: "
                                   "transposon and plasmid rates, genus count, "
                                   "synteny entropy and operon completeness, "
                                   "joined to the curated substrate class."),
    "etc_by_cluster.csv": P("etc_types.py", ["table:ro_etc"],
                            "Electron transport chain composition per RO type."),
    "evidence_by_cluster.csv": P("evidence_tiers.py", ["table:ro_evidence"],
                                 "Evidence tier distribution per RO type."),
    # Bu dosya ONCE pipeline disinda ELLE uretiliyordu; artik evidence_tiers.py
    # yaziyor (diamond all-vs-all) ve substrat adlarini chemistry.csv'nin
    # INGILIZCE alanlarindan aliyor. Elle uretilen surum Turkce etiketler
    # tasiyordu ve web'de yayinlaniyordu.
    "reference_pairs.csv": P("evidence_tiers.py",
                             ["ROs_71_Clean/", "chemistry.csv"],
                             "All-versus-all identity between the 71 curated "
                             "references with each pair's substrates; the "
                             "calibration behind the claim that no global "
                             "identity threshold implies a shared substrate."),
    "regulation_by_cluster.csv": P("analyze_regulation.py", ["table:ro_regulation"],
                                   "Regulatory architecture distribution per RO "
                                   "type."),
    "variant_signatures.csv": P("variant_signature.py", ["table:cluster_sdp",
                                                         "table:leaf_sdp"],
                                "Residue signature of every leaf at the columns "
                                "that discriminate the variants of its type."),
    "motif_stats.json": P("motif_stats.py", ["table:ro", "ROs_71_Clean/"],
                          "Conservation of the eight defining alignment columns "
                          "(Rieske ligands, bridging acid, catalytic triad) "
                          "across all confirmed ROs."),
    "operon_validation.json": P("operon_validation.py",
                                ["table:operon_gene", "table:neighbor_component",
                                 "table:neighbor"],
                                "Test of the operon convention against a "
                                "background of every gene in the same windows, "
                                "with a negative control and a 50-500 bp "
                                "threshold sensitivity sweep."),
    "redundancy.json": P("redundancy.py", ["table:ro"],
                         "Sequence redundancy: unique sequences versus entries, "
                         "copy-number distribution and the effect of "
                         "deduplication on every count."),
    "cooccurrence.json": P("cooccurrence.py", ["table:ro", "table:replicon"],
                           "Co-occurrence of RO types within 10 kb and on the "
                           "same replicon, against a permutation null that holds "
                           "per-replicon counts and type totals fixed."),
    "substrate_predictability.json": P("substrate_predictability.py",
                                       ["analysis_out/reference_pairs.csv",
                                        "cluster_ecology.csv"],
                                       "ROC and precision of the rule 'above X % "
                                       "identity means the same substrate', "
                                       "measured over all curated reference pairs."),
    "stats.json": P("stats_overview.py",
                    ["table:ro", "table:replicon", "table:gene_category",
                     "table:operon", "analysis_out/domain_by_cluster.csv",
                     "cluster_ecology.csv"],
                    "Hypothesis tests behind the headline findings, each run at "
                    "entry, unique-sequence and type-by-genus level with effect "
                    "sizes."),
    "habitat.json": P("isolation_source.py",
                      ["table:replicon_source", "table:ro", "table:replicon",
                       "cluster_ecology.csv", "chemistry.csv"],
                      "Habitat composition of the database and habitat crossed "
                      "with substrate class, reported both per entry and per "
                      "distinct organism name, plus which keyword decided each "
                      "habitat assignment."),
    "tree_all.nwk": P("build_phylogeny.py",
                      ["analysis_out/representatives_by_leaf.fasta",
                       "ROs_71_Clean/"],
                      "FastTree phylogeny of the 71 references plus one "
                      "representative per leaf, built on the hmmalign-trimmed "
                      "catalytic core."),
    "trees/*.nwk": P("build_phylogeny.py",
                     ["analysis_out/representatives_by_leaf.fasta"],
                     "Per-type phylogeny of that type's leaf representatives.",
                     producer_pattern="trees"),
    "ssn_edges.csv": P("build_phylogeny.py",
                       ["analysis_out/representatives_by_leaf.fasta"],
                       "Sequence similarity network edges from diamond "
                       "all-versus-all above the identity cutoff."),
    "ssn_nodes.csv": P("build_phylogeny.py",
                       ["analysis_out/representatives_by_leaf.fasta"],
                       "Sequence similarity network nodes with their type and "
                       "leaf annotation."),
    "cluster_identity_matrix.csv": P("build_phylogeny.py",
                                     ["analysis_out/representatives_by_leaf.fasta"],
                                     "Highest representative-to-representative "
                                     "identity for every pair of RO types."),
    # Bu iki FASTA'yi HICBIR mevcut script yazmiyor (olculdu). Haritada
    # tutuluyorlar ki uretici dogrulamasi onlari her kosuda isaretlesin;
    # Her ikisi de ARTIK evidence_tiers.py tarafindan uretiliyor. Onceden
    # hicbir script yazmiyordu ve make_hub.py bu dosyadan okudugu "155
    # yuksek-guven novel aday" sayisini yayinliyordu; yani yayinlanan sayi
    # mevcut kodla yeniden uretilemiyordu. Kural artik evidence_tiers.py
    # icinde acikca tanimli (dogrulanmis + merkezler tam + >=300 aa +
    # <%25 kimlik + >=3 uyeli varyant) ve 318 aday veriyor.
    "novel_high_confidence.fasta": P("evidence_tiers.py",
                                     ["table:ro_evidence", "table:ro_leaf", "table:leaf"],
                                     "Protein sequences of the high-confidence novel "
                                     "candidates, written under the rule documented in "
                                     "evidence_tiers.py."),
    "novel_high_confidence.csv": P("evidence_tiers.py",
                                   ["table:ro_evidence", "table:ro_leaf", "table:leaf"],
                                   "The same candidates as a table: nearest reference, "
                                   "identity to it, variant and organism."),
    "provenance.json": P("provenance.py", ["roar.sqlite", "analysis_out/"],
                         "This manifest: the origin, inputs, size and age of "
                         "every table and analysis file."),
}

# Pipeline tablosunda olup analysis_out/ icine artefakt yazmayan script'ler.
NON_ARTEFACT_SCRIPTS = {
    "run_ro_filter.py",        # ro_alpha.csv (kok dizin)
    "extract_genomic_context.py",  # genomic_context/
    "build_db.py", "annotate_ro.py",   # sadece tablolar
    "build_search_index.py",
    "make_report.py", "make_explorer.py", "make_hub.py",   # HTML
    "validate_curation.py",    # kapi; artefakt yazmaz
    "provenance.py",           # manifestoyu kendisi yazar, kendini sayma
}


# ---------------------------------------------------------------------------
# Yardimcilar
# ---------------------------------------------------------------------------

def iso(ts):
    return dt.datetime.fromtimestamp(ts).replace(microsecond=0).isoformat()


def newest_mtime(path):
    """Dosya icin mtime; dizin icin icindeki en yeni dosyanin mtime'i."""
    if not os.path.exists(path):
        return None
    if os.path.isfile(path):
        return os.path.getmtime(path)
    newest = os.path.getmtime(path)
    for root, _dirs, files in os.walk(path):
        for fn in files:
            try:
                newest = max(newest, os.path.getmtime(os.path.join(root, fn)))
            except OSError:
                pass
    return newest


def count_dir(path):
    n = 0
    for _root, _dirs, files in os.walk(path):
        n += len(files)
    return n


def count_records(path):
    """Dosya tipine gore kayit sayisi + sayilan seyin adi."""
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext == ".csv":
            with open(path, newline="", encoding="utf-8", errors="replace") as fh:
                rdr = csv.reader(fh)
                try:
                    next(rdr)
                except StopIteration:
                    return 0, "data rows (file is empty)"
                return sum(1 for row in rdr if row), "data rows"
        if ext == ".json":
            with open(path, encoding="utf-8") as fh:
                obj = json.load(fh)
            if isinstance(obj, list):
                return len(obj), "list items"
            if isinstance(obj, dict):
                return len(obj), "top-level keys"
            return 1, "scalar value"
        if ext in (".fasta", ".fa", ".faa"):
            with open(path, encoding="utf-8", errors="replace") as fh:
                return sum(1 for line in fh if line.startswith(">")), "sequences"
        if ext == ".nwk":
            text = open(path, encoding="utf-8", errors="replace").read()
            return len(re.findall(r"[(,]\s*([^(),:;]+):", text)), "tree tips"
    except Exception as exc:                      # noqa: BLE001
        return None, f"unreadable ({exc.__class__.__name__})"
    return None, "not counted (unknown file type)"


def read_script(name):
    path = os.path.join(ROOT, name)
    if not os.path.exists(path):
        return None
    return open(path, encoding="utf-8", errors="replace").read()


def code_only(text):
    """Docstring ve yorumlari atarak sadece KODU dondur.

    Neden gerekli: evidence_tiers.py'nin docstring'i
    'analysis_out/reference_pairs.csv uretir' diyor ama kodda o dosyayi yazan
    tek satir yok. Duz metin aramasi bu yalani dogrulamis gibi gorunur; o
    yuzden uretici dogrulamasi yalnizca kod govdesine bakar.
    """
    # Uc tirnakli bloklar (docstring ve coklu satir dizgiler)
    text = re.sub(r'"""(?:.|\n)*?"""', '""', text)
    text = re.sub(r"'''(?:.|\n)*?'''", "''", text)
    # Satir sonu yorumlari (dizgi icindeki '#' nadir, kabul edilebilir kayip)
    text = re.sub(r"(?m)#.*$", "", text)
    return text


def parse_readme_pipeline(path):
    """README_PIPELINE.md'deki pipeline tablosundan adim -> script eslemesi."""
    steps = {}
    if not os.path.exists(path):
        return steps
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 3:
            continue
        step = cells[0]
        if not re.match(r"^\d+[a-z0-9]*$", step):
            continue
        scripts = re.findall(r"`([A-Za-z0-9_]+\.py)`", cells[1])
        if not scripts:
            continue
        steps[step] = {"scripts": scripts, "io": cells[2]}
    return steps


# ---------------------------------------------------------------------------
# Manifesto kurulumu
# ---------------------------------------------------------------------------

class Manifest:
    def __init__(self, db_path, out_dir):
        self.db_path = db_path
        self.out_dir = out_dir
        self.warnings = []
        self.mtime_cache = {}

    def warn(self, kind, target, message):
        self.warnings.append({"kind": kind, "target": target, "message": message})

    # --- girdi adlarini diskteki yollara cevirme ---------------------------
    def resolve(self, name):
        """Girdi adi -> (concrete_path, is_source, source_keys).

        Girdi adi su bicimlerden biri olabilir:
            "table:X"              -> roar.sqlite (tablo basina mtime yok)
            SOURCE_DATA anahtari   -> kaynak veri
            "analysis_out/x.csv"   -> baska bir analiz dosyasi
            baska herhangi bir yol -> ara urun / kok dizindeki dosya
        """
        if name.startswith("table:"):
            return os.path.join(ROOT, os.path.basename(self.db_path)), False
        if name in SOURCE_DATA:
            return os.path.join(ROOT, name.rstrip("/")), True
        return os.path.join(ROOT, name.rstrip("/")), False

    def source_chain(self, inputs, seen=None):
        """Girdileri ozyinelemeli cozerek dayandigi KAYNAK verileri topla."""
        if seen is None:
            seen = set()
        sources = set()
        for name in inputs:
            if name in seen:
                continue
            seen.add(name)
            if name in SOURCE_DATA:
                sources.add(name)
                continue
            if name.startswith("table:"):
                tbl = name.split(":", 1)[1]
                prov = TABLE_PROVENANCE.get(tbl)
                if prov:
                    sources |= self.source_chain(prov["inputs"], seen)
                continue
            if name.startswith("analysis_out/"):
                rel = name[len("analysis_out/"):]
                prov = self.lookup_file(rel)
                if prov:
                    sources |= self.source_chain(prov["inputs"], seen)
                continue
            if name in INTERMEDIATE:
                sources |= self.source_chain(INTERMEDIATE[name][1], seen)
                continue
            # roar.sqlite gibi toplu girdiler: tum tablolarin zincirini al
            if os.path.basename(name) == os.path.basename(self.db_path):
                for prov in TABLE_PROVENANCE.values():
                    sources |= self.source_chain(prov["inputs"], seen)
                continue
            if name.rstrip("/") == "analysis_out":
                for prov in FILE_PROVENANCE.values():
                    sources |= self.source_chain(prov["inputs"], seen)
                continue
        return sources

    @staticmethod
    def lookup_file(rel):
        """Tam ad once, sonra glob kalibi."""
        if rel in FILE_PROVENANCE:
            return FILE_PROVENANCE[rel]
        for pattern, prov in FILE_PROVENANCE.items():
            if ("*" in pattern or "?" in pattern) and fnmatch.fnmatch(rel, pattern):
                return prov
        return None

    def mtime(self, path):
        if path not in self.mtime_cache:
            self.mtime_cache[path] = newest_mtime(path)
        return self.mtime_cache[path]

    def stale_inputs(self, target, artefact_mtime, inputs):
        """Kendisinden YENI olan girdileri listele (bayat zincir)."""
        if artefact_mtime is None or target == "analysis_out/provenance.json":
            return []
        stale = []
        for name in inputs:
            path, _is_source = self.resolve(name)
            if name.startswith("table:"):
                # Tablo basina mtime yok; DB dosyasinin zamani yaniltici olur.
                continue
            mt = self.mtime(path)
            if mt is None:
                self.warn("missing_input", target,
                          f"input '{name}' does not exist at {path}")
                continue
            if mt > artefact_mtime + 1:
                stale.append({"input": name, "input_mtime": iso(mt),
                              "newer_by_seconds": int(mt - artefact_mtime)})
        return stale

    def verify_producer(self, target, prov, basename):
        """Ureten script'in icinde cikti adi gerceten geciyor mu?"""
        script = prov.get("script")
        if script is None:
            self.warn("no_producer", target,
                      "no script in the current code base writes this artefact; "
                      "it is a left-over from an earlier run and will silently "
                      "go stale")
            return "none"
        text = read_script(script)
        if text is None:
            self.warn("missing_script", target,
                      f"producing script '{script}' not found in the repository")
            return "script_missing"
        needle = prov.get("producer_pattern") or basename
        body = code_only(text)
        if needle not in body:
            where = "only in a comment or docstring" if needle in text \
                    else "nowhere in the file"
            self.warn("producer_does_not_write", target,
                      f"'{script}' is recorded as the producer but '{needle}' "
                      f"appears {where}; either the map is wrong or the "
                      f"artefact is orphaned and will silently go stale")
            return "claimed_only"
        return "verified"

    # --- tablolar ----------------------------------------------------------
    def build_tables(self):
        out = {}
        if not os.path.exists(self.db_path):
            self.warn("missing_db", self.db_path, "database file not found")
            return out
        db_mtime = self.mtime(self.db_path)
        con = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True)
        names = [r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type IN ('table','view') "
            "ORDER BY name")]
        for tbl in names:
            try:
                rows = con.execute(f'SELECT COUNT(*) FROM "{tbl}"').fetchone()[0]
            except sqlite3.DatabaseError as exc:
                rows = None
                self.warn("unreadable_table", tbl, str(exc))
            if tbl in SHADOW_TABLES:
                out[tbl] = {
                    "script": "build_search_index.py" if tbl.startswith("ro_fts")
                              else "sqlite3 (internal)",
                    "inputs": [], "rows": rows,
                    "mtime": iso(db_mtime) if db_mtime else None,
                    "description": SHADOW_TABLES[tbl],
                    "producer_status": "internal", "source_data": [],
                    "stale_inputs": [],
                }
                continue
            prov = TABLE_PROVENANCE.get(tbl)
            if prov is None:
                self.warn("unmapped_table", tbl,
                          "table is not in TABLE_PROVENANCE; an unmapped "
                          "artefact is exactly the thing that silently goes "
                          "stale -- add it to the map")
                out[tbl] = {"script": None, "inputs": [], "rows": rows,
                            "mtime": iso(db_mtime) if db_mtime else None,
                            "description": "UNMAPPED -- origin unknown.",
                            "producer_status": "unmapped", "source_data": [],
                            "stale_inputs": []}
                continue
            status = self.verify_producer(f"table:{tbl}", prov, tbl)
            out[tbl] = {
                "script": prov["script"],
                "inputs": prov["inputs"],
                "rows": rows,
                "mtime": iso(db_mtime) if db_mtime else None,
                "mtime_note": "mtime of roar.sqlite; SQLite keeps no per-table "
                              "timestamp, so table-level staleness cannot be "
                              "detected from it",
                "description": prov["description"],
                "producer_status": status,
                "source_data": sorted(self.source_chain(prov["inputs"])),
                "stale_inputs": self.stale_inputs(f"table:{tbl}", db_mtime,
                                                  prov["inputs"]),
            }
        # Haritada olup veritabaninda olmayan tablolar
        for tbl in TABLE_PROVENANCE:
            if tbl not in names:
                self.warn("missing_table", tbl,
                          "table is in TABLE_PROVENANCE but not in the database; "
                          "the step that builds it has not been run")
        con.close()
        return out

    # --- analysis_out/ dosyalari -------------------------------------------
    def build_files(self):
        out = {}
        if not os.path.isdir(self.out_dir):
            self.warn("missing_out_dir", self.out_dir, "directory not found")
            return out
        found = []
        for root, _dirs, files in os.walk(self.out_dir):
            for fn in sorted(files):
                full = os.path.join(root, fn)
                rel = os.path.relpath(full, self.out_dir)
                found.append(rel)
        for rel in sorted(found):
            full = os.path.join(self.out_dir, rel)
            mt = self.mtime(full)
            n, kind = count_records(full)
            prov = self.lookup_file(rel)
            if prov is None:
                self.warn("unmapped_file", f"analysis_out/{rel}",
                          "file is not in FILE_PROVENANCE; an unmapped artefact "
                          "is exactly the thing that silently goes stale -- add "
                          "it to the map or delete the file")
                out[rel] = {"script": None, "inputs": [], "records": n,
                            "record_kind": kind, "mtime": iso(mt) if mt else None,
                            "size_bytes": os.path.getsize(full),
                            "description": "UNMAPPED -- origin unknown.",
                            "producer_status": "unmapped", "source_data": [],
                            "stale_inputs": []}
                continue
            status = self.verify_producer(f"analysis_out/{rel}", prov,
                                          os.path.basename(rel))
            out[rel] = {
                "script": prov["script"],
                "inputs": prov["inputs"],
                "records": n,
                "record_kind": kind,
                "mtime": iso(mt) if mt else None,
                "size_bytes": os.path.getsize(full),
                "description": prov["description"],
                "producer_status": status,
                "source_data": sorted(self.source_chain(prov["inputs"])),
                "stale_inputs": self.stale_inputs(f"analysis_out/{rel}", mt,
                                                  prov["inputs"]),
            }
        # Haritada olup diskte olmayan dosyalar
        for pattern in FILE_PROVENANCE:
            if "*" in pattern or "?" in pattern:
                if not any(fnmatch.fnmatch(r, pattern) for r in found):
                    self.warn("missing_file", f"analysis_out/{pattern}",
                              "nothing on disk matches this mapped pattern; the "
                              "step that produces it has not been run")
            elif pattern not in found and pattern != "provenance.json":
                self.warn("missing_file", f"analysis_out/{pattern}",
                          "file is in FILE_PROVENANCE but not on disk; the step "
                          "that produces it has not been run")
        return out

    # --- kaynak veri -------------------------------------------------------
    def build_sources(self):
        out = {}
        for name, meta in SOURCE_DATA.items():
            path = os.path.join(ROOT, name.rstrip("/"))
            if not os.path.exists(path):
                self.warn("missing_source", name,
                          "source data is referenced by the pipeline but not "
                          "present")
                out[name] = dict(meta, present=False)
                continue
            mt = self.mtime(path)
            entry = dict(meta, present=True, mtime=iso(mt) if mt else None)
            if meta["kind"] == "dir":
                entry["files"] = count_dir(path)
            else:
                entry["size_bytes"] = os.path.getsize(path)
                n, kind = count_records(path)
                if n is not None:
                    entry["records"] = n
                    entry["record_kind"] = kind
            out[name] = entry
        return out

    # --- README capraz kontrolu -------------------------------------------
    def cross_check_readme(self, readme):
        steps = parse_readme_pipeline(readme)
        if not steps:
            self.warn("readme_unparsed", readme,
                      "could not parse the pipeline table; the producing script "
                      "of each artefact could not be cross-checked")
            return {}
        listed = {s for step in steps.values() for s in step["scripts"]}
        mapped = {p["script"] for p in TABLE_PROVENANCE.values() if p["script"]}
        mapped |= {p["script"] for p in FILE_PROVENANCE.values() if p["script"]}
        mapped.discard("provenance.py")
        for script in sorted(mapped - listed):
            self.warn("script_not_in_readme", script,
                      "this script produces mapped artefacts but does not appear "
                      "in the README_PIPELINE.md pipeline table")
        for script in sorted(listed - mapped - NON_ARTEFACT_SCRIPTS):
            self.warn("script_produces_nothing_mapped", script,
                      "listed as a pipeline step but no table or analysis file "
                      "in the map names it as producer")
        return steps


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(ROOT, "roar.sqlite"))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "analysis_out"))
    ap.add_argument("--readme", default=os.path.join(ROOT, "README_PIPELINE.md"))
    ap.add_argument("--out", default=None,
                    help="default: <out-dir>/provenance.json")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if there is any warning")
    ap.add_argument("--quiet", action="store_true", help="print warnings only")
    args = ap.parse_args()

    out_path = args.out or os.path.join(args.out_dir, "provenance.json")

    man = Manifest(args.db, args.out_dir)
    readme_steps = man.cross_check_readme(args.readme)
    sources = man.build_sources()
    tables = man.build_tables()
    files = man.build_files()

    manifest = {
        "generated": dt.datetime.now().replace(microsecond=0).isoformat(),
        "root": ROOT,
        "database": {
            "path": os.path.relpath(args.db, ROOT),
            "mtime": iso(man.mtime(args.db)) if os.path.exists(args.db) else None,
            "size_bytes": os.path.getsize(args.db) if os.path.exists(args.db) else None,
            "note": "SQLite keeps no per-table modification time; every table "
                    "below reports the modification time of this file.",
        },
        "source_data": sources,
        "tables": tables,
        "files": files,
        "readme_pipeline_steps": readme_steps,
        "warnings": man.warnings,
        "counts": {
            "tables": len(tables),
            "files": len(files),
            "warnings": len(man.warnings),
        },
    }

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=True, sort_keys=False)
        fh.write("\n")

    if not args.quiet:
        print(f"[provenance] {len(tables)} tables, {len(files)} analysis files, "
              f"{len(sources)} source data roots -> {out_path}")
        unmapped = [k for k, v in list(tables.items()) + list(files.items())
                    if v["producer_status"] == "unmapped"]
        stale = [k for k, v in list(tables.items()) + list(files.items())
                 if v["stale_inputs"]]
        print(f"[provenance] unmapped artefacts: {len(unmapped)} | "
              f"artefacts older than an input: {len(stale)}")

    if man.warnings:
        print("", file=sys.stderr)
        print("=" * 72, file=sys.stderr)
        print(f"PROVENANCE WARNINGS: {len(man.warnings)}", file=sys.stderr)
        print("=" * 72, file=sys.stderr)
        by_kind = {}
        for w in man.warnings:
            by_kind.setdefault(w["kind"], []).append(w)
        for kind in sorted(by_kind):
            print(f"\n-- {kind} ({len(by_kind[kind])})", file=sys.stderr)
            for w in by_kind[kind]:
                print(f"   {w['target']}: {w['message']}", file=sys.stderr)
        print("", file=sys.stderr)
        if args.strict:
            return 1
    elif not args.quiet:
        print("[provenance] no warnings")
    return 0


if __name__ == "__main__":
    sys.exit(main())
