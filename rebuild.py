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
import hashlib
import json
import os
import subprocess
import sys
import time

STATE = "analysis_out/.rebuild_state.json"
OUT = "analysis_out"

# Cekirdek girdiler. Bunlardan biri degisirse asagidaki her sey bayatlar.
DB = "roar.sqlite"
CHEM = "chemistry.csv"
ECO = "cluster_ecology.csv"
REFS = "ROs_71_Clean/refs71.fasta"
MOTIF_HMM = "ROs_71_Clean/ROmotif71.hmm"

# Her kayit: (ad, ciktilar, komut, girdiler, agir mi)
#
# Girdilere BETIGIN KENDISI de dahildir: analiz kodu degistiginde ciktisi
# bayattir, ve bu oturumda tam olarak bu oldu -- stats_overview.py degisti,
# stats.json elle yeniden uretilmeseydi sayfa eski sayilari gosterecekti.
ARTEFACTS = [
    # --- cekirdek pipeline (agir olanlar isaretli)
    ("protein_filter", ["ro_alpha.csv", "ro_alpha.fasta", "ro_domtbl.out"],
     ["python3", "run_ro_filter.py", "--fasta", "combined_pfam.fasta",
      "--hmm", "RieskeDB71.hmm"],
     ["run_ro_filter.py", "ro_filter.py", "RieskeDB71.hmm"], True),
    ("genomic_context", ["genomic_context/neighbors.csv",
                         "genomic_context/ro_candidates.fasta"],
     ["python3", "extract_genomic_context.py", "--gbk-dir", "gbk_files",
      "--out-dir", "genomic_context"],
     ["extract_genomic_context.py"], True),
    ("database", [DB],
     ["python3", "build_db.py", "--context-dir", "genomic_context", "--db", DB],
     ["build_db.py", "genomic_context/neighbors.csv"], True),
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
     ["build_operons.py", DB, "component_hmms/components.hmm"], True),

    # --- veritabanina yazan ara adimlar
    ("evidence_tiers", [f"{OUT}/evidence_by_cluster.csv"],
     ["python3", "evidence_tiers.py", "--db", DB, "--out-dir", OUT],
     ["evidence_tiers.py", DB, CHEM], False),
    ("etc_types", [f"{OUT}/etc_by_cluster.csv"],
     ["python3", "etc_types.py", "--db", DB, "--out-dir", OUT],
     ["etc_types.py", DB], False),
    ("regulation", [],
     ["python3", "analyze_regulation.py", "--db", DB, "--out-dir", OUT],
     ["analyze_regulation.py", DB], False),

    # --- run_all.sh 5-12: grafikte eksik olan analiz ve DB-yazar adimlari
    ("analyze", [],
     ["python3", "analyze.py", "--db", DB, "--out-dir", OUT],
     ["analyze.py", DB], False),
    ("null_model", [f"{OUT}/null_model.csv"],
     ["python3", "null_model.py", "--db", DB, "--gbk-dir", "gbk_files",
      "--processes", "16", "--out", f"{OUT}/null_model.csv"],
     ["null_model.py", DB], True),
    ("variants", [],
     ["python3", "analyze_variants.py", "--db", DB, "--out-dir", OUT],
     ["analyze_variants.py", DB], False),
    ("subfamilies", [],
     ["python3", "discover_subfamilies.py", "--db", DB, "--out-dir", OUT],
     ["discover_subfamilies.py", DB], False),
    ("homogenize", [],
     ["python3", "recursive_homogenize.py", "--db", DB, "--out-dir", OUT],
     ["recursive_homogenize.py", DB], False),
    ("leaves", [],
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
    ("variant_signature", [],
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
      "--chemistry", CHEM, "--out", f"{OUT}/reference_redundancy.json"],
     ["reference_redundancy.py", DB, CHEM, REFS], False),
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
     ["isolation_source.py", DB, ECO, CHEM], False),
    ("threshold_sensitivity", [f"{OUT}/threshold_sensitivity.json"],
     ["python3", "threshold_sensitivity.py", "--db", DB,
      "--out", f"{OUT}/threshold_sensitivity.json"],
     ["threshold_sensitivity.py", DB, ECO], False),
    ("stratified_stats", [f"{OUT}/stratified_stats.json"],
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

    # --- ozet katmani: HER SEYDEN sonra
    ("phylogeny", [],
     ["python3", "build_phylogeny.py", "--db", DB, "--out-dir", OUT,
      "--fasttree", "FastTree"],
     ["build_phylogeny.py", DB, MOTIF_HMM], False),
    ("stats", [f"{OUT}/stats.json"],
     ["python3", "stats_overview.py", "--db", DB, "--out", f"{OUT}/stats.json"],
     ["stats_overview.py", DB, CHEM, ECO], False),
    ("search_index", [],
     ["python3", "build_search_index.py", "--db", DB, "--ecology", ECO,
      "--chemistry", CHEM],
     ["build_search_index.py", DB, ECO, CHEM], False),
    ("validation_gate", [],
     ["python3", "validate_curation.py", "--db", DB, "--ecology", ECO,
      "--chemistry", CHEM, "--out-dir", OUT],
     ["validate_curation.py", DB, ECO, CHEM, "ROs_71_Clean/ROmotif71_hmmaln.sto"], False),
    ("provenance", [f"{OUT}/provenance.json"],
     ["python3", "provenance.py", "--db", DB, "--out-dir", OUT],
     ["provenance.py", DB], False),
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


def digest(path):
    """Icerik ozeti. Buyuk dosyalarda bas + son + boyut yeter ve hizlidir."""
    if path.endswith(".sqlite"):
        return db_counts(path)
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
    missing_in = [p for p in inputs if not os.path.exists(p)]
    if missing_in:
        return f"girdi yok: {', '.join(missing_in)}"
    if outputs:
        oldest_out = min(os.path.getmtime(p) for p in outputs)
        # Veritabani mtime karsilastirmasinin DISINDA; yukaridaki nota bak.
        newer = [p for p in inputs
                 if not p.endswith(".sqlite")
                 and os.path.getmtime(p) > oldest_out]
        if newer:
            return f"girdi daha yeni: {', '.join(sorted(newer)[:3])}"
    prev = (state.get(name) or {}).get("inputs") or {}
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
    # Her adimin yazdigi tablolar: statik (kaynak) + dinamik (olculen).
    writes = {}
    for name, outputs, cmd, inputs, heavy in ARTEFACTS:
        w = static_writes(cmd[1]) if len(cmd) > 1 else set()
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

    # Bayatligi SIRAYLA ILERI YAY: bir adim bayatsa ciktilari ve yazdigi
    # tablolar "kirli" sayilir; onlari okuyan sonraki adimlar da bayattir.
    # Eskiden bayatlik bir kez hesaplaniyordu ve bir ust adim yeniden
    # uretilse bile, basta guncel gorunen alt adimlar calismiyordu.
    rows = []
    dirty_files, dirty_tables = set(), set()
    for name, outputs, cmd, inputs, heavy in ARTEFACTS:
        why = staleness(name, outputs, inputs, state,
                        ignore=ignore[name], counts=counts_now)
        if not why:
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
        if args.only and name not in args.only:
            continue
        rows.append((name, outputs, cmd, inputs, heavy, why))

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
        # --cpu kabul eden betiklere gecir
        try:
            helptext = subprocess.run(full[:2] + ["--help"], capture_output=True,
                                      text=True, timeout=60).stdout
            if "--cpu" in helptext:
                full += ["--cpu", str(args.cpu)]
            elif "--threads" in helptext:
                full += ["--threads", str(args.cpu)]
        except (subprocess.SubprocessError, OSError):
            pass
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
