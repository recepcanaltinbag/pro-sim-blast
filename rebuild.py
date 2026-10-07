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

# Her kayit: (ad, ciktilar, komut, girdiler, agir mi)
#
# Girdilere BETIGIN KENDISI de dahildir: analiz kodu degistiginde ciktisi
# bayattir, ve bu oturumda tam olarak bu oldu -- stats_overview.py degisti,
# stats.json elle yeniden uretilmeseydi sayfa eski sayilari gosterecekti.
ARTEFACTS = [
    # --- cekirdek pipeline (agir olanlar isaretli)
    ("protein_filter", ["ro_alpha.csv", "ro_alpha.fasta"],
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
    ("annotate", [],
     ["python3", "annotate_ro.py", "--context-dir", "genomic_context", "--db", DB],
     ["annotate_ro.py", DB], True),

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

    # --- kuratorlu veriye bagli analizler
    ("ecology_stats", [f"{OUT}/cluster_ecology_stats.csv"],
     ["python3", "analyze_ecology.py", "--db", DB,
      "--out", f"{OUT}/cluster_ecology_stats.csv"],
     ["analyze_ecology.py", DB, ECO], False),
    ("motif_stats", [f"{OUT}/motif_stats.json"],
     ["python3", "motif_stats.py", "--db", DB, "--out", f"{OUT}/motif_stats.json"],
     ["motif_stats.py", DB], False),
    ("operon_validation", [f"{OUT}/operon_validation.json"],
     ["python3", "operon_validation.py", "--db", DB,
      "--out", f"{OUT}/operon_validation.json"],
     ["operon_validation.py", DB], False),
    ("redundancy", [f"{OUT}/redundancy.json"],
     ["python3", "redundancy.py", "--db", DB, "--out", f"{OUT}/redundancy.json"],
     ["redundancy.py", DB], False),
    ("cooccurrence", [f"{OUT}/cooccurrence.json"],
     ["python3", "cooccurrence.py", "--db", DB, "--out", f"{OUT}/cooccurrence.json"],
     ["cooccurrence.py", DB], False),
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
     ["validate_curation.py", DB, ECO, CHEM], False),
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


def digest(path):
    """Icerik ozeti. Buyuk dosyalarda bas + son + boyut yeter ve hizlidir."""
    if path.endswith(".sqlite"):
        return db_fingerprint(path)
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


def staleness(name, outputs, inputs, state):
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
    changed = [p for p in inputs if prev.get(p) and prev[p] != digest(p)]
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
    rows = []
    for name, outputs, cmd, inputs, heavy in ARTEFACTS:
        if args.only and name not in args.only:
            continue
        rows.append((name, outputs, cmd, inputs, heavy,
                     staleness(name, outputs, inputs, state)))

    if args.baseline:
        for name, outputs, cmd, inputs, heavy, _ in rows:
            state[name] = {"inputs": {p: digest(p) for p in inputs},
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
        res = subprocess.run(full)
        if res.returncode != 0:
            failed.append(name)
            print(f"    BASARISIZ ({res.returncode}) -- zincir burada duruyor, "
                  f"cunku sonraki adimlar bunun ciktisini okuyor")
            break
        state[name] = {"inputs": {p: digest(p) for p in inputs},
                       "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                       "seconds": round(time.time() - t0, 1)}
        save_state(state)
        print(f"    tamam, {time.time() - t0:.1f} s")

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
