"""
ROAR-DB web application -- Rieske Oxygenase Annotation Resource.

Serves roar.sqlite (built by the pipeline in the parent directory) as a
browsable, searchable, downloadable database with an on-the-fly sequence
classifier (hmmsearch + hmmalign motif test, same code as the pipeline).

Run locally:
    ROAR_DB=../roar.sqlite uvicorn app:app --reload --port 8000

Environment:
    ROAR_DB        path to roar.sqlite               (default: ../roar.sqlite)
    ROAR_HMM       cluster HMM file                   (default: ../RieskeDB71.hmm)
    ROAR_MOTIF_HMM motif HMM file                     (default: ../ROs_71_Clean/ROmotif71.hmm)
    ROAR_ECOLOGY   cluster_ecology.csv                (default: ../cluster_ecology.csv)
    ROAR_BASE      URL prefix when hosted under a sub-path (static export)
"""

import csv
import datetime as _dt
import io
import functools
import json
import hashlib
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, OrderedDict
from typing import Optional

from fastapi import FastAPI, Form, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, Response, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
for p in (HERE, PARENT):
    if p not in sys.path:
        sys.path.insert(0, p)

import ro_taxonomy                                               # noqa: E402
from ro_filter import classify_target, parse_domtblout          # noqa: E402
from ro_motif import classify_motifs, read_stockholm_matchcols  # noqa: E402

DB_PATH = os.environ.get("ROAR_DB", os.path.join(PARENT, "roar.sqlite"))
HMM_PATH = os.environ.get("ROAR_HMM", os.path.join(PARENT, "RieskeDB71.hmm"))
MOTIF_HMM_PATH = os.environ.get("ROAR_MOTIF_HMM",
                                os.path.join(PARENT, "ROs_71_Clean", "ROmotif71.hmm"))
ECOLOGY_PATH = os.environ.get("ROAR_ECOLOGY", os.path.join(PARENT, "cluster_ecology.csv"))
CHEMISTRY_PATH = os.environ.get("ROAR_CHEMISTRY", os.path.join(PARENT, "chemistry.csv"))
BASE = os.environ.get("ROAR_BASE", "").rstrip("/")
STATIC_MODE = False          # set by freeze.py: links get ".html" suffixes

MAX_CLASSIFY_SEQS = 50
MAX_CLASSIFY_CHARS = 400_000
PAGE_SIZE = 100
NEIGHBOR_WINDOW = 10000

app = FastAPI(title="ROAR-DB", docs_url="/api/docs", redoc_url=None)
app.mount("/static", StaticFiles(directory=os.path.join(HERE, "static")), name="static")
templates = Jinja2Templates(directory=os.path.join(HERE, "templates"))


# ----------------------------------------------------------------- helpers
def connect():
    if not os.path.exists(DB_PATH):
        raise HTTPException(500, f"database not found: {DB_PATH}")
    con = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True, check_same_thread=False)
    con.row_factory = sqlite3.Row
    return con


def slug(candidate_id):
    """'NC_004999.1:16050-17400:1' -> 'NC_004999.1_16050-17400_1' (URL/file safe)."""
    return candidate_id.replace(":", "_")


_UNSLUG = re.compile(r"^(.*)_(\d+)-(\d+)_(-?[01])$")


def unslug(s):
    m = _UNSLUG.match(s)
    if not m:
        return s.replace("_", ":")
    return f"{m.group(1)}:{m.group(2)}-{m.group(3)}:{m.group(4)}"


def leaf_slug(leaf_id):
    """'3_303_NBDO#0' -> '3_303_NBDO-0' ('#' would be a URL fragment)."""
    return leaf_id.replace("#", "-")


def leaf_unslug(s):
    return s.replace("-", "#", 1) if "#" not in s else s


def u(path, **params):
    """URL builder aware of base prefix and static export mode."""
    if STATIC_MODE:
        if path == "/":
            path = "/index.html"
        elif not path.endswith(".html") and not path.startswith(("/download/", "/static/", "/api/")):
            path = path + ".html"
    if params:
        from urllib.parse import urlencode
        path += "?" + urlencode({k: v for k, v in params.items() if v not in (None, "")})
    return BASE + path


def fmt(n):
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return "–"


def pct(a, b):
    return f"{100.0 * a / b:.1f}%" if b else "–"


def build_info():
    """Sayfa altinda gosterilen surum bilgisi: commit, tarih, icerik sayilari.

    Bilimsel bir kaynakta "hangi surumu kullandim" sorusu cevaplanabilir olmali;
    statik site de sunucu da ayni bilgiyi basar.
    """
    info = {"commit": None, "commit_date": None,
            "built": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d")}
    try:
        out = subprocess.run(["git", "-C", PARENT, "log", "-1", "--format=%h|%cs"],
                             capture_output=True, text=True, timeout=10)
        if out.returncode == 0 and "|" in out.stdout:
            info["commit"], info["commit_date"] = out.stdout.strip().split("|", 1)
    except (subprocess.SubprocessError, OSError):
        pass
    #  VERI SURUMU, kod surumunden AYRI bir seydir ve ikisini tek sayiya
    #  indirmek yanlis olur: ayni commit, yeniden kosturulmus bir pipeline
    #  ile baska sayilar yayinlayabilir. `db_revision.py` her veri surumunu
    #  DEGISMEZ bir kayit olarak tutuyor; atif edilecek olan budur.
    info["data_revision"] = None
    info["data_revision_date"] = None
    path = os.path.join(PARENT, "analysis_out", "revisions.jsonl")
    try:
        with open(path, encoding="utf-8") as handle:
            lines = [ln for ln in handle if ln.strip()]
        if lines:
            last = json.loads(lines[-1])
            info["data_revision"] = last.get("revision")
            info["data_revision_date"] = (last.get("recorded_utc") or "")[:10]
            info["data_revision_entries"] = last.get("confirmed_entries")
    except (OSError, ValueError):
        pass
    return info


BUILD = build_info()
# Kisa tip adi tek yerden gelir: uc ad cakisiyor ve cakisanlara numara
# eklenmesi her sayfada ayni sekilde olmak zorunda. Tembel yuklenir, cunku
# atlas modulu bu dosyanin sonunda import ediliyor.
@functools.lru_cache(maxsize=1)
def short_names():
    import atlas as _atlas
    return _atlas.short_names(os.path.join(PARENT, "chemistry.csv"))


def gene(cluster):
    """Tip id'sinin gosterilecek kisa adi."""
    if not cluster:
        return ""
    return short_names().get(cluster, cluster.split("_", 2)[-1])


templates.env.globals.update(u=u, slug=slug, leaf_slug=leaf_slug, fmt=fmt, pct=pct,
                             build=BUILD, gene=gene)
templates.env.filters["fmt"] = fmt


def _cite(text):
    """Kaynak metnini baglantili hale getirir; atlas gec import edildigi
    icin tembel baglanir."""
    import atlas as _atlas
    return _atlas.linkify_citation(text)


templates.env.filters["cite"] = _cite


def load_ecology():
    eco = {}
    if os.path.exists(ECOLOGY_PATH):
        with open(ECOLOGY_PATH) as fh:
            for row in csv.DictReader(fh):
                eco[row["cluster"]] = row
    return eco


ECOLOGY = load_ecology()


def load_chemistry():
    rows = {}
    if os.path.exists(CHEMISTRY_PATH):
        with open(CHEMISTRY_PATH) as fh:
            for r in csv.DictReader(fh):
                rows[r["cluster"]] = r
    return rows


CHEMISTRY = load_chemistry()

# O2 OKSIJENININ KADERI. "Rieske oksijenaz" bir aile adidir ve ailenin icinde
# urune hic oksijen sokmayan uyeler var (desaturasyon, C-C bagi kurma): O2
# orada yalnizca yukseltgen. Site "dioksijenaz" diye suzen birine onlari da
# veriyordu. Tablo oxygen_fate.py tarafindan types/ ve reaction_classes.yaml'dan
# TURETILIR; chemistry.csv'ye sutun EKLENMEDI cunku o dosyayi yirmiden fazla
# betik okuyor. Dosya yoksa sutun bos kalir, sayfa calisir.
OXYGEN_FATE_PATH = os.environ.get(
    "ROAR_OXYGEN_FATE", os.path.join(PARENT, "oxygen_fate.csv"))

# Yayimlanan etiketler. Alan adlari makine icin, bu metinler insan icin.
# O2'nin rolu -- oksijen kaderinden turetilir. Ayri sutun olarak duruyor
# cunku "hangileri yalnizca yukseltgen" sorusu okuyucunun iki degeri kafasinda
# birlestirmesine birakilmamali.
O2_ROLE_LABEL = {
    "donor": "oxidant + O donor",
    "oxidant_only": "oxidant only",
    "unknown": "unknown",
}

OXYGEN_FATE_LABEL = {
    "product": "into the product",
    "leaving_group": "into the leaving group",
    "oxidant_only": "not incorporated",
    "unknown": "unknown",
}


def load_oxygen_fate():
    rows = {}
    if os.path.exists(OXYGEN_FATE_PATH):
        with open(OXYGEN_FATE_PATH) as fh:
            for r in csv.DictReader(fh):
                rows[r["cluster"]] = r
    return rows


OXYGEN_FATE = load_oxygen_fate()

# METABOLIK ROL. Oksijen kaderinden AYRI bir eksen: SxtT ve PrnD oksijeni
# urune sokar ama biyodegradasyon yapmaz, vanillat/kafein demetilazlarinda ise
# oksijen kopan formaldehite gider ama bunlar klasik biyodegradasyondur. Rol,
# enzimin Rieske olup olmadigi hakkinda HICBIR SEY soylemez.
METABOLIC_ROLE_PATH = os.environ.get(
    "ROAR_METABOLIC_ROLE", os.path.join(PARENT, "metabolic_role.csv"))

METABOLIC_ROLE_LABEL = {
    "catabolic": "catabolic (degradation)",
    "biosynthetic": "biosynthetic",
    "unknown": "unknown",
}


def load_metabolic_role():
    rows = {}
    if os.path.exists(METABOLIC_ROLE_PATH):
        with open(METABOLIC_ROLE_PATH) as fh:
            for r in csv.DictReader(fh):
                rows[r["cluster"]] = r
    return rows


METABOLIC_ROLE = load_metabolic_role()


def _oxygen_fate_view(cluster):
    """Tip sayfasi icin O2 oksijeninin kaderi, etiketi ve GEREKCESI.

    Gerekce birlikte gonderilir cunku bu alan bir KARARDIR, veri degil: bazi
    tiplerde deger reaksiyon sinifindan gelir, bazilarinda tipin kendi beyani
    sinifi ezer. Sayfa degeri gosterip nereden geldigini soyleyemezse okur
    onu veri sanar.
    """
    r = OXYGEN_FATE.get(cluster)
    if not r:
        return None
    fate = r.get("oxygen_fate") or "unknown"
    return {"fate": fate,
            "label": OXYGEN_FATE_LABEL.get(fate, fate),
            "basis": r.get("basis") or "",
            "from_type": r.get("resolved_from") == "type",
            "transfers": fate in ("product", "leaving_group")}


@functools.lru_cache(maxsize=1)
def _substrate_chebi(db_mtime):
    """cluster -> [satir]. `type_substrate_chebi`'den okur, tahmin etmez.

    KIMLIK ile ADAY ayri tutulur (`state`): bir ad eslesmesi, ayni adi
    tasiyan baska bir molekule denk gelebilir, ve sayfa bunu okuyucuya
    SOYLEMEK zorunda (ROADMAP 135.3).
    """
    con = connect()
    try:
        out = {}
        for row in con.execute("""
                SELECT cluster, chebi, chebi_label, formula, state, rhea,
                       is_curated_substrate
                  FROM type_substrate_chebi
                 ORDER BY state DESC, is_curated_substrate DESC, chebi"""):
            out.setdefault(row["cluster"], []).append(dict(row))
        return out
    except Exception:
        #  Tablo henuz uretilmemisse sayfa calismaya devam eder; eksik
        #  bir kimlik, kirik bir sayfadan iyidir.
        return {}
    finally:
        con.close()


def substrate_chebi(cluster):
    return _substrate_chebi(os.path.getmtime(DB_PATH)).get(cluster) or []


def cluster_parts(cluster):
    parts = (cluster or "").split("_", 2)
    return {"group": parts[0] if parts else "", "id": parts[1] if len(parts) > 1 else "",
            "gene": parts[2] if len(parts) > 2 else cluster}


#  CINS: dizgiden kesilmiyor, NCBI rank'indan okunuyor. Burada
#  "Candidatus X -> X" kurali vardi; `split()[0]`'in sahte-cins hatasini
#  dogru teshis etmis ama dizgi icinde cozmeye calismisti ve kalanlari
#  yakalamiyordu: `Mycolicibacterium` bugun bir ALTCINS (612 giris),
#  `'Nostoc`/`[Enterobacter]` tirnakli-parantezli adlar, ve yeniden
#  adlandirmalar (Sinorhizobium -> Ensifer). Tek sahip: ro_taxonomy
#  (ROADMAP 135.1). Burada ONBELLEKLENIR, cunku her istekte 17 bin satir
#  okumak sayfa cizimine giriyor.
@functools.lru_cache(maxsize=1)
def _genus_cache(db_mtime):
    with sqlite3.connect(DB_PATH) as con:
        return ro_taxonomy.genus_map(con)


def genus_map():
    """nucleotide_id -> Genus. DB dosyasi degisince onbellek yenilenir."""
    try:
        stamp = os.path.getmtime(DB_PATH)
    except OSError:
        stamp = 0
    return _genus_cache(stamp)


def render(request, name, **ctx):
    ctx.setdefault("request", request)
    ctx.setdefault("static_mode", STATIC_MODE)
    #  `gtax` her sablonda hazir: `_genus_basis.html` cins sayan HER
    #  sayfadan include ediliyor ve her birinin route'una elle eklemek,
    #  bir sayfanin unutulmasi demekti -- bu projenin tekrar eden hatasi
    #  tam olarak bir duzeltmeyi yayimlandigi her yere tasimamak.
    ctx.setdefault("gtax", atlas.read_json(apath("taxonomy_ids.json")))
    try:                                   # Starlette >= 0.29: (request, name, context)
        return templates.TemplateResponse(request, name, ctx)
    except TypeError:                      # older Starlette: (name, context-with-request)
        return templates.TemplateResponse(name, ctx)


# ------------------------------------------------------------------ stats
def totals(con):
    s = lambda q: con.execute(q).fetchone()[0]  # noqa: E731
    t = {
        "confirmed": s("SELECT COUNT(*) FROM ro WHERE is_confirmed=1"),
        "candidates": s("SELECT COUNT(*) FROM ro"),
        "clusters": s("SELECT COUNT(DISTINCT ro_cluster) FROM ro WHERE is_confirmed=1"),
        "replicons": s("SELECT COUNT(*) FROM replicon WHERE status='ok'"),
        "genomes_with_ro": s("SELECT COUNT(DISTINCT nucleotide_id) FROM ro WHERE is_confirmed=1"),
        "neighbors": s("SELECT COUNT(*) FROM neighbor n JOIN ro r USING(candidate_id) WHERE r.is_confirmed=1"),
        "plasmid": s("SELECT COUNT(*) FROM ro r JOIN replicon p USING(nucleotide_id) "
                     "WHERE r.is_confirmed=1 AND p.is_plasmid=1"),
    }
    for table, key, q in [
            # KOHORTA GORE YAYGINLIK (2026-11-02). "N replikon tarandi,
            # M tanesi en az bir tane tasiyor" bir YAYGINLIK iddiasi ve
            # paydasi "taranan genomlar". A5'te hedeflenerek cekilen kayitlar
            # tam da ICINDE RO OLDUGU BILINDIGI ICIN secildi; ayni paydaya
            # girerse oran ~%100'e cikar ve bir kesif gibi gorunur -- oysa
            # degisen sey olcum degil ORNEKLEME. Bu yuzden payda ve pay
            # YALNIZCA survey kohortundan okunuyor.
            ("replicon_cohort", "survey_replicons",
             "SELECT COUNT(*) FROM replicon p JOIN replicon_cohort c "
             "USING(nucleotide_id) WHERE p.status='ok' AND c.cohort='survey'"),
            ("replicon_cohort", "survey_with_ro",
             "SELECT COUNT(DISTINCT r.nucleotide_id) FROM ro r "
             "JOIN replicon_cohort c USING(nucleotide_id) "
             "WHERE r.is_confirmed=1 AND c.cohort='survey'"),
            # Hedeflenerek cekilenler ayri sayilir: gizlenmiyor, yalnizca
            # yaygınlık oranina sokulmuyor.
            ("replicon_cohort", "recovered_replicons",
             "SELECT COUNT(*) FROM replicon p JOIN replicon_cohort c "
             "USING(nucleotide_id) WHERE p.status='ok' "
             "AND c.cohort='protein_recovery'"),
            ("leaf", "leaves", "SELECT COUNT(*) FROM leaf"),
            ("operon", "operon_complete", "SELECT COUNT(*) FROM operon WHERE completeness=3"),
            ("operon", "operon_any", "SELECT COUNT(*) FROM operon WHERE completeness>0"),
            ("ro_subfamily", "novel", "SELECT COUNT(*) FROM ro_subfamily WHERE assignment_class='novel_candidate'"),
            ("ro_domain", "eukaryote", "SELECT COUNT(*) FROM ro_domain WHERE domain='Eukaryota'"),
            # 2026-10-30, 3/3 triad kapisi (ro_triad). Bu sayilar about.html'de
            # gosteriliyor ve ELLE YAZILMAMALI: triad okumasi hala degisiyor
            # (His capalamasi yeni), elle yazilan bir sayi bu projede zaten
            # bir kere sessizce bayatladi.
            ("ro_triad", "triad_complete",
             "SELECT COUNT(*) FROM ro_triad WHERE triad_complete=1"),
            # Kanonik olmayan His-His araligi: ligandlar homolog kolonlarda
            # ama ara literaturun 3-5'i degil. "Tam" sayiliyorlar; sayinin
            # icinde KAYBOLMAMALARI icin ayrica sayiliyor.
            ("ro_triad", "triad_non_canonical_spacing",
             "SELECT COUNT(*) FROM ro_triad WHERE "
             "triad_status='complete_non_canonical_his_spacing'"),
            ("ro_triad", "triad_incomplete",
             "SELECT COUNT(*) FROM ro_triad WHERE triad_complete=0"),
            ("ro_triad", "triad_absent",
             "SELECT COUNT(*) FROM ro_triad WHERE triad_status='carboxylate_absent'"),
            # His'i okunamayan girisler. BU SAYI ELLE YAZILIYDU (about.html'de
            # "20"), ve 2026-10-31'de yanlis cikti: capalama yeniden kostuktan
            # sonra 31 oldu, cumle de 19.114'e toplaniyordu ama 19.125'ten
            # soz ediyordu. Bu yuzden buraya alindi -- triad cumlesindeki bes
            # sayinin dordu zaten canliydi, elle yazilan tek sayi bu idi.
            ("ro_triad", "triad_his_absent",
             "SELECT COUNT(*) FROM ro_triad WHERE triad_status='his_absent'"),
            ("ro_triad", "triad_undecided",
             "SELECT COUNT(*) FROM ro_triad WHERE triad_status IN "
             "('his_needs_anchoring','carboxylate_untested')"),
            # 2026-10-31, uyelik TEMELI (ro_membership). 3/3 kapisi artik
            # uyelikte karar veriyor -- ama silme olarak degil, temel olarak:
            # triad'in istisnasiz oldugunu yazan bir kaynak bulunamadi, o
            # yuzden eksik triad'li giris atilmaz, hangi temelle uye oldugu
            # yayinlanir. Paydalar bu tablodan okunur.
            ("ro_membership", "basis_triad_complete",
             "SELECT COUNT(*) FROM ro_membership WHERE basis='triad_complete'"),
            ("ro_membership", "basis_documented_exception",
             "SELECT COUNT(*) FROM ro_membership WHERE basis='documented_exception'"),
            ("ro_membership", "basis_triad_incomplete",
             "SELECT COUNT(*) FROM ro_membership "
             "WHERE basis='triad_incomplete_retained'"),
            # 2026-11-01, BIRLESIM (ro_union). `confirmed` bir KAPSAM degil
            # KESISIM'dir: hem genomik baglami olan, hem kapiyi gecen
            # girisler. Bu sayilar o kesisimin disinda kalani gorunur kilar
            # ve ELLE YAZILMAMALI -- birlesim her kaynak eklendiginde buyur.
            ("ro_union", "union_sequences",
             "SELECT COUNT(*) FROM ro_union"),
            ("ro_union", "union_no_context",
             "SELECT COUNT(*) FROM ro_union WHERE has_genomic_context=0"),
            ("ro_union", "union_ungated",
             "SELECT COUNT(*) FROM ro_union WHERE passes_gate=0"),
            ("ro_union_source", "union_source_ids",
             "SELECT COUNT(*) FROM ro_union_source"),
            ("ro_union_source", "union_sources_dead",
             "SELECT COUNT(*) FROM ro_union_source WHERE status='deleted'"),
            # TEK KOPYA: genomik baglami olmayan VE tum kaynak kimlikleri
            # yukari akista silinmis diziler. "47.988 accession silindi,
            # yani 47.988 dizinin tek kopyasi bizde" demek YANLIS olurdu --
            # silinenlerin bir kismi genomik tarafta ya da canli bir baska
            # accession'da duruyor. Kosul bu yuzden diziye bakiyor,
            # accession'a degil.
            ("ro_union", "union_only_copy",
             "SELECT COUNT(*) FROM ro_union u WHERE u.has_genomic_context=0 "
             "AND NOT EXISTS (SELECT 1 FROM ro_union_source s "
             "WHERE s.seq_md5=u.seq_md5 "
             "AND s.status NOT IN ('deleted','demerged'))")]:
        if con.execute("SELECT name FROM sqlite_master WHERE name=?", (table,)).fetchone():
            t[key] = s(q)
    return t


def cluster_table(con):
    rows = con.execute("""
        SELECT r.ro_cluster cluster, COUNT(*) n,
               COUNT(DISTINCT r.nucleotide_id) replicons,
               SUM(p.is_plasmid) plasmid,
               SUM(o.has_beta) beta, SUM(o.has_ferredoxin) ferredoxin,
               SUM(o.has_reductase) reductase, SUM(o.completeness=3) complete
        FROM ro r
        JOIN replicon p ON p.nucleotide_id = r.nucleotide_id
        LEFT JOIN operon o ON o.candidate_id = r.candidate_id
        WHERE r.is_confirmed=1
        GROUP BY r.ro_cluster ORDER BY n DESC""").fetchall()
    leaves = dict(con.execute("SELECT cluster, COUNT(*) FROM leaf GROUP BY cluster"))
    euk = dict(con.execute(
        "SELECT cluster, SUM(domain='Eukaryota') FROM ro_domain GROUP BY cluster"))
    gmap = genus_map()
    genera = {}
    for cl, nuc in con.execute(
            "SELECT ro_cluster, nucleotide_id FROM ro WHERE is_confirmed=1"):
        got = gmap.get(nuc)
        if got is not None and got.is_ranked:
            genera.setdefault(cl, set()).add(got.key)
    out = []
    for r in rows:
        out.append(_cluster_row(dict(r), leaves, euk, genera))
    # Uyesi olmayan tipler ro tablosunda hic satir ACMIYOR, dolayisiyla bu
    # listede gorunmuyorlardi: 71 kuratorlu tipten 61'i listelenip geri kalan
    # 10'u sessizce kayboluyordu. Simdi sifir uyeyle listeleniyor ve sebebi
    # empty_types.json'dan okunup satira isaretleniyor.
    listed = {d["cluster"] for d in out}
    for cluster in sorted(set(CHEMISTRY) | set(ECOLOGY)):
        if cluster in listed:
            continue
        out.append(_cluster_row(
            {"cluster": cluster, "n": 0, "replicons": 0, "plasmid": 0, "beta": 0,
             "ferredoxin": 0, "reductase": 0, "complete": 0},
            leaves, euk, genera))
    return out


def _cluster_row(d, leaves, euk, genera):
    """cluster_table satirini kimya/ekoloji kuratorlu alanlariyla tamamlar."""
    cluster = d["cluster"]
    eco = ECOLOGY.get(cluster, {})
    chem = CHEMISTRY.get(cluster, {})
    d.update(cluster_parts(cluster))
    ofate = OXYGEN_FATE.get(cluster, {})
    mrole = METABOLIC_ROLE.get(cluster, {})
    d.update({"family": chem.get("family", ""),
              "reaction_class": chem.get("reaction_class", ""),
              "oxygen_fate": ofate.get("oxygen_fate", ""),
              "oxygen_fate_label": OXYGEN_FATE_LABEL.get(
                  ofate.get("oxygen_fate", ""), ""),
              "oxygen_fate_basis": ofate.get("basis", ""),
              "o2_role": ofate.get("o2_role", ""),
              "o2_role_label": O2_ROLE_LABEL.get(ofate.get("o2_role", ""), ""),
              "metabolic_role": mrole.get("metabolic_role", ""),
              "metabolic_role_label": METABOLIC_ROLE_LABEL.get(
                  mrole.get("metabolic_role", ""), ""),
              "metabolic_role_basis": mrole.get("basis", ""),
              "reaction": chem.get("reaction", ""),
              "pdb": chem.get("pdb", ""),
              "substrate": chem.get("substrate_en") or eco.get("substrate", ""),
              "sclass": eco.get("substrate_class", ""),
              "confidence": eco.get("confidence", ""), "leaves": leaves.get(cluster, 0),
              "euk": euk.get(cluster, 0) or 0, "genera": len(genera.get(cluster, ())),
              "empty_type": atlas.empty_type_view(cluster, apath("empty_types.json"))})
    return d


# ------------------------------------------------------------------ pages
@app.get("/", response_class=HTMLResponse)
def home(request: Request):
    con = connect()
    try:
        counts = dict(con.execute("SELECT ro_cluster, COUNT(*) FROM ro WHERE is_confirmed=1 "
                                  "GROUP BY ro_cluster"))
        groups = atlas.chemistry_groups(CHEMISTRY_PATH, counts, ECOLOGY)
        return render(request, "index.html", totals=totals(con), groups=groups,
                      reactions=atlas.reaction_summary(groups),
                      significance=atlas.significance_view(
                          CHEMISTRY_PATH, counts,
                          os.path.join(PARENT, "significance.yaml"), con),
                      motif=atlas.read_json(apath("motif_stats.json")))
    finally:
        con.close()


@app.get("/architecture", response_class=HTMLResponse)
def architecture(request: Request):
    """Anasayfanin 1. bolumu 2026-11-02'de BURAYA tasindi.

    Kullanici: "anasayfada cok yazi var ... ozellikle 1.2 ve ilk giris
    kismi, hic iyi bi database gibi durmuyo", sonra "kalabilirler ama daha
    ozet gibi olmalila". Olculdu: anasayfa 3.933 kelime, bu bolum 1.477'si.
    Icerik /atlas/structure, /atlas/evidence ve /about'ta GECMIYOR (olculdu),
    yani silinebilir degil. Anasayfa artik semayi, dort tamamlanma oranini
    ve bir paragraflik bulguyu tutuyor; gerekcenin tamami burada.
    """
    con = connect()
    try:
        return render(request, "architecture.html",
                      totals=totals(con),
                      motif=atlas.read_json(apath("motif_stats.json")))
    finally:
        con.close()


@app.get("/clusters", response_class=HTMLResponse)
def clusters(request: Request):
    con = connect()
    try:
        return render(request, "clusters.html", clusters=cluster_table(con))
    finally:
        con.close()


def cluster_detail(con, cluster):
    n = con.execute("SELECT COUNT(*) FROM ro WHERE is_confirmed=1 AND ro_cluster=?",
                    (cluster,)).fetchone()[0]
    if not n:
        # Uyesi olmayan tip bir HATA degildir: referans enzim kuratorlu sette
        # var ama taranan genomlarda ona atanan dogrulanmis protein cikmamis.
        # Ana sayfadaki kimya kartlari chemistry.csv'den geliyor (71 tip),
        # oysa ro tablosunda 61 tip var; aradaki fark 404 veriyordu.
        if cluster not in CHEMISTRY and cluster not in ECOLOGY:
            return None
        eco = ECOLOGY.get(cluster, {})
        return {"cluster": cluster, "n": 0, **cluster_parts(cluster),
                "substrate": CHEMISTRY.get(cluster, {}).get("substrate_en")
                             or eco.get("substrate", ""),
                "sclass": eco.get("substrate_class", ""),
                "note": eco.get("ecology_note", ""),
                "confidence": eco.get("confidence", ""),
                "replicons": 0, "plasmid": 0, "score": None, "coverage": None,
                "operon": {"beta": 0, "ferredoxin": 0, "reductase": 0, "complete": 0,
                           "nb_beta": 0, "nb_ferredoxin": 0, "nb_reductase": 0, "total": 0},
                "domains": [], "genera": [], "top_genera": [], "genera_n": 0,
                "neighborhood": [], "regulators": [], "layouts": [], "layouts_distant": [],
                "leaves": [],
                "classes": [], "empty": True}
    eco = ECOLOGY.get(cluster, {})
    info = {"cluster": cluster, "n": n, **cluster_parts(cluster),
            "substrate": eco.get("substrate", ""), "sclass": eco.get("substrate_class", ""),
            "note": eco.get("ecology_note", ""), "confidence": eco.get("confidence", "")}
    row = con.execute("""
        SELECT COUNT(DISTINCT r.nucleotide_id) replicons, SUM(p.is_plasmid) plasmid,
               AVG(r.hmm_score) score, AVG(r.model_coverage) coverage
        FROM ro r JOIN replicon p USING(nucleotide_id)
        WHERE r.is_confirmed=1 AND r.ro_cluster=?""", (cluster,)).fetchone()
    info.update(dict(row))
    info["operon"] = dict(con.execute("""
        SELECT SUM(has_beta) beta, SUM(has_ferredoxin) ferredoxin, SUM(has_reductase) reductase,
               SUM(completeness=3) complete, SUM(nearby_beta) nb_beta,
               SUM(nearby_ferredoxin) nb_ferredoxin, SUM(nearby_reductase) nb_reductase,
               COUNT(*) total
        FROM operon o JOIN ro r USING(candidate_id)
        WHERE r.is_confirmed=1 AND r.ro_cluster=?""", (cluster,)).fetchone())
    info["domains"] = con.execute(
        "SELECT domain, euk_group, COUNT(*) n FROM ro_domain WHERE cluster=? "
        "GROUP BY domain, euk_group ORDER BY n DESC", (cluster,)).fetchall()
    info["genera"] = con.execute("""
        SELECT organism, COUNT(*) n FROM ro r JOIN replicon p USING(nucleotide_id)
        WHERE r.is_confirmed=1 AND r.ro_cluster=? GROUP BY organism ORDER BY n DESC LIMIT 400""",
        (cluster,)).fetchall()
    #  Cins sayimi organizma dizgisinden degil replikonun taxid'inden
    #  gelir; cinsi NCBI'da atanmamis girisler bir cins SAYILMAZ.
    gmap = genus_map()
    gcount = Counter()
    for nuc, n in con.execute(
            "SELECT nucleotide_id, COUNT(*) FROM ro WHERE is_confirmed=1 "
            "AND ro_cluster=? GROUP BY 1", (cluster,)):
        got = gmap.get(nuc)
        if got is not None and got.is_ranked:
            gcount[got.key] += n
    info["top_genera"] = gcount.most_common(12)
    info["genera_n"] = len(gcount)
    # neighborhood: categories within +-4 genes, per-RO rate
    info["neighborhood"] = con.execute("""
        SELECT c.category, COUNT(DISTINCT nb.candidate_id) ros, COUNT(*) genes
        FROM neighbor nb
        JOIN ro r ON r.candidate_id = nb.candidate_id AND r.is_confirmed=1 AND r.ro_cluster=?
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        WHERE ABS(nb.gene_offset) <= 4
        GROUP BY c.category ORDER BY ros DESC""", (cluster,)).fetchall()
    info["regulators"] = con.execute("""
        SELECT c.category family, COUNT(*) n
        FROM neighbor nb
        JOIN ro r ON r.candidate_id = nb.candidate_id AND r.is_confirmed=1 AND r.ro_cluster=?
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regulator_family_v1'
        WHERE ABS(nb.gene_offset) <= 4
        GROUP BY c.category ORDER BY n DESC""", (cluster,)).fetchall()
    # Gen sirasi YAKINLIGA gore ayrilir. Kullanici hakli: uzak bir akrabanin
    # operon sirasi bu tipin kanonik sirasi hakkinda bir sey soylemiyor, ama
    # ayni listede durunca oyleymis gibi okunuyordu. Yakin uyeler (karakterize
    # ve yakin homolog) once, gerisi ayri bir blokta.
    info["layouts"] = con.execute("""
        SELECT o.layout, COUNT(*) n FROM operon o JOIN ro r USING(candidate_id)
        LEFT JOIN ro_evidence e USING(candidate_id)
        WHERE r.is_confirmed=1 AND r.ro_cluster=?
          AND e.tier IN ('characterized', 'close_homolog')
        GROUP BY o.layout ORDER BY n DESC LIMIT 8""", (cluster,)).fetchall()
    info["layouts_distant"] = con.execute("""
        SELECT o.layout, COUNT(*) n FROM operon o JOIN ro r USING(candidate_id)
        LEFT JOIN ro_evidence e USING(candidate_id)
        WHERE r.is_confirmed=1 AND r.ro_cluster=?
          AND (e.tier IS NULL OR e.tier NOT IN ('characterized', 'close_homolog'))
        GROUP BY o.layout ORDER BY n DESC LIMIT 8""", (cluster,)).fetchall()
    info["leaves"] = con.execute("""
        SELECT l.leaf_id, l.size, l.median_identity, l.is_homogeneous, l.depth,
               l.top_genera,
               l.representative, lp.label, lp.plasmid_rate, lp.neighbor_signature
        FROM leaf l LEFT JOIN leaf_profile lp USING(leaf_id)
        WHERE l.cluster=? ORDER BY l.size DESC""", (cluster,)).fetchall()
    info["classes"] = con.execute(
        "SELECT assignment_class, COUNT(*) n FROM ro_subfamily WHERE cluster=? "
        "GROUP BY assignment_class ORDER BY n DESC", (cluster,)).fetchall()
    return info


def members_query(con, cluster=None, leaf=None, q=None, page=1, size=PAGE_SIZE):
    where, params = ["r.is_confirmed=1"], []
    if cluster:
        where.append("r.ro_cluster=?"); params.append(cluster)
    if leaf:
        where.append("rl.leaf_id=?"); params.append(leaf)
    if q:
        where.append("(p.organism LIKE ? OR r.product LIKE ? OR r.protein_id LIKE ? "
                     "OR r.locus_tag LIKE ? OR r.candidate_id LIKE ?)")
        params += [f"%{q}%"] * 5
    sql = f"""
        FROM ro r JOIN replicon p USING(nucleotide_id)
        LEFT JOIN ro_leaf rl ON rl.candidate_id = r.candidate_id
        LEFT JOIN ro_subfamily rs ON rs.candidate_id = r.candidate_id
        LEFT JOIN ro_evidence ev ON ev.candidate_id = r.candidate_id
        LEFT JOIN operon o ON o.candidate_id = r.candidate_id
        WHERE {' AND '.join(where)}"""
    total = con.execute("SELECT COUNT(*) " + sql, params).fetchone()[0]
    rows = con.execute(f"""
        SELECT r.candidate_id, r.protein_id, r.locus_tag, r.product, r.ro_cluster, r.hmm_score,
               r.model_coverage, p.organism, p.is_plasmid, p.nucleotide_id,
               rl.leaf_id, rs.assignment_class, rs.core_identity,
               o.completeness, o.layout, ev.tier, ev.ref_identity, ev.nearest_ref
        {sql}
        -- EN YAKIN uyeler once. Kullanici hakli olarak "enzim sayfalarinda
        -- once cok yakin homologlara bakalim, uzak akrabalari ayri
        -- degerlendirelim" dedi; organizma adina gore siralamak en bilgilendirici
        -- uyeleri sayfalar arasina dagitiyordu.
        ORDER BY CASE ev.tier WHEN 'characterized' THEN 0 WHEN 'close_homolog' THEN 1
                              WHEN 'family_member' THEN 2 WHEN 'distant' THEN 3
                              WHEN 'novel' THEN 4 ELSE 5 END,
                 ev.ref_identity DESC, p.organism, r.candidate_id
        LIMIT ? OFFSET ?""",
        params + [size, (page - 1) * size]).fetchall()
    return total, rows


# Emekliye ayrilmis tip kimlikleri -> yerine gececek tip. Disaridan verilmis
# bir baglanti kirilmasin diye 301 ile yonlendirilir; statik yayinda
# freeze.py ayni tablodan meta-refresh sayfalari yazar.
#
# OMO, NDO(3_315) ve NDO(3_314) kaldirildi cunku dizileri sirasiyla OxoO,
# NahAc ve NarAa ile BIREBIR ayniydi: iki referans ayni diziye sahipse
# yarisma beraberlikle biter ve biri keyfi olarak butun uyeleri alir.
# NdmC'nin yonlendirilecegi bir yer YOK -- o bir kopya degil, Rieske motifi
# tasimayan bir protein; types/retired/ icinde gerekcesiyle duruyor.
RENAMED_TYPES = {
    "5_505_OxyA": "5_505_qxyA",
    "3_304_OMO": "1_101_OxoO",
    "3_315_NDO": "3_309_NahAc",
    "3_314_NDO": "3_316_NarAa",
}


@app.get("/cluster/{cluster}", response_class=HTMLResponse)
def cluster_page(request: Request, cluster: str, page: int = 1, q: Optional[str] = None):
    con = connect()
    try:
        # Yeniden adlandirilmis bir tipin ESKI adresi: 404 yerine yonlendir,
        # cunku disarida verilmis baglantilar ve yer imleri var.
        if cluster in RENAMED_TYPES:
            return RedirectResponse(u("/cluster/" + RENAMED_TYPES[cluster]), status_code=301)
        # Yonlendirilecek yeri OLMAYAN emekli tip: bos sayfa yerine gerekce.
        gone = atlas.retired_type_view(cluster)
        if gone:
            return render(request, "retired_type.html", gone=gone)
        info = cluster_detail(con, cluster)
        if not info:
            raise HTTPException(404, "cluster not found")
        total, rows = members_query(con, cluster=cluster, q=q, page=page)
        return render(request, "cluster.html", c=info, members=rows, total=total, page=page,
                      pages=(total + PAGE_SIZE - 1) // PAGE_SIZE, q=q or "",
                      chem=CHEMISTRY.get(cluster), x=cluster_extras(con, cluster),
                      ofate=_oxygen_fate_view(cluster),
                      chebi=substrate_chebi(cluster),
                      mrole=METABOLIC_ROLE.get(cluster))
    finally:
        con.close()


@app.get("/leaf/{leaf_id}", response_class=HTMLResponse)
def leaf_page(request: Request, leaf_id: str, page: int = 1):
    leaf_id = leaf_unslug(leaf_id)
    # Yeniden adlandirilan tipin varyant adresleri de yonlendirilir: ad
    # yaprak kimliginin ICINDE geciyor ve disarida verilmis baglantilar var.
    for old_name, new_name in RENAMED_TYPES.items():
        if leaf_id.startswith(old_name + "#"):
            return RedirectResponse(
                u("/leaf/" + leaf_slug(leaf_id.replace(old_name, new_name, 1))),
                status_code=301)
    con = connect()
    try:
        leaf = con.execute("""
            SELECT l.*, lp.label, lp.plasmid_rate, lp.transposon_rate, lp.neighbor_signature,
                   lp.enriched_categories, lp.genus_n
            FROM leaf l LEFT JOIN leaf_profile lp USING(leaf_id) WHERE l.leaf_id=?""",
            (leaf_id,)).fetchone()
        if not leaf:
            raise HTTPException(404, "leaf not found")
        total, rows = members_query(con, leaf=leaf_id, page=page)
        sdp = None
        if atlas.table_exists(con, "leaf_sdp"):
            r = con.execute("SELECT residues, signature, n_differ FROM leaf_sdp WHERE leaf_id=?",
                            (leaf_id,)).fetchone()
            c2 = con.execute("SELECT columns FROM cluster_sdp WHERE cluster=?",
                             (leaf["cluster"],)).fetchone()
            if r and c2:
                import json as _json
                sdp = {"residues": _json.loads(r[0]), "signature": r[1], "n_differ": r[2],
                       "columns": _json.loads(c2[0])}
        return render(request, "leaf.html", leaf=leaf, c=cluster_parts(leaf["cluster"]),
                      members=rows, total=total, page=page, sdp=sdp,
                      pages=(total + PAGE_SIZE - 1) // PAGE_SIZE)
    finally:
        con.close()


# ------------------------------------------------------------------ RO entry
COMPONENT_COLORS = {
    "alpha": "#c0392b", "alpha_other": "#e67e22", "beta": "#2980b9", "ferredoxin": "#27ae60",
    "reductase": "#f39c12", "rieske_other": "#16a085",
    "regulator": "#8e44ad", "transposon": "#2c3e50", "ring_cleavage": "#d35400",
    "transporter": "#7f8c8d", "dehydrogenase": "#95a5a6", "hydrolase": "#95a5a6",
    "hypothetical": "#dfe6e9", "other": "#b2bec3", "none": "#b2bec3",
}


def neighborhood_svg(ro, neighbors, operon_ids):
    """Gene arrow diagram of the +-10 kb window. Returns SVG markup."""
    left = max(0, ro["start"] - NEIGHBOR_WINDOW)
    right = ro["end"] + NEIGHBOR_WINDOW
    genes = [dict(n) for n in neighbors if not n["spans_origin"]]
    if genes:
        left = min(left, min(g["start"] for g in genes))
        right = max(right, max(g["end"] for g in genes))
    span = max(1, right - left)
    W, H, Y = 1000, 110, 55
    sx = lambda x: (x - left) / span * W  # noqa: E731

    def arrow(start, end, strand, color, label, title, bold=False, box=False):
        x1, x2 = sx(start), sx(end)
        if x2 - x1 < 4:
            x2 = x1 + 4
        h = 14 if not bold else 18
        head = min(8, (x2 - x1) / 2)
        if strand == -1:
            pts = [(x2, Y - h), (x1 + head, Y - h), (x1, Y), (x1 + head, Y + h), (x2, Y + h)]
        else:
            pts = [(x1, Y - h), (x2 - head, Y - h), (x2, Y), (x2 - head, Y + h), (x1, Y + h)]
        poly = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        stroke = "#111" if bold else "#555"
        out = f'<polygon points="{poly}" fill="{color}" stroke="{stroke}" stroke-width="{1.5 if bold else 0.8}"><title>{title}</title></polygon>'
        if box:
            out += f'<rect x="{x1-1:.1f}" y="{Y-h-5}" width="{x2-x1+2:.1f}" height="{2*h+10}" fill="none" stroke="#2d3436" stroke-dasharray="3,2" stroke-width="0.8"/>'
        if label and x2 - x1 > 24:
            ty = Y + h + 14 if strand != -1 else Y - h - 6
            out += f'<text x="{(x1+x2)/2:.1f}" y="{ty}" font-size="9" text-anchor="middle" fill="#2d3436">{label[:14]}</text>'
        return out

    import html as _h
    parts = [f'<svg viewBox="0 0 {W} {H}" width="100%" preserveAspectRatio="none" '
             f'xmlns="http://www.w3.org/2000/svg" class="genes">',
             f'<line x1="0" y1="{Y}" x2="{W}" y2="{Y}" stroke="#999" stroke-width="1"/>']
    for g in genes:
        comp = g.get("component") or "none"
        key = comp if comp != "none" else (g.get("category") or "other")
        color = COMPONENT_COLORS.get(key, "#b2bec3")
        label = g["gene"] or (comp if comp != "none" else "")
        title = _h.escape(f"{g['product']} | {g['start']}-{g['end']} ({'+' if g['strand']==1 else '-'}) | "
                          f"component: {comp} | category: {g.get('category')} | distance {g['distance']} bp")
        parts.append(arrow(g["start"], g["end"], g["strand"], color, label, title,
                           box=g["neighbor_id"] in operon_ids))
    parts.append(arrow(ro["start"], ro["end"], ro["strand"], COMPONENT_COLORS["alpha"],
                       ro["gene"] or "RO α", _h.escape(f"{ro['product']} (this entry)"),
                       bold=True, box=True))
    # scale bar
    kb = 2000
    parts.append(f'<line x1="{W-10-kb/span*W:.1f}" y1="{H-6}" x2="{W-10}" y2="{H-6}" stroke="#333" stroke-width="2"/>'
                 f'<text x="{W-10-kb/span*W/2:.1f}" y="{H-9}" font-size="9" text-anchor="middle">2 kb</text>')
    parts.append("</svg>")
    return "".join(parts)


def ro_detail(con, candidate_id):
    ro = con.execute("""
        SELECT r.*, p.organism, p.taxonomy, p.description, p.is_plasmid, p.is_circular,
               p.length replicon_length, p.cds_count,
               rl.leaf_id, rl.leaf_size, rs.assignment_class, rs.core_identity, rs.subfamily_id,
               d.domain, d.euk_group,
               t.taxid_current taxid, g.genus_name, g.genus_state,
               g.proxy_rank, g.proxy_name
        FROM ro r JOIN replicon p USING(nucleotide_id)
        LEFT JOIN replicon_taxon t ON t.nucleotide_id = r.nucleotide_id
        LEFT JOIN taxon_genus g ON g.nucleotide_id = r.nucleotide_id
        LEFT JOIN ro_leaf rl ON rl.candidate_id = r.candidate_id
        LEFT JOIN ro_subfamily rs ON rs.candidate_id = r.candidate_id
        LEFT JOIN ro_domain d ON d.candidate_id = r.candidate_id
        WHERE r.candidate_id=?""", (candidate_id,)).fetchone()
    if not ro:
        return None
    neighbors = con.execute("""
        SELECT nb.*, c.category, nc.component, nc.evidence
        FROM neighbor nb
        LEFT JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        LEFT JOIN neighbor_component nc
               ON nc.protein_key = nb.nucleotide_id||':'||nb.start||'-'||nb.end||':'||nb.strand
        WHERE nb.candidate_id=? ORDER BY nb.gene_offset""", (candidate_id,)).fetchall()
    operon = con.execute("SELECT * FROM operon WHERE candidate_id=?", (candidate_id,)).fetchone()
    operon_genes = con.execute(
        "SELECT * FROM operon_gene WHERE candidate_id=? ORDER BY position", (candidate_id,)).fetchall()
    operon_ids = {g["neighbor_id"] for g in operon_genes if g["neighbor_id"] is not None}
    # Triad okumasi ve UYELIK TEMELI. Sayfa eskiden kosulsuz "carries an intact
    # Rieske centre and catalytic triad" diyordu; ro_triad'a gore bu, ligandi
    # gercekten eksik olan girislerde YANLIS bir iddiaydi. Tablolar yoksa None
    # doner ve sayfa cumleyi hic kurmaz -- varsayilan iddia etmek degil.
    triad = membership = None
    if con.execute("SELECT name FROM sqlite_master "
                   "WHERE name='ro_triad'").fetchone():
        triad = con.execute("SELECT * FROM ro_triad WHERE candidate_id=?",
                            (candidate_id,)).fetchone()
    if con.execute("SELECT name FROM sqlite_master "
                   "WHERE name='ro_membership'").fetchone():
        membership = con.execute("SELECT * FROM ro_membership WHERE candidate_id=?",
                                 (candidate_id,)).fetchone()
    # motif residues on the stored sequence are not persisted; report flags only
    return {"ro": ro, "neighbors": neighbors, "operon": operon, "operon_genes": operon_genes,
            "svg": neighborhood_svg(ro, neighbors, operon_ids), "operon_ids": operon_ids,
            "cluster": cluster_parts(ro["ro_cluster"]),
            "triad": triad, "membership": membership,
            "eco": ECOLOGY.get(ro["ro_cluster"], {})}


@app.get("/entry", response_class=HTMLResponse)
def entry_page(request: Request):
    """Tek giris sayfasi kabugu; icerigi tarayici data/ parcalarindan cizer.

    Hem gelistirme sunucusunda hem statik ihracatta AYNI yol kullanilir ki iki
    render yolu ayrisamasin. Eski /ro/<slug> rotasi, sunucu tarafinda
    karsilastirma yapan verify_entry_render.py icin duruyor.
    """
    return render(request, "entry.html")


def entry_href(candidate_id):
    """Bir girisin adresi. Tek yerde tanimli, sablonlar bunu cagirir."""
    return u("/entry", id=slug(candidate_id))


@app.get("/ro/{ro_slug}", response_class=HTMLResponse)
def ro_page(request: Request, ro_slug: str):
    con = connect()
    try:
        d = ro_detail(con, unslug(ro_slug))
        if not d:
            raise HTTPException(404, "entry not found")
        return render(request, "ro.html",
                      x=ro_extras(con, d["ro"]["candidate_id"]),
                      xrefs=cross_refs(con, d["ro"]["sequence"]),
                      ecrev=ec_review_note(con, d["ro"]["sequence"]), **d)
    finally:
        con.close()


# ------------------------------------------------------------------ search
FACETS = [("tier", "Evidence level"), ("domain", "Domain of life"),
          ("family", "Chemical family"), ("cluster", "Enzyme type")]
FLAGS = [("plasmid", "is_plasmid", "on a plasmid"),
         ("partner", "has_partner", "operon partner verified"),
         ("regulator", "has_regulator", "divergent regulator upstream")]


# Tam eslesme bu sayidan az sonuc verdiyse onek eslesemesi denenir. Dusuk
# tutuluyor: amac bos ya da neredeyse bos bir sonucu kurtarmak, yuzlerce
# dogru sonucu sulandirmak degil.
WIDEN_BELOW = 25


def fts_expression(text):
    """Turn user input into a safe FTS5 expression.

    Raw input is tried first so that power users keep AND / OR / NEAR and prefix
    search. If SQLite rejects it, every token is quoted instead, which can never
    be a syntax error.
    """
    tokens = [t for t in re.split(r"[^\w*]+", text) if t]
    fallback = " ".join('"%s"' % t.replace('"', "") for t in tokens)
    return text.strip(), fallback


def fts_prefix_expression(text):
    """Her tokene sonek yildizi ekler; zaten yildizli olan birakilir.

    NEDEN. FTS5 TOKEN esler, bu yuzden "toluene" sorgusu
    "4-toluenesulfonate" substratini BULMUYOR: o metin "4" ve
    "toluenesulfonate" olarak tokenleniyor. Bir kimyasal adi arayan okuyucu
    turevlerini de gormek ister ve statik sitedeki arama (alt dizgi esleseme
    yapiyor) onlari zaten goruyordu. Olculdu: "toluene" uygulamada 4,
    yayinlanan sitede 660 giris donuyordu -- iki arama ayni veriden farkli
    cevaplar veriyordu. Onek sorgusu bu farki kapatiyor: 659.
    """
    tokens = [t for t in re.split(r"[^\w*]+", text) if t]
    if not tokens:
        return ""
    return " ".join(t if t.endswith("*") else '"%s"*' % t.replace('"', "")
                    for t in tokens)


def search_query(con, q, filters, flags, page=1, size=PAGE_SIZE):
    """Full-text search with facet filters. Returns (total, rows, facet_counts)."""
    where, params = [], []
    for field, _ in FACETS:
        value = filters.get(field)
        if value:
            where.append(f"s.{field} = ?")
            params.append(value)
    for name, column, _ in FLAGS:
        if flags.get(name):
            where.append(f"s.{column} = 1")
    clause = (" AND " + " AND ".join(where)) if where else ""

    text_query = bool(q and q.strip() and q.strip() != "*")
    joins = ("FROM ro_fts f JOIN ro_search s ON s.rowid = f.rowid "
             if text_query else "FROM ro_search s ") + (
        "JOIN ro r ON r.candidate_id = s.candidate_id "
        "LEFT JOIN ro_subfamily rs ON rs.candidate_id = s.candidate_id "
        "LEFT JOIN operon o ON o.candidate_id = s.candidate_id")
    cond = ("WHERE ro_fts MATCH ?" + clause) if text_query else (
        ("WHERE " + " AND ".join(where)) if where else "")
    order = "ORDER BY bm25(ro_fts)" if text_query else "ORDER BY s.organism, s.candidate_id"

    def run(expr):
        head = ([expr] + params) if text_query else list(params)
        total = con.execute(f"SELECT COUNT(*) {joins} {cond}", head).fetchone()[0]
        rows = con.execute(f"""
            SELECT s.candidate_id, s.protein_id, s.organism, s.product, s.cluster ro_cluster,
                   s.leaf_id, s.tier, s.domain, s.is_plasmid, s.family,
                   r.locus_tag, r.hmm_score, r.model_coverage,
                   rs.assignment_class, rs.core_identity, o.completeness, o.layout
            {joins} {cond} {order} LIMIT ? OFFSET ?""",
            head + [size, (page - 1) * size]).fetchall()
        facets = {}
        for field, _ in FACETS:
            glue = "AND" if cond else "WHERE"
            facets[field] = con.execute(
                f"SELECT s.{field} v, COUNT(*) n {joins} {cond} "
                f"{glue} s.{field} IS NOT NULL AND s.{field} != '' "
                f"GROUP BY 1 ORDER BY n DESC LIMIT 12", head).fetchall()
        return total, rows, facets

    raw, fallback = fts_expression(q)
    result = None
    try:
        result = run(raw)
    except sqlite3.OperationalError:
        if fallback:
            try:
                result = run(fallback)
            except sqlite3.OperationalError:
                result = None
    if result is None:
        return 0, [], {}, False

    # Tam token eslesemesi AZ sonuc verdiyse onek eslesemesiyle genisletilir ve
    # bu durum cagirana BILDIRILIR: sayfa "genisletildi" demek zorunda, yoksa
    # okuyucu neden "toluene" aramasinda "toluenesulfonate" gordugunu bilemez.
    if text_query and result[0] < WIDEN_BELOW and "*" not in q:
        prefix = fts_prefix_expression(q)
        if prefix:
            try:
                widened = run(prefix)
            except sqlite3.OperationalError:
                widened = None
            if widened and widened[0] > result[0]:
                return widened[0], widened[1], widened[2], True
    return result[0], result[1], result[2], False


@app.get("/search", response_class=HTMLResponse)
def search(request: Request, q: str = "", page: int = 1, tier: Optional[str] = None,
           domain: Optional[str] = None, family: Optional[str] = None,
           cluster: Optional[str] = None, plasmid: Optional[str] = None,
           partner: Optional[str] = None, regulator: Optional[str] = None):
    con = connect()
    try:
        q = (q or "").strip()
        filters = {"tier": tier, "domain": domain, "family": family, "cluster": cluster}
        flags = {"plasmid": plasmid, "partner": partner, "regulator": regulator}
        active = {k: v for k, v in filters.items() if v}
        active.update({k: "1" for k, v in flags.items() if v})
        if not q and not active:
            return render(request, "search.html", q="", members=[], total=0, page=1, pages=0,
                          clusters=[], facets={}, filters=filters, flags=flags, active=active)
        if not atlas.table_exists(con, "ro_fts"):
            total, rows = members_query(con, q=q, page=page)
            return render(request, "search.html", q=q, members=rows, total=total, page=page,
                          pages=(total + PAGE_SIZE - 1) // PAGE_SIZE, clusters=[], facets={},
                          filters=filters, flags=flags, active=active, widened=False)
        total, rows, facets, widened = search_query(con, q, filters, flags, page)
        clusters = [c for c in cluster_table(con)
                    if q and (q.lower() in c["cluster"].lower()
                              or q.lower() in (c["substrate"] or "").lower())] if q else []
        return render(request, "search.html", q=q, members=rows, total=total, page=page,
                      pages=(total + PAGE_SIZE - 1) // PAGE_SIZE, clusters=clusters,
                      facets=facets, filters=filters, flags=flags, active=active,
                      widened=widened)
    finally:
        con.close()


# ------------------------------------------------------------------ classify
FASTA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUO*\-]+$", re.I)


def parse_fasta(text):
    seqs = OrderedDict()
    name = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0] if line[1:].strip() else f"seq{len(seqs)+1}"
            base, k = name, 1
            while name in seqs:
                k += 1
                name = f"{base}_{k}"
            seqs[name] = []
        else:
            if name is None:
                name = "seq1"
                seqs[name] = []
            seqs[name].append(line.replace(" ", ""))
    out = OrderedDict()
    for k, chunks in seqs.items():
        s = "".join(chunks).upper().replace("*", "").replace("-", "")
        if not s:
            continue
        if not FASTA_RE.match(s):
            raise ValueError(f"'{k}' contains characters that are not amino acids")
        out[k] = s
    return out


def classify_sequences(seqs, cpu=2):
    """Run the pipeline's two tests on user sequences. Returns list of dicts."""
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "q.fa")
        with open(fa, "w") as fh:
            for k, s in seqs.items():
                fh.write(f">{k}\n{s}\n")
        dom = os.path.join(tmp, "q.dom")
        sto = os.path.join(tmp, "q.sto")
        subprocess.run(["hmmsearch", "--domtblout", dom, "--noali", "--cpu", str(cpu),
                        "-E", "1e-5", HMM_PATH, fa], check=True, stdout=subprocess.DEVNULL,
                       timeout=300)
        subprocess.run(["hmmalign", "--trim", "--amino", "-o", sto, MOTIF_HMM_PATH, fa],
                       check=True, stdout=subprocess.DEVNULL, timeout=300)
        hits = parse_domtblout(dom)
        aligned = read_stockholm_matchcols(sto)
    results = []
    for name, seq in seqs.items():
        res, status = classify_target(hits.get(name, [])) if name in hits else (None, "no_hit")
        motifs = classify_motifs(aligned[name]) if name in aligned else {}
        # top 3 alternative models for context
        alts = sorted(hits.get(name, []), key=lambda h: -h["score"])[:3]
        verdict = "not a Rieske oxygenase alpha"
        if status == "RO_alpha" and motifs.get("is_RO_alpha_motif"):
            verdict = "Rieske oxygenase alpha subunit"
        elif status == "RO_alpha":
            verdict = "coverage passed, catalytic motif incomplete"
        elif status == "fragment":
            verdict = "fragment (too short)"
        elif motifs.get("rieske_intact") and not motifs.get("catalytic_intact"):
            verdict = "Rieske domain only (ferredoxin / ISP type)"
        results.append({
            "name": name, "length": len(seq), "status": status, "verdict": verdict,
            "cluster": res["PredictedCluster"] if res and status == "RO_alpha" else (
                res["BestModelIfRejected"] if res else None),
            "score": res["Score"] if res else None,
            "evalue": res["E-value"] if res else None,
            "coverage": res["ModelCoverage"] if res else None,
            "rieske": f"{motifs.get('rieske_sites_found', 0)}/4",
            "catalytic": f"{motifs.get('catalytic_sites_found', 0)}/3",
            "bridging": motifs.get("bridging_asp"),
            "missing": motifs.get("missing_sites", ""),
            "alts": [(h["query_name"], h["score"], round(h["model_coverage"], 2)) for h in alts],
            "confirmed": status == "RO_alpha" and bool(motifs.get("is_RO_alpha_motif")),
        })
    return results


@app.get("/classify", response_class=HTMLResponse)
def classify_form(request: Request):
    return render(request, "classify.html", results=None, error=None, text="")


@app.post("/classify", response_class=HTMLResponse)
def classify_submit(request: Request, sequences: str = Form("")):
    text = sequences or ""
    error, results = None, None
    try:
        if len(text) > MAX_CLASSIFY_CHARS:
            raise ValueError(f"input too large (max {MAX_CLASSIFY_CHARS:,} characters)")
        seqs = parse_fasta(text)
        if not seqs:
            raise ValueError("no sequences found")
        if len(seqs) > MAX_CLASSIFY_SEQS:
            raise ValueError(f"too many sequences (max {MAX_CLASSIFY_SEQS})")
        results = classify_sequences(seqs)
    except ValueError as exc:
        error = str(exc)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as exc:
        error = f"classifier failed: {exc}"
    return render(request, "classify.html", results=results, error=error, text=text)


@app.get("/about", response_class=HTMLResponse)
def about(request: Request):
    con = connect()
    try:
        # Anlik goruntunun BUGUNKU durumu olcumden gelir. Metin "redundant
        # entries were deleted in between" diyordu ama kac tanesinin
        # silindigini soylemiyordu; 2026-11-01'de 73.699 accession'in
        # tamami UniProt'a soruldu.
        import json as _json                           # noqa: PLC0415
        fresh_path = apath("protein_freshness.json")
        fresh = (_json.load(open(fresh_path))
                 if os.path.exists(fresh_path) else None)
        link = None
        if fresh:
            link = dict(con.execute(
                """SELECT CASE WHEN s.entry_type = 'Inactive' THEN 'dead'
                                ELSE 'alive' END || '_' ||
                          CASE WHEN p.linked_candidate_id IS NULL
                                 OR p.linked_candidate_id = '' THEN 'unlinked'
                               ELSE 'linked' END,
                          COUNT(*)
                     FROM ro_protein p
                     JOIN ro_protein_status s USING(accession)
                    GROUP BY 1"""))
        return render(request, "about.html", totals=totals(con),
                      fresh=fresh, freshlink=link,
                      #  Model sabitlemesi: hangi HMM surumu yayimlanan
                      #  domain cagrilarini uretti (ROADMAP 135.5).
                      pins=atlas.read_json(apath("hmm_pins.json")),
                      #  Kapsama acigi: dis evrene gore ne kadarini
                      #  kapsiyoruz, takson ekseninde (ROADMAP 135.9).
                      coverage=atlas.read_json(apath("coverage_audit.json")))
    finally:
        con.close()


# ------------------------------------------------------------------ downloads
def fasta_response(rows, filename):
    buf = io.StringIO()
    for r in rows:
        buf.write(f">{r['candidate_id']} {r['protein_id'] or ''} cluster={r['ro_cluster']} "
                  f"organism=\"{r['organism']}\"\n{r['sequence']}\n")
    return PlainTextResponse(buf.getvalue(), headers={
        "Content-Disposition": f'attachment; filename="{filename}"'})


def csv_response(rows, filename, fieldnames=None):
    """Satir yoksa BILE baslik satiri yazilir.

    NEDEN. Eskiden `if rows:` korumasi yuzunden bos sonuc BOS BIR DOSYA
    uretiyordu. 2026-10-31'de olculdu: uyesi olmayan 8 tipin CSV'si
    (CarAa, IpBAa, BpdC1, XylC1, NahAc, DntAc, NinA, PbaA1) hic baslik
    tasimiyordu, yani indiren kisi sifir bayt bir dosya aliyor ve bunun
    "uye yok" mu yoksa "indirme bozuldu" mu oldugunu ayirt edemiyordu.
    Bos bir tablo bir SONUCTUR; bos bir dosya bir ARIZA gibi okunur.

    fieldnames verilmezse ilk satirdan turetilir; satir da yoksa
    ValueError yerine bos dosya doner (cagiran taraf kolon adlarini
    bilmiyorsa yapilabilecek baska bir sey yok).
    """
    buf = io.StringIO()
    names = fieldnames or (list(rows[0].keys()) if rows else None)
    if names:
        w = csv.DictWriter(buf, fieldnames=names)
        w.writeheader()
        for r in rows:
            w.writerow(dict(r))
    return Response(buf.getvalue(), media_type="text/csv", headers={
        "Content-Disposition": f'attachment; filename="{filename}"'})


def query_csv(con, sql, params=(), filename="rows.csv"):
    """Sorguyu kosar ve kolon adlarini CURSOR'DAN alir.

    Boylece sonuc bos olsa bile baslik satiri yazilabilir -- satirin
    kendisinden turetilen bir kolon listesi bos sonucta yoktur.
    """
    cur = con.execute(sql, params)
    names = [d[0] for d in cur.description]
    return csv_response(cur.fetchall(), filename, fieldnames=names)


MEMBER_COLS = """r.candidate_id, r.protein_id, r.locus_tag, r.gene, r.product, r.ro_cluster, r.ro_group,
       r.hmm_score, r.model_coverage, r.hmm_evalue, r.nucleotide_id, r.start, r.end, r.strand,
       p.organism, p.taxonomy, p.is_plasmid, rl.leaf_id, rs.assignment_class, rs.core_identity,
       o.completeness operon_completeness, o.has_beta, o.has_ferredoxin, o.has_reductase, o.layout"""
MEMBER_FROM = """FROM ro r JOIN replicon p USING(nucleotide_id)
       LEFT JOIN ro_leaf rl ON rl.candidate_id=r.candidate_id
       LEFT JOIN ro_subfamily rs ON rs.candidate_id=r.candidate_id
       LEFT JOIN operon o ON o.candidate_id=r.candidate_id"""


@app.get("/download/cluster/{cluster}.fasta")
def dl_cluster_fasta(cluster: str):
    con = connect()
    try:
        rows = con.execute(f"SELECT r.candidate_id, r.protein_id, r.ro_cluster, r.sequence, p.organism "
                           f"FROM ro r JOIN replicon p USING(nucleotide_id) "
                           f"WHERE r.is_confirmed=1 AND r.ro_cluster=?", (cluster,)).fetchall()
        return fasta_response(rows, f"{cluster}.fasta")
    finally:
        con.close()


@app.get("/download/cluster/{cluster}.csv")
def dl_cluster_csv(cluster: str):
    con = connect()
    try:
        return query_csv(
            con,
            f"SELECT {MEMBER_COLS} {MEMBER_FROM} "
            f"WHERE r.is_confirmed=1 AND r.ro_cluster=?",
            (cluster,), f"{cluster}.csv")
    finally:
        con.close()


@app.get("/download/leaf/{leaf_id}.fasta")
def dl_leaf_fasta(leaf_id: str):
    leaf_id = leaf_unslug(leaf_id)
    con = connect()
    try:
        rows = con.execute("SELECT r.candidate_id, r.protein_id, r.ro_cluster, r.sequence, p.organism "
                           "FROM ro r JOIN replicon p USING(nucleotide_id) JOIN ro_leaf rl USING(candidate_id) "
                           "WHERE r.is_confirmed=1 AND rl.leaf_id=?", (leaf_id,)).fetchall()
        return fasta_response(rows, f"{leaf_id.replace('#', '_')}.fasta")
    finally:
        con.close()


@app.get("/download/all_confirmed.fasta")
def dl_all_fasta():
    con = connect()
    try:
        rows = con.execute("SELECT r.candidate_id, r.protein_id, r.ro_cluster, r.sequence, p.organism "
                           "FROM ro r JOIN replicon p USING(nucleotide_id) WHERE r.is_confirmed=1").fetchall()
        return fasta_response(rows, "roar_confirmed.fasta")
    finally:
        con.close()


@app.get("/download/all_confirmed.csv")
def dl_all_csv():
    con = connect()
    try:
        rows = con.execute(f"SELECT {MEMBER_COLS} {MEMBER_FROM} WHERE r.is_confirmed=1").fetchall()
        return csv_response(rows, "roar_confirmed.csv")
    finally:
        con.close()


def union_rows(con):
    """`ro_union` + her dizinin TUM kaynak kimlikleri, tek satirda.

    Birim DIZI oldugu icin bir satir birden fazla accession tasiyabilir;
    kaynaklar `kaynak:accession:durum` bicimiyle boru isaretiyle birlesir.
    Indiren kisi boylece ayni dizinin UniProt'ta silinmis, NCBI genomunda
    onaylanmis olabilecegini TEK satirda gorur -- iki ayri dosyayi
    elle eslestirmek zorunda kalmaz.
    """
    sources = {}
    for h, src, acc, st in con.execute(
            "SELECT seq_md5, source, accession, status FROM ro_union_source "
            "ORDER BY source, accession"):
        sources.setdefault(h, []).append(f"{src}:{acc}:{st}")
    rows = []
    for r in con.execute(
            "SELECT seq_md5, length, n_sources, has_genomic_context, "
            "passes_gate, rieske_intact, triad_sites_found, ro_cluster, "
            "ro_group, genomic_candidate_id, organism, first_source, sequence "
            "FROM ro_union ORDER BY has_genomic_context DESC, passes_gate DESC, "
            "ro_cluster, seq_md5"):
        d = dict(r)
        d["sources"] = "|".join(sources.get(d["seq_md5"], ()))
        rows.append(d)
    return rows


@app.get("/download/all_union.fasta")
def dl_union_fasta():
    """Birlesimin TAMAMI: kapiyi gecmeyen adaylar DAHIL.

    `all_confirmed.fasta` kesisimi verir (genomik baglami olan VE kapiyi
    gecen). Bu dosya birlesimi verir. Ikisi ayri sorunun yanitidir ve
    ikisi de yayinda kalir; hangisini indirdigini basliktan ayirt
    edebilmek icin her baslik ctx ve gate alanlarini tasiyor.
    """
    con = connect()
    try:
        buf = io.StringIO()
        for d in union_rows(con):
            buf.write(
                '>%s cluster=%s group=%s ctx=%d gate=%d sources=%s '
                'organism="%s"\n%s\n' % (
                    d["seq_md5"], d["ro_cluster"] or "-", d["ro_group"] or "-",
                    d["has_genomic_context"], d["passes_gate"],
                    d["sources"] or "-", d["organism"] or "", d["sequence"]))
        return PlainTextResponse(buf.getvalue(), headers={
            "Content-Disposition": 'attachment; filename="roar_union.fasta"'})
    finally:
        con.close()


@app.get("/download/all_union.csv")
def dl_union_csv():
    con = connect()
    try:
        rows = union_rows(con)
        for d in rows:
            d.pop("sequence", None)
        return csv_response(rows, "roar_union.csv", fieldnames=[
            "seq_md5", "length", "n_sources", "has_genomic_context",
            "passes_gate", "rieske_intact", "triad_sites_found",
            "ro_cluster", "ro_group", "genomic_candidate_id",
            "organism", "first_source", "sources"])
    finally:
        con.close()


@app.get("/download/predictions/{name}")
def dl_prediction_table(name: str):
    """Tahmin tablolari OLDUGU GIBI servis edilir.

    Bunlar veritabanindan uretilmez: predictions_build.py'nin yazdigi dosyalar
    dogrudan verilir. Boylece sayfadaki sayilarla indirilen dosya ayni kaynaktan
    gelir ve ayrisamaz.
    """
    if name not in ("archaeal.tsv", "fungal.tsv", "predictions.json"):
        raise HTTPException(404, "unknown prediction table")
    path = os.path.join(PARENT, "predictions", name)
    if not os.path.exists(path):
        raise HTTPException(404, f"{name} not found (run predictions_build.py)")
    media = ("application/json" if name.endswith(".json")
             else "text/tab-separated-values")
    with open(path, "rb") as fh:
        body = fh.read()
    return Response(body, media_type=media, headers={
        "Content-Disposition": f'attachment; filename="roar_{name}"'})


@app.get("/download/operons.csv")
def dl_operons():
    con = connect()
    try:
        rows = con.execute("""SELECT og.candidate_id, r.ro_cluster, p.organism, og.position, og.component,
                                     og.category, og.product, og.gap_bp, og.protein_key
                              FROM operon_gene og JOIN ro r USING(candidate_id)
                              JOIN replicon p ON p.nucleotide_id=r.nucleotide_id
                              ORDER BY og.candidate_id, og.position""").fetchall()
        return csv_response(rows, "roar_operons.csv")
    finally:
        con.close()


@app.get("/download/ro/{ro_slug}.fasta")
def dl_ro_fasta(ro_slug: str):
    con = connect()
    try:
        cid = unslug(ro_slug)
        rows = con.execute("SELECT r.candidate_id, r.protein_id, r.ro_cluster, r.sequence, p.organism "
                           "FROM ro r JOIN replicon p USING(nucleotide_id) WHERE r.candidate_id=?",
                           (cid,)).fetchall()
        return fasta_response(rows, f"{ro_slug}.fasta")
    finally:
        con.close()


@app.get("/download/neighborhood/{ro_slug}.fasta")
def dl_ro_neighborhood(ro_slug: str):
    con = connect()
    try:
        cid = unslug(ro_slug)
        rows = con.execute("""
            SELECT nb.gene_offset, nb.product, nb.protein_id, np.protein_key, np.translation,
                   nc.component
            FROM neighbor nb
            JOIN neighbor_protein np ON np.protein_key = nb.nucleotide_id||':'||nb.start||'-'||nb.end||':'||nb.strand
            LEFT JOIN neighbor_component nc USING(protein_key)
            WHERE nb.candidate_id=? ORDER BY nb.gene_offset""", (cid,)).fetchall()
        buf = io.StringIO()
        for r in rows:
            buf.write(f">{r['protein_key']} offset={r['gene_offset']} {r['protein_id'] or ''} "
                      f"component={r['component'] or 'none'} \"{r['product']}\"\n{r['translation']}\n")
        return PlainTextResponse(buf.getvalue(), headers={
            "Content-Disposition": f'attachment; filename="{ro_slug}_neighborhood.fasta"'})
    finally:
        con.close()


# ------------------------------------------------------------------ JSON API
@app.get("/api/stats")
def api_stats():
    con = connect()
    try:
        return totals(con)
    finally:
        con.close()


@app.get("/api/clusters")
def api_clusters():
    con = connect()
    try:
        return cluster_table(con)
    finally:
        con.close()


@app.get("/api/cluster/{cluster}")
def api_cluster(cluster: str):
    con = connect()
    try:
        info = cluster_detail(con, cluster)
        if not info:
            raise HTTPException(404, "cluster not found")
        for k in ("domains", "genera", "neighborhood", "regulators", "layouts", "leaves", "classes"):
            info[k] = [dict(r) for r in info[k]]
        return info
    finally:
        con.close()


@app.get("/api/ro/{ro_slug}")
def api_ro(ro_slug: str):
    con = connect()
    try:
        d = ro_detail(con, unslug(ro_slug))
        if not d:
            raise HTTPException(404, "entry not found")
        return {"ro": dict(d["ro"]), "neighbors": [dict(n) for n in d["neighbors"]],
                "operon": dict(d["operon"]) if d["operon"] else None,
                "operon_genes": [dict(g) for g in d["operon_genes"]]}
    finally:
        con.close()


@app.get("/api/search")
def api_search(q: str = Query(..., min_length=2), page: int = 1):
    con = connect()
    try:
        total, rows = members_query(con, q=q, page=page)
        return {"total": total, "page": page, "results": [dict(r) for r in rows]}
    finally:
        con.close()


@app.post("/api/classify")
async def api_classify(request: Request):
    body = await request.body()
    text = body.decode("utf-8", "replace")
    if request.headers.get("content-type", "").startswith("application/json"):
        import json
        text = json.loads(text).get("sequences", "")
    try:
        seqs = parse_fasta(text)
        if not seqs or len(seqs) > MAX_CLASSIFY_SEQS or len(text) > MAX_CLASSIFY_CHARS:
            raise ValueError("provide 1-50 protein sequences in FASTA format")
        return {"results": classify_sequences(seqs)}
    except ValueError as exc:
        return JSONResponse({"error": str(exc)}, status_code=400)


# ===================================================================== ATLAS
import atlas  # noqa: E402
# Uyelik temeli aciklamalari KAYNAGINDA duruyor; sayfa kopyalamaz.
import triad_gate  # noqa: E402

ANALYSIS_DIR = os.environ.get("ROAR_ANALYSIS", os.path.join(PARENT, "analysis_out"))
templates.env.globals.update(
    leaf_stop_reason=atlas.leaf_stop_reason,
    LEAF_STOP_REASONS=atlas.LEAF_STOP_REASONS,
    MEMBERSHIP_BASES=triad_gate.MEMBERSHIP_BASES,
    element_slug=atlas.element_slug,
    hab_label=atlas.HABITAT_LABEL,
    layout_svg=atlas.layout_svg, operon_regulator_svg=atlas.operon_regulator_svg,
    dominant_arch=atlas.dominant_arch,
    reaction_scheme_svg=atlas.reaction_scheme_svg,
    chip_text=atlas.chip_text,
    # Sema aciklamalari SVG'nin icinde degil: orada sabit x konumlarina
    # yazilinca kirpiliyorlardi. Sayfada HTML olarak yaziliyorlar.
    REACTION_SCHEME_NOTE=atlas.REACTION_SCHEME_NOTE,
    gap_histogram=atlas.gap_histogram,
    # Grup ETIKETI tek yerden gelir: 6 ve 7 literaturde yok, numara olarak
    # yayimlanmazlar (bkz. atlas.group_label).
    group_label=atlas.group_label, is_literature_group=atlas.is_literature_group,
    row_label=atlas.row_label,
    LITERATURE_GROUPS=atlas.LITERATURE_GROUPS,
    amedian=atlas.median, GROUP_COLORS=atlas.GROUP_COLORS, TIER_COLORS=atlas.TIER_COLORS, TIERS=atlas.TIERS,
    TIER_LABEL=atlas.TIER_LABEL, TIER_MEANING=atlas.TIER_MEANING,
    COMPONENT_COLORS=COMPONENT_COLORS,
    entry_href=entry_href,
    # Ham taban onek. u('') KULLANILAMAZ: statik kipte u() sonuna ".html"
    # ekliyor ve taban "/pro-sim-blast.html" oluyordu -- giris sayfasi bu
    # yuzden butun veri isteklerinde 404 aldi (2026-10-19'da olculdu).
    BASE_PREFIX=BASE,
    # Kuratorlu referans sayisi VERIDEN gelir, sablonlara elle yazilmaz.
    # 2026-10-18: on iki sablonda "71" sabit yaziliydi ve referans seti 75'e
    # cikinca hepsi yanlisa dustu. Ayni hata motif_stats.py'de de olmustu
    # (orada 71 ve 426 sabitti). Sabit sayi, degisen veriyle sessizce ayrisir.
    N_REFS=len(CHEMISTRY))


def apath(*parts):
    return os.path.join(ANALYSIS_DIR, *parts)


def _overview(con):
    return atlas.per_cluster_overview(con, ECOLOGY)


def _carboxylate_spread():
    """Katalitik karboksilatin yapisal referanslarda hangi kolonlara dustugu.

    BU SAYILAR SAYFADA ELLE YAZILIYDU (2026-11-03'e kadar): "18 of them",
    "432 in 2 more", "the remaining 6", "407-411 and 449-451". Olcum iki
    dosyada ZATEN vardi; yani sayfa, var olan bir olcumu kopyalayarak
    yayimliyordu ve referans seti ya da capalama degistiginde sessizce
    bayatlayacakti -- bu projede ayni bicimde bir sayi daha once bayatladi.

    Duzyazi AYRISTIRILMIYOR. Denetimin `note` alani bantlari Turkce bir
    cumlede tasiyor; buradaki hicbir deger o cumleden okunmuyor:
      testable    : result.testable
      within_one  : references_at_431 + at_432  (431 +-1 penceresi)
      outliers    : testable - within_one
      bands       : column_bands'in LISTE degerli anahtarlari -- "431-432"
                    anahtari bir sayi dizesi tasidigi icin kendiliginden
                    disarida kalir
    Herhangi bir alan eksikse None doner ve sablon o cumleyi yazmaz; eksik
    olcumu tahminle doldurmak, olculmemis bir seyi yayimlamak olurdu.
    """
    audit = atlas.read_json(apath("carboxylate_column_audit.json")) or {}
    clade = atlas.read_json(apath("carboxylate_clade_shift.json")) or {}
    real = audit.get("the_real_column") or {}
    result = audit.get("result") or {}
    bands = [k for k, v in (clade.get("column_bands") or {}).items()
             if isinstance(v, list)]
    try:
        testable = int(result["testable"])
        at_col = int(real["references_at_431"])
        at_next = int(real["at_432"])
        column = int(real["column"])
    except (KeyError, TypeError, ValueError):
        return None
    if not bands:
        return None
    within_one = at_col + at_next

    # MESAFE ARALIGI: sayfa "in all 26 testable structures" diyordu, OLCUM
    # 20 (2026-11-03). Ham dosyada 6 yapi `column_call_distance = -1`
    # tasiyor -- bu bir mesafe DEGIL, "kolon cagrisinin yapida karsiligi
    # bulunamadi" nobetci degeri; ve o 6 yapi tam da bant aykirilarindan
    # (DdmC, NdmA, NdmB, BphA1, NagGH, ROCH34). Bir yapida (HcaE/8K0A)
    # mononukleer demir hic yok ve alan bile yazilmamis.
    #
    # Nobetci degeri ayiklamak ZORUNLU: naif bir min() sayfaya "-1 A"
    # yazardi. Ayiklanan sayi da GIZLENMIYOR, cunku 20'de olculen bir
    # araligi 26 yapiya atfetmek tam olarak bu projenin tekrar eden hatasi.
    raw = atlas.read_json(apath("carboxylate_column_audit_raw.json")) or []
    dists = [r["column_call_distance"] for r in raw
             if isinstance(r.get("column_call_distance"), (int, float))
             and r["column_call_distance"] > 0]
    sentinel = sum(1 for r in raw
                   if isinstance(r.get("column_call_distance"), (int, float))
                   and r["column_call_distance"] <= 0)
    no_field = sum(1 for r in raw if "column_call_distance" not in r)

    out = {"column": column, "testable": testable, "at_column": at_col,
           "at_next": at_next, "within_one": within_one,
           "outliers": testable - within_one,
           "bands": " and ".join(sorted(bands))}
    if dists:
        out.update({"dist_n": len(dists), "dist_min": min(dists),
                    "dist_max": max(dists), "dist_unmeasured": sentinel,
                    "dist_no_iron": no_field})
    return out


def _bridging_glu_by_type(con, group="5", min_n=20, top=4):
    """Grup 5 tiplerinde koprudeki karboksilatin Glu payi -- OLCULEREK.

    SAYFA BU IKI SAYIYI ELLE YAZIYORDU: "reaching 100 % in GbcA and above
    89 % in Stc2, CntA and BmoA". Olcum (2026-11-03): GbcA %100,0 dogru, ama
    Stc2 %93,4, CntA %90,7 ve **BmoA %87,4** -- yani BmoA 89'un ALTINDA ve
    iddia onun icin yanlis. Ustelik bu paylar uyelik her buyudugunde kayar
    (A5 kurtarmasi tam da grup 5'i buyutuyor), yani elle yazilan bir sayi
    burada kalici olarak yanlis olmaya mahkum.

    min_n: tek haneli uye sayilariyla hesaplanan bir yuzdeyi yayimlamamak
    icin taban; etc_types.py'nin MIN_TYPE_N olcutuyle ayni gerekce.
    """
    rows = con.execute("""
        SELECT m.ro_cluster, COUNT(*) n,
               SUM(CASE WHEN c.bridging_residue='E' THEN 1 ELSE 0 END) e
        FROM ro_carboxylate c JOIN ro_membership m USING(candidate_id)
        -- LIKE yerine on ek karsilastirmasi: '5_%' kalibinda alt tire
        -- LIKE'in joker karakteri ve ESCAPE kacisi gereksiz kirilganlik.
        WHERE substr(m.ro_cluster, 1, length(? ) + 1) = ? || '_'
        GROUP BY 1 HAVING COUNT(*) >= ?
        ORDER BY 1.0 * SUM(CASE WHEN c.bridging_residue='E' THEN 1 ELSE 0 END)
                 / COUNT(*) DESC
    """, (group, group, min_n)).fetchall()
    return [{"cluster": r[0], "name": r[0].split("_", 2)[-1],
             "n": r[1], "share": r[2] / r[1]} for r in rows[:top]]


@app.get("/atlas", response_class=HTMLResponse)
def atlas_home(request: Request):
    con = connect()
    try:
        return render(request, "atlas_home.html", over=_overview(con), totals=totals(con),
                      stats=atlas.read_json(apath("stats.json")),
                      carbox_spread=_carboxylate_spread(),
                      glu_top=_bridging_glu_by_type(con),
                      motif=atlas.read_json(apath("motif_stats.json")))
    finally:
        con.close()


@app.get("/atlas/groups", response_class=HTMLResponse)
def atlas_groups(request: Request):
    """Grup numaralarinin NE OLDUGUNU soyleyen sayfa.

    Bu sayfa 2026-11-01'de eklendi, cunku site "group 3" yaziyor ve grup 3'un
    ne oldugunu soyleyen hicbir sayfa YOKTU -- groups.yaml hic yayinlanmiyordu.
    Batie'nin bes sinifiyla karistirilma riski tam buradan doguyor: tanimsiz
    bir numara, okuyucunun bildigi baska bir bes'li semaya yapisir.
    """
    con = connect()
    try:
        gv = atlas.group_view(con, os.path.join(PARENT, "groups.yaml"))
        if not gv:
            raise HTTPException(404, "groups.yaml not readable")
        return render(request, "atlas_groups.html", gv=gv, page_id="groups")
    finally:
        con.close()


@app.get("/atlas/phylogeny", response_class=HTMLResponse)
def atlas_phylogeny(request: Request):
    if not os.path.exists(apath("tree_all.nwk")):
        raise HTTPException(404, "tree not built (run build_phylogeny.py)")
    con = connect()
    try:
        newick = open(apath("tree_all.nwk")).read()
        svg, n = atlas.render_tree_svg(newick,
                                       atlas.tip_info_from_db(con, u, ECOLOGY, short_names()),
                                       width=1040, row_h=11, label_w=380)
        # Yontem beyani OLCUMDEN gelir. Altyazi bir donem "FastTree under the
        # JTT model" diye sabit yaziyordu; fiilen kosan model LG idi.
        meta_path = apath("tree_meta.json")
        import json as _json
        tree_meta = (_json.load(open(meta_path)) if os.path.exists(meta_path) else None)
        # Gruplarin dal olusturup olusturmadigi da OLCUMDEN gelir. Sayfa bir
        # donem "this tree recovers the historical RO groups from sequence
        # alone" diye yaziyordu; olcum bunu yalnizca grup 5 ve 6 icin
        # dogruladi, grup 1-4 hicbir kesimde dal olusturmuyor.
        mono_path = apath("group_monophyly.json")
        mono = (_json.load(open(mono_path)) if os.path.exists(mono_path) else None)
        # BOOTSTRAP ayri bir olcumdur: SH-benzeri yerel destek hizalama
        # kolonlarini yeniden orneklemez, bootstrap orneklier. Bir grubun
        # klad OLUP OLMADIGI sorusu ancak bununla yanitlanir.
        boot = {}
        for which in ("references", "confident"):
            path = apath(f"bootstrap_{which}.json")
            if os.path.exists(path):
                boot[which] = _json.load(open(path))
        return render(request, "atlas_phylogeny.html", svg=svg, n_tips=n,
                      boot=boot or None,
                      cover=atlas.representative_coverage(con),
                      tree_meta=tree_meta, mono=mono,
                      mono_sum=atlas.monophyly_summary(mono),
                      support=atlas.newick_support_summary(newick))
    finally:
        con.close()


@functools.lru_cache(maxsize=1)
def _rep_coverage():
    """Temsilci kapsamasi; ag rotasi kendi baglantisini acmiyor."""
    con = connect()
    try:
        return atlas.representative_coverage(con)
    finally:
        con.close()


@app.get("/atlas/network", response_class=HTMLResponse)
def atlas_network(request: Request, min_identity: float = 35.0):
    return render(request, "atlas_network.html",
                  data=atlas.ssn(apath("ssn_nodes.csv"), apath("ssn_edges.csv"), min_identity),
                  matrix=atlas.identity_matrix(apath("cluster_identity_matrix.csv")),
                  refpairs=atlas.read_csv(apath("reference_pairs.csv"))[:40],
                  names=short_names(), min_identity=min_identity,
                  cover=_rep_coverage())


@app.get("/atlas/taxonomy", response_class=HTMLResponse)
def atlas_taxonomy(request: Request):
    con = connect()
    try:
        over = _overview(con)
        counted = Counter()
        for d in over:
            counted.update(d["phyla"])
        return render(request, "atlas_taxonomy.html", over=over,
                      top_phyla=[p for p, _ in counted.most_common(12)],
                      all_phyla=counted.most_common(20),
                      tree=atlas.taxonomy_tree(con),
                      gtax=atlas.read_json(apath("taxonomy_ids.json")),
                      gsens=atlas.read_json(apath("promoter_architecture.json")))
    finally:
        con.close()


@app.get("/atlas/ecology", response_class=HTMLResponse)
def atlas_ecology(request: Request):
    con = connect()
    try:
        hab = atlas.read_json(apath("habitat.json"))
        return render(request, "atlas_ecology.html", over=_overview(con),
                      stats=atlas.read_json(apath("stats.json")), habitat=hab,
                      hab_rows=atlas.habitat_rows(hab),
                      hab_label=atlas.HABITAT_LABEL,
                      kingdom_rows=atlas.host_kingdom_rows(hab),
                      enrich_family=atlas.habitat_enrichment(
                          hab, "chemical_family_habitat_profile", "contaminated_industrial"),
                      enrich_class=atlas.habitat_enrichment(
                          hab, "substrate_class_habitat_profile", "contaminated_industrial"),
                      enrich_mining=atlas.habitat_enrichment(
                          hab, "chemical_family_habitat_profile", "mining_acid_drainage"))
    finally:
        con.close()


@app.get("/atlas/operons", response_class=HTMLResponse)
def atlas_operons(request: Request):
    con = connect()
    try:
        profiles = con.execute("SELECT etc_profile, COUNT(*) n FROM ro_etc GROUP BY 1 "
                               "ORDER BY n DESC LIMIT 14").fetchall() \
            if atlas.table_exists(con, "ro_etc") else []
        by_group = con.execute(
            "SELECT r.ro_group grp, t.reductase_type red, t.ferredoxin_type fd, t.has_beta beta, "
            "COUNT(*) n FROM ro_etc t JOIN ro r USING(candidate_id) GROUP BY 1,2,3,4").fetchall() \
            if atlas.table_exists(con, "ro_etc") else []
        layouts = {}
        for cl, lay, n in con.execute(
                "SELECT r.ro_cluster, o.layout, COUNT(*) n FROM operon o JOIN ro r USING(candidate_id) "
                "WHERE r.is_confirmed=1 GROUP BY 1,2 ORDER BY 1, n DESC"):
            layouts.setdefault(cl, []).append((lay, n))
        reg_rows = {r["cluster"]: r for r in atlas.read_csv(apath("regulation_by_cluster.csv"))}
        return render(request, "atlas_operons.html", over=_overview(con), profiles=profiles,
                      by_group=[list(r) for r in by_group], layouts=layouts, reg_rows=reg_rows,
                      validation=atlas.read_json(apath("operon_validation.json")),
                      stats=atlas.read_json(apath("stats.json")))
    finally:
        con.close()


@app.get("/atlas/regulation", response_class=HTMLResponse)
def atlas_regulation(request: Request):
    con = connect()
    try:
        if not atlas.table_exists(con, "ro_regulation"):
            raise HTTPException(404, "regulation table not built (run analyze_regulation.py)")
        arch = con.execute("SELECT architecture, COUNT(*) n FROM ro_regulation "
                           "GROUP BY 1 ORDER BY n DESC").fetchall()
        fams = con.execute("SELECT upstream_family f, COUNT(*) n FROM ro_regulation "
                           "WHERE architecture='divergent_regulator' AND upstream_family IS NOT NULL "
                           "GROUP BY 1 ORDER BY n DESC").fetchall()
        gaps = {a: [r[0] for r in con.execute(
            "SELECT intergenic_bp FROM ro_regulation WHERE architecture=? AND intergenic_bp "
            "IS NOT NULL AND intergenic_bp <= 2000", (a,))] for a in
            ("divergent_regulator", "codirectional_regulator", "divergent_other", "codirectional_other")}
        products = con.execute(
            "SELECT upstream_product p, COUNT(*) n FROM ro_regulation "
            "WHERE architecture='divergent_regulator' GROUP BY 1 ORDER BY n DESC LIMIT 15").fetchall()
        per_cluster = atlas.read_csv(apath("regulation_by_cluster.csv"))
        # UZLASMA OLCULUR, ANLATILMAZ. Sayfada iki tablo ayni sayiyi veriyor:
        # mimari tablosu ro_regulation'i dogrudan sayiyor, tip tablosu ise
        # regulation_by_cluster.csv'den geliyor. Okuyucu 68 satiri toplamadan
        # ikisinin ayni oldugunu goremiyordu ve toplamayi kendi yapmak zorunda
        # kaliyordu -- yanlis bir ucuncu sayi tam boyle dogar. Toplam satiri
        # artik basiliyor ve validate_curation.py esitligi zorluyor.
        reg_totals = {k: sum(int(r.get(k) or 0) for r in per_cluster) for k in
                      ("n", "divergent_regulator", "codirectional_regulator",
                       "divergent_other", "codirectional_other", "unknown")}
        reg_totals["db_divergent_regulator"] = dict(arch).get("divergent_regulator", 0)
        return render(request, "atlas_regulation.html", arch=arch, fams=fams, gaps=gaps,
                      products=products, per_cluster=per_cluster, reg_totals=reg_totals,
                      over=_overview(con))
    finally:
        con.close()


@app.get("/atlas/evidence", response_class=HTMLResponse)
def atlas_evidence(request: Request):
    con = connect()
    try:
        over = _overview(con)
        tot = Counter()
        for d in over:
            tot.update(d["tiers"])
        mism = con.execute(
            "SELECT r.ro_cluster hmm, e.nearest_ref nearest, COUNT(*) n, AVG(e.ref_identity) ident "
            "FROM ro_evidence e JOIN ro r USING(candidate_id) "
            "WHERE e.same_as_hmm=0 AND e.nearest_ref IS NOT NULL GROUP BY 1,2 "
            "ORDER BY n DESC LIMIT 15").fetchall() if atlas.table_exists(con, "ro_evidence") else []
        idents = [r[0] or 0.0 for r in con.execute("SELECT ref_identity FROM ro_evidence")] \
            if atlas.table_exists(con, "ro_evidence") else []
        sweep = [(t, sum(1 for i in idents if i >= t), len(idents))
                 for t in (25, 30, 35, 40, 50, 60, 70, 80, 90, 95)]
        pairs = atlas.read_csv(apath("reference_pairs.csv"))
        diff = [p for p in pairs if p["same_substrate_label"] == "0"]
        same = [p for p in pairs if p["same_substrate_label"] == "1"]
        calib = {
            "max_diff": max((float(p["identity"]) for p in diff), default=None),
            "max_diff_pair": max(diff, key=lambda p: float(p["identity"])) if diff else None,
            "min_same": min((float(p["identity"]) for p in same), default=None),
            "min_same_pair": min(same, key=lambda p: float(p["identity"])) if same else None,
            "n_diff_over_60": sum(1 for p in diff if float(p["identity"]) >= 60),
            "n_diff": len(diff), "n_same": len(same)}
        return render(request, "atlas_evidence.html", over=over, tot=tot, mism=mism,
                      sweep=sweep, calib=calib,
                      xref=atlas.read_json(apath("cross_references.json")),
                      genes={c: (CHEMISTRY.get(c) or {}).get("gene") or c
                             for c in CHEMISTRY},
                      pred=atlas.read_json(apath("substrate_predictability.json")))
    finally:
        con.close()


NOVEL_WHERE = """
    FROM ro r
    JOIN ro_evidence e ON e.candidate_id = r.candidate_id
    JOIN replicon p ON p.nucleotide_id = r.nucleotide_id
    LEFT JOIN ro_leaf rl ON rl.candidate_id = r.candidate_id
    LEFT JOIN leaf l ON l.leaf_id = rl.leaf_id
    LEFT JOIN ro_domain d ON d.candidate_id = r.candidate_id
    LEFT JOIN operon o ON o.candidate_id = r.candidate_id
    LEFT JOIN ro_subfamily s ON s.candidate_id = r.candidate_id
    WHERE r.is_confirmed = 1 AND e.tier = 'novel' AND r.rieske_intact = 1
      AND r.catalytic_intact = 1 AND LENGTH(r.sequence) >= 300 AND l.size >= 3"""


@app.get("/atlas/novel", response_class=HTMLResponse)
def atlas_novel(request: Request):
    con = connect()
    try:
        if not atlas.table_exists(con, "ro_evidence"):
            raise HTTPException(404, "evidence table not built (run evidence_tiers.py)")
        rows = con.execute(f"""
            SELECT r.candidate_id, r.protein_id, r.ro_cluster, r.product, p.organism,
                   p.is_plasmid, e.ref_identity, e.nearest_ref, rl.leaf_id, l.size leaf_size,
                   d.domain, d.euk_group, o.has_beta, o.has_ferredoxin, o.has_reductase,
                   o.layout, s.assignment_class
            {NOVEL_WHERE} ORDER BY l.size DESC, e.ref_identity ASC""").fetchall()
        by_variant = {}
        for r in rows:
            by_variant.setdefault(r["leaf_id"], []).append(r)
        variants = sorted(by_variant.items(), key=lambda kv: -len(kv[1]))
        agreement = con.execute("""
            SELECT e.tier, s.assignment_class, COUNT(*) n
            FROM ro_evidence e JOIN ro_subfamily s USING(candidate_id)
            GROUP BY 1, 2""").fetchall()
        tiers = atlas.TIERS
        classes = ["core", "divergent", "alt_type", "novel_candidate"]
        matrix = {(t, c): 0 for t in tiers for c in classes}
        for r in agreement:
            if (r["tier"], r["assignment_class"]) in matrix:
                matrix[(r["tier"], r["assignment_class"])] = r["n"]
        domains = Counter(r["domain"] or "unassigned" for r in rows)
        live = {r[0] for r in con.execute(
            "SELECT DISTINCT ro_cluster FROM ro WHERE is_confirmed=1")}
        live |= {c for c in CHEMISTRY} | {c for c in ECOLOGY}   # artik hepsinin sayfasi var
        return render(request, "atlas_novel.html", rows=rows, variants=variants,
                      matrix=matrix, tiers=tiers, classes=classes, domains=domains,
                      n_total=len(rows), live_clusters=live,
                      budget=atlas.novelty_budget(con))
    finally:
        con.close()


@app.get("/download/novel_candidates.fasta")
def dl_novel():
    con = connect()
    try:
        rows = con.execute(f"""
            SELECT r.candidate_id, r.protein_id, r.ro_cluster, r.sequence, p.organism,
                   e.ref_identity, e.nearest_ref, rl.leaf_id
            {NOVEL_WHERE} ORDER BY l.size DESC""").fetchall()
        buf = io.StringIO()
        for r in rows:
            buf.write(f">{r['candidate_id']} {r['protein_id'] or ''} nearest={r['nearest_ref']} "
                      f"identity={r['ref_identity']:.1f} variant={r['leaf_id']} "
                      f"organism=\"{r['organism']}\"\n{r['sequence']}\n")
        return PlainTextResponse(buf.getvalue(), headers={
            "Content-Disposition": 'attachment; filename="roar_novel_candidates.fasta"'})
    finally:
        con.close()


@app.get("/atlas/cooccurrence", response_class=HTMLResponse)
def atlas_cooccurrence(request: Request):
    data = atlas.read_json(apath("cooccurrence.json"))
    if not data:
        raise HTTPException(404, "cooccurrence.json not found (run cooccurrence.py)")
    con = connect()
    try:
        chem = {k: v for k, v in CHEMISTRY.items()}
        return render(request, "atlas_cooccurrence.html", data=data, chem=chem,
                      totals=totals(con))
    finally:
        con.close()


@app.get("/atlas/control", response_class=HTMLResponse)
def atlas_control(request: Request):
    """Regulatorler, regulator ayrismasi ve hareketlilik."""
    return render(request, "atlas_control.html",
                  rel=atlas.operon_relations_view(apath("operon_relations.json")),
                  vo=atlas.variant_operon_view(apath("variant_and_operon.json")))


@app.get("/atlas/learning", response_class=HTMLResponse)
def atlas_learning(request: Request):
    """Diziden ve baglamdan ne ongorulebiliyor, ne ongorulemiyor."""
    return render(request, "atlas_learning.html",
                  learn=atlas.learning_view(apath("learning.json")))


@app.get("/atlas/origin", response_class=HTMLResponse)
def atlas_origin(request: Request):
    """Yakin ve uzak akrabalarin nerelerden geldigi, ve kirletici-oncesi ayrilma."""
    return render(request, "atlas_origin.html",
                  eco=atlas.ecological_origin_view(apath("ecological_origin.json")),
                  anc=atlas.ancestry_view(apath("ancestry.json")))


@app.get("/atlas/disagreements", response_class=HTMLResponse)
def atlas_disagreements(request: Request):
    """Yayimlanmis iddialarla olculen arasindaki farklar."""
    dis = atlas.read_json(apath("disagreements.json"))
    if not dis:
        raise HTTPException(404, "disagreements.json not found "
                                 "(run disagreements.py)")
    return render(request, "atlas_disagreements.html", dis=dis)


@app.get("/atlas/predictions", response_class=HTMLResponse)
def atlas_predictions(request: Request):
    """Arkeal ve fungal TAHMINLER -- kuratorlu setin DISINDA bir katman.

    Veri predictions/predictions.json'dan gelir; bu dosya predictions_build.py
    tarafindan uretilir ve types/, chemistry.csv, HMM kutuphanesi ve
    istatistiklerin hicbirine dokunmaz (bkz predictions/README.md).
    """
    path = os.path.join(PARENT, "predictions", "predictions.json")
    pr = atlas.read_json(path)
    if not pr:
        raise HTTPException(404, "predictions/predictions.json not found "
                                 "(run predictions_build.py)")
    return render(request, "atlas_predictions.html", pr=pr)


@app.get("/atlas/elements", response_class=HTMLResponse)
def atlas_elements(request: Request):
    """Duzenleyici aileleri ve IS aileleri: dizin, Sankey, ozgulluk testi."""
    el = atlas.control_elements_view(apath("control_elements.json"))
    if not el:
        raise HTTPException(404, "control_elements.json not found "
                                 "(run control_elements.py)")
    return render(request, "atlas_elements.html", el=el)


@app.get("/element/{family}", response_class=HTMLResponse)
def control_element_page(request: Request, family: str):
    """Tek bir duzenleyici ya da IS ailesi."""
    view = atlas.control_element_view(apath("control_elements.json"), family)
    if not view:
        raise HTTPException(404, f"no regulator or IS family named {family}")
    return render(request, "control_element.html", **view)


@app.get("/atlas/geography", response_class=HTMLResponse)
def atlas_geography(request: Request):
    """Enzimler nerelerde goruldu -- ve bu sorunun ne kadar sorulabilir oldugu."""
    geo = atlas.geography_view(apath("geography.json"))
    if not geo:
        raise HTTPException(404, "geography.json not found (run geography.py)")
    return render(request, "atlas_geography.html", geo=geo)


@app.get("/atlas/structure", response_class=HTMLResponse)
def atlas_structure(request: Request):
    """Aktif bolge: kristal yapilardan ve tahmin modellerinden."""
    return render(request, "atlas_structure.html",
                  site=atlas.active_site_view(apath("active_site.json")),
                  vo=atlas.variant_operon_view(apath("variant_and_operon.json")))


@app.get("/atlas/closeness", response_class=HTMLResponse)
def atlas_closeness(request: Request):
    """Yakinliga gore tabakalanmis sonuclar.

    Siteyi okuyan birinin en mesru sorusu su: bu sayilarin ne kadari kimyasi
    BILINEN uyelerden geliyor. Bu sayfa her iddiayi o eksende yeniden olcuyor.
    """
    strat = atlas.read_json(apath("stratified_stats.json")) or {}
    rows_order = [
        ("all", "every entry", "the single number the rest of the site reports"),
        ("close", "at least 60 % identical to a curated enzyme",
         "the members whose substrate label means something"),
        ("distant", "below 60 % identical", "reaction unknown for all of these"),
        ("tier:characterized", "characterised enzymes only",
         "the reference enzymes and their near-identical relatives"),
        ("tier:close_homolog", "close homologues", "where the rule is at its cleanest"),
        ("tier:family_member", "family members", ""),
        ("tier:distant", "distant relatives", ""),
        ("tier:novel", "no close characterised relative", "the most divergent members"),
        ("class:core", "core of its own subfamily", ""),
        ("class:divergent", "divergent within its subfamily", ""),
        ("class:novel_candidate", "novel candidates", ""),
    ]
    return render(request, "atlas_closeness.html", strat=strat,
                  rows_order=rows_order,
                  origin=strat.get("contaminated_origin"))


@app.get("/atlas/quality", response_class=HTMLResponse)
def atlas_quality(request: Request):
    con = connect()
    try:
        red = atlas.read_json(apath("redundancy.json"))
        return render(request, "atlas_quality.html", red=red,
                      stats=atlas.read_json(apath("stats.json")),
                      motif=atlas.read_json(apath("motif_stats.json")),
                      validation=atlas.read_json(apath("operon_validation.json")),
                      sens=atlas.read_json(apath("threshold_sensitivity.json")),
                      refdup=reference_redundancy(),
                      #  Dis dogrulama: kendi Rieske cagrimiz InterPro'nun
                      #  kendi profillerine karsi (ROADMAP 135.19).
                      ipr=atlas.read_json(apath("interpro_agreement.json")),
                      #  Ligand ligand denetim: asimetrik 139 dizide hangi
                      #  kume ligandi yok, ve motif HAM dizide var mi.
                      ligand=atlas.read_json(
                          apath("rieske_ligand_audit.json")),
                      #  124.3'un kalan yarisi: kutu zenginlesmesi, iki
                      #  null ve uc duzey (ROADMAP 135.22).
                      sigma=atlas.read_json(
                          apath("sigma70_enrichment.json")),
                      downloads=atlas.download_table(ANALYSIS_DIR), totals=totals(con))
    finally:
        con.close()


@app.get("/atlas/statistics", response_class=HTMLResponse)
def atlas_statistics(request: Request):
    data = atlas.read_json(apath("stats.json"))
    if not data:
        raise HTTPException(404, "stats.json not found (run stats_overview.py)")
    # Ogrenme sonuclari da buraya gelir: istatistik sayfasi "bu iliskilerden ne
    # ongorulebilir" sorusunu KISA bir bolumde cevaplar ve ayrintiyi kendi
    # sayfasina birakir.
    # Disaridan dogrulama: yayimlanmis kristal yapilardan gelen DENEYSEL alt
    # birim mimarisi, sitenin genomik baglamdan CIKARDIGI mimariyle yan yana
    # konur. Dosya yoksa bolum sessizce dusur.
    return render(request, "atlas_statistics.html", stats=data,
                  learn=atlas.learning_view(apath("learning.json")),
                  refstruct=atlas.read_json(apath("reference_structures.json")))


@app.get("/download/tree_all.nwk")
def dl_tree():
    if not os.path.exists(apath("tree_all.nwk")):
        raise HTTPException(404, "tree not built")
    return PlainTextResponse(open(apath("tree_all.nwk")).read(), headers={
        "Content-Disposition": 'attachment; filename="roar_tree.nwk"'})


@app.get("/download/analysis/{name}")
def dl_analysis(name: str):
    # Izin listesi atlas.PROVENANCE'tan turetilir; bkz. oradaki not.
    if name not in atlas.ANALYSIS_FILES or not os.path.exists(apath(name)):
        raise HTTPException(404, "not available")
    #  IKILI DOSYA METIN OLARAK OKUNAMAZ. Rota `open(...).read()` ile
    #  metin modunda okuyordu; bir .gz indirmesi UnicodeDecodeError ile
    #  500 verirdi. Tur uzantidan secilir ve icerik her zaman ikili
    #  okunur -- metin dosyalari da ikili olarak dogru gider.
    media = ("application/json" if name.endswith(".json")
             else "application/gzip" if name.endswith(".gz")
             else "text/plain" if name.endswith((".txt", ".gff3"))
             else "text/csv")
    with open(apath(name), "rb") as handle:
        payload = handle.read()
    return Response(payload, media_type=media,
                    headers={"Content-Disposition": f'attachment; filename="{name}"'})


# --- per-page extras -------------------------------------------------------
@functools.lru_cache(maxsize=1)
def reference_redundancy():
    """Referans fazlaligi raporu; tip sayfalarinda uyari olarak gosteriliyor."""
    return atlas.read_json(apath("reference_redundancy.json")) or {}


def cluster_extras(con, cluster):
    x = {"tiers": {}, "etc": [], "tree_svg": None, "nearest": [], "phyla": [],
         "transposon": 0, "regulation": None, "reg_families": []}
    # Bu tipin referansi bir baskasiyla ayni ya da neredeyse ayni mi? Eger
    # oyleyse uyelerin hangi profile dustugu kismen rastgele ve sayfadaki
    # sayinin bunu SOYLEMESI gerekiyor.
    x["redundancy"] = (reference_redundancy().get("affected") or {}).get(cluster)
    # Bu tipin BIRLESIMDEKI uyeligi. Sayfadaki "Members" sayisi kesisimdir
    # (genomik baglami olan VE kapiyi gecen); birlesim bundan buyuk olabilir
    # ve okuyucunun bunu TIP SAYFASINDA gormesi gerekiyor -- yoksa birlesim
    # yalnizca anasayfada duran bir toplam olarak kalir ve gezilemez.
    x["union"] = None
    if con.execute("SELECT name FROM sqlite_master "
                   "WHERE name='ro_union'").fetchone():
        row = con.execute(
            "SELECT COUNT(*), "
            "       SUM(CASE WHEN has_genomic_context=0 THEN 1 ELSE 0 END), "
            "       SUM(COALESCE(passes_gate,0)) "
            "FROM ro_union WHERE ro_cluster = ?", (cluster,)).fetchone()
        if row and row[0]:
            x["union"] = {"total": row[0], "no_context": row[1] or 0,
                          "gated": row[2] or 0}
    # BIRIM: bir giris bir GENOMIK LOKUSTUR, bir protein degil. Kullanici
    # 2026-11-02'de HcaE sayfasinda ayni protein kimliginin tekrar tekrar
    # gectigini gorup "entry duplication var, bug mi" diye sordu. Olculdu:
    # bug yok -- 259 girisin 258'i AYRI genom kaydi ve ayni genom+koordinat
    # ikilisi SIFIR kez tekrarliyor. Tekrar eden sey dizi ve accession:
    # HcaE E. coli enzimi, tek bir dizi 46 susta ayni; ve RefSeq WP_
    # kayitlari TASARIMI GEREGI ayni proteini tasiyan butun genomlar icin
    # TEK accession kullanir. Yani tekrar bir hata degil, birimin sonucu.
    # Sayfa bunu yazmiyordu, o yuzden okuyucu hata sonucuna vardi -- bu
    # projenin tekrar eden sinifi: BIRIM YAZILMAZSA OKUYUCU KENDI SONUCUNU
    # URETIR. Sayilar burada OLCULUYOR, sablonda yazilmiyor.
    row = con.execute(
        "SELECT COUNT(*), COUNT(DISTINCT nucleotide_id), "
        "       COUNT(DISTINCT sequence), "
        "       SUM(CASE WHEN protein_id LIKE 'WP_%' THEN 1 ELSE 0 END) "
        "FROM ro WHERE is_confirmed=1 AND ro_cluster=?", (cluster,)).fetchone()
    top = con.execute(
        "SELECT COUNT(*) n FROM ro WHERE is_confirmed=1 AND ro_cluster=? "
        "GROUP BY sequence ORDER BY n DESC LIMIT 1", (cluster,)).fetchone()
    overlap = con.execute(
        "SELECT COUNT(*) FROM (SELECT nucleotide_id, start, end, COUNT(*) n "
        "FROM ro WHERE is_confirmed=1 AND ro_cluster=? "
        "GROUP BY 1,2,3 HAVING n > 1)", (cluster,)).fetchone()[0]
    x["unit"] = ({"entries": row[0], "records": row[1], "sequences": row[2],
                  "wp": row[3] or 0, "max_copies": (top or [0])[0],
                  "same_locus_twice": overlap} if row and row[0] else None)
    # Bu tipin kendine ait uyesi yoksa SEBEBI yazilmali. Sayfa onceki halinde
    # sessizce bos kaliyordu ve kullanici bunu "sayfa eksik" diye okudu; eksik
    # olan sayfa degil, uye. empty_types.py bu kaydi uretiyor, dosya yoksa
    # None gelir ve sayfa eskisi gibi calisir.
    x["empty_type"] = atlas.empty_type_view(cluster, apath("empty_types.json"))
    # Ters yon: bu tipin altinda sayilan, kendi uyesi olmayan kuratorlu adlar.
    x["absorbs"] = atlas.absorbs_view(cluster, apath("empty_types.json"))
    # Atama ile dizi kimliginin UYUSUP uyusmadigi: tip sayfasindaki uye sayisi
    # bu baglam olmadan yaniltici olabiliyor.
    x["agreement"] = atlas.assignment_agreement(con).get(cluster)
    # Aktif bolge goruntusu: yalnizca yapisi olan tiplerde.
    x["site"] = atlas.active_site_for_type(apath("active_site.json"), cluster)
    if atlas.table_exists(con, "ro_evidence"):
        counted = dict(con.execute(
            "SELECT e.tier, COUNT(*) FROM ro_evidence e JOIN ro r USING(candidate_id) "
            "WHERE r.ro_cluster=? GROUP BY 1", (cluster,)))
        x["tiers"] = {t: counted.get(t, 0) for t in atlas.TIERS}
    if atlas.table_exists(con, "ro_etc"):
        x["etc"] = con.execute(
            "SELECT t.etc_profile p, COUNT(*) n FROM ro_etc t JOIN ro r USING(candidate_id) "
            "WHERE r.ro_cluster=? GROUP BY 1 ORDER BY n DESC LIMIT 8", (cluster,)).fetchall()
    if atlas.table_exists(con, "ro_regulation"):
        x["regulation"] = con.execute(
            "SELECT architecture a, COUNT(*) n FROM ro_regulation g JOIN ro r USING(candidate_id) "
            "WHERE r.ro_cluster=? GROUP BY 1 ORDER BY n DESC", (cluster,)).fetchall()
        x["reg_families"] = con.execute(
            "SELECT upstream_family f, COUNT(*) n FROM ro_regulation g JOIN ro r USING(candidate_id) "
            "WHERE r.ro_cluster=? AND g.architecture='divergent_regulator' AND upstream_family IS NOT NULL "
            "GROUP BY 1 ORDER BY n DESC LIMIT 8", (cluster,)).fetchall()
    x["phyla"] = Counter(atlas.phylum_of(t) for (t,) in con.execute(
        "SELECT p.taxonomy FROM ro r JOIN replicon p USING(nucleotide_id) "
        "WHERE r.is_confirmed=1 AND r.ro_cluster=?", (cluster,))).most_common(8)
    x["transposon"] = con.execute("""
        SELECT COUNT(DISTINCT r.candidate_id) FROM ro r
        JOIN neighbor nb ON nb.candidate_id=r.candidate_id
        JOIN gene_category c ON c.neighbor_id=nb.neighbor_id AND c.method='regex_v1'
             AND c.category='transposon'
        WHERE r.is_confirmed=1 AND r.ro_cluster=?""", (cluster,)).fetchone()[0]
    # "En yakin tipler" REFERANS-REFERANS kimliginden gelir.
    #
    # NEDEN DEGISTI. Eskiden cluster_identity_matrix.csv'den geliyordu, ama o
    # dosya iki tipin benzerligini DEGIL, iki tipe atanmis herhangi iki
    # temsilci arasindaki EN YUKSEK kimligi tasiyor. 2026-10-31'de olculdu:
    # qxyA sayfasi BedC1'e %91,4 benzer diyordu, oysa iki kuratorlu referans
    # birbirine %23 (ortak hizalama) / %27,4 (ikili mafft) benzer. 91,4'u
    # ureten cift iki TEK UYELI yaprak ve ikisi de `distant` kademesinde.
    # O sayi bir benzerlik olcusu degil, atama sinirinda bir CATISMA sinyali;
    # ayri dosyada, ayri adla ve ciftin kademesiyle birlikte yayimlaniyor.
    x["nearest"] = atlas.nearest_types(
        atlas.identity_matrix(apath("cluster_reference_identity.csv")), cluster)
    x["conflicts"] = atlas.assignment_conflicts(
        apath("cluster_assignment_conflicts.csv"), cluster)
    if atlas.table_exists(con, "replicon_source"):
        # "Distinct organisms" kolonu TUR sayisi diye sunuluyor ve basligi
        # "yogun dizilenen bir habitat ekoloji biriktirmeden giris biriktirir"
        # diyerek tam bu artefakti onledigini soyluyor. Ama COUNT(DISTINCT
        # p.organism) GenBank organizma dizgisini sayiyordu, yani SUSLARI:
        # ayni turun on izolati on "organizma" olarak gorunuyordu -- onledigini
        # soyledigi artefaktin kendisi. 2026-10-31'de olculdu, 4_411_XylX
        # rizosfer topragi: 28 yayimlanmis, gercek tur 8; 1_113_CdnA insan
        # klinik: 60'a karsi 32. Projenin kendi tur tanimi
        # isolation_source.species_of() (ilk iki kelime) kullaniliyordu.
        #
        # 2026-11-11: o tanim da BIRAKILDI. "Ilk iki kelime" GenBank
        # dizgisine bagliydi ve "Candidatus Pelagibacter ubique" gibi
        # adlarda tur yerine cins+adlandirma durumu veriyordu. Tur artik
        # NCBI rank='species' atasindan okunuyor (ro_taxonomy.species_map).
        # taxid'i olmayan kayit icin `atlas.species_of` son care olarak
        # kaliyor; su an boyle bir kayit YOK (17.073/17.073).
        hab_rows = con.execute("""
            SELECT s.habitat h, p.organism org, r.nucleotide_id nuc
            FROM ro r JOIN replicon p USING(nucleotide_id)
            JOIN replicon_source s ON s.nucleotide_id = r.nucleotide_id
            WHERE r.is_confirmed=1 AND r.ro_cluster=? AND s.habitat != 'unknown'
            """, (cluster,)).fetchall()
        spmap = ro_taxonomy.species_map(con)
        agg = {}
        for r in hab_rows:
            cell = agg.setdefault(r["h"], {"h": r["h"], "n": 0, "_sp": set()})
            cell["n"] += 1
            sp = spmap.get(r["nuc"]) or atlas.species_of(r["org"])
            if sp:
                cell["_sp"].add(sp)
        x["habitats"] = sorted(
            ({"h": c["h"], "n": c["n"], "species": len(c["_sp"])}
             for c in agg.values()), key=lambda c: -c["n"])[:8]
        x["habitat_coverage"] = con.execute("""
            SELECT SUM(CASE WHEN s.habitat != 'unknown' THEN 1 ELSE 0 END), COUNT(*)
            FROM ro r LEFT JOIN replicon_source s ON s.nucleotide_id = r.nucleotide_id
            WHERE r.is_confirmed=1 AND r.ro_cluster=?""", (cluster,)).fetchone()
    motif = atlas.read_json(apath("motif_stats.json")) or {}
    carbox = {}
    for site, entry in (motif.get("carboxylate_identity") or {}).items():
        row = (entry.get("per_cluster") or {}).get(cluster)
        if row:
            carbox[site] = {**row, "column": entry["column"]}
    x["carboxylate"] = carbox
    reg_rows = {r["cluster"]: r for r in atlas.read_csv(apath("regulation_by_cluster.csv"))}
    x["reg_row"] = reg_rows.get(cluster)
    # En sik operon dizilişi -- REFERANSA BENZEYEN uyeler uzerinden.
    #
    # NEDEN KISITLI. Bir tipin uyeleri profile gore atanir ve cogu UZAK olur:
    # qxyA'nin 164 uyesinden 104'u referansa %25-38 kimlikte. Onlarin operonu
    # "bu enzimin operonu" degil. Butun uyeler alindiginda qxyA'nin en sik
    # dizilişi [alpha] ve payi 60/164 (%37) -- cogunluk bile degil. Yalnizca
    # referansa benzeyenler alindiginda ayni diziliş 51/53 (%96) cikiyor, yani
    # ifade hem dogru hem anlamli hale geliyor. Butun tiplerde olculdu
    # (2026-10-28): ortalama mod payi %25 -> %33, gercek cogunluk saglayan tip
    # 5 -> 9 (57 tip icinde). Kullanici bunu istedi.
    #
    # Yakin uye sayisi cok azsa butun uyelere DUSULUR ve sayfa hangi temele
    # dayandigini SOYLER; sessizce taban degistirmek daha kotu olurdu.
    CLOSE = ("characterized", "close_homolog", "family_member")
    def _layout_rows(close_only):
        q = ("SELECT o.layout, COUNT(*) n FROM operon o JOIN ro r USING(candidate_id) "
             + ("JOIN ro_evidence e USING(candidate_id) " if close_only else "")
             + "WHERE r.is_confirmed=1 AND r.ro_cluster=? "
             + ("AND e.tier IN ('characterized','close_homolog','family_member') "
                if close_only else "")
             + "GROUP BY o.layout ORDER BY n DESC")
        return con.execute(q, (cluster,)).fetchall()

    MIN_CLOSE = 5
    rows = _layout_rows(True)
    total = sum(r["n"] for r in rows)
    x["layout_basis"] = "close"
    if total < MIN_CLOSE:
        rows = _layout_rows(False)
        total = sum(r["n"] for r in rows)
        x["layout_basis"] = "all"
    top_layout = rows[0] if rows else None
    x["layout_total"] = total
    x["layout_distinct"] = len(rows)
    if top_layout:
        x["layout_n"] = top_layout["n"]
        x["layout_share"] = 100.0 * top_layout["n"] / total if total else 0
        x["layout_alpha_only"] = top_layout["layout"] == "[alpha]"
    # Tip geneli sekil, LAYOUT ile AYNI TEMELDE hesaplanir.
    #
    # NEDEN. Eskiden layout yakin uyelerden (characterized/close/family),
    # mimari ise regulation_by_cluster.csv'den yani TUM uyelerden geliyordu ve
    # sekil ikisini birlestiriyordu -- iki ayri paydadan bir resim. qxyA'da
    # olculdu (2026-10-30): 164 uyede divergent_regulator 14, ama yakin 53
    # uyede YALNIZCA 1. Ustelik tum uyelerde en sik aile
    # "regulator_unclassified" oldugu icin sekil TetR'i hic YAZAMIYORDU,
    # oysa yakin uyedeki tek duzenleyici tam olarak TetR (78 bp).
    reg_basis = con.execute(
        "SELECT g.architecture a, g.upstream_family f, g.intergenic_bp bp "
        "FROM ro_regulation g JOIN ro r USING(candidate_id) "
        "JOIN ro_evidence e USING(candidate_id) "
        "WHERE r.is_confirmed=1 AND r.ro_cluster=? "
        "AND e.tier IN ('characterized','close_homolog','family_member')",
        (cluster,)).fetchall() if x["layout_basis"] == "close" else []
    if not reg_basis:
        reg_basis = con.execute(
            "SELECT g.architecture a, g.upstream_family f, g.intergenic_bp bp "
            "FROM ro_regulation g JOIN ro r USING(candidate_id) "
            "WHERE r.is_confirmed=1 AND r.ro_cluster=?", (cluster,)).fetchall()
    if top_layout and reg_basis:
        arch_n = Counter(r["a"] for r in reg_basis if r["a"])
        dominant = arch_n.most_common(1)[0][0] if arch_n else "unknown"
        fam_n = Counter(r["f"] for r in reg_basis
                        if r["a"] == dominant and r["f"]
                        and r["f"] != "regulator_unclassified")
        bps = sorted(r["bp"] for r in reg_basis
                     if r["a"] == dominant and r["bp"] is not None)
        x["canonical_arch"] = dominant
        x["canonical_reg_n"] = arch_n.get(dominant, 0)
        x["canonical_reg_total"] = sum(arch_n.values())
        x["canonical_svg"] = atlas.operon_regulator_svg(top_layout["layout"], {
            "architecture": dominant,
            "upstream_family": fam_n.most_common(1)[0][0] if fam_n else "",
            "intergenic_bp": bps[len(bps) // 2] if bps else None,
            # Bos birakilir: basligi SVG mimariye gore secer, cunku *_other
            # mimaride "regulator" demek olcumun kurmadigi bir iddia olurdu.
            "upstream_product": ""})
        x["canonical_arch_text"] = atlas.ARCH_TEXT.get(
            dominant, dominant.replace("_", " "))
        x["canonical_has_regulator"] = dominant in (
            "divergent_regulator", "codirectional_regulator")

    # --- Referansin KENDI mimarisi, tip geneli istatistikten AYRI ----------
    #
    # NEDEN AYRI. Tip geneli sayi ile literaturun tarif ettigi yapi ayni sey
    # DEGIL. qxyA'da yakin 53 uyenin yalnizca 1'inde divergent duzenleyici var
    # -- karakterize girisin kendisi, TetR, 78 bp -- yani tip geneli hicbir
    # korunmus duzenleyici BILDIRMIYOR. Buna karsin qxyAR operonu referansta
    # gercek ve belgeli (PMID:27792326): ters iplikteki TetR, operonla
    # promoter bolgesi arasinda duruyor, yani tek bir duzenleyici bolge ikisini
    # birlikte kontrol edebiliyor. Bu iki IDDIA ayri, o yuzden ayri gosterilir;
    # birini otekinin yerine koymak ya referansi silmek ya da tipe ait olmayan
    # bir korunmusluk iddia etmek olurdu.
    ref = con.execute("""
        SELECT r.candidate_id cid, e.tier tier, g.architecture a,
               g.upstream_family f, g.upstream_product prod, g.intergenic_bp bp,
               o.layout layout
        FROM ro r JOIN ro_evidence e USING(candidate_id)
        LEFT JOIN ro_regulation g USING(candidate_id)
        LEFT JOIN operon o USING(candidate_id)
        WHERE r.is_confirmed=1 AND r.ro_cluster=? AND e.tier='characterized'
        ORDER BY r.candidate_id LIMIT 1""", (cluster,)).fetchone()
    if ref and ref["a"]:
        x["ref_arch"] = {
            "cid": ref["cid"], "tier": ref["tier"], "architecture": ref["a"],
            "family": ref["f"] or "", "product": ref["prod"] or "",
            "intergenic_bp": ref["bp"], "layout": ref["layout"] or "[alpha]"}
        # Referans sekli TEK bir girisin sekli oldugu icin komsu genlerin
        # gercek urun adlari BILINIYOR; oklara baslik olarak onlar baglanir.
        # Oncesinde SCP2 gibi bir gen ekranda yalnizca "other" diye goruluyordu.
        ref_titles = [r[0] for r in con.execute(
            "SELECT COALESCE(NULLIF(product,''), category, component) "
            "FROM operon_gene WHERE candidate_id=? ORDER BY position",
            (ref["cid"],))]
        x["ref_arch_svg"] = atlas.operon_regulator_svg(
            ref["layout"] or "[alpha]",
            {"architecture": ref["a"], "upstream_family": ref["f"] or "",
             "intergenic_bp": ref["bp"], "upstream_product": ref["prod"] or ""},
            titles=ref_titles)
    if atlas.table_exists(con, "cluster_sdp"):
        row = con.execute("SELECT columns FROM cluster_sdp WHERE cluster=?", (cluster,)).fetchone()
        if row:
            import json as _json
            x["sdp_columns"] = _json.loads(row[0])
            x["sdp_leaves"] = [
                {"leaf_id": r[0], "size": r[1], "residues": _json.loads(r[2]), "n_differ": r[3]}
                for r in con.execute(
                    "SELECT leaf_id, size, residues, n_differ FROM leaf_sdp "
                    "WHERE cluster=? ORDER BY size DESC LIMIT 20", (cluster,))]
    tpath = apath("trees", f"{cluster}.nwk")
    if os.path.exists(tpath) and os.path.getsize(tpath) > 0:
        try:
            _nwk = open(tpath).read()
            x["tree_svg"], _ = atlas.render_tree_svg(
                _nwk, atlas.tip_info_from_db(con, u, ECOLOGY, short_names()),
                width=900, row_h=12, label_w=360)
            # Zayif dal isaretleri bu agacta da ciziliyor; sayilari sayfada
            # yazilmazsa okuyucu aciklanmamis daireler gorur.
            x["tree_support"] = atlas.newick_support_summary(_nwk)
        except Exception:
            x["tree_svg"] = None
    return x


#  Capraz referans URL kaliplari. Tek sahipli: cross_references.py ayni
#  kaliplari JSON'a yaziyor, dogrulayici dort adin hepsinin burada bir
#  karsiligi olmasini SINIYOR. Yeni bir veritabani eklenirse ikisi birlikte
#  degisir, yoksa sayfa URL'si olmayan bir ad basar.
XREF_DBS = [
    ("ec", "EC", "https://enzyme.expasy.org/EC/{}"),
    ("uniprot", "UniProt", "https://www.uniprot.org/uniprotkb/{}"),
    ("refseq_protein", "RefSeq protein",
     "https://www.ncbi.nlm.nih.gov/protein/{}"),
    ("embl_nucleotide", "ENA nucleotide",
     "https://www.ebi.ac.uk/ena/browser/view/{}"),
    #  2026-11-13'te eklenen uc tanimlayici. UCU DE BIRER IDDIA, salt
    #  adres degil: InterPro bir domain iddiasi, KEGG ve MetaCyc birer yol
    #  iddiasi, ve ucu de EC gibi otomatik anotasyonla yayiliyor. Disa
    #  baglanti olarak basilirlar, bu veritabaninin kendi ifadesi olarak
    #  degil.
    ("interpro", "InterPro",
     "https://www.ebi.ac.uk/interpro/entry/InterPro/{}/"),
    ("kegg", "KEGG gene", "https://www.genome.jp/entry/{}"),
    ("metacyc", "MetaCyc", "https://biocyc.org/getid?id={}"),
]


def cross_refs(con, sequence):
    """Bir dizinin dis tanimlayicilari: [(ad, etiket, [(kimlik, url)])].

    Birim DIZI (madde A1): ayni protein birkac accession altinda gecebilir,
    o yuzden anahtar `seq_md5`. Dizisi olmayan ya da `ro_xref` henuz
    kurulmamis bir kurulumda bos liste doner -- sayfa bolumu o zaman hic
    cizilmez, "yok" yazmaz, cunku olculmemis ile bulunmamis ayri seylerdir.
    """
    if not sequence or not atlas.table_exists(con, "ro_xref"):
        return []
    digest = hashlib.md5(sequence.strip().upper().encode()).hexdigest()
    rows = con.execute(
        "SELECT db, identifier FROM ro_xref WHERE seq_md5=? ORDER BY db, identifier",
        (digest,)).fetchall()
    if not rows:
        return []
    by_db = {}
    for db, ident in rows:
        by_db.setdefault(db, []).append(ident)
    out = []
    for name, label, template in XREF_DBS:
        idents = by_db.get(name)
        if idents:
            out.append((name, label,
                        [(i, template.format(i)) for i in idents]))
    return out


def ec_review_note(con, sequence):
    """EC'yi tasiyan UniProt kaydi incelenmis mi -- EC bir IDDIA.

    EC'lerin neredeyse tamami TrEMBL'den geliyor, yani otomatik anotasyonun
    yaydigi bir tahmin. Sayfa bunu yazmazsa okuyucu EC'yi kuratorlu substrat
    etiketiyle ayni kefeye koyar.
    """
    if not sequence or not atlas.table_exists(con, "ro_xref"):
        return None
    digest = hashlib.md5(sequence.strip().upper().encode()).hexdigest()
    row = con.execute(
        """SELECT SUM(CASE WHEN s.entry_type LIKE '%reviewed%'
                            AND s.entry_type NOT LIKE '%unreviewed%'
                           THEN 1 ELSE 0 END),
                  COUNT(*)
             FROM ro_union_source u
             JOIN ro_protein_status s ON s.accession = u.accession
            WHERE u.seq_md5 = ? AND u.source = 'uniprot'
              AND s.ec IS NOT NULL AND s.ec <> ''""", (digest,)).fetchone()
    if not row or not row[1]:
        return None
    return {"reviewed": row[0] or 0, "entries": row[1]}


def ro_extras(con, candidate_id):
    x = {}
    for table in ("ro_evidence", "ro_etc", "ro_regulation"):
        x[table] = con.execute(f"SELECT * FROM {table} WHERE candidate_id=?",
                               (candidate_id,)).fetchone() if atlas.table_exists(con, table) else None
    return x
