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
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, OrderedDict
from typing import Optional

from fastapi import FastAPI, Form, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
for p in (HERE, PARENT):
    if p not in sys.path:
        sys.path.insert(0, p)

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
    return info


BUILD = build_info()
templates.env.globals.update(u=u, slug=slug, leaf_slug=leaf_slug, fmt=fmt, pct=pct, build=BUILD)
templates.env.filters["fmt"] = fmt


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


def cluster_parts(cluster):
    parts = (cluster or "").split("_", 2)
    return {"group": parts[0] if parts else "", "id": parts[1] if len(parts) > 1 else "",
            "gene": parts[2] if len(parts) > 2 else cluster}


def genus(organism):
    if not organism:
        return ""
    w = organism.split()
    if w[0] in ("Candidatus", "uncultured") and len(w) > 1:
        return w[1]
    return w[0]


def render(request, name, **ctx):
    ctx.setdefault("request", request)
    ctx.setdefault("static_mode", STATIC_MODE)
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
            ("leaf", "leaves", "SELECT COUNT(*) FROM leaf"),
            ("operon", "operon_complete", "SELECT COUNT(*) FROM operon WHERE completeness=3"),
            ("operon", "operon_any", "SELECT COUNT(*) FROM operon WHERE completeness>0"),
            ("ro_subfamily", "novel", "SELECT COUNT(*) FROM ro_subfamily WHERE assignment_class='novel_candidate'"),
            ("ro_domain", "eukaryote", "SELECT COUNT(*) FROM ro_domain WHERE domain='Eukaryota'")]:
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
    genera = {}
    for cl, org in con.execute("SELECT ro_cluster, organism FROM ro r JOIN replicon p "
                               "USING(nucleotide_id) WHERE is_confirmed=1"):
        genera.setdefault(cl, set()).add(genus(org))
    out = []
    for r in rows:
        eco = ECOLOGY.get(r["cluster"], {})
        d = dict(r)
        d.update(cluster_parts(r["cluster"]))
        d.update({"substrate": eco.get("substrate", ""), "sclass": eco.get("substrate_class", ""),
                  "confidence": eco.get("confidence", ""), "leaves": leaves.get(r["cluster"], 0),
                  "euk": euk.get(r["cluster"], 0) or 0, "genera": len(genera.get(r["cluster"], ()))})
        out.append(d)
    return out


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
        return None
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
    gcount = Counter()
    for r in info["genera"]:
        gcount[genus(r["organism"])] += r["n"]
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
    info["layouts"] = con.execute("""
        SELECT layout, COUNT(*) n FROM operon o JOIN ro r USING(candidate_id)
        WHERE r.is_confirmed=1 AND r.ro_cluster=? GROUP BY layout ORDER BY n DESC LIMIT 10""",
        (cluster,)).fetchall()
    info["leaves"] = con.execute("""
        SELECT l.leaf_id, l.size, l.median_identity, l.is_homogeneous, l.top_genera,
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
        LEFT JOIN operon o ON o.candidate_id = r.candidate_id
        WHERE {' AND '.join(where)}"""
    total = con.execute("SELECT COUNT(*) " + sql, params).fetchone()[0]
    rows = con.execute(f"""
        SELECT r.candidate_id, r.protein_id, r.locus_tag, r.product, r.ro_cluster, r.hmm_score,
               r.model_coverage, p.organism, p.is_plasmid, p.nucleotide_id,
               rl.leaf_id, rs.assignment_class, rs.core_identity,
               o.completeness, o.layout
        {sql} ORDER BY p.organism, r.candidate_id LIMIT ? OFFSET ?""",
        params + [size, (page - 1) * size]).fetchall()
    return total, rows


@app.get("/cluster/{cluster}", response_class=HTMLResponse)
def cluster_page(request: Request, cluster: str, page: int = 1, q: Optional[str] = None):
    con = connect()
    try:
        info = cluster_detail(con, cluster)
        if not info:
            raise HTTPException(404, "cluster not found")
        total, rows = members_query(con, cluster=cluster, q=q, page=page)
        return render(request, "cluster.html", c=info, members=rows, total=total, page=page,
                      pages=(total + PAGE_SIZE - 1) // PAGE_SIZE, q=q or "",
                      chem=CHEMISTRY.get(cluster), x=cluster_extras(con, cluster))
    finally:
        con.close()


@app.get("/leaf/{leaf_id}", response_class=HTMLResponse)
def leaf_page(request: Request, leaf_id: str, page: int = 1):
    leaf_id = leaf_unslug(leaf_id)
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
               d.domain, d.euk_group
        FROM ro r JOIN replicon p USING(nucleotide_id)
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
    # motif residues on the stored sequence are not persisted; report flags only
    return {"ro": ro, "neighbors": neighbors, "operon": operon, "operon_genes": operon_genes,
            "svg": neighborhood_svg(ro, neighbors, operon_ids), "operon_ids": operon_ids,
            "cluster": cluster_parts(ro["ro_cluster"]),
            "eco": ECOLOGY.get(ro["ro_cluster"], {})}


@app.get("/ro/{ro_slug}", response_class=HTMLResponse)
def ro_page(request: Request, ro_slug: str):
    con = connect()
    try:
        d = ro_detail(con, unslug(ro_slug))
        if not d:
            raise HTTPException(404, "entry not found")
        return render(request, "ro.html", x=ro_extras(con, d["ro"]["candidate_id"]), **d)
    finally:
        con.close()


# ------------------------------------------------------------------ search
FACETS = [("tier", "Evidence level"), ("domain", "Domain of life"),
          ("family", "Chemical family"), ("cluster", "Enzyme type")]
FLAGS = [("plasmid", "is_plasmid", "on a plasmid"),
         ("partner", "has_partner", "operon partner verified"),
         ("regulator", "has_regulator", "divergent regulator upstream")]


def fts_expression(text):
    """Turn user input into a safe FTS5 expression.

    Raw input is tried first so that power users keep AND / OR / NEAR and prefix
    search. If SQLite rejects it, every token is quoted instead, which can never
    be a syntax error.
    """
    tokens = [t for t in re.split(r"[^\w*]+", text) if t]
    fallback = " ".join('"%s"' % t.replace('"', "") for t in tokens)
    return text.strip(), fallback


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
    try:
        return run(raw)
    except sqlite3.OperationalError:
        if not fallback:
            return 0, [], {}
        try:
            return run(fallback)
        except sqlite3.OperationalError:
            return 0, [], {}


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
                          filters=filters, flags=flags, active=active)
        total, rows, facets = search_query(con, q, filters, flags, page)
        clusters = [c for c in cluster_table(con)
                    if q and (q.lower() in c["cluster"].lower()
                              or q.lower() in (c["substrate"] or "").lower())] if q else []
        return render(request, "search.html", q=q, members=rows, total=total, page=page,
                      pages=(total + PAGE_SIZE - 1) // PAGE_SIZE, clusters=clusters,
                      facets=facets, filters=filters, flags=flags, active=active)
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
        return render(request, "about.html", totals=totals(con))
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


def csv_response(rows, filename):
    buf = io.StringIO()
    if rows:
        w = csv.DictWriter(buf, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(dict(r))
    return Response(buf.getvalue(), media_type="text/csv", headers={
        "Content-Disposition": f'attachment; filename="{filename}"'})


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
        rows = con.execute(f"SELECT {MEMBER_COLS} {MEMBER_FROM} WHERE r.is_confirmed=1 AND r.ro_cluster=?",
                           (cluster,)).fetchall()
        return csv_response(rows, f"{cluster}.csv")
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

ANALYSIS_DIR = os.environ.get("ROAR_ANALYSIS", os.path.join(PARENT, "analysis_out"))
templates.env.globals.update(
    layout_svg=atlas.layout_svg, operon_regulator_svg=atlas.operon_regulator_svg,
    reaction_scheme_svg=atlas.reaction_scheme_svg,
    gap_histogram=atlas.gap_histogram,
    amedian=atlas.median, GROUP_COLORS=atlas.GROUP_COLORS, TIER_COLORS=atlas.TIER_COLORS, TIERS=atlas.TIERS,
    TIER_LABEL=atlas.TIER_LABEL, TIER_MEANING=atlas.TIER_MEANING,
    COMPONENT_COLORS=COMPONENT_COLORS)


def apath(*parts):
    return os.path.join(ANALYSIS_DIR, *parts)


def _overview(con):
    return atlas.per_cluster_overview(con, ECOLOGY)


@app.get("/atlas", response_class=HTMLResponse)
def atlas_home(request: Request):
    con = connect()
    try:
        return render(request, "atlas_home.html", over=_overview(con), totals=totals(con),
                      stats=atlas.read_json(apath("stats.json")))
    finally:
        con.close()


@app.get("/atlas/phylogeny", response_class=HTMLResponse)
def atlas_phylogeny(request: Request):
    if not os.path.exists(apath("tree_all.nwk")):
        raise HTTPException(404, "tree not built (run build_phylogeny.py)")
    con = connect()
    try:
        svg, n = atlas.render_tree_svg(open(apath("tree_all.nwk")).read(),
                                       atlas.tip_info_from_db(con, u, ECOLOGY),
                                       width=1040, row_h=11, label_w=380)
        return render(request, "atlas_phylogeny.html", svg=svg, n_tips=n)
    finally:
        con.close()


@app.get("/atlas/network", response_class=HTMLResponse)
def atlas_network(request: Request, min_identity: float = 35.0):
    return render(request, "atlas_network.html",
                  data=atlas.ssn(apath("ssn_nodes.csv"), apath("ssn_edges.csv"), min_identity),
                  matrix=atlas.identity_matrix(apath("cluster_identity_matrix.csv")),
                  refpairs=atlas.read_csv(apath("reference_pairs.csv"))[:40],
                  min_identity=min_identity)


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
                      tree=atlas.taxonomy_tree(con))
    finally:
        con.close()


@app.get("/atlas/ecology", response_class=HTMLResponse)
def atlas_ecology(request: Request):
    con = connect()
    try:
        return render(request, "atlas_ecology.html", over=_overview(con),
                      stats=atlas.read_json(apath("stats.json")))
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
        return render(request, "atlas_regulation.html", arch=arch, fams=fams, gaps=gaps,
                      products=products, per_cluster=atlas.read_csv(apath("regulation_by_cluster.csv")),
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
        return render(request, "atlas_novel.html", rows=rows, variants=variants,
                      matrix=matrix, tiers=tiers, classes=classes, domains=domains,
                      n_total=len(rows))
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


@app.get("/atlas/quality", response_class=HTMLResponse)
def atlas_quality(request: Request):
    con = connect()
    try:
        red = atlas.read_json(apath("redundancy.json"))
        return render(request, "atlas_quality.html", red=red,
                      stats=atlas.read_json(apath("stats.json")),
                      motif=atlas.read_json(apath("motif_stats.json")),
                      validation=atlas.read_json(apath("operon_validation.json")),
                      downloads=atlas.download_table(ANALYSIS_DIR), totals=totals(con))
    finally:
        con.close()


@app.get("/atlas/statistics", response_class=HTMLResponse)
def atlas_statistics(request: Request):
    data = atlas.read_json(apath("stats.json"))
    if not data:
        raise HTTPException(404, "stats.json not found (run stats_overview.py)")
    return render(request, "atlas_statistics.html", stats=data)


@app.get("/download/tree_all.nwk")
def dl_tree():
    if not os.path.exists(apath("tree_all.nwk")):
        raise HTTPException(404, "tree not built")
    return PlainTextResponse(open(apath("tree_all.nwk")).read(), headers={
        "Content-Disposition": 'attachment; filename="roar_tree.nwk"'})


@app.get("/download/analysis/{name}")
def dl_analysis(name: str):
    allowed = {"ssn_edges.csv", "ssn_nodes.csv", "cluster_identity_matrix.csv",
               "reference_pairs.csv", "regulation_by_cluster.csv", "evidence_by_cluster.csv",
               "etc_by_cluster.csv", "leaf_profiles.csv", "cluster_ecology_stats.csv",
               "variant_signatures.csv", "motif_stats.json", "operon_validation.json",
               "redundancy.json", "cooccurrence.json", "substrate_predictability.json",
               "null_model.csv", "sdp_positions.csv", "stats.json"}
    if name not in allowed or not os.path.exists(apath(name)):
        raise HTTPException(404, "not available")
    media = "application/json" if name.endswith(".json") else "text/csv"
    return Response(open(apath(name)).read(), media_type=media,
                    headers={"Content-Disposition": f'attachment; filename="{name}"'})


# --- per-page extras -------------------------------------------------------
def cluster_extras(con, cluster):
    x = {"tiers": {}, "etc": [], "tree_svg": None, "nearest": [], "phyla": [],
         "transposon": 0, "regulation": None, "reg_families": []}
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
    x["nearest"] = atlas.nearest_types(
        atlas.identity_matrix(apath("cluster_identity_matrix.csv")), cluster)
    reg_rows = {r["cluster"]: r for r in atlas.read_csv(apath("regulation_by_cluster.csv"))}
    x["reg_row"] = reg_rows.get(cluster)
    top_layout = con.execute(
        "SELECT o.layout FROM operon o JOIN ro r USING(candidate_id) "
        "WHERE r.is_confirmed=1 AND r.ro_cluster=? GROUP BY o.layout "
        "ORDER BY COUNT(*) DESC LIMIT 1", (cluster,)).fetchone()
    if top_layout and x["reg_row"]:
        row = x["reg_row"]
        dominant = max(("divergent_regulator", "codirectional_regulator", "divergent_other",
                        "codirectional_other", "unknown"),
                       key=lambda k: int(row.get(k) or 0))
        x["canonical_arch"] = dominant
        x["canonical_svg"] = atlas.operon_regulator_svg(top_layout[0], {
            "architecture": dominant,
            "upstream_family": row.get("top_regulator_family") or "",
            "intergenic_bp": int(row["median_intergenic_bp"]) if row.get("median_intergenic_bp") else None,
            "upstream_product": "typical upstream regulator of this type"})
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
            x["tree_svg"], _ = atlas.render_tree_svg(
                open(tpath).read(), atlas.tip_info_from_db(con, u, ECOLOGY),
                width=900, row_h=12, label_w=360)
        except Exception:
            x["tree_svg"] = None
    return x


def ro_extras(con, candidate_id):
    x = {}
    for table in ("ro_evidence", "ro_etc", "ro_regulation"):
        x[table] = con.execute(f"SELECT * FROM {table} WHERE candidate_id=?",
                               (candidate_id,)).fetchone() if atlas.table_exists(con, table) else None
    return x
