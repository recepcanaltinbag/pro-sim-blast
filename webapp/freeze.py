"""
Static export of ROAR-DB for free static hosting (GitHub Pages, Netlify, Cloudflare Pages).

Renders every page through the FastAPI app with STATIC_MODE on (links get .html),
writes a client-side search index, and copies downloads. The sequence classifier
needs a server and is omitted from the static site (the Docker/HF Space has it).

    python3 freeze.py --out site --base /pro-sim-blast      # base = repo name for GitHub project pages
    python3 freeze.py --out site --base ""             # user/organisation pages or custom domain

Per-entry and per-variant download files are deliberately NOT exported: they would
add about 13,000 files and 50 MB while the sequence is already shown on the page
and the bulk FASTA covers every entry. GitHub Pages builds a site of this size
slowly, so file count is worth spending carefully.
"""

import argparse
import json
import os
import shutil
import sqlite3
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(HERE, "site"))
    ap.add_argument("--base", default="", help="URL prefix, e.g. /pro-sim-blast for a GitHub project page")
    ap.add_argument("--limit", type=int, default=0, help="only N entry pages (for testing)")
    ap.add_argument("--pages-only", action="store_true",
                    help="refresh the overview, type and variant pages plus static files and "
                         "bulk downloads; leave the per-entry pages in place (they are the slow "
                         "part and change only when the database itself changes)")
    args = ap.parse_args()

    os.environ["ROAR_BASE"] = args.base
    import app as A
    A.STATIC_MODE = True
    A.PAGE_SIZE = 10 ** 6      # static pages carry the full member list (no ?page= files)
    A.BASE = args.base.rstrip("/")
    from fastapi.testclient import TestClient
    client = TestClient(A.app)

    out = args.out
    if os.path.exists(out) and not args.pages_only:
        shutil.rmtree(out)
    os.makedirs(out, exist_ok=True)
    static_dst = os.path.join(out, "static")
    if os.path.exists(static_dst):
        shutil.rmtree(static_dst)
    shutil.copytree(os.path.join(HERE, "static"), static_dst)

    def save(url, path, binary=False):
        r = client.get(url)
        if r.status_code != 200:
            print(f"  !! {r.status_code} {url}")
            return
        full = os.path.join(out, path.lstrip("/"))
        os.makedirs(os.path.dirname(full), exist_ok=True)
        with open(full, "wb") as fh:
            fh.write(r.content)

    con = sqlite3.connect(f"file:{A.DB_PATH}?mode=ro", uri=True)
    clusters = [r[0] for r in con.execute(
        "SELECT DISTINCT ro_cluster FROM ro WHERE is_confirmed=1")]
    leaves = [r[0] for r in con.execute("SELECT leaf_id FROM leaf")]
    entries = [r[0] for r in con.execute(
        "SELECT candidate_id FROM ro WHERE is_confirmed=1 ORDER BY candidate_id")]
    if args.limit:
        entries = entries[:args.limit]

    print("[pages] fixed")
    save("/", "index.html"); save("/clusters", "clusters.html"); save("/about", "about.html")
    print("[pages] atlas")
    for p in ("atlas", "atlas/phylogeny", "atlas/network", "atlas/taxonomy", "atlas/operons",
              "atlas/regulation", "atlas/ecology", "atlas/evidence", "atlas/novel",
              "atlas/cooccurrence", "atlas/statistics", "atlas/quality"):
        save("/" + p, p + ".html")
    save("/download/tree_all.nwk", "download/tree_all.nwk")
    save("/download/novel_candidates.fasta", "download/novel_candidates.fasta")
    for name in ("ssn_edges.csv", "ssn_nodes.csv", "cluster_identity_matrix.csv",
                 "reference_pairs.csv", "regulation_by_cluster.csv", "evidence_by_cluster.csv",
                 "etc_by_cluster.csv", "leaf_profiles.csv", "cluster_ecology_stats.csv",
                 "variant_signatures.csv", "motif_stats.json", "operon_validation.json",
                 "redundancy.json", "cooccurrence.json", "substrate_predictability.json",
                 "habitat.json",
                 "null_model.csv", "sdp_positions.csv", "stats.json"):
        save("/download/analysis/" + name, "download/analysis/" + name)
    print(f"[pages] {len(clusters)} clusters")
    for c in clusters:
        save(f"/cluster/{c}", f"cluster/{c}.html")
        save(f"/download/cluster/{c}.fasta", f"download/cluster/{c}.fasta")
        save(f"/download/cluster/{c}.csv", f"download/cluster/{c}.csv")
    print(f"[pages] {len(leaves)} leaves")
    for l in leaves:
        save(f"/leaf/{A.leaf_slug(l)}", f"leaf/{A.leaf_slug(l)}.html")
    if args.pages_only:
        print("[search] index")
        write_search_index(con, out, A)
        con.close()
        print(f"[done] {out} (entry pages left untouched)")
        return

    print(f"[pages] {len(entries)} entries")
    for i, e in enumerate(entries, 1):
        s = A.slug(e)
        save(f"/ro/{s}", f"ro/{s}.html")
        if i % 1000 == 0:
            print(f"   {i}/{len(entries)}")
    print("[downloads] bulk")
    save("/download/all_confirmed.fasta", "download/all_confirmed.fasta")
    save("/download/all_confirmed.csv", "download/all_confirmed.csv")
    save("/download/operons.csv", "download/operons.csv")

    # client-side search index
    print("[search] index")
    write_search_index(con, out, A)
    con.close()
    print(f"[done] {out}")


def write_search_index(con, out, A):
    """Client-side index plus the static search page.

    The static site has no server, so the facets of the hosted application are
    reproduced in the browser. Every field the facets filter on is written into
    the index, which keeps the two versions of the search consistent.
    """
    has_search = con.execute(
        "SELECT name FROM sqlite_master WHERE name='ro_search'").fetchone()
    if has_search:
        rows = con.execute("""
            SELECT candidate_id, protein_id, organism, product, cluster, leaf_id,
                   tier, domain, family, is_plasmid, has_partner, has_regulator
            FROM ro_search""").fetchall()
        data = [[A.slug(r[0]), r[1] or "", r[2] or "", r[3] or "", r[4] or "", r[5] or "",
                 r[6] or "", r[7] or "", r[8] or "", r[9] or 0, r[10] or 0, r[11] or 0]
                for r in rows]
    else:
        rows = con.execute("""
            SELECT r.candidate_id, r.protein_id, p.organism, r.product, r.ro_cluster,
                   '', '', '', '', p.is_plasmid, 0, 0
            FROM ro r JOIN replicon p USING(nucleotide_id) WHERE r.is_confirmed=1""").fetchall()
        data = [[A.slug(r[0])] + list(r[1:]) for r in rows]
    with open(os.path.join(out, "search_index.json"), "w") as fh:
        json.dump(data, fh, separators=(",", ":"))
    with open(os.path.join(out, "search.html"), "w") as fh:
        fh.write(SEARCH_PAGE.replace("__BASE__", A.BASE))
    print(f"   {len(data)} entries indexed")


SEARCH_PAGE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>Search · ROAR-DB</title>
<link rel="stylesheet" href="__BASE__/static/style.css"></head><body>
<header class="top"><a class="brand" href="__BASE__/index.html"><span class="brand__mark">RO</span> ROAR-DB</a>
<nav><a href="__BASE__/atlas.html">Atlas</a><a href="__BASE__/clusters.html">Types</a>
<a href="__BASE__/search.html">Search</a><a href="__BASE__/about.html">Methods</a></nav></header>
<main>
<h1>Search</h1>
<form class="bigsearch" onsubmit="return false">
  <input id="q" type="search" placeholder="organism, product, protein identifier, locus tag, enzyme type or substrate" autofocus>
  <button onclick="run()">Search</button>
</form>
<p class="note">Several words are combined with AND. A trailing asterisk matches a prefix.
Filters on the left narrow the result without a new search.
<span id="status">loading index…</span></p>
<div class="searchlayout">
  <aside class="facets" id="facets"></aside>
  <div class="results">
    <h2>Entries <span class="note" id="count"></span></h2>
    <div class="tablewrap"><table class="data">
      <thead><tr><th>Entry</th><th>Organism</th><th>Product</th><th>Type</th>
      <th>Evidence</th><th>Variant</th></tr></thead>
      <tbody id="rows"></tbody></table></div>
    <nav class="pager" id="pager"></nav>
  </div>
</div>
</main>
<script>
const BASE = "__BASE__";
const TIER_COLORS = {characterized:'#1f7a4d', close_homolog:'#7fb069', family_member:'#f2c14e',
                     distant:'#f78154', novel:'#8e44ad'};
const F = {q:'', tier:'', domain:'', family:'', cluster:'', plasmid:0, partner:0, regulator:0, page:1};
const PAGE = 100;
let IDX = [], VIEW = [];
const esc = s => String(s==null?'':s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
// index columns: 0 slug 1 protein 2 organism 3 product 4 cluster 5 leaf 6 tier 7 domain 8 family 9 plasmid 10 partner 11 regulator
const TEXT = r => (r[0]+' '+r[1]+' '+r[2]+' '+r[3]+' '+r[4]+' '+r[5]+' '+r[8]).toLowerCase();

fetch(BASE + '/search_index.json').then(r => r.json()).then(d => {
  IDX = d.map(r => { r.push(TEXT(r)); return r; });
  document.getElementById('status').textContent = d.length.toLocaleString() + ' entries indexed.';
  const p = new URLSearchParams(location.search);
  ['q','tier','domain','family','cluster'].forEach(k => { if (p.get(k)) F[k] = p.get(k); });
  ['plasmid','partner','regulator'].forEach(k => { if (p.get(k)) F[k] = 1; });
  document.getElementById('q').value = F.q;
  run();
});

function matches(r) {
  if (F.tier && r[6] !== F.tier) return false;
  if (F.domain && r[7] !== F.domain) return false;
  if (F.family && r[8] !== F.family) return false;
  if (F.cluster && r[4] !== F.cluster) return false;
  if (F.plasmid && !r[9]) return false;
  if (F.partner && !r[10]) return false;
  if (F.regulator && !r[11]) return false;
  if (!F.q) return true;
  const hay = r[12];
  return F.q.toLowerCase().split(/\s+/).filter(Boolean).every(tok =>
    tok.endsWith('*') ? hay.includes(tok.slice(0, -1)) : hay.includes(tok));
}

function run() {
  F.q = document.getElementById('q').value.trim();
  F.page = 1;
  render();
}

function setFilter(k, v) { F[k] = (F[k] === v) ? '' : v; F.page = 1; render(); }
function toggle(k) { F[k] = F[k] ? 0 : 1; F.page = 1; render(); }
function goPage(n) { F.page = n; render(); window.scrollTo(0, 0); }

function render() {
  VIEW = IDX.filter(matches);
  document.getElementById('count').textContent = VIEW.length.toLocaleString();
  const start = (F.page - 1) * PAGE;
  document.getElementById('rows').innerHTML = VIEW.slice(start, start + PAGE).map(r => {
    const tier = r[6] ? `<span class="chip" style="background:${TIER_COLORS[r[6]]||'#888'};color:#fff;border:0">${esc(r[6].replace('_',' '))}</span>` : '';
    const leaf = r[5] ? `<a href="${BASE}/leaf/${esc(r[5].replace('#','-'))}.html">${esc(r[5].split('#').pop())}</a>` : '–';
    const dom = (r[7] && r[7] !== 'Bacteria') ? ` <span class="chip">${esc(r[7])}</span>` : '';
    const pl = r[9] ? ' <span class="chip chip--plasmid">plasmid</span>' : '';
    return `<tr><td><a href="${BASE}/ro/${esc(r[0])}.html">${esc(r[1] || r[0])}</a>${pl}</td>` +
           `<td><i>${esc(r[2])}</i>${dom}</td><td class="small">${esc(r[3])}</td>` +
           `<td><a href="${BASE}/cluster/${esc(r[4])}.html">${esc(r[4].split('_').slice(2).join('_'))}</a></td>` +
           `<td>${tier}</td><td>${leaf}</td></tr>`;
  }).join('') || '<tr><td colspan="6" class="empty">No entry matches. Try fewer words or remove a filter.</td></tr>';

  const pages = Math.ceil(VIEW.length / PAGE);
  document.getElementById('pager').innerHTML = pages > 1
    ? (F.page > 1 ? `<a href="#" onclick="goPage(${F.page-1});return false">previous</a>` : '') +
      `<span>page ${F.page} of ${pages}</span>` +
      (F.page < pages ? `<a href="#" onclick="goPage(${F.page+1});return false">next</a>` : '')
    : '';
  renderFacets();
}

function renderFacets() {
  const defs = [['tier', 'Evidence level', 6], ['domain', 'Domain of life', 7],
                ['family', 'Chemical family', 8], ['cluster', 'Enzyme type', 4]];
  let html = '';
  for (const [key, label, col] of defs) {
    const counts = new Map();
    for (const r of IDX) {
      const saved = F[key]; F[key] = '';
      if (matches(r) && r[col]) counts.set(r[col], (counts.get(r[col]) || 0) + 1);
      F[key] = saved;
    }
    const items = [...counts.entries()].sort((a, b) => b[1] - a[1]).slice(0, 12);
    if (!items.length) continue;
    html += `<div class="facet"><h3>${label}</h3><ul>` + items.map(([v, n]) =>
      `<li><a href="#" onclick="setFilter('${key}', ${JSON.stringify(v)});return false">` +
      (F[key] === v ? `<b>${esc(v)}</b>` : esc(v)) + `</a><span class="note">${n.toLocaleString()}</span></li>`).join('') + '</ul></div>';
  }
  html += '<div class="facet"><h3>Genomic context</h3><ul>' +
    [['plasmid', 'on a plasmid'], ['partner', 'verified operon partner'], ['regulator', 'divergent regulator']]
      .map(([k, l]) => `<li><a href="#" onclick="toggle('${k}');return false">${F[k] ? '<b>' + l + '</b>' : l}</a></li>`).join('') +
    '</ul></div>';
  document.getElementById('facets').innerHTML = html;
}
document.getElementById('q').addEventListener('keydown', e => { if (e.key === 'Enter') run(); });
</script></body></html>"""


if __name__ == "__main__":
    main()
