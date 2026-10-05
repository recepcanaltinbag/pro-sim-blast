"""
Static export of ROAR-DB for free static hosting (GitHub Pages, Netlify, Cloudflare Pages).

Renders every page through the FastAPI app with STATIC_MODE on (links get .html),
writes a client-side search index, and copies downloads. The sequence classifier
needs a server and is omitted from the static site (the Docker/HF Space has it).

    python3 freeze.py --out site --base /pro-sim-blast      # base = repo name for GitHub project pages
    python3 freeze.py --out site --base ""             # user/organisation pages or custom domain

Output size: roughly 11k entry pages * ~25 kB + downloads; fits GitHub Pages limits.
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
    args = ap.parse_args()

    os.environ["ROAR_BASE"] = args.base
    import app as A
    A.STATIC_MODE = True
    A.PAGE_SIZE = 10 ** 6      # static pages carry the full member list (no ?page= files)
    A.BASE = args.base.rstrip("/")
    from fastapi.testclient import TestClient
    client = TestClient(A.app)

    out = args.out
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)
    shutil.copytree(os.path.join(HERE, "static"), os.path.join(out, "static"))

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
              "atlas/regulation", "atlas/ecology", "atlas/evidence", "atlas/statistics"):
        save("/" + p, p + ".html")
    save("/download/tree_all.nwk", "download/tree_all.nwk")
    for name in ("ssn_edges.csv", "ssn_nodes.csv", "cluster_identity_matrix.csv",
                 "reference_pairs.csv", "regulation_by_cluster.csv", "evidence_by_cluster.csv",
                 "etc_by_cluster.csv", "leaf_profiles.csv", "cluster_ecology_stats.csv",
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
        save(f"/download/leaf/{A.leaf_slug(l)}.fasta", f"download/leaf/{A.leaf_slug(l)}.fasta")
    print(f"[pages] {len(entries)} entries")
    for i, e in enumerate(entries, 1):
        s = A.slug(e)
        save(f"/ro/{s}", f"ro/{s}.html")
        save(f"/download/ro/{s}.fasta", f"download/ro/{s}.fasta")
        if i % 1000 == 0:
            print(f"   {i}/{len(entries)}")
    print("[downloads] bulk")
    save("/download/all_confirmed.fasta", "download/all_confirmed.fasta")
    save("/download/all_confirmed.csv", "download/all_confirmed.csv")
    save("/download/operons.csv", "download/operons.csv")

    # client-side search index
    print("[search] index")
    rows = con.execute("""
        SELECT r.candidate_id, r.protein_id, r.locus_tag, r.product, r.ro_cluster, p.organism, p.is_plasmid
        FROM ro r JOIN replicon p USING(nucleotide_id) WHERE r.is_confirmed=1""").fetchall()
    with open(os.path.join(out, "search_index.json"), "w") as fh:
        json.dump([[A.slug(r[0]), r[1], r[2], r[3], r[4], r[5], r[6]] for r in rows], fh)
    with open(os.path.join(out, "search.html"), "w") as fh:
        fh.write(SEARCH_PAGE.replace("__BASE__", A.BASE))
    con.close()
    print(f"[done] {out}")


SEARCH_PAGE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>Search · ROAR-DB</title>
<link rel="stylesheet" href="__BASE__/static/style.css"></head><body>
<header class="top"><a class="brand" href="__BASE__/index.html"><span class="brand__mark">RO</span> ROAR-DB</a>
<nav><a href="__BASE__/clusters.html">Clusters</a><a href="__BASE__/search.html">Search</a><a href="__BASE__/about.html">Methods</a></nav></header>
<main><h1>Search</h1>
<form class="bigsearch" onsubmit="return false"><input id="q" type="search" placeholder="organism, product, protein ID, locus tag, entry ID or reference type" autofocus><button onclick="run()">Search</button></form>
<p class="muted" id="status">Loading index…</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Entry</th><th>Organism</th><th>Product</th><th>Type</th></tr></thead><tbody id="rows"></tbody></table></div>
</main>
<script>
let IDX=[];
fetch('__BASE__/search_index.json').then(r=>r.json()).then(d=>{IDX=d;document.getElementById('status').textContent=d.length.toLocaleString()+' entries indexed';
  const q=new URLSearchParams(location.search).get('q'); if(q){document.getElementById('q').value=q;run();}});
const esc=s=>String(s==null?'':s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function run(){const q=document.getElementById('q').value.trim().toLowerCase(); if(!q)return;
  const hits=IDX.filter(r=>r.some(v=>typeof v==='string'&&v.toLowerCase().includes(q))).slice(0,500);
  document.getElementById('status').textContent=hits.length+(hits.length==500?'+':'')+' matches';
  document.getElementById('rows').innerHTML=hits.map(r=>`<tr><td><a href="__BASE__/ro/${esc(r[0])}.html">${esc(r[1]||r[2]||r[0])}</a>${r[6]?' <span class="chip chip--plasmid">plasmid</span>':''}</td><td>${esc(r[5])}</td><td class="small">${esc(r[3])}</td><td><a href="__BASE__/cluster/${esc(r[4])}.html">${esc(r[4])}</a></td></tr>`).join('');}
document.getElementById('q').addEventListener('keydown',e=>{if(e.key==='Enter')run();});
</script></body></html>"""


if __name__ == "__main__":
    main()
