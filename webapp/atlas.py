"""
Atlas data layer for the overview pages (phylogeny, network, taxonomy, ecology,
operons/ETC, evidence, statistics) plus a small Newick -> SVG tree renderer.
Reads the SQLite DB and analysis_out/ files; no extra dependencies.
"""

import csv
import html
import json
import os
from collections import Counter, defaultdict

GROUP_COLORS = {"1": "#c0392b", "2": "#8e44ad", "3": "#2980b9", "4": "#27ae60",
                "5": "#f39c12", "?": "#95a5a6"}
TIERS = ["characterized", "close_homolog", "family_member", "distant", "novel"]
TIER_COLORS = {"characterized": "#1f7a4d", "close_homolog": "#7fb069",
               "family_member": "#f2c14e", "distant": "#f78154", "novel": "#8e44ad"}
TIER_LABEL = {
    "characterized": "characterized reference (≥95 % identity)",
    "close_homolog": "close homolog (≥60 %)",
    "family_member": "family member (40–60 %)",
    "distant": "distant (25–40 %): type by homology only",
    "novel": "novel (<25 %): no close reference",
}
TIER_MEANING = {
    "characterized": "Essentially the experimentally characterized reference enzyme "
                     "(or a strain variant). The substrate label applies.",
    "close_homolog": "Close homolog of a characterized enzyme. Same reaction is likely, "
                     "but substrate range may differ.",
    "family_member": "Same enzyme family as the reference type. Substrate cannot be "
                     "transferred with confidence; the label is a hypothesis.",
    "distant": "Rieske oxygenase alpha subunit assigned to its nearest type by profile "
               "score only. Function unknown; treat the substrate label as the nearest "
               "characterized relative, not as an annotation.",
    "novel": "Confirmed RO alpha subunit with no close characterized relative. "
             "Candidate for a new enzyme type.",
}


def group_of(cluster):
    return (cluster or "?").split("_")[0]


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh))


def read_json(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def phylum_of(taxonomy):
    parts = (taxonomy or "").split("; ")
    if not parts or not parts[0]:
        return "?"
    if parts[0] == "Eukaryota":
        return parts[1] if len(parts) > 1 else "Eukaryota"
    return parts[2] if len(parts) > 2 else parts[-1]


# ----------------------------------------------------------------- aggregates
def per_cluster_overview(con, ecology):
    """One row per cluster with taxonomy, ecology, ETC, evidence and operon summaries."""
    rows = {}
    for r in con.execute("""
        SELECT r.ro_cluster cluster, COUNT(*) n, SUM(p.is_plasmid) plasmid,
               SUM(o.has_beta) beta, SUM(o.has_ferredoxin) fd, SUM(o.has_reductase) red,
               SUM(o.completeness=3) complete
        FROM ro r JOIN replicon p USING(nucleotide_id)
        LEFT JOIN operon o ON o.candidate_id=r.candidate_id
        WHERE r.is_confirmed=1 GROUP BY r.ro_cluster"""):
        d = dict(r)
        d["group"] = group_of(d["cluster"])
        d["gene"] = d["cluster"].split("_", 2)[-1]
        e = ecology.get(d["cluster"], {})
        d["substrate"] = e.get("substrate", "")
        d["sclass"] = e.get("substrate_class", "unknown")
        d["confidence"] = e.get("confidence", "")
        d["phyla"] = Counter()
        d["genera"] = set()
        d["tiers"] = Counter()
        d["etc"] = Counter()
        d["transposon"] = 0
        rows[d["cluster"]] = d
    for cl, tax, org in con.execute(
            "SELECT r.ro_cluster, p.taxonomy, p.organism FROM ro r JOIN replicon p USING(nucleotide_id) "
            "WHERE r.is_confirmed=1"):
        rows[cl]["phyla"][phylum_of(tax)] += 1
        rows[cl]["genera"].add((org or "?").split()[0])
    if table_exists(con, "ro_evidence"):
        for cl, tier, n in con.execute(
                "SELECT r.ro_cluster, e.tier, COUNT(*) FROM ro_evidence e JOIN ro r USING(candidate_id) "
                "GROUP BY 1, 2"):
            rows[cl]["tiers"][tier] += n
    if table_exists(con, "ro_etc"):
        for cl, prof, n in con.execute(
                "SELECT r.ro_cluster, t.etc_profile, COUNT(*) FROM ro_etc t JOIN ro r USING(candidate_id) "
                "GROUP BY 1, 2"):
            rows[cl]["etc"][prof] += n
    for cl, n in con.execute("""
        SELECT r.ro_cluster, COUNT(DISTINCT r.candidate_id) FROM ro r
        JOIN neighbor nb ON nb.candidate_id=r.candidate_id
        JOIN gene_category c ON c.neighbor_id=nb.neighbor_id AND c.method='regex_v1' AND c.category='transposon'
        WHERE r.is_confirmed=1 GROUP BY 1"""):
        rows[cl]["transposon"] = n
    for d in rows.values():
        d["genera_n"] = len(d["genera"])
        d["genera"] = sorted(d["genera"])[:50]
        d["phyla"] = dict(d["phyla"].most_common())
        d["tiers"] = {t: d["tiers"].get(t, 0) for t in TIERS}
        d["etc"] = dict(d["etc"].most_common(6))
    return sorted(rows.values(), key=lambda d: -d["n"])


def table_exists(con, name):
    return con.execute("SELECT name FROM sqlite_master WHERE name=?", (name,)).fetchone() is not None


def identity_matrix(path):
    rows = read_csv(path)
    if not rows:
        return None
    clusters = [r["cluster"] for r in rows]
    mat = [[float(r[c]) for c in clusters] for r in rows]
    return {"clusters": clusters, "matrix": mat}


def nearest_types(matrix, cluster, k=5):
    if not matrix or cluster not in matrix["clusters"]:
        return []
    i = matrix["clusters"].index(cluster)
    pairs = [(matrix["clusters"][j], matrix["matrix"][i][j])
             for j in range(len(matrix["clusters"])) if j != i]
    return sorted(pairs, key=lambda p: -p[1])[:k]


def ssn(nodes_path, edges_path, min_identity=30.0):
    nodes = read_csv(nodes_path)
    edges = [e for e in read_csv(edges_path) if float(e["identity"]) >= min_identity]
    return {"nodes": [{"id": n["id"], "kind": n["kind"], "cluster": n["cluster"],
                       "group": group_of(n["cluster"]), "leaf": n["leaf_id"],
                       "size": int(n["leaf_size"] or 1), "organism": n["organism"],
                       "tier": n["tier"]} for n in nodes],
            "edges": [{"source": e["source"], "target": e["target"],
                       "identity": float(e["identity"])} for e in edges]}


# ----------------------------------------------------------------- newick
class Node:
    __slots__ = ("name", "length", "children", "x", "y")

    def __init__(self):
        self.name, self.length, self.children, self.x, self.y = "", 0.0, [], 0.0, 0.0


def parse_newick(text):
    """Minimal Newick parser supporting quoted labels and branch lengths."""
    text = text.strip().rstrip(";")
    pos = 0

    def read_label():
        nonlocal pos
        if pos < len(text) and text[pos] == "'":
            end = text.index("'", pos + 1)
            label = text[pos + 1:end]
            pos = end + 1
            return label
        start = pos
        while pos < len(text) and text[pos] not in ",:)(;":
            pos += 1
        return text[start:pos].strip()

    def read_node():
        nonlocal pos
        node = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                node.children.append(read_node())
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
        node.name = read_label()
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",)(;":
                pos += 1
            try:
                node.length = float(text[start:pos])
            except ValueError:
                node.length = 0.0
        return node

    return read_node()


def tips(node):
    if not node.children:
        return [node]
    out = []
    for c in node.children:
        out.extend(tips(c))
    return out


def ladderize(node):
    for c in node.children:
        ladderize(c)
    node.children.sort(key=lambda c: len(tips(c)))


def layout(root, row_h):
    """Rectangular layout: y per tip row, x = cumulative branch length."""
    counter = [0]

    def place(node, x):
        node.x = x
        if not node.children:
            node.y = counter[0] * row_h
            counter[0] += 1
        else:
            for c in node.children:
                place(c, x + max(c.length, 0.0))
            node.y = (node.children[0].y + node.children[-1].y) / 2

    place(root, 0.0)
    return counter[0]


def render_tree_svg(newick, tip_info, width=900, row_h=11, label_w=330):
    """Return SVG for a Newick tree.

    tip_info: name -> dict(label, color, href, title, marker)
    """
    root = parse_newick(newick)
    ladderize(root)
    n = layout(root, row_h)
    xmax = max(t.x for t in tips(root)) or 1.0
    scale = (width - label_w - 40) / xmax
    H = n * row_h + 30
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {H}" width="100%" '
             f'class="tree" style="height:auto">']

    def draw(node):
        x0 = 20 + node.x * scale
        for c in node.children:
            x1 = 20 + c.x * scale
            parts.append(f'<path d="M{x0:.1f},{node.y + 15:.1f} V{c.y + 15:.1f} H{x1:.1f}" '
                         f'fill="none" stroke="#777" stroke-width="1"/>')
            draw(c)
        if not node.children:
            info = tip_info.get(node.name, {})
            color = info.get("color", "#95a5a6")
            label = html.escape(info.get("label", node.name)[:60])
            title = html.escape(info.get("title", node.name))
            y = node.y + 15
            marker = info.get("marker")
            if marker == "ref":
                parts.append(f'<rect x="{x0 - 4:.1f}" y="{y - 4:.1f}" width="8" height="8" '
                             f'fill="{color}" stroke="#111" stroke-width="1"><title>{title}</title></rect>')
            else:
                parts.append(f'<circle cx="{x0:.1f}" cy="{y:.1f}" r="2.8" fill="{color}"><title>{title}</title></circle>')
            weight = ' font-weight="700"' if marker == "ref" else ""
            text = (f'<text x="{x0 + 7:.1f}" y="{y + 3.5:.1f}" font-size="9" fill="currentColor"'
                    f'{weight}>{label}<title>{title}</title></text>')
            href = info.get("href")
            parts.append(f'<a href="{html.escape(href)}">{text}</a>' if href else text)

    draw(root)
    # scale bar: 0.5 substitutions/site
    bar = 0.5 * scale
    parts.append(f'<line x1="20" y1="{H - 8}" x2="{20 + bar:.1f}" y2="{H - 8}" stroke="currentColor" stroke-width="2"/>'
                 f'<text x="{20 + bar + 5:.1f}" y="{H - 5}" font-size="9" fill="currentColor">0.5 subst./site</text>')
    parts.append("</svg>")
    return "".join(parts), n


def tip_info_from_db(con, url, ecology=None):
    """Build the tip_info map for all representatives and references."""
    info = {}
    for cid, cluster, org, leaf in con.execute("""
        SELECT r.candidate_id, r.ro_cluster, p.organism, rl.leaf_id
        FROM ro r JOIN replicon p USING(nucleotide_id)
        LEFT JOIN ro_leaf rl ON rl.candidate_id=r.candidate_id WHERE r.is_confirmed=1"""):
        info[cid] = {"label": f"{cluster.split('_', 2)[-1]} · {org or ''}",
                     "color": GROUP_COLORS.get(group_of(cluster), "#95a5a6"),
                     "href": url("/ro/" + cid.replace(":", "_")),
                     "title": f"{cid} | {cluster} | {org} | variant {leaf}"}
    if table_exists(con, "ro_evidence"):
        for cid, tier in con.execute("SELECT candidate_id, tier FROM ro_evidence"):
            if cid in info:
                info[cid]["title"] += f" | {tier}"
    for (cluster,) in con.execute("SELECT DISTINCT ro_cluster FROM ro WHERE is_confirmed=1"):
        name = "REF|" + cluster
        sub = (ecology or {}).get(cluster, {}).get("substrate", "")
        info[name] = {"label": f"★ {cluster} reference" + (f" ({sub})" if sub else ""),
                      "color": GROUP_COLORS.get(group_of(cluster), "#95a5a6"),
                      "href": url("/cluster/" + cluster), "marker": "ref",
                      "title": f"curated reference enzyme for {cluster}"}
    # references whose cluster has no confirmed members still appear in the tree
    return info
