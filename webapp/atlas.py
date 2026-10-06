"""
Atlas data layer for the overview pages (phylogeny, network, taxonomy, ecology,
operons/ETC, evidence, statistics) plus a small Newick -> SVG tree renderer.
Reads the SQLite DB and analysis_out/ files; no extra dependencies.
"""

import csv
import functools
import html
import json
import os
import re
from collections import Counter, defaultdict

GROUP_COLORS = {"1": "#c0392b", "2": "#8e44ad", "3": "#2980b9", "4": "#27ae60",
                "5": "#f39c12", "?": "#95a5a6"}
def _relative_luminance(colour):
    """WCAG bagil parlaklik. Girdi "#rrggbb"."""
    value = colour.lstrip("#")
    if len(value) != 6:
        return 0.0
    channels = []
    for index in (0, 2, 4):
        c = int(value[index:index + 2], 16) / 255
        channels.append(c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4)
    red, green, blue = channels
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


@functools.lru_cache(maxsize=1)
def _active_site_raw(path):
    return read_json(path) or {}


def active_site_for_type(path, cluster, radius="8.0"):
    """Bir tipin aktif bolgesi: demir, cepheyi saran kalintilar, baglar.

    Tip sayfasindaki uc boyutlu goruntuleyici butun proteini gosteriyordu;
    kullanici asil ilgilenilen yerin YALNIZCA aktif bolge ve oradaki
    etkilesimler oldugunu soyledi. Bu fonksiyon o goruntuyu kurmak icin gereken
    en kucuk veriyi dondurur: demirin zinciri ve numarasi, cepheyi saran
    kalintilarin numaralari, ve demire baglanan atomlarin mesafeleri.
    """
    raw = _active_site_raw(path)
    for entry in raw.get("structures") or []:
        if entry.get("type") != cluster or entry.get("status") != "ok":
            continue
        iron = entry.get("catalytic_iron") or {}
        pocket = ((entry.get("pocket") or {}).get(radius) or {})
        residues = []
        for r in pocket.get("residues") or []:
            residues.append({
                "chain": r.get("chain"),
                "number": r.get("author_number"),
                "name": r.get("residue_3letter"),
                "one": r.get("residue_1letter"),
                "distance": r.get("min_distance_to_iron_A"),
                "column": r.get("alignment_column"),
            })
        residues.sort(key=lambda r: (r["distance"] is None, r["distance"]))
        contacts = []
        for kind, pairs in (iron.get("coordinating_contacts") or {}).items():
            for label, distance in pairs:
                contacts.append({"kind": kind.replace("_", " "),
                                 "atom": label, "distance": distance})
        contacts.sort(key=lambda c: c["distance"])
        return {
            "pdb": entry.get("pdb_id"),
            "iron_chain": iron.get("chain"),
            "iron_residue": iron.get("residue"),
            "radius": float(radius),
            "residues": residues,
            "contacts": contacts,
            "ligands": iron.get("ligating_residues") or [],
        }
    return None


def operon_relations_view(path):
    """Operon iliskileri raporunu sayfanin ihtiyaci kadarina indirir."""
    raw = read_json(path)
    if not raw:
        return None
    reg = raw.get("regulation") or {}

    assoc = []
    for item in reg.get("associations") or []:
        v, null = item.get("cramers_v"), item.get("permuted_cramers_v_mean")
        if v is None or null is None:
            continue
        assoc.append({
            "label": item.get("question") or item.get("id", "").replace("_", " "),
            "unit": item.get("unit", ""),
            "v": v, "null": null,
            "excess": item.get("cramers_v_above_chance", round(v - null, 4)),
        })
    assoc.sort(key=lambda r: -r["excess"])

    mobility = raw.get("mobility_by_substrate_class") or {}
    rows = []
    for t in mobility.get("tests") or []:
        if t.get("ratio") is None:
            continue
        rows.append(t)
    # Yapisal kanit once, sonra metinden gelen, sonra ikisinin birlesimi.
    order = {"plasmid_only": 0, "transposon_regex_only": 1, "plasmid_or_transposon": 2}
    rows.sort(key=lambda t: (order.get(t.get("evidence"), 9), t.get("level") != "entry"))

    div = raw.get("regulator_vs_enzyme_divergence") or {}
    same = (div.get("same_regulator_family_only") or {}).get("by_enzyme_identity_band") or []
    bands = []
    for b in same:
        bands.append({
            "label": b.get("label"), "n_pairs": b.get("n_pairs"), "n_types": b.get("n_types"),
            "enzyme": b.get("enzyme_identity_median"),
            "regulator": b.get("regulator_identity_median"),
            "gap": round((b.get("enzyme_identity_median") or 0)
                         - (b.get("regulator_identity_median") or 0), 1),
            "share": b.get("share_of_pairs_regulator_less_conserved"),
        })

    bias_block = raw.get("annotation_bias") or {}
    bias = []
    for key, label in (("hard_floor", "records with no annotated neighbour at all"),
                       ("by_window_occupancy", "how many genes the submitter annotated nearby"),
                       ("by_replicon_cds_count", "size of the replicon"),
                       ("structural_evidence_also_tracks_annotation_density",
                        "the same test applied to plasmid evidence")):
        block = bias_block.get(key)
        if isinstance(block, dict):
            text = (block.get("finding") or block.get("what_it_shows")
                    or block.get("note") or block.get("summary")
                    or block.get("interpretation"))
            if text:
                bias.append({"label": label, "text": text})

    return {
        "method": raw.get("method") or {},
        "totals": raw.get("totals") or {},
        "regulation": reg,
        "assoc_rows": assoc,
        "divergence": {"bands": bands,
                       "retention": (div.get("regulator_family_retention") or {})
                       .get("by_enzyme_identity_band") or []},
        "mobility": {"rows": rows,
                     "regex_only": (mobility.get("composition") or {})
                     .get("mobile_only_because_of_regex")},
        "bias": bias,
        "bias_verdict": (bias_block.get("verdict") or {}).get("size_of_the_bias", ""),
        "limits": raw.get("limits") or [],
    }


def active_site_view(path):
    """2,6 MB'lik yapi raporunu SAYFANIN ihtiyaci kadarina indirir.

    Dosyanin tamami tarayiciya gonderilmez: her yapinin tam kalinti listesi ve
    her sutunun uye dagilimi onun buyuk kismini olusturuyor ve sayfada
    gosterilmiyor. Burada yalnizca ozetler cikariliyor.
    """
    raw = read_json(path)
    if not raw:
        return None
    cover = raw.get("coverage") or {}
    validation = (raw.get("predicted_pocket_transfer_validation") or {}).get("summary") or {}

    universal = []
    for radius in ("5.0", "8.0"):
        block = (raw.get("pocket_columns") or {}).get(radius) or {}
        for column in block.get("columns", []):
            universal.append({
                "radius": float(radius),
                "column": column.get("column"),
                "n_types": column.get("n_types"),
                "residues": column.get("residues") or column.get("residue") or {},
            })

    # Uye degiskenligi: hangi sutun korunuyor, hangisi oynuyor.
    variation = []
    merged = {}
    for source, store in (("crystal", raw.get("pocket_variation_across_members") or {}),
                          ("predicted", raw.get("predicted_pocket_variation_across_members") or {})):
        for cluster, info in store.items():
            if info.get("status") != "measured":
                continue
            # Kristal olcumu VARSA onu tercih et; tahmin ikincil.
            if cluster in merged and merged[cluster]["source"] == "crystal":
                continue
            merged[cluster] = {
                "cluster": cluster, "source": source,
                "members": info.get("n_members"),
                "columns": info.get("n_pocket_columns"),
                "invariant": len(info.get("invariant_columns") or []),
                "variable": len(info.get("variable_columns") or []),
            }
    for row in merged.values():
        if row["columns"]:
            row["invariant_share"] = round(row["invariant"] / row["columns"], 3)
            variation.append(row)
    variation.sort(key=lambda r: (r["invariant_share"], -(r["members"] or 0)))

    charge = []
    per_type = (raw.get("electrostatics") or {}).get("per_type") or []
    rows = per_type.values() if isinstance(per_type, dict) else per_type
    for item in rows:
        comp = item.get("composition") or {}
        charge.append({
            "cluster": item.get("type"), "pdb": item.get("pdb_id"),
            "family": (item.get("family") or "").replace("_", " "),
            "substrate": item.get("substrate"),
            "residues": comp.get("n_residues"),
            "acidic": comp.get("acidic_D_E"), "basic": comp.get("basic_K_R"),
            "aromatic": comp.get("aromatic_F_W_Y") or comp.get("aromatic"),
            "net": comp.get("net_formal_charge"),
            "net_excl": comp.get("net_formal_charge_excluding_iron_ligands"),
        })
    charge.sort(key=lambda r: (r["net"] if r["net"] is not None else 0))

    return {
        "coverage": cover,
        "validation": validation,
        "universal": universal,
        "variation": variation,
        "charge": charge,
        "not_electrostatics": (raw.get("electrostatics") or {}).get("explicitly_not_done"),
        "limitations": raw.get("limitations") or [],
        "n_predicted": len(raw.get("predicted_structures") or []),
        "n_predicted_failed": len(raw.get("predicted_undetermined") or []),
        "parameters": raw.get("parameters") or {},
    }


def assignment_agreement(con):
    """Her tip icin: uyelerin kaci GERCEKTEN bu tipin referansina en yakin?

    NEDEN. Atama HMM bit skoruna gore yapiliyor, "en yakin referans" ise dizi
    kimligine gore olculuyor. Ikisi ayri sorulardir ve AYRILABILIRLER. Kullanici
    bunu qxyA'da fark etti: tipe atanan 95 uyenin kuratorlu referansa kimligi
    ortanca %33, ve referansa en yakin 127 girisin yalnizca 81'i o tipe
    atanmis; kalani ayni grubun baska tiplerine dagilmis.

    Olculdu, genel resim: 43 tipin 16'sinda uyelerin COGU baska bir referansa
    daha yakin. Ama bu, profilin "tanimlayan proteini disarida biraktigi"
    anlamina gelmiyor: kimligin YUKSEK oldugu yerde (>=%60) uyusmazlik 2.868
    girisin yalnizca 121'inde (%4,2) ve bunlarin cogu zaten belgelenmis ikiz
    referans ciftleri. Uyusmazlik dusuk kimlikte yogunlasiyor, yani profiller
    uzak akrabalari ayirt edemedigi icin.

    Bu sayi tip sayfasinda gosteriliyor, cunku "bu tipin 610 uyesi var"
    cumlesi, uyelerin %70'i baska bir referansa daha yakinken yaniltici olur.
    """
    if not table_exists(con, "ro_evidence"):
        return {}
    rows = con.execute("""
        SELECT r.ro_cluster,
               COUNT(*),
               SUM(CASE WHEN e.nearest_ref = r.ro_cluster THEN 1 ELSE 0 END),
               AVG(e.ref_identity)
        FROM ro r JOIN ro_evidence e USING(candidate_id)
        WHERE r.is_confirmed = 1 AND r.ro_cluster IS NOT NULL
        GROUP BY 1""").fetchall()
    out = {}
    for cluster, total, agree, mean_identity in rows:
        if not total:
            continue
        out[cluster] = {
            "members": total,
            "nearest_is_own": agree,
            "agreement": round(agree / total, 4),
            "mean_identity_to_nearest": round(mean_identity or 0, 1),
        }
    return out


# Kaynak alanindaki tanimlayicilar BAGLANTIYA cevrilir. Alan serbest metin:
# icinde PDB kimligi, PMID, DOI, PMC numarasi ve bazen duz adres geciyor.
# Okuyucunun bir referansi takip edebilmesi icin bunlari elle aramasi
# gerekiyordu; desenler burada tek yerde tanimli, boylece her sayfada ayni
# sekilde calisiyor.
_CITATION_PATTERNS = (
    (re.compile(r"\bPMID[:\s]+(\d{4,9})\b", re.I),
     lambda m: ("https://pubmed.ncbi.nlm.nih.gov/%s/" % m.group(1), m.group(0))),
    (re.compile(r"\b(PMC\d{5,9})\b"),
     lambda m: ("https://pmc.ncbi.nlm.nih.gov/articles/%s/" % m.group(1), m.group(0))),
    (re.compile(r"\bdoi[:\s]+(10\.\d{4,9}/[^\s,;)\]]+)", re.I),
     lambda m: ("https://doi.org/%s" % m.group(1), m.group(0))),
    # PDB kimligi: rakamla baslayan dort karakter. Desen DAR tutuldu, cunku
    # dort harfli her kelimeyi yapi sanmak yanlis baglantilar uretirdi.
    (re.compile(r"\bPDB\s+([0-9][A-Za-z0-9]{3})\b"),
     lambda m: ("https://www.rcsb.org/structure/%s" % m.group(1).upper(), m.group(0))),
    (re.compile(r"(https?://[^\s,;)\]]+)"),
     lambda m: (m.group(1), m.group(1))),
)


def linkify_citation(text):
    """Kaynak metnindeki tanimlayicilari tiklanabilir yapar.

    Metin once KACISLANIR, sonra baglantilar eklenir; ters sirada yapmak
    uretilen etiketleri de kacislardi.
    """
    if not text:
        return ""
    escaped = html.escape(str(text))
    spans = []
    for pattern, build in _CITATION_PATTERNS:
        for match in pattern.finditer(escaped):
            if any(start < match.end() and match.start() < end for start, end, _ in spans):
                continue          # ust uste binen eslesmeyi atla
            url, label = build(match)
            spans.append((match.start(), match.end(), (url, label)))
    if not spans:
        return escaped
    spans.sort()
    out, cursor = [], 0
    for start, end, (url, label) in spans:
        out.append(escaped[cursor:start])
        out.append(f'<a href="{html.escape(url, quote=True)}" rel="noopener noreferrer" '
                   f'target="_blank">{label}</a>')
        cursor = end
    out.append(escaped[cursor:])
    return "".join(out)


def chip_text(background):
    """Arka plana gore OKUNUR metin rengi dondurur.

    NEDEN. Etiketlerin arka plani VERIDEN geliyor (TIER_COLORS,
    REACTION_COLORS) ama metin rengi sablonlarda "color:#fff" olarak SABIT
    yaziliydi. Olctum: alti etiket WCAG oraninda basarisiz, en kotusu
    "family member" 1,68:1 ve "close homolog" 2,53:1 -- yani okunmuyordu, ve
    en cok Evidence sayfasinda, yani konusu tam o etiketler olan sayfada.

    Renkleri DEGISTIRMEK yerine metni cevirmek secildi: renkler sayfalar
    arasinda anlam tasiyor ve kullanici onlari taniyor. Secim parlaklıktan
    HESAPLANIYOR, bir renk listesine bakilmiyor; boylece yarin bir ton
    degistirilirse kural kendiliginden dogru kaliyor. Sabit bir hex listesine
    bakan bir cozum sessizce bozulurdu.
    """
    dark_text, light_text = "#10161c", "#ffffff"
    background_luminance = _relative_luminance(background)

    def contrast(foreground):
        a, b = _relative_luminance(foreground), background_luminance
        high, low = max(a, b), min(a, b)
        return (high + 0.05) / (low + 0.05)

    return dark_text if contrast(dark_text) > contrast(light_text) else light_text


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
        d["euk"] = 0
        d["leaves"] = 0
        rows[d["cluster"]] = d
    for cluster, n in con.execute("SELECT cluster, COUNT(*) FROM leaf GROUP BY cluster"):
        if cluster in rows:
            rows[cluster]["leaves"] = n
    if table_exists(con, "ro_domain"):
        for cluster, n in con.execute(
                "SELECT cluster, SUM(domain = 'Eukaryota') FROM ro_domain GROUP BY cluster"):
            if cluster in rows:
                rows[cluster]["euk"] = n or 0
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


def representative_coverage(con):
    """Agacta ve agda KAC varyantin temsil edildigi, ve bedeli.

    Temsilci secimi iki uclu bir kural: tekil varyant kendisi temsilcidir, 5
    ve uzeri uye tasiyan varyant medoid verir, arada kalan 2-4 uyeli varyantlar
    ATLANIR. Sayfalar bir donem "her varyant icin bir temsilci" diyordu, oysa
    1.809 varyantin 1.205'i temsil ediliyor. Atlanan 604 varyant girislerin
    %14'unu tasiyor, dolayisiyla bu susulacak bir ayrinti degil.
    """
    if not table_exists(con, "leaf"):
        return None
    bands = {}
    for band, leaves, entries in con.execute("""
            SELECT CASE WHEN size = 1 THEN 'singleton'
                        WHEN size < 5 THEN 'small'
                        ELSE 'large' END,
                   COUNT(*), SUM(size)
            FROM leaf GROUP BY 1"""):
        bands[band] = {"variants": leaves, "entries": entries}
    total_variants = sum(b["variants"] for b in bands.values())
    total_entries = sum(b["entries"] for b in bands.values())
    represented = (bands.get("singleton", {}).get("variants", 0)
                   + bands.get("large", {}).get("variants", 0))
    skipped = bands.get("small", {"variants": 0, "entries": 0})
    missing_types = [row[0] for row in con.execute("""
        SELECT DISTINCT cluster FROM leaf
        WHERE cluster NOT IN (SELECT cluster FROM leaf WHERE size = 1 OR size >= 5)
        ORDER BY 1""")]
    return {
        "variants_total": total_variants,
        "variants_represented": represented,
        "variants_skipped": skipped["variants"],
        "entries_total": total_entries,
        "entries_skipped": skipped["entries"],
        "entries_skipped_share": (skipped["entries"] / total_entries) if total_entries else 0,
        "min_members_for_medoid": 5,
        "types_without_representative": missing_types,
        "bands": bands,
    }


def tip_info_from_db(con, url, ecology=None, names=None):
    """Build the tip_info map for all representatives and references.

    `names` cakisan kisa adlari ayristirilmis haliyle verir; verilmezse sade
    kisa ad kullanilir, boylece modul tek basina da calisir.
    """
    names = names or {}

    def gene(cluster):
        return names.get(cluster, cluster.split("_", 2)[-1])

    info = {}
    for cid, cluster, org, leaf in con.execute("""
        SELECT r.candidate_id, r.ro_cluster, p.organism, rl.leaf_id
        FROM ro r JOIN replicon p USING(nucleotide_id)
        LEFT JOIN ro_leaf rl ON rl.candidate_id=r.candidate_id WHERE r.is_confirmed=1"""):
        info[cid] = {"label": f"{gene(cluster)} · {org or ''}",
                     "color": GROUP_COLORS.get(group_of(cluster), "#95a5a6"),
                     "href": url("/ro/" + cid.replace(":", "_")),
                     "title": f"{cid} | {cluster} | {org} | variant {leaf}"}
    if table_exists(con, "ro_evidence"):
        for cid, tier in con.execute("SELECT candidate_id, tier FROM ro_evidence"):
            if cid in info:
                info[cid]["title"] += f" | {tier}"
    # agac/ag ucu adlari YAPRAK kimligidir ("3_313_KshA15#0"); temsilci girise baglanir
    for leaf_id, rep, cluster, size, ident in con.execute(
            "SELECT leaf_id, representative, cluster, size, median_identity FROM leaf"):
        org = con.execute("SELECT p.organism FROM ro r JOIN replicon p USING(nucleotide_id) "
                          "WHERE r.candidate_id=?", (rep,)).fetchone()
        org = org[0] if org else ""
        info[leaf_id] = {
            "label": f"{gene(cluster)}#{leaf_id.split('#')[-1]} · {org or ''} (n={size})",
            "color": GROUP_COLORS.get(group_of(cluster), "#95a5a6"),
            "href": url("/leaf/" + leaf_id.replace("#", "-")),
            "title": f"variant {leaf_id} | {size} members | median identity "
                     f"{ident:.2f} | representative {org}"}
    for (cluster,) in con.execute("SELECT DISTINCT ro_cluster FROM ro WHERE is_confirmed=1"):
        name = "REF|" + cluster
        sub = (ecology or {}).get(cluster, {}).get("substrate", "")
        info[name] = {"label": f"★ {cluster} reference" + (f" ({sub})" if sub else ""),
                      "color": GROUP_COLORS.get(group_of(cluster), "#95a5a6"),
                      "href": url("/cluster/" + cluster), "marker": "ref",
                      "title": f"curated reference enzyme for {cluster}"}
    # references whose cluster has no confirmed members still appear in the tree
    return info


# ----------------------------------------------------------- schematic SVGs
COMPONENT_FILL = {
    "alpha": "#c0392b", "[alpha]": "#c0392b", "alpha_other": "#e67e22", "beta": "#2980b9",
    "ferredoxin": "#27ae60", "reductase": "#f39c12", "rieske_other": "#16a085",
    "regulator": "#8e44ad", "transposon": "#2c3e50", "ring_cleavage": "#d35400",
    "transporter": "#7f8c8d", "dehydrogenase": "#95a5a6", "hydrolase": "#9aa5a8",
    "hypothetical": "#d8dcdd", "other": "#b2bec3", "none": "#b2bec3",
}
SHORT = {"[alpha]": "α", "alpha": "α", "alpha_other": "α′", "beta": "β", "ferredoxin": "Fd",
         "reductase": "Red", "rieske_other": "Rsk", "regulator": "Reg", "transposon": "IS",
         "ring_cleavage": "ring", "transporter": "tra", "dehydrogenase": "dh",
         "hydrolase": "hyd", "hypothetical": "?", "other": "·", "none": "·"}


def layout_svg(layout, height=34, gene_w=46, gap=3):
    """Render an operon layout string ('reductase > [alpha] > beta') as gene arrows, 5'->3'."""
    if not layout:
        return ""
    tokens = [t.strip() for t in layout.split(">") if t.strip()]
    width = len(tokens) * (gene_w + gap) + 10
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
             f'height="{height}" width="{width}" class="oplayout">']
    x = 5
    for token in tokens:
        fill = COMPONENT_FILL.get(token, "#b2bec3")
        is_alpha = token in ("[alpha]", "alpha")
        head = 9
        y0, y1 = (4, height - 10) if is_alpha else (7, height - 13)
        pts = (f"{x},{y0} {x + gene_w - head},{y0} {x + gene_w},{(y0 + y1) / 2:.1f} "
               f"{x + gene_w - head},{y1} {x},{y1}")
        stroke = "#111" if is_alpha else "#5a5a5a"
        label = SHORT.get(token, token[:3])
        parts.append(
            f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{1.4 if is_alpha else 0.7}"><title>{html.escape(token)}</title></polygon>'
            f'<text x="{x + (gene_w - head) / 2:.1f}" y="{(y0 + y1) / 2 + 3.5:.1f}" font-size="10" '
            f'text-anchor="middle" fill="#fff" font-weight="600">{html.escape(label)}</text>')
        x += gene_w + gap
    parts.append("</svg>")
    return "".join(parts)


def gap_histogram(values, bins=(0, 25, 50, 100, 150, 200, 300, 400, 600, 1000, 2000)):
    """Bucket intergenic distances for a compact bar chart."""
    counts = [0] * (len(bins) - 1)
    for v in values:
        for i in range(len(bins) - 1):
            if bins[i] <= v < bins[i + 1]:
                counts[i] += 1
                break
    labels = [f"{bins[i]}–{bins[i + 1]}" for i in range(len(bins) - 1)]
    return labels, counts


def median(values):
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    return float(s[n // 2]) if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def _field(row, key, default=None):
    """Read a key from a sqlite3.Row or a plain dict."""
    try:
        value = row[key]
    except (KeyError, IndexError, TypeError):
        return default
    return default if value is None else value


REG_SHORT = {"LysR": "LysR", "TetR": "TetR", "AraC": "AraC", "MarR": "MarR", "IclR": "IclR",
             "GntR": "GntR", "LuxR": "LuxR", "ArsR": "ArsR", "Crp_Fnr": "Crp", "sigma54": "s54",
             "two_component": "2comp", "sigma_factor": "sig", "regulator_unclassified": "Reg"}


def operon_regulator_svg(layout, reg=None, gene_w=46, gap=3):
    """Operon arrows preceded by the upstream regulator and the promoter region.

    The regulator is drawn pointing away from the operon when the architecture is
    divergent, the arrangement in which both share one intergenic promoter region.
    """
    tokens = [t.strip() for t in (layout or "").split(">") if t.strip()]
    if not tokens:
        return ""
    arch = _field(reg, "architecture", "") if reg is not None else ""
    family = _field(reg, "upstream_family", "") if reg is not None else ""
    bp = _field(reg, "intergenic_bp", None) if reg is not None else None
    show_reg = arch in ("divergent_regulator", "codirectional_regulator")
    divergent = arch == "divergent_regulator"

    pad, prom_w, row_h, height = 6, 54, 24, 72
    y0 = 20
    y1 = y0 + row_h
    mid = (y0 + y1) / 2.0
    width = pad * 2 + ((prom_w + gap) if (show_reg or bp) else 0) \
        + ((gene_w + gap) if show_reg else 0) + len(tokens) * (gene_w + gap)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
             f'height="{height}" width="{min(width, 980)}" class="oplayout">']
    x = pad
    head = 9

    if show_reg:
        label = REG_SHORT.get(family, "Reg")
        if divergent:
            pts = f"{x + gene_w},{y0} {x + head},{y0} {x},{mid:.1f} {x + head},{y1} {x + gene_w},{y1}"
        else:
            pts = f"{x},{y0} {x + gene_w - head},{y0} {x + gene_w},{mid:.1f} {x + gene_w - head},{y1} {x},{y1}"
        product = html.escape(str(_field(reg, "upstream_product", "regulator")))
        parts.append(
            f'<polygon points="{pts}" fill="#8e44ad" stroke="#4a235a" stroke-width="0.9">'
            f'<title>{product}</title></polygon>'
            f'<text x="{x + gene_w / 2:.1f}" y="{mid + 3.5:.1f}" font-size="9" text-anchor="middle" '
            f'fill="#fff" font-weight="600">{html.escape(label)}</text>'
            f'<text x="{x + gene_w / 2:.1f}" y="{y0 - 6:.1f}" font-size="9" text-anchor="middle" '
            f'fill="currentColor">{"divergent" if divergent else "same strand"}</text>')
        x += gene_w + gap

    if show_reg or bp:
        gap_label = (str(int(bp)) + " bp") if bp is not None else "promoter"
        parts.append(
            f'<rect x="{x}" y="{y0 - 2}" width="{prom_w}" height="{row_h + 4}" fill="none" '
            f'stroke="#8c2f22" stroke-dasharray="4,3" stroke-width="0.9"/>'
            f'<text x="{x + prom_w / 2:.1f}" y="{y1 + 14:.1f}" font-size="9" text-anchor="middle" '
            f'fill="currentColor">{gap_label}</text>')
        x += prom_w + gap

    for token in tokens:
        fill = COMPONENT_FILL.get(token, "#b2bec3")
        is_alpha = token in ("[alpha]", "alpha")
        ay0, ay1 = (y0 - 2, y1 + 2) if is_alpha else (y0, y1)
        amid = (ay0 + ay1) / 2.0
        pts = (f"{x},{ay0} {x + gene_w - head},{ay0} {x + gene_w},{amid:.1f} "
               f"{x + gene_w - head},{ay1} {x},{ay1}")
        parts.append(
            f'<polygon points="{pts}" fill="{fill}" '
            f'stroke="{"#111" if is_alpha else "#5a5a5a"}" '
            f'stroke-width="{1.4 if is_alpha else 0.7}">'
            f'<title>{html.escape(token)}</title></polygon>'
            f'<text x="{x + (gene_w - head) / 2:.1f}" y="{amid + 3.5:.1f}" font-size="10" '
            f'text-anchor="middle" fill="#fff" font-weight="600">'
            f'{html.escape(SHORT.get(token, token[:3]))}</text>')
        x += gene_w + gap

    parts.append(f'<text x="{pad}" y="{height - 4}" font-size="9" fill="currentColor">'
                 f'operon drawn 5′ to 3′; dashed box is the intergenic promoter region</text>')
    parts.append("</svg>")
    return "".join(parts)


# ------------------------------------------------------------- chemistry
FAMILY_ORDER = ["alkylbenzenes", "pah", "biaryls_ethers", "nitroaromatics", "haloaromatics",
                "sulfoaromatics", "aromatic_acids", "anilines", "quaternary_amines",
                "alkaloids", "terpenoids_steroids", "unknown"]
FAMILY_LABEL = {
    "alkylbenzenes": "Benzene and alkylbenzenes",
    "pah": "Polycyclic aromatic hydrocarbons",
    "biaryls_ethers": "Biaryls, diaryl ethers and heterocycles",
    "nitroaromatics": "Nitroaromatics",
    "haloaromatics": "Halogenated aromatics",
    "sulfoaromatics": "Sulfonated aromatics",
    "aromatic_acids": "Aromatic acids and lignin-derived compounds",
    "anilines": "Anilines",
    "quaternary_amines": "Quaternary amines and osmolytes",
    "alkaloids": "Alkaloids and purines",
    "terpenoids_steroids": "Terpenoids and steroids",
    "unknown": "Substrate not established",
}
FAMILY_NOTE = {
    "alkylbenzenes": "Fuel components and industrial solvents. These are the substrates on which "
                     "the chemistry of the family was first described, and toluene dioxygenase "
                     "remains its structural prototype.",
    "pah": "Combustion and tar residues, several of them priority pollutants. Their fused rings "
           "resist attack until two hydroxyls are installed on adjacent carbons.",
    "biaryls_ethers": "Biphenyls, dioxins, carbazole and lignin-derived biaryls. Some of these are "
                      "attacked at the carbon joining the two rings rather than on a ring edge, "
                      "which is what makes the skeleton fall apart.",
    "nitroaromatics": "Explosives and dye intermediates. The oxygenase installs the diol and the "
                      "nitro group leaves as nitrite, so one step both activates and detoxifies.",
    "haloaromatics": "Herbicides, solvents and their residues. The halogen departs as halide, "
                     "which is the step that makes these compounds biodegradable at all.",
    "sulfoaromatics": "Surfactant and dye intermediates. The sulfonate leaves as sulfite during "
                      "dihydroxylation.",
    "aromatic_acids": "Carboxylated aromatics from plants and from industry, including the "
                      "polyethylene terephthalate monomer. Most funnel into the beta-ketoadipate "
                      "pathway and then into central metabolism.",
    "anilines": "Dye and pesticide feedstocks. The amino group leaves as ammonia, giving catechol.",
    "quaternary_amines": "Choline, the betaines and carnitine. These enzymes carry the same two "
                         "metal centres as the ring-hydroxylating ones but act on open-chain "
                         "substrates with no aromatic ring at all, which is the clearest evidence "
                         "that the family is defined by its chemistry rather than by its substrate.",
    "alkaloids": "Caffeine and intermediates of saxitoxin biosynthesis. Here the enzymes strip "
                 "methyl groups from nitrogen or install single hydroxyls on complex scaffolds.",
    "terpenoids_steroids": "Conifer resin acids and 3-ketosteroids. Hydroxylation at a specific "
                           "ring position opens polycyclic skeletons that are otherwise inert.",
    "unknown": "Reference enzymes whose original publications could not be traced in this work. "
               "Their members are confirmed Rieske oxygenases, but this database claims no "
               "substrate for them.",
}
REACTION_LABEL = {
    "cis_dihydroxylation": "cis-dihydroxylation",
    "angular_dioxygenation": "angular dioxygenation",
    "dioxygenation_with_release": "dioxygenation with substituent release",
    "O_demethylation": "O-demethylation",
    "N_demethylation": "N-demethylation",
    "hydroxylation": "hydroxylation",
    "C_N_cleavage": "carbon-nitrogen bond cleavage",
    "unknown": "not established",
}
REACTION_NOTE = {
    "cis_dihydroxylation": "Both atoms of O2 are added to adjacent ring carbons, giving a cis-diol. "
                           "This is the reaction the family is named for.",
    "angular_dioxygenation": "Attack at the carbon that joins two rings. The product is unstable "
                             "and the ring system opens without further enzymes.",
    "dioxygenation_with_release": "Dihydroxylation at a substituted carbon, so the substituent "
                                  "leaves as nitrite, halide, sulfite, ammonia or carbon dioxide.",
    "O_demethylation": "One oxygen atom is used to oxidise an aryl methyl ether, releasing "
                       "formaldehyde and the free phenol.",
    "N_demethylation": "A methyl group on nitrogen is oxidised and released as formaldehyde.",
    "hydroxylation": "A single oxygen atom is inserted at one carbon; the second is reduced to "
                     "water.",
    "C_N_cleavage": "Oxygenation at the carbon next to a quaternary nitrogen, cleaving the "
                    "carbon-nitrogen bond.",
    "unknown": "No reaction is claimed for these types.",
}
REACTION_COLORS = {
    "cis_dihydroxylation": "#1f4e79", "angular_dioxygenation": "#2a8f87",
    "dioxygenation_with_release": "#8c2f22", "O_demethylation": "#b5731f",
    "N_demethylation": "#6d4c9a", "hydroxylation": "#1f7a4d",
    "C_N_cleavage": "#b03a6e", "unknown": "#8a8a8a",
}


def short_names(path):
    """Tip id'sinden gosterilecek kisa ad. Cakisanlara grup-numara eklenir.

    Uc kisa ad referans setinde IKI kez geciyor: BphA1 (2_201 ve 2_218),
    NDO (3_314 ve 3_315) ve NidA (3_317 ve 3_318). Her ikisini de ayni
    etiketle gostermek, iki ayri sayfaya giden iki satiri ayirt edilemez
    hale getiriyordu. Cakisma yoksa sade kisa ad kullanilir, cunku 68 tipin
    hepsine numara eklemek okunurlugu bosa dusurur.
    """
    clusters = [r["cluster"] for r in read_csv(path)]
    counts = Counter(c.split("_", 2)[-1] for c in clusters)
    names = {}
    for cluster in clusters:
        parts = cluster.split("_", 2)
        gene = parts[-1]
        names[cluster] = gene if counts[gene] == 1 else f"{gene} ({parts[1]})"
    return names


def chemistry_groups(path, counts, ecology=None):
    """Group the curated chemistry table by chemical family, newest counts attached."""
    rows = read_csv(path)
    names = short_names(path)
    by_family = defaultdict(list)
    for r in rows:
        r = dict(r)
        r["n"] = counts.get(r["cluster"], 0)
        r["gene"] = names.get(r["cluster"], r["cluster"].split("_", 2)[-1])
        r["group"] = group_of(r["cluster"])
        r["reaction_label"] = REACTION_LABEL.get(r["reaction_class"], r["reaction_class"])
        r["color"] = REACTION_COLORS.get(r["reaction_class"], "#8a8a8a")
        if ecology:
            r["sclass"] = ecology.get(r["cluster"], {}).get("substrate_class", "")
        by_family[r["family"]].append(r)
    groups = []
    for family in FAMILY_ORDER:
        members = sorted(by_family.get(family, []), key=lambda r: -r["n"])
        if not members:
            continue
        groups.append({
            "key": family, "label": FAMILY_LABEL.get(family, family),
            "note": FAMILY_NOTE.get(family, ""), "types": members,
            "n_types": len(members), "n_members": sum(m["n"] for m in members),
            "reactions": sorted({m["reaction_class"] for m in members}),
        })
    return groups


def reaction_summary(groups):
    """Members and types per reaction class, for the overview strip."""
    per = defaultdict(lambda: {"types": 0, "members": 0})
    for g in groups:
        for t in g["types"]:
            per[t["reaction_class"]]["types"] += 1
            per[t["reaction_class"]]["members"] += t["n"]
    out = []
    for key in ["cis_dihydroxylation", "dioxygenation_with_release", "angular_dioxygenation",
                "hydroxylation", "N_demethylation", "O_demethylation", "C_N_cleavage", "unknown"]:
        if key in per:
            out.append({"key": key, "label": REACTION_LABEL[key], "note": REACTION_NOTE[key],
                        "color": REACTION_COLORS[key], **per[key]})
    return out


# ------------------------------------------------- generic reaction schemes
#
# Substrat basina URUN YAPISI cizilmiyor: bu, her tipte hangi halka
# pozisyonunun saldiriya ugradigini varsaymayi gerektirir ve yanlis regiokimya
# gostermek hic gostermemekten kotudur. Onun yerine her REAKSIYON SINIFI icin
# genel mekanizma cizilir; bunlar substrattan bagimsiz ve kesindir.

def _hex_points(cx, cy, r):
    import math
    return [(cx + r * math.cos(math.radians(90 + 60 * i)),
             cy - r * math.sin(math.radians(90 + 60 * i))) for i in range(6)]


def _arene(cx, cy, r, aromatic=True, saturated=(), subs=(), fused=False):
    """Bir halka ciz. subs: [(vertex_index, 'OH')], saturated: sp3 kose indeksleri."""
    pts = _hex_points(cx, cy, r)
    poly = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    out = [f'<polygon points="{poly}" fill="none" stroke="currentColor" stroke-width="1.3"/>']
    if aromatic:
        out.append(f'<circle cx="{cx}" cy="{cy}" r="{r * 0.55:.1f}" fill="none" '
                   f'stroke="currentColor" stroke-width="1.1"/>')
    else:
        # kalan cift baglari ic cizgi olarak ciz (sp3 kosesi olmayan kenarlar)
        for i in range(6):
            j = (i + 1) % 6
            if i in saturated or j in saturated:
                continue
            if i % 2:
                continue
            x1, y1 = pts[i]
            x2, y2 = pts[j]
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            dx, dy = (cx - mx) * 0.22, (cy - my) * 0.22
            out.append(f'<line x1="{x1 + dx:.1f}" y1="{y1 + dy:.1f}" x2="{x2 + dx:.1f}" '
                       f'y2="{y2 + dy:.1f}" stroke="currentColor" stroke-width="1.1"/>')
    for index, label in subs:
        x, y = pts[index % 6]
        ox = (x - cx) * 0.55
        oy = (y - cy) * 0.55
        anchor = "start" if ox > 1 else ("end" if ox < -1 else "middle")
        out.append(f'<text x="{x + ox:.1f}" y="{y + oy + 3.5:.1f}" font-size="10" '
                   f'text-anchor="{anchor}" fill="currentColor">{html.escape(label)}</text>')
    return "".join(out)


def _arrow(x, y, width=56, above="", below="", marker_id="rxarrow"):
    out = [f'<line x1="{x}" y1="{y}" x2="{x + width}" y2="{y}" stroke="currentColor" '
           f'stroke-width="1.3" marker-end="url(#{marker_id})"/>']
    if above:
        out.append(f'<text x="{x + width / 2:.1f}" y="{y - 6}" font-size="9" '
                   f'text-anchor="middle" fill="currentColor">{above}</text>')
    if below:
        out.append(f'<text x="{x + width / 2:.1f}" y="{y + 13}" font-size="9" '
                   f'text-anchor="middle" fill="currentColor">{below}</text>')
    return "".join(out)


def _defs(marker_id="rxarrow"):
    """Ok isaretcisi tanimi. Marker id'si CAGRI BASINA tekil olmali: ana sayfada
    yedi reaksiyon semasi var ve hepsi ayni id'yi yazarsa belge gecersiz olur.
    Tarayici ilk tanima cozdugu icin gorsel sonuc dogru gorunur, bu yuzden
    gozle farkedilmez; dogrulama araci yakalar."""
    return ('<defs><marker id="' + marker_id + '" viewBox="0 0 10 10" refX="9" refY="5" '
            'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            '<path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>')


# Reaksiyon semalarinin ACIKLAMA cumleleri. SVG'nin ICINDE DEGIL burada
# duruyorlar ve sayfada HTML olarak yaziliyorlar. Sebep olculdu: cumleler SVG
# icinde sabit x konumlarina yaziliyordu ve ornegin "cis-diol, ring no longer
# aromatic" etiketi x=196'dan baslayip yaklasik 340 piksele uzaniyordu, oysa
# viewBox 300 piksel genisti. Sonuc: metin ya kirpiliyor ya da alttaki satirla
# ust uste biniyordu. HTML metni kendiliginden satir kirar, olceklenir ve
# secilebilir; SVG'de her satirin yerini elle hesaplamak gerekir. Cizim SVG'nin
# isi, kelime HTML'in isi.
REACTION_SCHEME_NOTE = {
    "cis_dihydroxylation":
        "Both oxygen atoms are added to the same face of the ring, giving a "
        "cis-diol. The ring is no longer aromatic.",
    "angular_dioxygenation":
        "The attack is at the ring junction, so the fused ring system comes "
        "apart rather than simply gaining hydroxyls.",
    "dioxygenation_with_release":
        "The substituent leaves during the reaction, as nitrite, a halide, "
        "sulfite, ammonia or carbon dioxide.",
    "O_demethylation":
        "The methyl group is released as formaldehyde, leaving a free phenol.",
    "N_demethylation":
        "One methyl group is removed per catalytic cycle and released as "
        "formaldehyde.",
    "hydroxylation":
        "One oxygen atom is inserted into the substrate; the other is reduced "
        "to water.",
    "C_N_cleavage":
        "The carbon-nitrogen bond is broken, releasing the amine and leaving "
        "an aldehyde.",
}


def reaction_scheme_svg(kind, width=262, height=86):
    """Reaksiyon sinifi icin genel mekanizma semasi."""
    r, cy = 19, 46
    # Tekil marker id: ayni sayfadaki yedi sema ayni id'yi paylasmasin.
    marker = "rxarrow-" + kind.replace("_", "-")
    arrow = functools.partial(_arrow, marker_id=marker)
    body = ""
    if kind == "cis_dihydroxylation":
        body = (_arene(34, cy, r) + arrow(62, cy, 58, "O₂, NAD(P)H")
                + _arene(150, cy, r, aromatic=False, saturated=(4, 5),
                         subs=[(4, "OH"), (5, "OH")]))
    elif kind == "angular_dioxygenation":
        body = (_arene(30, cy, r) + _arene(30 + r * 1.73, cy, r)
                + f'<circle cx="{30 + r * 0.87:.1f}" cy="{cy}" r="3" fill="#8c2f22"/>'
                + arrow(96, cy, 52, "O₂, NAD(P)H")
                + _arene(176, cy, r, subs=[(1, "OH")])
                + _arene(176 + r * 1.9, cy, r, subs=[(4, "OH")]))
    elif kind == "dioxygenation_with_release":
        body = (_arene(34, cy, r, subs=[(0, "X")]) + arrow(62, cy, 58, "O₂, NAD(P)H", "− X")
                + _arene(150, cy, r, aromatic=True, subs=[(0, "OH"), (1, "OH")]))
    elif kind == "O_demethylation":
        body = (_arene(34, cy, r, subs=[(0, "OCH₃")]) + arrow(76, cy, 54, "O₂, NAD(P)H", "− HCHO")
                + _arene(164, cy, r, subs=[(0, "OH")]))
    elif kind == "N_demethylation":
        body = ('<text x="18" y="52" font-size="13" fill="currentColor">R₂N–CH₃</text>'
                + arrow(86, cy, 54, "O₂, NAD(P)H", "− HCHO")
                + '<text x="150" y="52" font-size="13" fill="currentColor">R₂N–H</text>')
    elif kind == "hydroxylation":
        body = ('<text x="22" y="52" font-size="13" fill="currentColor">R–CH</text>'
                + arrow(74, cy, 54, "O₂, NAD(P)H", "+ H₂O")
                + '<text x="138" y="52" font-size="13" fill="currentColor">R–C–OH</text>')
    elif kind == "C_N_cleavage":
        body = ('<text x="12" y="52" font-size="13" fill="currentColor">R–CH₂–N⁺(CH₃)₃</text>'
                + arrow(116, cy, 50, "O₂, NAD(P)H")
                + '<text x="172" y="46" font-size="11" fill="currentColor">R–CHO</text>'
                + '<text x="172" y="60" font-size="11" fill="currentColor">+ N(CH₃)₃</text>')
    else:
        return ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
            f'width="100%" style="max-width:{width}px;height:auto" class="rxnscheme">'
            + _defs(marker) + body + "</svg>")


# ----------------------------------------------------------- taxonomy tree
def taxonomy_tree(con, max_children=25, min_count=2, max_depth=7):
    """NCBI soy dizgilerinden ic ice gecmis sayim agaci.

    Yigili cubuk grafigi "hangi filum ne kadar" sorusunu cevapliyor ama
    "bu enzim hangi dalda oturuyor" sorusunu cevaplamiyor. Agac onu gosterir.
    Kalabalik dugumlerde cocuk sayisi kirpilir ve kirpilan miktar ayrica
    yazilir, boylece eksik gosterim gizli kalmaz.
    """
    counts = Counter()
    for (taxonomy,) in con.execute("""
            SELECT p.taxonomy FROM ro r JOIN replicon p USING(nucleotide_id)
            WHERE r.is_confirmed = 1"""):
        parts = [t.strip() for t in (taxonomy or "").split(";") if t.strip()][:max_depth]
        if not parts:
            parts = ["unassigned"]
        for depth in range(1, len(parts) + 1):
            counts[tuple(parts[:depth])] += 1

    children = defaultdict(list)
    for path in counts:
        children[path[:-1]].append(path)

    def build(path, depth=0):
        kids = sorted(children.get(path, []), key=lambda p: -counts[p])
        kept = [k for k in kids if counts[k] >= min_count or depth < 3][:max_children]
        hidden = sum(counts[k] for k in kids if k not in kept)
        return {
            "name": path[-1] if path else "all entries",
            "n": counts[path] if path else sum(counts[p] for p in children[()]),
            "depth": depth,
            "children": [build(k, depth + 1) for k in kept],
            "hidden_taxa": len(kids) - len(kept),
            "hidden_entries": hidden,
        }

    root = build(())
    root["name"] = "all confirmed entries"
    return root


# ---------------------------------------------------- provenance of downloads
PROVENANCE = [
    ("all_confirmed.fasta", "route", "Amino acid sequence of every confirmed alpha subunit.",
     "extract_genomic_context.py, build_operons.py"),
    ("all_confirmed.csv", "route", "One row per confirmed entry: coordinates, type, variant, "
     "evidence level and operon summary.", "the whole pipeline"),
    ("operons.csv", "route", "Every gene of every predicted operon, with its component call.",
     "build_operons.py"),
    ("novel_candidates.fasta", "route", "Sequences with no close characterised relative.",
     "evidence_tiers.py"),
    ("tree_all.nwk", "file", "Maximum-likelihood tree of references and variant representatives.",
     "build_phylogeny.py"),
    ("reference_pairs.csv", "file", "Pairwise identity between the 71 curated references, with "
     "their substrates; the calibration for every identity threshold.", "evidence_tiers.py"),
    ("evidence_by_cluster.csv", "file", "Evidence level composition of each type.",
     "evidence_tiers.py"),
    ("etc_by_cluster.csv", "file", "Electron-transport configurations per type.", "etc_types.py"),
    ("regulation_by_cluster.csv", "file", "Upstream architecture and promoter region per type.",
     "analyze_regulation.py"),
    ("operon_validation.json", "file", "Position and strand tests of the operon rule, with the "
     "negative control and threshold sensitivity.", "operon_validation.py"),
    ("variant_signatures.csv", "file", "Residue signature of every variant at the columns that "
     "separate variants.", "variant_signature.py"),
    ("motif_stats.json", "file", "Conservation of the eight defining columns.", "motif_stats.py"),
    ("redundancy.json", "file", "Sequence redundancy and its effect on every count.",
     "redundancy.py"),
    ("stats.json", "file", "All hypothesis tests with effect sizes at three levels.",
     "stats_overview.py"),
    ("threshold_sensitivity.json", "file", "Every headline figure recomputed across eight "
     "inclusion thresholds, so the reader can see what the choice of threshold costs.",
     "threshold_sensitivity.py"),
    ("reference_redundancy.json", "file", "Where the curated reference set repeats itself, and "
     "how many members sit under a reference that is not uniquely defined.",
     "reference_redundancy.py"),
    ("stratified_stats.json", "file", "The bridging-carboxylate association recomputed within "
     "each evidence level and assignment class, so a reader can see how much of it rests on "
     "members whose chemistry is actually known.", "stratified_stats.py"),
    ("carboxylate_rows.json", "file", "One compact row per entry (group, type, residue, evidence "
     "level, assignment class) so the association can be recomputed in the browser.",
     "stratified_stats.py"),
    ("active_site.json", "file", "Catalytic-iron pocket of every type with a structure, from "
     "crystals where they exist and from predicted models elsewhere, with the validation of "
     "that transfer and how each pocket column varies across members.", "active_site.py"),
    ("operon_relations.json", "file", "Regulator, operon architecture and transposon "
     "associations, with the annotation bias measured.", "operon_relations.py"),
    ("learning.json", "file", "What can and cannot be predicted from sequence and from genomic "
     "context, with grouped cross-validation and the leaked figure shown alongside.",
     "learn_from_data.py"),
    ("ecological_origin.json", "file", "Habitat of each type's close and distant members, the "
     "dominant habitat of its main variant, and the wastewater picture, with a genus control.",
     "ecological_origin.py"),
    # Bu uc dosya indirme rotasindan SUNULUYORDU ama bu tabloda yoktu, yani
    # indirilebilir olup belgelenmemislerdi. Izin listesi artik bu tablodan
    # turetildigi icin eksiklik sessiz kalamaz.
    ("cooccurrence.json", "file", "Permutation test of which enzyme types share a replicon "
     "more often than their abundance explains.", "cooccurrence.py"),
    ("substrate_predictability.json", "file", "How well sequence identity predicts a shared "
     "substrate label, measured on the curated references.", "substrate_predictability.py"),
    ("habitat.json", "file", "Isolation source, host kingdom and geography per replicon, "
     "with the habitat vocabulary and every keyword that decided an assignment.",
     "isolation_source.py"),
    ("ssn_edges.csv", "file", "Similarity network edges above 30 % identity.",
     "build_phylogeny.py"),
    ("ssn_nodes.csv", "file", "Similarity network nodes.", "build_phylogeny.py"),
    ("cluster_identity_matrix.csv", "file", "Highest identity between representatives of each "
     "pair of types.", "build_phylogeny.py"),
    ("leaf_profiles.csv", "file", "Variant profiles: genus composition, plasmid rate, "
     "distinguishing neighbours.", "characterize_leaves.py"),
    ("cluster_ecology_stats.csv", "file", "Ecology and mobility per type.", "analyze_ecology.py"),
    ("null_model.csv", "file", "Neighbourhood enrichment against random windows.",
     "null_model.py"),
    ("sdp_positions.csv", "file", "Specificity-determining positions at type level.",
     "analyze_variants.py"),
]


# Yayinlanan analiz dosyalari TEK yerden turetilir. Daha once ayni liste
# hem `app.py` indirme rotasinda hem `freeze.py` icinde ayri ayri duruyordu;
# yeni bir dosya eklenince biri guncellenip oteki unutuluyordu ve dosya
# indirme tablosunda GORUNUP indirilemiyordu. PROVENANCE zaten her dosyayi
# ureticisiyle birlikte sayiyor, dolayisiyla dogru kaynak odur.
ANALYSIS_FILES = frozenset(name for name, kind, _desc, _script in PROVENANCE
                           if kind == "file")


def download_table(analysis_dir):
    """Her indirilebilir dosya: ne oldugu, hangi scriptin urettigi, kac satir."""
    rows = []
    for name, kind, description, script in PROVENANCE:
        entry = {"name": name, "kind": kind, "description": description, "script": script,
                 "size": None, "rows": None, "exists": kind == "route"}
        if kind == "file":
            path = os.path.join(analysis_dir, name)
            if os.path.exists(path):
                entry["exists"] = True
                entry["size"] = os.path.getsize(path)
                if name.endswith(".csv"):
                    with open(path) as fh:
                        entry["rows"] = max(0, sum(1 for _ in fh) - 1)
        rows.append(entry)
    return rows


# ------------------------------------------------------------------ habitat
KINGDOM_LABEL = {
    "human": "human", "animal": "animal", "plant": "plant",
    "fungus": "fungus", "alga": "alga",
    "ambiguous": "name is ambiguous", "not_a_host": "not an organism",
    "unresolved": "not yet mapped",
}


def host_kingdom_rows(habitat):
    """Konak alemi tablosu: kac replikon, bunlarin kaci habitatsiz.

    Son sutun asil gerekce: habitat sozlugunun sessiz kaldigi yerde konak
    alaninin konustugu kayit sayisi.
    """
    block = (habitat or {}).get("host_kingdom") or {}
    counts = block.get("counts") or {}
    crosstab = (habitat or {}).get("host_kingdom_by_habitat") or {}
    rows = []
    for kingdom, n in counts.items():
        if kingdom in ("<no host>", "none"):
            continue
        per_habitat = crosstab.get(kingdom, {})
        silent = per_habitat.get("unknown", 0) + per_habitat.get("other", 0)
        rows.append({
            "key": kingdom,
            "label": KINGDOM_LABEL.get(kingdom, kingdom.replace("_", " ")),
            "replicons": n,
            "habitat_silent": silent,
            "silent_share": (silent / n) if n else 0.0,
            "informative": kingdom not in ("not_a_host", "unresolved", "ambiguous"),
        })
    return sorted(rows, key=lambda r: (r["informative"] is False, -r["replicons"]))


HABITAT_LABEL = {
    "soil": "soil", "rhizosphere_plant": "plant and rhizosphere",
    "rhizosphere_soil": "rhizosphere soil",
    "plant_tissue": "plant tissue and endophytes",
    "plant_associated": "plant associated, position unstated",
    "built_environment": "built environment",
    "freshwater": "freshwater", "marine": "marine", "sediment": "sediment",
    "marine_sediment": "marine sediment", "freshwater_sediment": "freshwater sediment",
    "wastewater_sludge": "wastewater and sludge",
    "contaminated_industrial": "contaminated or industrial site",
    "mining_acid_drainage": "mine and acid drainage",
    "human_clinical": "human clinical", "animal_host": "animal host",
    "gut_faecal": "gut and faecal", "food_fermented": "food and fermentation",
    "air_dust": "air and dust", "extreme_thermal": "hot spring and thermal",
    "laboratory_strain": "laboratory strain",
    "water_unspecified": "water, compartment unstated",
    "engineered_system": "engineered system", "cave": "cave",
    "subsurface_deep": "deep subsurface", "fungal_associated": "fungus associated",
    "hypersaline": "hypersaline", "other": "not classifiable", "unknown": "no source recorded",
}


def habitat_rows(habitat, min_species=1):
    """Habitat dagilimi: giris ve TEKIL TUR sayisi yan yana.

    Ikisini birlikte vermek bu sayfanin butun noktasi: giris sayisi cok
    dizilenmis suslari sayar, tur sayisi saymaz. Aradaki fark yanliligin
    kendisidir, bu yuzden gizlenmez.
    """
    if not habitat:
        return []
    total_e = sum(v["entries"] for v in habitat["habitats"].values())
    total_s = sum(v["species"] for v in habitat["habitats"].values())
    rows = []
    for key, v in habitat["habitats"].items():
        if v["species"] < min_species:
            continue
        rows.append({
            "key": key, "label": HABITAT_LABEL.get(key, key.replace("_", " ")),
            "entries": v["entries"], "replicons": v["replicons"], "species": v["species"],
            "entry_share": v["entries"] / total_e if total_e else 0,
            "species_share": v["species"] / total_s if total_s else 0,
            "bias": (v["entries"] / total_e - v["species"] / total_s)
                    if total_e and total_s else 0,
            "entries_per_species": v["entries"] / v["species"] if v["species"] else None,
            "informative": key not in ("unknown", "other"),
        })
    # Gercek habitatlar once, "kaynak yok" ve "siniflanamaz" en sona: ikisi
    # habitat degil ve tur sayisi en yuksek olanlar oldugu icin tabloyu
    # tepeden isgal ediyorlardi.
    return sorted(rows, key=lambda r: (r["informative"] is False, -r["species"]))


def habitat_enrichment(habitat, profile_key, habitat_key, min_pairs=10):
    """Bir habitatin TUR BAZINDA zenginlesmesi, arka plana gore kat olarak.

    Arka plan tum tur-habitat ciftleri uzerinden hesaplanir. Giris sayisi degil
    tur sayisi kullanilir, cunku 400 dizilenmis Pseudomonas susu tek bir
    ekolojik gozlemdir.
    """
    if not habitat:
        return None
    prof = habitat.get(profile_key) or {}
    all_h = habitat.get("habitats") or {}
    bg_species = (all_h.get(habitat_key) or {}).get("species", 0)
    bg_total = sum(v.get("species", 0) for v in all_h.values())
    if not bg_total or not bg_species:
        return None
    background = bg_species / bg_total
    rows = []
    for key, entry in prof.items():
        pairs = entry.get("species_habitat_pairs") or 0
        if pairs < min_pairs:
            continue
        h = (entry.get("habitats") or {}).get(habitat_key) or {}
        share = h.get("species_fraction", 0.0)
        rows.append({
            "key": key, "label": key.replace("_", " "),
            "species": h.get("species", 0), "pairs": pairs,
            "entries_total": entry.get("entries_total", 0),
            "species_share": share, "entry_share": h.get("entry_fraction", 0.0),
            "fold": share / background if background else None,
        })
    return {"background": background, "background_species": bg_species,
            "background_pairs": bg_total,
            "rows": sorted(rows, key=lambda r: -(r["fold"] or 0))}
