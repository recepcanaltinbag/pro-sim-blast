"""
Interaktif ROAR-DB gezgini uret (explorer.html).

Tek dosya, disaridan hicbir sey yuklemez: veri sayfaya JSON olarak gomulur.
71 kume icin istatistik, taksonomi, komsuluk profili ve sinteni ornekleri.

Kullanici: kume arar, gruba/substrat sinifina gore filtreler, siralar, tiklar
ve detay panelinde o kumenin genomik baglamini gorur.
"""

import argparse
import csv
import json
import os
import sqlite3
from collections import Counter, defaultdict

# Null modelden gelen arka plan oranlari (kategori basina kontrol penceresi)
NULL_PATH = "analysis_out/null_model.csv"

CATEGORY_LABELS = {
    "ro_beta": "RO beta alt birim", "ro_alpha": "RO alpha alt birim",
    "ferredoxin": "ferredoksin", "reductase": "reduktaz",
    "ring_cleavage": "halka acilimi", "regulator": "regulator",
    "transposon": "mobil element", "transporter": "tasiyici",
    "dehydrogenase": "dehidrogenaz", "hydrolase": "hidrolaz",
    "hypothetical": "hipotetik", "other": "diger",
}

CATEGORY_SHORT = {
    "ro_alpha": "alpha", "ro_beta": "beta", "ferredoxin": "fdx",
    "reductase": "red", "ring_cleavage": "ring", "regulator": "reg",
    "transposon": "Tn", "transporter": "trn", "dehydrogenase": "dh",
    "hydrolase": "hyd", "hypothetical": "hyp", "other": "-",
}


def load_ecology(path):
    if not os.path.exists(path):
        return {}
    with open(path) as handle:
        return {row["cluster"]: row for row in csv.DictReader(handle)}


def load_variance(path="analysis_out/cluster_variance.csv"):
    """Kume ici kimlik ve alt-aile sayisi."""
    if not os.path.exists(path):
        return {}
    with open(path) as handle:
        return {row["cluster"]: row for row in csv.DictReader(handle)}


def load_subfamily(db_path):
    """Kume basina alt-aile sayisi ve novellik sinifi dagilimi."""
    import sqlite3
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    out = {}
    for row in con.execute("""
        SELECT cluster,
               COUNT(*) subfams,
               SUM(CASE WHEN size>=10 THEN 1 ELSE 0 END) ref_types
        FROM subfamily GROUP BY cluster"""):
        out.setdefault(row["cluster"], {})["subfams"] = row["subfams"]
        out[row["cluster"]]["ref_types"] = row["ref_types"]
    for row in con.execute("""
        SELECT cluster, assignment_class cls, COUNT(*) n
        FROM ro_subfamily GROUP BY cluster, assignment_class"""):
        out.setdefault(row["cluster"], {}).setdefault("cls", {})[row["cls"]] = row["n"]
    con.close()
    return out


def load_domains(db_path):
    """Kume basina yasam alani kirilimi: [bakteri, okaryot, arke, baskin_ok_grup]."""
    import sqlite3
    from collections import Counter
    con = sqlite3.connect(db_path); con.row_factory = sqlite3.Row
    if not con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='ro_domain'").fetchone():
        con.close(); return {}
    out = {}
    for row in con.execute("""SELECT cluster,
            SUM(domain='Bacteria') b, SUM(domain='Eukaryota') e,
            SUM(domain='Archaea') a FROM ro_domain GROUP BY cluster"""):
        out[row["cluster"]] = [row["b"], row["e"], row["a"], ""]
    for row in con.execute("""SELECT cluster, euk_group, COUNT(*) n FROM ro_domain
            WHERE domain='Eukaryota' GROUP BY cluster, euk_group
            ORDER BY cluster, n DESC"""):
        # ORDER BY n DESC: her kume icin ILK gorulen = en kalabalik okaryot grup
        if row["cluster"] in out and (not out[row["cluster"]][3]):
            out[row["cluster"]][3] = row["euk_group"]
    con.close()
    return out


def load_representatives(path="analysis_out/representatives.csv"):
    """leaf_id -> temsilci protein_id."""
    import csv as _csv
    if not os.path.exists(path):
        return {}
    with open(path) as handle:
        return {row["leaf_id"]: row["rep_protein_id"]
                for row in _csv.DictReader(handle) if row.get("rep_protein_id")}


def load_leaf_profiles(db_path, reps):
    """Kume basina varyant (yaprak) profilleri -- imza + etiket + mobilite + temsilci."""
    import sqlite3
    con = sqlite3.connect(db_path); con.row_factory = sqlite3.Row
    out = {}
    if con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='leaf_profile'").fetchone():
        for row in con.execute("""SELECT * FROM leaf_profile WHERE size>=5
                                  ORDER BY cluster, size DESC"""):
            idx = row["leaf_id"].split("#")[-1]
            out.setdefault(row["cluster"], []).append([
                idx, row["size"], round(row["median_identity"], 2),
                row["label"] or "", row["neighbor_signature"] or "",
                round(row["plasmid_rate"], 2), round(row["transposon_rate"], 2),
                row["genus_n"], row["enriched_categories"] or "",
                reps.get(row["leaf_id"], ""),
            ])
    con.close()
    return out


def load_leaves(db_path):
    """Kume basina yaprak sayisi + uye->yaprak eslemesi (ozyinelemeli homojenizasyon)."""
    import sqlite3
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    per_cluster, member_leaf = {}, {}
    if con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='leaf'").fetchone():
        for row in con.execute("""SELECT cluster, COUNT(*) leaves,
                SUM(CASE WHEN size>=10 THEN 1 ELSE 0 END) big,
                MAX(size) biggest, MAX(depth) depth
                FROM leaf GROUP BY cluster"""):
            per_cluster[row["cluster"]] = dict(leaves=row["leaves"], big=row["big"],
                                               biggest=row["biggest"], depth=row["depth"])
        for row in con.execute("SELECT candidate_id, leaf_id, median_identity "
                               "FROM ro_leaf JOIN leaf USING(leaf_id)"):
            member_leaf[row["candidate_id"]] = (row["leaf_id"], round(row["median_identity"], 2))
    con.close()
    return per_cluster, member_leaf


def load_null(path):
    if not os.path.exists(path):
        return {}
    with open(path) as handle:
        return {row["category"]: float(row["per_null_window"] or 0)
                for row in csv.DictReader(handle)}


def synteny_string(connection, candidate_id, cluster):
    """Bir RO'nun cevresindeki gen dizilimini yon oklariyla kur."""
    neighbors = connection.execute("""
        SELECT nb.gene_offset, nb.same_strand, nb.product, c.category
        FROM neighbor nb
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        WHERE nb.candidate_id = ? AND ABS(nb.gene_offset) <= 4
        ORDER BY nb.gene_offset
    """, (candidate_id,)).fetchall()

    genes, placed = [], False
    for neighbor in neighbors:
        if neighbor["gene_offset"] > 0 and not placed:
            genes.append({"label": cluster, "dir": "fwd", "self": True,
                          "product": "RO alpha subunit"})
            placed = True
        genes.append({
            "label": CATEGORY_SHORT.get(neighbor["category"], "?"),
            "dir": "fwd" if neighbor["same_strand"] else "rev",
            "self": False,
            "cat": neighbor["category"],
            "product": neighbor["product"],
        })
    if not placed:
        genes.append({"label": cluster, "dir": "fwd", "self": True,
                      "product": "RO alpha subunit"})
    return genes


def build_payload(connection, ecology, null_rates, variance, subfam_data,
                  leaf_data, member_leaf, member_class, leaf_profiles, domains):
    clusters = {}

    # --- Temel istatistikler
    rows = connection.execute("""
        SELECT r.ro_cluster cluster, r.candidate_id, r.ro_group, rep.organism,
               rep.is_plasmid, r.model_coverage, r.protein_id,
               MAX(CASE WHEN c.category='transposon' AND ABS(nb.distance)<=5000
                        THEN 1 ELSE 0 END) has_tn
        FROM ro r
        JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        LEFT JOIN neighbor nb ON nb.candidate_id = r.candidate_id
        LEFT JOIN gene_category c ON c.neighbor_id = nb.neighbor_id
             AND c.method='regex_v1'
        WHERE r.is_confirmed=1 AND rep.status='ok'
        GROUP BY r.candidate_id
    """).fetchall()

    for row in rows:
        name = row["cluster"]
        if not name or name == "N/A":
            continue
        entry = clusters.setdefault(name, {
            "cluster": name, "group": row["ro_group"], "n": 0, "tn": 0,
            "plasmid": 0, "genera": Counter(), "coverage": [], "examples": [],
        })
        entry["n"] += 1
        entry["tn"] += row["has_tn"] or 0
        entry["plasmid"] += row["is_plasmid"] or 0
        organism = row["organism"] or ""
        if organism.split():
            entry["genera"][organism.split()[0]] += 1
        if row["model_coverage"]:
            entry["coverage"].append(row["model_coverage"])
        if len(entry["examples"]) < 40:
            entry["examples"].append((row["candidate_id"], organism))
        # TUM uyeler -- indirilebilir liste icin
        leaf_id, leaf_ident = member_leaf.get(row["candidate_id"], ("", None))
        entry.setdefault("members", []).append([
            row["candidate_id"], row["protein_id"] or "", organism,
            1 if row["is_plasmid"] else 0,
            leaf_id.split("#")[-1] if leaf_id else "",
            member_class.get(row["candidate_id"], ""),
        ])

    # --- Komsuluk profili
    for row in connection.execute("""
        SELECT r.ro_cluster cluster, c.category, COUNT(*) n,
               SUM(nb.same_strand) same
        FROM ro r
        JOIN neighbor nb ON nb.candidate_id = r.candidate_id
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id AND c.method='regex_v1'
        WHERE r.is_confirmed=1
        GROUP BY r.ro_cluster, c.category
    """):
        entry = clusters.get(row["cluster"])
        if entry is not None:
            entry.setdefault("cats", {})[row["category"]] = {
                "n": row["n"], "same": row["same"] or 0}

    # --- Regulator aileleri
    for row in connection.execute("""
        SELECT r.ro_cluster cluster, c.category family, COUNT(*) n
        FROM ro r
        JOIN neighbor nb ON nb.candidate_id = r.candidate_id
        JOIN gene_category c ON c.neighbor_id = nb.neighbor_id
             AND c.method='regulator_family_v1'
        WHERE r.is_confirmed=1
        GROUP BY r.ro_cluster, c.category
    """):
        entry = clusters.get(row["cluster"])
        if entry is not None:
            entry.setdefault("fams", {})[row["family"]] = row["n"]

    # --- Paketle
    payload = []
    for name, entry in clusters.items():
        info = ecology.get(name, {})
        total_n = entry["n"]
        cats = entry.get("cats", {})
        per_ro = {k: v["n"] / total_n for k, v in cats.items()}

        category_profile = []
        for category, rate in sorted(per_ro.items(), key=lambda x: -x[1]):
            null_rate = null_rates.get(category, 0)
            category_profile.append({
                "cat": category,
                "label": CATEGORY_LABELS.get(category, category),
                "n": cats[category]["n"],
                "rate": round(rate, 3),
                "same": round(cats[category]["same"] / cats[category]["n"], 3)
                        if cats[category]["n"] else 0,
                "enr": round(rate / null_rate, 2) if null_rate else None,
            })

        # Sinteni ornekleri -- farkli organizmalardan sec
        seen_organisms, synteny = set(), []
        for candidate_id, organism in entry["examples"]:
            genus = organism.split()[0] if organism.split() else "?"
            if genus in seen_organisms:
                continue
            seen_organisms.add(genus)
            synteny.append({"organism": organism,
                            "genes": synteny_string(connection, candidate_id, name)})
            if len(synteny) >= 4:
                break

        sf = subfam_data.get(name, {})
        cls = sf.get("cls", {})
        lf = leaf_data.get(name, {})
        var = variance.get(name, {})
        median_id = float(var["median_identity"]) if var.get("median_identity") else None
        subfam = int(var["subfamilies"]) if var.get("subfamilies") else None
        payload.append({
            "ident": round(median_id, 3) if median_id is not None else None,
            "subfam": subfam,
            "subfams": sf.get("subfams"),
            "ref_types": sf.get("ref_types"),
            "n_core": cls.get("core", 0),
            "n_divergent": cls.get("divergent", 0),
            "n_novel": cls.get("novel_candidate", 0),
            "leaves": lf.get("leaves"),
            "leaves_big": lf.get("big"),
            "leaf_biggest": lf.get("biggest"),
            "leaf_depth": lf.get("depth"),
            "members": entry.get("members", []),
            "variants": leaf_profiles.get(name, []),
            "domain": domains.get(name, [0,0,0,""]),
            "sampled": int(var["sampled"]) if var.get("sampled") else None,
            "cluster": name,
            "group": int(entry["group"]) if entry["group"] else 0,
            "gene": name.split("_")[-1],
            "n": total_n,
            "tn_rate": round(entry["tn"] / total_n, 4),
            "pl_rate": round(entry["plasmid"] / total_n, 4),
            "genera_n": len(entry["genera"]),
            "top_genera": entry["genera"].most_common(8),
            "coverage": round(sum(entry["coverage"]) / len(entry["coverage"]), 3)
                        if entry["coverage"] else 0,
            "substrate": info.get("substrate", "bilinmiyor"),
            "sclass": info.get("substrate_class", "unknown"),
            "note": info.get("ecology_note", ""),
            "conf": info.get("confidence", "low"),
            "cats": category_profile,
            "fams": sorted(entry.get("fams", {}).items(), key=lambda x: -x[1])[:6],
            "synteny": synteny,
        })

    payload.sort(key=lambda c: -c["n"])
    return payload


def totals(connection):
    def scalar(query):
        return connection.execute(query).fetchone()[0]
    return {
        "confirmed": scalar("SELECT COUNT(*) FROM ro WHERE is_confirmed=1"),
        "candidates": scalar("SELECT COUNT(*) FROM ro"),
        "replicons": scalar("SELECT COUNT(*) FROM replicon WHERE status='ok'"),
        "neighbors": scalar("SELECT COUNT(*) FROM neighbor"),
        "plasmids": scalar("SELECT COUNT(*) FROM replicon WHERE is_plasmid=1"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--ecology", default="cluster_ecology.csv")
    parser.add_argument("--null", default=NULL_PATH)
    parser.add_argument("-o", "--output", default="explorer.html")
    args = parser.parse_args()

    connection = sqlite3.connect(args.db)
    connection.row_factory = sqlite3.Row
    leaf_data, member_leaf = load_leaves(args.db)
    member_class = {}
    for row in connection.execute("SELECT candidate_id, assignment_class FROM ro_subfamily"):
        member_class[row["candidate_id"]] = row["assignment_class"]
    payload = build_payload(connection, load_ecology(args.ecology),
                            load_null(args.null), load_variance(),
                            load_subfamily(args.db), leaf_data, member_leaf, member_class,
                            load_leaf_profiles(args.db, load_representatives()),
                            load_domains(args.db))
    summary = totals(connection)
    connection.close()

    data = json.dumps({"clusters": payload, "totals": summary},
                      ensure_ascii=False, separators=(",", ":"))
    with open(args.output, "w") as handle:
        handle.write(PAGE.replace("__DATA__", data))
    print(f"[yazildi] {args.output}  ({os.path.getsize(args.output):,} bayt, "
          f"{len(payload)} kume)")


PAGE = r"""<title>ROAR-DB Gezgini &middot; Rieske oksijenaz kumeleri</title>
<style>
:root {
  --ground:#F7F8F9; --panel:#FFFFFF; --sunken:#EDEFF1;
  --ink:#1A1F24; --ink-2:#4A555E; --ink-3:#78848D; --rule:#DCE1E4;
  --accent:#A6402A; --accent-soft:#F3E3DF;
  --teal:#2F6E7A; --teal-soft:#E0EDEF;
  --ok:#3F7A4E; --warn:#96701A;
  --sans:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  --serif:ui-serif,Georgia,"Iowan Old Style",Palatino,serif;
  --mono:ui-monospace,SFMono-Regular,"SF Mono",Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){:root{
  --ground:#13171A; --panel:#1A1F23; --sunken:#222A2F;
  --ink:#E4E8EA; --ink-2:#A6B1B8; --ink-3:#78858D; --rule:#2C353B;
  --accent:#E08163; --accent-soft:#33211C;
  --teal:#6FB3BF; --teal-soft:#16292D; --ok:#74B584; --warn:#D6A648;}}
:root[data-theme="dark"]{
  --ground:#13171A; --panel:#1A1F23; --sunken:#222A2F;
  --ink:#E4E8EA; --ink-2:#A6B1B8; --ink-3:#78858D; --rule:#2C353B;
  --accent:#E08163; --accent-soft:#33211C;
  --teal:#6FB3BF; --teal-soft:#16292D; --ok:#74B584; --warn:#D6A648;}
:root[data-theme="light"]{
  --ground:#F7F8F9; --panel:#FFFFFF; --sunken:#EDEFF1;
  --ink:#1A1F24; --ink-2:#4A555E; --ink-3:#78848D; --rule:#DCE1E4;
  --accent:#A6402A; --accent-soft:#F3E3DF;
  --teal:#2F6E7A; --teal-soft:#E0EDEF; --ok:#3F7A4E; --warn:#96701A;}

body{background:var(--ground);color:var(--ink);font-family:var(--sans);
  font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased;}
.num,.mono{font-family:var(--mono);font-variant-numeric:tabular-nums;}

/* ---------- ust bar ---------- */
.top{border-bottom:1px solid var(--rule);background:var(--panel);
  padding:1.5rem clamp(1rem,3vw,2rem) 0;}
.top__inner{max-width:1500px;margin:0 auto;}
.top h1{font-family:var(--serif);font-size:1.45rem;font-weight:600;
  letter-spacing:-.015em;margin-bottom:.25rem;}
.top__sub{color:var(--ink-3);font-size:.85rem;margin-bottom:1.1rem;}
.totals{display:flex;flex-wrap:wrap;gap:1.6rem;padding-bottom:1.2rem;}
.total__n{font-family:var(--mono);font-size:1.15rem;font-weight:600;
  letter-spacing:-.01em;}
.total__l{font-size:.73rem;color:var(--ink-3);text-transform:uppercase;
  letter-spacing:.06em;}

/* ---------- kontroller ---------- */
.controls{display:flex;flex-wrap:wrap;gap:.6rem;align-items:center;
  padding:0 0 1.1rem;}
.search{flex:1 1 240px;min-width:200px;padding:.5rem .75rem;font-size:.9rem;
  font-family:var(--sans);background:var(--ground);color:var(--ink);
  border:1px solid var(--rule);border-radius:4px;}
.search:focus{outline:2px solid var(--teal);outline-offset:-1px;}
.filters{display:flex;gap:.3rem;flex-wrap:wrap;}
.fbtn{padding:.4rem .7rem;font-size:.8rem;font-family:var(--sans);
  background:var(--ground);color:var(--ink-2);border:1px solid var(--rule);
  border-radius:4px;cursor:pointer;}
.fbtn:hover{border-color:var(--ink-3);}
.fbtn[aria-pressed="true"]{background:var(--ink);color:var(--ground);
  border-color:var(--ink);}
.fbtn:focus-visible{outline:2px solid var(--teal);outline-offset:1px;}

/* ---------- yerlesim ---------- */
.layout{max-width:1500px;margin:0 auto;padding:1.4rem clamp(1rem,3vw,2rem) 4rem;
  display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1.05fr);gap:1.4rem;
  align-items:start;}
@media (max-width:960px){.layout{grid-template-columns:1fr;}}

/* ---------- liste ---------- */
.list{border:1px solid var(--rule);background:var(--panel);
  max-height:78vh;overflow-y:auto;}
table{width:100%;border-collapse:collapse;font-size:.85rem;}
thead th{position:sticky;top:0;background:var(--panel);text-align:left;
  font-size:.7rem;text-transform:uppercase;letter-spacing:.06em;
  color:var(--ink-3);padding:.6rem .7rem;border-bottom:1px solid var(--rule);
  cursor:pointer;user-select:none;white-space:nowrap;}
thead th:hover{color:var(--ink);}
thead th.num{text-align:right;}
thead th[data-dir]::after{content:attr(data-dir);margin-left:.3rem;
  color:var(--accent);}
tbody tr{cursor:pointer;border-bottom:1px solid var(--rule);}
tbody tr:hover{background:var(--sunken);}
tbody tr[aria-selected="true"]{background:var(--accent-soft);}
tbody tr[aria-selected="true"] td{color:var(--ink);font-weight:600;}
td{padding:.5rem .7rem;color:var(--ink-2);}
td.num{text-align:right;font-family:var(--mono);font-variant-numeric:tabular-nums;}
.gcell{display:inline-flex;align-items:center;justify-content:center;
  width:1.35rem;height:1.35rem;border-radius:3px;font-family:var(--mono);
  font-size:.72rem;font-weight:600;}
.g1{background:#8E6FA8;color:#fff}.g2{background:#3D7EA6;color:#fff}
.g3{background:#3F8F73;color:#fff}.g4{background:#B5843A;color:#fff}
.g5{background:#A6402A;color:#fff}.g0{background:var(--sunken);color:var(--ink-3)}

/* ---------- detay ---------- */
.detail{border:1px solid var(--rule);background:var(--panel);
  max-height:78vh;overflow-y:auto;padding:1.4rem 1.4rem 2rem;}
.detail__empty{color:var(--ink-3);font-size:.9rem;padding:3rem 0;
  text-align:center;}
.detail h2{font-family:var(--mono);font-size:1.25rem;font-weight:600;
  letter-spacing:-.01em;margin-bottom:.2rem;}
.detail__sub{color:var(--ink-2);font-size:.9rem;margin-bottom:1.1rem;}
.detail h3{font-size:.78rem;text-transform:uppercase;letter-spacing:.07em;
  color:var(--ink-3);margin:1.6rem 0 .6rem;font-weight:650;}
.dstats{display:grid;grid-template-columns:repeat(auto-fit,minmax(88px,1fr));
  gap:1px;background:var(--rule);border:1px solid var(--rule);}
.dstat{background:var(--panel);padding:.7rem .75rem;}
.dstat__n{font-family:var(--mono);font-size:1.1rem;font-weight:600;
  letter-spacing:-.01em;}
.dstat__l{font-size:.68rem;color:var(--ink-3);text-transform:uppercase;
  letter-spacing:.05em;margin-top:.1rem;}

.chip{display:inline-block;font-size:.72rem;padding:.14em .5em;border-radius:2px;
  background:var(--sunken);color:var(--ink-2);}
.chip--xeno{background:var(--accent-soft);color:var(--accent);}
.chip--nat{background:var(--teal-soft);color:var(--teal);}
.chip--lowconf{background:var(--sunken);color:var(--warn);}

/* komsuluk profili */
.cat{display:grid;grid-template-columns:6.5rem 1fr 3.4rem 3.6rem;
  align-items:center;gap:.5rem;padding:.22rem 0;font-size:.79rem;}
.cat__l{color:var(--ink-2);white-space:nowrap;overflow:hidden;
  text-overflow:ellipsis;}
.cat__track{height:7px;background:var(--sunken);position:relative;}
.cat__bar{height:100%;background:var(--teal);}
.cat__mid{position:absolute;left:50%;top:-2px;bottom:-2px;width:1px;
  background:var(--ink-3);opacity:.5;}
.cat__v{font-family:var(--mono);font-size:.75rem;text-align:right;
  color:var(--ink-3);}
.cat__e{font-family:var(--mono);font-size:.75rem;text-align:right;font-weight:600;}
.e-up{color:var(--ok)}.e-down{color:var(--accent)}.e-flat{color:var(--ink-3)}

/* sinteni */
.syn{margin-bottom:.95rem;}
.syn__org{font-size:.78rem;color:var(--ink-3);font-style:italic;
  margin-bottom:.28rem;}
.syn__row{display:flex;flex-wrap:wrap;gap:3px;}
.gene{font-family:var(--mono);font-size:.68rem;padding:.18em .4em;
  background:var(--sunken);color:var(--ink-2);white-space:nowrap;
  border-left:2px solid transparent;}
.gene--self{background:var(--accent);color:#fff;font-weight:600;}
.gene--rev{border-left:2px solid var(--ink-3);}
.gene--fwd{border-right:2px solid var(--ink-3);}
.gene[data-c="regulator"]{background:var(--teal-soft);color:var(--teal);}
.gene[data-c="transposon"]{background:var(--accent-soft);color:var(--accent);}
.gene[data-c="ro_beta"],.gene[data-c="ferredoxin"],.gene[data-c="reductase"]
  {background:var(--teal-soft);color:var(--teal);}

.stack{display:flex;height:12px;border:1px solid var(--rule);overflow:hidden;margin:.5rem 0 .35rem;}
.seg{height:100%;}.seg--core{background:var(--teal);}.seg--div{background:var(--warn);}.seg--nov{background:var(--accent);}
.key{display:inline-block;width:.7rem;height:.7rem;border-radius:1px;margin:0 .2rem 0 .6rem;vertical-align:middle;}
.key--core{background:var(--teal);}.key--div{background:var(--warn);}.key--nov{background:var(--accent);}
.novelcount{font-family:var(--mono);font-weight:600;color:var(--accent);}
.domainbar{margin:.6rem 0 .2rem;padding:.55rem .65rem;background:var(--sunken);border-left:2px solid var(--teal);}
.domainbar--flag{border-left-color:var(--warn);background:var(--warn-soft);}
.domainbar__title{font-size:.78rem;font-weight:600;color:var(--ink);margin-bottom:.4rem;}
.flag{font-weight:400;color:var(--warn);font-size:.74rem;}
.dstack{display:flex;height:11px;border:1px solid var(--rule);overflow:hidden;}
.ds{height:100%;}.ds--b{background:var(--teal);}.ds--e{background:var(--accent);}.ds--a{background:var(--ink-3);}
.dlbl{display:flex;gap:.9rem;flex-wrap:wrap;font-size:.74rem;color:var(--ink-2);margin-top:.35rem;}
.key--b{background:var(--teal);}.key--e{background:var(--accent);}.key--a{background:var(--ink-3);}
.leafline{font-size:.8rem;color:var(--ink-2);margin-top:.6rem;padding:.5rem .65rem;
  background:var(--sunken);border-left:2px solid var(--teal);}
.leafline b{font-family:var(--mono);color:var(--ink);}.muted{color:var(--ink-3);}
.dlbtn{font-family:var(--sans);font-size:.72rem;padding:.2em .6em;margin-left:.6rem;
  background:var(--teal);color:#fff;border:none;border-radius:3px;cursor:pointer;
  vertical-align:middle;font-weight:600;}
.dlbtn:hover{filter:brightness(1.08);}.dlbtn:focus-visible{outline:2px solid var(--ink);}
.mtbl-wrap{border:1px solid var(--rule);max-height:340px;overflow-y:auto;margin-top:.5rem;}
.mtbl{width:100%;border-collapse:collapse;font-size:.78rem;}
.mtbl th{position:sticky;top:0;background:var(--panel);text-align:left;
  font-size:.66rem;text-transform:uppercase;letter-spacing:.05em;color:var(--ink-3);
  padding:.4rem .55rem;border-bottom:1px solid var(--rule);}
.mtbl td{padding:.32rem .55rem;border-bottom:1px solid var(--rule);color:var(--ink-2);}
.pl{color:var(--accent);font-size:.72rem;}
.mc{font-size:.68rem;padding:.1em .4em;border-radius:2px;background:var(--sunken);}
.mc--core{background:var(--teal-soft);color:var(--teal);}
.mc--novel_candidate{background:var(--accent-soft);color:var(--accent);}
.more{padding:.5rem .6rem;font-size:.76rem;color:var(--ink-3);}
.vhint{font-size:.74rem;color:var(--ink-3);font-weight:400;margin-left:.5rem;}
.variants{display:flex;flex-direction:column;gap:.5rem;margin-top:.6rem;}
.variant{border:1px solid var(--rule);border-left:3px solid var(--teal);
  padding:.55rem .7rem;background:var(--panel);}
.variant__head{display:flex;align-items:center;gap:.55rem;flex-wrap:wrap;}
.variant__id{font-size:.82rem;font-weight:600;color:var(--ink);}
.variant__n{font-family:var(--mono);font-size:.76rem;color:var(--ink-2);}
.variant__label{font-size:.82rem;color:var(--ink-2);margin-top:.28rem;}
.variant__rep{font-size:.74rem;color:var(--ink-3);margin-top:.2rem;}
.variant__rep .mono{color:var(--teal);}
.variant__sig{margin-top:.35rem;display:flex;flex-wrap:wrap;gap:.25rem;align-items:baseline;}
.siglbl{font-size:.72rem;color:var(--ink-3);}
.sigchip{font-size:.72rem;padding:.1em .45em;background:var(--teal-soft);color:var(--teal);border-radius:2px;}
.variant__meta{margin-top:.35rem;display:flex;flex-wrap:wrap;gap:.35rem;}
.vm{font-size:.7rem;padding:.1em .45em;background:var(--sunken);color:var(--ink-3);border-radius:2px;}
.vm--pl{background:var(--accent-soft);color:var(--accent);}
.vm--tn{background:var(--accent-soft);color:var(--accent);}
.vm--enr{background:var(--teal-soft);color:var(--teal);}
.dlbtn--sm{margin-left:auto;padding:.12em .5em;font-size:.68rem;}
.hom{font-family:var(--mono);font-size:.78rem;font-weight:600;}
.hom--ok{color:var(--ok)}.hom--mid{color:var(--warn)}.hom--bad{color:var(--accent)}
.hombox{border:1px solid var(--rule);border-left-width:3px;padding:.7rem .85rem;background:var(--sunken);}
.hombox--ok{border-left-color:var(--ok)}.hombox--mid{border-left-color:var(--warn)}
.hombox--bad{border-left-color:var(--accent)}
.hombox__top{display:flex;gap:.6rem;align-items:baseline;flex-wrap:wrap;}
.hom__verdict{font-size:.8rem;color:var(--ink-2);}
.hombox__sub{font-size:.78rem;color:var(--ink-3);margin-top:.25rem;}
.hombox__warn{font-size:.78rem;color:var(--ink-2);margin-top:.45rem;
  padding-top:.45rem;border-top:1px solid var(--rule);}
.genera{display:flex;flex-wrap:wrap;gap:.3rem;}
.genus{font-size:.76rem;padding:.16em .5em;background:var(--sunken);
  border-radius:2px;color:var(--ink-2);}
.genus b{font-family:var(--mono);color:var(--ink);font-weight:600;}
.note{font-size:.83rem;color:var(--ink-2);border-left:2px solid var(--rule);
  padding-left:.75rem;margin-top:.5rem;}
.empty{color:var(--ink-3);font-size:.85rem;padding:1.5rem 0;text-align:center;}
.legend{font-size:.74rem;color:var(--ink-3);margin-top:.5rem;}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important;}}
</style>

<div class="top"><div class="top__inner">
  <h1>ROAR-DB Gezgini</h1>
  <div class="top__sub">Dogrulanmis Rieske oksijenaz alpha alt birimleri &mdash;
    genomik baglam, taksonomi ve mobilite</div>
  <div class="totals" id="totals"></div>
  <div class="controls">
    <input class="search" id="q" type="search" placeholder="Kume, gen veya substrat ara&hellip;"
      aria-label="Kume ara">
    <div class="filters" id="gfilters" role="group" aria-label="Grup filtresi"></div>
    <div class="filters" id="sfilters" role="group" aria-label="Substrat filtresi"></div>
  </div>
</div></div>

<div class="layout">
  <div class="list">
    <table>
      <thead><tr>
        <th data-k="group">G</th>
        <th data-k="cluster">kume</th>
        <th data-k="substrate">substrat</th>
        <th data-k="n" class="num">RO</th>
        <th data-k="tn_rate" class="num">Tn</th>
        <th data-k="pl_rate" class="num">plaz</th>
        <th data-k="genera_n" class="num">cins</th>
        <th data-k="ident" class="num">kimlik</th>
        <th data-k="n_novel" class="num">novel</th>
      </tr></thead>
      <tbody id="rows"></tbody>
    </table>
  </div>
  <div class="detail" id="detail"></div>
</div>

<script>
const DATA = __DATA__;
const C = DATA.clusters;
const T = DATA.totals;
const fmt = n => n.toLocaleString('tr-TR');
const pct = v => (100*v).toFixed(1).replace('.',',') + '%';

document.getElementById('totals').innerHTML = [
  [T.confirmed, 'dogrulanmis RO'], [C.length, 'kume'],
  [T.replicons, 'replikon'], [T.neighbors, 'komsu kaydi'],
  [T.plasmids, 'plazmit'],
].map(([n,l]) => `<div><div class="total__n">${fmt(n)}</div>
  <div class="total__l">${l}</div></div>`).join('');

const SCLASS = {xenobiotic:['ksenobiyotik','xeno'],
  natural_aromatic:['dogal-aromatik','nat'],
  natural_specialized:['dogal-ozel','nat'], unknown:['bilinmiyor','']};

let state = {q:'', groups:new Set(), sclasses:new Set(), sort:'n', dir:-1, sel:null};

// --- filtre butonlari
const gf = document.getElementById('gfilters');
[1,2,3,4,5].forEach(g => {
  const b = document.createElement('button');
  b.className='fbtn'; b.textContent='G'+g; b.setAttribute('aria-pressed','false');
  b.onclick = () => { state.groups.has(g) ? state.groups.delete(g) : state.groups.add(g);
    b.setAttribute('aria-pressed', state.groups.has(g)); render(); };
  gf.appendChild(b);
});
const sf = document.getElementById('sfilters');
Object.entries(SCLASS).forEach(([k,[label]]) => {
  const b = document.createElement('button');
  b.className='fbtn'; b.textContent=label; b.setAttribute('aria-pressed','false');
  b.onclick = () => { state.sclasses.has(k) ? state.sclasses.delete(k) : state.sclasses.add(k);
    b.setAttribute('aria-pressed', state.sclasses.has(k)); render(); };
  sf.appendChild(b);
});

document.getElementById('q').addEventListener('input', e => {
  state.q = e.target.value.toLowerCase(); render();
});

document.querySelectorAll('thead th').forEach(th => {
  th.onclick = () => {
    const k = th.dataset.k;
    state.dir = (state.sort === k) ? -state.dir : -1;
    state.sort = k;
    document.querySelectorAll('thead th').forEach(x => x.removeAttribute('data-dir'));
    th.setAttribute('data-dir', state.dir < 0 ? '↓' : '↑');
    render();
  };
});

function visible() {
  return C.filter(c => {
    if (state.groups.size && !state.groups.has(c.group)) return false;
    if (state.sclasses.size && !state.sclasses.has(c.sclass)) return false;
    if (state.q) {
      const hay = (c.cluster+' '+c.substrate+' '+c.gene+' '+c.note).toLowerCase();
      if (!hay.includes(state.q)) return false;
    }
    return true;
  }).sort((a,b) => {
    const k = state.sort;
    const va = a[k] == null ? -Infinity : a[k], vb = b[k] == null ? -Infinity : b[k];
    if (typeof va === 'string') return state.dir * va.localeCompare(vb);
    return state.dir * (va - vb);
  });
}

function csvq(x){ return '"' + String(x==null?'':x).replace(/"/g,'""') + '"'; }
function esc(s){ return String(s==null?'':s).replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch])); }

function render() {
  const rows = visible();
  const tb = document.getElementById('rows');
  if (!rows.length) {
    tb.innerHTML = '<tr><td colspan="9" class="empty">Eslesen kume yok</td></tr>';
    return;
  }
  tb.innerHTML = rows.map(c => `<tr data-id="${c.cluster}"
      aria-selected="${state.sel===c.cluster}">
    <td><span class="gcell g${c.group}">${c.group||'?'}</span></td>
    <td class="mono">${c.cluster}</td>
    <td>${esc(c.substrate)}</td>
    <td class="num">${fmt(c.n)}</td>
    <td class="num">${pct(c.tn_rate)}</td>
    <td class="num">${pct(c.pl_rate)}</td>
    <td class="num">${c.genera_n}</td>
    <td class="num">${c.ident!=null ? `<span class="hom hom--${homClass(c.ident)}">${pct(c.ident)}</span>` : '–'}</td>
    <td class="num">${c.n_novel!=null && c.n_novel>0 ? `<span class="novelcount">${c.n_novel}</span>` : '–'}</td>
  </tr>`).join('');
  tb.querySelectorAll('tr[data-id]').forEach(tr => {
    tr.onclick = () => { state.sel = tr.dataset.id; render(); detail(state.sel); };
  });
}

function stackBar(a,b,c){
  const t=a+b+c||1;
  return `<span class="seg seg--core" style="width:${100*a/t}%"></span>`+
         `<span class="seg seg--div" style="width:${100*b/t}%"></span>`+
         `<span class="seg seg--nov" style="width:${100*c/t}%"></span>`;
}
function homClass(v){ return v>=0.6 ? 'ok' : v>=0.35 ? 'mid' : 'bad'; }
function homLabel(v){ return v>=0.6 ? 'homojen' : v>=0.35 ? 'karisik' : 'heterojen — tek enzim degil'; }
function enrClass(e) { return e == null ? 'e-flat' : e >= 1.5 ? 'e-up' : e < 0.8 ? 'e-down' : 'e-flat'; }

function detail(id) {
  const c = C.find(x => x.cluster === id);
  const d = document.getElementById('detail');
  if (!c) { d.innerHTML = '<div class="detail__empty">Bir kume secin</div>'; return; }

  const [slabel, stone] = SCLASS[c.sclass] || ['bilinmiyor',''];
  const maxRate = Math.max(...c.cats.map(x => x.rate), 0.01);

  d.innerHTML = `
    <h2>${c.cluster}</h2>
    <div class="detail__sub">
      <strong>${esc(c.substrate)}</strong>
      <span class="chip ${stone?'chip--'+stone:''}">${slabel}</span>
      ${c.conf==='dusuk' ? '<span class="chip chip--lowconf">substrat atamasi dogrulanmali</span>' : ''}
    </div>
    ${c.note ? `<div class="note">${esc(c.note)}</div>` : ''}

    <h3>Ozet</h3>
    <div class="dstats">
      <div class="dstat"><div class="dstat__n">${fmt(c.n)}</div><div class="dstat__l">RO</div></div>
      <div class="dstat"><div class="dstat__n">${pct(c.tn_rate)}</div><div class="dstat__l">Tn yakin</div></div>
      <div class="dstat"><div class="dstat__n">${pct(c.pl_rate)}</div><div class="dstat__l">plazmit</div></div>
      <div class="dstat"><div class="dstat__n">${c.genera_n}</div><div class="dstat__l">cins</div></div>
      <div class="dstat"><div class="dstat__n">${c.coverage.toFixed(2).replace('.',',')}</div><div class="dstat__l">coverage</div></div>
    </div>

    ${c.domain && (c.domain[0]+c.domain[1]+c.domain[2])>0 ? (()=>{
      const [b,e,a]=c.domain, t=b+e+a, ef=e/t;
      return `<div class="domainbar ${ef>=0.25?'domainbar--flag':''}">
        <div class="domainbar__title">Yasam alani${ef>=0.25?` <span class="flag">okaryot-agirlikli &mdash; bakteriyel yikim DB'sine ait olmayabilir</span>`:''}</div>
        <div class="dstack">
          <span class="ds ds--b" style="width:${100*b/t}%" title="Bakteri ${b}"></span>
          <span class="ds ds--e" style="width:${100*e/t}%" title="Okaryot ${e}"></span>
          <span class="ds ds--a" style="width:${100*a/t}%" title="Arke ${a}"></span>
        </div>
        <div class="dlbl"><span><span class="key key--b"></span>bakteri ${b}</span>
          ${e>0?`<span><span class="key key--e"></span>okaryot ${e}${c.domain[3]?' ('+c.domain[3]+')':''}</span>`:''}
          ${a>0?`<span><span class="key key--a"></span>arke ${a}</span>`:''}</div>
      </div>`;})() : ''}

    ${c.ident!=null ? `<h3>Kume ici homojenlik</h3>
      <div class="hombox hombox--${homClass(c.ident)}">
        <div class="hombox__top">
          <span class="hom hom--${homClass(c.ident)}">${pct(c.ident)} medyan kimlik</span>
          <span class="hom__verdict">${homLabel(c.ident)}</span>
        </div>
        ${c.subfam!=null ? `<div class="hombox__sub">${c.sampled} orneklenen uyeden
          <b>${c.subfam}</b> alt-aile (CD-HIT %70)</div>` : ''}
        ${c.ident<0.35 ? `<div class="hombox__warn">Bu kumeye atanan substrat
          tum uyeler icin gecerli sayilmamali &mdash; atama goreli en iyi modeli secer,
          mutlak bir benzerlik tabani yoktur.</div>` : ''}
      </div>` : ''}

    ${c.subfams!=null ? `<h3>Alt-aile yapisi ve novellik</h3>
      <div class="dstats">
        <div class="dstat"><div class="dstat__n">${c.subfams}</div><div class="dstat__l">alt-aile</div></div>
        <div class="dstat"><div class="dstat__n">${c.ref_types}</div><div class="dstat__l">referans tip</div></div>
        <div class="dstat"><div class="dstat__n" style="color:var(--accent)">${c.n_novel}</div><div class="dstat__l">novel aday</div></div>
      </div>
      <div class="stack" title="cekirdek / cevre / novel">
        ${stackBar(c.n_core, c.n_divergent, c.n_novel)}
      </div>
      <div class="legend"><span class="key key--core"></span>cekirdek
        <span class="key key--div"></span>cevre
        <span class="key key--nov"></span>novel aday (hicbir yerlesik tipe uymuyor)</div>
      ${c.leaves!=null ? `<div class="leafline">Ozyinelemeli homojenizasyon:
        <b>${c.leaves}</b> yaprak, <b>${c.leaves_big}</b> tanesi &ge;10 uye,
        en buyuk yaprak <b>${fmt(c.leaf_biggest)}</b> uye, derinlik ${c.leaf_depth}.
        <span class="muted">Her yaprak tek bir homojen enzim tipi.</span></div>` : ''}` : ''}

    <h3>Komsuluk profili</h3>
    ${c.cats.map(x => `<div class="cat">
      <span class="cat__l" title="${x.label}">${x.label}</span>
      <span class="cat__track"><span class="cat__bar"
        style="width:${(100*x.rate/maxRate).toFixed(0)}%"></span></span>
      <span class="cat__v">${x.rate.toFixed(2).replace('.',',')}</span>
      <span class="cat__e ${enrClass(x.enr)}">${x.enr!=null ? x.enr.toFixed(2).replace('.',',')+'×' : '–'}</span>
    </div>`).join('')}
    <div class="legend">RO basina ortalama komsu sayisi &middot; sagdaki deger
      null modele gore zenginlesme (1&times; = arka plan)</div>

    ${c.fams.length ? `<h3>Regulator aileleri</h3>
      <div class="genera">${c.fams.map(([f,n]) =>
        `<span class="genus">${f} <b>${n}</b></span>`).join('')}</div>` : ''}

    <h3>Taksonomi</h3>
    <div class="genera">${c.top_genera.map(([g,n]) =>
      `<span class="genus"><i>${g}</i> <b>${n}</b></span>`).join('')}</div>

    ${c.synteny.length ? `<h3>Sinteni ornekleri</h3>
      ${c.synteny.map(s => `<div class="syn">
        <div class="syn__org">${esc(s.organism)}</div>
        <div class="syn__row">${s.genes.map(g =>
          `<span class="gene ${g.self?'gene--self':''} gene--${g.dir}"
            ${g.cat?`data-c="${g.cat}"`:''} title="${esc(g.product)}"
            >${g.dir==='rev'?'◀':''}${g.label}${g.dir==='fwd'?'▶':''}</span>`
          ).join('')}</div>
      </div>`).join('')}
      <div class="legend">&#9654; RO ile ayni yon &middot; &#9664; ters yon
        &middot; uzerine gelince tam urun adi</div>` : ''}

    ${c.variants && c.variants.length ? `<h3>Varyantlar &middot; ${c.variants.length} ana yaprak
      <span class="vhint">ayni kume, farkli genomik baglam</span></h3>
      <div class="legend">Her yaprak homojen bir varyant. Komsuluk imzasi farkliysa
        muhtemelen farkli substrat &mdash; ayni enzim adi altinda ayri seyler.</div>
      <div class="variants">${c.variants.slice(0,14).map(v => `
        <div class="variant">
          <div class="variant__head">
            <span class="variant__id mono">${c.cluster}#${v[0]}</span>
            <span class="variant__n">${fmt(v[1])} uye</span>
            <span class="hom hom--${homClass(v[2])}">${pct(v[2])} ic kimlik</span>
            <button class="dlbtn dlbtn--sm" onclick="dlLeaf('${c.cluster}','${v[0]}')">CSV</button>
          </div>
          <div class="variant__label">${esc(v[3])}</div>
          ${v[9] ? `<div class="variant__rep">temsilci: <span class="mono">${v[9]}</span></div>` : ''}
          ${v[4] ? `<div class="variant__sig"><span class="siglbl">komsu:</span>
            ${v[4].split(' | ').slice(0,3).map(x=>`<span class="sigchip">${esc(x.slice(0,38))}</span>`).join('')}</div>` : ''}
          <div class="variant__meta">
            ${v[5]>=0.15?`<span class="vm vm--pl">plazmit %${Math.round(v[5]*100)}</span>`:''}
            ${v[6]>=0.15?`<span class="vm vm--tn">Tn %${Math.round(v[6]*100)}</span>`:''}
            <span class="vm">${v[7]} cins</span>
            ${v[8]?`<span class="vm vm--enr">zengin: ${v[8].split(';')[0]}</span>`:''}
          </div>
        </div>`).join('')}
        ${c.variants.length>14?`<div class="more">&hellip; ve ${c.variants.length-14} varyant daha</div>`:''}
      </div>` : ''}

    ${c.members && c.members.length ? `<h3>Uyeler &middot; ${fmt(c.members.length)} protein
      <button class="dlbtn" onclick="dlMembers('${c.cluster}')">CSV indir</button></h3>
      <div class="legend">protein ID (NCBI), organizma, replikon, homojen yaprak, sinif.
        Substrat: <b>${esc(c.substrate)}</b>. CSV tum uyeleri icerir.</div>
      <div class="mtbl-wrap"><table class="mtbl">
        <thead><tr><th>protein ID</th><th>organizma</th><th>rep</th>
        <th>yaprak</th><th>sinif</th></tr></thead>
        <tbody>${c.members.slice(0,60).map(m => `<tr>
          <td class="mono">${m[1]||m[0].slice(0,18)}</td>
          <td><i>${(m[2]||'?').slice(0,30)}</i></td>
          <td>${m[3]?'<span class="pl">plazmit</span>':'krom'}</td>
          <td class="mono">${m[4]!==''?'#'+m[4]:'–'}</td>
          <td><span class="mc mc--${m[5]||'na'}">${clsShort(m[5])}</span></td>
        </tr>`).join('')}</tbody>
      </table>${c.members.length>60?`<div class="more">&hellip; ve
        ${fmt(c.members.length-60)} tane daha &mdash; tamami CSV'de</div>`:''}</div>` : ''}
  `;
}

function clsShort(x){ return x==='core'?'cekirdek':x==='divergent'?'cevre':
  x==='novel_candidate'?'novel':x||'–'; }

function dlLeaf(cluster, leafIdx){
  const c = C.find(x => x.cluster === cluster);
  if (!c || !c.members) return;
  const rows = c.members.filter(m => String(m[4]) === String(leafIdx));
  const head = ['protein_id','candidate_id','organism','replicon','leaf',
                'assignment_class','cluster','substrate'];
  const lines = [head.join(',')];
  for (const m of rows){
    lines.push([m[1], m[0], m[2]||'', m[3]?'plasmid':'chromosome',
                cluster+'#'+m[4], m[5], cluster, c.substrate||'']
               .map(csvq).join(','));
  }
  triggerDownload(cluster + '_leaf' + leafIdx + '_members.csv', lines.join('\n'));
}

function dlMembers(cluster){
  const c = C.find(x => x.cluster === cluster);
  if (!c || !c.members) return;
  const head = ['protein_id','candidate_id','organism','replicon','leaf',
                'assignment_class','cluster','substrate','substrate_class'];
  const lines = [head.join(',')];
  for (const m of c.members){
    lines.push([m[1], m[0], m[2]||'', m[3]?'plasmid':'chromosome',
                m[4]!==''?c.cluster+'#'+m[4]:'', m[5], c.cluster,
                c.substrate||'', c.sclass]
               .map(csvq).join(','));
  }
  triggerDownload(cluster + '_members.csv', lines.join('\n'));
}

function triggerDownload(filename, text){
  const blob = new Blob([text], {type:'text/csv;charset=utf-8'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

render();
detail(null);
</script>
"""


if __name__ == "__main__":
    main()
