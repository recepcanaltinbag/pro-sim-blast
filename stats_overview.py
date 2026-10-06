"""
Istatistiksel cikarimlar -- veritabaninin "Statistics" sayfasinin kaynagi.

Her test icin: soru, test, n, etki buyuklugu, p, ve BAGIMSIZLIK uyarisi.
RO girisleri bagimsiz gozlem degildir (ayni tur/sus genomlari tekrar eder);
bu yuzden her test iki duzeyde yapilir:
    entry  : tum dogrulanmis RO'lar (n=11k)  -- p degerleri iyimser
    genus  : kume x cins basina tek gozlem (tekrar eden suslar cokertilir) -- muhafazakar
Sonuc ikisinde de ayni yonde ve genus duzeyinde p<0.05 ise "saglam" sayilir.

Cikti: analysis_out/stats.json  (web sayfasi dogrudan okur)
"""

import argparse
import csv
import json
import os
import sqlite3
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

SUBSTRATE_CLASSES = ["xenobiotic", "natural_aromatic", "natural_specialized"]


def trim_table(rows, cols, table):
    """Sifir toplamli satir/sutunlari at (ki-kare beklenen sifir hatasi)."""
    keep_c = [j for j in range(len(cols)) if sum(row[j] for row in table) > 0]
    keep_r = [i for i in range(len(rows)) if sum(table[i][j] for j in keep_c) > 0]
    return ([rows[i] for i in keep_r], [cols[j] for j in keep_c],
            [[table[i][j] for j in keep_c] for i in keep_r])


def cramers_v(table):
    table = np.asarray(table, dtype=float)
    chi2 = stats.chi2_contingency(table, correction=False)[0]
    n = table.sum()
    k = min(table.shape) - 1
    return float(np.sqrt(chi2 / (n * k))) if n and k else 0.0


def two_by_two(a_yes, a_no, b_yes, b_no):
    table = [[a_yes, a_no], [b_yes, b_no]]
    odds, p = stats.fisher_exact(table)
    ra = a_yes / max(1, a_yes + a_no)
    rb = b_yes / max(1, b_yes + b_no)
    return {"table": table, "rate_a": ra, "rate_b": rb,
            "ratio": (ra / rb) if rb else None, "odds_ratio": odds, "p": p}


def load(con, ecology_path):
    eco = {}
    with open(ecology_path) as fh:
        for r in csv.DictReader(fh):
            eco[r["cluster"]] = r
    rows = con.execute("""
        SELECT r.candidate_id, r.sequence, r.ro_cluster cluster, r.ro_group grp, p.organism,
               p.taxonomy, p.is_plasmid, p.cds_count,
               o.has_beta, o.has_ferredoxin, o.has_reductase, o.completeness,
               e.tier, e.ref_identity, t.reductase_type, t.ferredoxin_type, d.domain,
               s.host_kingdom,
               (SELECT COUNT(*) FROM neighbor nb JOIN gene_category c ON c.neighbor_id=nb.neighbor_id
                 AND c.method='regex_v1' WHERE nb.candidate_id=r.candidate_id AND c.category='transposon') transposons
        FROM ro r JOIN replicon p USING(nucleotide_id)
        LEFT JOIN operon o ON o.candidate_id=r.candidate_id
        LEFT JOIN ro_evidence e ON e.candidate_id=r.candidate_id
        LEFT JOIN ro_etc t ON t.candidate_id=r.candidate_id
        LEFT JOIN ro_domain d ON d.candidate_id=r.candidate_id
        LEFT JOIN replicon_source s ON s.nucleotide_id=r.nucleotide_id
        WHERE r.is_confirmed=1""").fetchall()
    cols = ["candidate_id", "sequence", "cluster", "group", "organism", "taxonomy",
            "is_plasmid", "cds_count",
            "has_beta", "has_ferredoxin", "has_reductase", "completeness", "tier", "ref_identity",
            "reductase_type", "ferredoxin_type", "domain", "host_kingdom", "transposons"]
    data = []
    for r in rows:
        d = dict(zip(cols, r))
        e = eco.get(d["cluster"], {})
        d["sclass"] = e.get("substrate_class", "unknown")
        d["confidence"] = e.get("confidence", "")
        tax = (d["taxonomy"] or "").split("; ")
        d["phylum"] = tax[2] if len(tax) > 2 else (tax[-1] if tax else "?")
        d["genus"] = (d["organism"] or "?").split()[0]
        d["mobile"] = int(bool(d["is_plasmid"]) or (d["transposons"] or 0) > 0)
        data.append(d)
    return data


def collapse_sequence(data):
    """Ayni amino asit dizisi basina tek gozlem.

    Genom veritabanlarinda ayni protein yuzlerce susta tekrarlanir; her kopyayi
    bagimsiz gozlem saymak etki buyuklugunu degil guven araligini sisirir.
    Ozellik degerleri kopyalar arasinda cogunluk oyuyla belirlenir.
    """
    groups = defaultdict(list)
    for d in data:
        groups[d.get("sequence") or d["candidate_id"]].append(d)
    out = []
    for _, items in groups.items():
        base = dict(items[0])
        for k in ("is_plasmid", "mobile", "has_beta", "has_ferredoxin", "has_reductase"):
            base[k] = int(np.mean([bool(i[k]) for i in items]) >= 0.5)
        base["n_entries"] = len(items)
        out.append(base)
    return out


def collapse_genus(data):
    """kume x cins basina tek temsil: oranlar cins icinde ortalanir, sonra yuvarlanir."""
    groups = defaultdict(list)
    for d in data:
        groups[(d["cluster"], d["genus"])].append(d)
    out = []
    for (cl, g), items in groups.items():
        base = dict(items[0])
        for k in ("is_plasmid", "mobile", "has_beta", "has_ferredoxin", "has_reductase"):
            base[k] = int(np.mean([bool(i[k]) for i in items]) >= 0.5)
        base["n_entries"] = len(items)
        out.append(base)
    return out


def test_binary_by_class(data, field, label, exclude_euk_heavy):
    """substrat sinifi (xenobiotic vs natural) x ikili ozellik, iki duzeyde."""
    res = {"question": label, "levels": {}}
    for level, rows in (("entry", data), ("sequence", collapse_sequence(data)),
                        ("genus", collapse_genus(data))):
        rows = [r for r in rows if r["sclass"] in SUBSTRATE_CLASSES
                and r["cluster"] not in exclude_euk_heavy and r["domain"] == "Bacteria"]
        x = [r for r in rows if r["sclass"] == "xenobiotic"]
        n = [r for r in rows if r["sclass"] != "xenobiotic"]
        res["levels"][level] = two_by_two(
            sum(bool(r[field]) for r in x), sum(not r[field] for r in x),
            sum(bool(r[field]) for r in n), sum(not r[field] for r in n))
        res["levels"][level]["n"] = len(rows)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--ecology", default="cluster_ecology.csv")
    ap.add_argument("--domain-csv", default="analysis_out/domain_by_cluster.csv")
    ap.add_argument("--chemistry", default="chemistry.csv")
    ap.add_argument("--out", default="analysis_out/stats.json")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    data = load(con, args.ecology)
    # okaryot-agirlikli kumeler bakteriyel mobilite testinden disarida
    euk_heavy = set()
    if os.path.exists(args.domain_csv):
        with open(args.domain_csv) as fh:
            for r in csv.DictReader(fh):
                try:
                    if float(r.get("eukaryota_rate", 0) or 0) >= 0.2:
                        euk_heavy.add(r["cluster"])
                except ValueError:
                    pass
    out = {"n_entries": len(data), "excluded_euk_heavy_clusters": sorted(euk_heavy), "tests": []}

    # 1. ksenobiyotik x plazmit / mobil element
    out["tests"].append({"id": "plasmid_by_class",
                         **test_binary_by_class(data, "is_plasmid",
                                                "Are xenobiotic-degrading RO types more often plasmid-borne than natural-substrate types?", euk_heavy)})
    out["tests"].append({"id": "mobile_by_class",
                         **test_binary_by_class(data, "mobile",
                                                "Are xenobiotic RO types more often associated with mobile elements (plasmid or transposase within 10 kb)?", euk_heavy)})

    # 2. taksonomik yayilim: kume basina cins/filum sayisi, sinifa gore (Mann-Whitney)
    per_cluster = defaultdict(lambda: {"genera": set(), "phyla": set(), "n": 0})
    cls_of = {}
    for d in data:
        pc = per_cluster[d["cluster"]]
        pc["genera"].add(d["genus"]); pc["phyla"].add(d["phylum"]); pc["n"] += 1
        cls_of[d["cluster"]] = d["sclass"]
    xen = [len(v["genera"]) / v["n"] for c, v in per_cluster.items() if cls_of[c] == "xenobiotic" and v["n"] >= 10]
    nat = [len(v["genera"]) / v["n"] for c, v in per_cluster.items()
           if cls_of[c] in ("natural_aromatic", "natural_specialized") and v["n"] >= 10]
    u, p = stats.mannwhitneyu(xen, nat, alternative="two-sided") if xen and nat else (None, None)
    out["tests"].append({
        "id": "taxonomic_breadth",
        "question": "Do natural-substrate RO types span more genera per member (older, vertically inherited) than xenobiotic types?",
        "unit": "cluster (n>=10 members)", "metric": "genera per member",
        "xenobiotic": {"n": len(xen), "median": float(np.median(xen)) if xen else None},
        "natural": {"n": len(nat), "median": float(np.median(nat)) if nat else None},
        "mannwhitney_U": float(u) if u is not None else None, "p": float(p) if p is not None else None})

    # 3. ETC profili x RO grubu (ki-kare + Cramér's V); beta x grup
    groups = sorted({d["group"] for d in data})
    red_types = ["FNR", "GR", "FNR+GR", "none"]
    table = [[sum(1 for d in data if d["group"] == g and d["reductase_type"] == t) for t in red_types] for g in groups]
    groups_r, red_types, table = trim_table(groups, red_types, table)
    chi2, p, dof, _ = stats.chi2_contingency(table)
    out["tests"].append({
        "id": "reductase_by_group",
        "question": "Is the reductase type found near the alpha subunit associated with the RO group?",
        "rows": groups_r, "cols": red_types, "table": table,
        "chi2": float(chi2), "dof": int(dof), "p": float(p), "cramers_v": cramers_v(table)})
    fd_types = ["rieske", "plant", "rieske+plant", "none"]
    table = [[sum(1 for d in data if d["group"] == g and d["ferredoxin_type"] == t) for t in fd_types] for g in groups]
    groups_f, fd_types, table = trim_table(groups, fd_types, table)
    chi2, p, dof, _ = stats.chi2_contingency(table)
    out["tests"].append({
        "id": "ferredoxin_by_group", "question": "Is the ferredoxin type associated with the RO group?",
        "rows": groups_f, "cols": fd_types, "table": table,
        "chi2": float(chi2), "dof": int(dof), "p": float(p), "cramers_v": cramers_v(table)})
    # Konak alemi x RO grubu. Yeni bir boyut: /host alani habitat'tan AYRI
    # tutuluyor ve bu, "hangi enzim hangi canliyla yasayan bakteride bulunuyor"
    # sorusunu sorulabilir kiliyor. Yalnizca alemi COZULMUS kayitlar girer;
    # "not_a_host" ve "ambiguous" disarida kalir cunku ekoloji tasimiyorlar.
    # Hem giris hem CINS duzeyinde olculur, cunku orneklem yanliligi burada
    # asiri: tek bir Arabidopsis ya da klinik projesi yuzlerce giris uretiyor.
    HOST_KINGDOMS = ["human", "animal", "plant"]
    for level, rows in (("entry", data), ("genus", collapse_genus(data))):
        hosted = [d for d in rows if d.get("host_kingdom") in HOST_KINGDOMS]
        if len(hosted) < 30:
            continue
        keys = sorted({d["group"] for d in hosted})
        table = [[sum(1 for d in hosted if d["group"] == k and d["host_kingdom"] == c)
                  for c in HOST_KINGDOMS] for k in keys]
        keys_t, cols_t, table_t = trim_table(keys, HOST_KINGDOMS, table)
        if len(keys_t) < 2 or len(cols_t) < 2:
            continue
        chi2, p, dof, _ = stats.chi2_contingency(table_t)
        out["tests"].append({
            "id": f"host_kingdom_by_group_{level}",
            "question": ("Do the RO groups differ in the kind of organism their carrier "
                         "was isolated from?" if level == "entry" else
                         "Does that difference survive collapsing strains to one "
                         "observation per type and genus?"),
            "rows": keys_t, "cols": cols_t, "table": table_t,
            "chi2": float(chi2), "dof": int(dof), "p": float(p),
            "cramers_v": cramers_v(table_t),
            "unit": "alpha subunit" if level == "entry" else "type and genus",
            "row_header": "RO group", "row_prefix": "group ",
            "n": sum(sum(r) for r in table_t)})

    # Ayni soru substrat SINIFI icin: ksenobiyotik kimya bitkiyle mi insanla mi
    # yasayan bakterilerde yogunlasiyor?
    for level, rows in (("entry", data), ("genus", collapse_genus(data))):
        hosted = [d for d in rows if d.get("host_kingdom") in HOST_KINGDOMS
                  and d["sclass"] in SUBSTRATE_CLASSES]
        if len(hosted) < 30:
            continue
        keys = sorted({d["sclass"] for d in hosted})
        table = [[sum(1 for d in hosted if d["sclass"] == k and d["host_kingdom"] == c)
                  for c in HOST_KINGDOMS] for k in keys]
        keys_t, cols_t, table_t = trim_table(keys, HOST_KINGDOMS, table)
        if len(keys_t) < 2 or len(cols_t) < 2:
            continue
        chi2, p, dof, _ = stats.chi2_contingency(table_t)
        out["tests"].append({
            "id": f"host_kingdom_by_class_{level}",
            "question": ("Is the substrate class of the enzyme associated with the kind of "
                         "organism its carrier was isolated from?" if level == "entry" else
                         "Does that association survive collapsing strains to one "
                         "observation per type and genus?"),
            "rows": keys_t, "cols": cols_t, "table": table_t,
            "chi2": float(chi2), "dof": int(dof), "p": float(p),
            "cramers_v": cramers_v(table_t),
            "unit": "alpha subunit" if level == "entry" else "type and genus",
            "row_header": "Substrate class", "row_prefix": "",
            "n": sum(sum(r) for r in table_t)})

    table = [[sum(1 for d in data if d["group"] == g and d["has_beta"] == b) for b in (1, 0)] for g in groups]
    groups_b, beta_cols, table = trim_table(groups, ["beta", "no beta"], table)
    chi2, p, dof, _ = stats.chi2_contingency(table)
    out["tests"].append({
        "id": "beta_by_group", "question": "Is a beta subunit in the operon (α3β3 architecture) group-specific?",
        "rows": groups_b, "cols": beta_cols, "table": table,
        "chi2": float(chi2), "dof": int(dof), "p": float(p), "cramers_v": cramers_v(table)})

    # 4. kanit duzeyi: substrat etiketi kac uye icin savunulabilir
    tiers = ["characterized", "close_homolog", "family_member", "distant", "novel"]
    tier_counts = Counter(d["tier"] for d in data)
    by_cluster = defaultdict(Counter)
    for d in data:
        by_cluster[d["cluster"]][d["tier"]] += 1
    out["evidence"] = {
        "tiers": tiers, "counts": [tier_counts[t] for t in tiers],
        "fraction_substrate_transferable": (tier_counts["characterized"] + tier_counts["close_homolog"]) / len(data),
        "clusters_mostly_distant": sorted(
            [c for c, cnt in by_cluster.items()
             if (cnt["distant"] + cnt["novel"]) / sum(cnt.values()) >= 0.8 and sum(cnt.values()) >= 20]),
    }

    # 5. operon tamligi x plazmit (mobil operonlar daha eksiksiz mi?)
    bact = [d for d in data if d["domain"] == "Bacteria" and d["cds_count"] and d["cds_count"] >= 20]
    pl = [d for d in bact if d["is_plasmid"]]
    ch = [d for d in bact if not d["is_plasmid"]]
    out["tests"].append({
        "id": "operon_completeness_by_plasmid",
        "question": "Do plasmid-borne alpha subunits carry more complete operons (β/Fd/reductase) than chromosomal ones? (replicons with >=20 CDS)",
        "plasmid": {"n": len(pl), "mean_components": float(np.mean([d["completeness"] or 0 for d in pl])) if pl else None},
        "chromosome": {"n": len(ch), "mean_components": float(np.mean([d["completeness"] or 0 for d in ch])) if ch else None},
        "mannwhitney_p": float(stats.mannwhitneyu([d["completeness"] or 0 for d in pl],
                                                  [d["completeness"] or 0 for d in ch]).pvalue) if pl and ch else None})

    # 5b. Kopru karboksilati Asp mi Glu mu -- grup ve reaksiyon sinifiyla iliskisi
    #
    # Olculen sey: iki karboksilatin DAVRANISI farkli. Katalitik demir
    # karboksilati neredeyse istisnasiz Asp (476:1), yani bu pozisyon ailenin
    # degismeyen parcasi. Kopru karboksilati ise 2,5:1 bolunuyor ve Glu dagilimi
    # rastgele degil: kuaterner amin dalini isaretliyor. Test bu izlenimi
    # olcuye donusturur.
    carbox = con.execute("""
        SELECT c.bridging_residue, c.catalytic_residue, r.ro_group, r.ro_cluster
        FROM ro_carboxylate c JOIN ro r USING(candidate_id)
        WHERE r.is_confirmed = 1""").fetchall() \
        if con.execute("SELECT name FROM sqlite_master WHERE name='ro_carboxylate'").fetchone() \
        else []
    if carbox:
        chem_class = {}
        with open(args.chemistry) as fh:
            for row in csv.DictReader(fh):
                chem_class[row["cluster"]] = row.get("reaction_class", "unknown")
        for test_id, key_fn, label in (
                ("bridging_residue_by_group", lambda r: r[2], "RO group"),
                ("bridging_residue_by_reaction",
                 lambda r: chem_class.get(r[3], "unknown"), "reaction class")):
            carbox_only = [r for r in carbox if r[0] in ("D", "E")]
            keys = sorted({key_fn(r) for r in carbox_only})
            cols = ["Asp", "Glu"]

            # YALNIZCA karboksilat tasiyan girisler. "other" kategorisi hizalama
            # boslugu demek ve orani gruplar arasinda cok degisiyor (grup 1'de
            # %27, grup 2'de %1); onu tabloya katmak kalinti kimligi testini
            # kaplama testine cevirir. Olculdu: karisik tabloda katalitik
            # kontrol V=0,20 cikiyordu, oysa Asp/Glu'ya daraltinca sinyal yok.
            table = [[sum(1 for r in carbox_only if key_fn(r) == k
                          and ("Asp" if r[0] == "D" else "Glu") == c)
                      for c in cols] for k in keys]
            keys_t, cols_t, table = trim_table(keys, cols, table)
            chi2, p, dof, _ = stats.chi2_contingency(table)
            out["tests"].append({
                "id": test_id,
                "question": f"Among entries that have a bridging carboxylate, is its identity "
                            f"(Asp or Glu) associated with the {label}?",
                "rows": keys_t, "cols": cols_t, "table": table,
                "chi2": float(chi2), "dof": int(dof), "p": float(p),
                "cramers_v": cramers_v(table)})
        # Katalitik karboksilat: KONTROL testi. Yine yalnizca karboksilat tasiyan
        # girisler; bu pozisyon degismezse burada sinyal cikmamali.
        cat_only = [r for r in carbox if r[1] in ("D", "E")]
        keys = sorted({r[2] for r in cat_only})
        cols = ["Asp", "Glu"]
        table = [[sum(1 for r in cat_only if r[2] == k
                      and ("Asp" if r[1] == "D" else "Glu") == c)
                  for c in cols] for k in keys]
        keys_t, cols_t, table = trim_table(keys, cols, table)
        chi2, p, dof, _ = stats.chi2_contingency(table)
        out["tests"].append({
            "id": "catalytic_residue_by_group",
            "question": "Among entries that have a catalytic carboxylate at all, is its "
                        "identity (Asp or Glu) associated with the RO group? This is the "
                        "control for the test above: the position is nearly invariant, so "
                        "a strong association here would mean the method is picking up "
                        "alignment artefacts rather than chemistry.",
            "rows": keys_t, "cols": cols_t, "table": table,
            "chi2": float(chi2), "dof": int(dof), "p": float(p),
            "cramers_v": cramers_v(table)})

    # 6. grup x filum (dagilim tablosu + Cramér's V)
    phyla = [p for p, _ in Counter(d["phylum"] for d in data).most_common(8)]
    table = [[sum(1 for d in data if d["group"] == g and d["phylum"] == ph) for ph in phyla] for g in groups]
    groups_p, phyla, table = trim_table(groups, phyla, table)
    chi2, p, dof, _ = stats.chi2_contingency(table)
    out["tests"].append({
        "id": "phylum_by_group", "question": "How are RO groups distributed across phyla? (top 8 phyla)",
        "rows": groups_p, "cols": phyla, "table": table,
        "chi2": float(chi2), "dof": int(dof), "p": float(p), "cramers_v": cramers_v(table)})

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"[yazildi] {args.out}")
    for t in out["tests"]:
        if "levels" in t:
            print(f"  {t['id']}")
            for level, r in t["levels"].items():
                print(f"      {level:9s} {r['rate_a']:.3f} vs {r['rate_b']:.3f} "
                      f"({r['ratio']:.2f}x, p={r['p']:.1e}, n={r['n']})")
        elif "cramers_v" in t:
            print(f"  {t['id']:32s} chi2={t['chi2']:.0f} dof={t['dof']} p={t['p']:.1e} V={t['cramers_v']:.2f}")
        else:
            print(f"  {t['id']:32s} {json.dumps({k: v for k, v in t.items() if k not in ('question', 'id')}, default=float)[:160]}")
    print("  evidence:", dict(zip(tiers, out["evidence"]["counts"])),
          f"transferable={out['evidence']['fraction_substrate_transferable']:.2f}",
          "mostly_distant:", out["evidence"]["clusters_mostly_distant"])


if __name__ == "__main__":
    main()
