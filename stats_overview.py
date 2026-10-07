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
import random
import re
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


def bh_adjust(pvals):
    """Benjamini-Hochberg: p listesi -> q listesi, ayni sirada.

    Bir tabloda 43 satir varsa duzeltmesiz bakildiginda birkac satirin sans
    eseri p<0,05 vermesi BEKLENIR. Asagidaki hucre ve satir testleri bir
    AILE olarak dusunulur ve yanlis kesif orani bu aile icinde kontrol edilir.
    """
    n = len(pvals)
    if not n:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, pvals[i] * n / rank)
        q[i] = running
    return q


def signal_map(rows, cols, table, genus_table=None, min_expected=5.0, alpha=0.05,
               min_lift=0.10):
    """Tablodaki sinyal NEREDE yasiyor?

    Toplu ki-kare tek bir soru sorar: "bu tabloda herhangi bir iliski var mi?"
    Bu soru iki yonde de yanlis yonlendirir. Topluca anlamsiz cikan bir tablo
    tek bir satirda gercek ve guclu bir iliski barindirabilir -- ornegin
    reduktaz tipi RO GRUBU ile iliskili olmayabilirken belli bir enzim tipinde,
    ya da belli bir reaksiyon sinifinda, iliski belirgindir. Tersi de olur:
    anlamli cikan bir tablonun butun sinyali tek bir satirdan gelebilir ve
    "gruplar farklidir" ifadesi geri kalan satirlar icin yanlis olur.

    Bu yuzden her kontenjans testi iki ek duzeyde cozulur:

    hucre duzeyi  standartlastirilmis Pearson artigi
                  z = (G - B) / sqrt(B (1 - satir payi)(1 - sutun payi))
                  buyuk orneklemde yaklasik N(0,1); beklenen sayi
                  min_expected altindaysa yaklasim gecerli olmadigi icin
                  hucre atlanir.
    satir duzeyi  o satir ile KALAN satirlarin toplami arasinda 2xC testi;
                  "bu tip/grup digerlerinden farkli mi" sorusunun cevabi.

    11 bin giriste p degeri neredeyse her hucrede kucuk cikar; bu yuzden
    istatistiksel isaret tek basina yeterli sayilmaz. Bir hucre ancak hem
    duzeltilmis q < alpha ise HEM DE o satirin sutun payi genel sutun payindan
    en az min_lift kadar (varsayilan 10 puan) sapiyorsa "maddi" sayilir.

    Iki aile ayri ayri Benjamini-Hochberg ile duzeltilir. genus_table verilirse
    (ayni satir ve sutun etiketleriyle, suslar cins basina cokertilerek kurulmus
    tablo) ayni hucre testi orada da yapilir ve q_genus olarak eklenir: ornekleme
    yanliligindan dogan hucreler bu ikinci gecisi gecemez.
    """
    t = np.asarray(table, dtype=float)
    if t.ndim != 2 or t.shape[0] < 2 or t.shape[1] < 2 or t.sum() <= 0:
        return None
    n = t.sum()
    rsum = t.sum(axis=1)
    csum = t.sum(axis=0)
    exp = np.outer(rsum, csum) / n

    def cell_z(arr):
        m = arr.sum()
        if m <= 0:
            return None
        rs, cs = arr.sum(axis=1), arr.sum(axis=0)
        e = np.outer(rs, cs) / m
        with np.errstate(divide="ignore", invalid="ignore"):
            den = np.sqrt(e * (1 - rs[:, None] / m) * (1 - cs[None, :] / m))
            z = np.where(den > 0, (arr - e) / den, 0.0)
        return e, z

    exp, z = cell_z(t)
    idx, raw = [], []
    for i in range(t.shape[0]):
        for j in range(t.shape[1]):
            if exp[i, j] < min_expected:
                continue
            idx.append((i, j))
            raw.append(float(2 * stats.norm.sf(abs(z[i, j]))))
    qs = bh_adjust(raw)

    gq = {}
    if genus_table is not None:
        g = np.asarray(genus_table, dtype=float)
        if g.shape == t.shape and g.sum() > 0:
            ge, gz = cell_z(g)
            gidx, graw = [], []
            for i in range(g.shape[0]):
                for j in range(g.shape[1]):
                    if ge[i, j] < min_expected:
                        continue
                    gidx.append((i, j))
                    graw.append(float(2 * stats.norm.sf(abs(gz[i, j]))))
            for (i, j), q in zip(gidx, bh_adjust(graw)):
                gq[(i, j)] = {"q": q, "z": float(gz[i, j])}

    cells = []
    for (i, j), p, q in zip(idx, raw, qs):
        share = float(t[i, j] / rsum[i]) if rsum[i] else 0.0
        base = float(csum[j] / n)
        entry = {"row": rows[i], "col": cols[j], "obs": int(t[i, j]),
                 "expected": round(float(exp[i, j]), 1), "z": round(float(z[i, j]), 2),
                 "p": p, "q": q, "share": round(share, 4), "baseline": round(base, 4),
                 "lift": round(share - base, 4),
                 "direction": "over" if z[i, j] > 0 else "under",
                 "flagged": bool(q < alpha),
                 "material": bool(q < alpha and abs(share - base) >= min_lift)}
        if (i, j) in gq:
            entry["q_genus"] = gq[(i, j)]["q"]
            entry["z_genus"] = round(gq[(i, j)]["z"], 2)
            entry["genus_tested"] = True
            entry["survives_genus"] = bool(gq[(i, j)]["q"] < alpha
                                           and gq[(i, j)]["z"] * z[i, j] > 0)
        elif gq:
            # Cins tablosu var ama bu hucrenin beklenen sayisi orada esigin
            # altinda: "gecemedi" degil, "sinanamadi".
            entry["genus_tested"] = False
        # "holds": maddi fark var ve cins duzeyinde sinanabildiyse orada da ayni
        # yonde duruyor. Cins tablosu hic kurulmadiysa olcut yalnizca maddiliktir
        # ve sayfa bunu ayrica belirtir.
        entry["holds"] = bool(entry["material"]
                              and entry.get("survives_genus", True))
        cells.append(entry)

    # satir duzeyi: bir satir vs kalan satirlarin toplami
    rowp, rowinfo = [], []
    for i in range(t.shape[0]):
        rest = t.sum(axis=0) - t[i]
        sub = np.vstack([t[i], rest])
        sub = sub[:, sub.sum(axis=0) > 0]
        info = {"row": rows[i], "n": int(rsum[i])}
        if sub.shape[1] < 2 or sub[0].sum() < 10 or sub[1].sum() < 10:
            info["p"] = None
            rowinfo.append(info)
            continue
        if sub.shape == (2, 2):
            p = float(stats.fisher_exact(sub)[1])
        else:
            p = float(stats.chi2_contingency(sub)[1])
        info["p"] = p
        info["cramers_v"] = round(cramers_v(sub), 3)
        rowp.append((len(rowinfo), p))
        rowinfo.append(info)
    if rowp:
        for (pos, _), q in zip(rowp, bh_adjust([p for _, p in rowp])):
            rowinfo[pos]["q"] = q
            rowinfo[pos]["flagged"] = bool(q < alpha)

    flagged = [c for c in cells if c["flagged"]]
    material = [c for c in cells if c["material"]]
    # Siralama: cins duzeyini de geçenler once, sonra etki buyuklugune gore.
    material.sort(key=lambda c: (0 if c["holds"] else 1, -abs(c["lift"])))
    return {"cells": cells, "rows": rowinfo,
            "top": material[:8],
            "n_cells_tested": len(cells),
            "n_cells_material": len(material),
            "n_cells_material_genus": sum(1 for c in material if c.get("holds")),
            "min_lift": min_lift,
            "n_cells_flagged": len(flagged),
            "n_rows_tested": sum(1 for r in rowinfo if r.get("p") is not None),
            "n_rows_flagged": sum(1 for r in rowinfo if r.get("flagged")),
            "has_genus": bool(gq),
            "n_cells_survive_genus": sum(1 for c in cells if c.get("survives_genus")),
            "alpha": alpha, "min_expected": min_expected}


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
               s.host_kingdom, s.habitat, s.geo,
               (SELECT COUNT(*) FROM neighbor nb JOIN gene_category c ON c.neighbor_id=nb.neighbor_id
                 AND c.method='regex_v1' WHERE nb.candidate_id=r.candidate_id AND c.category='transposon') transposons,
               (SELECT COUNT(*) FROM neighbor nb2 WHERE nb2.candidate_id=r.candidate_id) n_neighbors
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
            "reductase_type", "ferredoxin_type", "domain", "host_kingdom", "habitat",
            "geo", "transposons", "n_neighbors"]
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
        # Elektron tasima SISTEMI: klasik siniflandirmanin asil ekseni. Bazi
        # RO'lar uc bilesenli (reduktaz + ferredoksin + oksijenaz), bazilari iki
        # bilesenli (reduktaz-ferredoksin kaynasmis ya da ayri ferredoksin yok).
        # Bu ayrim reduktaz ve ferredoksin tiplerine AYRI AYRI bakildiginda
        # gorunmez; birlikte bakmak gerekir.
        red = (d.get("reductase_type") or "none") != "none"
        fdx = (d.get("ferredoxin_type") or "none") != "none"
        # BAGLAM YOKLUGU, YOKLUK DEGILDIR. 945 dogrulanmis giris (%8,3) icin
        # pencere icinde HIC komsu yok -- kisa kontig, parcali derleme ya da
        # cekilmemis kayit. Bu girisler simdiye kadar "beta yok" ve "yakinda
        # ortak yok" olarak sayiliyordu, yani kanit yoklugu yokluk kaniti
        # gibi islem goruyordu.
        #
        # Olculdu ve disaridan yakalandi: Miao & Schmidt 2025'in Tablo 1'i 21
        # tip icin DENEYSEL alt birim mimarisi veriyor ve cikarim 20'sinde
        # tutuyor. Tek uyusmazlik (2_201_BphA1) tam olarak bu durum: tipin tek
        # uyesinin hicbir komsusu yok, bu yuzden "beta yok" diye sayiliyor,
        # oysa kristal yapida beta alt birimi var. Yani uyusmazlik cikarimi
        # degil, bu kodlamayi curutuyor.
        d["has_context"] = (d.get("n_neighbors") or 0) > 0
        if not d["has_context"]:
            d["beta_label"] = "no context"
            d["etc_system"] = "no_context"
        else:
            d["beta_label"] = "beta" if d.get("has_beta") else "no beta"
            d["etc_system"] = ("three_component" if red and fdx else
                               "two_component" if red else
                               "ferredoxin_only" if fdx else "none_nearby")
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
    # Reaksiyon sinifi her girise yazilir: "hangi reaksiyon tarzi" sorusu
    # grup ve tip sorulariyla ayni araclarla sorulabilsin.
    reaction_of = {}
    if os.path.exists(args.chemistry):
        with open(args.chemistry) as fh:
            for row in csv.DictReader(fh):
                reaction_of[row["cluster"]] = row.get("reaction_class") or "unknown"
    for d in data:
        d["reaction_class"] = reaction_of.get(d["cluster"], "unknown")
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

    ETC_SYSTEMS = ["three_component", "two_component", "ferredoxin_only", "none_nearby"]

    def contingency(test_id, question, subset, key_field, col_field, keys, cols,
                    row_header, row_prefix="", unit="alpha subunit",
                    restriction=None, with_null=False, with_genus=True, note_key=None):
        """Bir kontenjans testi + sinyal haritasi + istege bagli null ve cins tablosu."""
        table = [[sum(1 for d in subset if d[key_field] == k and d[col_field] == c)
                  for c in cols] for k in keys]
        keys_t, cols_t, table_t = trim_table(keys, cols, table)
        if len(keys_t) < 2 or len(cols_t) < 2:
            return None
        chi2, p, dof, _ = stats.chi2_contingency(table_t)
        gtab = None
        if with_genus:
            gsub = collapse_genus(subset)
            gtab = [[sum(1 for d in gsub if d[key_field] == k and d[col_field] == c)
                     for c in cols_t] for k in keys_t]
        t = {"id": test_id, "question": question,
             "rows": keys_t, "cols": cols_t, "table": table_t,
             "chi2": float(chi2), "dof": int(dof), "p": float(p),
             "cramers_v": cramers_v(table_t),
             "row_header": row_header, "row_prefix": row_prefix, "unit": unit,
             "n": sum(sum(r) for r in table_t)}
        if restriction:
            t["restriction"] = restriction
        if with_null:
            t["null"] = permutation_null(keys_t, cols_t, subset, key_field, col_field)
        t["signals"] = signal_map(keys_t, cols_t, table_t, gtab)
        if gtab is not None:
            t["genus_table"] = gtab
        out["tests"].append(t)
        return t


    # 1. ksenobiyotik x plazmit / mobil element
    out["tests"].append({"id": "plasmid_by_class",
                         **test_binary_by_class(data, "is_plasmid",
                                                "Are xenobiotic-degrading RO types more often plasmid-borne than natural-substrate types?", euk_heavy)})
    out["tests"].append({"id": "mobile_by_class",
                         **test_binary_by_class(data, "mobile",
                                                "Are xenobiotic RO types more often associated with mobile elements (plasmid or transposase within 10 kb)?", euk_heavy)})

    # 2. taksonomik yayilim. ESKI metrik "uye basina cins sayisi" idi ve
    #    OLCULDU: tip buyuklugu ile Spearman korelasyonu -0,766 (p=9e-11),
    #    yani metrik yayilimi degil TIP BUYUKLUGUNU olcuyordu. Sebep mekanik:
    #    yeni uye eklendikce yeni cins bulma olasiligi duser (doygunluk), bu
    #    yuzden kalabalik tipler kendiliginden "dar" gorunur.
    #
    #    Dogrusu SEYRELTME (rarefaction): her tip AYNI sayida uyeye indirilir
    #    ve cins sayisi o ortak buyuklukte sayilir. Seyreltme sonrasi buyukluk
    #    korelasyonu +0,254'e (p=0,10) duser, yani yapisal etki kalkar. Sonuc
    #    da degisir: eski metrik ksenobiyotik tipleri daha GENIS gosteriyordu
    #    (0,316'ya 0,217), seyreltilmis olcumde yon TERSINE doner (12,2'ye
    #    13,3) ve fark zaten anlamli degildir. Yani eski egilim bir artefakti.
    RAREFY_TO = 20          # ortak buyukluk: 43 tip bu esigi gecıyor
    RAREFY_DRAWS = 300      # ortalama bu kadar cekilisten alinir
    rng = random.Random(7)

    def rarefied_genera(genera, size, draws=RAREFY_DRAWS):
        """Tipi `size` uyeye indirip kac cins kaldigini ortalar."""
        if len(genera) < size:
            return None
        return sum(len(set(rng.sample(genera, size))) for _ in range(draws)) / draws

    genus_lists = defaultdict(list)
    sclass_of = {}
    for d in data:
        genus_lists[d["cluster"]].append(d["genus"])
        sclass_of[d["cluster"]] = d["sclass"]
    rare = {}
    for cluster, genera in genus_lists.items():
        value = rarefied_genera(genera, RAREFY_TO)
        if value is not None:
            rare[cluster] = value
    rx = [v for c, v in rare.items() if sclass_of.get(c) == "xenobiotic"]
    rn = [v for c, v in rare.items()
          if sclass_of.get(c) in ("natural_aromatic", "natural_specialized")]
    if rx and rn:
        ru, rp = stats.mannwhitneyu(rx, rn, alternative="two-sided")
        sizes = [len(genus_lists[c]) for c in rare]
        rho, rho_p = stats.spearmanr(sizes, [rare[c] for c in rare])
        out["tests"].append({
            "id": "taxonomic_breadth_rarefied",
            "question": ("Do types acting on man-made substrates occupy more bacterial "
                         "genera than types acting on natural ones, once every type is "
                         "compared at the same sample size?"),
            "unit": f"enzyme type, rarefied to {RAREFY_TO} members",
            "metric": f"distinct genera among {RAREFY_TO} members, mean of "
                      f"{RAREFY_DRAWS} draws",
            "xenobiotic": {"n": len(rx), "median": float(np.median(rx))},
            "natural": {"n": len(rn), "median": float(np.median(rn))},
            "mannwhitney_U": float(ru), "p": float(rp),
            "size_correlation_after_rarefaction": {
                "spearman_rho": float(rho), "p": float(rho_p),
                "note": ("near zero is the point of rarefaction; the unrarefied metric "
                         "correlated with type size at rho = -0.77, p = 9e-11, so it was "
                         "measuring how many members a type has")},
            "restriction": (f"types with at least {RAREFY_TO} members: "
                            f"{len(rare)} of {len(genus_lists)}"),
            "n": len(rare)})

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

    # 3. ETC profili x RO grubu; beta x grup. Hepsi ayni yardimci uzerinden
    # kurulur, boylece her biri hem cins duzeyinde ikinci bir tabloya hem de
    # hucre bazinda sinyal haritasina sahip olur.
    groups = sorted({d["group"] for d in data})
    with_context = [d for d in data if d["has_context"]]
    no_context_n = len(data) - len(with_context)
    context_restriction = (
        f"entries whose genomic neighbourhood was retrieved: {len(with_context)} of "
        f"{len(data)}. The {no_context_n} entries with no neighbours at all are "
        "excluded rather than counted as having no partner, because for them the "
        "question was never asked")
    out["context_coverage"] = {
        "entries": len(data), "with_neighbours": len(with_context),
        "without_neighbours": no_context_n,
        "why_it_matters": (
            "an entry with no retrieved neighbourhood scores as having no beta "
            "subunit, no ferredoxin and no reductase. Counted that way it is "
            "indistinguishable from an enzyme that genuinely works alone, and it "
            "measures how fragmented the assembly was rather than how the enzyme is "
            "organised. The external architecture check against Miao & Schmidt 2025 "
            "turned on exactly this: the single type where the inference disagreed "
            "with a crystal structure was a type whose only member has no neighbours"),
    }
    contingency(
        "reductase_by_group",
        "Is the reductase type found near the alpha subunit associated with the RO group?",
        with_context, "group", "reductase_type", groups,
        ["FNR", "GR", "FNR+GR", "none"], "RO group", row_prefix="group ",
        restriction=context_restriction)
    contingency(
        "ferredoxin_by_group", "Is the ferredoxin type associated with the RO group?",
        with_context, "group", "ferredoxin_type", groups,
        ["rieske", "plant", "rieske+plant", "none"], "RO group", row_prefix="group ",
        restriction=context_restriction)
    # Cramer's V, satir sayisi artinca KENDILIGINDEN yukselir. Tip bazinda
    # tablo 5 satir yerine ~43 satir oldugu icin "tipte daha guclu" demek,
    # once bu yapisal etkiyi dislamayi gerektirir. Etiketler karistirilarak
    # AYNI SEKILDEKI tablo icin bir null dagilimi kurulur; gercek V ancak bu
    # null'in cok uzerindeyse iliski gercektir.
    def permutation_null(keys, values, subset, key_field, value_field, repeats=200):
        import random as _random
        rng = _random.Random(7)
        labels = [d[value_field] for d in subset]
        nulls = []
        for _ in range(repeats):
            rng.shuffle(labels)
            table = [[0] * len(values) for _ in keys]
            index = {k: i for i, k in enumerate(keys)}
            vindex = {v: j for j, v in enumerate(values)}
            for d, label in zip(subset, labels):
                i = index.get(d[key_field])
                j = vindex.get(label)
                if i is not None and j is not None:
                    table[i][j] += 1
            trimmed = trim_table(keys, values, table)[2]
            v = cramers_v(trimmed)
            if v is not None:
                nulls.append(v)
        if not nulls:
            return None
        mean = sum(nulls) / len(nulls)
        var = sum((x - mean) ** 2 for x in nulls) / len(nulls)
        return {"mean": round(mean, 4), "sd": round(var ** 0.5, 4),
                "repeats": len(nulls),
                "max": round(max(nulls), 4)}

    # ETC bileseni x ENZIM TIPI. Grup bazinda sorulan ayni soru, ama tip
    # bazinda: kullanici hakli olarak "belki enzim ozelinde olabilir, grup
    # degil" dedi. Grup bes kategoriye indiriyor ve tip icindeki farki
    # gizleyebilir. Seyreklik gercek bir risk: 71 tip x 4 ferredoksin tipi
    # cogu hucreyi bos birakir, bu yuzden yalnizca MIN_TYPE_N uyesi olan
    # tipler girer ve kac tipin girdigi raporlanir.
    MIN_TYPE_N = 20
    type_counts = Counter(d["cluster"] for d in data)
    big_types = sorted(t for t, n in type_counts.items() if n >= MIN_TYPE_N)
    big_subset = [d for d in with_context if d["cluster"] in big_types]
    type_restriction = (f"types with at least {MIN_TYPE_N} members: "
                        f"{len(big_types)} of {len(type_counts)} types, "
                        f"{len(big_subset)} entries")
    for field, values, label in (
            ("reductase_type", ["FNR", "GR", "FNR+GR", "none"], "reductase type"),
            ("ferredoxin_type", ["rieske", "plant", "rieske+plant", "none"], "ferredoxin type"),
    ):
        contingency(
            f"{field}_by_type",
            f"Is the {label} associated with the enzyme TYPE, rather than only with "
            f"the broad RO group?",
            big_subset, "cluster", field, big_types, values, "Enzyme type",
            restriction=type_restriction, with_null=True)

    # Ayni soru beta alt birimi icin: mimari tip bazinda mi grup bazinda mi
    # belirleniyor?
    contingency(
        "beta_by_type",
        "Is the presence of a beta subunit in the operon decided at the level of the "
        "enzyme type rather than the RO group?",
        big_subset, "cluster", "beta_label", big_types, ["beta", "no beta"],
        "Enzyme type", restriction=type_restriction, with_null=True)

    # Elektron tasima SISTEMI. Kullanicinin isaret ettigi ayrim: bazi RO'lar
    # uc bilesenli sistemdir (ayri reduktaz VE ayri ferredoksin), bazilari iki
    # bilesenli. Reduktaz tipi ve ferredoksin tipi AYRI AYRI test edildiginde bu
    # ayrim gorunmez, cunku her iki testte de "none" hucresi iki farkli biyolojik
    # duruma karsilik gelir: ortak gercekten yok, ya da ortak var ama obur
    # bilesen eksik. Birlesik degisken bu karisikligi kaldirir.
    #
    # Ayni soru UC duzeyde sorulur -- grup, enzim tipi, reaksiyon sinifi --
    # cunku bir duzeyde kaybolan iliski bir digerinde gercek olabilir.
    contingency(
        "etc_system_by_group",
        "Is the electron-transport SYSTEM near the alpha subunit -- three-component "
        "(separate reductase and ferredoxin), two-component (reductase only), or "
        "neither -- associated with the RO group?",
        with_context, "group", "etc_system", groups, ETC_SYSTEMS,
        "RO group", row_prefix="group ", restriction=context_restriction)

    contingency(
        "etc_system_by_type",
        "Is the electron-transport system a property of the individual enzyme type "
        "rather than of the broad RO group?",
        [d for d in with_context if d["cluster"] in big_types], "cluster",
        "etc_system", big_types, ETC_SYSTEMS, "Enzyme type",
        restriction=(f"types with at least {MIN_TYPE_N} members, and only entries whose "
                     f"neighbourhood was retrieved: {len(big_types)} of "
                     f"{len(type_counts)} types"),
        with_null=True)

    reaction_keys = sorted({d["reaction_class"] for d in data
                            if d["reaction_class"] != "unknown"})
    reaction_subset = [d for d in with_context if d["reaction_class"] != "unknown"]
    contingency(
        "etc_system_by_reaction",
        "Is the electron-transport system associated with the kind of reaction the "
        "enzyme performs?",
        reaction_subset, "reaction_class", "etc_system", reaction_keys, ETC_SYSTEMS,
        "Reaction class",
        restriction=("entries whose type carries a curated reaction class: "
                     f"{len(reaction_subset)} of {len(data)}"))

    # Ortaklar reaksiyon sinifiyla da sorulur: "bazi reaksiyon tarzlarinda
    # iliski olabilir" sorusu grup duzeyinde sorulan testin cevabini
    # degistirebilir.
    for field, values, label in (
            ("reductase_type", ["FNR", "GR", "FNR+GR", "none"], "reductase type"),
            ("ferredoxin_type", ["rieske", "plant", "rieske+plant", "none"],
             "ferredoxin type")):
        contingency(
            f"{field}_by_reaction",
            f"Is the {label} associated with the kind of reaction the enzyme performs?",
            reaction_subset, "reaction_class", field, reaction_keys, values,
            "Reaction class",
            restriction=("entries whose type carries a curated reaction class: "
                         f"{len(reaction_subset)} of {len(data)}"))

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

    contingency(
        "beta_by_group",
        "Is a beta subunit in the operon (α₃β₃ architecture) group-specific?",
        with_context, "group", "beta_label", groups, ["beta", "no beta"],
        "RO group", row_prefix="group ", restriction=context_restriction)

    # 3c. KANIT DUZEYI x HABITAT. Kullanicinin sorusu: "okyanusta olanlar uzak
    # akraba mi?" Olculebilir ve cevabi carpici. Ama ne olctugu konusunda
    # dikkatli olmak gerekiyor ve bu ciktida yaziliyor:
    #
    # Kanit duzeyi, girisin KURATORLU bir referansa ne kadar benzedigidir.
    # Kuratorlu referanslar ise kultur koleksiyonlarindan gelir ve o
    # koleksiyonlar toprak, klinik ve endustriyel kokenlidir. Dolayisiyla
    # "deniz = uzak" sonucu iki sekilde okunabilir ve test ikisini AYIRAMAZ:
    #   (a) deniz Rieske oksijenazlari gercekten ayri bir dal,
    #   (b) hicbir deniz enzimi karakterize edilmemis.
    # Ikisi de ilginc, ikisi de ayni sayiyi uretir, ve sayfa bunu soyluyor.
    HABITAT_MIN = 60
    tiers_order = ["characterized", "close_homolog", "family_member",
                   "distant", "novel"]
    hab_rows = [d for d in data
                if d.get("habitat") and d["habitat"] not in ("unknown", "other")
                and d.get("tier")]
    hab_counts = Counter(d["habitat"] for d in hab_rows)
    big_habitats = sorted(h for h, n in hab_counts.items() if n >= HABITAT_MIN)
    contingency(
        "evidence_tier_by_habitat",
        "Are the enzymes we cannot interpret -- the distant and novel members -- "
        "concentrated in particular habitats?",
        [d for d in hab_rows if d["habitat"] in big_habitats],
        "habitat", "tier", big_habitats, tiers_order,
        "Habitat", unit="alpha subunit",
        restriction=(f"habitats with at least {HABITAT_MIN} entries carrying both a "
                     f"habitat and an evidence level: {len(big_habitats)} of "
                     f"{len(hab_counts)} habitats"))

    # 3d. Ayni soru DOGRUDAN: "yorumlayamadiklarimiz" tek bir kategori olarak.
    # Bes duzeyli tabloda deniz girislerinin "distant" hucresi 10 puanlik
    # maddilik esiginin hemen altinda kaliyor (+9,5), ama uzak ve novel
    # BIRLIKTE sorulunca fark aciliyor. Soru zaten birlesik olani soruyor:
    # "bu habitatta bulduklarimizin ne kadarini yorumlayamiyoruz?"
    #
    # Her habitat icin iki-iki Fisher, hem giris hem CINS duzeyinde; cins
    # duzeyi sart, cunku tek bir derin orneklenmis proje bir habitati tek
    # basina tasiyabilir.
    uninterpretable = {"distant", "novel"}
    hab_level = {}
    for level, rows_at_level in (("entry", hab_rows),
                                 ("genus", collapse_genus(hab_rows))):
        per_hab = {}
        praw, keys = [], []
        total_un = sum(1 for d in rows_at_level if d["tier"] in uninterpretable)
        total_n = len(rows_at_level)
        for habitat in big_habitats:
            rows_h = [d for d in rows_at_level if d["habitat"] == habitat]
            if len(rows_h) < 20:
                continue
            a_yes = sum(1 for d in rows_h if d["tier"] in uninterpretable)
            a_no = len(rows_h) - a_yes
            b_yes = total_un - a_yes
            b_no = (total_n - len(rows_h)) - b_yes
            res = two_by_two(a_yes, a_no, b_yes, b_no)
            res["n"] = len(rows_h)
            res["genera"] = len({d["genus"] for d in rows_h})
            per_hab[habitat] = res
            praw.append(res["p"])
            keys.append(habitat)
        for habitat, q in zip(keys, bh_adjust(praw)):
            per_hab[habitat]["q"] = q
            per_hab[habitat]["flagged"] = bool(
                q < 0.05 and abs(per_hab[habitat]["rate_a"]
                                 - per_hab[habitat]["rate_b"]) >= 0.08)
        hab_level[level] = {
            "overall_rate": total_un / total_n if total_n else None,
            "n": total_n, "habitats": per_hab}

    # Bagimsiz bir ikinci dilim: habitat sozlugu yerine kaydin KENDI yer adi.
    # "Pacific Ocean", "Baltic Sea" gibi degerler country alaninda duruyor ve
    # habitat siniflandiricisindan tamamen ayri uretiliyor. Iki alan ayni yone
    # isaret ediyorsa giris duzeyindeki gozlem en azindan tek bir etiketleme
    # kuralinin artefakti degildir.
    # Kelime siniri SART. Basit alt dize eslesmesi "USA:Seattle"i ve
    # "Research Centre"i (icinde "sea" geciyor) deniz sayiyordu. Kurum adlari
    # da ayrica disarida: "Bigelow Laboratory for Ocean Sciences" bir laboratuvar.
    OCEAN_RE = re.compile(r"\b(ocean|sea)\b", re.I)
    NOT_A_PLACE_RE = re.compile(
        r"laborator|institut|universit|centre|center|science|museum|collection", re.I)
    ocean_rows = [d for d in data
                  if d.get("geo") and d["tier"]
                  and OCEAN_RE.search(d["geo"])
                  and not NOT_A_PLACE_RE.search(d["geo"])]
    all_tiered = [d for d in data if d["tier"]]
    cross_check = None
    if len(ocean_rows) >= 50 and all_tiered:
        cross_check = {
            "field": "the record's own place name (country), not the habitat vocabulary",
            "matches": "values containing 'Ocean' or 'Sea'",
            "n": len(ocean_rows),
            "rate": sum(1 for d in ocean_rows
                        if d["tier"] in uninterpretable) / len(ocean_rows),
            "overall_rate": sum(1 for d in all_tiered
                                if d["tier"] in uninterpretable) / len(all_tiered),
            "genera": len({d["genus"] for d in ocean_rows}),
        }

    out["tests"].append({
        "ocean_cross_check": cross_check,
        "id": "uninterpretable_by_habitat",
        "question": ("Which habitats return enzymes we cannot interpret? "
                     "A member is counted as uninterpretable when its nearest "
                     "curated reference is distant or absent, so no substrate "
                     "can be transferred to it."),
        "unit": "alpha subunit",
        "levels_by_habitat": hab_level,
        "what_this_cannot_separate": (
            "An evidence level measures distance to a CURATED reference, and the "
            "curated set comes from culture collections that are overwhelmingly "
            "soil, clinical and industrial in origin. A habitat that scores high "
            "here is therefore either genuinely full of divergent enzymes, or "
            "simply a habitat nobody has characterised an enzyme from. This test "
            "cannot tell those apart, and both are worth knowing."),
        "restriction": (f"habitats with at least 20 entries at the level shown, "
                        f"out of {len(hab_counts)} habitats"),
    })

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

    def carbox_genus_table(rows, residue_index, key_fn, keys, cols):
        """Karboksilat satirlarini kume x cins basina tek gozleme cokert.

        Ayni cinsin yuzlerce susu ayni kalintiyi tasir; onlari bagimsiz gozlem
        saymak hucre testlerini otomatik olarak anlamli yapar. Cins icinde
        cogunluk kalintisi alinir.
        """
        buckets = defaultdict(Counter)
        key_of = {}
        for r in rows:
            genus = (r[4] or "?").split()[0]
            bucket = (r[3], genus)
            buckets[bucket]["Asp" if r[residue_index] == "D" else "Glu"] += 1
            key_of[bucket] = key_fn(r)
        index = {k: i for i, k in enumerate(keys)}
        cindex = {c: j for j, c in enumerate(cols)}
        table = [[0] * len(cols) for _ in keys]
        for bucket, counts in buckets.items():
            i = index.get(key_of[bucket])
            j = cindex.get(counts.most_common(1)[0][0])
            if i is not None and j is not None:
                table[i][j] += 1
        return table

    # 5b. Kopru karboksilati Asp mi Glu mu -- grup ve reaksiyon sinifiyla iliskisi
    #
    # Olculen sey: iki karboksilatin DAVRANISI farkli. Katalitik demir
    # karboksilati neredeyse istisnasiz Asp (476:1), yani bu pozisyon ailenin
    # degismeyen parcasi. Kopru karboksilati ise 2,5:1 bolunuyor ve Glu dagilimi
    # rastgele degil: kuaterner amin dalini isaretliyor. Test bu izlenimi
    # olcuye donusturur.
    carbox = con.execute("""
        SELECT c.bridging_residue, c.catalytic_residue, r.ro_group, r.ro_cluster,
               p.organism
        FROM ro_carboxylate c JOIN ro r USING(candidate_id)
                              JOIN replicon p USING(nucleotide_id)
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
            gtab = carbox_genus_table(carbox_only, 0, key_fn, keys_t, cols_t)
            out["tests"].append({
                "id": test_id,
                "question": f"Among entries that have a bridging carboxylate, is its identity "
                            f"(Asp or Glu) associated with the {label}?",
                "rows": keys_t, "cols": cols_t, "table": table,
                "chi2": float(chi2), "dof": int(dof), "p": float(p),
                "cramers_v": cramers_v(table),
                "genus_table": gtab,
                "signals": signal_map(keys_t, cols_t, table, gtab)})
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
        cat_gtab = carbox_genus_table(cat_only, 1, lambda r: r[2], keys_t, cols_t)
        out["tests"].append({
            "genus_table": cat_gtab,
            "signals": signal_map(keys_t, cols_t, table, cat_gtab),
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
    contingency(
        "phylum_by_group",
        "How are RO groups distributed across phyla? (top 8 phyla)",
        [d for d in data if d["phylum"] in phyla], "group", "phylum", groups, phyla,
        "RO group", row_prefix="group ",
        restriction=f"the eight most frequent phyla: {', '.join(phyla)}")

    # Sinyal haritasi HER kontenjans testine eklenir. Yukarida contingency()
    # ile kurulan testler kendi haritasini zaten tasiyor; burada elle kurulmus
    # olanlar tamamlanir. Boylece sayfadaki her tablo icin "iliski tabloda
    # NEREDE" sorusu ayni araclarla cevaplanir ve toplu p degerinin anlamsiz
    # cikmasi tek basina "iliski yok" demeye yetmez.
    for t in out["tests"]:
        if t.get("signals") or "table" not in t or not isinstance(t.get("rows"), list):
            continue
        if not t["rows"] or not isinstance(t["table"][0], list):
            continue
        t["signals"] = signal_map(t["rows"], t["cols"], t["table"])

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
            sg = t.get("signals") or {}
            extra = ""
            if sg:
                extra = (f"  cells {sg['n_cells_material']}/{sg['n_cells_tested']}"
                         f" material, {sg['n_cells_material_genus']} hold"
                         f"{'' if sg.get('has_genus') else ' (no genus check)'}"
                         f", rows {sg['n_rows_flagged']}/{sg['n_rows_tested']}")
            print(f"  {t['id']:32s} chi2={t['chi2']:.0f} dof={t['dof']} "
                  f"p={t['p']:.1e} V={t['cramers_v']:.2f}{extra}")
        else:
            print(f"  {t['id']:32s} {json.dumps({k: v for k, v in t.items() if k not in ('question', 'id')}, default=float)[:160]}")
    print("  evidence:", dict(zip(tiers, out["evidence"]["counts"])),
          f"transferable={out['evidence']['fraction_substrate_transferable']:.2f}",
          "mostly_distant:", out["evidence"]["clusters_mostly_distant"])


if __name__ == "__main__":
    main()
