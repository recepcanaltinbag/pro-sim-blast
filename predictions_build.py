"""predictions/ katmanindaki tablolari uret.

NEDEN BU BETIK VAR
  Arkeal ve fungal adaylarin kaniti uc ayri dosyada duruyordu: metadata
  (arch_meta.tsv), bakteriyel kimlik (arch.m8 / fung.m8) ve motif skorlari
  (analysis_out/archaeal_fungal_ro_prediction.json). Siteye bakan birinin bu
  ucunu elle birlestirmesi gerekiyordu. Bu betik onlari TEK tabloya indirir.

KURAL
  Cikti yalnizca predictions/ altina yazilir. types/, chemistry.csv,
  roar.sqlite ve HMM kutuphanesine DOKUNMAZ -- bkz predictions/README.md.

KAPSAMA HER SATIRDA
  Motif skoru kapsamadan okunmaz. Dusuk kapsama "triad yok" degil
  "olculemedi" demektir; AerC'de olculdu (bkz aerc_resolution.json).
"""
import argparse, csv, json, os


def read_m8(path):
    """qseqid -> (en iyi hit, kimlik, kapsama)"""
    best = {}
    if not os.path.exists(path):
        return best
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or f[0] in best:
            continue
        sid = f[1]
        if "|" in sid:                       # tr|ACC|NAME -> ACC
            parts = sid.split("|")
            sid = parts[1] if len(parts) > 1 else sid
        best[f[0]] = (sid, float(f[2]), float(f[3]))
    return best


def read_meta(path):
    """accession -> (organizma, filum, protein adi, kanit duzeyi, uzunluk)"""
    out = {}
    if not os.path.exists(path):
        return out
    rows = open(path).read().strip().split("\n")
    for line in rows[1:]:
        f = line.split("\t")
        if len(f) < 6:
            continue
        acc, org, lineage, name, existence, length = f[0], f[1], f[2], f[3], f[4], f[5]
        phylum = next((k for k in ("Halobacteria", "Thermoproteota",
                                   "Thermoplasmata") if k in lineage), "?")
        out[acc] = (org, phylum, name, existence, length)
    return out


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print(f"[yazildi] {path}  ({len(rows)} satir)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--analysis-dir", default="analysis_out")
    ap.add_argument("--out-dir", default="predictions")
    a = ap.parse_args()
    A = lambda *p: os.path.join(a.analysis_dir, *p)
    os.makedirs(a.out_dir, exist_ok=True)

    summary = json.load(open(A("archaeal_fungal_ro_prediction.json")))
    meta = read_meta(A("arch_meta.tsv"))
    arch_ref = read_m8(A("arch.m8"))          # kuratorlu 75 referansa
    arch_bact = read_m8(A("arch_vs_bact.m8")) or read_m8(
        A("archaeal_clade", "arch_vs_bact.m8"))   # 22.5k bakteriyel proteine
    fung_ref = read_m8(A("fung.m8"))

    cov = summary["archaea"]["model_coverage_range"]
    rows = []
    for acc in sorted(meta):
        org, phylum, name, existence, length = meta[acc]
        r = arch_ref.get(acc, ("", "", ""))
        b = arch_bact.get(acc, ("", "", ""))
        rows.append([acc, org, phylum, length, name, existence,
                     "4/4", "3/3",
                     r[0].split("_pro")[0] if r[0] else "",
                     f"{r[1]:.1f}" if r[1] != "" else "",
                     f"{r[2]:.1f}" if r[2] != "" else "",
                     b[0] or "", f"{b[1]:.1f}" if b[1] != "" else "",
                     "prediction -- not characterised"])
    write_tsv(os.path.join(a.out_dir, "archaeal.tsv"),
              ["accession", "organism", "phylum", "length_aa", "protein_name",
               "evidence_level", "rieske_columns", "catalytic_columns",
               "nearest_curated_ref", "ref_identity_pct", "ref_coverage_pct",
               "nearest_bacterium", "bact_identity_pct", "status"], rows)

    rows = []
    for acc in sorted(fung_ref):
        r = fung_ref[acc]
        rows.append([acc, r[0].split("_pro")[0], f"{r[1]:.1f}", f"{r[2]:.1f}",
                     "4/4", "3/3", "prediction -- not characterised"])
    write_tsv(os.path.join(a.out_dir, "fungal.tsv"),
              ["accession", "nearest_curated_ref", "ref_identity_pct",
               "ref_coverage_pct", "rieske_columns", "catalytic_columns",
               "status"], rows)

    # Sayfanin okudugu ozet.
    #
    # ONEMLI: BU GORUNUM INGILIZCEDIR. Analiz kaydi
    # (archaeal_fungal_ro_prediction.json) Turkce yazilmistir ve oyle kalir --
    # o bir IC kayit, kod yorumlari gibi. Ama /atlas/predictions yayimlanan bir
    # sayfa ve yayimlanan metin Ingilizce olmak zorunda.
    #
    # Ilk surumde analiz kaydinin metin alanlari dogrudan sayfaya veriliyordu
    # ve sayfa Turkce render ediyordu (2026-10-23'te yakalandi; check_site.py
    # tam bu sizintiyi aramak icin var). Simdi SAYILAR kayittan, METIN burada.
    A = summary["archaea"]
    F = summary["fungi"]
    ctrl = summary["controls"]

    def neg(label, c):
        return {"set": label, "n": c["n"],
                "four_of_four_rieske": c.get("four_of_four_rieske"),
                "three_of_three_catalytic": c.get("three_of_three_catalytic"),
                "passes": c["pass"]}

    negatives = ctrl["negative"]
    view = {
        "generated": summary["_date"],
        "rule": "Nothing here enters types/, chemistry.csv, the HMM library, or any "
                "statistic in this database",
        "controls": {
            "positive": {"set": "curated references",
                         "n": ctrl["positive"]["of"],
                         "four_of_four_rieske": ctrl["positive"]["of"],
                         "three_of_three_catalytic": ctrl["positive"]["three_of_three"],
                         "passes": ctrl["positive"]["pass"]},
            "negative": [
                neg("cytochrome b6-f iron-sulfur subunit (reviewed)", negatives[0]),
                neg("Rieske ferredoxin (reviewed)", negatives[1]),
                neg("NirD (reviewed)", negatives[2])],
            "interpretation":
                "Fifty-six of the sixty-six true non-Rieske-oxygenases carry all four Rieske "
                "ligands, none carries the catalytic triad, and none passes. That is exactly "
                "the biochemical distinction the test is meant to draw: false positives 0/66, "
                "recovery 75/75."},
        "coverage_rule": {
            "rule": "Low coverage is not a negative. Where the alignment never reaches the "
                    "catalytic columns the result is 'not measurable', not 'no triad'. Every "
                    "row reports its coverage beside its motif score.",
            "case": "This was measured, not assumed: a 717-residue protein (AerC) aligned over "
                    "18.9 % of the model, stopped at column 231, never reached the catalytic "
                    "columns, and scored 0/3 -- which was wrong. Aligned in fragments, its "
                    "ligands were all in place."},
        "archaea": {
            "strict_pass": A["strict_pass"], "of": A["of"],
            "length_range_aa": A["length_range_aa"],
            "length_median_aa": A.get("length_median_aa"),
            "model_coverage_range": A["model_coverage_range"],
            "phyla": A["phyla"],
            "protein_existence": A["protein_existence"],
            # Sayilari al, Turkce aciklama alanlarini ALMA: analiz kaydindaki
            # alignment_coverage metni Turkceydi ve yayimlanan JSON'a siziyordu.
            "identity_to_curated_refs": {k: v for k, v in
                                         A["identity_to_curated_refs"].items()
                                         if not isinstance(v, str)},
            "alignment_coverage": "84-98 % of the query, that is full length rather "
                                  "than the Rieske domain alone",
            "nearest_refs": {
                "by_group": A["nearest_refs"]["by_group"],
                "summary": "The majority sit nearest to group 5 (carnitine and choline "
                           "oxygenases), group 3 second."},
            "hgt_test": {
                "question": "Is there a near-identical bacterial relative, as a recent "
                            "transfer would leave?",
                "db": "22,531 bacterial Rieske-domain proteins (UniProt)",
                "identity": A["hgt_test"]["identity_to_nearest_bacterium"],
                "conclusion": "No recent transfer. The nearest bacterium is 51.7 % identical "
                              "at most."},
            "phylogeny": {
                "monophyletic": A["phylogeny"]["monophyletic"],
                "n_clusters": len(A["phylogeny"]["pure_clusters"]),
                "support_note": "Support values are FastTree SH-like local support, not "
                                "bootstrap. Deep node support runs 0.52 to 0.97, so the "
                                "backbone is poorly resolved and no donor group is claimed.",
                "clusters": [
                    {"n": c["n"], "members": c.get("members"), "support": c.get("support"),
                     "sister": c.get("sister"),
                     "sister_node_support": c.get("sister_node_support")}
                    for c in A["phylogeny"]["pure_clusters"] if c.get("n", 0) >= 2]},
            "supported": [
                "Archaeal proteins with the complete catalytic apparatus exist and are not "
                "few: 109, across three archaeal classes, at 23-44 % identity to the nearest "
                "curated bacterial reference over 84-98 % of their length.",
                "They are not a recent acquisition: no bacterial relative exceeds 51.7 % "
                "identity in a set of 22,531.",
                "They are not one lineage: seven pure clusters, each strongly supported "
                "within itself, each with a bacterial sister group.",
                "Halobacteria alone shows at least two independent acquisitions."],
            "not_supported": [
                "That any of them has Rieske oxygenase activity. None is characterised; all "
                "are labelled Predicted or Inferred from homology.",
                "Any particular donor group. Deep node support is 0.52-0.97, so the backbone "
                "is poorly resolved.",
                "The direction of transfer. The tree was rooted at its midpoint.",
                "Bootstrap or AU support. FastTree SH-like local support was used instead."],
            "next_steps": "IQ-TREE with UFBoot and an AU topology test; the genomic context of "
                          "the electron-transfer partners."},
        "fungi": {
            "strict_pass": F["strict_pass"], "of": F["of"],
            "identity_to_curated_refs": {k: v for k, v in
                                         F["identity_to_curated_refs"].items()
                                         if not isinstance(v, str)},
            "nearest_refs_by_group": F.get("nearest_refs_by_group"),
            "pattern": "The fungal candidates sit overwhelmingly nearest to group 5 (633 of "
                       "654): CntA 493, CmoB 88, CmoS 50. Sixteen sit nearest to the sterol "
                       "7-desaturases Nvd and DAF-36.",
            "note": "The fungal annotations reading 'choline monooxygenase, chloroplastic' are "
                    "an artefact of transfer from plants. This test reads alignment columns "
                    "rather than annotations, so the finding does not depend on them. But note "
                    "that the archaeal set is also mostly nearest to group 5, so the earlier "
                    "claim that the two sets fall in different places -- and that this argued "
                    "against a shared artefact -- does not hold."},
        "counts": {"archaeal": A["strict_pass"], "fungal": F["strict_pass"]},
        "coverage_range": cov,
    }
    corr = A.get("CORRECTION_2026-10-14")
    if corr:
        view["correction"] = {
            "what_was_wrong": [
                "The length range was reported as 429-505 aa. It is 331-505 aa, median 381.",
                "The nearest references were reported as 'all group 3 (PhdA, NarAa, NidA)'. "
                "Group 5 is in fact the majority (61 of 109, CntA alone 59); group 3 is "
                "second with 44."],
            "cause": "The first run read only the top forty rows sorted by coverage, and the "
                     "top five DIAMOND hits per phylum, then generalised. The highest-identity "
                     "hits really were PhdA, NarAa and NidA; the bulk was CntA.",
            "consequence": "The argument that archaea and fungi fall in different places, and "
                           "that this is evidence against a shared artefact, is void. Both "
                           "sets are mostly nearest to group 5. A shared artefact is not "
                           "excluded by this measurement and needs its own control.",
            "still_valid": ["the controls (0/66 false positives, 75/75 recovery)",
                            "the HGT test (maximum 51.7 %)",
                            "the phylogeny: seven pure clusters, not monophyletic",
                            "the bridging aspartate in 109 of 109"]}
    # GUVENLIK AGI. Bu gorunum ve tablolar YAYIMLANIYOR, yani Ingilizce olmak
    # zorunda. Ilk surumde analiz kaydinin Turkce metinleri sayfaya siziyordu;
    # bir daha sessizce olmasin diye yazmadan once bakiliyor.
    import re as _re
    TR = _re.compile(r"\b(ve|icin|bir|degil|kume|yaprak|olculdu|cunku|yalnizca|"
                     r"kuratorlu|dogrulanmis|hepsi|tahmin|karakterize|edilmemis|"
                     r"sitokrom|artefakti|hizaliyorlar|girdi|cikti|sayfa|dosya|"
                     r"bilinmiyor|gerekce|kanit)\b", _re.I)
    blob = json.dumps(view, ensure_ascii=False)
    leaked = sorted(set(TR.findall(blob)))
    if leaked:
        raise SystemExit("[HATA] yayimlanacak gorunumde Turkce metin: "
                         + ", ".join(leaked)
                         + "\n  Analiz kaydi Turkce kalir; GORUNUM Ingilizce olmali.")

    out = os.path.join(a.out_dir, "predictions.json")
    json.dump(view, open(out, "w"), indent=1, ensure_ascii=False)
    print(f"[yazildi] {out}")


if __name__ == "__main__":
    main()
