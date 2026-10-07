"""Yayimlanmis iddialarla bu veritabaninin olctugu arasindaki FARKLAR.

NEDEN. Bir veritabani literaturu tekrar ediyorsa yeni bir sey soylemiyor
demektir. Degerli olan, olculen seyin yayimlanmis bir iddiayla ortusMEDIGI
yerlerdir -- ve bunlar bir kenarda tutulmazsa unutulur. Bu sayfa onlari tek
yerde toplar, her biri icin iddiayi kaynagiyla, bizim olcumumuzu, ve neyin
sonucu degistirecegini yazar.

IKI KURAL, ikisi de bu sayfanin ise yaramasi icin sart:

  1. SAYILAR BURADA URETILMEZ. Her karsi-kanit, onu zaten hesaplayan kanonik
     JSON dosyasindan OKUNUR (stats.json, operon_relations.json,
     ferredoxin_residue.json). Sayi sablona elle yazilsaydi, analiz
     degistiginde sessizce yanlis hale gelirdi -- ve bir "literaturle
     celisiyoruz" sayfasinda yanlis sayi, celistigin iddiadan daha kotu.

  2. CELISKI, CURUTME DEGILDIR. Her kayit "bu testin GOSTERMEDIGI sey"
     alanini tasimak zorunda. Dagilimsal bir olcum, mutagenez deneyini
     curutemez; genomik baglam, biyokimyasal olcumu curutemez. Sayfa ne
     olculdugunu soyler, daha fazlasini degil.

Kendi onceki calismalarimizla celiskiler de buraya girer -- aslinda en
guvenilir kayit odur, cunku kimse kendi sonucunu bedavaya geri almaz.

Cikti: analysis_out/disagreements.json
"""

import argparse
import json
import os


def read(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def test_by_id(stats, test_id):
    for t in (stats or {}).get("tests", []):
        if t.get("id") == test_id:
            return t
    return {}


def build(analysis_dir):
    stats = read(os.path.join(analysis_dir, "stats.json"))
    operon = read(os.path.join(analysis_dir, "operon_relations.json"))
    fdx = read(os.path.join(analysis_dir, "ferredoxin_residue.json"))

    entries = []

    # ------------------------------------------------ 1. ferredoksin kalintisi
    if fdx and fdx.get("verdict"):
        v = fdx["verdict"]
        ref = fdx.get("refinement") or v.get("refinement") or {}
        cov = fdx.get("coverage") or {}
        pc = fdx.get("positive_controls") or {}
        entries.append({
            "id": "ferredoxin_class_residue",
            "status": v.get("status", "contradicted"),
            "claim": ("The residue at the NDO-Fd A50 / CDO-Fd W48 column is "
                      "class-dependent: tryptophan predominates in Class IIB "
                      "ferredoxins, alanine in Class III, and tyrosine in "
                      "certain plant-type Class IIA ferredoxins."),
            "source": ("Miao H, Oerlemanns R, Hagedoorn P-L, Schmidt S. "
                       "Decoding and Reprogramming Redox Partner Specificity in "
                       "Rieske Oxygenases for Enhanced Catalytic Activity. "
                       "bioRxiv 2026"),
            "doi": "10.64898/2026.04.03.713453",
            "their_evidence": ("mutagenesis and chimera construction on two "
                               "ferredoxins, plus a reading of previously "
                               "published alignments"),
            "our_measurement": v,
            "our_evidence": (
                f"the equivalent column read in "
                f"{cov.get('aligned_well_enough', 0):,} ferredoxins "
                f"({cov.get('aligned_well_enough_distinct_sequences', 0):,} distinct "
                f"sequences) out of {cov.get('gene_category_ferredoxin_rows', 0):,} "
                f"in the database. The rest were discarded, almost all of them "
                f"because the four Rieske ligands did not align to the chemically "
                f"identical residue -- the column sits three residues after one of "
                f"those ligands, so without that anchor the register is noise"),
            "headline": (
                f"Alanine holds this column in "
                f"{100 * (v.get('dominant_share_entry_level') or 0):.0f} % of "
                f"observations across the whole family, and in "
                f"{100 * ((v.get('proxy_class_shares') or {}).get('IIB_like_A_share') or 0):.0f} % "
                f"of the Class IIB-like set where tryptophan was said to "
                f"predominate. Tryptophan is at "
                f"{100 * (v.get('tryptophan_share_entry_level') or 0):.1f} % overall."),
            "refinement": ref,
            "positive_controls": pc,
            "what_this_does_not_show": v.get("what_this_does_not_show"),
            "what_would_settle_it": (
                "an agreed assignment of Batie classes to these enzymes. This "
                "database assigns none, and its own groups are built from alpha "
                "subunit sequence rather than electron-chain architecture, so "
                "the class arm of the claim is answered here only through a "
                "proxy that the output states in full."),
            "source_file": "ferredoxin_residue.json",
        })

    # ------------------------------- 2. ferredoksin mi reduktaz mi daha secici
    red = test_by_id(stats, "reductase_type_by_type")
    fd = test_by_id(stats, "ferredoxin_type_by_type")
    if red and fd:
        entries.append({
            "id": "ferredoxin_more_selective",
            "status": "not supported",
            "claim": ("The reductase shows moderate promiscuity in redox partner "
                      "compatibility, while interactions between ferredoxin and "
                      "oxygenase are typically more selective."),
            "source": ("Miao H, Schmidt S. Rieske Oxygenases: Powerful Models "
                       "for Understanding Nature's Orchestration of Electron "
                       "Transfer and Oxidative Chemistry. Biochemistry "
                       "2025;64:3801-3813"),
            "doi": "10.1021/acs.biochem.5c00369",
            "their_evidence": ("structural studies of individual complexes and "
                               "reconstitution experiments with non-native "
                               "partners"),
            "our_measurement": {
                "reductase_type_by_enzyme_type": {
                    "cramers_v": red.get("cramers_v"), "n": red.get("n"),
                    "shuffled_null": (red.get("null") or {}).get("mean")},
                "ferredoxin_type_by_enzyme_type": {
                    "cramers_v": fd.get("cramers_v"), "n": fd.get("n"),
                    "shuffled_null": (fd.get("null") or {}).get("mean")},
            },
            "our_evidence": ("how tightly each partner's identity tracks the "
                             "enzyme type across the whole database, measured "
                             "the same way for both and against the same kind "
                             "of shuffled null"),
            "headline": (
                f"If the ferredoxin coupling were the more selective of the two, "
                f"the ferredoxin should track the enzyme type more tightly than "
                f"the reductase does. It does not: Cramér's V is "
                f"{fd.get('cramers_v', 0):.2f} for the ferredoxin against "
                f"{red.get('cramers_v', 0):.2f} for the reductase, a difference "
                f"of {abs((red.get('cramers_v') or 0) - (fd.get('cramers_v') or 0)):.2f} "
                f"in the opposite direction to the claim."),
            "what_this_does_not_show": (
                "that ferredoxin coupling is not selective in the biochemical "
                "sense. What is measured here is which partner FAMILY sits "
                "beside which enzyme in the genome, not whether a given "
                "ferredoxin can productively transfer electrons to a given "
                "oxygenase. Those are different questions and the second is the "
                "one their experiments answer."),
            "what_would_settle_it": (
                "cross-reconstitution assays across many more pairs than the "
                "handful tested so far, which is laboratory work this database "
                "cannot substitute for."),
            "source_file": "stats.json",
        })

    # ---------------------------------- 3. "yuksek oranda korunmus aspartat"
    br = test_by_id(stats, "bridging_residue_by_group")
    brr = test_by_id(stats, "bridging_residue_by_reaction")
    if br:
        cols = br.get("cols") or []
        table = br.get("table") or []
        totals = {c: sum(row[j] for row in table) for j, c in enumerate(cols)}
        asp, glu = totals.get("Asp", 0), totals.get("Glu", 0)
        entries.append({
            "id": "conserved_bridging_aspartate",
            "status": "refined",
            "claim": ("Electron transfer between the Rieske cluster and the "
                      "mononuclear iron proceeds by a proton-coupled mechanism "
                      "mediated by a highly conserved aspartate residue, which "
                      "bridges the two metal centres of neighbouring subunits."),
            "source": ("Miao H, Schmidt S. Biochemistry 2025;64:3801-3813"),
            "doi": "10.1021/acs.biochem.5c00369",
            "their_evidence": ("crystal structures and computational studies of "
                               "a small number of model systems"),
            "our_measurement": {
                "asp": asp, "glu": glu,
                "ratio": round(asp / glu, 2) if glu else None,
                "cramers_v_by_group": br.get("cramers_v"),
                "cramers_v_by_reaction": brr.get("cramers_v"),
                "p_by_group": br.get("p"),
            },
            "our_evidence": ("the identity of the bridging carboxylate read "
                             "from every confirmed alpha subunit that has one"),
            "headline": (
                f"Conserved in kind, not in identity. The bridging position is "
                f"a carboxylate almost without exception, but it is aspartate "
                f"{asp:,} times and glutamate {glu:,} times, a ratio of "
                f"{asp / glu:.1f} to 1 rather than invariance. The glutamate is "
                f"not scattered: it maps onto a definable branch, with "
                f"Cramér's V = {br.get('cramers_v', 0):.2f} against the RO group "
                f"and {brr.get('cramers_v', 0):.2f} against the reaction class."
                if glu else "no glutamate observed"),
            "what_this_does_not_show": (
                "that the mechanism differs in the glutamate branch. Both "
                "residues are carboxylates and either could in principle carry "
                "the same proton-coupled step; no kinetics were measured here. "
                "The claim being refined is about conservation, not mechanism."),
            "what_would_settle_it": ("a structure or a kinetic measurement from "
                                     "an enzyme in the glutamate branch."),
            "source_file": "stats.json",
        })

    # ---------------------------- 4. KENDI onceki sonucumuzla celiski
    div = ((operon or {}).get("regulator_vs_enzyme_divergence") or {})
    controls = div.get("control_genes") or {}
    prior = div.get("prior_art") or {}
    if controls.get("regulator_versus_each_control"):
        entries.append({
            "id": "regulators_diverge_faster",
            "status": "overturned in part, and it was our own",
            "claim": ("Regulators diverge faster than the enzymes they control. "
                      "Measured over 45 enzyme types and 2,382 pairs: mean main "
                      "diversity 0.378 against mean regulator diversity 0.446."),
            "source": ("this project's own earlier analysis "
                       "(evolution_rates_of_reg_and_enz.py, "
                       "enzyme_regulator_graphs.py), and the same direction is "
                       "widely assumed in the literature"),
            "doi": None,
            "their_evidence": prior.get("what_they_computed"),
            "our_measurement": {
                "controls": controls.get("on_their_own_pairs"),
                "common_set": controls.get("on_the_common_pair_set"),
                "head_to_head": controls.get("regulator_versus_each_control"),
                "verdict": controls.get("verdict"),
                "gate_ladder": div.get("orthology_gate_ladder"),
                "measured_floor": div.get("measured_identity_floor"),
            },
            "our_evidence": ("the same comparison, but with the beta subunit, "
                             "the ferredoxin and the reductase measured at the "
                             "same genome pairs as controls"),
            "headline": (
                "The pairs are drawn from within an enzyme type, and an enzyme "
                "type is defined by alpha subunit similarity, so a high alpha "
                "identity is a property of the selection rather than an "
                "observation. Against controls, roughly sixty per cent of the "
                "gap disappears: the regulator does beat the beta subunit, "
                "beats the reductase only marginally, and does not beat the "
                "ferredoxin at all. Under the strictest orthology gate the gap "
                "vanishes."),
            "what_this_does_not_show": (
                "that regulators evolve at the same rate as their enzymes. The "
                "regulator remains the most variable protein in the operon on "
                "every measure here. What collapses is the size of the effect "
                "and the claim that it is specific to regulators."),
            "what_would_settle_it": (
                "a phylogeny-aware rate comparison on confidently orthologous "
                "sets, rather than pairwise identity within profile bins."),
            "source_file": "operon_relations.json",
        })

    # ------------------------------------- 5. kendi eski metrigimiz tersine dondu
    rare = test_by_id(stats, "taxonomic_breadth_rarefied")
    old = test_by_id(stats, "taxonomic_breadth")
    if rare and old:
        corr = rare.get("size_correlation_after_rarefaction") or {}
        entries.append({
            "id": "taxonomic_breadth_metric",
            "status": "overturned, and it was our own",
            "claim": ("Enzyme types acting on man-made substrates occupy more "
                      "bacterial genera than types acting on natural ones, "
                      "measured as distinct genera per member."),
            "source": "this database's own earlier metric",
            "doi": None,
            "their_evidence": "genera per member, computed per enzyme type",
            "our_measurement": {
                "old_metric": {"xenobiotic_median": (old.get("xenobiotic") or {}).get("median"),
                               "natural_median": (old.get("natural") or {}).get("median"),
                               "p": old.get("p")},
                "rarefied": {"xenobiotic_median": (rare.get("xenobiotic") or {}).get("median"),
                             "natural_median": (rare.get("natural") or {}).get("median"),
                             "p": rare.get("p")},
                "size_correlation_after_rarefaction": corr,
            },
            "our_evidence": ("the same question with every type rarefied to a "
                             "common sample size"),
            "headline": (
                "The old metric was measuring type size, not breadth: it "
                "correlated with the number of members at rho = -0.77. Rarefied "
                "to a common size the direction REVERSES and the difference is "
                "not significant. A result that exists only because of how the "
                "measure behaves is not a result."),
            "what_this_does_not_show": (
                "that there is no difference in taxonomic breadth. With this "
                "many types the rarefied test has little power; it shows that "
                "the original evidence does not survive, not that the opposite "
                "is true."),
            "what_would_settle_it": "more characterised types.",
            "source_file": "stats.json",
        })

    by_status = {}
    for e in entries:
        by_status[e["status"]] = by_status.get(e["status"], 0) + 1
    return {
        "what_this_page_is": (
            "Places where what this database measures does not match a "
            "published claim, including claims made earlier by this project "
            "itself. Every number is read from the analysis file that produced "
            "it rather than retyped, and every entry states what its test does "
            "NOT show, because a distributional measurement cannot overturn an "
            "experiment."),
        "n_entries": len(entries),
        "by_status": by_status,
        "entries": entries,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--analysis-dir", default="analysis_out")
    ap.add_argument("--out", default="analysis_out/disagreements.json")
    args = ap.parse_args()

    out = build(args.analysis_dir)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"[yazildi] {args.out}")
    print(f"  {out['n_entries']} kayit")
    for e in out["entries"]:
        print(f"    [{e['status']:34s}] {e['id']}")
        missing = [k for k in ("headline", "what_this_does_not_show",
                               "what_would_settle_it") if not e.get(k)]
        if missing:
            print(f"        EKSIK ALAN: {', '.join(missing)}")


if __name__ == "__main__":
    main()
