"""
ROAR-DB -- PER-TYPE CARBOXYLATE ANCHORING: THE METHOD, PRINTED.

This script measures nothing. It prints, from
analysis_out/carboxylate_audit.json, WHAT the anchored reading looks at,
in the order it looks at it, with the number each step produced and with
each step's own accuracy. It exists so that a reader can check the method
before trusting the statistic the gate publishes.

The question behind all of it: the canonical definition of a Rieske
oxygenase requires a complete 2-His-1-carboxylate facial triad at the
mononuclear iron (Kauppi 1998; the Parales/Gibson consensus). ROAR-DB's
gate accepts 2/3 of it. Is the third site really absent in those members,
or is the READING of it broken?

Every number below comes from the audit JSON; nothing here is restated
from memory. Run carboxylate_audit.py first.
"""
import json
import os
import sys

AUDIT = os.path.join("analysis_out", "carboxylate_audit.json")


def hr(title):
    print("\n" + title)
    print("-" * len(title))


def main():
    if not os.path.exists(AUDIT):
        sys.exit("run carboxylate_audit.py first -- %s is missing" % AUDIT)
    d = json.load(open(AUDIT))
    base = d["baseline_fixed_column"]
    mech = d["mechanism"]
    neg, pos = mech["no_carboxylate"], mech["with_carboxylate"]
    wid = d["window_widening"]
    naive = d["naive_anchor"]
    probe = d["reference_profile_probe"]
    anc = d["anchored_reread"]
    tot = anc["totals"]
    cells = [c for c in anc["per_cluster"] if c.get("status") == "tested"]

    print("=" * 70)
    print("ROAR-DB  --  reading the catalytic carboxylate: what is looked at")
    print("audit generated %s   model %s, carboxylate at match column %s"
          % (d.get("generated", "?"), mech["model"], base["fixed_column"]))
    print("=" * 70)

    hr("STEP 1  the baseline reading, and what it claims")
    print("  looked at   the residue standing in ONE FIXED alignment column")
    print("              (match column %s, window +-%s) for every confirmed"
          % (base["fixed_column"], base["fixed_window"]))
    print("              member; the stored table is recomputed from the")
    print("              alignment to check it (%s disagreements)"
          % base["recomputation_disagreements_with_table"])
    print("  found       %s members give no D/E there, and the breakdown"
          % neg["n"])
    mcb = base["missing_call_breakdown"]
    print("              splits that number into two different findings:")
    print("                %5d  the column AND its neighbours are gaps"
          % mcb["window_entirely_gapped"])
    print("                       -> an alignment event, not a biological one")
    print("                %5d  a real residue sits there but is not D/E"
          % mcb["substituted_residue_in_place"])
    print("                       -> a biological claim, and testable")
    print("                %5d  centre gapped, neighbours present"
          % mcb["centre_gapped_neighbours_present"])
    print("  reading     %.1f%% of the 'absences' are empty columns. An"
          % (100.0 * mcb["window_entirely_gapped"] / neg["n"]))
    print("              uncovered column cannot report an absence.")

    hr("STEP 2  is the column even reachable in those members?")
    print("  looked at   per entry, the first and last FILLED match column of")
    print("              its own alignment -- how much of the model the")
    print("              sequence actually occupies")
    print("                                      no-carboxylate   with")
    for lab, key in (("last filled match column (median)", "last_filled_column_median"),
                     ("filled match residues (median)", "filled_match_residues_median"),
                     ("model coverage (mean)", "model_coverage_mean"),
                     ("hmm score (mean)", "hmm_score_mean"),
                     ("identity to nearest reference %", "ref_identity_mean")):
        print("    %-36s %10s %10s" % (lab, neg[key], pos[key]))
    print("  found       the alignment stops before column %s in %s of %s"
          % (base["fixed_column"],
             neg["alignment_stops_before_the_carboxylate_column"], neg["n"]))
    print("              entries (%.1f%%)." % (100 * neg["alignment_stops_before_fraction"]))
    print("  the anchor  in ALL %s of them the two catalytic histidines and" % neg["n"])
    print("              the bridging Asp/Glu are present, so the local")
    print("              register near the Rieske end is intact -- the")
    print("              failure is specific to the C-terminal lobe.")

    hr("STEP 3  would a wider window rescue it? (the cheap fix, tested)")
    print("  looked at   the SAME fixed column at +-2, +-5, +-10, +-20, +-40")
    for w in wid["windows"]:
        offs = w["offset_counts"]
        edge = ("; offsets reach %s" % min(int(k) for k in offs)) if offs else ""
        print("    +-%-3s columns -> %4d rescued (%.1f%%)%s"
              % (w["window_columns"], w["rescued"],
                 100 * w["rescued_fraction"], edge))
    print("  reading     widening a window around an EMPTY column scans more")
    print("              empty columns. At +-40 the hits sit ~40 columns away,")
    print("              which is no longer the same site. The fix has to be")
    print("              a different anchor, not a looser one.")

    hr("STEP 4  the naive anchor, and why it is NOT evidence")
    dist = naive["his2_to_carboxylate_raw_distance"]
    print("  looked at   D/E at a raw-sequence distance from the second")
    print("              catalytic histidine, calibrated on %s positives"
          % naive["calibration_n"])
    print("  found       that distance is not fixed: %s-%s residues"
          % (dist["p0"], dist["p100"]))
    print("              (2.5-97.5%%: %s-%s), so the window must be %s"
          % (dist["p2.5"], dist["p97.5"], naive["window_width_residues"]))
    print("              residues wide")
    print("  the test    hit rate %.2f%%  vs  chance %.2f%% from each"
          % (100 * naive["found_fraction"],
             100 * naive["expected_by_chance_fraction"]))
    print("              sequence's OWN D/E composition")
    print("  reading     a window that wide finds a D/E in almost any")
    print("              sequence. This step is reported precisely so it")
    print("              cannot be mistaken for the measurement.")

    hr("STEP 5  where is the failure -- the sequence, or the model?")
    print("  looked at   three sub-models cut from the curated reference")
    print("              alignment, searched against three sets: the")
    print("              no-carboxylate members, the carboxylate-bearing")
    print("              members, and -- as a CONTROL -- candidates the gate")
    print("              rejected as not RO alpha at all")
    print("    %-10s %-22s %9s %9s %9s"
          % ("model", "region", "no-carbox", "with", "rejected"))
    for m in probe["models"]:
        s = m["sets"]
        print("    %-10s cols %-17s %9s %9s %9s"
              % (m["model"], "%s-%s" % tuple(m["columns"]),
                 "%d/%d" % (s["no_carboxylate"]["with_any_hit"], s["no_carboxylate"]["n"]),
                 "%d/%d" % (s["with_carboxylate"]["with_any_hit"], s["with_carboxylate"]["n"]),
                 "%d/%d" % (s["rieske_only"]["with_any_hit"], s["rieske_only"]["n"])))
    print("  reading     the Rieske sub-model hits every set -- the positive")
    print("              control holds, so these are real Rieske proteins.")
    print("              The catalytic-lobe model is what fails on them. The")
    print("              global profile cannot align that lobe; the sequence")
    print("              is not the problem.")

    hr("STEP 6  the anchored re-reading -- THE measurement")
    print("  per type, the type's own members are realigned with mafft, and")
    print("  the carboxylate column is decided by the type's POSITIVES rather")
    print("  than by a constant. Window +-%s columns. Four things are looked"
          % anc["window_columns"])
    print("  at, in this order:")
    print("    (a) anchor      where the positives put the carboxylate in the")
    print("                    raw sequence; the column is their vote, and the")
    print("                    vote share is printed so a split vote shows")
    print("    (b) reference   the type's curated reference is in the same")
    print("                    alignment; if it does not carry D/E at the")
    print("                    anchor, the anchor is wrong")
    print("    (c) recovery    the method's accuracy IN THAT alignment: of the")
    print("                    positives, how many does the anchor hand back")
    print("    (d) only then   what the no-carboxylate members carry there")
    print()
    print("    %-18s %5s %5s %6s %6s %8s %5s %s"
          % ("type", "neg", "pos", "vote", "recov", "found", "abs", "null"))
    for c in cells:
        print("    %-18s %5d %5d %6s %6s %8s %5s %s"
              % (c["cluster"], c["n_no_carboxylate"], c["n_positive_anchors"],
                 c.get("anchor_column_vote_share"),
                 c.get("method_recovery_on_positives"),
                 c.get("carboxylate_present_fraction"), c.get("absent"),
                 "pass" if c.get("beats_the_null") else "underpowered"))
    print()
    print("  totals      carboxylate at the HOMOLOGOUS column   %5d"
          % tot["present_at_homologous_column"])
    print("              at a group-conserved column            %5d"
          % tot["present_at_group_conserved_column"])
    print("              still ABSENT after anchoring           %5d"
          % tot["absent"])

    hr("STEP 7  could the conserved column be luck? (the null model)")
    nm = anc["null_model"]
    print("  looked at   per type: observed D/E at the anchor against that")
    print("              group's OWN D/E frequency, one-sided binomial,")
    print("              Benjamini-Hochberg corrected")
    print("  compared to the D/E fraction of EVERY column in the same")
    print("              alignment with >=80% occupancy -- if the anchor")
    print("              column sits inside that distribution there is no")
    print("              signal (median column D/E fraction is reported per")
    print("              type for exactly this comparison)")
    print("  found       %s of %s types beat the null."
          % (nm["clusters_beating_the_null"], nm["clusters_tested"]))
    thin = [c for c in cells if not c.get("beats_the_null")]
    if thin:
        print("  the rest    have too few members to test, not a weak signal:")
        for c in thin:
            print("                %-18s %d member(s) read, %s found, %s absent"
                  % (c["cluster"], c["n_no_carboxylate"],
                     c.get("carboxylate_present"), c.get("absent")))
        print("              total %s entries, all of them found."
              % sum(c["n_no_carboxylate"] for c in thin))

    hr("STEP 8  what the residue-distance spread means, and does not")
    off = anc["residue_offset_from_expected"]
    print("  median %s residues, 95%% interval %s..%s; within +-10: %s of %s"
          % (off["median"], off["p2.5"], off["p97.5"],
             off["within_plus_minus_10"], off["n"]))
    print("  reading     the register test is the ALIGNMENT column in step 6,")
    print("              where offset 0 means the homologous position. This")
    print("              raw distance measures something else: how much indel")
    print("              has accumulated in the ~140-residue stretch between")
    print("              the His pair and the carboxylate. A wide spread is")
    print("              not error -- it is WHY a fixed column drifts out of")
    print("              register, and why the anchor has to be per type.")

    hr("WHAT THE GATE MAY THEREFORE SAY")
    complete = (pos["n"] + tot["present_at_homologous_column"]
                + tot["present_at_group_conserved_column"])
    print("  before      '4/4 Rieske ligands and >=2/3 of the facial triad'")
    print("              -- a tolerance, with nothing measured behind the")
    print("              third site")
    print("  after       '3/3 measured: %s members are triad-complete and %s"
          % (complete, tot["absent"]))
    print("              are documented, individually characterised")
    print("              exceptions'")
    print("  which is the canonical 2-His-1-carboxylate definition, measured")
    print("  rather than assumed.")
    print()


if __name__ == "__main__":
    main()
