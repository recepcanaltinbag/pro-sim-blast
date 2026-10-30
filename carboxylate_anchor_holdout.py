"""
ROAR-DB -- IS THE ANCHORED READING CIRCULAR? A HELD-OUT TEST.

WHY THIS EXISTS
  carboxylate_audit.py step 6 decides the carboxylate column from the vote
  of a type's OWN positives, and then reports the method's accuracy
  (method_recovery_on_positives) on THOSE SAME positives. That number is
  therefore optimistic by construction: the positives helped choose the
  column they are then scored against.

  This module removes the circularity. Per type the anchored positives are
  split deterministically (even index = FIT, odd index = HELD OUT). The
  column is voted by the FIT half only, and accuracy is measured on the
  HELD-OUT half, which had no say in it. It also reports whether the
  fit-half column is the same column the full vote picks -- if a half of
  the data lands on a different column, the anchor is not stable.

  Nothing here is written back to the database; the module is read-only and
  reproducible (mafft --thread 1).
"""
import argparse
import json
import os
import sqlite3
import sys
import tempfile
from collections import Counter

import carboxylate_audit as CA

LOW_VOTE = ["4_406_PDO", "1_108_PobA", "1_112_NdmA", "3_317_NidA"]
LARGEST = ["3_313_KshA15", "1_105_DdmC", "1_106_TsaM", "1_104_VanA",
           "1_111_GxtA", "5_501_CmoS"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    ap.add_argument("--ref-fasta", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--out", default="analysis_out/carboxylate_anchor_holdout.json")
    ap.add_argument("--window", type=int, default=3)
    ap.add_argument("--pos-sample", type=int, default=120)
    ap.add_argument("--clusters", default="default",
                    help="'default' (lowest vote shares + largest types), "
                         "'all', or a comma separated list")
    args = ap.parse_args()

    con = sqlite3.connect("file:%s?mode=ro" % args.db, uri=True)
    entries = CA.load_entries(con)
    con.close()
    rows, rf = CA.parse_stockholm(args.alignment, set(entries))
    match_idx = [i for i, ch in enumerate(rf) if ch != "."]
    index = {}
    for cid, row in rows.items():
        idx = CA.AlignmentIndex(row, match_idx, entries[cid]["sequence"])
        if idx.offset >= 0:
            index[cid] = idx
    refs = CA.read_fasta(args.ref_fasta)

    negs_by, pos_by = {}, {}
    for cid, rec in entries.items():
        bucket = (pos_by if (rec["stored_residue"] or "") in CA.CARBOX_OK
                  else negs_by)
        bucket.setdefault(rec["cluster"], []).append(cid)

    if args.clusters == "all":
        wanted = sorted(negs_by)
    elif args.clusters == "default":
        wanted = [c for c in LOW_VOTE + LARGEST if c in negs_by]
    else:
        wanted = [c.strip() for c in args.clusters.split(",") if c.strip()]

    workdir = tempfile.mkdtemp(prefix="carbox_holdout_")
    out = {"module": os.path.basename(__file__),
           "question": "does the anchor still read the carboxylate in "
                       "positives that had no vote in choosing it",
           "split": "deterministic: positives sorted by id, even index fits "
                    "the anchor, odd index is held out",
           "window_columns": args.window,
           "read_only": True, "reproducible": "mafft --thread 1",
           "per_cluster": []}

    for cluster in wanted:
        negs = sorted(negs_by.get(cluster, []))
        poss = sorted(pos_by.get(cluster, []))
        anchored = []
        for cid in poss:
            idx = index.get(cid)
            if idx is None:
                continue
            raw = idx.find(entries[cid]["sequence"], CA.CARBOX_COL,
                           CA.CARBOX_OK, args.window)[0]
            if raw is not None:
                anchored.append((cid, raw))
        if len(anchored) > args.pos_sample:
            step = max(1, len(anchored) // args.pos_sample)
            anchored = anchored[::step][:args.pos_sample]
        fit = anchored[0::2]
        held = anchored[1::2]
        if len(fit) < 3 or len(held) < 3:
            out["per_cluster"].append(
                {"cluster": cluster, "status": "too_few_positives_to_split",
                 "n_anchored_positives": len(anchored)})
            print("  %-18s too few positives to split (%d)"
                  % (cluster, len(anchored)))
            continue

        records = [(c, entries[c]["sequence"]) for c in negs]
        records += [(c, entries[c]["sequence"]) for c, _ in anchored]
        records += [("REF|" + k, refs[k]) for k in refs
                    if k.startswith(cluster + "_")]
        aln = CA.run_mafft(records, workdir, cluster + "_holdout", 1)
        maps = {k: CA.column_map(v) for k, v in aln.items()}
        width = len(next(iter(aln.values())))

        def vote(pairs):
            v = Counter()
            for cid, raw in pairs:
                col = maps[cid].get(raw)
                if col is not None:
                    v[col] += 1
            return v

        def reads_at(cid, col):
            row = aln[cid]
            for delta in sorted(range(-args.window, args.window + 1),
                                key=lambda x: (abs(x), x)):
                c = col + delta
                if 0 <= c < width and row[c].upper() in CA.CARBOX_OK:
                    return True
            return False

        full_v, fit_v = vote(anchored), vote(fit)
        if not fit_v:
            out["per_cluster"].append({"cluster": cluster,
                                       "status": "no_fit_vote"})
            continue
        full_col, full_top = full_v.most_common(1)[0]
        fit_col, fit_top = fit_v.most_common(1)[0]
        runner_up = (fit_v.most_common(2)[1][1] if len(fit_v) > 1 else 0)
        held_ok = sum(1 for cid, _ in held if reads_at(cid, fit_col))
        neg_found = sum(1 for cid in negs if reads_at(cid, fit_col))

        rec = {
            "cluster": cluster,
            "status": "tested",
            "n_anchored_positives": len(anchored),
            "n_fit": len(fit), "n_held_out": len(held),
            "n_no_carboxylate_members": len(negs),
            "full_vote_column": full_col,
            "full_vote_share": round(full_top / len(anchored), 3),
            "fit_vote_column": fit_col,
            "fit_vote_share": round(fit_top / len(fit), 3),
            "fit_margin_over_runner_up": fit_top - runner_up,
            "same_column_as_full_vote": fit_col == full_col,
            "held_out_recovery": round(held_ok / len(held), 3),
            "held_out_recovered": held_ok,
            "missing_call_members_found_with_fit_anchor": neg_found,
            "missing_call_fraction": (round(neg_found / len(negs), 3)
                                      if negs else None),
        }
        out["per_cluster"].append(rec)
        print("  %-18s fit=%3d held=%3d  held-out recovery %-6s  "
              "same column as full vote: %-5s  margin %d"
              % (cluster, len(fit), len(held), rec["held_out_recovery"],
                 rec["same_column_as_full_vote"],
                 rec["fit_margin_over_runner_up"]))

    tested = [c for c in out["per_cluster"] if c.get("status") == "tested"]
    if tested:
        nh = sum(c["n_held_out"] for c in tested)
        ok = sum(c["held_out_recovered"] for c in tested)
        out["summary"] = {
            "clusters_tested": len(tested),
            "held_out_positives": nh,
            "held_out_recovered": ok,
            "held_out_recovery_pooled": round(ok / nh, 4) if nh else None,
            "clusters_where_half_the_data_picks_the_same_column":
                sum(1 for c in tested if c["same_column_as_full_vote"]),
            "reading": "pooled held-out recovery is the method's accuracy on "
                       "positives that did not help choose the column; it is "
                       "the number to quote instead of "
                       "method_recovery_on_positives",
        }
        print("\n  POOLED held-out recovery: %d/%d = %s"
              % (ok, nh, out["summary"]["held_out_recovery_pooled"]))
        print("  same column as full vote: %d/%d clusters"
              % (out["summary"]["clusters_where_half_the_data_picks_the_same_column"],
                 len(tested)))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(out, open(args.out, "w"), indent=1, ensure_ascii=False)
    print("  -> %s" % args.out)


if __name__ == "__main__":
    sys.exit(main())
