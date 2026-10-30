"""Manifestodaki `rhea` capasindan oksijen kaderini HESAPLA ve kayitli degerle
karsilastir.

NEDEN AYRI BIR BETIK. Hesaplama once reaction_verification.json'dan okunan
Rhea kimlikleriyle yapiliyordu, yani capa bir ARA DOSYADA duruyordu. Capa
VERIDE durmali: tipin kendi manifestosunda. Boylece bir capa elle eklendiginde
(EC numarasi metinden cikarilarak, bkz. anchor_reactions.py) hesaplama onu
kendiliginden gorur.
"""
import argparse
import csv
import glob
import json
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reaction_balance as rb  # noqa: E402

OUT = "analysis_out/oxygen_fate_computed.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    resolved = {}
    fpath = os.path.join(a.root, "oxygen_fate.csv")
    if os.path.exists(fpath):
        with open(fpath) as fh:
            resolved = {r["cluster"]: r["oxygen_fate"]
                        for r in csv.DictReader(fh)}

    rows = []
    for f in sorted(glob.glob(os.path.join(a.root, "types", "*.yaml"))):
        m = yaml.safe_load(open(f))
        raw = str(m.get("rhea") or "").strip()
        if not raw:
            continue
        ids = [x.strip() for x in raw.split(",") if x.strip()]
        res = [rb.derive_fate(x) for x in ids]
        fates = sorted({d["fate"] for d in res if d["fate"]})
        rows.append({"cluster": m["cluster"],
                     "resolved_in_db": resolved.get(m["cluster"]),
                     "computed": fates, "details": res})
        print(f"[{m['cluster']}] {fates}", file=sys.stderr)
    rb.save_cache()

    decided = [r for r in rows if r["computed"] and r["computed"] != ["ambiguous"]]
    amb = [r for r in rows if r["computed"] == ["ambiguous"]]
    agree = [r for r in decided if r["resolved_in_db"] in r["computed"]]
    dis = [r for r in decided if r["resolved_in_db"] not in r["computed"]]
    json.dump({"_about": "The fate of the two oxygen atoms COMPUTED from the "
                         "Rhea anchor in the manifest, set against the value "
                         "resolved in the database. A reaction whose atoms the "
                         "count cannot follow is returned as ambiguous.",
               "_counts": {"anchored": len(rows), "computed": len(decided),
                           "agree": len(agree), "disagree": len(dis),
                           "ambiguous": len(amb)},
               "rows": rows}, open(os.path.join(a.root, a.out), "w"),
              indent=1, ensure_ascii=False)
    print(f"\ncapali tip: {len(rows)} | HESAPLANDI {len(decided)} "
          f"(uyusan {len(agree)}, FARK {len(dis)}) | sayim cozemiyor "
          f"{len(amb)}")
    for r in dis:
        print(f"  FARK {r['cluster']}: db={r['resolved_in_db']} "
              f"hesap={r['computed']}")


if __name__ == "__main__":
    main()
