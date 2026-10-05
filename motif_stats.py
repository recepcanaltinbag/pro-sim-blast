"""
Korunmus merkez istatistigi -- "bu veritabanindaki tum RO'lar ayni yapiyi mi
tasiyor" sorusunun ölçülmüş cevabi.

Dogrulanmis her RO icin ROmotif71 hizalamasindaki 8 tanimlayici kolon tek tek
sayilir:
    Rieske [2Fe-2S] merkezi   C-x-H ... C-x-x-H   (4 ligand, 4/4 zorunlu)
    Mononukleer Fe(II) merkezi 2-His-1-karboksilat (3 ligand, >=2/3 zorunlu)
    Alt birimler arasi kopru  Asp/Glu              (zorunlu degil, raporlanir)

Cikti: analysis_out/motif_stats.json -- web ana sayfasi bu dosyayi okur, boylece
"ortak yapi" iddiasi metin degil olcum olur.
"""

import argparse
import json
import os
import sqlite3
from collections import Counter

from ro_motif import (BRIDGING_SITE, CATALYTIC_SITES, RIESKE_SITES,
                      check_sites, read_stockholm_matchcols)

LABELS_EN = {
    "Rieske Cys-1": "Rieske cysteine 1",
    "Rieske His-1": "Rieske histidine 1",
    "Rieske Cys-2": "Rieske cysteine 2",
    "Rieske His-2": "Rieske histidine 2",
    "mononukleer Fe His-1": "mononuclear iron histidine 1",
    "mononukleer Fe His-2": "mononuclear iron histidine 2",
    "mononukleer Fe karboksilat": "mononuclear iron carboxylate (Asp or Glu)",
    "alt-birimler arasi elektron transfer Asp": "inter-subunit bridging aspartate",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    ap.add_argument("--out", default="analysis_out/motif_stats.json")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    confirmed = [r[0] for r in con.execute(
        "SELECT candidate_id FROM ro WHERE is_confirmed=1")]
    con.close()
    print(f"[okunuyor] hizalama: {args.alignment}")
    aligned = read_stockholm_matchcols(args.alignment)
    seqs = {c: aligned[c] for c in confirmed if c in aligned}
    print(f"[bilgi] {len(seqs)} / {len(confirmed)} dogrulanmis RO hizalamada")

    per_site = {}
    for sites in (RIESKE_SITES, CATALYTIC_SITES, [BRIDGING_SITE]):
        for column, accepted, label in sites:
            found = sum(1 for s in seqs.values()
                        if check_sites(s, [(column, accepted, label)])[0] == 1)
            per_site[label] = {
                "label": LABELS_EN.get(label, label),
                "column": column,
                "expected": "/".join(accepted),
                "found": found,
                "fraction": round(found / max(1, len(seqs)), 4),
            }

    rieske_counts, triad_counts, bridging = Counter(), Counter(), 0
    for s in seqs.values():
        r = check_sites(s, RIESKE_SITES)[0]
        c = check_sites(s, CATALYTIC_SITES)[0]
        rieske_counts[r] += 1
        triad_counts[c] += 1
        if check_sites(s, [BRIDGING_SITE])[0]:
            bridging += 1

    out = {
        "n": len(seqs),
        "sites": [per_site[k] for k in
                  [s[2] for s in RIESKE_SITES] + [s[2] for s in CATALYTIC_SITES]
                  + [BRIDGING_SITE[2]]],
        "rieske_ligands": {str(k): v for k, v in sorted(rieske_counts.items(), reverse=True)},
        "catalytic_triad": {str(k): v for k, v in sorted(triad_counts.items(), reverse=True)},
        "bridging_present": bridging,
        "model": {"name": "ROmotif71", "length": 426, "references": 71},
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"[yazildi] {args.out}")
    for s in out["sites"]:
        print(f"   kolon {s['column']:>4}  {s['expected']:4s} {s['label']:42s} "
              f"{s['found']:6d}  {100 * s['fraction']:.1f}%")
    print("   Rieske ligand sayisi:", dict(out["rieske_ligands"]))
    print("   Katalitik triad sayisi:", dict(out["catalytic_triad"]))
    print(f"   Kopru Asp/Glu: {bridging} ({100 * bridging / max(1, len(seqs)):.1f}%)")


if __name__ == "__main__":
    main()
