"""Her onaylanmis uyenin HANGI katalitik kolonlarda isabet verdigini cikar.

NEDEN BU BETIK VAR. Kapi (ro_motif.classify_motifs) katalitik triad icin
">=2/3" diyor ve DB yalnizca `catalytic_intact` ikili degerini sakliyor: 2/3
ile 3/3 ayirt EDILEMIYOR. Bu iki sorunu gizliyor:

  1. Karboksilat kolonu YANLIS (madde 73: 26 kristal yapinin 26'sinda cagirdigi
     kalinti demirden 10,9-14,7 A uzakta). Dolayisiyla 2/3 olan bir uyenin iki
     isabetinden biri o kolon olabilir -- yani uye GERCEKTE tek ligandla kabul
     edilmis olur.
  2. Hizalama kolona HIC VARMAMIS olabilir. check_sites bu durumda kalintiyi
     "-" okuyup eksik sayiyor, yani "hizalanmadi" ile "hizalandi ama kalinti
     yok" ayni gorunuyor. Kapi acisindan bu GUVENLI taraftir (eksik sayilan
     dislanir), ama kayitta ikisi ayni yazildigi icin neyin olculdugu
     bilinmiyor.

Bu betik SINIFLANDIRMAYI DEGISTIRMEZ, yalnizca olcer ve raporlar.
"""
import argparse
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ro_motif  # noqa: E402


def covered(aligned, column, window):
    """Dizinin hizalamasi bu kolonu KAPSIYOR mu (bosluk olsa bile)?

    Kapsamak ile isabet etmek ayri sey. Hizalama kolondan once bitmisse o
    kolon hakkinda hicbir sey OLCULMEMISTIR.
    """
    for off in range(-window, window + 1):
        i = column - 1 + off
        if 0 <= i < len(aligned) and aligned[i] != ".":
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", default="genomic_context/cand_aln.sto")
    ap.add_argument("--out", default="analysis_out/motif_site_audit.json")
    ap.add_argument("--db", default="roar.sqlite")
    a = ap.parse_args()

    confirmed = None
    if os.path.exists(a.db):
        import sqlite3
        con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
        confirmed = {r[0] for r in con.execute(
            "select candidate_id from ro where is_confirmed=1")}
        con.close()

    aligned = ro_motif.read_stockholm_matchcols(a.aln)
    w = ro_motif.COLUMN_WINDOW
    cat = ro_motif.CATALYTIC_SITES
    rie = ro_motif.RIESKE_SITES
    labels = [s[2] for s in cat]

    pattern = collections.Counter()
    n_total = n_conf = 0
    carboxylate_only_pairs = []
    uncovered = collections.Counter()
    for name, seq in aligned.items():
        key = name.split()[0]
        is_conf = confirmed is None or key in confirmed
        n_total += 1
        if not is_conf:
            continue
        n_conf += 1
        hits = tuple(
            ro_motif.check_sites(seq, [s], w)[0] == 1 for s in cat)
        rh = sum(ro_motif.check_sites(seq, [s], w)[0] for s in rie)
        pattern[(rh, hits)] += 1
        for s in cat:
            if not covered(seq, s[0], w):
                uncovered[s[2]] += 1
        # Iki His'ten biri yok ve karboksilat "var": isabetin yarisi artefakt
        if hits[2] and (not hits[0] or not hits[1]):
            carboxylate_only_pairs.append(key)

    rows = []
    for (rh, hits), n in sorted(pattern.items(), key=lambda kv: -kv[1]):
        rows.append({"rieske_found": rh,
                     "catalytic": {labels[i]: bool(hits[i])
                                   for i in range(len(labels))},
                     "catalytic_found": sum(hits), "members": n})
    out = {
        "_about": "Her onaylanmis uyenin katalitik kolon isabet deseni. "
                  "Siniflandirma DEGISTIRILMEDI, yalnizca olculdu.",
        "_alignment": a.aln,
        "_window": w,
        "_columns": [{"column": s[0], "accepted": s[1], "label": s[2]}
                     for s in cat],
        "aligned_sequences": n_total,
        "confirmed_scored": n_conf,
        "patterns": rows,
        "columns_not_covered_by_alignment": dict(uncovered),
        "carboxylate_compensating_for_a_missing_his": {
            "members": len(carboxylate_only_pairs),
            "note": "Bu uyelerde iki His'ten en az biri bulunamadi ama "
                    "karboksilat kolonu isabet verdi. Karboksilat kolonu "
                    "yanlis oldugu icin (madde 73) bu isabet metal "
                    "koordinasyonunun KANITI DEGIL: bu uyeler gercekte tek "
                    "gercek ligandla >=2/3 saymis olabilir.",
            "examples": sorted(carboxylate_only_pairs)[:20],
        },
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1, ensure_ascii=False)
    print(f"hizalanan: {n_total} | puanlanan onaylanmis: {n_conf}")
    for r in rows[:8]:
        print(f"  Rieske {r['rieske_found']}/4  katalitik "
              f"{r['catalytic_found']}/3  {r['catalytic']}  -> {r['members']}")
    print("hizalamanin KAPSAMADIGI kolonlar:", dict(uncovered))
    print("karboksilat eksik His'i telafi eden uye:",
          len(carboxylate_only_pairs))
    print("->", a.out)


if __name__ == "__main__":
    main()
