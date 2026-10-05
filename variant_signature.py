"""
Varyant duzeyinde kalinti imzasi -- "benzerlik var ama ayni enzim degil" iddiasinin
dogrudan kaniti.

GEREKCE
  analysis_out/reference_pairs.csv gosterdi ki kuratorlu referanslar arasinda FARKLI
  substratli ciftler %99,8 kimlige kadar cikiyor. Demek ki bu ailede substrat secimi
  global dizi benzerligiyle degil, az sayida pozisyonla belirleniyor. O halde dogru
  gorsel "iki protein ne kadar benzer" degil, "AYIRT EDICI POZISYONLARDA ne farkli".

NE HESAPLANIR
  Her kume icin, uyeleri varyantlara (leaf) bolunmus halde:
    within  = varyant ici kalinti entropisinin ortalamasi (varyant ici tutarlilik)
    between = varyant konsensuslarinin entropisi (varyantlar arasi ayrisma)
    skor    = between - within
  Skor yuksek kolonlar, varyantlari AYIRAN pozisyonlardir. Her varyant icin bu
  kolonlardaki konsensus kalintisi bir imza dizgisi olarak saklanir.

NE HESAPLANMAZ
  Bunlar yapisal olarak dogrulanmis substrat baglama cebi kalintilari DEGILDIR.
  Veritabaninda yapi yok; olculen sey "varyantlari ayiran hizalama kolonlari".
  Katalitik merkez kolonlari (metal ligandlari) ayrica isaretlenir, cunku onlar
  neredeyse degismez ve ayirt edici olmalari beklenmez -- tersi bir sinyal
  hizalama hatasina isaret eder.

Cikti: cluster_sdp + leaf_sdp tablolari, analysis_out/variant_signatures.csv
"""

import argparse
import csv
import json
import math
import os
import sqlite3
from collections import Counter, defaultdict

from ro_motif import MODEL_COLUMNS, DEFAULT_MODEL, read_stockholm_matchcols

MIN_CLUSTER_MEMBERS = 20     # altinda varyant karsilastirmasi anlamsiz
MIN_LEAF_MEMBERS = 3         # bir varyantin konsensusu icin en az uye
MIN_LEAVES = 2               # karsilastirma icin en az varyant
TOP_COLUMNS = 12             # varyant sayfalarinda gosterilecek kolon sayisi

# --- Neden bu kadar sert filtre ---
# Ilk denemede olcut sadece "between - within" entropisiydi. Sonuc: heterojen
# kumelerde (ornek CntA, 235 varyant) secilen kolonlarda her varyant BASKA bir
# kalinti gosteriyordu (skor ~3,5 bit = ~11 farkli durum). Bunlar spesifisite
# pozisyonu degil, hizalamanin en guvenilmez oldugu ilmek bolgeleridir: varyant
# ici entropi dusuk (uyeler neredeyse ayni), varyantlar arasi entropi azami.
# Gercek bir spesifisite pozisyonu AZ SAYIDA alternatif durum alir ve her
# durum birden fazla varyantta tekrarlanir. Asagidaki uc kapi bunu zorlar.
MIN_OCCUPANCY = 0.9          # kolon kume uyelerinin en az bu kadarinda doldurulmus olmali
MAX_STATES = 4               # varyant konsensuslari en fazla bu kadar farkli kalinti olabilir
MIN_TOP2_SHARE = 0.7         # en sik iki durum varyantlarin en az bu kadarini kapsamali
MAX_WITHIN_ENTROPY = 0.5     # kalinti varyant ICINDE korunmus olmali (bit)

SCHEMA = """
CREATE TABLE IF NOT EXISTS cluster_sdp (
    cluster    TEXT PRIMARY KEY,
    n_members  INTEGER,
    n_leaves   INTEGER,
    columns    TEXT,   -- JSON: [{column, score, within, between, consensus, is_ligand, domain}]
    method     TEXT
);
CREATE TABLE IF NOT EXISTS leaf_sdp (
    leaf_id    TEXT PRIMARY KEY,
    cluster    TEXT,
    size       INTEGER,
    residues   TEXT,   -- JSON: {column: residue}
    signature  TEXT,   -- okunabilir: "212H 217H 355D ..."
    n_differ   INTEGER -- kume konsensusundan kac kolonda ayriliyor
);
CREATE INDEX IF NOT EXISTS idx_leafsdp_cluster ON leaf_sdp(cluster);
"""


def entropy(counter):
    total = sum(counter.values())
    if total <= 1:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counter.values() if n)


def ligand_columns():
    model = MODEL_COLUMNS[DEFAULT_MODEL]
    cols = {c for c, _, _ in model["rieske"]} | {c for c, _, _ in model["catalytic"]}
    cols.add(model["bridging"][0])
    return cols


# Domain sinirlari YAKLASIKTIR. Model icinde Rieske ligandlari kolon 108'de
# bitiyor, ilk katalitik ligand 212'de basliyor; arada ikisine de atanamayan bir
# bolge var. Yapisal dogrulama yapilmadigi icin bolgeler konvansiyon olarak
# isaretlenir ve arayuzde oyle sunulur.
RIESKE_END = 120
CATALYTIC_START = 191


def zone_of(column):
    if column <= RIESKE_END:
        return "rieske"
    if column < CATALYTIC_START:
        return "linker"
    return "catalytic"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    ap.add_argument("--out-dir", default="analysis_out")
    ap.add_argument("--top", type=int, default=TOP_COLUMNS)
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    con.executescript(SCHEMA)
    con.execute("DELETE FROM cluster_sdp")
    con.execute("DELETE FROM leaf_sdp")

    print("[okunuyor] hizalama")
    aligned = read_stockholm_matchcols(args.alignment)
    length = len(next(iter(aligned.values())))
    ligands = ligand_columns()

    leaves = defaultdict(list)          # leaf_id -> [candidate_id]
    cluster_of_leaf = {}
    for candidate_id, cluster, leaf_id in con.execute(
            "SELECT candidate_id, cluster, leaf_id FROM ro_leaf"):
        if candidate_id in aligned:
            leaves[leaf_id].append(candidate_id)
            cluster_of_leaf[leaf_id] = cluster
    by_cluster = defaultdict(list)
    for leaf_id, members in leaves.items():
        by_cluster[cluster_of_leaf[leaf_id]].append(leaf_id)
    print(f"[bilgi] {len(leaves)} varyant, {len(by_cluster)} kume")

    cluster_rows, leaf_rows, csv_rows = [], [], []
    skipped_no_column, skipped_small = [], []
    for cluster, leaf_ids in sorted(by_cluster.items(),
                                    key=lambda kv: -sum(len(leaves[l]) for l in kv[1])):
        usable = [l for l in leaf_ids if len(leaves[l]) >= MIN_LEAF_MEMBERS]
        n_members = sum(len(leaves[l]) for l in leaf_ids)
        if n_members < MIN_CLUSTER_MEMBERS or len(usable) < MIN_LEAVES:
            skipped_small.append(cluster)
            continue

        # kolon basina varyant ici / varyantlar arasi entropi
        per_leaf_columns = {}
        for leaf_id in usable:
            seqs = [aligned[m] for m in leaves[leaf_id]]
            per_leaf_columns[leaf_id] = [
                Counter(s[i] for s in seqs if s[i] != "-") for i in range(length)]

        all_seqs = [aligned[m] for l in leaf_ids for m in leaves[l]]
        scored = []
        for i in range(length):
            occupancy = sum(1 for s in all_seqs if s[i] != "-") / len(all_seqs)
            if occupancy < MIN_OCCUPANCY:
                continue
            withins, states = [], Counter()
            for leaf_id in usable:
                counter = per_leaf_columns[leaf_id][i]
                if sum(counter.values()) < MIN_LEAF_MEMBERS:
                    continue
                withins.append(entropy(counter))
                states[counter.most_common(1)[0][0]] += 1
            if len(states) < 2 or len(states) > MAX_STATES or len(withins) < MIN_LEAVES:
                continue
            voted = sum(states.values())
            top2 = sum(n for _, n in states.most_common(2))
            if top2 / voted < MIN_TOP2_SHARE:
                continue
            within = sum(withins) / len(withins)
            if within > MAX_WITHIN_ENTROPY:
                continue
            between = entropy(states)
            scored.append((between - within, i + 1, within, between,
                           states.most_common(1)[0][0],
                           "".join(f"{r}{n}" for r, n in states.most_common())))
        scored.sort(reverse=True)
        picked = scored[:args.top]
        if not picked:
            skipped_no_column.append(cluster)
            continue

        columns_meta = [{
            "column": col, "score": round(score, 4), "within": round(w, 4),
            "between": round(b, 4), "consensus": cons, "states": st,
            "is_ligand": col in ligands,
            "domain": zone_of(col),
        } for score, col, w, b, cons, st in picked]
        cluster_rows.append((cluster, n_members, len(leaf_ids),
                             json.dumps(columns_meta), "variant_entropy_v1"))

        for leaf_id in leaf_ids:
            seqs = [aligned[m] for m in leaves[leaf_id]]
            residues, differ = {}, 0
            for meta in columns_meta:
                i = meta["column"] - 1
                counter = Counter(s[i] for s in seqs if s[i] != "-")
                residue = counter.most_common(1)[0][0] if counter else "-"
                residues[str(meta["column"])] = residue
                if residue != meta["consensus"]:
                    differ += 1
            signature = " ".join(f"{m['column']}{residues[str(m['column'])]}"
                                 for m in columns_meta)
            leaf_rows.append((leaf_id, cluster, len(seqs), json.dumps(residues),
                              signature, differ))
            csv_rows.append({
                "leaf_id": leaf_id, "cluster": cluster, "size": len(seqs),
                "n_differ_from_cluster_consensus": differ,
                "signature": signature,
                "columns": ";".join(str(m["column"]) for m in columns_meta),
            })

    con.executemany("INSERT INTO cluster_sdp VALUES (?,?,?,?,?)", cluster_rows)
    con.executemany("INSERT INTO leaf_sdp VALUES (?,?,?,?,?,?)", leaf_rows)
    con.commit()

    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, "variant_signatures.csv")
    with open(path, "w", newline="") as fh:
        if csv_rows:
            w = csv.DictWriter(fh, fieldnames=list(csv_rows[0].keys()))
            w.writeheader()
            w.writerows(csv_rows)

    print(f"[cluster_sdp] {len(cluster_rows)} kume")
    print(f"[leaf_sdp]    {len(leaf_rows)} varyant")
    ligand_hits = sum(1 for r in cluster_rows
                      for m in json.loads(r[3]) if m["is_ligand"])
    catalytic = sum(1 for r in cluster_rows
                    for m in json.loads(r[3]) if m["domain"] == "catalytic")
    rieske_zone = sum(1 for r in cluster_rows
                      for m in json.loads(r[3]) if m["domain"] == "rieske")
    total_cols = sum(len(json.loads(r[3])) for r in cluster_rows)
    print(f"[kontrol] ayirt edici kolonlar: katalitik bolge "
          f"{100 * catalytic / max(1, total_cols):.1f}%, Rieske bolgesi "
          f"{100 * rieske_zone / max(1, total_cols):.1f}%, geri kalani baglanti bolgesi")
    print(f"[kontrol] {ligand_hits} kolon metal ligandi pozisyonu. Motif testi +-2 kolon "
          f"penceresi kullandigi icin ligandin TAM kolonu kayabilir; bu yuzden ligand "
          f"pozisyonunun burada cikmasi eksik merkez anlamina gelmez.")
    if leaf_rows:
        differ = [r[5] for r in leaf_rows]
        print(f"[kontrol] varyantlar kume konsensusundan ortalama "
              f"{sum(differ) / len(differ):.1f} / {args.top} kolonda ayriliyor")
    print(f"[atlandi] {len(skipped_small)} kume kucuk veya tek varyantli; "
          f"{len(skipped_no_column)} kumede filtreyi gecen kolon yok "
          f"(bu kumelerde hizalama varyant ayrimi icin yeterince guvenilir degil)")
    print(f"[yazildi] {path}")
    con.close()


if __name__ == "__main__":
    main()
