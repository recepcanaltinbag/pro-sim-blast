"""Kanit duzeyi, BIRLESIM (`ro_union`) tarafi icin.

SORU. Tip sayfalarindaki kademe tablosu yalnizca GENOMIK tablodan
(`ro_evidence`) hesaplaniyor, cunku bir kademe ancak bir genomda
tuttugumuz dizi icin hesaplanabiliyordu. Sonuc, 2026-11 oturumunda
dogrudan yanlis okumaya yol acti: `5_505_qxyA` sayfasinda
`close_homolog` satiri 0 goruntu ve bu "bu tipin yakin homologu yok"
diye okundu. Oysa birlikte o tipe atanmis, kapiyi gecen binden fazla
dizi vardi -- yalnizca genomik baglamlari yoktu. Eksik olan sey veri
degil, O DIZILER ICIN HESAPLANMIS BIR KADEMEYDI.

Bu betik ayni kademeleri birlik dizileri icin de hesaplar.

OLCUM TEK YERDEN GELIR. Esikler, diamond'in cagrilma bicimi ve
kaplamanin hangi alandan okundugu `evidence_tiers.py`'den IMPORT
EDILIR, yeniden yazilmaz. Ikinci bir uygulama esikler ayni kalsa bile
sessizce ayrisirdi; bu projede tam bu hata birkac kez olctu
(hmmalign kolonlari, sigma70 null'u, DOI cozumleyicisi).

IKI SAYI AYNI BIRIM DEGIL -- VE BU RAPORDA YAZILI. `ro` tablosu
GIRISLERI sayar: ayni protein yuz genomda bulunursa yuz satir olur.
`ro_union` ise FARKLI DIZILERI sayar, yani tekilleştirilmiş. Dolayisiyla
"genomik 21.046 giris" ile "birlik 69.245 dizi" toplamlari
KARSILASTIRILAMAZ; oranlari da degil. Bu betigin cevabi her iki havuz
icin AYRI AYRI dagilim vermektir, bir fark degil. ROADMAP'te daha once
"%71 dizi kapsamasi" iddiasi tam bu yuzden geri cekilmisti.

BAGLAM SUTUNU KORUNUR. Her satir `has_genomic_context` tasir, boylece
sayfa "baglamiyla bilinen" ile "yalnizca protein duzeyinde bilinen"
kademeleri ayri gosterebilir. Ikisini tek sayida toplamak, kademe
tablosunun hangi paydayi saydigini yine belirsiz birakirdi.

BU SAYI NE DEGIL
  * Bir dizinin `characterized` kademesinde olmasi, o dizinin referans
    enzimin substratini cevirdigini GOSTERMEZ; kademe enzim TIPI
    guvenini olcer (bkz. evidence_tiers.py'nin esik kalibrasyonu).
  * Birlik kademeleri genomik kademelerin yerine GECMEZ: baglamsiz bir
    dizinin operonu, komsulugu ve replikonu yoktur, yani konum gerektiren
    hicbir istatistige giremez.

Cikti: `ro_union_evidence` tablosu + analysis_out/union_evidence.json
"""

import argparse
import json
import os
import sqlite3
import time
from collections import Counter, defaultdict

import evidence_tiers

SCHEMA = """
DROP TABLE IF EXISTS ro_union_evidence;
CREATE TABLE ro_union_evidence (
    seq_md5             TEXT PRIMARY KEY,
    nearest_ref         TEXT,     -- en yakin kuratorlu referans (model adi)
    ref_identity        REAL,     -- % kimlik (diamond)
    ref_qcov            REAL,     -- sorgu kaplamasi %
    same_as_assigned    INTEGER,  -- en yakin referans, atanan kume ile ayni mi
    tier                TEXT,
    has_genomic_context INTEGER   -- paydayi ayirmak icin KORUNUR
);
CREATE INDEX idx_rue_tier ON ro_union_evidence(tier);
CREATE INDEX idx_rue_ctx  ON ro_union_evidence(has_genomic_context);
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", default="analysis_out/union_evidence.json")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    rows = con.execute("""
        SELECT seq_md5, sequence, ro_cluster, has_genomic_context
        FROM ro_union
        WHERE passes_gate = 1 AND sequence IS NOT NULL AND LENGTH(sequence) > 0
    """).fetchall()
    if not rows:
        raise SystemExit("ro_union bos ya da yok -- once union adimi")

    meta = {r["seq_md5"]: (r["ro_cluster"], r["has_genomic_context"] or 0)
            for r in rows}
    print(f"[birlik] {len(rows):,} dizi kapiyi geciyor, kademe hesaplanacak")

    #  OLCUM evidence_tiers'in kendisinden: diamond cagrisi da, esikler de.
    best = evidence_tiers.best_hits(
        ((r["seq_md5"], r["sequence"]) for r in rows),
        args.refs, args.threads)

    out = []
    for md5, (cluster, ctx) in meta.items():
        sid, pid, qcov = best.get(md5, (None, None, 0.0))
        model = evidence_tiers.ref_model(sid)
        out.append((md5, model, pid, qcov,
                    int(model == cluster) if model and cluster else 0,
                    evidence_tiers.tier_of(pid, qcov), ctx))

    con.executescript(SCHEMA)
    con.executemany("INSERT INTO ro_union_evidence VALUES (?,?,?,?,?,?,?)", out)
    con.commit()

    #  Dagilim UC AYRI payda icin: tum birlik, baglamli, baglamsiz. Tek
    #  sayi vermek, kademe tablosunun hangi havuzu gosterdigi sorusunu
    #  yeniden belirsiz birakirdi -- bu betigin var olma sebebi o soru.
    def dist(predicate):
        c = Counter(r[5] for r in out if predicate(r))
        total = sum(c.values())
        return {"total": total,
                "tiers": {t: c[t] for t in evidence_tiers.TIER_ORDER}}

    by_cluster = defaultdict(lambda: Counter())
    for md5, model, pid, qcov, same, tier, ctx in out:
        cluster = meta[md5][0]
        if cluster:
            by_cluster[cluster][tier] += 1
            by_cluster[cluster]["_ctx" if ctx else "_noctx"] += 1

    report = {
        "module": "union_evidence.py",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "what_was_measured":
            "The same evidence tiers as the genomic table, computed for the "
            "union sequences: every gate-passing sequence in ro_union, "
            "whether or not it has genomic context. Thresholds and the "
            "DIAMOND call are imported from evidence_tiers.py, not restated.",
        "what_this_is_not":
            "These counts are NOT comparable with the genomic table's "
            "totals. The ro table counts ENTRIES, so one protein found in a "
            "hundred genomes is a hundred rows; ro_union counts DISTINCT "
            "SEQUENCES. Two counting rules cannot be divided by one another. "
            "A tier also says nothing about substrate: it measures "
            "confidence in the enzyme TYPE.",
        "pools": {
            "all_gate_passing": dist(lambda r: True),
            "with_genomic_context": dist(lambda r: r[6]),
            "without_genomic_context": dist(lambda r: not r[6]),
        },
        "by_cluster": {c: dict(v) for c, v in sorted(by_cluster.items())},
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    tmp = args.out + ".partial"
    with open(tmp, "w") as fh:
        json.dump(report, fh, indent=1, ensure_ascii=False)
    os.replace(tmp, args.out)

    for name, pool in report["pools"].items():
        line = ", ".join(f"{t} {pool['tiers'][t]:,}"
                         for t in evidence_tiers.TIER_ORDER)
        print(f"   {name:26s} n={pool['total']:>6,}  {line}")
    print("[yazildi]", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
