"""
Yasam alani (domain) siniflandirmasi: DB'nin bakteriyel mi, karisik mi oldugunu olc.

Temsilci analizinde bazi kume "tip" dizilerinin BITKI proteini oldugu ortaya
cikti (DdmC, CndA -> Gossypium/pamuk). Bu tek tuk degil: dogrulanmis RO'larin
%8,3'u okaryot (bitki, mantar, hayvan). PF00355 seti bakteriyel halka-hidroksileyen
oksijenazlarla okaryot Rieske proteinlerini (kloroplast CMO/CAO gibi) karistiriyor.

Onemli cikarim: katalitik-triad filtresi ferredoksini iyi eliyor ama "bakteriyel
yikim enzimi"ni "Rieske+katalitik merkez tasiyan herhangi bir okaryot protein"den
AYIRMIYOR. Bu bir veri-kompozisyonu gercegi -- gizlenmemeli, isaretlenmeli.

Bu modul her uyeyi domain'e (Bacteria/Eukaryota/Archaea) ve okaryotlari alt-gruba
(Viridiplantae=bitki/yesil alg, Fungi, Metazoa, diger) ayirir; kume bazinda
kirilim cikarir; okaryot-agirlikli kumeleri isaretler.

Cikti:
    roar.sqlite   yeni tablo: ro_domain
    analysis_out/domain_by_cluster.csv
"""

import argparse
import csv
import os
import sqlite3
from collections import Counter, defaultdict

# Okaryot alt-gruplarini kaba grupla. Degerler INGILIZCE, cunku bu sutun
# dogrudan web arayuzunde gosteriliyor: etiketi burada Turkce yazmak onu
# sayfaya tasiyor. Arayuz dili ile kod yorumlarinin dili ayri seylerdir.
EUK_GROUPS = {
    "Viridiplantae": "plant/alga",
    "Fungi": "fungus",
    "Metazoa": "animal",
    "Rhodophyta": "red alga",
    "Sar": "SAR",
    "Haptista": "other eukaryote",
    "Amoebozoa": "other eukaryote",
}


def classify(taxonomy):
    """Donen: (domain, alt_grup). taxonomy bos ise ('?', '?')."""
    if not taxonomy:
        return "?", "?"
    parts = [p.strip() for p in taxonomy.split(";")]
    domain = parts[0]
    if domain == "Eukaryota":
        second = parts[1] if len(parts) > 1 else "?"
        return "Eukaryota", EUK_GROUPS.get(second, "other eukaryote")
    return domain, domain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default="analysis_out")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    connection = sqlite3.connect(args.db)
    connection.row_factory = sqlite3.Row

    connection.executescript("""
        DROP TABLE IF EXISTS ro_domain;
        CREATE TABLE ro_domain (
            candidate_id TEXT PRIMARY KEY,
            cluster      TEXT,
            domain       TEXT,   -- Bacteria | Eukaryota | Archaea | Viruses
            euk_group    TEXT    -- okaryotsa plant/alga, fungus, animal...
        );
        CREATE INDEX idx_rodom_cluster ON ro_domain(cluster);
        CREATE INDEX idx_rodom_domain  ON ro_domain(domain);
    """)

    rows = connection.execute("""
        SELECT r.candidate_id, r.ro_cluster, rep.taxonomy
        FROM ro r JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        WHERE r.is_confirmed = 1 AND r.ro_cluster IS NOT NULL AND r.ro_cluster != 'N/A'
    """).fetchall()

    inserts, cluster_domain = [], defaultdict(Counter)
    for row in rows:
        domain, group = classify(row["taxonomy"])
        inserts.append((row["candidate_id"], row["ro_cluster"], domain, group))
        cluster_domain[row["ro_cluster"]][domain] += 1

    connection.executemany("INSERT INTO ro_domain VALUES (?,?,?,?)", inserts)
    connection.commit()

    # --- Genel dagilim
    overall = Counter(i[2] for i in inserts)
    total = sum(overall.values())
    print("=" * 66)
    print("DOGRULANMIS RO -- YASAM ALANI DAGILIMI")
    print("=" * 66)
    for domain, count in overall.most_common():
        print(f"  {domain:12s} {count:>6,}  (%{100*count/total:.1f})")

    euk_groups = Counter(i[3] for i in inserts if i[2] == "Eukaryota")
    print("\n  Okaryotlarin alt-grubu:")
    for group, count in euk_groups.most_common():
        print(f"    {group:14s} {count:>5,}")

    # --- Okaryot-agirlikli kumeler
    print("\n" + "=" * 66)
    print("OKARYOT-AGIRLIKLI KUMELER (bakteriyel yikim DB'sine ait olmayabilir)")
    print("=" * 66)
    print(f"  {'kume':16s} {'toplam':>7s} {'okaryot':>8s} {'%':>6s}  baskin grup")
    print("  " + "-" * 58)
    flagged = []
    for cluster, counts in cluster_domain.items():
        n = sum(counts.values())
        euk = counts.get("Eukaryota", 0)
        if n >= 10 and euk / n >= 0.25:
            groups = Counter(i[3] for i in inserts
                             if i[1] == cluster and i[2] == "Eukaryota")
            top = groups.most_common(1)[0][0] if groups else "?"
            flagged.append((cluster, n, euk, euk / n, top))
    for cluster, n, euk, frac, top in sorted(flagged, key=lambda x: -x[3]):
        print(f"  {cluster:16s} {n:>7,} {euk:>8,} {frac:>5.0%}  {top}")
    if not flagged:
        print("  (yok)")

    # --- CSV
    path = os.path.join(args.out_dir, "domain_by_cluster.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cluster", "total", "bacteria", "eukaryota", "archaea",
                         "eukaryota_rate", "top_euk_group"])
        for cluster, counts in sorted(cluster_domain.items()):
            n = sum(counts.values())
            euk = counts.get("Eukaryota", 0)
            groups = Counter(i[3] for i in inserts
                             if i[1] == cluster and i[2] == "Eukaryota")
            top = groups.most_common(1)[0][0] if groups else ""
            writer.writerow([cluster, n, counts.get("Bacteria", 0), euk,
                             counts.get("Archaea", 0), round(euk / n, 4) if n else 0, top])
    print(f"\n[yazildi] {path}")
    print(f"[bilgi] {len(flagged)} kume okaryot-agirlikli (>=%25 okaryot, n>=10)")
    connection.close()


if __name__ == "__main__":
    main()
