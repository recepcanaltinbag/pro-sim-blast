"""
Veri fazlaligi -- her yuzdenin altinda yatan sayim sorunu.

Genom veritabanlari ayni proteini cok kez icerir: ayni turun yuzlerce susu
dizilenmistir ve aralarinda ozdes RO alfa alt birimleri bulunur. Bu, "uyelerin
%X'i plazmit uzerinde" turu her ifadeyi etkiler, cunku X aslinda "dizilenmis
suslarin %X'i" demektir.

BU SCRIPT fazlaligi olcer ve her sayimin tekillestirilmis karsiligini verir:
    - ayni protein_id'yi (RefSeq WP_ numarasi) paylasan girisler
    - ayni AMINO ASIT DIZISINI paylasan girisler (protein_id'den bagimsiz)
    - kume boyutlarinin tekillestirme ile nasil degistigi
    - fazlaligin hangi cinslerde biriktigi

Girisler SILINMEZ: ayni protein farkli genomlarda farkli komsulukta bulunabilir
ve genomik baglam veritabaninin asil konusu oldugu icin her kopya ayri bir
gozlemdir. Olculen sey, bu kopyaların istatistige ne kadar etki ettigidir.

Cikti: analysis_out/redundancy.json
"""

import argparse
import json
import os
import sqlite3
from collections import Counter


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/redundancy.json")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    scalar = lambda sql: con.execute(sql).fetchone()[0]   # noqa: E731

    total = scalar("SELECT COUNT(*) FROM ro WHERE is_confirmed=1")
    distinct_seq = scalar("SELECT COUNT(DISTINCT sequence) FROM ro WHERE is_confirmed=1")
    dup_seq_entries = scalar("""
        SELECT COUNT(*) FROM ro WHERE is_confirmed=1 AND sequence IN
          (SELECT sequence FROM ro WHERE is_confirmed=1 GROUP BY 1 HAVING COUNT(*) > 1)""")
    distinct_pid = scalar("SELECT COUNT(DISTINCT protein_id) FROM ro "
                          "WHERE is_confirmed=1 AND protein_id != ''")
    dup_pid_entries = scalar("""
        SELECT COUNT(*) FROM ro WHERE is_confirmed=1 AND protein_id != '' AND protein_id IN
          (SELECT protein_id FROM ro WHERE is_confirmed=1 AND protein_id != ''
           GROUP BY 1 HAVING COUNT(*) > 1)""")
    max_copies = scalar("SELECT MAX(n) FROM (SELECT COUNT(*) n FROM ro "
                        "WHERE is_confirmed=1 GROUP BY sequence)")

    # Ayni dizinin kac FARKLI replikonda gorundugu: ayni genomda iki kopya
    # (gercek gen duplikasyonu) ile iki susta ayni gen farkli seylerdir.
    multi_replicon = scalar("""
        SELECT COUNT(*) FROM (
          SELECT sequence FROM ro WHERE is_confirmed=1
          GROUP BY sequence HAVING COUNT(DISTINCT nucleotide_id) > 1)""")
    same_replicon = scalar("""
        SELECT COUNT(*) FROM (
          SELECT sequence, nucleotide_id FROM ro WHERE is_confirmed=1
          GROUP BY sequence, nucleotide_id HAVING COUNT(*) > 1)""")

    clusters = [{"cluster": r[0], "entries": r[1], "unique_sequences": r[2],
                 "redundancy": round(1 - r[2] / r[1], 4)}
                for r in con.execute("""
                    SELECT ro_cluster, COUNT(*), COUNT(DISTINCT sequence) FROM ro
                    WHERE is_confirmed=1 GROUP BY 1 ORDER BY COUNT(*) DESC""")]

    genera = Counter()
    for (organism,) in con.execute("""
            SELECT p.organism FROM ro r JOIN replicon p USING(nucleotide_id)
            WHERE r.is_confirmed=1 AND r.sequence IN
              (SELECT sequence FROM ro WHERE is_confirmed=1 GROUP BY 1 HAVING COUNT(*) > 1)"""):
        genera[(organism or "?").split()[0]] += 1

    # Plazmit orani: her kopya tek tek sayildiginda ve dizi basina bir kez
    raw_plasmid = scalar("""
        SELECT COUNT(*) FROM ro r JOIN replicon p USING(nucleotide_id)
        WHERE r.is_confirmed=1 AND p.is_plasmid=1""")
    dedup_plasmid = scalar("""
        SELECT COUNT(*) FROM (
          SELECT r.sequence, MAX(p.is_plasmid) pl FROM ro r JOIN replicon p USING(nucleotide_id)
          WHERE r.is_confirmed=1 GROUP BY r.sequence) WHERE pl=1""")

    out = {
        "entries": total,
        "unique_sequences": distinct_seq,
        "redundancy": round(1 - distinct_seq / total, 4),
        "entries_with_a_duplicated_sequence": dup_seq_entries,
        "distinct_protein_ids": distinct_pid,
        "entries_sharing_a_protein_id": dup_pid_entries,
        "max_copies_of_one_sequence": max_copies,
        "sequences_in_more_than_one_replicon": multi_replicon,
        "sequences_twice_in_one_replicon": same_replicon,
        "plasmid_rate_per_entry": round(raw_plasmid / total, 4),
        "plasmid_rate_per_unique_sequence": round(dedup_plasmid / distinct_seq, 4),
        "clusters": clusters,
        "top_redundant_genera": genera.most_common(12),
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"[giris]            {total:,}")
    print(f"[tekil dizi]       {distinct_seq:,}  (fazlalik %{100 * out['redundancy']:.1f})")
    print(f"[kopyali giris]    {dup_seq_entries:,}; bir dizinin en fazla kopyasi {max_copies}")
    print(f"[protein_id]       {distinct_pid:,} tekil; {dup_pid_entries:,} giris paylasiyor")
    print(f"[ayni replikonda]  {same_replicon:,} dizi iki kez (gercek gen duplikasyonu olabilir)")
    print(f"[plazmit orani]    giris basina %{100 * out['plasmid_rate_per_entry']:.2f}, "
          f"tekil dizi basina %{100 * out['plasmid_rate_per_unique_sequence']:.2f}")
    worst = [c for c in clusters if c["entries"] >= 50][:5]
    print("[kume etkisi]      " + ", ".join(
        f"{c['cluster']} -%{100 * c['redundancy']:.0f}" for c in
        sorted(worst, key=lambda c: -c["redundancy"])))
    print(f"[yazildi] {args.out}")
    con.close()


if __name__ == "__main__":
    main()
