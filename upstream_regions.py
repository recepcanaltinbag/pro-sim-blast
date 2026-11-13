"""Operonlarin YUKARI AKIS intergenik DNA'sini cikarir.

NEDEN BU BETIK VAR. ROADMAP 124.3'un kalan yarisi sigma70 kutusu
zenginlesmesini eslenmis bir null'a karsi olcmeyi istiyor, ve o olcum
icin intergenik DIZI gerekiyor -- veritabaninda komsu genlerin protein
cevirileri var, DNA yok. Bu betik DNA'yi cikarir; zenginlesmeyi OLCMEZ.
Payda once (135.21), iddia sonra.

KENDI CIKTISINI DOGRULAR. Bolgenin uzunlugu zaten biliniyor:
`ro_regulation.intergenic_bp`, baska bir betik tarafindan koordinatlardan
hesaplanmis. Dolayisiyla cikarilan dizinin uzunlugu o sayiya ESIT OLMAK
ZORUNDA. Esit degilse koordinat okumam yanlis demektir, ve bu sessizce
gecmemeli: uyusmazlik sayilir ve bolge YAZILMAZ. Bu projede bir cikti
"gorunuyor diye dogru sayilamaz" dersi tekrar tekrar cikti.

IPLIK. Promotor, operonu SUREN ipliktedir. Operon eksi iplikteyse yukari
akis bolgesi ters tumleyen alinir; yoksa motif aramasi yanlis iplikte
yapilir ve "kutu yok" sonucu bir yon hatasindan gelir.

Cikti: analysis_out/upstream_regions.fasta.gz + upstream_regions.csv
"""

import argparse
import csv
import gzip
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_origin(gbk_dir, nucleotide_id):
    """GenBank kaydindan nukleotid dizisi. Dizi yoksa None.

    ORIGIN'den sonraki satirlarda konum numarasi ve bosluklu bloklar var;
    harf olmayan her sey atilir.
    """
    plain = os.path.join(gbk_dir, nucleotide_id + ".gbk")
    if os.path.exists(plain):
        handle = open(plain, "rt", encoding="utf-8", errors="replace")
    elif os.path.exists(plain + ".gz"):
        handle = gzip.open(plain + ".gz", "rt", encoding="utf-8",
                           errors="replace")
    else:
        return None
    chunks, seen = [], False
    with handle:
        for line in handle:
            if not seen:
                if line.startswith("ORIGIN"):
                    seen = True
                continue
            if line.startswith("//"):
                break
            chunks.append(re.sub(r"[^A-Za-z]", "", line))
    return "".join(chunks).upper() if seen else None


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def regions(con, readable):
    """Cikarilacak bolgeler: (cid, nuc, baslangic, bitis, strand, gap, arch)."""
    placeholders = ",".join("?" * len(readable))
    return con.execute(f"""
        SELECT r.candidate_id, r.nucleotide_id, g.operon_5p, n.start, n.end,
               o.strand, g.intergenic_bp, g.architecture
          FROM ro_regulation g
          JOIN ro r USING(candidate_id)
          JOIN operon o USING(candidate_id)
          JOIN neighbor n ON n.neighbor_id = g.upstream_gene_id
          JOIN replicon_kind k ON k.nucleotide_id = r.nucleotide_id
         WHERE r.is_confirmed = 1
           AND g.intergenic_bp IS NOT NULL
           AND k.record_kind = 'genome'
           AND r.nucleotide_id IN ({placeholders})
         ORDER BY r.nucleotide_id""", tuple(readable)).fetchall()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--gbk", default="gbk_files")
    parser.add_argument("--census",
                        default="analysis_out/replicon_sequence.csv")
    parser.add_argument("--out-fasta",
                        default="analysis_out/upstream_regions.fasta.gz")
    parser.add_argument("--out-csv",
                        default="analysis_out/upstream_regions.csv")
    args = parser.parse_args()

    with open(args.census, newline="", encoding="utf-8") as handle:
        readable = {row["nucleotide_id"] for row in csv.DictReader(handle)
                    if row["has_sequence"] == "1"}
    if not readable:
        sys.exit("dizisi okunabilir replikon yok -- sayim kosuldu mu?")

    con = sqlite3.connect(args.db)
    rows = regions(con, readable)
    print(f"[bolge] {len(rows)} aday bolge, "
          f"{len({r[1] for r in rows})} replikonda", flush=True)

    by_replicon = defaultdict(list)
    for row in rows:
        by_replicon[row[1]].append(row)

    tally = Counter()
    records, meta = [], []
    for n, (nuc, group) in enumerate(sorted(by_replicon.items()), 1):
        sequence = read_origin(args.gbk, nuc)
        if not sequence:
            tally["replicon_sequence_unreadable"] += len(group)
            continue
        for cid, _nuc, five_prime, up_start, up_end, strand, gap, arch in group:
            #  Bolge, yukari akis genin BITISI ile operonun 5' ucu
            #  arasindaki aralik. Hangisinin kucuk oldugu iplige gore
            #  degisir, o yuzden siralanir.
            low, high = sorted((up_end, five_prime)) if strand == 1 \
                else sorted((five_prime, up_start))
            low, high = max(0, low), min(len(sequence), high)
            piece = sequence[low:high]
            if len(piece) != gap:
                #  SESSIZCE YAZILMAZ: uzunluk zaten baska bir betik
                #  tarafindan hesaplanmis, uyusmazlik koordinat okumamin
                #  yanlis oldugunu soyler.
                tally["length_mismatch"] += 1
                continue
            if not piece:
                tally["empty"] += 1
                continue
            if strand != 1:
                piece = revcomp(piece)
            gc = sum(piece.count(base) for base in "GC")
            records.append((cid, piece))
            meta.append({"candidate_id": cid, "nucleotide_id": nuc,
                         "length": len(piece), "architecture": arch,
                         "gc_fraction": round(gc / len(piece), 4),
                         "strand": strand})
            tally["written"] += 1
        if n % 500 == 0:
            print(f"   {n}/{len(by_replicon)} replikon", flush=True)

    os.makedirs(os.path.dirname(args.out_fasta) or ".", exist_ok=True)
    with gzip.open(args.out_fasta, "wt", encoding="utf-8") as handle:
        for cid, piece in records:
            handle.write(f">{cid}\n{piece}\n")
    with open(args.out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "candidate_id", "nucleotide_id", "length", "architecture",
            "gc_fraction", "strand"])
        writer.writeheader()
        writer.writerows(meta)

    print(f"[sonuc] {dict(tally)}")
    print(f"[yazildi] {args.out_fasta}  ve  {args.out_csv}")
    #  Uyusmazlik VARSA bu bir hata: cikti yazilir ama kod 1 doner, cunku
    #  kismi bir bolge kumesi uzerinde olculen zenginlesme sessizce
    #  temsil edici gorunur.
    return 1 if tally["length_mismatch"] else 0


if __name__ == "__main__":
    sys.exit(main())
