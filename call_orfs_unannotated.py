"""
Anotasyonsuz ama nukleotidi olan kayitlar icin ORF cagrimi.

SORUN (kullanici tespiti, 2026-11-04). Bir kaydin CDS anotasyonu yoksa
`extract_genomic_context.py` ona `no_annotation` der ve oradaki dizi
genomik baglam KAZANMAZ. Ama anotasyon yokken NUKLEOTID VAR olabilir --
ve o zaman "genomik baglami yok" bir BIYOLOJI ifadesi degil, bir ARAC
ifadesidir. Ikisi ayni cumlede birlesmemeli.

OLCUM: 47 anotasyonsuz replikonun 46'sinda dizi var (255,5 Mb) ve bu
kayitlara 50 hedef dizi dusuyor. Kayitlar GERCEKTEN bos: ornekleme ile
bakildiginda yalnizca `source` ozelligi var, yani ayristiricimiz bir sey
dusurmuyor -- anotasyon NCBI tarafinda yok.

NE YAPIYOR. pyrodigal ile ORF cagirir ve cagrilan CDS'leri GenBank olarak
yazar. Sonra bu dizin AYNEN `extract_genomic_context.py`e verilir. Boylece
aday secimi, Rieske on filtresi ve komsuluk penceresi KENDI kurallariyla
calisir; burada hicbiri yeniden uygulanmaz. (Projenin tekrar eden hata
bicimi, karar veren kuralin iddia edilen kural olmamasidir; ikinci bir
cikarim yolu yazmak tam o hatadir.)

NE YAPMIYOR. Kapiyi gevsetmiyor. Cagrilan bir ORF'un cevirisi hedef diziyle
`a5_match_guard.py` kuralina gore BIREBIR tutmak zorunda; tutmazsa baglam
yine kurulmaz. ORF cagrimi bir bahane degil, yalnizca bakilmamis bir yere
bakmaktir.

PROVENANS. Cikti AYRI bir dizine yazilir ve her CDS `/note="pyrodigal"`
tasir, cunku cagrilmis bir komsulugun gen sinirlari OLCULMUS degil TAHMIN
edilmistir. Konuma bagli her olcum (operon, komsuluk) kaynagina gore
ayrilabilir olmali.

Kullanim:
  python3 call_orfs_unannotated.py \\
      --replicons /media/lin-bio/back2/roar_gc_a5/replicons.csv \\
      --replicons /media/lin-bio/back2/roar_gc_wide/replicons.csv \\
      --gbk-dir /media/lin-bio/back2/roar_gbk_a5_new --gbk-dir gbk_files \\
      --out-dir /media/lin-bio/back2/roar_gbk_prodigal
"""

import argparse
import csv
import gzip
import os
import sys

MIN_SEQ_LEN = 200     # altinda ORF cagrimi anlamsiz
META_BELOW = 100_000   # bu uzunlugun altinda meta modu (tek kayittan model egitilemez)


def unannotated(paths):
    """Cikarim replicons.csv dosyalarindan `no_annotation` satirlarini toplar."""
    seen, out = set(), []
    for path in paths:
        if not os.path.exists(path):
            print(f"[atlandi] {path} yok", file=sys.stderr)
            continue
        for row in csv.DictReader(open(path)):
            if row.get("status") != "no_annotation":
                continue
            key = row["nucleotide_id"]
            if key and key not in seen:
                seen.add(key)
                out.append(row)
    return out


def locate(filename, dirs):
    for d in dirs:
        p = os.path.join(d, filename)
        if os.path.exists(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--replicons", action="append", required=True,
                    help="cikarim replicons.csv; birden fazla verilebilir")
    ap.add_argument("--gbk-dir", action="append", required=True,
                    help="kaynak gbk dizini; birden fazla verilebilir")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    import pyrodigal
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    rows = unannotated(args.replicons)
    print(f"[hedef] {len(rows)} anotasyonsuz replikon")
    os.makedirs(args.out_dir, exist_ok=True)

    wrote = skipped = no_seq = 0
    total_orfs = 0
    for row in rows:
        src = locate(row["file"], args.gbk_dir)
        if not src:
            skipped += 1
            continue
        try:
            handle = gzip.open(src, "rt") if src.endswith(".gz") else open(src)
            records = list(SeqIO.parse(handle, "genbank"))
        except Exception as exc:
            print(f"[hata] {row['nucleotide_id']}: {exc}", file=sys.stderr)
            skipped += 1
            continue

        out_records, n_orfs = [], 0
        for rec in records:
            seq = str(rec.seq).upper()
            if len(seq.replace("N", "")) < MIN_SEQ_LEN:
                continue
            # meta modu: tek bir kisa kayittan gen modeli egitilemez, o yuzden
            # on-egitilmis modeller kullanilir. Uzun kayitlarda kendi modelini
            # egitmek daha dogru sonuc verir.
            finder = pyrodigal.GeneFinder(meta=len(seq) < META_BELOW)
            if len(seq) >= META_BELOW:
                try:
                    finder.train(seq)
                except Exception:
                    finder = pyrodigal.GeneFinder(meta=True)
            genes = finder.find_genes(seq)

            new = rec[:]
            new.features = [f for f in rec.features if f.type == "source"]
            for i, gene in enumerate(genes, 1):
                loc = SimpleLocation(gene.begin - 1, gene.end,
                                     strand=1 if gene.strand > 0 else -1)
                new.features.append(SeqFeature(loc, type="CDS", qualifiers={
                    "locus_tag": [f"PRODIGAL_{rec.id.replace('.', '_')}_{i:05d}"],
                    "product": ["hypothetical protein"],
                    # PROVENANS: her CDS kendi kaynagini tasir.
                    "note": ["pyrodigal"],
                    "translation": [gene.translate().rstrip("*")],
                    "codon_start": ["1"],
                    "transl_table": ["11"],
                }))
                n_orfs += 1
            out_records.append(new)

        if not out_records:
            no_seq += 1
            continue
        dest = os.path.join(args.out_dir, row["file"])
        if not dest.endswith(".gz"):
            dest += ".gz"
        with gzip.open(dest, "wt") as fh:
            SeqIO.write(out_records, fh, "genbank")
        wrote += 1
        total_orfs += n_orfs
        print(f"  {row['nucleotide_id']:<22} {n_orfs:>7,} ORF")

    print(f"\n[bitti] {wrote} dosya yazildi, {total_orfs:,} ORF cagrildi")
    print(f"        dizisi yok / kisa : {no_seq}")
    print(f"        kaynak bulunamadi : {skipped}")
    print(f"[dizin] {args.out_dir}")
    print("\nSonraki adim -- AYNI cikarim kurali bu dizine uygulanir:")
    print(f"  python3 extract_genomic_context.py --gbk-dir {args.out_dir} \\")
    print("      --out-dir /media/lin-bio/back2/roar_gc_prodigal")
    print("Ardindan a5_match_guard.py: cagrilan ORF'un cevirisi hedef diziyle")
    print("BIREBIR tutmazsa baglam kurulmaz.")


if __name__ == "__main__":
    main()
