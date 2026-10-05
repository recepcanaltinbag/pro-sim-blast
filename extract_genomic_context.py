"""
GenBank dosyalarindan RO alpha'lari ve genomik baglamlarini cikar.

NEDEN BOYLE: Mevcut pipeline RO'nun konumunu tblastn ile buluyor
(mongo_add_nuc_loc_and_related_genes.py) ve sadece start/end kaydediyor.
Bunun uc sorunu var:
  1. STRAND kaydedilmiyor -- orientation analizi imkansiz hale geliyor.
  2. tblastn koordinati anotasyondaki CDS sinirlariyla birebir ortusmeyebilir.
  3. Konum mongo'da; taksonomi baska yerde; komsu genler ucuncu yerde.

Bunun yerine RO'yu dogrudan ANOTASYONUN ICINDE buluyoruz. Boylece koordinat,
strand, komsular, locus_tag ve organizma bilgisi tek kaynaktan, tutarli gelir.

UC ASAMALI TARAMA (17.074 gbk, ~10M CDS icin gerekli):
  1. Tum CDS'leri cek (translation + koordinat + strand + product)
  2. UCUZ ON FILTRE: Rieske ligand motifi regex'i  C-x-H-x(14,18)-C-x-x-H
     Milyonlarca CDS'i binlere indirir, HMM'e sadece adaylar gider.
  3. Hayatta kalanlari hmmsearch + motif testinden gecir

KOMSU KAYDI: filtreleme YAPILMAZ. Mevcut kod is_matching_product() ile sadece
anahtar kelimeye uyan genleri sakliyor -- "hypothetical protein" komsular
tamamen kayboluyor. Oysa README'nin hedefi tam da bilinmeyen genlere fonksiyon
atfetmek. Ham komsu saklanir, siniflandirma sorgu aninda yapilir.
"""

import argparse
import csv
import gzip
import os
import re
import sys
from multiprocessing import Pool

from Bio import SeqIO

# Rieske [2Fe-2S] ligand imzasi: C-x-H ... C-x-x-H
# Referanslarda olculen aralik 14-18 kalinti; guvenli tarafta genis tutuldu.
RIESKE_REGEX = re.compile(r"C.H.{12,22}C..H")

# RO alpha uzunluk taban degeri -- fragmentleri ucuz asamada eler.
MIN_CDS_LENGTH = 250

NEIGHBOR_WINDOW = 10000   # RO'nun her iki yaninda kac bp taranacak


def feature_spans(feature):
    """Bir feature'in gercek parcalarini dondur: [(start, end), ...]

    NEDEN GEREKLI: dairesel replikonlarda orijini asan genler GenBank'ta
    join(...) ile yazilir, ornegin  join{[116527:116580](+), [0:355](+)}.
    BioPython'un feature.location.start / .end degerleri bu durumda parcalarin
    MIN ve MAX'ini verir -- yani 408 bp'lik bir gen 116.580 bp'lik bir aralik
    gibi gorunur ve replikondaki HER genin komsusu olur.

    Mevcut mongo_gbk_analysis.py bu tuzaga dusuyor (int(feature.location.start)).
    Plazmitler cogunlukla dairesel ve mobil RO kumelerinin bulundugu yer
    oldugu icin hata tam da ilgilendigimiz bolgede birikiyor.
    """
    if len(feature.location.parts) > 1:
        return [(int(part.start), int(part.end)) for part in feature.location.parts]
    return [(int(feature.location.start), int(feature.location.end))]


def signed_distance(ro_start, ro_end, spans, replicon_length=0, is_circular=False):
    """RO ile bir feature arasindaki isaretli en kisa mesafe.

    Negatif = RO'nun solunda, pozitif = saginda, 0 = ortusuyor.
    Dairesel replikonda orijin uzerinden gecen kisa yol da degerlendirilir.
    """
    candidates = []
    for start, end in spans:
        if end <= ro_start:
            candidates.append(end - ro_start)          # negatif
        elif start >= ro_end:
            candidates.append(start - ro_end)          # pozitif
        else:
            return 0                                    # ortusuyor

        if is_circular and replicon_length:
            # Orijin uzerinden dolasan alternatif mesafe
            if end <= ro_start:
                candidates.append(replicon_length - ro_end + start)
            elif start >= ro_end:
                candidates.append(-(replicon_length - end + ro_start))

    return min(candidates, key=abs) if candidates else 0


def wrap_offset(offset, n_genes, is_circular):
    """Dairesel replikonda gen sirasi farkini en kisa yone sar."""
    if not is_circular or n_genes < 2:
        return offset
    if offset > n_genes // 2:
        return offset - n_genes
    if offset < -(n_genes // 2):
        return offset + n_genes
    return offset


def parse_genbank(path):
    """Bir gbk dosyasindan tum CDS'leri ve replicon meta verisini cikar.

    Donen: (replicon_dict, [cds_dict, ...]) veya anotasyon yoksa (meta, [])
    """
    records = []
    try:
        handle = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
        with handle:
            records = list(SeqIO.parse(handle, "genbank"))
    except Exception as exc:
        return {"file": os.path.basename(path), "error": str(exc)}, []

    all_cds = []
    replicon = None

    for record in records:
        organism, taxonomy, mol_type, is_plasmid = "", "", "", False
        for feature in record.features:
            if feature.type == "source":
                organism = feature.qualifiers.get("organism", [""])[0]
                mol_type = feature.qualifiers.get("mol_type", [""])[0]
                if "plasmid" in feature.qualifiers:
                    is_plasmid = True
                break
        if record.annotations.get("taxonomy"):
            taxonomy = "; ".join(record.annotations["taxonomy"])
        if "plasmid" in record.description.lower():
            is_plasmid = True

        if replicon is None:
            replicon = {
                "nucleotide_id": record.id,
                "file": os.path.basename(path),
                "organism": organism,
                "taxonomy": taxonomy,
                "description": record.description,
                "length": len(record.seq),
                "mol_type": mol_type,
                "is_plasmid": is_plasmid,
            }

        # SADECE CDS -- "gene" feature'lari da almak her geni iki kez sayar.
        # Mevcut mongo_gbk_analysis.py bu hatayi yapiyor (feature.type in
        # ["gene","CDS"]); prokaryot gbk'da gene ve CDS sayilari esittir.
        is_circular = record.annotations.get("topology", "").lower() == "circular"
        if replicon is not None and replicon["nucleotide_id"] == record.id:
            replicon["is_circular"] = is_circular

        for feature in record.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation", [""])[0]
            spans = feature_spans(feature)
            all_cds.append({
                "nucleotide_id": record.id,
                "spans": spans,
                "spans_origin": len(spans) > 1,
                "replicon_length": len(record.seq),
                "is_circular": is_circular,
                # start/end raporlama icin; mesafe hesabi spans uzerinden yapilir
                "start": min(s for s, _ in spans),
                "end": max(e for _, e in spans),
                "strand": feature.location.strand or 0,
                "product": feature.qualifiers.get("product", ["hypothetical protein"])[0],
                "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
                "protein_id": feature.qualifiers.get("protein_id", [""])[0],
                "gene": feature.qualifiers.get("gene", [""])[0],
                "translation": translation,
            })

    if replicon is None:
        replicon = {"nucleotide_id": "", "file": os.path.basename(path),
                    "organism": "", "taxonomy": "", "description": "",
                    "length": 0, "mol_type": "", "is_plasmid": False}
    replicon["cds_count"] = len(all_cds)
    return replicon, all_cds


def find_ro_candidates(all_cds):
    """Ucuz on filtre: Rieske ligand motifi tasiyan, yeterince uzun CDS'ler."""
    candidates = []
    for index, cds in enumerate(all_cds):
        translation = cds["translation"]
        if len(translation) < MIN_CDS_LENGTH:
            continue
        if RIESKE_REGEX.search(translation):
            candidates.append(index)
    return candidates


def process_file(path):
    """Tek bir gbk dosyasini isle. Donen: (replicon, candidates, neighbors)"""
    replicon, all_cds = parse_genbank(path)
    if not all_cds:
        # Anotasyonsuz dosya. Sessizce bos donmek yerine isaretlenir --
        # orneklemde dosyalarin %17'si sifir CDS iceriyor ve bunlar
        # "komsusuz RO" gibi degil, "verisi eksik" olarak sayilmali.
        replicon["status"] = "no_annotation"
        return replicon, [], []

    replicon["status"] = "ok"
    candidate_indices = find_ro_candidates(all_cds)
    if not candidate_indices:
        return replicon, [], []

    candidates, neighbors = [], []
    for index in candidate_indices:
        cds = all_cds[index]
        candidate_id = f"{cds['nucleotide_id']}:{cds['start']}-{cds['end']}:{cds['strand']}"
        entry = dict(cds)
        entry["candidate_id"] = candidate_id
        candidates.append(entry)

        # Komsulari yakala -- HIC FILTRE YOK, hypothetical'lar dahil
        n_on_replicon = sum(1 for c in all_cds if c["nucleotide_id"] == cds["nucleotide_id"])
        for other_index, other in enumerate(all_cds):
            if other_index == index:
                continue
            # Cok kayitli dosyada (gbff) baska kontigin geni komsu sayilmasin
            if other["nucleotide_id"] != cds["nucleotide_id"]:
                continue
            distance = signed_distance(
                cds["start"], cds["end"], other["spans"],
                other["replicon_length"], other["is_circular"])
            if abs(distance) > NEIGHBOR_WINDOW:
                continue
            neighbors.append({
                "candidate_id": candidate_id,
                "nucleotide_id": other["nucleotide_id"],
                "start": other["start"],
                "end": other["end"],
                "strand": other["strand"],
                "distance": distance,
                # RO ile ayni yonde mi -- operon ve divergent regulator
                # analizinin temel sinyali. Mevcut pipeline'da hic yok.
                "same_strand": other["strand"] == cds["strand"],
                # kac gen uzakta; dairesel replikonda orijin uzerinden gecen
                # kisa yol icin sarilir (+3000 yerine -3 gibi)
                "gene_offset": wrap_offset(other_index - index, n_on_replicon,
                                           other["is_circular"]),
                "translation": other["translation"],
                "spans_origin": other["spans_origin"],
                "product": other["product"],
                "locus_tag": other["locus_tag"],
                "protein_id": other["protein_id"],
                "gene": other["gene"],
            })
    return replicon, candidates, neighbors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gbk-dir", default="gbk_files")
    parser.add_argument("--out-dir", default="genomic_context")
    parser.add_argument("--processes", type=int, default=6)
    parser.add_argument("--limit", type=int, default=0, help="test icin dosya sayisi")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    files = sorted(f for f in os.listdir(args.gbk_dir)
                   if f.endswith((".gbk", ".gb", ".gbk.gz", ".gbff")))
    if args.limit:
        files = files[:args.limit]
    paths = [os.path.join(args.gbk_dir, f) for f in files]
    print(f"[bilgi] {len(paths)} gbk dosyasi, {args.processes} surec")

    replicon_path = os.path.join(args.out_dir, "replicons.csv")
    candidate_path = os.path.join(args.out_dir, "ro_candidates.csv")
    neighbor_path = os.path.join(args.out_dir, "neighbors.csv")
    fasta_path = os.path.join(args.out_dir, "ro_candidates.fasta")
    nb_fasta_path = os.path.join(args.out_dir, "neighbor_proteins.fasta")

    replicon_fields = ["nucleotide_id", "file", "organism", "taxonomy", "description",
                       "length", "mol_type", "is_plasmid", "is_circular", "cds_count", "status"]
    candidate_fields = ["candidate_id", "nucleotide_id", "start", "end", "strand",
                        "product", "locus_tag", "protein_id", "gene"]
    neighbor_fields = ["candidate_id", "nucleotide_id", "start", "end", "strand",
                       "distance", "same_strand", "gene_offset", "spans_origin",
                       "product", "locus_tag", "protein_id", "gene"]

    counts = {"files": 0, "no_annotation": 0, "candidates": 0, "neighbors": 0}

    with open(replicon_path, "w", newline="") as rep_fh, \
         open(candidate_path, "w", newline="") as cand_fh, \
         open(neighbor_path, "w", newline="") as nb_fh, \
         open(fasta_path, "w") as fa_fh, \
         open(nb_fasta_path, "w") as nb_fa_fh:

        rep_writer = csv.DictWriter(rep_fh, fieldnames=replicon_fields, extrasaction="ignore")
        cand_writer = csv.DictWriter(cand_fh, fieldnames=candidate_fields, extrasaction="ignore")
        nb_writer = csv.DictWriter(nb_fh, fieldnames=neighbor_fields, extrasaction="ignore")
        rep_writer.writeheader(); cand_writer.writeheader(); nb_writer.writeheader()

        with Pool(args.processes) as pool:
            for replicon, candidates, neighbors in pool.imap_unordered(
                    process_file, paths, chunksize=4):
                counts["files"] += 1
                if replicon.get("status") == "no_annotation":
                    counts["no_annotation"] += 1
                rep_writer.writerow(replicon)
                for candidate in candidates:
                    cand_writer.writerow(candidate)
                    fa_fh.write(f">{candidate['candidate_id']}\n{candidate['translation']}\n")
                nb_writer.writerows(neighbors)
                # Komsu protein dizileri (operon bileseni tespiti icin, build_operons.py)
                seen = set()
                for nb in neighbors:
                    key = f"{nb['nucleotide_id']}:{nb['start']}-{nb['end']}:{nb['strand']}"
                    if nb["translation"] and key not in seen:
                        seen.add(key)
                        nb_fa_fh.write(f">{key}\n{nb['translation']}\n")
                counts["candidates"] += len(candidates)
                counts["neighbors"] += len(neighbors)

                if counts["files"] % 500 == 0:
                    sys.stdout.write(
                        f"\r  {counts['files']}/{len(paths)} dosya | "
                        f"{counts['candidates']} aday | {counts['neighbors']} komsu | "
                        f"{counts['no_annotation']} anotasyonsuz")
                    sys.stdout.flush()

    print(f"\n[bitti] {counts['files']} dosya islendi")
    print(f"        anotasyonsuz    : {counts['no_annotation']} "
          f"({100*counts['no_annotation']/max(1,counts['files']):.1f}%)")
    print(f"        RO adayi        : {counts['candidates']}")
    print(f"        komsu kaydi     : {counts['neighbors']}")
    print(f"[yazildi] {replicon_path}, {candidate_path}, {neighbor_path}, {fasta_path}")
    print("\nSonraki adim: adaylari HMM + motif testinden gecir")
    print(f"  hmmsearch --domtblout {args.out_dir}/cand_dom.out --noali --cpu 22 \\")
    print(f"            ROs_71_Clean/ROmotif71.hmm {fasta_path}")


if __name__ == "__main__":
    main()
