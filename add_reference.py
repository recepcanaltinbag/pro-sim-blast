"""
Literaturden yeni bir referans RO ekle -- tek komut, dogrulamali.

NEDEN: referans seti kapali bir liste degil. Literaturde karakterize edilmis ve
bu sette olmayan RO'lar var; birini eklemek su dosyalarin HEPSINI tutarli
tutmayi gerektiriyordu:
    ROs_71_Clean/refs71.fasta   dizi
    cluster_ecology.csv         substrat sinifi, ekoloji notu, guven
    chemistry.csv               substrat adi, SMILES, urun, reaksiyon, PDB, kaynak
    ROs_71_Clean/*.hmm          modeller
    ROs_71_Clean/motif_columns.json   motif kolonlari
Elle yapildiginda biri atlanirsa hata sessiz kaliyordu: ornegin chemistry.csv'ye
eklenmeyen bir tip web arayuzunde "substrat yok" diye gorunuyor ama hicbir
dogrulama uyarmiyordu.

BU SCRIPT eklemeyi atomik hale getirir: once her seyi DOGRULAR, sonra yazar.
Yazdiktan sonra ne yapilmasi gerektigini de soyler (model yeniden kurulumu ve
pipeline'in hangi adimindan itibaren tekrar kosulmasi gerektigi).

Kullanim:
    python3 add_reference.py --id 2_219_NewA --fasta yeni.fasta \\
        --substrate "4-chlorobiphenyl" --smiles "Clc1ccc(cc1)-c1ccccc1" \\
        --product "cis-dihydrodiol" --reaction "cis-dihydroxylation" \\
        --reaction-class cis_dihydroxylation --family biaryls_ethers \\
        --substrate-class xenobiotic --confidence high \\
        --source "PMID:12345678" --source-kind paper \\
        --note "Entry enzyme of the PCB pathway in strain X"

    python3 add_reference.py --list-vocabularies      # izinli degerleri yazdir

ID BICIMI: <grup>_<numara>_<GenAdi>. Grup 1-5 arasi, numara grup icinde tekil.
Mevcut en buyuk numarayi gormek icin --list-vocabularies kullan.
"""

import argparse
import csv
import os
import re
import shutil
import sys
from datetime import date

REFS_FASTA = "ROs_71_Clean/refs71.fasta"
ECOLOGY_CSV = "cluster_ecology.csv"
CHEMISTRY_CSV = "chemistry.csv"

VALID_SUBSTRATE_CLASS = ["xenobiotic", "natural_aromatic", "natural_specialized", "unknown"]
VALID_CONFIDENCE = ["low", "medium", "high"]
VALID_REACTION_CLASS = ["cis_dihydroxylation", "angular_dioxygenation",
                        "dioxygenation_with_release", "O_demethylation", "N_demethylation",
                        "hydroxylation", "C_N_cleavage", "unknown"]
VALID_FAMILY = ["alkylbenzenes", "pah", "biaryls_ethers", "nitroaromatics", "haloaromatics",
                "sulfoaromatics", "aromatic_acids", "anilines", "quaternary_amines",
                "alkaloids", "terpenoids_steroids", "unknown"]
VALID_SOURCE_KIND = ["paper", "structure", "curator", "paper+structure", "structure+paper",
                     "paper+curator", "general", "none"]

ID_PATTERN = re.compile(r"^([1-5])_(\d{3})_([A-Za-z0-9]+)$")
SMILES_OK = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)=#$%:/\\.*]+$")
AMINO = set("ACDEFGHIKLMNPQRSTVWYXBZUO")


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def smiles_plausible(smiles):
    """Kaba ama yararli SMILES kontrolu: RDKit yok, yazim hatasini yakalamak yeter."""
    if not smiles:
        return True, ""
    if not SMILES_OK.match(smiles):
        return False, "izinli olmayan karakter var"
    for open_ch, close_ch in (("(", ")"), ("[", "]")):
        if smiles.count(open_ch) != smiles.count(close_ch):
            return False, f"{open_ch}{close_ch} dengesiz"
    if not re.search(r"[A-Za-z]", smiles):
        return False, "hic atom yok"
    return True, ""


def list_vocabularies():
    existing = read_fasta(REFS_FASTA) if os.path.exists(REFS_FASTA) else {}
    ids = sorted(h.split("_")[0] + "_" + h.split("_")[1] for h in existing)
    per_group = {}
    for i in ids:
        g, n = i.split("_")
        per_group.setdefault(g, []).append(int(n))
    print(f"Mevcut referans sayisi: {len(existing)}")
    for g in sorted(per_group):
        nums = sorted(per_group[g])
        print(f"  grup {g}: {len(nums)} referans, en buyuk numara {nums[-1]}, "
              f"sonraki musait {nums[-1] + 1}")
    print("\nsubstrate_class :", ", ".join(VALID_SUBSTRATE_CLASS))
    print("confidence      :", ", ".join(VALID_CONFIDENCE))
    print("reaction_class  :", ", ".join(VALID_REACTION_CLASS))
    print("family          :", ", ".join(VALID_FAMILY))
    print("source_kind     :", ", ".join(VALID_SOURCE_KIND))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list-vocabularies", action="store_true")
    ap.add_argument("--id", help="<grup>_<numara>_<GenAdi>, ornek 2_219_NewA")
    ap.add_argument("--fasta", help="tek dizi iceren FASTA dosyasi")
    ap.add_argument("--sequence", help="dizi dogrudan (FASTA yerine)")
    ap.add_argument("--substrate", help="substrat adi, INGILIZCE")
    ap.add_argument("--smiles", default="")
    ap.add_argument("--product", default="")
    ap.add_argument("--reaction", default="")
    ap.add_argument("--reaction-class", default="unknown", choices=VALID_REACTION_CLASS)
    ap.add_argument("--family", default="unknown", choices=VALID_FAMILY)
    ap.add_argument("--substrate-class", default="unknown", choices=VALID_SUBSTRATE_CLASS)
    ap.add_argument("--confidence", default="medium", choices=VALID_CONFIDENCE)
    ap.add_argument("--pdb", default="")
    ap.add_argument("--source", default="")
    ap.add_argument("--source-kind", default="paper", choices=VALID_SOURCE_KIND)
    ap.add_argument("--note", default="")
    ap.add_argument("--description", default="", help="FASTA basligindaki serbest metin")
    ap.add_argument("--dry-run", action="store_true", help="dogrula, yazma")
    args = ap.parse_args()

    if args.list_vocabularies:
        list_vocabularies()
        return 0

    problems = []
    if not args.id or not ID_PATTERN.match(args.id):
        problems.append("--id <grup 1-5>_<uc haneli numara>_<GenAdi> biciminde olmali")
    if not args.substrate:
        problems.append("--substrate zorunlu (bilinmiyorsa 'unknown' yaz)")
    if not args.source and args.source_kind != "none":
        problems.append("--source zorunlu: bu veri nereden geliyor? (PMID, DOI, PDB, curator)")
    if args.pdb and not re.match(r"^[0-9A-Za-z]{4}$", args.pdb):
        problems.append("--pdb tam dort alfanumerik karakter olmali")
    ok, why = smiles_plausible(args.smiles)
    if not ok:
        problems.append(f"--smiles gecersiz gorunuyor: {why}")

    sequence = (args.sequence or "").strip().upper()
    if args.fasta:
        if not os.path.exists(args.fasta):
            problems.append(f"--fasta bulunamadi: {args.fasta}")
        else:
            entries = read_fasta(args.fasta)
            if len(entries) != 1:
                problems.append(f"--fasta tam bir dizi icermeli, {len(entries)} bulundu")
            else:
                sequence = next(iter(entries.values())).upper()
    if not sequence:
        problems.append("--fasta ya da --sequence ile dizi vermelisin")
    else:
        illegal = sorted(set(sequence) - AMINO)
        if illegal:
            problems.append(f"dizide amino asit olmayan karakter: {''.join(illegal)}")
        if len(sequence) < 300:
            problems.append(f"dizi {len(sequence)} kalinti: RO alpha icin kisa (>=300 beklenir)")

    existing = read_fasta(REFS_FASTA) if os.path.exists(REFS_FASTA) else {}
    existing_ids = {h.rsplit("_", 1)[0] if h.count("_") > 2 else h for h in existing}
    short_ids = {"_".join(h.split("_")[:3]) for h in existing}
    if args.id in short_ids:
        problems.append(f"--id zaten var: {args.id}")
    for header, seq in existing.items():
        if seq.upper() == sequence:
            problems.append(f"bu dizi zaten var: {header}")
            break

    eco = {r["cluster"] for r in csv.DictReader(open(ECOLOGY_CSV))} \
        if os.path.exists(ECOLOGY_CSV) else set()
    chem = {r["cluster"] for r in csv.DictReader(open(CHEMISTRY_CSV))} \
        if os.path.exists(CHEMISTRY_CSV) else set()
    if args.id in eco or args.id in chem:
        problems.append(f"{args.id} kuratorlu tablolarda zaten var")

    if problems:
        print("DOGRULAMA BASARISIZ:")
        for p in problems:
            print("  -", p)
        return 1

    print(f"[dogrulandi] {args.id}, {len(sequence)} kalinti, substrat '{args.substrate}'")
    if args.dry_run:
        print("[dry-run] hicbir dosya degistirilmedi")
        return 0

    stamp = date.today().isoformat()
    for path in (REFS_FASTA, ECOLOGY_CSV, CHEMISTRY_CSV):
        shutil.copy(path, f"{path}.bak-{stamp}")
    print(f"[yedek] *.bak-{stamp}")

    header = f"{args.id}_{args.description}" if args.description else args.id
    with open(REFS_FASTA, "a") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), 60):
            fh.write(sequence[i:i + 60] + "\n")

    def append_row(path, row):
        with open(path) as fh:
            fields = next(csv.reader(fh))
        with open(path, "a", newline="") as fh:
            csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore").writerow(row)

    append_row(ECOLOGY_CSV, {"cluster": args.id, "substrate": args.substrate,
                             "substrate_class": args.substrate_class,
                             "ecology_note": args.note, "confidence": args.confidence})
    append_row(CHEMISTRY_CSV, {
        "cluster": args.id, "substrate_en": args.substrate, "substrate_smiles": args.smiles,
        "product_en": args.product, "reaction": args.reaction,
        "reaction_class": args.reaction_class, "family": args.family, "pdb": args.pdb,
        "source": args.source, "source_kind": args.source_kind,
        "curation_confidence": args.confidence, "notes": args.note})
    print(f"[eklendi] {REFS_FASTA}, {ECOLOGY_CSV}, {CHEMISTRY_CSV}")

    n = len(existing) + 1
    print(f"""
SONRAKI ADIMLAR -- referans sayisi {len(existing)} -> {n}

1. Modelleri ve motif kolonlarini yeniden kur (hizalama degisti):
     python3 build_reference_models.py --refs {REFS_FASTA} \\
         --out-dir ROs_71_Clean --name ROmotif{n}
   Kolonlar beklenenden farkli cikarsa durup hizalamaya bak; sayilar ekrana
   basilir, sessizce kullanilmaz.

2. Kume atama modelini yeniden kur. RieskeDB{n}.hmm her referansin BLAST
   homologlarindan kuruluyor, bu adim bu depoda otomatik degil: yeni referansin
   homologlarini toplayip profilini ekle, sonra hmmpress.

3. Pipeline'i 1. adimdan itibaren kos: atama ve motif testi degisti, yani
   ro_alpha.csv'den itibaren her sey yeniden uretilmeli:
     bash run_all.sh 16

4. Yalnizca istatistikleri tazelemek istiyorsan (referans degismediyse) her
   adim bagimsiz kosulabilir, hepsi --db alir:
     python3 evidence_tiers.py && python3 cooccurrence.py && \\
     python3 stats_overview.py && python3 substrate_predictability.py && \\
     python3 validate_curation.py

5. Dogrulayiciyi kos: python3 validate_curation.py
""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
