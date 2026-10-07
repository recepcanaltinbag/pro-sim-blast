"""Kuratorluk duzeltmelerini BILDIRIMLE uygular, elle duzenlemeyle degil.

NEDEN. `chemistry.csv` elle duzenlenirse uc sey kaybolur: duzeltmenin NEDEN
yapildigi, KAYNAGI, ve eski degerin ne oldugu. Bu dosyadaki her duzeltme
eski degeri de yazar ve uygulamadan once DOGRULAR -- dosya bu arada baska
bir sekilde degistiyse arac durur ve hicbir sey yazmaz. Boylece ayni script
iki kez calistirilabilir, yarim uygulanmis bir durum olusmaz, ve her degisen
hucrenin yaninda kaynagi durur.

KAYNAK. Duzeltmelerin tamami iki belgeden ve bir yapi bankasindan geliyor:
  Runda M. (2024) Rieske oxygenases outside a nutshell. PhD tezi,
    Groningen Universitesi. doi:10.33612/diss.1175898687
  Miao H, Schmidt S. Biochemistry 2025;64:3801-3813.
    doi:10.1021/acs.biochem.5c00369
  RCSB PDB (her kimlik ayri ayri dogrulandi, bkz. reference_structures.csv)

EN ONEMLI DUZELTME. `3_307_NagGH` yanlis enzimdi. Uc bagimsiz kaynak ayni
seyi soyluyor -- tez Tablo 3 dipnotu onu "Salicylate 5-Monooxygenase" diye
adlandiriyor, tez metni "salicylate 5-hydroxylase ... salisilatin C5 halka
hidroksilasyonunu katalizler" diyor, ve referans dizimizle %100 ayni olan
7C8Z yapisinin RCSB basligi "Crystal structure of salicylate 5-hydroxylase
NagGH". Bizdeki kayit naftalin-2-sulfonat dioksijenasyonu diyordu ve
`curation_confidence: high` isaretliydi. Dizi dogru, ANOTASYON yanlisti --
174 uye bu yuzden yanlis reaksiyon sinifinda ve yanlis kimyasal ailede
sayiliyordu.

STEREOKIMYA. SxtT ve GxtA icin substrat ve urun ADLARI tezden kesin, ama
SMILES tezde YOK; PubChem'den alindi ve iki kayitta da bir stereomerkez
belirsiz. Bunlar yazilir ama `curation_confidence` medium'a cekilir ve
belirsizlik `notes` alanina acikca girer. Yanlis cizilmis bir stereomerkez,
cizilmemis olandan daha kotudur; ama kaydedilmis bir belirsizlik ikisinden
de iyidir.

Kullanim:
    python3 curate.py                 # yalnizca rapor
    python3 curate.py --apply         # chemistry.csv'ye yaz
"""

import argparse
import csv
import os
import shutil
import sys

CHEMISTRY = "chemistry.csv"

THESIS = "Runda M. 2024, PhD thesis, Groningen, doi:10.33612/diss.1175898687"
REVIEW = "Miao H, Schmidt S. Biochemistry 2025;64:3801-3813, doi:10.1021/acs.biochem.5c00369"

# Her duzeltme: (kume, alan, beklenen_eski_deger_ya_da_None, yeni_deger, gerekce)
# beklenen_eski None ise alanin bos olmasi beklenir.
FIXES = [
    # ---------------------------------------------------- 1. yanlis enzim
    ("3_307_NagGH", "substrate_en", "naphthalene-2-sulfonate",
     "salicylate (2-hydroxybenzoate)",
     f"three independent sources call this salicylate 5-hydroxylase: {THESIS} "
     "Table 3 footnote [e] and p.58; and RCSB titles 7C8Z, which is 100 % "
     "identical to this type's curated reference, 'Crystal structure of "
     "salicylate 5-hydroxylase NagGH'"),
    ("3_307_NagGH", "substrate_smiles", "OS(=O)(=O)c1ccc2ccccc2c1",
     "OC(=O)c1ccccc1O",
     "the old SMILES is naphthalene-2-sulfonate, the wrong substrate. "
     "Salicylate drawn from the compound name; neither document gives a SMILES"),
    ("3_307_NagGH", "product_en", "1,2-dihydroxynaphthalene and sulfite",
     "gentisate (2,5-dihydroxybenzoate)",
     "C5 ring hydroxylation of salicylate gives gentisate, the entry point of "
     f"the gentisate pathway in the naphthalene degrader Ralstonia sp. U2 ({THESIS} p.58)"),
    ("3_307_NagGH", "reaction",
     "dioxygenation with release of the sulfonate",
     "hydroxylation at ring position 5 of a hydroxylated benzoate",
     "follows from the substrate and product correction"),
    ("3_307_NagGH", "reaction_class", "dioxygenation_with_release", "hydroxylation",
     "the enzyme hydroxylates; it does not release a substituent. This moves "
     "174 members out of dioxygenation_with_release, where they were 14 % of "
     "the class"),
    ("3_307_NagGH", "family", "sulfoaromatics", "aromatic_acids",
     "salicylate is a hydroxybenzoate, not a sulfonated aromatic. These 174 "
     "members were 21 % of the sulfoaromatics family"),
    ("3_307_NagGH", "source_kind", "general", "paper+structure",
     "it had no individual citation at all; now it has a thesis chapter and a "
     "verified structure"),
    ("3_307_NagGH", "source",
     "established literature on the reference enzyme, no individual citation recorded",
     "PDB 7C8Z (salicylate 5-hydroxylase NagGH, Ralstonia sp. U2); "
     "doi:10.33612/diss.1175898687 Table 3 and p.58",
     "closes a citation gap at the same time as the correction"),
    ("3_307_NagGH", "pdb", "", "7C8Z",
     "verified at RCSB and 100 % identical to this type's curated reference; "
     "it is the evidence the correction rests on"),
    ("3_307_NagGH", "notes", None,
     "Salicylate 5-hydroxylase. Corrected 2026-10-07: this row previously read "
     "naphthalene-2-sulfonate dioxygenation, which three independent sources "
     "contradict. The curated reference sequence is unchanged and is 100 % "
     "identical to PDB 7C8Z. On 5-methylsalicylate the enzyme switches to "
     "methyl hydroxylation.",
     "records the correction in the row itself"),

    # ------------------------------------------- 2. iki bilinen bosluk kapandi
    ("1_110_SxtT", "substrate_en", "saxitoxin pathway intermediate",
     "beta-saxitoxinol",
     f"{THESIS} p.66: 'SxtT catalyzes the oxyfunctionalization of "
     "beta-saxitoxinol into saxitoxin'. p.57 adds that it also accepts "
     "dideoxysaxitoxin but that beta-saxitoxinol is the better and likely "
     "native substrate"),
    ("1_110_SxtT", "product_en", "hydroxylated alkaloid", "saxitoxin",
     f"{THESIS} p.66"),
    ("1_110_SxtT", "reaction", "tailoring hydroxylation on the saxitoxin scaffold",
     "hydroxylation at C12 of the saxitoxin scaffold, giving the geminal diol",
     "the position is what distinguishes SxtT from GxtA"),
    ("1_110_SxtT", "curation_confidence", "high", "medium",
     "the substrate and product names are certain, but the SMILES carries an "
     "unresolved stereocentre"),
    ("1_110_SxtT", "notes",
     "structure not drawn: the scaffold is stereochemically complex and the exact intermediate is uncertain",
     "Substrate and product identified 2026-10-07 from Runda 2024 (thesis) "
     "p.66 and p.57; primary literature Lukowski et al. JACS 2018;140:11863. "
     "SMILES NOT DRAWN: PubChem CID 49789083 'saxitoxinol' does not state "
     "which C12 epimer it is, and the thesis specifies beta. A curator must "
     "set that centre before a structure is drawn for this row.",
     "the gap is now named rather than unexplained, and the remaining "
     "uncertainty is stated precisely"),

    ("1_111_GxtA", "substrate_en", "gonyautoxin pathway intermediate",
     "saxitoxin (also accepts beta-saxitoxinol)",
     f"{THESIS} p.66: 'GxtA can accept both, beta-saxitoxinol and saxitoxin, "
     "catalyzing the incorporation of an OH-group at C11'"),
    ("1_111_GxtA", "product_en", "hydroxylated alkaloid", "11-beta-hydroxysaxitoxin",
     f"{THESIS} p.57: GxtA hydroxylates saxitoxin 'exclusively forming "
     "11-beta-hydroxysaxitoxin, which contains the required hydroxyl group "
     "for the biosynthesis of gonyautoxin'"),
    ("1_111_GxtA", "reaction", "tailoring hydroxylation on the saxitoxin scaffold",
     "hydroxylation at C11 of the saxitoxin scaffold, on the beta face",
     "the position and face are what distinguish GxtA from SxtT"),
    ("1_111_GxtA", "curation_confidence", "high", "medium",
     "as for SxtT: names certain, stereocentre unresolved"),
    ("1_111_GxtA", "notes", "structure not drawn, for the same reason as SxtT",
     "Substrate and product identified 2026-10-07 from Runda 2024 (thesis) "
     "p.66 and p.57; primary literature Lukowski et al. Nat Commun "
     "2020;11:2991. SMILES NOT DRAWN: PubChem CID 165362767 "
     "'11-hydroxysaxitoxin' leaves C11 unspecified and the thesis specifies "
     "beta. A curator must set that centre first.",
     "same treatment as SxtT"),

    # ------------------------------------------------- 3. eksik yapi kimlikleri
    ("3_316_NarAa", "pdb", "", "2B1X",
     f"{THESIS} Table 3, naphthalene dioxygenase of Rhodococcus sp. NCIMB "
     "12038. The only identifier in the thesis that the companion review's "
     "Table 1 does not also carry"),
    ("1_110_SxtT", "pdb", "", "6WN3", f"{THESIS} Table 3; verified at RCSB"),
    ("1_111_GxtA", "pdb", "", "6WNC", f"{THESIS} Table 3; verified at RCSB"),
    ("1_112_NdmA", "pdb", "", "6ICK", f"{THESIS} Table 3; verified at RCSB"),
    ("1_114_NdmB", "pdb", "", "6ICL", f"{THESIS} Table 3; verified at RCSB"),
    ("3_306_3NTDO", "pdb", "", "5XBP",
     f"{THESIS} Table 3 and {REVIEW} Table 1; 5XBP supersedes 5BRC"),
    ("4_406_PDO", "pdb", "", "7FJL", f"{THESIS} Table 3; verified at RCSB"),
    ("4_409_TPDO", "pdb", "", "7VJU",
     f"{THESIS} Table 3. 7Q04 is the same enzyme from Comamonas sp. E6 at "
     "99.5 % identity; 7VJU is the KF-1 enzyme this row is curated from"),
    ("2_217_HcaE", "pdb", "", "8K0A",
     f"{REVIEW} Table 1; verified at RCSB. Its primary citation is 'to be "
     "published', so it is citable only by identifier"),
    ("3_301_PhnA1a", "pdb", "", "2CKF",
     f"{REVIEW} Table 1, PAHDO of Sphingomonas CHY-1; mapped to this type by "
     "sequence at 100 % identity"),
    ("2_213_BphA1A2", "pdb", "", "1ULI",
     f"{REVIEW} Table 1, biphenyl dioxygenase of Rhodococcus jostii RHA1; "
     "mapped by sequence at 100 % identity"),
]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chemistry", default=CHEMISTRY)
    ap.add_argument("--apply", action="store_true",
                    help="dosyaya YAZ. Varsayilan yalnizca rapordur")
    args = ap.parse_args()

    with open(args.chemistry) as fh:
        reader = csv.DictReader(fh)
        fields = list(reader.fieldnames)
        rows = {r["cluster"]: r for r in reader}

    applied, already, refused = [], [], []
    for cluster, field, old, new, why in FIXES:
        row = rows.get(cluster)
        if row is None:
            refused.append((cluster, field, f"type not in {args.chemistry}"))
            continue
        if field not in fields:
            refused.append((cluster, field, "column does not exist"))
            continue
        current = (row.get(field) or "").strip()
        if current == new.strip():
            already.append((cluster, field))
            continue
        expected = "" if old is None else old.strip()
        if current != expected:
            refused.append((cluster, field,
                            f"expected {expected[:48]!r} but found {current[:48]!r}"))
            continue
        applied.append((cluster, field, current, new, why))

    print(f"{len(FIXES)} corrections declared")
    print(f"  {len(applied):3d} to apply")
    print(f"  {len(already):3d} already in place")
    print(f"  {len(refused):3d} REFUSED (the file does not match what was expected)")
    print()
    for cluster, field, cur, new, why in applied:
        print(f"  {cluster}  {field}")
        print(f"      from: {cur[:90] or '(empty)'}")
        print(f"      to:   {new[:90]}")
        print(f"      why:  {why[:150]}")
    if refused:
        print()
        print("REFUSED -- nothing was written for these:")
        for cluster, field, reason in refused:
            print(f"  {cluster}  {field}: {reason}")

    if not args.apply:
        print()
        print("Report only. Add --apply to write.")
        print("Afterwards: python3 rebuild.py --run   (chemistry.csv feeds many "
              "analyses and they will be stale)")
        return 1 if refused else 0

    if refused:
        print()
        print("Refusals present -- nothing written. Resolve them first.")
        return 1
    if not applied:
        print("\nNothing to do.")
        return 0

    shutil.copy2(args.chemistry, args.chemistry + ".bak")
    for cluster, field, _cur, new, _why in applied:
        rows[cluster][field] = new
    with open(args.chemistry, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for cluster in sorted(rows):
            w.writerow(rows[cluster])
    print()
    print(f"[yazildi] {args.chemistry} ({len(applied)} hucre)")
    print(f"[yedek]   {args.chemistry}.bak")
    print()
    print("Simdi: python3 rebuild.py --run")
    return 0


if __name__ == "__main__":
    sys.exit(main())
