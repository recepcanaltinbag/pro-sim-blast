# ROAR-DB literature review: reference-set redundancy and missing Rieske oxygenase types

Prepared for curator action. Nothing in the database was changed by this review.

Scope and inputs:

- `chemistry.csv` (71 curated reference enzymes)
- `analysis_out/reference_redundancy.json` (computational redundancy report)
- `ROs_71_Clean/refs71.fasta` (the 71 reference sequences actually used to build the profiles)

Method note. Every claim about *which protein a reference record is* was established by exact
sequence comparison of the ROAR-DB reference sequence against UniProtKB entries and against
Protein Data Bank deposited sequences, not by name matching. "EXACT" below means byte-identical
amino acid sequence. Literature claims are cited to primary papers with PMID/DOI.

---

## JOB 1 — Are the redundant reference pairs the same enzyme in the literature?

### Summary table

| Pair | Verdict | Confidence |
|---|---|---|
| `1_101_OxoO` / `3_304_OMO` | **Same protein, two names.** Both are the oxygenase component of 2-oxo-1,2-dihydroquinoline 8-monooxygenase of *Pseudomonas putida* 86. "OMO" is an abbreviation, not a separate enzyme. | Certain |
| `3_309_NahAc` / `3_315_NDO` | **Same protein, two names.** Both are naphthalene 1,2-dioxygenase alpha subunit, UniProt P0A110 (= P0A111). UniProt itself lists `nahAC`, `ndoB`, `nahA3`, `ndoC2` as synonyms of one gene. | Certain |
| `3_314_NDO` / `3_316_NarAa` | **Same protein, two names** — but *not* the protein the metadata claims. Both are NarAa of *Rhodococcus* sp. NCIMB 12038, **not** the *Pseudomonas* enzyme of PDB 1NDO that `3_314` is annotated with. | Certain |
| `1_114_NdmB` / `1_115_NdmC` | **Both records are NdmB.** `1_114` is the His6-TEV-tagged NdmB crystallography construct (PDB 6ICL); `1_115` is untagged NdmB (UniProt H9N290). **NdmC is a genuinely separate gene but is absent from the reference set.** | Certain |
| `2_202_EdoA1` / `2_203_cumA1` | **Two near-identical orthologues from two strains; no evidence they differ in substrate specificity.** The substrate labels reflect the compound each strain was isolated on, not a comparative assay. The single differing residue is not in the active site. | High |
| `1_102_CARDO` / `1_103_CarAa` | **Two orthologues from two organisms** (*Janthinobacterium* sp. J3 and *Pseudomonas resinovorans* CA10), 99.2 % identical. The literature explicitly reports that these two *car* operons are nearly identical. | Certain |

---

### 1. `1_101_OxoO` and `3_304_OMO` — the same enzyme under two names

**Sequence evidence.** Both reference records are 446 aa and are EXACT matches to PDB **1Z03**
chain A, "2-oxo-1,2-dihydroquinoline 8-monooxygenase, oxygenase component".

**Organism.** *Pseudomonas putida* 86 (soil isolate). Gene `oxoO`; the oxygenase component is
446 residues, consistent with the reference length.

**Literature.** The enzyme was purified and the genes cloned by Rosche and colleagues:

- Rosche B, Tshisuaka B, Hauer B, Lingens F, Fetzner S. *2-Oxo-1,2-dihydroquinoline
  8-monooxygenase: phylogenetic relationship to other multicomponent nonheme iron oxygenases.*
  J Bacteriol. 1997;179(11):3549-3554. PMID 9171399. doi:10.1128/jb.179.11.3549-3554.1997
- Martins BM, Svetlitchnaia T, Dobbek H. *2-Oxoquinoline 8-monooxygenase oxygenase component:
  active site modulation by Rieske-[2Fe-2S] center oxidation/reduction.* Structure.
  2005;13(5):817-824. PMID 15893671. doi:10.1016/j.str.2005.03.008 (PDB 1Z01, 1Z02, 1Z03)

**Verdict.** This is one protein. "OxoO" is the gene-product name; "OMO" is simply the acronym
used in reviews for *oxoquinoline monooxygenase*. There is no second enzyme in the literature.
The two records should be merged, keeping `oxoO` as the name. **Confidence: certain.**

Side observation of curatorial interest: the two names were assigned to *different groups*
(group 1 and group 3) by the clustering, so the same sequence currently anchors two different
type labels. 13 database entries are distributed between these two identical profiles (2 vs 11).

---

### 2. `3_309_NahAc` and `3_315_NDO` — the same protein, and the literature already treats several gene names as synonyms

**Sequence evidence.** Both records are 449 aa and are EXACT matches to:

- UniProt **P0A110** (`NDOB_PSEPU`), *Pseudomonas putida*, gene `ndoB`
- UniProt **P0A111** (`NDOB_PSEU8`), *Pseudomonas* sp. strain C18, gene `doxB`
- PDB **1NDO** chain A

**This is the strongest documentation case in the whole set.** UniProt P0A110 carries the gene
name `ndoB` with the explicit synonyms `nahA3`, **`nahAC`**, `ndoC2`. P0A110 and P0A111 are two
separate accessions holding *byte-identical* 449-residue sequences, from two separately published
gene clusters in two differently named strains. In other words, the primary literature named the
same protein `nahAc` (NAH7/classical *nah* operon), `ndoB` (*P. putida* NCIB 9816), `ndoC2`, and
`doxB` (*Pseudomonas* sp. C18, dibenzothiophene/naphthalene upper pathway), and the sequence
databases have since recognised them as one protein.

**Literature.**

- Kurkela S, Lehväslaiho H, Palva ET, Teeri TH. *Cloning, nucleotide sequence and
  characterization of genes encoding naphthalene dioxygenase of Pseudomonas putida strain
  NCIB9816.* Gene. 1988;73(2):355-362. PMID 3243438. doi:10.1016/0378-1119(88)90500-8 (`ndoB`)
- Simon MJ, Osslund TD, Saunders R, Ensley BD, Suggs S, Harcourt AA, Suen WC, Cruden DL,
  Gibson DT, Zylstra GJ. *Sequences of genes encoding naphthalene dioxygenase in Pseudomonas
  putida strains G7 and NCIB 9816-4.* Gene. 1993;127(1):31-37. PMID 8486285.
  doi:10.1016/0378-1119(93)90613-8 (`nahAc`)
- Denome SA, Stanley DC, Olson ES, Young KD. *Metabolism of dibenzothiophene and naphthalene in
  Pseudomonas strains: complete DNA sequence of an upper naphthalene catabolic pathway.*
  J Bacteriol. 1993;175(21):6890-6901. PMID 8226631. doi:10.1128/jb.175.21.6890-6901.1993 (`doxB`)
- Kauppi B, Lee K, Carredano E, Parales RE, Gibson DT, Eklund H, Ramaswamy S. *Structure of an
  aromatic-ring-hydroxylating dioxygenase — naphthalene 1,2-dioxygenase.* Structure.
  1998;6(5):571-586. (PDB 1NDO)

**Verdict.** One protein, four historical gene names. Merge. **Confidence: certain.**
Note that both profiles currently recruit zero database entries, so merging costs nothing.

---

### 3. `3_314_NDO` and `3_316_NarAa` — the same protein, and `3_314`'s PDB annotation is wrong

**Sequence evidence.** Both records are 470 aa and are EXACT matches to UniProt **Q9X3R9**,
*Rhodococcus* sp. NCIMB 12038, gene **`narAa`** ("aromatic ring-hydroxylating dioxygenase subunit
alpha" / "naphthalene dioxygenase large subunit"). The next closest relative is Q0PET4
(*Rhodococcus opacus* `narAa`) at 99.15 % over the same length.

**Important curation error found.** `chemistry.csv` gives `pdb = 1NDO` for `3_314_NDO`. That is
incorrect. PDB 1NDO is the 449-residue *Pseudomonas* enzyme (the one that `3_309`/`3_315` are).
The actinobacterial NarAa is a structurally and phylogenetically distinct naphthalene dioxygenase;
its own structure is PDB **2B1X** / **2B24** (Gakhar et al. 2005, below). Assigning 1NDO to a
Rhodococcus NarAa record conflates two unrelated naphthalene dioxygenase lineages.

**Literature.**

- Larkin MJ, Allen CC, Kulakov LA, Lipscomb DA. *Purification and characterization of a novel
  naphthalene dioxygenase from Rhodococcus sp. strain NCIMB12038.* J Bacteriol.
  1999;181(19):6200-6204. PMID 10498739. doi:10.1128/jb.181.19.6200-6204.1999
- Gakhar L, Malik ZA, Allen CC, Lipscomb DA, Larkin MJ, Ramaswamy S. *Structure and increased
  thermostability of Rhodococcus sp. naphthalene 1,2-dioxygenase.* J Bacteriol.
  2005;187(21):7222-7231. PMID 16237006. doi:10.1128/jb.187.21.7222-7231.2005

**Verdict.** One protein under two names; keep `NarAa` (the name used in the primary literature)
and drop the generic "NDO" record. Correct the PDB cross-reference. This matters practically:
30 of the 32 entries under this identical pair currently sit under `3_316_NarAa`, so the merge is
nearly free, but the wrong PDB id is propagating a misleading structural claim.
**Confidence: certain.**

---

### 4. `1_114_NdmB` and `1_115_NdmC` — both records are NdmB; NdmC is missing from the set

This was the pair the computational report flagged as "contained". The resolution is unambiguous.

**Sequence evidence.**

| Record | Length | Identity |
|---|---|---|
| `1_114_NdmB` | 373 aa | **EXACT** match to PDB **6ICL** chain A, "Methylxanthine N3-demethylase NdmB", *Pseudomonas putida*. The first 18 residues are the expression tag `MGSSHHHHHHENLYFQGS` (His6 + TEV site). |
| `1_115_NdmC` | 355 aa | **EXACT** match to UniProt **H9N290** (`NDMB_PSEPU`), "Methylxanthine N3-demethylase **NdmB**", gene `ndmB`, 355 aa. |

So `1_115_NdmC` is NdmB, mislabelled. `1_114_NdmB` is the same NdmB protein carrying a
crystallography tag. The 18-residue "containment" the redundancy module detected is exactly the
tag. (For completeness: `1_112_NdmA`, 369 aa, is likewise the tagged construct of UniProt
**H9N289** `NDMA_PSEPU` NdmA, 351 aa — correctly named, but tagged.)

**Is NdmC a genuinely separate gene?** Yes, and it is a strikingly different protein.

- UniProt **M1EY73**, "Methylxanthine N7-demethylase", gene `ndmC`, *Pseudomonas putida* CBB5,
  GenBank **JQ061129** / protein **AFD03118.1**, **284 aa**.
- Its only recognised domain is the "vanillate O-demethylase oxygenase-like C-terminal catalytic"
  module (Pfam PF19112), residues 105-274. **It has no Rieske [2Fe-2S] domain of its own.**
  This is why NdmC cannot work alone: it operates as a complex with NdmD (which supplies the
  Rieske centre) and requires the glutathione S-transferase-like NdmE.

**Literature.**

- Summers RM, Louie TM, Yu CL, Gakhar L, Louie KC, Subramanian M. *Novel, highly specific
  N-demethylases enable bacteria to live on caffeine and related purine alkaloids.* J Bacteriol.
  2012;194(8):2041-2049. PMID 22328667. doi:10.1128/JB.06637-11
  (GenBank: `ndmA` JQ061127, `ndmB` JQ061128, **`ndmC` JQ061129**, `ndmD` JQ061130)
- Summers RM, Mohanty SK, Gopishetty S, Subramanian M. *Genetic characterization of caffeine
  degradation by bacteria and its potential applications.* Microb Biotechnol. 2015;8(3):369-378.
  doi:10.1111/1751-7915.12262
- Kim JH, Kim BH, Brooks S, Kang SY, Summers RM, Song HK. *Structural and mechanistic insights
  into caffeine degradation by the bacterial N-demethylase complex.* J Mol Biol.
  2019;431(19):3647-3661. PMID 31412262. doi:10.1016/j.jmb.2019.08.004 (PDB 6ICL = NdmB,
  6ICO = NdmA with theophylline)

**Verdict and recommended action.** Delete one of the two NdmB records (or keep the untagged
H9N290 sequence and drop the tagged one), and **add the real NdmC (M1EY73 / AFD03118.1) as a new
reference**, flagged as a Rieske-domain-less partial oxygenase that functions only in the
NdmCD(E) complex. Note that the reference set currently contains *no* example of this
architecture, which is itself scientifically interesting: a catalytic-domain-only oxygenase that
borrows its Rieske centre in trans. Also note that the stated substrate for `1_115`
(7-methylxanthine, N7-demethylation) describes NdmC, not the sequence stored under that name, so
the chemistry row and the sequence currently disagree. **Confidence: certain** (direct sequence
identity to named UniProt and PDB records).

A further housekeeping point: all three Ndm records, plus `3_307_NagGH` and `1_102_CARDO`, carry
expression tags. Tags add 18-20 non-biological residues to profile training and should be
stripped before HMM building.

---

### 5. `2_202_EdoA1` and `2_203_cumA1` — the scientifically interesting pair

This is the one the owner was right to single out. The short answer: **these are two orthologues
from two different strains that differ at a single surface residue, and there is no published
experiment that compares their substrate preferences. The "ethylbenzene versus cumene" label is an
artefact of which compound each strain was isolated on.**

**Sequence evidence.**

| Record | Length | Identity |
|---|---|---|
| `2_202_EdoA1` | 459 aa | **EXACT** match to UniProt **Q9R956**, "Ethylbenzene dioxygenase large subunit", gene `edoA1`, *Pseudomonas fluorescens* (strain CA-4), EMBL **AF049851** / **AAD12763.1** |
| `2_203_cumA1` | 459 aa | **EXACT** match to UniProt **Q51743**, "Iron-sulfur protein large subunit of cumene dioxygenase", gene `cumA1`, *Pseudomonas fluorescens* IP01, EMBL **D37828** / **BAA07074.1**, = PDB **1WQL** chain A |

The two sequences differ at **exactly one position: 353 (Leu in EdoA1, Trp in CumA1).** A third
close relative, UniProt **P95566** (`ipbA1`, isopropylbenzene 2,3-dioxygenase of *Pseudomonas* sp.
JR1), is 98.5-98.7 % identical to both and also carries Trp353.

**The single difference is not in the active site.** Using the deposited 1WQL coordinates I
measured, for chain A, the minimum distance from each residue to the mononuclear Fe(II):

- The substrate pocket (within 10 Å of Fe) comprises His234, His240, Asp388 (the iron ligands)
  plus Gln227, Phe378, Tyr384, Phe228, Trp392, Leu333, Met232, His323, Ile336, Trp342 and others.
- **Trp353 is 15.4 Å from the mononuclear iron** (Val352 is 11.9 Å; Ser354 is 13.3 Å). For
  comparison, the canonical specificity residue Phe352 of naphthalene dioxygenase (PDB 1NDO) sits
  **5.3 Å** from the iron.

Residue 353 therefore lies outside the first and second shells of the substrate-binding pocket.
A Leu/Trp swap there has no structural basis for switching preference between ethylbenzene and
isopropylbenzene, two substrates that differ by one methyl group.

**What the literature actually says.** The abstract of the paper that deposited `edoA1` is
decisive, and it says the opposite of what the database records imply:

> "*Pseudomonas fluorescens* strain CA-4 is a bioreactor isolate previously characterised by the
> presence of a **side chain oxidation** pathway for ethylbenzene breakdown. In this report a
> second pathway involving ethylbenzene ring dioxygenation has been identified in this strain...
> The genes of the ring-dioxygenation have been cloned and sequenced. **They exhibit near identity
> to the gene clusters encoding the aromatic ring dioxygenase enzymes of two previously described
> isopropyl[benzene] degrading strains, *Pseudomonas* sp. strain JR1 and *P. fluorescens* IP01.**"

— Corkery DM, Dobson AD. *Reverse transcription-PCR analysis of the regulation of ethylbenzene
dioxygenase gene expression in Pseudomonas fluorescens CA-4.* FEMS Microbiol Lett.
1998;166(2):171-176. PMID 9770272. doi:10.1111/j.1574-6968.1998.tb13886.x

Three things follow. (a) The CA-4 paper is a transcriptional regulation study; UniProt records its
evidence as "NUCLEOTIDE SEQUENCE" only — **no purified-enzyme substrate assay was performed on
EdoA1**. (b) The authors themselves state the near-identity to the two isopropylbenzene/cumene
systems. (c) CA-4's originally described ethylbenzene route was side-chain oxidation to
2-phenylethanol (Corkery DM, O'Connor KE, Buckley CM, Dobson AD. *Ethylbenzene degradation by
Pseudomonas fluorescens strain CA-4.* FEMS Microbiol Lett. 1994;124(1):23-27. PMID 8001765), so
"ethylbenzene dioxygenase" in this strain is a name for a ring-attacking cluster found later, not
a demonstration of narrow ethylbenzene specificity.

The cumene enzyme, by contrast, has real enzymology and a structure:

- Habe H, Kimura T, Nojiri H, Yamane H, Omori T. *Cloning and nucleotide sequence of the genes
  involved in the meta-cleavage pathway of cumene degradation in Pseudomonas fluorescens IP01.*
  J Ferment Bioeng. 1996;81:247-254. (source of `cumA1`, EMBL D37828)
- Dong X, Fushinobu S, Fukuda E, Terada T, Nakamura S, Shimizu K, Nojiri H, Omori T, Shoun H,
  Wakagi T. *Crystal structure of the terminal oxygenase component of cumene dioxygenase from
  Pseudomonas fluorescens IP01.* J Bacteriol. 2005;187(7):2483-2490. PMID 15774891.
  doi:10.1128/jb.187.7.2483-2490.2005 (PDB 1WQL)

**Do both enzymes act on both substrates?** The alkylbenzene dioxygenases of this subfamily are
broad. Two independent lines of evidence:

- Chablain PA, Zgoda AL, Sarde CO, Truffaut N. *Genetic and molecular organization of the
  alkylbenzene catabolism operon in the psychrotrophic strain Pseudomonas putida 01G3.* Appl
  Environ Microbiol. 2001;67(1):453-458. PMID 11133479. doi:10.1128/aem.67.1.453-458.2001 — groups
  IpbA1 (IP01/JR1 lineage) and EdoA1 (CA-4) together as one ~99 %-identical cluster of
  "alkylbenzene" dioxygenases, and distinguishes strain groups by *growth* substrate rather than by
  purified-enzyme specificity.
- Berger JB, Marques SM, Wissner JL, Schelle JT, Beekwilder J, Damborsky J, Hauer B.
  *Regioselective sesquiterpene hydroxylation directed by tunnel remodeling in Rieske oxygenases.*
  JACS Au. 2026. PMID 41755873. doi:10.1021/jacsau.5c01130 — shows that wild-type cumene
  dioxygenase from IP01 hydroxylates the sesquiterpene beta-bisabolene, i.e. the enzyme's real
  substrate range extends far beyond cumene. A one-methyl difference between ethylbenzene and
  cumene is well inside that range.

**Verdict.** Two distinct database records (two strains, two accessions, one amino acid apart) but,
on the evidence available, **one functional enzyme type**. The differing substrate annotations
should be treated as *reported growth/induction substrates*, not as a demonstrated specificity
difference. I recommend keeping one reference (CumA1/Q51743, which has enzymology and a structure)
and recording EdoA1 as a strain variant, with the substrate field listing both ethylbenzene and
isopropylbenzene (cumene).

**Confidence: high** for "no published evidence of a specificity difference" and for the
structural placement of residue 353. **I could not find** a head-to-head kinetic comparison of
EdoA1 and CumA1 on both substrates; such an experiment does not appear to have been published, so
the possibility of a real but undocumented difference cannot be formally excluded. Stating this
openly is the honest position.

**Related naming redundancy that the curator should know about.** `2_203_cumA1` is annotated with
substrate *cumene* and `2_206_IpBAa` with substrate *isopropylbenzene*. **These are the same
chemical compound** (cumene = isopropylbenzene; the SMILES in `chemistry.csv` are identical,
`CC(C)c1ccccc1`). `2_206_IpBAa` is an EXACT match to UniProt **O51848** (`ipbAa`,
*Pseudomonas putida* RE204, Eaton RW & Timmis KN, J Bacteriol. 1986;168(1):123-131, PMID 3019995);
it is a genuinely distinct protein (90.6 % identical to CumA1, so a legitimate separate reference)
but the database presents one substrate under two names, which will fragment any
substrate-based analysis.

---

### 6. `1_102_CARDO` and `1_103_CarAa` — two orthologues, recognised as near-identical in the literature

**Sequence evidence.**

| Record | Length | Identity |
|---|---|---|
| `1_102_CARDO` | 392 aa | **EXACT** match to PDB **1WW9** chain A, terminal oxygenase component of carbazole 1,9a-dioxygenase, ***Janthinobacterium* sp. J3**. Carries a C-terminal `LEHHHHHH` tag; the untagged 384 residues are UniProt **Q84II6** (`carAa`, *Janthinobacterium* sp. J3). |
| `1_103_CarAa` | 384 aa | **EXACT** match to UniProt **Q8G8B6** (`CARAA_METRE`), CarAa of ***Pseudomonas (Metapseudomonas) resinovorans* CA10**. |

The two differ at three positions (39 N/D, 312 K/N, 358 S/V) — 99.22 % identity. So these are
**two orthologues from two different organisms**, not one protein under two names.

**The literature explicitly reports this near-identity.**

- Inoue K, Widada J, Nakai S, Endoh T, Urata M, Ashikawa Y, Shintani M, Saiki Y, Yoshida T,
  Habe H, Omori T, Nojiri H. *Divergent structures of carbazole degradative car operons isolated
  from gram-negative bacteria.* Biosci Biotechnol Biochem. 2004;68(7):1467-1480. PMID 15277751.
  doi:10.1271/bbb.68.1467 — reports that the *car* operons of CA10 and J3 have **nearly identical
  nucleotide sequences in their structural and intergenic regions but not in their flanking
  regions**, and defines a "*Pseudomonas*-type *car* gene cluster" shared by *Pseudomonas*,
  *Burkholderia* and *Janthinobacterium* isolates. The pattern (identical core, divergent flanks)
  is the signature of horizontal transfer of a mobile catabolic module.
- Nojiri H, Ashikawa Y, Noguchi H, Nam JW, Urata M, Fujimoto Z, Uchimura H, Terada T, Nakamura S,
  Shimizu K, Yoshida T, Habe H. *Structure of the terminal oxygenase component of angular
  dioxygenase, carbazole 1,9a-dioxygenase.* J Mol Biol. 2005;351(2):355-370. (PDB **1WW9**,
  *Janthinobacterium* sp. J3)
- Nam JW, Noguchi H, Fujimoto Z, Mizuno H, Ashikawa Y, Urata M, Terada T, Nakamura S, Shimizu K,
  Yoshida T, Habe H, Nojiri H. *Crystal structure of the ferredoxin component of carbazole
  1,9a-dioxygenase of Pseudomonas resinovorans strain CA10.* Proteins. 2005;58(4):779-789.
  PMID 15645447.
- Nojiri H et al. *Purification and characterization of carbazole 1,9a-dioxygenase, a
  three-component dioxygenase system of Pseudomonas resinovorans strain CA10.* Appl Environ
  Microbiol. 2002;68(12):5882-5890. (enzymology of the CA10 system)

**Curation error found.** `chemistry.csv` assigns `pdb = 1WW9` to **both** records. 1WW9 is the J3
structure only; the CA10 CarAa oxygenase has no structure of its own under that id. The CA10
record should either have no PDB id or cite the CA10-specific depositions.

**Verdict.** Genuinely two proteins, but functionally one type: same reaction (angular
dioxygenation at C-1/C-9a of carbazole), same product, 3 residues apart, and the literature
describes them as the same mobile *car* cluster. Keeping both is defensible for strain coverage,
but they should be marked as orthologues of one type rather than independent references —
particularly because `1_102_CARDO` recruits all 55 entries and `1_103_CarAa` recruits none, which
is a pure profile-competition artefact. **Confidence: certain** for the identity relationships,
**high** for the horizontal-transfer interpretation.

---

## Cases where the literature or the sequence databases treat separately named enzymes as one protein

The owner suggested that documenting this could be a contribution to the literature. These are
the concrete, citable cases found:

1. **`nahAc` = `ndoB` = `nahA3` = `ndoC2` = `doxB`.** UniProt P0A110 lists the first four as
   synonyms of one gene; P0A111 (`doxB`) holds a byte-identical 449-residue sequence under a
   separate accession and a separate organism name. Four independent primary publications
   (Kurkela 1988; Boronin 1989; Simon 1993; Denome 1993) named the same protein differently
   because they were working on separately isolated plasmids from differently named *Pseudomonas*
   strains. This is a textbook example of nomenclature inflation in this family.
2. **2-oxoquinoline 8-monooxygenase = "OxoO" = "OMO".** One gene product, two acronyms in common
   use; the second is a review abbreviation that has acquired a life of its own in databases.
3. **Cumene = isopropylbenzene.** `cumA1`, `ipbA1`, `ipbAa` and `edoA1` all describe
   alkylbenzene 2,3-dioxygenases of one subfamily (90-99.8 % identity) whose names derive from the
   isolation substrate of each strain, not from measured specificity. Chablain et al. 2001
   (PMID 11133479) already grouped them.
4. **The *car* cluster of CA10 and J3.** Inoue et al. 2004 (PMID 15277751) report near-identical
   structural regions; "CARDO" and "CarAa" are the same enzyme type in two hosts.
5. **`narAa` and `nidA` are the same Rhodococcus/Mycobacterium lineage** under two naming
   conventions (naphthalene dioxygenase vs naphthalene-inducible dioxygenase). ROAR-DB's
   `3_317_NidA` is an EXACT match to UniProt Q9X593 (`nidA`, *Rhodococcus* sp. I24) and is
   98.1-98.9 % identical to Q6TML2/Q6TMM0 (`narAa`, *Rhodococcus* sp. P200/P400) and 98.7 %
   identical to Q2WG94 (`nidA`, *Rhodococcus opacus*). These sit below the 99 % cut-off so the
   redundancy module did not flag them, but they are one type under two gene-name traditions.

A short methods/discussion paragraph making points 1-4 explicit, with these citations, would be a
legitimate and useful contribution: it documents that substrate labels in this enzyme family are
frequently *isolation-substrate labels* rather than specificity measurements.

---

## Additional curation problems found while verifying (not part of the original six pairs)

These were discovered by the same exact-sequence provenance check and are offered because they
affect the same kind of claim.

### `3_307_NagGH` is salicylate 5-hydroxylase, not a naphthalene-2-sulfonate dioxygenase

`chemistry.csv` records substrate "naphthalene-2-sulfonate", product
"1,2-dihydroxynaphthalene and sulfite", reaction "dioxygenation with release of the sulfonate".

The 442-residue reference sequence is an **EXACT match to PDB 7C8Z chain A**, "Salicylate
5-hydroxylase, large oxygenase component", including the 19-residue tag
`MGSHHHHHHSSGLVPRGSH`; the untagged 423 residues are UniProt **O52379** (`NAGG_RALSP`),
salicylate 5-hydroxylase large oxygenase component of *Ralstonia* sp. strain U2,
**EC 1.14.13.172**.

The characterised reaction of this protein is **salicylate (2-hydroxybenzoate) -> gentisate
(2,5-dihydroxybenzoate)**:

- Zhou NY, Al-Dulayymi J, Baird MS, Williams PA. *Salicylate 5-hydroxylase from Ralstonia sp.
  strain U2: a monooxygenase with close relationships to and shared electron transport proteins
  with naphthalene dioxygenase.* J Bacteriol. 2002;184(6):1547-1555.
- Hou YJ, Guo Y, Li DF, Zhou NY. *Structural and biochemical analysis reveals a distinct catalytic
  site of salicylate 5-monooxygenase NagGH from Rieske dioxygenases.* Appl Environ Microbiol.
  2021;87(6):e01629-20. PMID 33452034. doi:10.1128/aem.01629-20 (PDB **7C8Z**)
- Fuenmayor SL, Wild M, Boyes AL, Williams PA. *A gene cluster encoding steps in conversion of
  naphthalene to gentisate in Pseudomonas sp. strain U2.* J Bacteriol. 1998;180(9):2522-2530.

I found **no primary evidence** that NagGH acts on naphthalene-2-sulfonate. The reported
alternative substrates are substituted salicylates (2,4- and 2,6-dihydroxybenzoate). The curated
substrate should be changed to salicylate, the PDB id 7C8Z added, and — if naphthalenesulfonate
chemistry is wanted in the database — a genuine naphthalenesulfonate dioxygenase sought
separately. Note that this also means the set currently has **no** verified salicylate
5-hydroxylase entry even though it holds the protein. **Confidence: certain** on the sequence
identity; **high** on the substrate being wrong (absence of supporting literature).

### `1_113_CdnA`

The 356-residue sequence is an NdmA-type methylxanthine N1-demethylase orthologue (roughly
75-80 % identical to NdmA, UniProt H9N289) from a different organism. It is a legitimate separate
reference, but its product field ("demethylated xanthine") is unspecific and its organism is not
recorded. Worth resolving to a named accession.

### Expression tags in reference sequences

A scan of `ROs_71_Clean/refs71.fasta` for poly-histidine runs finds **exactly four** tagged
records, all taken from PDB/construct sequences:

| Record | Length | Tag |
|---|---|---|
| `1_112_NdmA` | 369 | N-terminal `MGSSHHHHHHENLYFQGS` (+18) |
| `1_114_NdmB` | 373 | N-terminal `MGSSHHHHHHENLYFQGS` (+18) |
| `3_307_NagGH` | 442 | N-terminal `MGSHHHHHHSSGLVPRGSH` (+19) |
| `1_102_CARDO` | 392 | C-terminal `LEHHHHHH` (+8) |

Tags contribute non-biological columns to the alignment and to the HMMs, and in two of these four
cases the tag is what produced a spurious redundancy relationship (the NdmB "containment" and the
CARDO/CarAa length difference). Stripping them before profile building is recommended.

---

## JOB 2 — Characterised Rieske oxygenase types missing from the 71

Selection criteria applied: (i) a Rieske non-heme iron oxygenase (the Rieske-cluster +
mononuclear-iron architecture, including members where the Rieske domain is supplied in trans);
(ii) an experimentally determined substrate, in vitro or by heterologous reconstitution;
(iii) chemistry or substrate class **not already represented** among the 71; (iv) preference for
work from 2020 onward. Further naphthalene/biphenyl/PAH dioxygenases were deliberately excluded.

I checked each candidate against the 71 reference substrates and reaction classes in
`chemistry.csv` before including it.

### Table of candidate new types

| # | Protein / gene | Organism | Substrate | Reaction type | Accession | PDB | Citation | Conf. |
|---|---|---|---|---|---|---|---|---|
| 1 | **TamC**; homologue ***Pt*TamC** | *Pseudoalteromonas citrea*; *P. tunicata* | tambjamine YP1 (linear bipyrrole); also non-native BE-18591 | **Oxidative carbocyclization at an unactivated primary (1°) C–H**, forming a C–C bond. First Rieske-catalysed 1° C–H activation reported. | not located (see Unverified) | none | Ramachandra M, Innis JLM, Yu J, Howe GW, Sauriol F, Oleschuk RD, Ross AC. J Am Chem Soc. 2025. PMID 39870577. doi:10.1021/jacs.4c17468 | High (chemistry), accession unresolved |
| 2 | **RedG**; **McpG** | *Streptomyces coelicolor* A3(2); *Streptomyces longispororuber* | undecylprodigiosin | **Regio- and stereodivergent oxidative carbocyclization** → streptorubin B (RedG) / metacycloprodigiosin (McpG) | RedG **O54095** (SCO5897, 395 aa); McpG **AEL16995** | none | Sydor PK, Barry SM, Odulate OM, Barona-Gomez F, Haynes SW, Corre C, Song L, Challis GL. Nat Chem. 2011;3:388-392. PMID 21505498. doi:10.1038/nchem.1024 | High |
| 3 | **AerC** | aeruginosin-producing cyanobacterium | N-isopentenyl-agmatine (product of the prenyltransferase Aer3) | **Carbocyclization plus C–C bond cleavage** building the Δ3-pyrroline ring of Aeap | not located | none | Zhang W, Ushimaru R, Kanaida M, Abe I. J Am Chem Soc. 2025. PMID 40080531. doi:10.1021/jacs.5c01994 | Medium-high (very recent) |
| 4 | **PrnD** (aminopyrrolnitrin oxygenase) | *Pseudomonas fluorescens* BL915 | aminopyrrolnitrin; also other arylamines | **Arylamine N-oxygenation: Ar-NH2 → Ar-NO2.** Oxygens shown to derive exclusively from O2. | UniProt **P95483** (`PRND_PSEFL`, 363 aa) | none | Lee JK, Simurdiak M, Zhao H. J Biol Chem. 2005;280:36719-36727. PMID 16150698. doi:10.1074/jbc.m505334200 | High |
| 5 | **Nvd (Neverland)**; **DAF-36** | *Drosophila melanogaster*; *Caenorhabditis elegans* (orthologues in *Bombyx*, zebrafish, *Xenopus*) | cholesterol | **Sterol C7–C8 desaturation** → 7-dehydrocholesterol (EC 1.14.19.21). Desaturation, not hydroxylation; animal host. | **Q1JUZ1** (Nvd, 429 aa); **Q17938** (DAF-36, 428 aa) | none | Yoshiyama-Yanagawa T et al. J Biol Chem. 2011;286:25756-25762. PMID 21632547. doi:10.1074/jbc.m111.244384 | High |
| 6 | **CAO** (chlorophyllide a oxygenase / chlorophyll b synthase) | *Arabidopsis thaliana* (and all chlorophyll-b-containing phototrophs) | chlorophyllide a / chlorophyll a | **Two successive oxygenations of the C7 methyl → formyl** (EC 1.14.13.122). Tetrapyrrole substrate; chemistry absent from the set. | **Q9MBA1** (`CAO_ARATH`, 536 aa) | none | Tanaka A et al. Proc Natl Acad Sci USA. 1998;95:12719-12723. PMID 9770552. doi:10.1073/pnas.95.21.12719 — and, for in vitro Rieske characterisation: Liu J, Knapp M, Jo M, Dill Z, Bridwell-Rabb J. ACS Cent Sci. 2022;8:1393-1403. PMID 36313167. doi:10.1021/acscentsci.2c00058 | High |
| 7 | **PAO** (pheophorbide a oxygenase; ACD1 / LLS1) | *Arabidopsis thaliana* | pheophorbide a | **Oxygenolytic opening of the chlorin macrocycle** → red chlorophyll catabolite (EC 1.14.15.17). Ferredoxin-dependent. | **Q9FYC2** (`PAO_ARATH`, 537 aa) | none | Pruzinská A, Tanner G, Anders I, Roca M, Hörtensteiner S. Proc Natl Acad Sci USA. 2003;100:15259-15264. PMID 14657372. doi:10.1073/pnas.2036571100 | High |
| 8 | **TIC55** | *Arabidopsis thaliana* | pFCC (a phyllobilin, linear tetrapyrrole) | **C32 hydroxylation** of a chlorophyll-breakdown phyllobilin; ferredoxin-dependent, membrane-bound | **Q9SK50** (`TIC55_ARATH`, 539 aa) | none | Hauenstein M, Christ B, Das A, Aubry S, Hörtensteiner S. Plant Cell. 2016;28:2510-2527. PMID 27655840. doi:10.1105/tpc.16.00630 | High |
| 9 | **TcsAB** (triclosan dioxygenase) | *Sphingomonas* sp. RD1 | triclosan | Ring hydroxylation / ether-bond cleavage → 2,4-dichlorophenol. Only 35.5 % identical to the nearest group IA alpha subunit, i.e. a genuinely new family member. | GenBank **WYX07280** (TcsA, oxygenase), **WYX07281** (reductase) | none | Yin Y, Ren H, Wu H, Lu Z. Environ Sci Technol. 2024;58:...  PMID 39012163. doi:10.1021/acs.est.4c02845 | High |
| 10 | **BrhABC** (berberine 11-hydroxylase) | *Burkholderia* sp. strain CJ1 | berberine (also palmatine) | **C11 hydroxylation of a protoberberine alkaloid**; stated by the authors to be the first bacterial oxygenase acting on a polycyclic aromatic alkaloid. 38 % identical to DdmC, 33 % to CndA (both already in the set) — so it belongs to a covered family but with new chemistry. | GenBank **BBN97246** (BrhA, oxygenase), **BBN97247** (BrhB, reductase), **BBN97248** (BrhC, ferredoxin) | none | Yoshida H, Takeda H, Wakana D, Sato F, Hosoe T. Biosci Biotechnol Biochem. 2020;84:1130-1138. PMID 32013749. doi:10.1080/09168451.2020.1722056 | High |
| 11 | **PmmA1B1** (+ redundant isoenzymes PmmA2B2, PmmA3B3; ferredoxin Orf05169) | *Pseudomonas putida* NCIMB 9866 | 4-hydroxy-3-methylbenzoate (from the priority pollutant 2,4-xylenol) | **Successive oxidation of an ortho-methyl group all the way to a carboxyl** → 4-hydroxyisophthalate. Rare two-step methyl→carboxyl oxidation; structures reported. | not located (see Unverified) | reported in paper | Meng H, Pan S, Chai Z, Zhou L, Zheng T, Zhang H, Fu H, Zou L, Liu D, Dai J, Yan D, Chao HJ. Appl Environ Microbiol. 2026. PMID 42554493. doi:10.1128/aem.01293-26 | High (chemistry), accession unresolved |
| 12 | **PpoC11** (pyrazon oxygenase, PPO) | *Phenylobacterium immobile* E DSM 1986 | pyrazon (5-amino-4-chloro-2-phenyl-3(2H)-pyridazinone) | cis-2,3-dihydroxylation of the phenyl ring of a pyridazinone herbicide; reconstituted in *E. coli* | **WP_091738533.1** | none | Hunold A, Escobedo-Hinojosa W, Potoudis E, Resende D, Farr T, Syrén PO, Hauer B. Appl Microbiol Biotechnol. 2021;105:2003-2015. PMID 33582834. doi:10.1007/s00253-021-11129-w | High |
| 13 | **SnoT** | *Streptomyces nogalater* | the L-rhodosamine sugar of a nogalamycin intermediate | **2''-hydroxylation of an amino sugar on a glycosylated anthracycline** — a carbohydrate substrate, absent from the set | not located | none | Nji Wandi B, Siitonen V, Palmu K, Metsä-Ketelä M. ChemBioChem. 2020;21:3062-3066. PMID 32557994. doi:10.1002/cbic.202000229 | High (chemistry), accession unresolved |
| 14 | **JerL**, **JerP**, **AmbP** | *Sorangium cellulosum* | jerangolid A and ambruticin VS-3 biosynthetic intermediates | JerL: hydroxymethylpyrone formation. **JerP and AmbP: clean tetrahydropyran desaturation** without over-oxidation — desaturation chemistry on a complex polyketide | GenBank **ABK32295** (JerL), **ABK32293** (JerP), **ABK32265** (AmbP) | none | Guth FM, Lindner F, Rydzek S, Peil A, Friedrich S, Hauer B, Hahn F. ACS Chem Biol. 2023;18:2450-2459. PMID 37948749. doi:10.1021/acschembio.3c00498 | High |
| 15 | **MsmA** (methanesulfonate monooxygenase hydroxylase alpha subunit) | *Methylosulfonomonas methylovora* | methanesulfonate (an **aliphatic** sulfonate) | **C–S bond cleavage** → formaldehyde + sulfite (EC 1.14.13.111). Carries an unusual CXH-X26-CXXH Rieske motif (long spacer). | **Q9X404** (`MSMA_METHY`, 414 aa) | none | De Marco P, Moradas-Ferreira P, Higgins TP, McDonald I, Kenna EM, Murrell JC. J Bacteriol. 1999;181:2244-2251. PMID 10094704 | High |
| 16 | **Rieske-type alkane monooxygenase** (large + small subunit, with ferredoxin and NADH reductase) | *Pusillimonas* sp. strain T7-7 | n-alkanes C5–C24 (maximal on pentadecane, C15); also nitromethane and methanesulfonic acid; **no activity on aromatics** | **Terminal hydroxylation of unactivated aliphatic alkanes.** Stated by the authors to be the first alkane monooxygenase in the Rieske non-heme iron oxygenase family. | not located; genome **CP002735** | none | Li P, Wang L, Feng L. J Bacteriol. 2013;195:1892-1901. PMID 23417490. doi:10.1128/jb.02107-12 | High (chemistry), accession unresolved |
| 17 | **PyfB** | pyrazofurin producer (actinomycete) | 4,6-dihydroxypyridazine-3-carboxylic acid | Oxygenation of a **pyridazine** ring that triggers a non-enzymatic ring contraction to the pyrazole core of pyrazofurin — an oxygenation used to set up a skeletal rearrangement | not located | none | Zheng Z, Lee YH, Ren D, Liu HW. J Am Chem Soc. 2025. PMID 41160757. doi:10.1021/jacs.5c15879 | Medium-high (very recent) |

### Priority recommendation

If only a handful can be curated, I would take, in order: **#5 Nvd/DAF-36** (clean accessions,
desaturation chemistry, animal host, no comparable entry), **#6 CAO** and **#7 PAO** (clean
accessions, tetrapyrrole chemistry, and they complete the plant Rieske oxygenase family of which
the set already has only CmoA/B/S), **#4 PrnD** (clean accession, N-oxygenation), **#2 RedG**
(clean accession, C–C bond formation), **#15 MsmA** and **#12 PpoC11** (clean accessions, new
substrate classes), **#9 TcsAB** and **#10 BrhA** (clean GenBank accessions, 2020/2024, pollutant
and alkaloid chemistry). That is nine entries with verified accessions covering five reaction
types absent from the current set: oxidative C–C carbocyclization, arylamine N-oxygenation, sterol
desaturation, tetrapyrrole macrocycle oxygenation/opening, and aliphatic C–S / C–H hydroxylation.

### Reaction classes currently absent from `chemistry.csv`

For reference, the existing `reaction_class` vocabulary is: hydroxylation, angular_dioxygenation,
O_demethylation, N_demethylation, cis_dihydroxylation, dioxygenation_with_release, C_N_cleavage.
The candidates above would add at least: **oxidative_carbocyclization** (C–C bond formation),
**N_oxygenation** (arylamine → nitro), **desaturation**, **macrocycle_oxygenolysis**, and
**C_S_cleavage** (aliphatic).

### Deliberately excluded

- Further naphthalene, biphenyl, PAH, toluene and nitroarene dioxygenases: the set already has
  ~20 of these and the owner asked for novel chemistry.
- Engineered variants rather than new enzymes: e.g. the thermotolerant ancestral VanA variant
  AncVanA3, which does show enhanced O-demethylation of the lignin monomer **3-O-methylgallate**
  (Lima AR et al. ACS Synth Biol. 2026. PMID 42234957. doi:10.1021/acssynbio.6c00325). 3-O-
  methylgallate is a substrate not present in the set, so this paper is worth reading for the
  substrate even though the protein is an ancestral reconstruction, not a natural isolate.
- Cumene dioxygenase acting on sesquiterpenes (PMID 41755873) — same protein as `2_203_cumA1`,
  so a substrate-scope extension of an existing reference rather than a new type. Worth adding to
  the notes field of `2_203`.
- Enzymes without the Rieske architecture that are sometimes grouped with these (P450 guaiacol
  O-demethylase GcoAB; flavin-dependent LadA; tetrahydrofolate-dependent DesA).

### Enrichments for existing entries found along the way

| Existing record | Enrichment |
|---|---|
| `4_406_PDO` (phthalate dioxygenase) | Molecular insights into substrate recognition and catalysis by phthalate dioxygenase from *Comamonas testosteroni*. J Biol Chem. 2021;297:101416. PMID 34800435 |
| `4_409_TPDO` (terephthalate) | Structural insights into dihydroxylation of terephthalate, a product of polyethylene terephthalate degradation. J Bacteriol. 2022;204:e00543-21. PMID 35007143 |
| `5_504_CntA` (carnitine) | Structural basis of carnitine monooxygenase CntA substrate specificity, inhibition, and intersubunit electron transfer. J Biol Chem. 2021;296:100229. PMID 33158989. Also FEBS J. 2023;290:1907 on an unusual active-site cysteine pair (PMID 36617384) |
| `5_506_Stc2` (stachydrine) | Light-driven oxidative demethylation catalysed by Stc2. ACS Catal. 2022;12:11109. PMID 37168530 |
| `1_110_SxtT`, `1_111_GxtA` | Structural basis for divergent C–H hydroxylation selectivity in two Rieske oxygenases. Nat Commun. 2020;11:2991. PMID 32532989. Design principles for site-selective hydroxylation by a Rieske oxygenase. Nat Commun. 2022;13:255. PMID 35017498. These give the structures and the substrate identities (beta-saxitoxinol for SxtT, saxitoxin for GxtA) that `chemistry.csv` currently records as uncertain |
| `1_106_TsaM` | Leveraging a structural blueprint to rationally engineer the Rieske oxygenase TsaM. Biochemistry. 2023;62:1844. PMID 37188334 |
| `3_313_KshA15` | Further studies on the 3-ketosteroid 9alpha-hydroxylase of *Rhodococcus ruber* Chol-4. Microorganisms. 2021;9:1171. PMID 34072338 |

---

## What I could NOT verify

Stated plainly, as requested.

1. **No head-to-head substrate comparison of EdoA1 and CumA1 exists.** I searched for a kinetic or
   whole-cell assay testing both enzymes on both ethylbenzene and cumene and found none. My
   conclusion that they are functionally one type rests on (a) the authors' own statement of near
   identity, (b) the absence of any enzymology on EdoA1, and (c) the structural observation that
   the single differing residue is 15.4 Å from the catalytic iron. A real but undocumented
   difference cannot be formally excluded.
2. **Accessions I could not locate** despite UniProt and NCBI protein searches:
   **TamC / *Pt*TamC** (candidate 1), **AerC** (3), **PmmA1/B1** (11), **SnoT** (13),
   **PyfB** (17), and the **Pusillimonas sp. T7-7 alkane monooxygenase** subunits (16). For the
   2025/2026 papers this is most likely because the deposited sequences are not yet annotated with
   those gene symbols; a curator should obtain them from the papers' data-availability statements.
   For SnoT the nogalamycin cluster of *S. nogalater* is deposited but the `snoT` symbol did not
   resolve in my queries.
3. **Whether residue 353 of CumA1/EdoA1 lies on a substrate access tunnel.** I established its
   distance to the iron (15.4 Å, i.e. outside the pocket) but did not model the access channel, so
   I cannot exclude an indirect effect on substrate entry. The JACS Au 2026 work on CDO tunnel
   loops (L284, I288, N279, A321) shows such second-shell effects are real in this very enzyme —
   none of those positions is 353.
4. **No support found for naphthalene-2-sulfonate as a NagGH substrate** — but absence of evidence
   in my searches is not proof of absence. The reported alternative substrates are substituted
   salicylates. A curator should check the original ROAR-DB provenance for where the
   naphthalene-2-sulfonate assignment came from.
5. **The identity of `1_113_CdnA`** (organism, strain, accession) was not resolved; it is clearly
   an NdmA-type N1-demethylase orthologue by sequence, but I did not pin its source record.
6. **The 3-aa differences between CARDO (J3) and CarAa (CA10)** have not, as far as I could find,
   been functionally tested. Whether positions 39/312/358 affect activity is unknown.
7. **Exact page numbers** for a few 2024-2026 papers (Environ Sci Technol 2024, the AEM 2026 and
   JACS 2025/2026 items) were not available from the indexing services at the time of writing; the
   DOI and PMID given are authoritative.
8. I did **not** independently re-derive the pairwise identities in
   `analysis_out/reference_redundancy.json`; I recomputed only the specific pairs discussed, and
   they agreed (99.78 % for EdoA1/CumA1, 99.22 % for CARDO/CarAa, 100 % for the three identical
   groups).
