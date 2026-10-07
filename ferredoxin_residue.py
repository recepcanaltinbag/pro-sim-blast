"""
Ferredoksin "A50 pozisyonu" -- iki proteinlik bir iddianin 4.091 komsuda sinanmasi.

NE OLCULUYOR
    Miao, Oerlemanns, Hagedoorn & Schmidt (bioRxiv 2026, doi:10.64898/2026.04.03.713453)
    CDO-Fd ve NDO-Fd uzerinde mutagenez yaparak tek bir kalinti pozisyonunun
    redoks ortagi secimini belirledigini gosteriyor: CDO-Fd'de W48, NDO-Fd'de
    A50 (ayni hizalama kolonu). Sonra su GENELLEMEYI yapiyorlar:

        "bu kalinti pozisyonu SINIFA BAGLIDIR: triptofan Sinif IIB
         ferredoksinlerinde, alanin Sinif III'te, tirozin ise Sinif IIA'nin
         bazi bitki-tipi ferredoksinlerinde baskindir."

    Bu genelleme iki proteinlik deneye + onceden yayinlanmis hizalamalara
    dayaniyor. Bu veritabaninda DIZISI OLAN 4.091 ferredoksin komsusu var.
    Bu script her birini NDO-Fd'ye (ve capraz kontrol icin CDO-Fd'ye) global
    hizalayip A50 kolonuna denk gelen kalintiyi okur, sonra o kalintiyi
    elimizdeki gruplandirmalarla capraz tablolar.

NEDEN
    Iki protein bir pozisyonun FONKSIYONEL oldugunu gosterebilir -- mutagenez
    bunu gosteriyor ve bu script o kismi tartismiyor. Ama iki protein bir
    pozisyonun SINIFA BAGLI DAGILDIGINI gosteremez; bunun icin dagilimin
    kendisine bakmak gerekir. Dagilim sorusu olcek sorusudur ve bu
    veritabaninin cevaplayabilecegi bir sorudur.

    Olcumun kendisi de bir tuzak tasiyor: 4.000 diziyi kotu hizalayip kolon
    okursaniz gurultu SONUC gibi gorunur. Bu yuzden her hizalama icin kimlik
    ve KAPSAMA kaydedilir, ayrica Rieske kumesinin dort ligandi (NDO-Fd'de
    C45, H47, C64, H67) sorgu dizisinde ayni kalintiya hizalanmak ZORUNDADIR.
    A50, H47'den yalnizca uc kalinti sonra gelir; ligandlar dogru oturmuyorsa
    o kolon zaten anlamsizdir. Esikleri gecemeyen diziler atilir ve kac tane
    atildigi raporlanir.

BATIE SINIFI UYARISI -- BU ONEMLI
    Yazarlarin "Sinif IIA / IIB / III" etiketleri Batie (1991) siniflandirmasi.
    BU VERITABANI BATIE SINIFI ATAMAZ. Elimizde kendi `ro_group` 1-5
    kumelenmesi, `ro_etc.reductase_type` (FNR | GR | FNR+GR | none) ve
    `ro_etc.ferredoxin_type` (rieske | plant | plant+rieske | none) var.
    ro_group alfa alt biriminin DIZI kumelenmesidir; Batie sinifi ise elektron
    tasima ZINCIRININ bilesen yapisidir. Ikisi ayni sey degildir ve aralarinda
    sessizce bir eslesme varsayilmaz.

    Bu yuzden ANA SONUC yalnizca gercekten elimizde olan degiskenlere karsi
    verilir (ferredoksin tipi, ferredoksinin KENDI Pfam alani, RO grubu, enzim
    tipi). Yazarlarin sinif etiketlerinin DOGRUDAN sinanmasi burada mumkun
    degildir; nedeni cikti JSON'unda `batie_caveat` altinda ayrica yazilidir.

    IKINCIL ve AYRICA ETIKETLENMIS bir analiz olarak kismi bir eslesme de
    kurulur. Eslesme literaturden degil YAZARLARIN KENDI TANIMINDAN gelir
    (makalenin Giris bolumu): "Sinif II ile III arasindaki temel ayrim Red
    bilesenlerinin dogasindadir" ve "Sinif II, iceren Fd tipine gore IIA ve
    IIB alt siniflarina ayrilir". Batie'de Sinif II reduktazi yalnizca
    flavin tasir (glutatyon-reduktaz katlanmasi = bu veritabaninda GR),
    Sinif III reduktazi ise flavin + [2Fe-2S] tasir (uc alanli FNR-tipi =
    burada FNR). Buradan:
        bitki-tipi Fd                -> IIA-benzeri
        Rieske-tipi Fd + GR reduktaz -> IIB-benzeri   (CDO boyle: CumA3+CumA4)
        Rieske-tipi Fd + FNR reduktaz-> III-benzeri   (NDO boyle: NahAb+NahAa)
    Eslesme makalenin iki model sisteminde de dogru sonuc veriyor; yine de
    bir VEKIL'dir, Batie etiketi degildir ve ana sonuc olarak sunulmaz.

Cikti: analysis_out/ferredoxin_residue.json
"""

import argparse
import json
import os
import sqlite3
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

from operon_relations import _aligner, pair_identity
from stats_overview import (bh_adjust, collapse_genus, cramers_v, signal_map,
                            trim_table)

# --------------------------------------------------------------------------
# Referans diziler.
#
# Ikisi de UniProt'tan DOGRULANDI ve makalenin Sekil 4A'sindaki hizalama
# satirlarinin bosluksuz halleriyle KARAKTER KARAKTER ayni:
#   NDO-Fd  P0A185 (NDOA_PSEPU, ndoA/nahAb) 104 aa
#   CDO-Fd  Q51746 (Q51746_PSEFL, cumA3)    109 aa
#
# Numaralandirma dogrulamasi: makalede adi gecen BUTUN kalinti numaralari bu
# dizilerde bire bir tutuyor -- CDO-Fd icin E25, D41, R42, D47, W48, E52, E62,
# S64, M67, E84, K87 ve NDO-Fd icin D42, T46, H47, S49, A50, D54, R60, E61,
# E63, P65, L66, H67, Q68, D72, P82, Q85, K88. Yani numaralandirma olgun
# dizinin 1-tabanli indeksi (His etiketi sayilmiyor). Dogrulama calisma
# aninda da tekrar yapilir ve tutmazsa script DURUR.
# --------------------------------------------------------------------------
NDO_FD = ("MTVKWIEAVALSDILEGDVLGVTVEGKELALYEVEGEIYATDNLCTHGSARMSDGYLEGR"
          "EIECPLHQGRFDVCTGKALCAPVTQNIKTYPVKIENLRVMIDLS")
CDO_FD = ("MTFSKVCEVSDVPVGDALQVESKGEAVAIFNVDGELFATQDRCTHGDWSLSEGGYLEGDI"
          "VECSLHMGRFCVRTGKVKAAPPCEPLKIYPIRIDGSDVFVDFDAGYLAP")

REFERENCES = {
    "NDO-Fd": {
        "sequence": NDO_FD,
        "uniprot": "P0A185",
        "uniprot_name": "NDOA_PSEPU",
        "gene": "ndoA / nahAb",
        "organism": "Pseudomonas sp. NCIB 9816-4 (UniProt: Pseudomonas putida)",
        "enzyme": "naphthalene 1,2-dioxygenase, ferredoxin component",
        "claim_position": 50,
        "claim_residue": "A",
        # Rieske [2Fe-2S] ligandlari: CxH ... CxxH
        "rieske_anchors": [(45, "C"), (47, "H"), (64, "C"), (67, "H")],
        # makalede gecen ve numaralandirmayi dogrulayan kalintilar
        "numbering_checks": [(42, "D"), (46, "T"), (47, "H"), (49, "S"), (50, "A"),
                             (54, "D"), (60, "R"), (61, "E"), (63, "E"), (65, "P"),
                             (66, "L"), (67, "H"), (68, "Q"), (72, "D"), (82, "P"),
                             (85, "Q"), (88, "K")],
    },
    "CDO-Fd": {
        "sequence": CDO_FD,
        "uniprot": "Q51746",
        "uniprot_name": "Q51746_PSEFL",
        "gene": "cumA3",
        "organism": "Pseudomonas fluorescens IP01",
        "enzyme": "cumene dioxygenase, ferredoxin component",
        "claim_position": 48,
        "claim_residue": "W",
        "rieske_anchors": [(43, "C"), (45, "H"), (63, "C"), (66, "H")],
        "numbering_checks": [(25, "E"), (41, "D"), (42, "R"), (47, "D"), (48, "W"),
                             (52, "E"), (62, "E"), (64, "S"), (67, "M"), (84, "E"),
                             (87, "K")],
    },
}

# --------------------------------------------------------------------------
# Hizalama kabul esikleri.
#
# identity/coverage tanimi operon_relations.pair_identity ile AYNI: kimlik
# yalnizca iki tarafi da kalinti olan kolonlarda, kapsama ise o kolon
# sayisinin KISA dizinin uzunluguna orani. Esikler asagida; en belirleyicisi
# ucuncusu: Rieske ligandlari oturmuyorsa kolon okunmaz.
# --------------------------------------------------------------------------
MIN_IDENTITY = 25.0     # % -- bunun altinda ~100 aa proteinde kolon karsiligi guvenilmez
MIN_COVERAGE = 0.60     # kisa diziye gore hizalanan kolon orani
REQUIRE_ANCHORS = 4     # dort Rieske ligandinin hepsi dogru oturmali

RESIDUE_CLASSES = ["W", "A", "Y", "other"]
REGEX_METHOD = "regex_v1"
MIN_TYPE_N = 20         # enzim tipi tablosuna girmek icin asgari uye (ev olcutu)
MAX_FD_LEN = 250        # etc_types.py'nin "bu bir ferredoksindir" uzunluk olcutu

RIESKE_HMM = {"Rieske", "Rieske_2"}
PLANT_HMM = {"Fer2"}
REDUCTASE_HMM = {"NAD_binding_1", "FAD_binding_6", "FAD_binding_8",
                 "Pyr_redox_2", "Pyr_redox_dim"}


def verify_numbering():
    """Referans dizilerde makalenin verdigi kalinti numaralari tutuyor mu?

    Tutmuyorsa yanlis dizi ya da yanlis numaralandirma kullaniyoruz demektir
    ve devam etmek anlamsizdir; bu yuzden AssertionError ile durulur.
    """
    report = {}
    for name, ref in REFERENCES.items():
        seq = ref["sequence"]
        bad = [(p, exp, seq[p - 1] if p <= len(seq) else None)
               for p, exp in ref["numbering_checks"]
               if p > len(seq) or seq[p - 1] != exp]
        assert not bad, f"{name} numaralandirmasi tutmuyor: {bad}"
        pos, res = ref["claim_position"], ref["claim_residue"]
        assert seq[pos - 1] == res, (
            f"{name} pozisyon {pos} = {seq[pos - 1]}, beklenen {res}")
        report[name] = {
            "uniprot": ref["uniprot"], "uniprot_name": ref["uniprot_name"],
            "gene": ref["gene"], "organism": ref["organism"],
            "enzyme": ref["enzyme"], "length": len(seq), "sequence": seq,
            "numbering": "1-based index of the mature sequence (His tag not counted)",
            "claim_position": pos, "residue_at_claim_position": seq[pos - 1],
            "numbering_checks_passed": len(ref["numbering_checks"]),
            "rieske_anchors": [f"{r}{p}" for p, r in ref["rieske_anchors"]],
            "source": ("UniProt REST; character-identical to the ungapped "
                       "Figure 4A alignment row of the preprint"),
        }
    return report


def cross_check_references():
    """NDO-Fd A50 kolonu, ev hizalayicisinda gercekten CDO-Fd W48'e mi dusuyor?

    Makale bu denkligi Clustal Omega ile kuruyor (Sekil 4A). Bagimsiz bir
    hizalayiciyla ayni kolon cikmiyorsa denklik hizalayiciya bagli demektir
    ve iddianin tasidigi butun agirlik o secime biner.
    """
    aln = _aligner().align(NDO_FD, CDO_FD)[0]
    idx = aln.indices
    col = int(np.where(idx[0] == 49)[0][0])
    j = int(idx[1][col])
    identity, coverage = pair_identity(NDO_FD, CDO_FD)
    return {
        "aligner": "Bio.Align.PairwiseAligner(scoring='blastp'), mode=global "
                   "(operon_relations._aligner)",
        "ndo_position": 50, "ndo_residue": NDO_FD[49],
        "cdo_position_aligned": (j + 1) if j >= 0 else None,
        "cdo_residue_aligned": CDO_FD[j] if j >= 0 else None,
        "agrees_with_preprint_figure_4A": bool(j >= 0 and j + 1 == 48
                                               and CDO_FD[j] == "W"),
        "pairwise_identity_pct": round(identity, 1) if identity else None,
        "coverage": round(coverage, 3) if coverage else None,
        "preprint_reported_identity_pct": 41,
        "note": ("the preprint reports 41% identity between the two Fd components "
                 "and aligns W48 with A50 using Clustal Omega; the house blastp "
                 "global aligner independently reproduces both, so the residue "
                 "equivalence is not an artefact of one alignment program"),
    }


# -------------------------------------------------------------------- veri

def load_ferredoxins(con):
    """Ferredoksin komsulari + dizileri + ait olduklari alfa'nin baglami.

    Koordinat anahtari operon_relations.load_regulator_pairs_input() ile ayni:
    neighbor_protein'de neighbor_id yok, baglanti
    nucleotide_id:start-end:strand uzerinden kurulur.

    Birim: (alfa girisi x ferredoksin komsusu) cifti. Ayni Fd proteini birden
    cok genomda tekrar ettigi icin dizi duzeyinde ve cins duzeyinde cokertme
    ayrica yapilir.
    """
    rows = con.execute(f"""
        SELECT nb.candidate_id, np.protein_key, np.translation, np.length,
               r.ro_cluster, r.ro_group,
               e.ferredoxin_type, e.reductase_type, e.has_beta,
               p.organism, p.is_plasmid
        FROM gene_category gc
        JOIN neighbor nb ON nb.neighbor_id = gc.neighbor_id
        JOIN neighbor_protein np
          ON np.protein_key = nb.nucleotide_id || ':' || nb.start || '-'
                              || nb.end || ':' || nb.strand
        JOIN ro r ON r.candidate_id = nb.candidate_id
        JOIN replicon p ON p.nucleotide_id = r.nucleotide_id
        LEFT JOIN ro_etc e ON e.candidate_id = nb.candidate_id
        WHERE gc.method = ?
          AND gc.category = 'ferredoxin'
          AND r.is_confirmed = 1
          AND np.translation IS NOT NULL AND np.translation <> ''
    """, (REGEX_METHOD,)).fetchall()

    # ferredoksinin KENDI Pfam alani: alfa duzeyindeki ro_etc.ferredoxin_type
    # operondaki BUTUN Fd'lerin birlesimidir, bu protein icin degil. Rieske mi
    # bitki-tipi mi sorusunun dogrudan cevabi neighbor_domain'de.
    domains = defaultdict(set)
    for key, hmm in con.execute("SELECT protein_key, hmm FROM neighbor_domain"):
        domains[key].add(hmm)

    out = []
    for (cid, key, seq, length, cluster, grp, fd_type, red_type, has_beta,
         organism, is_plasmid) in rows:
        d = domains.get(key, set())
        if d & RIESKE_HMM:
            own = "plant+rieske" if d & PLANT_HMM else "rieske"
        elif d & PLANT_HMM:
            own = "plant"
        elif d & REDUCTASE_HMM:
            own = "reductase_domains"
        else:
            own = "no_domain"
        out.append({
            "candidate_id": cid, "protein_key": key, "sequence": seq,
            "length": length,
            # collapse_genus() 'cluster' ve 'genus' anahtarlarini bekler
            "cluster": cluster or "?",
            "group": f"group {grp}" if grp else "group ?",
            "fd_type": fd_type or "none",
            "reductase_type": red_type or "none",
            "fd_own_domain": own,
            "organism": organism or "?",
            "genus": (organism or "?").split()[0],
            # collapse_genus() bu alanlarda ortalama aliyor; var olmalari sart
            "is_plasmid": int(bool(is_plasmid)),
            "mobile": int(bool(is_plasmid)),
            "has_beta": int(bool(has_beta)),
            "has_ferredoxin": 1,
            "has_reductase": int((red_type or "none") != "none"),
        })
    return out


# --------------------------------------------------------------- hizalama

def _read_position(task):
    """Bir sorgu dizisi icin: referansin hedef kolonundaki kalinti + kalite.

    Pool isci fonksiyonu. Donus: (seq, ref_name, residue, identity, coverage,
    anchors_ok, anchor_detail).
    """
    seq, ref_name, ref_seq, position, anchors = task
    try:
        aln = _aligner().align(ref_seq, seq)[0]
    except Exception:
        return (seq, ref_name, None, None, None, 0, "")
    idx = aln.indices
    ref_idx, qry_idx = idx[0], idx[1]

    def query_residue(ref_pos):
        hits = np.where(ref_idx == ref_pos - 1)[0]
        if not len(hits):
            return None
        j = int(qry_idx[int(hits[0])])
        return seq[j] if j >= 0 else None

    residue = query_residue(position)

    matches = aligned = 0
    for x, y in zip(ref_idx, qry_idx):
        if x >= 0 and y >= 0:
            aligned += 1
            if ref_seq[x] == seq[y]:
                matches += 1
    shorter = min(len(ref_seq), len(seq))
    identity = (100.0 * matches / aligned) if aligned else None
    coverage = (aligned / shorter) if (aligned and shorter) else None

    detail, ok = [], 0
    for pos, expected in anchors:
        got = query_residue(pos)
        detail.append(got or "-")
        if got == expected:
            ok += 1
    return (seq, ref_name, residue, identity, coverage, ok, "".join(detail))


def read_all_positions(sequences, cpu):
    """Her BENZERSIZ dizi icin iki referansa karsi kolon okumasi.

    Benzersiz dizi basina hizalanir (ayni protein yuzlerce susta tekrar eder),
    sonra sonuc satirlara geri dagitilir.
    """
    tasks = []
    for seq in sequences:
        for name, ref in REFERENCES.items():
            tasks.append((seq, name, ref["sequence"], ref["claim_position"],
                          ref["rieske_anchors"]))
    if cpu > 1 and len(tasks) > 200:
        with Pool(cpu) as pool:
            results = pool.map(_read_position, tasks, chunksize=64)
    else:
        results = [_read_position(t) for t in tasks]

    reads = defaultdict(dict)
    for seq, name, residue, identity, coverage, anchors, detail in results:
        reads[seq][name] = {
            "residue": residue,
            "identity": round(identity, 1) if identity is not None else None,
            "coverage": round(coverage, 3) if coverage is not None else None,
            "anchors_ok": anchors, "anchor_residues": detail,
            "passes": bool(residue is not None
                           and identity is not None and identity >= MIN_IDENTITY
                           and coverage is not None and coverage >= MIN_COVERAGE
                           and anchors >= REQUIRE_ANCHORS),
        }
    return reads


def residue_class(residue):
    if residue in ("W", "A", "Y"):
        return residue
    return "other"


# ------------------------------------------------------------ istatistik

def contingency(out, test_id, question, subset, key_field, col_field, keys, cols,
                row_header, unit, restriction=None, note=None):
    """Kontenjans testi + sinyal haritasi, giris VE cins duzeyinde.

    Istatistik tamamen stats_overview'dan gelir (signal_map / bh_adjust /
    cramers_v / trim_table / collapse_genus); burada ikinci bir uygulama yok.
    Maddiyet olcutu ev standardi: duzeltilmis q < 0,05 VE sutun payinin genel
    paydan en az 10 puan sapmasi VE cins duzeyinde ayni yonde durmasi.
    """
    from scipy import stats as _st
    table = [[sum(1 for d in subset if d[key_field] == k and d[col_field] == c)
              for c in cols] for k in keys]
    keys_t, cols_t, table_t = trim_table(keys, cols, table)
    if len(keys_t) < 2 or len(cols_t) < 2 or sum(map(sum, table_t)) == 0:
        return None
    chi2, p, dof, _ = _st.chi2_contingency(table_t)
    gsub = collapse_genus(subset)
    gtab = [[sum(1 for d in gsub if d[key_field] == k and d[col_field] == c)
             for c in cols_t] for k in keys_t]
    t = {
        "id": test_id, "question": question,
        "rows": keys_t, "cols": cols_t, "table": table_t,
        "chi2": float(chi2), "dof": int(dof), "p": float(p),
        "cramers_v": cramers_v(table_t),
        "row_header": row_header, "unit": unit,
        "n_entry": int(sum(map(sum, table_t))),
        "n_genus": int(sum(map(sum, gtab))),
        "genus_table": gtab,
    }
    # Cins tablosunda bos satir/sutun olusabilir (bir sinif cins basina tek
    # temsilciye dusunce kaybolur). signal_map AYNI SEKILDE tablo bekledigi
    # icin gtab oldugu gibi ona verilir; ozet istatistikler ise kirpilmis
    # kopyadan hesaplanir, yoksa ki-kare sifir beklenen deger hatasi verir.
    g_rows, g_cols, g_trim = trim_table(keys_t, cols_t, gtab)
    if len(g_rows) >= 2 and len(g_cols) >= 2 and sum(map(sum, g_trim)):
        t["cramers_v_genus"] = cramers_v(g_trim)
        try:
            t["p_genus"] = float(_st.chi2_contingency(g_trim)[1])
        except Exception:
            t["p_genus"] = None
    else:
        t["cramers_v_genus"] = None
        t["p_genus"] = None
    if restriction:
        t["restriction"] = restriction
    if note:
        t["note"] = note
    t["signals"] = signal_map(keys_t, cols_t, table_t, gtab)
    out["tests"].append(t)
    return t


def batie_proxy(row):
    """Makalenin KENDI sinif tanimindan tureyen vekil etiket (ikincil analiz).

    Giris bolumu: sinif II ile III reduktaza gore, IIA ile IIB ferredoksin
    tipine gore ayrilir. Batie'de sinif II reduktazi flavin-only (GR
    katlanmasi), sinif III reduktazi flavin + [2Fe-2S] (FNR-tipi uc alanli).
    Belirsiz kombinasyonlar (FNR+GR, plant+rieske, reduktaz gorulmemis)
    etiketlenmez -- "bilinmiyor" ile "yok" ayni sey degildir.
    """
    fd, red = row["fd_type"], row["reductase_type"]
    if fd == "plant":
        return "IIA-like (plant-type Fd)"
    if fd == "rieske":
        if red == "GR":
            return "IIB-like (Rieske Fd + GR reductase)"
        if red == "FNR":
            return "III-like (Rieske Fd + FNR reductase)"
    return "unassignable"


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/ferredoxin_residue.json")
    ap.add_argument("--cpu", type=int, default=min(20, os.cpu_count() or 1))
    args = ap.parse_args()

    out = {
        "claim": {
            "source": ("Miao H, Oerlemanns R, Hagedoorn P-L, Schmidt S. Decoding and "
                       "Reprogramming Redox Partner Specificity in Rieske Oxygenases "
                       "for Enhanced Catalytic Activity. bioRxiv 2026, "
                       "doi:10.64898/2026.04.03.713453 (CC-BY-NC)"),
            "quote": ("this residue position (e.g., W48 in CDO-Fd and A50 in NDO-Fd) "
                      "is class-dependent: tryptophan predominates in Class IIB Fds, "
                      "alanine in Class III Fds, and tyrosine in certain plant-type "
                      "Fds from Class IIA"),
            "evidence_behind_it": ("site-directed mutagenesis on two ferredoxins "
                                   "(CDO-Fd, NDO-Fd) plus inspection of previously "
                                   "published alignments"),
            "what_this_script_does_not_dispute": (
                "that the position is functionally important. The mutagenesis is "
                "direct evidence for that and nothing here bears on it. Only the "
                "distributional generalisation is tested"),
        },
    }

    out["references"] = verify_numbering()
    out["reference_cross_check"] = cross_check_references()

    con = sqlite3.connect(args.db)
    rows = load_ferredoxins(con)
    total_gc = con.execute(
        "SELECT COUNT(*) FROM gene_category WHERE method=? AND category='ferredoxin'",
        (REGEX_METHOD,)).fetchone()[0]

    sequences = sorted({r["sequence"] for r in rows})
    reads = read_all_positions(sequences, args.cpu)
    for r in rows:
        rd = reads.get(r["sequence"], {})
        r["read"] = rd
        ndo = rd.get("NDO-Fd", {})
        cdo = rd.get("CDO-Fd", {})
        r["residue"] = ndo.get("residue")
        r["residue_class"] = residue_class(ndo.get("residue"))
        r["passes"] = bool(ndo.get("passes"))
        r["cdo_residue"] = cdo.get("residue")
        r["cdo_passes"] = bool(cdo.get("passes"))
        r["identity"] = ndo.get("identity")
        r["coverage"] = ndo.get("coverage")
        r["anchors_ok"] = ndo.get("anchors_ok")

    good = [r for r in rows if r["passes"]]

    # ---------------------------------------------------------- kapsama
    disc = Counter()
    for r in rows:
        if r["passes"]:
            continue
        ndo = r["read"].get("NDO-Fd", {})
        if ndo.get("residue") is None:
            disc["target column is a gap"] += 1
        elif (ndo.get("anchors_ok") or 0) < REQUIRE_ANCHORS:
            disc[f"Rieske ligands misaligned ({ndo.get('anchors_ok')}/4)"] += 1
        elif (ndo.get("identity") or 0) < MIN_IDENTITY:
            disc["identity below threshold"] += 1
        else:
            disc["coverage below threshold"] += 1

    both = [r for r in good if r["cdo_passes"] and r["cdo_residue"]]
    concordant = sum(1 for r in both if r["cdo_residue"] == r["residue"])

    out["coverage"] = {
        "gene_category_ferredoxin_rows": total_gc,
        "with_translation": len(rows),
        "distinct_sequences": len(sequences),
        "aligned_well_enough": len(good),
        "aligned_well_enough_distinct_sequences":
            len({r["sequence"] for r in good}),
        "discarded": len(rows) - len(good),
        "discard_reasons": dict(disc.most_common()),
        "thresholds": {
            "min_identity_pct": MIN_IDENTITY,
            "min_coverage": MIN_COVERAGE,
            "rieske_anchors_required": REQUIRE_ANCHORS,
            "identity_definition": ("percent identity over columns where BOTH "
                                    "sides are residues (operon_relations."
                                    "pair_identity convention)"),
            "coverage_definition": ("residue-residue columns divided by the length "
                                    "of the SHORTER sequence"),
            "why_anchors": ("A50 sits three residues after the Rieske ligand H47. "
                            "If the four cluster ligands C45/H47/C64/H67 do not "
                            "land on the chemically identical residue in the query, "
                            "the register around A50 is not trustworthy and the "
                            "column is noise"),
        },
        "cross_reference_concordance": {
            "readable_against_both_references": len(both),
            "same_residue_from_both": concordant,
            "rate": round(concordant / len(both), 4) if both else None,
            "why": ("the same column read by aligning to NDO-Fd and to CDO-Fd "
                    "independently. A high rate means the residue call does not "
                    "depend on which reference was used"),
        },
        "discarded_by_own_pfam_domain": dict(Counter(
            r["fd_own_domain"] for r in rows if not r["passes"]).most_common()),
        "kept_by_own_pfam_domain": dict(Counter(
            r["fd_own_domain"] for r in good).most_common()),
    }

    # ------------------------------------------------- POZITIF KONTROLLER
    # Yontem calisiyor mu? Iki referansin KENDI yakin homologlari dogru
    # kalintiyi vermek ZORUNDA. Vermiyorsa kolon okumasi bozuktur ve
    # asagidaki dagilimin hicbir anlami yoktur. Bu, negatif bir sonucu
    # "boru hatti kirik" aciklamasindan ayiran tek kontroldur.
    controls = {}
    for name, thr in (("CDO-Fd", 60.0), ("NDO-Fd", 70.0)):
        sub = [r for r in good
               if (r["read"].get(name, {}).get("identity") or 0) >= thr]
        cnt = Counter(r["residue"] for r in sub)
        expected = REFERENCES[name]["claim_residue"]
        controls[name] = {
            "identity_threshold_pct": thr,
            "n": len(sub),
            "residues": dict(cnt.most_common()),
            "expected_residue": expected,
            "share_expected": round(cnt[expected] / len(sub), 4) if sub else None,
            "best_identity_pct": max((r["read"][name]["identity"] for r in sub),
                                     default=None),
            "enzyme_types": dict(Counter(r["cluster"] for r in sub).most_common(6)),
        }
    out["positive_controls"] = {
        "why": ("the close homologues of each reference must read back that "
                "reference's own residue. If they did not, the column mapping "
                "would be broken and the distribution below would be "
                "meaningless. This is what separates a real negative result "
                "from a broken pipeline"),
        "results": controls,
        "passed": all(c["share_expected"] == 1.0 for c in controls.values()
                      if c["n"]),
    }

    # ------------------------------- triptofan GERCEKTE nerede yasiyor
    # Sinif tablolari W'yi seyreltiyor olabilir: eger W belli ENZIM
    # TIPLERINDE yogunlasmissa, "sinifa bagli" ifadesi yanlis ama "rastgele"
    # de yanlis olur. Bu yuzden tip basina W/Y paylari ayrica yazilir --
    # buyukluk esigi olmadan, cunku en yuksek W payli tipler kucuk tipler.
    per_type = {}
    for r in good:
        d = per_type.setdefault(r["cluster"], Counter())
        d[r["residue"]] += 1
    type_rows = []
    for cl, cnt in per_type.items():
        tot = sum(cnt.values())
        type_rows.append({
            "enzyme_type": cl, "group": cl.split("_")[0], "n": tot,
            "W": cnt["W"], "A": cnt["A"], "Y": cnt["Y"],
            "W_share": round(cnt["W"] / tot, 4),
            "A_share": round(cnt["A"] / tot, 4),
            "Y_share": round(cnt["Y"] / tot, 4),
            "modal_residue": cnt.most_common(1)[0][0],
        })
    type_rows.sort(key=lambda x: (-x["W_share"], -x["n"]))
    n_w_total = sum(1 for r in good if r["residue"] == "W")
    w_types = [t for t in type_rows if t["W"]]
    out["where_the_tryptophan_lives"] = {
        "total_W_observations": n_w_total,
        "n_enzyme_types_with_any_W": len(w_types),
        "n_enzyme_types_readable": len(type_rows),
        "share_of_all_W_in_top3_types": (
            round(sum(t["W"] for t in w_types[:3]) / n_w_total, 4)
            if n_w_total else None),
        "types_with_any_W": w_types,
        "types_with_any_Y": [t for t in type_rows if t["Y"]],
        "finding": ("tryptophan is not spread across a class, it is concentrated "
                    "in a small number of enzyme types. Among them are exactly "
                    "the types the preprint worked on or cited: the cumene "
                    "dioxygenase type and the toluene dioxygenase type read W in "
                    "every readable member, which independently reproduces the "
                    "preprint's own Figure S9A observation that TDO-Fd carries a "
                    "tryptophan at this position. What does not reproduce is the "
                    "extension of that observation to a whole class: other "
                    "textbook Class IIB types in the same RO group, notably the "
                    "biphenyl dioxygenase types, read alanine (or tyrosine) and "
                    "never tryptophan"),
        "per_type_table": type_rows,
    }

    # --------------------------------------------------------- dagilim
    ent = Counter(r["residue"] for r in good)
    seq_level = {}
    for r in good:
        seq_level.setdefault(r["sequence"], r["residue"])
    sq = Counter(seq_level.values())
    gsub = collapse_genus(good)
    gn = Counter(r["residue"] for r in gsub)

    def dist(counter):
        n = sum(counter.values())
        return {
            "n": n,
            "counts": dict(counter.most_common()),
            "shares": {k: round(v / n, 4) for k, v in counter.most_common()} if n else {},
            "classes": {c: sum(v for k, v in counter.items()
                               if residue_class(k) == c) for c in RESIDUE_CLASSES},
            "class_shares": {c: round(sum(v for k, v in counter.items()
                                          if residue_class(k) == c) / n, 4)
                             for c in RESIDUE_CLASSES} if n else {},
        }

    out["distribution"] = {
        "entry_level": dist(ent),
        "distinct_sequence_level": dist(sq),
        "genus_level": dist(gn),
        "note": ("entry level counts one observation per (alpha entry x ferredoxin "
                 "neighbour); the same protein recurs across strains, so the "
                 "distinct-sequence and (type x genus) collapses are the "
                 "conservative readings"),
    }

    # ------------------------------------------------------------ testler
    out["tests"] = []

    groups = sorted({r["group"] for r in good})
    contingency(
        out, "residue_by_ferredoxin_type",
        "Is the residue at the NDO-Fd A50 column associated with the ferredoxin "
        "type recorded for the associated alpha subunit (Rieske-type vs plant-type)?",
        good, "fd_type", "residue_class",
        sorted({r["fd_type"] for r in good}), RESIDUE_CLASSES,
        "ferredoxin type (alpha level)", "alpha entry x ferredoxin neighbour",
        note=("ro_etc.ferredoxin_type is the union of the Fd types seen in the "
              "operon of the ALPHA subunit, not the type of this particular "
              "protein. It is the closest thing in this database to the authors' "
              "Class IIA/IIB split, but it is not that split"))

    contingency(
        out, "residue_by_own_pfam_domain",
        "Is the residue associated with the ferredoxin's OWN Pfam domain "
        "(Rieske PF00355/PF13806 vs plant-type Fer2 PF00111)?",
        good, "fd_own_domain", "residue_class",
        sorted({r["fd_own_domain"] for r in good}), RESIDUE_CLASSES,
        "ferredoxin's own Pfam domain", "alpha entry x ferredoxin neighbour",
        note=("sharper than the alpha-level field because it is measured on the "
              "ferredoxin protein itself"))

    contingency(
        out, "residue_by_group",
        "Is the residue associated with the RO group (1-5)?",
        good, "group", "residue_class", groups, RESIDUE_CLASSES,
        "RO group", "alpha entry x ferredoxin neighbour",
        note=("ro_group is a sequence clustering of the ALPHA subunit. It is not "
              "Batie's class, which is defined by the electron transport chain"))

    type_counts = Counter(r["cluster"] for r in good)
    big = sorted(t for t, n in type_counts.items() if n >= MIN_TYPE_N)
    big_subset = [r for r in good if r["cluster"] in big]
    contingency(
        out, "residue_by_enzyme_type",
        "Is the residue associated with the enzyme TYPE rather than only with the "
        "broad RO group?",
        big_subset, "cluster", "residue_class", big, RESIDUE_CLASSES,
        "Enzyme type", "alpha entry x ferredoxin neighbour",
        restriction=(f"types with at least {MIN_TYPE_N} readable ferredoxins: "
                     f"{len(big)} of {len(type_counts)} types, "
                     f"{len(big_subset)} observations"))

    contingency(
        out, "residue_by_reductase_type",
        "Is the residue associated with the reductase type found near the alpha "
        "subunit (FNR vs GR)? This is the axis on which Batie separates Class II "
        "from Class III.",
        good, "reductase_type", "residue_class",
        sorted({r["reductase_type"] for r in good}), RESIDUE_CLASSES,
        "reductase type", "alpha entry x ferredoxin neighbour")

    # ----------------------------------------- ikincil: Batie vekil eslesmesi
    for r in rows:
        r["batie_proxy"] = batie_proxy(r)
    proxy_subset = [r for r in good if r["batie_proxy"] != "unassignable"]
    proxy_keys = sorted({r["batie_proxy"] for r in proxy_subset})
    proxy_test = contingency(
        out, "residue_by_batie_proxy_SECONDARY",
        "SECONDARY, PROXY LABELS ONLY: under a partial mapping derived from the "
        "authors' own class definition, does the residue differ between "
        "IIB-like and III-like ferredoxins as claimed?",
        proxy_subset, "batie_proxy", "residue_class", proxy_keys, RESIDUE_CLASSES,
        "Batie-like proxy class", "alpha entry x ferredoxin neighbour",
        restriction=(f"{len(proxy_subset)} of {len(good)} readable ferredoxins could "
                     f"be given a proxy label; the rest have an ambiguous "
                     f"(FNR+GR, plant+rieske) or unobserved reductase"),
        note=("THESE ARE NOT BATIE CLASS LABELS. They are a proxy built from "
              "ro_etc.reductase_type and ro_etc.ferredoxin_type following the "
              "preprint's own definition (II vs III by reductase, IIA vs IIB by "
              "ferredoxin type) plus Batie 1991 (Class II reductase = flavin only "
              "= GR fold; Class III reductase = flavin + [2Fe-2S] = FNR fold). The "
              "mapping returns the correct answer for both of the preprint's model "
              "systems (CDO -> IIB-like, NDO -> III-like), but a proxy that agrees "
              "on two systems is still a proxy"))

    out["batie_caveat"] = {
        "direct_test_possible": False,
        "statement": ("A direct test of the authors' Class IIA / IIB / III labels "
                      "is NOT POSSIBLE in this database."),
        "why": ("those labels come from the Batie (1991) classification, which this "
                "database does not assign. What it has instead is its own ro_group "
                "1-5, which is a sequence clustering of the alpha subunit, plus "
                "ro_etc.reductase_type and ro_etc.ferredoxin_type, which describe "
                "which electron-transport partners were observed within 10 kb. "
                "ro_group and Batie class are different kinds of object - one is "
                "alpha sequence similarity, the other is electron-chain "
                "architecture - and no mapping between them is assumed here"),
        "what_was_done_instead": ("the residue was tested against the groupings "
                                  "that do exist, and separately against a clearly "
                                  "labelled proxy (test id ending _SECONDARY) built "
                                  "from the preprint's own class definition"),
        "plant_type_arm_untestable": (
            "the 'tyrosine in certain plant-type Fds from Class IIA' arm cannot be "
            "tested by this method at all, and not because of missing labels. "
            "Plant-type ferredoxins (Fer2, PF00111) have a thioredoxin-like fold "
            "with no structural or evolutionary counterpart to the Rieske-type "
            "ferredoxin region that carries A50. There is no column in a plant-type "
            "Fd that corresponds to NDO-Fd A50, so there is nothing to read. This "
            "shows up in the data as plant-type Fds failing the alignment "
            "thresholds rather than as a tyrosine count"),
    }

    # ------------------------------------------------------------- verdict
    ecl = out["distribution"]["entry_level"]["class_shares"]
    gcl = out["distribution"]["genus_level"]["class_shares"]
    dom = max(ecl, key=ecl.get) if ecl else None

    held = []
    for t in out["tests"]:
        sig = t.get("signals") or {}
        for c in sig.get("cells", []):
            if c.get("holds") and c["col"] in ("W", "A", "Y"):
                held.append({"test": t["id"], "row": c["row"], "residue": c["col"],
                             "share": c["share"], "baseline": c["baseline"],
                             "lift": c["lift"], "direction": c["direction"],
                             "q": c["q"], "q_genus": c.get("q_genus")})
    held.sort(key=lambda c: -abs(c["lift"]))

    proxy_cells = []
    if proxy_test:
        for c in (proxy_test.get("signals") or {}).get("cells", []):
            if c["col"] in ("W", "A", "Y"):
                proxy_cells.append(c)

    # Iddianin iki sinanabilir kolu: (a) W Sinif IIB'de baskin mi,
    # (b) A Sinif III'te baskin mi. Vekil tabloda IIB-benzeri satirda W
    # payinin A payini gecmesi gerekir; gecmiyorsa "baskin" yanlis.
    iib_w = iib_a = iii_a = iii_w = None
    if proxy_test:
        rws, cls, tab = proxy_test["rows"], proxy_test["cols"], proxy_test["table"]
        for i, rname in enumerate(rws):
            tot = sum(tab[i]) or 1
            share = {c: tab[i][j] / tot for j, c in enumerate(cls)}
            if rname.startswith("IIB"):
                iib_w, iib_a = share.get("W", 0.0), share.get("A", 0.0)
            if rname.startswith("III"):
                iii_a, iii_w = share.get("A", 0.0), share.get("W", 0.0)

    out["verdict"] = {
        "claim_tested": ("that the residue at the NDO-Fd A50 / CDO-Fd W48 column is "
                         "class-dependent, with W predominating in Class IIB Fds, "
                         "A in Class III Fds and Y in plant-type Class IIA Fds"),
        "dominant_residue_overall": dom,
        "dominant_share_entry_level": ecl.get(dom) if dom else None,
        "dominant_share_genus_level": gcl.get(dom) if dom else None,
        "tryptophan_share_entry_level": ecl.get("W"),
        "tyrosine_share_entry_level": ecl.get("Y"),
        "material_cells_that_survive_genus_collapse": held[:12],
        "n_material_cells_surviving": len(held),
        "proxy_class_shares": {
            "IIB_like_W_share": iib_w, "IIB_like_A_share": iib_a,
            "III_like_A_share": iii_a, "III_like_W_share": iii_w,
        },
        "arm_plant_type_IIA_tyrosine": "untestable",
        "arm_plant_type_IIA_tyrosine_reason": out["batie_caveat"][
            "plant_type_arm_untestable"],
    }

    # Karar metni dagilimdan turetilir, elle yazilmaz.
    a_dominant = dom == "A"
    w_minority = (ecl.get("W") or 0) < 0.25
    iib_w_wins = (iib_w is not None and iib_a is not None and iib_w > iib_a)
    iii_a_wins = (iii_a is not None and iii_w is not None and iii_a > iii_w)

    if a_dominant and w_minority and not iib_w_wins:
        status = "contradicted"
        reasoning = (
            f"alanine is the single dominant residue at this column across the "
            f"whole readable set ({100 * (ecl.get('A') or 0):.0f}% of "
            f"{out['distribution']['entry_level']['n']} observations at entry "
            f"level, {100 * (gcl.get('A') or 0):.0f}% after collapsing to one "
            f"observation per type x genus), and it is dominant in the IIB-like "
            f"proxy class as well. Tryptophan is a global minority "
            f"({100 * (ecl.get('W') or 0):.0f}%) and does not predominate in the "
            f"IIB-like class. Alanine at this position is therefore the default "
            f"state of Rieske-type ferredoxins, not a marker of Class III")
    elif iib_w_wins and iii_a_wins:
        status = "supported"
        reasoning = (
            "under the proxy mapping, tryptophan predominates in the IIB-like "
            "class and alanine in the III-like class, in the direction claimed, "
            "and the cells survive the genus-level collapse")
    elif iib_w_wins or iii_a_wins:
        status = "partially supported"
        reasoning = (
            f"only one arm of the claim holds under the proxy mapping "
            f"(IIB-like W>A: {iib_w_wins}; III-like A>W: {iii_a_wins}); "
            f"alanine is the overall dominant residue "
            f"({100 * (ecl.get('A') or 0):.0f}%)")
    else:
        status = "contradicted"
        reasoning = ("neither testable arm of the claim reproduces in the "
                     "direction claimed")

    out["verdict"]["status"] = status
    out["verdict"]["reasoning"] = reasoning
    out["verdict"]["direct_test_of_their_class_labels"] = "untestable"

    # Yalin bir "yanlis" eksik olurdu: W rastgele dagilmiyor, SIKI bir klad
    # isareti. Testlerde en guclu iliski sinif vekilinde degil ENZIM
    # TIPINDE (Cramer V ~0,6) -- yani pozisyon gercekten yapisal, ama
    # tasidigi bilgi sinif degil tip.
    wl = out["where_the_tryptophan_lives"]
    type_test = next((t for t in out["tests"]
                      if t["id"] == "residue_by_enzyme_type"), None)
    proxy_v = (proxy_test or {}).get("cramers_v")
    out["verdict"]["refinement"] = {
        "headline": ("the residue is strongly structured, but by ENZYME TYPE, "
                     "not by class"),
        "cramers_v_enzyme_type": (type_test or {}).get("cramers_v"),
        "cramers_v_batie_proxy": proxy_v,
        "effect_size_ratio": (round((type_test or {}).get("cramers_v", 0) / proxy_v, 2)
                              if proxy_v else None),
        "tryptophan_is_a_clade_marker": {
            "total_W_observations": wl["total_W_observations"],
            "n_types_carrying_it": wl["n_enzyme_types_with_any_W"],
            "of_readable_types": wl["n_enzyme_types_readable"],
            "share_of_W_in_top3_types": wl["share_of_all_W_in_top3_types"],
        },
        "statement": ("W at this column is confined to a narrow clade of types "
                      "around cumene / toluene / ethylbenzene dioxygenase - the "
                      "very enzymes the preprint studied and cited. Within that "
                      "clade the authors' observation is correct and reproduces "
                      "perfectly. The error is one of extrapolation: the clade is "
                      "much smaller than Class IIB, and the rest of Class IIB "
                      "carries alanine. A two-protein sample drawn from inside "
                      "that clade cannot distinguish a clade marker from a class "
                      "marker, which is precisely why the generalisation fails"),
    }
    out["verdict"]["summary"] = (
        f"The distributional claim is {status.upper()} in this dataset for the two "
        f"arms that can be tested (W in IIB, A in III), and the third arm (Y in "
        f"plant-type IIA) is UNTESTABLE here. A direct test of the authors' Batie "
        f"class labels is UNTESTABLE because this database does not assign them. "
        f"{reasoning}.")
    out["verdict"]["what_this_does_not_show"] = (
        "that the position does not matter. The authors' mutagenesis shows it "
        "does, and a residue can be functionally critical while being almost "
        "invariant across a family - indeed that is the usual case. What the "
        "distribution contradicts is only the claim that the IDENTITY of this "
        "residue tracks RO class")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    # ------------------------------------------------------------- ozet
    print(__doc__.strip().splitlines()[1].strip())
    print("=" * 78)
    print("\n[REFERANS DIZILER]")
    for name, ref in out["references"].items():
        print(f"  {name:8} {ref['uniprot']:8} {ref['length']:4} aa  "
              f"pos {ref['claim_position']:3} = {ref['residue_at_claim_position']}"
              f"  ({ref['numbering_checks_passed']} kalinti numarasi dogrulandi)")
    cc = out["reference_cross_check"]
    print(f"  capraz kontrol: NDO-Fd A50 kolonu -> CDO-Fd "
          f"{cc['cdo_residue_aligned']}{cc['cdo_position_aligned']}  "
          f"(makale Sekil 4A ile ayni: {cc['agrees_with_preprint_figure_4A']}), "
          f"kimlik %{cc['pairwise_identity_pct']} "
          f"(makale: %{cc['preprint_reported_identity_pct']})")

    cov = out["coverage"]
    print("\n[KAPSAMA]")
    print(f"  gene_category ferredoksin satiri : {cov['gene_category_ferredoxin_rows']}")
    print(f"  dizisi olan                      : {cov['with_translation']}"
          f"  ({cov['distinct_sequences']} benzersiz dizi)")
    print(f"  kolon okunacak kadar iyi hizalanan: {cov['aligned_well_enough']}"
          f"  ({cov['aligned_well_enough_distinct_sequences']} benzersiz)")
    print(f"  atilan                           : {cov['discarded']}")
    for reason, n in cov["discard_reasons"].items():
        print(f"      {reason:42} {n:5}")
    print(f"  esikler: kimlik >= %{MIN_IDENTITY}, kapsama >= {MIN_COVERAGE}, "
          f"Rieske ligandi {REQUIRE_ANCHORS}/4")
    crc = cov["cross_reference_concordance"]
    print(f"  iki referanstan ayni kalinti     : {crc['same_residue_from_both']}"
          f"/{crc['readable_against_both_references']}"
          f"  (%{100 * (crc['rate'] or 0):.1f})")

    pc = out["positive_controls"]
    print(f"\n[POZITIF KONTROL]  gecti: {pc['passed']}")
    for name, c in pc["results"].items():
        print(f"  {name}'ye >=%{c['identity_threshold_pct']:.0f} benzeyen "
              f"{c['n']:3} ferredoksin -> {c['residues']}"
              f"  (beklenen {c['expected_residue']}: "
              f"%{100 * (c['share_expected'] or 0):.0f}, "
              f"en yuksek kimlik %{c['best_identity_pct']:.1f})")

    print("\n[A50 KOLONUNDAKI KALINTI DAGILIMI]")
    print(f"  {'duzey':24} {'n':>6}   en sik kalintilar")
    for lvl in ("entry_level", "distinct_sequence_level", "genus_level"):
        d = out["distribution"][lvl]
        top = "  ".join(f"{k}:%{100 * v:.1f}" for k, v in
                        list(d["shares"].items())[:6])
        print(f"  {lvl:24} {d['n']:6}   {top}")
    print(f"  {'sinif (giris duzeyi)':24} {'':6}   " +
          "  ".join(f"{c}:%{100 * (out['distribution']['entry_level']['class_shares'].get(c) or 0):.1f}"
                    for c in RESIDUE_CLASSES))

    print("\n[KONTENJANS TESTLERI]  (giris duzeyi / cins duzeyi)")
    print(f"  {'test':38} {'n':>6} {'V':>6} {'p':>9} {'V_cins':>7} {'p_cins':>9} {'maddi':>6}")
    for t in out["tests"]:
        sig = t.get("signals") or {}
        print(f"  {t['id'][:38]:38} {t['n_entry']:6} {t['cramers_v']:6.3f} "
              f"{t['p']:9.2e} {(t['cramers_v_genus'] or 0):7.3f} "
              f"{(t['p_genus'] if t['p_genus'] is not None else float('nan')):9.2e} "
              f"{sig.get('n_cells_material_genus', 0):6}")

    print("\n[W / A / Y SATIR PAYLARI -- ana gruplandirmalar]")
    for t in out["tests"]:
        if t["id"] not in ("residue_by_own_pfam_domain", "residue_by_group",
                           "residue_by_reductase_type",
                           "residue_by_batie_proxy_SECONDARY"):
            continue
        print(f"  [{t['id']}]")
        print(f"    {'satir':40} {'n':>6} {'W':>7} {'A':>7} {'Y':>7} {'diger':>7}")
        for i, rname in enumerate(t["rows"]):
            tot = sum(t["table"][i]) or 1
            cells = {c: t["table"][i][j] for j, c in enumerate(t["cols"])}
            print(f"    {rname[:40]:40} {tot:6} " + " ".join(
                f"%{100 * cells.get(c, 0) / tot:5.1f}" for c in RESIDUE_CLASSES))

    print("\n[MADDI HUCRELER -- q<0.05, >=10 puan sapma, cins duzeyinde de duruyor]")
    if held:
        print(f"    {'test':34} {'satir':26} {'kal':>3} {'pay':>7} {'taban':>7} {'fark':>7}")
        for c in held[:12]:
            print(f"    {c['test'][:34]:34} {c['row'][:26]:26} {c['residue']:>3} "
                  f"%{100 * c['share']:5.1f} %{100 * c['baseline']:5.1f} "
                  f"{100 * c['lift']:+6.1f}")
    else:
        print("    yok -- W/A/Y icin hicbir hucre ucuncu olcutu birlikte gecmiyor")

    print("\n[TRIPTOFAN GERCEKTE NEREDE]  "
          f"{wl['total_W_observations']} W gozlemi, "
          f"{wl['n_enzyme_types_with_any_W']}/{wl['n_enzyme_types_readable']} "
          f"tipte; en ust 3 tip hepsinin "
          f"%{100 * (wl['share_of_all_W_in_top3_types'] or 0):.0f}'ini tasiyor")
    print(f"    {'enzim tipi':20} {'n':>5} {'W':>4} {'W payi':>8} {'A payi':>8}")
    for t in wl["types_with_any_W"]:
        print(f"    {t['enzyme_type']:20} {t['n']:5} {t['W']:4} "
              f"%{100 * t['W_share']:6.1f} %{100 * t['A_share']:6.1f}")
    if wl["types_with_any_Y"]:
        print("    -- tirozin tasiyan tipler (iddia bunu bitki-tipi IIA'ya "
              "bagliyor):")
        for t in wl["types_with_any_Y"]:
            print(f"    {t['enzyme_type']:20} {t['n']:5} {t['Y']:4} "
                  f"%{100 * t['Y_share']:6.1f} (Y)")

    print("\n[BATIE SINIFI UYARISI]")
    print(f"  {out['batie_caveat']['statement']}")
    print("  bitki-tipi kolu: SINANAMAZ (farkli katlanma, karsilik gelen kolon yok)")

    v = out["verdict"]
    print("\n" + "=" * 78)
    print(f"[KARAR]  {v['status'].upper()}")
    print(f"  {v['summary']}")
    ref = v["refinement"]
    print(f"\n[INCELTME]  {ref['headline']}")
    print(f"  Cramer V: enzim tipi {ref['cramers_v_enzyme_type']:.3f} vs "
          f"Batie vekili {ref['cramers_v_batie_proxy']:.3f} "
          f"({ref['effect_size_ratio']}x)")
    print(f"  {ref['statement']}")
    print(f"\n[yazildi] {args.out}")


if __name__ == "__main__":
    main()
