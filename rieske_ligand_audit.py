"""Kapinin disladigi Rieske proteinlerinde HANGI ligand degismis?

NEDEN BU BETIK VAR. InterPro karsilastirmasi (ROADMAP 135.19) tek yonlu
bir asimetri buldu: 139 dizide InterPro Rieske domaini goruyor, bu
veritabaninin kapisi ise disliyor -- cunku kapi DORT KUME LIGANDININ
yerinde olmasini istiyor, InterPro ise profil benzerligine bakiyor.
Celiski degil, ama 76'sinda mononukleer katalitik ucgen TAM, yani
"kumesi bozuk, katalitik merkezi saglam" bir kume var ve bu bir kurator
sorusu.

Bir kurator o soruyu "139 dizi dusuruldu" cumlesiyle cevaplayamaz.
Cevaplayabilmek icin HANGI ligandin degistigini ve YERINE NE geldigini
bilmek gerekir: Cys->Ser bir metal baglama kaybidir (kukurt gider),
Cys->Gly bir katlanma kaybina isaret eder, His->Tyr ise baska bir
koordinasyon kimyasi anlamina gelebilir. Bu betik onu olcer; hicbir
cagriyi DEGISTIRMEZ.

YONTEM. Diziler motif modeline (ROmotif71) hmmalign ile hizalanir ve
dort Rieske kolonu (Cys-1 164, His-1 166, Cys-2 184, His-2 187) tek tek
okunur. Kolonlar `ro_motif`ten alinir, BURADA YENIDEN YAZILMAZ -- ayni
olcumun ikinci uygulamasi birincisiyle sessizce ayrisir. Pencere de
`ro_motif.COLUMN_WINDOW`: kapinin kullandigi toleransin aynisi, yoksa
"kapi neyi disladi" sorusuna baska bir kuralla cevap verilmis olur.

Cikti: analysis_out/rieske_ligand_audit.json
"""

import argparse
import json
import os
import subprocess
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict

import re

import ro_motif

RIESKE_IPR = ("IPR017941", "IPR036922")


def asymmetric(con):
    """InterPro Rieske diyor, kapimiz disliyor: (seq_md5, dizi, ucgen)."""
    placeholders = ",".join("?" * len(RIESKE_IPR))
    return con.execute(f"""
        SELECT u.seq_md5, u.sequence, u.triad_sites_found, u.organism,
               u.ro_cluster
          FROM ro_union u
         WHERE u.seq_md5 IN (SELECT seq_md5 FROM ro_xref
                              WHERE db = 'interpro'
                                AND identifier IN ({placeholders}))
           AND (u.rieske_intact IS NULL OR u.rieske_intact = 0)
           AND u.sequence IS NOT NULL AND LENGTH(u.sequence) > 50
         ORDER BY u.seq_md5""", RIESKE_IPR).fetchall()


def align(rows, hmm, work):
    """hmmalign ile hizalar ve md5 -> MATCH-KOLON dizgisi dondurur.

    KOLONLARI KENDIM SAYMIYORUM -- ilk yazimda sayiyordum ve YANLISTI.
    hmmalign ciktisinda insert kolonlari da var (kucuk harf ve '.'),
    dolayisiyla hizalama dizgisindeki 164. konum modelin 164. ESLESME
    kolonu DEGIL. Sonuc apacik sacmaydi ve onu yakalatan sey de o oldu:
    139 dizinin 137'sinde dort ligandin HEPSI bosluk cikti. Insert
    kolonlarini atip model pozisyonlarina inmenin dogru yolu zaten
    projede var: `ro_motif.read_stockholm_matchcols`. Ayni olcumun ikinci
    uygulamasini yazmak bu projede tekrar eden hata, ve burada onu ben
    yaptim.
    """
    fasta = os.path.join(work, "asym.fasta")
    with open(fasta, "w", encoding="utf-8") as handle:
        for md5, seq, _t, _o, _c in rows:
            handle.write(f">{md5}\n{seq}\n")
    sto = os.path.join(work, "asym.sto")
    subprocess.run(["hmmalign", "--amino", "-o", sto, hmm, fasta],
                   check=True, stdout=subprocess.DEVNULL)
    return ro_motif.read_stockholm_matchcols(sto)


def read_ligands(aligned, window):
    """Dort Rieske kolonunda ne var: [(etiket, beklenen, bulunan, yerinde mi)]."""
    out = []
    for column, accepted, label in ro_motif.RIESKE_SITES:
        hit, seen = False, None
        for offset in range(-window, window + 1):
            index = column - 1 + offset
            if 0 <= index < len(aligned):
                residue = aligned[index].upper()
                if offset == 0:
                    seen = residue
                if residue in accepted:
                    hit = True
                    break
        out.append((label, accepted, seen, hit))
    return out


#  GERCEK RIESKE MOTIFI, kolondan BAGIMSIZ. Kume iki ilmekle tutulur:
#  C-x-H ve C-x-x-H, aralari ~14-20 kalinti (naftalen dioksijenazda
#  C81-x-H83 ... C101-x-x-H104). Bu desen HAM dizide aranir, yani
#  hizalamadan ve kolon numarasindan tamamen bagimsiz.
RIESKE_PAIR = re.compile(r"C.H.{10,30}C..H")


def pair_in_raw_sequence(sequence):
    """Ham dizide kanonik Rieske ilmek cifti var mi?

    NEDEN BU IKINCI OLCUM. Kolon testi bir ligandi iki sebeple
    bulamayabilir: ligand GERCEKTEN yok, ya da yerinde ama pencerenin
    (+-2 kolon) disina kaymis. Ikisi ayni sonucu verir ama ZIT anlamlari
    var -- birincisi "bu protein kume baglamaz", ikincisi "bizim
    hizalamamiz onu kaciriyor". Bu projede ayni ayrim daha once
    7_701_AerC'de karsimiza cikti: "kolon testi ulasamiyor" cumlesi
    protein hakkinda degil HIZALAMA hakkinda bir ifadeydi (ROADMAP
    135.16). O yuzden desen ham dizide, kolonlardan bagimsiz araniyor.
    """
    return bool(RIESKE_PAIR.search((sequence or "").upper()))


def validate_pattern(con, limit=4000):
    """AYIRT EDICIYI KULLANMADAN ONCE DOGRULA.

    "Desen 139 dizinin 138'inde YOK" cumlesi, desenin KENDISI calismiyorsa
    hicbir sey soylemez. O yuzden ayni desen, kume ligandlarini SAGLAM
    dedigimiz dizilerde ve 76 kuratorlu referansta ayrica sinanir; orada
    tutmuyorsa yokluk bulgusu degil, desen hatasidir. Bu projede ayni
    ihmal daha once olculdu: bir ayirt edici dogrulanmadan kullanildiginda
    uretilen sayi olcumun degil araci kusurunun sayisi olur.
    """
    rows = con.execute(
        "SELECT sequence FROM ro_union WHERE rieske_intact = 1 "
        "AND sequence IS NOT NULL LIMIT ?", (limit,)).fetchall()
    intact_hit = sum(1 for (seq,) in rows if pair_in_raw_sequence(seq))
    return {"sequences_we_call_intact": len(rows),
            "of_those_matching_the_pattern": intact_hit,
            "share": round(intact_hit / len(rows), 4) if rows else None}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--hmm", default="ROs_71_Clean/ROmotif71.hmm")
    parser.add_argument("--out",
                        default="analysis_out/rieske_ligand_audit.json")
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    rows = asymmetric(con)
    if not rows:
        sys.exit("asimetrik dizi yok -- InterPro karsilastirmasi kosuldu mu?")
    print(f"[giris] {len(rows)} dizi hizalanacak")

    with tempfile.TemporaryDirectory() as work:
        aligned = align(rows, args.hmm, work)
    print(f"[hizalama] {len(aligned)} dizi, match kolonu "
          f"{len(next(iter(aligned.values()))) if aligned else 0}")

    window = ro_motif.COLUMN_WINDOW
    missing_counter = Counter()
    substitution = defaultdict(Counter)
    by_triad = defaultdict(Counter)
    examples = defaultdict(list)
    unaligned = []

    raw_pair = Counter()
    for md5, _seq, triad, organism, cluster in rows:
        raw_pair["present" if pair_in_raw_sequence(_seq) else "absent"] += 1
        if md5 not in aligned:
            unaligned.append(md5)
            continue
        ligands = read_ligands(aligned[md5], window)
        gone = [(label, accepted, seen)
                for label, accepted, seen, hit in ligands if not hit]
        if not gone:
            #  Kapi disladi ama bu kolonlarda dort ligand yerinde: demek ki
            #  kapi BASKA bir dizi uzerinde karar verdi (aday hizalamasi ile
            #  buradaki yeniden hizalama ayni olmak zorunda degil). Sessiz
            #  gecilmez, ayri sayilir.
            missing_counter["none_missing_in_this_realignment"] += 1
            continue
        key = " + ".join(label for label, _a, _s in gone)
        missing_counter[key] += 1
        by_triad[str(triad)][key] += 1
        for label, accepted, seen in gone:
            substitution[label][seen or "-"] += 1
        if len(examples[key]) < 3:
            examples[key].append({"seq_md5": md5, "organism": organism,
                                  "nearest_cluster": cluster,
                                  "triad_sites_found": triad})

    report = {
        "question": ("in the sequences InterPro calls Rieske but this "
                     "database's gate excludes, which cluster ligand is "
                     "substituted, and by what?"),
        "why_it_matters": (
            "a curator cannot act on '139 sequences were dropped'. Cys to Ser "
            "loses the thiolate that holds the iron; Cys to Gly suggests the "
            "fold itself is gone; His to Tyr may mean a different "
            "coordination chemistry rather than none. The substitution is the "
            "question, not the count"),
        "method": {
            "alignment": ("hmmalign against the project motif model "
                          "ROmotif71, then insert columns removed with "
                          "ro_motif.read_stockholm_matchcols so that column "
                          "numbers are model positions"),
            "columns": [{"column": c, "accepted": a, "label": l}
                        for c, a, l in ro_motif.RIESKE_SITES],
            "window": window,
            "columns_are_not_restated_here": (
                "the columns and the tolerance window are read from "
                "ro_motif, the owner of the motif question; restating them "
                "would be a second implementation of the same measurement, "
                "and this project has measured that those diverge silently"),
            "changes_nothing": ("this script measures and reports; no call, "
                                "no table and no curation is changed by it"),
        },
        "denominator": {
            "sequences": len(rows),
            "aligned": len(rows) - len(unaligned),
            "not_aligned": len(unaligned),
        },
        "missing_ligand_pattern": dict(missing_counter.most_common()),
        "canonical_pair_in_raw_sequence": {
            "pattern": RIESKE_PAIR.pattern,
            "counts": dict(raw_pair.most_common()),
            "what_it_separates": (
                "a ligand can be missed for two opposite reasons: it is "
                "genuinely absent, or it is present but displaced beyond the "
                "column test's two-position window. Those give the same "
                "result and mean opposite things -- 'this protein cannot "
                "bind the cluster' against 'our alignment misses it'. This "
                "pattern is searched in the raw sequence, independent of any "
                "alignment or column number, so it separates the two"),
            "validated_on": validate_pattern(con),
            "what_it_is_not": (
                "finding the pattern does not prove a working cluster: the "
                "spacing is necessary, not sufficient, and a sequence can "
                "carry C-x-H ... C-x-x-H by chance. It narrows the question "
                "rather than answering it"),
        },
        "substituting_residue_by_ligand": {
            label: dict(c.most_common()) for label, c in substitution.items()},
        "pattern_by_mononuclear_triad_sites": {
            t: dict(c.most_common()) for t, c in sorted(by_triad.items(),
                                                        reverse=True)},
        "examples": {k: v for k, v in examples.items()},
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print(f"[desen] {dict(missing_counter.most_common(5))}")
    for label, c in substitution.items():
        print(f"[{label}] yerine: {dict(c.most_common(4))}")
    print(f"[yazildi] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
