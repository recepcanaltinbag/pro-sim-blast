"""Her kaydin NE OLDUGUNU siniflandirir: genom mu, transkript mi.

NEDEN BU DOSYA VAR. `replicon` tablosunun adi her satirin bir replikon
oldugunu soyluyor, ama soylemiyor: dogrulanmis uyesi olan 11.610 kaydin
958'i `XM_`/`NM_` ile baslayan bir mRNA transkripti. Bir transkript
replikon DEGIL -- kromozom ya da plazmit degil, bir genin kopyasi. Bu,
bu projede daha once olculmus hata sinifinin aynisi: HAYALET REPLIKON
satirinda tablo, biriminin ne oldugu konusunda yanlis bir sey iddia
ediyordu (ROADMAP 135.2). Burada iddia daha sessiz ama daha pahali,
cunku bu satirlar REPLIKON BASINA paydalarda duruyor: bir transkriptin
"komsu geni" ya da "tam operonu" olamaz, dolayisiyla her baglam
istatistigi 958 kayit boyunca sifir toplayip paydada kaliyor.

SATIRLAR SILINMEZ. Bu 958 kaydin DIZILERI gercek ve dogrulanmis uye
uretiyorlar; silmek 958 girisi kaybetmek olurdu. Yapilan sey, kaydin
turunu YAZMAK ve baglam paydalarindan cikarmak.

AYIRT EDICI, KULLANILMADAN ONCE DOGRULANIR. Akla gelen ilk olcut
accession onekidir (`XM_` -> mRNA). Ama bir onek bir dizgi; bu projede
bir karar biriminin dizgi isleminden turetilmesi tam olarak kacinilan
hata (bkz. ro_taxonomy). O yuzden tur, kaydin KENDI beyanindan okunur:
GenBank LOCUS satiri molekul tipini yazar --

    LOCUS       NM_001002612            1714 bp    mRNA    linear ...
    LOCUS       NC_002033             184457 bp    DNA     circular ...

-- ve onek yalnizca CAPRAZ KONTROL olarak kullanilir. Ikisi ayrisirsa
satir `disagreement` olarak sayilir ve siniflandirilmaz; sessizce bir
tarafa yazilmaz.

Cikti: replicon_kind tablosu + analysis_out/replicon_kind.json
"""

import argparse
import gzip
import json
import os
import re
import sqlite3
from collections import Counter

LOCUS = re.compile(rb"^LOCUS\s+(\S+)\s+(\d+)\s+bp\s+(\S+)", re.MULTILINE)

#  NCBI'nin KENDI accession oneki sozlugu (RefSeq ve GenBank). Bu bir
#  tahmin degil, belgelenmis bir kontrollu sozluk; ama yine de yalnizca
#  capraz kontrol icin kullanilir.
PREFIX_KIND = {
    "NC": "genome", "NZ": "genome", "AC": "genome", "NG": "genome",
    "NT": "genome", "NW": "genome", "NS": "genome", "CP": "genome",
    "XM": "transcript", "NM": "transcript",
    "XR": "transcript", "NR": "transcript",
}

#  LOCUS'un molekul tipi alani. GenBank'in kendi kumesi.
MOLTYPE_KIND = {
    "DNA": "genome", "ss-DNA": "genome", "ds-DNA": "genome",
    "mRNA": "transcript", "RNA": "transcript", "ss-RNA": "transcript",
    "cRNA": "transcript", "tRNA": "transcript", "rRNA": "transcript",
}


def open_record(gbk_dir, nucleotide_id):
    """Kaydin dosyasini acar. DIZIN IKI BICIMI BIR ARADA TUTUYOR: olculdu,
    17.074 dosyanin bir kismi duz `.gbk`, bir kismi `.gbk.gz`. Yalnizca
    birine bakmak 2.982 kaydi "okunamadi" diye isaretliyordu -- yani
    dosya duzeninden gelen bir eksik, kaydin kendi ozelligi gibi
    gorunuyordu."""
    plain = os.path.join(gbk_dir, nucleotide_id + ".gbk")
    if os.path.exists(plain):
        return open(plain, "rb")
    packed = plain + ".gz"
    if os.path.exists(packed):
        return gzip.open(packed, "rb")
    return None


def read_locus(handle):
    """LOCUS satirindan (accession, uzunluk, moltype). Dosyanin basi yeter."""
    head = handle.read(4096)
    match = LOCUS.search(head)
    if not match:
        return None
    return (match.group(1).decode(), int(match.group(2)),
            match.group(3).decode())


def classify(nucleotide_id, gbk_dir):
    base = nucleotide_id.split(".")[0]
    prefix = base.split("_")[0] if "_" in base else base[:2]
    by_prefix = PREFIX_KIND.get(prefix)

    handle = open_record(gbk_dir, nucleotide_id)
    if handle is None:
        locus = None
    else:
        with handle:
            locus = read_locus(handle)
    moltype = locus[2] if locus else None
    by_moltype = MOLTYPE_KIND.get(moltype) if moltype else None

    if by_moltype is None:
        #  Kaydin kendi beyani okunamadi: ONEKLE KARAR VERILMEZ. Bir
        #  dizgiden turetilen bir karar, olculmemis bir karardir.
        return (prefix, moltype, "unreadable", by_prefix)
    if by_prefix is not None and by_prefix != by_moltype:
        return (prefix, moltype, "disagreement", by_prefix)
    return (prefix, moltype, by_moltype, by_prefix)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--gbk", default="gbk_files")
    parser.add_argument("--out", default="analysis_out/replicon_kind.json")
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    nucs = [r[0] for r in con.execute(
        "SELECT nucleotide_id FROM replicon ORDER BY nucleotide_id")]

    rows, kinds, moltypes, unexpected = [], Counter(), Counter(), []
    for nuc in nucs:
        prefix, moltype, kind, by_prefix = classify(nuc, args.gbk)
        rows.append((nuc, prefix, moltype, kind))
        kinds[kind] += 1
        moltypes[str(moltype)] += 1
        if kind in ("disagreement", "unreadable"):
            unexpected.append({"nucleotide_id": nuc, "prefix": prefix,
                               "moltype": moltype, "prefix_says": by_prefix,
                               "outcome": kind})

    con.execute("DROP TABLE IF EXISTS replicon_kind")
    con.execute("""CREATE TABLE replicon_kind (
                       nucleotide_id TEXT PRIMARY KEY,
                       accession_prefix TEXT,
                       moltype TEXT,
                       record_kind TEXT)""")
    con.executemany("INSERT INTO replicon_kind VALUES (?,?,?,?)", rows)
    con.commit()

    #  BEDEL: transkript kayitlari baglam paydalarindan cikarilinca payda
    #  ne kadar kuculuyor, ve o kayitlarin kac dogrulanmis uyesi var.
    #  Ikisi ayri sayilir: payda duzeltmesi bir KAZANC, kaybedilen giris
    #  olmadiginin kaniti.
    effect = dict(con.execute("""
        SELECT k.record_kind, COUNT(DISTINCT r.nucleotide_id)
          FROM replicon_kind k JOIN ro r USING(nucleotide_id)
         WHERE r.is_confirmed = 1
         GROUP BY 1""").fetchall())
    entries = dict(con.execute("""
        SELECT k.record_kind, COUNT(*)
          FROM replicon_kind k JOIN ro r USING(nucleotide_id)
         WHERE r.is_confirmed = 1
         GROUP BY 1""").fetchall())

    report = {
        "question": "Is every row of the replicon table actually a replicon?",
        "answer": ("No: records whose own LOCUS line declares an mRNA molecule "
                   "are transcripts, not replicons. They are kept, because "
                   "their sequences are real and recruit confirmed members, "
                   "but they are not a per-replicon denominator."),
        "method": {
            "discriminator": ("the record's own declaration, the GenBank LOCUS "
                              "molecule type field"),
            "why_not_the_accession_prefix": (
                "an accession prefix is a string, and in this project a "
                "decision unit derived from a string operation is an unmeasured "
                "decision (see ro_taxonomy). The prefix is used only as a "
                "cross-check; where the two disagree the row is left "
                "unclassified rather than silently assigned."),
            "cross_check": ("NCBI's own documented accession prefix vocabulary "
                            "compared against the LOCUS molecule type"),
        },
        "record_kinds": dict(kinds.most_common()),
        "moltypes": dict(moltypes.most_common()),
        "prefix_and_moltype_disagree": len(unexpected),
        "unclassified_rows": unexpected[:20],
        "confirmed_members_by_record_kind": {
            "records": effect, "entries": entries},
        "what_changes": ("per-replicon context statistics take their "
                         "denominator from record_kind='genome'; a transcript "
                         "cannot have a neighbouring gene or a complete "
                         "operon, so counting it in the denominator reports a "
                         "measurement failure as a biological absence."),
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print(f"[tur]   {dict(kinds.most_common())}")
    print(f"[uyusmazlik] {len(unexpected)}")
    print(f"[dogrulanmis] kayit={effect} giris={entries}")
    print(f"[yazildi] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
