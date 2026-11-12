"""Kendi domain cagrimizi InterPro'nun KENDI modellerine karsi sinar.

NEDEN BU BETIK VAR. Bu veritabaninin en cok yuk tasiyan cagrisi "bu dizi
bir Rieske proteini"dir: aday havuzuna girmenin kosulu o, ve ondan sonraki
her sey (tip atamasi, operon, kimya) o cagrinin uzerine kuruluyor. Cagri
projenin KENDI modellerinden geliyor. 2026-11-13'te InterPro
tanimlayicilari cekilince (ROADMAP 135.18) ilk kez BAGIMSIZ bir ikinci
gorus elde edildi: InterPro ayni dizileri kendi profilleriyle anotluyor.

BAGIMSIZLIK TAM OLARAK NEREYE KADAR -- ve bu cumle onemli.
  BAGIMSIZ OLAN: bir dizi VERILDIGINDE, InterPro'nun ona Rieske domaini
    koyup koymadigi. O profiller baska insanlar tarafindan, baska
    verilerle, baska esiklerle kuruldu.
  BAGIMSIZ OLMAYAN: hangi dizilerin SORULDUGU. Aday havuzu bu projenin
    kendi Rieske/motif modelleriyle secildi, yani InterPro hic
    koymadigimiz bir dizi hakkinda bize itiraz EDEMEZ. Dolayisiyla
    buradaki uyum oranı KESINLIK (precision) hakkinda bir sey soyler,
    DUYARLILIK (recall) hakkinda HICBIR SEY soylemez: kacirdigimiz
    diziler bu olcumun disinda.
  OLCULMEYEN: InterPro'nun anotlamadigi diziler. 73.299 birlesik dizinin
    yalnizca 25.313'u InterPro'da bir kayit tasiyor; geri kalan 47.986
    CELISMIYOR, ANOTSUZ. Bir kanit yoklugu, yokluk kaniti degildir.

Cikti: analysis_out/interpro_agreement.json
"""

import argparse
import json
import os
import sqlite3
import sys
from collections import Counter, defaultdict

#  InterPro'nun Rieske kayitlari. Adlar kaynaktan alindi (InterPro API,
#  2026-11-13), tahmin edilmedi.
RIESKE = {
    "IPR017941": "Rieske [2Fe-2S] iron-sulphur domain",
    "IPR036922": "Rieske [2Fe-2S] iron-sulphur domain superfamily",
}

#  Alfa alt birimi / katalitik domain kayitlari. Bunlar Rieske'den AYRI
#  bir soru: Rieske domaini bir elektron aktarim merkezi, bu kayitlar
#  halka hidroksilleyen aile uyeligi.
ALPHA = {
    "IPR001663": "Aromatic-ring-hydroxylating dioxygenase, alpha subunit",
    "IPR015879": ("Aromatic-ring-hydroxylating dioxygenase, alpha subunit, "
                  "C-terminal domain"),
    "IPR015881": ("Aromatic-ring-hydroxylating dioxygenase ARHD, "
                  "Rieske 2Fe-2S-binding site"),
    "IPR043266": ("Ring-hydroxylating dioxygenase, alpha subunit NdoB-like, "
                  "C-terminal"),
}

#  Islevsel olarak AYRI aileler. Bunlarin varligi bir hata degil: bu
#  veritabani Rieske oksijenazlarin TAMAMINI topluyor, yalnizca halka
#  hidroksilleyenleri degil.
OTHER_FAMILIES = {
    "IPR050584": "Cholesterol 7-desaturase",
    "IPR044043": "Vanillate O-demethylase oxygenase-like, C-terminal catalytic",
    "IPR045605": ("3-ketosteroid-9-alpha-monooxygenase, oxygenase "
                  "component-like, C-terminal"),
    "IPR013626": "Pheophorbide a oxygenase",
    "IPR045623": "LigXa-like, C-terminal domain",
}


def load(con):
    """seq_md5 -> (bizim kaydimiz, InterPro kayitlari)."""
    ours = {}
    for md5, rieske, group, cluster, gate in con.execute("""
            SELECT seq_md5, rieske_intact, ro_group, ro_cluster, passes_gate
              FROM ro_union"""):
        ours[md5] = {"rieske_intact": rieske, "group": group,
                     "cluster": cluster, "passes_gate": gate}
    theirs = defaultdict(set)
    for md5, ident in con.execute(
            "SELECT seq_md5, identifier FROM ro_xref WHERE db = 'interpro'"):
        theirs[md5].add(ident)
    return ours, theirs


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out", default="analysis_out/interpro_agreement.json")
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    ours, theirs = load(con)
    annotated = [m for m in ours if m in theirs]

    #  1. RIESKE UYUMU. Birim dizi; payda "InterPro'da bir kaydi olan
    #  dizi", TUM diziler degil.
    with_rieske = [m for m in annotated
                   if theirs[m] & set(RIESKE)]
    ours_rieske = [m for m in annotated if ours[m]["rieske_intact"]]
    both = [m for m in ours_rieske if m in set(with_rieske)]
    ours_only = [m for m in ours_rieske if m not in set(with_rieske)]
    theirs_only = [m for m in with_rieske if not ours[m]["rieske_intact"]]

    #  2. AILE DAGILIMI. Bizim grubumuz ile InterPro'nun aile kaydi.
    #  Bu bir "dogru/yanlis" tablosu DEGIL: bu veritabani Rieske
    #  oksijenazlarin tamamini topluyor ve InterPro'nun ayri aile adlari
    #  (kolesterol 7-desaturaz, feoforbid a oksijenaz) tam olarak o
    #  cesitliligi gosteriyor.
    by_group = defaultdict(Counter)
    for md5 in annotated:
        group = ours[md5]["group"] or "?"
        for ident in theirs[md5]:
            if ident in ALPHA:
                by_group[group]["alpha_subunit_family"] += 1
            elif ident in OTHER_FAMILIES:
                by_group[group][OTHER_FAMILIES[ident]] += 1
        by_group[group]["_sequences"] += 1

    #  3. ASIMETRININ KENDISI. InterPro'nun Rieske dedigi ama bizim
    #  "saglam" demedigimiz diziler. Bu bir CELISKI DEGIL: InterPro'nun
    #  olcusu profil benzerligi (domain var mi), bizimki DORT LIGANDIN
    #  yerinde olmasi. Kumesiz bir Rieske katlanmasi kimyayi yapamaz, o
    #  yuzden kapi ligandlari istiyor. Ama bu diziler SESSIZCE dusmemeli:
    #  76'sinda mononukleer ucgenin ucu de tam, yani kuratorun bakmasi
    #  gereken bir kume (kume ligandi degismis, katalitik ucgen saglam).
    detail = con.execute("""
        SELECT u.triad_sites_found, COUNT(*)
          FROM ro_union u
         WHERE u.seq_md5 IN (SELECT seq_md5 FROM ro_xref
                              WHERE db = 'interpro' AND identifier IN (?, ?))
           AND (u.rieske_intact IS NULL OR u.rieske_intact = 0)
         GROUP BY 1 ORDER BY 1 DESC""", tuple(RIESKE)).fetchall()

    report = {
        "question": ("does InterPro, annotating the same sequences with its "
                     "own profiles, agree that they carry a Rieske domain?"),
        "unit": "one amino acid sequence (seq_md5)",
        "independence": {
            "independent": ("given a sequence, whether InterPro places a "
                            "Rieske domain on it. Those profiles were built "
                            "by other people, from other data, at other "
                            "thresholds"),
            "not_independent": ("which sequences get asked. The candidate "
                                "pool was selected with this project's own "
                                "Rieske and motif models, so InterPro cannot "
                                "object to a sequence we never included"),
            "therefore": ("this measures precision, and says nothing at all "
                          "about recall: the sequences we missed are outside "
                          "it"),
        },
        "denominators": {
            "union_sequences": len(ours),
            "with_any_interpro_entry": len(annotated),
            "without_any_interpro_entry": len(ours) - len(annotated),
            "why_the_rest_are_excluded": (
                "a sequence InterPro has not annotated does not contradict "
                "us, it is unannotated -- most are recent unreviewed entries. "
                "An absence of evidence is not evidence of absence"),
        },
        "rieske_agreement": {
            "interpro_entries_used": RIESKE,
            "our_rieske_intact": len(ours_rieske),
            "interpro_rieske": len(with_rieske),
            "both": len(both),
            "ours_only": len(ours_only),
            "interpro_only": len(theirs_only),
            "agreement_share": (round(len(both) / len(ours_rieske), 4)
                                if ours_rieske else None),
        },
        "interpro_sees_rieske_we_do_not": {
            "sequences": len(theirs_only),
            "all_are_excluded_by_our_gate": all(
                not ours[m]["passes_gate"] for m in theirs_only),
            "by_mononuclear_triad_sites_found": {
                str(t): n for t, n in detail},
            "why_this_is_not_a_contradiction": (
                "InterPro's measure is profile similarity -- is a Rieske "
                "domain present -- while ours is whether the four cluster "
                "ligands are in place. A Rieske fold without its cluster "
                "cannot do the chemistry, so the gate asks for the ligands. "
                "The two statements are different, not opposed"),
            "why_it_is_published_anyway": (
                "these sequences are dropped by the gate, and a drop that "
                "nobody can see is the failure mode this database keeps "
                "finding. 76 of the 139 carry all three mononuclear iron "
                "ligands, which makes them a set a curator should look at: a "
                "substituted cluster ligand with an intact catalytic triad"),
        },
        "family_entries_by_our_group": {
            g: dict(c.most_common()) for g, c in sorted(by_group.items())},
        "family_entries_are_not_errors": (
            "InterPro's separate family names -- cholesterol 7-desaturase, "
            "pheophorbide a oxygenase, vanillate O-demethylase, "
            "3-ketosteroid 9-alpha-monooxygenase -- are not disagreements. "
            "This database collects Rieske oxygenases as a whole rather than "
            "ring-hydroxylating dioxygenases alone, and those names are what "
            "that breadth looks like from the outside"),
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    r = report["rieske_agreement"]
    print(f"[payda] {len(annotated)}/{len(ours)} dizinin InterPro kaydi var")
    print(f"[rieske] bizde {r['our_rieske_intact']}, InterPro'da "
          f"{r['interpro_rieske']}, ikisinde {r['both']}  "
          f"(yalniz bizde {r['ours_only']}, yalniz InterPro'da "
          f"{r['interpro_only']})")
    print(f"[yazildi] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
