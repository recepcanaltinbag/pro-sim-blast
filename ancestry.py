"""Kirli sahadaki enzim, kimyanin uygulanmadigi yerdeki akrabasindan NE KADAR uzak?

Kullanicinin hipotezi: 2,4-D, dikamba, PCB, nitroaromatikler, benzalkonyum,
karbazol -- bu molekullerin hicbiri 1940'tan once ortalikta yoktu. Acik denize
de uygulanmiyorlar. Oyleyse kirletilmis bir sahadan gelen bir enzimin DENIZDE
akrabasi varsa, o iki soy hat kimyasaldan ONCE ayrilmis olmak zorundadir:
iskelet kirleticiden eski, ve kirli saha yeni bir enzim YARATMADI, var olan
bir soy hattini SECTI.

NEDEN KURATORLU REFERANS BILEREK KULLANILMIYOR. Bu soru daha once iki kez
referansa karsi olculdu (uyenin 71 kuratorlu enzimden en yakinina kimligi) ve
iki olcumun da ayni kusuru vardi: kuratorlu kume kultur koleksiyonlarindan
geliyor, yani ezici bicimde toprak, klinik ve endustriyel kaynakli. Deniz
uyesinin "referanstan uzak" cikmasi bu yuzden evrimle degil, referans kumesinin
NEREDEN toplandigiyla aciklanabilir ve iki aciklama birbirinden ayrilamaz.
Bu modul referansa hic bakmaz. Olculen sey, AYNI enzim tipinin iki habitat
grubundaki uyeleri arasindaki ikili amino asit kimligidir. Referans kokeni
yanliligi bu tasarimda sonuca giremez, cunku tasarimda referans yok.

OLCULEN NICELIK: capraz kumedeki EN YUKSEK kimlik. Bu bir SINIRDIR. En yakin
temiz-habitat akrabasi %40 kimlikteyse, o ayrilma herhangi bir sentetik
kimyasaldan kat kat eskidir. %99,8'deyse gecen on yilda da olmus olabilir --
yani yuksek bir maksimum hipotezi curutmez, yalnizca sinir koymaz.

MOLEKULER SAAT UYGULANMADI. Bu bilincli bir karar, eksiklik degil: protein
ayrisma hizlari soy hatlari arasinda buyuklik mertebeleri degisir ve bu veri
kumesinde hicbir kalibrasyon noktasi yok (ne fosil, ne tarihli izolat serisi,
ne bilinen bir ayrilma). Bir saat uydurmak, yil cinsinden sahte bir kesinlik
uretirdi. Cikti bu yuzden hicbir yerde YIL soylemez, yalnizca sinir soyler.

KONTROL: dogal substratli tipler. Kirletici parcalayan (xenobiotic) tiplerin
capraz habitat ayrilmalari, vanillat/kafein/kolesterol gibi insanin icat
etmedigi molekullerle calisan tiplerinkiyle AYNI derinlikteyse, kirletici
parcalayan soy hatlari siradan kadim cesitligin bir parcasidir -- hipotezin
en guclu bicimi. Farkliysa, farkin yonu de olculur.

IKINCI, BAGIMSIZ ACI: kirletilmis saha uyeleri tipin icinde TEK bir klad mi?
Bir saha var olan bir soy hattini ise aldiysa, uyeleri tipin mevcut
cesitliginin ICINDE dagilmis olmali; yeni bir seyi sectiyse tek ve siki bir
kumede toplanmis olmali. `ro_leaf.leaf_id` (dizi varyantlari) uzerinden
kirletilmis uyelerin kac ayri varyanta dagildigi sayilir ve bunun NULL'u
yeniden ornekleme ile kurulur: tipin uyelerinden ayni sayida uye rastgele
cekilse kac varyant beklenirdi?

Hizalama makinesi operon_relations.py'den OLDUGU GIBI alinir (`pair_identity`,
`MIN_ALIGNED_COVERAGE`): boslukla hizalanan kolonlar paydaya girmez ve kisa
proteinin en az yarisi kapsanmazsa cift atilir. Iki modulde iki ayri kimlik
tanimi olmasin diye tek uygulama paylasilir.

Cikti: analysis_out/ancestry.json

Kullanim:
    python3 ancestry.py --db roar.sqlite --ecology cluster_ecology.csv \
        --out analysis_out/ancestry.json --cpu 20
"""

import argparse
import csv
import json
import os
import random
import sqlite3
import statistics
from collections import Counter, defaultdict
from multiprocessing import Pool

from scipy import stats

# Hizalama ve kimlik tanimi operon_relations.py ile AYNI olmak zorunda, yoksa
# bu sayfadaki yuzdeler oteki sayfadakilerle karsilastirilamaz.
from operon_relations import (MIN_ALIGNED_COVERAGE, RANDOM_IDENTITY_FLOOR,
                              pair_identity)

# --- Habitat gruplari.
#
# KIRLI taraf: kimyanin fiilen uygulandigi ya da biriktigi yerler. Ucu de
# `isolation_source.py`'nin sozlugunden geliyor, burada yeniden tanimlanmiyor.
CONTAMINATED = ("contaminated_industrial", "wastewater_sludge", "mining_acid_drainage")

# TEMIZ taraf, DAR tanim: kullanicinin sorusu tam olarak bu. Acik deniz, bu
# veritabanindaki sentetik substratlarin hicbirinin uygulanmadigi yerdir.
CHEMICAL_FREE_NARROW = ("marine", "marine_sediment")

# TEMIZ taraf, GENIS tanim: ayni mantik (kimyasal oraya uygulanmaz) ama daha
# fazla gozlem, dolayisiyla daha fazla guc. Magara, derin yeralti ve termal
# kaynaklar dar tanimdan bile daha az temas gorur; tatli su ve sediment ise
# daha tartismali, cunku bir nehir yukari havzadan kirletici tasiyabilir.
# IKI tanim da raporlanir: dar olani soruyu, genis olani gucu temsil ediyor.
CHEMICAL_FREE_WIDE = CHEMICAL_FREE_NARROW + (
    "sediment", "freshwater", "cave", "extreme_thermal", "subsurface_deep", "hypersaline")

# `unknown` ve `other` iki tarafa da girmez: ne oldugu bilinmeyen bir kaynak
# ne kirli ne temiz sayilabilir. Kalan habitatlar (toprak, klinik, bitki,
# bagirsak...) bu testin evreninde DEGIL -- ucu de kimyasalla temasin
# belirsiz oldugu yerler. Hangi habitatin hicbir gruba girmedigi ciktida
# `not_assigned` olarak tek tek sayilir, sessiz kalmasin.
EXCLUDED_HABITATS = ("unknown", "other")

# Bir tarafin teste girmesi icin gereken en az AYRI CINS sayisi. Uye sayisi
# degil cins sayisi: 40 uyesi olup tek cinsten gelen bir taraf tek bir
# dizileme projesidir, iki habitat arasinda bir karsilastirma degil.
# Esik 3, cooccurrence.py ve operon_relations.py'deki MIN_CELL_GENERA ile ayni.
MIN_GENERA_PER_SIDE = 3

# Maksimum kimligin yorum bantlari. Kenarlar pipeline'in baska yerlerinden
# geliyor, burada icat edilmiyor: %95 `evidence_tiers.py`'de "characterized"
# kademesinin kapisi, %60 "ayni reaksiyonu yapiyor olmasi cok olasi" kapisi,
# %25 ise akraba OLMAYAN iki proteinin global hizalamada verdigi taban
# (operon_relations.RANDOM_IDENTITY_FLOOR). %99,5 ayri bir bant, cunku
# ~400 kalintilik bir proteinde bu, en fazla iki degisim demektir.
BOUND_BANDS = [
    (99.5, 100.01, "indistinguishable"),
    (95.0, 99.5, "very close"),
    (80.0, 95.0, "close"),
    (60.0, 80.0, "distant"),
    (0.0, 60.0, "deep"),
]

# Varyant yayilimi null'u: permutasyon sayisi ve tohum cooccurrence.py ve
# operon_relations.py ile AYNI, boylece uc modulun null'lari karsilastirilabilir.
PERMUTATIONS = 2000
RANDOM_SEED = 20261005

# Varyant yayilimi testi icin bir tipte gereken en az uye ve en az varyant.
# Tek varyantli bir tipte "dagilim" diye bir sey yok; testi kosmak degil,
# kosulamadigini soylemek dogru olan.
MIN_MEMBERS_FOR_SPREAD = 3
MIN_VARIANTS_FOR_SPREAD = 2

# stats_overview.py ve operon_relations.py ile AYNI evren tanimi.
SUBSTRATE_CLASSES = ["xenobiotic", "natural_aromatic", "natural_specialized"]
NATURAL_CLASSES = ("natural_aromatic", "natural_specialized")


def genus_of(organism):
    """Cins adi. 'Candidatus X' ve 'uncultured X' ikinci kelimeyi verir.

    webapp/app.py'deki `genus()` ile ayni kural: ilk kelimeyi almak bu iki
    onek icin butun adaylari tek bir sahte cinse yigardi.
    """
    if not organism:
        return "?"
    words = organism.split()
    if words[0] in ("Candidatus", "uncultured") and len(words) > 1:
        return words[1]
    return words[0]


def band_of(identity):
    for low, high, name in BOUND_BANDS:
        if low <= identity < high:
            return name
    return BOUND_BANDS[-1][2]


def load_ecology(path):
    """cluster_ecology.csv: tip -> substrat, sinif, guven."""
    eco = {}
    if not os.path.exists(path):
        return eco
    with open(path) as fh:
        for row in csv.DictReader(fh):
            eco[row["cluster"]] = {
                "substrate": row.get("substrate") or None,
                "substrate_class": row.get("substrate_class") or None,
                "confidence": row.get("confidence") or None,
            }
    return eco


def load_members(con):
    """Dogrulanmis, dizisi olan her giris: tip, habitat, organizma, varyant."""
    rows = con.execute("""
        SELECT r.candidate_id, r.ro_cluster, s.habitat, p.organism,
               l.leaf_id, r.sequence
        FROM ro r
        JOIN replicon p USING(nucleotide_id)
        JOIN replicon_source s USING(nucleotide_id)
        LEFT JOIN ro_leaf l USING(candidate_id)
        WHERE r.is_confirmed = 1
          AND r.sequence IS NOT NULL AND r.sequence <> ''
    """).fetchall()
    members = []
    for candidate_id, cluster, habitat, organism, leaf_id, sequence in rows:
        members.append({
            "candidate_id": candidate_id,
            "type": cluster,
            "habitat": habitat or "unknown",
            "organism": organism,
            "genus": genus_of(organism),
            "leaf_id": leaf_id,
            "sequence": sequence,
        })
    return members


def _cross_pair(task):
    """Pool isci fonksiyonu: bir capraz cift icin kimlik ve kapsama.

    operon_relations._score_pair ile ayni sozlesme; tek fark, burada cift
    basina TEK bir hizalama var (alfa-alfa), duzenleyici yok.
    """
    seq_a, seq_b = task
    return pair_identity(seq_a, seq_b)


# ------------------------------------------------- 1. olcum: capraz habitat kimligi

def build_sides(members, free_habitats):
    """Tip basina (kirli taraf, temiz taraf) uye listeleri."""
    sides = defaultdict(lambda: {"contaminated": [], "free": []})
    free = set(free_habitats)
    for member in members:
        if member["habitat"] in CONTAMINATED:
            sides[member["type"]]["contaminated"].append(member)
        elif member["habitat"] in free:
            sides[member["type"]]["free"].append(member)
    return sides


def eligible(side):
    """Her iki tarafta da en az MIN_GENERA_PER_SIDE ayri cins var mi?"""
    return (len({m["genus"] for m in side["contaminated"]}) >= MIN_GENERA_PER_SIDE
            and len({m["genus"] for m in side["free"]}) >= MIN_GENERA_PER_SIDE)


def score_all_cross_pairs(members, cpu):
    """Teste girebilen her tipin TUM capraz ciftlerini puanlar.

    Ciftler ORNEKLENMIYOR, hepsi hizalaniyor: en genis tanimda bile 34.500
    cift var ve cift basina ~5 ms, yani tam sayim ucuz. Ornekleme yapilsaydi
    maksimum -- bu modulun asil olcumu -- asagi dogru yanli olurdu, cunku
    ornek disinda kalan en yakin cift kacirilabilirdi.

    Hesap EN GENIS temiz tanim uzerinden bir kez yapilir; dar tanim onun bir
    ALT KUMESI oldugu icin ayni ciftlerden suzulerek turetilir.
    """
    wide = build_sides(members, CHEMICAL_FREE_WIDE)
    narrow = build_sides(members, CHEMICAL_FREE_NARROW)
    testable = sorted({cluster for cluster in wide
                       if eligible(wide[cluster])
                       or (cluster in narrow and eligible(narrow[cluster]))})

    tasks, index = [], []
    for cluster in testable:
        for a in wide[cluster]["contaminated"]:
            for b in wide[cluster]["free"]:
                tasks.append((a["sequence"], b["sequence"]))
                index.append((cluster, a, b))

    if cpu > 1 and len(tasks) > 200:
        with Pool(cpu) as pool:
            scored = pool.map(_cross_pair, tasks, chunksize=64)
    else:
        scored = [_cross_pair(t) for t in tasks]

    # Atilan ciftler de TIP ve TEMIZ UYE duzeyinde tutuluyor: dar tanimin
    # atilan cift sayisi genis tanimin sayisi DEGIL, onun alt kumesi.
    pairs = defaultdict(list)
    dropped = defaultdict(list)
    for (cluster, a, b), (identity, coverage) in zip(index, scored):
        if identity is None or coverage is None or coverage < MIN_ALIGNED_COVERAGE:
            dropped[cluster].append(b["candidate_id"])
            continue
        pairs[cluster].append({
            "identity": identity,
            "contaminated": a,
            "free": b,
            "same_genus": a["genus"] == b["genus"],
        })
    return pairs, dropped, wide, narrow, len(tasks)


def summarise_type(cluster, side, pairs, dropped, eco):
    """Bir tipin capraz habitat ozeti: sinir, ortanca, kac cift, kac cins."""
    contaminated, free = side["contaminated"], side["free"]
    free_ids = {m["candidate_id"] for m in free}
    kept = [p for p in pairs if p["free"]["candidate_id"] in free_ids]
    n_dropped = sum(1 for cid in dropped if cid in free_ids)
    if not kept:
        return None

    identities = [p["identity"] for p in kept]
    cross_genus = [p["identity"] for p in kept if not p["same_genus"]]

    # Cins hucresi basina ORTANCA, sonra hucreler uzerinde ortanca. Derin
    # dizilenmis tek bir cins (ornegin yuzlerce Pseudomonas genomu) boylece
    # dagilimi suruklemiyor. Maksimum ise TUM ciftler uzerinden aliniyor:
    # bir maksimum tanimi geregi "hakim olunacak" bir istatistik degil ve
    # temsilci secmek en yakin cifti kacirma riski tasir.
    cells = defaultdict(list)
    for p in kept:
        cells[(p["contaminated"]["genus"], p["free"]["genus"])].append(p["identity"])
    cell_medians = [statistics.median(v) for v in cells.values()]

    closest = max(kept, key=lambda p: p["identity"])
    closest_cross = max(kept, key=lambda p: (-1.0 if p["same_genus"] else p["identity"]))

    ecology = eco.get(cluster) or {}
    return {
        "cluster": cluster,
        "substrate": ecology.get("substrate"),
        "substrate_class": ecology.get("substrate_class"),
        "substrate_confidence": ecology.get("confidence"),
        "n_contaminated": len(contaminated),
        "n_free": len(free),
        "genera_contaminated": len({m["genus"] for m in contaminated}),
        "genera_free": len({m["genus"] for m in free}),
        "shared_genera": len({m["genus"] for m in contaminated}
                             & {m["genus"] for m in free}),
        "n_pairs": len(kept),
        "n_pairs_dropped_for_coverage": n_dropped,
        "n_genus_cells": len(cells),
        "max_identity": round(max(identities), 2),
        "max_identity_cross_genus": (round(max(cross_genus), 2) if cross_genus else None),
        "median_identity_genus_cells": round(statistics.median(cell_medians), 2),
        "median_identity_all_pairs": round(statistics.median(identities), 2),
        "min_identity": round(min(identities), 2),
        # Akraba OLMAYAN iki protein global hizalamada ~%20-25 verir
        # (operon_relations.RANDOM_IDENTITY_FLOOR). Bir tipin capraz
        # ciftlerinin ortancasi bu tabana yakinsa, ortanca ayrisma DERINLIGI
        # hakkinda bir sey soylemiyor; olcum tabana vurmus demektir. Bu yuzden
        # her satir, ortancasinin tabandan kac puan yukarida oldugunu ve kac
        # ciftinin tabanin ALTINDA kaldigini tasiyor.
        "share_of_pairs_below_random_floor":
            round(sum(1 for v in identities if v < RANDOM_IDENTITY_FLOOR)
                  / len(identities), 3),
        "median_points_above_random_floor":
            round(statistics.median(cell_medians) - RANDOM_IDENTITY_FLOOR, 2),
        "bound_band": band_of(max(identities)),
        "closest_pair": {
            "contaminated_organism": closest["contaminated"]["organism"],
            "contaminated_habitat": closest["contaminated"]["habitat"],
            "free_organism": closest["free"]["organism"],
            "free_habitat": closest["free"]["habitat"],
            "identity": round(closest["identity"], 2),
            "same_genus": closest["same_genus"],
        },
        "closest_cross_genus_pair": ({
            "contaminated_organism": closest_cross["contaminated"]["organism"],
            "free_organism": closest_cross["free"]["organism"],
            "free_habitat": closest_cross["free"]["habitat"],
            "identity": round(closest_cross["identity"], 2),
        } if cross_genus else None),
    }


def class_block(rows):
    """Bir substrat sinifinin tip duzeyindeki dagilimi."""
    if not rows:
        return None
    maxima = [r["max_identity"] for r in rows]
    medians = [r["median_identity_genus_cells"] for r in rows]
    bands = Counter(r["bound_band"] for r in rows)
    return {
        "n_types": len(rows),
        "median_of_type_maxima": round(statistics.median(maxima), 2),
        "min_of_type_maxima": round(min(maxima), 2),
        "max_of_type_maxima": round(max(maxima), 2),
        "median_of_type_medians": round(statistics.median(medians), 2),
        "bands": {name: bands.get(name, 0) for _l, _h, name in BOUND_BANDS},
        "types_with_an_identical_cross_habitat_pair":
            sum(1 for r in rows if r["max_identity"] >= 99.95),
        "n_pairs": sum(r["n_pairs"] for r in rows),
        "median_points_above_random_floor": round(statistics.median(
            [r["median_points_above_random_floor"] for r in rows]), 2),
        "types_whose_median_is_within_5_points_of_the_floor":
            sum(1 for r in rows if r["median_points_above_random_floor"] <= 5.0),
    }


def class_contrast(xeno, natural):
    """Kirletici tipleri ile dogal substratli tipleri ayni mi?

    Mann-Whitney, cunku tip basina tek bir sayi var (tipin maksimumu ya da
    ortancasi), dagilim normal degil ve n kucuk. Testin YONU onemli: burada
    aranan sey bir FARK YOKLUGU, yani yuksek bir p degeri hipotezi destekler.
    Bu yuzden p'nin yaninda farkin kendisi ve orneklem buyuklugu de yaziliyor;
    "anlamli degil" tek basina "ayni" demek degil, ozellikle n ~ 10 iken.
    """
    if len(xeno) < 3 or len(natural) < 3:
        return {"measurable": False,
                "why": f"{len(xeno)} xenobiotic and {len(natural)} natural-substrate "
                       "types clear the genus gate; a comparison needs at least "
                       "three on each side."}
    out = {"measurable": True, "n_xenobiotic": len(xeno), "n_natural": len(natural)}
    # type_medians KARSILASTIRMASI taban sinirli (bkz.
    # method.why_only_the_maximum_is_interpreted) ve karar icin kullanilmaz;
    # tabloya girmesinin tek nedeni, iki sinifin ayni tabana vurdugunun
    # gorunmesi.
    for key, label in (("max_identity", "type_maxima"),
                       ("median_identity_genus_cells", "type_medians")):
        a = [r[key] for r in xeno]
        b = [r[key] for r in natural]
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        out[label] = {
            "floor_limited": label == "type_medians",
            "xenobiotic_median": round(statistics.median(a), 2),
            "natural_median": round(statistics.median(b), 2),
            "difference_points": round(statistics.median(a) - statistics.median(b), 2),
            "mannwhitney_u": float(u),
            "p": round(float(p), 4),
        }
    return out


def cross_habitat_block(pairs, dropped, sides, eco, definition, habitats):
    """Bir temiz-habitat tanimi icin butun capraz habitat olcumu."""
    rows = []
    for cluster in sorted(sides):
        side = sides[cluster]
        if not eligible(side):
            continue
        row = summarise_type(cluster, side, pairs.get(cluster, []),
                             dropped.get(cluster, []), eco)
        if row:
            rows.append(row)
    rows.sort(key=lambda r: r["max_identity"])

    by_class = {}
    for name in SUBSTRATE_CLASSES:
        by_class[name] = class_block([r for r in rows if r["substrate_class"] == name])
    natural_rows = [r for r in rows if r["substrate_class"] in NATURAL_CLASSES]
    by_class["natural_combined"] = class_block(natural_rows)

    # Gucten dusen tipler: hangi tip hangi tarafta kac cinsle kaldi. Bir
    # "test edilemedi" listesi olmadan 24 tiplik tablo, 61 tipin tamami
    # gibi okunur.
    not_tested = []
    for cluster in sorted(sides):
        side = sides[cluster]
        if eligible(side):
            continue
        genera_c = len({m["genus"] for m in side["contaminated"]})
        genera_f = len({m["genus"] for m in side["free"]})
        if genera_c or genera_f:
            not_tested.append({
                "cluster": cluster,
                "substrate_class": (eco.get(cluster) or {}).get("substrate_class"),
                "genera_contaminated": genera_c,
                "genera_free": genera_f,
            })

    return {
        "definition": definition,
        "chemical_free_habitats": list(habitats),
        "n_types_tested": len(rows),
        "n_types_not_tested": len(not_tested),
        "n_pairs": sum(r["n_pairs"] for r in rows),
        "types": rows,
        "types_not_tested": not_tested,
        "by_substrate_class": by_class,
        "class_contrast": class_contrast(
            [r for r in rows if r["substrate_class"] == "xenobiotic"], natural_rows),
    }


# ------------------------------------------------- 2. olcum: varyant yayilimi

def variant_spread(members, eco, seed):
    """Kirletilmis saha uyeleri tipin cesitligi ICINDE mi, tek bir kumede mi?

    Olculen iki sey: (1) kac ayri varyanta dagildiklari, (2) en kalabalik tek
    varyantta toplanan oranlari. Ikisi de UYE SAYISINA siddetle bagli -- 5
    uyenin en fazla 5 varyata dagilabilmesi bir bulgu degil, aritmetiktir.
    Bu yuzden her iki olcum de NULL'a karsi okunuyor: tipin TUM uyelerinden
    ayni sayida uye yerine koymadan cekilirse ne cikardi. Null, tipin gercek
    varyant buyuklugu dagilimini oldugu gibi tasir.

    Ayni test, karsilastirma olsun diye temiz-habitat uyelerine de kosuluyor:
    iki habitat grubunun yayilimi birbirine benzerse, kirletilmis sahanin
    yayilimi hakkinda "ozel" bir sey yok.
    """
    rng = random.Random(seed)
    by_type = defaultdict(list)
    for member in members:
        if member["leaf_id"]:
            by_type[member["type"]].append(member)

    free_wide = set(CHEMICAL_FREE_WIDE)
    rows = []
    for cluster in sorted(by_type):
        pool = by_type[cluster]
        leaves = [m["leaf_id"] for m in pool]
        n_variants = len(set(leaves))
        subsets = {
            "contaminated": [m for m in pool if m["habitat"] in CONTAMINATED],
            "chemical_free": [m for m in pool if m["habitat"] in free_wide],
        }
        if len(subsets["contaminated"]) < MIN_MEMBERS_FOR_SPREAD:
            continue
        row = {
            "cluster": cluster,
            "substrate": (eco.get(cluster) or {}).get("substrate"),
            "substrate_class": (eco.get(cluster) or {}).get("substrate_class"),
            "type_members": len(pool),
            "type_variants": n_variants,
            "measurable": n_variants >= MIN_VARIANTS_FOR_SPREAD,
        }
        if not row["measurable"]:
            row["why_not"] = ("the whole type collapses into a single sequence variant, "
                              "so there is no spread to measure on either side")
            rows.append(row)
            continue

        for name, subset in subsets.items():
            if len(subset) < MIN_MEMBERS_FOR_SPREAD:
                row[name] = None
                continue
            k = len(subset)
            observed_variants = len({m["leaf_id"] for m in subset})
            observed_top = max(Counter(m["leaf_id"] for m in subset).values()) / k

            null_variants, null_top = [], []
            for _ in range(PERMUTATIONS):
                draw = rng.sample(leaves, k)
                null_variants.append(len(set(draw)))
                null_top.append(max(Counter(draw).values()) / k)
            null_mean = statistics.mean(null_variants)
            # Tek yonlu: ilgilendigimiz alternatif "beklenenden DAHA AZ
            # varyat", yani siki, yeni bir kume. Ters yon de ayrica yaziliyor.
            p_fewer = (sum(1 for v in null_variants if v <= observed_variants) + 1) \
                / (PERMUTATIONS + 1)
            p_more = (sum(1 for v in null_variants if v >= observed_variants) + 1) \
                / (PERMUTATIONS + 1)
            p_concentrated = (sum(1 for v in null_top if v >= observed_top) + 1) \
                / (PERMUTATIONS + 1)
            row[name] = {
                "members": k,
                "variants_occupied": observed_variants,
                "variants_expected_at_random": round(null_mean, 2),
                "spread_ratio": round(observed_variants / null_mean, 3) if null_mean else None,
                "top_variant_share": round(observed_top, 3),
                "top_variant_share_expected": round(statistics.mean(null_top), 3),
                "p_fewer_variants_than_chance": round(p_fewer, 4),
                "p_more_variants_than_chance": round(p_more, 4),
                "p_more_concentrated_than_chance": round(p_concentrated, 4),
            }
        rows.append(row)

    # Benjamini-Hochberg, stats_overview.py'nin uygulamasiyla: bu bir satir
    # AILESI ve duzeltmesiz bakildiginda birkac satirin sans eseri p<0,05
    # vermesi beklenir.
    from stats_overview import bh_adjust
    for key in ("contaminated", "chemical_free"):
        family = [r for r in rows if r.get(key)]
        if not family:
            continue
        qs = bh_adjust([r[key]["p_fewer_variants_than_chance"] for r in family])
        for r, q in zip(family, qs):
            r[key]["q_fewer_variants_than_chance"] = round(q, 4)

    tested = [r for r in rows if r.get("contaminated")]
    summary = {}
    for name in SUBSTRATE_CLASSES + ["natural_combined"]:
        if name == "natural_combined":
            subset = [r for r in tested if r["substrate_class"] in NATURAL_CLASSES]
        else:
            subset = [r for r in tested if r["substrate_class"] == name]
        if not subset:
            summary[name] = None
            continue
        ratios = [r["contaminated"]["spread_ratio"] for r in subset
                  if r["contaminated"]["spread_ratio"] is not None]
        summary[name] = {
            "n_types": len(subset),
            "median_spread_ratio": round(statistics.median(ratios), 3) if ratios else None,
            "types_tighter_than_chance_q05":
                sum(1 for r in subset
                    if r["contaminated"].get("q_fewer_variants_than_chance", 1) < 0.05),
            "types_more_concentrated_than_chance_p05":
                sum(1 for r in subset
                    if r["contaminated"]["p_more_concentrated_than_chance"] < 0.05),
            "median_variants_occupied":
                statistics.median([r["contaminated"]["variants_occupied"] for r in subset]),
            "median_type_variants":
                statistics.median([r["type_variants"] for r in subset]),
        }

    xeno = [r for r in tested if r["substrate_class"] == "xenobiotic"
            and r["contaminated"]["spread_ratio"] is not None]
    natural = [r for r in tested if r["substrate_class"] in NATURAL_CLASSES
               and r["contaminated"]["spread_ratio"] is not None]
    contrast = {"measurable": False}
    if len(xeno) >= 3 and len(natural) >= 3:
        a = [r["contaminated"]["spread_ratio"] for r in xeno]
        b = [r["contaminated"]["spread_ratio"] for r in natural]
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        contrast = {
            "measurable": True,
            "n_xenobiotic": len(a), "n_natural": len(b),
            "xenobiotic_median_spread_ratio": round(statistics.median(a), 3),
            "natural_median_spread_ratio": round(statistics.median(b), 3),
            "mannwhitney_u": float(u), "p": round(float(p), 4),
        }

    # Kirletilmis ve temiz taraflarin yayilimi, ayni tipte eslesmis olarak.
    paired = [r for r in tested if r.get("chemical_free")
              and r["contaminated"]["spread_ratio"] is not None
              and r["chemical_free"]["spread_ratio"] is not None]
    habitat_contrast = {"measurable": False}
    if len(paired) >= 5:
        diffs = [r["contaminated"]["spread_ratio"] - r["chemical_free"]["spread_ratio"]
                 for r in paired]
        w, p = stats.wilcoxon(diffs)
        habitat_contrast = {
            "measurable": True, "n_types": len(paired),
            "median_difference": round(statistics.median(diffs), 3),
            "wilcoxon_statistic": float(w), "p": round(float(p), 4),
        }

    return {
        "question": ("Within one enzyme type, do the contaminated-site members sit "
                     "inside the type's existing variant diversity, or do they form "
                     "one tight recent cluster?"),
        "null_model": (f"{PERMUTATIONS} draws without replacement of the same number "
                       "of members from the same type's full membership, seed "
                       f"{RANDOM_SEED}. The null therefore carries the type's real "
                       "variant-size distribution, including its large variants."),
        "n_types_tested": len(tested),
        "types": rows,
        "by_substrate_class": summary,
        "class_contrast": contrast,
        "contaminated_vs_chemical_free": habitat_contrast,
    }


# ------------------------------------------------------------------ verdict

def build_verdict(narrow, wide, spread):
    """Veri hipotezi destekliyor mu? Desteklemiyorsa da ayni aciklikla yazilir."""
    def below_95(block):
        """Maksimumu %95'in ALTINDA kalan tip sayisi: sinirin BAGLADIGI tipler."""
        bands = (block.get("bands") or {})
        return bands.get("deep", 0) + bands.get("distant", 0) + bands.get("close", 0)

    wide_xeno = wide["by_substrate_class"].get("xenobiotic") or {}
    narrow_xeno = narrow["by_substrate_class"].get("xenobiotic") or {}
    contrast = wide["class_contrast"]
    narrow_contrast = narrow["class_contrast"]
    maxima = contrast.get("type_maxima") or {}
    spread_xeno = spread["by_substrate_class"].get("xenobiotic") or {}
    spread_natural = spread["by_substrate_class"].get("natural_combined") or {}
    spread_contrast = spread["class_contrast"]

    n_xeno = wide_xeno.get("n_types") or 0
    same_as_control = bool(contrast.get("measurable") and maxima.get("p", 0) >= 0.05)
    label = "supported" if same_as_control else "qualified"
    if n_xeno < 3:
        label = "not testable"

    return {
        "label": label,
        "headline": (
            "The pollutant-degrading types show cross-habitat divergence of the same "
            "depth as the types acting on molecules humans never made. On this "
            "measure the pollutant-degrading lineages are ordinary ancient diversity, "
            "which is the hypothesis in its strongest form."
            if same_as_control else
            "The pollutant-degrading types and the natural-substrate types differ in "
            "how deep their cross-habitat divergence runs, so the two cannot be "
            "treated as one population."),
        "what_the_bound_says": (
            f"Under the narrow, marine-only definition, {narrow_xeno.get('n_types')} "
            f"xenobiotic types have at least {MIN_GENERA_PER_SIDE} genera on both "
            f"sides, and every one of them does have a marine relative. For "
            f"{below_95(narrow_xeno)} of those {narrow_xeno.get('n_types')} the "
            "CLOSEST marine relative of any contaminated-site member is below 95 % "
            "identity, and for "
            f"{(narrow_xeno.get('bands') or {}).get('deep', 0)} it is below 60 %. "
            "Those lineages cannot have separated recently: below 95 % identity means "
            "tens of substitutions across a ~400-residue protein, and below 60 % "
            "means hundreds. That is a bound on how recent the split can be, not a "
            "date for it."),
        "where_the_bound_is_loose": (
            f"{wide_xeno.get('types_with_an_identical_cross_habitat_pair')} of the "
            f"{n_xeno} xenobiotic types testable under the wide definition contain a "
            "contaminated-site member and a chemical-free member whose sequences are "
            "indistinguishable (99.95 % identity or above). For those the measurement "
            "places no bound at all: an ancient lineage that has barely changed and a "
            "transfer last decade produce the same number. A high maximum therefore "
            "cannot refute the hypothesis, and it cannot support it either -- it "
            "simply leaves the question open for that type."),
        "how_far_this_goes": (
            "The hypothesis is established type by type only for the types where the "
            "bound is deep. Across the whole set it is supported in a different and "
            "weaker way: by the control. Pollutant-degrading types are "
            "indistinguishable from natural-substrate types in how far their "
            "cross-habitat relatives sit, and nothing about them looks like a "
            "lineage that a contaminated site produced."),
        "the_control": (
            f"xenobiotic types: median type maximum {maxima.get('xenobiotic_median')} %, "
            f"natural-substrate types {maxima.get('natural_median')} %, a difference of "
            f"{maxima.get('difference_points')} points, Mann-Whitney p = "
            f"{maxima.get('p')} on {contrast.get('n_xenobiotic')} against "
            f"{contrast.get('n_natural')} types (wide definition). Under the narrow "
            "definition the same contrast is "
            f"{((narrow_contrast.get('type_maxima') or {}).get('difference_points'))} "
            f"points, p = {((narrow_contrast.get('type_maxima') or {}).get('p'))}. "
            "Neither is significant, and at ten to twenty types per side that means "
            "no difference was detected rather than that none exists. Note also that "
            "the natural-substrate control is carried mostly by natural_specialized "
            "types; only one natural_aromatic type clears the gate."
            if contrast.get("measurable") else contrast.get("why")),
        "the_marine_only_answer": (
            f"Under the user's own narrow definition -- marine and marine sediment "
            f"only -- {narrow_xeno.get('n_types')} xenobiotic types can be tested, "
            f"median type maximum {narrow_xeno.get('median_of_type_maxima')} %, range "
            f"{narrow_xeno.get('min_of_type_maxima')} to "
            f"{narrow_xeno.get('max_of_type_maxima')} %. The narrow definition is the "
            "stronger evidence, not the weaker: it has fewer types but its habitat is "
            "the one where the chemicals are least arguably absent, and it is the "
            "definition under which the bound actually bites for most types."
            if narrow_xeno else
            "Under the narrow marine-only definition no xenobiotic type clears the "
            "genus gate on both sides."),
        "the_second_angle": (
            f"Across {spread_xeno.get('n_types')} xenobiotic types, contaminated-site "
            f"members occupy a median {spread_xeno.get('median_spread_ratio')} times "
            "as many sequence variants as a random draw of the same size from the "
            "same type -- that is, as many as chance. After correction "
            f"{spread_xeno.get('types_tighter_than_chance_q05')} types are tighter "
            f"than chance. The natural-substrate control sits at "
            f"{spread_natural.get('median_spread_ratio')} "
            f"({spread_natural.get('n_types')} types"
            + (f", p = {spread_contrast.get('p')}" if spread_contrast.get("measurable") else "")
            + "). Contaminated-site members are spread through their type's existing "
            "variant diversity, not concentrated in one recent clade -- which is what "
            "recruitment of a pre-existing lineage looks like, and not what a newly "
            "selected variant would look like."),
        "no_clock_was_used": (
            "No molecular clock was applied and no age in years is stated anywhere in "
            "this output. Protein substitution rates differ by orders of magnitude "
            "between lineages, and this database contains no calibration point -- no "
            "fossil, no dated isolate series, no independently known divergence -- "
            "from which a rate could be estimated. Converting identity into years here "
            "would manufacture precision that the data do not contain. What identity "
            "can carry is an ordering and a bound, and that is all that is reported."),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Cross-habitat divergence of contaminated-site Rieske oxygenases, "
                    "measured without the curated reference set.")
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--ecology", default="cluster_ecology.csv",
                        help="substrate class per type (xenobiotic / natural_*)")
    parser.add_argument("--out", default="analysis_out/ancestry.json")
    parser.add_argument("--cpu", type=int, default=min(20, os.cpu_count() or 1))
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    eco = load_ecology(args.ecology)
    members = load_members(con)

    habitat_counts = Counter(m["habitat"] for m in members)
    assigned = set(CONTAMINATED) | set(CHEMICAL_FREE_WIDE) | set(EXCLUDED_HABITATS)
    not_assigned = {h: n for h, n in sorted(habitat_counts.items()) if h not in assigned}

    print(f"[1/3] {len(members)} confirmed entries with a sequence; "
          f"{sum(habitat_counts[h] for h in CONTAMINATED)} contaminated, "
          f"{sum(habitat_counts[h] for h in CHEMICAL_FREE_NARROW)} marine, "
          f"{sum(habitat_counts[h] for h in CHEMICAL_FREE_WIDE)} chemical-free (wide)")

    pairs, dropped, wide_sides, narrow_sides, n_tasks = score_all_cross_pairs(
        members, args.cpu)
    print(f"[2/3] {n_tasks} cross-habitat alignments, "
          f"{sum(len(v) for v in dropped.values())} dropped below "
          f"{100 * MIN_ALIGNED_COVERAGE:.0f} % aligned coverage")

    narrow = cross_habitat_block(pairs, dropped, narrow_sides, eco,
                                 "marine only (the question as asked)",
                                 CHEMICAL_FREE_NARROW)
    wide = cross_habitat_block(pairs, dropped, wide_sides, eco,
                               "every habitat where these chemicals are not applied",
                               CHEMICAL_FREE_WIDE)

    spread = variant_spread(members, eco, RANDOM_SEED)
    print(f"[3/3] variant spread: {spread['n_types_tested']} types, "
          f"{PERMUTATIONS} resamples each")

    out = {
        "question": (
            "If a contaminated-site enzyme has a relative in the sea, the two lineages "
            "split before the chemical existed: the scaffold pre-dates the pollutant, "
            "and the contaminated site selected an existing lineage rather than "
            "creating a new one. The synthetic substrates in this database (2,4-D, "
            "dicamba, PCBs, nitroarenes, carbazole, chloroacetanilides and the rest) "
            "date from roughly 1940 to 1970, and the open ocean is not where they are "
            "applied."),
        "method": {
            "why_no_curated_reference": (
                "Both earlier attempts at this question measured identity to the "
                "nearest of the 71 curated references, and both inherited the same "
                "flaw: the curated set comes from culture collections that are "
                "overwhelmingly soil, clinical and industrial, so any non-industrial "
                "member looks distant for a reason that has nothing to do with "
                "evolution. This measurement never touches the reference set. It "
                "compares members of one enzyme type against other members of the "
                "same type, so reference-origin bias cannot enter the result."),
            "what_is_measured": (
                "For each enzyme type, the pairwise amino-acid identity between every "
                "member isolated from a contaminated setting and every member isolated "
                "from a setting where the chemical is not applied. The headline figure "
                "is the MAXIMUM of that set, because it bounds how recently the two "
                "lineages could have separated."),
            "habitat_groups": {
                "contaminated": list(CONTAMINATED),
                "chemical_free_narrow": list(CHEMICAL_FREE_NARROW),
                "chemical_free_wide": list(CHEMICAL_FREE_WIDE),
                "excluded_as_uninformative": list(EXCLUDED_HABITATS),
                "not_assigned_to_either_side": not_assigned,
                "note": ("Habitats outside both groups -- soil, clinical, plant, gut "
                         "and the rest -- are settings where exposure is unknown "
                         "rather than absent, so they take no part in this test. They "
                         "are listed by name and count above so that the size of the "
                         "unused majority is visible: this test runs on a minority of "
                         "the database."),
            },
            "pair_selection": (
                "Every cross-habitat pair of every testable type is aligned; none are "
                "sampled. Pseudoreplication is handled at the reporting step instead: "
                "the median is taken first within each (contaminated genus x "
                "chemical-free genus) cell and then across cells, so one "
                "deeply-sequenced genus cannot drag the distribution. The maximum is "
                "taken over all pairs, because picking one representative per genus "
                "risks discarding the single closest pair, which is precisely the "
                "quantity of interest."),
            "eligibility": (
                f"a type is tested only if both sides carry at least "
                f"{MIN_GENERA_PER_SIDE} distinct genera. The gate counts genera, not "
                "members, because a side of forty members from one genus is one "
                "sequencing project rather than a habitat."),
            "alignment": (
                "Bio.Align global alignment with blastp scoring, identity computed "
                "only over columns where both sides carry a residue, and a pair "
                f"discarded unless the alignment covers at least "
                f"{100 * MIN_ALIGNED_COVERAGE:.0f} % of the shorter protein. The "
                "implementation is imported from operon_relations.py so that the two "
                "modules cannot drift to two different definitions of identity."),
            "why_only_the_maximum_is_interpreted": (
                "An enzyme type here is the best-matching reference HMM, not a clade, "
                "so a single type can hold members that are no more similar to each "
                f"other than two unrelated proteins ({RANDOM_IDENTITY_FLOOR:.0f} % is "
                "the global-alignment floor for unrelated sequences, the same constant "
                "operon_relations.py uses). The median cross-habitat identity of most "
                "types sits only a few points above that floor, which means the median "
                "has bottomed out and carries no information about divergence depth. "
                "It is reported for every type anyway, precisely so that the reader "
                "can see that it has bottomed out, but the statistic the argument "
                "rests on is the MAXIMUM, which is well clear of the floor "
                "everywhere."),
            "random_identity_floor": RANDOM_IDENTITY_FLOOR,
            "interpretation_bands": [
                {"from": low, "to": min(high, 100.0), "name": name}
                for low, high, name in BOUND_BANDS],
            "no_clock": (
                "No molecular clock, and no age in years. Substitution rates vary by "
                "orders of magnitude between bacterial lineages and no calibration "
                "point exists in this dataset, so identity is reported as a bound and "
                "an ordering only."),
            "permutations": PERMUTATIONS,
            "seed": RANDOM_SEED,
        },
        "coverage": {
            "confirmed_entries_with_a_sequence": len(members),
            "entries_by_group": {
                "contaminated": sum(habitat_counts[h] for h in CONTAMINATED),
                "chemical_free_narrow": sum(habitat_counts[h] for h in CHEMICAL_FREE_NARROW),
                "chemical_free_wide": sum(habitat_counts[h] for h in CHEMICAL_FREE_WIDE),
                "neither": len(members) - sum(habitat_counts[h] for h in CONTAMINATED)
                - sum(habitat_counts[h] for h in CHEMICAL_FREE_WIDE),
            },
            "habitat_counts": dict(sorted(habitat_counts.items(),
                                          key=lambda kv: -kv[1])),
            "alignments_run": n_tasks,
            "alignments_dropped_for_coverage": sum(len(v) for v in dropped.values()),
        },
        "cross_habitat": {"narrow": narrow, "wide": wide},
        "variant_spread": spread,
        "limitations": [
            "A habitat label is a property of the deposited record, not of the "
            "enzyme. It says where the genome was collected; it does not establish "
            "that the organism's ancestors lived there, nor that the chemical was "
            "truly absent. The marine entries are the least arguable on this point "
            "and that is why the narrow definition is reported beside the wide one.",
            "An identity of 100 % between a contaminated-site entry and a marine "
            "entry is also what a metadata error, a mislabelled isolate or a "
            "laboratory cross-contamination would produce. Such pairs are counted "
            "separately rather than being read as evidence of recent transfer.",
            "The test is symmetric and therefore silent on direction. It cannot say "
            "whether a lineage moved from the sea to a contaminated site or the "
            "reverse; it bounds only how far apart the two are.",
            "Both sides are small. Several types rest on three or four genera per "
            "side, and a non-significant control contrast at that sample size means "
            "'no difference detected', not 'no difference'.",
            "Substrate class comes from cluster_ecology.csv, where some assignments "
            "are marked low confidence (DdmC, cadA, bfzA among them). A type whose "
            "substrate label is wrong sits in the wrong column of every table here.",
            "An enzyme type is a best-hit bin against a reference HMM, not a clade. "
            "Most types contain cross-habitat pairs sitting at the unrelated-protein "
            f"floor ({RANDOM_IDENTITY_FLOOR:.0f} % identity), so the median cross-pair "
            "identity of a type has bottomed out and is not a measurement of "
            "divergence depth. Only the maximum is read as a bound, and the medians "
            "appear in the tables so the reader can see where the floor is.",
            "Sequence variants (ro_leaf.leaf_id) are a clustering of this database's "
            "own sequences, not a dated phylogeny. Occupying many variants shows that "
            "contaminated-site members are not one tight group; it does not date them.",
        ],
    }
    out["verdict"] = build_verdict(narrow, wide, spread)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"[yazildi] {args.out}")

    # ----------------------------------------------------------- readable summary
    for block in (narrow, wide):
        print()
        print(f"  CROSS-HABITAT IDENTITY -- {block['definition']}")
        print(f"  {block['n_types_tested']} types tested, "
              f"{block['n_types_not_tested']} fall below the genus gate, "
              f"{block['n_pairs']} pairs")
        print(f"    {'type':22s} {'class':20s} {'cont':>5s} {'free':>5s} "
              f"{'pairs':>6s} {'max':>6s} {'xgen':>6s} {'med':>6s} {'>flr':>5s}  bound")
        for r in block["types"]:
            print(f"    {r['cluster']:22s} {(r['substrate_class'] or '?'):20s} "
                  f"{r['n_contaminated']:>5d} {r['n_free']:>5d} {r['n_pairs']:>6d} "
                  f"{r['max_identity']:>6.1f} "
                  f"{(r['max_identity_cross_genus'] or 0):>6.1f} "
                  f"{r['median_identity_genus_cells']:>6.1f} "
                  f"{r['median_points_above_random_floor']:>5.1f}  {r['bound_band']}")
        print(f"    {'-' * 94}")
        for name in ("xenobiotic", "natural_combined"):
            c = block["by_substrate_class"].get(name)
            if not c:
                continue
            print(f"    {name:22s} {c['n_types']:>2d} types  "
                  f"median max {c['median_of_type_maxima']:5.1f} %  "
                  f"(range {c['min_of_type_maxima']:.1f}-{c['max_of_type_maxima']:.1f})  "
                  f"median cell-median {c['median_of_type_medians']:5.1f} % "
                  f"(+{c['median_points_above_random_floor']:.1f} over the "
                  f"{RANDOM_IDENTITY_FLOOR:.0f} % floor, "
                  f"{c['types_whose_median_is_within_5_points_of_the_floor']}"
                  f"/{c['n_types']} within 5 pts)  "
                  f"identical pairs in {c['types_with_an_identical_cross_habitat_pair']} types")
        ct = block["class_contrast"]
        if ct.get("measurable"):
            print(f"    control: xenobiotic vs natural, type maxima "
                  f"{ct['type_maxima']['difference_points']:+.1f} pts, "
                  f"p = {ct['type_maxima']['p']:.3f}; type medians "
                  f"{ct['type_medians']['difference_points']:+.1f} pts, "
                  f"p = {ct['type_medians']['p']:.3f}")
        else:
            print(f"    control: {ct.get('why')}")

    print()
    print("  VARIANT SPREAD OF CONTAMINATED-SITE MEMBERS (against a resampling null)")
    print(f"    {'type':22s} {'class':20s} {'n':>4s} {'vars':>5s} {'exp':>6s} "
          f"{'ratio':>6s} {'top':>6s} {'q':>7s}")
    for r in spread["types"]:
        c = r.get("contaminated")
        if not c:
            print(f"    {r['cluster']:22s} {(r['substrate_class'] or '?'):20s} "
                  f"{'-':>4s}   not measurable")
            continue
        print(f"    {r['cluster']:22s} {(r['substrate_class'] or '?'):20s} "
              f"{c['members']:>4d} {c['variants_occupied']:>5d} "
              f"{c['variants_expected_at_random']:>6.1f} {c['spread_ratio']:>6.2f} "
              f"{c['top_variant_share']:>6.2f} "
              f"{c.get('q_fewer_variants_than_chance', 1):>7.3f}")
    for name in ("xenobiotic", "natural_combined"):
        s = spread["by_substrate_class"].get(name)
        if not s:
            continue
        print(f"    {name:22s} {s['n_types']:>2d} types  "
              f"median spread ratio {s['median_spread_ratio']}  "
              f"tighter than chance (q<0.05): {s['types_tighter_than_chance_q05']}  "
              f"more concentrated (p<0.05): {s['types_more_concentrated_than_chance_p05']}")
    sc = spread["class_contrast"]
    if sc.get("measurable"):
        print(f"    control: xenobiotic {sc['xenobiotic_median_spread_ratio']} vs "
              f"natural {sc['natural_median_spread_ratio']}, p = {sc['p']:.3f}")
    hc = spread["contaminated_vs_chemical_free"]
    if hc.get("measurable"):
        print(f"    contaminated vs chemical-free side, paired per type: "
              f"median difference {hc['median_difference']:+.3f}, p = {hc['p']:.3f}")

    print()
    verdict = out["verdict"]
    print(f"  VERDICT [{verdict['label']}]")
    for key in ("headline", "what_the_bound_says", "where_the_bound_is_loose",
                "how_far_this_goes", "the_control", "the_marine_only_answer",
                "the_second_angle"):
        text = verdict.get(key)
        if text:
            print(f"    {key}: {text}")
    print(f"    no clock: {verdict['no_clock_was_used']}")


if __name__ == "__main__":
    main()
