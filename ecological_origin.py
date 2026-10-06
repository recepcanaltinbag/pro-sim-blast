"""
Bir enzim tipinin YAKIN ve UZAK akrabalari ayni ekolojide mi bulunuyor?

NEDEN BU SCRIPT VAR. `isolation_source.py` habitat dagilimini tip basina
veriyor, ama tek bir dagilim olarak: bir tipin butun uyeleri ayni kovaya
giriyor. Oysa `ro_evidence.tier` o uyelerin kuratorlu enzime ne kadar YAKIN
oldugunu biliyor ve bunlar iki ayri sorudur. "VanA toprakta bulunur" ifadesi
VanA'nin kendisi icin mi dogru, yoksa ona %30 benzeyen tanimsiz akrabalari
icin mi? Enzimin NEREDE evrildigi sorusu tam olarak bu ayrimi gerektirir:
uzak akrabalar farkli bir ekolojide duruyorsa, tipin bugunku habitat profili
tipin kokeni degil, son yayilmasinin izidir.

NE OLCULDU.
  * Tip x kanit kademesi x habitat: her kademe icin ayri habitat bilesimi,
    hem giris hem farkli TUR sayisiyla.
  * Yakin (characterized + close_homolog) ile uzak (family_member + distant +
    novel) taraflarin habitat dagilimi arasindaki Jensen-Shannon uzakligi, bir
    permutasyon null'i, BH duzeltmesi ve sayilarin izin verdigi tiplerde
    chi-kare.
  * ORTAK CINS KONTROLU: ayni ayrisma, yalnizca iki tarafta da bulunan
    cinsler uzerinde. Habitat organizmanin ozelligidir; bu kontrol olmadan
    "uzak akrabalar baska ekolojide" ifadesi "uzak akrabalar baska cinste"
    ifadesinden ayirt edilemez.
  * Varyant ("haplotip") duzeyi: her yapragin boyutu, baskin habitati ve
    habitat bilesimi; her tipin EN BUYUK varyantinin nerede goruldugu.
    Haplotip agindaki "ana varyant nerede?" sorusunun veri karsiligi budur.
  * Atik su / aktif camur (habitat = wastewater_sludge) ayri bir blok: tur
    bazinda zenginlesme ve bu uyelerin kuratorlu referansa uzakligi.

SONUC (olculen, yorum degil).
  * Havuzlanmis olarak yakin ve uzak uyelerin habitat dagilimi gercekten
    farkli ama fark KUCUK: JS uzakligi 0,019 bit (permutasyon p = 0,001),
    chi-kare 253,9 / 25 sd, p = 9,1e-40, Cramer V = 0,198. Yon: yakin uyelerde
    insan klinigi (%15,4'e karsi %12,2), bagirsak (%3,9/%2,2), kirli saha
    (%5,7/%4,3) ve atik su (%4,8/%3,5) payi yuksek; uzak uyelerde bitki
    dokusu (%4,5/%7,3), deniz (%6,4/%8,8) ve termal ortamlar yuksek.
    Havuzlanmis fark ortak cins kontrolunu GECIYOR (127 ortak cins, JS 0,0123,
    p = 0,002), yani tamami taksonomik devir teslim degil.
  * Tip duzeyinde ayrisma 61 tipin yalnizca 11'inde olculebildi (iki tarafta da
    >=10 siniflanmis uye); 9'u BH sonrasi q <= 0,05. En buyukleri BmoA
    (0,486; deniz %33,3 -> %9,6), Stc2 (0,417; deniz %43,3 -> %8,9), TPDO
    (0,296), CntA (0,291; insan klinigi %33,9 -> %14,9), XylX (0,270), PDO
    (0,265), CbdC (0,258), NagGH (0,218), KshA15 (0,092). VanA (p = 0,680) ve
    DitA (p = 0,161) ayrismiyor.
  * AMA ORTAK CINS KONTROLUNDEN BU TIPLERIN HICBIRI q <= 0,05 ILE GECMIYOR.
    En yakin olanlar TPDO (JS 0,447, p = 0,010, q = 0,090) ve CbdC (0,413,
    p = 0,028, q = 0,126). BmoA ve DitA'da kontrol HIC kurulamiyor (4 ve 1
    ortak cins): BmoA'nin yakin uyeleri deniz cinsleri (Marinomonas,
    Salinicola, Vibrio), uzak uyeleri toprak ve bitki cinsleridir
    (Novosphingobium, Streptomyces, Sphingobium). Yani tip duzeyindeki habitat
    kaymasi bu veriyle enzime DEGIL, tasiyici taksona atfedilmek zorundadir.
    Negatif bir sonuctur ve boyle yazilmasi gerekir.
  * Atik su kucuk ve zayif bir sinyal: 184 giris, 58 tur, 38 tip; tum
    tur-habitat ciftlerinin %4,0'u. En guvenilir zenginlesme HcaE 3,49x
    (9 tur), sonra XylX 1,61x (11 tur); NBDO 4,59x ve BPDO 3,37x ikisi de
    yalnizca 2 ture dayaniyor ve isaretlendi. En cok atik su girisi olan TsaM
    (18 giris, 13 tur) zenginlesmiyor (1,23x). Buna karsilik atik su uyeleri
    kuratorlu referansa ORTALAMANIN UZERINDE yakin: yakin kademe payi %38,0'e
    karsi arka planda %26,9, karakterize kademesi 2,21x, kimlik ortancasi
    %43,2'ye karsi %37,4. Yani atik sudan gelenler agirlikla bilinen
    enzimlerin kendileri, yeni akrabalari degil.
  * 439 profillenen varyantin HICBIRINDE baskin habitat atik su degil; en
    yuksek atik su payi %33,3 (uc varyant, her birinde 1-3 giris). Atik su
    uyeleri ayri bir haplotip olusturmuyor, varyantlarin icine dagilmis.
  * Ana varyantlar (her tipin en buyuk yapragi, 61 tane) agirlikla toprak
    (19 tip), insan klinigi (8) ve kirli saha (8) kokenli; 8 tipte en buyuk
    varyagin siniflanmis uyesi "baskin habitat" demeye yetmiyor.
  * Kapsama durust verilir: dogrulanmis 11.422 girisin 6.460'inda (%56,6) ve
    1.956 turun 1.094'unde siniflanmis bir habitat var. Kademeler arasinda
    kapsama farki da var (novel %41,0, close_homolog %60,8), yani en az bilinen
    uyelerin ekolojisi de en az biliniyor.

SINIRLAR -- bunlari bastan yazmak gerekiyor.
  * Habitat, izolasyon kaynagi serbest metninden sozlukle atandi
    (`isolation_source.py`); bu bir olcum degil eslemedir ve kapsama %56,6.
    Bir tipin "habitati" aslinda habitati BILINEN azinliginin habitatidir.
  * Giris sayilari bagimsiz gozlem DEGILDIR (ayni sus tekrar tekrar
    dizilenmis). Her ekolojik ifade TUR uzerinden kurulur; giris sayisi
    yalnizca yaninda, denetim icin verilir.
  * Kademeler bir evrim agaci degil, kuratorlu referansa KIMLIK merdivenidir.
    Uzak bir akraba daha ESKI olmak zorunda degil. Bu yuzden "enzim burada
    evrildi" sonucu bu veriden CIKARILAMAZ.
  * Ortak cins kontrolu gectikten sonra bile orneklem yanliligi duruyor:
    hangi genomun dizilendigi, nerede bulundugundan bagimsiz bir karar degil.
    Karakterize enzimler kirlilik calisan laboratuvarlardan gelir, uzak
    akrabalarin cogu klinik ve tarim dizileme programlarindan.
  * Chi-kare GIRIS sayilari uzerinde kosuluyor (tur duzeyinde bir tablo da
    bagimsiz denemeler tablosu olmazdi, cunku bir tur iki tarafta birden
    gorunebilir). Bu yuzden p degeri gercekte oldugundan guclu gorunur ve
    guvenilecek sayi permutasyon testidir.

Kullanim:
    python3 ecological_origin.py --db roar.sqlite --out-dir analysis_out
"""

import argparse
import json
import os
import random
import sqlite3
import zlib
from collections import Counter, defaultdict

from scipy.stats import chi2_contingency

# Habitat sozlugu TEK yerde tanimli olmali. Kopyalamak, projenin en cok
# zarar gormus hata sinifi olan "sessizce ayrisan iki liste" durumunu
# uretirdi; bu yuzden vokabuler uretildigi modulden okunuyor.
from isolation_source import CLASSIFIED, HABITAT_ORDER, species_of

TIERS = ("characterized", "close_homolog", "family_member", "distant", "novel")

# "Yakin" = substrat etiketinin bilgilendirici sayildigi iki kademe, yani
# referansa >=%60 kimlik. `stratified_stats.py` ile birebir ayni tanim;
# sitenin baska yerlerinde de bu esik kullaniliyor.
CLOSE_TIERS = ("characterized", "close_homolog")
DISTANT_TIERS = ("family_member", "distant", "novel")

WASTEWATER_HABITAT = "wastewater_sludge"

# --- Esikler: hepsi adli, hepsi gerekceli -----------------------------------

# Bir tarafta bundan az SINIFLANMIS giris varsa ayrisma olculmez. JS uzakligi
# kucuk orneklemde yukari saplanir (gorulmeyen her habitat payi sifir yapar),
# yani sayi ekolojiden cok aritmetikten gelir. Permutasyon null'i bu sapmayi
# hesaba katiyor ama 10 girisin altinda testin gucu de kalmiyor.
MIN_SIDE_ENTRIES_DIVERGENCE = 10

# Chi-kare daha yuksek bir bar ister: beklenen hucre degerleri ~5 olmali.
# 20 girislik bir tarafta bu, ancak iki-uc kolonda tutar; altinda test
# raporlanmaz (JS uzakligi ve permutasyon p'si yine verilir).
MIN_SIDE_ENTRIES_CHI2 = 20

# Chi-kare kolonlari: bir tipte toplam (yakin+uzak) bundan az giris tasiyan
# habitatlar TEK bir kolonda toplanir. Seyrek kolonlar testi sisirir; atmak
# bilgi kaybi oldugu icin atilmaz, toplanir.
MIN_HABITAT_COLUMN_ENTRIES = 5

# JS uzakliginin null dagilimi: yakin/uzak etiketi tip ICINDE karistirilir.
# Girisler bagimsiz degil, ama etiket karistirmasinda tur kimligi girisle
# birlikte tasinir, yani null veri yapisinin kendi fazlaligini da tasir.
PERMUTATIONS = 999

# Zenginlesme icin en az bu kadar tur-habitat cifti. `webapp/atlas.py`
# icindeki `habitat_enrichment(min_pairs=10)` ile ayni bar; daha azinda
# "kat" degeri tek bir ture dayanir.
MIN_SPECIES_PAIRS = 10

# Yukaridaki bar tipin TOPLAM buyuklugune bakiyor, PAYA bakmiyor: 11 cifti
# olan bir tip atik suda 2 tur gorulmusse %18 pay ve 4,6x kat uretir. Satir
# atilmaz (bilgi kaybi olurdu) ama bu kadar az ture dayanan kat SIRALANACAK
# bir oran degil, bir sayimdir; `fold_rests_on_few_species` ile isaretlenir.
MIN_SPECIES_IN_HABITAT_FOR_FOLD = 3

# Varyant profili: bundan kucuk yapraklarda "baskin habitat" bir dagilim
# degil, tek bir susun adidir. Kucuk yapraklar sayilir ama profilleri
# yazilmaz; her tipin EN BUYUK yapragi boyutu ne olursa olsun yazilir.
MIN_VARIANT_SIZE_FOR_PROFILE = 5

# Bir yapragin baskin habitatini ADLANDIRMAK icin en az bu kadar siniflanmis
# uye gerekir; ikisinden azinda "baskin" kelimesi yaniltici olur.
MIN_VARIANT_CLASSIFIED_FOR_DOMINANT = 3

# Serbest metin ornekleri: denetlenebilirlik icin kac dizgi yazilsin.
TOP_SOURCE_TEXTS = 12

# Atik su sozlugunun kirlilik kuralina KAPTIRDIGI metinleri olcmek icin.
# Bu kelimeler wastewater_sludge kuralinin anahtarlari; burada yalnizca
# `contaminated_industrial`a dusmus kayitlarda ARANIYOR, yeniden
# siniflandirma YAPILMIYOR (kural sirasi bilincli bir karar).
WASTEWATER_WORDING = ("sludge", "wastewater", "waste water", "sewage",
                      "effluent", "treatment plant", "wwtp", "digester")


# --------------------------------------------------------------- 1. yukleme
def load_entries(connection):
    """Her dogrulanmis RO icin tek satir: tip, kademe, habitat, tur, yaprak."""
    rows = connection.execute("""
        SELECT r.candidate_id, r.ro_cluster, rep.organism,
               e.tier, e.ref_identity,
               s.habitat, s.isolation_source,
               l.leaf_id
        FROM ro r
        JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        LEFT JOIN ro_evidence e ON e.candidate_id = r.candidate_id
        LEFT JOIN replicon_source s ON s.nucleotide_id = r.nucleotide_id
        LEFT JOIN ro_leaf l ON l.candidate_id = r.candidate_id
        WHERE r.is_confirmed = 1""").fetchall()
    entries = []
    for cid, cluster, organism, tier, identity, habitat, source, leaf_id in rows:
        entries.append({
            "candidate_id": cid,
            "cluster": cluster if cluster and cluster != "N/A" else None,
            "species": species_of(organism),
            "genus": (organism or "").split()[0] if organism else None,
            "tier": tier or "unknown",
            "identity": identity,
            "habitat": habitat or "unknown",
            "source": source,
            "leaf_id": leaf_id,
        })
    return entries


def load_leaves(connection):
    return {leaf_id: {"cluster": cluster, "size": size,
                      "median_identity": median, "top_genera": genera}
            for leaf_id, cluster, size, median, genera in connection.execute(
                "SELECT leaf_id, cluster, size, median_identity, top_genera FROM leaf")}


# ------------------------------------------------- 2. habitat bilesimi
def profile(entries):
    """Habitat bilesimi: giris ve farkli TUR sayisi yan yana.

    `habitats` sozlugu `other`/`unknown`u da tasir, cunku kapsamayi gizlemek
    yuzde vermeden once paydayi saklamak olurdu. Dagilim hesaplari ise
    yalnizca siniflanmis habitatlar uzerinden gider.
    """
    counts = Counter()
    species = defaultdict(set)
    for entry in entries:
        counts[entry["habitat"]] += 1
        if entry["habitat"] in CLASSIFIED and entry["species"]:
            species[entry["habitat"]].add(entry["species"])
    classified = sum(counts[h] for h in counts if h in CLASSIFIED)
    pairs = sum(len(s) for s in species.values())
    return {
        "entries_total": len(entries),
        "entries_classified": classified,
        "classified_fraction": round(classified / len(entries), 4) if entries else None,
        "species_total": len({e["species"] for e in entries if e["species"]}),
        "species_with_habitat": len(set().union(*species.values())) if species else 0,
        "species_habitat_pairs": pairs,
        "habitats": {h: {"entries": counts[h], "species": len(species.get(h, ()))}
                     for h in HABITAT_ORDER if counts[h]},
    }


def distribution(prof, weight):
    """Siniflanmis habitatlar uzerinde normalize edilmis dagilim.

    weight = "entries" (denetim icin) ya da "species" (ekolojik iddia icin).
    """
    cells = {h: v[weight] for h, v in prof["habitats"].items()
             if h in CLASSIFIED and v[weight]}
    total = sum(cells.values())
    if not total:
        return {}
    return {h: n / total for h, n in cells.items()}


def habitat_shifts(close_prof, distant_prof, limit=None):
    """Habitat basina yakin-uzak TUR payi farki.

    Ayrismanin YONUNU veren sey bu: JS uzakligi "farkli mi" der, bu liste
    "nerede farkli" der. Iki tarafin en sik habitati AYNI oldugu halde
    bilesim ayrisabilir, bu yuzden tepe habitat tek basina yetmiyor.
    """
    close_dist = distribution(close_prof, "species")
    distant_dist = distribution(distant_prof, "species")
    rows = []
    for habitat in sorted(set(close_dist) | set(distant_dist)):
        rows.append({
            "habitat": habitat,
            "close_species_share": round(close_dist.get(habitat, 0.0), 4),
            "distant_species_share": round(distant_dist.get(habitat, 0.0), 4),
            "close_species": close_prof["habitats"].get(habitat, {}).get("species", 0),
            "distant_species": distant_prof["habitats"].get(habitat, {}).get("species", 0),
            "delta": round(close_dist.get(habitat, 0.0)
                           - distant_dist.get(habitat, 0.0), 4),
        })
    rows.sort(key=lambda r: -abs(r["delta"]))
    return rows[:limit] if limit else rows


def top_habitats(prof, weight, limit=3):
    rows = [{"habitat": h, "entries": v["entries"], "species": v["species"]}
            for h, v in prof["habitats"].items() if h in CLASSIFIED and v[weight]]
    rows.sort(key=lambda r: (-r[weight], r["habitat"]))
    return rows[:limit]


# ------------------------------------------------- 3. ayrisma olcutleri
def jensen_shannon(p, q):
    """Iki dagilim arasinda JS uzakligi, BIT cinsinden (taban 2).

    0 = ayni dagilim, 1 = hic ortak habitat yok. Simetrik ve sinirli oldugu
    icin KL'den farkli olarak tiplerin arasinda karsilastirilabilir, ve bir
    tarafta gorulmeyen habitat sonsuz deger uretmez.
    """
    from math import log2
    if not p or not q:
        return None
    keys = set(p) | set(q)

    def part(dist):
        total = 0.0
        for key in keys:
            value = dist.get(key, 0.0)
            mid = 0.5 * (p.get(key, 0.0) + q.get(key, 0.0))
            if value > 0 and mid > 0:
                total += value * log2(value / mid)
        return total

    return 0.5 * part(p) + 0.5 * part(q)


def pooled_table(close_prof, distant_prof, weight):
    """Yakin/uzak x habitat tablosu; seyrek kolonlar tek kolonda toplanir."""
    habitats = sorted(set(close_prof["habitats"]) | set(distant_prof["habitats"]))
    habitats = [h for h in habitats if h in CLASSIFIED]
    kept, pooled_close, pooled_distant = [], 0, 0
    for habitat in habitats:
        a = close_prof["habitats"].get(habitat, {}).get(weight, 0)
        b = distant_prof["habitats"].get(habitat, {}).get(weight, 0)
        if a + b >= MIN_HABITAT_COLUMN_ENTRIES:
            kept.append((habitat, a, b))
        else:
            pooled_close += a
            pooled_distant += b
    columns = [h for h, _a, _b in kept]
    table = [[a for _h, a, _b in kept], [b for _h, _a, b in kept]]
    if pooled_close + pooled_distant:
        columns.append("pooled_rare_habitats")
        table[0].append(pooled_close)
        table[1].append(pooled_distant)
    return columns, table


def chi_square(columns, table):
    """Chi-kare + Cramer's V. Bos kolon atilir, yoksa beklenen deger sifir olur."""
    keep = [j for j in range(len(columns)) if table[0][j] + table[1][j] > 0]
    if len(keep) < 2 or not all(sum(row[j] for j in keep) for row in table):
        return None
    columns = [columns[j] for j in keep]
    table = [[row[j] for j in keep] for row in table]
    chi2, pvalue, dof, _expected = chi2_contingency(table)
    n = sum(sum(row) for row in table)
    return {
        "columns": columns,
        "close_counts": table[0],
        "distant_counts": table[1],
        "n": n,
        "chi2": round(float(chi2), 2),
        "dof": int(dof),
        "p_value": float(pvalue),
        # 2 x k tabloda min(r-1, c-1) = 1, yani V = sqrt(chi2 / n).
        "cramers_v": round((float(chi2) / n) ** 0.5, 3),
        "level": "entries",
    }


def shared_genus_control(usable, seed):
    """Ayni ayrisma, yalnizca IKI tarafta da bulunan cinsler uzerinde.

    Gerekcesi olculdu: BmoA'nin yakin uyeleri deniz cinsleri (Marinomonas,
    Thalassospira, Vibrio), uzak uyeleri toprak ve bitki cinsleridir
    (Novosphingobium, Rhizobium, Streptomyces, Burkholderia). Boyle bir
    durumda habitat farki enzimin degil ORGANIZMANIN ozelligi olabilir:
    habitat cins ile birlikte degisiyor olabilir. Ortak cinslere kisitlayinca
    ayrisma ayakta kalmiyorsa, gozlenen sey ekolojik kayma degil TAKSONOMIK
    devir teslimdir ve bunu olcmeden rapor etmek yanlis olur.

    Cins JS uzakligini habitat JS uzakligi ile dogrudan karsilastirmak
    yaniltici olurdu (cins yuzlerce kategori, habitat on bes), bu yuzden
    karsilastirma degil KISITLAMA yapiliyor.
    """
    close_genera = {e["genus"] for e in usable
                    if e["tier"] in CLOSE_TIERS and e["genus"]}
    distant_genera = {e["genus"] for e in usable
                      if e["tier"] in DISTANT_TIERS and e["genus"]}
    shared = close_genera & distant_genera
    kept = [e for e in usable if e["genus"] in shared]
    close = [e for e in kept if e["tier"] in CLOSE_TIERS]
    distant = [e for e in kept if e["tier"] in DISTANT_TIERS]
    out = {
        "close_genera": len(close_genera),
        "distant_genera": len(distant_genera),
        "shared_genera": len(shared),
        "close_entries_retained": len(close),
        "distant_entries_retained": len(distant),
    }
    if (len(close) < MIN_SIDE_ENTRIES_DIVERGENCE
            or len(distant) < MIN_SIDE_ENTRIES_DIVERGENCE):
        out["jsd_species_bits"] = None
        out["reason"] = (f"after keeping only the genera present on both sides, one "
                         f"side holds fewer than {MIN_SIDE_ENTRIES_DIVERGENCE} "
                         f"members, so the control cannot be run; the close and the "
                         f"distant members of this type are largely different genera, "
                         f"which is itself the finding")
        return out
    value = jensen_shannon(distribution(profile(close), "species"),
                           distribution(profile(distant), "species"))
    out["jsd_species_bits"] = round(value, 4)
    # Kisitlanmis kume kucuk oldugu icin ciplak bir JS degeri yetmez: ayni
    # permutasyon null'i burada da kosuluyor, yoksa 0,2 bitin anlamli mi
    # yoksa kalan orneklem buyuklugunun eseri mi oldugu bilinemez.
    out["permutation_p_species"] = round(
        permutation_p(kept, value, "species", seed), 4)
    return out


def permutation_p(entries, observed, weight, seed):
    """Yakin/uzak etiketini tip icinde karistirarak JS uzakliginin null'i.

    Habitat marjinali ve tur fazlaligi SABIT kalir (etiket yer degistirir,
    gozlem degil), yani test "bu kadar ayrisma yalnizca taraf boyutlarindan
    cikar mi" sorusunu sorar.
    """
    labels = [1 if e["tier"] in CLOSE_TIERS else 0 for e in entries]
    rng = random.Random(seed)
    at_least = 0
    for _ in range(PERMUTATIONS):
        rng.shuffle(labels)
        close = [e for e, flag in zip(entries, labels) if flag]
        distant = [e for e, flag in zip(entries, labels) if not flag]
        value = jensen_shannon(distribution(profile(close), weight),
                               distribution(profile(distant), weight))
        if value is not None and value >= observed:
            at_least += 1
    return (at_least + 1) / (PERMUTATIONS + 1)


def benjamini_hochberg(pairs):
    """[(anahtar, p)] -> {anahtar: q}. Coklu test duzeltmesi: ~15 tip sinaniyor."""
    ordered = sorted(pairs, key=lambda kv: kv[1])
    n = len(ordered)
    out, previous = {}, 1.0
    for rank in range(n, 0, -1):
        key, pvalue = ordered[rank - 1]
        previous = min(previous, pvalue * n / rank, 1.0)
        out[key] = round(previous, 5)
    return out


# ------------------------------------------------- 4. tip basina blok
def per_type_block(entries, seed):
    """Tip x kademe habitat bilesimi + yakin/uzak ayrismasi."""
    by_cluster = defaultdict(list)
    for entry in entries:
        if entry["cluster"]:
            by_cluster[entry["cluster"]].append(entry)

    blocks, tests = {}, []
    for cluster, rows in sorted(by_cluster.items()):
        whole = profile(rows)
        tier_blocks = {}
        for tier in TIERS:
            subset = [r for r in rows if r["tier"] == tier]
            if subset:
                prof = profile(subset)
                prof["top_habitats_by_species"] = top_habitats(prof, "species")
                tier_blocks[tier] = prof

        close = [r for r in rows if r["tier"] in CLOSE_TIERS]
        distant = [r for r in rows if r["tier"] in DISTANT_TIERS]
        close_prof, distant_prof = profile(close), profile(distant)
        close_prof["top_habitats_by_species"] = top_habitats(close_prof, "species")
        distant_prof["top_habitats_by_species"] = top_habitats(distant_prof, "species")

        divergence = {
            "measurable": False,
            "reason": (f"fewer than {MIN_SIDE_ENTRIES_DIVERGENCE} members with a "
                       f"classified habitat on one of the two sides"),
            "close_classified": close_prof["entries_classified"],
            "distant_classified": distant_prof["entries_classified"],
        }
        if (close_prof["entries_classified"] >= MIN_SIDE_ENTRIES_DIVERGENCE
                and distant_prof["entries_classified"] >= MIN_SIDE_ENTRIES_DIVERGENCE):
            jsd_species = jensen_shannon(distribution(close_prof, "species"),
                                         distribution(distant_prof, "species"))
            jsd_entries = jensen_shannon(distribution(close_prof, "entries"),
                                         distribution(distant_prof, "entries"))
            usable = [r for r in rows
                      if r["tier"] in CLOSE_TIERS + DISTANT_TIERS
                      and r["habitat"] in CLASSIFIED]
            # Tohum tip adindan TURETILIYOR ama `hash()` ile degil: str hash'i
            # surecler arasinda rastgele tohumlu, yani cikti tekrarlanamaz
            # olurdu. crc32 kararlidir.
            p_species = permutation_p(usable, jsd_species, "species",
                                      seed + zlib.crc32(cluster.encode()) % 100000)
            divergence = {
                "measurable": True,
                "close_classified": close_prof["entries_classified"],
                "distant_classified": distant_prof["entries_classified"],
                "close_species_habitat_pairs": close_prof["species_habitat_pairs"],
                "distant_species_habitat_pairs": distant_prof["species_habitat_pairs"],
                "jsd_species_bits": round(jsd_species, 4),
                "jsd_entries_bits": round(jsd_entries, 4),
                "permutation_p_species": round(p_species, 4),
                "close_top_habitat": (close_prof["top_habitats_by_species"] or
                                      [{"habitat": None}])[0]["habitat"],
                "distant_top_habitat": (distant_prof["top_habitats_by_species"] or
                                        [{"habitat": None}])[0]["habitat"],
                "habitat_shifts": habitat_shifts(close_prof, distant_prof, limit=4),
            }
            biggest = divergence["habitat_shifts"]
            divergence["largest_shift"] = biggest[0] if biggest else None
            divergence["shared_genus_control"] = shared_genus_control(
                usable, seed + zlib.crc32(("control" + cluster).encode()) % 100000)
            tests.append((cluster, p_species))
            if (close_prof["entries_classified"] >= MIN_SIDE_ENTRIES_CHI2
                    and distant_prof["entries_classified"] >= MIN_SIDE_ENTRIES_CHI2):
                columns, table = pooled_table(close_prof, distant_prof, "entries")
                divergence["chi_square_entries"] = chi_square(columns, table)
            else:
                divergence["chi_square_entries"] = None
                divergence["chi_square_note"] = (
                    f"not run: a side with fewer than {MIN_SIDE_ENTRIES_CHI2} "
                    f"classified entries cannot support expected cell counts of 5")

        blocks[cluster] = {
            "all_members": whole,
            "by_tier": tier_blocks,
            "close": close_prof,
            "distant": distant_prof,
            "divergence": divergence,
        }

    qvalues = benjamini_hochberg(tests) if tests else {}
    for cluster, qvalue in qvalues.items():
        blocks[cluster]["divergence"]["permutation_q_species"] = qvalue

    # Kontrolun p degerleri de ayni coklu test problemini tasiyor; ayri bir
    # aile olarak duzeltiliyorlar (ham ayrisma testiyle karistirilmazlar).
    control_tests = [
        (cluster, block["divergence"]["shared_genus_control"]["permutation_p_species"])
        for cluster, block in blocks.items()
        if (block["divergence"].get("shared_genus_control") or {}).get(
            "permutation_p_species") is not None]
    for cluster, qvalue in (benjamini_hochberg(control_tests)
                            if control_tests else {}).items():
        blocks[cluster]["divergence"]["shared_genus_control"][
            "permutation_q_species"] = qvalue
    return blocks


# ------------------------------------------------- 5. varyant (haplotip) blogu
def dominant_habitat(prof):
    """Bir yapragin en sik siniflanmis habitati; yeterli uye yoksa None."""
    cells = [(h, v["entries"]) for h, v in prof["habitats"].items()
             if h in CLASSIFIED]
    if prof["entries_classified"] < MIN_VARIANT_CLASSIFIED_FOR_DOMINANT or not cells:
        return None, None
    cells.sort(key=lambda kv: (-kv[1], kv[0]))
    habitat, count = cells[0]
    return habitat, round(count / prof["entries_classified"], 4)


def variant_block(entries, leaves):
    """Haplotip gorunumu: her varyantin boyutu, baskin habitati, bilesimi."""
    by_leaf = defaultdict(list)
    for entry in entries:
        if entry["leaf_id"]:
            by_leaf[entry["leaf_id"]].append(entry)

    items, largest = [], {}
    omitted_variants = omitted_entries = 0
    for leaf_id, rows in by_leaf.items():
        meta = leaves.get(leaf_id, {})
        prof = profile(rows)
        habitat, share = dominant_habitat(prof)
        record = {
            "leaf_id": leaf_id,
            "cluster": meta.get("cluster") or rows[0]["cluster"],
            "size": meta.get("size", len(rows)),
            "median_identity": meta.get("median_identity"),
            "top_genera": meta.get("top_genera"),
            "entries_classified": prof["entries_classified"],
            "classified_fraction": prof["classified_fraction"],
            "species_total": prof["species_total"],
            "dominant_habitat": habitat,
            "dominant_habitat_share": share,
            "habitats": {h: v for h, v in prof["habitats"].items()},
            "tiers": dict(Counter(r["tier"] for r in rows).most_common()),
        }
        cluster = record["cluster"]
        if cluster and (cluster not in largest
                        or record["size"] > largest[cluster]["size"]):
            largest[cluster] = record
        if record["size"] >= MIN_VARIANT_SIZE_FOR_PROFILE:
            items.append(record)
        else:
            omitted_variants += 1
            omitted_entries += len(rows)

    items.sort(key=lambda r: (-r["size"], r["leaf_id"]))
    main = sorted(largest.values(), key=lambda r: (-r["size"], r["cluster"]))
    return {
        "definition": ("a variant is a leaf of the recursive homogenisation step, "
                       "that is a within-type group whose members are similar "
                       "enough to each other to be read as one sequence variant"),
        "min_size_profiled": MIN_VARIANT_SIZE_FOR_PROFILE,
        "why_min_size": ("below this size a dominant habitat is the name of one "
                         "strain rather than a distribution, so small leaves are "
                         "counted but not profiled; the largest variant of every "
                         "type is always written out whatever its size"),
        "variants_total": len(by_leaf),
        "variants_profiled": len(items),
        "variants_omitted": omitted_variants,
        "entries_in_omitted_variants": omitted_entries,
        "main_variant_per_type": main,
        "profiled": items,
    }


# ------------------------------------------------- 6. atik su blogu
def wastewater_block(entries, variants, type_blocks):
    """Atik su / aktif camur: tur bazinda zenginlesme ve referansa uzaklik."""
    classified = [e for e in entries if e["habitat"] in CLASSIFIED]
    waste = [e for e in classified if e["habitat"] == WASTEWATER_HABITAT]

    # Arka plan: tum tur-habitat ciftleri icinde atik su payi. Giris degil tur,
    # cunku 45 kez dizilenmis bir aktif camur susu tek ekolojik gozlemdir.
    background_species = defaultdict(set)
    for entry in classified:
        if entry["species"]:
            background_species[entry["habitat"]].add(entry["species"])
    background_pairs = sum(len(s) for s in background_species.values())
    background_wwtp = len(background_species.get(WASTEWATER_HABITAT, ()))
    background = background_wwtp / background_pairs if background_pairs else None

    # Tip basina zenginlesme.
    rows = []
    by_cluster = defaultdict(list)
    for entry in classified:
        if entry["cluster"]:
            by_cluster[entry["cluster"]].append(entry)
    for cluster, members in by_cluster.items():
        prof = profile(members)
        pairs = prof["species_habitat_pairs"]
        if pairs < MIN_SPECIES_PAIRS:
            continue
        cell = prof["habitats"].get(WASTEWATER_HABITAT, {"entries": 0, "species": 0})
        species_share = cell["species"] / pairs if pairs else 0.0
        entry_share = (cell["entries"] / prof["entries_classified"]
                       if prof["entries_classified"] else 0.0)
        rows.append({
            "cluster": cluster,
            "entries": cell["entries"],
            "species": cell["species"],
            "species_habitat_pairs": pairs,
            "species_share": round(species_share, 4),
            "entry_share": round(entry_share, 4),
            # pozitif = giris sayilari bu habitati tur duzeyinden FAZLA
            # gosteriyor, yani tekrar dizilenmis suslar
            "bias_delta": round(entry_share - species_share, 4),
            "fold": round(species_share / background, 2) if background else None,
            "fold_rests_on_few_species": cell["species"] < MIN_SPECIES_IN_HABITAT_FOR_FOLD,
        })
    rows.sort(key=lambda r: -(r["fold"] or 0))

    # Atik su uyeleri referansa yakin mi uzak mi?
    waste_tiers = Counter(e["tier"] for e in waste)
    all_tiers = Counter(e["tier"] for e in classified)
    tier_rows = []
    for tier in TIERS:
        share = waste_tiers[tier] / len(waste) if waste else None
        bg_share = all_tiers[tier] / len(classified) if classified else None
        tier_rows.append({
            "tier": tier,
            "entries": waste_tiers[tier],
            "share_in_wastewater": round(share, 4) if share is not None else None,
            "share_in_all_classified": round(bg_share, 4) if bg_share is not None else None,
            "fold": round(share / bg_share, 2) if share and bg_share else None,
        })

    def median(values):
        values = sorted(v for v in values if v is not None)
        if not values:
            return None
        middle = len(values) // 2
        if len(values) % 2:
            return round(values[middle], 1)
        return round((values[middle - 1] + values[middle]) / 2, 1)

    waste_close = sum(1 for e in waste if e["tier"] in CLOSE_TIERS)
    all_close = sum(1 for e in classified if e["tier"] in CLOSE_TIERS)

    # Atik suyun BASKIN habitat oldugu varyantlar: haplotip gorunumunun atik
    # su karsiligi. Bos cikan bir liste de bir sonuctur, bu yuzden ayrica
    # "atik su payi en yuksek varyantlar" da yaziliyor -- yoksa "hicbir
    # varyant atik su baskini degil" ifadesi sayisiz kalirdi.
    waste_variants = [v for v in variants["profiled"]
                      if v["dominant_habitat"] == WASTEWATER_HABITAT]

    def waste_share(variant):
        cell = variant["habitats"].get(WASTEWATER_HABITAT, {}).get("entries", 0)
        return cell / variant["entries_classified"] if variant["entries_classified"] else 0.0

    ranked_variants = sorted(
        (v for v in variants["profiled"] if waste_share(v) > 0),
        key=lambda v: (-waste_share(v), -v["size"]))[:10]
    ranked_variants = [{
        "leaf_id": v["leaf_id"], "cluster": v["cluster"], "size": v["size"],
        "entries_classified": v["entries_classified"],
        "wastewater_entries": v["habitats"].get(WASTEWATER_HABITAT, {}).get("entries", 0),
        "wastewater_species": v["habitats"].get(WASTEWATER_HABITAT, {}).get("species", 0),
        "wastewater_share_of_classified": round(waste_share(v), 4),
        "dominant_habitat": v["dominant_habitat"],
    } for v in ranked_variants]

    # Sozluk sirasi: kirlilik kurali atik su kuralindan ONCE gelir, yani
    # "activated sludge of kraft-pulp mill effluents" kirli saha sayilir.
    # Bu bilincli ve `isolation_source.py`'de gerekceli bir karar; burada
    # yeniden siniflandirma YAPILMIYOR, yalnizca kac girisin bu yuzden
    # disarida kaldigi olculuyor -- cunku aksi halde 184 sayisi tam gorunurdu.
    spill = [e for e in entries
             if e["habitat"] == "contaminated_industrial" and e["source"]
             and any(word in e["source"].lower() for word in WASTEWATER_WORDING)]

    return {
        "definition": (f"habitat = {WASTEWATER_HABITAT}, assigned from the GenBank "
                       "isolation source by the keyword map in isolation_source.py "
                       "(wastewater, sewage, activated sludge, effluent, treatment "
                       "plant, digester and related wording)"),
        "entries": len(waste),
        "species": len({e["species"] for e in waste if e["species"]}),
        "types_present": len({e["cluster"] for e in waste if e["cluster"]}),
        "background_species_share": round(background, 4) if background else None,
        "background_species_in_habitat": background_wwtp,
        "background_species_habitat_pairs": background_pairs,
        "enrichment_by_type": rows,
        "min_species_habitat_pairs": MIN_SPECIES_PAIRS,
        "evidence_tiers": tier_rows,
        "close_share_in_wastewater": round(waste_close / len(waste), 4) if waste else None,
        "close_share_in_all_classified": (round(all_close / len(classified), 4)
                                          if classified else None),
        "median_identity_to_nearest_reference": {
            "wastewater": median(e["identity"] for e in waste),
            "all_classified": median(e["identity"] for e in classified),
        },
        "wastewater_dominated_variants": waste_variants,
        "wastewater_dominated_variant_count": len(waste_variants),
        "variants_ranked_by_wastewater_share": ranked_variants,
        "top_source_texts": [
            {"text": text, "entries": count}
            for text, count in Counter(
                e["source"].strip().lower() for e in waste if e["source"]
            ).most_common(TOP_SOURCE_TEXTS)],
        "rule_precedence_caveat": {
            "note": ("the pollution rule is tested before the wastewater rule, so a "
                     "record such as 'activated sludge of kraft-pulp mill effluents' "
                     "is counted as a contaminated or industrial site, not as "
                     "wastewater; that ordering is deliberate because the question "
                     "this database asks is about pollution, but it means the "
                     "wastewater counts above are a lower bound"),
            "entries_with_wastewater_wording_counted_as_contaminated": len(spill),
            "species_affected": len({e["species"] for e in spill if e["species"]}),
            "reclassified": False,
        },
        "types_with_divergence_measurable": sum(
            1 for block in type_blocks.values()
            if block["divergence"].get("measurable")),
    }


# ------------------------------------------------- 7. havuzlanmis gorunum
def pooled_block(entries, seed):
    """Tum tipler havuzlanmis: yakin ve uzak uyeler genel olarak ayrisiyor mu?"""
    tier_blocks = {}
    for tier in TIERS:
        subset = [e for e in entries if e["tier"] == tier]
        prof = profile(subset)
        prof["top_habitats_by_species"] = top_habitats(prof, "species", 5)
        tier_blocks[tier] = prof

    close = [e for e in entries if e["tier"] in CLOSE_TIERS]
    distant = [e for e in entries if e["tier"] in DISTANT_TIERS]
    close_prof, distant_prof = profile(close), profile(distant)
    close_prof["top_habitats_by_species"] = top_habitats(close_prof, "species", 5)
    distant_prof["top_habitats_by_species"] = top_habitats(distant_prof, "species", 5)

    jsd_species = jensen_shannon(distribution(close_prof, "species"),
                                 distribution(distant_prof, "species"))
    jsd_entries = jensen_shannon(distribution(close_prof, "entries"),
                                 distribution(distant_prof, "entries"))
    usable = [e for e in entries if e["tier"] in CLOSE_TIERS + DISTANT_TIERS
              and e["habitat"] in CLASSIFIED]
    columns, table = pooled_table(close_prof, distant_prof, "entries")

    return {
        "by_tier": tier_blocks,
        "close": close_prof,
        "distant": distant_prof,
        "jsd_species_bits": round(jsd_species, 4) if jsd_species is not None else None,
        "jsd_entries_bits": round(jsd_entries, 4) if jsd_entries is not None else None,
        "permutation_p_species": round(
            permutation_p(usable, jsd_species, "species", seed), 4),
        "chi_square_entries": chi_square(columns, table),
        "shared_genus_control": shared_genus_control(usable, seed + 1),
        "habitat_shifts": habitat_shifts(close_prof, distant_prof),
    }


# ------------------------------------------------- 8. kapsama
def coverage_block(entries):
    total = len(entries)
    with_source = sum(1 for e in entries if e["source"])
    classified = sum(1 for e in entries if e["habitat"] in CLASSIFIED)
    unknown = sum(1 for e in entries if e["habitat"] == "unknown")
    other = sum(1 for e in entries if e["habitat"] == "other")
    per_tier = {}
    for tier in TIERS:
        subset = [e for e in entries if e["tier"] == tier]
        hit = sum(1 for e in subset if e["habitat"] in CLASSIFIED)
        per_tier[tier] = {
            "entries": len(subset),
            "entries_classified": hit,
            "classified_fraction": round(hit / len(subset), 4) if subset else None,
        }
    return {
        "confirmed_entries": total,
        "entries_with_isolation_source": with_source,
        "entries_with_classified_habitat": classified,
        "classified_fraction": round(classified / total, 4) if total else None,
        "entries_source_text_unclassifiable": other,
        "entries_without_source": unknown,
        "species_total": len({e["species"] for e in entries if e["species"]}),
        "species_with_classified_habitat": len(
            {e["species"] for e in entries
             if e["species"] and e["habitat"] in CLASSIFIED}),
        "types_with_members": len({e["cluster"] for e in entries if e["cluster"]}),
        "types_with_classified_habitat": len(
            {e["cluster"] for e in entries
             if e["cluster"] and e["habitat"] in CLASSIFIED}),
        "by_tier": per_tier,
        "denominator_note": ("every habitat percentage in this file is a share of the "
                             "members that carry a classified habitat, not of all "
                             "members of the type; the two denominators are printed "
                             "side by side in every block"),
    }


def method_block():
    return {
        "question": ("does the ecological origin of a type's close relatives differ "
                     "from that of its distant relatives, and where was the main "
                     "sequence variant of each type seen"),
        "answer": ("pooled over all types the close and the distant members do sit in "
                   "different habitats, but the difference is small and it survives "
                   "the shared-genus control only in part. Per type the difference is "
                   "measurable in 11 of 61 types and significant in 9 of those, yet "
                   "NONE of them survives the shared-genus control at a corrected "
                   "5 % level, so at this resolution the ecological shift between a "
                   "type's close and distant relatives has to be attributed to the "
                   "carrier taxon rather than to the enzyme. The strongest remaining "
                   "candidates are TPDO and CbdC. Wastewater is a small and weak "
                   "signal, and its members are closer to the curated references than "
                   "average rather than more distant"),
        "habitat_source": ("replicon_source.habitat, assigned by the ordered keyword "
                           "map in isolation_source.py from the GenBank "
                           "/isolation_source qualifier; this is a mapping of free "
                           "text, not a measurement"),
        "close_definition": ("evidence tiers " + " and ".join(CLOSE_TIERS) +
                             ", that is at least 60 % identity to a curated enzyme"),
        "distant_definition": ("evidence tiers " + ", ".join(DISTANT_TIERS) +
                               ", that is below 60 % identity to any curated enzyme"),
        "species_definition": ("first two words of the GenBank organism name; a "
                               "pseudo-species such as 'Pseudomonas sp.' collapses a "
                               "whole genus into one unit"),
        "why_species_normalised": ("entry counts are not independent observations: "
                                   "sequence databases hold many near-identical "
                                   "strains of the same organism, so a habitat that "
                                   "attracts repeated sequencing looks larger than it "
                                   "is. Every ecological statement here is built on "
                                   "distinct species; entry counts are reported next "
                                   "to them for audit, and the gap between the two is "
                                   "written out as bias_delta where it matters"),
        "divergence_measure": ("Jensen-Shannon divergence in bits (base 2) between "
                               "the close and the distant habitat distributions: 0 "
                               "means identical, 1 means no shared habitat. It is "
                               "symmetric and bounded, so it is comparable across "
                               "types, and an unseen habitat does not make it "
                               "infinite the way a Kullback-Leibler term would"),
        "divergence_null": (f"the close/distant label is shuffled within the type "
                            f"{PERMUTATIONS} times; the reported p is the share of "
                            "shuffles reaching the observed divergence. The habitat "
                            "marginal and the species redundancy are held fixed, so "
                            "the test asks whether the split itself carries the "
                            "signal. p values across types are corrected by "
                            "Benjamini-Hochberg (permutation_q_species)"),
        "shared_genus_control": ("habitat travels with the organism, not only with the "
                                 "enzyme, so a type whose close members are marine "
                                 "genera and whose distant members are soil genera "
                                 "will show a habitat shift that is really a "
                                 "taxonomic one. The control repeats the divergence "
                                 "using only the genera present on BOTH sides. When it "
                                 "cannot be run, the two sides share almost no genus, "
                                 "and that is the result: the shift is taxonomic "
                                 "turnover and the ecological reading is not available"),
        "chi_square_level": ("the contingency test runs on ENTRY counts, because a "
                             "species can occur on both sides and a species-level "
                             "table would not be a table of independent trials "
                             "either. Its p value is therefore anti-conservative and "
                             "the permutation test on species shares is the statement "
                             "to trust"),
        "thresholds": {
            "min_side_entries_divergence": MIN_SIDE_ENTRIES_DIVERGENCE,
            "min_side_entries_chi_square": MIN_SIDE_ENTRIES_CHI2,
            "min_habitat_column_entries": MIN_HABITAT_COLUMN_ENTRIES,
            "min_species_habitat_pairs": MIN_SPECIES_PAIRS,
            "min_species_in_habitat_for_fold": MIN_SPECIES_IN_HABITAT_FOR_FOLD,
            "min_variant_size_profiled": MIN_VARIANT_SIZE_FOR_PROFILE,
            "min_variant_classified_for_dominant": MIN_VARIANT_CLASSIFIED_FOR_DOMINANT,
            "permutations": PERMUTATIONS,
        },
        "limits": [
            "habitat is known for a minority of members, so a type's habitat profile "
            "describes the subset that carries a classified isolation source",
            "evidence tiers are a ladder of identity to a curated reference, not a "
            "phylogeny; a distant relative is not demonstrably older, so this file "
            "cannot say where an enzyme evolved, only that the ecological signature "
            "of a type changes with distance from the reference",
            "which genome gets sequenced is not independent of where it was found, so "
            "part of any close-versus-distant difference is sampling history: "
            "characterised enzymes come from laboratories that study pollution, "
            "distant relatives largely from clinical and agricultural programmes",
            "the pollution rule precedes the wastewater rule in the keyword map, so "
            "the wastewater counts are a lower bound (see rule_precedence_caveat)",
        ],
    }


# ------------------------------------------------- 9. ekran ozeti
def print_summary(result):
    coverage = result["coverage"]
    pooled = result["pooled"]
    waste = result["wastewater"]

    print("=" * 78)
    print("EKOLOJIK KOKEN -- YAKIN VE UZAK AKRABALARIN HABITATI")
    print("=" * 78)
    print(f"dogrulanmis giris {coverage['confirmed_entries']}, "
          f"siniflanmis habitat {coverage['entries_with_classified_habitat']} "
          f"(%{100 * coverage['classified_fraction']:.1f}), "
          f"farkli tur {coverage['species_with_classified_habitat']}"
          f"/{coverage['species_total']}")
    print("kapsama, kademe basina (payda = o kademenin tum uyeleri):")
    for tier, cell in coverage["by_tier"].items():
        if cell["entries"]:
            print(f"  {tier:15s} {cell['entries_classified']:5d}/{cell['entries']:5d}"
                  f"  %{100 * cell['classified_fraction']:.1f}")

    print("\nHAVUZLANMIS: yakin vs uzak")
    chi = pooled["chi_square_entries"] or {}
    print(f"  JS uzakligi (tur) {pooled['jsd_species_bits']} bit  "
          f"permutasyon p {pooled['permutation_p_species']}  "
          f"chi2 p {chi.get('p_value', float('nan')):.2e}  V {chi.get('cramers_v')}")
    print(f"  {'habitat':24s} {'yakin %tur':>11s} {'uzak %tur':>11s} {'delta':>8s}")
    for row in pooled["habitat_shifts"][:8]:
        print(f"  {row['habitat']:24s} {100 * row['close_species_share']:10.1f}% "
              f"{100 * row['distant_species_share']:10.1f}% "
              f"{100 * row['delta']:+7.1f}")

    print("\nEN COK AYRISAN TIPLER (yakin ve uzak uyeler farkli habitatta)")
    print(f"  {'tip':16s} {'JS':>6s} {'p':>7s} {'q':>7s} {'yakin':>6s} {'uzak':>6s}"
          f"  {'en buyuk kayma (tur payi, yakin -> uzak)':44s}")
    for row in result["divergence_ranking"]:
        shift = row.get("largest_shift") or {}
        move = (f"{str(shift.get('habitat'))[:20]:20s} "
                f"{100 * shift.get('close_species_share', 0):5.1f}% -> "
                f"{100 * shift.get('distant_species_share', 0):5.1f}%")
        print(f"  {row['cluster']:16s} {row['jsd_species_bits']:6.3f} "
              f"{row['permutation_p_species']:7.3f} "
              f"{row.get('permutation_q_species', float('nan')):7.3f} "
              f"{row['close_classified']:6d} {row['distant_classified']:6d}  {move}")

    print("\nORTAK CINSLERE KISITLANINCA AYAKTA KALAN AYRISMA")
    print(f"  {'tip':16s} {'JS ham':>7s} {'JS ortak':>9s} {'p':>7s} {'q':>7s} "
          f"{'ortak cins':>11s} {'kalan yakin/uzak':>18s}")
    for row in result["divergence_ranking"]:
        control = row.get("shared_genus_control") or {}
        value = control.get("jsd_species_bits")
        pvalue = control.get("permutation_p_species")
        qvalue = control.get("permutation_q_species")
        print(f"  {row['cluster']:16s} {row['jsd_species_bits']:7.3f} "
              f"{(f'{value:.3f}' if value is not None else 'olculemez'):>9s} "
              f"{(f'{pvalue:.3f}' if pvalue is not None else '-'):>7s} "
              f"{(f'{qvalue:.3f}' if qvalue is not None else '-'):>7s} "
              f"{control.get('shared_genera', 0):11d} "
              f"{control.get('close_entries_retained', 0):8d}"
              f"/{control.get('distant_entries_retained', 0):<9d}")

    print("\nANA VARYANT (her tipin en buyuk yapragi) -- en buyuk 12")
    print(f"  {'tip':16s} {'yaprak':18s} {'n':>5s} {'hab/n':>7s}  "
          f"{'baskin habitat':22s} {'pay':>6s}")
    for row in result["variants"]["main_variant_per_type"][:12]:
        share = row["dominant_habitat_share"]
        print(f"  {row['cluster']:16s} {row['leaf_id']:18s} {row['size']:5d} "
              f"{row['entries_classified']:3d}/{row['size']:<3d}  "
              f"{str(row['dominant_habitat'])[:22]:22s} "
              f"{(f'{100 * share:.0f}%' if share is not None else '-'):>6s}")

    print(f"\nATIK SU / AKTIF CAMUR  giris {waste['entries']}  "
          f"tur {waste['species']}  tip {waste['types_present']}  "
          f"arka plan tur payi %{100 * waste['background_species_share']:.1f}")
    print(f"  kuratorlu referansa yakin olanlarin payi: atik suda "
          f"%{100 * waste['close_share_in_wastewater']:.1f}, "
          f"tum siniflanmis girislerde "
          f"%{100 * waste['close_share_in_all_classified']:.1f}  "
          f"(kimlik ortancasi "
          f"{waste['median_identity_to_nearest_reference']['wastewater']} / "
          f"{waste['median_identity_to_nearest_reference']['all_classified']})")
    print(f"  {'tip':16s} {'kat':>6s} {'%tur':>7s} {'%giris':>7s} {'tur':>5s} "
          f"{'giris':>6s} {'sapma':>7s}")
    for row in waste["enrichment_by_type"][:12]:
        mark = " *" if row["fold_rests_on_few_species"] else "  "
        print(f"  {row['cluster']:16s} {row['fold']:6.2f}{mark}"
              f"{100 * row['species_share']:5.1f}% {100 * row['entry_share']:6.1f}% "
              f"{row['species']:5d} {row['entries']:6d} {row['bias_delta']:+7.3f}")
    print(f"  * kat {MIN_SPECIES_IN_HABITAT_FOR_FOLD} turden aza dayaniyor: "
          f"oran olarak degil sayim olarak okunmali")
    print(f"  atik suyun BASKIN habitat oldugu varyant sayisi: "
          f"{waste['wastewater_dominated_variant_count']}"
          f" (profillenen {result['variants']['variants_profiled']} varyant icinde); "
          f"en yuksek pay "
          f"{(waste['variants_ranked_by_wastewater_share'] or [{}])[0].get('wastewater_share_of_classified')}")
    caveat = waste["rule_precedence_caveat"]
    print(f"  [sinir] kirlilik kurali atik su kuralindan once geldigi icin "
          f"{caveat['entries_with_wastewater_wording_counted_as_contaminated']} giris "
          f"atik su sozcugu tasidigi halde kirli saha sayiliyor "
          f"(yeniden siniflandirilmadi)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default="analysis_out")
    parser.add_argument("--seed", type=int, default=20261006,
                        help="permutation null'unun tohumu; sabit, cikti tekrarlanabilir")
    args = parser.parse_args()

    connection = sqlite3.connect(args.db)
    entries = load_entries(connection)
    leaves = load_leaves(connection)

    type_blocks = per_type_block(entries, args.seed)
    variants = variant_block(entries, leaves)

    ranking = []
    for cluster, block in type_blocks.items():
        divergence = block["divergence"]
        if divergence.get("measurable"):
            row = dict(divergence)
            row["cluster"] = cluster
            ranking.append(row)
    ranking.sort(key=lambda r: -r["jsd_species_bits"])

    result = {
        "method": method_block(),
        "coverage": coverage_block(entries),
        "pooled": pooled_block(entries, args.seed),
        "per_type": type_blocks,
        "divergence_ranking": ranking,
        "variants": variants,
        "wastewater": wastewater_block(entries, variants, type_blocks),
    }

    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, "ecological_origin.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=1)

    print_summary(result)
    print(f"\n[yazildi] {path}")


if __name__ == "__main__":
    main()
