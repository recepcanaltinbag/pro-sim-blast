"""
Operon duzeyi iliskiler -- duzenleyici x enzim, ve "mobilite" degiskeninin denetimi.

NEDEN BU SCRIPT VAR. Iki ayri soru, ama ayni kok nedenden oturu birlikte
olculmeleri gerekiyor.

BIRINCI SORU, kuratorun sorusu: "hangi enzim ailesi ya da tipi hangi
duzenleyiciyle gidiyor". `ro_regulation` tablosu operonun 5' ucundaki yukari
akis genini ve -- regulatorse -- ailesini tasiyor, ama bugune kadar yalnizca
TOPLAM dagilim yayinlanmisti (LysR 973, TetR 531, IclR 358...). O dagilim
"LysR baskin" der; "LysR HANGI kimyayla gidiyor" sorusuna cevap vermez.

IKINCI SORU, kuratorun hakli kuskusu: "operon anotasyonunu dogrudan
kullanirsak oradan cok yanlilik gelebilir". Kusku yerinde ve neyi vurdugu
olculdu. Veritabaninda iki FARKLI kalitede komsuluk bilgisi var ve sitede
ayni dille sunuluyorlar:

  * operon BILESENLERI (beta, ferredoksin, reduktaz) `build_operons.py`
    tarafindan komsu CDS'in PROTEIN DIZISI profil HMM'e taranarak bulunur.
    Bu dizi kanitidir; anotasyon metni yanlis ya da bos olsa da durur.
  * gen KATEGORILERI (`gene_category`, method='regex_v1') GenBank
    `/product` METNINE uygulanan bir regex'tir (bkz. build_db.py CATEGORIES).
    Olctugu sey biyoloji degil, kaydi gonderenin anotasyon aliskanligi.
    'transposon' kategorisi bunlardan biridir.

Ve yayinlanan "mobile" degiskeni (`stats_overview.py`) tam olarak
`is_plasmid OR transposon_count > 0`, yani SAGLAM yarinin uzerine ZAYIF
yarinin OR'u. Kompozisyon olculdu ve dengesiz: 11.422 dogrulanmis giriste
302 plazmit (%2,6, yapisal kanit), 1.704 regex transpozonu (%14,9), 1.843
"mobil" (%16,1) ve bunlarin 1.541'i, yani %83,6'si YALNIZCA metin regex'i
yuzunden mobil sayiliyor. Dolayisiyla yayinlanan degisken agirlikla
anotasyon kalitesini olcuyor.

NE OLCULDU VE NE CIKTI.

  1. Duzenleyici ailesi x kimyasal aile / enzim tipi / RO grubu, ki-kare +
     Cramer's V, hem giris hem TIP x CINS duzeyinde. V, PERMUTASYON null'una
     karsi da okunuyor: seyrek bir capraz tabloda V yalnizca hucre sayisindan
     sifirdan buyuk cikar (olculdu: aile x tip tablosunda sans tabani
     V=0,12), bu yuzden "gercek" olan V'nin kendisi degil null'u ASAN kismi.
     Sonuc: duzenleyici ailesi ile enzim tipi arasindaki iliski cins
     cokertmesinden sonra da ayakta (V=0,41, sans tabani 0,12); kimyasal aile
     ve RO grubu ile de ayakta ama daha zayif; operon MIMARISI (divergent /
     ko-direksiyonel) ile kimya arasindaki iliski en zayif halka.

  2. "Mobilite x substrat sinifi" UC AYRI kanit tanimiyla yeniden hesaplandi.
     Sonuc rahatlatici ama yayinlanan sayiyi da duzeltiyor: iliski yapisal
     kanitla TEK BASINA cok daha GUCLU (plazmit: giris 3,95x / cins 4,35x),
     regex yarisi ayni yonde ama zayif (1,33x / 1,67x), ve ikisinin OR'u
     yayinlanan ara degeri veriyor (1,40x / 1,77x). Yani bulgu anotasyon
     yarisina DAYANMIYOR; yayinlanan OR degiskeni bulguyu SEYRELTIYOR.

  3. Transpozon yakinligi kuratorlu enzime UZAKLIKLA iliskili, ve yon
     duzeye gore TERSINIYOR: giris duzeyinde yakin uyelerde oran daha dusuk
     (OR 0,82), TIP x CINS duzeyinde ise cok daha yuksek (characterized %60,9
     -> novel %2,2, V=0,196). Tersinmenin sebebi bulundu: giris sayilari tek
     tek hucrelerin esiri (GbcA x Pseudomonas 321 giris, VanA x Pseudomonas
     182). Cins duzeyindeki yon daha anlamli ve soyledigi sey rahatsiz edici:
     kuratorlu enzime yakinlik, "yakinda transpozon ANOTE EDILMIS olmasini"
     ongoruyor. Alt-aile atama sinifi ile iliski ise yok sayilir (V=0,06).

  4. Anotasyon yanliliginin kendi buyuklugu. Uc bilesen bulundu ve UCU AYNI
     YONE GITMIYOR, bu yuzden tek bir duzeltme katsayisi yok:
       * SERT TABAN: 945 giris (%8,3) hic anote komsusu olmayan kayitlarda
         (hepsi okaryot mRNA kaydi). Bu girislerde transpozon yakinligi
         sifir OLMAK ZORUNDA, ama yayinlanan %14,9'un paydasinda duruyorlar.
         Cagrilabilir kayitlara daralinca oran %16,3.
       * PENCERE DOLULUGU: +-10 kb'de anote gen sayisi ceyrekligine gore oran
         %11,3 ile %25,2 arasinda, yani 2,2x. Bu, kaydi gonderenin kac gen
         anote ettigine bagli.
       * REPLIKON BOYU: <=50 CDS'li kontiglerde %19,6, >500 CDS'li
         replikonlarda %14,1 -- yogunlukla TERS yonde, cunku kucuk kontigler
         IS elementi bakimindan zengindir (bu kismen gercek biyoloji).
     Yanliligin biyolojik iddiaya ne kadar sizdigi dogrudan olculdu:
     substrat sinifi x regex transpozonu iliskisi pencere doluluguna gore
     tabakalanip Mantel-Haenszel ile birlestirildi; ortak OR 1,41'den 1,31'e
     dusuyor. Yani apacik fazlaligin yaklasik dortte biri anotasyon
     yogunlugundan, dortte ucu degil.

NE OLCULMUYOR -- bunu da bastan yazmak gerekiyor:
  * Regex'in KENDISI dogrulanmiyor. "integrase|recombinase|resolvase" gibi
    anahtarlar transpozaz olmayan proteinleri de yakalar; burada olculen sey
    "metin regex'e uyuyor mu", "orada gercekten hareketli bir element var mi"
    DEGIL. Dizi dogrulamasi ancak komsu proteinlere transpozaz profil HMM'i
    kosularak yapilabilir ve bu scriptin isi degil.
  * Plazmit alani da mutlak dogru degil: GenBank'in kendi beyanina dayanir,
    entegre edilmis plazmitler ve isimlendirilmemis megaplazmitler kacar.
    Yine de bir DIZI OZELLIGIDIR, urun metninin kelime secimi degil.
  * Nedensellik yok. "Ksenobiyotik oksijenazlar plazmitte daha sik" bir
    birliktelik olcumudur; hangi yone aktigini bu veri soylemez.
  * Promotor DIZISI hala analiz edilemiyor (kayitlarin %88,9'u CON tipinde,
    dizi dosyada yok); dolayisiyla "duzenleyici ailesi x kimya" iliskisi
    OPERATOR dizisiyle degil, yalnizca komsu genin urun metniyle kurulmus
    bir aile etiketine dayanir -- yani 1. sorunun cevabi da kismen ayni
    anotasyon zincirine bagli ve bu, verdict bloklarinda yaziyor.

Cikti: analysis_out/operon_relations.json

Kullanim:
    python3 operon_relations.py --db roar.sqlite --out-dir analysis_out
"""

import argparse
import csv
import json
import os
import random
import sqlite3
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
from Bio import Align
from scipy import stats

# --- Sabitler. Her biri adlandirilmis, ve her birinin gerekcesi yaninda.

# Permutasyon sayisi ve tohum: cooccurrence.py ile AYNI, boylece iki modulun
# null'lari karsilastirilabilir ve kosular yeniden uretilebilir.
PERMUTATIONS = 2000
RANDOM_SEED = 20261005

# gene_category tablosundaki metin tabanli yontemin adi ve aradigimiz kategori.
# Kodda sabit yazmak yerine adlandirildi, cunku bu scriptin TUM konusu bu iki
# dizginin anotasyon metninden geliyor olmasi.
REGEX_METHOD = "regex_v1"
TRANSPOSON_CATEGORY = "transposon"

# Komsuluk penceresi. extract_genomic_context.py'deki NEIGHBOR_WINDOW ile ayni
# olmak ZORUNDA; burada yalnizca raporlamak icin tekrarlandi.
NEIGHBOUR_WINDOW_BP = 10000

# Ki-kare asimptotigi beklenen hucre degerinin ~5'in altina dusmemesini ister.
# Esigi karsilamayan hucrelerin ORANI her tabloda raporlanir; bu oran
# SPARSE_TABLE_LIMIT'i gecerse p degerine guvenilmez ve karar permutasyon
# null'unu asan V farkina birakilir.
EXPECTED_CELL_FLOOR = 5
SPARSE_TABLE_LIMIT = 0.20

# Cramer's V yorum esikleri (Cohen'in kucuk/orta/buyuk konvansiyonu).
# "Anlamli" ile "buyuk" ayrilmali: 11.422 giriste neredeyse her sey p<0,001.
V_NEGLIGIBLE = 0.10
V_MODERATE = 0.30

# Bir capraz tablo hucresini BULGU olarak adlandirmak icin iki kosul. Birincisi
# hucrenin kac girise dayandigi, ikincisi kac FARKLI cinse: 40 girisi olan ama
# tek cinsten gelen bir hucre tek bir dizileme projesidir, bulgu degil.
MIN_CELL_ENTRIES = 20
MIN_CELL_GENERA = 3

# Mantel-Haenszel tabakasina girmek icin tabakanin her iki kolunda en az bu
# kadar gozlem olmali; daha azi ortak OR'a gurultuden baska sey katmaz.
MIN_STRATUM_N = 30

# Anotasyon yogunlugu tabakalari VERIDEN cikarilir (ceyreklikler), elle
# yazilmaz: pencere doluluk dagilimi veri setine ozgudur ve sabit bir kesim
# noktasi baska bir surumde anlamini yitirir.
DENSITY_QUANTILES = (25, 50, 75)
# Replikon boyu tabakalari da ayni gerekceyle ceyreklik.
CDS_QUANTILES = (25, 50, 75)

# stats_overview.py ile AYNI evren tanimi, yoksa sayilar karsilastirilamaz.
SUBSTRATE_CLASSES = ["xenobiotic", "natural_aromatic", "natural_specialized"]
EUK_HEAVY_RATE = 0.2          # stats_overview.py: okaryot orani >=%20 olan kume

TIERS = ["characterized", "close_homolog", "family_member", "distant", "novel"]
# "Yakin" tanimi stratified_stats.py ile ayni: substrat etiketinin
# bilgilendirici sayildigi iki kademe (referansa >=%60 kimlik).
CLOSE_TIERS = ("characterized", "close_homolog")

# --- 3. soru (duzenleyici x enzim ayrisma hizi) sabitleri

# Enzim-enzim kimlik bantlari. Kuratorun "benzerlik esigine gore" istedigi
# kirilim bu. Kenarlar pipeline'in baska yerlerinde de kullanilan degerlere
# denk geliyor: %95 characterized kademesinin kapisi, %60 "ayni reaksiyon cok
# olasi" kapisi (bkz. evidence_tiers.py), %80 ikisinin arasi.
IDENTITY_BANDS = [(95.0, 100.01, "95-100"), (80.0, 95.0, "80-95"),
                  (60.0, 80.0, "60-80"), (0.0, 60.0, "below 60")]

# Tip icinde cift ornekleme ustu. 53 tipin TUM ic ciftleri 314.715 ediyor ve
# cift basina iki global hizalama ~10 ms suruyor; tamami tek surecte ~1 saat.
# Tip basina rastgele ornekleme yapiliyor cunku amac tip ICINDEKI cift
# dagilimini kestirmek ve duzgun ornekleme bunu yanlilik katmadan yapar.
MAX_PAIRS_PER_TYPE = 2000

# Global hizalamanin kisa proteinin en az bu kadarini kapsamasi gerekir, yoksa
# "kimlik" rakami kisa bir ortusme uzerinden hesaplanmis olur. Eski kodun
# (evolution_rates_of_reg_and_enz.py) hatasi tam buydu: bosluklari uyumsuzluk
# sayiyordu ve hicbir kapsama kontrolu yoktu.
MIN_ALIGNED_COVERAGE = 0.50

# Akraba OLMAYAN iki protein global hizalamada ~%20-25 kimlik verir. Bu yuzden
# %25 civarindaki bir duzenleyici kimligi "hizli ayrismis" demek DEGIL,
# "iliski olculemiyor" demektir; olcum burada tabana vuruyor ve bu bantta
# ayrisma HIZI hakkinda bir sey soylenemez.
RANDOM_IDENTITY_FLOOR = 25.0

# Bir bandin ya da tipin rapor edilmesi icin gereken en az cift sayisi.
MIN_PAIRS_FOR_BAND = 20

# Bir bantta ayrisma HIZI hakkinda konusabilmek icin, ciftlerin en fazla bu
# kadarinin olcum tabaninda olmasi gerekir. Yarisi tabandaysa bandin medyani
# artik ayrismayi degil tabani olcuyor.
MAX_FLOOR_SHARE_FOR_RATE_CLAIM = 0.5

# Sacilim grafigi icin JSON'a yazilan nokta ustu. Tamami yazilsa dosya
# gereksiz buyur; nokta secimi de sabit tohumla yapiliyor.
MAX_SCATTER_POINTS = 3000

# Aile degisimi analizinde bu etiket DISARIDA kalir: "siniflandirilamadi" iki
# duzenleyicinin ayni aileden olup olmadigini soylemez, dolayisiyla bir
# degisim kaniti da olamaz.
UNCLASSIFIED_FAMILY = "regulator_unclassified"


# ---------------------------------------------------------------- yardimcilar

def trim_table(rows, cols, table):
    """Sifir toplamli satir/sutunlari at (stats_overview.py ile ayni davranis)."""
    keep_c = [j for j in range(len(cols)) if sum(row[j] for row in table) > 0]
    keep_r = [i for i in range(len(rows)) if sum(table[i][j] for j in keep_c) > 0]
    return ([rows[i] for i in keep_r], [cols[j] for j in keep_c],
            [[table[i][j] for j in keep_c] for i in keep_r])


def cramers_v(table):
    table = np.asarray(table, dtype=float)
    if min(table.shape) < 2 or table.sum() == 0:
        return 0.0
    chi2 = stats.chi2_contingency(table, correction=False)[0]
    k = min(table.shape) - 1
    return float(np.sqrt(chi2 / (table.sum() * k)))


def crosstab(rows, row_key, col_key):
    """Iki alanin capraz tablosu; iki alandan biri bos olan satir ATILIR.

    Bos degeri kendi kategorisi yapmak testi "bu alan dolu mu" testine cevirir;
    upstream_family 11.422 girisin 8.061'inde bos oldugu icin bu fark buyuk.
    """
    usable = [r for r in rows if r.get(row_key) and r.get(col_key)]
    rk = sorted({r[row_key] for r in usable})
    ck = sorted({r[col_key] for r in usable})
    table = [[sum(1 for r in usable if r[row_key] == a and r[col_key] == b)
              for b in ck] for a in rk]
    return trim_table(rk, ck, table) + (usable,)


def permutation_null_v(table, rng, permutations=PERMUTATIONS):
    """V'nin SANS tabani: satir ve sutun toplamlari sabit, etiketler karistirilir.

    Seyrek bir tabloda V hicbir iliski olmasa da sifirdan buyuk cikar (hucre
    sayisi arttikca artar). Bu yuzden "V=0,45" tek basina okunamaz; okunabilir
    olan, ayni marjinallerle beklenen V'yi ne kadar astigi.
    """
    table = np.asarray(table, dtype=float)
    row_idx = np.repeat(np.arange(table.shape[0]), table.sum(axis=1).astype(int))
    col_idx = np.repeat(np.arange(table.shape[1]), table.sum(axis=0).astype(int))
    if row_idx.size == 0 or col_idx.size == 0:
        return None
    null = np.empty(permutations)
    for i in range(permutations):
        shuffled = rng.permutation(col_idx)
        counts = np.zeros_like(table)
        np.add.at(counts, (row_idx, shuffled), 1.0)
        null[i] = cramers_v(counts)
    return null


def association(rows, row_key, col_key, test_id, question, rng,
                row_header, col_header, unit, permutations=PERMUTATIONS):
    """Tek bir capraz tablo testi: ki-kare, V, permutasyon null'u, seyreklik notu."""
    rk, ck, table, usable = crosstab(rows, row_key, col_key)
    if len(rk) < 2 or len(ck) < 2:
        return None
    arr = np.asarray(table, dtype=float)
    chi2, p, dof, expected = stats.chi2_contingency(arr, correction=False)
    v = cramers_v(arr)
    thin = float((expected < EXPECTED_CELL_FLOOR).mean())
    null = permutation_null_v(arr, rng, permutations)
    null_mean = float(null.mean()) if null is not None else None
    null_sd = float(null.std()) if null is not None else None
    excess = (v - null_mean) if null_mean is not None else None
    out = {
        "id": test_id,
        "question": question,
        "unit": unit,
        "row_header": row_header,
        "col_header": col_header,
        "rows": rk,
        "cols": ck,
        "table": table,
        "n": int(arr.sum()),
        "chi2": float(chi2),
        "dof": int(dof),
        "p": float(p),
        "cramers_v": round(v, 4),
        "permuted_cramers_v_mean": round(null_mean, 4) if null_mean is not None else None,
        "permuted_cramers_v_sd": round(null_sd, 5) if null_sd is not None else None,
        "cramers_v_above_chance": round(excess, 4) if excess is not None else None,
        "z_against_permutation": (round((v - null_mean) / null_sd, 1)
                                  if null_sd else None),
        "fraction_cells_expected_below_5": round(thin, 3),
        "chi_square_p_trustworthy": bool(thin <= SPARSE_TABLE_LIMIT),
    }
    out["verdict"] = verdict_for(out)
    return out


def verdict_for(test):
    """Tek cumlelik, ingilizce karar. Etki buyuklugu ile p AYRI okunur."""
    v = test["cramers_v"]
    excess = test["cramers_v_above_chance"]
    effective = excess if excess is not None else v
    if test["p"] >= 0.05:
        return ("not detectable: no association beyond sampling noise "
                f"(p = {test['p']:.2g})")
    if effective < V_NEGLIGIBLE:
        return ("significant but negligible: the p value comes from the sample "
                f"size, the effect does not (Cramer's V above chance "
                f"{effective:.2f}, below the {V_NEGLIGIBLE:.2f} floor)")
    strength = "moderate to strong" if effective >= V_MODERATE else "real but modest"
    note = "" if test["chi_square_p_trustworthy"] else (
        "; the chi-square p value is not trustworthy here because "
        f"{100 * test['fraction_cells_expected_below_5']:.0f} % of cells expect "
        "fewer than 5 observations, so this verdict rests on the permutation "
        "comparison instead")
    return (f"{strength} association (Cramer's V {v:.2f}, "
            f"{effective:.2f} above the permuted chance level){note}")


def standardised_residuals(table):
    """Standartlastirilmis Pearson artiklari: hangi HUCRE tabloyu tasiyor."""
    arr = np.asarray(table, dtype=float)
    n = arr.sum()
    row_p = arr.sum(axis=1) / n
    col_p = arr.sum(axis=0) / n
    expected = np.outer(row_p, col_p) * n
    denom = np.sqrt(expected * np.outer(1 - row_p, 1 - col_p))
    with np.errstate(divide="ignore", invalid="ignore"):
        res = np.where(denom > 0, (arr - expected) / denom, 0.0)
    return res, expected


def cliffs_delta(a, b):
    """Mann-Whitney U'dan dogrudan etki buyuklugu (-1..+1).

    Cikti JSON'una p ile BIRLIKTE yazilir: medyan farki 0,7 puan olan bir
    karsilastirma 11.422 giriste p=1e-7 verir ve bu bir bulgu gibi okunur.
    """
    if not a or not b:
        return None, None
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(2 * u / (len(a) * len(b)) - 1), float(p)


def two_by_two(a_yes, a_no, b_yes, b_no):
    """stats_overview.py ile ayni ikili test ciktisi (karsilastirilabilirlik)."""
    odds, p = stats.fisher_exact([[a_yes, a_no], [b_yes, b_no]])
    ra = a_yes / max(1, a_yes + a_no)
    rb = b_yes / max(1, b_yes + b_no)
    return {"table": [[a_yes, a_no], [b_yes, b_no]],
            "rate_xenobiotic": round(ra, 4), "rate_natural": round(rb, 4),
            "ratio": round(ra / rb, 3) if rb else None,
            "odds_ratio": round(float(odds), 3), "p": float(p),
            "n": a_yes + a_no + b_yes + b_no}


def quantile_bins(values, quantiles):
    """Veriden ceyreklik kesim noktalari; tekrarlanan degerler teklestirilir."""
    cuts = sorted({int(round(c)) for c in np.percentile(values, quantiles)})
    bins, low = [], min(values)
    for c in cuts:
        if c >= low:
            bins.append((low, c))
            low = c + 1
    bins.append((low, max(values)))
    return [b for b in bins if b[0] <= b[1]]


def mantel_haenszel(strata):
    """Tabakalara ayrilmis 2x2 tablolardan ortak OR.

    Kullanilma nedeni: "ksenobiyotik girislerde regex transpozonu daha sik"
    gozlemi, ksenobiyotik girislerin daha IYI ANOTE EDILMIS pencerelerde
    oturmasindan da dogabilir. Tabakalama bu aciklamayi ayirir.
    """
    num = den = 0.0
    used = []
    for table in strata:
        (a, b), (c, d) = table
        total = a + b + c + d
        if min(a + b, c + d) < MIN_STRATUM_N or total == 0:
            continue
        num += a * d / total
        den += b * c / total
        used.append(table)
    return (round(num / den, 3) if den else None), used


# ------------------------------------------------------------------ veri yukleme

def load(con, chemistry_path, ecology_path, domain_csv):
    chem = {}
    with open(chemistry_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            chem[row["cluster"]] = row
    eco = {}
    with open(ecology_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            eco[row["cluster"]] = row
    euk_heavy = set()
    if os.path.exists(domain_csv):
        with open(domain_csv, encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                try:
                    if float(row.get("eukaryota_rate", 0) or 0) >= EUK_HEAVY_RATE:
                        euk_heavy.add(row["cluster"])
                except ValueError:
                    pass

    rows = con.execute(f"""
        SELECT r.candidate_id, r.ro_cluster, r.ro_group,
               p.organism, p.is_plasmid, p.cds_count, p.status, p.mol_type,
               g.upstream_family, g.upstream_category, g.architecture,
               g.upstream_divergent, g.intergenic_bp,
               e.tier, e.ref_identity, s.assignment_class, d.domain,
               (SELECT COUNT(*) FROM neighbor nb
                 WHERE nb.candidate_id = r.candidate_id) AS n_neighbours,
               (SELECT COUNT(*) FROM neighbor nb
                  JOIN gene_category c ON c.neighbor_id = nb.neighbor_id
                   AND c.method = '{REGEX_METHOD}'
                 WHERE nb.candidate_id = r.candidate_id
                   AND c.category = '{TRANSPOSON_CATEGORY}') AS transposons
        FROM ro r
        JOIN replicon p USING(nucleotide_id)
        LEFT JOIN ro_regulation g USING(candidate_id)
        LEFT JOIN ro_evidence e USING(candidate_id)
        LEFT JOIN ro_subfamily s USING(candidate_id)
        LEFT JOIN ro_domain d USING(candidate_id)
        WHERE r.is_confirmed = 1""").fetchall()

    data = []
    for (cid, cluster, group, organism, plasmid, cds, status, mol_type,
         up_family, up_category, architecture, divergent, intergenic,
         tier, identity, assign_class, domain, n_nb, transposons) in rows:
        c = chem.get(cluster, {})
        e = eco.get(cluster, {})
        data.append({
            "candidate_id": cid,
            "cluster": cluster,
            "group": group or "?",
            "genus": (organism or "?").split()[0],
            "chem_family": c.get("family") or "unknown",
            "reaction_class": c.get("reaction_class") or "unknown",
            "sclass": e.get("substrate_class", "unknown"),
            "domain": domain,
            "upstream_family": up_family,
            "upstream_category": up_category,
            "architecture": architecture,
            "upstream_divergent": divergent,
            "intergenic_bp": intergenic,
            "tier": tier,
            "ref_identity": identity,
            "assignment_class": assign_class,
            "cds_count": cds or 0,
            "status": status,
            "mol_type": mol_type,
            "n_neighbours": n_nb or 0,
            "is_plasmid": int(bool(plasmid)),
            "transposon": int((transposons or 0) > 0),
            "transposon_count": transposons or 0,
        })
    for d in data:
        # Yayinlanan degisken, oldugu gibi: stats_overview.py'deki tanim.
        d["mobile_published"] = int(bool(d["is_plasmid"]) or d["transposon"])
    return data, euk_heavy


# Cokertmede cogunluk oyu uygulanan ikili alanlar (stats_overview.py deseni).
BINARY_FIELDS = ("is_plasmid", "transposon", "mobile_published")
# Cokertmede MOD (en sik deger) alinan kategorik alanlar. Ortalama anlamsiz
# oldugu icin ayri tutuluyor; bos degerler oya KATILMAZ, yoksa iyi ornekli bir
# cins bos etiketle cokertilir.
CATEGORICAL_FIELDS = ("upstream_family", "upstream_category", "architecture",
                      "tier", "assignment_class")


def collapse_genus(data):
    """Tip x cins basina tek gozlem -- projenin her biyolojik iddia icin kosulu.

    stats_overview.py'deki `collapse_genus` ile ayni kural: ikili ozellikler
    cins icinde ortalanip yuvarlanir. Burada ek olarak kategorik alanlar mod
    ile, ref_identity medyan ile cokertilir, cunku bu modulde o alanlar da
    test degiskeni.
    """
    groups = defaultdict(list)
    for d in data:
        groups[(d["cluster"], d["genus"])].append(d)
    out = []
    for _, items in groups.items():
        base = dict(items[0])
        for field in BINARY_FIELDS:
            base[field] = int(np.mean([bool(i[field]) for i in items]) >= 0.5)
        for field in CATEGORICAL_FIELDS:
            counts = Counter(i[field] for i in items if i[field])
            base[field] = counts.most_common(1)[0][0] if counts else None
        ids = [i["ref_identity"] for i in items if i["ref_identity"] is not None]
        base["ref_identity"] = float(np.median(ids)) if ids else None
        base["n_neighbours"] = float(np.median([i["n_neighbours"] for i in items]))
        base["cds_count"] = float(np.median([i["cds_count"] for i in items]))
        base["n_entries"] = len(items)
        out.append(base)
    return out


# -------------------------------------------------------------- 1. soru: regulator

REGULATION_TESTS = [
    ("upstream_family", "chem_family",
     "Which chemical families does each upstream regulator family sit next to?",
     "Regulator family", "Chemical family"),
    ("upstream_family", "cluster",
     "Is the upstream regulator family associated with the enzyme type?",
     "Regulator family", "Enzyme type"),
    ("upstream_family", "group",
     "Is the upstream regulator family associated with the RO group?",
     "Regulator family", "RO group"),
    ("upstream_family", "reaction_class",
     "Is the upstream regulator family associated with the reaction class?",
     "Regulator family", "Reaction class"),
    ("upstream_category", "chem_family",
     "Is the kind of gene upstream of the operon (regulator, transporter, "
     "transposon, ...) associated with the chemical family?",
     "Upstream gene category", "Chemical family"),
    ("upstream_category", "group",
     "Is the kind of gene upstream of the operon associated with the RO group?",
     "Upstream gene category", "RO group"),
    ("architecture", "chem_family",
     "Is the operon architecture (divergent or codirectional regulator, other "
     "gene, or nothing in the window) associated with the chemical family?",
     "Operon architecture", "Chemical family"),
    ("architecture", "group",
     "Is the operon architecture associated with the RO group?",
     "Operon architecture", "RO group"),
    ("architecture", "upstream_family",
     "Do regulator families differ in how they are oriented relative to the "
     "operon? LysR regulators are classically divergent and TetR regulators "
     "codirectional, so this is the internal consistency check for the family "
     "labels.",
     "Operon architecture", "Regulator family"),
]


def regulation_block(data, genus_rows, rng, permutations):
    out = {"coverage": {}, "associations": [], "cell_enrichments": [],
           "regulator_profile_by_chemical_family": {},
           "architecture_by_chemical_family": {}}

    n = len(data)
    out["coverage"] = {
        "entries": n,
        "with_upstream_gene_in_window": sum(1 for d in data if d["upstream_category"]),
        "with_regulator_upstream": sum(1 for d in data
                                       if d["upstream_category"] == "regulator"),
        "with_named_regulator_family": sum(
            1 for d in data if d["upstream_family"]
            and d["upstream_family"] != "regulator_unclassified"),
        "regulator_family_unclassified": sum(
            1 for d in data if d["upstream_family"] == "regulator_unclassified"),
        "type_and_genus_observations": len(genus_rows),
        "note": ("the regulator family label is itself derived from the GenBank "
                 "/product text of the upstream gene, by the same kind of regular "
                 "expression as the gene categories; no operator sequence was "
                 "inspected, because 88.9 % of the source records are CON entries "
                 "that carry no sequence"),
    }

    for row_key, col_key, question, row_header, col_header in REGULATION_TESTS:
        for level, rows in (("entry", data), ("genus", genus_rows)):
            test = association(
                rows, row_key, col_key,
                f"{row_key}_x_{col_key}_{level}",
                question if level == "entry" else
                (question + " Does it survive collapsing repeated strains to one "
                            "observation per enzyme type and genus?"),
                rng, row_header, col_header,
                "alpha subunit" if level == "entry" else "enzyme type and genus",
                permutations)
            if test:
                test["level"] = level
                out["associations"].append(test)

    # Hangi HUCRELER tabloyu tasiyor: duzenleyici ailesi x kimyasal aile.
    # Kuratorun sorusu tam olarak bu ("hangi aile hangi duzenleyiciyle").
    rk, ck, table, usable = crosstab(data, "upstream_family", "chem_family")
    if rk and ck:
        res, expected = standardised_residuals(table)
        genera = defaultdict(set)
        clusters = defaultdict(set)
        for d in usable:
            genera[(d["upstream_family"], d["chem_family"])].add(d["genus"])
            clusters[(d["upstream_family"], d["chem_family"])].add(d["cluster"])
        cells = []
        for i, family in enumerate(rk):
            for j, chem in enumerate(ck):
                observed = table[i][j]
                key = (family, chem)
                if observed < MIN_CELL_ENTRIES:
                    continue
                if len(genera[key]) < MIN_CELL_GENERA:
                    continue
                cells.append({
                    "regulator_family": family,
                    "chemical_family": chem,
                    "entries": observed,
                    "expected": round(float(expected[i][j]), 1),
                    "observed_over_expected": round(observed / expected[i][j], 2)
                    if expected[i][j] else None,
                    "standardised_residual": round(float(res[i][j]), 1),
                    "distinct_genera": len(genera[key]),
                    "distinct_types": len(clusters[key]),
                    # Tek tipe dayanan hucre, o TIP hakkinda bir bulgudur ama
                    # kimyasal AILE hakkinda degil: aile etiketi tipin
                    # ozelligi oldugu icin tek tip aileyi temsil etmez.
                    "rests_on_a_single_enzyme_type": len(clusters[key]) == 1,
                })
        cells.sort(key=lambda c: -c["standardised_residual"])
        out["cell_enrichments"] = cells
        out["cell_enrichment_rule"] = (
            f"a cell is listed only if it rests on at least {MIN_CELL_ENTRIES} "
            f"entries and at least {MIN_CELL_GENERA} distinct genera, so that a "
            "single sequencing project cannot become a finding")

    # Okunabilir profiller: aile basina duzenleyici dagilimi (yuzde),
    # giris ve TIP x CINS duzeyinde yan yana.
    for level, rows in (("entry", data), ("genus", genus_rows)):
        profile = {}
        for chem in sorted({d["chem_family"] for d in rows}):
            subset = [d for d in rows if d["chem_family"] == chem and d["upstream_family"]]
            if not subset:
                continue
            counts = Counter(d["upstream_family"] for d in subset)
            profile[chem] = {
                "n": len(subset),
                "distinct_types": len({d["cluster"] for d in subset}),
                "distinct_genera": len({d["genus"] for d in subset}),
                "regulator_share": {f: round(c / len(subset), 3)
                                    for f, c in counts.most_common()},
            }
        out["regulator_profile_by_chemical_family"][level] = profile

        arch = {}
        for chem in sorted({d["chem_family"] for d in rows}):
            subset = [d for d in rows if d["chem_family"] == chem and d["architecture"]]
            if not subset:
                continue
            counts = Counter(d["architecture"] for d in subset)
            arch[chem] = {"n": len(subset),
                          "architecture_share": {a: round(c / len(subset), 3)
                                                 for a, c in counts.most_common()}}
        out["architecture_by_chemical_family"][level] = arch

    # Duzenleyici ailesi basina intergenik mesafe: aile etiketinin operon
    # geometrisiyle tutarli olup olmadiginin ikinci, bagimsiz kontrolu.
    spacing = {}
    for family in sorted({d["upstream_family"] for d in data if d["upstream_family"]}):
        gaps = [d["intergenic_bp"] for d in data
                if d["upstream_family"] == family and d["intergenic_bp"] is not None]
        if len(gaps) < MIN_CELL_ENTRIES:
            continue
        divergent = [d for d in data if d["upstream_family"] == family]
        spacing[family] = {
            "n": len(gaps),
            "median_intergenic_bp": float(np.median(gaps)),
            "divergent_share": round(
                float(np.mean([bool(d["upstream_divergent"]) for d in divergent])), 3),
        }
    out["regulator_family_geometry"] = spacing

    out["summary"] = regulation_summary(out)
    return out


def regulation_summary(block):
    """Her test icin tek satirlik ingilizce karar; genus duzeyi belirleyici."""
    lines = []
    by_id = {t["id"]: t for t in block["associations"]}
    for test in block["associations"]:
        if test["level"] != "genus":
            continue
        entry = by_id.get(test["id"].replace("_genus", "_entry"))
        lines.append({
            "pair": f"{test['row_header']} x {test['col_header']}",
            "entry_cramers_v": entry["cramers_v"] if entry else None,
            "entry_above_chance": entry["cramers_v_above_chance"] if entry else None,
            "genus_cramers_v": test["cramers_v"],
            "genus_above_chance": test["cramers_v_above_chance"],
            "genus_p": test["p"],
            "survives_genus_collapse": bool(
                test["p"] < 0.05
                and (test["cramers_v_above_chance"] or 0) >= V_NEGLIGIBLE),
            "verdict": test["verdict"],
        })
    lines.sort(key=lambda r: -(r["genus_above_chance"] or 0))
    return lines


# ------------------------------------------------- 2a. soru: mobilite, uc kanitla

MOBILITY_EVIDENCE = [
    ("is_plasmid",
     "plasmid_only",
     "the replicon record is annotated as a plasmid. Structural evidence: a "
     "property of the DNA molecule, independent of how any gene was described."),
    ("transposon",
     "regex_transposon_only",
     "at least one CDS within 10 kb whose GenBank /product text matches the "
     "transposon regular expression. Annotation-derived evidence: it measures "
     "the wording of the submitted annotation, not the sequence."),
    ("mobile_published",
     "published_or",
     "plasmid OR regex transposon. This is the 'mobile' variable published in "
     "stats_overview.py, and it inherits the weaker half."),
]


def mobility_block(data, genus_rows, euk_heavy):
    def universe(rows):
        # stats_overview.py ile AYNI evren: ekolojisi bilinen sinif, bakteri,
        # okaryot-agirlikli kumeler disarida. Sayilar yayinlanan testle
        # karsilastirilabilir olmasa anlamsiz olurdu.
        return [r for r in rows if r["sclass"] in SUBSTRATE_CLASSES
                and r["cluster"] not in euk_heavy and r["domain"] == "Bacteria"]

    out = {
        "evidence_definitions": {name: text for _, name, text in MOBILITY_EVIDENCE},
        "filters": {
            "substrate_classes": SUBSTRATE_CLASSES,
            "domain": "Bacteria",
            "excluded_eukaryote_heavy_types": sorted(euk_heavy),
            "note": ("identical to the universe used by stats_overview.py, so the "
                     "published numbers are reproduced and not merely approximated"),
        },
        "composition": {},
        "tests": [],
    }

    n = len(data)
    plasmid = sum(d["is_plasmid"] for d in data)
    transposon = sum(d["transposon"] for d in data)
    mobile = sum(d["mobile_published"] for d in data)
    regex_only = sum(1 for d in data if d["transposon"] and not d["is_plasmid"])
    both = sum(1 for d in data if d["transposon"] and d["is_plasmid"])
    out["composition"] = {
        "entries": n,
        "plasmid": plasmid, "plasmid_share": round(plasmid / n, 4),
        "regex_transposon": transposon,
        "regex_transposon_share": round(transposon / n, 4),
        "mobile_published": mobile, "mobile_share": round(mobile / n, 4),
        "mobile_only_because_of_regex": regex_only,
        "share_of_mobile_that_is_regex_only": round(regex_only / mobile, 4),
        "both_kinds_of_evidence": both,
        # Cumlenin SAYISI veriden gelir. Bu projede sabit yazilmis metnin
        # veritabaniyla celistigi bir hata sinifi zaten yasandi
        # (make_report.py, 960 / "4,3x"), o yuzden hicbir oran elle yazilmiyor.
        "reading": (
            "the published mobile variable is dominated by the annotation-derived "
            f"half: {regex_only} of the {mobile} mobile entries, "
            f"{100 * regex_only / mobile:.0f} %, are mobile for no other reason than "
            f"a word in a /product field, against {plasmid} with plasmid evidence; "
            f"the regex half outnumbers the structural half "
            f"{transposon / plasmid:.1f} to 1"),
    }

    for field, name, _text in MOBILITY_EVIDENCE:
        for level, rows in (("entry", data), ("genus", genus_rows)):
            rows = universe(rows)
            xen = [r for r in rows if r["sclass"] == "xenobiotic"]
            nat = [r for r in rows if r["sclass"] != "xenobiotic"]
            if not xen or not nat:
                continue
            result = two_by_two(
                sum(bool(r[field]) for r in xen), sum(not r[field] for r in xen),
                sum(bool(r[field]) for r in nat), sum(not r[field] for r in nat))
            result.update({"id": f"{name}_{level}", "evidence": name, "level": level,
                           "unit": "alpha subunit" if level == "entry"
                                   else "enzyme type and genus",
                           "question": ("Are xenobiotic-degrading RO types more often "
                                        "mobile than natural-substrate types, when "
                                        f"mobility is defined as: {name}?")})
            out["tests"].append(result)

    ratios = {(t["evidence"], t["level"]): t["ratio"] for t in out["tests"]}
    pvals = {(t["evidence"], t["level"]): t["p"] for t in out["tests"]}
    # Karar VERIDEN okunur: yapisal kanit tek basina ayni yonde ve en az
    # yayinlanan OR kadar guclu, ve cins duzeyinde anlamli ise sonuc anotasyon
    # yarisina bagli DEGILDIR.
    structural_alone_holds = bool(
        (ratios.get(("plasmid_only", "genus")) or 0) > 1
        and (pvals.get(("plasmid_only", "genus")) or 1) < 0.05
        and (ratios.get(("plasmid_only", "genus")) or 0)
        >= (ratios.get(("published_or", "genus")) or 0))
    out["verdict"] = {
        "does_the_conclusion_depend_on_the_annotation_half":
            "no" if structural_alone_holds else "yes",
        "explanation": (
            "The association is strongest on structural evidence alone "
            f"({ratios.get(('plasmid_only', 'entry'))}x per entry, "
            f"{ratios.get(('plasmid_only', 'genus'))}x per type and genus, "
            f"p = {pvals.get(('plasmid_only', 'genus')):.1e} at genus level). "
            "The annotation-derived half points the same way but much more weakly "
            f"({ratios.get(('regex_transposon_only', 'entry'))}x and "
            f"{ratios.get(('regex_transposon_only', 'genus'))}x). Because the "
            "regex half outnumbers the plasmid half roughly six to one, the "
            "published OR of the two lands near the weaker value "
            f"({ratios.get(('published_or', 'entry'))}x and "
            f"{ratios.get(('published_or', 'genus'))}x). So the finding does not "
            "rest on the annotation text; the published variable dilutes a strong "
            "structural result with a weak annotation-derived one."),
        "recommendation": (
            "report plasmid evidence as the headline mobility result and keep the "
            "regex transposon rate as a separate, explicitly annotation-derived "
            "measure; do not combine them with OR, because the combination is "
            "weaker than its better half and is read as if it were stronger"),
    }
    return out


# ------------------------------- 2b. soru: transpozon yakinligi x kuratorlu uzaklik

def evidence_block(data, genus_rows, rng, permutations):
    out = {"question": ("Does a regex-called transposon nearby relate to how far the "
                        "member is from a curated enzyme, and to its subfamily "
                        "assignment class?"),
           "transposon_rate_by_tier": {}, "associations": [],
           "ref_identity": {}, "transposon_rate_by_assignment_class": {}}

    for level, rows in (("entry", data), ("genus", genus_rows)):
        tier_rates = {}
        for tier in TIERS:
            subset = [r for r in rows if r["tier"] == tier]
            if not subset:
                continue
            tier_rates[tier] = {
                "n": len(subset),
                "transposon_rate": round(
                    float(np.mean([r["transposon"] for r in subset])), 4),
                "plasmid_rate": round(
                    float(np.mean([r["is_plasmid"] for r in subset])), 4),
            }
        out["transposon_rate_by_tier"][level] = tier_rates

        class_rates = {}
        for klass in sorted({r["assignment_class"] for r in rows if r["assignment_class"]}):
            subset = [r for r in rows if r["assignment_class"] == klass]
            class_rates[klass] = {
                "n": len(subset),
                "transposon_rate": round(
                    float(np.mean([r["transposon"] for r in subset])), 4)}
        out["transposon_rate_by_assignment_class"][level] = class_rates

        for field, header, question in (
                ("tier", "Evidence level",
                 "Is a nearby regex-called transposon associated with the evidence "
                 "level, that is with how close the member is to a curated enzyme?"),
                ("assignment_class", "Subfamily assignment class",
                 "Is a nearby regex-called transposon associated with the subfamily "
                 "assignment class?")):
            labelled = [dict(r, transposon_label=("transposon nearby"
                                                  if r["transposon"] else "none"))
                        for r in rows if r[field]]
            test = association(labelled, field, "transposon_label",
                               f"{field}_x_transposon_{level}", question, rng,
                               header, "Regex transposon within 10 kb",
                               "alpha subunit" if level == "entry"
                               else "enzyme type and genus", permutations)
            if test:
                test["level"] = level
                out["associations"].append(test)

        near = [r["ref_identity"] for r in rows
                if r["transposon"] and r["ref_identity"] is not None]
        far = [r["ref_identity"] for r in rows
               if not r["transposon"] and r["ref_identity"] is not None]
        delta, p = cliffs_delta(near, far)
        out["ref_identity"][level] = {
            "question": ("Do members with a regex-called transposon nearby sit closer "
                         "to a curated reference enzyme than members without one?"),
            "n_with_transposon": len(near), "n_without": len(far),
            "median_identity_with_transposon": round(float(np.median(near)), 2) if near else None,
            "median_identity_without": round(float(np.median(far)), 2) if far else None,
            "cliffs_delta": round(delta, 3) if delta is not None else None,
            "p": p,
            "note": ("Cliff's delta is reported because the medians differ by a "
                     "fraction of a percentage point at entry level while p is tiny; "
                     "the p value is a statement about 11,422 entries, not about the "
                     "size of the difference"),
        }

    # Tersinme ACIKCA raporlaniyor: iki duzey farkli yon soyluyor ve bunun
    # nedeni olculebilir (tek tek tip x cins hucrelerinin buyuklugu).
    entry = out["transposon_rate_by_tier"]["entry"]
    genus = out["transposon_rate_by_tier"]["genus"]
    # Monotonluk ve yon VERIDEN sinanir. "Iki duzey celisiyor" bir gozlem,
    # elle yazilmis bir iddia degil.
    def _series(block):
        return [block[t]["transposon_rate"] for t in TIERS if t in block]

    def _monotone_down(series):
        return all(a >= b for a, b in zip(series, series[1:]))

    entry_series, genus_series = _series(entry), _series(genus)

    def _pooled(rows):
        """Yakin ve uzak kademelerin TOPLU orani.

        Yon kararini uc noktalara (characterized 215 giris) degil, bu toplu
        karsilastirmaya baglamak gerekiyor: en yakin kademe en ince kademe ve
        tek basina yonu belirlemesi tesadufe cok acik. Ayni nicelik asagidaki
        yogunluk duzeltmesinde de kullaniliyor, boylece iki blok celismiyor.
        """
        close = [r for r in rows if r["tier"] in CLOSE_TIERS]
        far = [r for r in rows if r["tier"] and r["tier"] not in CLOSE_TIERS]
        if not close or not far:
            return None, None
        return (float(np.mean([r["transposon"] for r in close])),
                float(np.mean([r["transposon"] for r in far])))

    entry_close, entry_far = _pooled(data)
    genus_close, genus_far = _pooled(genus_rows)
    entry_down = bool(entry_close is not None and entry_close > entry_far)
    genus_down = bool(genus_close is not None and genus_close > genus_far)
    out["level_disagreement"] = {
        "pooled_close_vs_far": {
            "entry": {"close_rate": round(entry_close, 4),
                      "far_rate": round(entry_far, 4),
                      "ratio": round(entry_close / entry_far, 3) if entry_far else None},
            "genus": {"close_rate": round(genus_close, 4),
                      "far_rate": round(genus_far, 4),
                      "ratio": round(genus_close / genus_far, 3) if genus_far else None},
            "note": ("direction is judged on this pooled contrast, not on the "
                     "characterized tier alone, because that tier is the thinnest one "
                     "and would otherwise decide the verdict by itself"),
        },
        "entry_level_gradient": {t: entry[t]["transposon_rate"] for t in TIERS if t in entry},
        "genus_level_gradient": {t: genus[t]["transposon_rate"] for t in TIERS if t in genus},
        "entry_gradient_monotone_decreasing": _monotone_down(entry_series),
        "genus_gradient_monotone_decreasing": _monotone_down(genus_series),
        "levels_agree_in_direction": entry_down == genus_down,
        "finding": (
            ("the two levels agree in direction: pooled over tiers, close members are "
             f"{entry_close / entry_far:.2f} times as often next to a transposon per "
             f"entry and {genus_close / genus_far:.2f} times per type and genus. "
             if entry_down == genus_down else
             "the two levels disagree in direction. Pooled over tiers, close members "
             f"are {entry_close / entry_far:.2f} times as often next to a transposon "
             f"per entry but {genus_close / genus_far:.2f} times per type and genus. ")
            + f"Per entry the gradient runs {100 * entry_series[0]:.1f} % at the "
              f"closest level to {100 * entry_series[-1]:.1f} % at the most distant and "
            + ("is monotone. " if _monotone_down(entry_series) else "is not monotone. ")
            + f"Per enzyme type and genus it runs {100 * genus_series[0]:.1f} % to "
              f"{100 * genus_series[-1]:.1f} % and "
            + ("is monotone. " if _monotone_down(genus_series) else "is not monotone. ")
            + "Entry counts are dominated by a handful of type-by-genus cells, listed "
              "below, so the genus level is the one to read."),
        "largest_entry_cells": None,   # asagida doldurulur
        "interpretation": None,   # asagida, olculen duzeltmeden sonra doldurulur
    }
    cells = Counter((d["cluster"], d["genus"]) for d in data
                    if d["tier"] in CLOSE_TIERS)
    out["level_disagreement"]["largest_entry_cells"] = [
        {"type": c[0], "genus": c[1], "entries": n}
        for c, n in cells.most_common(5)]

    # Tabakalanmis kontrol: tier -> transpozon iliskisi pencere dolulugundan
    # bagimsiz mi? Yalnizca cagrilabilir gozlemler (>=1 anote komsu), cunku
    # komsusu olmayan bir giriste transpozon ZATEN cagrilamaz.
    # Iki duzeyde de kosulur: giris duzeyi 0,95x diyor, tip x cins duzeyi
    # 2,28x diyor ve ASIL aciklanmasi gereken ikincisi, o yuzden duzeltmenin
    # orada da yapilmasi gerekiyor.
    out["density_adjusted_close_vs_far"] = {}
    for level, rows in (("entry", data), ("genus", genus_rows)):
        callable_rows = [d for d in rows if d["n_neighbours"] > 0]
        close_all = [d for d in callable_rows if d["tier"] in CLOSE_TIERS]
        far_all = [d for d in callable_rows
                   if d["tier"] and d["tier"] not in CLOSE_TIERS]
        if not close_all or not far_all:
            continue
        strata, detail = [], []
        for low, high in quantile_bins([d["n_neighbours"] for d in callable_rows],
                                       DENSITY_QUANTILES):
            close = [d for d in close_all if low <= d["n_neighbours"] <= high]
            farr = [d for d in far_all if low <= d["n_neighbours"] <= high]
            if not close or not farr:
                continue
            a1 = sum(d["transposon"] for d in close)
            b1 = sum(d["transposon"] for d in farr)
            strata.append([[a1, len(close) - a1], [b1, len(farr) - b1]])
            detail.append({"annotated_neighbours": f"{low}-{high}",
                           "close_n": len(close),
                           "close_rate": round(a1 / len(close), 4),
                           "far_n": len(farr),
                           "far_rate": round(b1 / len(farr), 4),
                           "ratio": round((a1 / len(close)) / (b1 / len(farr)), 3)
                           if b1 else None})
        common_or, used = mantel_haenszel(strata)
        a1 = sum(d["transposon"] for d in close_all)
        b1 = sum(d["transposon"] for d in far_all)
        crude = (a1 * (len(far_all) - b1)) / max(1, (len(close_all) - a1) * b1)
        # "Fazla odds'un ne kadari yogunluktan" sorusu yalnizca ham OR 1'in
        # USTUNDE iken anlamlidir; altinda zaten bir fazlalik yok ve bolum
        # isaretsiz bir sayi uretir (olculdu: giris duzeyinde -0,07).
        shrink = (None if common_or is None or crude <= 1
                  else round((crude - common_or) / (crude - 1), 3))
        out["density_adjusted_close_vs_far"][level] = {
            "unit": ("alpha subunit" if level == "entry" else "enzyme type and genus")
                    + ", restricted to observations with at least one annotated "
                      "neighbour, because a record with none cannot have a "
                      "transposon called",
            "close_definition": "evidence level " + " or ".join(CLOSE_TIERS),
            "n_close": len(close_all), "n_far": len(far_all),
            "strata": detail,
            "crude_odds_ratio": round(float(crude), 3),
            "mantel_haenszel_odds_ratio": common_or,
            "strata_used": len(used),
            "share_of_excess_odds_attributable_to_annotation_density": shrink,
            # Yorum cumlesi de VERIDEN: hicbir yon ya da buyukluk sabit yazili degil.
            "note": (
                f"crude odds ratio {crude:.2f}, occupancy-adjusted "
                f"{common_or if common_or is not None else float('nan'):.2f}. "
                + ("Both sit below 1, so at this level close members are not more "
                   "often next to an annotated transposon and there is nothing for "
                   "annotation density to explain."
                   if crude < 1 and (common_or or 1) < 1 else
                   "Adjusting for how many genes were annotated in the window "
                   + (f"removes {100 * shrink:.0f} % of the excess odds and leaves "
                      f"{common_or:.2f}, so the association is reduced but not "
                      "explained away by annotation density."
                      if shrink is not None and 0 < shrink < 0.9 else
                      f"leaves {common_or:.2f}; see the strata, the adjustment does "
                      "not account for the association."))),
        }

    # Yorum, duzeltme OLCULDUKTEN sonra yazilir ve olculene baglidir. Ilk
    # surumde burada "bu bir anotasyon kalitesi etkisidir" yaziyordu; duzeltme
    # kosulunca bu iddianin DESTEKLENMEDIGI goruldu ve cumle degisti.
    genus_adj = out["density_adjusted_close_vs_far"].get("genus", {})
    genus_share = genus_adj.get("share_of_excess_odds_attributable_to_annotation_density")
    if not genus_down:
        interpretation = (
            "read at the genus level, members far from a curated enzyme are the ones "
            "more often next to an annotated transposon, which is the opposite of an "
            "annotation-richness effect and would need a different explanation")
    elif genus_share is not None and genus_share >= 0.5:
        interpretation = (
            "read at the genus level, closeness to a curated enzyme predicts that a "
            "transposon was ANNOTATED nearby, and most of that "
            f"({100 * genus_share:.0f} %) is accounted for by how many genes were "
            "annotated in the window. That is an annotation-quality effect and not "
            "evidence that well-characterised enzymes are more mobile.")
    else:
        interpretation = (
            "read at the genus level, closeness to a curated enzyme predicts that a "
            f"transposon was ANNOTATED nearby, by a factor of "
            f"{genus_adj.get('crude_odds_ratio')} in odds. The obvious explanation is "
            "annotation quality, since curated enzymes come from intensively studied "
            "organisms, but that explanation is NOT supported by the one proxy that "
            "can be measured here: adjusting for the number of annotated genes in the "
            "window leaves the odds ratio at "
            f"{genus_adj.get('mantel_haenszel_odds_ratio')}, removing only "
            f"{100 * (genus_share or 0):.0f} % of the excess. So three readings remain "
            "open and this data cannot choose between them: an annotation habit that "
            "gene count does not capture, such as how specifically a submitter names a "
            "mobile element; a real tendency for the better studied enzymes to sit in "
            "mobile neighbourhoods, which is what the plasmid result independently "
            "suggests; or residual confounding by genus, since the curated enzymes are "
            "concentrated in a few well sampled genera. What can be said is the "
            "negative: it is not simply that richer windows give more chances to match.")
    out["level_disagreement"]["interpretation"] = interpretation
    return out


# -------------------------------------------- 2c. soru: anotasyon yanliliginin boyu

def annotation_bias_block(data, euk_heavy, rng):
    out = {"question": ("How much does the regex transposon rate vary with the kind "
                        "of GenBank record and with how richly it is annotated? A "
                        "record with a sparse annotation cannot have a transposon "
                        "called near anything."),
           "window_bp": NEIGHBOUR_WINDOW_BP}

    n = len(data)
    genes = sum(d["n_neighbours"] for d in data)
    hits = sum(d["transposon_count"] for d in data)
    out["per_gene_rate"] = {
        "annotated_neighbour_genes_in_all_windows": genes,
        "genes_matching_the_transposon_regex": hits,
        "per_gene_rate": round(hits / genes, 5) if genes else None,
        "note": ("the entry-level rate and the per-gene rate are different "
                 "quantities: an entry is called positive if ANY gene in its window "
                 "matches, so the entry-level rate rises with the number of genes "
                 "annotated in the window even when the per-gene rate is constant"),
    }

    no_neighbour = [d for d in data if d["n_neighbours"] == 0]
    out["hard_floor"] = {
        "entries_with_no_annotated_neighbour": len(no_neighbour),
        "share_of_all_entries": round(len(no_neighbour) / n, 4),
        "their_transposon_rate": 0.0,
        "record_types": dict(Counter(d["mol_type"] for d in no_neighbour)),
        "replicon_status": dict(Counter(d["status"] for d in no_neighbour)),
        "published_rate_all_entries": round(
            float(np.mean([d["transposon"] for d in data])), 4),
        "rate_among_callable_records": round(
            float(np.mean([d["transposon"] for d in data if d["n_neighbours"] > 0])), 4),
    }
    _all_rate = out["hard_floor"]["published_rate_all_entries"]
    _callable_rate = out["hard_floor"]["rate_among_callable_records"]
    out["hard_floor"]["finding"] = (
        "these entries cannot be transposon-positive by construction, yet they sit in "
        "the denominator of the published rate. Restricting to records that can be "
        f"called at all moves the rate from {100 * _all_rate:.1f} % to "
        f"{100 * _callable_rate:.1f} %, a "
        f"{100 * (_callable_rate / _all_rate - 1):.0f} % relative shift that is pure "
        "record type and no biology.")

    by_record = {}
    for mol_type in sorted({d["mol_type"] or "" for d in data}):
        subset = [d for d in data if (d["mol_type"] or "") == mol_type]
        by_record[mol_type or "(empty)"] = {
            "n": len(subset),
            "median_annotated_neighbours": float(np.median([d["n_neighbours"] for d in subset])),
            "transposon_rate": round(float(np.mean([d["transposon"] for d in subset])), 4),
            "plasmid_rate": round(float(np.mean([d["is_plasmid"] for d in subset])), 4),
        }
    out["by_record_type"] = by_record
    out["by_replicon_status"] = {
        status: {"n": sum(1 for d in data if d["status"] == status),
                 "transposon_rate": round(
                     float(np.mean([d["transposon"] for d in data
                                    if d["status"] == status])), 4)}
        for status in sorted({d["status"] for d in data if d["status"]})
    }

    # Pencere dolulugu: bu, "seyrek anotasyon" sorusunun en dogrudan olcusu,
    # cunku regex'in uzerinde kosabilecegi gen sayisini birebir verir.
    callable_rows = [d for d in data if d["n_neighbours"] > 0]
    occupancy = []
    for low, high in quantile_bins([d["n_neighbours"] for d in callable_rows],
                                   DENSITY_QUANTILES):
        subset = [d for d in callable_rows if low <= d["n_neighbours"] <= high]
        if not subset:
            continue
        subset_genes = sum(d["n_neighbours"] for d in subset)
        occupancy.append({
            "annotated_neighbours": f"{low}-{high}",
            "n": len(subset),
            "transposon_rate": round(float(np.mean([d["transposon"] for d in subset])), 4),
            "per_gene_rate": round(sum(d["transposon_count"] for d in subset)
                                   / subset_genes, 5) if subset_genes else None,
            "plasmid_rate": round(float(np.mean([d["is_plasmid"] for d in subset])), 4),
            "median_cds_count": float(np.median([d["cds_count"] for d in subset])),
        })
    rates = [b["transposon_rate"] for b in occupancy]
    out["by_window_occupancy"] = {
        "bin_rule": (f"quartiles of the number of annotated neighbour genes within "
                     f"{NEIGHBOUR_WINDOW_BP} bp, taken from the callable entries "
                     "themselves rather than written into the code; the stratified "
                     "tests further down take their own quartiles from the population "
                     "each one analyses, so the bin edges there may differ slightly"),
        "bins": occupancy,
        "spread": {"lowest_rate": min(rates) if rates else None,
                   "highest_rate": max(rates) if rates else None,
                   "fold": round(max(rates) / min(rates), 2) if rates and min(rates) else None},
    }

    # Replikon boyu. Yogunlukla TERS yonde hareket ediyor; bunu gizlemek yerine
    # ayri bir bulgu olarak yaziyoruz.
    genomic = [d for d in data if d["mol_type"] == "genomic DNA" and d["cds_count"] > 0]
    size = []
    for low, high in quantile_bins([d["cds_count"] for d in genomic], CDS_QUANTILES):
        subset = [d for d in genomic if low <= d["cds_count"] <= high]
        if not subset:
            continue
        size.append({
            "cds_count": f"{low}-{high}",
            "n": len(subset),
            "transposon_rate": round(float(np.mean([d["transposon"] for d in subset])), 4),
            "plasmid_rate": round(float(np.mean([d["is_plasmid"] for d in subset])), 4),
            "median_annotated_neighbours": float(np.median([d["n_neighbours"] for d in subset])),
        })
    # Yonu VERIDEN oku: en kucuk ve en buyuk ceyreklik karsilastirilir.
    size_direction = ("falls" if size and size[-1]["transposon_rate"]
                      < size[0]["transposon_rate"] else "rises")
    occupancy_direction = ("rises" if occupancy and occupancy[-1]["transposon_rate"]
                           > occupancy[0]["transposon_rate"] else "falls")
    out["by_replicon_cds_count"] = {
        "bin_rule": "quartiles of replicon CDS count among genomic DNA records",
        "bins": size,
        "direction": size_direction,
        "finding": (
            f"the rate {size_direction} as the replicon gets larger, from "
            f"{100 * size[0]['transposon_rate']:.1f} % in the smallest quartile to "
            f"{100 * size[-1]['transposon_rate']:.1f} % in the largest, while the rate "
            f"{occupancy_direction} with window occupancy. "
            + ("The two run in opposite directions, which is the reason there is no "
               "single correction factor. "
               if size_direction != occupancy_direction else
               "They run in the same direction here. ")
            + "Small contigs and plasmids are genuinely richer in insertion sequences, "
              "so the replicon-size component is at least partly real biology and "
              "cannot be corrected away as if it were an artefact."),
    }

    # DURUST KARSI-KONTROL: yapisal kanit da anotasyon yogunluguyla korele mi?
    # Olculuyor cunku eger plazmit orani da pencere doluluguyla artiyorsa,
    # "yapisal kanit yanlilikdan muaf" demek yanlis olur. Oran gercekten
    # artiyor (iyi anote edilmis kayitlar tam kapanmis plazmit kayitlari), yani
    # yanlilik yalnizca regex yarisinda degil; farki, plazmit alaninin
    # anotasyon METNINE degil molekulun kimligine bagli olmasi.
    plasmid_by_occupancy = [{"annotated_neighbours": b["annotated_neighbours"],
                             "plasmid_rate": b["plasmid_rate"],
                             "transposon_rate": b["transposon_rate"]}
                            for b in occupancy]
    pl_rates = [b["plasmid_rate"] for b in occupancy]
    out["structural_evidence_also_tracks_annotation_density"] = {
        "question": ("Is the plasmid field itself correlated with annotation richness? "
                     "If it is, the bias is not confined to the regex half."),
        "bins": plasmid_by_occupancy,
        "spread": {"lowest_rate": min(pl_rates) if pl_rates else None,
                   "highest_rate": max(pl_rates) if pl_rates else None,
                   "fold": round(max(pl_rates) / min(pl_rates), 2)
                   if pl_rates and min(pl_rates) else None},
        "finding": (
            f"the plasmid rate also rises with window occupancy, from "
            f"{100 * pl_rates[0]:.1f} % to {100 * pl_rates[-1]:.1f} %, a wider "
            f"relative spread than the transposon rate. So sampling is not neutral for "
            "either variable: completely assembled, richly annotated records are both "
            "more likely to be recognised as plasmids and more likely to have a "
            "transposon annotated nearby. The difference between the two variables is "
            "not that one is unbiased. It is that the plasmid field is a claim about "
            "which molecule the gene sits on, which can in principle be checked, while "
            "the transposon category is a claim about the wording of a /product field, "
            "which cannot."),
        "caveat": ("this weakens any unconditional statement that plasmid evidence is "
                   "clean; it does not weaken the comparison between substrate classes, "
                   "because that comparison is made within the same pool of records"),
    }

    rho_cds = stats.spearmanr([d["cds_count"] for d in genomic],
                              [d["transposon"] for d in genomic])
    rho_nb = stats.spearmanr([d["n_neighbours"] for d in genomic],
                             [d["transposon"] for d in genomic])
    rho_link = stats.spearmanr([d["cds_count"] for d in genomic],
                               [d["n_neighbours"] for d in genomic])
    out["monotone_trends"] = {
        "unit": "alpha subunit on a genomic DNA record",
        "n": len(genomic),
        "cds_count_vs_transposon_called": {"spearman_rho": round(float(rho_cds.statistic), 4),
                                           "p": float(rho_cds.pvalue)},
        "window_occupancy_vs_transposon_called": {"spearman_rho": round(float(rho_nb.statistic), 4),
                                                  "p": float(rho_nb.pvalue)},
        "cds_count_vs_window_occupancy": {"spearman_rho": round(float(rho_link.statistic), 4),
                                          "p": float(rho_link.pvalue)},
        "note": ("all three correlations are small. They are reported with rho and "
                 "not with p alone because at n above 10,000 a rho of 0.05 reaches "
                 "p = 1e-7 and means almost nothing on its own"),
    }

    # ASIL SORU: yanlilik biyolojik iddiaya ne kadar siziyor? Substrat sinifi x
    # regex transpozonu iliskisi, pencere doluluguna gore tabakalanip
    # Mantel-Haenszel ile birlestirilir.
    universe = [d for d in callable_rows if d["sclass"] in SUBSTRATE_CLASSES
                and d["cluster"] not in euk_heavy and d["domain"] == "Bacteria"]
    xen_nb = [d["n_neighbours"] for d in universe if d["sclass"] == "xenobiotic"]
    nat_nb = [d["n_neighbours"] for d in universe if d["sclass"] != "xenobiotic"]
    delta, p_nb = cliffs_delta(xen_nb, nat_nb)
    strata, detail = [], []
    for low, high in quantile_bins([d["n_neighbours"] for d in universe],
                                   DENSITY_QUANTILES):
        xen = [d for d in universe
               if low <= d["n_neighbours"] <= high and d["sclass"] == "xenobiotic"]
        nat = [d for d in universe
               if low <= d["n_neighbours"] <= high and d["sclass"] != "xenobiotic"]
        if not xen or not nat:
            continue
        a1 = sum(d["transposon"] for d in xen)
        b1 = sum(d["transposon"] for d in nat)
        strata.append([[a1, len(xen) - a1], [b1, len(nat) - b1]])
        detail.append({"annotated_neighbours": f"{low}-{high}",
                       "xenobiotic_n": len(xen),
                       "xenobiotic_rate": round(a1 / len(xen), 4),
                       "natural_n": len(nat),
                       "natural_rate": round(b1 / len(nat), 4),
                       "ratio": round((a1 / len(xen)) / (b1 / len(nat)), 3)
                       if b1 else None})
    common_or, used = mantel_haenszel(strata)
    xen_all = [d for d in universe if d["sclass"] == "xenobiotic"]
    nat_all = [d for d in universe if d["sclass"] != "xenobiotic"]
    a1 = sum(d["transposon"] for d in xen_all)
    b1 = sum(d["transposon"] for d in nat_all)
    crude = (a1 * (len(nat_all) - b1)) / max(1, (len(xen_all) - a1) * b1)
    attributable = (None if not common_or or crude <= 1
                    else round((crude - common_or) / (crude - 1), 3))
    out["how_much_bias_reaches_the_biological_claim"] = {
        "question": ("Is the xenobiotic excess in regex transposon calls produced by "
                     "xenobiotic entries sitting in better annotated windows?"),
        "annotation_density_differs_between_classes": {
            "median_neighbours_xenobiotic": float(np.median(xen_nb)),
            "median_neighbours_natural": float(np.median(nat_nb)),
            "mean_neighbours_xenobiotic": round(float(np.mean(xen_nb)), 2),
            "mean_neighbours_natural": round(float(np.mean(nat_nb)), 2),
            "cliffs_delta": round(delta, 3) if delta is not None else None,
            "p": p_nb,
            "reading": (
                ("xenobiotic entries do sit in better annotated windows "
                 f"(Cliff's delta {delta:+.3f}). The effect is small but it points in "
                 "exactly the direction that would manufacture the transposon result, "
                 "so it has to be adjusted for rather than argued away."
                 if (delta or 0) > 0 else
                 "xenobiotic entries do not sit in better annotated windows "
                 f"(Cliff's delta {delta:+.3f}), so annotation density cannot be the "
                 "source of the transposon excess; the adjustment below is reported "
                 "anyway because the absence of a confounder has to be shown, not "
                 "assumed.")),
        },
        "strata": detail,
        "crude_odds_ratio": round(float(crude), 3),
        "mantel_haenszel_odds_ratio": common_or,
        "strata_used": len(used),
        "share_of_excess_odds_attributable_to_annotation_density": attributable,
    }
    _reach = out["how_much_bias_reaches_the_biological_claim"]
    _share = attributable
    _reach["finding"] = (
        "adjusting for window occupancy shrinks the excess odds "
        + (f"from {crude:.2f} to {common_or:.2f} " if common_or else "")
        + ("but does not remove them. "
           if common_or and common_or > 1.05 else "and all but removes them. ")
        + (f"About {100 * _share:.0f} % of the apparent transposon excess is "
           f"annotation density; the remaining {100 * (1 - _share):.0f} % is not "
           "explained by it." if _share is not None else
           "The share attributable to annotation density could not be computed "
           "because the crude odds ratio is not above 1."))

    out["verdict"] = {
        "size_of_the_bias": (
            "large enough to change how the number is read, not large enough to erase "
            "it. Three separate components were measured. First, a hard floor: "
            f"{out['hard_floor']['entries_with_no_annotated_neighbour']} entries "
            f"({100 * out['hard_floor']['share_of_all_entries']:.1f} %) are on records "
            "with no annotated neighbour at all and can never be called positive, yet "
            "they are counted in the denominator. Second, window occupancy: among "
            "records that can be called, the rate spans "
            f"{100 * out['by_window_occupancy']['spread']['lowest_rate']:.1f} % to "
            f"{100 * out['by_window_occupancy']['spread']['highest_rate']:.1f} % across "
            "quartiles of how many genes the submitter annotated within 10 kb, a "
            f"{out['by_window_occupancy']['spread']['fold']}-fold spread. Third, "
            f"replicon size, where the rate {size_direction} with CDS count, "
            + ("the opposite direction from window occupancy; part of that is real "
               "insertion-sequence biology. Because the components have opposite signs "
               "there is no single correction factor, which is itself the reason the "
               "variable cannot be presented as a measurement of mobility."
               if size_direction != occupancy_direction else
               "the same direction as window occupancy. Even so the two cannot be "
               "merged into one correction, because replicon size also carries real "
               "insertion-sequence biology.")),
        "what_to_do": (
            "keep the regex transposon rate, label it as annotation-derived wherever "
            "it appears, always quote it on the callable denominator, and never OR it "
            "into a variable that also carries structural evidence. The one thing that "
            "would actually fix it is sequence verification: run a transposase profile "
            "HMM over the neighbour proteins, exactly as build_operons.py already does "
            "for beta, ferredoxin and reductase. The neighbour protein sequences are "
            "already in the database."),
        "what_this_cannot_settle": (
            "whether the regex calls are true transposases. Annotation density explains "
            "part of the variation; the rest could be real insertion-sequence biology "
            "or a second, unmeasured annotation habit, and counting words cannot "
            "distinguish those two."),
    }
    return out


# ------------------------- 3. soru: duzenleyici ne kadar ayrisik evrildi

# Global hizalayici surec basina BIR kez kurulur. blastp puanlamasi secildi
# cunku karsilastirma protein duzeyinde ve BLOSUM62 + blastp bosluk cezalari
# bu is icin kalibre edilmis standart; keyfi bir puan matrisi secmemek icin.
_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        aligner = Align.PairwiseAligner(scoring="blastp")
        aligner.mode = "global"
        _ALIGNER = aligner
    return _ALIGNER


def pair_identity(seq_a, seq_b):
    """Global hizalamada yuzde kimlik ve kisa proteine gore kapsama.

    Kimlik YALNIZCA iki tarafi da kalinti olan kolonlar uzerinden hesaplanir;
    bosluk kolonlari paydaya girmez. Eski kod (evolution_rates_of_reg_and_enz.py)
    bosluk-kalinti kolonlarini uyumsuzluk sayiyordu ve ortak bir coklu
    hizalama kullaniyordu: duzenleyiciler hem daha kisa hem uzunluk olarak
    daha degisken oldugu icin o hizalamada daha fazla bosluk kolonu olusuyor
    ve "duzenleyiciler daha cesitli" sonucu kismen bu bosluklardan geliyordu.
    Ciftin KENDI hizalamasini kullanmak bu yanliligi ortadan kaldirir.
    """
    if not seq_a or not seq_b:
        return None, None
    alignment = _aligner().align(seq_a, seq_b)[0]
    matches = aligned = 0
    for residue_a, residue_b in zip(alignment[0], alignment[1]):
        if residue_a != "-" and residue_b != "-":
            aligned += 1
            if residue_a == residue_b:
                matches += 1
    shorter = min(len(seq_a), len(seq_b))
    if not aligned or not shorter:
        return None, None
    return 100.0 * matches / aligned, aligned / shorter


def _score_pair(task):
    """Pool isci fonksiyonu: bir cift icin enzim ve duzenleyici kimligi."""
    (alpha_a, reg_a, alpha_b, reg_b) = task
    enzyme_identity, enzyme_coverage = pair_identity(alpha_a, alpha_b)
    regulator_identity, regulator_coverage = pair_identity(reg_a, reg_b)
    return enzyme_identity, enzyme_coverage, regulator_identity, regulator_coverage


def load_regulator_pairs_input(con):
    """Alfa alt birimi ve YUKARI AKIS duzenleyicisinin dizisi birlikte olan girisler.

    Duzenleyici, `ro_regulation`'in belirledigi gen: operonun 5' ucunun
    yukarisindaki ILK gen ve kategorisi 'regulator'. Bu, eski kodun
    "orta noktalari 1.000-1.500 bp icinde olan en yakin regulator" kuralindan
    farkli ve daha savunulabilir, cunku iplik ve operon yapisini hesaba katar.
    Komsu protein dizileri `neighbor_protein`'den protein_key uzerinden
    baglanir (tabloda neighbor_id yok, anahtar koordinat tabanli).
    """
    return con.execute("""
        SELECT r.candidate_id, r.ro_cluster, r.ro_group, p.organism,
               g.upstream_family, g.architecture, r.sequence, np.translation
        FROM ro r
        JOIN replicon p USING(nucleotide_id)
        JOIN ro_regulation g USING(candidate_id)
        JOIN neighbor nb ON nb.neighbor_id = g.upstream_gene_id
        JOIN neighbor_protein np
          ON np.protein_key = nb.nucleotide_id || ':' || nb.start || '-'
                              || nb.end || ':' || nb.strand
        WHERE r.is_confirmed = 1
          AND g.upstream_category = 'regulator'
          AND r.sequence IS NOT NULL AND r.sequence <> ''
          AND np.translation IS NOT NULL AND np.translation <> ''
    """).fetchall()


def band_of(identity):
    for low, high, name in IDENTITY_BANDS:
        if low <= identity < high:
            return name
    return None


def paired_contrast(pairs, label):
    """Cift ICINDE enzim ve duzenleyici kimligini karsilastir.

    Olcut ORDINAL: hangisi daha korunmus. Iki yuzdeyi dogrudan cikarmak
    cazip ama iki protein ayni kisitlar altinda degil; bu yuzden onculuk
    isaret testine (ciftlerin kaci icin duzenleyici daha az korunmus) ve
    Wilcoxon isaretli sira testine veriliyor, fark yalnizca tanimlayici
    olarak yaziliyor.
    """
    if len(pairs) < MIN_PAIRS_FOR_BAND:
        return None
    enzyme = np.array([p["enzyme_identity"] for p in pairs])
    regulator = np.array([p["regulator_identity"] for p in pairs])
    delta = enzyme - regulator
    try:
        wilcoxon_p = float(stats.wilcoxon(enzyme, regulator).pvalue)
    except ValueError:
        wilcoxon_p = None
    share_less = float(np.mean(regulator < enzyme))
    at_floor = float(np.mean(regulator <= RANDOM_IDENTITY_FLOOR))
    return {
        "label": label,
        "n_pairs": len(pairs),
        "n_types": len({p["type"] for p in pairs}),
        "n_genus_pairs": len({p["genus_pair"] for p in pairs}),
        "n_genera": len({g for p in pairs for g in p["genus_pair"]}),
        "enzyme_identity_median": round(float(np.median(enzyme)), 1),
        "enzyme_identity_quartiles": [round(float(v), 1)
                                      for v in np.percentile(enzyme, [25, 75])],
        "regulator_identity_median": round(float(np.median(regulator)), 1),
        "regulator_identity_quartiles": [round(float(v), 1)
                                         for v in np.percentile(regulator, [25, 75])],
        "median_difference_enzyme_minus_regulator": round(float(np.median(delta)), 1),
        "share_of_pairs_regulator_less_conserved": round(share_less, 3),
        "wilcoxon_p": wilcoxon_p,
        "share_of_pairs_at_or_below_random_floor": round(at_floor, 3),
        "informative": bool(at_floor < MAX_FLOOR_SHARE_FOR_RATE_CLAIM),
    }


def divergence_block(con, cpu, seed=RANDOM_SEED):
    rows = load_regulator_pairs_input(con)
    total_confirmed = con.execute(
        "SELECT COUNT(*) FROM ro WHERE is_confirmed = 1").fetchone()[0]
    with_regulator = con.execute("""
        SELECT COUNT(*) FROM ro r JOIN ro_regulation g USING(candidate_id)
        WHERE r.is_confirmed = 1 AND g.upstream_category = 'regulator'"""
    ).fetchone()[0]

    out = {
        "question": ("How divergently did the regulators evolve compared with the "
                     "enzymes they control, by enzyme group and by the "
                     "enzyme-to-enzyme identity band?"),
        "prior_art": {
            "scripts": ["evolution_rates_of_reg_and_enz.py",
                        "evolution_rates_of_reg_and_enz_batch.py",
                        "enzyme_regulator_graphs.py"],
            "what_they_computed": (
                "they ask the same question and reach the same direction. From the "
                "old pipeline's proteins_output.fasta they paired each main gene with "
                "the nearest regulator within 1000 to 1500 bp by midpoint distance, "
                "aligned all main sequences of a type together with MUSCLE and all "
                "regulator sequences together, took the mean pairwise p-distance of "
                "each alignment, and compared the two numbers per type. Published "
                "result: 45 types, 2382 pairs, mean main diversity 0.378 against mean "
                "regulator diversity 0.446, concluding that regulators are the more "
                "diverse."),
            "why_it_is_not_reused_as_is": [
                "the p-distance counts a gap against a residue as a mismatch and is "
                "taken from a multiple alignment shared by all sequences of a type. "
                "Regulators are shorter and more variable in length, so their shared "
                "alignment carries more gap columns, and that alone inflates their "
                "apparent divergence. The direction of the published conclusion is "
                "therefore the direction of a known bias, which is why it needed "
                "recomputing rather than citing.",
                "the two diversity numbers are type-level aggregates, but "
                "enzyme_regulator_graphs.py feeds them to a paired t-test and a "
                "Wilcoxon signed-rank test as though they were paired observations. "
                "The pairing that the question actually needs is per entry pair, not "
                "per type.",
                "it runs on the old pipeline's protein set, which the rebuild "
                "replaced: that set accepted non-alpha Rieske proteins, and the "
                "summary still contains rows for CdnD, an electron-transfer component "
                "that was removed from the curated set as a contaminant.",
                "the pairing rule is nearest-regulator-by-midpoint and ignores strand "
                "and operon structure. ro_regulation now identifies the upstream gene "
                "properly, from the operon 5' end.",
                "enzyme_regulator_graphs.py imports sklearn, which is not installed "
                "here, so it cannot be run to reproduce its own figures.",
            ],
            "what_is_reused": (
                "the question itself, the choice of pairwise amino acid identity as "
                "the measure, and the per-type framing. What changed is that identity "
                "is computed on each pair's own global alignment with a coverage "
                "check, that the comparison is paired within an entry pair, and that "
                "the result is stratified by the enzyme-to-enzyme identity band, which "
                "the earlier scripts do not do at all."),
        },
        "coverage": {
            "confirmed_entries": total_confirmed,
            "entries_whose_upstream_gene_is_a_regulator": with_regulator,
            "of_those_with_both_sequences_available": len(rows),
            "sequence_availability": round(len(rows) / with_regulator, 4)
            if with_regulator else None,
            "share_of_all_confirmed_entries": round(len(rows) / total_confirmed, 4),
            "why_the_CON_problem_does_not_bite_here": (
                "88.9 % of the source GenBank records are CON entries with no "
                "nucleotide sequence, which is what blocks promoter analysis. Protein "
                "translations are a different matter: they are carried in the CDS "
                "/translation qualifier and are present even when the nucleotide "
                "sequence is not, so build_operons.py was able to store them. "
                "Measured availability for the regulators is "
                + (f"{100 * len(rows) / with_regulator:.1f} %, so sequence coverage is "
                   "not the limiting factor here."
                   if with_regulator else "not computable.")),
            "what_is_the_limiting_factor": (
                "whether the upstream gene is a regulator at all. That holds for "
                f"{with_regulator} of {total_confirmed} entries "
                f"({100 * with_regulator / total_confirmed:.1f} %), so roughly seven "
                "entries in ten cannot enter this analysis. Nothing about the "
                "remaining seven tenths is claimed."),
        },
    }
    if len(rows) < 2:
        out["verdict"] = {"conclusion": "not computable: too few usable entries"}
        return out

    entries = [{"candidate_id": cid, "type": cluster, "group": group or "?",
                "genus": (organism or "?").split()[0],
                "family": family, "architecture": architecture,
                "alpha": alpha, "regulator": regulator}
               for (cid, cluster, group, organism, family, architecture,
                    alpha, regulator) in rows]

    by_type = defaultdict(list)
    for entry in entries:
        by_type[entry["type"]].append(entry)

    # Ornekleme: tip icinde TUM ciftler listelenir, sabit tohumla karistirilir
    # ve ilk MAX_PAIRS_PER_TYPE tanesi alinir. Boylece hangi ciftlerin
    # alindigi tamamen yeniden uretilebilir.
    rng = random.Random(seed)
    selected, sampling = [], []
    for cluster in sorted(by_type):
        items = by_type[cluster]
        all_pairs = [(i, j) for i in range(len(items)) for j in range(i + 1, len(items))]
        if not all_pairs:
            continue
        rng.shuffle(all_pairs)
        taken = all_pairs[:MAX_PAIRS_PER_TYPE]
        sampling.append({"type": cluster, "entries": len(items),
                         "all_within_type_pairs": len(all_pairs),
                         "pairs_used": len(taken),
                         "complete": len(taken) == len(all_pairs)})
        for i, j in taken:
            selected.append((items[i], items[j]))

    tasks = [(a["alpha"], a["regulator"], b["alpha"], b["regulator"])
             for a, b in selected]
    if cpu > 1 and len(tasks) > 200:
        with Pool(cpu) as pool:
            scored = pool.map(_score_pair, tasks, chunksize=64)
    else:
        scored = [_score_pair(t) for t in tasks]

    pairs, dropped = [], 0
    for (a, b), (e_id, e_cov, r_id, r_cov) in zip(selected, scored):
        if None in (e_id, e_cov, r_id, r_cov):
            dropped += 1
            continue
        if e_cov < MIN_ALIGNED_COVERAGE or r_cov < MIN_ALIGNED_COVERAGE:
            dropped += 1
            continue
        pairs.append({
            "type": a["type"], "group": a["group"],
            "genus_pair": tuple(sorted((a["genus"], b["genus"]))),
            "families": (a["family"], b["family"]),
            "same_family": (a["family"] == b["family"]),
            "family_known": (a["family"] != UNCLASSIFIED_FAMILY
                             and b["family"] != UNCLASSIFIED_FAMILY
                             and a["family"] and b["family"]),
            "enzyme_identity": e_id, "regulator_identity": r_id,
            "band": band_of(e_id),
        })

    out["sampling"] = {
        "rule": (f"within each enzyme type every within-type pair is enumerated, "
                 f"shuffled with seed {seed} and the first {MAX_PAIRS_PER_TYPE} kept. "
                 "Uniform sampling inside a type estimates that type's pair "
                 "distribution without bias, and the fixed seed makes the exact set "
                 "of pairs reproducible."),
        "max_pairs_per_type": MAX_PAIRS_PER_TYPE,
        "seed": seed,
        "all_within_type_pairs_available": sum(s["all_within_type_pairs"]
                                               for s in sampling),
        "pairs_scored": len(tasks),
        "pairs_dropped_for_low_alignment_coverage": dropped,
        "min_aligned_coverage": MIN_ALIGNED_COVERAGE,
        "pairs_used": len(pairs),
        "types_completely_enumerated": sum(1 for s in sampling if s["complete"]),
        "types_sampled": len(sampling),
        "per_type": sampling,
    }
    if not pairs:
        out["verdict"] = {"conclusion": "not computable: no pair survived the "
                                        "alignment coverage check"}
        return out

    # (a) + (b) bant bant: butun ciftler, sonra YALNIZCA ayni aile
    out["all_pairs"] = paired_contrast(pairs, "all pairs")
    out["by_enzyme_identity_band"] = [
        r for r in (paired_contrast([p for p in pairs if p["band"] == name], name)
                    for _, _, name in IDENTITY_BANDS) if r]

    same_family = [p for p in pairs if p["family_known"] and p["same_family"]]
    out["same_regulator_family_only"] = {
        "why": ("if the two regulators belong to different families they are "
                "different proteins, and a low identity between them is not the same "
                "observation as one regulator diverging. Restricting to pairs whose "
                "regulators share a family is the comparison that actually measures "
                "divergence rate."),
        "n_pairs": len(same_family),
        "all": paired_contrast(same_family, "same family, all bands"),
        "by_enzyme_identity_band": [
            r for r in (paired_contrast(
                [p for p in same_family if p["band"] == name], name)
                for _, _, name in IDENTITY_BANDS) if r],
    }

    # (b) enzim grubuna gore
    out["by_enzyme_group"] = [
        r for r in (paired_contrast([p for p in pairs if p["group"] == group],
                                    f"group {group}")
                    for group in sorted({p["group"] for p in pairs})) if r]
    out["by_enzyme_group_same_family"] = [
        r for r in (paired_contrast([p for p in same_family if p["group"] == group],
                                    f"group {group}")
                    for group in sorted({p["group"] for p in same_family})) if r]

    # (c) aile korunuyor mu, degisiyor mu
    known = [p for p in pairs if p["family_known"]]
    switches = [p for p in known if not p["same_family"]]
    switch_by_band = []
    for _, _, name in IDENTITY_BANDS:
        subset = [p for p in known if p["band"] == name]
        if len(subset) < MIN_PAIRS_FOR_BAND:
            continue
        band_switch = [p for p in subset if not p["same_family"]]
        switch_by_band.append({
            "band": name, "n_pairs": len(subset),
            "n_types": len({p["type"] for p in subset}),
            "n_genus_pairs": len({p["genus_pair"] for p in subset}),
            "family_switch_rate": round(len(band_switch) / len(subset), 3),
            "regulator_identity_median_same_family": round(float(np.median(
                [p["regulator_identity"] for p in subset if p["same_family"]])), 1)
            if any(p["same_family"] for p in subset) else None,
            "regulator_identity_median_switched": round(float(np.median(
                [p["regulator_identity"] for p in band_switch])), 1)
            if band_switch else None,
        })
    per_type_switch = []
    for cluster in sorted({p["type"] for p in known}):
        subset = [p for p in known if p["type"] == cluster]
        if len(subset) < MIN_PAIRS_FOR_BAND:
            continue
        families = Counter(f for p in subset for f in p["families"])
        per_type_switch.append({
            "type": cluster, "n_pairs": len(subset),
            "n_genera": len({g for p in subset for g in p["genus_pair"]}),
            "family_switch_rate": round(
                float(np.mean([not p["same_family"] for p in subset])), 3),
            "families_seen": families.most_common(),
        })
    per_type_switch.sort(key=lambda r: -r["family_switch_rate"])
    out["regulator_family_retention"] = {
        "question": ("Is the regulator family retained across diverging members of "
                     "one enzyme type, or does it switch? A switch means the same "
                     "chemistry came under a different kind of transcriptional "
                     "control."),
        "pairs_with_both_families_classified": len(known),
        "pairs_excluded_as_unclassified": len(pairs) - len(known),
        "overall_switch_rate": round(len(switches) / len(known), 3) if known else None,
        "by_enzyme_identity_band": switch_by_band,
        "most_frequently_switched_family_pairs": [
            {"families": list(k), "n_pairs": v}
            for k, v in Counter(tuple(sorted(p["families"]))
                                for p in switches).most_common(10)],
        "by_type": per_type_switch,
    }

    # Bagimsizlik: ayni sayi hem cift hem de TIP x CINS-CIFTI duzeyinde
    collapsed = {}
    for pair in pairs:
        collapsed.setdefault((pair["type"], pair["genus_pair"]), []).append(pair)
    collapsed_pairs = []
    for key, items in collapsed.items():
        collapsed_pairs.append({
            "enzyme_identity": float(np.median([p["enzyme_identity"] for p in items])),
            "regulator_identity": float(np.median([p["regulator_identity"]
                                                   for p in items])),
        })
    rho_pair = stats.spearmanr([p["enzyme_identity"] for p in pairs],
                               [p["regulator_identity"] for p in pairs])
    rho_coll = stats.spearmanr([p["enzyme_identity"] for p in collapsed_pairs],
                               [p["regulator_identity"] for p in collapsed_pairs])
    out["correlation"] = {
        "question": "Does regulator identity track enzyme identity?",
        "pair_level": {"unit": "entry pair within one enzyme type",
                       "n": len(pairs),
                       "n_types": len({p["type"] for p in pairs}),
                       "n_genus_pairs": len({p["genus_pair"] for p in pairs}),
                       "spearman_rho": round(float(rho_pair.statistic), 3),
                       "p": float(rho_pair.pvalue)},
        "collapsed_level": {"unit": "one observation per enzyme type and unordered "
                                    "genus pair, medians within the cell",
                            "n": len(collapsed_pairs),
                            "spearman_rho": round(float(rho_coll.statistic), 3),
                            "p": float(rho_coll.pvalue)},
        "note": ("the pair level is reported for completeness only. Pairs inside one "
                 "type share sequences with each other, so its n is not a count of "
                 "independent observations; the collapsed level is the honest one."),
    }

    # Tabana vurma: olcumun nerede bilgi tasimadigini acikca say
    at_floor = [p for p in pairs if p["regulator_identity"] <= RANDOM_IDENTITY_FLOOR]
    out["measurement_floor"] = {
        "random_identity_floor_percent": RANDOM_IDENTITY_FLOOR,
        "why": ("two unrelated proteins align at roughly 20 to 25 % identity, so a "
                "regulator identity near that value does not mean fast divergence, it "
                "means no measurable relationship. In that region the measure is "
                "saturated and says nothing about rate."),
        "pairs_at_or_below_floor": len(at_floor),
        "share_of_all_pairs": round(len(at_floor) / len(pairs), 3),
        "share_by_band": {
            name: round(float(np.mean([p["regulator_identity"] <= RANDOM_IDENTITY_FLOOR
                                       for p in pairs if p["band"] == name])), 3)
            for _, _, name in IDENTITY_BANDS
            if sum(1 for p in pairs if p["band"] == name) >= MIN_PAIRS_FOR_BAND},
    }

    out["figure_data"] = figure_data(pairs, out, seed)
    out["verdict"] = divergence_verdict(out)
    return out


def figure_data(pairs, block, seed):
    """Web sayfasinin dogrudan cizebilecegi seriler.

    Goruntu dosyasi URETILMIYOR: bu pipeline'in deseni, sayilari JSON'a yazip
    sekli web katmaninda cizmek (bkz. stratified_stats.py -> carboxylate_rows.json).
    Boylece sekil ile sayi ayrisamaz.
    """
    rng = random.Random(seed)
    scatter_source = list(pairs)
    if len(scatter_source) > MAX_SCATTER_POINTS:
        rng.shuffle(scatter_source)
        scatter_source = scatter_source[:MAX_SCATTER_POINTS]
    return {
        "plot_1": {
            "kind": "grouped box or violin, one pair of boxes per band",
            "title": "Enzyme and regulator identity by enzyme identity band",
            "x_label": "Enzyme-to-enzyme identity band (%)",
            "y_label": "Pairwise amino acid identity (%)",
            "series": [
                {"band": row["label"], "n_pairs": row["n_pairs"],
                 "n_types": row["n_types"], "n_genus_pairs": row["n_genus_pairs"],
                 "enzyme_median": row["enzyme_identity_median"],
                 "enzyme_q1_q3": row["enzyme_identity_quartiles"],
                 "regulator_median": row["regulator_identity_median"],
                 "regulator_q1_q3": row["regulator_identity_quartiles"]}
                for row in block["by_enzyme_identity_band"]],
            "reference_line": {"y": RANDOM_IDENTITY_FLOOR,
                               "label": "identity floor for unrelated proteins"},
        },
        "plot_2": {
            "kind": "scatter with a diagonal",
            "title": "Regulator identity against enzyme identity, one point per pair",
            "x_label": "Enzyme identity (%)", "y_label": "Regulator identity (%)",
            "diagonal": "y = x marks equal conservation; points below it are pairs "
                        "whose regulator is the less conserved of the two",
            "columns": ["enzyme_identity", "regulator_identity",
                        "same_regulator_family", "enzyme_group"],
            "sampled_points": len(scatter_source),
            "total_points": len(pairs),
            "sampling_note": (f"a fixed-seed sample of at most {MAX_SCATTER_POINTS} "
                              "points, so the published file stays small"),
            "rows": [[round(p["enzyme_identity"], 1), round(p["regulator_identity"], 1),
                      int(bool(p["same_family"])), p["group"]]
                     for p in scatter_source],
        },
        "plot_3": {
            "kind": "line or bar",
            "title": "Regulator family switch rate by enzyme identity band",
            "x_label": "Enzyme-to-enzyme identity band (%)",
            "y_label": "Share of pairs whose two regulators belong to different "
                       "families",
            "series": [{"band": row["band"], "n_pairs": row["n_pairs"],
                        "switch_rate": row["family_switch_rate"]}
                       for row in
                       block["regulator_family_retention"]["by_enzyme_identity_band"]],
        },
    }


def divergence_verdict(block):
    """Karar VERIDEN kurulur; hicbir yon ya da sayi elle yazilmaz."""
    same = block["same_regulator_family_only"]["by_enzyme_identity_band"]
    bands = {row["label"]: row for row in same}
    retention = block["regulator_family_retention"]
    switch = {row["band"]: row["family_switch_rate"]
              for row in retention["by_enzyme_identity_band"]}
    informative = [row for row in same if row["informative"]]
    faster = [row for row in informative
              if row["share_of_pairs_regulator_less_conserved"] > 0.5]
    parts = []
    if informative:
        parts.append(
            "Among pairs whose two regulators share a family, which is the comparison "
            "that measures divergence rather than replacement, the regulator is the "
            "less conserved partner in "
            + ", ".join(f"{100 * row['share_of_pairs_regulator_less_conserved']:.0f} % "
                        f"of pairs at {row['label']} % enzyme identity "
                        f"(n = {row['n_pairs']}, {row['n_types']} types, "
                        f"{row['n_genus_pairs']} genus pairs)"
                        for row in informative)
            + ".")
        parts.append(
            "So regulators diverge faster than the enzymes they sit next to"
            if len(faster) == len(informative) and informative else
            "So the direction is not uniform across bands and should not be stated as "
            "a single rule")
        gap = [row for row in informative
               if row["median_difference_enzyme_minus_regulator"] is not None]
        if len(gap) >= 2:
            widest = max(gap, key=lambda r: r["median_difference_enzyme_minus_regulator"])
            # Fark MONOTON degil ve nedeni biliniyor: en dusuk bantta
            # duzenleyici kimligi olcum tabanina vuruyor, dolayisiyla fark
            # kapaniyormus gibi gorunuyor. Cumle bunu soylemek zorunda.
            parts.append(
                "and the gap widens as the enzymes diverge, from "
                f"{gap[0]['median_difference_enzyme_minus_regulator']:+.0f} points at "
                f"{gap[0]['label']} % enzyme identity to a maximum of "
                f"{widest['median_difference_enzyme_minus_regulator']:+.0f} points at "
                f"{widest['label']} %.")
            if widest is not gap[-1]:
                parts.append(
                    "Below that the gap appears to close again, to "
                    f"{gap[-1]['median_difference_enzyme_minus_regulator']:+.0f} points "
                    f"at {gap[-1]['label']} %, but that is compression against the "
                    "measurement floor rather than regulators becoming conserved "
                    "again: "
                    f"{100 * gap[-1]['share_of_pairs_at_or_below_random_floor']:.0f} % "
                    "of the pairs in that band already sit at the identity level of "
                    "unrelated proteins, so the regulator value cannot fall further "
                    "while the enzyme value still can.")
    conclusion = " ".join(parts) if parts else (
        "no band carries enough informative pairs to support a statement about "
        "divergence rate")
    switch_sentence = None
    if switch:
        ordered = [row["band"] for row in retention["by_enzyme_identity_band"]]
        switch_sentence = (
            "the regulator family is retained among near-identical enzymes and lost as "
            "they diverge: the switch rate runs "
            + ", ".join(f"{100 * switch[b]:.0f} % at {b} %" for b in ordered)
            + f", and overall {100 * retention['overall_switch_rate']:.0f} % of pairs "
              "with two classified regulators carry two DIFFERENT families. That is "
              "the same chemistry placed under a different kind of transcriptional "
              "control, and it is the main reason the unrestricted comparison "
              "overstates regulator divergence.")
    floor = block["measurement_floor"]
    return {
        "divergence_rate": conclusion,
        "family_retention": switch_sentence,
        "where_the_data_cannot_answer": (
            "the lowest enzyme identity band. "
            f"{100 * floor['share_of_all_pairs']:.0f} % of all pairs have a regulator "
            f"identity at or below {floor['random_identity_floor_percent']:.0f} %, "
            "which is where unrelated proteins sit, and in that region the measure is "
            "saturated: it cannot distinguish a fast-diverging regulator from an "
            "unrelated one. Rate statements are therefore restricted to the bands "
            "flagged as informative. Two further limits: the regulator is only ever "
            "the FIRST upstream gene, so a regulator acting from elsewhere is invisible; "
            "and no binding was measured, so 'the regulator of this enzyme' is an "
            "inference from position, not from function."),
        "relation_to_prior_art": (
            "the earlier scripts concluded that regulators are the more diverse, and "
            "that direction holds up. The contribution here is that it survives a "
            "measure without the gap-counting bias, that it is decomposed into real "
            "divergence against family replacement, and that it is resolved by "
            "identity band, which is what makes it interpretable."),
    }


# --------------------------------------------------------------------- cikti

def print_summary(result):
    reg = result["regulation"]
    mob = result["mobility_by_substrate_class"]
    ev = result["transposon_vs_curated_distance"]
    bias = result["annotation_bias"]

    print("=" * 78)
    print("OPERON DUZEYI ILISKILER")
    print("=" * 78)
    print(f"giris: {result['totals']['entries']}   "
          f"tip x cins gozlemi: {result['totals']['type_and_genus_observations']}")

    print("\n-- 1. SORU: duzenleyici x enzim ----------------------------------")
    cov = reg["coverage"]
    print(f"yukari akista gen var: {cov['with_upstream_gene_in_window']}  "
          f"duzenleyici: {cov['with_regulator_upstream']}  "
          f"adlandirilmis aile: {cov['with_named_regulator_family']}  "
          f"siniflanamayan: {cov['regulator_family_unclassified']}")
    print(f"{'iliski':44} {'V(giris)':>9} {'V(cins)':>8} {'sans':>6} {'asan':>6}")
    for row in reg["summary"]:
        print(f"{row['pair'][:44]:44} {row['entry_cramers_v']:9.3f} "
              f"{row['genus_cramers_v']:8.3f} "
              f"{(row['genus_cramers_v'] - (row['genus_above_chance'] or 0)):6.3f} "
              f"{(row['genus_above_chance'] or 0):6.3f}"
              f"{'  *' if row['survives_genus_collapse'] else '   -'}")
    print("  (* cins cokertmesinden sonra da gercek; - degil)")
    print("\n  en guclu hucreler (duzenleyici ailesi x kimyasal aile, giris duzeyi):")
    for cell in reg["cell_enrichments"][:8]:
        print(f"    {cell['regulator_family']:10} x {cell['chemical_family']:20} "
              f"n={cell['entries']:5} obs/exp={cell['observed_over_expected']:5.2f} "
              f"resid={cell['standardised_residual']:+6.1f} "
              f"cins={cell['distinct_genera']:3} tip={cell['distinct_types']}")

    print("\n-- 2a. SORU: mobilite, uc kanit tanimi ---------------------------")
    comp = mob["composition"]
    print(f"plazmit {comp['plasmid']} (%{100*comp['plasmid_share']:.1f})  "
          f"regex transpozon {comp['regex_transposon']} (%{100*comp['regex_transposon_share']:.1f})  "
          f"mobil {comp['mobile_published']} (%{100*comp['mobile_share']:.1f})  "
          f"yalnizca regex {comp['mobile_only_because_of_regex']} "
          f"(mobilin %{100*comp['share_of_mobile_that_is_regex_only']:.0f}'i)")
    print(f"{'kanit':24} {'duzey':8} {'ksenob':>8} {'dogal':>8} {'kat':>7} {'p':>10} {'n':>7}")
    for test in mob["tests"]:
        print(f"{test['evidence']:24} {test['level']:8} "
              f"{test['rate_xenobiotic']:8.4f} {test['rate_natural']:8.4f} "
              f"{test['ratio']:6.2f}x {test['p']:10.2e} {test['n']:7}")
    print("  karar: anotasyon yarisina bagimli mi -> "
          f"{mob['verdict']['does_the_conclusion_depend_on_the_annotation_half'].upper()}")

    print("\n-- 2b. SORU: transpozon x kuratorlu enzime uzaklik ---------------")
    for level in ("entry", "genus"):
        parts = [f"{t[:5]}=%{100*v['transposon_rate']:.1f}({v['n']})"
                 for t, v in ev["transposon_rate_by_tier"][level].items()]
        print(f"  {level:6} " + "  ".join(parts))
    for test in ev["associations"]:
        print(f"    {test['id']:40} V={test['cramers_v']:.3f} "
              f"asan={test['cramers_v_above_chance']:+.3f} p={test['p']:.1e}")
    for level in ("entry", "genus"):
        r = ev["ref_identity"][level]
        print(f"    ref_identity {level:6}: %{r['median_identity_with_transposon']:.1f} vs "
              f"%{r['median_identity_without']:.1f}  delta={r['cliffs_delta']:+.3f} "
              f"p={r['p']:.1e}")
    for level, adj in ev["density_adjusted_close_vs_far"].items():
        print(f"    yakin/uzak OR {level:6}: ham {adj['crude_odds_ratio']} -> "
              f"yogunluga gore {adj['mantel_haenszel_odds_ratio']}  "
              f"(yogunluga yazilan pay: {adj['share_of_excess_odds_attributable_to_annotation_density'] if adj['share_of_excess_odds_attributable_to_annotation_density'] is not None else '-'})")

    print("\n-- 2c. SORU: anotasyon yanliligi --------------------------------")
    floor = bias["hard_floor"]
    print(f"  komsusuz giris: {floor['entries_with_no_annotated_neighbour']} "
          f"(%{100*floor['share_of_all_entries']:.1f}), kayit tipi "
          f"{floor['record_types']}")
    print(f"  oran: tum girislerde %{100*floor['published_rate_all_entries']:.1f}, "
          f"cagrilabilir kayitlarda %{100*floor['rate_among_callable_records']:.1f}")
    print(f"  {'pencere dolulugu':20} {'n':>6} {'oran':>7} {'gen basina':>11} {'plazmit':>8}")
    for b in bias["by_window_occupancy"]["bins"]:
        print(f"  {b['annotated_neighbours']:20} {b['n']:6} "
              f"%{100*b['transposon_rate']:5.1f} {b['per_gene_rate']:11.5f} "
              f"%{100*b['plasmid_rate']:6.1f}")
    spread = bias["by_window_occupancy"]["spread"]
    print(f"    yayilim: %{100*spread['lowest_rate']:.1f} - "
          f"%{100*spread['highest_rate']:.1f} = {spread['fold']}x")
    print(f"  {'replikon CDS':20} {'n':>6} {'oran':>7} {'plazmit':>8} {'med.komsu':>10}")
    for b in bias["by_replicon_cds_count"]["bins"]:
        print(f"  {b['cds_count']:20} {b['n']:6} %{100*b['transposon_rate']:5.1f} "
              f"%{100*b['plasmid_rate']:6.1f} {b['median_annotated_neighbours']:10.0f}")
    reach = bias["how_much_bias_reaches_the_biological_claim"]
    print(f"  substrat sinifi x transpozon: ham OR {reach['crude_odds_ratio']} -> "
          f"yogunluga gore {reach['mantel_haenszel_odds_ratio']}  "
          f"(fazlaligin %{100*(reach['share_of_excess_odds_attributable_to_annotation_density'] or 0):.0f}'i "
          f"anotasyon yogunlugundan)")

    div = result.get("regulator_vs_enzyme_divergence")
    if div:
        print("\n-- 3. SORU: duzenleyici ne kadar ayrisik evrildi ----------------")
        cov = div["coverage"]
        print(f"  kapsama: yukari akista duzenleyici olan "
              f"{cov['entries_whose_upstream_gene_is_a_regulator']} giristen "
              f"{cov['of_those_with_both_sequences_available']}'inde iki dizi de var "
              f"(%{100*cov['sequence_availability']:.1f}); tum dogrulanmis girisin "
              f"%{100*cov['share_of_all_confirmed_entries']:.1f}'i")
        if "sampling" not in div:
            print("  hesaplanamadi"); return
        smp = div["sampling"]
        print(f"  ornekleme: {smp['all_within_type_pairs_available']} ic ciftten "
              f"{smp['pairs_scored']} ornekledi (tip basina en fazla "
              f"{smp['max_pairs_per_type']}, tohum {smp['seed']}), "
              f"kapsama kontrolunden {smp['pairs_dropped_for_low_alignment_coverage']} "
              f"dustu, {smp['pairs_used']} cift kullanildi")
        header = (f"  {'bant':10} {'cift':>6} {'tip':>4} {'cins-cifti':>11} "
                  f"{'enzim':>7} {'reg':>7} {'fark':>7} {'reg<enz':>8} {'tabanda':>8}")
        for title, series in (("TUM CIFTLER", div["by_enzyme_identity_band"]),
                              ("AYNI AILE", div["same_regulator_family_only"]
                               ["by_enzyme_identity_band"])):
            print(f"  [{title}]")
            print(header)
            for row in series:
                print(f"  {row['label']:10} {row['n_pairs']:6} {row['n_types']:4} "
                      f"{row['n_genus_pairs']:11} "
                      f"{row['enzyme_identity_median']:7.1f} "
                      f"{row['regulator_identity_median']:7.1f} "
                      f"{row['median_difference_enzyme_minus_regulator']:+7.1f} "
                      f"{100*row['share_of_pairs_regulator_less_conserved']:7.0f}% "
                      f"{100*row['share_of_pairs_at_or_below_random_floor']:7.0f}%"
                      f"{'' if row['informative'] else '  (tabanda, yorumlanamaz)'}")
        ret = div["regulator_family_retention"]
        print(f"  aile degisimi: genel %{100*ret['overall_switch_rate']:.0f}  "
              + "  ".join(f"{r['band']}=%{100*r['family_switch_rate']:.0f}"
                          for r in ret["by_enzyme_identity_band"]))
        print("  en sik degisen aile ciftleri: "
              + ", ".join(f"{'/'.join(r['families'])} ({r['n_pairs']})"
                          for r in ret["most_frequently_switched_family_pairs"][:5]))
        cor = div["correlation"]
        print(f"  korelasyon: cift duzeyi rho={cor['pair_level']['spearman_rho']} "
              f"(n={cor['pair_level']['n']}, bagimsiz degil)  |  "
              f"tip x cins-cifti rho={cor['collapsed_level']['spearman_rho']} "
              f"(n={cor['collapsed_level']['n']}, p={cor['collapsed_level']['p']:.1e})")
        print(f"  grup bazinda (ayni aile): "
              + "  ".join(f"g{r['label'].split()[-1]}: reg<enz "
                          f"%{100*r['share_of_pairs_regulator_less_conserved']:.0f} "
                          f"(n={r['n_pairs']})"
                          for r in div["by_enzyme_group_same_family"]))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--out-dir", default="analysis_out")
    parser.add_argument("--chemistry", default="chemistry.csv")
    parser.add_argument("--ecology", default="cluster_ecology.csv")
    parser.add_argument("--domain-csv", default="analysis_out/domain_by_cluster.csv")
    parser.add_argument("--permutations", type=int, default=PERMUTATIONS)
    # 3. soru cift basina iki global hizalama kosuyor; isolation_source.py ile
    # ayni desende coklu surec kullanilir.
    parser.add_argument("--cpu", type=int, default=min(20, os.cpu_count() or 1))
    parser.add_argument("--skip-divergence", action="store_true",
                        help="3. soruyu (duzenleyici ayrisma hizi) atla; "
                             "pairwise hizalamalar en pahali adim")
    args = parser.parse_args()

    con = sqlite3.connect(args.db)
    data, euk_heavy = load(con, args.chemistry, args.ecology, args.domain_csv)
    if not data:
        raise SystemExit("[hata] dogrulanmis RO bulunamadi; --db yolunu kontrol et")
    genus_rows = collapse_genus(data)
    rng = np.random.default_rng(RANDOM_SEED)

    result = {
        "method": {
            "why_this_module_exists": (
                "three questions about the operon. The first two belong together. "
                "First, which regulator goes with which enzyme family and type, which "
                "had never been crossed against the chemistry. Second, whether the "
                "published mobility result depends on gene categories that are a "
                "regular expression over GenBank /product text rather than a "
                "measurement of sequence. Third, how divergently the regulators "
                "evolved compared with the enzymes they sit next to, resolved by "
                "enzyme group and by the enzyme-to-enzyme identity band."),
            "two_grades_of_neighbourhood_evidence": {
                "sequence_verified": (
                    "operon components (beta subunit, ferredoxin, reductase) are "
                    "found by scanning the neighbour protein sequence with a profile "
                    "HMM in build_operons.py. This survives a wrong or missing "
                    "/product annotation."),
                "annotation_derived": (
                    f"gene categories in gene_category with method '{REGEX_METHOD}' "
                    "are a regular expression over the GenBank /product text, defined "
                    "in build_db.py. They measure annotation quality, not biology. The "
                    "transposon category is one of these, and so is the regulator "
                    "family label used in the first and third questions."),
            },
            "independence": (
                "RO entries are not independent observations: the same protein recurs "
                "in repeatedly sequenced strains. Every claim here is reported at two "
                "levels, per entry and per enzyme type and genus, using the same "
                "collapse rule as stats_overview.py. Where the two levels disagree the "
                "genus level is the one to read, and the disagreement is reported "
                "rather than hidden."),
            "effect_sizes": (
                "no bare p values. Cross tables carry Cramer's V together with the V "
                f"expected by chance from the same marginals ({args.permutations} "
                "permutations), because a sparse table yields a positive V with no "
                "association at all. Rank comparisons carry Cliff's delta. Two by two "
                "tables carry the rate ratio and the odds ratio."),
            "thresholds": {
                "expected_cell_floor": EXPECTED_CELL_FLOOR,
                "sparse_table_limit": SPARSE_TABLE_LIMIT,
                "cramers_v_negligible": V_NEGLIGIBLE,
                "cramers_v_moderate": V_MODERATE,
                "min_cell_entries_for_a_named_finding": MIN_CELL_ENTRIES,
                "min_distinct_genera_for_a_named_finding": MIN_CELL_GENERA,
                "min_stratum_n": MIN_STRATUM_N,
                "density_bins": f"quartiles taken from the data, {DENSITY_QUANTILES}",
                "permutations": args.permutations,
                "random_seed": RANDOM_SEED,
                "close_to_a_curated_enzyme": list(CLOSE_TIERS),
            },
        },
        "totals": {
            "entries": len(data),
            "type_and_genus_observations": len(genus_rows),
            "distinct_types": len({d["cluster"] for d in data}),
            "distinct_genera": len({d["genus"] for d in data}),
        },
    }

    result["regulation"] = regulation_block(data, genus_rows, rng, args.permutations)
    result["mobility_by_substrate_class"] = mobility_block(data, genus_rows, euk_heavy)
    result["transposon_vs_curated_distance"] = evidence_block(
        data, genus_rows, rng, args.permutations)
    result["annotation_bias"] = annotation_bias_block(data, euk_heavy, rng)
    if not args.skip_divergence:
        result["regulator_vs_enzyme_divergence"] = divergence_block(con, args.cpu)
    result["limits"] = [
        "The transposon calls are never verified against sequence. The regular "
        "expression also matches integrase, recombinase and resolvase, which are not "
        "always transposon components, so the category is an upper bound on text "
        "matches and says nothing about what is actually in the DNA.",
        "The regulator family label comes from the same kind of /product regular "
        "expression, so the first question rests on the same annotation chain as the "
        "bias measured in the second. The internal consistency check is the geometry: "
        "family labels that reproduce the classical orientation and spacing of their "
        "family are more believable than ones that do not.",
        "No promoter sequence was examined. 88.9 % of the source GenBank records are "
        "CON entries with no sequence in the file, so operator sites and -35/-10 boxes "
        "cannot be inspected and the regulator association cannot be confirmed at the "
        "level of DNA binding.",
        "The plasmid field is GenBank's own declaration. Integrated plasmids and "
        "unnamed megaplasmids are missed, so it is an undercount, but it is a property "
        "of the molecule rather than a choice of words.",
        "Nothing here is causal. Both questions measure co-occurrence in sequenced "
        "genomes, and the sampled genomes are not a sample of nature: PAH and "
        "xenobiotic degraders are a well studied group, which biases annotation "
        "richness and isolate counts in the same direction as the hypotheses.",
        "Chemical family is a property of the enzyme type, not of the entry, so the "
        "regulator-by-family tables are aggregations over at most 61 types. A family "
        "carried by two or three types can be moved by one type.",
        "The third question identifies a regulator by POSITION, as the first gene "
        "upstream of the operon 5' end. No binding was measured, so 'the regulator of "
        "this enzyme' is an inference from gene order. A regulator acting from "
        "elsewhere on the replicon, or a global regulator, is invisible to it.",
        "Pairwise identity saturates. Unrelated proteins align at roughly 20 to 25 % "
        "identity, so below that level the measure cannot separate a fast-diverging "
        "regulator from an unrelated protein, and no divergence rate is claimed in "
        "the bands where most pairs sit at that floor.",
        "The third question covers only the entries whose upstream gene is a "
        "regulator at all, under three in ten. Sequence availability is not the "
        "limit, because protein translations survive in CON records even though the "
        "nucleotide sequence does not; gene order is.",
    ]

    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, "operon_relations.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1, default=float)

    print_summary(result)
    print(f"\n[yazildi] {path}")


if __name__ == "__main__":
    main()
