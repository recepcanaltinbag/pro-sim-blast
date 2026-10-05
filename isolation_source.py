"""
Izolasyon kaynagi ve ekolojik meta veri -- gbk'da duran ama hic kullanilmayan alan.

NEDEN: veritabanindaki her RO bir replikona, her replikon bir GenBank kaydina
bagli ve o kayitlarin `source` ozelliginde izolasyonun NEREDEN yapildigi yaziyor
(17.073 kaydin 9.591'inde /isolation_source, 4.911'inde /host, 11.485'inde
/geo_loc_name veya /country). Pipeline bugune kadar yalnizca organizma ve
taksonomiyi okuyordu; "ksenobiyotik oksijenazlar kirli bolgelerden gelir"
turunden ekolojik iddialar ise ancak bu alanla sinanabilir.

NE OLCULMUYOR -- bunu bastan yazmak gerekiyor:
  * /isolation_source SERBEST METINDIR. Kontrollu bir sozluk degil; ayni habitat
    onlarca farkli yazimla gecer ve bir kismi hic bilgi tasimaz ("environmental",
    "not provided"). Bu script metni bir sozluge esler, OLCMEZ.
  * Esleme kelime temellidir (asagidaki HABITAT_RULES). Siralidir: ilk eslesen
    kural kazanir. Dolayisiyla her atama tek bir anahtar kelimeye indirgenebilir
    ve denetlenebilir -- cikti JSON'unda hangi kelimenin kac kez karar verdigi
    de yaziyor (`habitat_keyword_hits`).
  * Guvenle siniflanamayan metin `other`, metin olmayan/bos olan `unknown`
    kalir. Kapsama ISTEYEREK sismiyor: `host` alani habitat kararina
    KARISTIRILMIYOR, cunku konak adi orneklenen canliyi soyler, ortami soylemez
    (homo sapiens bir habitat degil; akciger de, dista bir yara da olabilir).
  * Orneklem yanliligi duzeltilmiyor, SADECE gorunur kiliniyor. Giris sayilari
    cok dizilenmis suslarin (Pseudomonas, Mycobacterium) tekrarindan olusur;
    bu yuzden her capraz tablo hem GIRIS hem de FARKLI TUR (organizma ikili
    adinin ilk iki kelimesi) sayisiyla veriliyor. "Ilk iki kelime" bir tur
    tanimi degildir: "Pseudomonas sp." tek bir sozde-tur olarak sayilir.
  * Cografya ulke duzeyindedir ve orneklenen YERI degil, YAYINI yapan grubun
    bildirdigi yeri gosterir. Koordinat (/lat_lon) okunmuyor.
  * Toplama yili yalnizca yildir; /collection_date'teki ay/gun atiliyor.
    Alan araliklar da iceriyor ("1900/1969"); ilk yil alinir, yani 1960
    oncesi yillarin bir kismi gercek izolasyon tarihi degil aralik basidir.

KARAR KURALLARI (siralama bilincli, tartisilabilir yerleri isaretliyorum):
  * `contaminated_industrial` kirlilik bilgisini FIZIKSEL ortam bilgisinin
    onune koyar: "PCB-contaminated lagoon sludge" camur degil kirli saha
    sayilir. Bu veritabaninin sorusu kirlilik oldugu icin bilincli bir secim.
  * `sediment` deniz/tatli sudan ONCE gelir: "marine sediment" sediment olur.
  * `rhizosphere_plant` topraktan ONCE gelir: "rhizosphere soil" bitki olur,
    "forest soil" toprak kalir.
  * `extreme_thermal` kirlilikten SONRA ama ortamlardan ONCE gelir: "fly ash
    dumping site of thermal power" sanayi, "decaying wood from a thermal pond"
    termal olur.
  * Olumsuzlama: "non-axenic culture from freshwater" gibi dizgilerde "non"
    one eki atlanir (bkz. `_hit`).
  * Cok yayginligina ragmen ESLENMEYEN metinler var: ciplak "water" (ne tatli
    ne tuzlu oldugu belirsiz), "cave", "hypersaline"/"saltern", "desert sand",
    "bioreactor", konak tur adlari. Bunlar `other` icinde kaliyor ve ciktida
    tek tek listeleniyor; yanlis siniflamaktansa siniflamamak yegdir.

Adimlar:
    1. En az bir dogrulanmis RO tasiyan replikonlarin gbk'larini coklu surecle
       (build_operons.py ile ayni desen) yeniden parse et, `source` niteliklerini
       cek -> replicon_source
    2. isolation_source -> habitat (HABITAT_RULES, sirali)
    3. Habitat x substrat sinifi / kimyasal aile capraz tablolari, tip basina
       habitat profili, ulke/yil dagilimlari, kapsama blogu
       -> analysis_out/habitat.json

Kullanim:
    python3 isolation_source.py --db roar.sqlite --gbk-dir gbk_files \
        --out analysis_out/habitat.json --cpu 20

Tekrar calistirmak guvenli: replicon_source doluysa gbk parse atlanir ama
habitat atamasi her kosuda yeniden yapilir (kural degisikligi hemen yansir).
`--force` gbk parse'i da zorlar.
"""

import argparse
import csv
import json
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

from Bio import SeqIO

SCHEMA = """
CREATE TABLE IF NOT EXISTS replicon_source (
    nucleotide_id    TEXT PRIMARY KEY,
    isolation_source TEXT,
    host             TEXT,
    geo              TEXT,     -- /geo_loc_name veya /country, oldugu gibi
    country          TEXT,     -- geo'nun ilk parcasi (":" oncesi)
    collection_year  INTEGER,
    habitat          TEXT      -- HABITAT_RULES sonucu
);
CREATE INDEX IF NOT EXISTS idx_repsrc_habitat ON replicon_source(habitat);
"""

# ---------------------------------------------------------------- sozluk
# Metin once normalize edilir: kucuk harf, [a-z0-9] disindaki her sey bosluk.
# Yani "PCB-contaminated" -> "pcb contaminated", "ENVO:00002170" -> "envo 00002170".
# Anahtar kelime tam kelime olarak aranir, sonunda istege bagli "s"/"es" kabul
# edilir (soil/soils, nodule/nodules). Bu yuzden "mine" "mineral"i tutmaz,
# "oil" "soil"i tutmaz, "coral" "coralloid"u tutmaz.
#
# SIRA ONEMLIDIR: listedeki ilk eslesen kural kazanir.
HABITAT_RULES = (
    # 1. Once laboratuvar: ortam degil, kokeni belirsiz kultur.
    ("laboratory_strain", (
        "laboratory", "laboratory strain", "lab strain", "axenic culture",
        "pure culture", "genetically modified", "in vitro",
    )),
    # 2. Kirlilik/sanayi: fiziksel ortamdan ONCE (bkz. modul docstring).
    ("contaminated_industrial", (
        "contaminated", "contamination", "polluted", "pollution", "xenobiotic",
        "crude oil", "oil sludge", "oil field", "oil well", "oil shale",
        "fuel oil", "waste oil", "oil refinery", "refinery", "petroleum",
        "petrochemical", "diesel", "gasoline", "kerosene", "hydrocarbon",
        "tar", "coal tar", "creosote", "coke",
        "pcb", "pah", "hch", "hexachlorocyclohexane", "pentachlorophenol",
        "chlorophenol", "nitrophenol", "trichloroethylene", "dioxin",
        "phthalate", "benzene", "toluene", "naphthalene", "biphenyl",
        "atrazine", "isoproturon", "pesticide", "herbicide", "insecticide",
        "tnt", "explosive",
        "landfill", "slag", "smelter", "fly ash",
        "industrial", "pulp mill", "paper mill", "kraft", "tannery",
        "gas works", "power plant", "bioremediation",
        "phenol", "cresol", "aniline", "solvent",
    )),
    # 2b. Madencilik AYRI kategori. Eskiden contaminated_industrial icindeydi
    #     ama asit maden drenaji ve metal atigi ORGANIK kirlilik degildir; bu
    #     veritabaninin sorusu organik substrat yikimi oldugu icin ikisini ayni
    #     kovaya atmak ana bulguyu suluyordu. OLCULDU: maden girisleri
    #     cikarilinca PAH enzimlerinin kirli-saha zenginlesmesi 2,97x'ten
    #     3,65x'e cikti, yani madenler sinyali azaltiyordu.
    ("mining_acid_drainage", (
        "mine", "mining", "tailing", "ore deposit",
    )),
    # 3. Termal/volkanik imza: kirlilikten SONRA, cunku "thermal power plant"
    #    sanayi tesisidir, sicak kaynak degil ("fly ash ... thermal power").
    ("extreme_thermal", (
        "hot spring", "hot water", "thermal", "geothermal", "hydrothermal",
        "fumarole", "solfataric", "solfatara", "volcanic", "volcano", "geyser",
        "steam", "boiling",
    )),
    # 4. Muhendislik urunu su/camur sistemleri.
    ("wastewater_sludge", (
        "wastewater", "waste water", "sewage", "sewer", "activated sludge",
        "sludge", "effluent", "treatment plant", "wwtp", "digester",
        "anaerobic digestion", "septic", "waste treatment",
    )),
    # 5. Insan klinigi. Bilincli olarak DISARIDA tutulanlar: "tumor", "lesion",
    #    "tissue", "skin", "hip", "ear" -- bitki/hayvan orneklerinde de geciyor.
    ("human_clinical", (
        "sputum", "sputa", "blood", "urine", "urinary", "wound", "pus",
        "abscess", "patient", "clinical", "nosocomial", "hospital",
        "intensive care", "health care", "icu", "nicu", "catheter",
        "bronchoscopy", "bronchial", "bronchus", "trachea", "tracheal",
        "endotracheal", "lung", "pulmonary", "respiratory", "throat",
        "pharyngeal", "nasopharyngeal", "nasal", "sinus", "antral",
        "conjunctiva", "cornea", "cerebrospinal", "cerebral", "brain",
        "lymph node", "spleen", "biopsy", "vagina", "vaginal", "cervical",
        "perirectal", "rectal swab", "oral", "surgical", "dialysis",
        "peritoneal", "aids", "cystic fibrosis", "mycobacteriosis",
        "granuloma", "granulomatous", "ulcer", "celiac", "duodenal",
        "gastric", "mucosa", "sepsis", "septicemia", "bacteremia",
        "expectoration", "bodily fluid", "homo sapiens", "human",
    )),
    # 6. Sindirim sistemi / diski -- insan ya da hayvan olabilir, ayirmiyoruz.
    ("gut_faecal", (
        "feces", "faeces", "fecal", "faecal", "stool", "gut", "intestine",
        "intestinal", "caecum", "cecum", "cecal", "colon", "colonic",
        "rumen", "ruminal", "stomach", "manure", "dung", "excrement",
        "digestive tract",
    )),
    # 7. Gida / fermantasyon.
    ("food_fermented", (
        "milk", "cheese", "yogurt", "yoghurt", "butter", "cream", "curd",
        "whey", "meat", "beef", "pork", "ham", "sausage", "salami",
        "chicken breast", "retail", "seafood", "food", "beverage", "juice",
        "beer", "wine", "vinegar", "cider", "sake", "brewery", "winery",
        "kombucha", "kimchi", "doenjang", "miso", "soy sauce",
        "fermented", "fermentation", "pickle", "dough", "bread",
    )),
    # 8. Hayvan konagi (sindirim disi dokular ve simbiyozlar).
    ("animal_host", (
        "nidamental gland", "squid", "oyster", "mussel", "clam", "snail",
        "sponge", "coral", "fish", "seahorse", "crab",
        "insect", "beetle", "mosquito", "termite", "wasp", "bee", "honeybee",
        "aphid", "larva", "nematode", "tick", "louse", "silkworm",
        "earthworm", "egg", "nest", "mouse", "mice", "rat", "pig", "porcine",
        "cattle", "bovine", "cow", "sheep", "goat", "chicken", "poultry",
        "bird", "dog", "horse", "animal", "drosophila", "caenorhabditis",
        "bos taurus", "mus musculus", "sea cucumber", "sea urchin",
    )),
    # 9. Bitki ile iliskili her sey (rizosfer, nodul, filosfer, bitki dokusu).
    #    Topraktan ONCE: "rhizosphere soil" bitki, "forest soil" toprak.
    ("rhizosphere_plant", (
        "rhizosphere", "rhizospheric", "rhizoplane", "phyllosphere",
        "endophyte", "endophytic", "root", "nodule", "leaf", "leaves",
        "stem", "shoot", "seedling", "seed", "grain", "tuber", "flower",
        "pollen", "bark", "wood", "trunk", "gall", "phloem", "plant",
        "moss", "fruit", "apple", "grape", "grapevine", "peach",
        "strawberry", "tomato", "potato", "rice", "wheat", "maize", "corn",
        "soybean", "bean", "pea", "vegetable", "watercress", "straw",
        "leaf litter",
    )),
    # 10. Sediment, deniz/tatli sudan ONCE: "marine sediment" sediment olur.
    # Hipersalin: tatli su/deniz kurallarindan ONCE, cunku "hypersaline lake"
    # icinde "lake" gecer ve aksi halde tatli su sayilir (olculdu, duzeltildi).
    ("hypersaline", ("hypersaline", "saltern", "salt pan", "salt lake",
                     "brine", "solar salt")),
    ("sediment", ("sediment", "silt", "mud", "mudflat", "benthic")),
    ("marine", (
        "marine", "sea", "seawater", "sea water", "sea ice", "ocean",
        "oceanic", "coastal", "coast", "seashore", "intertidal",
        "tidal", "tide", "estuarine", "estuary", "brackish", "mangrove",
        "pelagic", "reef", "salt marsh", "atlantic", "pacific",
    )),
    ("freshwater", (
        "freshwater", "fresh water", "lake", "river", "pond", "stream",
        "creek", "brook", "canal", "reservoir", "groundwater", "ground water",
        "drinking water", "tap water", "well water", "mineral water",
        "spring water", "rain water", "wetland", "marsh", "waterfall",
        "aquifer", "hyporheic",
    )),
    ("soil", (
        "soil", "topsoil", "subsoil", "compost", "humus", "peat", "paddy",
        "farmland", "cropland", "pasture", "loam",
    )),
    ("air_dust", ("air", "aerosol", "dust", "atmospheric")),
)

# --- Sonradan eklenen kategoriler (ilk surumde `other` icinde kalmislardi) ---
# Listenin SONUNDA: daha belirli kurallar once calisir, yani "marine water"
# denizde kalir, yalnizca ciplak "water" asagiya duser.
HABITAT_RULES = HABITAT_RULES + (
    # Kompartmani belirsiz su: en buyuk siniflanmayan kume (161 giris).
    # "Su oldugunu biliyoruz, hangisi oldugunu bilmiyoruz" ile "hicbir sey
    # bilmiyoruz" ayri seylerdir; ikincisi `unknown`.
    ("water_unspecified", ("water", "aqueous")),
    # Muhendislik sistemleri: kirletici adi geciyorsa yukaridaki kirlilik
    # kurali zaten yakalar; buraya ortami belirsiz olanlar duser.
    ("engineered_system", ("bioreactor", "reactor", "enrichment culture",
                           "chemostat", "biofilter", "batch culture")),
    ("cave", ("cave", "speleothem", "karst")),
    ("subsurface_deep", ("subsurface", "borehole", "aquifer", "groundwater")),
    ("fungal_associated", ("mycosphere", "sporocarp", "fruiting body",
                           "mushroom")),
)

HABITAT_ORDER = tuple(name for name, _ in HABITAT_RULES) + ("other", "unknown")
CLASSIFIED = set(name for name, _ in HABITAT_RULES)

# Metin var ama bilgi tasimiyor -> unknown (other degil: "siniflanamadi" ile
# "soylenmemis" ayri seylerdir ve kapsama blogunda ayri sayiliyorlar).
UNINFORMATIVE_EXACT = {
    "", "na", "n a", "none", "nd", "unknown", "missing", "not known",
    "not available", "not provided", "not applicable", "not determined",
    "not collected", "unspecified", "environmental", "environment",
    "sample", "isolate", "specimen",
}
UNINFORMATIVE_PREFIX = (
    "not provided", "not available", "not applicable", "not known",
    "not collected", "not determined", "not specified", "missing",
    "unknown", "unspecified",
)

SOURCE_QUALIFIERS = ("isolation_source", "host", "geo_loc_name", "country",
                     "collection_date")
YEAR_RE = re.compile(r"\b(1[89]\d{2}|20\d{2})\b")
MIN_YEAR, MAX_YEAR = 1850, 2030


def normalise(value):
    """Serbest metni kelime temelli eslemeye hazirla."""
    if not value:
        return ""
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def _keyword_regex(keyword):
    return re.compile(r"\b" + re.escape(keyword) + r"(e?s)?\b")


_COMPILED_RULES = tuple(
    (habitat, tuple((kw, _keyword_regex(kw)) for kw in keywords))
    for habitat, keywords in HABITAT_RULES
)


def _hit(regex, text):
    """Kelime eslesmesi, "non ..." olumsuzlamasi atlanarak.

    Veride "non-axenic culture from freshwater", "mangroves; non-axenic
    culture" gibi dizgiler var; normalizasyon tireyi bosluga cevirdigi icin
    duz arama bunlari yanlis tarafa koyar. Bitisik olumsuzlamalar
    ("uncontaminated") kelime siniri sayesinde zaten eslesmiyor.
    """
    for match in regex.finditer(text):
        if not text[:match.start()].endswith("non "):
            return True
    return False


def classify_habitat(isolation_source):
    """isolation_source -> (habitat, karari veren anahtar kelime)."""
    text = normalise(isolation_source)
    if not text:
        return "unknown", None
    if text in UNINFORMATIVE_EXACT or text.startswith(UNINFORMATIVE_PREFIX):
        return "unknown", "<uninformative>"
    for habitat, keywords in _COMPILED_RULES:
        for keyword, regex in keywords:
            if _hit(regex, text):
                return habitat, keyword
    return "other", None


def parse_year(collection_date):
    """/collection_date -> yil (int) ya da None. Ay/gun atilir."""
    if not collection_date:
        return None
    match = YEAR_RE.search(collection_date)
    if not match:
        return None
    year = int(match.group(1))
    return year if MIN_YEAR <= year <= MAX_YEAR else None


def split_country(geo):
    """"Germany: Cologne" -> "Germany". Ulke adi oldugu gibi birakilir."""
    if not geo:
        return None
    return geo.split(":")[0].strip() or None


# ------------------------------------------------------- 1. gbk source parse
def _source_worker(job):
    """(path, wanted_ids) -> [(nucleotide_id, iso, host, geo, year, fallback)]

    Yalnizca `source` ozelligi okunuyor ama SeqIO kaydin tamamini parse eder;
    bu kabul edilen maliyet (build_operons.py ile ayni desen, ~9.300 dosya).
    """
    path, wanted = job
    out = []
    try:
        with open(path) as handle:
            records = list(SeqIO.parse(handle, "genbank"))
    except Exception as exc:          # tek dosya tum isi durdurmasin
        sys.stderr.write(f"[hata] {path}: {exc}\n")
        return out

    for record in records:
        source = None
        for feature in record.features:
            if feature.type == "source":
                source = feature.qualifiers
                break
        values = {}
        for key in SOURCE_QUALIFIERS:
            raw = (source or {}).get(key, [""])[0].strip()
            values[key] = raw or None
        # geo_loc_name yeni, country eski ad; ikisi ayni alan.
        geo = values["geo_loc_name"] or values["country"]
        row = (values["isolation_source"], values["host"], geo,
               parse_year(values["collection_date"]))
        if record.id in wanted:
            out.append((record.id,) + row + (0,))
        elif len(records) == 1 and len(wanted) == 1:
            # Surum/ad uyusmazligi: dosyada tek kayit ve tek beklenen id varsa
            # eslestir, ama bunu ayrica say (kapsama blogunda raporlanir).
            out.append((next(iter(wanted)),) + row + (1,))
    return out


def parse_sources(connection, gbk_dir, processes, force):
    existing = connection.execute("SELECT COUNT(*) FROM replicon_source").fetchone()[0]
    if existing and not force:
        print(f"[atlandi] replicon_source dolu ({existing} replikon)")
        return None        # parse edilmedi; kapsama blogunda null olarak gecer

    wanted = defaultdict(set)
    files = {}
    for nuc, path in connection.execute("""
            SELECT DISTINCT rep.nucleotide_id, rep.file
            FROM replicon rep JOIN ro r ON r.nucleotide_id = rep.nucleotide_id
            WHERE r.is_confirmed = 1"""):
        full = os.path.join(gbk_dir, path or f"{nuc}.gbk")
        wanted[full].add(nuc)
        files[nuc] = full

    jobs, missing = [], []
    for path, ids in wanted.items():
        if os.path.exists(path):
            jobs.append((path, ids))
        else:
            missing.extend(ids)
            sys.stderr.write(f"[uyari] gbk yok: {path}\n")
    print(f"[replicon_source] {len(jobs)} gbk dosyasi, {processes} surec")

    connection.execute("DELETE FROM replicon_source")
    seen, fallbacks = set(), 0
    with Pool(processes) as pool:
        for i, rows in enumerate(pool.imap_unordered(_source_worker, jobs,
                                                     chunksize=8), 1):
            for nuc, iso, host, geo, year, fallback in rows:
                connection.execute(
                    "INSERT OR REPLACE INTO replicon_source VALUES (?,?,?,?,?,?,NULL)",
                    (nuc, iso, host, geo, split_country(geo), year))
                seen.add(nuc)
                fallbacks += fallback
            if i % 500 == 0:
                connection.commit()
                sys.stdout.write(f"\r  {i}/{len(jobs)} dosya | {len(seen)} replikon")
                sys.stdout.flush()
    # Parse edilemeyen/bulunamayan replikonlar da satir alir: payda acik kalsin.
    for nuc in files:
        if nuc not in seen:
            connection.execute(
                "INSERT OR REPLACE INTO replicon_source VALUES (?,NULL,NULL,NULL,NULL,NULL,NULL)",
                (nuc,))
    connection.commit()
    print(f"\r[replicon_source] {len(seen)} replikon okundu, "
          f"{len(files) - len(seen)} bos satir (dosya yok / parse hatasi)"
          + " " * 10)
    if fallbacks:
        print(f"   [not] {fallbacks} kayitta record.id beklenen id ile ayni degildi, "
              f"tek kayitli dosya oldugu icin eslestirildi")
    return fallbacks


# -------------------------------------------------------- 2. habitat atamasi
def assign_habitats(connection):
    """Depolanan metinden habitat'i her kosuda yeniden hesapla."""
    rows = connection.execute(
        "SELECT nucleotide_id, isolation_source FROM replicon_source").fetchall()
    hits = defaultdict(Counter)
    updates = []
    for nuc, iso in rows:
        habitat, keyword = classify_habitat(iso)
        hits[habitat][keyword or "<no rule matched>"] += 1
        updates.append((habitat, nuc))
    connection.executemany(
        "UPDATE replicon_source SET habitat = ? WHERE nucleotide_id = ?", updates)
    connection.commit()
    print(f"[habitat] {len(updates)} replikon siniflandirildi")
    return {h: dict(c.most_common()) for h, c in hits.items()}


# ------------------------------------------------------------- 3. toplamalar
def read_csv_map(path, key="cluster"):
    if not os.path.exists(path):
        sys.stderr.write(f"[uyari] {path} yok, ilgili sutunlar bos kalir\n")
        return {}
    with open(path) as handle:
        return {row[key]: row for row in csv.DictReader(handle)}


def species_of(organism):
    """Organizma ikili adinin ilk iki kelimesi. Tur tanimi DEGIL (bkz. docstring)."""
    if not organism:
        return None
    parts = organism.split()
    return " ".join(parts[:2]) if parts else None


def load_entries(connection, ecology, chemistry):
    """Her dogrulanmis RO icin tek satir: tip, tur, habitat, kimya, cografya."""
    rows = connection.execute("""
        SELECT r.candidate_id, r.nucleotide_id, r.ro_cluster, rep.organism,
               s.isolation_source, s.host, s.country, s.collection_year, s.habitat
        FROM ro r
        JOIN replicon rep ON rep.nucleotide_id = r.nucleotide_id
        LEFT JOIN replicon_source s ON s.nucleotide_id = r.nucleotide_id
        WHERE r.is_confirmed = 1""").fetchall()
    entries = []
    for (cid, nuc, cluster, organism, iso, host, country, year, habitat) in rows:
        eco = ecology.get(cluster, {})
        chem = chemistry.get(cluster, {})
        entries.append({
            "candidate_id": cid,
            "nucleotide_id": nuc,
            "cluster": cluster if cluster and cluster != "N/A" else None,
            "species": species_of(organism),
            "isolation_source": iso,
            "host": host,
            "country": country,
            "year": year,
            "habitat": habitat or "unknown",
            "substrate_class": eco.get("substrate_class") or "unknown",
            "family": chem.get("family") or "unknown",
        })
    return entries


def crosstab(entries, axis):
    """axis x habitat: hem giris hem farkli tur sayisi."""
    counts = defaultdict(Counter)
    species = defaultdict(lambda: defaultdict(set))
    for entry in entries:
        key, habitat = entry[axis], entry["habitat"]
        counts[key][habitat] += 1
        if entry["species"]:
            species[key][habitat].add(entry["species"])
    return {
        key: {
            "entries": {h: counts[key][h] for h in HABITAT_ORDER if counts[key][h]},
            "species": {h: len(species[key][h]) for h in HABITAT_ORDER
                        if species[key].get(h)},
            "entries_total": sum(counts[key].values()),
            "species_total": len(set().union(*species[key].values()))
                             if species[key] else 0,
        }
        for key in counts
    }


def species_normalised_profile(entries, axis):
    """Isin asil noktasi: ayni dagilim GIRIS payiyla ve TUR payiyla yan yana.

    Payda yalnizca SINIFLANMIS habitatlar (other/unknown disarida) -- boylece
    iki pay da 1'e toplanir ve karsilastirilabilir. Bir tur iki habitatta
    gorunuyorsa ikisinde de sayilir; bu yuzden tur paydasi "tur x habitat
    cifti" sayisidir, `species_distinct` ise ayni turu bir kez sayar.
    """
    usable = [e for e in entries
              if e["habitat"] in CLASSIFIED and e["species"] and e[axis]]
    counts = defaultdict(Counter)
    species = defaultdict(lambda: defaultdict(set))
    for entry in usable:
        counts[entry[axis]][entry["habitat"]] += 1
        species[entry[axis]][entry["habitat"]].add(entry["species"])

    out = {}
    for key in counts:
        entry_total = sum(counts[key].values())
        pair_total = sum(len(s) for s in species[key].values())
        habitats = {}
        for habitat in HABITAT_ORDER:
            n_entries = counts[key].get(habitat, 0)
            n_species = len(species[key].get(habitat, ()))
            if not n_entries and not n_species:
                continue
            entry_fraction = n_entries / entry_total if entry_total else 0.0
            species_fraction = n_species / pair_total if pair_total else 0.0
            habitats[habitat] = {
                "entries": n_entries,
                "entry_fraction": round(entry_fraction, 4),
                "species": n_species,
                "species_fraction": round(species_fraction, 4),
                # pozitif = giris sayilari bu habitati tur duzeyinden FAZLA
                # gosteriyor (cok dizilenmis suslar)
                "bias_delta": round(entry_fraction - species_fraction, 4),
            }
        out[key] = {
            "entries_total": entry_total,
            "species_habitat_pairs": pair_total,
            "species_distinct": len(set().union(*species[key].values()))
                                if species[key] else 0,
            "habitats": habitats,
        }
    return out


def per_type_profile(entries, ecology, chemistry, top_n=5):
    by_cluster = defaultdict(list)
    for entry in entries:
        if entry["cluster"]:
            by_cluster[entry["cluster"]].append(entry)
    out = {}
    for cluster, rows in by_cluster.items():
        species = defaultdict(set)
        counts = Counter()
        for entry in rows:
            if entry["habitat"] in CLASSIFIED and entry["species"]:
                species[entry["habitat"]].add(entry["species"])
            counts[entry["habitat"]] += 1
        ranked = sorted(((h, len(s)) for h, s in species.items()),
                        key=lambda kv: (-kv[1], kv[0]))
        out[cluster] = {
            "entries": len(rows),
            "species_distinct": len({e["species"] for e in rows if e["species"]}),
            "habitats_distinct": len(species),
            "substrate_class": ecology.get(cluster, {}).get("substrate_class", "unknown"),
            "family": chemistry.get(cluster, {}).get("family", "unknown"),
            "entries_by_habitat": {h: counts[h] for h in HABITAT_ORDER if counts[h]},
            "top_habitats_by_species": [
                {"habitat": h, "species": n, "entries": counts[h]}
                for h, n in ranked[:top_n]],
        }
    return out


def coverage_block(connection, entries, fallbacks):
    total_replicons = connection.execute("""
        SELECT COUNT(DISTINCT rep.nucleotide_id)
        FROM replicon rep JOIN ro r ON r.nucleotide_id = rep.nucleotide_id
        WHERE r.is_confirmed = 1""").fetchone()[0]
    rows = connection.execute("""
        SELECT isolation_source, host, geo, country, collection_year, habitat
        FROM replicon_source""").fetchall()

    missing_iso = uninformative = 0
    for iso, _h, _g, _c, _y, _hab in rows:
        text = normalise(iso)
        if not text:
            missing_iso += 1
        elif text in UNINFORMATIVE_EXACT or text.startswith(UNINFORMATIVE_PREFIX):
            uninformative += 1

    unclassified = Counter()
    for iso, _h, _g, _c, _y, habitat in rows:
        if habitat == "other" and iso:
            unclassified[iso.strip().lower()] += 1

    habitat_of = {}
    for nuc, habitat, host in connection.execute(
            "SELECT nucleotide_id, habitat, host FROM replicon_source"):
        habitat_of[nuc] = (habitat, host)

    ro_iso = sum(1 for e in entries if e["isolation_source"])
    ro_classified = sum(1 for e in entries if e["habitat"] in CLASSIFIED)
    unknown_with_host = sum(
        1 for _nuc, (habitat, host) in habitat_of.items()
        if habitat in ("other", "unknown") and host)

    return {
        "replicons_with_confirmed_ro": total_replicons,
        "replicons_in_table": len(rows),
        "replicons_with_isolation_source": sum(1 for r in rows if r[0]),
        "replicons_with_host": sum(1 for r in rows if r[1]),
        "replicons_with_geo": sum(1 for r in rows if r[2]),
        "replicons_with_country": sum(1 for r in rows if r[3]),
        "replicons_with_collection_year": sum(1 for r in rows if r[4]),
        "replicons_with_classified_habitat":
            sum(1 for r in rows if r[5] in CLASSIFIED),
        "isolation_source_absent": missing_iso,
        "isolation_source_uninformative_text": uninformative,
        "confirmed_ros_total": len(entries),
        "confirmed_ros_with_isolation_source": ro_iso,
        "confirmed_ros_with_isolation_source_fraction":
            round(ro_iso / len(entries), 4) if entries else None,
        "confirmed_ros_with_classified_habitat": ro_classified,
        "confirmed_ros_with_classified_habitat_fraction":
            round(ro_classified / len(entries), 4) if entries else None,
        "confirmed_ros_habitat_other":
            sum(1 for e in entries if e["habitat"] == "other"),
        "confirmed_ros_habitat_unknown":
            sum(1 for e in entries if e["habitat"] == "unknown"),
        "unclassified_strings_distinct": len(unclassified),
        "unclassified_strings_entries": sum(unclassified.values()),
        "unclassified_top_50": [{"isolation_source": s, "replicons": n}
                                for s, n in unclassified.most_common(50)],
        "replicons_unclassified_but_host_present": unknown_with_host,
        "record_id_fallback_matches": fallbacks,
        "note": ("Every figure in this file must be read against "
                 "confirmed_ros_with_classified_habitat: isolation_source is a "
                 "free-text GenBank qualifier, it is absent for a large part of "
                 "the records, and the /host qualifier is deliberately NOT used "
                 "to assign a habitat. Entry counts are dominated by repeatedly "
                 "sequenced strains, so the distinct-species counts next to them "
                 "are the comparable numbers."),
    }


def habitat_summary(entries):
    counts = Counter()
    species = defaultdict(set)
    replicons = defaultdict(set)
    for entry in entries:
        counts[entry["habitat"]] += 1
        replicons[entry["habitat"]].add(entry["nucleotide_id"])
        if entry["species"]:
            species[entry["habitat"]].add(entry["species"])
    total = sum(counts.values())
    return {
        habitat: {
            "entries": counts[habitat],
            "entry_fraction": round(counts[habitat] / total, 4) if total else 0.0,
            "replicons": len(replicons[habitat]),
            "species": len(species[habitat]),
        }
        for habitat in HABITAT_ORDER if counts[habitat]
    }


def geography_block(entries):
    countries = defaultdict(lambda: {"entries": 0, "replicons": set(),
                                     "species": set()})
    years = defaultdict(lambda: {"entries": 0, "replicons": set()})
    for entry in entries:
        if entry["country"]:
            slot = countries[entry["country"]]
            slot["entries"] += 1
            slot["replicons"].add(entry["nucleotide_id"])
            if entry["species"]:
                slot["species"].add(entry["species"])
        if entry["year"]:
            slot = years[entry["year"]]
            slot["entries"] += 1
            slot["replicons"].add(entry["nucleotide_id"])
    country_out = {
        name: {"entries": slot["entries"], "replicons": len(slot["replicons"]),
               "species": len(slot["species"])}
        for name, slot in sorted(countries.items(),
                                 key=lambda kv: -kv[1]["entries"])
    }
    year_out = {
        str(year): {"entries": slot["entries"],
                    "replicons": len(slot["replicons"])}
        for year, slot in sorted(years.items())
    }
    return country_out, year_out


# ------------------------------------------------------------------- rapor
def print_summary(result):
    cov = result["coverage"]
    print("\n=== Kapsama (her sayinin paydasi) ===")
    print(f"  dogrulanmis RO                  : {cov['confirmed_ros_total']}")
    print(f"  replikon (>=1 dogrulanmis RO)   : {cov['replicons_with_confirmed_ro']}")
    print(f"  isolation_source olan replikon  : {cov['replicons_with_isolation_source']}"
          f"  ({100 * cov['replicons_with_isolation_source'] / max(1, cov['replicons_in_table']):.1f}%)")
    print(f"  host olan                       : {cov['replicons_with_host']}")
    print(f"  cografya olan                   : {cov['replicons_with_geo']}")
    print(f"  toplama yili olan               : {cov['replicons_with_collection_year']}")
    print(f"  SINIFLANMIS habitat (replikon)  : {cov['replicons_with_classified_habitat']}")
    print(f"  RO, isolation_source'li replikonda: {cov['confirmed_ros_with_isolation_source']}"
          f"  ({100 * cov['confirmed_ros_with_isolation_source_fraction']:.1f}%)")
    print(f"  RO, siniflanmis habitatta       : {cov['confirmed_ros_with_classified_habitat']}"
          f"  ({100 * cov['confirmed_ros_with_classified_habitat_fraction']:.1f}%)")
    print(f"  metin var ama bilgi yok         : {cov['isolation_source_uninformative_text']}")
    print(f"  metin var, kural tutmadi (other): {cov['unclassified_strings_entries']} replikon, "
          f"{cov['unclassified_strings_distinct']} farkli metin")
    print(f"  other/unknown ama host var      : {cov['replicons_unclassified_but_host_present']}")

    print("\n=== Habitat dagilimi (giris = dogrulanmis RO) ===")
    print(f"  {'habitat':24s} {'giris':>7s} {'%':>6s} {'replikon':>9s} {'tur':>6s}")
    for habitat, slot in result["habitats"].items():
        print(f"  {habitat:24s} {slot['entries']:7d} "
              f"{100 * slot['entry_fraction']:6.1f} {slot['replicons']:9d} "
              f"{slot['species']:6d}")

    print("\n=== Siniflanamayan en sik 15 metin (hepsi 'other') ===")
    for item in result["coverage"]["unclassified_top_50"][:15]:
        print(f"  {item['replicons']:5d}  {item['isolation_source'][:70]}")

    print("\n=== Substrat sinifi x habitat: giris payi vs TUR payi ===")
    print("    (bias_delta > 0 => giris sayilari bu habitati tur duzeyinden fazla gosteriyor)")
    for cls, slot in sorted(result["substrate_class_habitat_profile"].items()):
        print(f"\n  {cls}  (giris {slot['entries_total']}, "
              f"farkli tur {slot['species_distinct']})")
        ranked = sorted(slot["habitats"].items(),
                        key=lambda kv: -kv[1]["species_fraction"])
        for habitat, cell in ranked[:8]:
            print(f"    {habitat:24s} giris {100 * cell['entry_fraction']:5.1f}%  "
                  f"tur {100 * cell['species_fraction']:5.1f}%  "
                  f"(n_tur {cell['species']:4d})  delta {cell['bias_delta']:+.3f}")

    print("\n=== Ulke (ilk 10) ===")
    for name, slot in list(result["countries"].items())[:10]:
        print(f"  {name[:30]:30s} giris {slot['entries']:5d}  tur {slot['species']:4d}")

    years = result["collection_years"]
    if years:
        keys = list(years)
        print(f"\n=== Toplama yili: {keys[0]}-{keys[-1]}, "
              f"{sum(v['entries'] for v in years.values())} giris ===")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--gbk-dir", default="gbk_files")
    ap.add_argument("--out", default="analysis_out/habitat.json")
    ap.add_argument("--ecology", default="cluster_ecology.csv")
    ap.add_argument("--chemistry", default="chemistry.csv")
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--force", action="store_true",
                    help="replicon_source dolu olsa da gbk'lari yeniden parse et")
    args = ap.parse_args()

    connection = sqlite3.connect(args.db)
    connection.executescript(SCHEMA)

    fallbacks = parse_sources(connection, args.gbk_dir, args.cpu, args.force)
    keyword_hits = assign_habitats(connection)

    ecology = read_csv_map(args.ecology)
    chemistry = read_csv_map(args.chemistry)
    entries = load_entries(connection, ecology, chemistry)

    countries, years = geography_block(entries)
    result = {
        "method": {
            "habitat_vocabulary": list(HABITAT_ORDER),
            "assignment": ("ordered keyword map on the GenBank /isolation_source "
                           "qualifier only; first matching rule wins; /host is "
                           "stored but never used to assign a habitat"),
            "text_normalisation": ("lowercase, every non-alphanumeric character "
                                   "becomes a space, whole-word match with an "
                                   "optional trailing s/es"),
            "species_definition": ("first two words of the GenBank organism name; "
                                   "a pseudo-species such as 'Pseudomonas sp.' "
                                   "collapses a whole genus into one unit"),
            "collection_year": ("first four-digit year found in /collection_date; "
                                "month and day are dropped and a range such as "
                                "'1900/1969' yields its first year"),
            "rule_order": [{"habitat": h, "keywords": list(k)}
                           for h, k in HABITAT_RULES],
            "uninformative_text": sorted(UNINFORMATIVE_EXACT),
        },
        "coverage": coverage_block(connection, entries, fallbacks),
        "habitats": habitat_summary(entries),
        "habitat_keyword_hits": keyword_hits,
        "habitat_by_substrate_class": crosstab(entries, "substrate_class"),
        "habitat_by_chemical_family": crosstab(entries, "family"),
        "substrate_class_habitat_profile":
            species_normalised_profile(entries, "substrate_class"),
        "chemical_family_habitat_profile":
            species_normalised_profile(entries, "family"),
        "per_ro_type": per_type_profile(entries, ecology, chemistry),
        "countries": countries,
        "collection_years": years,
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(result, handle, indent=2, sort_keys=False)
    print_summary(result)
    print(f"\n[yazildi] {args.out}")
    connection.close()


if __name__ == "__main__":
    main()
