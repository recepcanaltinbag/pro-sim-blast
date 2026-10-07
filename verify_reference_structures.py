"""
Referans yapilar -- 2025 derlemesinin Tablo 1'ini RCSB'ye ve 71 tipe karsi dogrular.

NE DOGRULANIYOR
  Miao H, Schmidt S. Biochemistry 2025;64:3801-3813 (doi:10.1021/acs.biochem.5c00369,
  CC BY 4.0) derlemesinin Tablo 1'i, deneysel olarak COZULMUS Rieske oksijenaz
  yapilarinin atif verilebilir bir envanteridir: 26 satir, her biri bir kisa ad,
  organizma, ALT BIRIM BILESIMI (a3 / a3b3 / a3a'3 / a3a3) ve bir PDB kimligi.
  Bu script o 26 satiri reference_structures.csv'den okur ve

    1. her PDB kimligini RCSB REST API'sinde cozer: baslik, kaynak organizma,
       ayri polimer varliklari ve zincir sayilari, birincil atif DOI/PubMed,
       superseded/obsoleted durumu;
    2. her satiri chemistry.csv'deki bir tipe (cluster) baglar -- UC kanit
       sinifindan EN GUCLUSU ile: (a) pdb_match, kimlik zaten chemistry.csv'nin
       pdb sutununda; (b) organism_and_name, tip adi derlemedeki kisa ad ile ve
       organizma source/notes metni ile ortusuyor; (c) sequence_identity,
       RCSB'nin varlik dizisi kuratorlu referans dizisi ile >= %95 ayni.
       Gerekcelendirilemeyen satir BOS BIRAKILIR ve sebebi mapping_note'a yazilir
       -- tahmin URETILMEZ;
    3. derlemenin DENEYSEL bilesim etiketini, bu veritabaninin GENOMIK BAGLAMDAN
       CIKARDIGI beta alt birim varligi ile karsilastirir.

NEDEN
  Iki ayri acik var ve ikisi de bu tabloyla kapaniyor.

  (1) Kuratorlu set 71 tip tanimliyor; bunlarin 50'sinde atif verilebilir bir
      kaynak, 7'sinde hicbir yapi (ne deneysel ne AlphaFold) yok. Derleme hakemli
      ve acik erisimli bir envanter oldugu icin, eslesen her satir o tipe
      atif verilebilir bir DENEYSEL yapi kazandirir.

  (2) Site, beta alt biriminin varligini (ro_etc.has_beta) genomik baglamdan
      CIKARIYOR ve istatistik sayfasi bunu beta_by_type / Cramer V = 0.94 olarak
      DISARIDAN DOGRULANMAMIS bicimde sunuyor. Derlemenin bilesim sutunu ise
      kristalografik olarak OLCULMUS bir a3-karsi-a3b3 etiketidir. Yani eslesen
      her tip icin bagimsiz bir dogruluk olcutu veriyor. architecture_validation
      bolumu tam olarak bu karsilastirmadir: cikarimin disaridan ilk sinavi.

  Ayrica iki yonlu uyusmazlik aranir: derlemenin kimligi ile chemistry.csv'nin
  kimligi ayni proteinin FARKLI DURUMLARI mi (apo/holo, baska ligand, baska
  cozunurluk -- o zaman ikisi de dogru, fark ONEMSIZ), yoksa GERCEKTEN farkli
  proteinler mi (o zaman biri kurasyon hatasi)? Karar dizi kimligi ile verilir,
  baslik benzerligi ile degil.

OLCULEN SEYLER
  - PDB kimliginin cozulup cozulmedigi, deneysel yontem, cozunurluk
  - RCSB kaynak organizmasinin derlemedeki organizma ile cins/tur duzeyinde
    ortusmesi (yeniden siniflandirma ayri raporlanir, uyusmazlik sayilmaz)
  - AYRI polimer varligi sayisi: a3 -> 1 katalitik polimer, a3b3 -> 2,
    a3a'3 -> 2 farkli katalitik polimer, a3a3 -> 1. Kristalizasyon icin eklenen
    yabanci polimerler (ornegin lizozim) kaynak organizmasi farkli oldugu icin
    sayilmaz.
  - birincil atif DOI ve PubMed kimligi
  - RCSB varlik dizisi ile kuratorlu referans dizisi arasindaki kimlik yuzdesi
    (blastp; nident / referans uzunlugu -- yani referansin TAMAMI uzerinden)
  - eslesen her tip icin: derlemenin bilesimi vs ro_etc.has_beta cogunlugu

Girdi : reference_structures.csv, chemistry.csv, ROs_71_Clean/refs71.fasta, roar.sqlite
Cikti : reference_structures.csv (yerinde genisletilir), analysis_out/reference_structures.json
Onbellek: structures/rcsb_entry/<PDB>.json, structures/rcsb_entity/<PDB>_<N>.json
          --offline ile tekrar kosusta hic istek gitmez.
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sqlite3
import sys
import tempfile
import time
import urllib.error
import urllib.request
from collections import defaultdict

ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{0}"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{0}/{1}"
USER_AGENT = "ROAR-DB/1.0 (Rieske oxygenase database; reference structure verification)"
HTTP_TIMEOUT_S = 30
REQUEST_PAUSE_S = 0.5

# Derlemenin bilesim etiketi -> beklenen AYRI katalitik polimer sayisi ve
# beta alt biriminin varligi. a3a'3 iki FARKLI katalitik alt birimdir (alfa ve
# alfa-prime), beta degil; a3a3 ayni alfanin homoheksameridir.
COMPOSITION = {
    "a3":    {"distinct": 1, "beta": False, "label": "alpha3"},
    "a3b3":  {"distinct": 2, "beta": True,  "label": "alpha3beta3"},
    "a3a'3": {"distinct": 2, "beta": False, "label": "alpha3alpha'3"},
    "a3a3":  {"distinct": 1, "beta": False, "label": "alpha3alpha3"},
}

# CSV'ye eklenecek sutunlar -- mevcut sutunlar ve satir sirasi KORUNUR.
NEW_COLUMNS = ["mapped_cluster", "mapping_basis", "mapping_identity_pct",
               "rcsb_title", "rcsb_organism", "rcsb_entities",
               "rcsb_citation_doi", "agrees_with_review", "mapping_note"]

MIN_SEQ_IDENTITY = 95.0
_STOPWORDS = {"sp.", "sp", "strain", "subsp.", "str."}


# --------------------------------------------------------------------------- #
# Onbellekli HTTP
# --------------------------------------------------------------------------- #

def fetch_json(url, cache_path, offline=False, pause=REQUEST_PAUSE_S):
    """JSON uc noktasini onbellekli cek. Donen: (veri|None, durum).

    Onbellekte varsa AG'A CIKMAZ. HTTP hatasi da onbeleklenir, cunku
    "bu kimlik cozulmuyor" da bir sonuctur ve tekrar sorulmasi gerekmez.
    """
    if cache_path and os.path.exists(cache_path) and os.path.getsize(cache_path) > 0:
        try:
            with open(cache_path) as handle:
                return json.load(handle), "cached"
        except (ValueError, OSError):
            pass
    if offline:
        return None, "absent from cache and --offline was given"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_S) as response:
            payload = response.read()
    except urllib.error.HTTPError as exc:
        data = {"_http_error": exc.code}
        _write_cache(cache_path, data)
        time.sleep(pause)
        return data, "http %d" % exc.code
    except (urllib.error.URLError, OSError) as exc:
        return None, "request failed: %s" % exc
    try:
        data = json.loads(payload)
    except ValueError as exc:
        return None, "response is not valid JSON: %s" % exc
    _write_cache(cache_path, data)
    time.sleep(pause)
    return data, "downloaded"


def _write_cache(cache_path, data):
    if not cache_path:
        return
    directory = os.path.dirname(cache_path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    temporary = cache_path + ".part"
    with open(temporary, "w") as handle:
        json.dump(data, handle)
    os.replace(temporary, cache_path)


# --------------------------------------------------------------------------- #
# FASTA / dizi yardimcilari
# --------------------------------------------------------------------------- #

def read_fasta(path):
    """{baslik: dizi} -- baslik ilk bosluga kadar."""
    out, name, buf = {}, None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        out[name] = "".join(buf)
    return out


def cluster_of(header):
    """refs71.fasta basligindan tip kimligi: ilk uc alt cizgi alani."""
    return "_".join(header.split("_")[:3])


def name_of(cluster):
    """3_315_NDO -> NDO"""
    parts = cluster.split("_")
    return parts[2] if len(parts) > 2 else cluster


def organism_tokens(text):
    """Organizma adindan anlamli kelimeler: cins, tur, sus -- dolgu kelimeler atilir."""
    tokens = [t.strip(",.").lower() for t in (text or "").split()]
    return [t for t in tokens if t and t not in _STOPWORDS]


def species_epithet(text):
    """'Comamonas testosteroni' -> 'testosteroni'; 'Comamonas sp. E6' -> ''.

    Bir TUR IDDIASI yoksa bos doner. 'sp.' bir tur adi degildir ve bir sus kodu
    (E6, U2, RHA1, JS765) da degildir: rakam iceren ya da cok kisa belirtecler
    tur epiteti SAYILMAZ. Boylece 'Ralstonia sp. strain U2' ile 'Ralstonia sp.'
    arasinda uydurma bir tur uyusmazligi raporlanmaz.
    """
    raw = [t.strip(",.").lower() for t in (text or "").split()]
    if len(raw) < 2:
        return ""
    candidate = raw[1]
    if candidate in _STOPWORDS or not candidate.isalpha() or len(candidate) < 4:
        return ""
    return candidate


# --------------------------------------------------------------------------- #
# RCSB kayitlarini okuma
# --------------------------------------------------------------------------- #

def entity_sequence(entity):
    raw = (entity.get("entity_poly", {}) or {}).get("pdbx_seq_one_letter_code_can") or ""
    return re.sub(r"\s+", "", raw).upper()


def entity_organisms(entity):
    out = []
    for source in entity.get("rcsb_entity_source_organism", []) or []:
        name = source.get("ncbi_scientific_name") or source.get("scientific_name")
        if name and name not in out:
            out.append(name)
    return out


def read_entry(pdb_id, cache_dir, offline, pause):
    """Bir PDB kaydini ve TUM polimer varliklarini oku. Donen: dict."""
    entry, status = fetch_json(ENTRY_URL.format(pdb_id),
                               os.path.join(cache_dir, "rcsb_entry", pdb_id + ".json"),
                               offline, pause)
    record = {"pdb": pdb_id, "fetch_status": status, "resolves": False,
              "title": "", "organisms": [], "entities": [], "citation": {},
              "superseded": [], "method": "", "resolution": None}
    if entry is None:
        record["error"] = status
        return record
    if "_http_error" in entry:
        record["error"] = "RCSB returned HTTP %s -- the id does not resolve" % entry["_http_error"]
        return record
    record["resolves"] = True
    record["title"] = (entry.get("struct", {}) or {}).get("title", "") or ""
    info = entry.get("rcsb_entry_info", {}) or {}
    record["method"] = info.get("experimental_method", "") or ""
    resolution = info.get("resolution_combined") or []
    record["resolution"] = resolution[0] if resolution else None
    record["ligands"] = info.get("nonpolymer_bound_components") or []

    # Birincil atif: id == "primary". Yoksa ilk atif.
    citations = entry.get("citation", []) or []
    primary = next((c for c in citations if c.get("id") == "primary"), None) \
        or (citations[0] if citations else {})
    record["citation"] = {
        "doi": primary.get("pdbx_database_id_DOI") or "",
        "pubmed": str(primary.get("pdbx_database_id_PubMed") or ""),
        "title": primary.get("title") or "",
        "journal": primary.get("rcsb_journal_abbrev") or primary.get("journal_abbrev") or "",
        "year": primary.get("year"),
    }

    # Obsoleted / superseded bayraklari -- eski bir kaydin yerine gecmis olmak
    # sorun degil, ama KENDISI eskimis bir kaydi atif vermek sorundur.
    for obs in entry.get("pdbx_database_PDB_obs_spr", []) or []:
        record["superseded"].append({
            "id": obs.get("id"), "pdb_id": obs.get("pdb_id"),
            "replace_pdb_id": obs.get("replace_pdb_id"), "date": obs.get("date")})

    identifiers = entry.get("rcsb_entry_container_identifiers", {}) or {}
    for entity_id in identifiers.get("polymer_entity_ids", []) or []:
        entity, estatus = fetch_json(
            ENTITY_URL.format(pdb_id, entity_id),
            os.path.join(cache_dir, "rcsb_entity", "%s_%s.json" % (pdb_id, entity_id)),
            offline, pause)
        if entity is None or "_http_error" in entity:
            record["entities"].append({"entity_id": entity_id, "error": estatus})
            continue
        container = entity.get("rcsb_polymer_entity_container_identifiers", {}) or {}
        chains = container.get("auth_asym_ids", []) or []
        references = [
            "%s:%s" % (r.get("database_name"), r.get("database_accession"))
            for r in (container.get("reference_sequence_identifiers") or [])]
        sequence = entity_sequence(entity)
        record["entities"].append({
            "entity_id": entity_id,
            "description": (entity.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description", "") or "",
            "chains": chains,
            "chain_count": len(chains),
            "length": len(sequence),
            "organisms": entity_organisms(entity),
            "reference_ids": references,
            "sequence": sequence,
        })
    for ent in record["entities"]:
        for org in ent.get("organisms", []):
            if org not in record["organisms"]:
                record["organisms"].append(org)
    return record


# --------------------------------------------------------------------------- #
# Dizi kimligi
# --------------------------------------------------------------------------- #

def identity_table(entity_fasta, refs_path, threads=4):
    """{varlik_adi: [(cluster, kimlik_yuzdesi, nident, ref_len), ...]} -- blastp.

    Kimlik, ALIGNMENT uzerinden degil REFERANSIN TAMAMI uzerinden hesaplanir
    (nident / referans uzunlugu). Boylece kismi bir hizalanma yuksek bir yerel
    pident ile yanlisca "ayni protein" gibi gorunmez.
    """
    if not shutil.which("blastp"):
        return None
    with tempfile.TemporaryDirectory() as tmp:
        query = os.path.join(tmp, "entities.fasta")
        with open(query, "w") as handle:
            handle.write(entity_fasta)
        command = ["blastp", "-query", query, "-subject", refs_path,
                   "-outfmt", "6 qseqid sseqid pident nident length qlen slen bitscore",
                   "-max_target_seqs", "500", "-evalue", "1e-3",
                   "-num_threads", str(threads), "-seg", "no"]
        try:
            raw = subprocess.run(command, check=True, capture_output=True,
                                 text=True).stdout
        except (subprocess.CalledProcessError, OSError) as exc:
            print("[uyari] blastp basarisiz: %s" % exc, file=sys.stderr)
            return None
    table = defaultdict(list)
    for line in raw.splitlines():
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        query_id, subject = fields[0], fields[1]
        nident, ref_len = int(fields[3]), int(fields[6])
        if not ref_len:
            continue
        table[query_id].append((cluster_of(subject), 100.0 * nident / ref_len,
                                nident, ref_len, float(fields[7])))
    # Ayni tip icin birden fazla HSP olabilir: en iyisini tut.
    out = {}
    for query_id, hits in table.items():
        best = {}
        for cluster, pct, nident, ref_len, bits in hits:
            if cluster not in best or pct > best[cluster][1]:
                best[cluster] = (cluster, pct, nident, ref_len, bits)
        out[query_id] = sorted(best.values(), key=lambda h: -h[1])
    return out


# --------------------------------------------------------------------------- #
# Eslestirme
# --------------------------------------------------------------------------- #

def pdb_index(chemistry_rows):
    """{PDB kimligi: [cluster, ...]} -- chemistry.csv'nin pdb sutunundan."""
    index = defaultdict(list)
    for row in chemistry_rows:
        for token in re.findall(r"\b\d[A-Za-z0-9]{3}\b", (row.get("pdb") or "")):
            token = token.upper()
            if row["cluster"] not in index[token]:
                index[token].append(row["cluster"])
    return index


def map_row(row, record, chem_by_cluster, pdb_by_id, identities, refs_by_cluster):
    """Bir derleme satirini bir tipe bagla.

    Kanit siniflari guc sirasiyla denenir: pdb_match > organism_and_name >
    sequence_identity. Hicbiri gerekcelenemiyorsa mapped_cluster BOS kalir.
    Donen: dict(mapped_cluster, mapping_basis, mapping_identity_pct, notes[],
                alternatives[], best_identity)
    """
    notes, result = [], {"mapped_cluster": "", "mapping_basis": "unmapped",
                         "mapping_identity_pct": "", "alternatives": [],
                         "best_identity": None, "best_identity_cluster": ""}
    ro_name = (row.get("ro_name") or "").strip()
    pdb = (row.get("pdb") or "").strip().upper()

    # Bu kaydin hangi varligi "alfa"? refs71'e en iyi vuran varlik.
    best = None
    for entity in record.get("entities", []):
        key = "%s_%s" % (record["pdb"], entity["entity_id"])
        for hit in (identities or {}).get(key, []):
            if best is None or hit[1] > best[1]:
                best = (hit[0], hit[1], hit[2], hit[3], entity["entity_id"])
    if best:
        result["best_identity"] = round(best[1], 2)
        result["best_identity_cluster"] = best[0]

    seq_ranked = []
    if best:
        key = "%s_%s" % (record["pdb"], best[4])
        seq_ranked = [(h[0], h[1]) for h in (identities or {}).get(key, [])]

    # --- (a) pdb_match -----------------------------------------------------
    candidates = list(pdb_by_id.get(pdb, []))
    if candidates:
        # Birden fazla tip ayni kimligi tasiyabilir (kuratorlu sette ayni enzim
        # iki adla durabilir, ya da kimlik YANLIS bir tipe de yazilmis olabilir).
        # Once dizi kimligi, sonra ad ortusmesi ile ayirt et.
        scored = []
        by_pct = dict(seq_ranked)
        for cluster in candidates:
            pct = by_pct.get(cluster)
            scored.append((pct if pct is not None else -1.0,
                           1 if name_of(cluster).lower() == ro_name.lower() else 0,
                           cluster))
        scored.sort(reverse=True)
        top = scored[0]
        # Dizi verisi varsa ve EN IYI aday bile referansla ortusmuyorsa, kimlik
        # yanlis tipe yazilmis demektir: pdb_match'i KULLANMA.
        if top[0] >= MIN_SEQ_IDENTITY or (top[0] < 0 and identities is None):
            result["mapped_cluster"] = top[2]
            result["mapping_basis"] = "pdb_match"
            if top[0] >= 0:
                # mapping_identity_pct sutunu yalnizca sequence_identity icin
                # doldurulur; sayi yine de gorunur olsun diye nota yazilir.
                notes.append("sequence corroborates the pdb_match: %.1f%% identical "
                             "to this type's curated reference" % top[0])
            others = [c for _, _, c in scored[1:]]
            for pct, _, cluster in scored[1:]:
                if pct >= MIN_SEQ_IDENTITY:
                    result["alternatives"].append(
                        "%s (%.1f%% identical reference, duplicate curated name "
                        "carrying the same PDB id)" % (cluster, pct))
                elif pct >= 0:
                    result["alternatives"].append(
                        "%s (carries %s in chemistry.csv but its curated reference "
                        "is only %.1f%% identical to this entry -- CURATION ERROR)"
                        % (cluster, pdb, pct))
            if others:
                notes.append("chemistry.csv assigns %s to %d types: %s"
                             % (pdb, len(candidates), ", ".join(candidates)))
            return result, notes
        notes.append("%s appears in chemistry.csv for %s, but no curated reference "
                     "there reaches %.0f%% identity with this entry (best %.1f%%); "
                     "pdb_match was rejected"
                     % (pdb, ", ".join(candidates), MIN_SEQ_IDENTITY, max(top[0], 0.0)))

    # --- (b) organism_and_name ---------------------------------------------
    review_tokens = set(organism_tokens(row.get("organism")))
    name_hits = []
    for cluster, chem in chem_by_cluster.items():
        if name_of(cluster).lower() != ro_name.lower():
            continue
        haystack = ((chem.get("source") or "") + " " + (chem.get("notes") or "")).lower()
        shared = sorted(t for t in review_tokens if len(t) > 3 and t in haystack)
        if shared:
            name_hits.append((len(shared), cluster, shared))
    if name_hits:
        name_hits.sort(reverse=True)
        _, cluster, shared = name_hits[0]
        result["mapped_cluster"] = cluster
        result["mapping_basis"] = "organism_and_name"
        notes.append("type name '%s' matches the review's RO name and the organism "
                     "term(s) %s appear in this type's source/notes"
                     % (name_of(cluster), "/".join(shared)))
        pct = dict(seq_ranked).get(cluster)
        if pct is not None:
            notes.append("sequence corroborates: %.1f%% identical to the curated "
                         "reference" % pct)
        for other, other_pct in seq_ranked:
            if other != cluster and other_pct >= MIN_SEQ_IDENTITY:
                result["alternatives"].append(
                    "%s (%.1f%% identical reference -- duplicate in the curated set)"
                    % (other, other_pct))
        return result, notes

    # --- (c) sequence_identity ---------------------------------------------
    if seq_ranked and seq_ranked[0][1] >= MIN_SEQ_IDENTITY:
        # Ad ortusmesi varsa esit kimlikli adaylar arasinda onu sec.
        top_pct = seq_ranked[0][1]
        tied = [c for c, p in seq_ranked if p >= top_pct - 0.001]
        chosen = next((c for c in tied if name_of(c).lower() == ro_name.lower()), tied[0])
        result["mapped_cluster"] = chosen
        result["mapping_basis"] = "sequence_identity"
        result["mapping_identity_pct"] = "%.1f" % dict(seq_ranked)[chosen]
        for other, other_pct in seq_ranked:
            if other != chosen and other_pct >= MIN_SEQ_IDENTITY:
                result["alternatives"].append(
                    "%s (%.1f%% identical reference -- duplicate in the curated set)"
                    % (other, other_pct))
        return result, notes

    # --- eslesmedi ---------------------------------------------------------
    if identities is None:
        notes.append("blastp is not available, so the sequence basis could not be "
                     "tested; no PDB id or name/organism match either")
    elif seq_ranked:
        notes.append("no curated reference reaches %.0f%% identity (best: %s at "
                     "%.1f%%), and neither the PDB id nor the name/organism match a "
                     "type -- this enzyme has no counterpart among the 71 types"
                     % (MIN_SEQ_IDENTITY, seq_ranked[0][0], seq_ranked[0][1]))
    else:
        notes.append("no curated reference produced a blastp hit at all, and neither "
                     "the PDB id nor the name/organism match a type")
    return result, notes


# --------------------------------------------------------------------------- #
# Derleme ile RCSB arasindaki uyum
# --------------------------------------------------------------------------- #

def check_agreement(row, record):
    """Derlemenin satiri ile RCSB kaydi uyusuyor mu? Donen: (verdict, notes[])."""
    notes, verdict = [], "yes"
    if not record.get("resolves"):
        return "no", ["the PDB id does not resolve at RCSB: %s"
                      % record.get("error", "unknown reason")]

    # --- organizma ---------------------------------------------------------
    # Yabanci polimerlerin (kristalizasyon katkilari) organizmasi bu kontrole
    # KARISMAMALI: yalnizca derlemenin cinsini tasiyan adlar kiyaslanir.
    review = organism_tokens(row.get("organism"))
    review_genus = review[0] if review else ""
    all_names = [n for n in record.get("organisms", []) if n]
    own_names = [n for n in all_names if organism_tokens(n)[:1] == [review_genus]]
    if not own_names:
        verdict = "no"
        notes.append("organism disagrees: the review says '%s', RCSB says %s"
                     % (row.get("organism"), " / ".join(all_names) or "nothing"))
    else:
        review_species = species_epithet(row.get("organism"))
        rcsb_species = [species_epithet(n) for n in own_names]
        rcsb_species = [e for e in rcsb_species if e]
        if review_species and rcsb_species and review_species not in rcsb_species:
            verdict = "partial"
            notes.append("organism genus agrees but the species epithet disagrees: "
                         "the review says '%s', RCSB says %s"
                         % (row.get("organism"), " / ".join(own_names)))
        elif bool(review_species) != bool(rcsb_species):
            # Bir taraf tur adi veriyor, digeri 'sp.' diyor. Bu bir CELISKI
            # degil, farkli bir ozgulluk duzeyi: verdict dusurulmez.
            notes.append("the two sources name this organism at different levels of "
                         "specificity -- the review says '%s', RCSB says %s; the "
                         "genus agrees and this is nomenclature, not a different "
                         "protein" % (row.get("organism"), " / ".join(own_names)))

        # Birincil atif suşu BASKA bir cinsle adlandiriyorsa (yeniden
        # siniflandirma) bunu soyle: atif verilen makale ile depozit ayrisiyor.
        title = (record.get("citation", {}) or {}).get("title", "") or ""
        strain = [t for t in review[1:]
                  if t and (any(c.isdigit() for c in t) or len(t) <= 3)]
        lowered = title.lower()
        if title and review_genus and review_genus not in lowered \
                and any(t in lowered for t in strain):
            notes.append("the primary citation names this strain under a different "
                         "genus than both the deposition and the review (citation "
                         "title: \"%s\") -- a strain reclassification, not a "
                         "different protein" % title)

    # --- bilesim -----------------------------------------------------------
    spec = COMPOSITION.get((row.get("subunit_composition") or "").strip())
    counted, foreign = [], []
    for entity in record.get("entities", []):
        organisms = entity.get("organisms", [])
        if organisms and review_genus and \
                not any(organism_tokens(o)[:1] == [review_genus] for o in organisms):
            foreign.append(entity)
        else:
            counted.append(entity)
    if foreign:
        notes.append("%d polymer entit%s from a different organism (%s) %s not "
                     "counted toward the subunit composition -- a crystallization "
                     "additive, not part of the RO"
                     % (len(foreign), "y" if len(foreign) == 1 else "ies",
                        "; ".join("%s / %s" % (e["description"],
                                               "/".join(e["organisms"]))
                                  for e in foreign),
                        "is" if len(foreign) == 1 else "are"))
    if spec is None:
        verdict = "partial" if verdict == "yes" else verdict
        notes.append("the review's composition label '%s' is not one this script "
                     "knows how to check" % row.get("subunit_composition"))
    elif len(counted) != spec["distinct"]:
        if spec["label"] == "alpha3alpha'3" and len(counted) == 1:
            verdict = "partial"
            notes.append("the review's alpha3alpha'3 label describes the "
                         "NdmA3/NdmB3 heterocomplex, but this entry contains only "
                         "one of the two catalytic subunits (%d distinct polymer), "
                         "so the composition cannot be confirmed from this entry "
                         "alone" % len(counted))
        else:
            verdict = "no"
            notes.append("composition disagrees: the review says %s (expects %d "
                         "distinct catalytic polymer%s) but the entry has %d: %s"
                         % (spec["label"], spec["distinct"],
                            "" if spec["distinct"] == 1 else "s", len(counted),
                            "; ".join("%s (%d chains, %d aa)"
                                      % (e["description"], e["chain_count"],
                                         e["length"]) for e in counted)))
    else:
        chains = sum(e["chain_count"] for e in counted)
        if chains < 3:
            notes.append("the number of distinct polymers (%d) matches %s, but the "
                         "deposited asymmetric unit holds only %d chain%s -- the "
                         "trimer is generated by crystallographic symmetry and is "
                         "not verifiable from the entry record alone"
                         % (len(counted), spec["label"], chains,
                            "" if chains == 1 else "s"))

    # --- kaydin kendi durumu ----------------------------------------------
    for obs in record.get("superseded", []):
        if obs.get("id") == "SPRSDE":
            notes.append("this entry supersedes the earlier %s (%s) -- the review "
                         "cites the current entry, which is correct"
                         % (obs.get("replace_pdb_id"), (obs.get("date") or "")[:10]))
        else:
            verdict = "no"
            notes.append("RCSB flags this entry as %s (%s)"
                         % (obs.get("id"), obs.get("replace_pdb_id")))

    if not record.get("citation", {}).get("doi") and \
            not record.get("citation", {}).get("pubmed"):
        notes.append("RCSB has no DOI or PubMed id for the primary citation "
                     "(journal: %s) -- the structure is citable only by its PDB id"
                     % (record["citation"].get("journal") or "unknown"))
    return verdict, notes


# --------------------------------------------------------------------------- #
# Kurasyon farklari: ayni protein mi, farkli protein mi?
# --------------------------------------------------------------------------- #

def compare_curated_ids(rows, records, chem_by_cluster, identities):
    """Derlemenin kimligi ile chemistry.csv'nin kimligi ayni mi?

    Ayni degilse karar DIZIYE gore verilir: iki kayit ayni diziyi tasiyorsa
    fark ONEMSIZDIR (apo/holo, baska ligand, baska cozunurluk); dizi farkliysa
    biri KURASYON HATASIDIR.
    """
    out = []
    for row, record, mapping in rows:
        cluster = mapping["mapped_cluster"]
        if not cluster:
            continue
        chem = chem_by_cluster.get(cluster, {})
        curated = (chem.get("pdb") or "").strip().upper()
        review = (row.get("pdb") or "").strip().upper()
        if not curated or curated == review:
            continue
        other = records.get(curated)
        item = {"cluster": cluster, "review_pdb": review, "curated_pdb": curated,
                "material": None, "verdict": "", "detail": ""}
        if other is None or not other.get("resolves"):
            item["material"] = True
            item["verdict"] = "the curated id could not be resolved at RCSB"
            out.append(item)
            continue

        def alpha(rec):
            best = None
            for entity in rec.get("entities", []):
                key = "%s_%s" % (rec["pdb"], entity["entity_id"])
                for hit in (identities or {}).get(key, []):
                    if hit[0] == cluster and (best is None or hit[1] > best[0]):
                        best = (hit[1], entity)
            if best:
                return best[1], best[0]
            longest = max(rec.get("entities", []) or [{}],
                          key=lambda e: e.get("length", 0))
            return longest, None

        review_entity, review_pct = alpha(record)
        curated_entity, curated_pct = alpha(other)
        seq_a = review_entity.get("sequence", "")
        seq_b = curated_entity.get("sequence", "")
        same_organism = bool(set(record.get("organisms", [])) &
                             set(other.get("organisms", [])))
        if seq_a and seq_b and seq_a == seq_b:
            item["material"] = False
            item["verdict"] = "the same protein in different states"
            item["detail"] = (
                "identical deposited sequences (%d aa, 0 differences), same source "
                "organism; they differ only in state: %s has %s at %s A, %s has %s "
                "at %s A. Both ids are correct and the difference is immaterial."
                % (len(seq_a), review, _state(record), record.get("resolution"),
                   curated, _state(other), other.get("resolution")))
        elif seq_a and seq_b and len(seq_a) == len(seq_b):
            diffs = [(i + 1, p, q) for i, (p, q) in enumerate(zip(seq_a, seq_b))
                     if p != q]
            if len(diffs) <= 3 and same_organism:
                item["material"] = False
                item["verdict"] = ("the same protein; the curated id is a point "
                                   "variant of it")
                item["detail"] = (
                    "same length (%d aa) and same organism, differing at %d "
                    "position%s: %s. %s is the wild type (%s, %s A), %s is the "
                    "substituted form (%s, %s A). Same protein, so the difference "
                    "is immaterial as to identity -- but the curated id is an "
                    "engineered variant, which matters if the entry is used as the "
                    "wild-type active-site reference."
                    % (len(seq_a), len(diffs), "" if len(diffs) == 1 else "s",
                       ", ".join("%s%d%s" % (p, i, q) for i, p, q in diffs),
                       review, _state(record), record.get("resolution"),
                       curated, _state(other), other.get("resolution")))
            else:
                item["material"] = True
                item["verdict"] = "genuinely different proteins"
                item["detail"] = ("same length but %d sequence differences"
                                  % len(diffs))
        else:
            item["material"] = True
            item["verdict"] = "genuinely different proteins"
            item["detail"] = (
                "different proteins: %s is %s (%s, %d aa%s) and %s is %s (%s, %d "
                "aa%s). The curated reference of %s is %s identical to %s, so the "
                "curated id is the one that represents this type; the review's "
                "table simply does not list it."
                % (review, review_entity.get("description", "?"),
                   "/".join(record.get("organisms", [])) or "?",
                   review_entity.get("length", 0),
                   "" if review_pct is None else ", %.1f%% to the curated reference" % review_pct,
                   curated, curated_entity.get("description", "?"),
                   "/".join(other.get("organisms", [])) or "?",
                   curated_entity.get("length", 0),
                   "" if curated_pct is None else ", %.1f%% to the curated reference" % curated_pct,
                   cluster,
                   "more" if (curated_pct or 0) >= (review_pct or 0) else "less",
                   review))
        out.append(item)
    return out


def curated_ids_vs_review(chemistry_rows, review_rows, triples, records, identities):
    """chemistry.csv'nin TASIDIGI her PDB kimligini derlemenin listesine karsi koy.

    Uc sonuc mumkun:
      listed_and_matches       -- kimlik Tablo 1'de var ve tipin referans dizisi
                                  o kayda uyuyor;
      listed_but_mismatched    -- kimlik Tablo 1'de var ama tipin KENDI referansi
                                  o kayitla ortusmuyor: kimlik yanlis tipe
                                  yazilmis, KURASYON HATASI;
      absent_from_review       -- kimlik Tablo 1'de hic yok. Bu tek basina hata
                                  degildir: tipin referansi o kayitla %100
                                  ortusuyorsa kurasyon dogru, derlemenin tablosu
                                  EKSIKTIR. Karar yine diziye gore verilir.
    """
    review_ids = {(r.get("pdb") or "").strip().upper() for r in review_rows}
    mapped_clusters = {m["mapped_cluster"] for _, _, m in triples if m["mapped_cluster"]}

    def identity(pdb, cluster):
        best = None
        for entity in (records.get(pdb, {}) or {}).get("entities", []):
            key = "%s_%s" % (pdb, entity["entity_id"])
            for hit in (identities or {}).get(key, []):
                if hit[0] == cluster and (best is None or hit[1] > best):
                    best = hit[1]
        return best

    out = []
    for chem in chemistry_rows:
        cluster = chem["cluster"]
        ids = sorted({t.upper() for t in
                      re.findall(r"\b\d[A-Za-z0-9]{3}\b", chem.get("pdb") or "")})
        for pdb in ids:
            own = identity(pdb, cluster)
            item = {"cluster": cluster, "curated_pdb": pdb,
                    "listed_in_review": pdb in review_ids,
                    "identity_to_this_types_reference_pct":
                        None if own is None else round(own, 2),
                    "mapped_by_a_review_row": cluster in mapped_clusters}
            if own is not None and own >= MIN_SEQ_IDENTITY:
                if pdb in review_ids:
                    item["status"] = "listed_and_matches"
                    item["verdict"] = ("the review lists this id and the structure "
                                       "matches this type's curated reference "
                                       "(%.1f%%); nothing to fix" % own)
                else:
                    # Derlemenin Tablo 1'inde EN YAKIN kayit hangisi?
                    closest, closest_pct = "", -1.0
                    for other in sorted(review_ids):
                        pct = identity(other, cluster)
                        if pct is not None and pct > closest_pct:
                            closest, closest_pct = other, pct
                    item["status"] = "absent_from_review"
                    item["closest_entry_in_table_1"] = closest
                    item["closest_entry_identity_pct"] = None if closest_pct < 0 \
                        else round(closest_pct, 2)
                    if closest_pct >= MIN_SEQ_IDENTITY:
                        # Derleme AYNI proteini baska bir kimlikle listeliyor:
                        # iki kimlik ayni proteinin farkli durumlari. Ayrintili
                        # karar curated_pdb_disagreements bolumunde.
                        item["same_protein_listed_under"] = closest
                        item["verdict"] = (
                            "this exact id is not in Table 1, but the review lists "
                            "the same protein under %s (%.1f%% to the same curated "
                            "reference). Both ids are correct; see "
                            "curated_pdb_disagreements for which states they are."
                            % (closest, closest_pct))
                    else:
                        item["verdict"] = (
                            "this id is NOT in Table 1, and no entry the review "
                            "does list comes close to this type's curated reference "
                            "(closest: %s at %.1f%%, against %.1f%% for the curated "
                            "id). The curated id is the right one for this type and "
                            "the review's table simply omits this structure -- this "
                            "is not a curation error."
                            % (closest or "none", max(closest_pct, 0.0), own))
            elif own is not None:
                item["status"] = "listed_but_mismatched" if pdb in review_ids \
                    else "mismatched_and_absent_from_review"
                item["verdict"] = (
                    "CURATION ERROR: this type's curated reference is only %.1f%% "
                    "identical to %s, so %s does not represent this type's "
                    "reference protein." % (own, pdb, pdb))
            else:
                item["status"] = "no_sequence_comparison"
                item["verdict"] = ("no blastp hit between this type's reference and "
                                   "this entry, so the assignment could not be "
                                   "tested by sequence")
            out.append(item)
    return out


def _state(record):
    ligands = record.get("ligands") or []
    if not ligands:
        return "no bound heteroatoms"
    return "ligands " + "/".join(ligands)


# --------------------------------------------------------------------------- #
# Mimari dogrulama -- derlemenin DENEYSEL etiketi vs cikarilan has_beta
# --------------------------------------------------------------------------- #

def architecture_validation(rows, db_path):
    """Derlemenin bilesim etiketini ro_etc.has_beta cogunlugu ile karsilastir."""
    out = {"what_is_compared": (
        "The review's subunit_composition column is an experimentally determined "
        "alpha3 versus alpha3beta3 label. This database instead INFERS the presence "
        "of a beta subunit from genomic context (ro_etc.has_beta), and the "
        "statistics page reports that inference as beta_by_type with Cramer V = "
        "0.94 without any outside validation. For every review row that maps onto a "
        "type, the review therefore supplies independent ground truth for that "
        "inference. Agreement is scored per type, not per row, because two review "
        "rows can land on the same type."),
        "query": ("SELECT r.ro_cluster, e.has_beta FROM ro r JOIN ro_etc e "
                  "USING(candidate_id) WHERE r.ro_cluster = ? AND r.is_confirmed = 1"),
        "rows": [], "by_type": [], "disagreements": [], "not_testable": []}
    if not os.path.exists(db_path):
        out["error"] = "roar.sqlite not found at %s" % db_path
        return out
    connection = sqlite3.connect("file:%s?mode=ro" % db_path, uri=True)
    try:
        per_type = {}
        for row, record, mapping in rows:
            cluster = mapping["mapped_cluster"]
            spec = COMPOSITION.get((row.get("subunit_composition") or "").strip())
            if not cluster or spec is None:
                out["not_testable"].append({
                    "ro_name": row.get("ro_name"), "pdb": row.get("pdb"),
                    "reason": "the row is not mapped to a type" if not cluster
                    else "unknown composition label"})
                continue
            counts = connection.execute(
                "SELECT SUM(CASE WHEN e.has_beta=1 THEN 1 ELSE 0 END), "
                "       SUM(CASE WHEN e.has_beta=0 THEN 1 ELSE 0 END), COUNT(*) "
                "FROM ro r JOIN ro_etc e USING(candidate_id) "
                "WHERE r.ro_cluster = ? AND r.is_confirmed = 1", (cluster,)).fetchone()
            yes, no, total = (counts[0] or 0), (counts[1] or 0), (counts[2] or 0)
            entry = {"ro_name": row.get("ro_name"), "pdb": row.get("pdb"),
                     "mapped_cluster": cluster,
                     "review_composition": row.get("subunit_composition"),
                     "review_has_beta": spec["beta"],
                     "confirmed_members": total, "has_beta_1": yes, "has_beta_0": no}
            if total == 0 or (yes + no) == 0:
                entry["db_has_beta_majority"] = None
                entry["agreement"] = "not_testable"
                entry["note"] = ("this type has no confirmed members in the ro "
                                 "table, so there is no inferred label to compare")
                out["not_testable"].append({
                    "ro_name": row.get("ro_name"), "pdb": row.get("pdb"),
                    "mapped_cluster": cluster,
                    "reason": "the type has no confirmed members"})
            else:
                fraction = yes / float(yes + no)
                majority = fraction >= 0.5
                entry["beta_fraction"] = round(fraction, 4)
                entry["db_has_beta_majority"] = majority
                entry["agreement"] = "agree" if majority == spec["beta"] else "disagree"
                if entry["agreement"] == "disagree":
                    entry["direction"] = (
                        "the database infers a beta subunit where the review's "
                        "structure shows none (DB over-calls beta)" if majority else
                        "the review's structure contains a beta subunit but the "
                        "database infers none (DB under-calls beta)")
                    out["disagreements"].append(entry)
            out["rows"].append(entry)
            if entry.get("agreement") in ("agree", "disagree"):
                per_type.setdefault(cluster, entry)
        out["by_type"] = sorted(per_type.values(), key=lambda e: e["mapped_cluster"])
        agree = sum(1 for e in out["by_type"] if e["agreement"] == "agree")
        total_types = len(out["by_type"])
        out["types_tested"] = total_types
        out["types_agree"] = agree
        out["types_disagree"] = total_types - agree
        out["agreement_rate"] = round(agree / float(total_types), 4) if total_types else None
        out["rows_tested"] = sum(1 for e in out["rows"]
                                 if e.get("agreement") in ("agree", "disagree"))
        out["rows_agree"] = sum(1 for e in out["rows"] if e.get("agreement") == "agree")
    finally:
        connection.close()
    return out


# --------------------------------------------------------------------------- #

def _strip_sequences(record):
    """Kaydin JSON'a yazilacak hali: ham diziler cikarilir, boyut patlamasin."""
    out = {k: v for k, v in record.items() if k != "entities"}
    out["entities"] = [{k: v for k, v in e.items() if k != "sequence"}
                       for e in record.get("entities", [])]
    return out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--chemistry", default="chemistry.csv")
    parser.add_argument("--in", dest="infile", default="reference_structures.csv",
                        help="derlemenin Tablo 1'i (yerinde genisletilir)")
    parser.add_argument("--out", default="analysis_out/reference_structures.json")
    parser.add_argument("--out-csv", dest="out_csv", default=None,
                        help="varsayilan: --in ile ayni dosya (yerinde yazar)")
    parser.add_argument("--offline", action="store_true",
                        help="yalnizca onbellekten oku, hic istek gonderme")
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--refs", default="ROs_71_Clean/refs71.fasta")
    parser.add_argument("--cache", default="structures",
                        help="RCSB yanitlarinin onbellek dizini")
    parser.add_argument("--pause", type=float, default=REQUEST_PAUSE_S,
                        help="istekler arasi bekleme (saniye)")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--completeness", default="analysis_out/completeness_audit.json",
                        help="yapisi hic olmayan tipleri buradan okur; yoksa atlanir")
    args = parser.parse_args()
    out_csv = args.out_csv or args.infile

    with open(args.infile, newline="") as handle:
        reader = csv.DictReader(handle)
        source_columns = [c for c in reader.fieldnames if c not in NEW_COLUMNS]
        review_rows = list(reader)
    chemistry_rows = list(csv.DictReader(open(args.chemistry, newline="")))
    chem_by_cluster = {r["cluster"]: r for r in chemistry_rows}
    refs = read_fasta(args.refs)
    refs_by_cluster = {cluster_of(h): (h, s) for h, s in refs.items()}
    print("[girdi] %d derleme satiri, %d kuratorlu tip, %d referans dizisi"
          % (len(review_rows), len(chemistry_rows), len(refs_by_cluster)))

    # --- 1. RCSB dogrulamasi ----------------------------------------------
    wanted = []
    for row in review_rows:
        pid = (row.get("pdb") or "").strip().upper()
        if pid and pid not in wanted:
            wanted.append(pid)
    # chemistry.csv'nin tasidigi kimlikler de gerekiyor: farklari kiyaslamak icin.
    for pid in sorted(pdb_index(chemistry_rows)):
        if pid not in wanted:
            wanted.append(pid)
    records = {}
    for pid in wanted:
        records[pid] = read_entry(pid, args.cache, args.offline, args.pause)
        state = "ok" if records[pid]["resolves"] else records[pid].get("error", "?")
        print("[rcsb] %-5s %-11s %s" % (pid, records[pid]["fetch_status"], state))

    # --- 2. dizi kimligi tablosu ------------------------------------------
    blocks = []
    for pid, record in records.items():
        for entity in record.get("entities", []):
            if entity.get("sequence"):
                blocks.append(">%s_%s\n%s" % (pid, entity["entity_id"],
                                              entity["sequence"]))
    identities = identity_table("\n".join(blocks) + "\n", args.refs, args.threads) \
        if blocks else None
    if identities is None:
        print("[uyari] blastp yok ya da calismadi: dizi kanitiyla eslestirme atlandi")
    else:
        print("[dizi] %d varlik dizisi %d referansa karsi blastp ile kiyaslandi"
              % (len(blocks), len(refs_by_cluster)))

    pdb_by_id = pdb_index(chemistry_rows)
    out_rows, triples = [], []
    for row in review_rows:
        pid = (row.get("pdb") or "").strip().upper()
        record = records.get(pid, {"pdb": pid, "resolves": False,
                                   "error": "not fetched"})
        mapping, map_notes = map_row(row, record, chem_by_cluster, pdb_by_id,
                                     identities, refs_by_cluster)
        verdict, agree_notes = check_agreement(row, record)
        notes = agree_notes + map_notes
        for alternative in mapping["alternatives"]:
            notes.append("also matches " + alternative)

        entities_text = "; ".join(
            "%s x%d (%d aa%s)" % (e.get("description") or "?", e.get("chain_count", 0),
                                  e.get("length", 0),
                                  "" if not e.get("reference_ids")
                                  else ", " + "/".join(e["reference_ids"]))
            for e in record.get("entities", []))
        citation = record.get("citation", {}) or {}
        doi = citation.get("doi") or ""
        if not doi and citation.get("pubmed"):
            doi = "PMID:" + citation["pubmed"]

        out = dict(row)
        out.update({
            "mapped_cluster": mapping["mapped_cluster"],
            "mapping_basis": mapping["mapping_basis"],
            "mapping_identity_pct": mapping["mapping_identity_pct"],
            "rcsb_title": record.get("title", ""),
            "rcsb_organism": "; ".join(record.get("organisms", [])),
            "rcsb_entities": entities_text,
            "rcsb_citation_doi": doi,
            "agrees_with_review": verdict,
            "mapping_note": " | ".join(n for n in notes if n),
        })
        out_rows.append(out)
        triples.append((row, record, mapping))

    # --- 3. kurasyon farklari ---------------------------------------------
    differences = compare_curated_ids(triples, records, chem_by_cluster, identities)
    curated_audit = curated_ids_vs_review(chemistry_rows, review_rows, triples,
                                          records, identities)

    # --- 4. mimari dogrulama ----------------------------------------------
    architecture = architecture_validation(triples, args.db)

    # --- 5. ozet ----------------------------------------------------------
    mapped = [r for r in out_rows if r["mapped_cluster"]]
    gains = []
    for row in out_rows:
        cluster = row["mapped_cluster"]
        if not cluster:
            continue
        curated = (chem_by_cluster.get(cluster, {}).get("pdb") or "").strip()
        if not curated:
            gains.append({"cluster": cluster, "gains_pdb": row["pdb"],
                          "ro_name": row["ro_name"],
                          "organism": row["rcsb_organism"],
                          "citation": row["rcsb_citation_doi"],
                          "basis": row["mapping_basis"],
                          "identity_pct": row["mapping_identity_pct"]})
    gained_clusters = sorted({g["cluster"] for g in gains})

    summary = {
        "review": {
            "citation": "Miao H, Schmidt S. Rieske Oxygenases: Powerful Models for "
                        "Understanding Nature's Orchestration of Electron Transfer "
                        "and Oxidative Chemistry. Biochemistry 2025;64:3801-3813",
            "doi": "10.1021/acs.biochem.5c00369",
            "licence": "CC BY 4.0", "table": "Table 1"},
        "rows_total": len(out_rows),
        "rows_resolved_at_rcsb": sum(1 for p in (r["pdb"].strip().upper()
                                                 for r in out_rows)
                                     if records.get(p, {}).get("resolves")),
        "rows_mapped": len(mapped),
        "rows_unmapped": len(out_rows) - len(mapped),
        "mapping_basis_counts": {
            basis: sum(1 for r in out_rows if r["mapping_basis"] == basis)
            for basis in ("pdb_match", "organism_and_name", "sequence_identity",
                          "unmapped")},
        "agreement_counts": {
            verdict: sum(1 for r in out_rows if r["agrees_with_review"] == verdict)
            for verdict in ("yes", "partial", "no")},
        "distinct_types_mapped": len({r["mapped_cluster"] for r in mapped}),
        "types_gaining_an_experimental_structure": gained_clusters,
        "types_gaining_detail": gains,
        "curated_pdb_disagreements": differences,
        "curated_pdb_disagreements_material": [d for d in differences if d["material"]],
        "curated_pdb_ids_vs_table_1": curated_audit,
        "curated_pdb_ids_that_are_curation_errors": [
            c for c in curated_audit
            if c["status"] in ("listed_but_mismatched",
                               "mismatched_and_absent_from_review")],
        "curated_pdb_ids_the_review_omits": [
            c for c in curated_audit if c["status"] == "absent_from_review"],
    }

    # completeness_audit.json "hicbir yapisi olmayan" tipleri sayiyor (ne deneysel
    # ne AlphaFold). Bu tablonun en keskin katkisi tam o listeyi kisaltmasidir.
    try:
        with open(args.completeness) as handle:
            no_structure = json.load(handle)["summary"]["structure"][
                "types_with_no_structure_at_all"]
        summary["types_with_no_structure_at_all_before"] = sorted(no_structure)
        summary["of_those_now_gaining_an_experimental_structure"] = sorted(
            set(no_structure) & set(gained_clusters))
        summary["of_those_still_without_any_structure"] = sorted(
            set(no_structure) - set(gained_clusters))
    except (OSError, ValueError, KeyError) as exc:
        summary["types_with_no_structure_at_all_before"] = None
        summary["no_structure_cross_reference_note"] = (
            "could not read %s (%s), so the cross-reference against the types that "
            "had no structure of any kind was skipped" % (args.completeness, exc))

    payload = {
        "module": os.path.basename(__file__),
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "what_was_measured": (
            "Every PDB id in Table 1 of Miao & Schmidt 2025 was resolved against the "
            "RCSB data API; each row was mapped onto one of this database's 71 "
            "curated enzyme types by PDB id, by name plus organism, or by sequence "
            "identity to the curated reference; and the review's experimentally "
            "determined subunit composition was compared against this database's "
            "genomic-context inference of beta-subunit presence."),
        "parameters": {"min_sequence_identity_pct": MIN_SEQ_IDENTITY,
                       "identity_definition": "nident / curated reference length",
                       "chemistry": args.chemistry, "input_csv": args.infile,
                       "db": args.db, "refs": args.refs, "cache": args.cache,
                       "offline": args.offline,
                       "blastp_available": identities is not None},
        "summary": summary,
        "architecture_validation": architecture,
        "rows": out_rows,
        "rcsb_records": {pid: _strip_sequences(rec)
                         for pid, rec in records.items()},
        "limitations": [
            "The composition check counts DISTINCT polymer entities, which is what "
            "the review's alpha3 versus alpha3beta3 label encodes. It does not "
            "verify the trimeric assembly itself: for entries whose asymmetric unit "
            "holds fewer than three chains the trimer is generated by "
            "crystallographic symmetry, and that is recorded as a note rather than a "
            "disagreement.",
            "The alpha3alpha'3 label describes a heterocomplex of two different "
            "catalytic subunits, but the two subunits were deposited as separate "
            "entries. Neither entry alone can confirm that composition.",
            "Sequence identity is computed against the single curated reference "
            "sequence of a type, not against that type's members, so it certifies "
            "that the structure represents the REFERENCE protein -- not that it "
            "represents every member assigned to the type.",
            "Where two curated types carry identical reference sequences, a "
            "structure matches both at the same identity. The script then prefers "
            "the type whose name matches the review's RO name and records the other "
            "as a duplicate; it does not decide which of the two names is correct.",
            "ro_etc.has_beta is itself an inference, and for a type with very few "
            "confirmed members its majority rests on one or two genomes. "
            "Disagreements carry the member count so that low-power cases are "
            "visible.",
        ],
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False, default=str)

    columns = source_columns + NEW_COLUMNS
    temporary = out_csv + ".part"
    with open(temporary, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns,
                                quoting=csv.QUOTE_ALL, lineterminator="\n")
        writer.writeheader()
        for row in out_rows:
            writer.writerow({c: row.get(c, "") for c in columns})
    os.replace(temporary, out_csv)

    print("\n[cikti] %s" % out_csv)
    print("[cikti] %s" % args.out)
    print("[ozet] %d/%d satir RCSB'de cozuldu, %d satir tipe baglandi "
          "(%d ayri tip), %d satir baglanamadi"
          % (summary["rows_resolved_at_rcsb"], summary["rows_total"],
             summary["rows_mapped"], summary["distinct_types_mapped"],
             summary["rows_unmapped"]))
    print("[ozet] uyum: %s" % summary["agreement_counts"])
    print("[ozet] deneysel yapi KAZANAN tip: %d -- %s"
          % (len(gained_clusters), ", ".join(gained_clusters)))
    closed = summary.get("of_those_now_gaining_an_experimental_structure")
    if closed is not None:
        print("[ozet] bunlardan %d'si HIC yapisi olmayan %d tipten: %s"
              % (len(closed), len(summary["types_with_no_structure_at_all_before"]),
                 ", ".join(closed) or "-"))
    if architecture.get("agreement_rate") is not None:
        print("[mimari] %d tip sinandi, %d uyuyor, %d uyusmuyor -- oran %.1f%%"
              % (architecture["types_tested"], architecture["types_agree"],
                 architecture["types_disagree"],
                 100.0 * architecture["agreement_rate"]))
    for item in architecture.get("disagreements", []):
        print("[mimari/UYUSMAZLIK] %s (%s, %s): derleme=%s, veritabani beta=%s "
              "(%d/%d uye) -- %s"
              % (item["mapped_cluster"], item["ro_name"], item["pdb"],
                 item["review_composition"], item["db_has_beta_majority"],
                 item["has_beta_1"], item["confirmed_members"], item["direction"]))
    for item in differences:
        print("[fark] %s: derleme=%s, chemistry.csv=%s -> %s (%s)"
              % (item["cluster"], item["review_pdb"], item["curated_pdb"],
                 item["verdict"], "MATERIAL" if item["material"] else "immaterial"))
    for item in summary["curated_pdb_ids_that_are_curation_errors"]:
        print("[KURASYON HATASI] %s -> %s: tipin referansi bu kayitla yalnizca "
              "%s%% ortusuyor" % (item["cluster"], item["curated_pdb"],
                                  item["identity_to_this_types_reference_pct"]))
    for item in summary["curated_pdb_ids_the_review_omits"]:
        if item.get("same_protein_listed_under"):
            print("[baska kimlikle] %s -> %s: derleme ayni proteini %s olarak "
                  "listeliyor (%s%%)"
                  % (item["cluster"], item["curated_pdb"],
                     item["same_protein_listed_under"],
                     item.get("closest_entry_identity_pct")))
        else:
            print("[derlemede yok] %s -> %s: tip referansina %s%%, Tablo 1'deki en "
                  "yakin kayit %s yalnizca %s%% -- kurasyon DOGRU, tablo eksik"
                  % (item["cluster"], item["curated_pdb"],
                     item["identity_to_this_types_reference_pct"],
                     item.get("closest_entry_in_table_1"),
                     item.get("closest_entry_identity_pct")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
