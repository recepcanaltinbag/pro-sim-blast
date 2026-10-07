"""
ROAR-DB -- 71 kuratorlu tipin EKSIKLIK DENETIMI (completeness audit).

NEDEN BU MODUL VAR
  Veritabani 71 kuratorlu RO tipi yayinliyor ama her tip ayni olcude
  belgelenmis DEGIL. Bir tipin yaninda kristal yapi, okunabilir bir
  reaksiyon ve izlenebilir bir kaynak varken, bir diger tipin yaninda
  yalnizca "established literature on the reference enzyme, no individual
  citation" yazili. Bu fark bugune kadar HICBIR yerde tek tabloda
  toplanmamisti: yapi kapsamini active_site.py, kanit kademesini
  evidence_tiers.py, kaynak alanlarini ise yalnizca chemistry.csv'yi acan
  insan goruyordu. Sonuc olarak "neyi tamamlamak gerekiyor" sorusunun cevabi
  her seferinde elle arastirilmak zorundaydi ve eksik tipler sessizce
  yayinda kaliyordu.

  Bu modul denetimi OLCUME cevirir: tip basina tek satir, ucu birbirinden
  bagimsiz uc eksen (yapi / kimya / literatur) ve her eksende "ne yapilirsa
  duzelir" etiketi. Modul SALT OKUNUR: chemistry.csv ve roar.sqlite'a
  dokunmaz, yapı dosyasi indirmez.

NE OLCULDU
  1. YAPI
     a. chemistry.csv'deki `pdb` alani var mi; varsa bicimi gecerli mi
        (4 karakter, ilk karakter 1-9 arasi bir rakam, kalani alfanumerik --
        RCSB'nin klasik kimlik bicimi). Bozuk kimlik AYRICA raporlanir,
        cunku bos alan ile hatali alan ayni sey degildir.
     b. Sitenin YERELDE gercekten neye sahip oldugu: structures/ altindaki
        deneysel dosya (<PDBID>.cif/.pdb) ve structures/alphafold/ altindaki
        ongorulen model (AF-<ACC>-F1-model_v<N>.pdb). Her tip ucden birine
        dusurulur: experimental / predicted_only / none.
     c. roar.sqlite icinde yapi tutan bir tablo olup olmadigi da ARANIR
        (adinda struct/pdb/model gecen tablolar). Bulunmazsa bu da bir
        olcumdur: yapi bilgisi veritabaninda degil, dosya sisteminde ve
        analysis_out/active_site.json icinde yasiyor.
     d. Deneysel yapisi olmayan tipler icin AlphaFold modeli ALINABILIR MI:
        referans UniProt kimligi active_site.json'un
        prediction_accession_mapping bloklarindan okunur (o blok kimligi
        UniParc'ta CRC64 tam-dizi eslesmesiyle bulur, yani tahmin degil) ve
        AlphaFold DB API'sine sorulur. Model INDIRILMEZ; yalnizca var/yok,
        kimlik, surum, ortalama pLDDT ve URL kaydedilir.

        BURADA BIR HATA YAKALANDI: UniParc kimlikleri bazi kayitlarda surum
        sonekiyle doner ("T1YXQ1.1"). AlphaFold API'si surum sonekli kimlige
        HTTP 400 verir, yani active_site.py'nin "modeli yok" dedigi bazi
        tiplerin modeli aslinda VAR. Bu modul her kimligi hem oldugu gibi hem
        de surum soneki atilmis haliyle dener ve hangi bicimin calistigini
        yazar.

  2. KIMYA / REAKSIYON
     Her tipte substrate_en, substrate_smiles, product_en, reaction,
     reaction_class alanlarinin dolu olup olmadigi alan alan raporlanir.
     SMILES dizgeleri ayrica ayristirilir. RDKit varsa gercek kimyasal
     ayristirma yapilir; YOKSA (bu makinede yok) yalnizca HAFIF bir sozdizim
     denetimi yapilir -- parantez dengesi, koseli parantez dengesi, halka
     kapanis rakamlarinin eslesmesi, atom simgelerinin izinli kumede olmasi.
     Bu HAFIF denetim kimyasal gecerlilik DEGILDIR (valans, aromatiklik,
     yuk dengesi sinanmaz) ve ciktida acikca boyle etiketlenir.
     Uyesi olmayan tipler ve kimya satiri olmayan tipler ayri ayri listelenir.

  3. LITERATUR / ATIF
     source, source_kind, curation_confidence okunur ve `source` metni
     IZLENEBILIR bir kimlik tasiyor mu diye taranir: DOI, PMID, PMCID,
     dergi cilt:sayfa atifi, dergi makale numarasi, PDB kimligi. Buradan dort
     sinif cikar:
        resolvable_identifier  DOI/PMID/PMCID var -- okur dogrudan ulasir
        journal_citation       cilt:sayfa ya da makale no var -- ulasilabilir
        structure_only         tek kanit bir PDB kimligi, makale yok
        vague_free_text        serbest metin, okur takip edemez
     "paper" ile "structure" kanitlarinin ayrimi source_kind'dan DEGIL,
     metinden de dogrulanir; ikisi celisirse celiski raporlanir.

NE OLCULMEDI
  - Kaynagin DOGRULUGU. Bir DOI'nin var olmasi, o makalenin o enzimi
    gercekten tanimladigi anlamina gelmez; bu modul erisilebilirligi olcer,
    icerigi degil. DOI/PMID canli olarak sorgulanmaz (ag kullanimi yalnizca
    AlphaFold API'si ile sinirli tutuldu).
  - SMILES'in kimyasal dogrulugu (RDKit yoksa). Bkz. yukarisi.
  - product_smiles diye bir kolon chemistry.csv'de HIC YOK: urun yalnizca
    ingilizce metin olarak tutuluyor. Bu, tek bir tipin eksigi degil, semanin
    kendisinin eksigi oldugu icin ozet blogunda ayrica bildirilir.

Cikti: analysis_out/completeness_audit.json  (makine okunur, tip basina bir nesne)
       analysis_out/completeness_audit.md    (ne yapilirsa duzelir, oncelikli)
"""

import argparse
import csv
import json
import os
import re
import sqlite3
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone

ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{0}"

# chemistry.csv'de bir tipin "reaksiyonu tanimlanmis" sayilmasi icin dolu
# olmasi gereken alanlar. family/notes bilerek DISARIDA: family bir
# siniflandirma, notes ise serbest yorum alani.
REQUIRED_CHEMISTRY_FIELDS = ["substrate_en", "substrate_smiles", "product_en",
                             "reaction", "reaction_class"]

# RCSB klasik kimlik bicimi: 4 karakter, ilki 1-9 bir rakam.
PDB_ID_RE = re.compile(r"^[1-9][A-Za-z0-9]{3}$")

# Yapi dosyasi olarak kabul edilen uzantilar (structures/ altinda).
STRUCTURE_SUFFIXES = [".cif", ".pdb", ".ent", ".mmcif"]

# --- izlenebilir kimlik desenleri -------------------------------------------
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s,;)\]]+")
PMID_RE = re.compile(r"PMID[:\s]*(\d{4,9})", re.I)
PMCID_RE = re.compile(r"PMC(\d{4,9})", re.I)
# "J Bacteriol 1997;179:3549" ya da "2002;184:509-518"
JOURNAL_VOLPAGE_RE = re.compile(r"\b(19|20)\d{2}\s*;\s*\d+\s*:\s*\d+")
# "J Bacteriol 2025 jb.00221-25" -- cilt:sayfa yok, makale numarasi var
ARTICLE_NO_RE = re.compile(r"\b[a-z]{2,6}\.\d{3,6}-\d{2,4}\b")
PDB_IN_TEXT_RE = re.compile(r"\bPDB\s+([1-9][A-Za-z0-9]{3})", re.I)

# --- HAFIF SMILES denetimi ---------------------------------------------------
# Yalnizca sozdizim. Iki harfli simgeler once denenir, yoksa "Cl" -> "C"+"l".
TWO_LETTER_ORGANIC = ["Cl", "Br"]
ONE_LETTER_ORGANIC = ["B", "C", "N", "O", "P", "S", "F", "I"]
AROMATIC_LOWER = ["c", "n", "o", "p", "s", "b"]
AROMATIC_TWO_LOWER = ["se", "as"]
BOND_CHARS = set("-=#$:/\\~")
ELEMENTS = set("""H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe
Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe
Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb
Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr""".split())


def smiles_syntax_check(smiles):
    """SMILES'i HAFIF denetle: parantez/koseli parantez/halka dengesi + atomlar.

    Bu bir KIMYASAL dogrulama degildir. Valans, aromatik halka gecerliligi,
    yuk dengesi ve stereo tutarliligi SINANMAZ. RDKit varsa cagiran taraf
    gercek ayristirmayi tercih eder; bu fonksiyon yalnizca yedek yoldur.
    """
    problems = []
    atoms = []
    depth = 0
    ring_hits = Counter()
    i = 0
    n = len(smiles)
    prev_token = None
    while i < n:
        ch = smiles[i]
        if ch == "[":
            close = smiles.find("]", i)
            if close < 0:
                problems.append("unclosed square bracket at position %d" % i)
                break
            inner = smiles[i + 1:close]
            symbol = re.match(r"\d*([A-Za-z][a-z]?)", inner)
            if not symbol:
                problems.append("bracket atom [%s] carries no element symbol" % inner)
            else:
                sym = symbol.group(1)
                if sym not in ELEMENTS and sym.capitalize() not in ELEMENTS:
                    problems.append("unknown element symbol in [%s]" % inner)
            atoms.append(inner)
            prev_token = "atom"
            i = close + 1
            continue
        if ch == "(":
            depth += 1
            if i + 1 < n and smiles[i + 1] == ")":
                problems.append("empty branch () at position %d" % i)
            prev_token = "open"
            i += 1
            continue
        if ch == ")":
            depth -= 1
            if depth < 0:
                problems.append("closing parenthesis without an opening one at position %d" % i)
                depth = 0
            prev_token = "close"
            i += 1
            continue
        if ch == "%":
            num = smiles[i + 1:i + 3]
            if not num.isdigit() or len(num) != 2:
                problems.append("%%-ring bond at position %d is not followed by two digits" % i)
                i += 1
                continue
            ring_hits[num] += 1
            prev_token = "ring"
            i += 3
            continue
        if ch.isdigit():
            if prev_token is None:
                problems.append("ring-closure digit at position %d has no atom to attach to" % i)
            ring_hits[ch] += 1
            prev_token = "ring"
            i += 1
            continue
        if ch in BOND_CHARS:
            prev_token = "bond"
            i += 1
            continue
        if ch in ".@+-":
            prev_token = "misc"
            i += 1
            continue
        two = smiles[i:i + 2]
        if two in TWO_LETTER_ORGANIC or two in AROMATIC_TWO_LOWER:
            atoms.append(two)
            prev_token = "atom"
            i += 2
            continue
        if ch in ONE_LETTER_ORGANIC or ch in AROMATIC_LOWER or ch == "*":
            atoms.append(ch)
            prev_token = "atom"
            i += 1
            continue
        problems.append("character %r at position %d is not legal SMILES" % (ch, i))
        i += 1
    if depth != 0:
        problems.append("%d parenthesis/-es left unclosed" % depth)
    unpaired = sorted(k for k, v in ring_hits.items() if v % 2)
    if unpaired:
        problems.append("ring-closure label(s) %s appear an odd number of times"
                        % ", ".join(unpaired))
    if not atoms:
        problems.append("no atom found")
    if smiles and smiles[0] in BOND_CHARS:
        problems.append("string starts with a bond symbol")
    if smiles and smiles[-1] in BOND_CHARS:
        problems.append("string ends with a bond symbol")
    return {"ok": not problems, "n_atom_tokens": len(atoms), "problems": problems}


def smiles_check(smiles, rdkit_mol_from_smiles):
    """RDKit varsa gercek ayristirma, yoksa hafif sozdizim denetimi."""
    if not (smiles or "").strip():
        return {"present": False, "ok": False, "method": "none",
                "problems": ["substrate_smiles is empty"]}
    if rdkit_mol_from_smiles is not None:
        mol = rdkit_mol_from_smiles(smiles)
        if mol is None:
            return {"present": True, "ok": False, "method": "rdkit",
                    "problems": ["RDKit could not parse the SMILES"]}
        return {"present": True, "ok": True, "method": "rdkit",
                "n_atoms": mol.GetNumAtoms(), "problems": []}
    res = smiles_syntax_check(smiles)
    res.update(present=True, method="lightweight_syntax_only")
    return res


# --- yapi tarafi -------------------------------------------------------------
def pdb_id_report(raw):
    """chemistry.csv'deki pdb alanini degerlendirir."""
    value = (raw or "").strip()
    if not value:
        return {"present": False, "value": None, "wellformed": False, "problems": []}
    problems = []
    if len(value) != 4:
        problems.append("PDB identifier %r is %d characters, not 4" % (value, len(value)))
    if not PDB_ID_RE.match(value):
        problems.append("PDB identifier %r does not match the RCSB pattern "
                        "[1-9] followed by three alphanumerics" % value)
    return {"present": True, "value": value.upper(),
            "wellformed": not problems, "problems": problems}


def find_experimental_file(structures_dir, pdb_id):
    """structures/ altinda bu PDB kimligine ait yerel dosyayi arar."""
    if not pdb_id:
        return None
    for suffix in STRUCTURE_SUFFIXES:
        for name in (pdb_id.upper() + suffix, pdb_id.lower() + suffix):
            path = os.path.join(structures_dir, name)
            if os.path.exists(path):
                return path
    return None


def find_alphafold_file(structures_dir, accession, version):
    """structures/alphafold/ altinda indirilmis modeli arar."""
    if not accession:
        return None
    folder = os.path.join(structures_dir, "alphafold")
    if version:
        path = os.path.join(folder, "AF-%s-F1-model_v%s.pdb" % (accession, version))
        if os.path.exists(path):
            return path
    if not os.path.isdir(folder):
        return None
    prefix = "AF-%s-F1-model" % accession
    for name in sorted(os.listdir(folder)):
        if name.startswith(prefix) and name.endswith(".pdb"):
            return os.path.join(folder, name)
    return None


def structure_tables_in_db(con):
    """roar.sqlite'da adinda struct/pdb/model gecen tablolari arar."""
    rows = con.execute("SELECT name FROM sqlite_master WHERE type IN ('table','view')").fetchall()
    hits = []
    for (name,) in rows:
        low = name.lower()
        if "struct" in low or "pdb" in low or "model" in low:
            hits.append(name)
    return sorted(hits)


def strip_accession_version(accession):
    """UniParc bazi kimlikleri surum sonekiyle dondurur: 'T1YXQ1.1' -> 'T1YXQ1'."""
    return accession.split(".")[0] if accession else accession


def read_cached_alphafold(cache_dir, accession):
    """active_site.py'nin BIRAKTIGI onbellegi SALT OKUNUR kullanir."""
    path = os.path.join(cache_dir, accession + ".json")
    if not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (ValueError, OSError):
        return None


def query_alphafold(accession, timeout):
    """AlphaFold DB API'sine tek bir GET. Model INDIRILMEZ, meta veri okunur."""
    url = ALPHAFOLD_API_URL.format(accession)
    try:
        with urllib.request.urlopen(url, timeout=timeout) as fh:
            payload = json.loads(fh.read().decode("utf-8", "replace"))
            return payload, fh.getcode()
    except urllib.error.HTTPError as exc:
        return None, exc.code
    except Exception as exc:                                  # ag/DNS/zaman asimi
        return None, "error: %s" % exc


def alphafold_record(payload):
    if isinstance(payload, list) and payload:
        rec = payload[0]
        if isinstance(rec, dict):
            return rec
    if isinstance(payload, dict) and payload.get("uniprotAccession"):
        return payload
    return None


def probe_alphafold(accessions, cache_dir, timeout, sleep_s, probed):
    """Verilen kimlikleri hem OLDUGU GIBI hem SURUM SONEKI ATILMIS dene.

    `probed` cagrilar arasi paylasilan bir sozluk: ayni kimlik iki kez
    sorulmaz (hem kibarlik hem belirlenimlilik icin).
    """
    attempts = []
    for accession in accessions:
        for candidate in [accession, strip_accession_version(accession)]:
            if not candidate or any(a["accession"] == candidate for a in attempts):
                continue
            if candidate in probed:
                attempts.append(dict(probed[candidate]))
                continue
            cached = read_cached_alphafold(cache_dir, candidate)
            if cached is not None:
                rec = alphafold_record(cached)
                attempt = {"accession": candidate, "how": "local cache written by active_site.py",
                           "http_status": None, "model_exists": rec is not None}
            else:
                if sleep_s and probed:
                    time.sleep(sleep_s)
                payload, status = query_alphafold(candidate, timeout)
                rec = alphafold_record(payload)
                attempt = {"accession": candidate, "how": "AlphaFold DB API",
                           "http_status": status, "model_exists": rec is not None}
            if rec is not None:
                attempt.update({
                    "alphafold_version": rec.get("latestVersion"),
                    "global_mean_plddt": rec.get("globalMetricValue"),
                    "model_url": rec.get("pdbUrl"),
                    "model_sequence_length": len(rec.get("uniprotSequence") or "") or None,
                })
            probed[candidate] = attempt
            attempts.append(dict(attempt))
    return attempts


# --- literatur tarafi --------------------------------------------------------
def source_report(source, source_kind):
    """`source` metnini izlenebilir kimlikler icin tarar ve siniflandirir."""
    text = (source or "").strip()
    ids = {
        "doi": sorted(set(DOI_RE.findall(text))),
        "pmid": sorted(set(PMID_RE.findall(text))),
        "pmcid": ["PMC" + x for x in sorted(set(PMCID_RE.findall(text)))],
        "journal_volume_page": [m.group(0) for m in JOURNAL_VOLPAGE_RE.finditer(text)],
        "journal_article_number": sorted(set(ARTICLE_NO_RE.findall(text))),
        "pdb": sorted({x.upper() for x in PDB_IN_TEXT_RE.findall(text)}),
    }
    kind = (source_kind or "").strip().lower()
    kind_says_paper = "paper" in kind
    kind_says_structure = "structure" in kind
    text_has_paper = bool(ids["doi"] or ids["pmid"] or ids["pmcid"]
                          or ids["journal_volume_page"] or ids["journal_article_number"])
    text_has_structure = bool(ids["pdb"])

    if ids["doi"] or ids["pmid"] or ids["pmcid"]:
        classification = "resolvable_identifier"
    elif ids["journal_volume_page"] or ids["journal_article_number"]:
        classification = "journal_citation"
    elif text_has_structure:
        classification = "structure_only"
    else:
        classification = "vague_free_text"

    conflicts = []
    if kind_says_paper and not text_has_paper:
        conflicts.append("source_kind claims a paper but the source text carries no "
                         "DOI, PMID, PMCID or journal citation")
    if text_has_paper and not kind_says_paper:
        conflicts.append("the source text carries a paper citation but source_kind "
                         "does not say 'paper'")
    if kind_says_structure and not text_has_structure:
        conflicts.append("source_kind claims a structure but the source text names no PDB id")
    return {
        "source": source,
        "source_kind": source_kind,
        "identifiers": ids,
        "has_paper_evidence": text_has_paper,
        "has_structure_evidence": text_has_structure,
        "classification": classification,
        "citable": classification in ("resolvable_identifier", "journal_citation"),
        "source_kind_conflicts": conflicts,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--chemistry", default="chemistry.csv")
    ap.add_argument("--structures", default="structures",
                    help="yerel yapi onbellegi (SALT OKUNUR)")
    ap.add_argument("--active-site", default="analysis_out/active_site.json",
                    help="referans UniProt kimligi haritasinin kaynagi")
    ap.add_argument("--out", default="analysis_out/completeness_audit.json")
    ap.add_argument("--out-md", default="analysis_out/completeness_audit.md")
    ap.add_argument("--alphafold-probe", choices=["none", "needed", "all"], default="needed",
                    help="'needed': yalnizca deneysel yapisi olmayan ve kimligi "
                         "cozulememis tipler icin AlphaFold API'sine sorar")
    ap.add_argument("--alphafold-timeout", type=float, default=45.0)
    ap.add_argument("--alphafold-sleep", type=float, default=2.0,
                    help="API cagrilari arasi bekleme (saniye) -- kibarlik")
    args = ap.parse_args()

    try:
        from rdkit import Chem
        rdkit_parse = Chem.MolFromSmiles
        rdkit_available = True
    except Exception:
        rdkit_parse = None
        rdkit_available = False

    # --- kuratorlu kimya satirlari ---
    with open(args.chemistry, newline="") as fh:
        chem_rows = list(csv.DictReader(fh))
    chem = {r["cluster"]: r for r in chem_rows}

    # --- uyeler ---
    con = sqlite3.connect(args.db)
    members = {}
    for cluster, total, confirmed in con.execute(
            "SELECT ro_cluster, COUNT(*), SUM(CASE WHEN is_confirmed=1 THEN 1 ELSE 0 END) "
            "FROM ro WHERE ro_cluster IS NOT NULL AND ro_cluster NOT IN ('','N/A') "
            "GROUP BY ro_cluster"):
        members[cluster] = {"rows": total, "confirmed": confirmed or 0}
    struct_tables = structure_tables_in_db(con)
    con.close()

    # --- referans kimlik haritasi (active_site.py uretti) ---
    accession_map = {}
    active_site_present = os.path.exists(args.active_site)
    if active_site_present:
        with open(args.active_site) as fh:
            raw = json.load(fh)
        for entry in raw.get("prediction_accession_mapping") or []:
            accession_map[entry.get("type")] = entry

    all_types = sorted(set(chem) | set(members))
    af_cache_dir = os.path.join(args.structures, "alphafold")
    probed = {}
    rows = []

    for name in all_types:
        row = chem.get(name)
        mem = members.get(name, {"rows": 0, "confirmed": 0})
        amap = accession_map.get(name, {})

        # ---------------- 1. YAPI ----------------
        pdb = pdb_id_report(row.get("pdb") if row else None)
        exp_file = find_experimental_file(args.structures, pdb["value"])
        declared_accession = amap.get("uniprot_accession")
        tried = list(amap.get("accessions_tried") or [])
        af_version = amap.get("alphafold_version")
        af_file = find_alphafold_file(args.structures, declared_accession, af_version)

        af = {
            "reference_accession": declared_accession,
            "accession_source": ("analysis_out/active_site.json -- UniParc CRC64 exact "
                                 "sequence match, not a similarity guess"),
            "mapping_status": amap.get("status"),
            "mapping_problems": list(amap.get("mapping_problems") or amap.get("problems") or []),
            "uniparc_id": amap.get("uniparc_id"),
            "accessions_tried_by_active_site": tried,
            "alphafold_version": af_version,
            "global_mean_plddt": amap.get("global_mean_plddt"),
            "model_url": amap.get("model_url"),
            "model_downloaded_locally": af_file,
            "model_exists": bool(declared_accession),
            "probe": [],
        }

        needs_probe = (exp_file is None and not declared_accession and tried
                       and args.alphafold_probe in ("needed", "all"))
        if args.alphafold_probe == "all" and tried:
            needs_probe = True
        if needs_probe:
            af["probe"] = probe_alphafold(tried, af_cache_dir, args.alphafold_timeout,
                                          args.alphafold_sleep, probed)
            hit = next((a for a in af["probe"] if a.get("model_exists")), None)
            if hit:
                af.update({
                    "model_exists": True,
                    "reference_accession": hit["accession"],
                    "alphafold_version": hit.get("alphafold_version"),
                    "global_mean_plddt": hit.get("global_mean_plddt"),
                    "model_url": hit.get("model_url"),
                    "recovered_by_stripping_version_suffix":
                        hit["accession"] not in tried,
                })

        if exp_file:
            local_state = "experimental"
        elif af_file:
            local_state = "predicted_only"
        else:
            local_state = "none"

        structure = {
            "pdb_field": pdb,
            "experimental_file_local": exp_file,
            "has_experimental_structure": exp_file is not None,
            "pdb_declared_but_file_absent": pdb["present"] and exp_file is None,
            "alphafold": af,
            "local_state": local_state,
        }

        # ---------------- 2. KIMYA ----------------
        if row is None:
            chemistry = {"has_chemistry_row": False,
                         "missing_fields": list(REQUIRED_CHEMISTRY_FIELDS),
                         "smiles_check": {"present": False, "ok": False, "method": "none",
                                          "problems": ["no chemistry.csv row at all"]},
                         "fields": {}}
        else:
            missing = [f for f in REQUIRED_CHEMISTRY_FIELDS if not (row.get(f) or "").strip()]
            chemistry = {
                "has_chemistry_row": True,
                "missing_fields": missing,
                "fields": {f: (row.get(f) or "").strip() or None
                           for f in REQUIRED_CHEMISTRY_FIELDS},
                "family": (row.get("family") or "").strip() or None,
                "notes_present": bool((row.get("notes") or "").strip()),
                "smiles_check": smiles_check(row.get("substrate_smiles"), rdkit_parse),
            }

        # ---------------- 3. LITERATUR ----------------
        if row is None:
            literature = {"source": None, "source_kind": None, "curation_confidence": None,
                          "classification": "missing", "citable": False,
                          "identifiers": {}, "has_paper_evidence": False,
                          "has_structure_evidence": False, "source_kind_conflicts": []}
        else:
            literature = source_report(row.get("source"), row.get("source_kind"))
            literature["curation_confidence"] = (row.get("curation_confidence") or "").strip() or None

        # ---------------- ne yapilirsa duzelir ----------------
        actions = []
        if not chemistry["has_chemistry_row"]:
            actions.append("needs a chemistry.csv row")
        for field in chemistry["missing_fields"]:
            actions.append("needs a %s" % field)
        if chemistry["smiles_check"].get("present") and not chemistry["smiles_check"]["ok"]:
            actions.append("needs a corrected substrate_smiles")
        if structure["pdb_declared_but_file_absent"]:
            actions.append("needs the declared PDB entry cached locally")
        if pdb["present"] and not pdb["wellformed"]:
            actions.append("needs a corrected pdb identifier")
        if not structure["has_experimental_structure"]:
            if af.get("recovered_by_stripping_version_suffix"):
                actions.append("needs the AlphaFold accession fixed "
                               "(drop the UniParc version suffix) and the model cached")
            elif af_file:
                pass                      # ongorulen model zaten yerelde
            elif af.get("model_exists"):
                actions.append("needs the AlphaFold model cached")
            elif tried:
                actions.append("no AlphaFold model exists for the mapped accession -- "
                               "needs an experimental structure or another accession")
            else:
                actions.append("needs a UniProt accession before any AlphaFold "
                               "model can be reached")
        if literature["classification"] == "vague_free_text":
            actions.append("needs a citable source (DOI or PMID)")
        elif literature["classification"] == "structure_only":
            actions.append("needs a paper citation (only a PDB entry is cited)")
        if mem["confirmed"] == 0:
            actions.append("has no confirmed members in the database")
        if literature.get("source_kind_conflicts"):
            actions.append("needs source_kind reconciled with the source text")

        rows.append({
            "type": name,
            "in_chemistry_csv": row is not None,
            "members": {"rows_in_ro": mem["rows"], "confirmed": mem["confirmed"],
                        "has_confirmed_members": mem["confirmed"] > 0},
            "structure": structure,
            "chemistry": chemistry,
            "literature": literature,
            "complete": not actions,
            "actions": actions,
        })

    # ---------------- ozet ----------------
    def names(pred):
        return [r["type"] for r in rows if pred(r)]

    field_gaps = {f: names(lambda r, f=f: f in r["chemistry"]["missing_fields"])
                  for f in REQUIRED_CHEMISTRY_FIELDS}
    action_index = defaultdict(list)
    for r in rows:
        for a in r["actions"]:
            action_index[a].append(r["type"])

    summary = {
        "types_audited": len(rows),
        "types_in_chemistry_csv": len(chem),
        "types_with_members_in_db": len(members),
        "types_complete_on_all_three_axes": len(names(lambda r: r["complete"])),
        "structure": {
            "types_with_a_pdb_field": len(names(lambda r: r["structure"]["pdb_field"]["present"])),
            "distinct_pdb_entries": len({r["structure"]["pdb_field"]["value"] for r in rows
                                         if r["structure"]["pdb_field"]["value"]}),
            "malformed_pdb_identifiers": names(
                lambda r: r["structure"]["pdb_field"]["present"]
                and not r["structure"]["pdb_field"]["wellformed"]),
            "pdb_declared_but_file_absent": names(
                lambda r: r["structure"]["pdb_declared_but_file_absent"]),
            "local_state_counts": dict(Counter(r["structure"]["local_state"] for r in rows)),
            "types_with_experimental_structure": names(
                lambda r: r["structure"]["local_state"] == "experimental"),
            "types_with_only_a_predicted_model": names(
                lambda r: r["structure"]["local_state"] == "predicted_only"),
            "types_with_no_structure_at_all": names(
                lambda r: r["structure"]["local_state"] == "none"),
            "alphafold_model_recoverable_by_fixing_the_accession": names(
                lambda r: r["structure"]["alphafold"].get(
                    "recovered_by_stripping_version_suffix")),
            "no_alphafold_model_for_the_mapped_accession": names(
                lambda r: r["structure"]["local_state"] == "none"
                and r["structure"]["alphafold"]["accessions_tried_by_active_site"]
                and not r["structure"]["alphafold"]["model_exists"]),
            "no_uniprot_accession_at_all": names(
                lambda r: r["structure"]["local_state"] == "none"
                and not r["structure"]["alphafold"]["accessions_tried_by_active_site"]),
            "structure_tables_found_in_sqlite": struct_tables,
            "note": ("roar.sqlite carries no structure table; structure state lives in the "
                     "structures/ cache and in analysis_out/active_site.json"
                     if not struct_tables else
                     "structure-like tables were found and are listed above"),
        },
        "chemistry": {
            "types_with_a_chemistry_row": len(names(lambda r: r["in_chemistry_csv"])),
            "types_with_members_but_no_chemistry_row": names(
                lambda r: not r["in_chemistry_csv"] and r["members"]["confirmed"] > 0),
            "types_in_chemistry_csv_without_confirmed_members": names(
                lambda r: r["in_chemistry_csv"] and r["members"]["confirmed"] == 0),
            "per_field_gaps": field_gaps,
            "types_with_every_required_field": len(
                names(lambda r: r["in_chemistry_csv"] and not r["chemistry"]["missing_fields"])),
            "smiles_validation_method": ("rdkit" if rdkit_available
                                         else "lightweight_syntax_only"),
            "smiles_validation_caveat": (
                "RDKit is installed, SMILES were parsed chemically" if rdkit_available else
                "RDKit is NOT installed on this machine. The SMILES were only checked "
                "for syntax: balanced parentheses and brackets, paired ring-closure "
                "labels and legal atom symbols. This is NOT a chemical validation -- "
                "valence, aromaticity and charge balance were not tested."),
            "smiles_present": len(names(lambda r: r["chemistry"]["smiles_check"].get("present"))),
            "smiles_failing_the_check": names(
                lambda r: r["chemistry"]["smiles_check"].get("present")
                and not r["chemistry"]["smiles_check"]["ok"]),
            "schema_level_gap": ("chemistry.csv has no product_smiles column at all, so for "
                                 "all %d types the product is text only and no reaction can be "
                                 "balanced or drawn mechanically" % len(chem)),
        },
        "literature": {
            "classification_counts": dict(Counter(r["literature"]["classification"]
                                                  for r in rows)),
            "source_kind_counts": dict(Counter(r["literature"].get("source_kind")
                                               for r in rows)),
            "curation_confidence_counts": dict(Counter(r["literature"].get("curation_confidence")
                                                       for r in rows)),
            "types_with_a_paper": names(lambda r: r["literature"]["has_paper_evidence"]),
            "types_with_only_a_structure_as_evidence": names(
                lambda r: r["literature"]["classification"] == "structure_only"),
            "types_with_neither_a_paper_nor_a_structure": names(
                lambda r: r["literature"]["classification"] == "vague_free_text"),
            "types_with_a_resolvable_identifier": names(
                lambda r: r["literature"]["classification"] == "resolvable_identifier"),
            "types_with_a_journal_citation_but_no_doi_or_pmid": names(
                lambda r: r["literature"]["classification"] == "journal_citation"),
            "types_with_source_kind_conflicts": names(
                lambda r: bool(r["literature"].get("source_kind_conflicts"))),
        },
        "actions_by_what_would_fix_it": {k: sorted(v) for k, v in
                                         sorted(action_index.items(),
                                                key=lambda kv: (-len(kv[1]), kv[0]))},
    }

    out = {
        "module": "completeness_audit.py",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "what_was_measured": (
            "For every curated enzyme type -- the union of chemistry.csv cluster values and "
            "the distinct ro.ro_cluster values -- whether it has a 3D structure (experimental "
            "or a cached AlphaFold model), a complete reaction/chemistry description, and a "
            "literature source a reader could follow up."),
        "parameters": {
            "db": args.db, "chemistry": args.chemistry, "structures": args.structures,
            "active_site": args.active_site, "alphafold_probe": args.alphafold_probe,
            "alphafold_sleep_s": args.alphafold_sleep,
            "required_chemistry_fields": REQUIRED_CHEMISTRY_FIELDS,
            "rdkit_available": rdkit_available,
            "active_site_json_present": active_site_present,
        },
        "method_notes": [
            "Read-only: chemistry.csv and roar.sqlite are never written, and no structure "
            "file is downloaded. The AlphaFold DB API is queried for metadata only.",
            "The reference UniProt accession is not guessed from similarity. It comes from "
            "active_site.py, which finds it by CRC64 exact-sequence lookup in UniParc.",
            "UniParc returns some accessions with a version suffix (T1YXQ1.1). The AlphaFold "
            "API answers HTTP 400 for those, so every accession is tried both as given and "
            "with the suffix stripped, and the attempt log records which form worked.",
            "A type counts as having a structure the site can show only when the file is "
            "actually on disk under structures/; a PDB id in chemistry.csv is a claim, not "
            "a cached structure.",
        ],
        "limitations": [
            "Source accessibility is measured, not source correctness: a DOI that exists is "
            "not evidence that the paper characterises this enzyme.",
            "DOIs, PMIDs and PMCIDs were not resolved over the network.",
            ("SMILES were only syntax-checked because RDKit is unavailable here; this is not "
             "a chemical validation." if not rdkit_available else
             "SMILES were parsed with RDKit."),
            "'no AlphaFold model exists' means the AlphaFold DB has no entry for the mapped "
            "accession on the day of the run; the database grows, so this can change.",
        ],
        "summary": summary,
        "types": rows,
    }

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=False)
    print("[yazildi] %s  (%d tip)" % (args.out, len(rows)))

    write_markdown(args.out_md, out)
    print("[yazildi] %s" % args.out_md)
    print_summary(out)


def write_markdown(path, out):
    """Insan okunur ozet: ne yapilirsa duzelir, oncelik sirasina gore."""
    s = out["summary"]
    rows = out["types"]
    by_type = {r["type"]: r for r in rows}

    def listing(items):
        return ", ".join("`%s`" % x for x in items) if items else "_yok / none_"

    lines = []
    lines.append("# ROAR-DB completeness audit")
    lines.append("")
    lines.append("Uretildi / generated: %s -- `completeness_audit.py`" % out["generated_utc"])
    lines.append("Makine okunur surum: `analysis_out/completeness_audit.json`")
    lines.append("")
    lines.append("## Nerede duruyoruz")
    lines.append("")
    lines.append("| olcu | deger |")
    lines.append("| --- | --- |")
    lines.append("| denetlenen tip | %d |" % s["types_audited"])
    lines.append("| uc eksende de eksigi olmayan tip | %d |"
                 % s["types_complete_on_all_three_axes"])
    lines.append("| deneysel yapisi yerelde olan | %d |"
                 % len(s["structure"]["types_with_experimental_structure"]))
    lines.append("| yalnizca ongorulen (AlphaFold) modeli olan | %d |"
                 % len(s["structure"]["types_with_only_a_predicted_model"]))
    lines.append("| hic yapisi olmayan | %d |"
                 % len(s["structure"]["types_with_no_structure_at_all"]))
    lines.append("| izlenebilir kimlikli kaynagi olan (DOI/PMID/PMCID) | %d |"
                 % len(s["literature"]["types_with_a_resolvable_identifier"]))
    lines.append("| dergi atifi olan (DOI/PMID yok) | %d |"
                 % len(s["literature"]["types_with_a_journal_citation_but_no_doi_or_pmid"]))
    lines.append("| tek kanidi bir PDB kaydi olan | %d |"
                 % len(s["literature"]["types_with_only_a_structure_as_evidence"]))
    lines.append("| ne makale ne yapi atfi olan | %d |"
                 % len(s["literature"]["types_with_neither_a_paper_nor_a_structure"]))
    lines.append("| dogrulanmis uyesi olmayan | %d |"
                 % len(s["chemistry"]["types_in_chemistry_csv_without_confirmed_members"]))
    lines.append("")

    lines.append("## Oncelikli is listesi -- ne yapilirsa duzelir")
    lines.append("")
    priority = [
        ("1. Hemen duzelebilir: AlphaFold kimligi bozuk",
         "needs the AlphaFold accession fixed (drop the UniParc version suffix) "
         "and the model cached",
         "`active_site.py` bu tipler icin UniParc'in surum sonekli kimligini "
         "(`T1YXQ1.1`) AlphaFold API'sine gonderiyor; API surum sonekli kimlige "
         "HTTP 400 veriyor ve tip 'modeli yok' diye kaydediliyor. Sonek atilinca "
         "model donuyor ve dizi uzunlugu kuratorlu dizi ile birebir ayni. "
         "Tek satirlik bir duzeltme iki tipe yapi kazandirir."),
        ("2. Eksik SMILES", "needs a substrate_smiles",
         "Bu tiplerde substrat adi var ama makine okunur yapisi yok; arama, "
         "benzerlik ve cizim bu tipleri atliyor."),
        ("3. Izlenemeyen kaynak", "needs a citable source (DOI or PMID)",
         "`source` alani serbest metin ('established literature on the reference "
         "enzyme, no individual citation'). Okur iddianin kaynagina ulasamaz. "
         "Her tip icin tek bir DOI ya da PMID yeterli."),
        ("4. Yalnizca yapi atfi var, makale yok",
         "needs a paper citation (only a PDB entry is cited)",
         "Yapinin varligi substrat atamasini kanitlamaz; bu tiplerde substrati "
         "bildiren makale atfi eksik."),
        ("5. UniProt kimligi hic yok",
         "needs a UniProt accession before any AlphaFold model can be reached",
         "Kuratorlu dizinin CRC64 saglama toplami UniParc'ta ya hic kayit "
         "bulmuyor ya da bulunan kayit UniProtKB kimligi tasimiyor. Kimlik "
         "olmadan AlphaFold'a gidilemez; cozum kuratorlu diziyi UniProt'ta "
         "kayitli bir girisle degistirmek ya da kimligi elle eklemek."),
        ("6. AlphaFold modeli gercekten yok",
         "no AlphaFold model exists for the mapped accession -- "
         "needs an experimental structure or another accession",
         "Kimlik dogru ve surum soneki de denendi; AlphaFold DB'de kayit yok. "
         "Bu tipler yapi ekseninde ancak deneysel bir yapi ya da baska bir "
         "kimlikle kapanir."),
        ("7. Dogrulanmis uyesi olmayan tip", "has no confirmed members in the database",
         "Tip kuratorlu sette duruyor ama HMM taramasi ona hicbir giris "
         "atamadi. Bunlar cogu kez ayni enzimin baska bir adla ikinci kez "
         "girilmesinden kaynaklanir; birlestirilmeli ya da neden bos oldugu "
         "notlanmali."),
        ("8. source_kind metinle celisiyor", "needs source_kind reconciled with the source text",
         "`source_kind` bir makale ya da yapi iddia ediyor ama `source` metninde "
         "o kanit yok (ya da tersi)."),
    ]
    actions = s["actions_by_what_would_fix_it"]
    for title, key, why in priority:
        items = actions.get(key, [])
        lines.append("### %s (%d tip)" % (title, len(items)))
        lines.append("")
        lines.append(why)
        lines.append("")
        lines.append(listing(items))
        lines.append("")
        if key.startswith("needs the AlphaFold accession fixed"):
            for t in items:
                af = by_type[t]["structure"]["alphafold"]
                lines.append("- `%s`: %s -> **%s**, v%s, mean pLDDT %s, %s"
                             % (t, ", ".join(af["accessions_tried_by_active_site"]),
                                af["reference_accession"], af["alphafold_version"],
                                af["global_mean_plddt"], af["model_url"]))
            lines.append("")

    handled = {k for _, k, _ in priority}
    rest = {k: v for k, v in actions.items() if k not in handled}
    if rest:
        lines.append("### Diger eksikler")
        lines.append("")
        for k, v in sorted(rest.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            lines.append("- **%s** (%d): %s" % (k, len(v), listing(v)))
        lines.append("")

    lines.append("## Yapi ekseni -- tip tip")
    lines.append("")
    lines.append("Deneysel yapisi yerelde olan %d tip %d ayri PDB kaydina dusuyor "
                 "(ayni enzim sette birden fazla adla duruyor):"
                 % (len(s["structure"]["types_with_experimental_structure"]),
                    s["structure"]["distinct_pdb_entries"]))
    lines.append("")
    for t in s["structure"]["types_with_experimental_structure"]:
        st = by_type[t]["structure"]
        lines.append("- `%s` -- %s (`%s`)" % (t, st["pdb_field"]["value"],
                                              st["experimental_file_local"]))
    lines.append("")
    lines.append("Hic yapisi olmayan %d tip: %s"
                 % (len(s["structure"]["types_with_no_structure_at_all"]),
                    listing(s["structure"]["types_with_no_structure_at_all"])))
    lines.append("")
    lines.append("## Tip basina degil, SEMA duzeyinde eksikler")
    lines.append("")
    lines.append("- %s" % s["chemistry"]["schema_level_gap"])
    lines.append("- %s" % s["structure"]["note"])
    lines.append("- SMILES denetimi: %s" % s["chemistry"]["smiles_validation_caveat"])
    low = [r["type"] for r in rows
           if (r["literature"].get("curation_confidence") or "") == "low"]
    lines.append("- `curation_confidence` %d tipte `low`, kalaninda `high`: %s. "
                 "Yani guven alani pratikte ayrim yapmiyor; ara bir `medium` "
                 "kademesi hic kullanilmamis."
                 % (len(low), listing(low)))
    lines.append("- Yalnizca ongorulen (AlphaFold) modeli olan %d tip bir EKSIK "
                 "degil, bir SINIR: ongorulen modelde metal yok, bu yuzden aktif "
                 "bolge demire gore tanimlanamaz (bkz. `active_site.py`)."
                 % len(s["structure"]["types_with_only_a_predicted_model"]))
    lines.append("")
    lines.append("## Olculmeyen seyler")
    lines.append("")
    for lim in out["limitations"]:
        lines.append("- %s" % lim)
    lines.append("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def print_summary(out):
    s = out["summary"]
    rows = out["types"]
    print()
    print("=== ROAR-DB completeness audit ===")
    print("%-58s %s" % ("denetlenen tip", s["types_audited"]))
    print("%-58s %s" % ("uc eksende de eksigi olmayan", s["types_complete_on_all_three_axes"]))
    print()
    print("--- yapi ---")
    for label, key in [("deneysel yapi (yerelde)", "types_with_experimental_structure"),
                       ("yalnizca ongorulen model (yerelde)", "types_with_only_a_predicted_model"),
                       ("hic yapi yok", "types_with_no_structure_at_all")]:
        print("%-58s %s" % (label, len(s["structure"][key])))
    print("%-58s %s" % ("chemistry.csv'de pdb alani dolu", s["structure"]["types_with_a_pdb_field"]))
    print("%-58s %s" % ("bozuk pdb kimligi", len(s["structure"]["malformed_pdb_identifiers"])))
    print("%-58s %s" % ("pdb bildirilmis ama dosya yok",
                        len(s["structure"]["pdb_declared_but_file_absent"])))
    print("%-58s %s" % ("sqlite'da yapi tablosu",
                        s["structure"]["structure_tables_found_in_sqlite"] or "yok"))
    print()
    print("--- kimya ---")
    print("%-58s %s" % ("kimya satiri olan tip", s["chemistry"]["types_with_a_chemistry_row"]))
    print("%-58s %s" % ("butun zorunlu alanlari dolu",
                        s["chemistry"]["types_with_every_required_field"]))
    for field, items in s["chemistry"]["per_field_gaps"].items():
        print("%-58s %s" % ("  eksik: " + field, ", ".join(items) or "-"))
    print("%-58s %s" % ("SMILES denetim yontemi", s["chemistry"]["smiles_validation_method"]))
    print("%-58s %s" % ("SMILES denetimi basarisiz",
                        ", ".join(s["chemistry"]["smiles_failing_the_check"]) or "-"))
    print("%-58s %s" % ("uyesi olmayan tip",
                        len(s["chemistry"]["types_in_chemistry_csv_without_confirmed_members"])))
    print("%-58s %s" % ("uyesi var ama kimya satiri yok",
                        ", ".join(s["chemistry"]["types_with_members_but_no_chemistry_row"])
                        or "yok (dogrulandi)"))
    print()
    print("--- literatur ---")
    lit = s["literature"]
    for label, key in [("izlenebilir kimlik (DOI/PMID/PMCID)",
                        "types_with_a_resolvable_identifier"),
                       ("dergi atifi, DOI/PMID yok",
                        "types_with_a_journal_citation_but_no_doi_or_pmid"),
                       ("yalnizca yapi kaniti", "types_with_only_a_structure_as_evidence"),
                       ("ne makale ne yapi", "types_with_neither_a_paper_nor_a_structure")]:
        print("%-58s %s" % (label, len(lit[key])))
    print("%-58s %s" % ("source_kind celiskisi", len(lit["types_with_source_kind_conflicts"])))
    print()
    print("--- ne yapilirsa duzelir (tip sayisi) ---")
    for action, items in s["actions_by_what_would_fix_it"].items():
        print("%4d  %s" % (len(items), action))
    print()
    incomplete = [r["type"] for r in rows if not r["complete"]]
    print("eksigi olan tip: %d / %d" % (len(incomplete), len(rows)))


if __name__ == "__main__":
    main()
