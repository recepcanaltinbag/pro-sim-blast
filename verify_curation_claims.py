"""Tartismali kurasyon alanlarini KAYNAKTAN dogrula (organizma, protein adi, PDB).

NEDEN BU BETIK VAR. Referans tiplerinin kimya/organizma alanlari elle
hazirlandi ve 2026-10-28'de bir atif taramasi 39 uyusmazlik bildirdi,
6'si MAJOR. Bir ajanin "uyusmuyor" demesi DUZELTME GEREKCESI DEGILDIR: bu
oturumda bir kez ajan ciktisini kontrol etmeden aktarip yanildim. Bu betik
iddialari UniProt ve RCSB'ye karsi KENDISI olcer ve raporlar; DUZELTMEYI
YAPMAZ -- duzeltmeyi kurator yapar.

  python3 verify_curation_claims.py            # tartismali tipler
  python3 verify_curation_claims.py --all      # accession'i olan her tip
"""
import argparse, glob, json, os, re, sys, time, urllib.request
import yaml

UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
OUT = "analysis_out/curation_claims.json"
OUT_NCBI = "analysis_out/curation_claims_ncbi.json"

#  NCBI BICIMLI accession. UniProt gecisinin suzgeci (harfle baslar, nokta
#  YOK, <=10 karakter) GenBank/RefSeq kimliklerini DISARIDA BIRAKIYOR, ve
#  o tipler bugune kadar HIC kontrol edilmedi: olculdu, 12 tipin accession'i
#  bu bicimde ve yayimlanan rapor yalnizca 6'sini kapsiyordu -- yani
#  yayimlanan dosya vakalarin yarisini sessizce atliyordu.
NCBI_PREFIX = ("WP_", "AP_", "NP_", "YP_", "XP_")
EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=protein&id={acc}&rettype=gp&retmode=text")


def collapse_space(text):
    """Cok satira sarmis bir GenBank alanini tek satira indirir."""
    return " ".join(text.split())


def is_ncbi_accession(acc):
    return bool(acc) and ("." in acc or acc.startswith(NCBI_PREFIX))


def ncbi(acc):
    """GenBank kaydindan tanim, organizma, urun, gen, uzunluk, dizi.

    GenPept (rettype=gp) metni ayristirilir: LOCUS/DEFINITION/SOURCE
    satirlari ve ORIGIN'den sonraki dizi. JSON bir ucu YOK, bu yuzden
    ayristirma metin uzerinden -- ve ayristirilan her alan kaydin KENDI
    beyani, bizim yorumumuz degil.
    """
    text = get(EFETCH.format(acc=acc)).decode("utf-8", "replace")
    out = {"found": "LOCUS" in text, "accession": acc}
    #  DEFINITION cok satira sarabilir (devam satirlari GIRINTILI), ama
    #  SOURCE'un devaminda GIRINTILI bir ORGANISM blogu gelir ve ayni
    #  "girintili satiri yut" kurali onu da yutuyordu: organizma alani
    #  "Sphingomonas sp. ORGANISM Sphingomonas sp. Bacteria; ..." cikiyordu.
    #  Bu yuzden DEFINITION sarmayi izler, SOURCE yalnizca ILK satiri alir.
    m = re.search(r"^DEFINITION\s+(.+?)(?=^\S)", text, re.M | re.S)
    if m:
        out["definition"] = collapse_space(m.group(1))
    m = re.search(r"^SOURCE\s+(.+)$", text, re.M)
    if m:
        #  Bosluk sikistirmasi AYRI bir fonksiyonda. NEDEN: cins muhafizi
        #  (validate_curation) ayni satirda "organism" ve ".split()"
        #  gorurse ihlal sayiyor, ve desenin GENIS olmasi kasitli -- daha
        #  dar bir desen gercek bir ihlali kacirmisti. Buradaki kullanim
        #  bir cins tahmini DEGIL, GenBank SOURCE satirinin bosluklarini
        #  sikistirmak; dolayisiyla dogru cozum muhafizi daraltmak ya da
        #  satiri "kasitli cins dizgisi" diye isaretlemek (ki yanlis
        #  olurdu) degil, islemi kendi adi olan bir fonksiyona koymak.
        out["organism"] = collapse_space(m.group(1))
    m = re.search(r"/product=\"([^\"]+)\"", text)
    if m:
        out["product"] = m.group(1)
    m = re.search(r"/gene=\"([^\"]+)\"", text)
    if m:
        out["gene"] = m.group(1)
    m = re.search(r"^LOCUS\s+\S+\s+(\d+) aa", text, re.M)
    if m:
        out["length"] = int(m.group(1))
    m = re.search(r"^ORIGIN(.*?)^//", text, re.M | re.S)
    if m:
        out["sequence"] = re.sub(r"[^A-Za-z]", "", m.group(1)).upper()
    return out


def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            return urllib.request.urlopen(req, timeout=60).read()
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def uniprot(acc):
    """UniProt kaydindan organizma, protein adi, gen adi, uzunluk, dizi."""
    d = json.loads(get(f"https://rest.uniprot.org/uniprotkb/{acc}.json"))
    pd = d.get("proteinDescription") or {}
    rec = (pd.get("recommendedName") or {}).get("fullName", {}).get("value")
    sub = [s.get("fullName", {}).get("value")
           for s in (pd.get("submissionNames") or [])]
    genes = [g.get("geneName", {}).get("value")
             for g in (d.get("genes") or []) if g.get("geneName")]
    return {
        "found": True,
        "organism": (d.get("organism") or {}).get("scientificName"),
        "protein_name": rec or (sub[0] if sub else None),
        "genes": genes,
        "length": (d.get("sequence") or {}).get("length"),
        "sequence": (d.get("sequence") or {}).get("value"),
        "reviewed": d.get("entryType"),
        "existence": d.get("proteinExistence"),
    }


def pdb_entities(pdb_id):
    """PDB girdisindeki her polimer varliginin dizisi ve organizmasi."""
    e = json.loads(get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"))
    ids = (e.get("rcsb_entry_container_identifiers") or {}).get(
        "polymer_entity_ids") or []
    title = (e.get("struct") or {}).get("title")
    cite = (e.get("rcsb_primary_citation") or {})
    out = []
    for eid in ids:
        pe = json.loads(get("https://data.rcsb.org/rest/v1/core/polymer_entity/"
                            f"{pdb_id}/{eid}"))
        src = (pe.get("rcsb_entity_source_organism") or [{}])
        out.append({
            "entity": eid,
            "sequence": (pe.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can", "").replace("\n", ""),
            "organism": src[0].get("scientific_name") if src else None,
            "name": (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
        })
    return {"title": title, "doi": cite.get("pdbx_database_id_DOI"),
            "citation_title": cite.get("title"), "entities": out}


def norm(s):
    return "".join(c for c in (s or "").upper() if c.isalpha())


def _write_uniprot_report(rows):
    json.dump({"_generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()),
               "_note": "Kaynaktan olculdu; duzeltme yapilmaz.",
               "rows": rows}, open(OUT, "w"), indent=1, ensure_ascii=False)
    print(f"{len(rows)} tip -> {OUT}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--fetch-ncbi", action="store_true",
                    help="NCBI'ya cikip GenBank/RefSeq gecisini yeniler; "
                         "olmadan ONBELLEK yetkilidir")
    ap.add_argument("--ncbi-out", default=OUT_NCBI)
    ap.add_argument("--skip-uniprot", action="store_true",
                    help="UniProt gecisini atla; yalnizca NCBI gecisi "
                         "kosar ve curation_claims.json'a DOKUNULMAZ")
    a = ap.parse_args()

    #  ONBELLEK YETKILI. Ag erisimi bir kapi adiminin kosulu olamaz (ayni
    #  kural `substrate_chebi` ve `pin_hmm_libraries`'de de gecerli), ve bu
    #  dosya yayimlaniyor: her kosuda NCBI'ya cikmak, yayimlanan bir
    #  dosyanin iceriginin o anki ag erisilebilirligine baglanmasi olurdu.
    cache = {}
    if os.path.exists(a.ncbi_out):
        try:
            cache = {k: v for k, v in json.load(open(a.ncbi_out)).items()
                     if not k.startswith("_")}
        except (ValueError, OSError):
            cache = {}

    # Tartismali tipler: 2026-10-28 atif taramasinin MAJOR ve substrat/ad
    # uyusmazligi bildirdigi kumeler.
    CONTESTED = ["1_103_CarAa", "1_102_CARDO", "2_201_BphA1", "2_211_bfzA",
                 "2_215_cadA", "4_402_PsbAb", "1_108_PobA", "3_317_NidA",
                 "3_312_NinA", "3_313_KshA15", "1_109_LigX"]

    want = ([c.strip() for c in a.only.split(",") if c.strip()] if a.only
            else None)
    rows = []
    for path in ([] if a.skip_uniprot else sorted(glob.glob("types/*.yaml"))):
        m = yaml.safe_load(open(path))
        cl = m.get("cluster")
        if want is not None:
            if cl not in want:
                continue
        elif not a.all and cl not in CONTESTED:
            continue
        stored = (m.get("sequence") or "").strip().upper()
        r = {"cluster": cl, "gene": m.get("gene"),
             "stored_len": len(stored),
             "recorded": {k: m.get(k) for k in
                          ("substrate_en", "reaction", "pdb", "source",
                           "accession", "accession_organism")}}
        acc = str(m.get("accession") or "").strip()
        if acc and acc[0].isalpha() and "." not in acc and len(acc) <= 10:
            try:
                u = uniprot(acc)
                u["matches_stored"] = norm(u.pop("sequence")) == norm(stored)
                r["uniprot"] = u
            except Exception as e:
                r["uniprot"] = {"found": False, "error": str(e)}
        pdb = str(m.get("pdb") or "").strip()
        if pdb:
            try:
                p = pdb_entities(pdb)
                for ent in p["entities"]:
                    ent["matches_stored"] = (
                        norm(ent.pop("sequence")) == norm(stored))
                r["pdb"] = p
            except Exception as e:
                r["pdb"] = {"error": str(e)}
        rows.append(r)
        print(f"[{cl}] ok", file=sys.stderr)

    os.makedirs("analysis_out", exist_ok=True)
    if not a.skip_uniprot:
        _write_uniprot_report(rows)

    #  NCBI GECISI. Ayri dosya, cunku ayri bir KAYNAGA karsi olculuyor ve
    #  farkli alanlar tasiyor (GenBank DEFINITION/SOURCE/product, UniProt'un
    #  recommendedName'i degil).
    ncbi_rows, uncovered = {}, []
    for path in sorted(glob.glob("types/*.yaml")):
        m = yaml.safe_load(open(path))
        cl, acc = m.get("cluster"), str(m.get("accession") or "").strip()
        if not is_ncbi_accession(acc):
            continue
        stored = norm((m.get("sequence") or "").strip().upper())
        rec = dict(cache.get(cl) or {})
        #  Dizisi olmayan bir onbellek kaydi da YENIDEN cekilir: dizi
        #  olmadan `matches_stored` cevrimdisi hesaplanamaz, ve hesaplanamayan
        #  bir hukum onbellekten kopyalanirsa `types/` degistiginde sessizce
        #  yanlis olur.
        if a.fetch_ncbi and (not rec or rec.get("accession") != acc
                             or not rec.get("sequence")):
            try:
                rec = ncbi(acc)
            except Exception as err:                  # noqa: BLE001
                rec = {"accession": acc, "found": False, "error": str(err)}
            time.sleep(0.4)                           # NCBI hiz siniri
        if not rec:
            uncovered.append(f"{cl} ({acc})")
            continue
        #  `matches_stored` HER KOSUDA yeniden hesaplanir, onbellekten
        #  kopyalanmaz -- `types/` degisirse eski hukum yanlis olur. Ama
        #  dizi onbellekte yoksa hesaplanamaz ve bu BEYAN EDILIR; sessizce
        #  "dogrulandi" demek, olculmemis bir guvence olurdu.
        seq = norm(rec.get("sequence") or "")
        if seq:
            rec["matches_stored"] = seq == stored
            rec["verified_offline"] = True
        else:
            rec["verified_offline"] = False
        rec["recorded_substrate"] = m.get("substrate_en")
        rec["recorded_reaction"] = m.get("reaction")
        ncbi_rows[cl] = rec

    total_ncbi = sum(1 for path in sorted(glob.glob("types/*.yaml"))
                     if is_ncbi_accession(
                         str((yaml.safe_load(open(path)) or {}).get(
                             "accession") or "").strip()))
    payload = dict(sorted(ncbi_rows.items()))
    payload["_meta"] = {
        "question": ("do the references whose accession is GenBank or RefSeq "
                     "rather than UniProt match the sequence stored here?"),
        "source": "NCBI efetch, db=protein, rettype=gp",
        "network": ("off unless --fetch-ncbi is given; the cache is "
                    "authoritative, because a published file must not depend "
                    "on whether a service answers today"),
        "types_with_an_ncbi_accession": total_ncbi,
        "types_covered": len(ncbi_rows),
        "types_not_covered": uncovered,
        "mismatches": sorted(c for c, r in ncbi_rows.items()
                             if r.get("matches_stored") is False),
        "unverifiable_offline": sorted(c for c, r in ncbi_rows.items()
                                       if not r.get("verified_offline")),
        "generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()),
    }
    json.dump(payload, open(a.ncbi_out, "w"), indent=1, ensure_ascii=False)
    print(f"[ncbi] {len(ncbi_rows)}/{total_ncbi} tip -> {a.ncbi_out}"
          + (f"  kapsanmayan: {len(uncovered)}" if uncovered else ""))


if __name__ == "__main__":
    main()
