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
import argparse, glob, json, os, sys, time, urllib.request
import yaml

UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"
OUT = "analysis_out/curation_claims.json"


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--only", default="")
    a = ap.parse_args()

    # Tartismali tipler: 2026-10-28 atif taramasinin MAJOR ve substrat/ad
    # uyusmazligi bildirdigi kumeler.
    CONTESTED = ["1_103_CarAa", "1_102_CARDO", "2_201_BphA1", "2_211_bfzA",
                 "2_215_cadA", "4_402_PsbAb", "1_108_PobA", "3_317_NidA",
                 "3_312_NinA", "3_313_KshA15", "1_109_LigX"]

    want = ([c.strip() for c in a.only.split(",") if c.strip()] if a.only
            else None)
    rows = []
    for path in sorted(glob.glob("types/*.yaml")):
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
    json.dump({"_generated": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()),
               "_note": "Kaynaktan olculdu; duzeltme yapilmaz.",
               "rows": rows}, open(OUT, "w"), indent=1, ensure_ascii=False)
    print(f"{len(rows)} tip -> {OUT}")


if __name__ == "__main__":
    main()
