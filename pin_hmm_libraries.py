"""Her HMM kutuphanesini SABITLER: hangi modelin hangi SURUMU kullanildi.

SORU. Yayimlanmis bir HMM isabeti tekrar uretilebilir mi? Bir domain
cagrisi ("bu komsu bir ferredoksin") modelin kendisine baglidir, ve Pfam
modelleri surum atlar: `component_hmms/components.hmm` icinde PF00355.33
duruyor, InterPro bugun PF00355.34 veriyor. Iki surum ayni esikte ayni
diziye ayni skoru VERMEK ZORUNDA DEGIL. Dolayisiyla "Pfam modeli
kullanildi" cumlesi bir yontem tarifi degil.

NE OLCULUYOR. Diskteki her kutuphanenin icinde ne oldugu: model basina
NAME / ACC (surumle birlikte) / LENG, model blogunun ozeti, ve dosyanin
tamaminin ozeti. Olculen sey KULLANILAN kutuphanedir, indirilmesi gereken
kutuphane degil.

NE OLCULMUYOR -- ve bu kasitli. Pfam SURUM ETIKETI (37.0 gibi) model
dosyasinin icinde YOK; HMM kaydi yalnizca accession'in kendi surumunu
tasir. Bu yuzden burada bir surum etiketi UYDURULMUYOR: alinti birimi
`PF00355.33` gibi accession.surum ciftidir, ki o etiketten DAHA kesindir
(bir release icinde de modeller ayri ayri surum atar).

SURUKLENME. `--check-upstream` ile InterPro'ya her accession'in BUGUN
hangi surumde oldugu sorulur ve fark kaydedilir. Bu bir hata degil, bir
OLGU: sabitlenmis kutuphane bilerek geride olabilir. Hata, farkin
BILINMEMESIDIR. Ag erisimi gerektirdigi icin kapi bu kismi istemez;
sonuc onbellege yazilir ve raporda tarihiyle durur.

PROJENIN KENDI MODELLERI (ROmotif71, RieskeDB71) Pfam'dan gelmez; onlar
yalnizca OZET ile sabitlenir ve rapor bunu ayrica soyler -- bir accession
surumu olmayan modeli "surumlendi" diye gostermek, olculmemis bir
guvence olurdu.

Cikti: analysis_out/hmm_pins.json
"""

import argparse
import gzip
import hashlib
import json
import os
import re
import time
import urllib.request

#  Sabitlenen kutuphaneler. Her biri AYRI bir soruya cevap veriyor ve
#  bu yuzden ayri dosyada duruyor (bkz. fetch_family_hmms.py basligi).
LIBRARIES = [
    ("component_hmms/components.hmm",
     "operon components (beta, ferredoxin, reductase) -- input of operon derivation"),
    ("family_hmms/mobility.hmm", "transposase and IS families of neighbouring genes"),
    ("family_hmms/recombinase.hmm", "recombinase and integrase families"),
    ("family_hmms/regulator.hmm", "transcriptional regulator families"),
    ("ROs_71_Clean/ROmotif71.hmm",
     "project-built motif model of the Rieske alpha subunit"),
    ("RieskeDB71.hmm",
     "project-built library of the 76 curated reference types, used for type "
     "assignment by bit score"),
]

PFAM_ACC = re.compile(r"^PF\d{5}\.\d+$")
UPSTREAM_URL = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}/?annotation=hmm"
CACHE = "analysis_out/hmm_upstream_cache.json"


def read_models(path):
    """Bir HMMER3 kutuphanesindeki her modelin kimligi + blok ozeti.

    Kendi ayristirici: HMM kaydi `HMMER3/f` ile baslar ve `//` ile biter.
    Bir model blogunun ozeti, o modelin TAM metni uzerinden alinir, yani
    bir esik ya da bir olasilik degisse ozet degisir.
    """
    opener = gzip.open if path.endswith(".gz") else open
    models, block = [], []
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            block.append(line)
            if line.startswith("//"):
                text = "".join(block)
                got = {"digest": hashlib.sha256(text.encode()).hexdigest()[:16]}
                for key, field in (("NAME", "name"), ("ACC", "accession"),
                                   ("DESC", "description"), ("LENG", "length")):
                    match = re.search(rf"^{key}\s+(.+)$", text, re.MULTILINE)
                    if match:
                        got[field] = match.group(1).strip()
                if "length" in got:
                    got["length"] = int(got["length"])
                models.append(got)
                block = []
    return models


def file_digest(path):
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            sha.update(chunk)
    return sha.hexdigest()[:16]


def upstream_version(acc, cache):
    """InterPro'nun BUGUN verdigi accession.surum. Onbellekten okur."""
    base = acc.split(".")[0]
    if base in cache:
        return cache[base]
    try:
        with urllib.request.urlopen(UPSTREAM_URL.format(acc=base), timeout=60) as res:
            raw = res.read()
        try:
            text = gzip.decompress(raw).decode("utf-8", "replace")
        except OSError:
            text = raw.decode("utf-8", "replace")
        match = re.search(r"^ACC\s+(PF\d{5}\.\d+)$", text, re.MULTILINE)
        got = match.group(1) if match else None
    except Exception as exc:                      # ag hatasi bir olcum degil
        got = f"error: {type(exc).__name__}"
    cache[base] = got
    return got


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="analysis_out/hmm_pins.json")
    ap.add_argument("--check-upstream", action="store_true",
                    help="InterPro'ya sorup suruklenmeyi kaydet (ag gerekir)")
    args = ap.parse_args()

    libraries, pfam_pins, project_pins = [], [], []
    for path, purpose in LIBRARIES:
        if not os.path.exists(path):
            libraries.append({"library": path, "purpose": purpose,
                              "present": False})
            continue
        models = read_models(path)
        versioned = [m for m in models if PFAM_ACC.match(m.get("accession", ""))]
        libraries.append({
            "library": path,
            "purpose": purpose,
            "present": True,
            "file_digest": file_digest(path),
            "models": len(models),
            "pfam_versioned_models": len(versioned),
            "pinned_by": ("accession version and digest" if versioned
                          else "digest only (not a Pfam library)"),
        })
        for model in models:
            row = dict(model, library=path)
            (pfam_pins if PFAM_ACC.match(model.get("accession", ""))
             else project_pins).append(row)

    drift = {"checked": False}
    if args.check_upstream:
        cache = {}
        if os.path.exists(CACHE):
            with open(CACHE) as handle:
                cache = json.load(handle).get("versions", {})
        rows = []
        for pin in pfam_pins:
            now = upstream_version(pin["accession"], cache)
            if now and not str(now).startswith("error:") and now != pin["accession"]:
                rows.append({"library": pin["library"], "name": pin["name"],
                             "pinned": pin["accession"], "upstream_now": now})
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        with open(CACHE, "w") as handle:
            json.dump({"checked_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                    time.gmtime()),
                       "versions": cache}, handle, indent=1)
        drift = {
            "checked": True,
            "checked_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "pfam_models_checked": len(pfam_pins),
            "models_behind_upstream": len(rows),
            "detail": rows,
            "reading": "A pinned library that is behind upstream is not an "
                       "error: re-downloading would change every published HMM "
                       "call. What would be an error is not knowing the "
                       "difference, which is why it is recorded here.",
        }

    report = {
        "question": "Which exact model produced each published HMM call, and "
                    "can that call be reproduced?",
        "method": {
            "unit": "one HMM model inside one library file",
            "citable_pin": "the accession with its own version (e.g. PF00355.33) "
                           "plus a digest of the model block; models version "
                           "independently, so this is finer than a Pfam release "
                           "label",
            "no_release_label": "The Pfam release (e.g. 37.0) is not written "
                                "inside an HMM record, so no release label is "
                                "invented here.",
            "project_models": "ROmotif71 and RieskeDB71 are built by this "
                              "project, not fetched from Pfam. They carry no "
                              "accession version and are pinned by digest only; "
                              "they are listed separately rather than being "
                              "shown as if they were versioned.",
        },
        "libraries": libraries,
        "pfam_models": sorted(pfam_pins, key=lambda r: r["accession"]),
        "project_models": sorted(project_pins, key=lambda r: r.get("name", "")),
        "upstream_drift": drift,
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=1)

    for lib in libraries:
        if not lib["present"]:
            print(f"[yok]  {lib['library']}")
            continue
        print(f"[sabit] {lib['library']:34s} {lib['models']:3d} model, "
              f"{lib['pfam_versioned_models']:3d} Pfam surumlu, "
              f"ozet {lib['file_digest']}")
    if drift["checked"]:
        print(f"[yukari akis] {drift['pfam_models_checked']} Pfam modeli "
              f"soruldu, {drift['models_behind_upstream']} tanesi geride")
        for row in drift["detail"]:
            print(f"    {row['name']:18s} {row['pinned']} -> {row['upstream_now']}")
    print(f"[yazildi] {args.out}")


if __name__ == "__main__":
    main()
