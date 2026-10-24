"""citations.yaml icindeki her DOI'yi Crossref'e karsi dogrula.

NEDEN. Atiflar elle yazilinca yanlis cilt/sayfa/yil kolayca girer ve bir
veritabaninda yanlis atif, atifsizdan kotudur. Bu betik her DOI'yi cozer ve
kutukteki kunyeyle karsilastirir. Uyusmazligi RAPOR EDER, kendi basina
DUZELTMEZ -- duzeltmeyi kurator yapar.

Cikis kodu 1: en az bir DOI cozulemedi ya da baslik uyusmadi.
"""
import argparse, json, os, re, sys, time, urllib.parse, urllib.request
import yaml

UA = "ROAR-DB/1.0 (mailto:recocanextra@gmail.com)"


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    return json.load(urllib.request.urlopen(req, timeout=60))["message"]


def norm(s):
    """Karsilastirma icin sadelestir.

    Crossref basliklarda HTML etiketi dondurebiliyor: Pruzinska 2003'un
    basligi "Pheophorbide <i>a</i> oxygenase" olarak geliyor ve etiket
    cikarilmazsa "ia" diye okunup baslik UYUSMAZ sayiliyordu -- yanlis alarm.
    """
    s = re.sub(r"<[^>]+>", "", s or "")
    return "".join(ch for ch in s.lower() if ch.isalnum())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", default="citations.yaml")
    ap.add_argument("--types-dir", default="types",
                    help="manifestolardaki citations anahtarlari da dogrulanir")
    a = ap.parse_args()
    cites = yaml.safe_load(open(a.file))["citations"]

    ids = [c["id"] for c in cites]
    bad = 0

    # Manifestolardaki atif anahtarlari kutukte var mi?
    #
    # NEDEN. types/*.yaml icindeki `citations` bir ANAHTAR listesi; kunye
    # citations.yaml'da duruyor. Anahtar yanlis yazilirsa ya da kutukten bir
    # kayit silinirse, sayfa atifsiz kalir ve bunu hicbir sey soylemez.
    import glob                                       # noqa: PLC0415
    known = set(ids)
    missing = []
    for path in sorted(glob.glob(os.path.join(a.types_dir, "*.yaml"))):
        try:
            man = yaml.safe_load(open(path)) or {}
        except Exception:                             # noqa: BLE001
            continue
        for key in (man.get("citations") or []):
            if key not in known:
                missing.append(f"{os.path.basename(path)}: {key}")
    if missing:
        print(f"[HATA] kutukte olmayan atif anahtari ({len(missing)}):")
        for m in missing[:10]:
            print("   ", m)
        bad += 1
    else:
        n = sum(1 for p2 in glob.glob(os.path.join(a.types_dir, "*.yaml"))
                if (yaml.safe_load(open(p2)) or {}).get("citations"))
        print(f"[OK   ] manifesto atif anahtarlari: {n} tipte, hepsi kutukte")
    if len(ids) != len(set(ids)):
        print("[HATA] id cakismasi var")
        bad += 1

    for c in cites:
        for field in ("id", "authors", "title", "year", "doi", "verified_by", "used_for"):
            if field not in c:
                print(f"[HATA] {c.get('id','?')}: '{field}' alani yok")
                bad += 1
    print()

    for c in cites:
        doi = c["doi"]
        try:
            m = crossref(doi)
        except Exception as exc:                      # noqa: BLE001
            print(f"[COZULEMEDI] {c['id']}  {doi}  ({type(exc).__name__})")
            bad += 1
            time.sleep(0.5)
            continue
        title = (m.get("title") or [""])[0]
        if norm(title)[:40] not in norm(c["title"]) and norm(c["title"])[:40] not in norm(title):
            print(f"[BASLIK UYUSMADI] {c['id']}  {doi}")
            print(f"    kutuk   : {c['title'][:76]}")
            print(f"    crossref: {title[:76]}")
            bad += 1
            time.sleep(0.5)
            continue
        notes = []
        yr = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
        if yr and int(c["year"]) != int(yr):
            notes.append(f"yil kutuk={c['year']} crossref={yr}")
        for key, got in (("volume", m.get("volume")), ("pages", m.get("page"))):
            if got and c.get(key) and str(c[key]) != str(got):
                notes.append(f"{key} kutuk={c[key]} crossref={got}")
            elif got and not c.get(key):
                notes.append(f"{key} kutukte YOK, crossref={got}")
        flag = "UYARI" if notes else "OK   "
        print(f"[{flag}] {c['id']:26s} {doi}")
        for n in notes:
            print(f"         {n}")
        time.sleep(0.5)

    print()
    if bad:
        print(f"{bad} sorun var.")
        return 1
    print(f"{len(cites)} atifin hepsi cozuldu ve basliklari uyusuyor.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
