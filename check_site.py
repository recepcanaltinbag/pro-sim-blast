"""
Yayina cikmadan once sitenin kendisini denetler: her sayfa aciliyor mu, her
dahili baglanti 200 donuyor mu, arayuzde Turkce metin kalmis mi.

NEDEN BU SCRIPT VAR. Arayuzde Turkce etiketler bulundu ("bitki/alg"), ve ilk
elle yapilan tarama bunlari KACIRDI, cunku sadece islev sozcuklerine (ve, icin,
bir) bakiyordu; kacan sey bir ICERIK sozcuguydu ve veritabanina bir Python
sozlugunden geliyordu. Ayni sekilde uye sayisi sifir olan tipler anasayfada
listelenirken sayfalari 404 donuyordu. Iki hata da ayni sebepten gorunmez
kaldi: kimse sayfalari toptan, otomatik olarak denemiyordu. Bu script onu
yapar ve FAIL varsa sifirdan farkli cikar, boylece deploy oncesi durur.

Kullanim:
    python3 check_site.py              # hepsi
    python3 check_site.py --quick      # tip sayfalarini orneklem al
"""

import argparse
import glob
import json
import os
import re
import sqlite3
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "webapp"))

# Turkce arama listesi. Hem islev hem ICERIK sozcukleri; birincisi ilk taramada
# yeterli sanilmisti ve yetmedi.
TURKISH_WORDS = """
ve veya icin ile olan bir bu su daha cok az yok var gore sayisi degil kayit
kayitlar ama cunku iliskili tablo sonuc hepsi tum bitki bitkiler algler mantar
mantarlar hayvan hayvanlar diger digerleri kume kumeler tur turler yuksek
dusuk orta buyuk kucuk toplam oran orani sayi adet yeni eski okaryot
okaryotlar bakteri bakteriyel yikim enzim enzimler dizisi dizi diziler kaynak
kaynagi bilinmeyen bilinmiyor belirsiz supheli guven guvenli dogrulanan
dogrulanmis adaylar aday bulundu bulunan hata hatali eksik tamam secili secim
ornek ornekler ortalama topraktan toprak deniz tatli tuzlu sulak hava yuzey
derin kaya kirmizi yesil mavi sari siyah beyaz baskin ayri ortak ozel genel
yerel sadece ayrica ancak yani ise hem yazildi bilgi uyari sayfasi sayfa
baslik aciklama ozet detay ayrinti nerede nasil neden hangi kac kadar sonra
simdi zaman yil gun olarak olmasi olmak yapilan yapan eden edilen veren
verilen alinan icinde uzerinde altinda yaninda arasinda disinda kisi insan
klinik hasta hastane ciger beyin doku analizi hesap hesaplanan olcum olculen
sonuclar dagilim dagilimi grafik tablosu sekil resim harita kloroplast yaprak
kok govde tohum meyve cicek
""".split()

# Ingilizce ile ESYAZIMLI olanlar ve Latin ikili adlandirma kisaltmalari.
# Bunlar listede kalirsa her sayfa yanlisca isaretlenir.
HOMOGRAPHS = {"var", "sp", "subsp", "cf", "str", "ssp", "protein", "proteinler",
              "test", "testi", "once", "an", "su", "kan", "analiz", "alg", "sari",
              "ice", "tam", "son", "ton", "pan", "ban", "men", "ten"}

TURKISH_RE = re.compile(
    r"\b(" + "|".join(sorted(set(TURKISH_WORDS) - HOMOGRAPHS, key=len, reverse=True)) + r")\b",
    re.I)
TAGS_RE = re.compile(r"<script.*?</script>|<style.*?</style>", re.S)
HREF_RE = re.compile(r'href="([^"#?][^"]*)"')


def visible_text(html):
    """Gorunur metin: script ve style ICERIGI atilir, sonra etiketler silinir.

    Script icerigi atilmak zorunda, cunku orada href'ler JavaScript ile
    birlestiriliyor ("' + rcsb + '") ve tarayici onlari baglanti sanmaz.
    """
    return re.sub(r"<[^>]+>", " ", TAGS_RE.sub(" ", html))


def check_search_parity(root, client):
    """Yayinlanan aramanin bir KIMYASALI bulabildigini dogrular.

    Statik arama indeksi bir donem yalnizca giris basina alanlari tasiyordu
    (organizma, urun, tip, varyant, aile) ve kimya yoktu. Sonuc: uygulamada
    "benzalkonium" 95 giris donuyordu, yayinlanan sitede 0. Sitenin butun
    duzeni enzimleri kimyasallarina gore gruplamak oldugu icin bu kucuk bir
    eksik degildi. Bu kontrol her kuratorlu substratin yayinlanan indekste
    gercekten BULUNABILIR oldugunu sinar.
    """
    index_path = os.path.join(root, "search_index.json")
    page_path = os.path.join(root, "search.html")
    if not (os.path.exists(index_path) and os.path.exists(page_path)):
        return ["static search index or page missing; run webapp/freeze.py"], 0

    with open(index_path, encoding="utf-8") as fh:
        index = json.load(fh)
    page = open(page_path, encoding="utf-8").read()
    match = re.search(r"const CHEM = (\{.*?\});", page, re.S)
    if not match:
        return ["search.html carries no chemistry map, so a chemical name cannot be found"], 0
    chem = json.loads(match.group(1))
    if "renderTypes" not in page:
        return ["search.html does not show matching enzyme types, so a substrate with no "
                "confirmed member is unreachable"], 0

    # Sayfanin kendi arama metnini AYNI sekilde kurar.
    haystack = []
    for row in index:
        extra = chem.get(row[4], ["", "", ""])
        haystack.append(" ".join(str(x) for x in (
            row[0], row[1], row[2], row[3], row[4], row[5], row[8],
            extra[0], extra[1], extra[2])).lower())

    problems, tested = [], 0
    for cluster, values in chem.items():
        substrate = (values[1] or "").strip()
        if not substrate:
            continue
        # Substratin ilk anlamli kelimesi aranabilir olmali.
        word = next((w for w in re.split(r"[^A-Za-z]+", substrate.lower())
                     if len(w) > 4), None)
        if not word:
            continue
        tested += 1
        # Sayfa IKI yolda ariyor: giris indeksi ve tip sozlugu. Uye sayisi
        # sifir olan tiplerin substrati hicbir giriste gecmez, ama tip sonuclari
        # bolumunde bulunur; kontrol ikisini de saymak zorunda, yoksa dogru
        # calisan bir sayfayi hatali bildirir.
        in_entries = any(word in hay for hay in haystack)
        in_types = any(word in (key + " " + " ".join(values)).lower()
                       for key, values in chem.items())
        if not (in_entries or in_types):
            problems.append(f"'{word}' from substrate {substrate!r} ({cluster}) "
                            f"matches neither an entry nor a type")
    return problems, tested


def check_static_site(root):
    """Yayinlanan statik siteyi DISKTE gezer: her baglanti bir dosyaya denk mi?

    Dinamik uygulamayi denemek yetmiyor. Statik ihracat baglantilara `.html`
    ekliyor ve yollari goreli hale getiriyor; o cevirinin bir yerde bozulmasi
    uygulamada hic gorunmez. Kullanicinin "acilmayan sayfa" sikayetlerinin
    bir kismi tam olarak buydu.
    """
    from urllib.parse import unquote, urljoin

    pages = sorted(glob.glob(os.path.join(root, "**", "*.html"), recursive=True))
    if not pages:
        return ["static site not built; run webapp/freeze.py first"], 0, 0

    # Taban yolu BUILD'DEN okunur, elle verilmez. GitHub proje sayfasi siteyi
    # /pro-sim-blast/ altinda sunuyor; `freeze.py --base` verilmeden
    # uretilince butun baglantilar /about.html gibi KOKE gider ve 404 olur.
    # Bu bir kere gerceklesti ve yayinlanan sitenin gezinmesini kirdi, bu
    # yuzden taban artik olculen bir sey.
    index = os.path.join(root, "index.html")
    base = ""
    if os.path.exists(index):
        m = re.search(r'href="([^"]*)/static/style\.css"',
                      open(index, encoding="utf-8", errors="replace").read())
        if m:
            base = m.group(1)

    problems, links = [], 0
    for page in pages:
        rel = os.path.relpath(page, root)
        html = open(page, encoding="utf-8", errors="replace").read()
        for href in set(HREF_RE.findall(TAGS_RE.sub(" ", html))):
            if href.startswith(("http", "mailto:", "//", "javascript:", "data:")):
                continue
            links += 1
            if href.startswith("/"):
                # Mutlak yol taban onekiyle BASLAMAK zorunda, yoksa sunucuda
                # site kokunun disina dusuyor.
                if base and not (href == base or href.startswith(base + "/")):
                    problems.append(f"{rel} -> {href} (missing base prefix {base!r})")
                    continue
                target = unquote(href[len(base):].lstrip("/").split("#")[0].split("?")[0])
            else:
                # Goreli yol, sayfanin KENDI konumuna gore cozulur.
                target = unquote(urljoin(rel.replace(os.sep, "/"), href)
                                 .split("#")[0].split("?")[0])
            if not target or target.startswith(".."):
                problems.append(f"{rel} -> {href} (escapes the site root)")
                continue
            path = os.path.join(root, target.replace("/", os.sep))
            if os.path.isdir(path):
                path = os.path.join(path, "index.html")
            if not os.path.exists(path):
                problems.append(f"{rel} -> {href}")
    if not base:
        problems.insert(0, "index.html carries no base prefix; the build was made "
                           "without --base and every absolute link will 404 when the "
                           "site is served from a subdirectory")
    return problems, len(pages), links


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quick", action="store_true",
                        help="tip sayfalarinin onda birini dene")
    parser.add_argument("--static", default=os.path.join("webapp", "site"),
                        help="yayinlanan statik sitenin klasoru")
    parser.add_argument("--skip-static", action="store_true",
                        help="statik site denetimini atla")
    args = parser.parse_args()

    from fastapi.testclient import TestClient
    import app as appmod

    client = TestClient(appmod.app)
    connection = sqlite3.connect(appmod.DB_PATH)

    types = [r[0] for r in connection.execute(
        "SELECT DISTINCT ro_cluster FROM ro "
        "WHERE ro_cluster IS NOT NULL AND ro_cluster <> 'N/A' ORDER BY 1")]
    if args.quick:
        types = types[::10]

    pages = (["/", "/clusters", "/classify", "/about", "/atlas",
              "/search?q=naphthalene"]
             + ["/atlas/" + name for name in
                ("phylogeny", "network", "taxonomy", "ecology", "operons", "regulation",
                 "evidence", "novel", "cooccurrence", "quality", "statistics")]
             + ["/cluster/" + t for t in types])

    failures, turkish = [], []
    link_status, broken = {}, defaultdict(set)

    for page in pages:
        response = client.get(page)
        if response.status_code != 200:
            failures.append((page, response.status_code))
            continue
        body = response.text
        words = sorted({m.group(1).lower() for m in TURKISH_RE.finditer(visible_text(body))})
        if words:
            turkish.append((page, words))
        for href in sorted(set(HREF_RE.findall(TAGS_RE.sub(" ", body)))):
            if href.startswith(("http", "mailto:", "//", "javascript:")):
                continue
            target = href if href.startswith("/") else "/" + href
            if target not in link_status:
                link_status[target] = client.get(target).status_code
            if link_status[target] != 200:
                broken[page].add((target, link_status[target]))

    # Yayinlanan veri dosyalarinda da Turkce aranir: indirme baglantilari
    # sayfanin bir parcasi, icerikleri de kullaniciya gidiyor.
    # Liste atlas.PROVENANCE'tan ALINIR, freeze.py kaynagindan KAZINMAZ.
    # Kazima bir donem calisiyordu cunku dosya adlari orada duz dizgi olarak
    # yaziliydi; liste turetilmis hale gelince regex hicbir sey bulamadi ve bu
    # kontrol sessizce bos calisti. Olculdu: 21 dosya yerine 2.
    import atlas as _atlas
    published = set(_atlas.ANALYSIS_FILES)
    data_turkish = []
    for name in sorted(published):
        path = os.path.join("analysis_out", name)
        if not os.path.exists(path):
            continue
        text = open(path, encoding="utf-8", errors="replace").read()
        words = sorted({m.group(1).lower() for m in TURKISH_RE.finditer(text)})
        if words:
            data_turkish.append((name, words))

    print("=" * 72)
    print("SITE DENETIMI")
    print("=" * 72)
    print(f"  pages requested          {len(pages)}")
    print(f"  pages served             {len(pages) - len(failures)}")
    print(f"  internal links checked   {len(link_status)}")
    print(f"  published data files     {len(published)}")

    static_problems, search_problems = [], []
    if not args.skip_static:
        static_problems, static_pages, static_links = check_static_site(args.static)
        print(f"  static pages on disk     {static_pages}")
        print(f"  static links resolved    {static_links}")
        search_problems, searched = check_search_parity(args.static, client)
        print(f"  substrates searchable    {searched - len(search_problems)}/{searched}")

    problems = 0
    if failures:
        problems += len(failures)
        print(f"\n[FAIL] {len(failures)} page(s) did not return 200")
        for page, code in failures[:20]:
            print(f"        {code}  {page}")
    if broken:
        instances = sum(len(v) for v in broken.values())
        problems += instances
        print(f"\n[FAIL] {instances} broken internal link(s) on {len(broken)} page(s)")
        for page, items in list(broken.items())[:20]:
            for target, code in sorted(items)[:4]:
                print(f"        {code}  {target}   (linked from {page})")
    if turkish:
        problems += len(turkish)
        print(f"\n[FAIL] Turkish text visible on {len(turkish)} page(s)")
        for page, words in turkish[:20]:
            print(f"        {page}: {', '.join(words[:6])}")
    if data_turkish:
        problems += len(data_turkish)
        print(f"\n[FAIL] Turkish text in {len(data_turkish)} published data file(s)")
        for name, words in data_turkish[:20]:
            print(f"        {name}: {', '.join(words[:6])}")

    if search_problems:
        problems += len(search_problems)
        print(f"\n[FAIL] {len(search_problems)} substrate(s) cannot be found in the "
              f"published search")
        for item in search_problems[:12]:
            print(f"        {item}")

    if static_problems:
        problems += len(static_problems)
        print(f"\n[FAIL] {len(static_problems)} unresolved link(s) in the static export")
        for item in static_problems[:20]:
            print(f"        {item}")

    if problems:
        print(f"\n{problems} problem(s). Not ready to deploy.")
        return 1
    print("\nAll checks passed. Ready to deploy.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
