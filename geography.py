"""Cografya: her enzim NERELERDE goruldu, ve bu soru ne kadar sorulabilir.

Kullanicinin istegi acikti: enzim sayfalarinda "nerelerde goruldu" haritasi,
varyantlarin yayilimi. Veri `replicon_source.country` alaninda duruyor ve
121 ayri ulke/deniz adiyla temiz sayilir. Asil is veriyi cikarmak degil, onu
YANILTICI OLMAYACAK bicimde sunmak.

Olculen seyin ne OLMADIGI onemli. Ham sayim, enzimin nerede yasadigini degil,
hangi ulkenin kac genom DIZDIGINI olcer: ABD tek basina girislerin dortte
birini tasiyor ve bu, Rieske oksijenazlarin Amerika'da yogunlastigi anlamina
gelmez. Bu yuzden sayfa uc ayri sayi tasir:

  ham giris sayisi      -- en yaniltici olani, yalnizca eksiksizlik icin
  ayri CINS sayisi      -- ayni susun tekrari sayilmaz
  arka plana gore sapma -- asil olcu. Bir tipin ulke dagilimi, BUTUN
                           dogrulanmis girislerin ulke dagilimiyla
                           karsilastirilir. "Bu enzim, Rieske oksijenazlarin
                           genelinden farkli yerlerde mi bulunuyor?" sorusunun
                           cevabi budur ve orneklem yanliligi iki tarafta da
                           ayni olduğu icin buyuk olcude sadelesir.

Sapma, istatistik sayfasiyla AYNI araci kullanir (stats_overview.signal_map):
standartlastirilmis Pearson artigi, Benjamini-Hochberg duzeltmesi, ve maddi
fark esigi. Iki sayfada iki ayri yontem olmasin diye tek uygulama paylasilir.

Cikti: analysis_out/geography.json
"""

import argparse
import json
import os
import sqlite3
from collections import Counter, defaultdict

from stats_overview import bh_adjust, signal_map  # noqa: F401  (tek uygulama)

# ISO-3 kodlari: koroplet harita bunlarla cizilir. Deniz ve tarihsel devletler
# bilerek disarida: bir harita yuzeyine oturmuyorlar ve ayri listelenirler.
ISO3 = {
    "Afghanistan": "AFG", "Algeria": "DZA", "Antarctica": "ATA",
    "Argentina": "ARG", "Armenia": "ARM", "Australia": "AUS", "Austria": "AUT",
    "Bangladesh": "BGD", "Belarus": "BLR", "Belgium": "BEL", "Bolivia": "BOL",
    "Brazil": "BRA", "Cambodia": "KHM", "Cameroon": "CMR", "Canada": "CAN",
    "Chile": "CHL", "China": "CHN", "Colombia": "COL", "Costa Rica": "CRI",
    "Croatia": "HRV", "Czech Republic": "CZE", "Denmark": "DNK",
    "Dominican Republic": "DOM", "Ecuador": "ECU", "Egypt": "EGY",
    "Eritrea": "ERI", "Ethiopia": "ETH", "Finland": "FIN", "France": "FRA",
    "French Guiana": "GUF", "Georgia": "GEO", "Germany": "DEU", "Ghana": "GHA",
    "Greece": "GRC", "Guadeloupe": "GLP", "Guatemala": "GTM", "Haiti": "HTI",
    "Honduras": "HND", "Hong Kong": "HKG", "Hungary": "HUN", "Iceland": "ISL",
    "India": "IND", "Indonesia": "IDN", "Iran": "IRN", "Iraq": "IRQ",
    "Ireland": "IRL", "Israel": "ISR", "Italy": "ITA", "Jamaica": "JAM",
    "Japan": "JPN", "Jordan": "JOR", "Kazakhstan": "KAZ", "Kenya": "KEN",
    "Korea": "KOR", "Kuwait": "KWT", "Laos": "LAO", "Lebanon": "LBN",
    "Malaysia": "MYS", "Mayotte": "MYT", "Mexico": "MEX", "Mongolia": "MNG",
    "Morocco": "MAR", "Mozambique": "MOZ", "Myanmar": "MMR", "Namibia": "NAM",
    "Netherlands": "NLD", "New Zealand": "NZL", "Nigeria": "NGA",
    "Northern Mariana Islands": "MNP", "Norway": "NOR", "Oman": "OMN",
    "Pakistan": "PAK", "Palau": "PLW", "Panama": "PAN",
    "Papua New Guinea": "PNG", "Peru": "PER", "Philippines": "PHL",
    "Poland": "POL", "Portugal": "PRT", "Puerto Rico": "PRI",
    "Reunion": "REU", "Romania": "ROU", "Russia": "RUS",
    "Saudi Arabia": "SAU", "Senegal": "SEN", "Singapore": "SGP",
    "Slovakia": "SVK", "Somalia": "SOM", "South Africa": "ZAF",
    "South Korea": "KOR", "Spain": "ESP", "Sudan": "SDN", "Svalbard": "SJM",
    "Sweden": "SWE", "Switzerland": "CHE", "Taiwan": "TWN",
    "Tanzania": "TZA", "Thailand": "THA", "Trinidad and Tobago": "TTO",
    "Tunisia": "TUN", "Turkey": "TUR", "USA": "USA", "Uganda": "UGA",
    "Ukraine": "UKR", "United Arab Emirates": "ARE", "United Kingdom": "GBR",
    "Uruguay": "URY", "Uzbekistan": "UZB", "Venezuela": "VEN",
    "Viet Nam": "VNM", "Virgin Islands": "VIR", "Zimbabwe": "ZWE",
}

# Bir ulke yuzeyine oturmayan kayitlar. Atilmazlar -- ekolojik olarak en
# ilginc olanlar arasindalar -- ama haritada degil, kendi listelerinde durur.
NON_TERRITORIAL = {
    "Arctic Ocean": "ocean", "Atlantic Ocean": "ocean",
    "Baltic Sea": "sea", "Indian Ocean": "ocean",
    "Mediterranean Sea": "sea", "North Sea": "sea", "Pacific Ocean": "ocean",
}

# Artik var olmayan devletler: kaydin tarihi oldugu gibi korunur, cunku
# toplama yili da saklaniyor ve kaydi "duzeltmek" veriyi tahrif etmek olur.
HISTORICAL = {"USSR": "dissolved 1991", "Yugoslavia": "dissolved 1992"}

MIN_ENTRIES_FOR_TYPE = 20     # bunun altindaki tip icin ulke dagilimi anlamsiz
MIN_ENTRIES_FOR_COUNTRY = 15  # bunun altindaki ulke sapma testine girmez


def load(con):
    rows = con.execute("""
        SELECT r.candidate_id, r.ro_cluster, r.ro_group, p.organism,
               s.country, s.geo, s.collection_year, s.habitat,
               e.tier, l.leaf_id
        FROM ro r
        JOIN replicon p USING(nucleotide_id)
        JOIN replicon_source s ON s.nucleotide_id = r.nucleotide_id
        LEFT JOIN ro_evidence e ON e.candidate_id = r.candidate_id
        LEFT JOIN ro_leaf l ON l.candidate_id = r.candidate_id
        WHERE r.is_confirmed = 1
          AND s.country IS NOT NULL AND s.country <> ''
    """).fetchall()
    out = []
    for cid, cluster, group, organism, country, geo, year, habitat, tier, leaf in rows:
        out.append({
            "candidate_id": cid, "cluster": cluster, "group": group or "?",
            "genus": (organism or "?").split()[0],
            "country": country, "geo": geo, "year": year,
            "habitat": habitat or "unknown", "tier": tier or "unknown",
            "leaf": leaf,
        })
    return out


def place_kind(name):
    if name in NON_TERRITORIAL:
        return NON_TERRITORIAL[name]
    if name in HISTORICAL:
        return "historical"
    if name in ISO3:
        return "country"
    return "unmapped"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out", default="analysis_out/geography.json")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    data = load(con)
    total_confirmed = con.execute(
        "SELECT COUNT(*) FROM ro WHERE is_confirmed = 1").fetchone()[0]

    # ---------------------------------------------------------------- ulkeler
    per_country = defaultdict(lambda: {"entries": 0, "genera": set(),
                                       "types": set(), "habitats": Counter(),
                                       "years": []})
    for d in data:
        c = per_country[d["country"]]
        c["entries"] += 1
        c["genera"].add(d["genus"])
        c["types"].add(d["cluster"])
        c["habitats"][d["habitat"]] += 1
        if d["year"]:
            c["years"].append(d["year"])

    countries = []
    for name, c in sorted(per_country.items(), key=lambda kv: -kv[1]["entries"]):
        years = sorted(c["years"])
        countries.append({
            "name": name,
            "iso3": ISO3.get(name),
            "kind": place_kind(name),
            "entries": c["entries"],
            "genera": len(c["genera"]),
            "types": len(c["types"]),
            # Ayni cinsin tekrari sayilmasin diye ikinci bir olcu: cins basina
            # kac giris dustugu. Buyuk sayi, derin orneklenmis birkac susu
            # gosterir, cesitliligi degil.
            "entries_per_genus": round(c["entries"] / len(c["genera"]), 2),
            "top_habitat": c["habitats"].most_common(1)[0][0] if c["habitats"] else None,
            "year_range": [years[0], years[-1]] if years else None,
        })

    # --------------------------------------------------- tip x ulke sapmasi
    # Arka plan: BUTUN dogrulanmis girislerin ulke dagilimi. Bir tipin kendi
    # dagilimi buna gore olculur; dizileme cabasi iki tarafta da ayni oldugu
    # icin buyuk olcude sadelesir.
    type_counts = Counter(d["cluster"] for d in data)
    big_types = sorted(t for t, n in type_counts.items()
                       if n >= MIN_ENTRIES_FOR_TYPE)
    country_counts = Counter(d["country"] for d in data)
    big_countries = [c for c, n in country_counts.most_common()
                     if n >= MIN_ENTRIES_FOR_COUNTRY]

    subset = [d for d in data if d["cluster"] in big_types
              and d["country"] in big_countries]
    table = [[sum(1 for d in subset if d["cluster"] == t and d["country"] == c)
              for c in big_countries] for t in big_types]

    # Cins duzeyi: ayni tip + ayni cins + ayni ulke tek gozlem sayilir. Tek
    # bir projenin yuzlerce susu bir hucreyi tek basina anlamli yapmasin.
    seen = set()
    gsubset = []
    for d in subset:
        key = (d["cluster"], d["genus"], d["country"])
        if key in seen:
            continue
        seen.add(key)
        gsubset.append(d)
    gtable = [[sum(1 for d in gsubset if d["cluster"] == t and d["country"] == c)
               for c in big_countries] for t in big_types]

    geo_signals = signal_map(big_types, big_countries, table, gtable)

    by_type = {}
    for cluster in sorted(type_counts):
        rows = [d for d in data if d["cluster"] == cluster]
        counts = Counter(d["country"] for d in rows)
        genera_by_country = defaultdict(set)
        for d in rows:
            genera_by_country[d["country"]].add(d["genus"])
        by_type[cluster] = {
            "entries_with_a_country": len(rows),
            "n_countries": len(counts),
            "countries": [
                {"name": name, "iso3": ISO3.get(name), "kind": place_kind(name),
                 "entries": n, "genera": len(genera_by_country[name])}
                for name, n in counts.most_common()],
        }

    # ------------------------------------------------------------ varyantlar
    # "Ana varyant nerede goruldu, uzak akrabalari nerede?" Yaprak (leaf)
    # kimligi varyant yerine gecer: ayni yaprak = ayni dizi varyanti.
    variants = {}
    for cluster in big_types:
        rows = [d for d in data if d["cluster"] == cluster and d["leaf"]]
        if not rows:
            continue
        leaves = Counter(d["leaf"] for d in rows)
        main_leaf, main_n = leaves.most_common(1)[0]
        main_rows = [d for d in rows if d["leaf"] == main_leaf]
        other_rows = [d for d in rows if d["leaf"] != main_leaf]
        variants[cluster] = {
            "n_variants": len(leaves),
            "main_variant": main_leaf,
            "main_variant_entries": main_n,
            "main_variant_countries": [
                {"name": k, "iso3": ISO3.get(k), "entries": v}
                for k, v in Counter(d["country"] for d in main_rows).most_common()],
            "other_variant_countries": [
                {"name": k, "iso3": ISO3.get(k), "entries": v}
                for k, v in Counter(d["country"] for d in other_rows).most_common()],
            "countries_unique_to_other_variants": sorted(
                {d["country"] for d in other_rows} - {d["country"] for d in main_rows}),
        }

    out = {
        "coverage": {
            "confirmed_entries": total_confirmed,
            "entries_with_a_country": len(data),
            "share": round(len(data) / total_confirmed, 4) if total_confirmed else None,
            "distinct_places": len(per_country),
            "places_on_the_map": sum(1 for c in countries if c["kind"] == "country"),
            "places_off_the_map": [c["name"] for c in countries
                                   if c["kind"] != "country"],
            "what_is_missing": (
                "entries whose source record carries no /geo_loc_name or /country "
                "qualifier. Absence is not evidence of anything: it is a metadata "
                "gap, and it is not random -- older submissions carry it far more "
                "often than recent ones."),
        },
        "caveats": [
            "A country here is where the genome was collected and deposited, not "
            "where the enzyme lives. The map is in large part a map of sequencing "
            "capacity, which is why raw counts are reported beside the number of "
            "distinct genera and why the only inference drawn is the deviation "
            "from the overall distribution.",
            "The deviation test compares each type's country distribution against "
            "the country distribution of all confirmed entries, so sampling effort "
            "largely cancels. It is repeated after collapsing each type, genus and "
            "country to one observation, and only cells that survive both are "
            "reported.",
            "Oceans, seas and dissolved states are kept as recorded and listed "
            "separately rather than being forced onto a national boundary.",
        ],
        "countries": countries,
        "by_type": by_type,
        "variants": variants,
        "deviation": {
            "question": ("Is this enzyme type found in a different set of places "
                         "than Rieske oxygenases in general?"),
            "types_tested": big_types,
            "countries_tested": big_countries,
            "min_entries_for_type": MIN_ENTRIES_FOR_TYPE,
            "min_entries_for_country": MIN_ENTRIES_FOR_COUNTRY,
            "table": table,
            "genus_table": gtable,
            "signals": geo_signals,
        },
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"[yazildi] {args.out}")
    print(f"  {len(data)} of {total_confirmed} confirmed entries carry a country "
          f"({100 * len(data) / total_confirmed:.1f} %)")
    print(f"  {len(per_country)} distinct places, "
          f"{out['coverage']['places_on_the_map']} of them on the map")
    print(f"  off the map: {', '.join(out['coverage']['places_off_the_map']) or 'none'}")
    if geo_signals:
        print(f"  deviation: {len(big_types)} types x {len(big_countries)} countries, "
              f"{geo_signals['n_cells_material']} of {geo_signals['n_cells_tested']} "
              f"cells material, {geo_signals['n_cells_material_genus']} hold at "
              f"genus level")
        for c in geo_signals["top"][:8]:
            print(f"    {c['row']:24s} {c['col']:18s} "
                  f"{100 * c['share']:5.1f} % vs {100 * c['baseline']:5.1f} % "
                  f"({100 * c['lift']:+5.1f} pts){'  holds' if c['holds'] else ''}")


if __name__ == "__main__":
    main()
