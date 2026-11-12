"""Taksonomi KONTROLLU SOZLUK olur: NCBI taxid + rank'li soyagaci.

NEDEN BU DOSYA VAR
  Bu veritabaninin kendi istatistik kurali soyle diyor: bir iddia giris,
  dizi ve CINS duzeyinde test edilir, ve KARARI CINS DUZEYI VERIR. Yani
  yayimlanan her sinyalin birimi cins. Ama cinsin kendisi bugune kadar
  hicbir yerde TANIMLI degildi; bes ayri dosyada ayri ayri, su satirla
  tahmin ediliyordu:

      genus = organism.split()[0]

  Bu bir siniflandirma degil, bir dizgi islemi. Ve olculdu, hata uretiyor:

  1. "Candidatus" bir cins DEGIL, adlandirma durumudur. `split()[0]`
     onu cins sanip 17 AYRI gercek cinsi (Micropelagos, Rhodobacter,
     Pelagibacter, Puniceispirillum, ...) tek kovaya yiginiyor: 57 giris.
     Birlestirme yonu de en kotusu -- ayri cinsleri tek cins gibi
     gostermek, bir sinyalin tek cinste YOGUNLASTIGI izlenimini guclendirir,
     ki cins muhafizinin varlik sebebi tam olarak bunu engellemek.
  2. "uncultured" ayni sekilde cins sayiliyor.
  3. Cinsi dizgiden tahmin etmek, kaydin kendi soyagacini kullanmamak
     demek; GenBank dosyasinin icinde o soyagaci ZATEN duruyor.

  Ve iki tahminci birbiriyle de tutarsizdi: `analyze.py` eksik organizmaya
  "Unknown" deyip onu BIR CINS olarak sayiyor, `analyze_ecology.py` None
  deyip satiri DUSURUYOR. Yani ayni kumenin cins paydasi, hangi betigi
  kosturdugunuza gore degisiyordu.

NEDEN RANK, SON ALAN DEGIL
  Akla gelen ilk duzeltme, `replicon.taxonomy` soyagacinin SON alanini
  cins saymak. Olculdu: 17.074 replikonun 2.006'sinda (%11,7) o alan cins
  degil, bir TUR KOMPLEKSI ("Mycobacterium tuberculosis complex",
  "Bacillus cereus group", "Ralstonia solanacearum species complex").
  Soyagacinda rank etiketi YOK; alan sirasi rank anlamina gelmiyor.
  Dolayisiyla cinsi rank'tan okumak icin NCBI taxonomy dump'i gerekir:
  `nodes.dmp` her taxid'in rank'ini yazar.

BIRIM VE PAYDA
  `replicon_taxon` satiri = REPLIKON (bir nucleotide_id). Cins duzeyi
  istatistik bu tabloyu `ro.nucleotide_id` uzerinden baglar; bir girisin
  cinsi, uzerinde durdugu replikonun cinsidir.

KAYNAK
  taxid: her GenBank kaydinin `source` ozelliginde `/db_xref="taxon:N"`.
  Yerel, aga cikmaz, 17.074/17.074 dosyada var.
  rank + ad + soyagaci: NCBI taxonomy dump (nodes.dmp, names.dmp,
  merged.dmp, delnodes.dmp), md5'i dogrulanarak indirilir.

Cikti:
    sqlite  replicon_taxon(nucleotide_id, taxid, taxid_current, name, rank)
    sqlite  taxon_lineage(taxid, rank, ancestor_taxid, ancestor_name)
    analysis_out/taxonomy_ids.json
"""

import argparse
import glob
import gzip
import json
import os
import re
import sqlite3
from collections import Counter

TAXON_XREF = re.compile(rb'/db_xref="taxon:(\d+)"')

#  Yayimlanan soyagacinda tutulan rank'lar. NCBI'da 40'tan fazla rank var
#  ("subspecies", "clade", "no rank", ...); hepsini tasimak tabloyu
#  sisirir ve hicbiri yayimlanmiyor. Bunlar yayimlananlar.
PUBLISHED_RANKS = ("superkingdom", "domain", "phylum", "class",
                   "order", "family", "genus", "species")

#  CINS SAYILMAYACAK adlar. Bunlar cins degil, adlandirma durumu ya da
#  yer tutucu; bir tanesi bile cins kovasina girerse o kova birbiriyle
#  ilgisiz organizmalari toplar.
NOT_A_GENUS = ("candidatus", "uncultured", "unclassified", "candidate",
               "environmental samples", "incertae sedis")

#  CINSI OLMAYAN takson icin vekil rank sirasi. NCBI'nin guncel agacinda
#  bazi turlerin yolunda cins rank'li dugum HIC YOK: olculdu, 614
#  dogrulanmis giris boyle. Ornegi `Mycolicibacterium fortuitum`:
#      tur -> ALTCINS "Mycobacterium (subgen. Mycolicibacterium)" -> AILE
#  Yani `Mycolicibacterium` bugunun NCBI taksonomisinde bir cins degil,
#  bir altcins; dizgi tahmini onu cins sanip bu DB'nin 7. en kalabalik
#  "cinsini" uretiyordu. Bu satirlar DUSURULMEZ (payda sessizce kuculur)
#  ve cins UYDURULMAZ; `genus_state='unassigned'` ile tasinir ve en yakin
#  adlandirilmis rank VEKIL olarak yazilir. Bir test ya vekille gruplar
#  ya da disarida birakir; ikisini de BEYAN ETMEK zorunda.
PROXY_RANKS = ("subgenus", "family", "order", "class", "phylum")


def read_taxids(gbk_dir):
    """Her GenBank kaydinin basligindan taxid'i okur.

    Yalnizca BASLIK okunur: taxid `source` ozelliginde, ilk CDS'ten once
    gecer. 14.092 sikistirilmamis + 2.982 gz dosyada tam dosyayi acmak
    bedavaya yakin degil.
    """
    found, missing = {}, []
    paths = sorted(glob.glob(os.path.join(gbk_dir, "*.gbk"))
                   + glob.glob(os.path.join(gbk_dir, "*.gbk.gz")))
    for path in paths:
        base = os.path.basename(path)
        nuc = base[:-7] if base.endswith(".gbk.gz") else base[:-4]
        opener = gzip.open if path.endswith(".gz") else open
        head = b""
        try:
            with opener(path, "rb") as handle:
                for _ in range(400):
                    line = handle.readline()
                    if not line:
                        break
                    head += line
                    if line.startswith(b"     gene") or line.startswith(b"     CDS"):
                        break
        except Exception as exc:                      # bozuk/kesik dosya
            missing.append((nuc, f"okunamadi: {exc}"))
            continue
        match = TAXON_XREF.search(head)
        if match:
            found[nuc] = int(match.group(1))
        else:
            missing.append((nuc, "taxon db_xref yok"))
    return found, missing


def read_taxdump(dump_dir):
    """nodes.dmp + names.dmp + merged.dmp + delnodes.dmp okur."""
    nodes = {}
    with open(os.path.join(dump_dir, "nodes.dmp"), encoding="utf-8") as handle:
        for line in handle:
            parts = line.split("\t|\t")
            nodes[int(parts[0])] = (int(parts[1]), parts[2].strip())
    names = {}
    with open(os.path.join(dump_dir, "names.dmp"), encoding="utf-8") as handle:
        for line in handle:
            parts = line.split("\t|\t")
            if parts[3].startswith("scientific name"):
                names[int(parts[0])] = parts[1].strip()
    merged = {}
    path = os.path.join(dump_dir, "merged.dmp")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                parts = line.split("\t|")
                merged[int(parts[0])] = int(parts[1])
    deleted = set()
    path = os.path.join(dump_dir, "delnodes.dmp")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                deleted.add(int(line.split("\t|")[0]))
    return nodes, names, merged, deleted


def lineage_of(taxid, nodes, limit=80):
    """taxid'den koke dogru taxid listesi (kendisi dahil)."""
    chain, seen = [], set()
    cur = taxid
    while cur and cur not in seen and len(chain) < limit:
        chain.append(cur)
        seen.add(cur)
        parent = nodes.get(cur, (None, None))[0]
        if parent == cur:                             # kok kendi ebeveyni
            break
        cur = parent
    return chain


def resolve(taxid, nodes, merged, deleted):
    """Birlestirilmis taxid'i gunceline tasir. Dondurur: (guncel, durum)."""
    if taxid in nodes:
        return taxid, "current"
    if taxid in merged and merged[taxid] in nodes:
        return merged[taxid], "merged"
    if taxid in deleted:
        return None, "deleted"
    return None, "unknown"


def write_tables(con, rows, lineages, genus_rows):
    con.executescript("""
        DROP TABLE IF EXISTS replicon_taxon;
        DROP TABLE IF EXISTS taxon_lineage;
        CREATE TABLE replicon_taxon (
            nucleotide_id  TEXT PRIMARY KEY REFERENCES replicon(nucleotide_id),
            taxid          INTEGER,   -- GenBank kaydinda yazan
            taxid_current  INTEGER,   -- merged.dmp ile tasinmis hali
            status         TEXT,      -- current | merged | deleted | unknown
            name           TEXT,
            rank           TEXT
        );
        CREATE TABLE taxon_lineage (
            taxid          INTEGER,   -- replicon_taxon.taxid_current
            rank           TEXT,      -- PUBLISHED_RANKS icinden
            ancestor_taxid INTEGER,
            ancestor_name  TEXT,
            PRIMARY KEY (taxid, rank)
        );
        DROP TABLE IF EXISTS taxon_genus;
        CREATE TABLE taxon_genus (
            nucleotide_id TEXT PRIMARY KEY REFERENCES replicon(nucleotide_id),
            taxid_current INTEGER,
            genus_taxid   INTEGER,  -- rank'i tam olarak 'genus' olan ata
            genus_name    TEXT,
            genus_state   TEXT,     -- ranked | unassigned | no_taxid
            proxy_rank    TEXT,     -- unassigned ise en yakin adlandirilmis rank
            proxy_name    TEXT
        );
        CREATE INDEX idx_rt_taxid ON replicon_taxon(taxid_current);
        CREATE INDEX idx_tl_rank  ON taxon_lineage(rank, ancestor_name);
        CREATE INDEX idx_tg_genus ON taxon_genus(genus_state, genus_name);
    """)
    con.executemany("INSERT INTO replicon_taxon VALUES (?,?,?,?,?,?)", rows)
    con.executemany("INSERT INTO taxon_lineage VALUES (?,?,?,?)", lineages)
    con.executemany("INSERT INTO taxon_genus VALUES (?,?,?,?,?,?,?)", genus_rows)
    con.commit()


def nearest_ranks(taxid, nodes, names):
    """taxid -> {rank: (ata_taxid, ata_adi)}; her rank'in EN YAKIN atasi.

    NEDEN AYRI FONKSIYON. Cins kurali bu dosyanin sahipliginde, ama cins
    sorusunu yalnizca replikonlar sormuyor: kapsama denetimi (coverage_audit)
    DIS evrenin taxid'lerini de cinse cevirmek zorunda ve o taxid'ler bu
    DB'nin hicbir tablosunda yok. Ikinci bir cozucu yazmak, bu dosyanin
    varlik sebebi olan hatayi (ayni olcumun ikinci uygulamasi birincisiyle
    sessizce ayrisir) aynen geri getirirdi.
    """
    out = {}
    for anc in lineage_of(taxid, nodes):
        rank = nodes.get(anc, (None, None))[1]
        if rank in PUBLISHED_RANKS or rank in PROXY_RANKS:
            #  lineage_of koke dogru gidiyor: ilk gorulen EN YAKIN olandir.
            out.setdefault(rank, (anc, names.get(anc)))
    return out


def genus_from_ranks(ranks):
    """Cins KARARI. Dondurur: (state, genus_taxid, genus_name, proxy_rank,
    proxy_name). Kurali tek yerde tutan fonksiyon budur; cins soran her
    yol -- replikon tablosu da, dis evren karsilastirmasi da -- bundan gecer.
    """
    got = ranks.get("genus")
    if got:
        return "ranked", got[0], got[1], None, None
    for rank in PROXY_RANKS:
        if rank in ranks:
            return "unassigned", None, None, rank, ranks[rank][1]
    return "unassigned", None, None, None, None


def resolve_genus(rows, lineages):
    """Her replikon icin cinsi RANK'tan belirler; yoksa vekili yazar."""
    by_taxid = {}
    for taxid, rank, anc, anc_name in lineages:
        by_taxid.setdefault(taxid, {})[rank] = (anc, anc_name)
    out = []
    for nuc, _taxid, current, _status, _name, _rank in rows:
        if current is None:
            out.append((nuc, None, None, None, "no_taxid", None, None))
            continue
        state, gtaxid, gname, proxy_rank, proxy_name = genus_from_ranks(
            by_taxid.get(current, {}))
        out.append((nuc, current, gtaxid, gname, state, proxy_rank, proxy_name))
    return out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--gbk", default="gbk_files")
    parser.add_argument("--dump", default="/media/lin-bio/back2/roar_taxonomy",
                        help="nodes.dmp/names.dmp'in bulundugu dizin")
    parser.add_argument("--out", default="analysis_out/taxonomy_ids.json")
    args = parser.parse_args()

    found, missing = read_taxids(args.gbk)
    print(f"[taxid] {len(found)} replikon, {len(missing)} eksik")
    nodes, names, merged, deleted = read_taxdump(args.dump)
    print(f"[dump]  {len(nodes)} dugum, {len(names)} ad, {len(merged)} birlestirilmis")

    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    known = {r[0] for r in con.execute("SELECT nucleotide_id FROM replicon")}

    rows, status_count = [], Counter()
    lineage_rows, needed = [], {}
    for nuc, taxid in sorted(found.items()):
        if nuc not in known:                 # DB'de olmayan dosya: atlanir
            continue
        current, status = resolve(taxid, nodes, merged, deleted)
        status_count[status] += 1
        name = names.get(current) if current else None
        rank = nodes.get(current, (None, None))[1] if current else None
        rows.append((nuc, taxid, current, status, name, rank))
        if current:
            needed[current] = True

    for taxid in needed:
        for anc in lineage_of(taxid, nodes):
            rank = nodes.get(anc, (None, None))[1]
            if rank in PUBLISHED_RANKS or rank in PROXY_RANKS:
                lineage_rows.append((taxid, rank, anc, names.get(anc)))
    #  Ayni rank bir soyagacinda iki kere gecebilir (ornegin iki "clade"
    #  degil ama NCBI bazi dallarda "no rank" araya girince tekrar
    #  olabiliyor). En YAKIN olani tutulur; lineage_of koke dogru gittigi
    #  icin ilk gorulen en yakindir.
    seen, deduped = set(), []
    for taxid, rank, anc, anc_name in lineage_rows:
        if (taxid, rank) in seen:
            continue
        seen.add((taxid, rank))
        deduped.append((taxid, rank, anc, anc_name))

    genus_rows = resolve_genus(rows, deduped)
    write_tables(con, rows, deduped, genus_rows)

    #  OLCUM: rank'tan okunan cins, dizgiden tahmin edilen cinsle nerede
    #  ayrisiyor? Bu sayi duzeltmenin BEDELI degil, BOYUTU; ve iki yonu
    #  ayri sayilir, cunku ikisi farkli hatalar:
    #    changed     -> dizgi baska bir cins adi veriyordu (yeniden
    #                   adlandirma, ornegin Pseudomonas -> Stutzerimonas)
    #    unassigned  -> dizginin cins sandigi sey bugun cins DEGIL
    by_nuc = {r[0]: r for r in genus_rows}
    per_entry = con.execute("""
        SELECT ro.nucleotide_id AS nuc, rep.organism AS organism, COUNT(*) AS n
          FROM ro JOIN replicon rep USING (nucleotide_id)
         WHERE ro.is_confirmed = 1 GROUP BY 1, 2""").fetchall()
    total = changed = no_genus = 0
    fake = Counter()
    proxy_ranks = Counter()
    examples, unassigned_examples = [], []
    for row in per_entry:
        total += row["n"]
        split_guess = (row["organism"] or "").split()
        split_guess = split_guess[0] if split_guess else None
        got = by_nuc.get(row["nuc"])
        state = got[4] if got else "no_taxid"
        ranked = got[3] if got else None
        if split_guess and split_guess.lower() in NOT_A_GENUS:
            fake[split_guess] += row["n"]
        if state != "ranked":
            no_genus += row["n"]
            proxy_ranks[(got[5] if got else None) or "yok"] += row["n"]
            if len(unassigned_examples) < 10:
                unassigned_examples.append(
                    {"organism": row["organism"], "string_guess": split_guess,
                     "nearest_named_rank": got[5] if got else None,
                     "nearest_named_name": got[6] if got else None,
                     "entries": row["n"]})
        elif ranked != split_guess:
            changed += row["n"]
            if len(examples) < 15:
                examples.append({"organism": row["organism"],
                                 "string_guess": split_guess,
                                 "ranked_genus": ranked, "entries": row["n"]})

    report = {
        "method": {
            "question": "What changes when each entry's genus is decided by "
                        "NCBI rank instead of by guessing at the organism name?",
            "unit": "one replicon (replicon_taxon); a genus-level statistic "
                    "joins through ro.nucleotide_id",
            "genus_definition": "the nearest ancestor whose rank is exactly "
                                "'genus' in the NCBI lineage (nodes.dmp)",
            "not_a_genus": list(NOT_A_GENUS),
            "taxid_source": 'the GenBank source feature, /db_xref="taxon:N"',
            "taxdump_files": ["nodes.dmp", "names.dmp", "merged.dmp",
                              "delnodes.dmp"],
        },
        "denominators": {
            "replicons_in_db": len(known),
            "replicons_with_taxid": len(rows),
            "replicons_without_taxid": len(known) - len(rows),
            "confirmed_entries": total,
        },
        "taxid_status": dict(status_count),
        "distinct_taxids": len(needed),
        "distinct_ranked_genera": len({r[3] for r in genus_rows if r[3]}),
        "genus_change": {
            "entries_whose_genus_changes": changed,
            "share_changed": round(changed / total, 4) if total else None,
            "examples": examples,
        },
        "genus_unassigned": {
            "entries": no_genus,
            "share": round(no_genus / total, 4) if total else None,
            "handling": "Neither dropped nor invented: carried as "
                        "genus_state='unassigned' with the nearest named rank "
                        "recorded as a proxy. A genus-level test either groups "
                        "by that proxy or leaves the entry out, and either way "
                        "it is obliged to declare this number.",
            "nearest_named_rank": dict(proxy_ranks),
            "examples": unassigned_examples,
        },
        "fake_genus_buckets_removed": dict(fake),
        #  "Candidatus" kovasinin ICINDE kac GERCEK cins vardi. Bu sayi,
        #  dizgi tahmincisinin hatasinin YONUNU gosterir: ayri cinsleri
        #  tek cins gibi gostermek, bir sinyalin tek cinste yogunlastigi
        #  izlenimini guclendirir -- cins muhafizinin engellemek icin var
        #  oldugu seyin tam kendisi.
        "candidatus_genera": len({
            (r["organism"] or "").split()[1]
            for r in per_entry
            if (r["organism"] or "").startswith("Candidatus")
            and len((r["organism"] or "").split()) > 1}),
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)

    print(f"[tablo] replicon_taxon {len(rows)}, taxon_lineage {len(deduped)}")
    print(f"[cins]  rank'li tekil cins {len({r[3] for r in genus_rows if r[3]})}, "
          f"degisen giris {changed}/{total}, cinsi atanmamis {no_genus}")
    print(f"[sahte] {dict(fake)}")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
