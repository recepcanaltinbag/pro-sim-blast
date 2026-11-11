"""Cins sorusunun TEK SAHIBI.

NEDEN BU DOSYA VAR
  Bu DB'nin istatistik kurali "karari cins duzeyi verir" diyor, ama cins
  bes ayri dosyada ayri ayri `organism.split()[0]` ile tahmin ediliyordu
  ve ikisi eksik organizmada birbirine zit davraniyordu (biri "Unknown"u
  bir cins sayip paydaya katiyor, oteki satiri dusuruyordu). Yani ayni
  kumenin cins paydasi, hangi betigi kosturdugunuza gore degisiyordu.
  Ayni hata sinifi bu projede daha once olculmustu: bir olcumun IKINCI bir
  uygulamasi, birincisiyle sessizce ayrisir.

  `taxonomy_ids.py` cinsi NCBI rank'indan belirleyip `taxon_genus`'a
  yaziyor; bu dosya o tablonun TEK okuyucusudur. Cins soran her betik
  buradan gecer.

KULLANIM
    import ro_taxonomy
    genus = ro_taxonomy.genus_map(con)          # nucleotide_id -> Genus
    key   = genus[nuc].key                      # gruplama anahtari
    genus[nuc].name                             # gosterilecek ad
    genus[nuc].is_ranked                        # rank'tan mi geldi

NEDEN `key` ILE `name` AYRI
  614 dogrulanmis girisin (%3,2) NCBI'da cins rank'li atasi YOK -- en
  bilineni `Mycolicibacterium`, ki bugun bir altcins. Bu satirlar
  dusurulmez (payda sessizce kuculur) ve cinsleri uydurulmaz. En yakin
  adlandirilmis rank VEKIL olarak tasinir, ama `key` vekilin rank'ini
  tasidigi icin bir gruplama vekili GERCEK bir cinsle asla ayni kovaya
  koyamaz. Tahmin edilen cins bunu yapabiliyordu; bu yapamaz.
"""

from collections import namedtuple

Genus = namedtuple("Genus", "key name taxid is_ranked state proxy_rank")

#  Cinsi atanmamis satirlari dusurmek yerine GORUNUR kilan anahtar on eki.
#  Bir cins tablosunda "subgenus:Mycobacterium (subgen. Mycolicibacterium)"
#  gibi bir satir gorulurse, o satirin cins OLMADIGI okunur.
PROXY_PREFIX = "rank:"
UNKNOWN_KEY = "taxid yok"


def genus_map(con):
    """nucleotide_id -> Genus. `taxon_genus` tablosunu okur, tahmin etmez."""
    rows = con.execute("""
        SELECT nucleotide_id, genus_taxid, genus_name, genus_state,
               proxy_rank, proxy_name
          FROM taxon_genus""").fetchall()
    out = {}
    for nuc, taxid, name, state, proxy_rank, proxy_name in rows:
        if state == "ranked":
            out[nuc] = Genus(name, name, taxid, True, state, None)
        elif proxy_rank:
            key = f"{PROXY_PREFIX}{proxy_rank}:{proxy_name}"
            out[nuc] = Genus(key, proxy_name, None, False, state, proxy_rank)
        else:
            out[nuc] = Genus(UNKNOWN_KEY, None, None, False, state, None)
    return out


def genus_counts(genus, nucleotide_ids):
    """Bir giris kumesinin cins dagilimi + beyan edilmesi gereken sayilar.

    Dondurur: (ranked_counts, declaration). `declaration` cins duzeyi bir
    iddianin YANINDA yayimlanmak zorunda olan sayilardir; cagiran onu
    gormezden gelirse yuzde sessizce yanlis paydaya oturur.
    """
    from collections import Counter
    counts = Counter()
    unassigned = unknown = 0
    for nuc in nucleotide_ids:
        got = genus.get(nuc)
        if got is None or got.state == "no_taxid":
            unknown += 1
        elif got.is_ranked:
            counts[got.key] += 1
        else:
            unassigned += 1
    return counts, {
        "entries_with_ranked_genus": sum(counts.values()),
        "entries_genus_unassigned": unassigned,
        "entries_without_taxid": unknown,
        "distinct_ranked_genera": len(counts),
    }


def species_map(con):
    """nucleotide_id -> NCBI tur (rank='species') bilimsel adi.

    NEDEN. Tur de cins gibi dizgiden kesiliyordu: `organism.split()[:2]`,
    yani "ilk iki kelime". O kural GenBank dizgisine bagli ve dizgi sus
    adini da tasiyor; ustelik "Candidatus Pelagibacter ubique" gibi uc
    kelimeli adlarda ilk iki kelime TUR DEGIL cins+adlandirma durumu
    veriyor. Burada tur, cins gibi, RANK'tan okunur.
    """
    rows = con.execute("""
        SELECT t.nucleotide_id, l.ancestor_name
          FROM replicon_taxon t
          JOIN taxon_lineage l ON l.taxid = t.taxid_current
         WHERE l.rank = 'species'""").fetchall()
    return {nuc: name for nuc, name in rows}


def _by_candidate(con, per_nuc):
    """`nucleotide_id` anahtarli bir haritayi `candidate_id` anahtarina cevirir.

    NEDEN. Bazi betikler (alt-aile kesfi, yaprak homojenlestirme, ozellik
    ogrenme) replikonu hic gormez; elinde yalnizca aday kimligi ve bir
    organizma DIZGISI vardir, ve cinsi o dizgiden kesiyorlardi. Replikonu
    sorguya eklemek o betiklerin veri akisini degistirir; ayni haritayi
    aday kimligine cevirmek degistirmez. Ayni tablo, ayri anahtar.
    """
    rows = con.execute("SELECT candidate_id, nucleotide_id FROM ro").fetchall()
    return {cid: per_nuc[nuc] for cid, nuc in rows if nuc in per_nuc}


def genus_by_candidate(con):
    """candidate_id -> Genus. Eksik aday, cinsi atanmamis aday DEGILDIR:
    anahtar hic bulunmaz, boylece cagiran onu paydaya katmaz."""
    return _by_candidate(con, genus_map(con))


def species_by_candidate(con):
    """candidate_id -> NCBI tur adi (rank'tan)."""
    return _by_candidate(con, species_map(con))
