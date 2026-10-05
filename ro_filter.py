"""
RO alpha-subunit filtresi -- domain mimarisine dayali.

ESKI YONTEM (filter_hmm_out.py):
    hmmsearch --tblout + full-sequence E-value < 1e-10

    Sorun: ferredoksin, sitokrom bc1 Rieske ISP, NirD gibi Rieske domaini tasiyan
    ama katalitik domaini OLMAYAN proteinler bu esigi geciyor -- cunku Rieske
    domainleri gercekten homolog, E-value dogru olarak dusuk cikiyor.
    Olcum (300 kisa + 300 uzun protein): her iki grubun da %100'u E<1e-10 geciyor.
    Esigin ayirt etme gucu sifir. Tam sette 117.537 gecenin ~32.000'i alpha degil.

YENI YONTEM:
    hmmsearch --domtblout + HMM model coverage

    RieskeDB.hmm profilleri (mlen 338-587) N-terminal Rieske [2Fe-2S] domaini ile
    C-terminal katalitik mononukleer Fe(II) domainini BIRLIKTE modelliyor.
    Alpha olmayanlar profilin sadece N-terminal parcasina hizalaniyor:
        ferredoksin / NirD / bc1 ISP : medyan coverage 0.27
        gercek RO alpha              : medyan coverage 0.95
        72 referans seed             : minimum coverage 0.97

    Uzunluk filtresi bu isi YAPAMAZ: bc1 ISP (352 aa) ve b6f ISP (522 aa) gibi
    alpha-olmayan proteinler tam RO alpha uzunluk araliginda. Onlarda da katalitik
    domain yok, uzunluk transmembran ankraj gibi baska parcalardan geliyor.
"""

import os
from collections import defaultdict

# domtblout sutun indeksleri (0-tabanli)
_TARGET, _TLEN, _QUERY, _QLEN = 0, 2, 3, 5
_FULL_E, _FULL_SCORE = 6, 7
_HMM_FROM, _HMM_TO = 15, 16

# Varsayilanlar 2788 tam-boy RO alpha + 1570 kesin non-alpha uzerinde kalibre edildi.
#
#   yontem              tam-boy alpha tutulan    kesin non-alpha elenen
#   E<1e-10 (eski)              100.0%                   82.2%
#   coverage >= 0.40             99.9%                   95.8%
#   coverage >= 0.45             99.8%                   95.8%     <-- secilen
#   coverage >= 0.60             88.9%                   96.5%
#   coverage >= 0.80             79.9%                   97.4%
#
# 0.45 ustunde non-alpha elemesi neredeyse hic iyilesmiyor (95.8 -> 97.4) ama
# duyarlilik cokuyor (99.8 -> 79.9). Denge noktasi 0.45.
# 73 referans seed'in tamaminin coverage'i >= 0.97, hicbiri kaybedilmiyor.
DEFAULT_MIN_COVERAGE = 0.45
DEFAULT_MAX_EVALUE = 1e-10

# UniProt'ta "Naphthalene 1,2-dioxygenase" adiyla kayitli 110 aa'lik parcalar var.
# Bunlar yanlis pozitif DEGIL, eksik kayit -- ayri isaretlenmeli ki alpha
# olmayanlarla karistirilmasin ve downstream analize sessizce sizmasin.
DEFAULT_MIN_LENGTH = 300

# Referans setinden cikarilan modeller. Ikisi de RO alpha subunit DEGIL:
#   2_205_IsoMO  Rieske merkezi hic yok (CxH...CxxH motifi bulunamadi).
#                Cozunur di-demir monooksijenaz -- farkli bir enzim ailesi.
#                RieskeDB.hmm'de 133 BLAST hitiyle egitilmis bir modeli var.
#   1_113_CdnD   Rieske merkezi var ama katalitik triad yok. Adi zaten
#                "CdnD_electrontransfer" -- elektron transfer bileseni.
#
# Atama rekabetci (en yuksek bit skoru) oldugu icin bu modellerin hitlerini
# siniflandirma aninda elemek, HMM'i onlarsiz yeniden kosmakla ayni sonucu
# verir: protein bir sonraki en iyi modeline duser.
EXCLUDED_MODELS = frozenset({"2_205_IsoMO", "1_113_CdnD"})


def _merge_intervals(intervals):
    """Cakisan araliklari birlestir, toplam uzunlugu dondur.

    Bir protein modele birden fazla parca halinde hizalanabilir; coverage'i
    tek domainden degil, HMM ekseninde birlesik kaplamadan hesaplamak gerekir.
    """
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    return total + cur_end - cur_start + 1


def parse_domtblout(domtblout_file):
    """domtblout dosyasini oku, (target, query) cifti basina birlesik coverage hesapla.

    Donen: {target_name: [hit_dict, ...]} -- her hit bir HMM modeli icin
    coverage, skor ve E-value tasir.
    """
    envelopes = defaultdict(list)   # (target, query) -> [(hmm_from, hmm_to), ...]
    meta = {}                       # (target, query) -> (qlen, tlen, full_e, full_score)

    with open(domtblout_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 19:
                continue
            key = (fields[_TARGET], fields[_QUERY])
            envelopes[key].append((int(fields[_HMM_FROM]), int(fields[_HMM_TO])))
            meta[key] = (
                int(fields[_QLEN]),
                int(fields[_TLEN]),
                float(fields[_FULL_E]),
                float(fields[_FULL_SCORE]),
            )

    hits = defaultdict(list)
    for (target, query), spans in envelopes.items():
        qlen, tlen, full_e, full_score = meta[(target, query)]
        covered = _merge_intervals(spans)
        hits[target].append({
            "target_name": target,
            "query_name": query,
            "model_length": qlen,
            "protein_length": tlen,
            "e_value": full_e,
            "score": full_score,
            "model_coverage": covered / qlen if qlen else 0.0,
        })
    return hits


def classify_target(hits_for_target,
                    min_coverage=DEFAULT_MIN_COVERAGE,
                    max_evalue=DEFAULT_MAX_EVALUE,
                    min_length=DEFAULT_MIN_LENGTH,
                    excluded_models=EXCLUDED_MODELS):
    """Bir proteinin tum model hitlerinden RO atamasi yap.

    Coverage barajini gecen modeller arasindan EN YUKSEK SKORLU olani secer --
    en dusuk E-value'yu degil. E-value protein uzunluguna ve veritabani boyutuna
    duyarli; bit skoru kume atamasi icin daha kararli.

    Donen: (result_dict, status). status uc degerden biri:
        "RO_alpha"    filtreyi gecti
        "fragment"    coverage yeterli ama protein cok kisa (eksik kayit)
        "not_alpha"   Rieske domaini var, katalitik domain yok
    Ucunu ayirmak onemli: "fragment" bir yanlis pozitif degil, eksik veri.
    Tek bir "N/A" etiketine yikmak ikisini karistirir.
    """
    if excluded_models:
        hits_for_target = [h for h in hits_for_target
                           if h["query_name"] not in excluded_models]

    if not hits_for_target:
        return None, "no_hit"

    passing = [h for h in hits_for_target
               if h["model_coverage"] >= min_coverage and h["e_value"] < max_evalue]

    if passing:
        best = max(passing, key=lambda h: h["score"])
        status = "RO_alpha" if best["protein_length"] >= min_length else "fragment"
    else:
        best = max(hits_for_target, key=lambda h: h["score"])
        status = "not_alpha"

    is_ro = status == "RO_alpha"
    parts = best["query_name"].split("_")
    target = best["target_name"]

    result = {
        "Accession Name": target,
        "UniProt Accession": target.split("|")[1] if "|" in target else target,
        "PredictedCluster": best["query_name"] if is_ro else "N/A",
        "E-value": best["e_value"],
        "Score": best["score"],
        "ModelCoverage": round(best["model_coverage"], 4),
        "ProteinLength": best["protein_length"],
        "Status": status,
        "Group": parts[0] if is_ro and len(parts) > 0 else None,
        "Group and ID": parts[1] if is_ro and len(parts) > 1 else None,
        "Gene Abbreviation": parts[2] if is_ro and len(parts) > 2 else None,
        # Elenenlerde en iyi hitin kim oldugu -- eleme sebebini incelemek icin
        "BestModelIfRejected": None if is_ro else best["query_name"],
    }
    return result, status


def filter_domtblout(domtblout_file,
                     min_coverage=DEFAULT_MIN_COVERAGE,
                     max_evalue=DEFAULT_MAX_EVALUE,
                     min_length=DEFAULT_MIN_LENGTH,
                     output_file=None):
    """Bir domtblout dosyasindaki TUM proteinleri siniflandir.

    Tek hmmsearch cagrisinin ciktisi uzerinde calisir. Mevcut pipeline protein
    basina bir hmmsearch calistiriyor (189.657 subprocess, her biri 73 HMM'i
    bastan yukluyor); birlesik fasta uzerinde tek cagri ayni sonucu verir.

    Donen: {status: [result_dict, ...]}
    """
    hits = parse_domtblout(domtblout_file)
    by_status = defaultdict(list)
    for target_hits in hits.values():
        result, status = classify_target(target_hits, min_coverage, max_evalue, min_length)
        by_status[status].append(result)

    if output_file:
        import csv
        rows = [row for status_rows in by_status.values() for row in status_rows]
        if rows:
            write_header = not os.path.isfile(output_file)
            with open(output_file, "a", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
                if write_header:
                    writer.writeheader()
                writer.writerows(rows)

    return dict(by_status)
