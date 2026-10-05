"""
RO alpha subunit KATALITIK MERKEZ dogrulamasi -- dogrudan kalinti testi.

Coverage filtresi (ro_filter.py) mimariyi DOLAYLI olcer: "protein profilin ne
kadarini kapliyor". Bu modul dogrudan olcer: "Rieske ligandlari ve mononukleer
demir triad'i yerinde mi".

Konumlar 73 kuratorlu referansin ROmotif.hmm'e hmmalign'lanmasindan cikarildi
(match-state kolonlari, model LENG=412):

  Rieske [2Fe-2S] merkezi (N-terminal domain) -- C-X-H ... C-X-X-H
      kolon  80   C   %99
      kolon  82   H   %99
      kolon 101   C   %100
      kolon 104   H   %99

  Mononukleer Fe(II) katalitik merkezi (C-terminal domain) -- 2-His-1-karboksilat
      kolon 209   H   %99
      kolon 214   H   %93
      kolon 345   D   %97

  Ek olarak kolon 206 (D, %88): bir alt birimin Rieske merkezini komsu alt
  birimin mononukleer demirine baglayan elektron-transfer koprusu. Zorunlu
  tutulmuyor (referanslarda bile %88), ama raporlaniyor.

Ferredoksin, sitokrom bc1 Rieske ISP ve NirD'de Rieske ligandlari VAR ama
katalitik triad YOK -- ayrimin biyokimyasal temeli tam olarak bu.
"""

from collections import Counter

# Kolon numaralari MODELE OZGUdur -- hizalama degisince kayarlar.
# Iki model icin de cikarildi; ayni biyokimyasal yapi, farkli koordinatlar.
#
#                     ROmotif.hmm (73 ref)   ROmotif71.hmm (71 ref, temiz)
#   model uzunlugu         412                      426
#   Rieske Cys-1            80  (%99)                85  (%100)
#   Rieske His-1            82  (%99)                87  (%100)
#   Rieske Cys-2           101  (%100)              105  (%100)
#   Rieske His-2           104  (%99)               108  (%100)
#   kopru Asp              206  (%88)               209  (%89)
#   katalitik His-1        209  (%99)               212  (%100)
#   katalitik His-2        214  (%93)               217  (%94)
#   katalitik karboksilat  345  (%97)               355  (%99)
#
# Temiz modelde korunma oranlari her kolonda daha yuksek -- IsoMO ve CdnD
# cikinca hizalama kesinlesti. Bagil yapi ayni: C-x-H ... C-x-x-H, iki
# katalitik His arasi 5 kalinti, karboksilat C-terminalde.
#
# YENI BIR MODEL KURARSAN bu kolonlari yeniden cikarmalisin:
#   hmmalign --trim -o refs_aln.sto YENI.hmm refs.fasta
#   ardindan match-state kolonlarinda korunum sayimi (bkz. modul sonundaki
#   derive_motif_columns fonksiyonu).

MODEL_COLUMNS = {
    "ROmotif71": {   # temiz 71-referans modeli -- varsayilan
        "rieske":    [(85, "C", "Rieske Cys-1"), (87, "H", "Rieske His-1"),
                      (105, "C", "Rieske Cys-2"), (108, "H", "Rieske His-2")],
        "catalytic": [(212, "H", "mononukleer Fe His-1"),
                      (217, "H", "mononukleer Fe His-2"),
                      (355, "DE", "mononukleer Fe karboksilat")],
        "bridging":  (209, "DE", "alt-birimler arasi elektron transfer Asp"),
    },
    "ROmotif": {     # orijinal 73-referans modeli
        "rieske":    [(80, "C", "Rieske Cys-1"), (82, "H", "Rieske His-1"),
                      (101, "C", "Rieske Cys-2"), (104, "H", "Rieske His-2")],
        "catalytic": [(209, "H", "mononukleer Fe His-1"),
                      (214, "H", "mononukleer Fe His-2"),
                      (345, "DE", "mononukleer Fe karboksilat")],
        "bridging":  (206, "DE", "alt-birimler arasi elektron transfer Asp"),
    },
}

DEFAULT_MODEL = "ROmotif71"

RIESKE_SITES = MODEL_COLUMNS[DEFAULT_MODEL]["rieske"]
CATALYTIC_SITES = MODEL_COLUMNS[DEFAULT_MODEL]["catalytic"]
BRIDGING_SITE = MODEL_COLUMNS[DEFAULT_MODEL]["bridging"]

# --- Kriterin gevsetilmesi: neden ---
#
# Katalitik triad'in 3/3'unu tam kolonda aramak FAZLA RIJIT. 73 kuratorlu
# referansin 5'i tek bir kolonda takiliyor (OxoO/CARDO/CarAa/OMO kolon 214'te K,
# PobA kolon 345'te N) -- bunlarin Rieske motifleri elle dogrulandi, gercek RO'lar.
# Sorun subaile hizalama kaymasi, kalintinin yoklugu degil.
#
# Kalibrasyon (2788 tam-boy alpha / 1570 kesin non-alpha):
#     kriter                                  alpha tutulan   non-alpha elenen
#     Rieske 4/4 + katalitik 3/3, pencere 0        83.5%            97.5%
#     Rieske 4/4 + katalitik 2/3, pencere +-2      99.7%            95.3%   <-- secilen
#
# Rieske 4/4 sarti sert tutuldu: referanslar icinde bunu sadece IsoMO kaybediyor,
# o da zaten gercek bir negatif (Rieske merkezi hic yok).
#
# Elenemeyen %4.7'nin tamami incelendi: hepsi >=300 aa (medyan 373 aa),
# "2Fe-2S ferredoxin" etiketli ama Rieske ligandlari ve katalitik triad'i tam.
# Gercek ferredoksin ~110 aa'dir; bunlar yanlis anotasyonlu RO alpha'lar.
# Yani gercek yanlis-pozitif orani pratikte ~0, tavan etiket gurultusunden.
MIN_CATALYTIC_SITES = 2
COLUMN_WINDOW = 2


def read_stockholm_matchcols(sto_file):
    """hmmalign Stockholm ciktisini oku, sadece match-state kolonlarini dondur.

    hmmalign ciktisinda match state'ler BUYUK harf veya '-', insertion'lar
    kucuk harf veya '.' ile gosterilir. Insertion kolonlari atilinca geriye
    kalan kolonlar dogrudan HMM model pozisyonlarina karsilik gelir (1-tabanli).
    """
    blocks = {}
    order = []
    with open(sto_file) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            name, chunk = parts
            if name not in blocks:
                blocks[name] = []
                order.append(name)
            blocks[name].append(chunk)

    sequences = {name: "".join(chunks) for name, chunks in blocks.items()}
    if not sequences:
        return {}

    reference = sequences[order[0]]
    match_indices = [i for i, char in enumerate(reference) if char == "-" or char.isupper()]
    return {name: "".join(seq[i] for i in match_indices)
            for name, seq in sequences.items()}


def check_sites(aligned_sequence, sites, window=COLUMN_WINDOW):
    """Verilen kolonlarda beklenen kalinti var mi. Donen: (bulunan, toplam, eksikler)

    window > 0 ise kolonun +-window komsulugunda da arar. Subaile hizalama
    kaymalarini tolere etmek icin gerekli (bkz. modul basligindaki gerekce).
    """
    found, missing = 0, []
    for column, accepted, label in sites:
        hit = False
        for offset in range(-window, window + 1):
            index = column - 1 + offset
            if 0 <= index < len(aligned_sequence) \
                    and aligned_sequence[index].upper() in accepted:
                hit = True
                break
        if hit:
            found += 1
        else:
            index = column - 1
            residue = aligned_sequence[index] if index < len(aligned_sequence) else "-"
            missing.append(f"{label}@{column}={residue}")
    return found, len(sites), missing


def classify_motifs(aligned_sequence, min_catalytic=MIN_CATALYTIC_SITES,
                    window=COLUMN_WINDOW):
    """Bir hizalanmis proteinin motif durumunu cikar.

    Donen dict:
        rieske_intact      Rieske ligandlarinin 4/4'u yerinde mi (sert sart)
        catalytic_intact   katalitik triad'in >=min_catalytic'i yerinde mi
        is_RO_alpha_motif  ikisi de saglandi mi  <-- asil karar
    """
    r_found, r_total, r_missing = check_sites(aligned_sequence, RIESKE_SITES, window)
    c_found, c_total, c_missing = check_sites(aligned_sequence, CATALYTIC_SITES, window)
    b_found, _, _ = check_sites(aligned_sequence, [BRIDGING_SITE], window)

    rieske_intact = r_found == r_total
    catalytic_intact = c_found >= min_catalytic

    return {
        "rieske_sites_found": r_found,
        "rieske_sites_total": r_total,
        "catalytic_sites_found": c_found,
        "catalytic_sites_total": c_total,
        "rieske_intact": rieske_intact,
        "catalytic_intact": catalytic_intact,
        "bridging_asp": bool(b_found),
        "is_RO_alpha_motif": rieske_intact and catalytic_intact,
        "missing_sites": ";".join(r_missing + c_missing),
    }


def motif_report(sto_file):
    """Bir hmmalign ciktisindaki tum proteinleri motif acisindan siniflandir."""
    aligned = read_stockholm_matchcols(sto_file)
    return {name: classify_motifs(seq) for name, seq in aligned.items()}


def summarize(report):
    """Konsol ozeti -- kac protein hangi kategoride."""
    counter = Counter()
    for record in report.values():
        if record["is_RO_alpha_motif"]:
            counter["Rieske + katalitik (gercek RO alpha)"] += 1
        elif record["rieske_intact"]:
            counter["sadece Rieske (ferredoksin / ISP / NirD tipi)"] += 1
        elif record["catalytic_intact"]:
            counter["sadece katalitik (Rieske kaybolmus / kismi)"] += 1
        else:
            counter["ikisi de yok"] += 1
    return counter
