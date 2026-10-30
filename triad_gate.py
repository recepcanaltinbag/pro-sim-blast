"""
2-His-1-karboksilat triad'inin 3/3 kapisi -- literaturun tanimiyla.

NEDEN BU DOSYA VAR
------------------
Mononukleer Fe(II) merkezi literaturde istisnasiz "2-His-1-karboksilat facial
triad" olarak tanimlanir (Barry & Challis 2013, ACS Catal 3:2362; Ferraro 2005,
Biochem Biophys Res Commun 338:175). Karboksilat sus degil: NDO'da D362A mutanti
TAMAMEN inaktiftir. Yani uc liganddan ikisini isteyen bir kabul kurali
(MIN_CATALYTIC_SITES = 2) literaturle uyumlu DEGILDIR.

Buna ragmen kapi uzun sure 2/3'te kaldi, cunku 3/3'u SABIT kolonda uygulamak
2388 giris (%12.5) ve 45 tip kaybettiriyordu. 2026-10-30 denetimi
(analysis_out/carboxylate_audit.json) o 2388'in 2357'sinde karboksilatin
FIILEN VAR oldugunu, sadece sabit kolonda okunamadigini olctu: eksik
cagrilarin %81.5'i kolonun ve komsularinin hizalamada BOS olmasindan geliyor.
Yani sorun kuralda degil OLCUMDE idi.

BU DOSYANIN YAPTIGI
-------------------
Triad'in uc uyesini AYRI AYRI, her biri icin dogru okumayla cagirir:

    His-1, His-2   ROmotif71 sabit kolonu +-2    (iki kolon da iyi korunmus)
    karboksilat    (a) sabit kolon +-2, yoksa
                   (b) tip-ici capali okuma: ayni tipin karboksilati bilinen
                       uyelerine mafft ile hizalanip homolog kolondan okunur

(b) denetimin 6. adiminda olculdu ve capraz dogrulandi: tutulan yarida
bulunan kolon, tum oylamanin kolonuyla 10/10 tipte ayni; ayri tutulan yarida
geri okuma 490/490. Yontem bu dosyada tekrar kosmaz, denetimin per_entry_calls
blogundan okunur -- ve okumadan once TAZELIK kontrolu yapilir.

CIKTI
-----
ro_triad          giris basina uc ligand + 3/3 karari + karboksilatin kaynagi
ro_triad_meta     tazelik parmak izi (hangi kolonlar, hangi denetim, kac giris)

--apply verilmezse HICBIR SEY yazilmaz; sadece ne olacagi bildirilir.
--apply ile ro_triad yazilir. ro.is_confirmed'e DOKUNULMAZ: yayinlanmis
denominatorlari sessizce kaydirmamak icin karar once ayri tabloda durur
(bkz. ingest_protein_entries.py'deki ayni gerekce).
"""

import argparse
import json
import os
import re
import sqlite3
from collections import Counter, defaultdict
from datetime import datetime, timezone

from ro_motif import CATALYTIC_SITES, read_stockholm_matchcols

# Kume, dizgi DEGIL. Dizgi olarak yazilinca ("DE") bos okuma sessizce
# "var" sayiliyordu: Python'da `"" in "DE"` True'dur. Bu hata AerC'nin 7
# girisini "triad tam" diye isaretledi (2026-10-30'da yakalandi).
CARBOXYLATE_RESIDUES = frozenset("DE")
HIS_WIDE_WINDOW = 5

# Parales/Gibson konsensusu: karboksilat ve iki histidin ham dizide
# D-x(2)-H-x(3-5)-H olarak dizilir (Complete_thesis.pdf, bol. 3.1, atif [104]).
# Burada His CAGRISINI DOGRULAMAK icin kullanilir: sabit kolon +-2'de His
# bulunamayan girislerde pencereyi genisletmek His'i "buluyor" ama BULDUGU
# His yanlis olanidir -- bkz. asagidaki olcum.
TRIAD_CONSENSUS = re.compile(r"[DE].{2}H.{3,5}H")

# Denetimin capali okuma cagrilari. "present_*" = karboksilat bulundu.
PRESENT_CALLS = ("present_at_homologous_column",
                 "present_at_group_conserved_column",
                 "present_at_shifted_column",
                 "present_at_reference_anchored_column")


def read_site(row, column, accepted, window=2):
    """Kolonun +-window komsulugunda beklenen kalintiyi ara.

    Donen: (kalinti, kayma) -- bulunamazsa (kolondaki kalinti, None).
    """
    for delta in sorted(range(-window, window + 1), key=lambda x: (abs(x), x)):
        index = column - 1 + delta
        if 0 <= index < len(row) and row[index].upper() in accepted:
            return row[index].upper(), delta
    index = column - 1
    return (row[index].upper() if 0 <= index < len(row) else "-"), None


def his_pair_spacing(row, his1_col, his2_col, window):
    """Pencerede bulunan iki His arasindaki KALINTI sayisi (gapler atilir).

    Literatur H-x(3-5)-H diyor, yani bu sayi 3-5 olmalidir. Sabit kolonda
    His bulunamayan girislerde pencere genisletilince bu sayi 0 ya da 7
    cikiyor: pencere bitisik bir HH'yi ya da baska bir His'i yakaliyor,
    triad'in His'ini DEGIL. Bu yuzden genis pencere bir cagri uretmez,
    sadece "neden cagri yapilmadi"yi belgeler.
    """
    found = []
    for column in (his1_col, his2_col):
        for delta in sorted(range(-window, window + 1), key=lambda x: (abs(x), x)):
            index = column - 1 + delta
            if 0 <= index < len(row) and row[index].upper() == "H":
                found.append(index)
                break
        else:
            return None
    low, high = sorted(found)
    return sum(1 for char in row[low + 1:high] if char != "-")


def check_audit_freshness(audit, confirmed_ids, fixed_negatives):
    """Denetim bu DB ve bu kolonlar icin mi uretilmis.

    Yanlis kolonda uretilmis bir denetimi sessizce kullanmak bu projede bir
    kere oldu (kolon 355 / 471 / 431 vakasi); bir daha olmasin diye burada
    gurultulu basarisizlik var.
    """
    problems = []
    columns = audit.get("inputs", {}).get("columns", {}).get("catalytic")
    want = [list(site) for site in CATALYTIC_SITES]
    if columns != want:
        problems.append("denetim katalitik kolonlari %s, ro_motif %s diyor"
                        % (columns, want))

    n_confirmed = audit.get("inputs", {}).get("confirmed_entries")
    if n_confirmed != len(confirmed_ids):
        problems.append("denetim %s onaylanmis giris gormus, DB'de %d var"
                        % (n_confirmed, len(confirmed_ids)))

    calls = audit.get("per_entry_calls") or {}
    if not calls:
        problems.append("denetimde per_entry_calls blogu yok")
    else:
        missing = fixed_negatives - set(calls)
        extra = set(calls) - fixed_negatives
        if missing:
            problems.append("sabit kolonda karboksilatsiz %d giris denetimde yok "
                            "(ornek: %s)" % (len(missing), sorted(missing)[0]))
        if extra:
            problems.append("denetimde DB'nin negatifi olmayan %d giris var "
                            "(ornek: %s)" % (len(extra), sorted(extra)[0]))
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    ap.add_argument("--audit", default="analysis_out/carboxylate_audit.json")
    ap.add_argument("--his-anchor", default="analysis_out/triad_anchor_his.json",
                    help="triad_anchor_his.py ciktisi; yoksa His SADECE sabit "
                         "kolondan okunur ve karari verilemeyenler "
                         "'his_needs_anchoring' kalir")
    ap.add_argument("--out", default="analysis_out/triad_gate.json")
    ap.add_argument("--apply", action="store_true",
                    help="ro_triad tablosunu yaz (verilmezse sadece bildirir)")
    args = ap.parse_args()

    his1_col, his1_ok, his1_label = CATALYTIC_SITES[0]
    his2_col, his2_ok, his2_label = CATALYTIC_SITES[1]
    carbox_col, carbox_ok, carbox_label = CATALYTIC_SITES[2]
    print("[kolonlar] %s=%d  %s=%d  %s=%d"
          % (his1_label, his1_col, his2_label, his2_col, carbox_label, carbox_col))

    con = sqlite3.connect(args.db)
    rows = dict(con.execute("SELECT candidate_id, ro_cluster FROM ro "
                            "WHERE is_confirmed = 1"))
    stored = dict(con.execute(
        "SELECT c.candidate_id, c.catalytic_residue FROM ro_carboxylate c "
        "JOIN ro r ON r.candidate_id = c.candidate_id WHERE r.is_confirmed = 1"))
    print("[DB] onaylanmis giris: %d   ro_carboxylate satiri: %d"
          % (len(rows), len(stored)))

    fixed_negatives = {cid for cid in rows
                       if (stored.get(cid) or "") not in CARBOXYLATE_RESIDUES}
    print("[DB] sabit kolonda karboksilat YOK: %d" % len(fixed_negatives))

    with open(args.audit) as handle:
        audit = json.load(handle)
    problems = check_audit_freshness(audit, set(rows), fixed_negatives)
    if problems:
        print("\n[HATA] denetim bu DB/kolon seti icin TAZE DEGIL:")
        for item in problems:
            print("   - %s" % item)
        print("\ncarboxylate_audit.py yeniden kosulmali. Hicbir sey yazilmadi.")
        return 2
    print("[tazelik] denetim %s, kolonlar ve giris sayisi uyusuyor"
          % audit["generated"])

    print("[okunuyor] %s" % args.alignment)
    aligned = read_stockholm_matchcols(args.alignment)
    print("[bilgi] hizalamada %d dizi" % len(aligned))

    # His'in capali okumasi. Yoksa His sabit kolonda kalir ve karari
    # verilemeyen girisler "his_needs_anchoring" olarak isaretlenir --
    # hicbir durumda "His yok" diye yayinlanmaz.
    his_anchor, his_anchor_stamp = {}, None
    if os.path.exists(args.his_anchor):
        with open(args.his_anchor) as handle:
            doc = json.load(handle)
        want = {"his1": his1_col, "his2": his2_col, "carboxylate": carbox_col}
        if doc.get("inputs", {}).get("columns") != want:
            print("[HATA] %s baska kolonlar icin uretilmis: %s"
                  % (args.his_anchor, doc.get("inputs", {}).get("columns")))
            return 2
        his_anchor = doc["per_entry_calls"]["histidine"]
        his_anchor.update(doc["per_entry_calls"].get("carboxylate_his") or {})
        his_anchor_stamp = doc["generated"]
        print("[His capasi] %s, %d giris icin cagri var"
              % (his_anchor_stamp, len(his_anchor)))
    else:
        print("[His capasi] yok (%s); His sadece sabit kolondan okunacak"
              % args.his_anchor)

    calls = audit["per_entry_calls"]
    out_rows = []
    tally = Counter()
    source_tally = Counter()
    by_type_incomplete = Counter()
    status_by_type = defaultdict(Counter)
    his_spacing_seen = Counter()
    canonical_spacing = Counter()
    for cid, cluster in rows.items():
        row = aligned.get(cid)
        if row is None:
            tally["hizalamada yok"] += 1
            out_rows.append((cid, cluster, None, None, None, "unaligned",
                             None, None, 0, "unaligned",
                             "absent from the alignment"))
            continue

        h1, h1_off = read_site(row, his1_col, his1_ok)
        h2, h2_off = read_site(row, his2_col, his2_ok)
        his_source = "fixed_column"
        his_note = ""
        if h1_off is None or h2_off is None:
            # Sabit kolonda His yok. Pencereyi genisletmek His'i "buluyor" ama
            # bulunan cift literaturun aralik sartini tutmuyor (asagida olculur),
            # bu yuzden burada YOKLUK DA ILAN EDILMEZ: giris capali okumaya
            # devredilir. Karboksilatta yapilan hatanin ayni sinifi.
            spacing = his_pair_spacing(row, his1_col, his2_col, HIS_WIDE_WINDOW)
            his_spacing_seen[spacing] += 1
            anchored = his_anchor.get(cid)
            if anchored is None:
                his_source = "needs_anchoring"
                his_note = ("no histidine at the fixed column; the pair found "
                            "within +-%d columns is separated by %s residues "
                            "(literature: 3-5), so no call was made"
                            % (HIS_WIDE_WINDOW, spacing))
            elif anchored["call"] in PRESENT_CALLS:
                his_source = "anchored"
                h1, h2 = anchored["his1"], anchored["his2"]
                h1_off = h2_off = anchored.get("offset", (0, 0))[0]
                his_note = ("both histidines read at this type's own anchored "
                            "columns; the pair is separated by %s residues"
                            % anchored.get("spacing"))
                if not anchored.get("canonical_spacing", True):
                    his_note += (" -- longer than the usual 3-5, an insertion "
                                 "in the iron-binding loop, and the two columns "
                                 "carry a histidine in %s of this type's "
                                 "members with a fixed-column reading"
                                 % anchored.get("anchor_column_his_fraction"))
            else:
                his_source = "anchored"
                h1_off = h2_off = None
                # Capa kolonunda OTURAN kalinti. Pencerede bulunani yazmak
                # "his_absent" satirlarinda kalinti olarak H gosteriyordu.
                h1 = anchored.get("his1_at_column", anchored["his1"])
                h2 = anchored.get("his2_at_column", anchored["his2"])
                his_note = ("no histidine at this type's own anchored columns "
                            "(the type's positives carry one in %s of cases), "
                            "so the site is genuinely incomplete"
                            % anchored.get("anchor_column_his_fraction"))

        if his_source == "fixed_column":
            canonical_spacing[his_pair_spacing(row, his1_col, his2_col, 2)] += 1

        if cid in fixed_negatives:
            call = calls[cid]
            if call["anchored_call"] in PRESENT_CALLS:
                carbox = call["anchored_residue"]
                source = ("anchored_homologous"
                          if call["anchored_call"] == "present_at_homologous_column"
                          else "anchored_group_conserved")
                offset = call["anchored_offset"]
                note = ""
            elif call["anchored_call"] == "untested":
                carbox, source, offset = call["anchored_residue"], "untested", None
                note = ("fewer than 5 anchor positives in this type, and the "
                        "type's curated reference carries no carboxylate at "
                        "this column either, so the column is not homologous "
                        "here: the absence is NOT established")
            else:
                carbox, source, offset = call["anchored_residue"], "absent", None
                note = "not found even after the anchored re-reading"
        else:
            carbox, offset = read_site(row, carbox_col, carbox_ok)
            source, note = "fixed_column", ""

        has_carbox = (carbox or "") in CARBOXYLATE_RESIDUES
        his_ok = h1_off is not None and h2_off is not None
        complete = int(his_ok and has_carbox)

        if complete:
            status = "complete"
        elif not his_ok:
            status = ("his_absent" if his_source == "anchored"
                      else "his_needs_anchoring")
        elif source == "untested":
            status = "carboxylate_untested"
        else:
            status = "carboxylate_absent"
        tally[status] += 1
        source_tally[source] += 1
        if not complete:
            by_type_incomplete[cluster] += 1
            status_by_type[status][cluster] += 1
        note = "; ".join(x for x in (note, his_note) if x)
        out_rows.append((cid, cluster, h1, h2, carbox, source, offset,
                         his_source, complete, status, note))

    print("\n[sonuc] %s" % dict(tally.most_common()))
    print("[karboksilatin kaynagi] %s" % dict(source_tally.most_common()))
    print("[His1->His2 kalinti arasi, sabit kolonda okunanlar] %s"
          % dict(canonical_spacing.most_common(5)))
    print("[His1->His2 kalinti arasi, genis pencerede bulunanlar] %s   "
          "<- literatur 3-5 diyor, yani genis pencere YANLIS His'i buluyor; "
          "karar capali okumadan geliyor" % dict(his_spacing_seen.most_common()))
    print("\n3/3'u tamamlamayan tipler (durumuna gore):")
    for status in sorted(status_by_type, key=lambda s: -sum(status_by_type[s].values())):
        print("   %s  (%d giris)" % (status, sum(status_by_type[status].values())))
        for cluster, count in status_by_type[status].most_common():
            print("      %-16s %d" % (cluster, count))

    print("\n[dikkat] karboksilatin yoklugu %d giriste KANITLANMADI, His okumasi "
          "%d giriste hala capalanmamis. Bu %d girisin hicbiri 'triad yok' diye "
          "yayinlanamaz; kapi siki tutulmadan once cozulmeleri gerekir."
          % (tally["carboxylate_untested"], tally["his_needs_anchoring"],
             tally["carboxylate_untested"] + tally["his_needs_anchoring"]))

    summary = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "module": "triad_gate.py",
        "rule": "Rieske 4/4 (annotate_ro) + mononuclear Fe triad 3/3",
        "carboxylate_reading": "fixed column +-2; if not found there, the "
                               "within-type anchored re-reading",
        "his_reading": "fixed column +-2; if no histidine is there, the "
                       "within-type anchored reading, where the call needs "
                       "the exact anchored column (a shifted hit must also "
                       "match the literature's H-x(3-5)-H spacing)",
        "audit_used": {"path": args.audit, "generated": audit["generated"]},
        "his_anchor_used": ({"path": args.his_anchor,
                             "generated": his_anchor_stamp}
                            if his_anchor_stamp else None),
        "columns": {"his1": his1_col, "his2": his2_col, "carboxylate": carbox_col},
        "confirmed_entries": len(rows),
        "triad_complete": tally["complete"],
        "triad_incomplete": sum(v for k, v in tally.items() if k != "complete"),
        "status_counts": dict(tally.most_common()),
        "carboxylate_source_counts": dict(source_tally.most_common()),
        "incomplete_by_type": dict(by_type_incomplete.most_common()),
        "incomplete_by_status_and_type": {s: dict(c.most_common())
                                          for s, c in status_by_type.items()},
        "his_pair_spacing_fixed_column": dict(canonical_spacing.most_common()),
        "his_pair_spacing_wide_window": dict(his_spacing_seen.most_common()),
        "his_wide_window_verdict":
            "Widening the window to +-%d does find a histidine in the entries "
            "that have none at the fixed column, but the pair it finds is not "
            "separated by the 3-5 residues the literature requires "
            "(H-x(3-5)-H): the window catches an adjacent HH instead. These "
            "entries are therefore neither complete nor absent; they need the "
            "same within-type anchored reading the carboxylate got."
            % HIS_WIDE_WINDOW,
        "absence_not_proven": tally["carboxylate_untested"],
        "his_absent": tally["his_absent"],
        "his_undecided": tally["his_needs_anchoring"],
        "applied": bool(args.apply),
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(summary, handle, indent=1)
    print("\n[yazildi] %s" % args.out)

    if not args.apply:
        print("[kuru kosu] ro_triad YAZILMADI. --apply ile yaz.")
        con.close()
        return 0

    con.executescript("""
        DROP TABLE IF EXISTS ro_triad;
        CREATE TABLE ro_triad (
            candidate_id        TEXT PRIMARY KEY,
            ro_cluster          TEXT,
            his1_residue        TEXT,
            his2_residue        TEXT,
            carboxylate_residue TEXT,
            carboxylate_source  TEXT,  -- fixed_column | anchored_* | absent | untested
            carboxylate_offset  INTEGER,
            his_source          TEXT,  -- fixed_column | needs_anchoring
            triad_complete      INTEGER,
            triad_status        TEXT,  -- complete | carboxylate_absent |
                                       -- carboxylate_untested | his_needs_anchoring
            note                TEXT
        );
        DROP TABLE IF EXISTS ro_triad_meta;
        CREATE TABLE ro_triad_meta (key TEXT PRIMARY KEY, value TEXT);
    """)
    con.executemany("INSERT INTO ro_triad VALUES (?,?,?,?,?,?,?,?,?,?,?)", out_rows)
    con.execute("CREATE INDEX idx_triad_complete ON ro_triad(triad_complete)")
    con.execute("CREATE INDEX idx_triad_source ON ro_triad(carboxylate_source)")
    con.execute("CREATE INDEX idx_triad_status ON ro_triad(triad_status)")
    con.executemany("INSERT INTO ro_triad_meta VALUES (?,?)",
                    [(k, json.dumps(v)) for k, v in summary.items()])
    con.commit()
    con.close()
    print("[yazildi] ro_triad (%d satir) + ro_triad_meta" % len(out_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
