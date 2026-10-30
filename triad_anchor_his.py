"""
Triad'in geri kalan iki ligandini da capala: His-1, His-2 ve capasiz tipler.

NEDEN BU DOSYA VAR
------------------
Karboksilat 2026-10-30'da sabit kolondan capali okumaya gecirildi ve "uyelerin
%12,5'inde karboksilat yok" bulgusunun %98,7'si hizalama yanilgisi cikti. Ayni
sinifin kalan iki ornegi:

  (1) Iki His hala SADECE sabit kolondan okunuyor. 55 giris orada His
      bulamiyor. Pencereyi +-5'e cikarmak His'i "buluyor" ama buldugu cift
      literaturun H-x(3-5)-H sartini tutmuyor (arasi 0 ya da 7 kalinti),
      yani YANLIS histidini buluyor. Dogru yol karboksilattakiyle ayni:
      tip-ici yeniden hizalama, kolonu tipin KENDI pozitifleri belirler.

  (2) Tipin kendi pozitifi 5'ten az oldugunda capali okuma tipi ATLIYOR ve
      elindeki KURATORLU REFERANSI hic denemiyor. 7_701_AerC tam bu durumda:
      7 uyesinin hicbirinde sabit kolonda karboksilat yok, dolayisiyla
      capa pozitifi de yok -- ama 7_701_AerC_AerC_pro refs71.fasta'da duruyor
      ve referansin karboksilat pozisyonu referans hizalamasindan BILINIYOR.

METOT
-----
Her is icin ayni uc adim:
  a) CAPA: pozitiflerin (ya da (2)'de kuratorlu referansin) ligand kalintisinin
     HAM pozisyonu, tip-ici mafft hizalamasinda hangi kolona dusuyor -- oy
     cokluguyla secilir.
  b) OKUMA: karari verilemeyen girisler o kolonun +-window komsulugunda okunur.
  c) DOGRULAMA, uc ayri sart:
       - yontem pozitiflerin kendisinde geri okunabiliyor mu (recovery)
       - AYRI TUTULAN yaridan secilen kolon, tum oyun kolonuyla ayni mi
         (dongusellik testi; karboksilatta 10/10 tipte ayni cikti)
       - His icin EK sart: bulunan ciftin HAM dizideki arasi 3-5 kalinti mi
         (literaturun H-x(3-5)-H'si). Bu sart olmadan genis pencere yanlis
         His'i kabul ediyordu; burada cagriyi o sart reddedebilir.

Cikti: analysis_out/triad_anchor_his.json. DB'ye bu dosya YAZMAZ; karari
triad_gate.py okur ve ro_triad'a isler.
"""

import argparse
import json
import os
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict

from carboxylate_audit import (AlignmentIndex, column_map, parse_stockholm,
                               read_fasta, run_mafft, write_fasta)
from ro_motif import CATALYTIC_SITES

HIS1_COL, HIS2_COL = CATALYTIC_SITES[0][0], CATALYTIC_SITES[1][0]
CARBOX_COL = CATALYTIC_SITES[2][0]
CARBOXYLATE_RESIDUES = frozenset("DE")

# Literaturun ham-dizi araligi: D-x(2)-H-x(3-5)-H (Parales/Gibson konsensusu,
# Complete_thesis.pdf bol. 3.1, atif [104]). Olculen dagilim bunu dogruluyor:
# sabit kolonda okunan 19.059 giriste ara 4 (18.445) ya da 3 (607).
#
# AMA BU ARALIK KARAR VERMEZ, SADECE KANITTIR. 2026-10-30'da olculdu: 5_504_CntA'nin
# 28 kararsiz girisinin 27'sinde iki His de tipin KENDI kolonlarinda (371 ve 379)
# duruyor -- o kolonlar 120 pozitifin 120'sinde H -- ama aradaki ara 7 kalinti,
# cunku bu alt-varyantta demir baglayan ilmekte 3 kalintilik bir insersiyon var
# (GYH-FSLSGP-CH, pozitiflerde CYH---CGPAH). Aralik sarti KARAR verseydi 27
# gercek enzim "triad eksik" diye atilacakti. Literatur iki His'in demiri
# baglamasini ister, aralarindaki ilmegin uzunlugunu SARTA BAGLAMAZ.
# Karar: kolon homolojisi. Aralik raporlanir ve kanonik olmayanlar isaretlenir.
HIS_SPACING_CANONICAL = range(3, 6)


def align(records, workdir, tag, threads, reuse):
    """mafft -- reuse ise ONCEKI ciktiyi kullan, ama girdi BIREBIR ayniysa.

    Girdiyi karsilastirmadan yeniden kullanmak bu projedeki bayatlik
    hatalarinin kaynagi; burada .fa dosyasinin icerigi byte duzeyinde
    karsilastirilir, farkliysa mafft yeniden kosar.
    """
    fasta = os.path.join(workdir, "co_%s.fa" % tag)
    out = os.path.join(workdir, "co_%s.aln" % tag)
    if reuse and os.path.exists(fasta) and os.path.getsize(out or "") > 0:
        probe = fasta + ".probe"
        write_fasta(probe, records)
        same = open(probe, "rb").read() == open(fasta, "rb").read()
        os.unlink(probe)
        if same:
            return read_fasta(out)
    return run_mafft(records, workdir, tag, threads)


def gapless_distance(row, i, j):
    low, high = sorted((i, j))
    return sum(1 for char in row[low + 1:high] if char != "-")


def read_at(row, column, accepted, window):
    """Kolon +-window icinde kalintiyi ara. Donen: (kalinti, kayma, kolon)."""
    for delta in sorted(range(-window, window + 1), key=lambda x: (abs(x), x)):
        index = column + delta
        if 0 <= index < len(row) and row[index].upper() in accepted:
            return row[index].upper(), delta, index
    return (row[column].upper() if 0 <= column < len(row) else "?"), None, None


def vote_column(aligned, anchors, maps):
    """Capalarin ham pozisyonlari hangi kolona dusuyor -- oy cokluguyla."""
    votes = Counter()
    for cid, raw in anchors:
        column = maps[cid].get(raw)
        if column is not None:
            votes[column] += 1
    if not votes:
        return None, 0, 0.0
    column, top = votes.most_common(1)[0]
    return column, top, round(top / len(anchors), 3)


def anchor_his_for_cluster(cluster, undecided, positives, entries, index, refs,
                           workdir, threads, window, pos_sample, log, reuse):
    """Bir tipin His kolonlarini kendi pozitiflerine capala ve kararsizlari oku."""
    anchors = []
    for cid in positives:
        idx = index.get(cid)
        if idx is None:
            continue
        sequence = entries[cid]["sequence"]
        raw1 = idx.find(sequence, HIS1_COL, "H", 2)[0]
        raw2 = idx.find(sequence, HIS2_COL, "H", 2)[0]
        if raw1 is not None and raw2 is not None:
            anchors.append((cid, raw1, raw2))
    if len(anchors) > pos_sample:
        step = max(1, len(anchors) // pos_sample)
        anchors = anchors[::step][:pos_sample]

    out = {"cluster": cluster, "n_undecided": len(undecided),
           "n_positive_anchors": len(anchors)}
    if len(anchors) < 5:
        out["status"] = "no_positive_anchor_available"
        return out, {cid: {"call": "untested"} for cid in undecided}

    ref_names = [k for k in refs if k.startswith(cluster + "_")]
    records = [(c, entries[c]["sequence"]) for c in undecided]
    records += [(c, entries[c]["sequence"]) for c, _, _ in anchors]
    records += [("REF|" + k, refs[k]) for k in ref_names]
    aligned = align(records, workdir, "his_" + cluster, threads, reuse)
    maps = {k: column_map(v) for k, v in aligned.items()}

    col1, top1, share1 = vote_column(aligned, [(c, r) for c, r, _ in anchors], maps)
    col2, top2, share2 = vote_column(aligned, [(c, _r) for c, _, _r in anchors], maps)
    if col1 is None or col2 is None:
        out["status"] = "anchor_column_not_resolvable"
        return out, {cid: {"call": "untested"} for cid in undecided}

    # Dongusellik testi: capalarin yarisi oy verir, otekine bakilir.
    fit = anchors[0::2]
    held = anchors[1::2]
    fit1 = vote_column(aligned, [(c, r) for c, r, _ in fit], maps)[0]
    fit2 = vote_column(aligned, [(c, _r) for c, _, _r in fit], maps)[0]
    held_ok = sum(1 for cid, _, _ in held
                  if read_at(aligned[cid], fit1, "H", window)[1] is not None
                  and read_at(aligned[cid], fit2, "H", window)[1] is not None)
    recovered = sum(1 for cid, _, _ in anchors
                    if read_at(aligned[cid], col1, "H", window)[1] is not None
                    and read_at(aligned[cid], col2, "H", window)[1] is not None)

    # Capa kolonlarinin kendi kaniti: pozitiflerin kacinda o kolon His.
    # Cagri bunun uzerine kurulur. (Pozitifler kolonu KENDI His'lerine gore
    # oyladigi icin bu oran yuksek olmak zorunda; dusukse capa guvenilmez ve
    # bu da raporlanir.)
    column_support = tuple(
        round(sum(1 for cid, _, _ in anchors
                  if aligned[cid][column].upper() == "H") / len(anchors), 3)
        for column in (col1, col2))

    calls = {}
    spacings = Counter()
    for cid in undecided:
        row = aligned[cid]
        res1, off1, at1 = read_at(row, col1, "H", window)
        res2, off2, at2 = read_at(row, col2, "H", window)
        if off1 is None or off2 is None:
            calls[cid] = {"call": "absent", "his1": res1, "his2": res2,
                          "spacing": None}
            spacings["bulunamadi"] += 1
            continue
        spacing = gapless_distance(row, at1, at2)
        spacings[spacing] += 1
        canonical = spacing in HIS_SPACING_CANONICAL
        # Karar: His TAM capa kolonunda mi (kayma 0). Capa kolonu homolog
        # pozisyonun kendisidir; komsu bir kolondaki His bir KAYMADIR ve
        # tek basina yetmez -- orada pencere bitisik bir HH'yi yakalayabilir
        # (1_106_TsaM ve 6_601_PAO'da tam bu oldu: capa kolonunda N/D/R var,
        # pencere yanindaki His'i buluyor ve ara 0 ya da 1 cikiyor). Kaymis
        # bir okuma ancak literaturun H-x(3-5)-H araligi destekliyorsa cagri
        # uretir.
        exact = off1 == 0 and off2 == 0
        if exact:
            call = "present_at_homologous_column"
        elif canonical:
            call = "present_at_shifted_column"
        else:
            call = "absent"
        calls[cid] = {"call": call, "his1": res1, "his2": res2,
                      "spacing": spacing, "offset": (off1, off2),
                      "at_exact_anchor_column": exact,
                      "canonical_spacing": canonical,
                      # Kolonun KENDI kanidi: ayni kolonda kac pozitif His
                      # tasiyor. Cagri buna dayanir, araliga degil.
                      "anchor_column_his_fraction": [column_support[0],
                                                     column_support[1]]}

    out.update({
        "status": "tested",
        "curated_reference_in_alignment": ref_names,
        "anchor_columns": [col1, col2],
        "anchor_column_vote_share": [share1, share2],
        "method_recovery_on_positives": round(recovered / len(anchors), 3),
        "holdout": {"n_fit": len(fit), "n_held": len(held),
                    "fit_columns": [fit1, fit2],
                    "fit_equals_full_vote": [fit1 == col1, fit2 == col2],
                    "held_out_recovery": (round(held_ok / len(held), 3)
                                          if held else None)},
        "reference_residues_at_anchor": {
            k: (aligned["REF|" + k][col1], aligned["REF|" + k][col2])
            for k in ref_names},
        "anchor_column_his_fraction_in_positives": list(column_support),
        "spacing_counts": {str(k): v for k, v in spacings.most_common()},
        "called_at_exact_anchor_column": sum(
            1 for c in calls.values() if c.get("at_exact_anchor_column")),
        "non_canonical_spacing": sum(
            1 for c in calls.values()
            if c.get("spacing") is not None and not c["canonical_spacing"]),
        "calls": dict(Counter(c["call"] for c in calls.values()).most_common()),
    })
    log("  %-14s n=%-3d capa kol %d/%d (pozitifte His %s)  geri okuma %5.1f%%  "
        "tutulan yari %s  ara %s  -> %s"
        % (cluster, len(undecided), col1, col2, list(column_support),
           100 * out["method_recovery_on_positives"],
           out["holdout"]["held_out_recovery"],
           out["spacing_counts"], out["calls"]))
    return out, calls


def anchor_carboxylate_to_reference(cluster, undecided, entries, refs, ref_index,
                                    workdir, threads, window, log, reuse):
    """Pozitifi olmayan tipin karboksilatini KURATORLU REFERANSA capala."""
    ref_names = [k for k in refs if k.startswith(cluster + "_")]
    out = {"cluster": cluster, "n_undecided": len(undecided),
           "curated_reference": ref_names}
    usable = [(k, ref_index[k]) for k in ref_names if ref_index.get(k) is not None]
    if not usable:
        out["status"] = "reference_carboxylate_position_unknown"
        return out, {cid: {"call": "untested"} for cid in undecided}

    records = [(c, entries[c]["sequence"]) for c in undecided]
    records += [("REF|" + k, refs[k]) for k in ref_names]
    aligned = align(records, workdir, "ref_" + cluster, threads, reuse)
    maps = {k: column_map(v) for k, v in aligned.items()}

    votes = Counter()
    for name, raw in usable:
        column = maps["REF|" + name].get(raw)
        if column is not None:
            votes[column] += 1
    if not votes:
        out["status"] = "anchor_column_not_resolvable"
        return out, {cid: {"call": "untested"} for cid in undecided}
    anchor_col = votes.most_common(1)[0][0]

    calls, residues = {}, Counter()
    for cid in undecided:
        res, delta, _ = read_at(aligned[cid], anchor_col, CARBOXYLATE_RESIDUES,
                                window)
        residues[res] += 1
        calls[cid] = ({"call": "present_at_reference_anchored_column",
                       "residue": res, "offset": delta} if delta is not None
                      else {"call": "absent", "residue": res, "offset": None})
    out.update({
        "status": "tested",
        "anchor_source": "curated_reference",
        "reference_raw_positions": {k: r for k, r in usable},
        "anchor_column": anchor_col,
        "reference_residue_at_anchor": {
            k: aligned["REF|" + k][anchor_col] for k in ref_names},
        "residues_read": dict(residues.most_common()),
        "calls": dict(Counter(c["call"] for c in calls.values()).most_common()),
    })
    log("  %-14s n=%-3d REFERANS capasi kol %d  okunan %s  -> %s"
        % (cluster, len(undecided), anchor_col, out["residues_read"],
           out["calls"]))
    return out, calls


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--alignment", default="genomic_context/cand_aln.sto")
    ap.add_argument("--ref-alignment", default="ROs_71_Clean/refs71_hmmaln.sto")
    ap.add_argument("--ref-fasta", default="ROs_71_Clean/refs71.fasta")
    ap.add_argument("--out", default="analysis_out/triad_anchor_his.json")
    ap.add_argument("--window", type=int, default=3)
    ap.add_argument("--pos-sample", type=int, default=120)
    ap.add_argument("--mafft-threads", type=int, default=1,
                    help="VARSAYILAN 1: cok parcacikli mafft bit duzeyinde "
                         "tekrarlanabilir degil")
    ap.add_argument("--work", default=None)
    ap.add_argument("--reuse-alignments", action="store_true",
                    help="onceki mafft ciktisini kullan, ama girdi fasta'si "
                         "BIREBIR ayniysa")
    args = ap.parse_args()

    def log(message):
        print(message, flush=True)

    workdir = args.work or tempfile.mkdtemp(prefix="triad_anchor_")
    os.makedirs(workdir, exist_ok=True)
    log("[gecici dizin] %s" % workdir)

    con = sqlite3.connect("file:%s?mode=ro" % args.db, uri=True)
    entries = {r[0]: {"sequence": r[1], "cluster": r[2]} for r in con.execute(
        "SELECT candidate_id, sequence, ro_cluster FROM ro WHERE is_confirmed = 1")}
    try:
        triad = {r[0]: (r[1], r[2]) for r in con.execute(
            "SELECT candidate_id, triad_status, carboxylate_source FROM ro_triad")}
    except sqlite3.OperationalError:
        log("[HATA] ro_triad yok. Once: python3 triad_gate.py --apply")
        return 2
    con.close()
    log("[DB] onaylanmis giris %d, ro_triad satiri %d" % (len(entries), len(triad)))

    his_undecided = sorted(c for c, (s, _) in triad.items()
                           if s == "his_needs_anchoring")
    carbox_untested = sorted(c for c, (s, _) in triad.items()
                             if s == "carboxylate_untested")
    log("[is] His karari verilemeyen %d giris, karboksilati test edilemeyen %d"
        % (len(his_undecided), len(carbox_untested)))
    if not his_undecided and not carbox_untested:
        log("[bilgi] capalanacak giris yok.")
        return 0

    log("[okunuyor] %s" % args.alignment)
    rows, rf = parse_stockholm(args.alignment, set(entries))
    match_idx = [i for i, ch in enumerate(rf) if ch != "."]
    index = {}
    for cid, row in rows.items():
        idx = AlignmentIndex(row, match_idx, entries[cid]["sequence"])
        if idx.offset >= 0:
            index[cid] = idx
    log("[bilgi] ham diziye eslenen %d giris" % len(index))

    refs = read_fasta(args.ref_fasta)
    log("[okunuyor] referans hizalamasi %s" % args.ref_alignment)
    ref_rows, ref_rf = parse_stockholm(args.ref_alignment)
    ref_match = [i for i, ch in enumerate(ref_rf) if ch != "."]
    ref_carbox_raw = {}
    for name, row in ref_rows.items():
        if name not in refs:
            continue
        idx = AlignmentIndex(row, ref_match, refs[name])
        if idx.offset < 0:
            continue
        raw = idx.find(refs[name], CARBOX_COL, CARBOXYLATE_RESIDUES, 2)[0]
        if raw is not None:
            ref_carbox_raw[name] = raw
    log("[bilgi] karboksilat pozisyonu bilinen referans: %d / %d"
        % (len(ref_carbox_raw), len(refs)))

    his_positives = defaultdict(list)
    for cid, (status, _) in triad.items():
        if status != "his_needs_anchoring" and cid in index:
            his_positives[entries[cid]["cluster"]].append(cid)

    by_cluster = defaultdict(list)
    for cid in his_undecided:
        by_cluster[entries[cid]["cluster"]].append(cid)

    log("\n[adim 1] His kolonlarini tip-ici capala")
    his_reports, his_calls = [], {}
    for cluster in sorted(by_cluster, key=lambda c: -len(by_cluster[c])):
        report, calls = anchor_his_for_cluster(
            cluster, sorted(by_cluster[cluster]),
            sorted(his_positives.get(cluster, [])), entries, index, refs,
            workdir, args.mafft_threads, args.window, args.pos_sample, log,
            args.reuse_alignments)
        his_reports.append(report)
        his_calls.update(calls)

    log("\n[adim 2] capasiz tiplerin karboksilatini KURATORLU REFERANSA capala")
    ref_reports, ref_calls = [], {}
    ref_by_cluster = defaultdict(list)
    for cid in carbox_untested:
        ref_by_cluster[entries[cid]["cluster"]].append(cid)
    for cluster in sorted(ref_by_cluster):
        report, calls = anchor_carboxylate_to_reference(
            cluster, sorted(ref_by_cluster[cluster]), entries, refs,
            ref_carbox_raw, workdir, args.mafft_threads, args.window, log,
            args.reuse_alignments)
        ref_reports.append(report)
        ref_calls.update(calls)

    his_tally = Counter(c["call"] for c in his_calls.values())
    ref_tally = Counter(c["call"] for c in ref_calls.values())
    log("\n[sonuc] His capali okuma: %s" % dict(his_tally.most_common()))
    log("[sonuc] referans capali karboksilat: %s" % dict(ref_tally.most_common()))

    out = {
        "generated": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc).isoformat(timespec="seconds"),
        "module": "triad_anchor_his.py",
        "writes_to_db": False,
        "inputs": {"db": args.db, "alignment": args.alignment,
                   "reference_alignment": args.ref_alignment,
                   "columns": {"his1": HIS1_COL, "his2": HIS2_COL,
                               "carboxylate": CARBOX_COL},
                   "window_columns": args.window,
                   "his_spacing_canonical": [min(HIS_SPACING_CANONICAL),
                                             max(HIS_SPACING_CANONICAL)],
                   "decision_rule": "column homology decides; the His-His "
                                    "spacing is reported, not required",
                   "mafft_threads": args.mafft_threads},
        "histidine_anchoring": {"per_cluster": his_reports,
                                "call_counts": dict(his_tally.most_common())},
        "carboxylate_reference_anchoring": {
            "per_cluster": ref_reports,
            "call_counts": dict(ref_tally.most_common())},
        "per_entry_calls": {
            "histidine": {cid: his_calls[cid] for cid in sorted(his_calls)},
            "carboxylate": {cid: ref_calls[cid] for cid in sorted(ref_calls)}},
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(out, handle, indent=1)
    log("\n[yazildi] %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
