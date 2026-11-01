"""`ro_protein`'deki UniProt kayitlari BUGUN hala var mi -- olculur.

NEDEN BU DOSYA VAR
  `ro_protein` 73.699 satir tutuyor ve bu sayi `about.html`'de yayimlaniyor.
  Tablo, PF00355 setinin bir indirmesinden kuruldu (madde 95: indirme tarihi
  kayitli). 2026-11-01'de 300'luk rastgele bir ornek UniProt'a soruldu ve
  baglantisiz girislerin **%63'u** `Inactive` cikti; sebebi sorulan 38
  kaydin 38'i `DELETED` -- BIRLESTIRME DEGIL, yani yeni bir accession
  altinda yasamiyorlar. UniProt aradan gecen surede gereksiz proteomlari
  buyuk olcude sildi.

  Ornek bir tahmindir; bu betik TAMAMINI sorar. Silinmis kayitlar
  SESSIZCE ATILMAZ: "silinmis kayit" ile "veri yok" ayri seylerdir, ve
  34 bin satir aciklamasiz yok olursa yayimlanmis her protein-duzeyi sayi
  anlamini degistirir. O yuzden ayri bir tabloya DURUM yazilir, satirlar
  yerinde kalir.

NEDEN AYRI TABLO, `ro_protein`'e SUTUN DEGIL
  `ro_protein`'i ingest_protein_entries.py uretiyor ve o AGIR bir adim
  (hmmalign + motif taramasi). Durumu oraya sutun olarak eklemek, bir ag
  sorgusunu tazelemek icin butun motif olcumunu yeniden kosturmak demekti.
  Ayri tablo, tazeligi kendi girdisine baglar.

NE TOPLANIR
  entry_type      : "UniProtKB unreviewed (TrEMBL)" | "UniProtKB reviewed
                    (Swiss-Prot)" | "Inactive" | "missing" (hic donmedi)
  inactive_reason : DELETED | MERGED | DEMERGED (yalnizca Inactive icin)
  merged_into     : birlestirildiyse yeni accession -- bu satirlar
                    KAYIP DEGIL, yalnizca kimlik degistirmis
  embl_nuc        : EMBL/ENA nukleotit accession'i (varsa, ilk taneси)
                    -- genomik baglamin kurulabilecegi yolun kaniti
  refseq          : RefSeq protein accession'i (varsa)
  ec              : EC numarasi (varsa) -- C maddesi icin, literaturle
                    baglanti kuran standart tanimlayici

Kullanim:
    python3 protein_freshness.py --db roar.sqlite --out-dir analysis_out
    python3 protein_freshness.py --limit 400        # deneme kosusu
"""

import argparse
import collections
import json
import os
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

API = "https://rest.uniprot.org/uniprotkb/search"
# Yigin boyu OLCULDU (2026-11-01): 100 accession calisiyor, 200 HTTP 400
# donduruyor -- sinir URL uzunlugu degil (3.130 karakterde sorun yok), sorgu
# yan tumcesi SAYISI. 100'de kaliyoruz.
BATCH = 100
FIELDS = "accession,protein_name,ec,xref_embl,xref_refseq,organism_name"


def fetch(accessions, retries=4, pause=0.4):
    """Bir yigin accession icin UniProt kayitlarini dondur.

    Ag hatasinda ustel bekleyerek yeniden dener. Butun yigin kalici olarak
    basarisiz olursa None doner -- cagiran taraf bunu "olculemedi" diye
    kaydeder, "silinmis" diye DEGIL. Ikisini karistirmak, ag kesintisini
    veri kaybi gibi gostermek olurdu.
    """
    query = " OR ".join(f"accession:{a}" for a in accessions)
    url = (f"{API}?format=json&size=500&fields={FIELDS}"
           f"&query={urllib.parse.quote(query)}")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as handle:
                return json.load(handle).get("results", [])
        except (urllib.error.URLError, TimeoutError, ValueError) as exc:
            if attempt == retries - 1:
                print(f"  [hata] yigin atlandi: {exc}", file=sys.stderr)
                return None
            time.sleep(pause * (2 ** attempt))
    return None


def first_xref(record, database):
    for ref in record.get("uniProtKBCrossReferences") or []:
        if ref.get("database") == database:
            return ref.get("id")
    return None


def ec_numbers(record):
    """Kayittaki EC numaralari, virgulle. Hem onerilen hem alternatif adlar."""
    found = []
    desc = record.get("proteinDescription") or {}
    blocks = []
    if desc.get("recommendedName"):
        blocks.append(desc["recommendedName"])
    blocks.extend(desc.get("alternativeNames") or [])
    for block in blocks:
        for item in block.get("ecNumbers") or []:
            value = item.get("value")
            if value and value not in found:
                found.append(value)
    return ",".join(found)


def inactive_bits(record):
    reason = record.get("inactiveReason") or {}
    kind = reason.get("inactiveReasonType")
    merged = reason.get("mergeDemergedTo") or []
    return kind, (",".join(merged) if merged else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--out-dir", default="analysis_out")
    ap.add_argument("--limit", type=int, default=0,
                    help="yalnizca ilk N accession (deneme kosusu)")
    ap.add_argument("--pause", type=float, default=0.25,
                    help="yiginlar arasi bekleme, saniye")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    rows = [r[0] for r in con.execute(
        "SELECT accession FROM ro_protein "
        "WHERE accession IS NOT NULL AND accession <> '' ORDER BY accession")]
    if args.limit:
        rows = rows[:args.limit]
    print(f"[sorgulanacak] {len(rows):,} accession, {BATCH}'luk yiginlar")

    con.execute("DROP TABLE IF EXISTS ro_protein_status")
    con.execute("""
        CREATE TABLE ro_protein_status (
            accession       TEXT PRIMARY KEY,
            entry_type      TEXT,   -- ... | Inactive | missing | unmeasured
            inactive_reason TEXT,   -- DELETED | MERGED | DEMERGED
            merged_into     TEXT,   -- birlestirildiyse yeni accession(lar)
            embl_nuc        TEXT,   -- EMBL/ENA nukleotit kaydi
            refseq          TEXT,
            ec              TEXT,
            organism        TEXT
        )""")
    con.execute("CREATE INDEX idx_rps_type ON ro_protein_status(entry_type)")

    seen = 0
    batches = [rows[i:i + BATCH] for i in range(0, len(rows), BATCH)]
    for index, chunk in enumerate(batches, 1):
        results = fetch(chunk)
        if results is None:
            # Olculemedi. "missing" DEGIL -- o, UniProt'un bilmedigi demek.
            con.executemany(
                "INSERT OR REPLACE INTO ro_protein_status(accession, entry_type) "
                "VALUES (?, 'unmeasured')", [(a,) for a in chunk])
            con.commit()
            continue
        got = set()
        payload = []
        for record in results:
            acc = record.get("primaryAccession")
            if not acc:
                continue
            got.add(acc)
            kind, merged = inactive_bits(record)
            payload.append((
                acc, record.get("entryType"), kind, merged,
                first_xref(record, "EMBL"), first_xref(record, "RefSeq"),
                ec_numbers(record),
                (record.get("organism") or {}).get("scientificName"),
            ))
        # Hic donmeyenler: UniProt bu accession'i taniMIyor.
        for acc in chunk:
            if acc not in got:
                payload.append((acc, "missing", None, None, None, None, "", None))
        con.executemany(
            "INSERT OR REPLACE INTO ro_protein_status VALUES (?,?,?,?,?,?,?,?)",
            payload)
        con.commit()
        seen += len(chunk)
        if index % 25 == 0 or index == len(batches):
            print(f"  {index}/{len(batches)} yigin, {seen:,} accession")
        time.sleep(args.pause)

    # --- ozet
    by_type = dict(con.execute(
        "SELECT entry_type, COUNT(*) FROM ro_protein_status GROUP BY 1"))
    by_reason = dict(con.execute(
        "SELECT inactive_reason, COUNT(*) FROM ro_protein_status "
        "WHERE entry_type='Inactive' GROUP BY 1"))
    alive = con.execute(
        "SELECT COUNT(*) FROM ro_protein_status WHERE entry_type NOT IN "
        "('Inactive','missing','unmeasured')").fetchone()[0]
    recoverable = con.execute(
        """SELECT COUNT(*) FROM ro_protein_status s
           JOIN ro_protein p ON p.accession = s.accession
           WHERE s.entry_type NOT IN ('Inactive','missing','unmeasured')
             AND s.embl_nuc IS NOT NULL AND s.embl_nuc <> ''
             AND (p.linked_candidate_id IS NULL OR p.linked_candidate_id = '')
             AND p.passes_current_gate = 1""").fetchone()[0]
    with_ec = con.execute(
        "SELECT COUNT(*) FROM ro_protein_status WHERE ec <> ''").fetchone()[0]

    doc = {
        "source": "UniProtKB REST",
        "n_queried": len(rows),
        "by_entry_type": by_type,
        "inactive_by_reason": by_reason,
        "alive": alive,
        "ec_present": with_ec,
        "recoverable_genomic_context": recoverable,
        "note": ("recoverable_genomic_context counts entries that are still "
                 "live in UniProt, carry an EMBL nucleotide cross-reference, "
                 "pass the motif gate and are not already linked to a genomic "
                 "candidate -- that is, the entries for which a genomic record "
                 "could be fetched and an operon built."),
    }
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, "protein_freshness.json")
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(doc, handle, indent=1)

    print(f"\n[yazildi] {out_path}")
    for key, value in sorted(by_type.items(), key=lambda kv: -kv[1]):
        share = 100.0 * value / max(1, len(rows))
        print(f"  {key or '?':34s} {value:7,d}  ({share:.1f}%)")
    if by_reason:
        print("  Inactive sebebi:",
              ", ".join(f"{k}={v:,}" for k, v in sorted(by_reason.items())))
    print(f"  halen canli                      {alive:7,d}")
    print(f"  EC numarasi olan                 {with_ec:7,d}")
    print(f"  genomik baglami kurulabilir      {recoverable:7,d}")
    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
