# -*- coding: utf-8 -*-
"""`ro_union`: her kaynaktan gelen RO adaylarinin BIRLESIMI, birim DIZI.

NEDEN BU TABLO VAR. Bugun ayni alpha alt birimi uc yerde uc kere sayiliyor
(`ro`, `ro_protein`, `ro_protein_status`) ve sayfadan hangi sayinin ne
oldugu okunamiyor. 19.125 bir KAPSAM degil, bir KESISIM: hem genomik
baglami olan, hem kapiyi gecen girisler. Kullanicinin istedigi sey
(2026-11-01) kesisim degil BIRLESIM: "NCBI'daki, ordaki, burdaki tum RO
olabilecek her sey bu DB'den toparlanabilsin".

BIRIM NEDEN ACCESSION DEGIL, DIZI. UniProt anlik goruntusundeki
accession'larin %65,1'i artik silinmis durumda; ayni dizi bir baska
proteomda bambaska bir accession'la yasiyor olabilir. Accession'i birim
saymak, UniProt'un temizlik kararlarini bizim biyolojik sayimimiz gibi
gostermek olur. Birim dizidir; accession o dizinin bir KAYNAK KIMLIGIDIR
ve `ro_union_source` tablosunda, kaynagi ve durumu yaninda tutulur.

KAPI KARAR VERMEZ, ETIKETLER. `passes_gate=0` olan satirlar tabloda
KALIR -- cikarilmazlar. "RO olabilecek her sey" ancak boyle toparlanir;
kullanici hangi tiere bakacagina kendi karar verir.
"""
import argparse
import hashlib
import os
import sqlite3
import sys

FASTA = "ro_alpha.fasta"


def md5(seq):
    return hashlib.md5(seq.strip().upper().encode()).hexdigest()


def read_fasta(path):
    """{ad: dizi}. Basliktaki ilk alan ad olarak alinir."""
    out, name, chunks = {}, None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name:
                chunks.append(line.strip())
    if name:
        out[name] = "".join(chunks)
    return out


SCHEMA = """
DROP TABLE IF EXISTS ro_union;
DROP TABLE IF EXISTS ro_union_source;
CREATE TABLE ro_union (
    seq_md5             TEXT PRIMARY KEY,
    sequence            TEXT,
    length              INTEGER,
    n_sources           INTEGER,
    has_genomic_context INTEGER,
    passes_gate         INTEGER,
    rieske_intact       INTEGER,
    triad_sites_found   INTEGER,
    ro_cluster          TEXT,
    ro_group            TEXT,
    genomic_candidate_id TEXT,
    organism            TEXT,
    first_source        TEXT
);
CREATE INDEX idx_rou_cluster ON ro_union(ro_cluster);
CREATE INDEX idx_rou_gate    ON ro_union(passes_gate);
CREATE INDEX idx_rou_ctx     ON ro_union(has_genomic_context);
CREATE TABLE ro_union_source (
    seq_md5   TEXT,
    source    TEXT,   -- ncbi_genome | uniprot | refseq_protein
    accession TEXT,
    status    TEXT,   -- live | deleted | merged | demerged | confirmed
    PRIMARY KEY (seq_md5, source, accession)
);
CREATE INDEX idx_rous_md5 ON ro_union_source(seq_md5);
CREATE INDEX idx_rous_src ON ro_union_source(source, status);
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="roar.sqlite")
    ap.add_argument("--fasta", default=FASTA)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    con = sqlite3.connect(args.db)

    # ---- 1. kaynak: NCBI genom kayitlari (anotasyonlu, kapiyi gecmis) ----
    rec = {}          # seq_md5 -> alan sozlugu
    src = []          # (seq_md5, source, accession, status)
    n_gen = 0
    for cid, nuc, pid, cluster, group, seq in con.execute(
            "SELECT candidate_id, nucleotide_id, protein_id, ro_cluster, "
            "ro_group, sequence FROM ro "
            "WHERE is_confirmed=1 AND sequence IS NOT NULL"):
        n_gen += 1
        h = md5(seq)
        e = rec.setdefault(h, {
            "sequence": seq.strip().upper(), "length": len(seq.strip()),
            "has_genomic_context": 0, "passes_gate": 0,
            "rieske_intact": None, "triad_sites_found": None,
            "ro_cluster": None, "ro_group": None,
            "genomic_candidate_id": None, "organism": None,
            "first_source": "ncbi_genome",
        })
        e["has_genomic_context"] = 1
        e["passes_gate"] = 1          # is_confirmed kapiyi gecmis demektir
        e["ro_cluster"] = e["ro_cluster"] or cluster
        e["ro_group"] = e["ro_group"] or group
        e["genomic_candidate_id"] = e["genomic_candidate_id"] or cid
        src.append((h, "ncbi_genome", pid or cid, "confirmed"))
    print("[kaynak 1] NCBI genom: %s giris, %s tekil dizi"
          % (format(n_gen, ","), format(len(rec), ",")), flush=True)

    # Organizma adini replicon'dan tasi (genomik tarafin organizmasi orada).
    org = {}
    for cid, organism in con.execute(
            "SELECT r.candidate_id, p.organism FROM ro r "
            "JOIN replicon p USING(nucleotide_id) WHERE r.is_confirmed=1"):
        org[cid] = organism
    for h, e in rec.items():
        if e["genomic_candidate_id"]:
            e["organism"] = org.get(e["genomic_candidate_id"])

    # ---- 2. kaynak: UniProt anlik goruntusu (canli ve silinmis) ----
    if not os.path.exists(args.fasta):
        print("[hata] %s yok -- protein tarafi atlanamaz, durduruluyor"
              % args.fasta, file=sys.stderr)
        return 2
    seqs = read_fasta(args.fasta)
    print("[kaynak 2] %s dizisi okundu" % format(len(seqs), ","), flush=True)

    status = {a: (t, r) for a, t, r in con.execute(
        "SELECT accession, entry_type, inactive_reason FROM ro_protein_status")}
    org_up = dict(con.execute(
        "SELECT accession, organism FROM ro_protein_status"))

    n_prot = n_new = 0
    for pname, acc, cluster, rieske, triad, passes in con.execute(
            "SELECT protein_name, accession, ro_cluster, rieske_intact, "
            "triad_sites_found, passes_current_gate FROM ro_protein"):
        seq = seqs.get(pname)
        if not seq:
            continue
        n_prot += 1
        h = md5(seq)
        e = rec.get(h)
        if e is None:
            n_new += 1
            e = rec[h] = {
                "sequence": seq.strip().upper(), "length": len(seq.strip()),
                "has_genomic_context": 0, "passes_gate": 0,
                "rieske_intact": None, "triad_sites_found": None,
                "ro_cluster": None, "ro_group": None,
                "genomic_candidate_id": None, "organism": None,
                "first_source": "uniprot",
            }
        # Olcumu yalnizca BOS olani doldur: genomik taraf capali okumayla
        # olculuyor, protein tarafi sabit kolonla. Genomik olcum daha
        # guvenilir oldugu icin onun uzerine yazilmaz.
        if e["rieske_intact"] is None:
            e["rieske_intact"] = rieske
        if e["triad_sites_found"] is None:
            e["triad_sites_found"] = triad
        if not e["ro_cluster"]:
            e["ro_cluster"] = cluster
        e["passes_gate"] = max(e["passes_gate"], passes or 0)
        if not e["organism"]:
            e["organism"] = org_up.get(acc)
        etype, reason = status.get(acc, (None, None))
        if etype == "Inactive":
            st = (reason or "DELETED").lower()
        elif etype in (None, "missing", "unmeasured"):
            st = "unmeasured"
        else:
            st = "live"
        src.append((h, "uniprot", acc, st))
    print("[kaynak 2] UniProt: %s giris, %s YENI tekil dizi"
          % (format(n_prot, ","), format(n_new, ",")), flush=True)

    # ---- kaynak sayisi ----
    per = {}
    for h, source, acc, st in src:
        per.setdefault(h, set()).add((source, acc))
    for h, e in rec.items():
        e["n_sources"] = len(per.get(h, ()))

    print("\n[birlesim] %s tekil dizi" % format(len(rec), ","))
    ctx = sum(1 for e in rec.values() if e["has_genomic_context"])
    gate = sum(1 for e in rec.values() if e["passes_gate"])
    print("   genomik baglami olan        : %s" % format(ctx, ","))
    print("   genomik baglami OLMAYAN     : %s" % format(len(rec) - ctx, ","))
    print("   kapiyi gecen                : %s" % format(gate, ","))
    print("   kapiyi GECMEYEN (aday kalir): %s" % format(len(rec) - gate, ","))

    if args.dry_run:
        print("[dry-run] tablo yazilmadi")
        return 0

    con.executescript(SCHEMA)
    con.executemany(
        "INSERT INTO ro_union (seq_md5, sequence, length, n_sources, "
        "has_genomic_context, passes_gate, rieske_intact, triad_sites_found, "
        "ro_cluster, ro_group, genomic_candidate_id, organism, first_source) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
        [(h, e["sequence"], e["length"], e["n_sources"],
          e["has_genomic_context"], e["passes_gate"], e["rieske_intact"],
          e["triad_sites_found"], e["ro_cluster"], e["ro_group"],
          e["genomic_candidate_id"], e["organism"], e["first_source"])
         for h, e in rec.items()])
    con.executemany(
        "INSERT OR IGNORE INTO ro_union_source "
        "(seq_md5, source, accession, status) VALUES (?,?,?,?)", src)
    con.commit()
    print("[yazildi] ro_union %s satir, ro_union_source %s satir"
          % (format(len(rec), ","), format(len(src), ",")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
