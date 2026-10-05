"""
Operon bulucu -- dizi tabanli bilesen tespiti + ayni-iplik/kisa-aralik kurali.

NEDEN: README_PIPELINE'daki bilinen sinir "operon tamligi olculemiyor" idi;
beta/ferredoksin/reduktaz komsulari sik sik 'hypothetical protein' diye
anotasyonlu oldugu icin urun adi regex'i biyolojiyi degil anotasyon
kelime dagarcigini olcuyordu. Bu adim komsu PROTEIN DIZILERINI gbk'dan
cekip Pfam HMM'leriyle test eder; anotasyon metnine bagimli degildir.

Adimlar:
    1. ro.sequence sutununu doldur (genomic_context/ro_candidates.fasta)
    2. Dogrulanmis RO iceren replikonlarin gbk'larini yeniden parse et,
       neighbor tablosundaki CDS'lerin translation'larini al -> neighbor_protein
    3. hmmsearch component_hmms/components.hmm -> neighbor_domain
    4. Her komsuya bilesen etiketi: beta | ferredoxin | reductase | alpha_other
    5. Her dogrulanmis RO icin operon: RO'dan iki yone, ayni iplikte,
       genler arasi bosluk <= --max-gap (varsayilan 150 bp) oldukca yuru.
       -> operon, operon_gene tablolari

Kullanim:
    python3 build_operons.py --db roar.sqlite --gbk-dir gbk_files --cpu 20

Tekrar calistirmak guvenli: tablolar bastan uretilir, uzun adimlar
(gbk parse, hmmsearch) ciktilari varsa atlanir (--force ile zorlanir).
"""

import argparse
import csv
import gzip
import os
import sqlite3
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool

from Bio import SeqIO

csv.field_size_limit(10 ** 7)

COMPONENT_HMM = os.path.join("component_hmms", "components.hmm")

# Pfam modeli -> hangi bilesenin kaniti. Siralama onemli degil; karar
# assign_component() icinde uzunluk + kombinasyonla verilir.
BETA_MODELS = {"Ring_hydroxyl_B"}
FERREDOXIN_MODELS = {"Fer2", "Rieske", "Rieske_2"}
REDUCTASE_MODELS = {"NAD_binding_1", "FAD_binding_6", "FAD_binding_8",
                    "Pyr_redox_2", "Pyr_redox_dim"}
RIESKE_MODELS = {"Rieske", "Rieske_2"}

# Rieske tipi ferredoksin ~100-130 aa, bitki tipi [2Fe-2S] ferredoksin ~100 aa.
# RO alpha >= 300 aa. Arada kalan Rieske'li proteinler (ISP vb.) 'rieske_other'.
MAX_FERREDOXIN_LEN = 250
MIN_ALPHA_LEN = 300

DEFAULT_MAX_GAP = 150        # bp; prokaryot operon tahmininde yaygin esik
DEFAULT_EVALUE = 1e-5

SCHEMA = """
CREATE TABLE IF NOT EXISTS neighbor_protein (
    protein_key   TEXT PRIMARY KEY,   -- nucleotide_id:start-end:strand
    nucleotide_id TEXT,
    start         INTEGER,
    end           INTEGER,
    strand        INTEGER,
    length        INTEGER,
    translation   TEXT
);
CREATE TABLE IF NOT EXISTS neighbor_domain (
    protein_key TEXT REFERENCES neighbor_protein(protein_key),
    hmm         TEXT,
    score       REAL,
    evalue      REAL,
    coverage    REAL
);
CREATE TABLE IF NOT EXISTS neighbor_component (
    protein_key TEXT PRIMARY KEY,
    component   TEXT,    -- beta | ferredoxin | reductase | alpha_other | rieske_other | none
    evidence    TEXT     -- hangi modeller
);
CREATE TABLE IF NOT EXISTS operon (
    candidate_id   TEXT PRIMARY KEY REFERENCES ro(candidate_id),
    nucleotide_id  TEXT,
    strand         INTEGER,
    start          INTEGER,
    end            INTEGER,
    n_genes        INTEGER,
    has_beta       INTEGER,
    has_ferredoxin INTEGER,
    has_reductase  INTEGER,
    completeness   INTEGER,   -- 0..3: operon icindeki beta+ferredoksin+reduktaz
    nearby_beta       INTEGER,  -- +-10 kb icinde (operon disinda olsa da)
    nearby_ferredoxin INTEGER,
    nearby_reductase  INTEGER,
    layout         TEXT       -- 5'->3' bilesen dizisi, RO '[alpha]' ile
);
CREATE TABLE IF NOT EXISTS operon_gene (
    candidate_id TEXT,
    position     INTEGER,     -- 0 = RO, negatif = yukari (5'), pozitif = asagi (3')
    neighbor_id  INTEGER,     -- NULL ise RO'nun kendisi
    protein_key  TEXT,
    component    TEXT,
    category     TEXT,
    product      TEXT,
    gap_bp       INTEGER      -- bir onceki gene olan bosluk (negatif = ortusme)
);
CREATE INDEX IF NOT EXISTS idx_og_cand ON operon_gene(candidate_id);
CREATE INDEX IF NOT EXISTS idx_ro_cluster ON ro(ro_cluster, is_confirmed);
CREATE INDEX IF NOT EXISTS idx_nb_nuc ON neighbor(nucleotide_id, start, end, strand);
CREATE INDEX IF NOT EXISTS idx_roleaf_leaf ON ro_leaf(leaf_id);
CREATE INDEX IF NOT EXISTS idx_nd_key  ON neighbor_domain(protein_key);
"""


def protein_key(nucleotide_id, start, end, strand):
    return f"{nucleotide_id}:{start}-{end}:{strand}"


# --------------------------------------------------------------- 1. ro.sequence
def fill_ro_sequences(connection, fasta):
    cols = [r[1] for r in connection.execute("PRAGMA table_info(ro)")]
    if "sequence" not in cols:
        connection.execute("ALTER TABLE ro ADD COLUMN sequence TEXT")
    if connection.execute("SELECT COUNT(*) FROM ro WHERE sequence IS NULL").fetchone()[0] == 0:
        print("[atlandi] ro.sequence dolu")
        return
    if not os.path.exists(fasta):
        print(f"[uyari] {fasta} yok, ro.sequence doldurulamadi")
        return
    rows, name, chunks = [], None, []
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    rows.append(("".join(chunks), name))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name:
        rows.append(("".join(chunks), name))
    connection.executemany("UPDATE ro SET sequence=? WHERE candidate_id=?", rows)
    connection.commit()
    print(f"[ro.sequence] {len(rows)} dizi yazildi")


# ------------------------------------------------------- 2. neighbor_protein
def _extract_worker(job):
    path, wanted = job      # wanted: set of (start, end, strand)
    out = []
    try:
        handle = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
        with handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    parts = feature.location.parts
                    start = min(int(p.start) for p in parts)
                    end = max(int(p.end) for p in parts)
                    strand = feature.location.strand or 0
                    if (start, end, strand) not in wanted:
                        continue
                    translation = feature.qualifiers.get("translation", [""])[0]
                    if translation:
                        out.append((protein_key(record.id, start, end, strand),
                                    record.id, start, end, strand,
                                    len(translation), translation))
    except Exception as exc:        # tek dosya tum isi durdurmasin
        sys.stderr.write(f"[hata] {path}: {exc}\n")
    return out


def extract_neighbor_proteins(connection, gbk_dir, processes, force):
    existing = connection.execute("SELECT COUNT(*) FROM neighbor_protein").fetchone()[0]
    if existing and not force:
        print(f"[atlandi] neighbor_protein dolu ({existing})")
        return
    connection.execute("DELETE FROM neighbor_protein")

    wanted = defaultdict(set)
    files = dict(connection.execute("SELECT nucleotide_id, file FROM replicon"))
    for nuc, start, end, strand in connection.execute("""
            SELECT DISTINCT n.nucleotide_id, n.start, n.end, n.strand
            FROM neighbor n JOIN ro r ON r.candidate_id = n.candidate_id
            WHERE r.is_confirmed = 1"""):
        wanted[nuc].add((start, end, strand))

    jobs = []
    for nuc, keys in wanted.items():
        path = os.path.join(gbk_dir, files.get(nuc, f"{nuc}.gbk"))
        if os.path.exists(path):
            jobs.append((path, keys))
        else:
            sys.stderr.write(f"[uyari] gbk yok: {path}\n")
    print(f"[neighbor_protein] {len(jobs)} gbk dosyasi, {processes} surec")

    total = 0
    with Pool(processes) as pool:
        for i, rows in enumerate(pool.imap_unordered(_extract_worker, jobs, chunksize=8), 1):
            connection.executemany(
                "INSERT OR REPLACE INTO neighbor_protein VALUES (?,?,?,?,?,?,?)", rows)
            total += len(rows)
            if i % 500 == 0:
                connection.commit()
                sys.stdout.write(f"\r  {i}/{len(jobs)} dosya | {total} protein")
                sys.stdout.flush()
    connection.commit()
    print(f"\r[neighbor_protein] {total} protein yazildi" + " " * 20)


# -------------------------------------------------------- 3. neighbor_domain
def run_component_hmms(connection, hmm_file, work_dir, cpu, force):
    existing = connection.execute("SELECT COUNT(*) FROM neighbor_domain").fetchone()[0]
    if existing and not force:
        print(f"[atlandi] neighbor_domain dolu ({existing})")
        return
    os.makedirs(work_dir, exist_ok=True)
    fasta = os.path.join(work_dir, "neighbor_proteins.fasta")
    domtbl = os.path.join(work_dir, "neighbor_components.domtbl")

    with open(fasta, "w") as handle:
        for key, seq in connection.execute(
                "SELECT protein_key, translation FROM neighbor_protein"):
            handle.write(f">{key}\n{seq}\n")

    if force or not os.path.exists(domtbl) or os.path.getsize(domtbl) == 0:
        cmd = ["hmmsearch", "--domtblout", domtbl, "--noali", "--cpu", str(cpu),
               "-E", str(DEFAULT_EVALUE), hmm_file, fasta]
        print("[calisiyor]", " ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)

    # ro_filter.parse_domtblout ayni formati okur; coverage = HMM ekseninde kaplama
    from ro_filter import parse_domtblout
    hits = parse_domtblout(domtbl)
    rows = [(key, h["query_name"], h["score"], h["e_value"], round(h["model_coverage"], 3))
            for key, hs in hits.items() for h in hs]
    connection.execute("DELETE FROM neighbor_domain")
    connection.executemany("INSERT INTO neighbor_domain VALUES (?,?,?,?,?)", rows)
    connection.commit()
    print(f"[neighbor_domain] {len(rows)} domain hiti ({len(hits)} protein)")


# ----------------------------------------------------- 4. neighbor_component
def assign_component(models, length):
    """Pfam hit kumesi + uzunluk -> bilesen etiketi."""
    if models & BETA_MODELS:
        return "beta"
    if models & REDUCTASE_MODELS:
        return "reductase"
    if models & RIESKE_MODELS and length >= MIN_ALPHA_LEN:
        return "alpha_other"
    if models & FERREDOXIN_MODELS and length <= MAX_FERREDOXIN_LEN:
        return "ferredoxin"
    if models & RIESKE_MODELS:
        return "rieske_other"
    return "none"


def assign_components(connection):
    models_by_key = defaultdict(set)
    for key, hmm in connection.execute("SELECT protein_key, hmm FROM neighbor_domain"):
        models_by_key[key].add(hmm)
    rows = []
    for key, length in connection.execute("SELECT protein_key, length FROM neighbor_protein"):
        models = models_by_key.get(key, set())
        rows.append((key, assign_component(models, length), ";".join(sorted(models))))
    connection.execute("DELETE FROM neighbor_component")
    connection.executemany("INSERT INTO neighbor_component VALUES (?,?,?)", rows)
    connection.commit()
    counts = defaultdict(int)
    for _, comp, _ in rows:
        counts[comp] += 1
    print("[neighbor_component]", dict(counts))


# ----------------------------------------------------------------- 5. operon
def walk(ro, neighbors, max_gap):
    """RO'dan iki yone ayni iplikte yuru. Donen: [(neighbor, gap)] siralı liste."""
    left = sorted((n for n in neighbors if n["gene_offset"] < 0),
                  key=lambda n: -n["gene_offset"])     # -1, -2, ...
    right = sorted((n for n in neighbors if n["gene_offset"] > 0),
                   key=lambda n: n["gene_offset"])     # +1, +2, ...

    def chain(side, direction):
        prev_start, prev_end = ro["start"], ro["end"]
        out = []
        for n in side:
            if n["strand"] != ro["strand"]:
                break
            if n["spans_origin"]:
                gap = 0                             # koordinat guvenilmez, mesafe kucuk
            elif direction < 0:
                gap = prev_start - n["end"]
            else:
                gap = n["start"] - prev_end
            if gap > max_gap:
                break
            out.append((n, gap))
            prev_start, prev_end = n["start"], n["end"]
        return out

    return chain(left, -1), chain(right, +1)


def build_operons(connection, max_gap):
    connection.execute("DELETE FROM operon")
    connection.execute("DELETE FROM operon_gene")

    component = dict(connection.execute(
        "SELECT protein_key, component FROM neighbor_component"))
    category = {}
    for nid, cat in connection.execute(
            "SELECT neighbor_id, category FROM gene_category WHERE method='regex_v1'"):
        category[nid] = cat

    ros = [dict(zip(("candidate_id", "nucleotide_id", "start", "end", "strand", "product"), r))
           for r in connection.execute(
               "SELECT candidate_id, nucleotide_id, start, end, strand, product "
               "FROM ro WHERE is_confirmed=1")]
    neighbors_by_ro = defaultdict(list)
    cols = ("neighbor_id", "candidate_id", "nucleotide_id", "start", "end", "strand",
            "distance", "same_strand", "gene_offset", "spans_origin", "product")
    for row in connection.execute(
            "SELECT neighbor_id, candidate_id, nucleotide_id, start, end, strand, distance, "
            "same_strand, gene_offset, spans_origin, product FROM neighbor"):
        neighbors_by_ro[row[1]].append(dict(zip(cols, row)))

    op_rows, gene_rows = [], []
    for ro in ros:
        nbs = neighbors_by_ro.get(ro["candidate_id"], [])
        left, right = walk(ro, nbs, max_gap)

        def comp_of(n):
            return component.get(protein_key(n["nucleotide_id"], n["start"], n["end"],
                                             n["strand"]), "none")

        # 5'->3' sira: minus iplikte okuma yonu ters
        ordered = [(n, g, -i - 1) for i, (n, g) in enumerate(left)][::-1] \
            + [(None, 0, 0)] + [(n, g, i + 1) for i, (n, g) in enumerate(right)]
        if ro["strand"] == -1:
            ordered = [(n, g, -p) for n, g, p in ordered[::-1]]

        comps = []
        for n, gap, pos in ordered:
            if n is None:
                comps.append("[alpha]")
                gene_rows.append((ro["candidate_id"], 0, None, None, "alpha", "ro_alpha",
                                  ro["product"], 0))
            else:
                comp = comp_of(n)
                comps.append(comp if comp != "none" else category.get(n["neighbor_id"], "other"))
                gene_rows.append((ro["candidate_id"], pos, n["neighbor_id"],
                                  protein_key(n["nucleotide_id"], n["start"], n["end"], n["strand"]),
                                  comp, category.get(n["neighbor_id"], "other"),
                                  n["product"], gap))

        in_op = {comp_of(n) for n, _ in left + right}
        nearby = {comp_of(n) for n in nbs}
        genes = [n for n, _ in left + right]
        starts = [ro["start"]] + [n["start"] for n in genes]
        ends = [ro["end"]] + [n["end"] for n in genes]
        has = {c: int(c in in_op) for c in ("beta", "ferredoxin", "reductase")}
        op_rows.append((
            ro["candidate_id"], ro["nucleotide_id"], ro["strand"], min(starts), max(ends),
            len(genes) + 1, has["beta"], has["ferredoxin"], has["reductase"],
            sum(has.values()),
            int("beta" in nearby), int("ferredoxin" in nearby), int("reductase" in nearby),
            " > ".join(comps)))

    connection.executemany("INSERT INTO operon VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", op_rows)
    connection.executemany("INSERT INTO operon_gene VALUES (?,?,?,?,?,?,?,?)", gene_rows)
    connection.commit()

    print(f"[operon] {len(op_rows)} RO icin operon turetildi (max_gap={max_gap} bp)")
    for label, q in [
            ("operonda beta", "has_beta=1"), ("operonda ferredoksin", "has_ferredoxin=1"),
            ("operonda reduktaz", "has_reductase=1"), ("tam (3/3)", "completeness=3"),
            ("hic bilesen yok", "completeness=0"),
            ("+-10kb icinde beta", "nearby_beta=1")]:
        n = connection.execute(f"SELECT COUNT(*) FROM operon WHERE {q}").fetchone()[0]
        print(f"   {label:24s}: {n:6d}  ({100*n/max(1,len(op_rows)):.1f}%)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default="roar.sqlite")
    parser.add_argument("--gbk-dir", default="gbk_files")
    parser.add_argument("--context-dir", default="genomic_context")
    parser.add_argument("--hmm", default=COMPONENT_HMM)
    parser.add_argument("--work-dir", default="operon_work")
    parser.add_argument("--cpu", type=int, default=8)
    parser.add_argument("--max-gap", type=int, default=DEFAULT_MAX_GAP)
    parser.add_argument("--force", action="store_true",
                        help="gbk parse ve hmmsearch'u ciktilar olsa da yeniden yap")
    parser.add_argument("--operons-only", action="store_true",
                        help="sadece operon turetmeyi tekrarla (ornegin --max-gap degisince)")
    args = parser.parse_args()

    connection = sqlite3.connect(args.db)
    connection.executescript(SCHEMA)

    if not args.operons_only:
        fill_ro_sequences(connection, os.path.join(args.context_dir, "ro_candidates.fasta"))
        extract_neighbor_proteins(connection, args.gbk_dir, args.cpu, args.force)
        run_component_hmms(connection, args.hmm, args.work_dir, args.cpu, args.force)
        assign_components(connection)
    build_operons(connection, args.max_gap)
    connection.close()


if __name__ == "__main__":
    main()
