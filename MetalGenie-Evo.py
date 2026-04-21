#!/usr/bin/env python3
"""
MetalGenie-Evo  –  Modernised FeGenie-compatible HMM annotation pipeline
======================================================================
Improvements over the original FeGenie and the first reimplementation:

  1. GFF-based clustering  – uses actual bp coordinates + strand from Prodigal
     GFF output instead of fragile ORF-index heuristics. Falls back to index-
     based clustering when no GFF is available (--gff_dir not given).

  2. Strand-aware operons  – optional flag (--strand_aware): genes on opposite
     strands are split into separate clusters before operon-rule filtering.

  3. Parallel hmmsearch    – ProcessPoolExecutor dispatches one job per
     (genome, HMM) pair; --threads controls the pool size.

  4. Configurable operon rules  – loaded from operon_rules.json in the HMM
     library directory. Hardcoded FeGenie defaults are used as fallback so the
     tool works with an unmodified FeGenie iron/ directory. MetHMMDB metal-
     resistance categories are automatically assigned report_all semantics.

  5. Tblout result caching  – already-computed .tblout files are reused on
     subsequent runs (useful during incremental analysis / debugging).

  6. Normalised heatmap  – optional --norm flag outputs gene counts normalised
     by total predicted ORFs per genome (same as FeGenie --norm).

Input
-----
  --faa_dir      Directory of ORF FASTA files (.faa) from Prodigal
  --faa_ext      Extension of ORF files (default: faa)
  --gff_dir      Directory of Prodigal GFF files (same basename as .faa);
                 enables coordinate-based clustering. Optional.
  --hmm_dir      HMM library directory (FeGenie iron/ structure or output of
                 build_hmm_library.py)
  --out          Output directory

HMM library structure expected
-------------------------------
  hmm_dir/
    <category_1>/
        gene_A.hmm
        gene_B.hmm
    <category_2>/
        ...
    HMM-bitcutoffs.txt      HMM_stem<tab>bitscore
    FeGenie-map.txt         HMM_stem<tab>readable_gene_name   (optional)
    operon_rules.json       Configurable operon rules          (optional)

Outputs (in --out/)
-------------------
  MetalGenie-Evo-summary.csv
  MetalGenie-Evo-geneSummary-clusters.csv   (FeGenie-R-script compatible)
  MetalGenie-Evo-heatmap-data.csv
"""

import argparse
import json
import os
import re
import sys
import fnmatch
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


# ═══════════════════════════════════════════════════════════════════════════════
# CONSTANTS  –  FeGenie hardcoded operon rules (used when no JSON is present)
# ═══════════════════════════════════════════════════════════════════════════════

_DEFAULT_OPERON_RULES = [
    # ── Electron shuttle (FLEET) ─────────────────────────────────────────────
    {
        "name": "FLEET",
        "categories": ["iron_oxidation"],
        "genes": ["EetA", "EetB", "Ndh2", "FmnB", "FmnA", "DmkA", "DmkB", "PplA"],
        "rule": "require_n_of",
        "min_genes": 5,
        "on_fail": "passthrough_non_members",
    },
    # ── Magnetosome ───────────────────────────────────────────────────────────
    {
        "name": "MAM",
        "categories": ["magnetosome_formation"],
        "genes": ["MamA", "MamB", "MamE", "MamK", "MamP", "MamM",
                  "MamQ", "MamI", "MamL", "MamO"],
        "rule": "require_n_of",
        "min_genes": 5,
        "on_fail": "passthrough_non_members",
    },
    # ── FoxABC ────────────────────────────────────────────────────────────────
    {
        "name": "FOXABC",
        "categories": ["iron_oxidation"],
        "genes": ["FoxA", "FoxB", "FoxC"],
        "rule": "require_n_of",
        "min_genes": 2,
        "on_fail": "passthrough_non_members",
    },
    # ── FoxEYZ  (FoxE mandatory) ──────────────────────────────────────────────
    {
        "name": "FOXEYZ",
        "categories": ["iron_oxidation"],
        "genes": ["FoxE", "FoxY", "FoxZ"],
        "rule": "require_anchor",
        "anchor": "FoxE",
        "on_fail": "passthrough_non_members",
    },
    # ── DFE operons ───────────────────────────────────────────────────────────
    {
        "name": "DFE1",
        "categories": ["iron_reduction", "probable_iron_reduction"],
        "genes": ["DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"],
        "rule": "require_n_of",
        "min_genes": 3,
        "on_fail": "passthrough_non_members",
    },
    {
        "name": "DFE2",
        "categories": ["iron_reduction", "probable_iron_reduction"],
        "genes": ["DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"],
        "rule": "require_n_of",
        "min_genes": 3,
        "on_fail": "passthrough_non_members",
    },
    # ── Mtr/Mto disambiguation ────────────────────────────────────────────────
    {
        "name": "MtrMto",
        "categories": ["iron_oxidation", "iron_reduction",
                       "possible_iron_oxidation_and_possible_iron_reduction"],
        "genes": ["MtrA", "MtrB_TIGR03509", "MtrC_TIGR03507", "MtoA", "CymA"],
        "rule": "mtr_disambiguation",
        "on_fail": "keep_all",
    },
    # ── Siderophore transport (need ≥2 or lone trusted gene) ─────────────────
    {
        "name": "SIDERO_TRANSPORT",
        "categories": [
            "iron_aquisition-siderophore_transport_potential",
            "iron_aquisition-heme_transport",
            "iron_aquisition-siderophore_transport",
        ],
        "genes": [],    # all genes in these categories qualify
        "rule": "require_n_cat_or_lone_trusted",
        "min_genes": 2,
        "trusted_lone": [
            "FutA1-iron_ABC_transporter_iron-binding-rep",
            "FutA2-iron_ABC_transporter_iron-binding-rep",
            "FutC-iron_ABC_transporter_ATPase-rep",
            "LbtU-LvtA-PiuA-PirA-RhtA",
            "LbtU-LbtB-legiobactin_receptor",
            "LbtU_LbtB-legiobactin_receptor_2",
            "IroC-salmochelin_transport-rep",
        ],
        "on_fail": "drop",
    },
    # ── Siderophore synthesis (≥3 unique synth genes) ─────────────────────────
    {
        "name": "SIDERO_SYNTH",
        "categories": ["iron_aquisition-siderophore_synthesis"],
        "genes": [],
        "rule": "require_n_cat",
        "min_genes": 3,
        "on_fail": "drop",
    },
    # ── Iron transport / heme oxygenase (≥2 genes) ────────────────────────────
    {
        "name": "IRON_TRANSPORT",
        "categories": ["iron_aquisition-iron_transport",
                       "iron_aquisition-heme_oxygenase"],
        "genes": [],
        "rule": "require_n_cat",
        "min_genes": 2,
        "on_fail": "drop",
    },
]

# Categories matching these patterns skip all operon filtering (report_all)
_REPORT_ALL_PATTERNS = [
    "metal_resistance-*",
    "iron_storage",
]


# ═══════════════════════════════════════════════════════════════════════════════
# IO UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

def read_fasta(path):
    """Return dict: orf_id → sequence (no spaces in key)."""
    seqs, header, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(parts)
                header, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if header is not None:
        seqs[header] = "".join(parts)
    return seqs


def read_cutoffs(path):
    """Return dict: hmm_stem → float bitscore."""
    cutoffs = {}
    if not os.path.isfile(path):
        return cutoffs
    with open(path) as fh:
        for line in fh:
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                try:
                    cutoffs[ls[0]] = float(ls[1])
                except ValueError:
                    pass
    return cutoffs


def read_map(path):
    """Return dict: hmm_stem → readable gene name."""
    gmap = {}
    if not os.path.isfile(path):
        return gmap
    with open(path) as fh:
        for line in fh:
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                gmap[ls[0]] = ls[1]
    return gmap


def load_operon_rules(hmm_dir):
    """
    Load operon_rules.json from the HMM library directory.
    Falls back to hardcoded FeGenie defaults if the file is absent.
    """
    rules_path = Path(hmm_dir) / "operon_rules.json"
    if rules_path.exists():
        with open(rules_path) as fh:
            data = json.load(fh)
        rules = data.get("rules", [])
        report_all_pats = data.get("report_all_categories", _REPORT_ALL_PATTERNS)
        print(f"[INFO] Loaded {len(rules)} operon rules from {rules_path}")
    else:
        rules = _DEFAULT_OPERON_RULES
        report_all_pats = _REPORT_ALL_PATTERNS
        print("[INFO] Using built-in FeGenie operon rules (no operon_rules.json found)")
    return rules, report_all_pats


# ═══════════════════════════════════════════════════════════════════════════════
# GFF PARSING  (Prodigal output)
# ═══════════════════════════════════════════════════════════════════════════════

def load_prodigal_gff(gff_path):
    """
    Parse a Prodigal GFF file.
    Returns dict: orf_id → {"contig": str, "start": int, "end": int, "strand": str}

    Prodigal GFF example line:
      NODE_1  Prodigal_v2.6.3  CDS  1  1533  85.7  +  0  ID=NODE_1_1;...
    ORF ID in FASTA header = seqname + "_" + ordinal  (= NODE_1_1 here)
    """
    coords = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]
            # Extract ID from attributes  (ID=<orf_id>;...)
            m = re.search(r"ID=([^;]+)", parts[8])
            if m:
                orf_id = m.group(1).strip()
                coords[orf_id] = {
                    "contig": contig,
                    "start":  start,
                    "end":    end,
                    "strand": strand,
                }
    return coords


def load_gff_dir(gff_dir, faa_files):
    """
    Load GFF files matching the FAA basenames.
    Returns dict: genome_name → orf_coords_dict (from load_prodigal_gff)
    Accepted GFF extensions: .gff, .gff3, .prodigal.gff
    """
    gff_dir = Path(gff_dir)
    genome_coords = {}
    for faa in faa_files:
        stem = Path(faa).stem
        found = None
        for ext in (".gff", ".gff3", ".prodigal.gff"):
            candidate = gff_dir / (stem + ext)
            if candidate.exists():
                found = candidate
                break
        if found:
            genome_coords[Path(faa).name] = load_prodigal_gff(str(found))
        else:
            print(f"  [WARN] No GFF found for {faa.name}, will use index-based clustering",
                  file=sys.stderr)
    return genome_coords


# ═══════════════════════════════════════════════════════════════════════════════
# CLUSTERING
# ═══════════════════════════════════════════════════════════════════════════════

def _index_from_name(orf_name):
    """Extract (contig, int_index) from a Prodigal ORF name like 'scaffold_001_5'."""
    parts = orf_name.rsplit("_", 1)
    if len(parts) == 2:
        try:
            return parts[0], int(parts[1])
        except ValueError:
            pass
    return orf_name, 0


def cluster_by_index(orf_set, max_gap=5):
    """
    Group ORFs into proximity clusters using Prodigal ordinal index.
    orf_set: iterable of orf_id strings
    Returns list of lists (each inner list = one cluster of orf_ids).
    """
    # Group by contig
    by_contig = defaultdict(list)
    for orf in orf_set:
        contig, idx = _index_from_name(orf)
        by_contig[contig].append((idx, orf))

    clusters = []
    for contig, entries in by_contig.items():
        entries.sort(key=lambda x: x[0])
        group = [entries[0][1]]
        for i in range(1, len(entries)):
            gap = entries[i][0] - entries[i - 1][0]
            if gap <= max_gap:
                group.append(entries[i][1])
            else:
                clusters.append(group)
                group = [entries[i][1]]
        clusters.append(group)
    return clusters


def cluster_by_coordinates(orf_set, orf_coords, max_bp_gap=5000,
                            strand_aware=False):
    """
    Group ORFs into proximity clusters using actual genomic bp coordinates.
    orf_coords: dict from load_prodigal_gff
    Returns list of lists of orf_ids.
    """
    by_contig = defaultdict(list)
    for orf in orf_set:
        c = orf_coords.get(orf)
        if c is None:
            # Fallback to index for this ORF
            contig, idx = _index_from_name(orf)
            by_contig[contig].append((idx * 300, idx * 300, "+", orf))
        else:
            by_contig[c["contig"]].append(
                (c["start"], c["end"], c["strand"], orf)
            )

    clusters = []
    for contig, entries in by_contig.items():
        entries.sort(key=lambda x: x[0])

        # Optionally split by strand before clustering
        if strand_aware:
            by_strand = defaultdict(list)
            for e in entries:
                by_strand[e[2]].append(e)
            strand_groups = list(by_strand.values())
        else:
            strand_groups = [entries]

        for sg in strand_groups:
            if not sg:
                continue
            group = [sg[0][3]]
            for i in range(1, len(sg)):
                # Gap = start of current - end of previous
                gap = sg[i][0] - sg[i - 1][1]
                if gap <= max_bp_gap:
                    group.append(sg[i][3])
                else:
                    clusters.append(group)
                    group = [sg[i][3]]
            clusters.append(group)

    return clusters


def build_clusters(genome, orf_hits, orf_coords, max_gap, max_bp_gap,
                   strand_aware):
    """
    Dispatch to coordinate- or index-based clustering depending on availability
    of GFF data for this genome.
    Returns list of lists of orf_ids.
    """
    if orf_coords:
        return cluster_by_coordinates(
            orf_hits.keys(), orf_coords,
            max_bp_gap=max_bp_gap, strand_aware=strand_aware
        )
    else:
        return cluster_by_index(orf_hits.keys(), max_gap=max_gap)


# ═══════════════════════════════════════════════════════════════════════════════
# HMMER
# ═══════════════════════════════════════════════════════════════════════════════

def _hmmsearch_job(args_tuple):
    """Top-level function for ProcessPoolExecutor (must be picklable)."""
    hmm_file, faa_file, tblout_path, bitscore_cutoff, threads = args_tuple
    if Path(tblout_path).exists():
        return tblout_path, True, ""   # cached
    cmd = [
        "hmmsearch",
        "--cpu",    str(threads),
        "-T",       str(max(bitscore_cutoff, 0)),
        "--tblout", tblout_path,
        "--noali",
        "-o",       "/dev/null",
        hmm_file,
        faa_file,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    ok = result.returncode == 0
    return tblout_path, ok, result.stderr if not ok else ""


def parse_tblout(tblout_path):
    """Parse hmmsearch --tblout; return list of (orf, evalue, bitscore)."""
    hits = []
    if not os.path.isfile(tblout_path):
        return hits
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                evalue   = float(parts[4])
                bitscore = float(parts[5])
            except ValueError:
                continue
            if evalue < 0.1:
                hits.append((parts[0], evalue, bitscore))
    return hits


def run_all_hmmsearches(faa_files, cat_hmms, cutoffs, out_tmp,
                        threads_total, hmm_threads=1):
    """
    Launch hmmsearch for every (genome, HMM) pair using a process pool.

    threads_total  – total CPU budget
    hmm_threads    – CPUs given to each individual hmmsearch call (default 1)
    Pool workers   – threads_total // hmm_threads
    """
    jobs = []
    for faa in faa_files:
        for cat, hmm_list in cat_hmms.items():
            for stem, hmm_path in hmm_list:
                tblout = out_tmp / f"{faa.name}__{stem}.tblout"
                cutoff = cutoffs.get(stem, 0)
                jobs.append((str(hmm_path), str(faa), str(tblout),
                             cutoff, hmm_threads))

    n_workers = max(1, threads_total // hmm_threads)
    total     = len(jobs)
    done      = 0
    errors    = 0

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_hmmsearch_job, j): j for j in jobs}
        for fut in as_completed(futures):
            done += 1
            tblout, ok, err = fut.result()
            if not ok:
                errors += 1
                print(f"\n  [WARN] hmmsearch failed for {Path(tblout).name}: {err}",
                      file=sys.stderr)
            sys.stdout.write(f"\r[INFO] hmmsearch: {done}/{total}  "
                             f"({errors} errors)  ")
            sys.stdout.flush()
    print()
    return


# ═══════════════════════════════════════════════════════════════════════════════
# BEST-HIT COLLECTION  (after all hmmsearches are done)
# ═══════════════════════════════════════════════════════════════════════════════

def collect_best_hits(faa_files, cat_hmms, out_tmp):
    """
    Scan all tblout files and keep only the best-scoring hit per ORF per genome.
    Returns dict: genome → {orf → {hmm_stem, cat, evalue, bitscore, cutoff}}
    """
    best_hit = defaultdict(dict)
    for faa in faa_files:
        genome = faa.name
        for cat, hmm_list in cat_hmms.items():
            for stem, _ in hmm_list:
                tblout = out_tmp / f"{genome}__{stem}.tblout"
                for orf, evalue, bitscore in parse_tblout(str(tblout)):
                    prev = best_hit[genome].get(orf)
                    if prev is None or bitscore > prev["bitscore"]:
                        best_hit[genome][orf] = {
                            "hmm_stem": stem,
                            "cat":      cat,
                            "evalue":   evalue,
                            "bitscore": bitscore,
                        }
    return best_hit


# ═══════════════════════════════════════════════════════════════════════════════
# OPERON FILTERING ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

def _category_matches(cat, pattern_list):
    return any(fnmatch.fnmatch(cat, p) for p in pattern_list)


def _unique_genes_in(rows, gene_set):
    return len({r["hmm_stem"] for r in rows if r["hmm_stem"] in gene_set})


def _rows_in_cats(rows, cat_set):
    return [r for r in rows if r["cat"] in cat_set]


def _mtr_disambiguation(rows):
    """
    Mtr/Mto operon context rules (ported from FeGenie):
      MtoA + MtrB (no MtrC) → iron_oxidation
      MtrA + MtrB            → iron_reduction
      MtrC alone             → iron_reduction context
    Returns modified copy of rows.
    """
    stems = {r["hmm_stem"] for r in rows}
    updated = [dict(r) for r in rows]    # shallow copy

    if "MtoA" in stems and "MtrB_TIGR03509" in stems \
            and "MtrC_TIGR03507" not in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrB_TIGR03509", "MtoA", "CymA"}:
                r["cat"] = "iron_oxidation"

    elif "MtrA" in stems and "MtrB_TIGR03509" in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrA", "MtrB_TIGR03509"}:
                r["cat"] = "iron_reduction"

    elif "MtrC_TIGR03507" in stems:
        for r in updated:
            if r["hmm_stem"] in {"MtrA", "MtrB_TIGR03509"}:
                r["cat"] = "iron_reduction"

    return updated


def apply_rule(cluster_rows, rule_def):
    """
    Apply a single JSON-defined operon rule to a cluster.
    Returns the (possibly modified, possibly filtered) list of rows.
    """
    gene_set  = set(rule_def.get("genes", []))
    cat_set   = set(rule_def.get("categories", []))
    on_fail   = rule_def.get("on_fail", "keep_all")
    rule_type = rule_def["rule"]

    # Only activate if the cluster contains relevant genes/categories
    stems_present = {r["hmm_stem"] for r in cluster_rows}
    cats_present  = {r["cat"]      for r in cluster_rows}

    rule_triggered = (
        (gene_set and stems_present & gene_set) or
        (cat_set  and cats_present  & cat_set)
    )
    if not rule_triggered:
        return cluster_rows   # rule does not apply to this cluster

    # ── require_n_of ─────────────────────────────────────────────────────────
    if rule_type == "require_n_of":
        min_n    = rule_def.get("min_genes", 1)
        n_found  = _unique_genes_in(cluster_rows, gene_set)
        members  = [r for r in cluster_rows if r["hmm_stem"] in gene_set]
        non_mbrs = [r for r in cluster_rows if r["hmm_stem"] not in gene_set]
        if n_found >= min_n:
            return cluster_rows
        # Threshold not met
        if on_fail == "passthrough_non_members" and non_mbrs:
            return non_mbrs
        elif on_fail == "keep_all":
            return cluster_rows
        return []   # drop

    # ── require_anchor ────────────────────────────────────────────────────────
    if rule_type == "require_anchor":
        anchor   = rule_def.get("anchor", "")
        members  = [r for r in cluster_rows if r["hmm_stem"] in gene_set]
        non_mbrs = [r for r in cluster_rows if r["hmm_stem"] not in gene_set]
        if anchor in {r["hmm_stem"] for r in members}:
            return cluster_rows
        if on_fail == "passthrough_non_members" and non_mbrs:
            return non_mbrs
        elif on_fail == "keep_all":
            return cluster_rows
        return []

    # ── require_n_cat ─────────────────────────────────────────────────────────
    if rule_type == "require_n_cat":
        min_n      = rule_def.get("min_genes", 2)
        cat_members = _rows_in_cats(cluster_rows, cat_set)
        if len({r["hmm_stem"] for r in cat_members}) >= min_n:
            return cluster_rows
        if on_fail == "keep_all":
            return cluster_rows
        return []

    # ── require_n_cat_or_lone_trusted ─────────────────────────────────────────
    if rule_type == "require_n_cat_or_lone_trusted":
        min_n        = rule_def.get("min_genes", 2)
        trusted_set  = set(rule_def.get("trusted_lone", []))
        cat_members  = _rows_in_cats(cluster_rows, cat_set)
        unique_cat   = {r["hmm_stem"] for r in cat_members}
        if len(unique_cat) > 1:      # multiple distinct HMMs → keep
            return cluster_rows
        if trusted_set & stems_present:  # lone but trusted → keep
            return cluster_rows
        if len(unique_cat) >= min_n:
            return cluster_rows
        if on_fail == "keep_all":
            return cluster_rows
        return []

    # ── mtr_disambiguation ────────────────────────────────────────────────────
    if rule_type == "mtr_disambiguation":
        return _mtr_disambiguation(cluster_rows)

    # ── report_all / unknown rule type ────────────────────────────────────────
    return cluster_rows


def filter_cluster(cluster_rows, operon_rules, report_all_patterns,
                   all_results=False):
    """
    Run all applicable operon rules against a cluster in order.
    Each rule may reduce and/or relabel the rows.

    Returns final list of rows to keep.
    """
    if all_results:
        return cluster_rows

    # Check if entire cluster's category falls under report_all
    cats = {r["cat"] for r in cluster_rows}
    if all(_category_matches(c, report_all_patterns) for c in cats):
        return cluster_rows

    rows = cluster_rows
    for rule_def in operon_rules:
        rows = apply_rule(rows, rule_def)
        if not rows:
            return []   # cluster eliminated
    return rows


# ═══════════════════════════════════════════════════════════════════════════════
# SECOND-PASS  per-gene post-filters (Cyc1, Cyc2, regulation)
# ═══════════════════════════════════════════════════════════════════════════════

FE_REDOX_CATS = {"iron_reduction", "iron_oxidation"}


def count_heme_motifs(seq):
    """Count c-type heme-binding motifs (CXXCH variants) in a protein sequence."""
    if not seq:
        return 0
    # CX2CH, CX3CH, CX4CH, CX14CH, CX15CH  (FeGenie's set)
    patterns = [r"C(..)CH", r"C(...)CH", r"C(....)CH",
                r"C(.{14})CH", r"C(.{15})CH"]
    return sum(len(re.findall(p, seq)) for p in patterns)


def second_pass_filter(cluster_rows, gene_to_cat, seq_dict, all_results=False):
    """
    Per-gene checks applied after first-pass operon filtering.
    Handles Cyc1 co-occurrence, Cyc2 length+heme, regulation presence.
    """
    if all_results:
        return cluster_rows

    hmm_stems = [r["hmm_stem"] for r in cluster_rows]
    kept = []

    for r in cluster_rows:
        stem = r["hmm_stem"]
        cat  = r["cat"]

        # ── Cyc1: require ≥2 Fe-redox genes in same cluster ─────────────────
        if stem == "Cyc1":
            fe_count = len({h for h in set(hmm_stems)
                            if gene_to_cat.get(h, "") in FE_REDOX_CATS})
            if fe_count >= 2:
                kept.append(r)
            continue

        # ── Cyc2 variants: len ≥ 365 aa AND has CXXCH motif ─────────────────
        if re.match(r"Cyc2", stem):
            seq = seq_dict.get(r["genome"], {}).get(r["orf"], "")
            if len(seq) >= 365 and count_heme_motifs(seq) > 0:
                kept.append(r)
            continue

        # ── iron_gene_regulation: at least one regulation gene in cluster ────
        if cat == "iron_gene_regulation":
            if any("regulation" in gene_to_cat.get(h, "")
                   for h in hmm_stems):
                kept.append(r)
            continue

        kept.append(r)

    return kept


# ═══════════════════════════════════════════════════════════════════════════════
# OUTPUT WRITERS
# ═══════════════════════════════════════════════════════════════════════════════

def write_summary(path, final_rows):
    with open(path, "w") as fh:
        fh.write("category,genome/assembly,orf,gene,bitscore,"
                 "bitscore_cutoff,cluster_id,heme_c_motifs,protein_sequence\n")
        prev_cid = None
        for r in final_rows:
            if prev_cid is not None and r["cluster_id"] != prev_cid:
                fh.write("#,#,#,#,#,#,#,#\n")
            fh.write(
                f"{r['cat']},{r['genome']},{r['orf']},{r['gene_name']},"
                f"{r['bitscore']:.1f},{r['cutoff']},{r['cluster_id']},"
                f"{r['heme_motifs']},{r['sequence']}\n"
            )
            prev_cid = r["cluster_id"]


def write_gene_summary(path, final_rows):
    """FeGenie-geneSummary-clusters.csv compatible format (used by R scripts)."""
    with open(path, "w") as fh:
        fh.write("process,assembly,orf,gene,bitscore,cluster_id\n")
        prev_cid = None
        for r in final_rows:
            if prev_cid is not None and r["cluster_id"] != prev_cid:
                fh.write("#,#,#,#,#,#\n")
            fh.write(
                f"{r['cat']},{r['genome']},{r['orf']},{r['gene_name']},"
                f"{r['bitscore']:.1f},{r['cluster_id']}\n"
            )
            prev_cid = r["cluster_id"]


def write_heatmap(path, final_rows, all_genomes, norm_dict=None):
    """
    Gene-count matrix.  If norm_dict is provided (genome → total ORF count),
    values are normalised per genome (gene_count / total_orfs × 1000).
    """
    all_cats = sorted({r["cat"] for r in final_rows})

    # Count unique cluster IDs per genome per category
    count_matrix = defaultdict(lambda: defaultdict(set))
    for r in final_rows:
        count_matrix[r["cat"]][r["genome"]].add(r["cluster_id"])

    with open(path, "w") as fh:
        fh.write("X," + ",".join(all_genomes) + "\n")
        for cat in all_cats:
            row_vals = []
            for g in all_genomes:
                raw = len(count_matrix[cat].get(g, set()))
                if norm_dict and norm_dict.get(g, 0) > 0:
                    val = f"{raw / norm_dict[g] * 1000:.4f}"
                else:
                    val = str(raw)
                row_vals.append(val)
            fh.write(cat + "," + ",".join(row_vals) + "\n")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        prog="MetalGenie-Evo",
        description="FeGenie-compatible HMM annotation with improved clustering,"
                    " parallelism, and configurable operon rules",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Input
    parser.add_argument("--faa_dir",   required=True,
                        help="Directory of ORF .faa files (Prodigal output)")
    parser.add_argument("--faa_ext",   default="faa",
                        help="Extension of ORF files (without dot)")
    parser.add_argument("--gff_dir",
                        help="Directory of Prodigal .gff files (enables "
                             "coordinate-based clustering). Optional.")
    parser.add_argument("--hmm_dir",   required=True,
                        help="HMM library directory (FeGenie iron/ structure)")
    parser.add_argument("--out",       default="metalgenie_evo_out",
                        help="Output directory")
    # Clustering
    parser.add_argument("--max_gap",   type=int, default=5,
                        help="Max ORF-index gap for clustering (index mode)")
    parser.add_argument("--max_bp_gap", type=int, default=5000,
                        help="Max bp gap between gene ends/starts (GFF mode)")
    parser.add_argument("--strand_aware", action="store_true",
                        help="Split clusters at strand changes (GFF mode only)")
    # Execution
    parser.add_argument("--threads",   type=int, default=4,
                        help="Total CPU threads (split across parallel hmmsearch)")
    parser.add_argument("--hmm_threads", type=int, default=1,
                        help="Threads per individual hmmsearch call")
    # Filtering
    parser.add_argument("--all_results", action="store_true",
                        help="Report all HMM hits; skip operon-context filters")
    parser.add_argument("--norm",      action="store_true",
                        help="Normalise heatmap counts per total ORFs × 1000")
    # Misc
    parser.add_argument("--keep_tblout", action="store_true",
                        help="Keep per-genome hmmsearch .tblout files")

    args = parser.parse_args()

    faa_dir  = Path(args.faa_dir)
    hmm_dir  = Path(args.hmm_dir)
    out_dir  = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    tblout_dir = out_dir / "_tblout_cache"
    tblout_dir.mkdir(exist_ok=True)

    # ── Collect FAA files ────────────────────────────────────────────────────
    faa_files = sorted(faa_dir.glob(f"*.{args.faa_ext}"))
    if not faa_files:
        sys.exit(f"[ERROR] No .{args.faa_ext} files found in {faa_dir}")
    print(f"[INFO] {len(faa_files)} genome/bin FAA files found")

    # ── Load HMM library ─────────────────────────────────────────────────────
    cutoffs  = read_cutoffs(str(hmm_dir / "HMM-bitcutoffs.txt"))
    gene_map = read_map(str(hmm_dir / "FeGenie-map.txt"))

    cat_hmms       = defaultdict(list)   # cat → [(stem, path)]
    hmm_stem_to_cat = {}

    for entry in sorted(hmm_dir.iterdir()):
        if entry.is_dir() and not entry.name.startswith("."):
            cat = entry.name
            for hmm_file in sorted(entry.glob("*.hmm")):
                stem = hmm_file.stem
                cat_hmms[cat].append((stem, hmm_file))
                hmm_stem_to_cat[stem] = cat

    if not cat_hmms:
        sys.exit(f"[ERROR] No HMM sub-directories found in {hmm_dir}")

    total_hmms = sum(len(v) for v in cat_hmms.values())
    print(f"[INFO] {total_hmms} HMMs across {len(cat_hmms)} categories")

    operon_rules, report_all_patterns = load_operon_rules(hmm_dir)

    # ── Load GFF coordinates (optional) ──────────────────────────────────────
    genome_coords = {}
    if args.gff_dir:
        print(f"[INFO] Loading GFF coordinates from {args.gff_dir}…")
        genome_coords = load_gff_dir(args.gff_dir, faa_files)
        n_with_gff = sum(1 for f in faa_files if f.name in genome_coords)
        print(f"       {n_with_gff}/{len(faa_files)} genomes have GFF data")
    else:
        print("[INFO] No --gff_dir: using ORF-index-based clustering")

    # ── Load sequences ────────────────────────────────────────────────────────
    print("[INFO] Loading protein sequences…")
    seq_dict = {faa.name: read_fasta(str(faa)) for faa in faa_files}

    # ── Run all hmmsearches (parallel) ────────────────────────────────────────
    print(f"[INFO] Launching hmmsearch "
          f"({args.threads} total threads, {args.hmm_threads} per job)…")
    run_all_hmmsearches(
        faa_files, cat_hmms, cutoffs, tblout_dir,
        args.threads, args.hmm_threads
    )

    # ── Collect best hit per ORF ──────────────────────────────────────────────
    print("[INFO] Collecting best HMM hits per ORF…")
    best_hit = collect_best_hits(faa_files, cat_hmms, tblout_dir)

    # ── Cluster + filter ──────────────────────────────────────────────────────
    print("[INFO] Clustering and filtering…")
    cluster_id  = 0
    final_rows  = []

    for faa in faa_files:
        genome     = faa.name
        orf_hits   = best_hit.get(genome, {})
        if not orf_hits:
            continue
        orf_coords = genome_coords.get(genome, {})

        raw_clusters = build_clusters(
            genome, orf_hits, orf_coords,
            max_gap=args.max_gap,
            max_bp_gap=args.max_bp_gap,
            strand_aware=args.strand_aware,
        )

        for orf_group in raw_clusters:
            # Build row dicts for this cluster
            cluster_rows = []
            for orf in orf_group:
                hit = orf_hits.get(orf)
                if hit is None:
                    continue
                cluster_rows.append({
                    "cat":        hit["cat"],
                    "genome":     genome,
                    "orf":        orf,
                    "hmm_stem":   hit["hmm_stem"],
                    "bitscore":   hit["bitscore"],
                    "cutoff":     cutoffs.get(hit["hmm_stem"], 0),
                    "evalue":     hit["evalue"],
                    "cluster_id": cluster_id,
                })

            if not cluster_rows:
                cluster_id += 1
                continue

            # First-pass operon filter
            filtered = filter_cluster(
                cluster_rows, operon_rules, report_all_patterns,
                args.all_results
            )

            # Second-pass per-gene filter
            filtered = second_pass_filter(
                filtered, hmm_stem_to_cat, seq_dict, args.all_results
            )

            for r in filtered:
                # Attach readable gene name and sequence
                r["gene_name"]   = gene_map.get(r["hmm_stem"], r["hmm_stem"])
                r["sequence"]    = seq_dict.get(genome, {}).get(r["orf"], "")
                r["heme_motifs"] = count_heme_motifs(r["sequence"])
                final_rows.append(r)

            cluster_id += 1

    # Sort by cluster_id then orf
    final_rows.sort(key=lambda r: (r["cluster_id"], r["orf"]))

    # ── Normalisation dict ────────────────────────────────────────────────────
    norm_dict = None
    if args.norm:
        norm_dict = {faa.name: len(seq_dict.get(faa.name, {}))
                     for faa in faa_files}

    all_genomes = sorted(f.name for f in faa_files)

    # ── Write outputs ─────────────────────────────────────────────────────────
    summary_path  = out_dir / "MetalGenie-Evo-summary.csv"
    clusters_path = out_dir / "MetalGenie-Evo-geneSummary-clusters.csv"
    heatmap_path  = out_dir / "MetalGenie-Evo-heatmap-data.csv"

    print(f"[INFO] Writing {summary_path.name}…")
    write_summary(str(summary_path), final_rows)

    print(f"[INFO] Writing {clusters_path.name}…")
    write_gene_summary(str(clusters_path), final_rows)

    print(f"[INFO] Writing {heatmap_path.name}…")
    write_heatmap(str(heatmap_path), final_rows, all_genomes, norm_dict)

    # ── Cleanup tblout cache ──────────────────────────────────────────────────
    if not args.keep_tblout:
        shutil.rmtree(tblout_dir, ignore_errors=True)
    else:
        print(f"[INFO] tblout cache kept at {tblout_dir}/")

    # ── Summary stats ─────────────────────────────────────────────────────────
    n_genomes_hit = len({r["genome"] for r in final_rows})
    cat_counts    = defaultdict(int)
    for r in final_rows:
        cat_counts[r["cat"]] += 1

    print(f"\n{'─'*60}")
    print(f"  MetalGenie-Evo  –  run complete")
    print(f"  {len(final_rows)} ORFs reported across "
          f"{n_genomes_hit}/{len(faa_files)} genomes")
    print(f"\n  Hits per category:")
    for cat, n in sorted(cat_counts.items()):
        print(f"    {n:5d}  {cat}")
    print(f"\n  Outputs  →  {out_dir}/")
    print(f"{'─'*60}")


if __name__ == "__main__":
    main()
