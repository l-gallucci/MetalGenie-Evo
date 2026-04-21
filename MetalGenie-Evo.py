#!/usr/bin/env python3
"""
MetalGenie-Evo  —  HMM-based annotation of iron cycling and metal resistance genes
===================================================================================
Built on FeGenie (Garber et al. 2020) with the following extensions:

  1. GFF-based operon clustering  — uses actual bp coordinates from Prodigal GFF
     instead of ORF-index heuristics. Falls back to index-based clustering when
     GFF files are unavailable (--gff_dir not provided).

  2. Strand-aware clustering      — optional: genes on opposite strands are split
     into separate clusters before operon-rule filtering (--strand_aware).

  3. Parallel hmmsearch           — ProcessPoolExecutor dispatches one job per
     (genome, HMM) pair; --threads controls the pool size.

  4. Configurable operon rules    — loaded from operon_rules.json in the HMM
     library directory. FeGenie defaults are used as fallback.

  5. Tblout caching               — completed hmmsearch results are reused on
     subsequent runs.

  6. Normalised gene-count heatmap — --norm flag (same as FeGenie --norm).

  7. Coverage-based heatmap       — --bam / --bams for BAM-derived abundance,
     or --depth for pre-computed depth files (jgi / BBMap / samtools format).

Input
-----
  --faa_dir   Directory of ORF FASTA files (.faa) from Prodigal
  --faa_ext   Extension of ORF files (default: faa)
  --gff_dir   Directory of Prodigal GFF files (optional, same basename as .faa)
  --hmm_dir   HMM library directory
  --out       Output directory

Outputs (in --out/)
-------------------
  MetalGenie-Evo-summary.csv
  MetalGenie-Evo-geneSummary-clusters.csv   (FeGenie R-script compatible)
  MetalGenie-Evo-heatmap-data.csv           (gene counts)
  MetalGenie-Evo-coverage-heatmap.csv       (coverage-based, only with --bam/--depth)
"""

import argparse
import csv
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
# DEFAULT OPERON RULES  (used when hmm_library/operon_rules.json is absent)
# ═══════════════════════════════════════════════════════════════════════════════

_DEFAULT_OPERON_RULES = [
    {"name": "FLEET", "categories": ["iron_oxidation"],
     "genes": ["EetA","EetB","Ndh2","FmnB","FmnA","DmkA","DmkB","PplA"],
     "rule": "require_n_of", "min_genes": 5, "on_fail": "passthrough_non_members"},
    {"name": "MAM", "categories": ["magnetosome_formation"],
     "genes": ["MamA","MamB","MamE","MamK","MamP","MamM","MamQ","MamI","MamL","MamO"],
     "rule": "require_n_of", "min_genes": 5, "on_fail": "passthrough_non_members"},
    {"name": "FOXABC", "categories": ["iron_oxidation"],
     "genes": ["FoxA","FoxB","FoxC"],
     "rule": "require_n_of", "min_genes": 2, "on_fail": "passthrough_non_members"},
    {"name": "FOXEYZ", "categories": ["iron_oxidation"],
     "genes": ["FoxE","FoxY","FoxZ"],
     "rule": "require_anchor", "anchor": "FoxE", "on_fail": "passthrough_non_members"},
    {"name": "DFE1", "categories": ["iron_reduction","probable_iron_reduction"],
     "genes": ["DFE_0448","DFE_0449","DFE_0450","DFE_0451"],
     "rule": "require_n_of", "min_genes": 3, "on_fail": "passthrough_non_members"},
    {"name": "DFE2", "categories": ["iron_reduction","probable_iron_reduction"],
     "genes": ["DFE_0461","DFE_0462","DFE_0463","DFE_0464","DFE_0465"],
     "rule": "require_n_of", "min_genes": 3, "on_fail": "passthrough_non_members"},
    {"name": "MtrMto",
     "categories": ["iron_oxidation","iron_reduction",
                    "possible_iron_oxidation_and_possible_iron_reduction"],
     "genes": ["MtrA","MtrB_TIGR03509","MtrC_TIGR03507","MtoA","CymA"],
     "rule": "mtr_disambiguation", "on_fail": "keep_all"},
    {"name": "SIDERO_TRANSPORT",
     "categories": ["iron_aquisition-siderophore_transport_potential",
                    "iron_aquisition-heme_transport",
                    "iron_aquisition-siderophore_transport"],
     "genes": [],
     "rule": "require_n_cat_or_lone_trusted", "min_genes": 2,
     "trusted_lone": ["FutA1-iron_ABC_transporter_iron-binding-rep",
                      "FutA2-iron_ABC_transporter_iron-binding-rep",
                      "FutC-iron_ABC_transporter_ATPase-rep",
                      "LbtU-LvtA-PiuA-PirA-RhtA",
                      "LbtU-LbtB-legiobactin_receptor",
                      "LbtU_LbtB-legiobactin_receptor_2",
                      "IroC-salmochelin_transport-rep"],
     "on_fail": "drop"},
    {"name": "SIDERO_SYNTH",
     "categories": ["iron_aquisition-siderophore_synthesis"],
     "genes": [], "rule": "require_n_cat", "min_genes": 3, "on_fail": "drop"},
    {"name": "IRON_TRANSPORT",
     "categories": ["iron_aquisition-iron_transport","iron_aquisition-heme_oxygenase"],
     "genes": [], "rule": "require_n_cat", "min_genes": 2, "on_fail": "drop"},
]

_REPORT_ALL_PATTERNS = ["metal_resistance-*", "iron_storage"]


# ═══════════════════════════════════════════════════════════════════════════════
# IO UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

def read_fasta(path):
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
# GFF PARSING
# ═══════════════════════════════════════════════════════════════════════════════

def load_prodigal_gff(gff_path):
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
            m = re.search(r"ID=([^;]+)", parts[8])
            if m:
                orf_id = m.group(1).strip()
                coords[orf_id] = {"contig": contig, "start": start,
                                  "end": end, "strand": strand}
    return coords


def load_gff_dir(gff_dir, faa_files):
    gff_dir = Path(gff_dir)
    genome_coords = {}
    for faa in faa_files:
        stem  = Path(faa).stem
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
    parts = orf_name.rsplit("_", 1)
    if len(parts) == 2:
        try:
            return parts[0], int(parts[1])
        except ValueError:
            pass
    return orf_name, 0


def cluster_by_index(orf_set, max_gap=5):
    by_contig = defaultdict(list)
    for orf in orf_set:
        contig, idx = _index_from_name(orf)
        by_contig[contig].append((idx, orf))
    clusters = []
    for contig, entries in by_contig.items():
        entries.sort(key=lambda x: x[0])
        group = [entries[0][1]]
        for i in range(1, len(entries)):
            if entries[i][0] - entries[i-1][0] <= max_gap:
                group.append(entries[i][1])
            else:
                clusters.append(group)
                group = [entries[i][1]]
        clusters.append(group)
    return clusters


def cluster_by_coordinates(orf_set, orf_coords, max_bp_gap=5000, strand_aware=False):
    by_contig = defaultdict(list)
    for orf in orf_set:
        c = orf_coords.get(orf)
        if c is None:
            contig, idx = _index_from_name(orf)
            by_contig[contig].append((idx * 300, idx * 300, "+", orf))
        else:
            by_contig[c["contig"]].append((c["start"], c["end"], c["strand"], orf))
    clusters = []
    for contig, entries in by_contig.items():
        entries.sort(key=lambda x: x[0])
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
                if sg[i][0] - sg[i-1][1] <= max_bp_gap:
                    group.append(sg[i][3])
                else:
                    clusters.append(group)
                    group = [sg[i][3]]
            clusters.append(group)
    return clusters


def build_clusters(genome, orf_hits, orf_coords, max_gap, max_bp_gap, strand_aware):
    if orf_coords:
        return cluster_by_coordinates(orf_hits.keys(), orf_coords,
                                      max_bp_gap=max_bp_gap,
                                      strand_aware=strand_aware)
    else:
        return cluster_by_index(orf_hits.keys(), max_gap=max_gap)


# ═══════════════════════════════════════════════════════════════════════════════
# HMMER
# ═══════════════════════════════════════════════════════════════════════════════

def _hmmsearch_job(args_tuple):
    hmm_file, faa_file, tblout_path, bitscore_cutoff, threads = args_tuple
    if Path(tblout_path).exists():
        return tblout_path, True, ""
    cmd = ["hmmsearch", "--cpu", str(threads),
           "-T", str(max(bitscore_cutoff, 0)),
           "--tblout", tblout_path, "--noali",
           "-o", "/dev/null", hmm_file, faa_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return tblout_path, result.returncode == 0, result.stderr if result.returncode else ""


def parse_tblout(tblout_path):
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
    jobs = []
    for faa in faa_files:
        for cat, hmm_list in cat_hmms.items():
            for stem, hmm_path in hmm_list:
                tblout = out_tmp / f"{faa.name}__{stem}.tblout"
                jobs.append((str(hmm_path), str(faa), str(tblout),
                             cutoffs.get(stem, 0), hmm_threads))
    n_workers = max(1, threads_total // hmm_threads)
    total = len(jobs)
    done = errors = 0
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_hmmsearch_job, j): j for j in jobs}
        for fut in as_completed(futures):
            done += 1
            _, ok, err = fut.result()
            if not ok:
                errors += 1
                print(f"\n  [WARN] hmmsearch failed: {err[:80]}", file=sys.stderr)
            sys.stdout.write(f"\r[INFO] hmmsearch: {done}/{total}  ({errors} errors)  ")
            sys.stdout.flush()
    print()


def collect_best_hits(faa_files, cat_hmms, out_tmp):
    best_hit = defaultdict(dict)
    for faa in faa_files:
        genome = faa.name
        for cat, hmm_list in cat_hmms.items():
            for stem, _ in hmm_list:
                tblout = out_tmp / f"{genome}__{stem}.tblout"
                for orf, evalue, bitscore in parse_tblout(str(tblout)):
                    prev = best_hit[genome].get(orf)
                    if prev is None or bitscore > prev["bitscore"]:
                        best_hit[genome][orf] = {"hmm_stem": stem, "cat": cat,
                                                  "evalue": evalue, "bitscore": bitscore}
    return best_hit


# ═══════════════════════════════════════════════════════════════════════════════
# OPERON FILTERING
# ═══════════════════════════════════════════════════════════════════════════════

def _category_matches(cat, pattern_list):
    return any(fnmatch.fnmatch(cat, p) for p in pattern_list)


def _unique_genes_in(rows, gene_set):
    return len({r["hmm_stem"] for r in rows if r["hmm_stem"] in gene_set})


def _rows_in_cats(rows, cat_set):
    return [r for r in rows if r["cat"] in cat_set]


def _mtr_disambiguation(rows):
    stems   = {r["hmm_stem"] for r in rows}
    updated = [dict(r) for r in rows]
    if "MtoA" in stems and "MtrB_TIGR03509" in stems and "MtrC_TIGR03507" not in stems:
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
    gene_set  = set(rule_def.get("genes", []))
    cat_set   = set(rule_def.get("categories", []))
    on_fail   = rule_def.get("on_fail", "keep_all")
    rule_type = rule_def["rule"]
    stems_present = {r["hmm_stem"] for r in cluster_rows}
    cats_present  = {r["cat"]      for r in cluster_rows}
    rule_triggered = ((gene_set and stems_present & gene_set) or
                      (cat_set  and cats_present  & cat_set))
    if not rule_triggered:
        return cluster_rows

    if rule_type == "require_n_of":
        min_n    = rule_def.get("min_genes", 1)
        n_found  = _unique_genes_in(cluster_rows, gene_set)
        non_mbrs = [r for r in cluster_rows if r["hmm_stem"] not in gene_set]
        if n_found >= min_n:
            return cluster_rows
        if on_fail == "passthrough_non_members" and non_mbrs:
            return non_mbrs
        elif on_fail == "keep_all":
            return cluster_rows
        return []

    if rule_type == "require_anchor":
        anchor   = rule_def.get("anchor", "")
        non_mbrs = [r for r in cluster_rows if r["hmm_stem"] not in gene_set]
        if anchor in {r["hmm_stem"] for r in cluster_rows if r["hmm_stem"] in gene_set}:
            return cluster_rows
        if on_fail == "passthrough_non_members" and non_mbrs:
            return non_mbrs
        elif on_fail == "keep_all":
            return cluster_rows
        return []

    if rule_type == "require_n_cat":
        min_n = rule_def.get("min_genes", 2)
        if len({r["hmm_stem"] for r in _rows_in_cats(cluster_rows, cat_set)}) >= min_n:
            return cluster_rows
        return [] if on_fail != "keep_all" else cluster_rows

    if rule_type == "require_n_cat_or_lone_trusted":
        min_n       = rule_def.get("min_genes", 2)
        trusted_set = set(rule_def.get("trusted_lone", []))
        cat_members = _rows_in_cats(cluster_rows, cat_set)
        unique_cat  = {r["hmm_stem"] for r in cat_members}
        if len(unique_cat) > 1 or trusted_set & stems_present:
            return cluster_rows
        if len(unique_cat) >= min_n:
            return cluster_rows
        return [] if on_fail != "keep_all" else cluster_rows

    if rule_type == "mtr_disambiguation":
        return _mtr_disambiguation(cluster_rows)

    return cluster_rows


def filter_cluster(cluster_rows, operon_rules, report_all_patterns, all_results=False):
    if all_results:
        return cluster_rows
    cats = {r["cat"] for r in cluster_rows}
    if all(_category_matches(c, report_all_patterns) for c in cats):
        return cluster_rows
    rows = cluster_rows
    for rule_def in operon_rules:
        rows = apply_rule(rows, rule_def)
        if not rows:
            return []
    return rows


FE_REDOX_CATS = {"iron_reduction", "iron_oxidation"}


def count_heme_motifs(seq):
    if not seq:
        return 0
    patterns = [r"C(..)CH", r"C(...)CH", r"C(....)CH",
                r"C(.{14})CH", r"C(.{15})CH"]
    return sum(len(re.findall(p, seq)) for p in patterns)


def second_pass_filter(cluster_rows, gene_to_cat, seq_dict, all_results=False):
    if all_results:
        return cluster_rows
    hmm_stems = [r["hmm_stem"] for r in cluster_rows]
    kept = []
    for r in cluster_rows:
        stem = r["hmm_stem"]
        cat  = r["cat"]
        if stem == "Cyc1":
            fe_count = len({h for h in set(hmm_stems)
                            if gene_to_cat.get(h, "") in FE_REDOX_CATS})
            if fe_count >= 2:
                kept.append(r)
            continue
        if re.match(r"Cyc2", stem):
            seq = seq_dict.get(r["genome"], {}).get(r["orf"], "")
            if len(seq) >= 365 and count_heme_motifs(seq) > 0:
                kept.append(r)
            continue
        if cat == "iron_gene_regulation":
            if any("regulation" in gene_to_cat.get(h, "") for h in hmm_stems):
                kept.append(r)
            continue
        kept.append(r)
    return kept


# ═══════════════════════════════════════════════════════════════════════════════
# COVERAGE — BAM and pre-computed depth files
# ═══════════════════════════════════════════════════════════════════════════════

def _orf_to_contig(orf_name):
    """Derive contig name from a Prodigal ORF name (strip trailing _N)."""
    parts = orf_name.rsplit("_", 1)
    if len(parts) == 2:
        try:
            int(parts[1])
            return parts[0]
        except ValueError:
            pass
    return orf_name


def compute_coverage_from_bam(bam_path, out_dir, genome_label):
    """
    Run samtools coverage on a BAM file and return dict: contig → mean_depth.

    Requires samtools ≥ 1.10 (samtools coverage command).
    Falls back to samtools depth if samtools coverage is unavailable.
    """
    contig_depth = {}

    # Try samtools coverage first (samtools ≥ 1.10)
    depth_file = out_dir / f"{genome_label}.coverage"
    cmd = ["samtools", "coverage", "-H", "-o", str(depth_file), bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0 and depth_file.exists():
        # Columns: rname startpos endpos numreads covbases coverage meandepth ...
        with open(depth_file) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) >= 7:
                    try:
                        contig_depth[parts[0]] = float(parts[6])
                    except ValueError:
                        pass
        return contig_depth

    # Fallback: samtools depth (per-position, need to average)
    print(f"  [WARN] samtools coverage failed for {genome_label}, "
          f"falling back to samtools depth", file=sys.stderr)
    cmd = ["samtools", "depth", "-a", bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] samtools depth also failed for {genome_label}: "
              f"{result.stderr[:120]}", file=sys.stderr)
        return {}

    sums   = defaultdict(float)
    counts = defaultdict(int)
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            try:
                sums[parts[0]]   += float(parts[2])
                counts[parts[0]] += 1
            except ValueError:
                pass
    for contig in sums:
        contig_depth[contig] = sums[contig] / counts[contig] if counts[contig] else 0.0

    return contig_depth


def load_depth_file(depth_path):
    """
    Load a pre-computed depth file into dict: contig → mean_depth.

    Accepted formats:
      • jgi_summarize_bam_contig_depths:
        contigName  contigLen  totalAvgDepth  [sample1  sample1-var  ...]
        (header line expected; uses totalAvgDepth = column 3)

      • samtools coverage -H output:
        rname  startpos  endpos  numreads  covbases  coverage  meandepth ...
        (uses meandepth = column 7)

      • BBMap pileup.sh:
        #ID  Avg_fold  Length  Ref_GC  Covered_percent ...
        (uses Avg_fold = column 2)

      • Two-column plain:
        contig<TAB>depth   (no header)
    """
    contig_depth = {}
    if not os.path.isfile(depth_path):
        print(f"  [WARN] Depth file not found: {depth_path}", file=sys.stderr)
        return contig_depth

    with open(depth_path) as fh:
        lines = [l.rstrip() for l in fh if l.strip()]

    if not lines:
        return contig_depth

    header = lines[0]
    cols   = header.split("\t")

    # jgi format: first column header is "contigName"
    if cols[0].lower() in ("contigname", "#contigname"):
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= 3 and parts[0] != "contigName":
                try:
                    contig_depth[parts[0]] = float(parts[2])
                except ValueError:
                    pass
        return contig_depth

    # samtools coverage -H format: first col = rname, col 7 = meandepth
    if cols[0].lower() == "rname" or (len(cols) >= 7 and cols[6].lower() == "meandepth"):
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= 7:
                try:
                    contig_depth[parts[0]] = float(parts[6])
                except ValueError:
                    pass
        return contig_depth

    # BBMap pileup.sh: #ID  Avg_fold
    if cols[0].startswith("#"):
        for line in lines[1:]:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    contig_depth[parts[0]] = float(parts[1])
                except ValueError:
                    pass
        return contig_depth

    # Plain two-column: contig<TAB>depth
    for line in lines:
        parts = line.split("\t")
        if len(parts) >= 2:
            try:
                contig_depth[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return contig_depth


def load_bams_tsv(bams_path):
    """
    Load a BAM map TSV: genome_label<TAB>bam_path
    Returns dict: genome_label → bam_path
    """
    bam_map = {}
    with open(bams_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                bam_map[parts[0].strip()] = parts[1].strip()
    return bam_map


def build_contig_coverage(faa_files, bam_map, depth_map, out_dir):
    """
    Build coverage dict: genome → {contig → mean_depth}
    Works from BAM files (via samtools) or pre-computed depth files.
    bam_map  : genome_label → bam_path   (may be empty)
    depth_map: genome_label → depth_path (may be empty)
    """
    coverage_dir = out_dir / "coverage_tmp"
    coverage_dir.mkdir(exist_ok=True)

    genome_coverage = {}   # genome_name → {contig → depth}

    for faa in faa_files:
        label = faa.name    # e.g. "bin_001.faa"
        stem  = faa.stem    # e.g. "bin_001"

        # Prefer pre-computed depth file
        depth_path = depth_map.get(label) or depth_map.get(stem)
        if depth_path:
            print(f"  [INFO] Loading depth file for {stem}…")
            genome_coverage[label] = load_depth_file(depth_path)
            continue

        # Fall back to BAM
        bam_path = bam_map.get(label) or bam_map.get(stem)
        if bam_path:
            print(f"  [INFO] Computing coverage from BAM for {stem}…")
            genome_coverage[label] = compute_coverage_from_bam(
                bam_path, coverage_dir, stem
            )

    return genome_coverage


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
            fh.write(f"{r['cat']},{r['genome']},{r['orf']},{r['gene_name']},"
                     f"{r['bitscore']:.1f},{r['cutoff']},{r['cluster_id']},"
                     f"{r['heme_motifs']},{r['sequence']}\n")
            prev_cid = r["cluster_id"]


def write_gene_summary(path, final_rows):
    with open(path, "w") as fh:
        fh.write("process,assembly,orf,gene,bitscore,cluster_id\n")
        prev_cid = None
        for r in final_rows:
            if prev_cid is not None and r["cluster_id"] != prev_cid:
                fh.write("#,#,#,#,#,#\n")
            fh.write(f"{r['cat']},{r['genome']},{r['orf']},{r['gene_name']},"
                     f"{r['bitscore']:.1f},{r['cluster_id']}\n")
            prev_cid = r["cluster_id"]


def write_heatmap(path, final_rows, all_genomes, norm_dict=None):
    """Gene-count heatmap. With norm_dict: counts normalised per total ORFs × 1000."""
    all_cats     = sorted({r["cat"] for r in final_rows})
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


def write_coverage_heatmap(path, final_rows, all_genomes, genome_coverage):
    """
    Coverage-based heatmap.
    For each category × genome: sum of mean-depth values of ORF contigs.
    This mirrors FeGenie's BAM-based heatmap logic.
    """
    all_cats = sorted({r["cat"] for r in final_rows})

    # sum of contig depths per category per genome
    cov_matrix = defaultdict(lambda: defaultdict(float))
    for r in final_rows:
        genome  = r["genome"]
        contig  = _orf_to_contig(r["orf"])
        cat     = r["cat"]
        depth   = genome_coverage.get(genome, {}).get(contig, 0.0)
        cov_matrix[cat][genome] += depth

    with open(path, "w") as fh:
        fh.write("X," + ",".join(all_genomes) + "\n")
        for cat in all_cats:
            row_vals = [f"{cov_matrix[cat].get(g, 0.0):.4f}" for g in all_genomes]
            fh.write(cat + "," + ",".join(row_vals) + "\n")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        prog="MetalGenie-Evo",
        description="HMM-based annotation of iron cycling and metal resistance genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # ── Input ─────────────────────────────────────────────────────────────────
    parser.add_argument("--faa_dir",   required=True,
                        help="Directory of ORF .faa files (Prodigal output)")
    parser.add_argument("--faa_ext",   default="faa",
                        help="Extension of ORF files (without dot)")
    parser.add_argument("--gff_dir",
                        help="Directory of Prodigal .gff files (enables "
                             "coordinate-based clustering)")
    parser.add_argument("--hmm_dir",   required=True,
                        help="HMM library directory")
    parser.add_argument("--out",       default="metalgenie_evo_out",
                        help="Output directory")
    # ── Clustering ────────────────────────────────────────────────────────────
    parser.add_argument("--max_gap",    type=int, default=5,
                        help="Max ORF-index gap (index-mode clustering, "
                             "FeGenie-compatible)")
    parser.add_argument("--max_bp_gap", type=int, default=5000,
                        help="Max bp gap between gene ends (GFF mode)")
    parser.add_argument("--strand_aware", action="store_true",
                        help="Split clusters at strand changes (GFF mode only)")
    # ── Execution ─────────────────────────────────────────────────────────────
    parser.add_argument("--threads",     type=int, default=4,
                        help="Total CPU threads")
    parser.add_argument("--hmm_threads", type=int, default=1,
                        help="Threads per individual hmmsearch call")
    # ── Filtering ─────────────────────────────────────────────────────────────
    parser.add_argument("--all_results", action="store_true",
                        help="Report all HMM hits; skip operon-context filters")
    # ── Gene-count heatmap ────────────────────────────────────────────────────
    parser.add_argument("--norm", action="store_true",
                        help="Normalise gene-count heatmap per total ORFs × 1000")
    # ── Coverage heatmap (optional) ───────────────────────────────────────────
    parser.add_argument("--bam",
                        help="Single sorted BAM file for coverage heatmap. "
                             "The genome label is matched to FAA filenames by stem. "
                             "Requires samtools ≥ 1.10 in PATH.")
    parser.add_argument("--bams",
                        help="TSV file: genome_label<TAB>bam_path  (one per line). "
                             "Use when each genome/bin has its own BAM file. "
                             "Requires samtools ≥ 1.10 in PATH.")
    parser.add_argument("--depth",
                        help="Pre-computed depth file (single genome). "
                             "Accepts jgi_summarize_bam_contig_depths, "
                             "samtools coverage -H, BBMap pileup.sh, or "
                             "plain contig<TAB>depth format. "
                             "The genome label is matched to FAA filenames by stem.")
    parser.add_argument("--depths",
                        help="TSV file: genome_label<TAB>depth_path  (one per line). "
                             "Use when each genome/bin has its own depth file.")
    # ── Misc ──────────────────────────────────────────────────────────────────
    parser.add_argument("--keep_tblout", action="store_true",
                        help="Keep per-genome hmmsearch .tblout files")

    args = parser.parse_args()

    faa_dir  = Path(args.faa_dir)
    hmm_dir  = Path(args.hmm_dir)
    out_dir  = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    tblout_dir = out_dir / "_tblout_cache"
    tblout_dir.mkdir(exist_ok=True)

    # ── FAA files ─────────────────────────────────────────────────────────────
    faa_files = sorted(faa_dir.glob(f"*.{args.faa_ext}"))
    if not faa_files:
        sys.exit(f"[ERROR] No .{args.faa_ext} files found in {faa_dir}")
    print(f"[INFO] {len(faa_files)} genome/bin FAA files found")

    # ── HMM library ───────────────────────────────────────────────────────────
    # Try MetalGenie-map.txt first, fall back to FeGenie-map.txt for compatibility
    gene_map = read_map(str(hmm_dir / "MetalGenie-map.txt"))
    if not gene_map:
        gene_map = read_map(str(hmm_dir / "FeGenie-map.txt"))

    cutoffs        = read_cutoffs(str(hmm_dir / "HMM-bitcutoffs.txt"))
    cat_hmms       = defaultdict(list)
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

    # ── GFF coordinates (optional) ────────────────────────────────────────────
    genome_coords = {}
    if args.gff_dir:
        print(f"[INFO] Loading GFF coordinates from {args.gff_dir}…")
        genome_coords = load_gff_dir(args.gff_dir, faa_files)
        n_gff = sum(1 for f in faa_files if f.name in genome_coords)
        print(f"       {n_gff}/{len(faa_files)} genomes have GFF data")
    else:
        print("[INFO] No --gff_dir: using ORF-index-based clustering "
              "(FeGenie-compatible)")

    # ── Load sequences ────────────────────────────────────────────────────────
    print("[INFO] Loading protein sequences…")
    seq_dict = {faa.name: read_fasta(str(faa)) for faa in faa_files}

    # ── hmmsearch ─────────────────────────────────────────────────────────────
    print(f"[INFO] Launching hmmsearch "
          f"({args.threads} total threads, {args.hmm_threads} per job)…")
    run_all_hmmsearches(faa_files, cat_hmms, cutoffs, tblout_dir,
                        args.threads, args.hmm_threads)

    # ── Best hit per ORF ──────────────────────────────────────────────────────
    print("[INFO] Collecting best HMM hits per ORF…")
    best_hit = collect_best_hits(faa_files, cat_hmms, tblout_dir)

    # ── Cluster + filter ──────────────────────────────────────────────────────
    print("[INFO] Clustering and filtering…")
    cluster_id = 0
    final_rows = []

    for faa in faa_files:
        genome    = faa.name
        orf_hits  = best_hit.get(genome, {})
        if not orf_hits:
            continue
        orf_coords = genome_coords.get(genome, {})

        raw_clusters = build_clusters(
            genome, orf_hits, orf_coords,
            max_gap=args.max_gap, max_bp_gap=args.max_bp_gap,
            strand_aware=args.strand_aware,
        )

        for orf_group in raw_clusters:
            cluster_rows = []
            for orf in orf_group:
                hit = orf_hits.get(orf)
                if hit is None:
                    continue
                cluster_rows.append({
                    "cat": hit["cat"], "genome": genome, "orf": orf,
                    "hmm_stem": hit["hmm_stem"], "bitscore": hit["bitscore"],
                    "cutoff": cutoffs.get(hit["hmm_stem"], 0),
                    "evalue": hit["evalue"], "cluster_id": cluster_id,
                })
            if not cluster_rows:
                cluster_id += 1
                continue

            filtered = filter_cluster(cluster_rows, operon_rules,
                                       report_all_patterns, args.all_results)
            filtered = second_pass_filter(filtered, hmm_stem_to_cat,
                                          seq_dict, args.all_results)

            for r in filtered:
                r["gene_name"]   = gene_map.get(r["hmm_stem"], r["hmm_stem"])
                r["sequence"]    = seq_dict.get(genome, {}).get(r["orf"], "")
                r["heme_motifs"] = count_heme_motifs(r["sequence"])
                final_rows.append(r)

            cluster_id += 1

    final_rows.sort(key=lambda r: (r["cluster_id"], r["orf"]))

    # ── Normalisation dict (gene counts) ──────────────────────────────────────
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

    # ── Coverage heatmap (optional) ───────────────────────────────────────────
    bam_map   = {}
    depth_map = {}

    if args.bams:
        bam_map = load_bams_tsv(args.bams)
        print(f"[INFO] Loaded BAM map: {len(bam_map)} entries from {args.bams}")
    elif args.bam:
        # Single BAM: match by stem to all FAA files
        bam_stem = Path(args.bam).stem
        for faa in faa_files:
            bam_map[faa.name] = args.bam
            bam_map[faa.stem] = args.bam
        print(f"[INFO] Single BAM file: {args.bam}")

    if args.depths:
        depth_map = load_bams_tsv(args.depths)   # same two-column TSV format
        print(f"[INFO] Loaded depth map: {len(depth_map)} entries")
    elif args.depth:
        for faa in faa_files:
            depth_map[faa.name] = args.depth
            depth_map[faa.stem] = args.depth
        print(f"[INFO] Single depth file: {args.depth}")

    if bam_map or depth_map:
        print("[INFO] Computing coverage-based heatmap…")
        genome_coverage = build_contig_coverage(
            faa_files, bam_map, depth_map, out_dir
        )
        if genome_coverage:
            cov_path = out_dir / "MetalGenie-Evo-coverage-heatmap.csv"
            print(f"[INFO] Writing {cov_path.name}…")
            write_coverage_heatmap(str(cov_path), final_rows,
                                   all_genomes, genome_coverage)
        else:
            print("[WARN] No coverage data could be loaded.", file=sys.stderr)

    # ── Cleanup ───────────────────────────────────────────────────────────────
    if not args.keep_tblout:
        shutil.rmtree(tblout_dir, ignore_errors=True)
    else:
        print(f"[INFO] tblout cache kept at {tblout_dir}/")

    # ── Summary ───────────────────────────────────────────────────────────────
    n_genomes_hit = len({r["genome"] for r in final_rows})
    cat_counts    = defaultdict(int)
    for r in final_rows:
        cat_counts[r["cat"]] += 1

    print(f"\n{'─'*60}")
    print(f"  MetalGenie-Evo  —  run complete")
    print(f"  {len(final_rows)} ORFs reported across "
          f"{n_genomes_hit}/{len(faa_files)} genomes")
    print(f"\n  Hits per category:")
    for cat, n in sorted(cat_counts.items()):
        print(f"    {n:5d}  {cat}")
    print(f"\n  Outputs  →  {out_dir}/")
    print(f"{'─'*60}")


if __name__ == "__main__":
    main()
