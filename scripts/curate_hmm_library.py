#!/usr/bin/env python3
"""
curate_hmm_library.py
=====================
Standalone curation script for the MetalGenie-Evo HMM library.

PURPOSE
-------
This script is NOT part of MetalGenie-Evo's runtime. It is provided so that
anyone can:
  1. Inspect how the bundled hmm_library/ was built
  2. Reproduce it from scratch using the original sources
  3. Verify deduplication decisions
  4. Add new HMM sources and regenerate the library

It takes all raw HMM sources (FeGenie, MetHMMDB, Tabuteau/iron_cheese,
InterPro/NCBI NF*, custom), assigns them to FeGenie-style pathway categories,
deduplicates across sources, and writes the final hmm_library/ directory
together with a full provenance log (hmm_registry.tsv).

SOURCES
-------
  1. FeGenie original          already organised by category subdirectory
  2. new_iron_hmms             flat dir — Tabuteau et al. 2025
                               (KOfam + custom HMMs for iron acquisition)
  3. others                    flat dir — NCBI NF* and Pfam PF* downloads
  4. MetHMMDB individual/      flat dir — metal mobility resistance genes

OUTPUT STRUCTURE (mirrors FeGenie iron/ exactly)
-------------------------------------------------
  hmm_library/
    iron_reduction/
    iron_oxidation/
    iron_storage/
    iron_gene_regulation/
    magnetosome_formation/
    possible_iron_oxidation_and_possible_iron_reduction/
    probable_iron_reduction/
    iron_aquisition-iron_transport/
    iron_aquisition-heme_transport/
    iron_aquisition-heme_oxygenase/
    iron_aquisition-siderophore_synthesis/
    iron_aquisition-siderophore_transport/
    iron_aquisition-siderophore_transport_potential/
    metal_resistance-arsenic/
    metal_resistance-copper/
    metal_resistance-mercury/
    metal_resistance-cobalt_zinc_cadmium/
    metal_resistance-chromium/
    metal_resistance-lead/
    metal_resistance-nickel/
    metal_resistance-silver/
    metal_resistance-tellurite/
    metal_resistance-antimony/
    HMM-bitcutoffs.txt
    FeGenie-map.txt
    operon_rules.json       (copied from source repo if present)
    hmm_registry.tsv        (full provenance log)

USAGE
-----
# Full build from scratch (developer / first time)
python scripts/curate_hmm_library.py \\
    --fegenie_dir  /path/to/FeGenie/hmms/iron/ \\
    --flat_dir     /path/to/new_iron_hmms   iron_aquisition \\
    --flat_dir     /path/to/others          iron_aquisition \\
    --methmmdb_dir /path/to/MetHMMDB/individual/ \\
    --out_dir      hmm_library/ \\
    --log          curation_report.tsv

# Verify an existing library against its registry
python scripts/curate_hmm_library.py --verify hmm_library/

# Re-run after editing curation_report.tsv to fix REVIEW_NEEDED rows
python scripts/curate_hmm_library.py \\
    --from_report  curation_report.tsv \\
    --out_dir      hmm_library/
"""

import argparse
import csv
import os
import re
import shutil
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


# ═══════════════════════════════════════════════════════════════════════════════
# VALID CATEGORIES
# ═══════════════════════════════════════════════════════════════════════════════

VALID_CATEGORIES = [
    # ── FeGenie original ──────────────────────────────────────────────────────
    "iron_reduction",
    "iron_oxidation",
    "iron_storage",
    "iron_gene_regulation",
    "magnetosome_formation",
    "possible_iron_oxidation_and_possible_iron_reduction",
    "probable_iron_reduction",
    "iron_aquisition-iron_transport",
    "iron_aquisition-heme_transport",
    "iron_aquisition-heme_oxygenase",
    "iron_aquisition-siderophore_synthesis",
    "iron_aquisition-siderophore_transport",
    "iron_aquisition-siderophore_transport_potential",
    # ── MetHMMDB metal resistance ─────────────────────────────────────────────
    "metal_resistance-arsenic",
    "metal_resistance-copper",
    "metal_resistance-mercury",
    "metal_resistance-cobalt_zinc_cadmium",
    "metal_resistance-chromium",
    "metal_resistance-lead",
    "metal_resistance-nickel",
    "metal_resistance-silver",
    "metal_resistance-gold",
    "metal_resistance-tellurite",
    "metal_resistance-antimony",
    "metal_resistance-bismuth",
    "metal_resistance-unknown",
]

REVIEW_NEEDED = "REVIEW_NEEDED"


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY LOOKUP TABLES
# ═══════════════════════════════════════════════════════════════════════════════

# KOfam K-numbers → MetalGenie-Evo category
# Source: Tabuteau et al. 2025 Table S3 + KEGG annotations
KOFAM_CAT = {
    # Biosynthesis
    "K00216": "iron_aquisition-siderophore_synthesis",
    "K00257": "iron_aquisition-siderophore_synthesis",
    "K00824": "iron_aquisition-siderophore_synthesis",
    "K01252": "iron_aquisition-siderophore_synthesis",
    "K01582": "iron_aquisition-siderophore_synthesis",
    "K01780": "iron_aquisition-siderophore_synthesis",
    "K01851": "iron_aquisition-siderophore_synthesis",
    "K01909": "iron_aquisition-siderophore_synthesis",
    "K04759": "iron_aquisition-siderophore_synthesis",
    "K04788": "iron_aquisition-siderophore_synthesis",   # entA
    "K04789": "iron_aquisition-siderophore_synthesis",   # entB
    "K04790": "iron_aquisition-siderophore_synthesis",   # entC
    "K04791": "iron_aquisition-siderophore_synthesis",   # entD
    "K04792": "iron_aquisition-siderophore_synthesis",   # entE
    "K04793": "iron_aquisition-siderophore_synthesis",   # entF
    "K05375": "iron_aquisition-siderophore_synthesis",   # iucA aerobactin
    "K07116": "iron_aquisition-siderophore_synthesis",   # iucB
    "K07224": "iron_aquisition-siderophore_synthesis",   # iucC
    "K07230": "iron_aquisition-siderophore_synthesis",   # iucD
    "K10531": "iron_aquisition-siderophore_synthesis",   # dhbA
    "K11604": "iron_aquisition-siderophore_synthesis",
    "K12237": "iron_aquisition-siderophore_synthesis",   # pvdH
    "K12238": "iron_aquisition-siderophore_synthesis",   # pvdA
    "K12239": "iron_aquisition-siderophore_synthesis",
    "K12240": "iron_aquisition-siderophore_synthesis",
    "K12241": "iron_aquisition-siderophore_synthesis",
    "K12346": "iron_aquisition-siderophore_synthesis",
    "K13745": "iron_aquisition-siderophore_synthesis",
    "K14759": "iron_aquisition-siderophore_synthesis",
    "K15584": "iron_aquisition-siderophore_synthesis",
    "K15652": "iron_aquisition-siderophore_synthesis",
    "K15653": "iron_aquisition-siderophore_synthesis",
    "K15681": "iron_aquisition-siderophore_synthesis",
    "K15721": "iron_aquisition-siderophore_synthesis",
    "K16088": "iron_aquisition-siderophore_synthesis",
    "K16089": "iron_aquisition-siderophore_synthesis",
    "K16090": "iron_aquisition-siderophore_synthesis",
    "K16091": "iron_aquisition-siderophore_synthesis",
    "K19611": "iron_aquisition-siderophore_synthesis",
    "K19791": "iron_aquisition-siderophore_synthesis",
    "K19792": "iron_aquisition-siderophore_synthesis",
    "K21476": "iron_aquisition-siderophore_synthesis",
    "K21721": "iron_aquisition-siderophore_synthesis",
    "K21949": "iron_aquisition-siderophore_synthesis",
    "K21992": "iron_aquisition-siderophore_synthesis",
    "K22148": "iron_aquisition-siderophore_synthesis",
    "K22149": "iron_aquisition-siderophore_synthesis",
    "K22150": "iron_aquisition-siderophore_synthesis",
    "K22151": "iron_aquisition-siderophore_synthesis",
    "K22152": "iron_aquisition-siderophore_synthesis",
    "K22153": "iron_aquisition-siderophore_synthesis",
    "K23120": "iron_aquisition-siderophore_synthesis",
    "K23122": "iron_aquisition-siderophore_synthesis",
    "K23181": "iron_aquisition-siderophore_synthesis",
    "K23185": "iron_aquisition-siderophore_synthesis",
    "K23227": "iron_aquisition-siderophore_synthesis",
    "K23371": "iron_aquisition-siderophore_synthesis",
    "K23372": "iron_aquisition-siderophore_synthesis",
    "K23373": "iron_aquisition-siderophore_synthesis",
    "K23374": "iron_aquisition-siderophore_synthesis",
    "K23375": "iron_aquisition-siderophore_synthesis",
    "K23376": "iron_aquisition-siderophore_synthesis",
    "K23446": "iron_aquisition-siderophore_synthesis",
    "K23447": "iron_aquisition-siderophore_synthesis",
    "K23725": "iron_aquisition-siderophore_synthesis",
    "K24101": "iron_aquisition-siderophore_synthesis",
    "K24108": "iron_aquisition-siderophore_synthesis",
    "K24109": "iron_aquisition-siderophore_synthesis",
    "K24110": "iron_aquisition-siderophore_synthesis",
    "K24111": "iron_aquisition-siderophore_synthesis",
    "K24112": "iron_aquisition-siderophore_synthesis",
    "K25109": "iron_aquisition-siderophore_synthesis",
    "K25113": "iron_aquisition-siderophore_synthesis",
    "K25282": "iron_aquisition-siderophore_synthesis",
    "K25286": "iron_aquisition-siderophore_synthesis",
    "K25287": "iron_aquisition-siderophore_synthesis",
    "K25308": "iron_aquisition-siderophore_synthesis",
    # Transport / import
    "K02012": "iron_aquisition-siderophore_transport",   # fepA
    "K02014": "iron_aquisition-siderophore_transport",   # fepB
    "K02016": "iron_aquisition-siderophore_transport",   # fepC
    "K02361": "iron_aquisition-siderophore_transport",   # fepD
    "K02362": "iron_aquisition-siderophore_transport",   # fepE
    "K02363": "iron_aquisition-siderophore_transport",   # fepG
    "K02364": "iron_aquisition-siderophore_transport",   # fes
    "K03894": "iron_aquisition-siderophore_transport",   # fecA
    "K03895": "iron_aquisition-siderophore_transport",   # fecB
    "K03896": "iron_aquisition-siderophore_transport",   # fecC
    "K03897": "iron_aquisition-siderophore_transport",   # fecD
    "K04778": "iron_aquisition-siderophore_transport",   # iutA
    "K04780": "iron_aquisition-siderophore_transport",   # fhuA
    "K04781": "iron_aquisition-siderophore_transport",   # fhuB
    "K04782": "iron_aquisition-siderophore_transport",   # fhuC
    "K04783": "iron_aquisition-siderophore_transport",   # fhuD
    "K04784": "iron_aquisition-siderophore_transport",   # fhuE
    "K04785": "iron_aquisition-siderophore_transport",   # fhuF
    "K04786": "iron_aquisition-siderophore_transport",   # fiu
    "K04787": "iron_aquisition-siderophore_transport",   # cirA
    "K07243": "iron_aquisition-siderophore_transport",   # iutA/irgA
    "K08197": "iron_aquisition-siderophore_transport",   # viuA
}

# NCBI Protein Family Models
NCBI_NF_CAT = {
    "NF037942": "iron_aquisition-siderophore_synthesis",
    "NF040982": "iron_aquisition-siderophore_transport",
    "NF040984": "iron_aquisition-siderophore_transport",
    "NF040985": "iron_aquisition-siderophore_transport",
    "NF040986": "iron_aquisition-siderophore_transport",
    "NF040987": "iron_aquisition-siderophore_transport",
    "NF041018": "iron_aquisition-iron_transport",
    "NF041019": "iron_aquisition-iron_transport",
    "NF041020": "iron_aquisition-iron_transport",
}

PFAM_CAT = {
    "PF16525": "iron_aquisition-siderophore_synthesis",
}

# MetHMMDB metal keyword → category
METAL_KW_CAT = {
    "arsenic": "metal_resistance-arsenic",
    "arsb": "metal_resistance-arsenic",
    "arsc": "metal_resistance-arsenic",
    "arsd": "metal_resistance-arsenic",
    "arse": "metal_resistance-arsenic",
    "copper": "metal_resistance-copper",
    "czc": "metal_resistance-cobalt_zinc_cadmium",
    "cobalt": "metal_resistance-cobalt_zinc_cadmium",
    "zinc": "metal_resistance-cobalt_zinc_cadmium",
    "cadmium": "metal_resistance-cobalt_zinc_cadmium",
    "mercury": "metal_resistance-mercury",
    "mer": "metal_resistance-mercury",
    "chromate": "metal_resistance-chromium",
    "chromium": "metal_resistance-chromium",
    "lead": "metal_resistance-lead",
    "nickel": "metal_resistance-nickel",
    "silver": "metal_resistance-silver",
    "gold": "metal_resistance-gold",
    "tellurite": "metal_resistance-tellurite",
    "antimony": "metal_resistance-antimony",
    "bismuth": "metal_resistance-bismuth",
}

# Keyword heuristics for iron acquisition flat dirs
_SYNTH_KW = [
    "synthetase", "synthase", "biosynthesis", "nrps", "nonribosomal",
    "pvd", "pyoverdine", "fusarinine", "coprogen", "rhizobactin",
    "rhizoferrin", "desferrioxamine", "ferrichrome", "aerobactin",
    "staphyloferrin", "vibrioferrin", "iuc", "dhb", "entA", "entB",
    "entC", "entD", "entE", "entF", "fsl", "sid", "sidc",
    "adenylation", "condensation", "thiolation", "isochorismate",
    "rumb", "rhbf", "desd",
]
_TRANSPORT_KW = [
    "receptor", "tonb", "outer membrane", "abc transporter",
    "binding protein", "permease", "import", "uptake",
    "fhua", "fhub", "fhuc", "fhud", "fhue",
    "fepa", "fepb", "fepc", "fepd", "fepg",
    "feca", "fecb", "fecc", "fecd",
    "iuta", "irga", "cira", "fiu",
    "futa", "futb", "futc", "viua", "fata", "pvdr",
    "fcua", "fite", "fiua", "lbtu",
]
_HEME_KW = ["heme", "haemin", "haemophore", "hasa", "husa", "hmbr",
            "shua", "bhua"]


# ═══════════════════════════════════════════════════════════════════════════════
# HMMER3 PARSER
# ═══════════════════════════════════════════════════════════════════════════════

def parse_hmm_meta(hmm_path):
    """Parse NAME, ACC, DESC, NSEQ, GA, TC, DATE from a single .hmm file.

    DATE field in HMMER3 format: 'DATE  Fri Mar 18 09:00:00 2022'
    Stored as a datetime object (or None if absent/unparseable).
    """
    from datetime import datetime as _dt
    m = {"name": "", "acc": "", "desc": "", "nseq": 0,
         "ga_seq": None, "tc_seq": None, "date": None}
    with open(hmm_path, errors="replace") as fh:
        for line in fh:
            t = line.split()[0] if line.split() else ""
            if t == "NAME":
                m["name"] = line.split(None, 1)[1].strip()
            elif t == "ACC":
                m["acc"]  = line.split(None, 1)[1].strip()
            elif t == "DESC":
                m["desc"] = (line.split(None, 1)[1].strip()
                             if len(line.split(None, 1)) > 1 else "")
            elif t == "NSEQ":
                try:
                    m["nseq"] = int(line.split()[1])
                except (IndexError, ValueError):
                    pass
            elif t == "GA":
                parts = line.split()
                try:
                    m["ga_seq"] = float(parts[1].rstrip(";"))
                except (IndexError, ValueError):
                    pass
            elif t == "TC":
                parts = line.split()
                try:
                    m["tc_seq"] = float(parts[1].rstrip(";"))
                except (IndexError, ValueError):
                    pass
            elif t == "DATE":
                # e.g.  DATE  Fri Mar 18 09:00:00 2022
                date_str = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else ""
                for fmt in ("%a %b %d %H:%M:%S %Y",   # Fri Mar 18 09:00:00 2022
                            "%a %b  %d %H:%M:%S %Y",  # single-digit day with extra space
                            "%Y-%m-%d",                # ISO fallback
                            "%d/%m/%Y"):               # European fallback
                    try:
                        m["date"] = _dt.strptime(date_str, fmt)
                        break
                    except ValueError:
                        continue
            elif t == "//":
                break
    return m


# ═══════════════════════════════════════════════════════════════════════════════
# CUTOFF LOADER
# ═══════════════════════════════════════════════════════════════════════════════

def load_cutoff_tsv(tsv_path):
    """
    Load external bitscore cutoffs from a TSV file.

    Accepts two formats:
      1. Two-column, no header:      stem<TAB>cutoff
      2. With header row: any TSV that has columns named 'stem' (or 'name' or
         'hmm') and 'cutoff' (or 'bitscore' or 'threshold' or 'score').
         Extra columns (gene_name, category, description…) are silently ignored.

    Returns dict: stem → float cutoff.

    Typical sources:
      • Tabuteau et al. 2025  Table S3  (HMM_name / bitscore columns)
      • FeGenie HMM-bitcutoffs.txt      (already tab-separated, no header)
      • Any custom table you export from Excel
    """
    cutoffs = {}
    if not os.path.isfile(tsv_path):
        print(f"  [WARN] Cutoff file not found: {tsv_path}", file=sys.stderr)
        return cutoffs

    with open(tsv_path, errors="replace") as fh:
        first_line = fh.readline().rstrip()
        fh.seek(0)

        # Detect whether there is an explicit header row.
        # Heuristic: a header contains recognisable column-name tokens.
        # FeGenie-style files (stem<TAB>cutoff, no header) are two-column
        # where the second field is always numeric — detect that case first.
        cols = first_line.split("\t")
        second_is_numeric = len(cols) >= 2 and re.match(r"^[\d.]+$", cols[1].strip())
        known_header_words = {"stem", "name", "hmm", "cutoff", "bitscore",
                              "threshold", "score", "model", "hmm_name"}
        first_is_header_word = cols[0].strip().lower() in known_header_words
        has_header = first_is_header_word and not second_is_numeric

        if has_header:
            reader = csv.DictReader(fh, delimiter="\t")
            # Normalise column names: lower-case, strip whitespace
            rows = []
            for row in reader:
                rows.append({k.lower().strip(): (v or "").strip()
                             for k, v in row.items()
                             if k is not None})

            # Find stem column — exact match first, then substring
            avail_cols = rows[0].keys() if rows else {}
            _STEM_NAMES  = ["stem", "name", "hmm", "hmm_name", "model"]
            _CUT_NAMES   = ["cutoff", "bitscore", "threshold", "score"]

            stem_col = next(
                (c for c in _STEM_NAMES if c in avail_cols), None
            )
            if not stem_col:   # fallback: first column whose name contains a stem keyword
                stem_col = next(
                    (c for c in avail_cols
                     if any(k in c for k in _STEM_NAMES)), None
                )

            cutoff_col = next(
                (c for c in _CUT_NAMES if c in avail_cols), None
            )
            if not cutoff_col:  # fallback: any column whose name contains "cutoff" or "score"
                cutoff_col = next(
                    (c for c in avail_cols
                     if any(k in c for k in _CUT_NAMES)), None
                )

            if not stem_col or not cutoff_col:
                print(f"  [WARN] Could not identify stem/cutoff columns in {tsv_path}.",
                      file=sys.stderr)
                print(f"         Columns found: {list(rows[0].keys()) if rows else '(empty)'}",
                      file=sys.stderr)
                return cutoffs

            for row in rows:
                stem = row[stem_col].strip()
                try:
                    cutoffs[stem] = float(row[cutoff_col])
                except (ValueError, KeyError):
                    pass
        else:
            # Simple two-column: stem<TAB>cutoff  (like FeGenie-bitcutoffs.txt)
            for line in fh:
                ls = line.rstrip().split("\t")
                if len(ls) >= 2:
                    try:
                        cutoffs[ls[0].strip()] = float(ls[1].strip())
                    except ValueError:
                        pass

    return cutoffs


def apply_cutoffs(models, fegenie_cutoffs, extra_cutoffs=None):
    """
    Apply bitscore cutoffs to every model.

    Priority:
      1. extra_cutoffs   (user-supplied --cutoff_tsv, optional)
      2. fegenie_cutoffs (auto-loaded from HMM-bitcutoffs.txt — covers Tabuteau
                          names since they are already listed there)
      3. GA embedded in the HMM header
      4. TC embedded in the HMM header
      5. 0.0  (fallback — rely on E-value < 0.1 only)

    fegenie_cutoffs is always passed; extra_cutoffs may be empty/None.
    """
    if extra_cutoffs is None:
        extra_cutoffs = {}

    stats = {"extra": 0, "fegenie_txt": 0,
             "embedded_ga": 0, "embedded_tc": 0, "missing": 0}

    for m in models:
        stem = m["stem"]
        if stem in extra_cutoffs:
            m["cutoff"] = extra_cutoffs[stem]
            m["cutoff_source"] = "extra_tsv"
            stats["extra"] += 1
        elif stem in fegenie_cutoffs:
            m["cutoff"] = fegenie_cutoffs[stem]
            m["cutoff_source"] = "fegenie_bitcutoffs"
            stats["fegenie_txt"] += 1
        elif m.get("ga_seq") and float(m["ga_seq"]) > 0:
            m["cutoff"] = float(m["ga_seq"])
            m["cutoff_source"] = "embedded_GA"
            stats["embedded_ga"] += 1
        elif m.get("tc_seq") and float(m["tc_seq"]) > 0:
            m["cutoff"] = float(m["tc_seq"])
            m["cutoff_source"] = "embedded_TC"
            stats["embedded_tc"] += 1
        else:
            m["cutoff"] = 0.0
            m["cutoff_source"] = "none"
            stats["missing"] += 1

    return models, stats


# ═══════════════════════════════════════════════════════════════════════════════
# CATEGORY INFERENCE
# ═══════════════════════════════════════════════════════════════════════════════

def _kw(text, words):
    t = text.lower()
    return any(w.lower() in t for w in words)


def infer_category(stem, name, desc, source_hint="iron_aquisition",
                   pre_defined=None):
    """
    Infer MetalGenie-Evo category for a single HMM.
    Returns (category, confidence: 'high'|'low'|'review').
    """
    combined = f"{stem} {name} {desc}"

    # User-supplied map takes absolute precedence
    if pre_defined and stem in pre_defined:
        return pre_defined[stem], "high"
    ko = re.match(r"^(K\d{5})", stem)
    if ko and ko.group(1) in KOFAM_CAT:
        return KOFAM_CAT[ko.group(1)], "high"

    # NCBI NF
    nf = re.match(r"^(NF\d+)", stem)
    if nf and nf.group(1) in NCBI_NF_CAT:
        return NCBI_NF_CAT[nf.group(1)], "high"

    # Pfam PF
    pf = re.match(r"^(PF\d+)", stem)
    if pf and pf.group(1) in PFAM_CAT:
        return PFAM_CAT[pf.group(1)], "high"

    # MetHMMDB metal keywords (check name first, then desc)
    for kw, cat in METAL_KW_CAT.items():
        if kw in stem.lower() or kw in name.lower() or kw in desc.lower():
            return cat, "high"

    # Iron acquisition keyword heuristics
    if "iron_aquisition" in source_hint or "iron" in source_hint:
        if _kw(combined, _HEME_KW):
            return "iron_aquisition-heme_transport", "low"
        if _kw(combined, _SYNTH_KW):
            return "iron_aquisition-siderophore_synthesis", "low"
        if _kw(combined, _TRANSPORT_KW):
            return "iron_aquisition-siderophore_transport", "low"

    # Tabuteau custom names
    lower = stem.lower()
    synth_stems = ["fusarinine", "tafc", "nrps", "rhizoferrin", "copfam",
                   "fer3", "rhbf", "desd", "fso1", "rumb", "rhizobactin",
                   "sidc", "put_sid", "all_seq", "put_nrps", "put_cop",
                   "putrhiz", "put_b4eug", "put_b4euh", "put_d4au",
                   "put_crl", "put_fcua", "put_ferr", "put_tafc",
                   "put_rumb", "put_wp", "put_q6csn", "put_q5d6",
                   "put_i1rn", "put_rhb", "put_sid_a", "f9x9v1"]
    transport_stems = ["fcua", "fhua", "fite", "fiua", "put_fhua",
                       "put_fhue", "put_fitA", "put_fiua", "put_b4exx",
                       "put_b4exy", "put_q6csn3"]
    if any(s in lower for s in synth_stems):
        return "iron_aquisition-siderophore_synthesis", "low"
    if any(s in lower for s in transport_stems):
        return "iron_aquisition-siderophore_transport", "low"

    return REVIEW_NEEDED, "review"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE LOADERS
# ═══════════════════════════════════════════════════════════════════════════════

def load_fegenie(fegenie_dir):
    """Load FeGenie iron/ directory — already categorised."""
    base = Path(fegenie_dir)
    cutoffs, gene_map = {}, {}

    cp = base / "HMM-bitcutoffs.txt"
    if cp.exists():
        for line in open(cp):
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                try:
                    cutoffs[ls[0]] = float(ls[1])
                except ValueError:
                    pass

    mp = base / "FeGenie-map.txt"
    if mp.exists():
        for line in open(mp):
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                gene_map[ls[0]] = ls[1]

    models = []
    for cat_dir in sorted(base.iterdir()):
        if not cat_dir.is_dir() or cat_dir.name.startswith("."):
            continue
        for hmm_file in sorted(cat_dir.glob("*.hmm")):
            meta = parse_hmm_meta(str(hmm_file))
            stem = hmm_file.stem
            models.append({
                "stem":       stem,
                "hmm_file":   str(hmm_file),
                "name":       meta["name"] or stem,
                "acc":        meta["acc"],
                "desc":       meta["desc"],
                "nseq":       meta["nseq"],
                "ga_seq":     meta["ga_seq"],
                "tc_seq":     meta["tc_seq"],
                "date":       meta["date"],
                "category":   cat_dir.name,
                "confidence": "high",
                "gene_name":  gene_map.get(stem, stem),
                "cutoff":     cutoffs.get(stem, 0.0),
                "source":     "fegenie",
            })
    return models


def load_flat_dir(dir_path, hint, source_label):
    """Load a flat directory of .hmm files with category inference."""
    models = []
    for hmm_file in sorted(Path(dir_path).glob("*.hmm")):
        meta  = parse_hmm_meta(str(hmm_file))
        stem  = hmm_file.stem
        cat, conf = infer_category(stem, meta["name"] or stem, meta["desc"], hint)
        cutoff = meta.get("ga_seq") or meta.get("tc_seq") or 0.0
        models.append({
            "stem":       stem,
            "hmm_file":   str(hmm_file),
            "name":       meta["name"] or stem,
            "acc":        meta["acc"],
            "desc":       meta["desc"],
            "nseq":       meta["nseq"],
            "ga_seq":     meta["ga_seq"],
            "tc_seq":     meta["tc_seq"],
            "date":       meta["date"],
            "category":   cat,
            "confidence": conf,
            "gene_name":  re.sub(r"_\d+$", "", stem),
            "cutoff":     cutoff,
            "source":     source_label,
        })
    return models


def load_methmmdb_individual(dir_path):
    """Load MetHMMDB individual/ folder — keyword heuristics only (fallback)."""
    return load_flat_dir(dir_path, "metal_resistance", "methmmdb")


# MetHMMDB metal_type field → MetalGenie-Evo category
_METAL_TYPE_TO_CAT = {
    "arsenic":    "metal_resistance-arsenic",
    "copper":     "metal_resistance-copper",
    "mercury":    "metal_resistance-mercury",
    "zinc":       "metal_resistance-cobalt_zinc_cadmium",
    "cadmium":    "metal_resistance-cobalt_zinc_cadmium",
    "cobalt":     "metal_resistance-cobalt_zinc_cadmium",
    "chromium":   "metal_resistance-chromium",
    "chromate":   "metal_resistance-chromium",
    "lead":       "metal_resistance-lead",
    "nickel":     "metal_resistance-nickel",
    "silver":     "metal_resistance-silver",
    "gold":       "metal_resistance-gold",
    "tellurite":  "metal_resistance-tellurite",
    "tellurium":  "metal_resistance-tellurite",
    "antimony":   "metal_resistance-antimony",
    "bismuth":    "metal_resistance-bismuth",
    "iron":       "iron_resistance",
    "multimetal": "metal_resistance-multimetal",
}


def _metal_type_to_category(metal_type_list):
    """
    Convert MetHMMDB metal_type list (e.g. ["Arsenic"]) to a category string.
    Multiple metals → "metal_resistance-multimetal".
    Unknown → "metal_resistance-unknown".
    """
    if not metal_type_list:
        return "metal_resistance-unknown"
    metals = [m.lower().strip() for m in metal_type_list]
    if len(metals) > 1:
        return "metal_resistance-multimetal"
    return _METAL_TYPE_TO_CAT.get(metals[0], f"metal_resistance-{metals[0]}")


def load_methmmdb_json(json_path, hmm_base_dir):
    """
    Load MetHMMDB using its metadata JSON (the authoritative source).

    Parameters
    ----------
    json_path    : path to MetHMMDB metadata JSON
    hmm_base_dir : base directory of the MetHMMDB repo (individual/ folder will
                   be resolved relative to this, or absolute paths in hmm_file
                   are used directly)

    Returns
    -------
    models : list of model dicts (one per JSON entry)
    related_graph : dict  stem → set of related stems (from MetHMMDB itself)
    """
    import json as _json

    with open(json_path) as fh:
        entries = _json.load(fh)

    # Normalise: the JSON may be a list directly, or wrapped in a key
    if isinstance(entries, dict):
        # Try common wrapper keys
        for key in ("models", "hmms", "data", "entries"):
            if key in entries:
                entries = entries[key]
                break
        else:
            entries = list(entries.values())

    models = []
    related_graph = {}   # stem → set of related stems
    base = Path(hmm_base_dir)

    missing_files = []

    for entry in entries:
        stem = entry.get("id", "")
        if not stem:
            continue

        # Resolve HMM file path
        hmm_rel = entry.get("hmm_file", "")
        hmm_abs = None
        if hmm_rel:
            # Try absolute first, then relative to base dir
            candidate = Path(hmm_rel)
            if candidate.is_absolute() and candidate.exists():
                hmm_abs = str(candidate)
            else:
                # Strip leading path components until we find the file
                for suffix in [hmm_rel,
                                base / hmm_rel,
                                base / Path(hmm_rel).name,
                                base / "individual" / Path(hmm_rel).name]:
                    if Path(suffix).exists():
                        hmm_abs = str(suffix)
                        break

        if hmm_abs is None:
            missing_files.append(stem)
            continue

        # Parse embedded HMMER3 meta for nseq / GA / TC
        hmm_meta = parse_hmm_meta(hmm_abs)

        # Category from metal_type (authoritative)
        metal_types = entry.get("metal_type", [])
        category = _metal_type_to_category(metal_types)

        # Gene name from resistance_type (human-readable)
        resistance_type = entry.get("resistance_type", "")
        gene_name = resistance_type if resistance_type else stem

        # nseq: prefer JSON sequences_count, fallback to HMMER NSEQ
        nseq = entry.get("sequences_count", hmm_meta["nseq"] or 0)

        # Cutoff: prefer GA/TC embedded in HMM (MetHMMDB sets these properly)
        cutoff = hmm_meta.get("ga_seq") or hmm_meta.get("tc_seq") or 0.0

        models.append({
            "stem":       stem,
            "hmm_file":   hmm_abs,
            "name":       hmm_meta["name"] or stem,
            "acc":        hmm_meta["acc"],
            "desc":       resistance_type,
            "nseq":       nseq,
            "ga_seq":     hmm_meta["ga_seq"],
            "tc_seq":     hmm_meta["tc_seq"],
            "date":       hmm_meta["date"],
            "category":   category,
            "confidence": "high",
            "gene_name":  gene_name,
            "cutoff":     cutoff,
            "source":     "methmmdb",
            "_related":   set(entry.get("related_models", [])),
        })

        related_graph[stem] = set(entry.get("related_models", []))

    if missing_files:
        print(f"  [WARN] {len(missing_files)} MetHMMDB entries had no resolvable "
              f"HMM file and were skipped.", file=sys.stderr)
        if len(missing_files) <= 10:
            for s in missing_files:
                print(f"         missing: {s}", file=sys.stderr)

    return models, related_graph


# ═══════════════════════════════════════════════════════════════════════════════
# DEDUPLICATION — cross-source only
# ═══════════════════════════════════════════════════════════════════════════════

def _tokens(text):
    return set(re.findall(r"[a-zA-Z0-9]+", text.lower()))


def _jaccard(a, b):
    return len(a & b) / len(a | b) if (a | b) else 1.0


def _stem_base(name):
    """
    Normalise a HMM name to its bare gene identifier for comparison.

    Strips:
      • trailing variant suffixes:  _1, _2, _V2, _VF, _V3_corr, _rep, _corr
      • leading prefixes:           put_, Put_, All_seq_, all_seq_
      • trailing _TIGR / _Pfam / _NCBI accession cruft
    Then lowercases.

    Examples:
      fhuA_V3_corr  →  fhua
      put_FhuA      →  fhua
      All_seq_FslA  →  fsla
      arsB_1        →  arsb
      MtrB_TIGR03509 →  mtrb
    """
    s = name.strip()
    # Strip leading put_ / Put_ / All_seq_ / all_seq_
    s = re.sub(r"(?i)^(put_|all_seq_)", "", s)
    # Strip trailing _VF, _V2, _V3_corr, _V4, _rep, _corr, _TIGR\d+, _PF\d+
    s = re.sub(r"(?i)(_V\d+(_corr)?|_VF|_rep|_corr|_TIGR\d+|_PF\d+)$", "", s)
    # Strip any remaining trailing _N numeric variant
    s = re.sub(r"_\d+$", "", s)
    return s.lower()


def are_duplicates(m1, m2, thr=0.6):
    """
    Return True only when two models from DIFFERENT sources represent the same
    gene family.  Same-source models (e.g. two MetHMMDB variants) are NEVER
    considered duplicates — they are kept as distinct entries.

    Checks (in order, any True → duplicate):
      1. ACC identical  (non-empty, version-stripped)
      2. NAME identical (case-insensitive)
      3. Normalised stem identical  (strips put_, _VF, _V3_corr, _1… etc.)
      4. Jaccard(NAME tokens) ≥ thr  AND  Jaccard(DESC tokens) ≥ thr
    """
    # Same source → not a duplicate, always keep both
    if m1["source"] == m2["source"]:
        return False

    a1 = re.sub(r"\.\d+$", "", m1["acc"])
    a2 = re.sub(r"\.\d+$", "", m2["acc"])
    if a1 and a2 and a1 == a2:
        return True
    if m1["name"] and m2["name"] and m1["name"].lower() == m2["name"].lower():
        return True
    if _stem_base(m1["name"]) == _stem_base(m2["name"]) and m1["name"]:
        return True
    nj = _jaccard(_tokens(m1["name"]), _tokens(m2["name"]))
    dj = _jaccard(_tokens(m1["desc"]),  _tokens(m2["desc"]))
    if nj >= thr and dj >= thr:
        return True
    return False


def pick_best(m1, m2):
    """
    Choose which model to keep when two cross-source duplicates are detected.

    Rule: keep the more RECENT one (DATE field), because a newer build of the
    same HMM has been trained on more/better-curated sequences.
    Tie (same date or both None): keep the one with more training sequences (nseq).
    Second tie: keep the one with a non-zero cutoff.
    """
    d1 = m1.get("date")
    d2 = m2.get("date")

    if d1 and d2:
        if d1 > d2:
            return m1
        if d2 > d1:
            return m2
        # Same date → fall through to nseq

    # One has a date, other doesn't → prefer the one with a date (more metadata)
    if d1 and not d2:
        return m1
    if d2 and not d1:
        return m2

    # No dates → prefer higher nseq
    if m1["nseq"] != m2["nseq"]:
        return m1 if m1["nseq"] > m2["nseq"] else m2

    # Final tie-break: prefer the one with a cutoff set
    c1 = float(m1.get("cutoff") or 0)
    c2 = float(m2.get("cutoff") or 0)
    return m1 if c1 >= c2 else m2


def deduplicate(all_models, jaccard_thr):
    accepted, dup_log = [], []
    for candidate in all_models:
        found = False
        for i, existing in enumerate(accepted):
            if are_duplicates(candidate, existing, jaccard_thr):
                best = pick_best(candidate, existing)
                worst = existing if best is candidate else candidate
                dup_log.append({
                    "kept_stem":         best["stem"],
                    "kept_source":       best["source"],
                    "kept_cutoff":       str(best.get("cutoff", 0)),
                    "removed_stem":      worst["stem"],
                    "removed_source":    worst["source"],
                    "removed_category":  worst["category"],
                    "removed_cutoff":    str(worst.get("cutoff", 0)),
                    "cutoff_conflict":   (
                        "yes" if (best.get("cutoff") and worst.get("cutoff")
                                  and abs(float(best.get("cutoff", 0))
                                          - float(worst.get("cutoff", 0))) > 1.0)
                        else "no"
                    ),
                    "reason": (
                        "acc_match"  if best["acc"] and best["acc"] == worst["acc"]
                        else "name_match" if best["name"].lower() == worst["name"].lower()
                        else "stem_match" if _stem_base(best["name"]) == _stem_base(worst["name"])
                        else "jaccard"
                    ),
                })
                accepted[i] = best
                found = True
                break
        if not found:
            accepted.append(candidate)
    return accepted, dup_log


# ═══════════════════════════════════════════════════════════════════════════════
# REGISTRY + WRITER
# ═══════════════════════════════════════════════════════════════════════════════

REPORT_FIELDS = ["stem", "hmm_file", "name", "acc", "desc", "nseq",
                 "category", "confidence", "gene_name", "cutoff", "source"]

REG_FIELDS = ["stem", "name", "acc", "category", "gene_name", "source",
              "hmm_file", "nseq", "cutoff", "added_date", "status"]


def write_library(out_dir, accepted, dup_log, log_path=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    today = datetime.today().strftime("%Y-%m-%d")
    cutoff_lines, map_lines, registry_rows = [], [], []

    for m in accepted:
        if m["category"] == REVIEW_NEEDED:
            continue   # never write unreviewed HMMs

        cat_dir = out_dir / m["category"]
        cat_dir.mkdir(exist_ok=True)
        dest = cat_dir / f"{m['stem']}.hmm"

        src = m.get("hmm_file", "")
        if src and os.path.isfile(src):
            if str(dest.resolve()) != str(Path(src).resolve()):
                shutil.copy2(src, dest)
        else:
            print(f"  [WARN] HMM file not found for {m['stem']}: {src}",
                  file=sys.stderr)
            continue

        if m.get("cutoff"):
            cutoff_lines.append(f"{m['stem']}\t{float(m['cutoff']):.1f}")
        map_lines.append(f"{m['stem']}\t{m.get('gene_name', m['stem'])}")

        registry_rows.append({
            "stem":       m["stem"],
            "name":       m["name"],
            "acc":        m["acc"],
            "category":   m["category"],
            "gene_name":  m.get("gene_name", m["stem"]),
            "source":     m["source"],
            "hmm_file":   m.get("hmm_file", ""),
            "nseq":       str(m.get("nseq", 0)),
            "cutoff":     str(m.get("cutoff", 0)),
            "added_date": today,
            "status":     "active",
        })

    (out_dir / "HMM-bitcutoffs.txt").write_text("\n".join(cutoff_lines) + "\n")
    (out_dir / "FeGenie-map.txt").write_text("\n".join(map_lines) + "\n")

    with open(out_dir / "hmm_registry.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REG_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in sorted(registry_rows, key=lambda x: (x["category"], x["stem"])):
            w.writerow(r)

    # Write deduplication log
    if dup_log:
        dup_path = out_dir / "deduplication_log.tsv"
        with open(dup_path, "w", newline="") as fh:
            w = csv.DictWriter(fh,
                               fieldnames=["kept_stem", "kept_source",
                                           "kept_cutoff",
                                           "removed_stem", "removed_source",
                                           "removed_category", "removed_cutoff",
                                           "cutoff_conflict", "reason"],
                               delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for row in dup_log:
                w.writerow(row)
        print(f"[INFO] Deduplication log → {dup_path}")

    # Copy operon_rules.json if it exists in the repo but is not already in out_dir
    rules_src = Path(__file__).parent.parent / "hmm_library" / "operon_rules.json"
    rules_dst = out_dir / "operon_rules.json"
    if rules_src.exists() and rules_src.resolve() != rules_dst.resolve():
        shutil.copy2(rules_src, rules_dst)

    return len(registry_rows)


# ═══════════════════════════════════════════════════════════════════════════════
# VERIFY MODE
# ═══════════════════════════════════════════════════════════════════════════════

def verify_library(lib_dir):
    lib_dir = Path(lib_dir)
    reg_path = lib_dir / "hmm_registry.tsv"
    if not reg_path.exists():
        print(f"[ERROR] No hmm_registry.tsv found in {lib_dir}")
        return

    ok, missing, orphan = 0, [], []

    # Check every registered HMM exists on disk
    with open(reg_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            expected = lib_dir / row["category"] / f"{row['stem']}.hmm"
            if expected.exists():
                ok += 1
            else:
                missing.append(f"{row['category']}/{row['stem']}.hmm")

    # Check every .hmm file is in the registry
    registered = set()
    with open(reg_path) as fh:
        registered = {r["stem"] for r in csv.DictReader(fh, delimiter="\t")}
    for hmm_file in lib_dir.rglob("*.hmm"):
        if hmm_file.stem not in registered:
            orphan.append(str(hmm_file.relative_to(lib_dir)))

    print(f"[VERIFY] {ok} HMMs registered and present on disk")
    if missing:
        print(f"[VERIFY] {len(missing)} registered but MISSING from disk:")
        for f in missing[:20]:
            print(f"    {f}")
    if orphan:
        print(f"[VERIFY] {len(orphan)} .hmm files on disk but NOT in registry:")
        for f in orphan[:20]:
            print(f"    {f}")
    if not missing and not orphan:
        print("[VERIFY] Library is consistent ✓")


def _print_cutoff_stats(stats):
    total = sum(stats.values())
    print(f"[INFO] Bitscore cutoff coverage ({total} HMMs):")
    print(f"       extra --cutoff_tsv   : {stats['extra']:4d}")
    print(f"       FeGenie bitcutoffs   : {stats['fegenie_txt']:4d}  "
          f"(includes Tabuteau names already in that file)")
    print(f"       embedded GA          : {stats['embedded_ga']:4d}")
    print(f"       embedded TC          : {stats['embedded_tc']:4d}")
    print(f"       missing (→ 0)        : {stats['missing']:4d}"
          + (" ← E-value < 0.1 only" if stats["missing"] else " ✓"))


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        prog="curate_hmm_library.py",
        description="Build, inspect or verify the MetalGenie-Evo HMM library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
# Full build — MetHMMDB via JSON (authoritative, recommended)
python scripts/curate_hmm_library.py \\
    --fegenie_dir  /path/to/FeGenie/hmms/iron/ \\
    --flat_dir     /path/to/new_iron_hmms   iron_aquisition  tabuteau \\
    --flat_dir     /path/to/others          iron_aquisition  interpro \\
    --methmmdb_json /path/to/MetHMMDB/metadata.json \\
    --methmmdb_dir  /path/to/MetHMMDB/ \\
    --cutoff_tsv   tabuteau_TableS3.tsv \\
    --out_dir      hmm_library/ \\
    --log          curation_report.tsv

# MetHMMDB via individual/ folder only (keyword heuristics fallback)
python scripts/curate_hmm_library.py \\
    --fegenie_dir  /path/to/FeGenie/hmms/iron/ \\
    --methmmdb_dir /path/to/MetHMMDB/individual/ \\
    --out_dir      hmm_library/ \\
    --log          curation_report.tsv

# Fix REVIEW_NEEDED rows, then rebuild
python scripts/curate_hmm_library.py \\
    --from_report  curation_report.tsv \\
    --cutoff_tsv   tabuteau_TableS3.tsv \\
    --out_dir      hmm_library/

# Verify library consistency
python scripts/curate_hmm_library.py --verify hmm_library/
        """,
    )

    parser.add_argument("--fegenie_dir",
                        help="FeGenie hmms/iron/ directory (already categorised)")
    parser.add_argument("--flat_dir",  action="append", nargs=3,
                        metavar=("PATH", "HINT", "LABEL"),
                        help="Flat HMM dir: path  category_hint  source_label. "
                             "Repeat for multiple dirs. "
                             "Example: --flat_dir /path/new_iron_hmms iron_aquisition tabuteau")
    parser.add_argument("--methmmdb_dir",
                        help="MetHMMDB repo base directory (used with --methmmdb_json "
                             "to resolve hmm_file paths, or alone for keyword-based loading)")
    parser.add_argument("--methmmdb_json",
                        help="Path to MetHMMDB metadata JSON file "
                             "(e.g. MetHMMDB/metadata.json or models.json). "
                             "Provides authoritative metal_type → category mapping "
                             "and uses related_models for intra-database deduplication. "
                             "Requires --methmmdb_dir to resolve HMM file paths.")
    parser.add_argument("--out_dir",
                        help="Output library directory")
    parser.add_argument("--log",
                        help="Write per-HMM curation report TSV here "
                             "(edit REVIEW_NEEDED rows, then use --from_report)")
    parser.add_argument("--from_report",
                        help="Re-read a previously generated --log TSV "
                             "(after manual editing of REVIEW_NEEDED rows) "
                             "and write the final library")
    parser.add_argument("--cutoff_tsv", action="append", metavar="FILE",
                        help="TSV file with bitscore cutoffs (stem + cutoff columns). "
                             "Repeat for multiple files. "
                             "Overrides GA/TC values embedded in HMM headers. "
                             "Accepts: FeGenie-style 2-col no-header, or any TSV "
                             "with recognisable column names (stem/name/hmm, "
                             "cutoff/bitscore/threshold). "
                             "Example: --cutoff_tsv tabuteau_TableS3.tsv "
                             "--cutoff_tsv interpro_cutoffs.tsv")
    parser.add_argument("--verify", metavar="LIB_DIR",
                        help="Verify an existing library against its registry")
    parser.add_argument("--jaccard_thr", type=float, default=0.6,
                        help="Jaccard similarity threshold for deduplication")
    args = parser.parse_args()

    # ── Verify mode ───────────────────────────────────────────────────────────
    if args.verify:
        verify_library(args.verify)
        return

    # ── Auto-load FeGenie bitcutoffs (covers Tabuteau names too) ─────────────
    fegenie_cutoffs = {}
    if args.fegenie_dir:
        bcp = Path(args.fegenie_dir) / "HMM-bitcutoffs.txt"
        if bcp.exists():
            fegenie_cutoffs = load_cutoff_tsv(str(bcp))
            print(f"[INFO] Auto-loaded {len(fegenie_cutoffs)} cutoffs "
                  f"from FeGenie HMM-bitcutoffs.txt")

    # ── Optional extra cutoff files ───────────────────────────────────────────
    extra_cutoffs = {}
    for tsv_path in (args.cutoff_tsv or []):
        loaded = load_cutoff_tsv(tsv_path)
        print(f"[INFO] Loaded {len(loaded)} extra cutoffs from {tsv_path}")
        extra_cutoffs.update(loaded)

    # ── Rebuild from reviewed report ──────────────────────────────────────────
    if args.from_report:
        print(f"[INFO] Loading reviewed report from {args.from_report}…")
        all_models = []
        review_left = 0
        with open(args.from_report) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row["category"] == REVIEW_NEEDED:
                    review_left += 1
                    print(f"  [SKIP] {row['stem']}: still REVIEW_NEEDED", file=sys.stderr)
                    continue
                try:
                    cutoff = float(row.get("cutoff", 0) or 0)
                except ValueError:
                    cutoff = 0.0
                all_models.append({
                    "stem":       row["stem"],
                    "hmm_file":   row["hmm_file"],
                    "name":       row["name"],
                    "acc":        row.get("acc", ""),
                    "desc":       row.get("desc", ""),
                    "nseq":       int(row.get("nseq", 0) or 0),
                    "category":   row["category"],
                    "confidence": row.get("confidence", "low"),
                    "gene_name":  row.get("gene_name", row["stem"]),
                    "cutoff":     cutoff,
                    "source":     row.get("source", "flat_extra"),
                })
        if review_left:
            print(f"\n[WARN] {review_left} HMMs still have REVIEW_NEEDED — "
                  f"they will be skipped.", file=sys.stderr)

        accepted, dup_log = deduplicate(all_models, args.jaccard_thr)

        # In --from_report mode the cutoff column in the TSV already contains
        # the correct values from the original run. Only re-apply if the user
        # provided extra --cutoff_tsv files or --fegenie_dir.
        if fegenie_cutoffs or extra_cutoffs:
            accepted, cutoff_stats = apply_cutoffs(accepted, fegenie_cutoffs, extra_cutoffs)
            _print_cutoff_stats(cutoff_stats)
        else:
            # Report the cutoff coverage directly from the loaded TSV
            n_with = sum(1 for m in accepted if float(m.get("cutoff") or 0) > 0)
            n_without = len(accepted) - n_with
            print(f"[INFO] Bitscore cutoffs from report TSV: "
                  f"{n_with} set, {n_without} missing (→ E-value < 0.1 only)")
        n = write_library(args.out_dir, accepted, dup_log)
        print(f"[DONE] {n} HMMs written to {args.out_dir}/")
        return

    # ── Full build mode ───────────────────────────────────────────────────────
    all_models = []

    if args.fegenie_dir:
        print(f"[INFO] Loading FeGenie from {args.fegenie_dir}…")
        fe = load_fegenie(args.fegenie_dir)
        print(f"       {len(fe)} models in {len({m['category'] for m in fe})} categories")
        all_models.extend(fe)

    for dir_path, hint, label in (args.flat_dir or []):
        print(f"[INFO] Loading flat dir '{label}' from {dir_path} (hint: {hint})…")
        flat = load_flat_dir(dir_path, hint, label)
        print(f"       {len(flat)} models")
        all_models.extend(flat)

    if args.methmmdb_json:
        if not args.methmmdb_dir:
            sys.exit("[ERROR] --methmmdb_json requires --methmmdb_dir to resolve HMM paths.")
        print(f"[INFO] Loading MetHMMDB from JSON: {args.methmmdb_json}…")
        md, _ = load_methmmdb_json(args.methmmdb_json, args.methmmdb_dir)
        print(f"       {len(md)} models loaded (all variants kept)")
        all_models.extend(md)
    elif args.methmmdb_dir:
        print(f"[INFO] Loading MetHMMDB individual/ from {args.methmmdb_dir} "
              f"(keyword heuristics — consider adding --methmmdb_json)…")
        md = load_methmmdb_individual(args.methmmdb_dir)
        print(f"       {len(md)} models")
        all_models.extend(md)

    if not all_models:
        sys.exit("[ERROR] No models loaded. Provide at least one source.")

    print(f"\n[INFO] Total before deduplication: {len(all_models)}")

    # ── Deduplication ─────────────────────────────────────────────────────────
    print(f"[INFO] Deduplicating (Jaccard ≥ {args.jaccard_thr})…")
    accepted, dup_log = deduplicate(all_models, args.jaccard_thr)
    print(f"       {len(accepted)} unique  |  {len(dup_log)} duplicates removed")

    # ── Apply external cutoffs ────────────────────────────────────────────────
    accepted, cutoff_stats = apply_cutoffs(accepted, fegenie_cutoffs, extra_cutoffs)
    _print_cutoff_stats(cutoff_stats)

    # ── Category summary ──────────────────────────────────────────────────────
    cat_counts = defaultdict(int)
    review_count = 0
    for m in accepted:
        cat_counts[m["category"]] += 1
        if m["category"] == REVIEW_NEEDED:
            review_count += 1

    print(f"\n[INFO] Proposed assignment summary:")
    for cat, n in sorted(cat_counts.items()):
        flag = " ← NEEDS MANUAL REVIEW" if cat == REVIEW_NEEDED else ""
        print(f"  {n:4d}  {cat}{flag}")

    # ── Write report TSV ──────────────────────────────────────────────────────
    if args.log:
        with open(args.log, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=REPORT_FIELDS, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            for m in sorted(accepted, key=lambda x: (x["category"], x["stem"])):
                w.writerow(m)
        print(f"\n[INFO] Curation report → {args.log}")
        if review_count:
            print(f"       {review_count} rows need manual assignment.")
            print(f"       Open the TSV, fix REVIEW_NEEDED rows, save, then run:")
            print(f"       python scripts/curate_hmm_library.py \\")
            print(f"           --from_report {args.log} \\")
            print(f"           --out_dir     {args.out_dir or 'hmm_library/'}")
        else:
            print(f"       All assignments resolved. Proceeding to write library…")

    # ── Write library if no unresolved reviews ────────────────────────────────
    if not review_count and args.out_dir:
        n = write_library(args.out_dir, accepted, dup_log, args.log)
        cat_summary = defaultdict(int)
        for m in accepted:
            if m["category"] != REVIEW_NEEDED:
                cat_summary[m["category"]] += 1
        print(f"\n[INFO] Final library — {n} HMMs in {len(cat_summary)} categories:")
        for cat, n_cat in sorted(cat_summary.items()):
            print(f"  {n_cat:4d}  {cat}")
        print(f"\n[DONE] Library written to {args.out_dir}/")
    elif review_count and args.out_dir:
        print(f"\n[WARN] Library NOT written because {review_count} HMMs "
              f"still have REVIEW_NEEDED.")
        print(f"       Fix them in {args.log} and re-run with --from_report.")


if __name__ == "__main__":
    main()
