#!/usr/bin/env python3
"""
assign_flat_hmms.py  –  Propose category assignments for flat HMM directories
==============================================================================
Reads one or more flat directories of .hmm files (no sub-directory structure),
parses HMMER3 headers, and produces a tab-delimited TSV that you review and
edit before passing to build_hmm_library.py.

Sources handled:
  • new_iron_hmms  (Tabuteau et al. 2025 / iron_cheese_metagenomes)
      K* → KOfam iron genes  (biosynthesis vs import guessed from lookup table)
      put_*, All_seq_*, Fer3_*, etc. → custom-built HMMs (biosynthesis / import)
  • others          (NCBI NF* and Pfam PF* downloads from InterPro)
  • Any other flat directory

Category vocabulary (MetalGenie-Evo):
  iron_aquisition-siderophore_synthesis
  iron_aquisition-siderophore_transport
  iron_aquisition-siderophore_transport_potential
  iron_aquisition-iron_transport
  iron_aquisition-heme_transport
  iron_aquisition-heme_oxygenase
  iron_gene_regulation
  iron_oxidation
  iron_reduction
  iron_storage
  magnetosome_formation
  REVIEW_NEEDED          ← fill these in manually before running build_hmm_library.py

Usage
-----
# Inspect and produce draft TSV
python scripts/assign_flat_hmms.py \\
    --dirs  new_iron_hmms:iron_aquisition \\
            others:iron_aquisition \\
    --out   hmm_assignments.tsv

# After editing hmm_assignments.tsv, build the library
python scripts/build_hmm_library.py \\
    --fegenie_dir   /path/to/FeGenie/hmms/iron/ \\
    --flat_dirs     hmm_assignments.tsv \\
    --methmmdb_dir  MetHMMDB/individual/ \\
    --out_dir       hmm_library/

The --dirs argument accepts  path:hint  pairs where hint is a broad category
hint used to initialise guessing (iron_aquisition, iron_reduction, etc.).
"""

import argparse
import os
import re
import sys
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# KOfam ID → MetalGenie-Evo category lookup
# Based on KEGG annotations for iron-related KOs
# ─────────────────────────────────────────────────────────────────────────────

KOFAM_CAT = {
    # ── Siderophore biosynthesis ──────────────────────────────────────────────
    "K00216": "iron_aquisition-siderophore_synthesis",   # DHPS / menaquinone
    "K00257": "iron_aquisition-siderophore_synthesis",   # acsA / acyl-CoA
    "K00824": "iron_aquisition-siderophore_synthesis",   # 2,3-DHBA aminotransferase
    "K01252": "iron_aquisition-siderophore_synthesis",   # iucA (aerobactin synth)
    "K01582": "iron_aquisition-siderophore_synthesis",   # lysine decarboxylase
    "K01780": "iron_aquisition-siderophore_synthesis",
    "K01851": "iron_aquisition-siderophore_synthesis",
    "K01909": "iron_aquisition-siderophore_synthesis",   # EntD
    "K04759": "iron_aquisition-siderophore_synthesis",   # dhbF (bacillibactin)
    "K04788": "iron_aquisition-siderophore_synthesis",   # entA
    "K04789": "iron_aquisition-siderophore_synthesis",   # entB
    "K04790": "iron_aquisition-siderophore_synthesis",   # entC
    "K04791": "iron_aquisition-siderophore_synthesis",   # entD
    "K04792": "iron_aquisition-siderophore_synthesis",   # entE
    "K04793": "iron_aquisition-siderophore_synthesis",   # entF
    "K05375": "iron_aquisition-siderophore_synthesis",   # iucA (aerobactin)
    "K07116": "iron_aquisition-siderophore_synthesis",   # iucB
    "K07224": "iron_aquisition-siderophore_synthesis",   # iucC
    "K07230": "iron_aquisition-siderophore_synthesis",   # iucD
    "K10531": "iron_aquisition-siderophore_synthesis",   # dhbA
    "K11604": "iron_aquisition-siderophore_synthesis",   # vibB
    "K12237": "iron_aquisition-siderophore_synthesis",   # pvdH
    "K12238": "iron_aquisition-siderophore_synthesis",   # pvdA
    "K12239": "iron_aquisition-siderophore_synthesis",   # pvdF
    "K12240": "iron_aquisition-siderophore_synthesis",   # pvdQ
    "K12241": "iron_aquisition-siderophore_synthesis",   # pvdY
    "K12346": "iron_aquisition-siderophore_synthesis",   # viuB / vbs
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
    # ── Siderophore / iron transport (import) ─────────────────────────────────
    "K02012": "iron_aquisition-siderophore_transport",   # fepA
    "K02014": "iron_aquisition-siderophore_transport",   # fepB
    "K02016": "iron_aquisition-siderophore_transport",   # fepC
    "K02361": "iron_aquisition-siderophore_transport",   # fepD
    "K02362": "iron_aquisition-siderophore_transport",   # fepE
    "K02363": "iron_aquisition-siderophore_transport",   # fepG
    "K02364": "iron_aquisition-siderophore_transport",   # fes (enterobactin esterase)
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

# ─────────────────────────────────────────────────────────────────────────────
# NCBI Protein Family IDs → category
# NF03xxxx / NF04xxxx — known iron-related families
# ─────────────────────────────────────────────────────────────────────────────

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
    "PF16525": "iron_aquisition-siderophore_synthesis",   # SbnA-like
}

# ─────────────────────────────────────────────────────────────────────────────
# Keyword heuristics on NAME/DESC
# ─────────────────────────────────────────────────────────────────────────────

_SYNTH_KW = [
    "synthetase", "synthase", "biosynthesis", "NRPS", "nonribosomal",
    "entA", "entB", "entC", "entD", "entE", "entF",
    "pvd", "pyoverdine", "siderophore biosyn", "fusarinine", "coprogen",
    "rhizobactin", "rhizoferrin", "desferrioxamine", "ferrichrome",
    "aerobactin", "staphyloferrin", "vibrioferrin",
    "iucA", "iucB", "iucC", "iucD",
    "dhbA", "dhbB", "dhbC", "dhbE", "dhbF",
    "FslA", "FslC", "NRPS", "adenylation domain",
    "condensation domain", "thiolation", "isochorismate",
]

_TRANSPORT_KW = [
    "receptor", "TonB-dependent", "outer membrane receptor",
    "ABC transporter", "binding protein", "permease",
    "import", "uptake", "siderophore transport",
    "fepA", "fepB", "fepC", "fepD", "fepG",
    "fhuA", "fhuB", "fhuC", "fhuD", "fhuE",
    "fecA", "fecB", "fecC", "fecD",
    "iutA", "irgA", "cirA", "fiu",
    "FutA", "FutB", "FutC",
    "viuA", "fatA", "pvdR",
]

_HEME_KW  = ["heme", "haemin", "haemophore", "hemophore", "HasA", "HusA",
             "HmbR", "ShuA", "BhuA"]
_REDOX_KW = ["reductase", "ferritin", "bacterioferritin"]
_REG_KW   = ["fur", "regulator", "repressor", "response regulator", "sigma"]


def _kw_match(text, keywords):
    t = text.lower()
    return any(k.lower() in t for k in keywords)


def guess_category(stem, name, desc, hint="iron_aquisition"):
    """
    Guess MetalGenie-Evo category for a single HMM.
    Returns (category, confidence) where confidence is 'high'/'low'.
    """
    combined = f"{stem} {name} {desc}"

    # ── KOfam K-number lookup ─────────────────────────────────────────────────
    ko_match = re.match(r"^(K\d{5})", stem)
    if ko_match:
        ko = ko_match.group(1)
        if ko in KOFAM_CAT:
            return KOFAM_CAT[ko], "high"

    # ── NCBI NF-number lookup ─────────────────────────────────────────────────
    nf_match = re.match(r"^(NF\d+)", stem)
    if nf_match:
        nf = nf_match.group(1)
        if nf in NCBI_NF_CAT:
            return NCBI_NF_CAT[nf], "high"

    # ── Pfam PF lookup ────────────────────────────────────────────────────────
    pf_match = re.match(r"^(PF\d+)", stem)
    if pf_match:
        pf = pf_match.group(1)
        if pf in PFAM_CAT:
            return PFAM_CAT[pf], "high"

    # ── Keyword heuristics ────────────────────────────────────────────────────
    if _kw_match(combined, _SYNTH_KW):
        return "iron_aquisition-siderophore_synthesis", "low"
    if _kw_match(combined, _HEME_KW):
        return "iron_aquisition-heme_transport", "low"
    if _kw_match(combined, _TRANSPORT_KW):
        return "iron_aquisition-siderophore_transport", "low"
    if _kw_match(combined, _REDOX_KW):
        return "iron_storage", "low"
    if _kw_match(combined, _REG_KW):
        return "iron_gene_regulation", "low"

    # Tabuteau custom names (put_*, All_seq_*, Fer3_*, etc.)
    lower = stem.lower()
    if any(x in lower for x in ["fusarinine", "tafc", "sid", "sidc",
                                  "nrps", "rhizoferrin", "copfam",
                                  "fer3", "rhbf", "desd", "fso1",
                                  "rumb", "rhizobactin"]):
        return "iron_aquisition-siderophore_synthesis", "low"
    if any(x in lower for x in ["fcua", "fhua", "fitA", "fiua", "fut",
                                  "fhue", "lbtu", "pvdr", "iron"]):
        return "iron_aquisition-siderophore_transport", "low"

    return "REVIEW_NEEDED", "low"


# ─────────────────────────────────────────────────────────────────────────────
# HMM header parser (minimal — NAME + DESC only)
# ─────────────────────────────────────────────────────────────────────────────

def parse_hmm_name_desc(hmm_path):
    name, desc = "", ""
    with open(hmm_path, errors="replace") as fh:
        for line in fh:
            tag = line.split()[0] if line.split() else ""
            if tag == "NAME":
                name = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else ""
            elif tag == "DESC":
                desc = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else ""
            elif tag == "//":
                break
    return name, desc


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Propose category assignments for flat HMM directories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--dirs", nargs="+", required=True,
        metavar="PATH:HINT",
        help="One or more  path:hint  pairs. "
             "HINT is a broad initial category hint "
             "(iron_aquisition | iron_reduction | metal_resistance | …). "
             "Example:  --dirs  new_iron_hmms:iron_aquisition  others:iron_aquisition"
    )
    parser.add_argument(
        "--out", required=True,
        help="Output TSV file to review and edit"
    )
    args = parser.parse_args()

    rows = []
    n_high = 0
    n_review = 0

    for dir_hint in args.dirs:
        if ":" in dir_hint:
            dir_path, hint = dir_hint.rsplit(":", 1)
        else:
            dir_path, hint = dir_hint, "iron_aquisition"

        hmm_files = sorted(Path(dir_path).glob("*.hmm"))
        if not hmm_files:
            print(f"  [WARN] No .hmm files found in {dir_path}", file=sys.stderr)
            continue

        print(f"[INFO] {dir_path}: {len(hmm_files)} HMM files")

        for hmm_file in hmm_files:
            # Skip HMMER index files
            if hmm_file.suffix != ".hmm":
                continue
            stem = hmm_file.stem
            name, desc = parse_hmm_name_desc(str(hmm_file))
            cat, conf = guess_category(stem, name, desc, hint)

            if conf == "high":
                n_high += 1
            if cat == "REVIEW_NEEDED":
                n_review += 1

            rows.append({
                "stem":      stem,
                "hmm_file":  str(hmm_file.resolve()),
                "name":      name or stem,
                "desc":      desc,
                "category":  cat,
                "confidence": conf,
                "gene_name": re.sub(r"_\d+$", "", stem),
                "cutoff":    "",
                "notes":     "",
            })

    # Write TSV
    header = ["stem", "hmm_file", "name", "desc", "category",
              "confidence", "gene_name", "cutoff", "notes"]

    with open(args.out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[k]) for k in header) + "\n")

    # Summary
    cats = {}
    for r in rows:
        cats[r["category"]] = cats.get(r["category"], 0) + 1

    print(f"\n[INFO] Total HMMs processed: {len(rows)}")
    print(f"       High-confidence assignments: {n_high}")
    print(f"       Need manual review (REVIEW_NEEDED): {n_review}")
    print(f"\n[INFO] Proposed categories:")
    for cat, n in sorted(cats.items()):
        flag = " ← REVIEW THESE" if cat == "REVIEW_NEEDED" else ""
        print(f"  {n:4d}  {cat}{flag}")

    print(f"\n[INFO] Draft TSV written to: {args.out}")
    print(     "       Open it in a spreadsheet, fix REVIEW_NEEDED rows,")
    print(     "       then pass it to build_hmm_library.py via --flat_dirs")


if __name__ == "__main__":
    main()
