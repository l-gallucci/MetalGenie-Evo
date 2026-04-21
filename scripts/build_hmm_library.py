#!/usr/bin/env python3
"""
build_hmm_library.py  –  Build / update a MetalGenie-Evo HMM library
==================================================================
Merges HMM profiles from multiple sources into the MetalGenie-Evo
library format (one sub-directory per functional category) and handles:

  • FeGenie iron/ directory (original FeGenie HMM set)
  • MetHMMDB  – either as a concatenated .hmm file OR from the
    individual/ subfolder (one .hmm per gene family)
  • Any extra .hmm file you supply with name:category[:gene_name] syntax

Deduplication
─────────────
Two profiles are considered duplicates when:
  1. Their ACC fields are identical (non-empty), or
  2. Their NAME fields are identical, or
  3. NAME Jaccard similarity AND DESC Jaccard similarity both exceed
     --jaccard_thr (default 0.6), or
  4. Their NAME stems (trailing _N variant suffix stripped) match.

Ties are resolved by NSEQ count, then by source priority
(FeGenie > MetHMMDB individual > MetHMMDB concatenated > extra).

Registry (hmm_registry.tsv)
───────────────────────────
Every accepted model is tracked with:
  stem, name, acc, category, gene_name, library, source_file,
  nseq, cutoff, added_date, updated_date, status

Re-running with --update marks removed entries as status=removed instead
of deleting them, enabling reproducible version-aware pipelines.

Usage examples
──────────────
# First build: FeGenie + MetHMMDB (individual folder)
python scripts/build_hmm_library.py \\
    --fegenie_dir  /path/to/FeGenie/hmms/iron/ \\
    --methmmdb_dir /path/to/MetHMMDB/individual/ \\
    --out_dir      hmm_library/

# First build: FeGenie + MetHMMDB (concatenated file)
python scripts/build_hmm_library.py \\
    --fegenie_dir  /path/to/FeGenie/hmms/iron/ \\
    --methmmdb     /path/to/MetHMMDB/MetHMMDb.hmm \\
    --out_dir      hmm_library/

# Update to a newer MetHMMDB release (keeps FeGenie models untouched)
python scripts/build_hmm_library.py \\
    --methmmdb_dir /path/to/MetHMMDB_v2/individual/ \\
    --out_dir      hmm_library/ \\
    --update
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


# ─────────────────────────────────────────────────────────────────────────────
# HMMER3 .hmm file parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_hmm_headers(hmm_path):
    """
    Parse a (possibly multi-model) HMMER3 .hmm file.
    Returns list of dicts – one per model – with header metadata
    and byte offsets for later extraction.
    """
    models = []
    current = None
    byte_pos = 0

    with open(hmm_path, "rb") as fh:
        for raw_line in fh:
            line = raw_line.decode("utf-8", errors="replace").rstrip()
            byte_pos += len(raw_line)

            if line.startswith("HMMER3/f"):
                if current is not None:
                    models.append(current)
                current = {
                    "name": "", "acc": "", "desc": "", "nseq": 0,
                    "ga_seq": None, "tc_seq": None,
                    "start_byte": byte_pos - len(raw_line),
                    "end_byte": None,
                    "source": str(Path(hmm_path).name),
                }
            elif current is not None:
                tag = line.split()[0] if line.split() else ""
                if tag == "NAME":
                    current["name"] = line.split(None, 1)[1].strip()
                elif tag == "ACC":
                    current["acc"]  = line.split(None, 1)[1].strip()
                elif tag == "DESC":
                    current["desc"] = (line.split(None, 1)[1].strip()
                                       if len(line.split(None, 1)) > 1 else "")
                elif tag == "NSEQ":
                    try:
                        current["nseq"] = int(line.split()[1])
                    except (IndexError, ValueError):
                        pass
                elif tag == "GA":
                    parts = line.split()
                    try:
                        current["ga_seq"] = float(parts[1].rstrip(";"))
                    except (IndexError, ValueError):
                        pass
                elif tag == "TC":
                    parts = line.split()
                    try:
                        current["tc_seq"] = float(parts[1].rstrip(";"))
                    except (IndexError, ValueError):
                        pass
                elif tag == "//":
                    current["end_byte"] = byte_pos
                    models.append(current)
                    current = None

    if current is not None:
        current["end_byte"] = byte_pos
        models.append(current)

    return models


def extract_hmm_bytes(hmm_path, start_byte, end_byte):
    """Extract raw bytes of one model from a multi-model .hmm file."""
    with open(hmm_path, "rb") as fh:
        fh.seek(start_byte)
        return fh.read(end_byte - start_byte)


# ─────────────────────────────────────────────────────────────────────────────
# Deduplication helpers
# ─────────────────────────────────────────────────────────────────────────────

def _tokenize(text):
    return set(re.findall(r"[a-zA-Z0-9]+", text.lower()))


def _jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def _name_stem(name):
    """Strip trailing variant suffix  _N  from a HMM name."""
    return re.sub(r"_\d+$", "", name).lower()


def models_are_duplicate(m1, m2, jaccard_thr=0.6):
    acc1 = re.sub(r"\.\d+$", "", m1["acc"])
    acc2 = re.sub(r"\.\d+$", "", m2["acc"])
    if acc1 and acc2 and acc1 == acc2:
        return True
    if m1["name"] and m2["name"] and m1["name"].lower() == m2["name"].lower():
        return True
    if _name_stem(m1["name"]) == _name_stem(m2["name"]) \
            and m1["name"] and m2["name"]:
        return True
    nj = _jaccard(_tokenize(m1["name"]), _tokenize(m2["name"]))
    dj = _jaccard(_tokenize(m1["desc"]),  _tokenize(m2["desc"]))
    if nj >= jaccard_thr and dj >= jaccard_thr:
        return True
    return False


_SRC_PRIORITY = {"fegenie": 0, "methmmdb_individual": 1,
                 "methmmdb_concat": 2, "extra": 3}


def pick_best(m1, m2):
    p1 = _SRC_PRIORITY.get(m1.get("library", "extra"), 3)
    p2 = _SRC_PRIORITY.get(m2.get("library", "extra"), 3)
    if p1 != p2:
        return m1 if p1 < p2 else m2
    return m1 if m1["nseq"] >= m2["nseq"] else m2


# ─────────────────────────────────────────────────────────────────────────────
# MetHMMDB category inference
# ─────────────────────────────────────────────────────────────────────────────

_METAL_KW = {
    "arsenic": "metal_resistance-arsenic",
    "arsb": "metal_resistance-arsenic", "arsc": "metal_resistance-arsenic",
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


def infer_methmmdb_category(model, meta_table=None):
    stem = _name_stem(model["name"])
    if meta_table and stem in meta_table:
        row = meta_table[stem]
        return row["category"], row.get("gene_name", stem)
    desc = model["desc"].lower()
    for kw, cat in _METAL_KW.items():
        if kw in desc or kw in model["name"].lower():
            return cat, re.sub(r"_\d+$", "", model["name"])
    return f"metal_resistance-{stem.split('_')[0]}", \
           re.sub(r"_\d+$", "", model["name"])


def load_meta_tsv(path):
    table = {}
    if not path or not os.path.isfile(path):
        return table
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            stem = re.sub(r"_\d+$", "", row.get("name", "")).strip().lower()
            if stem:
                table[stem] = {
                    "category":  row.get("category", "metal_resistance-unknown"),
                    "gene_name": row.get("gene_name", stem),
                }
    return table


# ─────────────────────────────────────────────────────────────────────────────
# Loaders
# ─────────────────────────────────────────────────────────────────────────────

def load_fegenie(fegenie_dir):
    base = Path(fegenie_dir)
    cutoffs, gene_map = {}, {}

    cutoffs_path = base / "HMM-bitcutoffs.txt"
    if cutoffs_path.exists():
        for line in open(cutoffs_path):
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                try:
                    cutoffs[ls[0]] = float(ls[1])
                except ValueError:
                    pass

    map_path = base / "FeGenie-map.txt"
    if map_path.exists():
        for line in open(map_path):
            ls = line.rstrip().split("\t")
            if len(ls) >= 2:
                gene_map[ls[0]] = ls[1]

    models = []
    for cat_dir in sorted(base.iterdir()):
        if not cat_dir.is_dir() or cat_dir.name.startswith("."):
            continue
        for hmm_file in sorted(cat_dir.glob("*.hmm")):
            parsed = parse_hmm_headers(str(hmm_file))
            m = parsed[0] if parsed else {}
            stem = hmm_file.stem
            models.append({
                "name":      m.get("name", stem),
                "acc":       m.get("acc", ""),
                "desc":      m.get("desc", ""),
                "nseq":      m.get("nseq", 0),
                "ga_seq":    m.get("ga_seq"),
                "tc_seq":    m.get("tc_seq"),
                "category":  cat_dir.name,
                "library":   "fegenie",
                "hmm_file":  str(hmm_file),
                "gene_name": gene_map.get(stem, stem),
                "stem":      stem,
                "cutoff":    cutoffs.get(stem, 0.0),
            })
    return models


def load_methmmdb_individual(individual_dir, meta_table=None):
    """Load MetHMMDB from its individual/ folder (one .hmm per gene)."""
    models = []
    for hmm_file in sorted(Path(individual_dir).glob("*.hmm")):
        parsed = parse_hmm_headers(str(hmm_file))
        m = parsed[0] if parsed else {}
        if not m.get("name"):
            m["name"] = hmm_file.stem
        category, gene_name = infer_methmmdb_category(m, meta_table)
        stem = hmm_file.stem
        models.append({
            "name":      m["name"],
            "acc":       m.get("acc", ""),
            "desc":      m.get("desc", ""),
            "nseq":      m.get("nseq", 0),
            "ga_seq":    m.get("ga_seq"),
            "tc_seq":    m.get("tc_seq"),
            "category":  category,
            "library":   "methmmdb_individual",
            "hmm_file":  str(hmm_file),
            "gene_name": gene_name,
            "stem":      stem,
            "cutoff":    m.get("ga_seq") or m.get("tc_seq") or 0.0,
        })
    return models


def load_methmmdb_concat(concat_hmm, meta_table=None):
    """Load MetHMMDB from a single concatenated .hmm file."""
    models = []
    for m in parse_hmm_headers(concat_hmm):
        if not m.get("name"):
            continue
        category, gene_name = infer_methmmdb_category(m, meta_table)
        full_stem = m["name"].replace(" ", "_")
        models.append({
            "name":        m["name"],
            "acc":         m.get("acc", ""),
            "desc":        m.get("desc", ""),
            "nseq":        m.get("nseq", 0),
            "ga_seq":      m.get("ga_seq"),
            "tc_seq":      m.get("tc_seq"),
            "category":    category,
            "library":     "methmmdb_concat",
            "hmm_file":    None,
            "hmm_bytes":   extract_hmm_bytes(concat_hmm,
                               m["start_byte"], m["end_byte"]),
            "gene_name":   gene_name,
            "stem":        full_stem,
            "cutoff":      m.get("ga_seq") or m.get("tc_seq") or 0.0,
        })
    return models


# ─────────────────────────────────────────────────────────────────────────────
# Registry
# ─────────────────────────────────────────────────────────────────────────────

_REG_FIELDS = ["stem", "name", "acc", "category", "gene_name",
               "library", "source_file", "nseq", "cutoff",
               "added_date", "updated_date", "status"]


def load_registry(path):
    reg = {}
    if not os.path.isfile(path):
        return reg
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            reg[row["stem"]] = row
    return reg


def save_registry(path, reg):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_REG_FIELDS,
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in sorted(reg.values(),
                          key=lambda r: (r["category"], r["stem"])):
            w.writerow(row)


# ─────────────────────────────────────────────────────────────────────────────
# Library writer
# ─────────────────────────────────────────────────────────────────────────────

def write_library(out_dir, accepted, registry):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cutoff_lines, map_lines = [], []

    for m in accepted:
        cat_dir = out_dir / m["category"]
        cat_dir.mkdir(exist_ok=True)
        out_hmm = cat_dir / f"{m['stem']}.hmm"

        if m.get("hmm_file") and os.path.isfile(m["hmm_file"]):
            if str(out_hmm.resolve()) != str(Path(m["hmm_file"]).resolve()):
                shutil.copy2(m["hmm_file"], out_hmm)
        elif m.get("hmm_bytes"):
            with open(out_hmm, "wb") as fh:
                fh.write(m["hmm_bytes"])
        else:
            print(f"  [WARN] No data for {m['stem']}, skipping", file=sys.stderr)
            continue

        if m.get("cutoff"):
            cutoff_lines.append(f"{m['stem']}\t{float(m['cutoff']):.1f}")
        map_lines.append(f"{m['stem']}\t{m.get('gene_name', m['stem'])}")

    (out_dir / "HMM-bitcutoffs.txt").write_text("\n".join(cutoff_lines) + "\n")
    (out_dir / "FeGenie-map.txt").write_text("\n".join(map_lines) + "\n")
    save_registry(out_dir / "hmm_registry.tsv", registry)

    print(f"[INFO] Wrote {len(accepted)} models to {out_dir}/")


def load_flat_tsv(tsv_path):
    """
    Load a reviewed assignment TSV (output of assign_flat_hmms.py).
    Skips rows with category == REVIEW_NEEDED or empty category.
    """
    models = []
    skipped = 0
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cat = row.get("category", "").strip()
            if not cat or cat == "REVIEW_NEEDED":
                skipped += 1
                continue
            stem = row.get("stem", "").strip()
            hmm_file = row.get("hmm_file", "").strip()
            if not hmm_file or not os.path.isfile(hmm_file):
                print(f"  [WARN] File not found for {stem}: {hmm_file}", file=sys.stderr)
                skipped += 1
                continue
            try:
                cutoff = float(row.get("cutoff", 0) or 0)
            except ValueError:
                cutoff = 0.0
            models.append({
                "name":      row.get("name", stem),
                "acc":       "",
                "desc":      row.get("desc", ""),
                "nseq":      0,
                "ga_seq":    None, "tc_seq": None,
                "category":  cat,
                "library":   "flat_extra",
                "hmm_file":  hmm_file,
                "gene_name": row.get("gene_name", stem),
                "stem":      stem,
                "cutoff":    cutoff,
            })
    if skipped:
        print(f"  [WARN] {skipped} rows skipped (REVIEW_NEEDED or missing file)")
    return models


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Build/update a MetalGenie-Evo HMM library",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fegenie_dir",    help="FeGenie hmms/iron/ directory")
    p.add_argument("--methmmdb_dir",
                   help="MetHMMDB individual/ folder (preferred over --methmmdb)")
    p.add_argument("--methmmdb",       help="MetHMMDB concatenated .hmm file")
    p.add_argument("--methmmdb_meta",
                   help="MetHMMDB metadata TSV (columns: name, category, gene_name)")
    p.add_argument("--flat_dirs",
                   help="TSV produced by assign_flat_hmms.py (after manual review). "
                        "Columns: stem, hmm_file, name, desc, category, confidence, "
                        "gene_name, cutoff, notes")
    p.add_argument("--out_dir",        required=True, help="Output library directory")
    p.add_argument("--update",         action="store_true",
                   help="Update mode: add/replace models, keep existing ones")
    p.add_argument("--jaccard_thr",    type=float, default=0.6,
                   help="Jaccard threshold for NAME+DESC deduplication")
    args = p.parse_args()

    out_dir  = Path(args.out_dir)
    reg_path = out_dir / "hmm_registry.tsv"
    today    = datetime.today().strftime("%Y-%m-%d")

    registry = load_registry(str(reg_path)) if args.update else {}
    if args.update:
        print(f"[INFO] Update mode – registry has {len(registry)} entries")

    meta_table = load_meta_tsv(args.methmmdb_meta)
    all_models = []

    if args.fegenie_dir:
        print(f"[INFO] Loading FeGenie from {args.fegenie_dir}…")
        fe = load_fegenie(args.fegenie_dir)
        print(f"       {len(fe)} models")
        all_models.extend(fe)

    if args.methmmdb_dir:
        print(f"[INFO] Loading MetHMMDB individual/ from {args.methmmdb_dir}…")
        md = load_methmmdb_individual(args.methmmdb_dir, meta_table)
        print(f"       {len(md)} models")
        all_models.extend(md)
    elif args.methmmdb:
        print(f"[INFO] Loading MetHMMDB (concat) from {args.methmmdb}…")
        mc = load_methmmdb_concat(args.methmmdb, meta_table)
        print(f"       {len(mc)} models")
        all_models.extend(mc)

    if args.flat_dirs:
        print(f"[INFO] Loading reviewed flat HMMs from {args.flat_dirs}…")
        flat = load_flat_tsv(args.flat_dirs)
        print(f"       {len(flat)} models")
        all_models.extend(flat)

    if not all_models:
        sys.exit("[ERROR] No models loaded.")

    print(f"\n[INFO] Total before deduplication: {len(all_models)}")
    print(f"[INFO] Deduplicating (Jaccard ≥ {args.jaccard_thr})…")

    accepted, duplicates = [], []
    for candidate in all_models:
        found_dup = False
        for i, existing in enumerate(accepted):
            if models_are_duplicate(candidate, existing, args.jaccard_thr):
                best = pick_best(candidate, existing)
                worst_stem = existing["stem"] if best is candidate else candidate["stem"]
                duplicates.append((best["stem"], worst_stem))
                if best is candidate:
                    accepted[i] = candidate
                found_dup = True
                break
        if not found_dup:
            accepted.append(candidate)

    print(f"       {len(accepted)} unique  |  {len(duplicates)} duplicates removed")
    if duplicates[:15]:
        print("  Sample duplicates (kept ← removed):")
        for k, r in duplicates[:15]:
            print(f"    {k:40s} ← {r}")
        if len(duplicates) > 15:
            print(f"    … and {len(duplicates)-15} more")

    # Update registry
    for m in accepted:
        stem = m["stem"]
        if stem in registry:
            registry[stem].update({
                "updated_date": today,
                "category":     m["category"],
                "cutoff":       str(m.get("cutoff", 0)),
                "gene_name":    m.get("gene_name", stem),
                "status":       "active",
            })
        else:
            registry[stem] = {
                "stem":         stem,
                "name":         m["name"],
                "acc":          m["acc"],
                "category":     m["category"],
                "gene_name":    m.get("gene_name", m["name"]),
                "library":      m.get("library", "unknown"),
                "source_file":  m.get("hmm_file", ""),
                "nseq":         str(m.get("nseq", 0)),
                "cutoff":       str(m.get("cutoff", 0)),
                "added_date":   today,
                "updated_date": today,
                "status":       "active",
            }

    if args.update:
        active = {m["stem"] for m in accepted}
        for stem in registry:
            if stem not in active:
                registry[stem]["status"] = "removed"

    write_library(str(out_dir), accepted, registry)

    cat_counts = defaultdict(int)
    for m in accepted:
        cat_counts[m["category"]] += 1
    print("\n[INFO] Models per category:")
    for cat, n in sorted(cat_counts.items()):
        print(f"  {n:4d}  {cat}")
    print(f"\n[DONE] Library ready at {out_dir}/")


if __name__ == "__main__":
    main()
