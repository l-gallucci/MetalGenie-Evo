# MetalGenie-Evo

**MetalGenie-Evo** is a modernised reimplementation of [FeGenie](https://github.com/Arkadiy-Garber/FeGenie) — a tool for identifying iron-cycling genes in genome and metagenome assemblies. MetalGenie-Evo preserves all of FeGenie's biological logic while replacing its ageing internals with a faster, more flexible, and extensible architecture. It also introduces native support for **MetHMMDB** (metal mobility resistance genes), enabling simultaneous profiling of iron cycling and metal resistance in a single run.

> **Citation for the original tool:**
> Garber AI et al. (2020) *FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies.* Front. Microbiol. 11:37. [doi:10.3389/fmicb.2020.00037](https://www.frontiersin.org/articles/10.3389/fmicb.2020.00037/full)

---

## Contents

- [What changed from FeGenie](#what-changed-from-fegenie)
- [How it works](#how-it-works)
- [Installation](#installation)
- [Usage](#usage)
- [Output files](#output-files)
- [Operon rules](#operon-rules)
- [Extending the HMM library](#extending-the-hmm-library)
- [Differences from FeGenie at a glance](#differences-from-fegenie-at-a-glance)

---

## What changed from FeGenie

### 1 — Coordinate-based operon clustering (new)

FeGenie determines whether two HMM hits belong to the same operon by comparing the **numerical suffix** of Prodigal's ORF names (e.g. `scaffold_1_5` and `scaffold_1_9` differ by 4, so they cluster together with the default gap of 5). This works only when ORF ordinals reflect genomic proximity, which is true for Prodigal but fragile for other callers and completely wrong if ORFs have been renamed.

MetalGenie-Evo reads Prodigal's **GFF output** (`--gff_dir`) and clusters by actual base-pair distance between gene endpoints (`--max_bp_gap`, default 5 000 bp). The ordinal-index fallback is kept for cases where GFF files are unavailable.

### 2 — Strand-aware clustering (new)

With `--strand_aware`, MetalGenie-Evo separates hits on opposite DNA strands into distinct clusters before applying operon rules. Divergently transcribed genes rarely form a functional operon, so this reduces false-positive operon calls — especially relevant for magnetosome islands and iron-reduction outer-membrane complexes.

### 3 — Parallel execution (new)

FeGenie calls `hmmsearch` sequentially, one HMM × one genome at a time. MetalGenie-Evo dispatches all `(genome, HMM)` pairs to a `ProcessPoolExecutor`, using `--threads` CPUs in total. Already-computed `.tblout` files are cached and reused, so interrupted runs can be resumed cheaply.

### 4 — Configurable operon rules (new)

FeGenie's operon-context filters are scattered as hardcoded `if/elif` branches across ~600 lines of Python. MetalGenie-Evo externalises these rules into a **JSON file** (`hmm_library/operon_rules.json`) that travels with the HMM library. New databases (like MetHMMDB) can declare their own rules — or declare `report_all_categories` to bypass filtering entirely — without touching the code.

### 5 — MetHMMDB integration (new)

[MetHMMDB](https://github.com/Haelmorn/MetHMMDB) provides 254 HMM profiles for 121 metal mobility resistance genes. `scripts/build_hmm_library.py` merges MetHMMDB into the MetalGenie-Evo library, handling:
- Parsing of MetHMMDB's individual `.hmm` files **or** the concatenated file
- Automatic functional-category assignment from gene names and descriptions
- Cross-database deduplication (NAME / ACC exact match + Jaccard similarity on tokens)
- A versioned registry (`hmm_registry.tsv`) for reproducible updates

### 6 — Cleaner codebase

The original FeGenie is a single 3 100-line file. MetalGenie-Evo separates concerns into focused functions, uses `pathlib` and `subprocess` instead of `os.system`, and produces the same output CSV format so existing R plotting scripts work without modification.

---

## How it works

```
Input FAA files (Prodigal ORFs)
        │
        ▼
hmmsearch (parallel, per-HMM bitscore cutoffs from HMM-bitcutoffs.txt)
        │
        ▼
Best-hit collection  ─────  one hit per ORF, highest bitscore wins
        │
        ▼
Genomic clustering
  ├── GFF mode  →  bp-distance + optional strand separation
  └── Index mode  →  Prodigal ordinal gap (fallback)
        │
        ▼
Operon-context filtering  (operon_rules.json)
  ├── Per-rule checks  (FLEET, MAM, Fox, Mtr, DFE, siderophore …)
  └── report_all bypass for metal_resistance-* categories
        │
        ▼
Second-pass per-gene filters
  ├── Cyc1: co-occurrence with ≥2 Fe-redox genes
  ├── Cyc2: length ≥ 365 aa + CXXCH heme motif
  └── iron_gene_regulation: ≥1 regulation gene in cluster
        │
        ▼
Heme-motif counting  (CXXCH, CX3CH, CX4CH, CX14CH, CX15CH)
        │
        ▼
Output CSVs
```

All FeGenie operon rules are preserved exactly:

| System | Rule |
|---|---|
| FLEET (iron oxidation e⁻ shuttle) | ≥ 5 unique FLEET genes in cluster |
| Magnetosome (MAM) | ≥ 5 unique MAM genes in cluster |
| FoxABC | ≥ 2 of 3 subunits |
| FoxEYZ | FoxE (anchor) must be present |
| DFE operons | ≥ 3 of 4–5 subunits |
| Mtr/Mto disambiguation | MtoA+MtrB → oxidation; MtrA+MtrB → reduction |
| Siderophore transport | ≥ 2 distinct transport HMMs, or lone trusted receptor |
| Siderophore synthesis | ≥ 3 distinct synthesis HMMs |
| Iron/heme transport | ≥ 2 distinct transport HMMs |
| Cyc1 | Co-occurs with ≥ 2 Fe-redox genes |
| Cyc2 | Length ≥ 365 aa + ≥ 1 CXXCH motif |

---

## Installation

### With conda (recommended)

```bash
git clone https://github.com/your-username/MetalGenie-Evo.git
cd MetalGenie-Evo
bash setup.sh
conda activate metalgenie-evo
```

### Manual

Ensure the following tools are in your `$PATH`:

| Tool | Version | Purpose |
|---|---|---|
| Python | ≥ 3.8 | Core runtime |
| HMMER (`hmmsearch`) | ≥ 3.3 | HMM search |
| Prodigal | ≥ 2.6.3 | ORF prediction (if providing nucleotide input) |

```bash
git clone https://github.com/your-username/MetalGenie-Evo.git
cd MetalGenie-Evo
chmod +x MetalGenie-Evo.py scripts/build_hmm_library.py
```

---

## HMM library

The `hmm_library/` directory is **included in this repository** and is ready to use out of the box. It contains:

- **FeGenie iron HMMs** — all original profiles for iron cycling genes (iron reduction, iron oxidation, iron acquisition, iron storage, magnetosome formation, iron gene regulation)
- **MetHMMDB metal resistance HMMs** — 254 profiles for 121 metal mobility resistance genes (arsenic, copper, mercury, zinc/cobalt/cadmium, chromium, and more)
- `HMM-bitcutoffs.txt` — per-HMM bitscore thresholds
- `FeGenie-map.txt` — HMM stem → readable gene name mapping
- `hmm_registry.tsv` — full provenance record (source, NSEQ, cutoff, dates)
- `operon_rules.json` — configurable operon-context filtering rules

No setup step is required. Clone the repo and run.

## Usage

### Basic run (ORF FASTA input)

If you already have Prodigal `.faa` files:

```bash
MetalGenie-Evo.py \
    --faa_dir  orfs/ \
    --hmm_dir  hmm_library/ \
    --out      results/
```

### With GFF-based clustering (recommended)

Prodigal produces both `.faa` and `.gff` output. Pass the GFF directory to enable coordinate-based operon clustering:

```bash
# Run Prodigal (metagenome mode)
mkdir -p orfs/ gffs/
for f in genomes/*.fna; do
    base=$(basename "$f" .fna)
    prodigal -i "$f" -p meta \
             -a "orfs/${base}.faa" \
             -f gff -o "gffs/${base}.gff" -q
done

# Run MetalGenie-Evo
MetalGenie-Evo.py \
    --faa_dir    orfs/ \
    --gff_dir    gffs/ \
    --hmm_dir    hmm_library/ \
    --out        results/ \
    --threads    16 \
    --max_bp_gap 5000
```

### Strand-aware clustering

```bash
MetalGenie-Evo.py \
    --faa_dir      orfs/ \
    --gff_dir      gffs/ \
    --hmm_dir      hmm_library/ \
    --out          results/ \
    --threads      16 \
    --strand_aware
```

### Normalised heatmap output

```bash
MetalGenie-Evo.py \
    --faa_dir orfs/ \
    --hmm_dir hmm_library/ \
    --out     results/ \
    --norm                   # normalise counts by total ORFs × 1000
```

### Report all hits (skip operon filtering)

```bash
MetalGenie-Evo.py \
    --faa_dir    orfs/ \
    --hmm_dir    hmm_library/ \
    --out        results/ \
    --all_results
```

### Full argument reference

| Argument | Default | Description |
|---|---|---|
| `--faa_dir` | *(required)* | Directory of ORF `.faa` files |
| `--faa_ext` | `faa` | Extension of ORF files (without dot) |
| `--gff_dir` | — | Directory of Prodigal `.gff` files |
| `--hmm_dir` | *(required)* | HMM library directory |
| `--out` | `metalgenie_evo_out` | Output directory |
| `--threads` | `4` | Total CPU threads |
| `--hmm_threads` | `1` | CPUs per individual hmmsearch call |
| `--max_gap` | `5` | Max ORF-index gap (index-mode clustering) |
| `--max_bp_gap` | `5000` | Max bp gap between gene ends (GFF mode) |
| `--strand_aware` | off | Split clusters at strand changes (GFF only) |
| `--all_results` | off | Skip all operon-context filters |
| `--norm` | off | Normalise heatmap counts per total ORFs |
| `--keep_tblout` | off | Keep cached `.tblout` files after run |

---

## Output files

All output files are written to `--out/`.

### `MetalGenie-Evo-summary.csv`

One row per reported ORF. Cluster separators (`#,#,…`) mark operon boundaries.

| Column | Description |
|---|---|
| `category` | Functional category (e.g. `iron_reduction`) |
| `genome/assembly` | Source genome filename |
| `orf` | ORF identifier |
| `gene` | Readable gene name (from `FeGenie-map.txt`) |
| `bitscore` | hmmsearch bitscore |
| `bitscore_cutoff` | Per-HMM cutoff used |
| `cluster_id` | Genomic cluster index |
| `heme_c_motifs` | Count of CXXCH heme-binding motifs |
| `protein_sequence` | Amino acid sequence |

### `MetalGenie-Evo-geneSummary-clusters.csv`

Compact version compatible with FeGenie's R plotting scripts (`DotPlot.R`, `dendro-heatmap.R`).

### `MetalGenie-Evo-heatmap-data.csv`

Gene-count matrix: rows = functional categories, columns = genomes.
With `--norm`: values are `(gene_count / total_orfs) × 1000`.

---

## Operon rules

Operon rules are defined in `hmm_library/operon_rules.json`. Place a copy of this file in your `--hmm_dir` to customise filtering. If the file is absent, MetalGenie-Evo falls back to the built-in FeGenie defaults.

### Rule schema

```json
{
  "report_all_categories": ["metal_resistance-*", "iron_storage"],

  "rules": [
    {
      "name":       "FLEET",
      "categories": ["iron_oxidation"],
      "genes":      ["EetA", "EetB", "Ndh2", "FmnB", "FmnA", "DmkA", "DmkB", "PplA"],
      "rule":       "require_n_of",
      "min_genes":  5,
      "on_fail":    "passthrough_non_members"
    }
  ]
}
```

### `report_all_categories`

Glob patterns. Any cluster whose **entire** category set matches one of these patterns bypasses all rules and is reported as-is. MetHMMDB categories (`metal_resistance-*`) and `iron_storage` are included by default.

### Rule types

| `rule` | Behaviour |
|---|---|
| `require_n_of` | Cluster must contain ≥ `min_genes` unique members of `genes` |
| `require_anchor` | A specific `anchor` gene must be present |
| `require_n_cat` | Cluster must have ≥ `min_genes` distinct HMMs whose category is in `categories` |
| `require_n_cat_or_lone_trusted` | As above, but lone hits in `trusted_lone` are also accepted |
| `mtr_disambiguation` | Re-assigns category based on Mtr/Mto subunit co-presence |

### `on_fail` actions

| `on_fail` | What happens when a rule threshold is not met |
|---|---|
| `passthrough_non_members` | Keep only genes **not** in the rule's gene set |
| `drop` | Remove the entire cluster |
| `keep_all` | Keep everything (used by disambiguation rules) |

### Adding rules for a new database

```json
{
  "name":       "MY_OPERON",
  "categories": ["my_category"],
  "genes":      ["GeneA", "GeneB", "GeneC"],
  "rule":       "require_n_of",
  "min_genes":  2,
  "on_fail":    "drop"
}
```

Or to report every hit unconditionally:

```json
{
  "report_all_categories": ["metal_resistance-*", "iron_storage", "my_category"]
}
```

---

## Extending the HMM library

The bundled library covers iron cycling (FeGenie) and metal resistance (MetHMMDB). To add more gene families:

1. Place individual `.hmm` files in a new subdirectory of `hmm_library/` named after the functional category.
2. Add per-HMM bitscore cutoffs to `hmm_library/HMM-bitcutoffs.txt` (tab-delimited: `stem<tab>cutoff`). Omit the line to rely on E-value < 0.1 only.
3. Add readable gene names to `hmm_library/FeGenie-map.txt` (`stem<tab>gene_name`).
4. Optionally add operon rules for the new categories to `hmm_library/operon_rules.json`.

To rebuild the library from scratch or update MetHMMDB to a newer release, use `scripts/build_hmm_library.py`:

```bash
# Add a new MetHMMDB release without touching FeGenie models
python scripts/build_hmm_library.py \
    --methmmdb_dir MetHMMDB_v2/individual/ \
    --out_dir      hmm_library/ \
    --update
```

---

## Differences from FeGenie at a glance

| Feature | FeGenie | MetalGenie-Evo |
|---|---|---|
| Operon clustering | ORF ordinal index | **bp coordinates (GFF)** + ordinal fallback |
| Strand awareness | No | **Yes (`--strand_aware`)** |
| Parallelism | Sequential | **ProcessPoolExecutor** |
| Operon rules | Hardcoded Python | **JSON config file** |
| Result caching | No | **Yes (tblout cache)** |
| MetHMMDB support | No | **Yes (individual/ or concat)** |
| HMM deduplication | No | **Yes (ACC + NAME + Jaccard)** |
| Versioned HMM registry | No | **Yes (hmm_registry.tsv)** |
| Normalised heatmap | `--norm` | **`--norm`** (same) |
| R script compatibility | Yes | **Yes (same CSV format)** |
| FeGenie filter rules | All | **All (exact port)** |

---

## License

MIT License. See [LICENSE](LICENSE).

MetalGenie-Evo is not affiliated with the original FeGenie authors.
If you use MetalGenie-Evo in your research, please also cite the original FeGenie paper (see top of this README) and MetHMMDB if you use its profiles.
