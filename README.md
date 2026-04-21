# MetalGenie-Evo

**MetalGenie-Evo** is built on top of [FeGenie](https://github.com/Arkadiy-Garber/FeGenie), an established tool for identifying iron-cycling genes in genome and metagenome assemblies. MetalGenie-Evo extends FeGenie's biological logic with additional features — coordinate-based operon clustering, parallel execution, and configurable filtering rules — and broadens its scope to include metal mobility resistance genes through integration of **MetHMMDB** and additional curated HMM sources.

## Citations

If you use MetalGenie-Evo in your research, please cite all of the following that apply:

**FeGenie** (iron cycling HMMs and operon logic):
> Garber AI, Nealson KH, Okamoto A, McAllister SM, Chan CS, Barco RA, Merino N (2020) *FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies.* Front. Microbiol. 11:37. [doi:10.3389/fmicb.2020.00037](https://www.frontiersin.org/articles/10.3389/fmicb.2020.00037/full)

**MetHMMDB** (metal mobility resistance gene HMMs):
> Kciuchcinski K et al. (2025) *Fast and accurate detection of metal resistance genes using MetHMMDB.* bioRxiv. [doi:10.1101/2024.12.26.629440](https://www.biorxiv.org/content/10.1101/2024.12.26.629440v2)

**Tabuteau et al.** (additional iron acquisition HMMs for bacteria and fungi):
> Tabuteau S, Hervé V, Irlinger F, Monnet C (2025) *Metagenomic profiling and genome-centric analysis reveal iron acquisition systems in cheese-associated bacteria and fungi.* Environmental Microbiology 27(12):e70218. [doi:10.1111/1462-2920.70218](https://doi.org/10.1111/1462-2920.70218)

**KOfam** (KEGG Ortholog HMMs used in Tabuteau et al.):
> Aramaki T, Blanc-Mathieu R, Endo H, Ohkubo K, Kanehisa M, Goto S, Ogata H (2020) *KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold.* Bioinformatics 36:2251–2252. [doi:10.1093/bioinformatics/btz859](https://doi.org/10.1093/bioinformatics/btz859)

**NCBI Protein Family Models** (NF* HMMs):
> Li W, O'Neill KR, Haft DH, DiCuccio M, Chetvernin V, Badretdin A, Coulouris G, Chitsaz F, Derbyshire MK, Durkin AS, Gonzales NR, Gwadz M, Lanczycki CJ, Song JS, Thanki N, Wang J, Yamashita RA, Yang M, Zheng C, Marchler-Bauer A, Thibaud-Nissen F (2021) *RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation.* Nucleic Acids Research 49:D1020–D1028. [doi:10.1093/nar/gkaa1105](https://doi.org/10.1093/nar/gkaa1105)

---

## Contents

- [What's new compared to FeGenie](#whats-new-compared-to-fegenie)
- [How it works](#how-it-works)
- [Installation](#installation)
- [HMM library](#hmm-library)
- [Usage](#usage)
- [Output files](#output-files)
- [Operon rules](#operon-rules)
- [Extending the HMM library](#extending-the-hmm-library)
- [Differences from FeGenie at a glance](#differences-from-fegenie-at-a-glance)

---

## What's new compared to FeGenie

### 1 — Coordinate-based operon clustering

FeGenie clusters genes into putative operons by comparing the numerical suffix of Prodigal's ORF names (e.g. `scaffold_1_5` and `scaffold_1_9` differ by 4, within the default gap of 5). MetalGenie-Evo optionally reads Prodigal's **GFF output** (`--gff_dir`) and clusters by actual base-pair distance between gene endpoints (`--max_bp_gap`, default 5 000 bp), which is more robust across different gene callers and assembly formats. The original ordinal-index method is kept as a fallback when GFF files are not available, ensuring full backward compatibility.

> **Note for reproducibility:** FeGenie's published results use the ordinal-index method with `--max_gap 5`. To replicate those results exactly, run MetalGenie-Evo without `--gff_dir` and with `--max_gap 5`.

### 2 — Strand-aware clustering

With `--strand_aware`, MetalGenie-Evo separates hits on opposite DNA strands into distinct clusters before applying operon rules. This option reduces false-positive operon calls for divergently transcribed gene pairs, and is particularly relevant for magnetosome islands and iron-reduction outer-membrane complexes.

### 3 — Parallel execution

MetalGenie-Evo dispatches all `(genome, HMM)` pairs to a `ProcessPoolExecutor`, using `--threads` CPUs in total. Already-computed `.tblout` files are cached and reused on subsequent runs, so interrupted analyses can be resumed without re-running hmmsearch.

### 4 — Configurable operon rules

MetalGenie-Evo externalises FeGenie's operon-context filters into a **JSON file** (`hmm_library/operon_rules.json`) that travels with the HMM library. All original FeGenie rules are encoded exactly. New databases can declare their own rules — or declare `report_all_categories` to bypass filtering entirely — without modifying the source code.

### 5 — Expanded HMM library

The bundled library integrates four sources (see [Citations](#citations)):
- **FeGenie** — original iron cycling profiles
- **Tabuteau et al. 2025** — additional iron acquisition HMMs for bacteria and fungi, sourced from KOfam, FeGenie, and NCBI Protein Family Models
- **NCBI NF\* / Pfam** — selected profiles for iron transport and siderophore systems
- **MetHMMDB** — 254 profiles for 121 metal mobility resistance genes

All sources are merged, deduplicated, and tracked in a versioned registry (`hmm_registry.tsv`). The curation process is fully reproducible via `scripts/curate_hmm_library.py`.

### 6 — Modular codebase

MetalGenie-Evo separates concerns into focused, testable functions and uses `pathlib` and `subprocess` throughout. Output CSV formats are identical to FeGenie's, so existing R plotting scripts (`DotPlot.R`, `dendro-heatmap.R`) work without modification.

### 7 — Metagenome-specific features

When working with metagenomic assemblies rather than isolated genomes, MetalGenie-Evo provides:

- **Integrated Prodigal** (`--fna_dir --meta`) — runs ORF prediction internally in metagenomic mode, producing `.faa` and `.gff` automatically
- **Contig length filter** (`--min_contig_len`) — skips ORFs on contigs shorter than a given threshold, reducing noise from incomplete assemblies
- **Relaxed operon thresholds** (`--relaxed_operons`) — halves minimum gene-count requirements for clusters on short contigs where operons may be fragmented across assembly breaks
- **TPM-normalised coverage** (`--norm_coverage`) — normalises the coverage heatmap to TPM using contig lengths, enabling cross-sample comparisons
- **Contig column in all outputs** — every output file includes the source contig, enabling downstream linkage to binning or taxonomic classification results
- **Long-format output** (`MetalGenie-Evo-results-long.tsv`) — tidy TSV with one row per ORF, no cluster separators, no sequences; directly loadable in R (`read_tsv()`) or pandas

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
  └── Index mode  →  Prodigal ordinal gap (fallback, FeGenie-compatible)
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

MetalGenie-Evo is installed by cloning this repository and creating the conda environment. The `environment.yml` handles everything — it installs all external dependencies (HMMER, Prodigal, samtools) **and** registers the `MetalGenie-Evo` command via pip in one step.

### Dependencies

**Core pipeline (required):**

| Tool | Min version | Purpose |
|---|---|---|
| Python | 3.8 | Runtime |
| HMMER (`hmmsearch`) | 3.3 | HMM search |
| Prodigal | 2.6.3 | ORF prediction from nucleotide input |

**Coverage heatmap (optional, only needed for `--bam`):**

| Tool | Min version | Purpose |
|---|---|---|
| samtools | 1.10 | Per-contig coverage from BAM files |

> If you already have pre-computed depth files from MetaBAT2 (`jgi_summarize_bam_contig_depths`) or BBMap (`pileup.sh`), use `--depth` instead of `--bam` — samtools is not needed.

---

### Option 1 — Conda / Mamba (recommended)

```bash
# 1. Clone the repository
git clone https://github.com/your-username/MetalGenie-Evo.git
cd MetalGenie-Evo

# 2. Create the conda environment — installs all dependencies
#    AND registers the MetalGenie-Evo command automatically
conda env create -f environment.yml
# or faster:
mamba env create -f environment.yml

# 3. Activate
conda activate metalgenie-evo

# 4. Verify
MetalGenie-Evo --help
```

The environment is named `metalgenie-evo` automatically. No `chmod`, no PATH editing, no `setup.sh` required.

---

### Option 2 — Add to an existing conda environment

If you already have an active environment with HMMER and Prodigal:

```bash
git clone https://github.com/your-username/MetalGenie-Evo.git
cd MetalGenie-Evo

# Install any missing tools
conda install -c bioconda hmmer>=3.3 prodigal>=2.6.3 samtools>=1.10

# Register the MetalGenie-Evo command
pip install -e .

MetalGenie-Evo --help
```

---

### Option 3 — Manual (no conda)

```bash
git clone https://github.com/your-username/MetalGenie-Evo.git
cd MetalGenie-Evo

# Verify dependencies are in PATH
python3 --version    # must be >= 3.8
hmmsearch -h         # must be >= 3.3
prodigal -v          # must be >= 2.6.3
samtools --version   # must be >= 1.10  (optional)

# Register the command
pip install -e .

MetalGenie-Evo --help
```


## HMM library

The `hmm_library/` directory is **included in this repository** and is ready to use out of the box. It contains:

- **FeGenie iron HMMs** — original profiles for iron cycling genes (iron reduction, iron oxidation, iron acquisition, iron storage, magnetosome formation, iron gene regulation)
- **Tabuteau et al. HMMs** — additional iron acquisition profiles from KOfam, FeGenie, and NCBI NF* sources, curated for bacteria and fungi
- **MetHMMDB metal resistance HMMs** — profiles for metal mobility resistance genes (arsenic, copper, mercury, zinc/cobalt/cadmium, chromium, and more)
- `HMM-bitcutoffs.txt` — per-HMM bitscore thresholds
- `MetalGenie-map.txt` — HMM stem → readable gene name mapping
- `hmm_registry.tsv` — full provenance record (source, NSEQ, cutoff, dates)
- `deduplication_log.tsv` — log of all cross-source duplicate decisions
- `operon_rules.json` — configurable operon-context filtering rules

No setup step is required. Clone the repo and run.

---

## Usage

### Basic run (ORF FASTA input)

If you already have Prodigal `.faa` files:

```bash
MetalGenie-Evo \
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
MetalGenie-Evo \
    --faa_dir    orfs/ \
    --gff_dir    gffs/ \
    --hmm_dir    hmm_library/ \
    --out        results/ \
    --threads    16 \
    --max_bp_gap 5000
```

### Strand-aware clustering

```bash
MetalGenie-Evo \
    --faa_dir      orfs/ \
    --gff_dir      gffs/ \
    --hmm_dir      hmm_library/ \
    --out          results/ \
    --threads      16 \
    --strand_aware
```

### Normalised heatmap output

```bash
MetalGenie-Evo \
    --faa_dir orfs/ \
    --hmm_dir hmm_library/ \
    --out     results/ \
    --norm                   # normalise counts by total ORFs × 1000
```

### Report all hits (skip operon filtering)

```bash
MetalGenie-Evo \
    --faa_dir    orfs/ \
    --hmm_dir    hmm_library/ \
    --out        results/ \
    --all_results
```


### Coverage-based heatmap (BAM files)

When you have BAM files from read mapping, MetalGenie-Evo can produce an additional heatmap where values represent read coverage instead of gene counts.

**Single BAM file** (same BAM for all genomes/bins):
```bash
MetalGenie-Evo \
    --faa_dir orfs/ \
    --hmm_dir hmm_library/ \
    --out     results/ \
    --bam     mapping/sample.sorted.bam
```

**Per-genome BAM files** (different BAM per genome):
```bash
# Prepare a TSV: genome_label<TAB>bam_path
echo -e "bin_001.faa\tmapping/bin_001.sorted.bam" > bam_map.tsv
echo -e "bin_002.faa\tmapping/bin_002.sorted.bam" >> bam_map.tsv

MetalGenie-Evo \
    --faa_dir orfs/ \
    --hmm_dir hmm_library/ \
    --out     results/ \
    --bams    bam_map.tsv
```

**Pre-computed depth files** (from MetaBAT2, BBMap, or samtools):
```bash
# If you already ran jgi_summarize_bam_contig_depths or BBMap pileup.sh:
MetalGenie-Evo \
    --faa_dir orfs/ \
    --hmm_dir hmm_library/ \
    --out     results/ \
    --depth   mapping/contigs.depth       # single file

# Or per-genome:
MetalGenie-Evo \
    --faa_dir  orfs/ \
    --hmm_dir  hmm_library/ \
    --out      results/ \
    --depths   depth_map.tsv              # genome_label<TAB>depth_path
```

Accepted depth file formats: `jgi_summarize_bam_contig_depths`, `samtools coverage -H`, BBMap `pileup.sh`, or a plain two-column `contig<TAB>depth` file.

The coverage heatmap is written as `MetalGenie-Evo-coverage-heatmap.csv` alongside the standard gene-count heatmap.

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
| `--bam` | — | Single sorted BAM file for coverage heatmap (requires samtools ≥ 1.10) |
| `--bams` | — | TSV file: `genome_label<TAB>bam_path` for per-genome BAM files |
| `--depth` | — | Pre-computed depth file (jgi / BBMap / samtools format) |
| `--depths` | — | TSV file: `genome_label<TAB>depth_path` for per-genome depth files |
| `--norm_coverage` | off | Normalise coverage heatmap to TPM |
| **Metagenome options** | | |
| `--fna_dir` | — | Directory of nucleotide assemblies (mutually exclusive with `--faa_dir`); Prodigal run internally |
| `--fna_ext` | `fna` | Extension of assembly files |
| `--meta` | off | Prodigal metagenomic mode (`-p meta`) |
| `--min_contig_len` | `0` | Skip ORFs on contigs shorter than this (bp); 0 = no filter |
| `--relaxed_operons` | off | Halve operon min-gene thresholds for short contigs |
| `--relaxed_threshold` | `10000` | Contig length (bp) below which thresholds are relaxed |

---

## Output files

All output files are written to `--out/`.

### `MetalGenie-Evo-summary.csv`

One row per reported ORF. Cluster separators (`#,#,…`) mark operon boundaries.

| Column | Description |
|---|---|
| `category` | Functional category (e.g. `iron_reduction`) |
| `genome/assembly` | Source genome filename |
| `contig` | Source contig identifier |
| `orf` | ORF identifier |
| `gene` | Readable gene name (from `MetalGenie-map.txt`) |
| `bitscore` | hmmsearch bitscore |
| `bitscore_cutoff` | Per-HMM cutoff used |
| `cluster_id` | Genomic cluster index |
| `heme_c_motifs` | Count of CXXCH heme-binding motifs |
| `protein_sequence` | Amino acid sequence |

### `MetalGenie-Evo-geneSummary-clusters.csv`

Compact version compatible with FeGenie's R plotting scripts (`DotPlot.R`, `dendro-heatmap.R`). Includes `contig` column.

### `MetalGenie-Evo-heatmap-data.csv`

Gene-count matrix: rows = functional categories, columns = genomes.
With `--norm`: values are `(gene_count / total_orfs) × 1000`.

### `MetalGenie-Evo-results-long.tsv`

Tidy long-format TSV — one row per ORF, no cluster separators, no protein sequences. Includes `contig` and `contig_len` columns. Directly loadable in R or pandas for comparative analyses.

### `MetalGenie-Evo-coverage-heatmap.csv`

Coverage-based matrix (only written when `--bam`, `--bams`, `--depth`, or `--depths` is provided).
With `--norm_coverage`: values are TPM-normalised using contig lengths.

---

## Operon rules

Operon rules are defined in `hmm_library/operon_rules.json`. If the file is absent, MetalGenie-Evo falls back to the built-in FeGenie defaults.

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

To add new HMM sources, use `scripts/curate_hmm_library.py`. This script handles category assignment, cross-source deduplication, cutoff loading, and provenance tracking. A full curation report (`curation_report.tsv`) and deduplication log are written alongside the library for transparency and reproducibility.

```bash
# Full rebuild from all sources
python scripts/curate_hmm_library.py \
    --fegenie_dir   /path/to/FeGenie/hmms/iron/ \
    --flat_dir      /path/to/new_iron_hmms  iron_aquisition  tabuteau \
    --flat_dir      /path/to/others         iron_aquisition  interpro \
    --methmmdb_json /path/to/MetHMMDB/metadata.json \
    --methmmdb_dir  /path/to/MetHMMDB/ \
    --out_dir       hmm_library/ \
    --log           curation_report.tsv

# Verify library consistency
python scripts/curate_hmm_library.py --verify hmm_library/
```

To add individual HMMs manually without rebuilding:

1. Place `.hmm` files in a subdirectory of `hmm_library/` named after the functional category.
2. Add bitscore cutoffs to `hmm_library/HMM-bitcutoffs.txt` (`stem<tab>cutoff`).
3. Add readable gene names to `hmm_library/MetalGenie-map.txt` (`stem<tab>gene_name`).
4. Optionally add operon rules to `hmm_library/operon_rules.json`.

---

## Differences from FeGenie at a glance

| Feature | FeGenie | MetalGenie-Evo |
|---|---|---|
| Operon clustering | ORF ordinal index | **bp coordinates (GFF)** + ordinal fallback |
| Strand awareness | No | **Yes (`--strand_aware`)** |
| Parallelism | Sequential | **ProcessPoolExecutor** |
| Operon rules | Hardcoded Python | **JSON config file** |
| Result caching | No | **Yes (tblout cache)** |
| MetHMMDB support | No | **Yes** |
| Additional iron acquisition HMMs | No | **Yes (Tabuteau et al. 2025)** |
| HMM deduplication | No | **Yes (ACC + NAME + Jaccard, cross-source only)** |
| Versioned HMM registry | No | **Yes (hmm_registry.tsv)** |
| Normalised heatmap | `--norm` | **`--norm`** (same) |
| Coverage heatmap | `--bam` | **`--bam` / `--bams` / `--depth` / `--depths`** |
| Coverage normalisation | No | **TPM (`--norm_coverage`)** |
| Integrated Prodigal | No | **Yes (`--fna_dir --meta`)** |
| Contig length filter | No | **Yes (`--min_contig_len`)** |
| Relaxed operon thresholds | No | **Yes (`--relaxed_operons`)** |
| Contig column in outputs | No | **Yes** |
| Long-format output | No | **Yes (`-results-long.tsv`)** |
| R script compatibility | Yes | **Yes (same CSV format)** |
| FeGenie filter rules | All | **All (exact port)** |

---

## License

MIT License. See [LICENSE](LICENSE).

MetalGenie-Evo is not affiliated with the original FeGenie authors.
