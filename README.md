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
> Li W et al. (2021) *RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation.* Nucleic Acids Research 49:D1020–D1028. [doi:10.1093/nar/gkaa1105](https://doi.org/10.1093/nar/gkaa1105)

**UniOP** (operon prediction, if using `--operon_prediction`):
> Su H, Zhang R, Söding J (2024) *UniOP: a universal operon prediction for high-throughput prokaryotic (meta-)genomic data using intergenic distance.* bioRxiv. [doi:10.1101/2024.11.11.623000](https://www.biorxiv.org/content/10.1101/2024.11.11.623000)

---

## Contents

- [What's new compared to FeGenie](#whats-new-compared-to-fegenie)
- [How it works](#how-it-works)
- [Installation](#installation)
- [HMM library](#hmm-library)
- [Usage](#usage)
- [Operon prediction with UniOP](#operon-prediction-with-uniop)
- [Anvi'o integration](#anvio-integration)
- [Output files](#output-files)
- [Operon filtering logic](#operon-filtering-logic)
- [Extending the HMM library](#extending-the-hmm-library)
- [Visualisation](#visualisation)
- [Differences from FeGenie at a glance](#differences-from-fegenie-at-a-glance)

---

## What's new compared to FeGenie

### 1 — Coordinate-based operon clustering

FeGenie clusters genes by comparing the numerical suffix of Prodigal's ORF names. MetalGenie-Evo optionally reads Prodigal's GFF output (`--gff_dir`) and clusters by actual base-pair distance between gene endpoints (`--max_bp_gap`, default 5 000 bp). The original ordinal-index method is kept as a fallback.

> **Note for reproducibility:** FeGenie's published results use the ordinal-index method with `--max_gap 5`. To replicate exactly, run without `--gff_dir` and with `--max_gap 5`.

### 2 — Strand-aware clustering

With `--strand_aware`, MetalGenie-Evo separates hits on opposite DNA strands into distinct clusters before applying operon rules. Particularly relevant for magnetosome islands and iron-reduction outer-membrane complexes.

### 3 — Parallel execution

All `(genome, HMM)` pairs are dispatched to a `ProcessPoolExecutor`. Already-computed `.tblout` files are cached and reused on subsequent runs.

### 4 — FeGenie-exact operon filtering + configurable rules for model organisms

By default, MetalGenie-Evo uses an exact reimplementation of FeGenie's operon-context filtering logic, including the per-ORF pass/break behaviour for iron acquisition categories. Users working with well-characterised model organisms can override this with a custom `operon_rules.json` in the HMM library directory. See [Operon filtering logic](#operon-filtering-logic).

### 5 — Expanded HMM library

Integrates four sources: FeGenie, Tabuteau et al. 2025 (KOfam + custom), NCBI NF*/Pfam, and MetHMMDB. All sources are merged, deduplicated, and tracked in a versioned registry.

### 6 — Metagenome-specific features

- **Integrated Prodigal** (`--fna_dir --meta`) — ORF prediction internally in metagenomic mode
- **Contig length filter** (`--min_contig_len`) — skips ORFs on short contigs
- **Relaxed operon thresholds** (`--relaxed_operons`) — halves thresholds for short contigs
- **TPM-normalised coverage** (`--norm_coverage`)
- **Contig column** in all outputs
- **Long-format output** (`MetalGenie-Evo-results-long.tsv`)

### 7 — UniOP operon prediction (optional)

The `--operon_prediction` flag runs [UniOP](https://github.com/hongsua/UniOP) on each genome/MAG and produces `MetalGenie-Evo-OperonStructure.tsv`, linking each HMM hit to its predicted operon. UniOP uses only intergenic distances and is fully applicable to MAGs and novel organisms. The operon prediction is an independent, additional output — it does not filter or modify HMM hits.

### 8 — Anvi'o integration (optional)

The `--anvio` flag produces `MetalGenie-Evo-anvio-functions.tsv`, directly importable with `anvi-import-functions`. When used with `--bakta_gff_dir`, MetalGenie-Evo maps Prodigal ORF names to Bakta gene IDs via coordinate matching, enabling seamless integration with Anvi'o databases built from Bakta external gene calls.

---

## How it works

```
Input FAA files (Prodigal ORFs or --fna_dir for internal prediction)
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
Operon-context filtering
  ├── FeGenie exact port (default)  — per-ORF pass/break logic
  ├── JSON rule engine (if operon_rules.json present)
  └── --all_results  — report everything unfiltered
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
Output CSVs  +  optional UniOP operon structure  +  optional Anvi'o TSV
```

FeGenie operon rules implemented exactly:

| System | Rule |
|---|---|
| FLEET | ≥ 5 unique FLEET genes in cluster |
| Magnetosome (MAM) | ≥ 5 unique MAM genes in cluster |
| FoxABC | ≥ 2 of 3 subunits |
| FoxEYZ | FoxE anchor must be present |
| DFE operons | ≥ 3 of 4–5 subunits |
| Mtr/Mto disambiguation | MtoA+MtrB → oxidation; MtrA+MtrB → reduction |
| Siderophore transport | ≥ 2 distinct transport HMMs, or lone trusted receptor |
| Siderophore synthesis | ≥ 3 distinct synthesis HMMs |
| Iron/heme transport | ≥ 2 distinct transport HMMs |
| Cyc1 | Co-occurs with ≥ 2 Fe-redox genes |
| Cyc2 | Length ≥ 365 aa + ≥ 1 CXXCH motif |

---

## Installation

MetalGenie-Evo is installed by cloning this repository and creating the conda environment. The `environment.yml` installs all external dependencies and registers the `MetalGenie-Evo` command via pip in one step.

### Dependencies

**Core pipeline (required):**

| Tool | Min version | Purpose |
|---|---|---|
| Python | 3.8 | Runtime |
| HMMER (`hmmsearch`) | 3.3 | HMM search |
| Prodigal | 2.6.3 | ORF prediction from nucleotide input |

**Coverage heatmap (optional):**

| Tool | Min version | Purpose |
|---|---|---|
| samtools | 1.10 | Per-contig coverage from BAM files |

**Operon prediction (optional):**

| Tool | Notes |
|---|---|
| UniOP | Not on conda/pip — install separately (see below) |

> If you have pre-computed depth files from MetaBAT2 or BBMap, use `--depth` — samtools is not needed.

---

### Option 1 — Conda / Mamba (recommended)

```bash
git clone https://github.com/l-gallucci/MetalGenie-Evo.git
cd MetalGenie-Evo

conda env create -f environment.yml   # or: mamba env create -f environment.yml
conda activate metalgenie-evo

pip install -e .

MetalGenie-Evo --help
```

No `chmod`, no PATH editing, no `setup.sh` required.

---

### Option 2 — Add to an existing conda environment

```bash
git clone https://github.com/l-gallucci/MetalGenie-Evo.git
cd MetalGenie-Evo
conda install -c bioconda hmmer>=3.3 prodigal>=2.6.3 samtools>=1.10
pip install -e .
MetalGenie-Evo --help
```

---

### Option 3 — Manual (no conda)

```bash
git clone https://github.com/l-gallucci/MetalGenie-Evo.git
cd MetalGenie-Evo
pip install -e .
MetalGenie-Evo --help
```

---

### Installing UniOP (optional)

UniOP is required only for `--operon_prediction`. It is a standalone Python script, not distributed through conda or pip.

```bash
# Clone UniOP
git clone https://github.com/hongsua/UniOP.git

# Copy Prodigal into UniOP's working directory
# (required only when passing FNA files; not needed with Prodigal-format FAA headers)
cp $(which prodigal) UniOP/src/

# Install Python dependencies
pip install pandas numpy scikit-learn
```

Pass the path when running MetalGenie-Evo:
```bash
MetalGenie-Evo ... --operon_prediction --uniop_path /path/to/UniOP/src/UniOP
```

---

## HMM library

The `hmm_library/` directory is included in this repository and is ready to use. It contains:

- **FeGenie iron HMMs** — original profiles for iron cycling genes
- **Tabuteau et al. HMMs** — additional iron acquisition profiles from KOfam, FeGenie, and NCBI NF* sources
- **MetHMMDB metal resistance HMMs** — profiles for metal mobility resistance genes
- `HMM-bitcutoffs.txt` — per-HMM bitscore thresholds
- `MetalGenie-map.txt` — HMM stem → readable gene name mapping
- `hmm_registry.tsv` — full provenance record
- `deduplication_log.tsv` — log of all cross-source duplicate decisions
- `operon_rules.json` — configurable filtering rules (for model organisms; absent = FeGenie default)

---

## Usage

### Basic run (Prodigal FAA input)

```bash
MetalGenie-Evo \
    --faa_dir  orfs/ \
    --hmm_dir  hmm_library/ \
    --out      results/
```

### From nucleotide assemblies (Prodigal internal)

```bash
MetalGenie-Evo \
    --fna_dir  assemblies/ \
    --fna_ext  fna \
    --meta \
    --hmm_dir  hmm_library/ \
    --out      results/ \
    --threads  16
```

### With Bakta-annotated genomes

```bash
MetalGenie-Evo \
    --faa_dir  bakta_output/ \     # Bakta .faa files
    --gff_dir  bakta_output/ \     # Bakta .gff3 files (same directory is fine)
    --hmm_dir  hmm_library/ \
    --out      results/ \
    --threads  16
```

### Report all hits (skip operon filtering)

```bash
MetalGenie-Evo --faa_dir orfs/ --hmm_dir hmm_library/ --out results/ --all_results
```

### Coverage-based heatmap

```bash
# From BAM files
MetalGenie-Evo --faa_dir orfs/ --hmm_dir hmm_library/ --out results/ \
    --bams bam_map.tsv --norm_coverage

# From pre-computed depth files
MetalGenie-Evo --faa_dir orfs/ --hmm_dir hmm_library/ --out results/ \
    --depths depth_map.tsv
```

Where `bam_map.tsv` / `depth_map.tsv` is a two-column TSV: `genome_label<TAB>file_path`.

### Full argument reference

| Argument | Default | Description |
|---|---|---|
| `--faa_dir` | *(required\*)* | Directory of ORF `.faa` files |
| `--fna_dir` | *(required\*)* | Directory of nucleotide assemblies; Prodigal run internally |
| `--faa_ext` | `faa` | Extension of ORF files |
| `--fna_ext` | `fna` | Extension of assembly files |
| `--meta` | off | Prodigal metagenomic mode (`-p meta`) |
| `--gff_dir` | — | Prodigal/Bakta GFF files for bp-based clustering |
| `--hmm_dir` | *(required)* | HMM library directory |
| `--out` | `metalgenie_evo_out` | Output directory |
| `--threads` | `4` | Total CPU threads |
| `--hmm_threads` | `1` | CPUs per individual hmmsearch call |
| `--max_gap` | `5` | Max ORF-index gap (index-mode, FeGenie-compatible) |
| `--max_bp_gap` | `5000` | Max bp gap between gene ends (GFF mode) |
| `--strand_aware` | off | Split clusters at strand changes (GFF mode only) |
| `--all_results` | off | Skip all operon-context filters |
| `--norm` | off | Normalise gene-count heatmap per total ORFs × 1000 |
| `--keep_tblout` | off | Keep cached `.tblout` files |
| `--min_contig_len` | `0` | Skip ORFs on contigs shorter than this (bp) |
| `--relaxed_operons` | off | Halve operon min-gene thresholds for short contigs |
| `--relaxed_threshold` | `10000` | Contig length (bp) below which thresholds are relaxed |
| `--bam` | — | Single sorted BAM file (requires samtools ≥ 1.10) |
| `--bams` | — | TSV: `genome<TAB>bam_path` |
| `--depth` | — | Pre-computed depth file (jgi / BBMap / samtools / plain) |
| `--depths` | — | TSV: `genome<TAB>depth_path` |
| `--norm_coverage` | off | TPM-normalise coverage heatmap |
| `--operon_prediction` | off | Run UniOP and write `OperonStructure.tsv` |
| `--uniop_path` | `uniop` | Path to UniOP script |
| `--bakta_gff_dir` | — | Bakta GFF3 directory for Prodigal↔Bakta ID mapping |
| `--anvio` | off | Write Anvi'o-compatible functions TSV |

\* `--faa_dir` and `--fna_dir` are mutually exclusive; exactly one is required.

---

## Operon prediction with UniOP

UniOP predicts operons from intergenic distances alone — no RNA-seq, no functional annotations, no reference genomes required. It achieves AUC-PR of 0.95–0.99 across diverse prokaryotic genomes and has been validated on 3,269 MAGs across 15 phyla. Runtime averages 1.3 seconds per genome on CPU.

MetalGenie-Evo runs UniOP as an independent step after HMM annotation. Operon predictions are an additional output and do not filter or modify HMM hits.

### Input detection

MetalGenie-Evo detects the FAA header format automatically:

- **Prodigal FAA headers** contain embedded coordinates (`# start # end # strand`) → UniOP is called with `-a faa_file`
- **Bakta / NCBI FAA headers** have no coordinates → requires `--fna_dir`; MetalGenie-Evo passes the FNA files to UniOP with `-i fna_file`

### Usage

```bash
# Prodigal FAA (headers have coordinates)
MetalGenie-Evo \
    --faa_dir          orfs/ \
    --hmm_dir          hmm_library/ \
    --out              results/ \
    --operon_prediction \
    --uniop_path       /path/to/UniOP/src/UniOP

# Bakta FAA (needs FNA for UniOP)
MetalGenie-Evo \
    --fna_dir          assemblies/ \
    --bakta_gff_dir    bakta_output/ \
    --hmm_dir          hmm_library/ \
    --out              results/ \
    --operon_prediction \
    --uniop_path       /path/to/UniOP/src/UniOP \
    --anvio
```

UniOP results are cached in `results/_uniop/<genome>/` and reused on subsequent runs.

---

## Anvi'o integration

The `--anvio` flag produces `MetalGenie-Evo-anvio-functions.tsv` for direct import with `anvi-import-functions`.

### Without Bakta mapping

If you used Prodigal FAA files directly, `gene_callers_id` contains Prodigal ORF names. Map these to Anvi'o integer IDs before importing:

```bash
# Export gene calls from your contigs database
anvi-export-gene-calls -c CONTIGS.db -o gene_calls.tsv

# Join on coordinates to get integer IDs, then import
anvi-import-functions \
    -c CONTIGS.db \
    -i MetalGenie-Evo-anvio-functions.tsv \
    -p MetalGenie-Evo
```

### With Bakta mapping (`--bakta_gff_dir`) — recommended

When using `--fna_dir` + `--bakta_gff_dir`, MetalGenie-Evo runs Prodigal internally, then matches each Prodigal ORF to its Bakta gene ID via coordinate overlap (±3 bp tolerance). The `gene_callers_id` column in the Anvi'o TSV will contain Bakta IDs (e.g. `AMXMAG_00053`), which map directly to Anvi'o integer IDs when the contigs database was built from Bakta external gene calls.

```bash
# Step 1 — Build Anvi'o contigs database from Bakta external gene calls
anvi-script-process-genbank \
    --input-genbank bakta_output/genome.gbff \
    --output-dir    anvio_input/

anvi-gen-contigs-database \
    -f anvio_input/genome.fa \
    --external-gene-calls anvio_input/gene_calls.tsv \
    -o CONTIGS.db

# Step 2 — Run MetalGenie-Evo with Bakta mapping
MetalGenie-Evo \
    --fna_dir       assemblies/ \
    --bakta_gff_dir bakta_output/ \
    --hmm_dir       hmm_library/ \
    --out           results/ \
    --anvio

# Step 3 — Import directly (Bakta IDs match Anvi'o gene_callers_id)
anvi-import-functions \
    -c CONTIGS.db \
    -i results/MetalGenie-Evo-anvio-functions.tsv \
    -p MetalGenie-Evo
```

---

## Output files

### Standard outputs (always written)

| File | Description |
|---|---|
| `MetalGenie-Evo-summary.csv` | Per-ORF detailed results with cluster separators |
| `MetalGenie-Evo-geneSummary-clusters.csv` | FeGenie R-script compatible compact summary |
| `MetalGenie-Evo-heatmap-data.csv` | Gene-count matrix (categories × genomes) |
| `MetalGenie-Evo-results-long.tsv` | Tidy long-format TSV, one row per ORF |

### Optional outputs

| File | Requires |
|---|---|
| `MetalGenie-Evo-coverage-heatmap.csv` | `--bam` / `--bams` / `--depth` / `--depths` |
| `MetalGenie-Evo-OperonStructure.tsv` | `--operon_prediction` |
| `MetalGenie-Evo-anvio-functions.tsv` | `--anvio` |

### Column reference — `MetalGenie-Evo-summary.csv`

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

### Column reference — `MetalGenie-Evo-OperonStructure.tsv`

| Column | Description |
|---|---|
| `operon_id` | UniOP operon ID (e.g. `genome_OP0001`) or `singleton_...` |
| `genome` | Source genome filename |
| `contig` | Source contig |
| `orf` | Prodigal ORF identifier |
| `bakta_gene_id` | Bakta gene ID — only present with `--bakta_gff_dir` |
| `gene` | Readable gene name |
| `category` | Functional category |
| `hmm_stem` | HMM model name |
| `bitscore` | hmmsearch bitscore |
| `e_value` | hmmsearch E-value |
| `unioperon_members` | Other HMM-positive ORFs in the same UniOP operon |

---

## Operon filtering logic

### Default: FeGenie exact port

By default (no `operon_rules.json` in `--hmm_dir`), MetalGenie-Evo uses an exact reimplementation of FeGenie's operon-context filtering. Each cluster is routed to exactly one rule handler based on which special genes are present. For iron acquisition categories, filtering uses a **per-ORF pass/break pattern**: a gene that does not meet the co-occurrence threshold is silently skipped, but other genes in the same cluster are still reported. This replicates FeGenie's original behaviour.

Use `--all_results` to bypass all filtering entirely.

### For model organisms: JSON rule engine

If `operon_rules.json` is present in `--hmm_dir`, the JSON rule engine is used instead. Intended for users working with well-characterised organisms.

```json
{
  "report_all_categories": ["metal_resistance-*", "iron_storage"],
  "rules": [
    {
      "name":       "FLEET",
      "categories": ["iron_oxidation"],
      "genes":      ["EetA","EetB","Ndh2","FmnB","FmnA","DmkA","DmkB","PplA"],
      "rule":       "require_n_of",
      "min_genes":  5,
      "on_fail":    "passthrough_non_members"
    }
  ]
}
```

`report_all_categories` accepts glob patterns — any cluster whose entire category set matches bypasses all rules. MetHMMDB categories (`metal_resistance-*`) and `iron_storage` bypass rules by default.

| Rule type | Behaviour |
|---|---|
| `require_n_of` | ≥ `min_genes` unique members of `genes` |
| `require_anchor` | Specific `anchor` gene must be present |
| `require_n_cat` | ≥ `min_genes` distinct HMMs in `categories` |
| `require_n_cat_or_lone_trusted` | As above, or lone hit in `trusted_lone` |
| `mtr_disambiguation` | Re-assigns category based on Mtr/Mto co-presence |

| `on_fail` | Behaviour |
|---|---|
| `passthrough_non_members` | Keep genes not in the rule's gene set |
| `drop` | Remove the entire cluster |
| `keep_all` | Keep everything |

---

## Extending the HMM library

```bash
# Full rebuild from all sources
python scripts/curate_hmm_library.py \
    --fegenie_dir   /path/to/FeGenie/hmms/iron/ \
    --flat_dir      /path/to/new_iron_hmms  iron_aquisition  tabuteau \
    --methmmdb_json /path/to/MetHMMDB/metadata.json \
    --methmmdb_dir  /path/to/MetHMMDB/ \
    --out_dir       hmm_library/ \
    --log           curation_report.tsv

# Verify
python scripts/curate_hmm_library.py --verify hmm_library/
```

To add individual HMMs manually: place `.hmm` files in a subdirectory of `hmm_library/` named after the category, add cutoffs to `HMM-bitcutoffs.txt`, names to `MetalGenie-map.txt`, and optionally rules to `operon_rules.json`.

---

## Visualisation

MetalGenie-Evo includes `scripts/plot_heatmap.R` for publication-ready heatmaps.

**Features:** hierarchical clustering (Ward's method) on both axes; white→orange→red for gene counts, white→blue for coverage; static PDF/PNG via `pheatmap`; interactive self-contained HTML via `plotly`.

**Install R packages:**
```bash
conda install -c conda-forge r-pheatmap r-plotly r-htmlwidgets \
                              r-optparse r-rcolorbrewer r-scales
```

**Usage:**
```bash
# Process entire results directory
Rscript scripts/plot_heatmap.R --input results/

# Static PDF only
Rscript scripts/plot_heatmap.R --input results/ --type static --format pdf --out figures/

# Interactive HTML, filter low-count categories
Rscript scripts/plot_heatmap.R --input results/ --type interactive --min_count 1
```

| Argument | Default | Description |
|---|---|---|
| `--input` | *(required)* | CSV file or output directory |
| `--type` | `both` | `static` \| `interactive` \| `both` |
| `--format` | `both` | `pdf` \| `png` \| `both` |
| `--out` | same as input | Output directory |
| `--width` / `--height` | auto | Plot dimensions in inches |
| `--min_count` | `0` | Min value to include a category row |
| `--no_cluster_rows` / `--no_cluster_cols` | off | Disable clustering axes |

---

## Differences from FeGenie at a glance

| Feature | FeGenie | MetalGenie-Evo |
|---|---|---|
| Operon clustering | ORF ordinal index | **bp coordinates (GFF)** + ordinal fallback |
| Strand awareness | No | **Yes (`--strand_aware`)** |
| Parallelism | Sequential | **ProcessPoolExecutor** |
| Operon filter logic | Hardcoded Python | **Exact port + JSON override for model organisms** |
| Result caching | No | **Yes (tblout cache)** |
| MetHMMDB support | No | **Yes** |
| Additional iron HMMs | No | **Yes (Tabuteau et al. 2025)** |
| HMM deduplication | No | **Yes (ACC + NAME + Jaccard, cross-source only)** |
| Versioned HMM registry | No | **Yes (`hmm_registry.tsv`)** |
| Normalised heatmap | `--norm` | **`--norm`** (same) |
| Coverage heatmap | `--bam` | **`--bam` / `--bams` / `--depth` / `--depths`** |
| Coverage normalisation | No | **TPM (`--norm_coverage`)** |
| Integrated Prodigal | No | **Yes (`--fna_dir --meta`)** |
| Contig length filter | No | **Yes (`--min_contig_len`)** |
| Relaxed operon thresholds | No | **Yes (`--relaxed_operons`)** |
| Contig column in outputs | No | **Yes** |
| Long-format output | No | **Yes (`-results-long.tsv`)** |
| Operon prediction | No | **Yes (UniOP, `--operon_prediction`)** |
| Anvi'o integration | No | **Yes (`--anvio`, `--bakta_gff_dir`)** |
| R script compatibility | Yes | **Yes (same CSV format)** |
| FeGenie filter rules | All | **All (exact port)** |

---

## License

MIT License. See [LICENSE](LICENSE).

MetalGenie-Evo is not affiliated with the original FeGenie authors.