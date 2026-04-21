# Note on HMM Library Population

The `hmm_library/` directory is intentionally empty in the repository.
You need to populate it once before running MetalGenie-Evo.

## Option 1 – FeGenie HMMs only (iron genes)

```bash
python scripts/build_hmm_library.py \
    --fegenie_dir /path/to/FeGenie/hmms/iron/ \
    --out_dir     hmm_library/
```

## Option 2 – FeGenie + MetHMMDB (metal resistance genes)

MetHMMDB provides individual HMM files in its `individual/` subfolder.
Using that is preferred over the concatenated file.

```bash
# Clone MetHMMDB
git clone https://github.com/Haelmorn/MetHMMDB
# or: git clone https://github.com/kciuchcinski/MetHMMDB

python scripts/build_hmm_library.py \
    --fegenie_dir  /path/to/FeGenie/hmms/iron/ \
    --methmmdb_dir MetHMMDB/individual/ \
    --out_dir      hmm_library/
```

If you only have the concatenated file:

```bash
python scripts/build_hmm_library.py \
    --fegenie_dir /path/to/FeGenie/hmms/iron/ \
    --methmmdb    MetHMMDB/MetHMMDb.hmm \
    --out_dir     hmm_library/
```

## Option 3 – Update to a new MetHMMDB release

```bash
python scripts/build_hmm_library.py \
    --methmmdb_dir MetHMMDB_v2/individual/ \
    --out_dir      hmm_library/ \
    --update
```

The `--update` flag preserves existing models and only adds/replaces new ones.
Removed models are flagged as `status=removed` in `hmm_registry.tsv`.

## Adding custom operon rules

Copy `hmm_library/operon_rules.json` and add rules for your new categories.
See the main README for the rule schema.
