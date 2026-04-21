#!/usr/bin/env bash
# MetalGenie-Evo  –  quick smoke test using synthetic data
# Requires: prodigal, hmmsearch (from HMMER) in PATH
# Usage: bash test_data/test_run.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
FAA_DIR="$SCRIPT_DIR/faa"
GFF_DIR="$SCRIPT_DIR/gff"
OUT_DIR="$SCRIPT_DIR/test_output"
HMM_DIR="$REPO_DIR/hmm_library"

echo "========================================"
echo "  MetalGenie-Evo  –  test run"
echo "========================================"

mkdir -p "$OUT_DIR"

echo "[INFO] Input FAA:  $FAA_DIR"
echo "[INFO] Input GFF:  $GFF_DIR"
echo "[INFO] HMM lib:    $HMM_DIR"
echo "[INFO] Output:     $OUT_DIR"
echo ""

# Run MetalGenie-Evo
python3 "$REPO_DIR/MetalGenie-Evo.py" \
    --faa_dir    "$FAA_DIR" \
    --gff_dir    "$GFF_DIR" \
    --hmm_dir    "$HMM_DIR" \
    --out        "$OUT_DIR" \
    --threads    4 \
    --max_bp_gap 5000 \
    --keep_tblout

echo ""
echo "[INFO] Output files:"
ls -lh "$OUT_DIR"/*.csv 2>/dev/null || echo "  (no CSV outputs – check for errors above)"

echo ""
echo "[DONE] Test run complete."
