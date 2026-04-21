#!/usr/bin/env bash
# MetalGenie-Evo setup — verifies the installation
# Not required if you used: conda env create -f environment.yml
# (the environment.yml already runs pip install -e . automatically)
set -e

echo "========================================"
echo "  MetalGenie-Evo  –  setup"
echo "========================================"

# ── Check conda environment is active ────────────────────────────────────────
if [[ -z "$CONDA_PREFIX" ]]; then
    echo "[ERROR] No conda environment is active."
    echo ""
    echo "  Please run:"
    echo "    conda env create -f environment.yml"
    echo "    conda activate metalgenie-evo"
    exit 1
fi
echo "[INFO] Active environment : $CONDA_DEFAULT_ENV"

# ── Check required tools ──────────────────────────────────────────────────────
MISSING=()
for tool in python3 hmmsearch prodigal; do
    command -v "$tool" &>/dev/null || MISSING+=("$tool")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "[ERROR] Missing tools: ${MISSING[*]}"
    echo "        Run: conda env create -f environment.yml"
    exit 1
fi
echo "[INFO] All required tools found"

# ── Install package if not already installed ──────────────────────────────────
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! command -v MetalGenie-Evo &>/dev/null; then
    echo "[INFO] Installing MetalGenie-Evo..."
    pip install -e "$REPO_DIR" --quiet
fi

# ── Verify ────────────────────────────────────────────────────────────────────
if command -v MetalGenie-Evo &>/dev/null; then
    echo "[INFO] MetalGenie-Evo installed successfully"
else
    echo "[ERROR] Installation failed"
    exit 1
fi

echo ""
echo "========================================"
echo "  Setup complete!"
echo ""
echo "  Usage:"
echo "    conda activate $CONDA_DEFAULT_ENV"
echo "    MetalGenie-Evo --help"
echo "========================================"
