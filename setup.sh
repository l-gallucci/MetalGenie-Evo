#!/usr/bin/env bash
# FeGenie-Evo setup script
# Installs dependencies via conda/mamba if available, otherwise checks for system tools.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV="fegenie-evo"

echo "========================================"
echo "  FeGenie-Evo  –  setup"
echo "========================================"

# ── Check for conda/mamba ───────────────────────────────────────────────────
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
else
    echo "[WARN] Neither conda nor mamba found."
    echo "       Checking for system dependencies instead..."
    CONDA_CMD=""
fi

if [[ -n "$CONDA_CMD" ]]; then
    echo "[INFO] Using $CONDA_CMD to create environment '$CONDA_ENV'..."
    $CONDA_CMD env create -n "$CONDA_ENV" \
        -f "$SCRIPT_DIR/environment.yml" --yes 2>/dev/null || \
    $CONDA_CMD env update -n "$CONDA_ENV" \
        -f "$SCRIPT_DIR/environment.yml" --yes

    echo ""
    echo "[OK] Environment '$CONDA_ENV' ready."
    echo "     Activate with:  conda activate $CONDA_ENV"
    echo "     Then run:       FeGenie-Evo.py --help"
else
    # Manual dependency check
    MISSING=()
    for tool in python3 hmmsearch prodigal blastp; do
        if ! command -v "$tool" &>/dev/null; then
            MISSING+=("$tool")
        fi
    done

    if [[ ${#MISSING[@]} -gt 0 ]]; then
        echo "[ERROR] Missing required tools: ${MISSING[*]}"
        echo "        Please install them manually or use conda."
        exit 1
    fi

    echo "[OK] All required tools found in PATH."
fi

# ── Make main script executable ─────────────────────────────────────────────
chmod +x "$SCRIPT_DIR/FeGenie-Evo.py"
chmod +x "$SCRIPT_DIR/scripts/build_hmm_library.py"

# ── Optionally add to PATH ───────────────────────────────────────────────────
echo ""
read -r -p "Add FeGenie-Evo.py to PATH in ~/.bashrc? [y/N] " response
if [[ "$response" =~ ^[Yy]$ ]]; then
    echo "export PATH=\"\$PATH:$SCRIPT_DIR\"" >> "$HOME/.bashrc"
    echo "[OK] Added to ~/.bashrc. Run 'source ~/.bashrc' to apply."
fi

echo ""
echo "========================================"
echo "  Setup complete!"
echo "  Test the tool:  bash test_data/test_run.sh"
echo "========================================"
