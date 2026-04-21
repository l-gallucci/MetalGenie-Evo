#!/usr/bin/env bash
set -e

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL="$REPO_DIR/MetalGenie-Evo"

echo "========================================"
echo "  MetalGenie-Evo  –  setup"
echo "========================================"

# ── Check we are inside a conda environment ───────────────────────────────────
if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
    echo "[ERROR] No conda environment is active."
    echo ""
    echo "  Please create and activate the environment first:"
    echo "    conda env create -f environment.yml"
    echo "    conda activate metalgenie-evo"
    echo "  Then re-run:  bash setup.sh"
    exit 1
fi

echo "[INFO] Using conda environment: $CONDA_DEFAULT_ENV"

# ── Check required tools are available ───────────────────────────────────────
MISSING=()
for tool in python3 hmmsearch prodigal; do
    command -v "$tool" &>/dev/null || MISSING+=("$tool")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "[ERROR] Missing tools in current environment: ${MISSING[*]}"
    echo "        Make sure the environment was created from environment.yml"
    exit 1
fi

# ── Executable permissions ────────────────────────────────────────────────────
chmod +x "$TOOL"
echo "[INFO] $TOOL is executable"

# ── Install wrapper into ~/.local/bin ─────────────────────────────────────────
LOCAL_BIN="$HOME/.local/bin"
mkdir -p "$LOCAL_BIN"
WRAPPER="$LOCAL_BIN/MetalGenie-Evo"

# Locate conda init script
CONDA_SH=""
for candidate in "$HOME/miniforge3/etc/profile.d/conda.sh" \
                 "$HOME/mambaforge/etc/profile.d/conda.sh" \
                 "$HOME/anaconda3/etc/profile.d/conda.sh" \
                 "$HOME/miniconda3/etc/profile.d/conda.sh" \
                 "${CONDA_PREFIX}/../../etc/profile.d/conda.sh"; do
    if [[ -f "$candidate" ]]; then CONDA_SH="$candidate"; break; fi
done

cat > "$WRAPPER" << WRAPPER_EOF
#!/usr/bin/env bash
_CONDA_SH="$CONDA_SH"
if [[ -n "\$_CONDA_SH" && -f "\$_CONDA_SH" ]]; then
    source "\$_CONDA_SH"
    conda activate $CONDA_DEFAULT_ENV 2>/dev/null
fi
exec "$TOOL" "\$@"
WRAPPER_EOF

chmod +x "$WRAPPER"
echo "[INFO] Wrapper installed: $WRAPPER"

# ── Check if ~/.local/bin is in PATH ─────────────────────────────────────────
if echo "$PATH" | grep -q "$LOCAL_BIN"; then
    NEEDS_SOURCE=false
else
    if [[ -f "$HOME/.zshrc" ]]; then
        SHELL_RC="$HOME/.zshrc"
    else
        SHELL_RC="$HOME/.bashrc"
    fi
    if ! grep -qF ".local/bin" "$SHELL_RC" 2>/dev/null; then
        echo "" >> "$SHELL_RC"
        echo "# Added by MetalGenie-Evo setup" >> "$SHELL_RC"
        echo "export PATH=\"\$HOME/.local/bin:\$PATH\"" >> "$SHELL_RC"
    fi
    NEEDS_SOURCE=true
fi

echo ""
echo "========================================"
echo "  Setup complete!"
echo ""
if [[ "$NEEDS_SOURCE" == true ]]; then
    echo "  Run once to reload your shell:"
    echo "    source $SHELL_RC"
    echo ""
fi
echo "  Then from any directory:"
echo "    MetalGenie-Evo --help"
echo ""
echo "  The wrapper activates '$CONDA_DEFAULT_ENV' automatically."
echo "========================================"