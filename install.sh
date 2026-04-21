#!/usr/bin/env bash
# MRBIGR2 installer — creates or updates a conda env for running the FastMCP server.
# Usage:
#   ./install.sh                              # default env name, Python 3.10
#   ./install.sh my_env_name 3.11             # custom env name, Python 3.11
#   ./install.sh --with-java                  # also install openjdk into the conda env
#   ./install.sh my_env_name 3.11 --with-java # custom env + optional Java
#
# The script is idempotent: if the env already exists, it only installs missing
# Python dependencies. Java is optional and can be installed with --with-java
# for strict ClusterONE support in the network module. Perl is still only verified
# for ANNOVAR-related scripts.
# Finally it writes a .mcp.json in the parent project directory so that Claude Code
# picks up the MCP server on next launch — no manual restart required.

set -eo pipefail

INSTALL_JAVA=0
POSITIONAL_ARGS=()

for arg in "$@"; do
    case "$arg" in
        --with-java)
            INSTALL_JAVA=1
            ;;
        -h|--help)
            cat <<'EOF'
Usage:
  ./install.sh [env_name] [python_version] [--with-java]

Examples:
  ./install.sh
  ./install.sh my_env 3.11
  ./install.sh --with-java
  ./install.sh my_env 3.11 --with-java

Notes:
  - Java is not installed by default.
  - Pass --with-java to install openjdk into the conda environment so that
    net.module_identify can use the bundled ClusterONE jar instead of the
    NetworkX fallback implementation.
EOF
            exit 0
            ;;
        *)
            POSITIONAL_ARGS+=("$arg")
            ;;
    esac
done

ENV_NAME="${POSITIONAL_ARGS[0]:-mrbigr2}"
PY_VERSION="${POSITIONAL_ARGS[1]:-3.10}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

py_deps=(
    "pandas>=1.3.0"
    "numpy>=1.20.0"
    "scipy>=1.7.0"
    "scikit-learn>=1.0.0"
    "matplotlib>=3.4.0"
    "seaborn>=0.11.0"
    "statsmodels>=0.13.0"
    "gseapy>=0.10.0"
    "networkx>=2.6.0"
    "fastmcp>=2.0.0"
)

# Map install spec -> python import name, for "is it already installed" checks.
declare -A import_name=(
    ["pandas>=1.3.0"]="pandas"
    ["numpy>=1.20.0"]="numpy"
    ["scipy>=1.7.0"]="scipy"
    ["scikit-learn>=1.0.0"]="sklearn"
    ["matplotlib>=3.4.0"]="matplotlib"
    ["seaborn>=0.11.0"]="seaborn"
    ["statsmodels>=0.13.0"]="statsmodels"
    ["gseapy>=0.10.0"]="gseapy"
    ["networkx>=2.6.0"]="networkx"
    ["fastmcp>=2.0.0"]="fastmcp"
)

echo "=================================================="
echo "MRBIGR2 installer"
echo "  target env : $ENV_NAME"
echo "  python     : $PY_VERSION"
echo "  with java  : $INSTALL_JAVA"
echo "  root       : $SCRIPT_DIR"
echo "=================================================="

# --- conda setup --------------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
    echo "[ERROR] conda not found in PATH. Install miniconda/anaconda first." >&2
    exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "[info] conda env '$ENV_NAME' already exists — checking dependencies."
else
    echo "[info] creating conda env '$ENV_NAME' (python=$PY_VERSION)"
    conda create -y -n "$ENV_NAME" "python=$PY_VERSION"
fi

conda activate "$ENV_NAME"

# --- python deps --------------------------------------------------------------
missing=()
for spec in "${py_deps[@]}"; do
    mod="${import_name[$spec]}"
    if python -c "import $mod" >/dev/null 2>&1; then
        echo "[ok]   $mod already installed"
    else
        echo "[miss] $mod"
        missing+=("$spec")
    fi
done

if [ "${#missing[@]}" -gt 0 ]; then
    echo "[info] installing missing packages: ${missing[*]}"
    pip install --quiet "${missing[@]}"
fi

# --- bundled binaries ---------------------------------------------------------
UTILS="$SCRIPT_DIR/utils"
for bin in plink gemma.linux FastTree; do
    if [ -f "$UTILS/$bin" ]; then
        chmod +x "$UTILS/$bin" || true
        echo "[ok]   bundled: utils/$bin"
    else
        echo "[warn] missing bundled binary: utils/$bin"
    fi
done

# --- java (for cluster_one) ---------------------------------------------------
if command -v java >/dev/null 2>&1; then
    java_ver="$(java -version 2>&1 | head -n 1)"
    echo "[ok]   java present: $java_ver"
elif [ "$INSTALL_JAVA" -eq 1 ]; then
    echo "[info] installing openjdk into conda env '$ENV_NAME'"
    conda install -y -n "$ENV_NAME" -c conda-forge openjdk
    if command -v java >/dev/null 2>&1; then
        java_ver="$(java -version 2>&1 | head -n 1)"
        echo "[ok]   java installed: $java_ver"
    else
        echo "[warn] requested Java installation, but java is still not in PATH after conda install." >&2
    fi
else
    cat <<'EOF' >&2
[warn] Java not found.
       The network module (net.module_identify) can still run with a
       NetworkX fallback, but strict ClusterONE results require Java.
       To install Java into the conda env, rerun:
         ./install.sh [env_name] [python_version] --with-java
       Or install one of:
         sudo apt-get install default-jre
         conda install -c conda-forge openjdk
EOF
fi

# --- perl (for ANNOVAR) -------------------------------------------------------
if command -v perl >/dev/null 2>&1; then
    perl_ver="$(perl -e 'print $^V')"
    echo "[ok]   perl present: $perl_ver"
else
    echo "[warn] Perl not found. The annotation module (anno.annotate_with_annovar)" >&2
    echo "       needs perl to run utils/*.pl. Install via your OS package manager." >&2
fi

# --- editable install ---------------------------------------------------------
# Expose `mrbigr` and `mrbigr-skill` console scripts and put the `mrbigr`
# package on sys.path so the import-verification step below can succeed.
echo "[info] installing the mrbigr package in editable mode"
pip install --quiet -e "$SCRIPT_DIR"

# --- final check --------------------------------------------------------------
echo "=================================================="
echo "Verifying import chain..."
python - <<'PY'
try:
    from mrbigr.core import geno, pheno, gwas, vis, anno, qtl, mr, go, net, peak
    import fastmcp
    print("[ok]   all modules import")
except Exception as e:
    print(f"[FAIL] {type(e).__name__}: {e}")
    raise SystemExit(1)
PY

# --- deploy Agent Skills-compatible workflow skills --------------------------
echo ""
echo "[info] Installing MRBIGR2 workflow skills for supported agents..."
MRBIGR_ROOT="$SCRIPT_DIR" mrbigr-skill install-all --target all
echo "[ok]   skills installed under supported Agent Skills directories"

# --- configure .mcp.json for Claude Code --------------------------------------
PYTHON_BIN="$(which python)"
SERVER_PY="$SCRIPT_DIR/src/server.py"

# Walk up to find a Claude Code project root (.claude/ dir or git root),
# otherwise default to the parent of SCRIPT_DIR.
find_project_root() {
    local dir="$SCRIPT_DIR"
    while [ "$dir" != "/" ]; do
        if [ -d "$dir/.claude" ] || [ -d "$dir/.git" ]; then
            echo "$dir"
            return
        fi
        dir="$(dirname "$dir")"
    done
    echo "$(dirname "$SCRIPT_DIR")"
}

PROJECT_ROOT="$(find_project_root)"
MCP_JSON="$PROJECT_ROOT/.mcp.json"

echo ""
echo "[info] Configuring MCP server for Claude Code..."
echo "  project root : $PROJECT_ROOT"
echo "  .mcp.json    : $MCP_JSON"

# Build JSON with the detected python path and absolute server.py path
cat > "$MCP_JSON" <<MCPEOF
{
  "mcpServers": {
    "mrbigr2": {
      "command": "$PYTHON_BIN",
      "args": ["$SERVER_PY"],
      "cwd": "$SCRIPT_DIR",
      "env": {
        "MRBIGR_ROOT": "$SCRIPT_DIR"
      }
    }
  }
}
MCPEOF

echo "[ok]   .mcp.json written"

# --- CLI wrapper already bundled in package -----------------------------------
WRAPPER="$SCRIPT_DIR/mrbigr"
chmod +x "$WRAPPER" 2>/dev/null || true
echo "[ok]   CLI wrapper: $WRAPPER"

echo "=================================================="
echo "Install complete."
echo ""
echo "Immediate use (current session, no restart needed):"
echo "  $WRAPPER list                                # show all tools"
echo "  $WRAPPER gwas_lmm --phe phe.csv --geno geno  # run GWAS"
echo "  $WRAPPER get_top_snps --gwas_file result.txt  # top SNPs"
echo ""
echo "Next session (MCP auto-loaded, no CLI needed):"
echo "  Restart Claude Code in $PROJECT_ROOT"
echo "  The MCP server 'mrbigr2' will be auto-detected."
echo "  MRBIGR2 skills will be available globally in supported agents."
echo "  Verify in Claude Code/Codex/OpenCode by listing available skills."
echo "=================================================="
