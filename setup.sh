#!/usr/bin/env bash
# Provision dependencies, install a copy of the package/models, then verify it.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME=magiclamp
SPEC=environment.yml
DRY_RUN=0
INSTALL_ONLY=0
CHECK_ARGS=(doctor --smoke-test)
while (($#)); do
    case "$1" in
        --name) [[ $# -ge 2 ]] || { echo '--name needs a value' >&2; exit 2; }; ENV_NAME=$2; shift 2 ;;
        --plots) SPEC=environment-plots.yml; CHECK_ARGS+=(--plots); shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --install-only|--register-only) INSTALL_ONLY=1; shift ;;
        -h|--help) echo 'Usage: bash setup.sh [--name NAME] [--plots] [--dry-run] [--install-only]'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ "$ENV_NAME" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo 'Invalid environment name' >&2; exit 2; }
if ((DRY_RUN)); then
    printf 'Environment: %s\nSpecification: %s/%s\nInstall into existing environment: %s\n' "$ENV_NAME" "$SCRIPT_DIR" "$SPEC" "$INSTALL_ONLY"
    exit 0
fi
command -v conda >/dev/null || { echo 'Install Miniforge/conda first; see docs/INSTALL.md.' >&2; exit 1; }
if ((!INSTALL_ONLY)); then
    # Refuse an existing environment; never remove it or change global channels.
    # Older conda-env releases reject --yes. The process-local setting works
    # with both older and current conda without changing the user's config.
    CONDA_ALWAYS_YES=true conda env create --name "$ENV_NAME" --file "$SCRIPT_DIR/$SPEC"
fi
ENV_PREFIX="$(conda run --no-capture-output -n "$ENV_NAME" python -c 'import sys; print(sys.prefix)')"
[[ -d "$ENV_PREFIX/conda-meta" ]] || { echo 'Could not locate the conda environment.' >&2; exit 1; }
# Refuse unrelated executables; allow replacing the old managed launcher or our installed entry point.
for COMMAND in MagicLamp.py magiclamp; do
    LAUNCHER="$ENV_PREFIX/bin/$COMMAND"
    if [[ -e "$LAUNCHER" || -L "$LAUNCHER" ]]; then
        if ! head -n 12 "$LAUNCHER" | grep -Eq '# MagicLamp managed launcher|from MagicLamp import main|from latest.cli import main'; then
            echo "Refusing to overwrite an unrelated executable: $LAUNCHER" >&2
            exit 1
        fi
    fi
done
# Python dependencies/build tools are supplied by conda; do not change them with pip.
conda run --no-capture-output -n "$ENV_NAME" python -m pip install --no-deps --no-build-isolation "$SCRIPT_DIR"
conda run --no-capture-output -n "$ENV_NAME" magiclamp "${CHECK_ARGS[@]}"
printf 'Ready. Run: conda activate %s\nThen: magiclamp GramGenie --help\n' "$ENV_NAME"
