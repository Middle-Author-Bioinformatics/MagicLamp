#!/usr/bin/env bash
# Create a dedicated environment without removing existing environments or changing channels.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME=magiclamp
SPEC=environment.yml
DRY_RUN=0
REGISTER_ONLY=0
while (($#)); do
    case "$1" in
        --name) [[ $# -ge 2 ]] || { echo '--name needs a value' >&2; exit 2; }; ENV_NAME=$2; shift 2 ;;
        --plots) SPEC=environment-plots.yml; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --register-only) REGISTER_ONLY=1; shift ;;
        -h|--help) echo 'Usage: bash setup.sh [--name NAME] [--plots] [--dry-run] [--register-only]'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ "$ENV_NAME" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo 'Invalid environment name' >&2; exit 2; }
if ((DRY_RUN)); then
    printf 'Environment: %s\nSpecification: %s/%s\nRegister only: %s\n' "$ENV_NAME" "$SCRIPT_DIR" "$SPEC" "$REGISTER_ONLY"
    exit 0
fi
command -v conda >/dev/null || { echo 'Install conda/mamba before running setup.sh.' >&2; exit 1; }
if ((!REGISTER_ONLY)); then
    # conda refuses an existing prefix; it is never removed automatically.
    conda env create --name "$ENV_NAME" --file "$SCRIPT_DIR/$SPEC" --yes
fi
ENV_PREFIX="$(conda run --no-capture-output -n "$ENV_NAME" python -c 'import sys; print(sys.prefix)')"
[[ -d "$ENV_PREFIX/conda-meta" ]] || { echo 'Could not locate the conda environment.' >&2; exit 1; }
mkdir -p "$ENV_PREFIX/bin"
LAUNCHER="$ENV_PREFIX/bin/MagicLamp.py"
if [[ -e "$LAUNCHER" ]] && ! head -n 3 "$LAUNCHER" | grep -qF '# MagicLamp managed launcher'; then
    echo "Refusing to overwrite an unrelated executable: $LAUNCHER" >&2
    exit 1
fi
{
    echo '#!/usr/bin/env bash'
    echo '# MagicLamp managed launcher'
    printf 'exec %q %q "$@"\n' "$ENV_PREFIX/bin/python" "$SCRIPT_DIR/MagicLamp.py"
} > "$LAUNCHER"
chmod +x "$LAUNCHER"
printf 'Ready. Run: conda activate %s\nThen: MagicLamp.py help\n' "$ENV_NAME"
