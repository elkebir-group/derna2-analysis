#!/usr/bin/env bash
# Run the derna2 paper pipeline (pipelines/Snakefile) with the flags every run needs.
#
#   pipelines/run.sh [--dry-run] [--cores N] [-- TARGET ...]
#
# A target is a rule name (e1_figures, fig1_runs, ...) or an output path relative to the
# repository root, such as data/derived/figures/fig1_glpe_hbb.pdf. No target means `all`:
# every figure, from the deposited results.
#
# Loads pipelines/config/host.yaml after config.yaml when it exists. Passes
# --rerun-triggers mtime, so a code or config change never reruns an existing output by
# itself: delete the output, or add --forcerun <rule>, to redo it.
#
# Written for bash 3.2 (macOS) and 4.2 (CentOS 7) as well as current bash: under `set -u`
# those treat an empty array as unset, so arrays are only expanded through ${a[@]+...}.
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
snakemake=${SNAKEMAKE:-snakemake}
if ! command -v "$snakemake" > /dev/null 2>&1; then
    if [ -x "$root/.venv/bin/snakemake" ]; then
        snakemake=$root/.venv/bin/snakemake
    else
        echo "snakemake not found: install the pipeline extra, or set SNAKEMAKE" >&2
        exit 1
    fi
fi

cores=$(($(getconf _NPROCESSORS_ONLN) > 2 ? $(getconf _NPROCESSORS_ONLN) - 2 : 1))
extra=()
while [ $# -gt 0 ]; do
    case $1 in
        --dry-run) extra+=(--dry-run); shift ;;
        --cores) cores=$2; shift 2 ;;
        --) shift; break ;;
        *) extra+=("$1"); shift ;;
    esac
done

configs=("$root/pipelines/config/config.yaml")
if [ -f "$root/pipelines/config/host.yaml" ]; then
    configs+=("$root/pipelines/config/host.yaml")
fi

cd "$root"
exec "$snakemake" \
    --snakefile pipelines/Snakefile \
    --configfile "${configs[@]}" \
    --cores "$cores" \
    --rerun-triggers mtime \
    --printshellcmds \
    ${extra[@]+"${extra[@]}"} \
    "$@"
