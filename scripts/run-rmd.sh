#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=sturcotte_rmd
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=endikaprieto@vhio.net
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=100gb
#SBATCH --time=10:00:00
#SBATCH --output=logs/%x_%j.log

# Runs all chunks of an .Rmd (except eval=FALSE ones) in the project environment, without rendering a report.
# Usage (from the project root, sturcotte/). Run it directly: it creates logs/ and submits itself with sbatch
# (any option before the .Rmd is passed to sbatch):
# cd /mnt/bioinfnas/immuno/eprieto/analyses/sturcotte
# scripts/run-rmd.sh -J deconv-pt3385 scripts/deconvolution/deconv-analysis.Rmd
# scripts/run-rmd.sh -J tcr-viewer --mem=80gb scripts/reports/tcr-umap-viewer.Rmd
# (sbatch -J <name> scripts/run-rmd.sh <file>.Rmd also works, but only if logs/ already exists)

set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$PWD}"
if [[ ! -x "$REPO/scripts/interactive-R.sh" ]]; then
  echo "Submit from the project root (sturcotte/): $REPO" >&2
  exit 1
fi

# Submission mode (outside a Slurm job): Slurm opens --output before the job starts and
# does not create its directory, so logs/ must exist before calling sbatch
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  if [[ $# -lt 1 || ! -f "${!#}" ]]; then
    echo "Usage: scripts/run-rmd.sh [-J <name>] scripts/<step>/<file>.Rmd" >&2
    exit 1
  fi
  mkdir -p "$REPO/logs"
  exec sbatch "${@:1:$#-1}" "$0" "${!#}"
fi

if [[ $# -ne 1 || ! -f "$1" ]]; then
  echo "Usage: scripts/run-rmd.sh [-J <name>] scripts/<step>/<file>.Rmd" >&2
  exit 1
fi
RMD="$(realpath "$1")"

exec "$REPO/scripts/interactive-R.sh" --no-save --no-restore -s -e "
  code <- tempfile(fileext = '.R')
  invisible(knitr::purl('$RMD', output = code, quiet = TRUE, documentation = 0))
  options(device = function(...) grDevices::pdf(file = NULL, width = 24, height = 12), warn = 1)
  source(code, echo = TRUE, max.deparse.length = Inf)
"