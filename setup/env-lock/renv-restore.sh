#!/bin/bash
# Installs the packages recorded in renv-sturcotte/renv.lock into the project library (from source).
# Usage, from the project root on a compute node:  sbatch -p long -c 24 --mem=64G -t 48:00:00 setup/env-lock/renv-restore.sh
# On a slow shared disk one package with thousands of files (BH) can take more than renv's default hour,
# so the per-package time limit is raised and finished packages are kept if another one fails.
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
export MAKEFLAGS="-j4"
export RENV_CONFIG_INSTALL_JOBS=4
export RENV_CONFIG_INSTALL_TRANSACTIONAL=FALSE
scripts/interactive-R.sh --no-save --no-restore -e '
options(Ncpus = 4, timeout = 3600, renv.install.timeout = 12 * 3600)
renv::restore(prompt = FALSE)
print(renv::status())
'
