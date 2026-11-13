# Identification and prioritization of candidate antigen-responsive TCR clonotypes following APC co-culture

## CHUM-MED x VHIO-AG

This repository contains the code used to find, in five patients (**pt3385, pt600280, pt600311, pt3328 and pt3513**), the T-cell receptor (TCR) clonotypes that become activated when expanded tumor-infiltrating lymphocytes (TILs) are co-cultured with antigen-presenting cells (APCs) loaded with tumor antigens, and to relate them to the T cells found in the tumor.

## Workflow

1. **Preprocessing** (section 4). TCR sequences are reconstructed with MiXCR from the single-cell V(D)J FASTQ files of each experimental condition (`Unpulsed`, `CD3CD28`, `Pept_pool` and individual `R` conditions), using IMGT-derived references, and restricted to the cells of the Seurat objects.
2. **Deconvolution analysis** (section 5). For each patient, cells are filtered, CD4 and CD8 T cells are separated, an activation signature is derived from the stimulated conditions, and every clonotype is tested for activation in each condition against the Unpulsed background. Responding clonotypes are ranked for functional validation (Tier and Level), with putative negative controls.
3. **Reports** (section 6). The prioritized clonotypes are matched to the single-cell data of the TILs, and two interactive HTML reports are built: a viewer of every clonotype in the experimental conditions and in the TILs, and an analysis of the expansion of the TIL clonotypes of one cluster.

```
1. preprocessing    deconv-cell-barcode ─► script_fastq_concat ─► check_samplesheet_vs_fastq
                    ─► script_mixcr_imgt ─► script_extract_cell_alleles ──┐
                                                                         │
2. deconvolution    deconv-analysis (one run per patient and compartment) ◄┘
                    └─► patients/<patient>/deconv/tcr/prioritized_<compartment>_tcrs_for_validation.tsv
                                                                         │
3. reports          til-frequency ◄── TIL single-cell object ────────────┤
                    └─► tcr-umap-viewer ─► patients/prioritized_tcrs_umap_viewer.html
                    cxcl13-expansion ───► patients/cd4-cxcl13-expansion/CD4_CXCL13_expansion_report.html
```

The computational environment is pinned:

| Layer | Pinned by |
| --- | --- |
| R 4.5.2 + system libraries | `setup/env-lock/conda_renv_sturcotte.explicit.txt` (conda, linux-64) |
| 300 R packages (CRAN, Bioconductor 3.22, GitHub commits) | `renv-sturcotte/renv.lock` (renv 1.3.1) |
| Python 3.11 for MAGIC (`magic-impute` 3.0.0) | `setup/env-lock/conda_magic_env.explicit.txt` + `magic_env.pip.txt` |
| MiXCR (Docker image, licence required) | `scripts/preprocessing/script_mixcr_imgt.sh` |

## Repository layout

```
sturcotte/
├── scripts/                          # interactive-R.sh, run-rmd.sh (SLURM launcher, section 3)
│   ├── preprocessing/                # FASTQ → MiXCR TCR assignments per cell (section 4)
│   ├── deconvolution/                # deconv-analysis.Rmd: QC, activation, TCR prioritization (section 5)
│   └── reports/                      # TIL frequencies and the two HTML reports, with their page templates (section 6)
├── setup/                            # setup-m1.R (libraries, functions, plot style), install-Rmagic-venv.Rmd
│   ├── env-lock/                     # conda/pip locks and how renv.lock was derived (section 1)
│   └── refs/                         # ProjecTILs reference maps (not tracked; downloaded by deconv-analysis.Rmd)
├── renv-sturcotte/                   # renv project: renv.lock, activate.R, settings.json
├── patients/                         # outputs per patient and reports (not tracked)
├── til-data/                         # single-cell object of the TILs (not tracked, section 2)
└── logs/                             # SLURM logs (not tracked; created by run-rmd.sh)
```

## 1. Recreate the environment

Requirements: Linux x86_64 with SLURM, [miniforge/conda](https://github.com/conda-forge/miniforge) and git. Loading the largest Seurat objects needs about 80 GB of RAM.

```bash
git clone https://github.com/eprieto012/sturcotte-vhio.git sturcotte
cd sturcotte

# R 4.5.2 and system libraries
conda create -p ./envs/renv_sturcotte --file setup/env-lock/conda_renv_sturcotte.explicit.txt

# Python environment for MAGIC
conda create -p ./envs/magic_env --file setup/env-lock/conda_magic_env.explicit.txt
./envs/magic_env/bin/pip install -r setup/env-lock/magic_env.pip.txt

# Tell the launchers and setup-m1.R where the environments are
export ENV_REPRO="$PWD/envs/renv_sturcotte"
export MAGIC_PYTHON="$PWD/envs/magic_env/bin/python"

# Install the exact R package versions (from source: several hours)
sbatch -p long -c 24 --mem=64G -t 48:00:00 setup/env-lock/renv-restore.sh
```

`scripts/interactive-R.sh` puts the conda R first on `PATH`, unsets `R_LIBS*` and starts R from `renv-sturcotte/`, whose `.Rprofile` activates renv, so only the renv library is used. Always start R through it. The scripts find the project folder by walking up from the working directory, so they also work when run from their own folder (e.g. in RStudio). To check the environment:

```r
renv::project()   # .../sturcotte/renv-sturcotte
renv::status()    # should report: "No issues found"
```

### How the lockfile was derived

The analyses were first run in an RStudio Server Docker image (R 4.5.2, Bioconductor 3.22) with a user library. `renv.lock` records the versions installed there of every package the scripts use, including optional packages that change results when present (`presto`, `glmGamPoi`); the recommended packages have the versions shipped with R 4.5.2. The files that document this are in `setup/env-lock/`: `docker-installed-packages.tsv` and `docker-session.txt` (the original environment), `dump-docker-env.R` and `make-lock.R` (how they were read).

## 2. Input data

| Data | Where in the scripts | Tracked |
| --- | --- | --- |
| Single-cell V(D)J FASTQ files per patient and condition | `RAW_DATA=/mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/raw_data` | no |
| Seurat objects of the co-culture experiments, one per patient | `SEURAT=/mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/seurat` | no |
| Deep TCRβ sequencing: productive clonotype frequencies in Tumor, Blood, PE, Bulk and TIDAL01 | `/mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/tcrb/clonotypes_adaptive_bulk_productivefreqs_update.rds` (copy in `patients/`) | no |
| Single-cell CITE-seq object of the TILs (13 patients, PD1+ and PD1dim sorted fractions) | `til-data/cite1_2_3_4_5_6_8_withoutPat400304_TCR_ADT_virus_singlets.rds` | no |

The absolute paths are those of the VHIO cluster; edit them at the top of each script to run elsewhere.

## 3. Running the R Markdown scripts

Interactively, on a compute node (not on the front end):

```bash
srun --partition=highmem --nodes=1 --ntasks=1 --cpus-per-task=25 --mem=100G --time=0-8:00:00 --pty /bin/bash -i
scripts/interactive-R.sh
```

As a SLURM job, from the project root. `run-rmd.sh` runs every chunk of the script (except those marked `eval=FALSE`) without rendering a document, and writes the log to `logs/`:

```bash
scripts/run-rmd.sh -J tcr-viewer scripts/reports/tcr-umap-viewer.Rmd
```

## 4. Preprocessing

Scripts in `scripts/preprocessing/`. Examples for `pt3385`.

1. **MiXCR sample sheets from the Seurat cell barcodes**: `deconv-cell-barcode.Rmd`. Writes one `10x-sample-sheet.tsv` per condition in `SEURAT/pt3385/10x_sample_sheets/<CONDITION>/`.
2. **Concatenate FASTQ files of different lanes**, only when a condition has several `Lane*` directories (it overwrites the merged files in `Data`):
   ```bash
   sbatch scripts/preprocessing/script_fastq_concat.sh $RAW_DATA/pt3385
   ```
3. **Check the sample-sheet / FASTQ condition mapping** before MiXCR:
   ```bash
   cd $RAW_DATA && python3 /path/to/sturcotte/scripts/preprocessing/check_samplesheet_vs_fastq.py --max-reads 1000000
   ```
4. **MiXCR and TCR allele inference**, once per condition (`10x-sc-xcr-vdj` preset, IMGT-derived reference library; Docker and a MiXCR licence are required):
   ```bash
   sbatch scripts/preprocessing/script_mixcr_imgt.sh \
     $RAW_DATA/pt3385/Unpulsed \
     $SEURAT/pt3385/10x_sample_sheets/Unpulsed/10x-sample-sheet.tsv \
     --no-skip --no-validate-only --check-imgt-updates
   ```
   Outputs: `RAW_DATA/pt3385/<CONDITION>/imgt/` and `.../validation/`.
5. **Keep the MiXCR records of the Seurat cells** (`--dedup no` keeps separate TRA/TRB rows):
   ```bash
   sbatch scripts/preprocessing/script_extract_cell_alleles.sh --root $RAW_DATA --patients pt3385 --dedup no
   ```
   Output: `RAW_DATA/pt3385/<CONDITION>/imgt/*.cell_alleles_filtered.tsv`.

## 5. Deconvolution analysis

`scripts/deconvolution/deconv-analysis.Rmd`, run once per patient and compartment (set `sample`, and `compartment` to `cd4` or `cd8`, at the top of the corresponding sections). It combines Cell Ranger and MiXCR annotations, performs QC, doublet detection, MAGIC imputation and CD4/CD8 assignment, derives the activation signature, scores each cell (UCell), tests each clonotype in each condition, adds the bulk TCRβ repertoire context and selects candidates and negative controls. Libraries and functions are loaded from `setup/setup-m1.R`; if the MAGIC Python environment does not exist yet, see `setup/install-Rmagic-venv.Rmd`.

Outputs for each patient, in `patients/<patient>/deconv/`:

- `tcr/prioritized_<compartment>_tcrs_for_validation.tsv`: prioritization table.
- `plots/`: QC, cell-type and UMAP figures (`IC01` to `IC11`).
- `markers/`: DEGs of each condition against Unpulsed and the activation signature.
- `seurat/`: filtered object (`obj_filt.rds`) and one object per compartment (`obj_<compartment>_TCR.rds`).

### Interpreting the prioritisation tables and UMAPs

#### Table columns

| Column                                           | Interpretation                                                                                                                                                                                                          |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Prioritize`                                   | `TRUE` marks a positive candidate selected for  functional validation. `FALSE` means it was not selected as a positive candidate; this also includes putative negative controls.                                   |
| `Tier` / `Level`                             | Summarise the evidence for antigen responsiveness and the selection stringency, or identify a putative negative-control class.                                                                                          |
| `Condition`                                    | Stimulation condition with significant activation, such as an individual `R` condition or `Pept_pool`. Blank for negative controls.                                                                                 |
| `Unpulsed`, `CD3CD28`, `Pept_pool`, `R*` | Activated cells / total clonotype cells, with the activation percentage. A cell is classified as activated when its gene-expression-based UCell activation score exceeds the 99th percentile of scores in `Unpulsed`. |
| `*_padj`                                       | Benjamini–Hochberg-adjusted P value for increased activation in that condition; values below 0.01 are significant.                                                                                                     |
| `*_cluster`, `*_phenotype`                   | Most frequent Seurat cluster or functional cell-state annotation for that clonotype and condition, with the proportion of cells.                                                                                        |
| Repertoire frequency / rank                      | Clonotype abundance and rank in the available bulk repertoires (e.g. Tumor, Blood, PE, Bulk and TIDAL01). These data provide context but do not prove antigen specificity.                                             |

#### Tier 1 vs Tier 2: where is the response?

| Tier                                           | Definition                                                                                                                                 |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| **Tier 1 — deconvoluted response**      | Significant activation (`padj < 0.01`) in at least one individual `R` condition, with or without a significant `Pept_pool` response. |
| **Tier 2 — peptide-pool response only** | Significant activation in `Pept_pool`, with no significant activation in any individual `R` condition.                                 |

#### L1 vs L2: which positive candidates are selected?

All prioritised positive clonotypes must have at least one cell in `Unpulsed` and less than 5% basal activation. The condition used to meet an activation-percentage threshold must contain at least six clonotype cells.

| Category            | Selection criteria                                                                                                                                          |
| ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Tier 1-L1** | At least one `R` condition has `padj < 0.01` and more than 60% activated cells.                                                                         |
| **Tier 1-L2** | At least one `R` condition is significant, but none meets L1; `Pept_pool` has at least 50% activated cells. `Pept_pool` significance is not required. |
| **Tier 2-L1** | `Pept_pool` is significant and has more than 60% activated cells; no `R` condition is significant.                                                      |
| **Tier 2-L2** | `Pept_pool` is significant and has 50–60% activated cells; no `R` condition is significant. Used only if there are no Tier 2–L1 candidates.           |

L2 is a rescue category, not automatically a weak response. A Tier 1–L2 clonotype can respond strongly to `Pept_pool` but fail the stricter L1 requirement in an individual `R` condition.

#### Putative negative controls (C1–C3)

Negative controls have no significant antigen-induced activation and are marked `Tier = Negative control`, `Prioritize = FALSE`, with a blank
`Condition`.

| Control                                                | Selection criteria                                                                                                                                                                                    |
| ------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **C1 — strongest control evidence**             | Detected in `Unpulsed`; more than 50% activation among more than five `CD3CD28` cells; less than 2% activation in every detected `Pept_pool`/`R` condition, with at least one such condition. |
| **C2 — lower CD3CD28 cell count**               | Same requirements as C1, but with 1–5 `CD3CD28` cells.                                                                                                                                             |
| **C3 — not detected after antigen stimulation** | Absent from all `Pept_pool`/`R` conditions, with at least two cells in `Unpulsed` and/or `CD3CD28`.                                                                                           |

Clonotypes detected in both control conditions are preferred. Controls are selected in the order C1 → C2 → C3 until three are obtained. These are putative controls, not experimentally verified non-reactive TCRs; C3 may reflect limited sampling.

#### UMAP plots

IC11 displays up to 12 prioritised positive TCRs, following the table order; the complete candidate list remains in the table. Each plot
highlights cells belonging to one clonotype, coloured by Seurat cluster, while other cells are grey. Panels correspond to `Unpulsed`, `CD3CD28`,
`Pept_pool` and individual `R` conditions.

## 6. Reports

Scripts in `scripts/reports/`. They read the outputs of section 5 and the TIL object, so they need a compute node (about 80 GB of RAM).

### 6.1 Frequency of the prioritized clonotypes in the TILs

`til-frequency.Rmd` adds to every prioritization table the number of TIL cells of the same patient carrying each clonotype, as `cells / TIL cells with a TCR (%)`, and its rank among the TIL clonotypes, for both sorted fractions together (`TIL_sc`) and separately (`TIL_sc_PD1p`, `TIL_sc_PD1d`). `TIL_sc_matched_by` and `TIL_sc_match_key` record how each clonotype was matched:

- by the CDR3β nucleotide sequence and TRBJ gene, when the deconvolution object contains V(D)J sequences (pt3385, pt600280, pt600311);
- by the clonotype identifier otherwise (pt3328, pt3513), and only for identifiers present in the bulk TCRβ table: identifiers outside it were assigned independently in each dataset, so they are left as `NA`.

The existing columns of the tables are not modified, and the script can be run again (it replaces its own columns).

### 6.2 Clonotype viewer

`tcr-umap-viewer.Rmd` builds `patients/prioritized_tcrs_umap_viewer.html`, a self-contained page (no server or internet needed). Run `til-frequency.Rmd` first. For the selected patient, compartment and clonotype it shows:

- the cells of the clonotype on the UMAP of each experimental condition, coloured by Seurat cluster, with its activation, adjusted P value, main cluster and phenotype per condition, Tier and Level;
- its frequency and rank in the TILs, how it was matched, and its cells on the TIL UMAP next to the TIL clusters, the patients and two rows of markers (protein and RNA);
- its frequency and rank in the deep TCRβ-seq repertoires.

Besides the clonotypes of the prioritization tables, every TIL clonotype of the patient can be selected (those with at least two TIL cells are listed; all can be found with the filter box). For the ones found in the experimental conditions, the activation statistics are recomputed with the rule of section 5; the script checks that this reproduces the values of the tables. These entries are labelled "Non prioritized" and exist only in the viewer.

### 6.3 Expansion of the TIL clonotypes of a cluster

`cxcl13-expansion.Rmd` takes the clonotypes with cells in one cluster of the TIL object (`CD4_CXCL13` of `clusters_old`; set `cluster_column` and `cluster_name` to use another one) in the five patients, and classifies them with the deep TCRβ-seq data: **Expanded** when the frequency in PE (primary expansion) is higher than in Bulk. It then compares the cells of expanded and non-expanded clonotypes (Wilcoxon test, Seurat `FindMarkers`) for genes (RNA) and proteins (CITE-seq antibodies), with all the patients together and for each patient. Outputs in `patients/cd4-cxcl13-expansion/`: clonotype list, counts per patient, DEG/DEP tables and `CD4_CXCL13_expansion_report.html`.

## Reproducibility notes

- Seeds are fixed in the scripts (`set.seed(123)`) and in Python (NumPy, `random`) by `setup/setup-m1.R`.
- The analyses were run in the Docker image described in section 1; `renv.lock` reproduces its package versions on a conda R 4.5.2. Small numerical differences (e.g. UMAP coordinates) are possible between the two because of the BLAS library and the locale.
- The Docker image was rebuilt after the deconvolution analyses were run; the packages that come from the image rather than from the user library (e.g. `readr`, `forcats`, `ape`) are recorded with the versions of the rebuilt image.

## Contact

Questions about the code: Endika Prieto (endikaprieto@vhio.net).

## License

The code is released under the [MIT License](LICENSE). Data are not covered by it.
