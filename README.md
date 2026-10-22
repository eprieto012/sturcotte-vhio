# Identification and prioritization of candidate antigen-responsive TCR clonotypes following APC co-culture

## CHUM-MED x VHIO-AG

This README summarises the execution order for identifying and prioritising TCR clonotypes with increased activation after APC co-culture. The complete workflow applies to **pt3385, pt600280, pt600311, pt3328 and pt3513**.

## Main paths

```Shell
RAW_DATA=/mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/raw_data
SEURAT=/mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/seurat
SCRIPTS=/home/endikaprieto@vhio.org/sturcotte
PATIENT_OUTPUTS=/home/endikaprieto@vhio.org/sturcotte/patients
```

## Execution order

### 1. Generate MiXCR sample sheets from Seurat cell barcodes

**Script:**
`/home/endikaprieto@vhio.org/sturcotte/deconv-cell-barcode.Rmd`

Run the R Markdown file in RStudio on the VHIO OOD environment. It extracts Seurat cell barcodes and creates one `10x-sample-sheet.tsv` per
condition.

**Output example for `pt3385`:**

`/mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/seurat/pt3385/10x_sample_sheets/Unpulsed/10x-sample-sheet.tsv`

### 2. Concatenate FASTQ files from different sequencing lanes, when needed

**Script:**
`/home/endikaprieto@vhio.org/sturcotte/script_fastq_concat.sh`

Use this when a condition contains multiple `Lane*` directories. The script concatenates R1/R2 and, when present, I1/I2 files into the condition’s `Data` directory.

**Example for `pt3385`:**

```bash
sbatch /home/endikaprieto@vhio.org/sturcotte/script_fastq_concat.sh \
  /mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/raw_data/pt3385
```

The script overwrites the merged FASTQ files in `Data`; run it only when concatenation is needed.

### 3. Check the sample-sheet and FASTQ condition mapping

**Script:**
`/home/endikaprieto@vhio.org/sturcotte/check_samplesheet_vs_fastq.py`

Run the check before MiXCR to verify the condition/sample mapping. The example checks up to one million reads.

```bash
cd /mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/raw_data
python3 /home/endikaprieto@vhio.org/sturcotte/check_samplesheet_vs_fastq.py \
  --max-reads 1000000
```

### 4. Run MiXCR and infer TCR alleles

**Script:** `/home/endikaprieto@vhio.org/sturcotte/script_mixcr_imgt.sh`

Run the script once per condition. It uses the `10x-sc-xcr-vdj` preset and builds or updates an IMGT-derived reference library from current IMGT reference annotations. These updated suballele annotations support high-resolution TCR gene assignment and reconstruction. Docker and a working MiXCR licence are required.

**Example for `pt3385`, condition `Unpulsed`:**

```bash
sbatch /home/endikaprieto@vhio.org/sturcotte/script_mixcr_imgt.sh \
  /mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/raw_data/pt3385/Unpulsed \
  /mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/seurat/pt3385/10x_sample_sheets/Unpulsed/10x-sample-sheet.tsv \
  --no-skip \
  --no-validate-only \
  --check-imgt-updates
```

Repeat for `CD3CD28`, `Pept_pool` and each available `R` condition, replacing the condition in both paths. Outputs are written under `RAW_DATA/pt3385/<CONDITION>/imgt/` and `RAW_DATA/pt3385/<CONDITION>/validation/`.

### 5. Filter MiXCR assignments to the Seurat cell barcodes

**Script:**
`/home/endikaprieto@vhio.org/sturcotte/script_extract_cell_alleles.sh`

This script retains MiXCR records whose cell barcodes occur in the condition’s sample sheet. Use `--dedup no` to preserve all matching rows, including separate TRA/TRB entries.

**Example for `pt3385`:**

```bash
sbatch /home/endikaprieto@vhio.org/sturcotte/script_extract_cell_alleles.sh \
  --root /mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/raw_data \
  --patients pt3385 \
  --dedup no
```

**Output example:**

`RAW_DATA/pt3385/<CONDITION>/imgt/*.cell_alleles_filtered.tsv`

### 6. Combine Cell Ranger and MiXCR data and run the TCR analysis

Run the relevant R Markdown file in RStudio on the VHIO OOD environment, after the filtered MiXCR files have been generated. The analysis combines Cell Ranger and MiXCR annotations, performs QC and clonotype analysis, prioritises candidate TCRs, adds available bulk-repertoire context and generates plots.

**Script:** `/home/endikaprieto@vhio.org/sturcotte/deconv-analysis.Rmd`

Libraries and custom functions are loaded through `~/sturcotte/setup/setup-m1.R`. If the Rmagic Python environment has not yet been created, run
`~/sturcotte/setup/install-Rmagic-venv.Rmd` first.

**Main outputs for `pt3385`:**

- Prioritisation tables: `~/sturcotte/patients/pt3385/deconv/tcr/`
- UMAP and other plots: `~/sturcotte/patients/pt3385/deconv/plots/`
- Main candidate tables: `~/sturcotte/patients/pt3385/deconv/tcr/prioritized_cd4_tcrs_for_validation.tsv` and
  `/sturcotte/patients/pt3385/deconv/tcr/prioritized_cd8_tcrs_for_validation.tsv`

## Interpreting the prioritisation tables and UMAPs

### Table columns

| Column                                           | Interpretation                                                                                                                                                                                                          |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Prioritize`                                   | `TRUE` marks a positive candidate selected for  functional validation. `FALSE` means it was not selected as a positive candidate; this also includes putative negative controls.                                   |
| `Tier` / `Level`                             | Summarise the evidence for antigen responsiveness and the selection stringency, or identify a putative negative-control class.                                                                                          |
| `Condition`                                    | Stimulation condition with significant activation, such as an individual `R` condition or `Pept_pool`. Blank for negative controls.                                                                                 |
| `Unpulsed`, `CD3CD28`, `Pept_pool`, `R*` | Activated cells / total clonotype cells, with the activation percentage. A cell is classified as activated when its gene-expression-based UCell activation score exceeds the 99th percentile of scores in `Unpulsed`. |
| `*_padj`                                       | Benjamini–Hochberg-adjusted P value for increased activation in that condition; values below 0.01 are significant.                                                                                                     |
| `*_cluster`, `*_phenotype`                   | Most frequent Seurat cluster or functional cell-state annotation for that clonotype and condition, with the proportion of cells.                                                                                        |
| Repertoire frequency / rank                      | Clonotype abundance and rank in the available bulk repertoires (e.g. Tumor, Blood, PE, Bulk and TIDAL01). These data provide context but do not prove antigen specificity.                                             |

### Tier 1 vs Tier 2: where is the response?

| Tier                                           | Definition                                                                                                                                 |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| **Tier 1 — deconvoluted response**      | Significant activation (`padj < 0.01`) in at least one individual `R` condition, with or without a significant `Pept_pool` response. |
| **Tier 2 — peptide-pool response only** | Significant activation in `Pept_pool`, with no significant activation in any individual `R` condition.                                 |

### L1 vs L2: which positive candidates are selected?

All prioritised positive clonotypes must have at least one cell in `Unpulsed` and less than 5% basal activation. The condition used to meet an activation-percentage threshold must contain at least six clonotype cells.

| Category            | Selection criteria                                                                                                                                          |
| ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Tier 1-L1** | At least one `R` condition has `padj < 0.01` and more than 60% activated cells.                                                                         |
| **Tier 1-L2** | At least one `R` condition is significant, but none meets L1; `Pept_pool` has at least 50% activated cells. `Pept_pool` significance is not required. |
| **Tier 2-L1** | `Pept_pool` is significant and has more than 60% activated cells; no `R` condition is significant.                                                      |
| **Tier 2-L2** | `Pept_pool` is significant and has 50–60% activated cells; no `R` condition is significant. Used only if there are no Tier 2–L1 candidates.           |

L2 is a rescue category, not automatically a weak response. A Tier 1–L2 clonotype can respond strongly to `Pept_pool` but fail the stricter L1 requirement in an individual `R` condition.

### Putative negative controls (C1–C3)

Negative controls have no significant antigen-induced activation and are marked `Tier = Negative control`, `Prioritize = FALSE`, with a blank
`Condition`.

| Control                                                | Selection criteria                                                                                                                                                                                    |
| ------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **C1 — strongest control evidence**             | Detected in `Unpulsed`; more than 50% activation among more than five `CD3CD28` cells; less than 2% activation in every detected `Pept_pool`/`R` condition, with at least one such condition. |
| **C2 — lower CD3CD28 cell count**               | Same requirements as C1, but with 1–5 `CD3CD28` cells.                                                                                                                                             |
| **C3 — not detected after antigen stimulation** | Absent from all `Pept_pool`/`R` conditions, with at least two cells in `Unpulsed` and/or `CD3CD28`.                                                                                           |

Clonotypes detected in both control conditions are preferred. Controls are selected in the order C1 → C2 → C3 until three are obtained. These are putative controls, not experimentally verified non-reactive TCRs; C3 may reflect limited sampling.

### UMAP plots

IC11 displays up to 12 prioritised positive TCRs, following the table order; the complete candidate list remains in the table. Each plot
highlights cells belonging to one clonotype, coloured by Seurat cluster, while other cells are grey. Panels correspond to `Unpulsed`, `CD3CD28`,
`Pept_pool` and individual `R` conditions.
