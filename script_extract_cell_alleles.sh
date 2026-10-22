#!/usr/bin/env bash
#SBATCH -p short
#SBATCH --job-name=extract_cell_alleles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=5:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -euo pipefail

ROOT="$(pwd)"
DEDUP="no"
PATIENTS=()

show_usage() {
    cat <<'USAGE'

Usage:
  Single patient:
        sbatch script_extract_cell_alleles.sh --root /path/to/deconv/raw_data --patients pt3385 --dedup no
        e.g. sbatch script_extract_cell_alleles.sh --root /mnt/petasan_immuno/raw_data/itag/collabs/sturcotte/scSeq/deconv/raw_data --patients pt3385 --dedup no
  Multiple patients:
        sbatch script_extract_cell_alleles.sh --root /path/to/deconv/raw_data --patients "pt600280 pt600311" --dedup no
  All pt* patient directories:
        sbatch script_extract_cell_alleles.sh --root /path/to/deconv/raw_data --dedup no

Options:
  --root DIR        Directory containing pt* patient folders (default: current directory).
  --patients LIST   Space-separated patient IDs, quoted or unquoted.
                    If omitted, process all pt* directories in --root.
  --dedup yes|no    yes: keep only the first matching row per cell barcode.
                    no: retain all matching rows (default; preserves TRA/TRB entries).
  -h, --help        Display this help message.

Input structure:
  ROOT/ptXXXX/CONDITION/10x-sample-sheet.tsv
  ROOT/ptXXXX/CONDITION/imgt/*.cell_alleles.tsv

Output:
  ROOT/ptXXXX/CONDITION/imgt/*.cell_alleles_filtered.tsv
  The name 'PeptidePool' is changed to 'Pept_pool' in output filenames.

USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --root)
            if [[ $# -lt 2 || "$2" == --* ]]; then
                echo "[ERROR] --root requires a directory." >&2
                show_usage >&2
                exit 1
            fi
            ROOT="$2"
            shift 2
            ;;
        --patients)
            shift
            if [[ $# -eq 0 || "$1" == --* ]]; then
                echo "[ERROR] --patients requires at least one patient ID." >&2
                show_usage >&2
                exit 1
            fi
            while [[ $# -gt 0 && "$1" != --* ]]; do
                read -r -a parsed_patients <<< "$1"
                PATIENTS+=("${parsed_patients[@]}")
                shift
            done
            ;;
        --dedup)
            if [[ $# -lt 2 || "$2" == --* ]]; then
                echo "[ERROR] --dedup requires yes or no." >&2
                show_usage >&2
                exit 1
            fi
            DEDUP="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            show_usage >&2
            exit 1
            ;;
    esac
done

if [[ ! -d "$ROOT" ]]; then
    echo "[ERROR] Root directory does not exist: $ROOT" >&2
    exit 1
fi
ROOT="$(cd "$ROOT" && pwd -P)"

if [[ "$DEDUP" != "yes" && "$DEDUP" != "no" ]]; then
    echo "[ERROR] --dedup must be yes or no." >&2
    exit 1
fi

if [[ ${#PATIENTS[@]} -eq 0 ]]; then
    mapfile -t PATIENTS < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d -name 'pt*' -printf '%f\n' | sort)
fi

if [[ ${#PATIENTS[@]} -eq 0 ]]; then
    echo "[ERROR] No patient directories found under $ROOT" >&2
    exit 1
fi

echo "[INFO] ROOT=$ROOT"
echo "[INFO] PATIENTS=${PATIENTS[*]}"
echo "[INFO] DEDUP=$DEDUP"

processed=0
warnings=0

for patient in "${PATIENTS[@]}"; do
    patient_dir="$ROOT/$patient"
    if [[ ! -d "$patient_dir" ]]; then
        echo "[WARN] Patient directory does not exist: $patient_dir" >&2
        warnings=$((warnings + 1))
        continue
    fi

    found_sheet=0
    while IFS= read -r sheet; do
        found_sheet=1
        cond_dir="$(dirname "$sheet")"
        cond="$(basename "$cond_dir")"

        if [[ ! -d "$cond_dir/imgt" ]]; then
            echo "[WARN] $patient/$cond: missing imgt directory" >&2
            warnings=$((warnings + 1))
            continue
        fi

        # --- Verify the sample sheet contains a CELL barcode column ---
        sheet_col="$(awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i ~ /^CELL/) {print i; exit}}' "$sheet")"
        if [[ -z "$sheet_col" ]]; then
            echo "[WARN] $patient/$cond: no CELL column in sample sheet" >&2
            warnings=$((warnings + 1))
            continue
        fi

        mapfile -t alleles < <(find "$cond_dir/imgt" -maxdepth 1 -type f -name '*.cell_alleles.tsv' | sort)
        if [[ ${#alleles[@]} -eq 0 ]]; then
            echo "[WARN] $patient/$cond: no *.cell_alleles.tsv files found" >&2
            warnings=$((warnings + 1))
            continue
        fi

        sheet_unique="$(awk -F '\t' -v c="$sheet_col" 'NR>1 && $c!="" {print $c}' "$sheet" | sort -u | wc -l)"

        for allele in "${alleles[@]}"; do
            allele_base="$(basename "$allele" .cell_alleles.tsv)"
            if [[ "$allele_base" == *"PeptidePool"* ]]; then
                allele_base="${allele_base/PeptidePool/Pept_pool}"
            fi
            out="$cond_dir/imgt/${allele_base}.cell_alleles_filtered.tsv"

            # --- Filter allele records by barcodes from the sample sheet ---
            if [[ "$DEDUP" == "yes" ]]; then
                awk -F '\t' -v OFS='\t' -v c="$sheet_col" '
                    NR==FNR {if (FNR>1 && $c!="") keep[$c]=1; next}
                    FNR==1 {print; next}
                    ($1 in keep) && !seen[$1]++ {print}
                ' "$sheet" "$allele" > "$out"
            else
                awk -F '\t' -v OFS='\t' -v c="$sheet_col" '
                    NR==FNR {if (FNR>1 && $c!="") keep[$c]=1; next}
                    FNR==1 {print; next}
                    ($1 in keep) {print}
                ' "$sheet" "$allele" > "$out"
            fi

            in_rows=$(( $(wc -l < "$allele") - 1 ))
            out_rows=$(( $(wc -l < "$out") - 1 ))
            echo "[OK] $patient/$cond: sheet_unique=$sheet_unique, out_rows=$out_rows, in_rows=$in_rows -> $out"
            processed=$((processed + 1))
        done
    done < <(find "$patient_dir" -mindepth 2 -maxdepth 2 -type f -name '10x-sample-sheet.tsv' | sort)

    if [[ $found_sheet -eq 0 ]]; then
        echo "[WARN] $patient: no 10x-sample-sheet.tsv files found" >&2
        warnings=$((warnings + 1))
    fi
done

echo "[INFO] Processed: $processed .cell_alleles.tsv files"
echo "[INFO] Warnings: $warnings"
