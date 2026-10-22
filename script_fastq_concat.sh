#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=fastq_concat
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=endikaprieto@vhio.net
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=5:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

# Usage:
# sbatch script_fastq_concat.sh /path/to/folder

set -euo pipefail
ROOT_DIR="${1:-}"
if [[ -z "${ROOT_DIR}" ]]; then
	echo "Usage: $0 /path/to/pt600280" >&2
	exit 1
fi
if [[ ! -d "${ROOT_DIR}" ]]; then
	echo "Root directory does not exist: ${ROOT_DIR}" >&2
	exit 1
fi
processed=0
skipped=0
for BASE_DIR in "${ROOT_DIR}"/*; do
	[[ -d "${BASE_DIR}" ]] || continue
	SAMPLE_NAME="$(basename "${BASE_DIR}")"
	OUT_DIR="${BASE_DIR}/Data"
	R1_OUT="${OUT_DIR}/${SAMPLE_NAME}_R1.fastq.gz"
	R2_OUT="${OUT_DIR}/${SAMPLE_NAME}_R2.fastq.gz"
	I1_OUT="${OUT_DIR}/${SAMPLE_NAME}_I1.fastq.gz"
	I2_OUT="${OUT_DIR}/${SAMPLE_NAME}_I2.fastq.gz"
	mapfile -t R1_FILES < <(find "${BASE_DIR}"/Lane* -maxdepth 2 -type f -name "*.fastq.gz" 2>/dev/null | grep -E "_R1(_|\.)" | sort)
	mapfile -t R2_FILES < <(find "${BASE_DIR}"/Lane* -maxdepth 2 -type f -name "*.fastq.gz" 2>/dev/null | grep -E "_R2(_|\.)" | sort)
	mapfile -t I1_FILES < <(find "${BASE_DIR}"/Lane* -maxdepth 2 -type f -name "*.fastq.gz" 2>/dev/null | grep -E "_I1(_|\.)" | sort)
	mapfile -t I2_FILES < <(find "${BASE_DIR}"/Lane* -maxdepth 2 -type f -name "*.fastq.gz" 2>/dev/null | grep -E "_I2(_|\.)" | sort)
	if [[ ${#R1_FILES[@]} -eq 0 && ${#R2_FILES[@]} -eq 0 ]]; then
		echo "Skipping ${SAMPLE_NAME}: no FASTQ files found in Lane*"
		((skipped+=1))
		continue
	fi
	if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
		echo "Skipping ${SAMPLE_NAME}: missing R1 or R2 files" >&2
		((skipped+=1))
		continue
	fi
	mkdir -p "${OUT_DIR}"
	echo "Processing ${SAMPLE_NAME}"
	echo "  Concatenating ${#R1_FILES[@]} R1 files -> ${R1_OUT}"
	cat "${R1_FILES[@]}" > "${R1_OUT}"
	echo "  Concatenating ${#R2_FILES[@]} R2 files -> ${R2_OUT}"
	cat "${R2_FILES[@]}" > "${R2_OUT}"
	if [[ ${#I1_FILES[@]} -gt 0 ]]; then
		echo "  Concatenating ${#I1_FILES[@]} I1 files -> ${I1_OUT}"
		cat "${I1_FILES[@]}" > "${I1_OUT}"
	else
		echo "  I1 not found, skipping I1 output"
	fi
	if [[ ${#I2_FILES[@]} -gt 0 ]]; then
		echo "  Concatenating ${#I2_FILES[@]} I2 files -> ${I2_OUT}"
		cat "${I2_FILES[@]}" > "${I2_OUT}"
	else
		echo "  I2 not found, skipping I2 output"
	fi
	echo "  Done: ${OUT_DIR}"
	((processed+=1))
done
echo "Summary:"
echo "  Processed folders: ${processed}"
echo "  Skipped folders: ${skipped}"