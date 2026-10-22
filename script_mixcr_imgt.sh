#!/bin/bash
#SBATCH -p highmem
#SBATCH --job-name=mixcr_imgt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=endikaprieto@vhio.net
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=300gb
#SBATCH --time=96:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "$0")"

show_usage() {
    cat <<'EOF'

Usage:
  Single patient folder:
        sbatch script_mixcr_imgt.sh /path/to/patient_folder /path/to/10x-sample-sheet.tsv_or_condition_root --no-skip --no-validate-only --check-imgt-updates

EOF
}

run_batch_mode() {
    local root_dir="$1"
    local input_sample_sheet="$2"
    local skip_align="$3"
    local validate_only="$4"
    local imgt_policy="$5"
    local processed=0
    local skipped=0
    local failed=0
    local sample_dir=""
    local sample_name=""

    if [ ! -d "${root_dir}" ]; then
        echo "ERROR: Root directory does not exist: ${root_dir}" >&2
        return 1
    fi

    if [ ! -e "${input_sample_sheet}" ]; then
        echo "ERROR: Sample sheet path does not exist: ${input_sample_sheet}" >&2
        return 1
    fi

    if [ -f "${input_sample_sheet}" ] && [ ! -s "${input_sample_sheet}" ]; then
        echo "ERROR: Sample sheet file is empty: ${input_sample_sheet}" >&2
        return 1
    fi

    echo "Running batch mode"
    echo "Root directory: ${root_dir}"
    echo "Sample sheet: ${input_sample_sheet}"
    echo "IMGT policy: ${imgt_policy}"

    for sample_dir in "${root_dir}"/*; do
        [ -d "${sample_dir}" ] || continue
        sample_name="$(basename "${sample_dir}")"

        if [ ! -d "${sample_dir}/Data" ] && [ ! -d "${sample_dir}/data" ]; then
            echo "Skipping ${sample_name}: missing Data/data folder"
            ((skipped+=1))
            continue
        fi

        echo ""
        echo "--- Running patient ${sample_name} ---"
        if bash "${SCRIPT_PATH}" --single "${sample_dir}" "${input_sample_sheet}" "${skip_align}" "${validate_only}" "${imgt_policy}"; then
            ((processed+=1))
        else
            echo "FAILED ${sample_name}" >&2
            ((failed+=1))
        fi
    done

    echo ""
    echo "Batch summary:"
    echo "  Processed: ${processed}"
    echo "  Skipped: ${skipped}"
    echo "  Failed: ${failed}"

    if [ "${failed}" -gt 0 ]; then
        return 1
    fi
}

# configuration
LICENSE="E-VRKEGFHFHNBQIFCUJZOTZPNQISODMMGKRROTCOLBAOTTAOKA"
DEFAULT_IMAGE="ghcr.io/milaboratory/mixcr/mixcr:4.7.0-235-develop"
IMAGE="${MIXCR_IMAGE_OVERRIDE:-${DEFAULT_IMAGE}}"

if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
    show_usage
    exit 0
fi

if [ "${1:-}" = "--batch" ]; then
    ROOT_DIR="${2:-}"
    INPUT_SAMPLE_SHEET="${3:-}"
    SKIP_ALIGN="${4:---no-skip}"
    VALIDATE_ONLY="${5:---no-validate-only}"
    IMGT_POLICY="${6:---check-imgt-updates}"

    if [ -z "${ROOT_DIR}" ] || [ -z "${INPUT_SAMPLE_SHEET}" ]; then
        show_usage
        exit 1
    fi

    run_batch_mode "${ROOT_DIR}" "${INPUT_SAMPLE_SHEET}" "${SKIP_ALIGN}" "${VALIDATE_ONLY}" "${IMGT_POLICY}"
    exit $?
fi

if [ "${1:-}" = "--single" ]; then
    shift
fi

WORK_DIR="${1:-}"
INPUT_SAMPLE_SHEET="${2:-}"
SKIP_ALIGN="${3:---no-skip}"
VALIDATE_ONLY="${4:---no-validate-only}"
IMGT_POLICY="${5:---check-imgt-updates}"

if [ -z "${WORK_DIR}" ] || [ -z "${INPUT_SAMPLE_SHEET}" ]; then
    show_usage
    exit 1
fi

case "${IMGT_POLICY}" in
    --if-missing-imgt|--check-imgt-updates|--refresh-imgt)
        ;;
    *)
        echo "ERROR: Invalid IMGT policy: ${IMGT_POLICY}" >&2
        show_usage
        exit 1
        ;;
esac

if [ ! -d "${WORK_DIR}" ]; then
    echo "ERROR: WORK_DIR does not exist: ${WORK_DIR}" >&2
    exit 1
fi

if [ ! -e "${INPUT_SAMPLE_SHEET}" ]; then
    echo "ERROR: Sample sheet path does not exist: ${INPUT_SAMPLE_SHEET}" >&2
    exit 1
fi

if [ -f "${INPUT_SAMPLE_SHEET}" ] && [ ! -s "${INPUT_SAMPLE_SHEET}" ]; then
    echo "ERROR: Sample sheet file is empty: ${INPUT_SAMPLE_SHEET}" >&2
    exit 1
fi

LIB_DIR="${WORK_DIR}/libraries"
LIB_NAME="imgt_hsa_tr_latest"
IMGT_ZIP_URL="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory.zip"
SAMPLE_SHEET="${WORK_DIR}/10x-sample-sheet.tsv"
SAMPLES=""
PRESET_DEFAULT="10x-sc-xcr-vdj"
PRESET="${PRESET_DEFAULT}"
CONSENSUS_STAGE="${MIXCR_CONSENSUS_STAGE:-4}"
ALIGN_INPUT_MODE="${MIXCR_ALIGN_INPUT:-consensus}"
LIBRARY_MODE="${MIXCR_LIBRARY_MODE:-custom}"
USE_SAMPLE_SHEET="${MIXCR_USE_SAMPLE_SHEET:-1}"
SKIP_REFINETAGS="${MIXCR_SKIP_REFINETAGS:-0}"
DISABLE_REFINE_POSTFILTER="${MIXCR_DISABLE_REFINE_POSTFILTER:-0}"
MIXCR_THREADS="${MIXCR_THREADS:-${SLURM_CPUS_PER_TASK:-40}}"

case "${CONSENSUS_STAGE}" in
    1|2|3|4)
        ;;
    *)
        echo "ERROR: MIXCR_CONSENSUS_STAGE must be one of: 1,2,3,4 (got '${CONSENSUS_STAGE}')" >&2
        exit 1
        ;;
esac

case "${ALIGN_INPUT_MODE}" in
    consensus|refined)
        ;;
    *)
        echo "ERROR: MIXCR_ALIGN_INPUT must be 'consensus' or 'refined' (got '${ALIGN_INPUT_MODE}')" >&2
        exit 1
        ;;
esac

case "${LIBRARY_MODE}" in
    custom|default)
        ;;
    *)
        echo "ERROR: MIXCR_LIBRARY_MODE must be 'custom' or 'default' (got '${LIBRARY_MODE}')" >&2
        exit 1
        ;;
esac

case "${USE_SAMPLE_SHEET}" in
    0|1)
        ;;
    *)
        echo "ERROR: MIXCR_USE_SAMPLE_SHEET must be 0 or 1 (got '${USE_SAMPLE_SHEET}')" >&2
        exit 1
        ;;
esac

case "${SKIP_REFINETAGS}" in
    0|1)
        ;;
    *)
        echo "ERROR: MIXCR_SKIP_REFINETAGS must be 0 or 1 (got '${SKIP_REFINETAGS}')" >&2
        exit 1
        ;;
esac

case "${DISABLE_REFINE_POSTFILTER}" in
    0|1)
        ;;
    *)
        echo "ERROR: MIXCR_DISABLE_REFINE_POSTFILTER must be 0 or 1 (got '${DISABLE_REFINE_POSTFILTER}')" >&2
        exit 1
        ;;
esac
# functions
log_section() {
    echo ""
    echo "============================================================================="
    echo "$1"
    echo "============================================================================="
}
run_mixcr() {
    docker run --rm \
      -e MI_LICENSE="${LICENSE}" \
      -v "${LIB_DIR}":/root/.mixcr/libraries:rw \
      -v "${WORK_DIR}":/work \
      "${IMAGE}" \
      mixcr "$@"
}

preset_exists() {
    local preset_name="$1"
    run_mixcr listPresets | awk 'NR > 1 {print $1}' | grep -Fxq "${preset_name}"
}

normalize_preset() {
    local requested="$1"

    if preset_exists "${requested}"; then
        echo "${requested}"
        return 0
    fi

    if [ "${requested}" = "10x-sc-xcr-vdj-gemx-v2" ]; then
        # Some MiXCR builds do not ship gemx-v2 preset; fall back to generic 10x VDJ.
        if preset_exists "10x-sc-xcr-vdj"; then
            echo "10x-sc-xcr-vdj"
            return 0
        fi
    fi

    echo "ERROR: Requested preset '${requested}' is not available in MiXCR image ${IMAGE}" >&2
    echo "Available related presets:" >&2
    run_mixcr listPresets | grep -E '10x-sc-xcr-vdj|gemx' >&2 || true
    return 1
}

utc_now() {
    date -u +"%Y-%m-%dT%H:%M:%SZ"
}

file_mtime_utc() {
    local path="$1"
    if [ ! -e "${path}" ]; then
        echo "missing"
        return 0
    fi
    date -u -d "@$(stat -c %Y "${path}")" +"%Y-%m-%dT%H:%M:%SZ"
}

fetch_imgt_last_modified_utc() {
    local header=""
    header="$(curl -fsSI "${IMGT_ZIP_URL}" | tr -d '\r' | awk 'BEGIN{IGNORECASE=1} /^Last-Modified:/{sub(/^Last-Modified:[[:space:]]*/, ""); print; exit}')" || return 1
    [ -n "${header}" ] || return 1
    date -u -d "${header}" +"%Y-%m-%dT%H:%M:%SZ" 2>/dev/null || echo "${header}"
}

read_metadata_value() {
    local metadata_file="$1"
    local key="$2"
    if [ ! -s "${metadata_file}" ]; then
        return 0
    fi
    awk -F'\t' -v k="${key}" '$1 == k {print $2; exit}' "${metadata_file}"
}

write_imgt_metadata() {
    local metadata_file="$1"
    local built_at_utc="$2"
    local remote_last_modified_utc="$3"
    local policy="$4"

    cat > "${metadata_file}" <<EOF
key	value
built_at_utc	${built_at_utc}
remote_last_modified_utc	${remote_last_modified_utc}
imgt_policy	${policy}
mixcr_image	${IMAGE}
imgt_zip_url	${IMGT_ZIP_URL}
EOF
}

write_imgt_status_report() {
    local status_file="$1"
    local action="$2"
    local reason="$3"
    local remote_last_modified_utc="$4"
    local remote_check_status="$5"
    local previous_remote_last_modified_utc="$6"

    cat > "${status_file}" <<EOF
key	value
checked_at_utc	$(utc_now)
update_action	${action}
update_reason	${reason}
imgt_policy	${IMGT_POLICY}
remote_last_modified_utc	${remote_last_modified_utc}
remote_check_status	${remote_check_status}
previous_remote_last_modified_utc	${previous_remote_last_modified_utc}
library_file	${LIB_DIR}/${LIB_NAME}.json.gz
library_file_mtime_utc	$(file_mtime_utc "${LIB_DIR}/${LIB_NAME}.json.gz")
library_metadata_file	${LIB_DIR}/${LIB_NAME}.metadata.tsv
mixcr_image	${IMAGE}
imgt_zip_url	${IMGT_ZIP_URL}
EOF
}

resolve_sample_sheet_path() {
    local work_dir="$1"
    local sample_sheet_input="$2"
    local condition_name=""
    local direct_candidate=""
    local filename=""
    local all_sheet_files=()
    local resolved=""

    condition_name="$(basename "${work_dir}")"

    if [ -f "${sample_sheet_input}" ] && [ -s "${sample_sheet_input}" ]; then
        echo "${sample_sheet_input}"
        return 0
    fi

    if [ -d "${sample_sheet_input}" ]; then
        for filename in "10x-sample-sheet.tsv" "10x_sample_sheet.tsv"; do
            direct_candidate="${sample_sheet_input}/${condition_name}/${filename}"
            if [ -s "${direct_candidate}" ]; then
                echo "${direct_candidate}"
                return 0
            fi
        done

        for filename in "10x-sample-sheet.tsv" "10x_sample_sheet.tsv"; do
            direct_candidate="${sample_sheet_input}/${filename}"
            if [ -s "${direct_candidate}" ]; then
                echo "${direct_candidate}"
                return 0
            fi
        done

        mapfile -t all_sheet_files < <(find "${sample_sheet_input}" -mindepth 2 -maxdepth 3 -type f \( -name "10x-sample-sheet.tsv" -o -name "10x_sample_sheet.tsv" \) | sort)

        if [ ${#all_sheet_files[@]} -eq 0 ]; then
            echo "ERROR: No sample sheet files found under ${sample_sheet_input}" >&2
            return 1
        fi

        resolved="$(python3 - "${condition_name}" "${all_sheet_files[@]}" <<'PY'
import difflib
import os
import re
import sys

condition = sys.argv[1]
paths = sys.argv[2:]

def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())

target = norm(condition)
folder_names = [os.path.basename(os.path.dirname(path)) for path in paths]
folder_norms = [norm(name) for name in folder_names]

for path, folder_norm in zip(paths, folder_norms):
    if folder_norm == target:
        print(path)
        raise SystemExit(0)

substring_matches = [path for path, folder_norm in zip(paths, folder_norms) if target in folder_norm or folder_norm in target]
if len(substring_matches) == 1:
    print(substring_matches[0])
    raise SystemExit(0)

scores = [difflib.SequenceMatcher(None, target, folder_norm).ratio() for folder_norm in folder_norms]
best_score = max(scores)
best_indexes = [idx for idx, score in enumerate(scores) if score == best_score]

if len(best_indexes) == 1 and best_score >= 0.55:
    print(paths[best_indexes[0]])
    raise SystemExit(0)

raise SystemExit(1)
PY
)" || true

        if [ -n "${resolved}" ] && [ -s "${resolved}" ]; then
            echo "${resolved}"
            return 0
        fi

        echo "ERROR: Could not resolve a unique sample sheet for condition '${condition_name}' under ${sample_sheet_input}" >&2
        echo "Available sample sheet candidates:" >&2
        printf '  %s\n' "${all_sheet_files[@]}" >&2
        return 1
    fi

    echo "ERROR: Sample sheet input must be a file or a directory: ${sample_sheet_input}" >&2
    return 1
}

initialize_samples_from_sheet() {
    local resolved_sheet_source=""
    local resolved_abs=""
    local target_abs=""

    resolved_sheet_source="$(resolve_sample_sheet_path "${WORK_DIR}" "${INPUT_SAMPLE_SHEET}")"
    resolved_abs="$(readlink -f "${resolved_sheet_source}")"
    target_abs="$(readlink -f "${SAMPLE_SHEET}" 2>/dev/null || true)"

    if [ "${resolved_abs}" != "${target_abs}" ]; then
        cp -f "${resolved_sheet_source}" "${SAMPLE_SHEET}"
    else
        echo "Sample sheet already in place: ${SAMPLE_SHEET}"
    fi

    echo "Resolved sample sheet for $(basename "${WORK_DIR}"): ${resolved_sheet_source}"

    if [ ! -s "${SAMPLE_SHEET}" ]; then
        echo "ERROR: Missing or empty sample sheet after resolution: ${SAMPLE_SHEET}" >&2
        return 1
    fi

    SAMPLES="$(awk 'NR > 1 && $1 != "" {print $1}' "${SAMPLE_SHEET}" | sort -u | xargs)"

    if [ -z "${SAMPLES}" ]; then
        echo "ERROR: No sample names found in ${SAMPLE_SHEET}" >&2
        return 1
    fi
}

select_fastq_file() {
    local read_label="$1"
    shift
    local files=("$@")
    local preferred=()
    local f

    if [ ${#files[@]} -eq 0 ]; then
        echo "ERROR: No ${read_label} FASTQ files found" >&2
        return 1
    fi

    for f in "${files[@]}"; do
        if [[ "$(basename "${f}")" =~ _${read_label}\.fastq\.gz$ ]]; then
            preferred+=("${f}")
        fi
    done

    if [ ${#preferred[@]} -eq 1 ]; then
        echo "${preferred[0]}"
        return 0
    fi

    if [ ${#files[@]} -eq 1 ]; then
        echo "${files[0]}"
        return 0
    fi

    echo "ERROR: Multiple ${read_label} FASTQ candidates found; keep only one merged file in Data/data" >&2
    printf '  %s\n' "${files[@]}" >&2
    return 1
}

prepare_consensus_from_raw() {
    local consensus_file="${WORK_DIR}/samples/sample_result.consensus.${CONSENSUS_STAGE}.mic"
    local consensus_exists=0
    local raw_data_dir_host=""
    local raw_data_dir_container=""
    local r1_files=()
    local r2_files=()
    local raw_r1_host=""
    local raw_r2_host=""
    local raw_r1_container=""
    local raw_r2_container=""

    if [ -s "${consensus_file}" ]; then
        consensus_exists=1
    fi

    if [ -d "${WORK_DIR}/Data" ]; then
        raw_data_dir_host="${WORK_DIR}/Data"
        raw_data_dir_container="/work/Data"
    elif [ -d "${WORK_DIR}/data" ]; then
        raw_data_dir_host="${WORK_DIR}/data"
        raw_data_dir_container="/work/data"
    else
        echo "Warning: Neither ${WORK_DIR}/Data nor ${WORK_DIR}/data exists; using default preset ${PRESET}" >&2
    fi

    if [ -n "${raw_data_dir_host}" ]; then
        mapfile -t r1_files < <(find "${raw_data_dir_host}" -maxdepth 1 -type f -name "*.fastq.gz" | grep -E "_R1(_|\.)" | sort || true)
        mapfile -t r2_files < <(find "${raw_data_dir_host}" -maxdepth 1 -type f -name "*.fastq.gz" | grep -E "_R2(_|\.)" | sort || true)

        if raw_r1_host="$(select_fastq_file "R1" "${r1_files[@]}")" && raw_r2_host="$(select_fastq_file "R2" "${r2_files[@]}")"; then
            raw_r1_container="${raw_data_dir_container}/$(basename "${raw_r1_host}")"
            raw_r2_container="${raw_data_dir_container}/$(basename "${raw_r2_host}")"
        else
            if [ "${consensus_exists}" -eq 1 ]; then
                echo "Warning: Could not resolve R1/R2 FASTQ files for chemistry autodetection; keeping preset ${PRESET}" >&2
            else
                echo "ERROR: Could not resolve unique R1/R2 FASTQ files in ${raw_data_dir_host}" >&2
                return 1
            fi
        fi

        if [ -n "${MIXCR_PRESET_OVERRIDE:-}" ]; then
            PRESET="${MIXCR_PRESET_OVERRIDE}"
            echo "Using MixCR preset override from MIXCR_PRESET_OVERRIDE: ${PRESET}"
        else
            PRESET="${PRESET_DEFAULT}"
            echo "Using default MixCR preset: ${PRESET}"
            echo "Tip: set MIXCR_PRESET_OVERRIDE to force another preset when needed"
        fi
        PRESET="$(normalize_preset "${PRESET}")"
        echo "Using MixCR preset: ${PRESET}"
    fi

    if [ "${consensus_exists}" -eq 1 ]; then
        echo "Existing consensus file found: ${consensus_file}"
        echo "Using consensus stage: ${CONSENSUS_STAGE}"
        echo "Using MixCR preset for align: ${PRESET}"
        return 0
    fi

    if [ -z "${raw_data_dir_host}" ] || [ -z "${raw_r1_host}" ] || [ -z "${raw_r2_host}" ]; then
        echo "ERROR: Missing raw FASTQ inputs and consensus file does not exist; cannot run mixcr analyze" >&2
        return 1
    fi

    log_section "STEP PRE: GENERATE sample_result.consensus.${CONSENSUS_STAGE}.mic FROM RAW FASTQ"
    echo "Raw data folder: ${raw_data_dir_host}"
    echo "Using R1: ${raw_r1_host}"
    echo "Using R2: ${raw_r2_host}"
    echo "Using consensus stage: ${CONSENSUS_STAGE}"
        if [ "${USE_SAMPLE_SHEET}" = "1" ]; then
                ANALYZE_SAMPLE_SHEET_ARGS=(--sample-sheet /work/10x-sample-sheet.tsv)
                echo "Using sample-sheet during analyze: yes"
        else
                ANALYZE_SAMPLE_SHEET_ARGS=()
                echo "Using sample-sheet during analyze: no"
        fi

    echo "Using MiXCR threads: ${MIXCR_THREADS}"
    run_mixcr analyze "${PRESET}" \
      --species hsa \
      --threads "${MIXCR_THREADS}" \
            "${ANALYZE_SAMPLE_SHEET_ARGS[@]}" \
      "${raw_r1_container}" \
      "${raw_r2_container}" \
      /work/samples/sample_result

    if [ ! -s "${consensus_file}" ]; then
        echo "ERROR: MiXCR analyze finished but selected consensus file was not created: ${consensus_file}" >&2
        echo "Tip: try MIXCR_CONSENSUS_STAGE=1 or MIXCR_CONSENSUS_STAGE=2" >&2
        return 1
    fi
}

build_imgt_library() {
    local lib_file="${LIB_DIR}/${LIB_NAME}.json.gz"
    local metadata_file="${LIB_DIR}/${LIB_NAME}.metadata.tsv"
    local status_file="${WORK_DIR}/validation/${LIB_NAME}.status.tsv"
    local remote_last_modified_utc="unknown"
    local previous_remote_last_modified_utc="unknown"
    local remote_check_status="not_checked"
    local update_action="reused"
    local update_reason="library_present"
    local should_rebuild=0
    local built_at_utc=""

    mkdir -p "${WORK_DIR}/imgt_ref/raw" "${WORK_DIR}/imgt_ref/prepared" "${LIB_DIR}" "${WORK_DIR}/validation"

    previous_remote_last_modified_utc="$(read_metadata_value "${metadata_file}" "remote_last_modified_utc")"
    if [ -z "${previous_remote_last_modified_utc}" ]; then
        previous_remote_last_modified_utc="unknown"
    fi

    if [ "${IMGT_POLICY}" != "--if-missing-imgt" ] || [ ! -s "${lib_file}" ]; then
        if remote_last_modified_utc="$(fetch_imgt_last_modified_utc)"; then
            remote_check_status="ok"
        else
            remote_last_modified_utc="unknown"
            remote_check_status="failed"
            echo "WARNING: Could not retrieve IMGT Last-Modified header from ${IMGT_ZIP_URL}" >&2
        fi
    fi

    if [ ! -s "${lib_file}" ]; then
        should_rebuild=1
        update_reason="library_missing"
    elif [ "${IMGT_POLICY}" = "--refresh-imgt" ]; then
        should_rebuild=1
        update_reason="forced_refresh"
    elif [ "${IMGT_POLICY}" = "--check-imgt-updates" ] && [ "${remote_check_status}" = "ok" ] && [ "${previous_remote_last_modified_utc}" != "${remote_last_modified_utc}" ]; then
        should_rebuild=1
        update_reason="remote_last_modified_changed"
    fi

    if [ "${should_rebuild}" -eq 0 ]; then
        echo "Using existing MiXCR IMGT library: ${lib_file}"
        if [ ! -s "${metadata_file}" ]; then
            write_imgt_metadata "${metadata_file}" "$(file_mtime_utc "${lib_file}")" "${remote_last_modified_utc}" "${IMGT_POLICY}"
        fi
        write_imgt_status_report "${status_file}" "${update_action}" "${update_reason}" "${remote_last_modified_utc}" "${remote_check_status}" "${previous_remote_last_modified_utc}"
        echo "IMGT status report: ${status_file}"
        return 0
    fi

    update_action="rebuilt"
    log_section "STEP 0: DOWNLOAD IMGT/V-QUEST REFERENCES AND BUILD MIXCR LIBRARY"
    echo "IMGT policy: ${IMGT_POLICY}"
    echo "Update reason: ${update_reason}"
    if [ "${remote_check_status}" = "ok" ]; then
        echo "Remote IMGT Last-Modified (UTC): ${remote_last_modified_utc}"
    fi

    curl -fsSL "${IMGT_ZIP_URL}" -o "${WORK_DIR}/imgt_ref/IMGT_V-QUEST_reference_directory.zip"
    rm -rf "${WORK_DIR}/imgt_ref/raw/IMGT_V-QUEST_reference_directory"
    unzip -q -o "${WORK_DIR}/imgt_ref/IMGT_V-QUEST_reference_directory.zip" -d "${WORK_DIR}/imgt_ref/raw"
    RAW_TR_DIR="$(find "${WORK_DIR}/imgt_ref/raw" -type d -path "*/Homo_sapiens/TR" | head -n 1)"
    if [ -z "${RAW_TR_DIR}" ]; then
        echo "ERROR: Could not find Homo_sapiens/TR inside IMGT ZIP"
        exit 1
    fi
    python3 - "${RAW_TR_DIR}" "${WORK_DIR}/imgt_ref/prepared" <<'PY'
import re
import sys
from pathlib import Path
raw = Path(sys.argv[1])
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
wanted = ["TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ"]
functional_values = {"F", "(F)", "[F]"}
def iter_fasta(path):
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            yield header, "".join(seq)
for stem in wanted:
    inp = raw / f"{stem}.fasta"
    outp = out / f"{stem}.fasta"
    if not inp.exists():
        raise SystemExit(f"Missing expected IMGT file: {inp}")
    n_in = 0
    n_out = 0
    seen = set()
    with open(outp, "w") as wh:
        for header, seq in iter_fasta(inp):
            n_in += 1
            fields = header.split("|")
            if len(fields) < 5:
                raise SystemExit(f"Unexpected IMGT header in {inp}: {header}")
            gene = fields[1].strip()
            func = fields[3].strip()
            partial = any("partial" in x.lower() for x in fields)
            if func not in functional_values:
                continue
            if partial:
                continue
            seq = re.sub(r"[\s\.\-]", "", seq).upper().replace("U", "T")
            bad = re.sub(r"[ACGT]", "", seq)
            if bad:
                raise SystemExit(f"Unexpected nucleotide(s) in {gene} from {inp}: {sorted(set(bad))}")
            if not seq or gene in seen:
                continue
            seen.add(gene)
            wh.write(f">{gene}\n")
            for i in range(0, len(seq), 80):
                wh.write(seq[i:i+80] + "\n")
            n_out += 1
    print(f"{stem}: wrote {n_out}/{n_in} functional non-partial records -> {outp}")
PY
    run_mixcr buildLibrary \
      --debug \
      --force-overwrite \
      --species hsa \
      --taxon-id 9606 \
      --chain TRA \
      --v-genes-from-fasta /work/imgt_ref/prepared/TRAV.fasta \
      --v-gene-feature VRegion \
      --j-genes-from-fasta /work/imgt_ref/prepared/TRAJ.fasta \
      --j-gene-feature JRegion \
      /root/.mixcr/libraries/${LIB_NAME}.TRA.json.gz \
      > "${WORK_DIR}/validation/${LIB_NAME}.TRA.buildLibrary.log" 2>&1
    run_mixcr buildLibrary \
      --debug \
      --force-overwrite \
      --species hsa \
      --taxon-id 9606 \
      --chain TRB \
      --v-genes-from-fasta /work/imgt_ref/prepared/TRBV.fasta \
      --v-gene-feature VRegion \
      --d-genes-from-fasta /work/imgt_ref/prepared/TRBD.fasta \
      --d-gene-feature DRegion \
      --j-genes-from-fasta /work/imgt_ref/prepared/TRBJ.fasta \
      --j-gene-feature JRegion \
      /root/.mixcr/libraries/${LIB_NAME}.TRB.json.gz \
      > "${WORK_DIR}/validation/${LIB_NAME}.TRB.buildLibrary.log" 2>&1
    run_mixcr mergeLibrary \
      --force-overwrite \
      /root/.mixcr/libraries/${LIB_NAME}.TRA.json.gz \
      /root/.mixcr/libraries/${LIB_NAME}.TRB.json.gz \
      /root/.mixcr/libraries/${LIB_NAME}.json.gz
    run_mixcr debugLibrary \
      /root/.mixcr/libraries/${LIB_NAME}.json.gz \
      > "${WORK_DIR}/validation/${LIB_NAME}.debugLibrary.txt" 2>&1

        built_at_utc="$(utc_now)"
        write_imgt_metadata "${metadata_file}" "${built_at_utc}" "${remote_last_modified_utc}" "${IMGT_POLICY}"
        write_imgt_status_report "${status_file}" "${update_action}" "${update_reason}" "${remote_last_modified_utc}" "${remote_check_status}" "${previous_remote_last_modified_utc}"

        echo "Built library: ${LIB_DIR}/${LIB_NAME}.json.gz"
        echo "Debug report: ${WORK_DIR}/validation/${LIB_NAME}.debugLibrary.txt"
        echo "Metadata: ${metadata_file}"
        echo "IMGT status report: ${status_file}"
}
# workflow
mkdir -p "${WORK_DIR}/imgt" "${WORK_DIR}/samples" "${WORK_DIR}/validation" "${LIB_DIR}"
initialize_samples_from_sheet
build_imgt_library
if [ "$VALIDATE_ONLY" == "--validate-only" ]; then
    log_section "VALIDATION ONLY MODE"
    echo "Skipping MixCR pipeline after library build/debug validation."
    exit 0
fi
if [ "$SKIP_ALIGN" != "--skip-align" ]; then
    prepare_consensus_from_raw
    if [ "${USE_SAMPLE_SHEET}" = "0" ] && [ "$(echo "${SAMPLES}" | wc -w)" -ne 1 ]; then
        echo "ERROR: MIXCR_USE_SAMPLE_SHEET=0 requires exactly one sample name in sample sheet" >&2
        echo "Detected samples: ${SAMPLES}" >&2
        exit 1
    fi
    ALIGN_INPUT_FILE=""
    if [ "${ALIGN_INPUT_MODE}" = "refined" ]; then
        ALIGN_INPUT_FILE="/work/samples/sample_result.refined.mic"
    else
        ALIGN_INPUT_FILE="/work/samples/sample_result.consensus.${CONSENSUS_STAGE}.mic"
    fi

    if [ ! -s "${WORK_DIR}${ALIGN_INPUT_FILE#/work}" ]; then
        echo "ERROR: Align input file not found or empty: ${WORK_DIR}${ALIGN_INPUT_FILE#/work}" >&2
        exit 1
    fi

    log_section "STEP 1: ALIGN with custom IMGT/V-QUEST-derived library"
    if [ "${LIBRARY_MODE}" = "default" ]; then
        ALIGN_LIBRARY_ARGS=()
        echo "Using MiXCR default built-in library"
    else
        ALIGN_LIBRARY_ARGS=(--library "${LIB_NAME}")
        echo "Using custom library: ${LIB_NAME}"
    fi
    if [ "${USE_SAMPLE_SHEET}" = "1" ]; then
        ALIGN_SAMPLE_SHEET_ARGS=(--sample-sheet /work/10x-sample-sheet.tsv)
        echo "Using sample-sheet during align: yes"
    else
        ALIGN_SAMPLE_SHEET_ARGS=()
        echo "Using sample-sheet during align: no"
    fi
    if [ "${DISABLE_REFINE_POSTFILTER}" = "1" ]; then
        ALIGN_PRESET_OVERRIDES=(-M refineTagsAndSort.parameters.postFilter=null)
        echo "Disabling refineTagsAndSort postFilter via align preset override"
    else
        ALIGN_PRESET_OVERRIDES=()
    fi
    echo "Using align input mode: ${ALIGN_INPUT_MODE}"
    echo "Using align input file: ${ALIGN_INPUT_FILE}"
        run_mixcr align \
            --force-overwrite \
            "${ALIGN_LIBRARY_ARGS[@]}" \
            --preset "${PRESET}" \
            "${ALIGN_PRESET_OVERRIDES[@]}" \
      --species hsa \
      --report /work/imgt/imgt_result.align.report.txt \
      --json-report /work/imgt/imgt_result.align.report.json \
      --save-output-file-names /work/imgt/imgt_result.align.list.tsv \
      "${ALIGN_SAMPLE_SHEET_ARGS[@]}" \
    "${ALIGN_INPUT_FILE}" \
      /work/imgt/imgt_result.alignments.vdjca
    echo "✓ Alignment complete"
    log_section "STEPS 2-6: PER-SAMPLE PROCESSING (refine → extend → assemble → assembleContigs → assembleCells)"
    for SAMPLE in ${SAMPLES}; do
            SAMPLE_ALIGN_INPUT=""
            SAMPLE_EXTEND_INPUT=""
            if [ "${USE_SAMPLE_SHEET}" = "1" ]; then
                SAMPLE_ALIGN_INPUT="/work/imgt/imgt_result.${SAMPLE}.alignments.vdjca"
            else
                SAMPLE_ALIGN_INPUT="/work/imgt/imgt_result.alignments.vdjca"
            fi
      echo ""
      echo "--- Processing ${SAMPLE} ---"
            if [ "${SKIP_REFINETAGS}" = "1" ]; then
                echo "Skipping refineTagsAndSort for ${SAMPLE} (MIXCR_SKIP_REFINETAGS=1)"
                SAMPLE_EXTEND_INPUT="${SAMPLE_ALIGN_INPUT}"
            else
                run_mixcr refineTagsAndSort \
                    --force-overwrite \
                    --report /work/imgt/imgt_result.${SAMPLE}.refine.report.txt \
                    --json-report /work/imgt/imgt_result.${SAMPLE}.refine.report.json \
                    "${SAMPLE_ALIGN_INPUT}" \
                    /work/imgt/imgt_result.${SAMPLE}.refined.vdjca
                SAMPLE_EXTEND_INPUT="/work/imgt/imgt_result.${SAMPLE}.refined.vdjca"
            fi
      run_mixcr extend \
        --force-overwrite \
        --report /work/imgt/imgt_result.${SAMPLE}.extend.report.txt \
        --json-report /work/imgt/imgt_result.${SAMPLE}.extend.report.json \
                "${SAMPLE_EXTEND_INPUT}" \
        /work/imgt/imgt_result.${SAMPLE}.extended.vdjca
      run_mixcr assemble \
        --force-overwrite \
        --report /work/imgt/imgt_result.${SAMPLE}.assemble.report.txt \
        --json-report /work/imgt/imgt_result.${SAMPLE}.assemble.report.json \
        /work/imgt/imgt_result.${SAMPLE}.extended.vdjca \
        /work/imgt/imgt_result.${SAMPLE}.clna
      run_mixcr assembleContigs \
        --force-overwrite \
        --report /work/imgt/imgt_result.${SAMPLE}.assembleContigs.report.txt \
        --json-report /work/imgt/imgt_result.${SAMPLE}.assembleContigs.report.json \
        /work/imgt/imgt_result.${SAMPLE}.clna \
        /work/imgt/imgt_result.${SAMPLE}.contigs.clns
      run_mixcr assembleCells \
        --force-overwrite \
        --report /work/imgt/imgt_result.${SAMPLE}.assembleCells.report.txt \
        --json-report /work/imgt/imgt_result.${SAMPLE}.assembleCells.report.json \
        /work/imgt/imgt_result.${SAMPLE}.contigs.clns \
        /work/imgt/imgt_result.${SAMPLE}.assembledCells.clns
      echo "✓ ${SAMPLE} complete"
    done
    echo "✓ All samples processed"
fi


# STEP 7: findAlleles (always run to ensure latest results)
log_section "STEP 7: FIND ALLELES with MixCR"

FIND_ALLELES_INPUTS=()
for SAMPLE in ${SAMPLES}; do
  FIND_ALLELES_INPUTS+=("/work/imgt/imgt_result.${SAMPLE}.contigs.clns")
done

docker run --rm \
  -e MI_LICENSE="${LICENSE}" \
  -v "${LIB_DIR}":/root/.mixcr/libraries:ro \
  -v "${WORK_DIR}":/work \
  "${IMAGE}" \
  mixcr findAlleles \
    --force-overwrite \
    --output-template '/work/imgt/{file_name}_alleles.clns' \
    --export-alleles-mutations /work/imgt/imgt_result.alleles_stats.tsv \
    "${FIND_ALLELES_INPUTS[@]}"

echo "✓ MixCR findAlleles complete"

# STEP 8: exportClones (for all samples)
log_section "STEP 8: EXPORT CLONES"

for SAMPLE in ${SAMPLES}; do
  echo "Exporting ${SAMPLE}..."
  docker run --rm \
    -e MI_LICENSE="${LICENSE}" \
    -v "${LIB_DIR}":/root/.mixcr/libraries:ro \
    -v "${WORK_DIR}":/work \
    "${IMAGE}" \
    mixcr exportClones \
      --force-overwrite \
      --drop-default-fields \
      --dont-split-files \
      --split-by-tags Cell \
      -tags Cell \
      -cellId \
      -cellGroup \
      -cloneId \
      -chains \
      -topChains \
      -readCount \
      -uniqueTagCount Molecule \
      -vHit \
      -dHit \
      -jHit \
      -cHit \
      -vHitsWithScore \
      -dHitsWithScore \
      -jHitsWithScore \
      -cHitsWithScore \
      -vBestIdentityPercent \
      -dBestIdentityPercent \
      -jBestIdentityPercent \
      -cBestIdentityPercent \
      -nFeature CDR3 \
      -aaFeature CDR3 \
      /work/imgt/imgt_result.${SAMPLE}.contigs_alleles.clns \
      /work/imgt/imgt_result.${SAMPLE}.cell_alleles.tsv
done

echo "✓ Export complete"

# validation
log_section "STEP 9-11: CREDIBILITY VALIDATION AND FILTERING"

WORK_DIR_ENV="${WORK_DIR}" python3 << 'VALIDATION_SCRIPT'
import csv
from collections import defaultdict
import os
import re

work_dir = os.environ["WORK_DIR_ENV"]
analysis_dir = f"{work_dir}/validation"
os.makedirs(analysis_dir, exist_ok=True)

sample_sheet = f"{work_dir}/10x-sample-sheet.tsv"
with open(sample_sheet, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    samples = sorted({row.get('Sample', '').strip() for row in reader if row.get('Sample', '').strip()})

if not samples:
    raise SystemExit(f"No samples found in sample sheet: {sample_sheet}")

# Segment-aware validation.
# findAlleles is used for V/D/J allele support when present in alleles_stats.tsv.
# Constant-region C alleles/suballeles are not expected to be represented reliably
# in alleles_stats.tsv, so TRAC/TRBC suballeles are summarized separately from
# bestCHit/allCHitsWithScore exported clone/cell assignments.
segment_hit_fields = {
    'V': ['bestVHit', 'allVHitsWithScore'],
    'D': ['bestDHit', 'allDHitsWithScore'],
    'J': ['bestJHit', 'allJHitsWithScore'],
    'C': ['bestCHit', 'allCHitsWithScore'],
}

best_hit_fields = {
    'V': ['bestVHit'],
    'D': ['bestDHit'],
    'J': ['bestJHit'],
    'C': ['bestCHit'],
}

all_hit_fields = {
    'V': ['allVHitsWithScore'],
    'D': ['allDHitsWithScore'],
    'J': ['allJHitsWithScore'],
    'C': ['allCHitsWithScore'],
}

tcr_gene_prefix_to_segment = {
    'TRAV': 'V', 'TRBV': 'V', 'TRDV': 'V', 'TRGV': 'V',
    'TRBD': 'D', 'TRDD': 'D',
    'TRAJ': 'J', 'TRBJ': 'J', 'TRDJ': 'J', 'TRGJ': 'J',
    'TRAC': 'C', 'TRBC': 'C', 'TRDC': 'C', 'TRGC': 'C',
}

missing_values = {'', '-', '.', 'NA', 'N/A', 'NaN', 'nan', 'None', 'null'}
removed_status = 'REMOVED_BECAUSE_NOT_REPRESENTED_IN_SOURCE_FILES'


def safe_int(value, default=0):
    try:
        if value is None:
            return default
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default


def classify_gene_segment(gene_name):
    gene_name = (gene_name or '').strip()
    for prefix, segment in tcr_gene_prefix_to_segment.items():
        if gene_name.startswith(prefix):
            return segment
    return 'OTHER'


def strip_hit_annotations(hit):
    hit = (hit or '').strip()
    if hit in missing_values:
        return ''

    # Remove MixCR score/identity annotations, e.g. TRBC1*01(123.4), keeping the allele.
    hit = re.split(r'\s*\(', hit, maxsplit=1)[0].strip()
    hit = hit.split()[0].strip()

    # Remove rare key=value wrappers if they appear in exported compound fields.
    if '=' in hit:
        hit = hit.split('=', 1)[-1].strip()

    return '' if hit in missing_values else hit


def clean_hit_to_allele(hit):
    return strip_hit_annotations(hit)


def clean_hit_to_gene(hit):
    allele = clean_hit_to_allele(hit)
    if not allele:
        return ''
    gene = allele.split('*', 1)[0].strip()
    return '' if gene in missing_values else gene


def split_hit_field_to_alleles(value):
    value = (value or '').strip()
    if value in missing_values:
        return []

    # MixCR multi-hit fields are commonly delimited by comma, semicolon or pipe.
    # Do not split on '/', because some TCR gene names contain '/', e.g. TRAV23/DV6.
    parts = re.split(r'[;,|]', value)
    alleles = []
    for part in parts:
        allele = clean_hit_to_allele(part)
        if allele:
            alleles.append(allele)
    return alleles


def split_hit_field_to_genes(value):
    genes = []
    for allele in split_hit_field_to_alleles(value):
        gene = clean_hit_to_gene(allele)
        if gene:
            genes.append(gene)
    return genes


def get_first_nonempty_field(row, field_candidates):
    for field in field_candidates:
        value = row.get(field, '')
        if value and value.strip() not in missing_values:
            return field, value
    return None, ''


def get_cell_key(row, row_index):
    # Prefer true cell identifiers when present; otherwise fall back to row index.
    for field in ['cellId', 'Cell', 'cell', 'tagValueCELL', 'tagValueCell', 'cellGroup']:
        value = row.get(field, '')
        if value and value.strip() not in missing_values:
            return value.strip()
    return f"row_{row_index}"


def confidence_from_pct(pct):
    return "HIGH" if pct >= 80 else ("MEDIUM" if pct >= 50 else "LOW")


def validation_status_for_gene(gene, credible_genes, marginal_genes):
    if gene in credible_genes:
        return 'YES'
    if gene in marginal_genes:
        return 'MAYBE'
    return 'NO'


print("STEP 9: Loading data...")

# Load alleles_stats generated by mixcr findAlleles.
alleles_data = {}
alleles_file = f"{work_dir}/imgt/imgt_result.alleles_stats.tsv"
with open(alleles_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene_name = row.get('geneName', '')
        allele_name = row.get('alleleName', '')

        if not gene_name or not allele_name:
            continue

        if gene_name not in alleles_data:
            alleles_data[gene_name] = {}

        alleles_data[gene_name][allele_name] = {
            'enough_info': row.get('enoughInfo', 'false').lower() == 'true',
            'status': row.get('status', ''),
            'clones': safe_int(row.get('clonesCount', 0)),
            'naives': safe_int(row.get('naivesCount', 0)),
            'total': safe_int(row.get('totalClonesCountForGene', 0)),
            'segment': classify_gene_segment(gene_name),
        }

# Extract real cell-level assignments from V, D, J and C hits.
# This is gene-level support used for credibility filtering.
print("STEP 10: Extracting real V/D/J/C gene assignments...")

real_data = defaultdict(lambda: {'total': 0, 'segments': set(), 'source_fields': set()})
counted_cell_gene_pairs = set()

# This independent C-level structure preserves allele/suballele suffixes, e.g. TRBC1*01.
c_allele_data = defaultdict(lambda: {
    'cells': set(),
    'rows': 0,
    'read_count': 0,
    'molecule_count': 0,
    'samples': set(),
    'source_fields': set(),
})
c_gene_cells = defaultdict(set)
c_gene_rows = defaultdict(int)
c_gene_source_fields = defaultdict(set)

for sample in samples:
    cell_file = f"{work_dir}/imgt/imgt_result.{sample}.cell_alleles.tsv"
    if os.path.exists(cell_file):
        with open(cell_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row_index, row in enumerate(reader, start=1):
                cell_key = get_cell_key(row, row_index)
                read_count = safe_int(row.get('readCount', 0))
                molecule_count = safe_int(
                    row.get('uniqueTagCount', row.get('uniqueTagCount Molecule', row.get('uniqueTagCount_Molecule', 0)))
                )

                for fallback_segment, field_candidates in segment_hit_fields.items():
                    source_field, value = get_first_nonempty_field(row, field_candidates)
                    if not source_field:
                        continue

                    for gene in split_hit_field_to_genes(value):
                        segment = classify_gene_segment(gene)
                        if segment == 'OTHER':
                            segment = fallback_segment

                        dedup_key = (sample, cell_key, gene)
                        if dedup_key not in counted_cell_gene_pairs:
                            real_data[gene]['total'] += 1
                            counted_cell_gene_pairs.add(dedup_key)

                        real_data[gene]['segments'].add(segment)
                        real_data[gene]['source_fields'].add(source_field)

                # Constant-region allele/suballele report.
                # Prefer bestCHit for quantitative counting to avoid overcounting all candidate hits.
                c_source_field, c_value = get_first_nonempty_field(row, best_hit_fields['C'])
                if not c_source_field:
                    # Fallback: if bestCHit is absent, use allCHitsWithScore but mark the source field.
                    c_source_field, c_value = get_first_nonempty_field(row, all_hit_fields['C'])

                if c_source_field:
                    for c_allele in split_hit_field_to_alleles(c_value):
                        c_gene = clean_hit_to_gene(c_allele)
                        if not c_gene or classify_gene_segment(c_gene) != 'C':
                            continue

                        c_cell_key = (sample, cell_key, c_gene)
                        c_allele_key = (c_gene, c_allele)
                        c_gene_cells[c_gene].add(c_cell_key)
                        c_gene_rows[c_gene] += 1
                        c_gene_source_fields[c_gene].add(c_source_field)

                        c_allele_data[c_allele_key]['cells'].add((sample, cell_key, c_gene, c_allele))
                        c_allele_data[c_allele_key]['rows'] += 1
                        c_allele_data[c_allele_key]['read_count'] += read_count
                        c_allele_data[c_allele_key]['molecule_count'] += molecule_count
                        c_allele_data[c_allele_key]['samples'].add(sample)
                        c_allele_data[c_allele_key]['source_fields'].add(c_source_field)
    else:
        print(f"WARNING: Missing exported clone file: {cell_file}")

# Categorize genes in alleles_stats by MixCR enoughInfo.
print("STEP 11: Generating credibility report...")

genes_with_info = {}
genes_without_info = {}

for gene, alleles in alleles_data.items():
    has_true = any(a['enough_info'] for a in alleles.values())
    if has_true:
        genes_with_info[gene] = alleles
    else:
        genes_without_info[gene] = alleles

# Generate credibility report.
credible_genes = set()
marginal_genes = set()
report_file = f"{analysis_dir}/credibility_report.tsv"
all_report_genes = sorted(set(alleles_data.keys()) | set(real_data.keys()))

with open(report_file, 'w') as f:
    f.write("Gene\tSegmentType\tCredible\tReason\tCells_in_Data\tValidation_Status\tEnoughInfo\tMixCR_Allele_Statuses\tSourceHitFields\n")

    for gene in all_report_genes:
        cell_count = real_data[gene]['total'] if gene in real_data else 0
        segment_values = set()
        if gene in real_data:
            segment_values.update(real_data[gene]['segments'])
        if gene in alleles_data:
            segment_values.update(a['segment'] for a in alleles_data[gene].values())
        segment_type = ",".join(sorted(segment_values)) if segment_values else classify_gene_segment(gene)

        mixcr_statuses = sorted(
            set(a['status'] for a in alleles_data.get(gene, {}).values() if a['status'])
        )
        mixcr_status = ",".join(mixcr_statuses) if mixcr_statuses else "NA"

        source_fields = sorted(real_data[gene]['source_fields']) if gene in real_data else []
        source_fields_text = ",".join(source_fields) if source_fields else "NA"

        if gene in genes_with_info:
            credible = "YES"
            reason = "Has enoughInfo=true"
            validation_status = "RELIABLE"
            credible_genes.add(gene)
        else:
            if cell_count < 10:
                credible = "NO"
                reason = "Too few cells (<10)"
                validation_status = "UNRELIABLE"
            elif cell_count < 50:
                credible = "MAYBE"
                reason = "Few cells (10-50), use with caution"
                validation_status = "MARGINAL"
                marginal_genes.add(gene)
            else:
                credible = "YES"
                reason = f"Sufficient cells ({cell_count}), likely reliable despite MixCR warning"
                validation_status = "PROBABLY_OK"
                credible_genes.add(gene)

        enough_info = "true" if gene in genes_with_info else "false"
        f.write(
            f"{gene}\t{segment_type}\t{credible}\t{reason}\t{cell_count}\t"
            f"{validation_status}\t{enough_info}\t{mixcr_status}\t{source_fields_text}\n"
        )

print(f"✓ Report: {report_file}")
print(f"  Credible genes: {len(credible_genes)}")
print(f"  Marginal genes: {len(marginal_genes)}")
print(f"  Genes with real V/D/J/C assignments: {len(real_data)}")
print(f"  TRAJ/TRBJ/TRAC/TRBC genes in real assignments: {sum(1 for g in real_data if g.startswith(('TRAJ', 'TRBJ', 'TRAC', 'TRBC')))}")

# Filter alleles_stats. Note: C alleles are usually absent from findAlleles alleles_stats.tsv,
# so C suballeles are added through the independent constant-region reports below.
credible_file = f"{analysis_dir}/alleles_stats_CREDIBLE_ONLY.tsv"
comprehensive_file = f"{analysis_dir}/alleles_stats_CREDIBLE_PLUS_MARGINAL.tsv"

with open(alleles_file, 'r') as f_in, open(credible_file, 'w') as f_out:
    reader = csv.DictReader(f_in, delimiter='\t')
    writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter='\t')
    writer.writeheader()
    for row in reader:
        if row.get('geneName', '') in credible_genes:
            writer.writerow(row)

with open(alleles_file, 'r') as f_in, open(comprehensive_file, 'w') as f_out:
    reader = csv.DictReader(f_in, delimiter='\t')
    writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter='\t')
    writer.writeheader()
    for row in reader:
        gene = row.get('geneName', '')
        if gene in credible_genes or gene in marginal_genes:
            writer.writerow(row)

print(f"✓ Filtered files created:")
print(f"  - {credible_file}")
print(f"  - {comprehensive_file}")

# Build findAlleles-only dominant summary for V/D/J genes present in credible_file.
print("STEP 12: Generating findAlleles-only dominant allele summary...")

findalleles_grouped = defaultdict(list)
with open(credible_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row.get('geneName', '')
        allele = row.get('alleleName', '')
        clones = safe_int(row.get('clonesCount', 0))
        total = safe_int(row.get('totalClonesCountForGene', 0))
        status = row.get('status', '')

        if gene and allele and total > 0:
            findalleles_grouped[gene].append({
                'allele': allele,
                'support': clones,
                'total': total,
                'status': status,
                'segment': classify_gene_segment(gene),
                'support_unit': 'clones',
                'evidence_source': 'findAlleles',
            })

findalleles_summary_file = f"{analysis_dir}/dominant_alleles_summary_findAlleles_only.tsv"
findalleles_dominant_rows = []

for gene, alleles in sorted(findalleles_grouped.items()):
    active = [a for a in alleles if a['status'] != removed_status]
    if not active:
        continue

    dominant = max(active, key=lambda x: x['support'])
    total = max(a['total'] for a in active)
    pct = round(dominant['support'] / total * 100, 1) if total > 0 else 0
    confidence = confidence_from_pct(pct)
    segment_type = classify_gene_segment(gene)

    if pct < 95:
        all_options_pct = "; ".join(
            f"{a['allele']}:{a['support']}/{total} ({round(a['support'] / total * 100, 1)}%)"
            for a in sorted(active, key=lambda x: x['support'], reverse=True)
        )
    else:
        all_options_pct = ""

    findalleles_dominant_rows.append({
        'geneName': gene,
        'segmentType': segment_type,
        'dominantAllele': dominant['allele'],
        'dominantSupportCount': dominant['support'],
        'totalSupportForGene': total,
        'supportUnit': 'clones',
        'dominantPct': f"{pct}%",
        'nAllelesTestedOrObserved': len(active),
        'confidence': confidence,
        'evidenceSource': 'findAlleles',
        'validationStatus': validation_status_for_gene(gene, credible_genes, marginal_genes),
        'allAlleleOptionsPct_ifDominantLt95': all_options_pct,
    })

with open(findalleles_summary_file, 'w', newline='') as f:
    fieldnames = [
        'geneName',
        'segmentType',
        'dominantAllele',
        'dominantSupportCount',
        'totalSupportForGene',
        'supportUnit',
        'dominantPct',
        'nAllelesTestedOrObserved',
        'confidence',
        'evidenceSource',
        'validationStatus',
        'allAlleleOptionsPct_ifDominantLt95',
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(findalleles_dominant_rows)

# Constant-region C allele/suballele summaries from exported bestCHit assignments.
print("STEP 13: Generating constant-region C allele/suballele summaries from bestCHit...")

c_alleles_summary_file = f"{analysis_dir}/constant_region_alleles_summary.tsv"
c_dominant_summary_file = f"{analysis_dir}/constant_region_dominant_alleles_summary.tsv"

c_rows = []
for (gene, allele), data in sorted(c_allele_data.items()):
    total_cells_for_gene = len(c_gene_cells[gene])
    cells_count = len(data['cells'])
    pct = round(cells_count / total_cells_for_gene * 100, 1) if total_cells_for_gene > 0 else 0

    c_rows.append({
        'geneName': gene,
        'segmentType': 'C',
        'alleleName': allele,
        'cellsCount': cells_count,
        'totalCellsForGene': total_cells_for_gene,
        'allelePct': f"{pct}%",
        'rowsCount': data['rows'],
        'readCountSum': data['read_count'],
        'uniqueMoleculeCountSum': data['molecule_count'],
        'nSamples': len(data['samples']),
        'samples': ",".join(sorted(data['samples'])),
        'sourceFields': ",".join(sorted(data['source_fields'])),
        'validationStatus': validation_status_for_gene(gene, credible_genes, marginal_genes),
    })

with open(c_alleles_summary_file, 'w', newline='') as f:
    fieldnames = [
        'geneName',
        'segmentType',
        'alleleName',
        'cellsCount',
        'totalCellsForGene',
        'allelePct',
        'rowsCount',
        'readCountSum',
        'uniqueMoleculeCountSum',
        'nSamples',
        'samples',
        'sourceFields',
        'validationStatus',
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(c_rows)

c_dominant_rows = []
for gene in sorted(c_gene_cells.keys()):
    gene_alleles = [row for row in c_rows if row['geneName'] == gene]
    if not gene_alleles:
        continue

    dominant = max(gene_alleles, key=lambda x: safe_int(x['cellsCount']))
    total = safe_int(dominant['totalCellsForGene'])
    dominant_support = safe_int(dominant['cellsCount'])
    pct = round(dominant_support / total * 100, 1) if total > 0 else 0
    confidence = confidence_from_pct(pct)

    if pct < 95:
        all_options_pct = "; ".join(
            f"{row['alleleName']}:{row['cellsCount']}/{total} ({row['allelePct']})"
            for row in sorted(gene_alleles, key=lambda x: safe_int(x['cellsCount']), reverse=True)
        )
    else:
        all_options_pct = ""

    c_dominant_rows.append({
        'geneName': gene,
        'segmentType': 'C',
        'dominantAllele': dominant['alleleName'],
        'dominantSupportCount': dominant_support,
        'totalSupportForGene': total,
        'supportUnit': 'cells',
        'dominantPct': f"{pct}%",
        'nAllelesTestedOrObserved': len(gene_alleles),
        'confidence': confidence,
        'evidenceSource': 'bestCHit',
        'validationStatus': validation_status_for_gene(gene, credible_genes, marginal_genes),
        'allAlleleOptionsPct_ifDominantLt95': all_options_pct,
    })

with open(c_dominant_summary_file, 'w', newline='') as f:
    fieldnames = [
        'geneName',
        'segmentType',
        'dominantAllele',
        'dominantSupportCount',
        'totalSupportForGene',
        'supportUnit',
        'dominantPct',
        'nAllelesTestedOrObserved',
        'confidence',
        'evidenceSource',
        'validationStatus',
        'allAlleleOptionsPct_ifDominantLt95',
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(c_dominant_rows)

# Combined final report: findAlleles V/D/J plus C suballeles from bestCHit.
print("STEP 14: Generating combined dominant allele summary including C suballeles...")

summary_file = f"{analysis_dir}/dominant_alleles_summary.tsv"
combined_rows = []
combined_rows.extend(findalleles_dominant_rows)
combined_rows.extend(row for row in c_dominant_rows if row['validationStatus'] in {'YES', 'MAYBE'})

segment_order = {'V': 0, 'D': 1, 'J': 2, 'C': 3, 'OTHER': 4}
combined_rows.sort(key=lambda x: (segment_order.get(x['segmentType'], 99), x['geneName'], x['dominantAllele']))

with open(summary_file, 'w', newline='') as f:
    fieldnames = [
        'geneName',
        'segmentType',
        'dominantAllele',
        'dominantSupportCount',
        'totalSupportForGene',
        'supportUnit',
        'dominantPct',
        'nAllelesTestedOrObserved',
        'confidence',
        'evidenceSource',
        'validationStatus',
        'allAlleleOptionsPct_ifDominantLt95',
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(combined_rows)

print(f"✓ Dominant allele summary: {summary_file}")
print(f"✓ findAlleles-only dominant summary: {findalleles_summary_file}")
print(f"✓ Constant-region allele summary: {c_alleles_summary_file}")
print(f"✓ Constant-region dominant summary: {c_dominant_summary_file}")
print(f"  C genes with observed allele/suballele assignments: {len(c_gene_cells)}")
print(f"  C allele/suballele rows: {len(c_rows)}")
print(f"  TRAC/TRBC genes in combined summary: {sum(1 for row in combined_rows if row['geneName'].startswith(('TRAC', 'TRBC')))}")

VALIDATION_SCRIPT

#summary
log_section "FINAL SUMMARY"

echo ""
echo "✓ Pipeline complete!"
echo ""
echo "Results in: ${WORK_DIR}/"
echo "  - imgt/imgt_result.alleles_stats.tsv (original MixCR findAlleles output; mainly V/J alleles)"
echo "  - imgt/imgt_result.<Sample>.cell_alleles.tsv (exported clone/cell hits including C hits)"
echo "  - libraries/${LIB_NAME}.json.gz (patient-specific IMGT library used in this run)"
echo "  - libraries/${LIB_NAME}.metadata.tsv (library build timestamp and IMGT source timestamp)"
echo ""
echo "Validation in: ${WORK_DIR}/validation/"
echo "  - ${LIB_NAME}.status.tsv (last IMGT check date, remote Last-Modified, update action)"
echo "  - credibility_report.tsv (segment-aware V/D/J/C credibility breakdown)"
echo "  - alleles_stats_CREDIBLE_ONLY.tsv (credible genes present in alleles_stats.tsv)"
echo "  - alleles_stats_CREDIBLE_PLUS_MARGINAL.tsv (credible + marginal genes present in alleles_stats.tsv)"
echo "  - dominant_alleles_summary_findAlleles_only.tsv (V/D/J dominant alleles from findAlleles)"
echo "  - constant_region_alleles_summary.tsv (TRAC/TRBC allele/suballele support from bestCHit)"
echo "  - constant_region_dominant_alleles_summary.tsv (dominant TRAC/TRBC suballeles from bestCHit)"
echo "  - dominant_alleles_summary.tsv (combined V/D/J + C dominant allele summary)"
echo ""