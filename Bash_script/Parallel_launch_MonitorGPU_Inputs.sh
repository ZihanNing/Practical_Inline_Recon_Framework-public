#!/usr/bin/env bash
# This is a bash script used to call the custom recon launched in the background by a MATLAB handle
# It can optionally wait for specific input files to appear before launching.
#
# Input:
#   -l LOG_FILE             Path to the GPU status log (to log GPU usage) [required]
#   -f MATLAB_FUNCTION      Name of the MATLAB function to run [required]
#   -r RECON_ID             Recon identifier to pass to the MATLAB function [required]
#   -w MAX_WAIT_TIME        (in sec) after the maximum queue time, the GPU status will be double-checked [default: 900]
#   -i INTERVAL             (in sec) during the queue-up, interval to check the log for updates [default: 10]
#   -c FLAG_CHECK_INPUTS    0 = do not check for input files (default)
#                           1 = check that all input files exist before launching
#   -u INPUT_LIST           (only if -c 1) A quoted, space-separated list of full paths or patterns to check
#
# Behavior:
#   • If FLAG_CHECK_INPUTS=1, the script will wait until every entry in INPUT_LIST exists
#     before proceeding to check for a free GPU. If any input is missing, it sleeps for INTERVAL
#     seconds and retries, accumulating WAIT_TIME.
#   • Once inputs are ready (or if FLAG_CHECK_INPUTS=0), it checks the GPU log every INTERVAL
#     seconds. When it finds a free GPU, it launches MATLAB in the background with the specified function.
#   • If WAIT_TIME ≥ MAX_WAIT_TIME at any point, it will run verify_GPU_status.sh in the background
#     and reset WAIT_TIME to zero.
#
# by Zihan Ning
# @ King's College London
# Last modified: 3 June 2025

usage() {
  cat <<EOF
Usage: $0 \
  -l LOG_FILE \
  -f MATLAB_FUNCTION \
  -r RECON_ID \
  [-w MAX_WAIT_TIME] \
  [-i INTERVAL] \
  [-c FLAG_CHECK_INPUTS] \
  [-u INPUT_LIST]

  -l LOG_FILE             Path to the GPU status log (required)
  -f MATLAB_FUNCTION      Name of the MATLAB function to run (required)
  -r RECON_ID             Recon identifier to pass to the MATLAB function (required)
  -w MAX_WAIT_TIME        Max wait time in seconds before running verify_GPU_status.sh (default: 900)
  -i INTERVAL             Polling interval in seconds (default: 10)
  -c FLAG_CHECK_INPUTS    0 = do not check for input files (default)
                          1 = check that all input files exist before launching
  -u INPUT_LIST           Quoted, space-separated list of full paths or filename patterns to check.
                          Only used if -c 1. For example:
                            -u "/path/to/file1 /path/to/file2 /path/to/ExterREF_CMS_*"
EOF
  exit 1
}

# --- Defaults ---
MAX_WAIT_TIME=900
INTERVAL=10
FLAG_CHECK_INPUTS=0
INPUT_LIST=""

# --- Parse options ---
while getopts ":l:f:r:w:i:c:u:" opt; do
  case $opt in
    l) LOG_FILE="$OPTARG" ;;
    f) MATLAB_FUNCTION="$OPTARG" ;;
    r) RECON_ID="$OPTARG" ;;
    w) MAX_WAIT_TIME="$OPTARG" ;;
    i) INTERVAL="$OPTARG" ;;
    c) FLAG_CHECK_INPUTS="$OPTARG" ;;
    u) INPUT_LIST="$OPTARG" ;;
    *) usage ;;
  esac
done

# --- Validate required args ---
if [[ -z "$LOG_FILE" || -z "$MATLAB_FUNCTION" || -z "$RECON_ID" ]]; then
  usage
fi

# If input checking is enabled, ensure INPUT_LIST is provided
if [[ "$FLAG_CHECK_INPUTS" -eq 1 && -z "$INPUT_LIST" ]]; then
  echo "Error: -u INPUT_LIST must be provided when -c 1" >&2
  usage
fi

WAIT_TIME=0

get_latest_gpu_status() {
  grep -v '^$' "$LOG_FILE" | tail -n1
}

check_gpu_free() {
  local line="${1#*: }" parts status gpu_id state
  IFS=';' read -ra parts <<<"$line"
  for status in "${parts[@]}"; do
    status="$(echo "$status" | xargs)"        # trim whitespace
    gpu_id="${status%%_*}"
    state="${status#*_}"
    state="${state%%_*}"
    if [[ "$state" == "free" ]]; then
      echo "$gpu_id"
      return 0
    fi
  done
  return 1
}

check_inputs_exist() {
  local missing=0
  # Split INPUT_LIST into an array by whitespace
  read -ra files <<<"$INPUT_LIST"
  for item in "${files[@]}"; do
    # Use globbing: if the pattern matches at least one file, consider it present
    shopt -s nullglob
    matches=( $item )
    shopt -u nullglob
    if [[ ${#matches[@]} -eq 0 ]]; then
      echo "$item"
      missing=1
    fi
  done
  return $missing
}

# --- Main loop ---
while true; do
  if [[ ! -f "$LOG_FILE" ]]; then
    echo "Error: log file not found: $LOG_FILE" >&2
    exit 1
  fi

  # If input-checking is enabled, verify all inputs are present
  if [[ "$FLAG_CHECK_INPUTS" -eq 1 ]]; then
    if missing_file=$(check_inputs_exist); then
      # Some files/patterns are missing
      echo "Waiting for input(s): $missing_file not found. Waited $WAIT_TIME/$MAX_WAIT_TIME sec. Retrying in $INTERVAL sec..."
      sleep "$INTERVAL"
      WAIT_TIME=$(( WAIT_TIME + INTERVAL ))
      if (( WAIT_TIME >= MAX_WAIT_TIME )); then
        echo "Max wait exceeded while waiting for inputs. Running verify_GPU_status.sh in background..."
        bash verify_GPU_status.sh &
        WAIT_TIME=0
      fi
      continue
    fi
    # All inputs exist; proceed to GPU-check
  fi

  # Check GPU status
  latest=$(get_latest_gpu_status)
  if [[ -z "$latest" ]]; then
    echo "Error: no valid lines in $LOG_FILE" >&2
    exit 1
  fi

  if free_gpu=$(check_gpu_free "$latest"); then
    echo "GPU $free_gpu is free. Launching MATLAB..."
    ulimit -c unlimited
    ulimit -n 4096
    # Call the MATLAB function with the Recon_ID argument
    matlab -r "try, ${MATLAB_FUNCTION}('${RECON_ID}'), catch e, disp(getReport(e)), exit(1), end, exit(0);" &
    exit 0
  fi

  echo "No free GPU; waited $WAIT_TIME/$MAX_WAIT_TIME sec. Retrying in $INTERVAL sec..."
  sleep "$INTERVAL"
  WAIT_TIME=$(( WAIT_TIME + INTERVAL ))

  if (( WAIT_TIME >= MAX_WAIT_TIME )); then
    echo "Max wait exceeded. Running verify_GPU_status.sh in background..."
    bash verify_GPU_status.sh &
    WAIT_TIME=0
  fi
done

