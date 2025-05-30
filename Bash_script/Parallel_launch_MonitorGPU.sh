#!/usr/bin/env bash
# This is a bash script used to call the custom recon launched in the background by a matlab handle
# Input:
# 	LOG_FILE: the path of the log file (to log GPU usage)
# 	MATLAB_FUNCTION: the matlab handle to be called with the custom recon integrated
# 	MAX_WAIT_TIME: (in sec) after the maximum queuing time, the GPU status will be double check
# 	INTERVAL: (in sec) during the queue-up, the interval to check the log for update
#
# by Zihan Ning
# @ King's college london
# Last modifed: 29 May 2025

usage() {
  cat <<EOF
Usage: $0 -l LOG_FILE -f MATLAB_FUNCTION -r Recon_ID [-w MAX_WAIT_TIME] [-i INTERVAL]

  -l LOG_FILE          Path to the GPU status log (required)
  -f MATLAB_FUNCTION   Name of the MATLAB function to run (required)
  -r Recon_ID          Recon identifier to pass to the MATLAB function (required)
  -w MAX_WAIT_TIME     Max wait time in seconds before running verify_GPU_status.sh (default: 900)
  -i INTERVAL          Polling interval in seconds (default: 10)
EOF
  exit 1
}

# --- Defaults ---
MAX_WAIT_TIME=900
INTERVAL=10

# --- Parse options ---
while getopts ":l:f:r:w:i:" opt; do
  case $opt in
    l) LOG_FILE="$OPTARG" ;;
    f) MATLAB_FUNCTION="$OPTARG" ;;
    r) RECON_ID="$OPTARG" ;;
    w) MAX_WAIT_TIME="$OPTARG" ;;
    i) INTERVAL="$OPTARG" ;;
    *) usage ;;
  esac
done

# --- Validate required args ---
[[ -z "$LOG_FILE" || -z "$MATLAB_FUNCTION" || -z "$RECON_ID" ]] && usage

WAIT_TIME=0

get_latest_gpu_status() {
  grep -v '^$' "$LOG_FILE" | tail -n1
}

check_gpu_free() {
  local line="${1#*: }" parts status gpu_id state
  IFS=';' read -ra parts <<<"$line"
  for status in "${parts[@]}"; do
    status="$(echo "$status" | xargs)"        # trim
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

# --- Main loop ---
while true; do
  if [[ ! -f "$LOG_FILE" ]]; then
    echo "Error: log file not found: $LOG_FILE" >&2
    exit 1
  fi

  latest=$(get_latest_gpu_status)
  if [[ -z "$latest" ]]; then
    echo "Error: no valid lines in $LOG_FILE" >&2
    exit 1
  fi

  if free_gpu=$(check_gpu_free "$latest"); then
    echo "GPU $free_gpu is free. Launching MATLAB..."
    ulimit -c unlimited
    ulimit -n 4096
    # Call your MATLAB function with the Recon_ID argument
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

