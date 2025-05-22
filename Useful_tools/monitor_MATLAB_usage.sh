#!/bin/bash

# Check for PID
if [ -z "$1" ]; then
  echo "Usage: $0 <PID>"
  exit 1
fi

PID=$1
LOGFILE="monitor_log_PID${PID}.txt"

# Validate PID
if ! ps -p $PID > /dev/null; then
  echo "Process with PID $PID not found."
  exit 1
fi

# System info
TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
TOTAL_MEM_GB=$(awk -v mem_kb=$TOTAL_MEM_KB 'BEGIN {printf "%.2f", mem_kb / 1024 / 1024}')
TOTAL_CORES=$(nproc)

echo "Monitoring PID: $PID"
echo "System RAM: $TOTAL_MEM_GB GB | Logical CPUs: $TOTAL_CORES"
echo "Logging to: $LOGFILE"
echo

# Header
HEADER="Timestamp | CPU % of Total | CPU Used (cores) | %MEM | RAM Used (GB) | GPU Util % | GPU MEM Used (MB)"
echo "$HEADER"
echo "$HEADER" > "$LOGFILE"

# Loop
while ps -p $PID > /dev/null; do
  TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")

  read CPU MEM <<< $(ps -p $PID -o %cpu,%mem --no-headers)
  RSS_KB=$(ps -p $PID -o rss=)
  RAM_GB=$(awk -v rss_kb=$RSS_KB 'BEGIN {printf "%.2f", rss_kb / 1024 / 1024}')
  USED_CORES=$(awk -v cpu=$CPU 'BEGIN {printf "%.2f", cpu / 100}')
  PERCENT_OF_TOTAL_CPU=$(awk -v cpu=$CPU -v total=$TOTAL_CORES 'BEGIN {printf "%.2f", cpu / (total * 100) * 100}')

  # GPU info
  GPU_MEM_USED="N/A"
  GPU_UTIL="N/A"
  if command -v nvidia-smi &> /dev/null; then
    GPU_MEM_USED=$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits | grep $PID | awk '{print $2}' | head -n1)
    GPU_UTIL=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits | head -n1)
  fi

  LINE="$TIMESTAMP | $PERCENT_OF_TOTAL_CPU% | $USED_CORES | $MEM | $RAM_GB | ${GPU_UTIL:-N/A}% | ${GPU_MEM_USED:-N/A}"
  echo "$LINE"
  echo "$LINE" >> "$LOGFILE"

  sleep 1
done

echo "Process $PID has exited." | tee -a "$LOGFILE"

