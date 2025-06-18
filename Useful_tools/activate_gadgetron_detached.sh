#!/bin/bash

# Check if any program is running on port 9002 (the default port of gadgetron)
if lsof -Pi :9002 -sTCP:LISTEN -t >/dev/null; then
    echo "A program is running on port 9002. Killing the program..."
    kill $(lsof -t -i:9002)
    echo "A previous Gadgetron client or other program on port 9002 has been killed."
else
    echo "No program is running on port 9002."
fi

# Source the Conda initialization script
source /home/gadgetron/anaconda3/etc/profile.d/conda.sh

# Activate the conda environment
conda activate gadgetron

# Run Gadgetron in the background with nohup, and log to a timestamped file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$HOME/detached_gadgetron_log_${TIMESTAMP}.txt"

echo "Starting Gadgetron client in detached mode. Logs will be saved to $LOG_FILE"

nohup gadgetron > "$LOG_FILE" 2>&1 &


