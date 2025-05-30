#!/bin/bash

# Path to the GPU status log file
LOG_FILE="$HOME/DISORDER_gadgetron/log/GPU_status_log.txt"
CHECK_INTERVAL=1  # in second
MAX_CHECKS=5      # Number of consecutive checks required
MIN_UNUSED_MEM=95  # Minimum percentage of unused memory to consider the GPU as free
MAX_LAST_UPDATE_TIME=20  # maximum time of the last updated timestamp to be verified

# Function to read the latest non-empty line from the log file
get_latest_gpu_status() {
    # Read the last non-empty line from the file
    latest_line=$(grep -v '^$' "$LOG_FILE" | tail -n 1)
    echo "$latest_line"
}

# Function to parse GPU status from the latest line and return an associative array
parse_gpu_status() {
    local status_line="$1"
    declare -A gpu_status_map

    # Remove everything before the first ": "
    status_line=$(echo "$status_line" | sed 's/.*\(GPU1\)/\1/')

    # Split status by GPU
    IFS=';' read -r -a gpu_info <<< "$status_line"

    for gpu in "${gpu_info[@]}"; do
        # Trim leading/trailing whitespace
        gpu=$(echo "$gpu" | xargs)

        # Parse GPU status and timestamp, keeping original format
        gpu_id=$(echo "$gpu" | awk -F'_' '{print $1}')
        gpu_status=$(echo "$gpu" | awk -F'_' '{print $2}')
        gpu_time=$(echo "$gpu" | awk -F'_' '{print $3 "_" $4}')

        # Store GPU status and timestamp in the associative array (preserve the original format)
        gpu_status_map["$gpu_id"]="${gpu_status}_${gpu_time}"
    done

    # Return the associative array
    echo "$(declare -p gpu_status_map)"
}

# Function to convert the GPU timestamp to a format bash can recognize
convert_gpu_timestamp() {
    local gpu_time="$1"
    # Replace the underscore between date and time with a space
    echo "$gpu_time" | sed 's/_/ /'
}

# Function to calculate the time difference between the GPU timestamp and the current time
calculate_time_difference() {
    local gpu_time="$1"

    # Convert the timestamp to a format recognized by bash
    gpu_time=$(convert_gpu_timestamp "$gpu_time")

    # Get the current time in the same format
    current_time=$(date "+%Y-%m-%d %H:%M:%S")

    # Convert both timestamps to epoch time (seconds since 1970-01-01)
    gpu_epoch=$(date -d "$gpu_time" +%s)
    current_epoch=$(date -d "$current_time" +%s)

    # Calculate the difference in minutes
    time_diff=$(( (current_epoch - gpu_epoch) / 60 ))
    echo "$time_diff"
}

# Function to check current GPU memory usage using nvidia-smi
check_gpu_memory_usage() {
    local gpu_id="$1"
    local gpu_index="$2"  # GPU index as required by nvidia-smi
    local all_above_threshold=true
    local memory_threshold=$MIN_UNUSED_MEM

    echo "Checking memory usage for $gpu_id (index $gpu_index)..."

    for ((i=1; i<=MAX_CHECKS; i++)); do
        # Run nvidia-smi and get the memory usage (total and used) for the specific GPU
        gpu_memory_info=$(nvidia-smi --id="$gpu_index" --query-gpu=memory.total,memory.used --format=csv,noheader,nounits)
        
        # Extract total and used memory
        total_memory=$(echo "$gpu_memory_info" | awk -F',' '{print $1}')
        used_memory=$(echo "$gpu_memory_info" | awk -F',' '{print $2}')
        
        # Calculate the unused memory percentage
        unused_memory=$((total_memory - used_memory))
        unused_percentage=$((100 * unused_memory / total_memory))
        
        # Print memory usage details
        echo "Check $i: GPU $gpu_id -> Total: ${total_memory}MiB, Used: ${used_memory}MiB, Unused: ${unused_memory}MiB (${unused_percentage}%)"

        # If any check has a lower unused memory percentage than the threshold, set flag to false
        if (( unused_percentage < memory_threshold )); then
            all_above_threshold=false
        fi

        # Sleep for the specified interval before the next check
        sleep $CHECK_INTERVAL
    done

    # Return the final result
    if [ "$all_above_threshold" = true ]; then
        echo "true"
    else
        echo "false"
    fi
}


# Function to update the GPU status log file if a GPU becomes free
update_gpu_status_log() {
    local gpu_id="$1"

    # Read the latest non-empty log line again to regenerate gpu_status_map
    latest_status=$(get_latest_gpu_status)

    # Regenerate the gpu_status_map by calling parse_gpu_status
    eval "$(parse_gpu_status "$latest_status")"

    local current_time=$(date "+%Y-%m-%d_%H:%M:%S")
    local new_log_entry="${current_time}: "

    # Loop through all GPUs and update only the specified GPU, keeping others unchanged
    for id in GPU1 GPU2 GPU3; do
        if [[ "$id" == "$gpu_id" ]]; then
            # For the GPU being updated, change status to 'free' and update the timestamp
            new_log_entry+="${id}_free_${current_time}; "
        else
            # For the other GPUs, retrieve their current status and timestamp from gpu_status_map
            gpu_info="${gpu_status_map[$id]}"
            if [ -n "$gpu_info" ]; then
                new_log_entry+="${id}_${gpu_info}; "
            else
                # If the GPU info isn't found, log an error message and skip this GPU
                echo "Error: Unable to retrieve the status for $id"
                new_log_entry+="${id}_unknown_0000-00-00_00:00:00; "
            fi
        fi
    done

    # Remove the trailing semicolon and whitespace
    new_log_entry=$(echo "$new_log_entry" | sed 's/; $//')

    # Append the new log entry to the log file
    echo "$new_log_entry" >> "$LOG_FILE"
    echo "New log entry added: $new_log_entry"
}

# Main function to monitor GPU status and update the log file if necessary
monitor_gpus() {
    latest_status=$(get_latest_gpu_status)
    echo "Latest status: $latest_status"

    # Call parse_gpu_status and capture the associative array it returns
    eval "$(parse_gpu_status "$latest_status")"  # Capture the associative array generated by parse_gpu_status

    # Print the gpu_status_map before updating the log
    echo "gpu_status_map just after parsing the log:"
    declare -p gpu_status_map  # Print the entire gpu_status_map

    # Now we can access the gpu_status_map associative array
    echo "Checking GPU status, timestamp, and time intervals:"
    for gpu_id in GPU1 GPU2 GPU3; do
        # Extract status and timestamp from gpu_status_map
        gpu_info="${gpu_status_map[$gpu_id]}"
        if [ -n "$gpu_info" ]; then
            gpu_status=$(echo "$gpu_info" | awk -F'_' '{print $1}')
            gpu_time=$(echo "$gpu_info" | awk -F'_' '{print $2 "_" $3}')
            
            # Convert GPU timestamp to a recognizable format
            converted_time=$(convert_gpu_timestamp "$gpu_time")
            echo "$gpu_id status: $gpu_status, time: $converted_time"

            # Calculate the time difference from the current time
            time_difference=$(calculate_time_difference "$gpu_time")
            echo "Time interval between now and $gpu_id's last status: $time_difference minutes"

            # Check if the time difference is more than 30 minutes
            if [ "$gpu_status" == "occupied" ] && (( time_difference > $MAX_LAST_UPDATE_TIME )); then
                echo "$gpu_id's status is older than 30 minutes!"
                
                # Determine the GPU index (assuming GPU1 -> index 0, GPU2 -> index 1, etc.)
                gpu_index=$(( ${gpu_id:3:1} - 1 ))
                
                # Check the current memory usage of this GPU 3 times
                all_checks_passed=$(check_gpu_memory_usage "$gpu_id" "$gpu_index" | tee /dev/tty | tail -n 1)

                # If all checks passed the threshold, update the log file
                if [ "$all_checks_passed" = "true" ]; then
                    echo "$gpu_id meets the condition to be marked as free!"
                    
                    # Update the log for this specific GPU
                    update_gpu_status_log "$gpu_id"
                fi
            fi
        else
            echo "$gpu_id not found in log"
        fi
    done
}


# Start the GPU monitoring
monitor_gpus
