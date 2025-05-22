#!/bin/bash

# Check if nvidia-smi is installed
if ! command -v nvidia-smi &> /dev/null
then
    echo "nvidia-smi could not be found. Make sure NVIDIA drivers are installed."
    exit 1
fi

# Check if mpstat (from sysstat package) is installed
if ! command -v mpstat &> /dev/null
then
    echo "mpstat could not be found. Please install the sysstat package."
    exit 1
fi

# Function to get CPU usage
get_cpu_usage() {
    echo "------ CPU Usage ------"
    # Get the idle percentage from mpstat
    local idle=$(mpstat 1 1 | awk '/Average:/ {print $12}')
    # Calculate CPU usage as 100 - idle
    local usage=$(echo "100 - $idle" | bc)
    local available=$idle
    printf "Used: %s%% | Available: %s%%\n\n" "$usage" "$available"
}

# Function to normalize memory sizes from free -h
normalize_size() {
    local size=$1
    local value=$(echo "$size" | grep -o -E '[0-9.]+')       # Extract the numeric part
    local unit=$(echo "$size" | grep -o -E '[a-zA-Z]+')    # Extract the unit

    # Convert based on the unit
    case "$unit" in
        G|Gi)
            echo "$value"  # GiB, no conversion needed
            ;;
        M|Mi)
            # Convert MiB to GiB with two decimal places
            echo "$(echo "scale=2; $value / 1024" | bc)"
            ;;
        K|Ki)
            # Convert KiB to GiB with two decimal places
            echo "$(echo "scale=2; $value / 1048576" | bc)"
            ;;
        *)
            # If no unit or unrecognized unit, assume the value is in GiB
            echo "$value"
            ;;
    esac
}

# Function to get RAM and swap usage using free -h
get_ram_usage() {
    echo "------ RAM Usage ------"
    
    # Get memory details from free -h command (total, used, free for Mem and Swap)
    total_mem=$(free -h | grep Mem | awk '{print $2}')
    used_mem=$(free -h | grep Mem | awk '{print $3}')
    free_mem=$(free -h | grep Mem | awk '{print $4}')
    
    total_swap=$(free -h | grep Swap | awk '{print $2}')
    used_swap=$(free -h | grep Swap | awk '{print $3}')
    free_swap=$(free -h | grep Swap | awk '{print $4}')
    
    # Normalize the units (GiB)
    total_mem_gb=$(normalize_size "$total_mem")
    used_mem_gb=$(normalize_size "$used_mem")
    free_mem_gb=$(normalize_size "$free_mem")

    total_swap_gb=$(normalize_size "$total_swap")
    used_swap_gb=$(normalize_size "$used_swap")
    free_swap_gb=$(normalize_size "$free_swap")

    # Calculate free percentages
    # Prevent division by zero
    if (( $(echo "$total_mem_gb > 0" | bc -l) )); then
        free_mem_percent=$(echo "scale=2; ($free_mem_gb / $total_mem_gb) * 100" | bc)
    else
        free_mem_percent="N/A"
    fi

    if (( $(echo "$total_swap_gb > 0" | bc -l) )); then
        free_swap_percent=$(echo "scale=2; ($free_swap_gb / $total_swap_gb) * 100" | bc)
    else
        free_swap_percent="N/A"
    fi
    
    # Display memory and swap usage with free percentages
    printf "Mem: Total: %sGi | Used: %sGi | Free percentage: %s%%\n" "$total_mem_gb" "$used_mem_gb" "$free_mem_percent"
    printf "Swap: Total: %sGi | Used: %sGi | Free percentage: %s%%\n\n" "$total_swap_gb" "$used_swap_gb" "$free_swap_percent"
}

# Function to get GPU usage using nvidia-smi
get_gpu_usage() {
    echo "------ GPU Usage ------"
    # Print table header
    printf "%-5s %-25s %-18s %-17s %-20s\n" "Index" "GPU Type" "Total Memory (MB)" "Used Memory (MB)" "Available Memory (MB)"
    echo "-----------------------------------------------------------------------------------------------"
    # Fetch and format GPU data
    nvidia-smi --query-gpu=index,name,memory.total,memory.used,memory.free --format=csv,noheader,nounits | while IFS=',' read -r index name total used free
    do
        printf "%-5s %-25s %-18s %-17s %-20s\n" "$index" "$name" "$total" "$used" "$free"
    done
    echo ""
}

# Main loop to monitor CPU, RAM, and GPU usage every 5 seconds
while true; do
    clear  # Clear the terminal for clean output
    get_cpu_usage  # Get and display CPU usage
    get_ram_usage  # Get and display RAM usage (Mem and Swap with free percentages)
    get_gpu_usage  # Get and display GPU usage
    sleep 5        # Update every 5 seconds
done

