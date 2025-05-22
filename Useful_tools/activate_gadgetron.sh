#!/bin/bash
# A short bash script to activate gadgetron client on the default port (9002) quickly
# 
# Zihan Ning <zihan.1.ning@kcl.ac.uk>
# King's college london, 22-May-2025


# Check if any program is running on port 9002 (the defualt port of gadgetron)
if lsof -Pi :9002 -sTCP:LISTEN -t >/dev/null; then
    # If a program is found, kill it
    echo "A program is running on port 9002. Killing the program..."
    kill $(lsof -t -i:9002)
    echo "A previous Gadgetron client/or other programme on port 9002 has been killed."
else
    echo "No program is running on port 9002."
fi

# Source the Conda initialization script
# Modify the pathway that suits your installation
source /home/zn23/anaconda3/etc/profile.d/conda.sh

# Activate the conda environment
conda activate gadgetron

# Execute the gadgetron command
gadgetron









