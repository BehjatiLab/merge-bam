#!/bin/bash

module load cellgen/conda

# Define the environment name and path
ENV_NAME="merge-bam"
ENV_PATH="/software/cellgen/team274/miniconda3/envs/$ENV_NAME"

# Check if the environment exists
if [ -d "$ENV_PATH" ]; then
    echo "Environment '$ENV_NAME' exists at $ENV_PATH."
else
    echo "Environment '$ENV_NAME' does not exist. Creating it at $ENV_PATH."

    # Create the new conda environment with Python (specify the version if needed)
    conda create -p "$ENV_PATH" python=3.11.9 -y

    # Activate the newly created environment
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate "$ENV_PATH"

    # Install pip in the conda environment
    conda install pip -y

    # Install the required packages using pip from the requirements.txt file
    pip install -r requirements.txt

    echo "Environment '$ENV_NAME' created and packages installed."
fi
