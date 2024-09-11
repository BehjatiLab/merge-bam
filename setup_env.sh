#!/bin/bash

module load cellgen/conda

# Define the environment name and path
ENV_NAME="merge-bam"
ENV_PATH="/software/cellgen/team274/miniconda3/envs/$ENV_NAME"
REPO_URL="https://github.com/BehjatiLab/merge-bam"

# Create the new conda environment with Python (specify the version if needed)
echo "Creating environment '$ENV_NAME' at $ENV_PATH."
conda create -p "$ENV_PATH" python=3.11.9 -y

# Activate the newly created environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_PATH"

# Install pip in the conda environment
conda install pip -y

# Clone the GitHub repository
git clone "$REPO_URL"
REPO_DIR=$(basename "$REPO_URL" .git)

# Navigate to the cloned directory
cd "$REPO_DIR" || exit

# Install the package from the cloned repository
pip install .

# Navigate back to the parent directory and remove the cloned repository
cd ..
rm -rf "$REPO_DIR"

echo "Environment '$ENV_NAME' created, and package installed from $REPO_URL."
