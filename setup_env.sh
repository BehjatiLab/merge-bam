#!/bin/bash

module load cellgen/conda

# Define the environment name and path
ENV_NAME="merge-bam"
ENV_PATH="/software/cellgen/team274/miniconda3/envs/$ENV_NAME"
REPO_URL="https://github.com/BehjatiLab/merge-bam"

# Create the new conda environment with Python (specify the version if needed)
echo "Creating environment '$ENV_NAME' at $ENV_PATH."
conda create -p "$ENV_PATH" python=3.11.9 -y

if [ $? -ne 0 ]; then
  echo "Failed to create the environment. Exiting."
  exit 1
fi

# Activate the newly created environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_PATH"
if [ $? -ne 0 ]; then
  echo "Failed to activate the environment. Exiting."
  exit 1
fi

# Install pip in the conda environment
conda install pip -y

if [ $? -ne 0 ]; then
  echo "Failed to install pip. Exiting."
  exit 1
fi

# Clone the GitHub repository
git clone "$REPO_URL"

if [ $? -ne 0 ]; then
  echo "Failed to clone the repository. Exiting."
  exit 1
fi

REPO_DIR=$(basename "$REPO_URL" .git)

# Navigate to the cloned directory
cd "$REPO_DIR" || exit

# Install the package from the cloned repository
pip install .

if [ $? -ne 0 ]; then
  echo "Failed to install the package. Exiting."
  exit 1
fi

# Navigate back to the parent directory and remove the cloned repository
cd ..
rm -rf "$REPO_DIR"

echo "Environment '$ENV_NAME' created, and package installed from $REPO_URL."
