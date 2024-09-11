#!/bin/bash

module load cellgen/conda

# Define the environment name and path
ENV_NAME="merge-bam"
ENV_PATH="/software/cellgen/team274/miniconda3/envs/$ENV_NAME"
REPO_URL="https://github.com/BehjatiLab/merge-bam"

# Create the new conda environment with Python
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
echo "Installing pip in the environment."
conda install pip -y
if [ $? -ne 0 ]; then
  echo "Failed to install pip. Exiting."
  exit 1
fi

# Clone the GitHub repository
echo "Cloning the repository from $REPO_URL."
git clone "$REPO_URL"
if [ $? -ne 0 ]; then
  echo "Failed to clone the repository. Exiting."
  exit 1
fi

REPO_DIR=$(basename "$REPO_URL" .git)

# Navigate to the cloned directory
cd "$REPO_DIR" || { echo "Failed to navigate to the repository directory. Exiting."; exit 1; }

# Install the package from the cloned repository
echo "Installing the package from the cloned repository."
pip install . --verbose
if [ $? -ne 0 ]; then
  echo "Failed to install the package. Exiting."
  exit 1
fi

# Set executable permissions for the installed scripts
BIN_DIR=$(python -c 'import sysconfig; print(sysconfig.get_path("scripts"))')
chmod +x "$BIN_DIR/create-test-bed" "$BIN_DIR/merge-bam" "$BIN_DIR/bed-from-txt"

# Verify installation by checking if the entry points are installed
echo "Checking if the commands are installed correctly."
if ! command -v create-test-bed &> /dev/null || ! command -v merge-bam &> /dev/null || ! command -v bed-from-txt &> /dev/null; then
  echo "Error: Command line tools were not installed correctly."
  exit 1
fi

# Navigate back to the parent directory and remove the cloned repository
cd ..
rm -rf "$REPO_DIR"

echo "Environment '$ENV_NAME' created, and package installed from $REPO_URL."

# List installed packages and check that 'my_bam_tools' is among them
echo "Checking installed packages:"
pip list | grep 'my_bam_tools'

# Check if modules can be imported
echo "Verifying module imports."
python -c "import bed_from_mut_txt; import create_position_bed; import merge_bam; print('Modules imported successfully.')"

echo "Installation and setup completed successfully."
