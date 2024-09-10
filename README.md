# BAM Filtering and Merging Module

This module provides scripts for filtering and merging BAM files based on intervals specified in BED files. The module includes three scripts:

- `setup_env.sh`: Sets up a Conda environment with the necessary dependencies.
- `merge_bam.py`: Filters BAM files based on intervals from a BED file and merges the results into a single output BAM file.
- `create_test_bed.py`: Generates a test BED file with specified genomic positions and optionally compresses and indexes the file using bgzip and Tabix.

## Scripts Overview

### 1. `setup_env.sh`

This script creates a Conda environment named `merge-bam` and installs the necessary Python packages from a `requirements.txt` file using pip.

**Functionality:**
- Checks if the environment already exists in the specified directory `/software/cellgen/team274/miniconda3/envs/`.
- If the environment exists, it prints a message indicating its existence.
- If the environment does not exist, it creates the environment, installs pip, and then installs the required packages.

**Usage:**
Run the script from the command line:
```bash
bash setup_env.sh
```

### 2. `merge_bam.py`

This script filters reads from one or more BAM files based on intervals from a BED file and merges the filtered reads into a single output BAM file. It supports multithreading to speed up processing and can resume from existing outputs if interrupted.

**Functionality:**
- Loads intervals from a BED file, which can be plain or compressed with bgzip and indexed with Tabix.
- Filters BAM files by the loaded intervals and writes the filtered reads to temporary BAM files.
- Merges all temporary BAM files into a final output BAM file.
- Keeps track of processed intervals to avoid redundant processing.

**Usage:**
```
python merge_bam.py -b /path/to/bam1.bam /path/to/bam2.bam \
                    -e /path/to/intervals.bed.gz \
                    -o /path/to/output.bam \
                    -t 4 \
                    --temp-dir /path/to/temp \
                    --overwrite
```

**Arguments:**
- `-b, --bam-files`: Paths to input BAM files to be filtered and merged.
- `-l`, `--bam-list`: Path to a text file containing a list of input BAM files to be filtered and merged.
- `-e, --bed-file`: Path to the BED file containing genomic intervals for filtering.
- `-o, --output-file`: Path to the output BAM file.
- `-t, --threads`: Number of threads to use for processing (default: 4).
- `--temp-dir`: Directory to store temporary files.
- `--overwrite`: Overwrite existing temporary files and output file if they exist.

### 3. `create_test_bed.py`

This script generates a BED file with a specified number of positions. It can optionally compress the file using bgzip and create a Tabix index.

**Functionality:**
- Generates random genomic positions within a specified range on a specified chromosome.
- Writes these positions to a BED file.
- Optionally compresses the BED file with bgzip and indexes it with Tabix.

**Usage:**
```
python create_test_bed.py -n 1000 \
                          -o /path/to/output.bed \
                          -c chr1 \
                          --bgzip
```

**Arguments:**
- `-n, --num-positions`: Number of genomic positions to generate.
- `-o, --output-file`: Path to the output BED file.
- `-c, --chromosome`: Chromosome name for the generated positions.
- `--bgzip`: Compress the output BED file with bgzip and create a Tabix index.

## Example workflow
1. Set up the Conda environment:
```bash
bash setup_env.sh
```

2. Generate a test BED file with random genomic positions:
```bash
python create_test_bed.py -n 1000 -o test.bed -c chr1 --bgzip
```

3. Filter and merge BAM files based on the generated BED file:
```bash
python merge_bam.py -b /path/to/bam1.bam /path/to/bam2.bam -e test.bed.gz -o merged.bam -t 4 --temp-dir temp --overwrite
```

This workflow demonstrates how to set up the environment, generate a test BED file, and filter BAM files based on the generated BED intervals, merging them into a single output file.

## Dependencies
- `Python==3.11.9`
- `pysam==0.22.1`
