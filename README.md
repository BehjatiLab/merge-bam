# BAM Filtering and Merging Command-Line Tools

This module provides command-line tools for filtering and merging BAM files based on intervals specified in BED files. It includes the following scripts:

- `setup_env.sh`: Sets up a Conda environment with the necessary dependencies.
- `merge-bam`: Filters BAM files based on intervals from a BED file and merges the results into a single output BAM file.
- `create-test-bed`: Generates a test BED file with specified genomic positions and optionally compresses and indexes the file using bgzip and Tabix.
- `bed-from-txt`: Converts a text file with mutation data into a BGZIPPED and TABIX indexed BED file.

----

## Command-Line Tools Overview

### 1. Setting Up the Environment with `setup_env.sh`

This script creates a Conda environment named `merge-bam` and installs the necessary Python packages, including the BAM tools, from a GitHub repository. This script is inteded for group members using the Sanger compute farm. Once this is run once, the package should exist for everyone in the group's conda environment directory.

**Usage:**
Run the script from the command line:
```bash
bash setup_env.sh
```

### 2. Filtering and Merging BAM Files with `merge-bam`

This command filters reads from one or more BAM files based on intervals from a BED file and merges the filtered reads into a single output BAM file. It supports multithreading to speed up processing and can resume from existing outputs if interrupted.

**Functionality:**

- Loads intervals from a BED file, which can be plain or compressed with bgzip and indexed with Tabix.
- Filters BAM files by the loaded intervals and writes the filtered reads to temporary BAM files.
- Merges all temporary BAM files into a final output BAM file.
- Keeps track of processed intervals to avoid redundant processing.

**Usage:**
```bash
merge-bam -b /path/to/bam1.bam /path/to/bam2.bam \
          -e /path/to/intervals.bed.gz \
          -o /path/to/output.bam \
          -t 4 \
          --temp-dir /path/to/temp \
          --overwrite
```

**Arguments:**
- `-b`, `--bam-files`: Paths to input BAM files to be filtered and merged.
- `-l`, `--bam-list`: Path to a text file containing a list of input BAM files to be filtered and merged.
- `-e`, `--bed-file`: Path to the BED file containing genomic intervals for filtering.
- `-o`, `--output-file`: Path to the output BAM file.
- `-t`, `--threads`: Number of threads to use for processing (default: 4).
- `--temp-dir`: Directory to store temporary files.
- `--overwrite`: Overwrite existing temporary files and output file if they exist.

### 3. Creating a Test BED File with `create-test-bed`

This command generates a BED file with a specified number of positions. It can optionally compress the file using bgzip and create a Tabix index.

**Functionality:**

- Generates random genomic positions within a specified range on a specified chromosome.
- Writes these positions to a BED file.
- Optionally compresses the BED file with bgzip and indexes it with Tabix.

**Usage:**
```bash
create-test-bed -n 1000 \
                -o /path/to/output.bed \
                -c chr1 \
                --bgzip
```

**Arguments:**

- `-n`, `--num-positions`: Number of genomic positions to generate.
- `-o`, `--output-file`: Path to the output BED file.
- `-c`, `--chromosome`: Chromosome to generate positions on.
- `--bgzip`: Compress the output BED file with bgzip and create a Tabix index.

### 4. Converting a Text File to a BED File with `bed-from-txt`

This command converts a text file with mutation data into a BED file that can be used with the `merge-bam` command. It can optionally compress the file using bgzip and create a Tabix index. This was made so that R dataframes using chr1_1000_A_T mutation IDs can be rapidly turned into bed files, specifically lists of germline positions.

**Functionality:**

- Parses the input text file to extract positions and mutations.
- Writes these as intervals to a BED file.
- Supports different indexing conventions (1-based or 0-based) based on the VCF_indexed_input flag.
- Compresses and indexes the BED file using bgzip and Tabix.

**Usage:**
```bash
bed-from-txt -i /path/to/mutations.txt \
             -o /path/to/output.bed \
             --vcf-indexed-input \
             --bgzip
```

**Arguments:**

- `-i`, `--input-file`: Path to the input text file.
- `-o`, `--output-file`: Path to the output BED file.
- `--vcf-indexed-input`: Use 1-based indexing for the input positions. If not specified, 0-based indexing is used.
- `--bgzip`: Compress the output BED file with bgzip and create a Tabix index.

----

## Dependencies
- `python==3.11.9`
- `pysam==0.22.1`




