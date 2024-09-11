#!/usr/bin/env python

import pysam
import argparse
from concurrent.futures import ProcessPoolExecutor
import os
import json


def filter_bam_by_interval(
        bam_file, intervals, temp_bam_path, processed_intervals
):
    """
    Filter reads from a single BAM file that intersect with regions in the
    provided intervals and append to output BAM.

    :param bam_file:
        Path to the BAM file to be filtered.
    :param intervals:
        List of intervals (chrom, start, end) to filter reads by.
    :param temp_bam_path:
        Path to the temporary output BAM file.
    :param processed_intervals:
        A set of intervals that have already been processed.
    """
    with pysam.AlignmentFile(bam_file, 'rb') as bam, pysam.AlignmentFile(
            temp_bam_path, 'wb', header=bam.header
    ) as out_bam:
        for interval in intervals:
            chrom, start, end = interval
            if interval in processed_intervals:
                continue
            for read in bam.fetch(chrom, start, end):
                out_bam.write(read)


# TODO add support for other types of tabix indexed files such as VCF.
def load_intervals(file):
    """
    Load intervals from a BED file or a tabix-indexed gzipped BED file (.bed.gz).

    :param file:
        Path to the BED file (either plain .bed or .bed.gz with .tbi index).
    :return:
        List of intervals (chrom, start, end).
    """
    intervals = []

    if file.endswith('.gz'):
        # Check if the .tbi index exists for the gzipped BED file
        if not os.path.exists(file + '.tbi'):
            raise FileNotFoundError(
                f"Tabix index (.tbi) not found for {file}"
            )

        # Use pysam.TabixFile to read the gzipped BED file
        with pysam.TabixFile(file) as tbx:
            for line in tbx.fetch():
                chrom, start, end = line.strip().split()[:3]
                intervals.append((chrom, int(start), int(end)))
    else:
        # Read plain BED file
        with open(file, 'r') as bed:
            for line in bed:
                chrom, start, end = line.strip().split()[:3]
                intervals.append((chrom, int(start), int(end)))

    return intervals


def load_processed_intervals(output_bam):
    """
    Load processed intervals from an existing output BAM file or tracking file
    if available.

    :param output_bam:
        Path to the output BAM file.
    :return:
        Set of processed intervals.
    """
    tracking_file = output_bam + '.processed.json'
    if os.path.exists(tracking_file):
        with open(tracking_file, 'r') as file:
            return set(tuple(interval) for interval in json.load(file))
    return set()


def save_processed_intervals(output_bam, processed_intervals):
    """
    Save processed intervals to a tracking file.

    :param output_bam:
        Path to the output BAM file.
    :param processed_intervals:
        Set of processed intervals.
    """
    tracking_file = output_bam + '.processed.json'
    with open(tracking_file, 'w') as file:
        json.dump(list(processed_intervals), file)


def process_bam_files(
        bam_files, bed_file, output_bam, num_threads, overwrite, temp_dir
):
    intervals = load_intervals(bed_file)
    processed_intervals = set()

    # Check if output exists and handle accordingly
    if os.path.exists(output_bam) and not overwrite:
        print(f"Resuming from existing output: {output_bam}")
        processed_intervals = load_processed_intervals(output_bam)
    elif os.path.exists(output_bam) and overwrite:
        print(f"Overwriting existing output: {output_bam}")
        os.remove(output_bam)
        if os.path.exists(output_bam + '.bai'):
            os.remove(output_bam + '.bai')
        if os.path.exists(output_bam + '.processed.json'):
            os.remove(output_bam + '.processed.json')
    else:
        print(f"Creating new output: {output_bam}")

    # Create the temporary directory if it does not exist
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Create paths for temporary BAM files with sequential naming
    temp_bam_paths = [
        os.path.join(temp_dir, f"{i}.temp.bam") for i in range(len(bam_files))
    ]

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = []

        # Submit jobs for parallel processing
        for i, (bam_file, temp_bam_path) in enumerate(
                zip(bam_files, temp_bam_paths)
        ):
            if os.path.exists(temp_bam_path):
                print(
                    f"Found existing temporary BAM: {temp_bam_path}, "
                    f"skipping processing for this BAM file."
                )
                continue

            print(f"Processing {temp_bam_path}")
            future = executor.submit(
                filter_bam_by_interval, bam_file, intervals,
                temp_bam_path, processed_intervals
            )
            futures.append(future)

        # Wait for all processes to finish
        for future in futures:
            future.result()

    # Merge temporary BAM files into final output
    pysam.merge(output_bam, *temp_bam_paths)
    pysam.index(output_bam)

    # Update processed intervals
    save_processed_intervals(output_bam, set(intervals))


def main():
    parser = argparse.ArgumentParser(
        description="Filter BAM files by BED intervals and create a new BAM "
                    "file."
    )
    parser.add_argument(
        '-b', '--bam-files', nargs='+',
        help="Paths to BAM files, or a text file containing paths.",
        required=False
    )
    parser.add_argument(
        '-l', '--bam-list',
        help="Text file containing paths to BAM files, one per line.",
        required=False
    )
    parser.add_argument(
        '-e', '--bed-file',
        help="Path to BED file containing regions of interest.", required=True
    )
    parser.add_argument(
        '-o', '--output-bam', help="Path to the output BAM file.",
        required=True
    )
    parser.add_argument(
        '-t', '--threads', type=int, default=4,
        help="Number of threads to use for processing."
    )
    parser.add_argument(
        '--overwrite', action='store_true',
        help="Overwrite the existing output BAM file if it exists."
    )
    parser.add_argument(
        '--temp-dir',
        help="Directory for temporary BAM files. This directory will not be "
             "automatically deleted.",
        required=True
    )

    args = parser.parse_args()

    # Load BAM files from command line or list file
    bam_files = args.bam_files if args.bam_files else []
    if args.bam_list:
        with open(args.bam_list, 'r') as bam_list_file:
            bam_files.extend([line.strip() for line in bam_list_file])

    if not bam_files:
        parser.error(
            "No BAM files provided. Use --bam-files or --bam-list to specify "
            "input BAM files."
        )

    # Process BAM files
    process_bam_files(
        bam_files, args.bed_file, args.output_bam, args.threads, args.overwrite,
        args.temp_dir
    )


if __name__ == "__main__":
    main()

