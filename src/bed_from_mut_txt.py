#!/usr/bin/env python

import argparse
import pysam
import os
import heapq
from tempfile import NamedTemporaryFile


def create_bed_from_text(
        input_file, output_file, vcf_indexed_input=False, bgzip=False,
        chunk_size=100000
):
    """
    Create a sorted BED file from a text file with format "chr1_10144_T_C".

    :param input_file:
        Path to the input text file.
    :param output_file:
        Path to the output BED file.
    :param vcf_indexed_input:
        Boolean indicating whether the input is VCF-indexed.
    :param bgzip:
        Boolean indicating whether to bgzip and index the BED file.
    :param chunk_size:
        Number of lines to process at once (for sorting in chunks).
    """

    # Function to write sorted chunks to temporary files
    def write_sorted_chunk(lines, temp_files):
        lines.sort(key=lambda x: (x[0], x[1]))
        with NamedTemporaryFile(
                delete=False, mode='w',
                dir=os.path.dirname(output_file)
        ) as temp_file:
            for line in lines:
                temp_file.write('\t'.join(map(str, line)) + '\n')
            temp_files.append(temp_file.name)

    # Read input, process in chunks, and write sorted chunks
    temp_files = []
    current_chunk = []

    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            chrom, pos, ref, alt = line.split('_')
            pos = int(pos)

            if vcf_indexed_input:
                start = pos - 1  # VCF is 1-indexed, BED is 0-indexed
                end = pos
            else:
                start = pos  # BED start position
                end = pos + 1

            current_chunk.append((chrom, start, end, ref, alt))

            if len(current_chunk) >= chunk_size:
                write_sorted_chunk(current_chunk, temp_files)
                current_chunk = []

        # Write remaining lines in the last chunk
        if current_chunk:
            write_sorted_chunk(current_chunk, temp_files)

    # Merge sorted chunks into the final output file
    with open(output_file, 'w') as bed_file:
        files = [open(temp, 'r') for temp in temp_files]
        merged = heapq.merge(
            *files,
            key=lambda x: (x.split('\t')[0], int(x.split('\t')[1]))
        )
        for line in merged:
            bed_file.write(line)

    # Clean up temporary files
    for file in files:
        file.close()
    for temp in temp_files:
        os.remove(temp)

    print(f"Generated and sorted BED file: {output_file}")

    # Optionally bgzip and tabix index the BED file
    if bgzip:
        bgzip_file_path = output_file + '.gz'
        # Compress the BED file
        pysam.tabix_compress(output_file, bgzip_file_path, force=True)
        # Create Tabix index
        pysam.tabix_index(bgzip_file_path, preset="bed", force=True)
        # Remove the uncompressed BED file
        os.remove(output_file)
        print(f"Bgzipped and indexed BED file: {bgzip_file_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Create a sorted BGZIPPED and TABIX-indexed BED file from a"
                    " text file with lines like 'chr1_10144_T_C'."
    )
    parser.add_argument(
        '-i', '--input-file', type=str, required=True,
        help="Path to the input text file containing lines like 'chr1_10144_T_C'."
    )
    parser.add_argument(
        '-o', '--output-file', type=str, required=True,
        help="Path to the output BED file."
    )
    parser.add_argument(
        '--vcf-indexed-input', action='store_true',
        help="Treat the input positions as VCF-indexed (1-based)."
    )
    parser.add_argument(
        '--bgzip', action='store_true',
        help="Bgzip and tabix index the BED file."
    )
    parser.add_argument(
        '--chunk-size', type=int, default=100000,
        help="Number of lines to process per chunk for sorting."
    )

    args = parser.parse_args()

    # Create BED file from the input text file
    create_bed_from_text(
        args.input_file, args.output_file, args.vcf_indexed_input, args.bgzip,
        args.chunk_size
    )


if __name__ == "__main__":
    main()
