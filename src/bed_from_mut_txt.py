#!/usr/bin/env python

import argparse
import pysam
import os


def create_bed_from_text(
        input_file, output_file, vcf_indexed_input=False, bgzip=False
):
    """
    Create a BED file from a text file with format "chr1_10144_T_C".

    :param input_file:
        Path to the input text file.
    :param output_file:
        Path to the output BED file.
    :param vcf_indexed_input:
        Boolean indicating whether the input is VCF-indexed.
    :param bgzip:
        Boolean indicating whether to bgzip and index the BED file.
    """
    # Write positions directly to the BED file as they are processed
    with open(output_file, 'w') as bed_file:
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

                # Write each entry to the BED file immediately
                bed_file.write(f"{chrom}\t{start}\t{end}\t{ref}\t{alt}\n")

    print(f"Generated BED file: {output_file}")

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
        description="Create a BGZIPPED and TABIX-indexed BED file from a text "
                    "file with lines like 'chr1_10144_T_C'."
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

    args = parser.parse_args()

    # Create BED file from the input text file
    create_bed_from_text(
        args.input_file, args.output_file, args.vcf_indexed_input, args.bgzip
    )


if __name__ == "__main__":
    main()
