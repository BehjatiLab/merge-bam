import argparse
import random
import pysam
import os


def generate_bed_file(
        output_file, num_positions=1, chromosome='chr1', position_list=None, bgzip=False
):
    """
    Generate a BED file with a specified number of positions, optionally
    bgzip and tabix index the file.

    :param output_file:
        Path to the output BED file.
    :param num_positions:
        Number of positions to generate.
    :param chromosome:
        Chromosome name to use (default is 'chr1').
    :param position_list:
        List of positions to use instead of generating random positions.
    :param bgzip:
        Boolean indicating whether to bgzip and index the BED file.
    """
    if position_list:
        num_positions = len(position_list)

    positions = set()

    # Generate unique start positions
    while len(positions) < num_positions:
        if position_list:
            start = position_list.pop()
        else:
            start = random.randint(1000, 1_000_000)
        positions.add(start)

    # Sort positions
    positions = sorted(positions)

    # Write positions to the BED file
    bed_file_path = output_file
    with open(bed_file_path, 'w') as bed_file:
        for start in positions:
            end = start + 1
            bed_file.write(f"{chromosome}\t{start}\t{end}\n")

    print(f"Generated BED file with {num_positions} positions: {bed_file_path}")

    # Optionally bgzip and tabix index the BED file
    if bgzip:
        bgzip_file_path = bed_file_path + '.gz'
        # Compress the BED file
        pysam.tabix_compress(bed_file_path, bgzip_file_path, force=True)
        # Create Tabix index
        pysam.tabix_index(bgzip_file_path, preset="bed", force=True)
        # Remove the uncompressed BED file
        os.remove(bed_file_path)
        print(f"Bgzipped and indexed BED file: {bgzip_file_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a test BED file with specified number of "
                    "positions."
    )
    parser.add_argument(
        '-n', '--num-positions', type=int, required=True,
        help="Number of positions to generate."
    )
    parser.add_argument(
        '-o', '--output-file', type=str, required=True,
        help="Path to the output BED file."
    )
    parser.add_argument(
        '-c', '--chromosome', type=str, default='chr1',
        help="Chromosome name (default: chr1)."
    )
    parser.add_argument(
        '--bgzip', action='store_true',
        help="Bgzip and tabix index the BED file."
    )

    args = parser.parse_args()

    # Generate BED file
    generate_bed_file(
        args.output_file, args.num_positions, args.chromosome, args.bgzip
    )


if __name__ == "__main__":
    main()
