import os
import pysam
import shutil
from src.merge_bam import process_bam_files
from src.create_position_bed import generate_bed_file


def create_dummy_bam(file_path, reads_at_positions):
    """
    Create a dummy BAM file with reads at specified positions.

    :param file_path: Path to the output BAM file.
    :param reads_at_positions: A dictionary with positions as keys and number of reads as values.
    """
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }

    with pysam.AlignmentFile(file_path, "wb", header=header) as outf:
        for pos, num_reads in reads_at_positions.items():
            for i in range(num_reads):
                read = pysam.AlignedSegment()
                read.query_name = f"read_{pos}_{i}"
                read.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
                read.flag = 0
                read.reference_id = 0  # 'chr1'
                read.reference_start = pos
                read.mapping_quality = 60
                read.cigar = [(0, len(read.query_sequence))]  # 36M
                read.query_qualities = pysam.qualitystring_to_array(
                    "F" * len(read.query_sequence)
                )
                outf.write(read)

            # Add a read that stops just short of the position
            read = pysam.AlignedSegment()
            read.query_name = f"read_{pos}_before"
            read.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
            read.flag = 0
            read.reference_id = 0  # 'chr1'
            read.reference_start = pos - len(read.query_sequence)
            read.mapping_quality = 60
            read.cigar = [(0, len(read.query_sequence))]  # 36M
            read.query_qualities = pysam.qualitystring_to_array(
                "F" * len(read.query_sequence)
            )
            outf.write(read)

            # Add a read that starts just after the position
            read = pysam.AlignedSegment()
            read.query_name = f"read_{pos}_after"
            read.query_sequence = "GCTTAGCTAGCTACCTATATCTTGGTCTTGGCCGA"
            read.flag = 0
            read.reference_id = 0  # 'chr1'
            read.reference_start = pos + 1
            read.mapping_quality = 60
            read.cigar = [(0, len(read.query_sequence))]  # 36M
            read.query_qualities = pysam.qualitystring_to_array(
                "F" * len(read.query_sequence)
            )
            outf.write(read)

    pysam.sort("-o", file_path, file_path)

    pysam.index(file_path)


def test_merge_bams():
    # Define the positions and read counts for each BAM file
    bam1_positions = {20001:  7, 50000: 15, 70000: 10}
    bam2_positions = {20001: 20, 50000:  5, 70000:  8}

    # Create the dummy BAM files
    bam1_path = "test_bam1.bam"
    bam2_path = "test_bam2.bam"
    create_dummy_bam(bam1_path, bam1_positions)
    create_dummy_bam(bam2_path, bam2_positions)

    # Create the test BED file
    generate_bed_file(
        "test_intervals.bed", position_list=[20001, 50000]
    )

    # Define the output BAM path
    merged_bam_path = "merged.bam"

    # Define the temporary directory for processing
    temp_dir = "tmp_test_dir"
    os.makedirs(temp_dir, exist_ok=True)

    # Use process_bam_files function from merge_bam module to merge the BAM
    # files
    process_bam_files(
        bam_files=[bam1_path, bam2_path],
        bed_file="test_intervals.bed",
        output_bam=merged_bam_path,
        num_threads=4,
        overwrite=True,
        temp_dir=temp_dir
    )

    # Check the merged BAM file to ensure correctness
    expected_counts = {20001: 27, 50000: 20, 70000: 0}
    observed_counts = {20001: 0, 50000: 0, 70000: 0}

    with pysam.AlignmentFile(merged_bam_path, "rb") as merged_bam:
        for read in merged_bam:
            pos = read.reference_start
            if pos in observed_counts:
                observed_counts[pos] += 1

    # Verify that the merge is correct
    assert observed_counts == expected_counts, (
        f"Test failed: {observed_counts} != {expected_counts}"
    )
    print("Merge test passed!")

    # Clean up test files
    for path in [bam1_path, bam2_path, merged_bam_path]:
        os.remove(path)
        os.remove(path + '.bai')
    os.remove("test_intervals.bed")
    os.remove("merged.bam.processed.json")
    shutil.rmtree(temp_dir)


if __name__ == "__main__":
    test_merge_bams()
