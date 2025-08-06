#!/usr/bin/env python3

"""
Process BAM files for X chromosome genotyping analysis
Extracts X chromosome reads and renames cell barcodes from -1 to -RNA
"""

import pysam
import sys

def main():
    # Read input/output BAM paths from command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python process_bam.py <input_bam> <output_bam>")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]

    print(f"Processing {input_bam} -> {output_bam}")

    # Open the input and output BAM files
    try:
        in_bam = pysam.AlignmentFile(input_bam, 'rb')
        out_bam = pysam.AlignmentFile(output_bam, 'wb', template=in_bam)
    except Exception as e:
        print(f"ERROR: Failed to open BAM files: {e}")
        sys.exit(1)

    # Process reads for RNA-seq - extract chrX and rename barcodes
    processed_reads = 0
    for aln in in_bam.fetch("chrX"):
        if aln.has_tag("CB"):  # CB is the cell barcode tag
            old_bc = aln.get_tag('CB')
            # Rename the barcode with '-RNA' suffix
            new_bc = old_bc.replace("-1", "-RNA")
            aln.set_tag('CB', new_bc)
            out_bam.write(aln)
            processed_reads += 1

    # Close files
    in_bam.close()
    out_bam.close()

    print(f"Processing completed. Processed {processed_reads} reads from chrX")

if __name__ == "__main__":
    main()
