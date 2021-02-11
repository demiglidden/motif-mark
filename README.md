# motif-mark
Updated 10 February 2021, by Demi Glidden

This script visualizes nucleotide motifs on genes with a single exon. It currently supports any number of genes, of any length, and 5 motifs.

Please ensure that the genes are in fasta format with exons in uppercase and introns in lowercase. The motifs should be passed as a text file with each motif on its own line.

This is a Python3 script. Please install the following packages before executing script:
    pycairo
    argparse
    re

General Notes:
All features and genes are to scale (both inter-gene and intra-gene)
Exons are marked as black boxes.
Motifs are slightly transparent so that you can see where they overlap if they do.


