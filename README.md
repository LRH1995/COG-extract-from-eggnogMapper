Some tips about search-and-extract-allCOG-for-conserved-single-copy-gene.py
aims and attentions： This script was designed to extract protein sequences that belong to a COG number from the user-offered ".faa" and ".emapper.annotations" files. For multicopy gene families, this script will compare their evalue and the one with the smallest evalue will be retained.

Step 1: Opening the script, and slide down to the last line. The following parameters must be changed before every run. This script will change our input files, so we must copy them before every run.
inputfile1=genome.faa
inputfile2=genome.emapper.annotations
output=/the/path/that/you/hope/the/results/store/

Step 2: Performing this script, you will find many new files named by COG number.
