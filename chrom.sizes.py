#!/usr/bin/env python

import sys
import time
from operator import itemgetter


def main():


	file_dir = sys.argv[1]
	ref_genome = sys.argv[2]

	outfilename = "/" + file_dir + “/chrom.sizes.txt"

	if ref_genome = "hg":
		list_dir = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
				"chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX","chrY","chrM"]
	elif ref_genome == "hgdm":
		list_dir = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
						"chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM", "dm_chrX",
						"dm_chr2L", "dm_chr2R", "dm_chr3L", "dm_chr3R", "dm_chr4","dm_chrM"]



	INFILE = 0
	OUTFILE = open(outfilename, "w")

	abspos = 0
	for x in list_dir:
		infilename = x + ".fa"
		print infilename
		if INFILE:
			INFILE.close()
		INFILE = open(infilename, "r")
		INFILE.readline()
		chromlength = 0
		for y in INFILE:
			aline = y.rstrip()
			chromlength += len(aline)
		OUTFILE.write(x)
		OUTFILE.write("\t")
		OUTFILE.write(str(chromlength))
		OUTFILE.write("\t")
		OUTFILE.write(str(abspos))
		OUTFILE.write("\n")
		abspos += chromlength
		
	INFILE.close()
	OUTFILE.close()


if __name__ == "__main__":
	main()
