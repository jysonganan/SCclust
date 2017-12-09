#!/usr/bin/env python

import re
import sys

def main():
	read_length = int(sys.argv[1])
        ref_genome = sys.argv[2]
	if ref_genome == "hg":
		list_dir = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
				"chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
	elif ref_genome == "hgdm":
		list_dir = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
						"chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "dm_chrX",
						"dm_chr2L", "dm_chr2R", "dm_chr3L", "dm_chr3R", "dm_chr4"]
	
	for li in list_dir:
			thisChrom = li
			INFILE = open(thisChrom + ".fa", "r")
			if thisChrom == "chrY":
				INFILE.close()
				INFILE = open("chrY.psr.fa", "r")
			chr = []
			x = ""
			y = ""
			INFILE.readline()
			for x in INFILE:
				chr.append(x.rstrip())
			x = "".join(chr)
			y = x.upper()

			for i in range(len(y) - read_length + 1):
				print ">" + li + "." + str(i)
				print y[i:(i+read_length)]

			INFILE.close()

if __name__ == "__main__":
	main()
