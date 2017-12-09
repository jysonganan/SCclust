import re
import sys
import random
import string

def varbin_GC(chromFa_dir, data_dir, Nk):

	outfilename = data_dir + "/varbin.gc.content." + Nk + ".txt"
	bins = fileToArray(data_dir+ "/bin.boundaries." + Nk + ".sorted.txt", False)

	prevChrom = ""
	INFILE = False
	OUTFILE = open(outfilename, "w")
	OUTFILE.write("bin.chrom\tbin.start.chrompos\tbin.start.abspos\tbin.end.chrompos\tbin.length\tmappable.positions\tgc.content\n")

	for arow in bins:
		thisChrom = arow[0].strip()
		if thisChrom == "23":
			thisChrom = "X"
		if thisChrom == "24":
			thisChrom = "Y"
		thisChrom = thisChrom
		thisStart = int(arow[1].strip())
		thisEnd = int(arow[3].strip())

		if thisChrom == prevChrom:
			pass
		else:
			if INFILE:
				INFILE.close()
			INFILE = open(chromFa_dir + "/" + thisChrom + ".fa", "r")
			if thisChrom == "chrY":
				INFILE.close()
				INFILE = open(chromFa_dir + "/chrY.psr.fa", "r")
			chr = []
			x = ""
			y = ""
			INFILE.readline()
			for x in INFILE:
				chr.append(x.rstrip())
			x = "".join(chr)
			y = x.upper()
			print "after read " + thisChrom
			prevChrom = thisChrom

		gcContent = float(len(re.findall("[CG]", y[(thisStart):(thisEnd+1)]))) / float(len(re.findall("[ACGT]", y[(thisStart):(thisEnd+1)])))
		OUTFILE.write("\t".join(arow))
		OUTFILE.write("\t")
		OUTFILE.write(str(gcContent))
		OUTFILE.write("\n")
		OUTFILE.flush()

	INFILE.close()
	OUTFILE.close()


def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")

	ra = []
	if skipFirst:
		input.readline()
	for x in input:
		arow = x.rstrip().split("\t")
		ra.append(arow)

	input.close()
	return(ra)
