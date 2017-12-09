import re
import sys
import random
import string

def main():
	thisState = "out"
	prevChrom = ""
	prevStart = ""
	prevEnd = ""
	
	for aline in sys.stdin:
		arow = aline.rstrip().split("\t")
		if len(arow) < 11:
			continue
			
		if thisState == "out":
			if arow[1] == "0":
				prevChrom = arow[2]
				prevStart = arow[3]
				prevEnd = str(int(arow[3]) + 1)
				thisState = "in"
		else:
			if arow[1] == "0":
				if arow[2] == prevChrom:
					prevEnd = str(int(arow[3]) + 1)
				else:
					print prevChrom, "\t", prevStart, "\t", prevEnd
					prevChrom = arow[2]
					prevStart = arow[3]
					prevEnd = str(int(arow[3]) + 1)
			else:
				print prevChrom, "\t", prevStart, "\t", prevEnd
				thisState = "out"
		
	if thisState == "in":
		print prevChrom, "\t", prevStart, "\t", prevEnd


if __name__ == "__main__":
	main()
