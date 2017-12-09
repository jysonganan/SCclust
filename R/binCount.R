#'make bin count files
#'
#'Make bin count files from SAM files
#'@param SAM_dir The directory that contains the sorted SAM files with duplicates removed. The SAM file names: cellname.rmdup.sam.
#'e.g. CJA1023.rmdup.sam
#'The bin count files will also be saved in SAM_dir. The output bin count file names: cellname.varbin.Nk.txt, cellname.varbin.Nk.stats.txt
#'e.g. CJA1023.varbin.20k.txt, CJA1023.varbin.20k.stats.txt
#'@param data_dir The directory that contains the bin boundary file; also contains mappable region files such as
#' chrom.mappable.txt, mappable.regions.sorted.txt, chrom.sizes.txt
#'@param cell_name The cell name.
#'@param Nk Default value: 5k. It denotes the number of bins. 5000 = 5k, 10000 = 10k, ...
#'@export


bin_counts <- function(SAM_dir, data_dir, cell_name, Nk = "5k"){

  rPython::python.load(system.file("bin_count.py", package = "SCclust"))
  rPython::python.call("bin_count", SAM_dir, data_dir, cell_name, Nk)

  file1 <- paste(cell_name, ".varbin.", Nk,".txt", sep = "")
  file2 <- paste(cell_name, ".varbin.", Nk, ".stats.txt", sep = "")


  cmd1 <- paste("mv ", SAM_dir, "/varbin.txt ", SAM_dir, "/", file1, sep = "")
  cmd2 <- paste("mv ", SAM_dir, "/varbin.stats.txt ", SAM_dir, "/", file2, sep = "")

  system(cmd1)
  system(cmd2)
}
