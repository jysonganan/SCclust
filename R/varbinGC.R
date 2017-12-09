#'Compute GC content.
#'
#'Given the bin boundaries, compute the GC content in each bin.
#'@param chromFa_dir File directory that contains the reference genome.
#'@param data_dir File directory. The directory that contains the files of mappable regions, bin boundaries etc.
#'The GC content file (varbin.gc.content.Nk.txt) will also be output to this directory.
#'@param Nk Default value: 5k. It denotes the number of bins. 5000 = 5k, 10000 = 10k, ...
#'@export

varbinGC <- function(chromFa_dir, data_dir, Nk = "5k"){

  rPython::python.load(system.file("varbin_GC.py", package = "SCclust"))
  rPython::python.call("varbin_GC", chromFa_dir, data_dir, Nk)
}
