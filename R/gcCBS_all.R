#'Generate the segmented profile for each cell.
#'
#'Generate the segmented profile for each cell in the input directory with CBS. (Including GC correction)
#'@param input_file_dir The directory that contains the bin count files for all cells. For example, one bin count file for a cell is CJA1023.varbin.20k.txt
#'@param Nk Default value: 5k. It denotes the number of bins. 5000 = 5k, 10000 = 10k, ...
#'@param gc The GC content table.
#'@param alpha Parameter alpha defined in DNAcopy package.
#'@param nperm Parameter nperm defined in DNAcopy package.
#'@param undo.SD Parameter undo.SD defined in DNAcopy package.
#'@param min.width Parameter min.width defined in DNAcopy package.
#'@param method Default: "multiplier" to transform the ratio data to integer CN state. When genome is hg-dm hybrid, set method as "dmploidies" to use dm plodies as reference.
#'@param genome Default: "hg" (Human genome). hg-dm hybrid genome: "hgdm".
#'@return The list containing seg.quantal and ratio.quantal matrix for all cells in the input_file_dir.
#'@export

cbs.segment_all <- function(input_file_dir, Nk = "5k", gc, alpha, nperm, undo.SD, min.width, method = "multiplier", genome = "hg") {

  if ((genome == "hg") & (method = "dmploidies")) {
    stop("method dmploidies only works when genome is hgdm!")
  }

  chrom.numeric <- substring(gc$bin.chrom, 4)
  chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
  chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
  if (genome == "hgdm"){
    chrom.numeric[which(gc$bin.chrom == "dm_chrX")] <- "25"
    chrom.numeric[which(gc$bin.chrom == "dm_chr2L")] <- "26"
    chrom.numeric[which(gc$bin.chrom == "dm_chr2R")] <- "27"
    chrom.numeric[which(gc$bin.chrom == "dm_chr3L")] <- "28"
    chrom.numeric[which(gc$bin.chrom == "dm_chr3R")] <- "29"
    chrom.numeric[which(gc$bin.chrom == "dm_chr4")] <- "30"
  }
  chrom.numeric <- as.numeric(chrom.numeric)
  pattern.file <- paste("\\.varbin.", Nk, ".txt$", sep = "")
  file.names <- list.files(path = input_file_dir, pattern = pattern.file)

  input_list <- list()
  for (j in 1:length(file.names)){
    input_list[[j]] <- read.table(paste(input_file_dir, "/", file.names[j], sep=""), header = T)
    input_list[[j]]$chrom <- chrom.numeric
    input_list[[j]]$gc.content <- gc$gc.content
  }

  output_list <- lapply(input_list, cbs.segment_1, gc, alpha, nperm, undo.SD, min.width)
  output_mat = output_mat_1 <- output_list[[1]][,c("chrom","chrompos","abspos")]

  for (j in 1:length(file.names)){
    output_mat <- cbind(output_mat, output_list[[j]][,"seg.quantal"])
    output_mat_1 <- cbind(output_mat_1, output_list[[j]][,"ratio.quantal"])
  }

  colnames(output_mat) = colnames(output_mat_1) <-c("chrom","chrompos","abspos",sub("\\..*","",file.names))
  return(list(seg.quantal = output_mat, ratio.quantal = output_mat_1))
}



cbs.segment_1 <- function(bin_mat, gc, alpha, nperm, undo.SD, min.width){
  a <- bin_mat$bincount + 1
  bin_mat$ratio <- a / mean(a)
  bin_mat$lowratio <- lowess.gc(bin_mat$gc.content, bin_mat$ratio)
  #thisRatioNobig <- bin_mat[-bad[, 1], ]
  set.seed(25)
  CNA.object <- DNAcopy::CNA(log(bin_mat$lowratio, base=2), bin_mat$chrom,
                    bin_mat$chrompos, data.type="logratio")
  smoothed.CNA.object <- DNAcopy::smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- DNAcopy::segment(smoothed.CNA.object, alpha=alpha,
                                                  nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2)
  thisShort <- segment.smoothed.CNA.object[[2]]


  m <- matrix(data=0, nrow=nrow(bin_mat), ncol=1)
  prevEnd <- 0
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.mark[i]
    m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
    prevEnd = thisEnd
  }
  cbs.long.nobad <- m[, 1]

  #####  NEW STUFF  also check min.width=2 above

  workShort <- thisShort
  workShort$segnum <- 0
  workShort$seg.start <- 0
  workShort$seg.end <- 0
  prevEnd <- 0
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.mark[i]
    workShort$seg.start[i] <- thisStart
    workShort$seg.end[i] <- thisEnd
    workShort$segnum[i] <- i
    prevEnd = thisEnd
  }

  discardSegments <- TRUE
  while (discardSegments) {
    orderShort <- workShort[order(workShort$num.mark, abs(workShort$seg.mean)), ]
    if (orderShort[1, "num.mark"] < min.width) {
      workShort <- remove.segment(workShort, orderShort[1, "segnum"], bin_mat, undo.SD)
    } else {
      discardSegments <- FALSE
    }
  }

  workShort <- sdundo.all(workShort, bin_mat, undo.SD)
  thisShort <- workShort

  #####  END NEW STUFF

  m <- matrix(data=0, nrow=nrow(bin_mat), ncol=1)
  prevEnd <- 0
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.mark[i]
    m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
    prevEnd = thisEnd
  }

  bin_mat$seg.mean.LOWESS <- m[, 1]

  thisGrid <- seq(1.5, 5.5, by=0.05)
  thisOuter <- bin_mat$seg.mean.LOWESS %o% thisGrid
  thisOuterRound <- round(thisOuter)
  thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
  thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
  thisMultiplier <- thisGrid[which.min(thisOuterColsums)]

  bin_mat$seg.quantal <- bin_mat$seg.mean.LOWESS * thisMultiplier
  bin_mat$ratio.quantal <- bin_mat$lowratio * thisMultiplier ###

  if (method == "dmploidies"){
    dmbinID <- which(bin_mat$chrom%in%(25:30))
    median_ratio <- median(bin_mat$lowratio[dmbinID], na.rm = TRUE)
    theMultiplier <- 2/median_ratio
    bin_mat$seg.quantal <- bin_mat$seg.mean.LOWESS * theMultiplier
    bin_mat$ratio.quantal <- bin_mat$lowratio * theMultiplier}

  return(bin_mat)
}
