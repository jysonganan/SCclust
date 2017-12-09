
preprocess_annot <- function(annot){
  annot[annot[,"bin.chrom"] == "chrX", "bin.chrom"] <- "chr23"
  annot[annot[,"bin.chrom"] == "chrY", "bin.chrom"] <- "chr24"

  annot[,"bin.chrom"] <- as.numeric(substring(annot[,"bin.chrom"], first=4,
                                                  last=nchar(annot[,"bin.chrom"])))

  return(annot)
}


#Identify the centromere areas

#
# drop_areas <- function(cyto){
#   centromere = c("p11","q11")
#   cyto[cyto[,1] == "chrX",1] <- "chr23"
#   cyto[cyto[,1] == "chrY",1] <- "chr24"
#   cyto[,1] <- as.numeric(substring(cyto[,1], first = 4, last = nchar(cyto[,1])))
#   cyto <- cyto[order(cyto[,1]),]
#   centroleft <- cyto[grep(centromere[1], cyto[,4]),]
#   centroright <- cyto[grep(centromere[2], cyto[,4]),]
#   centroleft <- centroleft[match(unique(centroleft[,1]), centroleft[,1]),]
#   centroright <- centroright[nrow(centroright):1,]
#   centroright <- centroright[match(unique(centroright[,1]), centroright[,1]),]
#   centroright <- centroright[nrow(centroright):1,]
#   dropareas <- cbind(centroleft[,c(1,2)], centroright[,3])
#   dimnames(dropareas)[[2]] <- c("chrom", "from", "to")
#   return(dropareas)
# }


#'Preprocess the segmented copy number profiles
#'
#'Input the segmented copy number profiles of some single cells,
#'generate the breakpoint table which provides the genomic location
#'and segmented copy number values for all cells.
#'@param segfile The segmented copy number profiles of some single cells.
#'@param gc The GC content table.
#'@param eviltwins Bad cells. NULL or a vector of cell names.
#'@param ploidies Logical. If TRUE (Default), the ploidies and homoloss information are taken into consideration.
#'@return A list for two tables: breakpoint_table and ploidies_table. Each row in breakpoint_table correspond
#'to one segment divided by breakpoints in one cell. The columns in breakpoint_table: cell names (profid),
#'chromosome information(chrom, chromstart, chromend), absolute genomic location (abstart, absend), the bin location
#'(segstarts, segends), count of bins in the segement (segbins) and the copy number value for the segment (segvals, cvals).
#'Note: cvals are derived by incoporating ploidies.
#'@export


preprocess_segfile <- function(segfile, gc, eviltwins = NULL, ploidies = TRUE){
  annot <- preprocess_annot(gc)

  ## generate annotation file with info: chrom, chromstart, chromend, abstart, absend
  #make sure the location of 5000 bins in segfile/annot are same
  annotXk <- cbind(annot[match(4e8*segfile[,"chrom"] + segfile[,"chrompos"],
                               4e8*annot[,"bin.chrom"] + annot[,"bin.start"]),
                         c("bin.chrom", "bin.start", "bin.end"),], segfile[,"abspos"])
  dimnames(annotXk)[[2]] <- c("chrom", "chromstart", "chromend", "abstart")
  annotXk <- cbind(annotXk, annotXk[,"chromend"] +
                     annotXk[,"abstart"] - annotXk[,"chromstart"])
  dimnames(annotXk)[[2]][ncol(annotXk)] <- "absend"


  ## generate the breakpoint matrix (logic matrix, True for a breakpoint)
  # the rows that are the start of a chr, are all denoted as TRUE(breakpoint)
  a <- round(as.matrix(segfile[,-(1:3)]))
  b <- a - rbind(matrix(nrow=1, ncol=ncol(a), data=0), a[-nrow(a),])
  breakpoint_mat <- (b != 0) # True for a breakpoint in each cell
  breakpoint_mat[match(unique(segfile[,"chrom"]), segfile[,"chrom"]),] <- TRUE


  ## breakpoint_table: each row corresponds to a break point
  # cols: (1)the cell name: profid (2)location: chrom chromstart abstart chromend absend
  #segstarts segends segbins (3)the CN value: segvals
  segstarts <- row(a)[breakpoint_mat]
  segvals <- a[breakpoint_mat]
  profid <- dimnames(a)[[2]][cumsum(segstarts == 1)]
  segends <- c(segstarts[-1] - 1, nrow(a))
  segends[segends == 0] <- nrow(a)
  segbins <- segends - segstarts + 1
  segloc <- cbind(annotXk[segstarts, c("chrom", "chromstart", "abstart")],
                  annotXk[segends, c("chromend", "absend")])
  breakpoint_table <- data.frame(I(profid), segloc, segstarts, segends, segbins, segvals)


  ## ploidies table
  #modeofmodes.R
  ploidymod <- apply(a[annotXk[,"chrom"] < 23,], 2, getmode)
  ploidymed <- apply(a[annotXk[,"chrom"] < 23,], 2, median)
  ploidychromod <- apply(a[annotXk[,"chrom"] < 23,], 2, modeofmodes,
                         otherlabel = annotXk[annotXk[,"chrom"] < 23, "chrom"],
                         tiebreaker = 2, tiebreakerside = "greater")
  homoloss <- colSums(!a[annotXk[,"chrom"] < 23,])/sum(annotXk[,"chrom"] < 23)
  ploidies_table <- cbind(ploidymed, ploidymod, ploidychromod, homoloss)

  breakpoint_table <- cbind(breakpoint_table, breakpoint_table[,"segvals"] -
                              ploidies_table[breakpoint_table[,"profid"], "ploidychromod"])
  dimnames(breakpoint_table)[[2]][ncol(breakpoint_table)] <- "cvals"

  ## Ploidies taken into consideration
  if (ploidies == TRUE){
    breakpoint_table <-
      breakpoint_table[breakpoint_table[,"profid"]%in%dimnames(ploidies_table)[[1]][ploidies_table[,"homoloss"] < 0.01],]
  }

  ## delete bad cells
  goodCells <- setdiff(dimnames(ploidies_table)[[1]], eviltwins)
  breakpoint_table <- breakpoint_table[breakpoint_table[,"profid"]%in%goodCells,]

  return(list(breakpoint_table = breakpoint_table, ploidies_table = ploidies_table))
}

## tshort-->breakpoint_table
## ploidies-->ploidies_table
