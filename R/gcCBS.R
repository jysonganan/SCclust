lowess.gc <- function(jtkx, jtky) {
  jtklow <- lowess(jtkx, log(jtky), f=0.05)
  jtkz <- approx(jtklow$x, jtklow$y, jtkx)
  return(exp(log(jtky) - jtkz$y))
}

remove.segment <- function( rsShort, rsSegnum, ratioData, sd.undo ) {

  ############################################
  ## deciding the appending location (left or right)
  ############################################

  appendLeft <- TRUE
  checkSdundo <- FALSE

  if (rsSegnum == 1) {
    appendLeft <- FALSE
  } else {
    if (rsSegnum == nrow(rsShort)) {
      appendLeft <- TRUE
    } else {
      rightIndex <- rsSegnum + 1
      leftIndex <- rsSegnum - 1

      if (rsShort[rightIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
        appendLeft <- TRUE
      } else {
        if (rsShort[leftIndex, "chrom"] != rsShort[rsSegnum, "chrom"]) {
          appendLeft <- FALSE
        } else {
          if (abs(rsShort[leftIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"]) < abs(rsShort[rightIndex, "seg.mean"] - rsShort[rsSegnum, "seg.mean"])) {
            appendLeft <- TRUE
            checkSdundo <- TRUE
          } else {
            appendLeft <- FALSE
            checkSdundo <- TRUE
          }}}
    }}

  appendIndex <- 99999999
  if (appendLeft) {
    appendIndex <- rsSegnum - 1
  } else {
    appendIndex <- rsSegnum + 1
  }

  ## in the short table(each row is one segment),append the short segment to its left/right segment.

  ############################################
  ## make the new short table after appending the short segment.
  ############################################

  tempShort <- rsShort
  newLocStart <- -1
  newLocEnd <- -1
  if (appendLeft) {
    tempShort[appendIndex, "loc.end"] <- tempShort[rsSegnum, "loc.end"]
    tempShort[appendIndex, "seg.end"] <- tempShort[rsSegnum, "seg.end"]
  } else {
    tempShort[appendIndex, "loc.start"] <- tempShort[rsSegnum, "loc.start"]
    tempShort[appendIndex, "seg.start"] <- tempShort[rsSegnum, "seg.start"]
  }

  tempShort[appendIndex, "num.mark"] <- tempShort[appendIndex, "num.mark"] + tempShort[rsSegnum, "num.mark"]
  ## updating seg mean
  tempShort[appendIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[appendIndex, "seg.start"]:tempShort[appendIndex, "seg.end"]], base=2))

  # cat("append", tempShort[appendIndex, "chrom"], tempShort[appendIndex, "loc.start"], tempShort[appendIndex, "loc.end"],
  #     tempShort[appendIndex, "num.mark"], tempShort[appendIndex, "seg.mean"], tempShort[appendIndex, "seg.start"], tempShort[appendIndex, "seg.end"], "\n")

  tempShort <- tempShort[-rsSegnum, ]
  tempShort$segnum <- seq(1:nrow(tempShort))



  if (checkSdundo) {
    thisSd <- -1
    if (appendLeft) {
      leftIndex <- appendIndex
      rightIndex <- appendIndex + 1
    } else {
      leftIndex <- appendIndex - 2
      rightIndex <- appendIndex - 1
    }
    #thisSd <- sd(ratioData[tempShort$seg.start[leftIndex]:tempShort$seg.start[rightIndex], "lowratio"])

    ## sd of log2 (reason to use sqrt(2) here)GC normalized bin count data (the data before segmentation)
    thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)

    ############################################
    ## after appending, check the new splitting.
    ############################################

    ## check whether the new splitting a real change point

    if (abs(tempShort$seg.mean[leftIndex] - tempShort$seg.mean[rightIndex]) < (sd.undo * thisSd) ) {

      # cat("left", tempShort[leftIndex, "chrom"], tempShort[leftIndex, "loc.start"],
      #     tempShort[leftIndex, "loc.end"], tempShort[leftIndex, "num.mark"], tempShort[leftIndex, "seg.mean"],
      #     tempShort[leftIndex, "seg.start"], tempShort[leftIndex, "seg.end"], "\n")
      # cat("right", tempShort[rightIndex, "chrom"], tempShort[rightIndex, "loc.start"],
      #     tempShort[rightIndex, "loc.end"], tempShort[rightIndex, "num.mark"], tempShort[rightIndex, "seg.mean"],
      #     tempShort[rightIndex, "seg.start"], tempShort[rightIndex, "seg.end"], "\n")

      ##  remove changepoint (combine)
      tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
      tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
      tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
      ## updating seg mean
      tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base=2))
      tempShort <- tempShort[-rightIndex, ]
      tempShort$segnum <- seq(1:nrow(tempShort))
    }
  }




  return(tempShort)
}











sdundo.all <- function (sdShort, ratioData, sd.undo) {

  tempShort <- sdShort
  thisSd <- mad(diff(ratioData[, "lowratio"])) / sqrt(2)

  while ( TRUE ) {

    chrom <- tempShort$chrom
    chrom.shift <- c(tempShort$chrom[-1], tempShort$chrom[1])
    breakpoints <- which(chrom == chrom.shift)
    #breakpoints intra(inside) chrom

    if (length(breakpoints) < 1) {
      break
    }

    breakpoints.shift <- breakpoints + 1
    undo.breakpoints <- breakpoints[which(abs(tempShort$seg.mean[breakpoints] - tempShort$seg.mean[breakpoints.shift]) < thisSd * sd.undo)]

    if (length(undo.breakpoints) < 1) {
      break
    }

    undo.breakpoints.shift <- undo.breakpoints + 1
    undo.df <- tempShort[undo.breakpoints, ]
    undo.df$seg.mean.diff <- abs(tempShort$seg.mean[undo.breakpoints] - tempShort$seg.mean[undo.breakpoints.shift])
    min.index <- which.min(undo.df$seg.mean.diff)
    leftIndex <- undo.df$segnum[min.index]
    rightIndex <- leftIndex + 1
    tempShort[leftIndex, "loc.end"] <- tempShort[rightIndex, "loc.end"]
    tempShort[leftIndex, "seg.end"] <- tempShort[rightIndex, "seg.end"]
    tempShort[leftIndex, "num.mark"] <- tempShort[leftIndex, "num.mark"] + tempShort[rightIndex, "num.mark"]
    tempShort[leftIndex, "seg.mean"] <- mean(log(ratioData$lowratio[tempShort[leftIndex, "seg.start"]:tempShort[rightIndex, "seg.end"]], base=2))
    tempShort <- tempShort[-rightIndex, ]
    tempShort$segnum <- seq(1:nrow(tempShort))
  }

  return(tempShort)
}


