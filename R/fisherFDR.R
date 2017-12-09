#'Compute FDRs for Fisher's test p-values.
#'
#'Linear fit to the tail of empirical null distribution of Fisher p-values;
#'FDR computation: compare true to simulated CDF(empirical null).
#'@param true_fisherPV The Fisher's test p-values for the observation.
#'@param sim_fisherPV The Fisher's test p-values for the permutations.
#'@param cell_names A character vector. The names of cells.
#'@param lm_max Numeric value. Default: 0.001. The threshold parameter for the linear fit.
#'@param graphic Logical. If TRUE (Default), generate the CDF plots of p-values and FDR values.
#'@return A list containing the matrix of the FDR values (mat_fdr) and distance matrix based on Fisher's test p-values (mat_dist).
#'@export


fdr_fisherPV <- function(true_fisherPV, sim_fisherPV, cell_names, lm_max = 0.001, graphic = TRUE){
  #cell_names = NULL
  ## Sort the true and the null data and get counts for each unique fisher pv value
  sim_sort <- sort(sim_fisherPV)
  sim_unique <- unique(sim_sort)
  sim_count <- tapply(match(sim_sort, sim_unique), match(sim_sort, sim_unique), length)
  true_sort <- sort(true_fisherPV)
  true_unique <- unique(true_sort)
  true_count <- tapply(match(true_sort, true_unique), match(true_sort, true_unique), length)


  ## linear fit to the tail of empirical null distribution of Fisher p-values
  ##############################################################################

  # lm_max -- A parameter in a linear fit
  #Observed p-values are often far lower than any p-value sampled from the null; to determine FDR in
  #such cases get a power-law fit to the low-p tail of the null CDF and use it to extrapolate to very
  #low p-values. Use the actual null CDF to estimate FDR for higher p-values.
  lmfit <- lm(log(cumsum(sim_count[(cumsum(sim_count)/sum(sim_count)) < lm_max])/sum(sim_count))~
             log(sim_unique[(cumsum(sim_count)/sum(sim_count)) < lm_max]))

  ## generate the plot CDF vs FisherPV
  if(graphic == TRUE){
    curve(exp(lmfit$coefficients[1] + lmfit$coefficients[2]*log(x)),
          from = min(true_unique), to = max(true_unique), log = "xy")
    points(sim_unique, cumsum(sim_count)/sum(sim_count))
    points(true_unique, cumsum(true_count)/sum(true_count), col = "red")
  }


  ## FDR computation: compare true to simulated CDF(empirical null)
  ##############################################################################
  dat_TrueSim <- cbind(c(log(true_unique), log(sim_unique)),
                       c(log(cumsum(true_count)/sum(true_count)), log(cumsum(sim_count)/sum(sim_count))),
                       c(rep(0, length(true_unique)), rep(1, length(sim_unique))))
  dat_TrueSim <- dat_TrueSim[order(dat_TrueSim[,1]),]

  ## (1) estimate FDR for higher p-values
  #(linear interpolation)
  num_sim <- cumsum(dat_TrueSim[,3])[dat_TrueSim[,3] == 0]  #num of sim with less pv than the specific true value
  x1pos <- match(num_sim, cumsum(dat_TrueSim[,3]))[num_sim > 0]
  x2pos <- match(num_sim + 1, cumsum(dat_TrueSim[,3]))[num_sim > 0]
  logpv1 <- dat_TrueSim[x1pos,1]
  logpv2 <- dat_TrueSim[x2pos,1]
  logcdf1 <- dat_TrueSim[x1pos,2]
  logcdf2 <- dat_TrueSim[x2pos,2]

  logfdr_interp <- rep(0, length(true_unique))
  logfdr_interp[num_sim > 0] <- (logcdf2 - logcdf1)*log(true_unique)[num_sim > 0]/(logpv2 - logpv1) +
    (logcdf1*logpv2 - logcdf2*logpv1)/(logpv2 - logpv1) - log(cumsum(true_count)/sum(true_count))[num_sim > 0]
  logfdr_interp[true_unique > max(sim_unique)] <- 0

  ## (2) estimate FDR for low p-values
  #(a power-law fit to the low-p tail of the null CDF, use it to extrapolate to very low p-values)
  logfdr <- logfdr_interp
  lmu <- max(sim_unique[(cumsum(sim_count)/sum(sim_count)) < lm_max])

  logfdr_lmfit <-
    lmfit$coefficients[2]*log(true_unique) + lmfit$coefficients[1] - log(cumsum(true_count)/sum(true_count))

  if(is.finite(lmu) & min(true_unique) < min(sim_unique)){

    logfdr[true_unique < lmu & true_unique > min(sim_unique)]<-
      (logfdr_lmfit[true_unique < lmu & true_unique > min(sim_unique)]*
         (log(lmu) - log(true_unique[true_unique < lmu & true_unique > min(sim_unique)])) -
         logfdr_interp[true_unique < lmu & true_unique > min(sim_unique)]*(log(min(sim_unique)) -
            log(true_unique[true_unique < lmu & true_unique > min(sim_unique)])))/(log(lmu) - log(min(sim_unique)))

    logfdr[true_unique < min(sim_unique)] <- logfdr_lmfit[true_unique < min(sim_unique)]
  }

  logfdr <- cummax(logfdr)
  logfdr[logfdr > 0] <- 0
  if(graphic){
    plot(true_unique, exp(logfdr), log = "xy")}

  logfdr_all <- logfdr[match(true_fisherPV, true_unique)]

  # mat_fisherPV <- matrix(ncol = (1 + sqrt(1 + 8*length(true_fisherPV)))/2,
  #                 nrow = (1 + sqrt(1 + 8*length(true_fisherPV)))/2, data = 0)
  # mat_fisherPV[upper.tri(mat_fisherPV)] <- log10(true_fisherPV)
  # mat_fisherPV <- pmin(mat_fisherPV, t(mat_fisherPV))
  #dimnames(mat_fisherPV) <- list(cellnames, cellnames)


  mat_fdr <- matrix(ncol = (1 + sqrt(1 + 8*length(true_fisherPV)))/2,
                 nrow = (1 + sqrt(1 + 8*length(true_fisherPV)))/2,data = 0)
  mat_fdr[upper.tri(mat_fdr)] <- logfdr_all/log(10)     ##############################################??log(10)
  mat_fdr <- pmin(mat_fdr,t(mat_fdr))
  colnames(mat_fdr) = rownames(mat_fdr) <- cell_names

  mat_dist <- matrix(ncol = (1 + sqrt(1 + 8*length(true_fisherPV)))/2,
                     nrow = (1 + sqrt(1 + 8*length(true_fisherPV)))/2,data = 0)
  mat_dist[upper.tri(mat_dist)] <- log10(true_fisherPV)
  mat_dist <- pmin(mat_dist,t(mat_dist))
  colnames(mat_dist) = rownames(mat_dist) <- cell_names


  # if(is.null(cell_names) == TRUE){
  #   colnames(mat_fdr) = rownames(mat_fdr) <- 1:dim(mat_fdr)[[1]]
  # }else{
  #   colnames(mat_fdr) = rownames(mat_fdr) <- cell_names
  # }
  return(list(mat_fdr = mat_fdr, mat_dist = mat_dist))
}
