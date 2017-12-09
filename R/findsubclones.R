#'Identify subclones in hierarchical tree.
#'
#'Based on hierarchical clustering, identify the hard/soft clones.
#'@param hc The hclust objects with clones identified.
#'@param pinmat The pinmat.
#'@param pins The pins.
#'@param min_node_size An integer. Default: 6. The minimum node size for a subclone.
#'@param sim_round The number of permutation simulations for subclone identification. Default: 500.
#'@param lm_max Numeric value. Default: 0.001. The threshold parameter for the linear fit to identify subclones.
#'@param hc_method Default: average
#'@param base_share An integer. Default: 3. A balance parameter for controlling minimal number of shared features in a subclone node.
#'@param fdr_thresh FDR criterion for subclone nodes. Default: -2.
#'@param share_min A feature is considered shared if present in share_min fraction of leaves in a node.Default: 0.90.
#'@param bymax Logical. If TRUE (Default), use maximal of mean FDR for the node to find subclones.
#'@param climb_from_size An integer. Default: 2.
#'@param climb_to_share An integer. Default: 3.
#'@param graphic Logical. If TRUE (Default), generate the hierarchical tree plot with hard/soft clones labeled.
#'@return A list of hclust objects for clones.
#'@export


find_subclone <- function(hc, pinmat, pins, min_node_size = 6, sim_round = 500, lm_max = 0.001, hc_method = "average",
                          base_share = 3, fdr_thresh = -2, share_min = 0.90, bymax = TRUE, climb_from_size = 2,
                          climb_to_share = 3, graphic = TRUE){
  clonedef <- "soft"

  bigclones<- unique(hc$softclones[clonedef,])[hc$nodesize[unique(hc$softclones[clonedef,])] >= min_node_size]

  sub_hc_clone <- list()

  if(length(bigclones) > 0){

    for (i in 1:length(bigclones)){
      clonemat <- pinmat[,dimnames(pinmat)[[2]]%in%hc$leaflabel_list[[bigclones[i]]], drop = F]
      clonepins <- pins[(rowSums(clonemat) > 0) & (rowSums(clonemat) < ncol(clonemat)),,drop = F]
      clonemat <- clonemat[(rowSums(clonemat) > 0)&(rowSums(clonemat) < ncol(clonemat)),,drop = F]

      n_share <- base_share + sum(rowSums(clonemat) > (share_min * ncol(clonemat)))
      cellnames <- dimnames(clonemat)[[2]]

      res <- simFisher_parallel(clonepins, clonemat, sim_round = sim_round)
      vtrue <- res$true_fisherPV
      msim <- res$sim_fisherPV

      res_1 <- fdr_fisherPV(vtrue, msim, cellnames, lm_max = lm_max, graphic = FALSE)

      mfdr <- res_1$mat_fdr
      mdist <- res_1$mat_dist

      subhc <- hclust_tree(clonemat, mfdr, mdist, hc_method = hc_method)
      subhc_clone <- find_clone(subhc, fdr_thresh = fdr_thresh, share_min = share_min, n_share = n_share, bymax = bymax,
                             climb_from_size = climb_from_size, climb_to_share = climb_to_share, graphic = FALSE)
      sub_hc_clone[[i]] <- subhc_clone
    }

  }

  return(sub_hc_clone)
}





