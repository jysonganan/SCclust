
#'Generate files for the downstream visualization using viewer.
#'@param output_file_dir The directory made to save the output files for the viewer.
#'@param seg.quantal The seg.quantal profiles.
#'@param ratio.quantal The ratio.quantal profiles
#'@param pins The pins.
#'@param pinmat The pinmat.
#'@param mat_dist The dissmilarity based on Fisher's test p-values for hierarchical clustering.
#'@param hc_clone The hc object with clones identified.
#'@param sub_hc_clone The list of hc objects in clones with subclones identified.
#'@param subcloneTooBig Default: 0.8
#'@param smear The smear parameter. Default value: 2
#'@param study The name of your study.
#'@export


output_viewer <- function(output_file_dir, seg.quantal, ratio.quantal, pins, pinmat, mat_dist, hc_clone, sub_hc_clone,
                          subcloneTooBig = 0.8, smear = 2, study){


  # system(paste("mkdir ", paste("../", output_file_dir, sep = ""), sep = ""))

  pytableP <- TreePy(data = as.dist(mat_dist), method = "average")
  pytableP <- cbind(pytableP, hc_clone$maxfdr)
  dimnames(pytableP)[[2]][ncol(pytableP)] <- "log10fdr"


  clonetype <- "soft"
  clonetable <- data.frame(hc_clone$labels, rep(0, length(hc_clone$labels)),
                           rep(0, length(hc_clone$labels)), stringsAsFactors = F)
  dimnames(clonetable)[[2]] <- c("ID", "clone", "subclone")


  ##clones
  for (nodes in unique(hc_clone$softclones[clonetype,])){
    clonetable[hc_clone$leafID_list[[nodes]], "clone"] <- nodes
  }

  ##subclones
  if (length(sub_hc_clone) > 0){

    for (i in 1:length(sub_hc_clone)){
      subhc <- sub_hc_clone[[i]]

      clunique <- unique(subhc$softclones[clonetype,])
      if (length(clunique) > 1){
        for(nodes in clunique){
          clonetable[match(subhc$leaflabel_list[[nodes]], clonetable[,1]), "subclone"] <- nodes
        }
      }

      if(length(clunique) == 1){
        if(subhc$nodesize[clunique] < (subcloneTooBig * max(subhc$nodesize))){
          clonetable[match(subhc$leaflabel_list[[clunique]], clonetable[,1]), "subclone"] <- clunique}
      }
    }
  }



  if (!is.null(output_file_dir)){
    write.table(seg.quantal, paste(output_file_dir, "/uber.hg19.", study, ".seg.quantal.R.seg.txt", sep = ""),
                col.names = T, row.names = F, sep = "\t", quote = F)
    write.csv(ratio.quantal,paste(output_file_dir, "/uber.hg19.", study, ".lowratio.quantal.R.ratio.csv", sep = ""))
    write.table(pins, paste(output_file_dir, "/", study, "smear", smear, "bpPins.pins.txt", sep = ""),
                col.names = T, row.names = F, sep = "\t", quote = F)
    write.table(pinmat, paste(output_file_dir, "/", study, "smear", smear, "bpPinMat.pinmat.txt", sep = ""),
                col.names = T, row.names = F, sep = "\t", quote = F)
    write.table(pytableP, paste(output_file_dir, "/", study, "smear", smear, "bpFisherTreePyP.tree.txt", sep = ""),
                col.names = T, row.names = F, sep = "\t", quote = F)
    write.table(clonetable, paste(output_file_dir, "/", study, "smear", smear, "bpFisherPcloneTracks.clone.txt", sep = ""),
                col.names = T,row.names = F, sep = "\t", quote = F)
  }

}
