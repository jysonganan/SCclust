## ----echo=FALSE,results='asis', message=FALSE, warning=FALSE, fig.width=10, fig.height=2.5----

library(knitr)
library(SCclust)
tab2<-example1.varbin.20k
knitr::kable(tab2[1:6,],caption ="bin_mat")

gc <- varbin.gc.20k
tab3<-gc_one(tab2, gc)

tab4<-cbs.segment_one(tab3,alpha = 0.05, nperm = 1000, undo.SD = 1.0, min.width = 5, method = "multiplier", genome = "hg", graphic = TRUE)

knitr::kable(tab3[1:6,],caption ="bin_mat_normalized")
knitr::kable(tab4[1:6,],caption ="bin_mat_segmented")



