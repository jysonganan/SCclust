TreePy<-function(data, method, metric="euclidean"){
	#source("/Volumes/Wigler-1/user/gsun/hclust/TBEST/R/TreeStat.R")
	if(data.class(data)=="dist"){
		Clustersize<-TBEST:::TreeStat(data,mystat="fldc",method,metric)[,4]
		hc<-hclust(data,method)
	}
	else{
		Clustersize<-TBEST:::TreeStat(data.matrix(data),mystat="fldc",method,metric)[,4]
		hc<-hclust(dist(data.matrix(data),method=metric),method)
	}
	Py<-cbind(hc$merge,hc$height,Clustersize)
	d<-Py[,1:2]
	d[Py[,1:2]<0]<- -d[Py[,1:2]<0]-1
	d[Py[,1:2]>0]<- d[Py[,1:2]>0]+nrow(Py)
	Py[,1:2]<-d
	dimnames(Py)[[2]]<-c("index1","index2","height","clustersize")
	return(Py)
}
