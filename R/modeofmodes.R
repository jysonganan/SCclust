#v and otherlabel are integer vectors of equal length. Find the mode of v for
#each value of otherlabel and find the mode of modes. If the latter 
#computation results in ties, find the values nearest tiebreaker (integer)
#and of those choose the highest if tiebreakerside (character) is "greater" or
#the lowest.
modeofmodes<-function(v,otherlabel,tiebreaker,tiebreakerside){
  countmat<-as.matrix(tapply(X=v,INDEX=list(as.factor(v),as.factor(otherlabel)),
                             FUN=length))
  modes<-dimnames(countmat)[[1]][apply(countmat,2,which.max)]
  modecounts<-tapply(X=modes,INDEX=as.factor(modes),FUN=length)
  topmodes<-as.numeric(names(modecounts))[modecounts==max(modecounts)]
  untie<-topmodes[abs(topmodes-tiebreaker)==min(abs(topmodes-tiebreaker))]
  if(tiebreakerside=="greater")return(max(untie))
  return(min(untie))
}


getmode<-function(x){as.numeric(names(which.max(tapply(X=x, INDEX=as.factor(x), FUN=length))))}