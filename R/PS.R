PS<-function(pi,cov=0){
  a=strsplit(pi[,3],",")
  tag=a[[which.max(pi[,2])]]
  ll=lapply(a,function(a){ return(length(setdiff(a,tag))) })         
  g=pi[unlist(ll)<=cov,]
  if(!is.null(dim(g))){ 
    g=g[order(as.numeric(g[,2]),decreasing=T),]
  }
  pi=pi[unlist(ll)>cov,]
  if(ncol(as.matrix(pi))==1){ g=c(list(g),list(pi)) }else if(nrow(pi)==0){ return(list(g)) }else{ g=c(list(g),PS(pi)) }
  return(g)
}