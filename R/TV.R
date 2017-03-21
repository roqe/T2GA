TV<-function(z,S){
  if(length(z)==length(which(z==0))){ return(0) }
  kik=which(z!=0)
  return(t(z[kik])%*%ginv(S[kik,kik])%*%z[kik])
}