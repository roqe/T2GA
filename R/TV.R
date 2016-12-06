TV<-function(z,S,dgv=0.4){
  if(length(z)==length(which(z==0))){ return(0) }
  kik=which(z!=0)
  diag(S)=dgv
  return(t(z[kik])%*%compSS(S[kik,kik])%*%z[kik])
}