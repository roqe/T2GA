##### T-square
TS=function(pathway,ppi,purb,dgv=0.4){
  pathway=matrix(pathway[order(pathway[,2]),],ncol=3)
  z=matrix(as.numeric(pathway[,3]))
  m=which(esb_ID[,2]%in%pathway[,2])
  S=diag(x=dgv,ncol=length(z),nrow=length(z))
  if(dim(S)[1]!=1){
    for(i in 2:nrow(S)){
      for(j in 1:(i-1)){
        x1=pathway[,2][i]
        x2=pathway[,2][j]
        if(names(ppi)[1]=="Uniprot1"){
          s1=x1
          s2=x2
        }else{
          s1=esb_ID[m,1][which(esb_ID[m,2]%in%x1)]
          s2=esb_ID[m,1][which(esb_ID[m,2]%in%x2)]
        }
        if(length(s1)!=0 & length(s2)!=0){
          p=ppi[which(ppi[,1]%in%s1),]
          p=p[which(p[,2]%in%s2),3]
          if(length(p)>0){
            if(z[which(pathway[,2]%in%x1)]*z[which(pathway[,2]%in%x2)]<0){
              S[i,j]=-mean(p)
              S[j,i]=-mean(p)
            }else{
              S[i,j]=mean(p)
              S[j,i]=mean(p)
            }
          } 
        }
      }
    }
  }
  S=make.positive.definite(S)
  r=rankMatrix(S,tol=1e-10)[1]
  T2=TV(z,S)
  I=diag(x=dgv,ncol=length(z),nrow=length(z))
  T2I=TV(z,I)
  return(c(as.character(pathway[1,1]),
           paste0(pathway[,2],collapse=","),nrow(pathway),r,
           T2,pchisq(T2,r,lower.tail = F),
           T2I,pchisq(T2I,r,lower.tail = F)))
}