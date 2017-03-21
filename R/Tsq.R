##### T-square
Tsq=function(pathway,ppi,per,purb){
  pathway=matrix(pathway[order(pathway[,2]),],ncol=3)
  z=matrix(as.numeric(pathway[,3]))
  m=which(esb_ID[,2]%in%pathway[,2])
  S=diag(x=1,ncol=length(z),nrow=length(z))
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
  T2=TV(z,S)
  I=diag(x=dgv,ncol=length(z),nrow=length(z))
  T2I=TV(z,I)
  H0=sapply(1:per,function(i){
    tz=rnorm(length(z))
    tz[which(-purb<tz&tz<purb)]=0
    return(c(TV(tz,S),TV(tz,I)))
  })
  pv=(length(which(H0[1,]>=as.numeric(T2)))+1)/per
  pvI=(length(which(H0[2,]>=as.numeric(T2I)))+1)/per
  return(c(as.character(pathway[1,1]),nrow(pathway),paste0(pathway[,2],collapse=","),T2,pv,T2I,pvI,pchisq(T2,length(z),lower.tail=FALSE)))
}