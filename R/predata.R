#' @import plyr

predata=function(data,outth=100){
  print(paste("    #(input site/probe):",nrow(data)))
  ### remove missing
  rm=data[,2]=="NAN"|data[,2]=="NaN"|data[,2]=="NA"|data[,2]=="na"|data[,2]=="-"|data[,1]==""|data[,2]==""|is.na(data[,2])|is.na(data[,1])
  data=data[!rm,]
  ### multiple ids one value
  data1=data[nchar(data[,1])<11,]
  data2=data[nchar(data[,1])>10,]
  if(nrow(data2)!=0){
    data2=do.call(rbind,apply(data2,1,function(v){
      a=unlist(strsplit(as.character(v[1]),"|"))
      return(cbind(a,v[2]))
    }))
  }
  colnames(data2)=colnames(data1)
  data=rbind(data1,data2)
  colnames(data)=c("id","exp")
  rownames(data)=c()
  data[,"id"]=substr(data[,"id"],1,6)
  data=as.data.frame(data,stringsAsFactors=F)
  ### one id mutiple values
  class(data[,2])="numeric"
  data=ddply(data,.(id),summarise,M=median(exp))
  ### normalization
  if(min(data[,2])>=0){
    data[data[,2]==0,2]=min(data[data[,2]!=0,2])
    data[,2]=log2(data[,2])
  }
  ### outliers replacement
  out=data[,2]>=outth|data[,2]<=-outth
  data[out,2]=data[out,2]/abs(data[out,2])*max(abs(data[,2]))
  ### standardization
  data[,2]=data[,2]/sd(data[,2])
  return(data)
}
