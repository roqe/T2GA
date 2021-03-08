#' Compute T^2 and its p-value for each pathway
#'
#' This function computes the T^2 score and its significance level.
#'  
#' @param data Processed data using importdata function.
#' @param purb Perturbance threshold, default is 1.5 (after normalization).
#' @param pathDB Pathway database: "KEGG" or "Reactome".
#' @param ppi Protein-protein interaction database: "STRING_ppi" or "HitPredict_ppi".
#' @param intg Apply pathway integration or not, default is TRUE.
#' @param alpha Significance level, default is 0.05.
#' @param ncore Number of parallel computing cores, default is 1.
#' @keywords compute $T^2$
#' @export
#' @importFrom parallel mclapply
#' @importFrom corpcor make.positive.definite
#' @importFrom Matrix rankMatrix
#' @importFrom MASS ginv
#' @importFrom data.table rbindlist
#' @examples
#' dat1=importdata(TCR_5min)
#' res1=computeT2(dat1)
#' dat2=importdata(TCR_5min,TCR_15min)
#' res2=computeT2(dat2,pathDB="Reactome",ppi=HitPredict_v4)

computeT2=function(data,purb=1.5,pathDB="KEGG",ppi=STRING_v11,intg=TRUE,alpha=0.05,ncore=1){
  if(pathDB=="Reactome"){
    vex=Reactome_vex
    pid=Reactome_pid
  }else{
    vex=KEGG_vex
    pid=KEGG_pid
  }
  print("=================================================")
  print(paste0(" Using pathway database:  ",pathDB))
  print(paste0(" Using ppi database:      ",ifelse(nchar(ppi[1,1])>8,"STRING","HitPredict")))
  print("-------------------------------------------------")
  # data[-purb<data[,2]&data[,2]<purb,2]=0  
  data=data[-purb>data[,2]|data[,2]>purb,]
  ### data mapping
  n=which(vex[,2]%in%data[,1])
  print(paste("    #(mapped entries):   ",length(unique(vex[vex[,2]%in%data[,1],2]))))
  vexData=cbind(vex[n,],as.array(apply(vex[n,],1,function(v){ return(data[data[,1]%in%v[2],2]) })))
  colnames(vexData)[3]="exp"
  mp=table(vexData[,1])
  print(paste("    #(mapped pathways):  ",length(mp[mp!=0])))
  ### pathway integration
  pi=t(sapply(array(names(mp[mp!=0])),function(cp){
    pathway=as.matrix(vexData[which(vexData[,1]%in%cp),])
    if(ncol(pathway)==1){ pathway=matrix(pathway,ncol=3) }
    return(c(as.character(pathway[1,1]),nrow(pathway),paste0(pathway[,2],collapse=",")))
  }))
  ps=PS(pi)
  print(paste("    #(summary pathways): ",length(ps)))
  ### compute T2
  rplist=unlist(lapply(ps,function(ps){return(as.matrix(ps)[1,1])}))
  input=unlist(ifelse(intg,list(rplist),list(pi[,1])))
  r=do.call(rbind,mclapply(input,function(cp){
    pathway=as.matrix(vexData[which(vexData[,1]%in%cp),])
    if(ncol(pathway)==1){pathway=matrix(pathway,ncol=3)}
    return(TS(pathway,ppi,purb))
  },mc.cores=ncore))#
  ##### Result output ----------------------------------------
  rr=apply(r,1,function(r){
    rm=strsplit(r[2],",")[[1]]
    cm=matrix(pi[which(pi[,2]%in%r[3]),],ncol=3)
    sm=c()
    for(i in 1:nrow(cm)){
      if(length(setdiff(strsplit(cm[i,3],",")[[1]],rm))==0){ sm=rbind(sm,c(cm[i,1],r[-1])) }
    }
    ttl=pid[pid[,1]%in%sm[,1],]
    if(!is.null(nrow(ttl))){
      ttl=as.character(ttl[order(ttl[,1]),2])
    }else{
      ttl=as.character(ttl[2])
    }
    return(list(cbind(ttl,sm)))
  })
  rrr=do.call(rbind,unlist(rr,recursive=F))[,1:7]
  colnames(rrr)=c("Pathway title","Pathway ID","Uniprot IDs","#Mapped","df","T-square","p-value")
  rrr=as.data.table(unique(rrr[as.numeric(rrr[,7])<=alpha,]))
  print(paste("    #(enriched pathways):",nrow(rrr)))
  print("=================================================")
  return(rrr)
}