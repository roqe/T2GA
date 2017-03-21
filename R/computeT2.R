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
#' @param ncore Number of parallel computing cores, default is 7.
#' @param per Number of permutation to construct the null distribution.
#' @keywords compute $T^2$
#' @export
#' @importFrom parallel mclapply
#' @examples
#' dat1=importdata(tcr_05)
#' res1=computeT2(dat1)
#' dat2=importdata(fileName1=tcr_05,fileName2=tcr_15)
#' res2=computeT2(dat2,pathDB="Reactome",ppi=HitPredict_ppi)

computeT2=function(data,purb=1.5,pathDB="KEGG",ppi=STRING_ppi,intg=TRUE,alpha=0.05,ncore=3,per=10000){
  if(pathDB=="Reactome"){
    vex=Reactome_vex
    pid=Reactome_pid
  }else{
    vex=KEGG_vex
    pid=KEGG_pid
  }
  data[-purb<data[,2]&data[,2]<purb,2]=0  
  ### data mapping
  n=which(vex[,2]%in%data[,1])
  print(paste("    #(mapped entries):  ",length(unique(vex[vex[,2]%in%data[,1],2]))))
  vexData=cbind(vex[n,],as.array(apply(vex[n,],1,function(v){ return(data[data[,1]%in%v[2],2]) })))
  colnames(vexData)[3]="exp"
  mp=table(vexData[,1])
  print(paste("    #(mapped pathways): ",length(mp[mp!=0])))
  ### pathway integration
  pi=t(sapply(array(names(mp[mp!=0])),function(cp){
    pathway=as.matrix(vexData[which(vexData[,1]%in%cp),])
    if(ncol(pathway)==1){ pathway=matrix(pathway,ncol=3) }
    return(c(as.character(pathway[1,1]),nrow(pathway),paste0(pathway[,2],collapse=",")))
  }))
  ps=PS(pi)
  print(paste("    #(summary pathways):",length(ps)))
  ### compute T2
  rplist=unlist(lapply(ps,function(ps){return(as.matrix(ps)[1,1])}))
  if(intg==TRUE){
    input=rplist
  }else{ input=pi[,1] }
  r=do.call(rbind,mclapply(input,function(cp){
    pathway=as.matrix(vexData[which(vexData[,1]%in%cp),])
    if(ncol(pathway)==1){pathway=matrix(pathway,ncol=3)}
    return(TS(pathway,ppi,per,purb))
  },mc.cores=ncore,mc.preschedule=FALSE))#
  ##### Result output ----------------------------------------
  rr=apply(r,1,function(r){
    rm=strsplit(r[3],",")[[1]]
    cm=matrix(pi[which(pi[,2]%in%r[2]),],ncol=3)
    sm=c()
    for(i in 1:nrow(cm)){
      if(length(setdiff(strsplit(cm[i,3],",")[[1]],rm))==0){ sm=rbind(sm,c(cm[i,1],r[2:8])) }
    }
    ttl=pid[pid[,1]%in%sm[,1],]
    if(!is.null(nrow(ttl))){
      ttl=as.character(ttl[order(ttl[,1]),2])
    }else{
      ttl=as.character(ttl[2])
    }
    return(list(cbind(ttl,sm)))
  })
  rrr=do.call(rbind,unlist(rr,recursive=F))[,1:6]
  colnames(rrr)=c("Pathway title","Pathway ID","#Mapped","Uniprot IDs","T-square","p-value")
  rrr=rrr[as.numeric(rrr[,6])<=alpha,]
  print(paste("    #(enriched pathways): ",nrow(rrr)))
  return(rrr)
}