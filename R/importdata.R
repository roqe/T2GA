#' Data import and pre-process
#'
#' This function removes NA, replaces extreme values, and standardizes the data. It also map the data with uniprot identifiers to ensp idnetifiers.
#' @param fileName1 Expression data with Uniprot identifiers.
#' @param fileName2 Expression data with Uniprot identifiers, optional for time-course data
#' @param outth Outlier threshold, default is 10.
#' @param type data type, options: "exp" for expression values, "pv" for p-values, default is "exp"
#' @keywords pre-process
#' @export
#' @examples
#' dat1=importdata(TCR_5min)
#' dat2=importdata(TCR_5min,TCR_15min)

importdata=function(fileName1,fileName2=data.frame(),outth=10,typ="exp"){
  print("=================================================")
  print(" Dataset summary:")
  print("-------------------------------------------------")
  if(nrow(fileName2)==0){
    data<-as.data.frame(predata(as.matrix(fileName1),outth,typ),stringsAsFactors=F)
  }else{
    data1<-predata(as.matrix(fileName1),outth,typ)
    data2<-predata(as.matrix(fileName2),outth,typ)
    #time series, divide time 1 and time 2
    data12<-setdiff(data1$id,data2$id)
    data12<-cbind(data12,data1[data1$id%in%data12,2]-data2$M[which.min(abs(data2$M))])
    data21<-setdiff(data2$id,data1$id)
    data21<-cbind(data21,data2[data2$id%in%data21,2]-data1$M[which.min(abs(data1$M))])
    data<-intersect(data1$id,data2$id)
    data<-cbind(data,data1[complete.cases(data1[match(data1$id,data),]),2]-data2[complete.cases(data2[match(data2$id,data),]),2])
    data<-as.data.frame(rbind(data,data12,data21),stringsAsFactors=F)#union set
  }
  colnames(data)<-c("id","exp")
  class(data$exp)<-"numeric"
  print(paste("    #(input proteins):   ",nrow(data)))
  print("=================================================")
  return(data)
}