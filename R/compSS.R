#' @importFrom MASS ginv

compSS<-function(S){
  SS=try(solve(S),silent=T)
  if(!is.matrix(SS)) {
    SS=ginv(S)
  }
  return(SS)
}