#' @importFrom corpcor make.positive.definite
#' @importFrom MASS ginv

compSS<-function(S){
  SS=try(solve(S),silent=T)
  if(!is.matrix(SS)) {
    SS=try(ginv(S),silent=T)
  }
  return(make.positive.definite(SS))
}