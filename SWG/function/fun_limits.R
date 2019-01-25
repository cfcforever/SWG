fun_limits <- function(data, min0 = F){
  if (min0==T){
    return(c(0, max(abs(data))))
  }else{
    return(c(-max(abs(data)), max(abs(data))))
  }
}