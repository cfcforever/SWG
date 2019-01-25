FAR <- function(R_nat, R_real, threshold){
  
  n = length(threshold)
  
  res = rep(NA, n)
  for (k in 1:n){
    p0 = sum(R_nat>=threshold[k]) / length(R_nat)
    p1 = sum(R_real>=threshold[k])/ length(R_real)
    
    res[k] = 1 - p0/p1
  }
    
  return (res)
  
}