##weighting
pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond,NLON)

W_HGT = array(NaN,dim=dim(HGT.Atl))

for(pix in 1:940){
  W_HGT[,pix] = HGT.Atl[,pix]/scale[pix]
}

pcaW=prcomp(W_HGT)


##############################################################################
search.PCAlevel <- function(PCAGL,p){
  p.var=((PCAGL$sdev)*(PCAGL$sdev)/sum((PCAGL$sdev)*(PCAGL$sdev)))
  s=0
  i=1
  t=0
  while( (i<=(length(PCAGL$sdev))) && (t==0) ){
    s=s+p.var[i]
    #cat((s*100))
    if(p>=(s*100)){
      i=i+1
    }
    else{
      t=1
    }
  }
  
  if(t==0){
    stop("ERROR IN search.PCAlevel\n")
  }
  
  return(i)
}
