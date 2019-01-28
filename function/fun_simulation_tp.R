DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

DATE_proj = DATE_ERAI[DATE_ERAI$Y>=y1 & DATE_ERAI$Y<=y2,]
LENGTH_proj = dim(DATE_proj)[1]

P1 = array(NaN,c(LENGTH_proj,1)) # matrice de probas d'occurrence
shape = array(NaN,c(LENGTH_proj,1))
rate = array(NaN,c(LENGTH_proj,1))

DATE_OBS= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

for (seas in 1:4){
  # seas = 1
  cat(season[seas],'\n')
  
  DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  idxdates=which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  DATE_ERAI=DATE_ERAI[idxdates,]
  
  idx_proj_dates = which( (DATE_ERAI['Y']>=y1) & (DATE_ERAI['Y']<=y2) )    # 2005 - 2017
  
  load(paste0('~/Documents/LSCE/SWG/slp_SD_diagnosis/predictor/predictor_', SEAS[seas], '_1.RData'))
  
  PCS=PCS[idx_proj_dates,]
  NB_pred=dim(PCS)[2]
  
  
  # Nom du fichier sauvegardé contenant les estimations
  filename = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SWG_ERAI_ESD_tp_", season[seas], "_cross-val_", year.begin, "_", year.end, "_1.RData")
  load(filename)
  
  tempP1 = array(NaN,c(length(idx_proj_dates),length(fit_stations_ro)))
  tempshape = array(NaN,c(length(idx_proj_dates),length(fit_stations_ro)))
  temprate = array(NaN,c(length(idx_proj_dates),length(fit_stations_ro)))
  
  slot=array(NA,length(fit_stations_rf))
  
  for (i in 1:length(fit_stations_rf)){
    slot[i]=length(slotNames(fit_stations_rf[[i]]))
  }
  idx = which(slot>0)
  
  for (k in idx){
    
    ## occurrence
    if(statio_pix_ro[[k]]=='statio'){
      par_ro=fit_stations_ro[[k]]
      temp=par_ro[1]
      tempP1[,k] = replicate(length(idx_proj_dates),temp)
    }
    
    if(statio_pix_ro[[k]]!='statio'){
      par_ro = coef(fit_stations_ro[[k]],matrix=T)
      temp = par_ro[1] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_ro[2:(NB_pred+1)])
      tempP1[,k] = exp(temp)/(1+exp(temp))
    }
    
    
    #intensity
    par_rf=coef(fit_stations_rf[[k]],matrix=T)
    if(statio_pix_rf[[k]]=='statio'){
      temp2 =par_rf
      tempshape[,k] =replicate(length(idx_proj_dates),exp (temp2[2]))
      temprate[,k] = replicate(length(idx_proj_dates),exp (temp2[1]))
    }
    
    if(statio_pix_rf[[k]]!='statio'){
      tempshape[,k] = exp(par_rf[1,2] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_rf[2:(NB_pred+1),2]))
      temprate[,k] = exp(par_rf[1,1] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_rf[2:(NB_pred+1),1]))
    }
    
  } # end for k
  
  td = which(DATE_OBS$Y==y1)[1]-1
  
  P1[(idxdates[idx_proj_dates] - td),]    = tempP1    ; cat(range(tempP1),    "\n")
  shape[(idxdates[idx_proj_dates] - td),] = tempshape ; cat(range(tempshape), "\n")
  rate[(idxdates[idx_proj_dates] - td),]  = temprate  ; cat(range(temprate),  "\n")
  
} # end for seas
save(P1, shape, rate, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/precip_P1_shape_rate_", y1, "_", y2, "_", type, "_1.RData"))
