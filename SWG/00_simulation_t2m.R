fun_simulation <- function(predictor, NUM, parameter, output, type, 
                           y1 = 2005, y2 = 2017, 
                           year.begin = 1979, year.end = 2004){
  
  #### case_1: simulation of mean and sd from 1998 to 2017 for counter-factual ####
  # period for projection (1998 - 2017)
  DATE_proj = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))
  LENGTH_proj = dim(atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day')))[1]
  
  # for temperature, we use gaussian distribution, so mean and sd involve
  mean = array(NaN,c(LENGTH_proj,1)) # matrice de probas d'occurrence
  sd = array(NaN,c(LENGTH_proj,1))
  
  DATE_OBS= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  
  MON=c(12,1:11)
  MON=matrix(MON,nrow=4,ncol=3,byrow=T)
  SEAS=c('DJF','MAM','JJA','SON')
  season=c('Winter','Spring','Summer','Fall')
  
  for (seas in 1:4){
    # seas = 2
    cat(season[seas],'\n')
    
    DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))
    idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
    DATE_ERAI=DATE_ERAI[idxdates,]
    
    idx_proj_dates = which( (DATE_ERAI['Y']>=y1) & (DATE_ERAI['Y']<=y2) ) 
    
    load(paste0(predictor, SEAS[seas], "_", NUM, ".RData"))
    PCS=PCS[idx_proj_dates,] # Pour T
    NB_pred=dim(PCS)[2]
    
    
    # Nom du fichier sauvegardé contenant les estimations
    filename = paste0(parameter, season[seas], "_cross-val_", year.begin, "_", year.end, "_", NUM, ".RData")
    load(filename)
    
    tempmean = array(NaN,c(length(idx_proj_dates),length(fit_stations_tt)))
    tempsd = array(NaN,c(length(idx_proj_dates),length(fit_stations_tt)))
    
    slot=array(NA,length(fit_stations_tt))
    for (i in 1:length(fit_stations_tt)){
      slot[i]=length(slotNames(fit_stations_tt[[i]]))
    }
    
    idx=which(slot>0)
    
    for (k in idx){
      #cat(k,'\n')
      par_tt=coef(fit_stations_tt[[k]],matrix=T)
      tempsd[,k] = exp(par_tt[1,2] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_tt[2:(NB_pred+1),2]))
      tempmean[,k] = par_tt[1,1] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_tt[2:(NB_pred+1),1])
    } # end for k
    
    td = which(DATE_OBS$Y==y1)[1]-1
    
    mean[(idxdates[idx_proj_dates] - td),] = tempmean   ; cat(range(tempmean), "\n")
    sd[(idxdates[idx_proj_dates] - td),] = tempsd       ; cat(range(tempsd)  , "\n")
    
  } # end for seas
  
  save(sd, mean, file = paste0(output, "t2m_mean_sd_", y1, "_", y2, "_", type, "_", NUM, ".RData"))
}

for (k in 1:63){
  fun_simulation(predictor = '~/Documents/LSCE/SWG/NA_1.5/pca/predictor_',
                 parameter = '~/Documents/LSCE/SWG/NA_1.5/tmean/1979_1998/SWG_ERAI_ESD_tmean_',
                 NUM = k, 
                 type = "nat",
                 output = '~/Documents/LSCE/SWG/NA_1.5/tmean/nat/', 
                 y1 = 1999, y2 = 2017,
                 year.begin = 1979, year.end = 1998)
}

for (k in 1:63){
  fun_simulation(predictor = '~/Documents/LSCE/SWG/NA_1.5/pca/predictor_',
                 parameter = '~/Documents/LSCE/SWG/NA_1.5/tmean/1999_2017/SWG_ERAI_ESD_tmean_',
                 NUM = k, 
                 type = "rea",
                 output = '~/Documents/LSCE/SWG/NA_1.5/tmean/rea/', 
                 y1 = 1999, y2 = 2017,
                 year.begin = 1999, year.end = 2017)
}
