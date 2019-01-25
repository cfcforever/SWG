nPC = 2
#### case_2: only SysDyn from 1979 to 1997 ####
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

# Estimation of parameters in 4 seaons
for (seas in 1:4){
  # seas = 1
  cat(season[seas],'\n')
  
  # DATE_ERAI as the same as DATE_OBS
  DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  
  # choose for the right months for season 
  idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  
  # DATE_ERAI only in season
  DATE_ERAI = DATE_ERAI[idxdates,]
  
  # /!\ load PCS of season - VERY IMPORTANT
  load(paste0("~/Documents/LSCE/SWG/t2m/case_2/geop_",SEAS[seas],"_SysDyn.RData"))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  
  # load tmean - this is main data
  load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
  dat = tmean
  
  # number of station - here only 1
  NB_STATION = 1
  
  idx_control_dates = which((DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']<=1997))
  
  # PCS in estimation years
  PCS=PCS[idx_control_dates,]
  pc_name=names(PCS)
  
  # cherche les indices de jour dans DATE_ECA qui correspondent aux dates ERAI pour la période de calibration
  idx_dates=array(0,dim=length(idx_control_dates))
  for (i in 1:length(idx_control_dates)){
    idx_dates[i]=which(((DATE_OBS['Y'])[,1]==(DATE_ERAI['Y'])[idx_control_dates[i],1])
                       & ((DATE_OBS['m'])[,1]==(DATE_ERAI['m'])[idx_control_dates[i],1])
                       & ((DATE_OBS['d'])[,1]==(DATE_ERAI['d'])[idx_control_dates[i],1]))
  } # end for i
  
  fmlatt=as.formula(paste("TT ~ ",paste(pc_name, collapse= "+")))
  fit_stations_tt <- vector("list", NB_STATION)
  
  for (i in 1:NB_STATION){
    
    temp=dat[idx_dates,i]
    
    if(length(which(is.na(temp)!=TRUE))==0){
      fit_stations_tt[[i]]=NaN
    }
    if(length(which(is.na(temp)!=TRUE))!=0){
      tt = data.frame(TT=temp,PCS)   # ncol = 1+20
      fit_stations_tt[[(i)]] = try(vgam(fmlatt ,uninormal(zero=NULL), data=tt,maxit=1000, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=T)
    }                          
    
  } # end for i NB_STATION
  
  output = paste('~/Documents/LSCE/SWG/t2m/case_2/SWG_ERAI_ESD_tmean_',season[seas],'_cross-val_1979_1997.RData',sep="")
  cat(output,'\n')
  save(fit_stations_tt,file = output)
} # end for seas


#### case_2: simulation of mean and sd from 1998 to 2017 for counter-factual ####
# period for projection (1998 - 2017)
DATE_proj = atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day'))
LENGTH_proj = dim(atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day')))[1]

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
  
  DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  DATE_ERAI=DATE_ERAI[idxdates,]
  
  idx_proj_dates = which( (DATE_ERAI['Y']>=1998) & (DATE_ERAI['Y']<=2017) ) 
  
  load(paste("~/Documents/LSCE/SWG/t2m/case_2/geop_",SEAS[seas],"_SysDyn.RData",sep=""))
  
  PCS=PCS[idx_proj_dates,] # Pour T
  NB_pred=dim(PCS)[2]
  
  
  # Nom du fichier sauvegardé contenant les estimations
  filename=paste('~/Documents/LSCE/SWG/t2m/case_2/SWG_ERAI_ESD_tmean_',season[seas],'_cross-val_1979_1997.RData',sep="")
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
  
  td = which(DATE_OBS$Y==1998)[1]-1
  
  mean[(idxdates[idx_proj_dates] - td),] = tempmean   ; cat(range(tempmean), "\n")
  sd[(idxdates[idx_proj_dates] - td),] = tempsd       ; cat(range(tempsd)  , "\n")
  
} # end for seas

save(sd, mean, file = "~/Documents/LSCE/SWG/t2m/case_2/tmean_mean_sd_1998_2017_nat.RData")

#### case_2: estimation of model-2 from 1998 to 2017 for factual ####
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

# Estimation of parameters in 4 seaons
for (seas in 1:4){
  # seas = 1
  cat(season[seas],'\n')
  
  # DATE_ERAI as the same as DATE_OBS
  DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  
  # choose for the right months for season 
  idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  
  # DATE_ERAI only in season
  DATE_ERAI = DATE_ERAI[idxdates,]
  
  # /!\ load PCS of season - VERY IMPORTANT
  load(paste0("~/Documents/LSCE/SWG/t2m/case_2/geop_",SEAS[seas],"_SysDyn.RData"))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  
  # load tmean - this is main data
  load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
  dat = tmean
  
  # number of station - here only 1
  NB_STATION = 1
  
  idx_control_dates = which((DATE_ERAI['Y']>=1998) & (DATE_ERAI['Y']<=2017))
  
  # PCS in estimation years
  PCS=PCS[idx_control_dates,]
  pc_name=names(PCS)
  
  # cherche les indices de jour dans DATE_ECA qui correspondent aux dates ERAI pour la période de calibration
  idx_dates=array(0,dim=length(idx_control_dates))
  for (i in 1:length(idx_control_dates)){
    idx_dates[i]=which(((DATE_OBS['Y'])[,1]==(DATE_ERAI['Y'])[idx_control_dates[i],1])
                       & ((DATE_OBS['m'])[,1]==(DATE_ERAI['m'])[idx_control_dates[i],1])
                       & ((DATE_OBS['d'])[,1]==(DATE_ERAI['d'])[idx_control_dates[i],1]))
  } # end for i
  
  fmlatt=as.formula(paste("TT ~ ",paste(pc_name, collapse= "+")))
  fit_stations_tt <- vector("list", NB_STATION)
  
  for (i in 1:NB_STATION){
    
    temp=dat[idx_dates,i]
    
    if(length(which(is.na(temp)!=TRUE))==0){
      fit_stations_tt[[i]]=NaN
    }
    if(length(which(is.na(temp)!=TRUE))!=0){
      tt = data.frame(TT=temp,PCS)   # ncol = 1+20
      fit_stations_tt[[(i)]] = try(vgam(fmlatt ,uninormal(zero=NULL), data=tt,maxit=1000, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=T)
    }                          
    
  } # end for i NB_STATION
  
  output = paste('~/Documents/LSCE/SWG/t2m/case_2/SWG_ERAI_ESD_tmean_',season[seas],'_cross-val_1998_2017.RData',sep="")
  cat(output,'\n')
  save(fit_stations_tt,file = output)
  
} # end for seas


#### case_2: estimation of mean and sd from 1998 to 2017 for factual ####
# period for projection (1998 - 2017)
DATE_proj = atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day'))
LENGTH_proj = dim(atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day')))[1]

# for temperature, we use gaussian distribution, so mean and sd involve
mean = array(NaN,c(LENGTH_proj,1)) # matrice de probas d'occurrence
sd = array(NaN,c(LENGTH_proj,1))

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

for (seas in 1:4){
  # seas = 1
  cat(season[seas],'\n')
  
  DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  
  # choose for the right months for season 
  idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  
  # DATE_ERAI only in season
  DATE_ERAI = DATE_ERAI[idxdates,]
  
  # number of station - here only 1
  NB_STATION = 1
  
  idx_proj_dates = which((DATE_ERAI['Y']>=1998) & (DATE_ERAI['Y']<=2017))
  
  # /!\ load PCS of season - VERY IMPORTANT
  load(paste0("~/Documents/LSCE/SWG/t2m/case_2/geop_",SEAS[seas],"_SysDyn.RData"))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS=PCS[idx_proj_dates,] # Pour T
  NB_pred=dim(PCS)[2]
  
  filename=paste('~/Documents/LSCE/SWG/t2m/case_2/SWG_ERAI_ESD_tmean_',season[seas],'_cross-val_1998_2017.RData',sep="")
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
  
  td = which(DATE_OBS$Y==1998)[1]-1
  
  mean[(idxdates[idx_proj_dates] - td),] = tempmean   ; cat(range(tempmean), "\n")
  sd[(idxdates[idx_proj_dates] - td),] = tempsd       ; cat(range(tempsd)  , "\n")
  
}
save(mean, sd, file = "~/Documents/LSCE/SWG/t2m/case_2/tmean_mean_sd_1998_2017_rea.RData")
