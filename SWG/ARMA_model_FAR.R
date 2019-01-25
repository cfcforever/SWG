library(ggplot2)
library(forecast)
library(tseries)

#### load t2m data - obs ####
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData") 
range(tmean)

#### estimation 1979-1997 ####
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

for (seas in 1:4){
  DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  
  # choose for the right months for season 
  idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  
  # DATE_ERAI only in season
  DATE_ERAI = DATE_ERAI[idxdates,]
  
  dat = tmean
  NB_STATION = 1
  
  idx_control_dates = which((DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']<=1997))
  
  # cherche les indices de jour dans DATE_ECA qui correspondent aux dates ERAI pour la période de calibration
  idx_dates=array(0,dim=length(idx_control_dates))
  for (i in 1:length(idx_control_dates)){
    idx_dates[i]=which(((DATE_OBS['Y'])[,1]==(DATE_ERAI['Y'])[idx_control_dates[i],1])
                       & ((DATE_OBS['m'])[,1]==(DATE_ERAI['m'])[idx_control_dates[i],1])
                       & ((DATE_OBS['d'])[,1]==(DATE_ERAI['d'])[idx_control_dates[i],1]))
  } # end for i
  
  temp=dat[idx_dates,1]
  
  model <- auto.arima(temp, start.p = 0, max.p = 10, d = 0, start.q = 0, max.q = 10)
  
  save(model, file = paste0("~/Documents/LSCE/SWG/t2m/case_0/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1979_1997.RData"))
} 

#### simulation 1998-2017 nat ####
DATE_proj = atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day'))
LENGTH_proj = dim(atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day')))[1]

DATE = atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day'))

for(i in 1:100){
  
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  cat('Run', NUM, '\n')
  
  Sample = array(NaN,c(LENGTH_proj,1))
  for (seas in 1:4){
    
    idxdates = which(DATE_proj['m']==MON[seas,1] | DATE_proj['m']==MON[seas,2] | DATE_proj['m']==MON[seas,3])
    
    load(paste0("~/Documents/LSCE/SWG/t2m/case_0/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1979_1997.RData"))
    
    Sample[idxdates, ] = arima.sim(n=length(idxdates), list(ar = model$coef[1:model$arma[1]], ma = model$coef[(model$arma[1]+1):(model$arma[1]+model$arma[2])]),sd = sqrt(model$sigma2)) + tail(model$coef,1)
  }
  
  Sample=round(Sample,digits=2)
  
  filesimu = paste('~/Documents/LSCE/SWG/t2m/case_0/SIMU_FAR/SIMU_SWG_tmean_1979_1997', '_run_', NUM, '.RData', sep="")
  save(Sample,DATE,file = filesimu)
  
}

#### estimation 1998-2017 ####
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

for (seas in 1:4){
  DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  
  # choose for the right months for season 
  idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  
  # DATE_ERAI only in season
  DATE_ERAI = DATE_ERAI[idxdates,]
  
  dat = tmean
  NB_STATION = 1
  
  idx_control_dates = which((DATE_ERAI['Y']>=1998) & (DATE_ERAI['Y']<=2017))
  
  # cherche les indices de jour dans DATE_ECA qui correspondent aux dates ERAI pour la période de calibration
  idx_dates=array(0,dim=length(idx_control_dates))
  for (i in 1:length(idx_control_dates)){
    idx_dates[i]=which(((DATE_OBS['Y'])[,1]==(DATE_ERAI['Y'])[idx_control_dates[i],1])
                       & ((DATE_OBS['m'])[,1]==(DATE_ERAI['m'])[idx_control_dates[i],1])
                       & ((DATE_OBS['d'])[,1]==(DATE_ERAI['d'])[idx_control_dates[i],1]))
  } # end for i
  
  temp=dat[idx_dates,1]
  
  model <- auto.arima(temp, start.p = 0, max.p = 10, d = 0, start.q = 0, max.q = 10)
  
  save(model, file = paste0("~/Documents/LSCE/SWG/t2m/case_0/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1998_2017.RData"))
} 


#### simulation 1998-2017 rea ####
DATE_proj = atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day'))
LENGTH_proj = dim(atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day')))[1]

DATE = atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day'))

for(i in 1:100){
  
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  cat('Run', NUM, '\n')
  
  Sample = array(NaN,c(LENGTH_proj,1))
  for (seas in 1:4){
    
    idxdates = which(DATE_proj['m']==MON[seas,1] | DATE_proj['m']==MON[seas,2] | DATE_proj['m']==MON[seas,3])
    
    load(paste0("~/Documents/LSCE/SWG/t2m/case_0/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1998_2017.RData"))
    
    Sample[idxdates, ] = arima.sim(n=length(idxdates), list(ar = model$coef[1:model$arma[1]], ma = model$coef[(model$arma[1]+1):(model$arma[1]+model$arma[2])]),sd = sqrt(model$sigma2)) + tail(model$coef,1)
  }
  
  Sample=round(Sample,digits=2)
  
  filesimu = paste('~/Documents/LSCE/SWG/t2m/case_0/SIMU_FAR/SIMU_SWG_tmean_1998_2017', '_run_', NUM, '.RData', sep="")
  save(Sample,DATE,file = filesimu)
  
}


#### FAR ####
Sample.list = vector("list", 100)
names(Sample.list) = paste0("Sample_", 1:100)
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(paste0("~/Documents/LSCE/SWG/t2m/case_0/SIMU_FAR/SIMU_SWG_tmean_1979_1997_run_", NUM, ".RData"))
  Sample.list[[i]] = Sample
}
Sample.nat = do.call(unlist(cbind.data.frame), Sample.list)
save(Sample.nat, file = paste0("~/Documents/LSCE/SWG/t2m/case_0/Sample_1998_2017_nat.RData"))

Sample.list = vector("list", 100)
names(Sample.list) = paste0("Sample_", 1:100)
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(paste0("~/Documents/LSCE/SWG/t2m/case_0/SIMU_FAR/SIMU_SWG_tmean_1998_2017_run_", NUM, ".RData"))
  Sample.list[[i]] = Sample
}
Sample.rea = do.call(unlist(cbind.data.frame), Sample.list)
save(Sample.rea, file = paste0("~/Documents/LSCE/SWG/t2m/case_0/Sample_1998_2017_rea.RData"))

thres = 10
apply(Sample.data, MARGIN = 1, FUN = function(x){sum(x>=thres)/length(x)})