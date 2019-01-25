#### load ncdf4 package ####
library(ncdf4)
tail(as.Date(1:24957, origin = "1949-12-31")) # "2018-04-30"
tail(as.Date(1:25202, origin = "1949-12-31")) # "2018-12-31"

#### precip_48_50_N_01_04_E.RData ####
filenames = list.files("~/Documents/LSCE/precip_48_50_N_01_04_E/")
nc <- nc_open("~/Documents/LSCE/precip_48_50_N_01_04_E/igridensembles_05_tg_0001.25_048.25_n.nc")
print(nc)
data = ncvar_get(nc, "tg")
sum(is.na(data))
range(data[!is.na(data)])

time = as.Date(1:24957, origin = "1949-12-31")
nt = length(time)
nc_close(nc)

lon = seq(from = 1.25, to = 3.75, by = 0.5); nx = length(lon)
lat = seq(from = 48.25, to =  49.75, by = 0.5); ny = length(lat)
coord_grid = as.data.frame(matrix(NA, nrow = nx*ny, ncol = 2))
colnames(coord_grid) = c("lon", "lat")
coord_grid$lon = rep(lon, each = ny)
coord_grid$lat = rep(lat, n = nx)

data = as.data.frame(matrix(NA, nrow = nt, ncol = nx*ny))
rownames(data) = time
for (k in 1:(nx*ny)){
  nc <- nc_open(paste0("~/Documents/LSCE/precip_48_50_N_01_04_E/igridensembles_05_tg_000", coord_grid[k,1], "_0", coord_grid[k,2], "_n.nc"))
  data[,k] = ncvar_get(nc, "tg")[1:nt]
  nc_close(nc)
}

precip = as.data.frame(apply(data, MARGIN = 1, FUN = mean))
colnames(precip) = c("precip")
range(precip)
save(precip, file = "~/Documents/LSCE/SWG/precip_48_50_N_01_04_E.RData")

which(rownames(precip) == "1979-01-01")
which(rownames(precip) == "2017-12-31")
precip = precip[which(rownames(precip) == "1979-01-01"):which(rownames(precip) == "2017-12-31"),]
precip = as.data.frame(precip)
rownames(precip) = time[which(time=="1979-01-01"):which(time=="2017-12-31")]
save(precip, file = "~/Documents/LSCE/SWG/precip_48_50_N_01_04_E_1979_2017.RData")


#### precip_47_51_N_10_15_E.RData ####
filenames = list.files("~/Documents/LSCE/precip_47_51_N_10_15_E/")
nc <- nc_open("~/Documents/LSCE/precip_47_51_N_10_15_E/igridensembles_05_rr_0010.25_047.25_n.nc")
print(nc)
data = ncvar_get(nc, "rr")
sum(is.na(data))
range(data[!is.na(data)])
time = ncvar_get(nc, "time")

time = as.Date(1:24957, origin = "1949-12-31")
nt = length(time)
nc_close(nc)

lon = seq(from = 10.25, to = 14.75, by = 0.5); nx = length(lon)
lat = seq(from = 47.25, to = 50.75, by = 0.5); ny = length(lat)
coord_grid = as.data.frame(matrix(NA, nrow = nx*ny, ncol = 2))
colnames(coord_grid) = c("lon", "lat")
coord_grid$lon = rep(lon, each = ny)
coord_grid$lat = rep(lat, n = nx)

data = as.data.frame(matrix(NA, nrow = nt, ncol = nx*ny))
rownames(data) = time
for (k in 1:(nx*ny)){
  nc <- nc_open(paste0("~/Documents/LSCE/precip_47_51_N_10_15_E/igridensembles_05_rr_00", coord_grid[k,1], "_0", coord_grid[k,2], "_n.nc"))
  data[,k] = ncvar_get(nc, "rr")[1:nt]
  nc_close(nc)
}

precip = as.data.frame(apply(data, MARGIN = 1, FUN = mean))
colnames(precip) = c("precip")
range(precip)
save(precip, file = "~/Documents/LSCE/SWG/precip_47_51_N_10_15_E.RData")

which(rownames(precip) == "1979-01-01")
which(rownames(precip) == "2017-12-31")
precip = precip[which(rownames(precip) == "1979-01-01"):which(rownames(precip) == "2017-12-31"),]
precip = as.data.frame(precip)
rownames(precip) = time[which(time=="1979-01-01"):which(time=="2017-12-31")]
save(precip, file = "~/Documents/LSCE/SWG/precip_47_51_N_10_15_E_1979_2017.RData")


#### Estimation parameters ####
rm(list = ls())

library(stats4)
library(VGAM)
library(timeDate)

threshocc=0
threshfall=0
thresh_name=0
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

for (seas in 1:4){
  # seas = 1
  cat(season[seas],'\n')
  
  DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  idxdates=which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  DATE_ERAI=DATE_ERAI[idxdates,]
  
  load(paste("~/Documents/LSCE/SWG/tp/slp_",SEAS[seas],"_PCS.RData",sep=""))
  PCS = PCS[,1:20]
  
  load("~/Documents/LSCE/SWG/precip_47_51_N_10_15_E.RData")
  dat = precip
  
  NB_STATION = 1
  
  idx_control_dates = which((DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']<=2004))  # from 1979 to 2004 total 26 years
  idx_proj_dates    = which((DATE_ERAI['Y']>=2005) & (DATE_ERAI['Y']<=2017))  # 2005 - 2017
  
  DATECALIB = DATE_ERAI[idx_control_dates,]
  DATEPROJ  = DATE_ERAI[idx_proj_dates,]
  
  PCS=PCS[idx_control_dates,]
  pc_name=names(PCS)
  
  # cherche les indices de jour dans DATE_ECA qui correspondent aux dates ERAI pour la période de calibration
  idx_dates=array(0,dim=length(idx_control_dates))
  for (i in 1:length(idx_control_dates)){
    idx_dates[i]=which(((DATE_OBS['Y'])[,1]==(DATE_ERAI['Y'])[idx_control_dates[i],1])
                       & ((DATE_OBS['m'])[,1]==(DATE_ERAI['m'])[idx_control_dates[i],1])
                       & ((DATE_OBS['d'])[,1]==(DATE_ERAI['d'])[idx_control_dates[i],1]))
  } # end for i
  
  
  ######## Rain occurence, logistic regression,vglm ######
  fmlaro = as.formula(paste("RO ~ ", paste(pc_name, collapse= "+")))
  fit_stations_ro = vector("list", NB_STATION)
  statio_pix_ro = array(NaN, NB_STATION)
  
  for (i in 1:NB_STATION){
    
    idx1 = which(dat[idx_dates,i] > threshocc) ### & stations_eca_precip[idx_dates,i]<1000
    if(length(idx1)<2){
      fit_stations_ro[[i]]=NaN
    }
    
    if(length(idx1)>=2){    
      rain_occur=array(0,dim=length(idx_dates))
      rain_occur[idx1]=1
      ro=data.frame(RO=rain_occur,PCS)
      
      fit_stations_ro[[i]]=try(vglm(fmlaro ,binomialff(multiple.responses=TRUE), data=ro, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=F)
    
      if(length(slotNames(fit_stations_ro[[i]]))==0){
        fit_stations_ro[[i]] = length(which(stations_claris_precip[idx_dates,i] > threshocc))/length(stations_claris_precip[idx_dates,i])
        statio_pix_ro[i]='statio'
      }
      
      if(length(slotNames(fit_stations_ro[[i]]))!=0){
        if(length(slotNames(fit_stations_ro[[i]]))>0 & is.finite(fit_stations_ro[[i]]@criterion$loglikelihood)!=TRUE){
          fit_stations_ro[[i]] = length(which(dat[idx_dates,i] > threshocc))/length(dat[idx_dates,i])
          statio_pix_ro[i]='statio'
        } # end if
      } # end if length(slotNames...) != 0
      
    } # end if length idx1 >= 2
    
  } # end for i NB__STATION

  
  ############## Rainfall amounts, gamma distribution, vglm
  
  fmlarf = as.formula(paste("RF ~ ",paste(pc_name, collapse= "+")))
  fit_stations_rf <- vector("list", NB_STATION)
  statio_pix_rf<- array(NaN, NB_STATION)
  
  for (i in 1:NB_STATION){
    
    idx1 = which(dat[idx_dates,i] > threshfall) ### & stations_claris_precip[idx_dates,i]<1000
    
    if(length(idx1)<2){
      fit_stations_rf[[i]]=NaN
    }
    
    if(length(idx1)>=2){    
      
      rain=dat[idx_dates[idx1],i]
      rf=data.frame(RF=rain,PCS[idx1,])
      
      fit_stations_rf[[i]] = try(vgam(fmlarf ,gammaR(zero=NULL), data=rf,maxit=1000, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=T)#, trace = TRUE, crit="c"
      
      if(length(slotNames(fit_stations_rf[[i]]))==0){
        fit_stations_rf[[i]] = try(vgam(rain~1 ,gammaR(zero=NULL), data=rf,maxit=1000, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=T)
        statio_pix_rf[i]='statio'
      }
      if(length(slotNames(fit_stations_rf[[i]]))!=0){
        if(length(slotNames(fit_stations_rf[[i]]))>0 & is.finite(fit_stations_rf[[i]]@criterion$loglikelihood)!=TRUE){
          fit_stations_rf[[i]] = try(vgam(rain~1 ,gammaR(zero=NULL), data=rf,maxit=1000, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=T)
          statio_pix_rf[i]='statio'
        }
      }	
    }
    
    cat("station ",i,": \n")
    warnings()
    
    
  } # end for i
  
  cat(which(statio_pix_ro=='statio'),'\n')
  cat(which(statio_pix_rf=='statio'),'\n')
  a=date()
  cat(a,'\n')
  
  # Nom du fichier de sortie en RData
  filename=paste('~/Documents/LSCE/SWG/SWG_ERAI_ESD_precip_',season[seas],'_cross-val.RData',sep="")
  cat(filename,'\n')
  save(statio_pix_ro,fit_stations_ro,statio_pix_rf,fit_stations_rf, file = filename)
  
} # end for seas


#### Simulation ####
rm(list = ls())

library(stats4)
library(VGAM)
library(timeDate)

DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

DATE_proj = atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day'))
LENGTH_proj = dim(atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day')))[1]

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
  
  idx_proj_dates = which( (DATE_ERAI['Y']>=2005) & (DATE_ERAI['Y']<=2017) )    # 2005 - 2017
  
  load(paste("~/Documents/LSCE/SWG/t2m/geop_",SEAS[seas],"_PCS.RData",sep=""))
  
  PCS=PCS[idx_proj_dates,1:20]
  NB_pred=dim(PCS)[2]
  
  
  # Nom du fichier sauvegardé contenant les estimations
  filename=paste('~/Documents/LSCE/SWG/SWG_ERAI_ESD_precip_',season[seas],'_cross-val.RData',sep="")
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
  
  td = which(DATE_OBS$Y==2005)[1]-1
  
  P1[(idxdates[idx_proj_dates] - td),] = tempP1       ; cat(range(tempP1),    "\n")
  shape[(idxdates[idx_proj_dates] - td),] = tempshape ; cat(range(tempshape), "\n")
  rate[(idxdates[idx_proj_dates] - td),] = temprate   ; cat(range(temprate),  "\n")
  
} # end for seas

save(P1, shape, rate, file = "~/Documents/LSCE/SWG/precip_P1_shape_rate.RData")
# load("~/Documents/LSCE/SWG/precip_P1_shape_rate.RData")

#### Sampling ####
rm(list=ls())

DATE = atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day'))

Sample_binomocc_gammaity<- function(par_ro,par_rf1,par_rf2,thresh=0){
  a=date()
  cat(a,'\n')
  
  Nb_stations=dim(par_ro)[2]
  echant=array(NaN,dim=dim(par_ro))
  for (station in 1:Nb_stations){
    simul_pluie=runif(dim(par_ro)[1],0,1)
    for (j in 1:(dim(par_ro)[1])){
      if (simul_pluie[j]<par_ro[j,station]){
        p=runif(1)
        C=pgamma(thresh,shape=par_rf1[j,station],scale=1/par_rf2[j,station])
        echant[j,station]=qgamma(p*(1-C)+C,shape=par_rf1[j,station],scale=1/par_rf2[j,station])
      }
      else {echant[j,station]=0}
    }
  }
  a=date()
  cat(a,'\n')
  return(echant)
}

idx = 1
for(i in 1:100){
  
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  cat('Run',NUM,'\n')
  
  Sample_temp = Sample_binomocc_gammaity(P1,shape,rate)
  Sample=array(NaN,dim=dim(P1))
  Sample[,idx]=Sample_temp
  
  Sample=round(Sample,digits=2)
  
  filesimu = paste('~/Documents/LSCE/SWG/SIMU_precip/SIMU_SWG_precip', '_run_', NUM, '.RData', sep="")
  save(Sample,DATE,file = filesimu)
  
}


#### Diagnosis - distribution ####
rm(list=ls())

load("~/Documents/LSCE/SWG/SIMU_precip/SIMU_SWG_precip_run_006.RData")
Sample = Sample[,1]
plot(Sample)
range(Sample)

load("~/Documents/LSCE/SWG/precip_47_51_N_10_15_E_1979_2017.RData")
Obs = precip[which(rownames(precip)=="2005-01-01"):which(rownames(precip)=="2017-12-31"),]
plot(Obs)
range(Obs) 

hist(Obs)
hist(Sample)

library(goftest)


library(CDFt)
res <- CramerVonMisesTwoSamples(Sample, Obs)
pvalue = 1/6*exp(-res)
pvalue

# library(RVAideMemoire)

# qqplot
for (k in 1:100){
  
  NUM=cbind(formatC(k, digits = 0, width = 3, format = "f", flag = "0"))
  
  load(paste0("~/Documents/LSCE/SWG/SIMU_precip/SIMU_SWG_precip_run_", NUM, ".RData"))
  Sample = Sample[,1]
  
  pmax = max(c(Obs, Sample))
  
  qqplot(Obs, Sample, xlab = "obs", ylab = "sim", pch = 20,
         main = paste0("precip - Obs vs Sample_", NUM), cex.main = 1.5, cex.lab = 1, xlim = c(0,pmax), ylim = c(0,pmax))
  abline(0,1,col="red")
  
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/SIMU_precip/qqplot_precip_run_", NUM, ".pdf"), width = 7, height = 7)
}



#### Diagnosis - autoregression ####
library(ggplot2)

load("~/Documents/LSCE/SWG/precip_47_51_N_10_15_E_1979_2017.RData")
Obs = precip[which(rownames(precip)=="2005-01-01"):which(rownames(precip)=="2017-12-31"),]

for (k in 1:100){
  
  NUM=cbind(formatC(k, digits = 0, width = 3, format = "f", flag = "0"))
  
  load(paste0("~/Documents/LSCE/SWG/SIMU_precip/SIMU_SWG_precip_run_", NUM, ".RData"))
  Sample = Sample[,1]
  
  nlag = 31
  data_acf = as.data.frame(matrix(data = NA, nrow = 2*nlag, ncol = 3))
  
  acf_obs = acf(Obs, plot = F, lag.max = nlag-1)
  acf_obs = with(acf_obs, data.frame(lag, acf))
  acf_obs = acf_obs[-1,]
  acf_sim = acf(Sample, plot = F, lag.max = nlag-1)
  acf_sim = with(acf_sim, data.frame(lag, acf))
  acf_sim = acf_sim[-1,]
  data_acf = rbind(acf_obs, acf_sim)
  data_acf$type = c(rep("obs", nlag-1), rep("sim", nlag-1))
  
  color = c("black", "red")
  {
    p = ggplot(data = data_acf, mapping = aes(x = lag, y = acf))
    p = p + geom_bar(stat = "identity", aes(fill = type), position = "dodge")
    p = p + scale_fill_manual(values = color)
    p = p + xlab("day")
    p = p + ggtitle(paste0("acf - Obs vs Sample_", k))
    print(p)
  }
  
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/SIMU_precip/acf_precip_run_", NUM, ".pdf"), width = 9, height = 7)
}



#### FAR ####
source("fun_FAR.R")

load("~/Documents/LSCE/SWG/precip_47_51_N_10_15_E_1979_2017.RData")
dat = precip[,1]

DATE = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

Seq_nat  = dat[DATE$Y>=1979 & DATE$Y<=1997]; range(Seq_nat)
Seq_real = dat[DATE$Y>=1998 & DATE$Y<=2017]; range(Seq_real)

FAR(R_nat = Seq_nat, R_real = Seq_real, threshold = 20)

x = seq(from = 30, to = min(max(Seq_nat), max(Seq_real)), by = 0.1)
plot(x, FAR(Seq_nat, Seq_real, x), type = "l")

