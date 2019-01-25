#### load ncdf4 package ####
library(ncdf4)
tail(as.Date(1:24957, origin = "1949-12-31")) # "2018-04-30"
tail(as.Date(1:25202, origin = "1949-12-31")) # "2018-12-31"

#### tmean_48_50_N_01_04_E.RData ####
file.dir = "~/Documents/LSCE/"
filenames = list.files(paste0(file.dir, "tmean_48_50_N_01_04_E"))
nc <- nc_open("~/Documents/LSCE/tmean_48_50_N_01_04_E/igridensembles_05_tg_0001.25_048.25_n.nc")
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
  nc <- nc_open(paste0("~/Documents/LSCE/tmean_48_50_N_01_04_E/igridensembles_05_tg_000", coord_grid[k,1], "_0", coord_grid[k,2], "_n.nc"))
  data[,k] = ncvar_get(nc, "tg")[1:nt]
  nc_close(nc)
}

tmean = as.data.frame(apply(data, MARGIN = 1, FUN = mean))
colnames(tmean) = c("tmean")
range(tmean)
save(tmean, file = "~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E.RData")

which(rownames(tmean) == "1979-01-01")
which(rownames(tmean) == "2017-12-31")
tmean = tmean[which(rownames(tmean) == "1979-01-01"):which(rownames(tmean) == "2017-12-31"),]
tmean = as.data.frame(tmean)
rownames(tmean) = time[which(time=="1979-01-01"):which(time=="2017-12-31")]
save(tmean, file = "~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")


#### tmean_47_51_N_10_15_E.RData ####
file.dir = "~/Documents/LSCE/"
filenames = list.files(paste0(file.dir, "tmean_47_51_N_10_15_E"))
nc <- nc_open("~/Documents/LSCE/tmean_47_51_N_10_15_E/igridensembles_05_tg_0010.25_047.25_n.nc")
print(nc)
data = ncvar_get(nc, "tg")
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
  nc <- nc_open(paste0("~/Documents/LSCE/tmean_47_51_N_10_15_E/igridensembles_05_tg_00", coord_grid[k,1], "_0", coord_grid[k,2], "_n.nc"))
  data[,k] = ncvar_get(nc, "tg")[1:nt]
  nc_close(nc)
}

tmean = as.data.frame(apply(data, MARGIN = 1, FUN = mean))
colnames(tmean) = c("tmean")
range(tmean)
save(tmean, file = "~/Documents/LSCE/SWG/tmean_47_51_N_10_15_E.RData")


#### Estimation parameters ####
rm(list = ls())

library(stats4)
library(VGAM)
library(timeDate)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

for (seas in 1:4){
  # seas = 1
  cat(season[seas],'\n')
  
  DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  DATE_ERAI = DATE_ERAI[idxdates,]
  
  load(paste("~/Documents/LSCE/SWG/t2m/geop_",SEAS[seas],"_PCS.RData",sep=""))
  PCS = PCS[,1:20]
  
  load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E.RData")
  dat = tmean
  
  NB_STATION = 1
  
  idx_control_dates = which((DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']<=2004))  # from 1979 to 2004 total 26 years
  
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
      tt = data.frame(TT=temp,PCS)
      fit_stations_tt[[(i)]] = try(vgam(fmlatt ,uninormal(zero=NULL), data=tt,maxit=1000, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=T)
      ## cat(fit_stations_tt[[i]]@criterion$loglikelihood,'\n')
    }                          
    
  } # end for i NB_STATION
  
  filename=paste('~/Documents/LSCE/SWG/SWG_ERAI_ESD_tmean_',season[seas],'_cross-val.RData',sep="")
  cat(filename,'\n')
  save(fit_stations_tt,file = filename)

} # end for seas


#### Simulation ####
rm(list = ls())

DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

DATE_proj = atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day'))
LENGTH_proj = dim(atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day')))[1]

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
  idxdates=which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  DATE_ERAI=DATE_ERAI[idxdates,]
  
  idx_proj_dates = which( (DATE_ERAI['Y']>=2005) & (DATE_ERAI['Y']<=2017) )    # 2005 - 2017
  
  load(paste("~/Documents/LSCE/SWG/t2m/geop_",SEAS[seas],"_PCS.RData",sep=""))
  
  PCS=PCS[idx_proj_dates,1:20] # Pour T
  NB_pred=dim(PCS)[2]
  
  
  # Nom du fichier sauvegardé contenant les estimations
  filename=paste('~/Documents/LSCE/SWG/SWG_ERAI_ESD_tmean_',season[seas],'_cross-val.RData',sep="")
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
  
  td = which(DATE_OBS$Y==2005)[1]-1
  
  mean[(idxdates[idx_proj_dates] - td),] = tempmean   ; cat(range(tempmean), "\n")
  sd[(idxdates[idx_proj_dates] - td),] = tempsd       ; cat(range(tempsd)  , "\n")
  
} # end for seas

save(sd, mean, file = "~/Documents/LSCE/SWG/tmean_mean_sd.RData")
load("~/Documents/LSCE/SWG/tmean_mean_sd.RData")

#### Sampling ####

DATE = atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day'))

Sample_gaussian_temp<- function(par_tt1,par_tt2){
  a=date()
  cat(a,'\n')
  Nb_stations=dim(par_tt1)[2]
  echant=array(NaN,dim=dim(par_tt1))
  for (station in 1:Nb_stations){
    for (j in 1:(dim(par_tt1)[1])){
      echant[j,station]=rnorm(1,mean=par_tt1[j,station],sd=par_tt2[j,station])
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
  
  Sample_temp = Sample_gaussian_temp(mean,sd)
  Sample=array(NaN,dim=dim(sd))
  Sample[,idx]=Sample_temp
  
  Sample=round(Sample,digits=2)
  
  #  filesimu = paste('/home/users/mvrac/R/ESD/Output/SWG/SIMU_SWG_',VAR,'_ERAI_ESD_Claris_',experiment,'_run_',NUM,'.RData',sep="")
  filesimu = paste('~/Documents/LSCE/SWG//SIMU_tmean/SIMU_SWG_tmean', '_run_', NUM, '.RData', sep="")
  save(Sample,DATE,file = filesimu)
  
}


#### Diagnosis - distribution ####
rm(list=ls())

load("~/Documents/LSCE/SWG/SIMU_tmean/SIMU_SWG_tmean_run_001.RData")
Sample = Sample[,1]
plot(Sample)
range(Sample) # -8.69 27.67

load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
Obs = tmean[which(rownames(tmean)=="2005-01-01"):which(rownames(tmean)=="2017-12-31"),]
plot(Obs)
range(Obs) # -7.878333 27.243333

library(goftest)
hist(Obs)
hist(Sample)

library(CDFt)
res <- CramerVonMisesTwoSamples(Sample, Obs)
pvalue = 1/6*exp(-res)

library(RVAideMemoire)
CvM.test(Sample, Obs)


# qqplot
qqplot(Obs, Sample, xlab = "obs", ylab = "sim", pch = 20,
       main = "tmean", cex.main = 1.5, cex.lab = 1)
abline(0,1,col="red")


#### Diagnosis - autoregression ####
library(ggplot2)

load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
Obs = tmean[which(rownames(tmean)=="2005-01-01"):which(rownames(tmean)=="2017-12-31"),]

for (k in 1:100){
  
  NUM=cbind(formatC(k, digits = 0, width = 3, format = "f", flag = "0"))
  
  load(paste0("~/Documents/LSCE/SWG/SIMU_tmean/SIMU_SWG_tmean_run_", NUM, ".RData"))
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
  
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/SIMU_tmean/acf_tmean_run_", NUM, ".pdf"), width = 9, height = 7)
}



#### FAR ####
source("fun_FAR.R")

load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
dat = tmean[,1]

DATE = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

Seq_nat  = dat[DATE$Y>=1979 & DATE$Y<=1997]
Seq_real = dat[DATE$Y>=1998 & DATE$Y<=2017]
  
FAR(R_nat = Seq_nat, R_real = Seq_real, threshold = 20)

x = seq(from = 20, to = min(max(Seq_nat), max(Seq_real)), by = 0.01)
plot(x, FAR(Seq_nat, Seq_real, x), type = "l")





