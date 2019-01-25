# datadir = "~/Documents/LSCE/SWG/"

library(ncdf4)
library(timeDate)

nc <- nc_open("~/Documents/LSCE/nc/temp_2m.1979-2017_FR.nc")
print(nc)
data = ncvar_get(nc, "t2m") 
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)

nt = length(time)
dat = data[,,seq(from = 1, by = 4, to = nt)] + data[,,seq(from = 2, by = 4, to = nt)] +
      data[,,seq(from = 3, by = 4, to = nt)] + data[,,seq(from = 4, by = 4, to = nt)]
dat = dat/4 - 273.15

nt = nt/4
dim(dat) = c(length(lon)*length(lat), nt)
dat = t(dat)
dat = apply(dat, 1, mean)
dat = as.data.frame(dat)
t2m = dat[,1]


######## Estimation des parametres de SWG pour t2m ################
library(stats4)
library(VGAM)
library(timeDate)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_OBS  = atoms(timeSequence(from="1979-01-01",to="2004-12-31",by='day'))

for(seas in 1:4){
  # seas = 1
  cat(season[seas],'\n')
  
  DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2004-12-31",by='day'))
  idxdates=which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  DATE_ERAI=DATE_ERAI[idxdates,]
  
  # load PCA
  load(paste("~/Documents/LSCE/SWG/t2m/",SEAS[seas],"_PCS.RData",sep=""))
  # load data
  # load("~/Documents/LSCE/SWG/t2m/gridobs.RData")
  
  NB_STATION = 1
  
  idx_control_dates = which((DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']>=2004))
  
  PCS=PCS[idx_control_dates,1:10]
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
      tt=data.frame(TT=temp,PCS)
      fit_stations_tt[[(i)]] = try(vgam(fmlatt ,uninormal(zero=NULL), data=tt,maxit=1000, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=T)
      ## cat(fit_stations_tt[[i]]@criterion$loglikelihood,'\n')
    }                          
    
  } # end for i NB_STATION
  
  
  filename=paste('~/Documents/LSCE/SWG/t2m/SWG_t2M_ESD','_',season[seas],'.RData',sep="")
  cat(filename,'\n')
  save(fit_stations_tt,file = filename)
  
} # end for seas

load("~/Documents/LSCE/SWG/t2m/SWG_t2M_ESD_Winter.RData")
fit_stations_tt

VAR = "t2m"

library(stats4)
library(VGAM)
library(timeDate)

DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

mean = array(NaN,c(dim(DATE_ERAI)[1],1))
sd = array(NaN,c(dim(DATE_ERAI)[1],1))

###save=paste('PARAM_SWG_ERAI_VALUE_ECA_TT.RData',sep="") # Nom du fichier qui sera sauvegardé

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')


for(seas in 1:4){
  
  cat(season[seas],'\n')
  
  #    DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2008-12-31",by='day'))###limite à 2008 mem si donnees jusque 2012
  DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  idxdates=which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
  DATE_ERAI=DATE_ERAI[idxdates,]
  
  load(paste("~/Documents/LSCE/SWG/t2m/",SEAS[seas],"_PCS.RData",sep=""))
  
  
  
  idx_proj_dates = which((DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']<=2017))
  
  
  
  PCS=PCS[idx_proj_dates,] # Pour T
  NB_pred=dim(PCS)[2]
  
  
  # Nom du fichier sauvegardé contenant les estimations
  filename=paste('~/Documents/LSCE/SWG/t2m/SWG_t2M_ESD','_',season[seas],'.RData',sep="")
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
  
  mean[idxdates[idx_proj_dates],] = tempmean
  sd[idxdates[idx_proj_dates],] = tempsd
  
}#end seas


save(sd,mean,file = paste('~/Documents/LSCE/SWG/t2m/FieldParameters_SWG_ERAI_ESD_',VAR,'.RData',sep=""))

load(paste('~/Documents/LSCE/SWG/t2m/FieldParameters_SWG_ERAI_ESD_',VAR,'.RData',sep=""))



##################
#####sampling########
##################

DATE = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

#datadir="/home/estimr1/pvait/RDATA/"

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

idx=which(sd[1,]!='NaN')

for(i in 2:2){
  
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  cat('Run',NUM,'\n')
  
  Sample_temp = Sample_gaussian_temp(par_tt1 = mean, par_tt2 = sd)
  Sample=array(NaN,dim=dim(sd))
  Sample[,idx]=Sample_temp
  
  Sample=round(Sample,digits=2)
  
  filesimu = paste('~/Documents/LSCE/SWG/t2m/SIMU_SWG_',VAR,'_run_',NUM,'.RData',sep="")
  save(Sample,DATE,file = filesimu)
  
}
## q(save='no')

Sample = echant[,1]
range(Sample)
plot(Sample)
plot(t2m)



