# datadir = "~/Documents/LSCE/SWG/"

library(ncdf4)
library(timeDate)

nc <- nc_open("~/Documents/LSCE/nc/precip_total.1979-2017_FR.nc")
print(nc)
data = ncvar_get(nc, "tp")
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)

nt = length(time)
dat = data[,,seq(from = 1, by = 2, to = nt)] + data[,,seq(from = 2, by = 2, to = nt)]
dat = dat*10^3
range(dat)

nt = nt/2
dim(dat) = c(length(lon)*length(lat), nt)
dat = t(dat)
dat = apply(dat, 1, mean)

load("~/Documents/LSCE/SWG/tp/gridobs_tp.RData")

## weighting
pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = grids_tp*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = grids_tp[,pix]/scale[pix]
}

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = T)
PCS = as.data.frame(pca$x[,1:12])
save(PCS, file = paste0("~/Documents/LSCE/SWG/tp/", SEAS[3], "_PCS.RData"))


#### grids_total_precip ####
grids_total_precip = as.data.frame(apply(grids_tp, 1, mean))
dim(grids_total_precip)
save(grids_total_precip, file = "~/Documents/LSCE/SWG/tp/grids_total_precip.RData")


######## Estimation des parametres de SWG pour tp ################
experiment = "kfold"

library(stats4)
library(VGAM)
library(timeDate)

threshocc=0.1
threshfall=0.1
thresh_name=0
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
DATE_OBS  = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

for(ifold in 1:2){
  # ifold = 1
  a=date()
  cat("fold",ifold,a,'\n')
  
  for(seas in 1:4){
    # seas = 1
    cat(season[seas],'\n')
    
    
    DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2012-12-31",by='day'))
    DATE_OBS  = atoms(timeSequence(from="1979-01-01",to="2012-12-10",by='day'))
    
    
    idxdates=which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
    DATE_ERAI=DATE_ERAI[idxdates,]
    
    load(paste("~/Documents/LSCE/SWG/tp/",SEAS[seas],"_PCS.RData",sep=""))
    
    load("~/Documents/LSCE/SWG/tp/gridobs_tp.RData")
    #load(paste(datadir,"ECA_stations_obs/ECA_stationobs.RData",sep=""))
    
    load("~/Documents/LSCE/SWG/tp/grids_total_precip.RData")
    
    NB_STATION = ncol(grids_total_precip)
    #rm(LAT,LON)
    
    
    ############################################
    ########### in the following, ifold corresponds to the calibration period ############
    ############################################
    #    if(ifold==1){idx_control_dates=which(DATE_ERAI['Y']>1984)}
    #    if(ifold==2){idx_control_dates=which((DATE_ERAI['Y']<1985) | (DATE_ERAI['Y']>1990))}
    #    if(ifold==3){idx_control_dates=which((DATE_ERAI['Y']<1991) | (DATE_ERAI['Y']>1996))}
    #    if(ifold==4){idx_control_dates=which((DATE_ERAI['Y']<1997) | (DATE_ERAI['Y']>2002))}
    #    if(ifold==5){idx_control_dates=which(DATE_ERAI['Y']<2003)}
    
    
    if(ifold==1){
      idx_control_dates=which( (DATE_ERAI['Y']>=1989) & (DATE_ERAI['Y']<=1998) ) # 1984 - 2003
      idx_proj_dates=which(DATE_ERAI['Y']<=1988) # 1979 - 1983
    }
    if(ifold==2){
      idx_control_dates=which(DATE_ERAI['Y']<=1988)
      idx_proj_dates=which( (DATE_ERAI['Y']>=1989) & (DATE_ERAI['Y']<=1998) ) 
    }
    
    
    pc_name=names(PCS)
    PCS=PCS[idx_control_dates,]
    
    # cherche les indices de jour dans DATE_ECA qui correspondent aux dates ERAI pour la période de calibration
    idx_dates=array(0,dim=length(idx_control_dates))
    for (i in 1:length(idx_control_dates)){
      idx_dates[i]=which(((DATE_OBS['Y'])[,1]==(DATE_ERAI['Y'])[idx_control_dates[i],1])
                         & ((DATE_OBS['m'])[,1]==(DATE_ERAI['m'])[idx_control_dates[i],1])
                         & ((DATE_OBS['d'])[,1]==(DATE_ERAI['d'])[idx_control_dates[i],1]))
    } # end for i
    
    
    
    
    ## Rain occurence, logistic regression,vglm
    fmlaro = as.formula(paste("RO ~ ", paste(pc_name, collapse= "+")))
    fit_stations_ro = vector("list", NB_STATION)
    statio_pix_ro = array("Non-statio", NB_STATION)
    
    for (i in 1:NB_STATION){
      
      idx1 = which(grids_total_precip[idx_dates,i] > threshocc) ### & stations_eca_precip[idx_dates,i]<1000
      if(length(idx1)<2){
        fit_stations_ro[[i]]=NaN
      }
      
      if(length(idx1)>=2){    
        rain_occur=array(0,dim=length(idx_dates))
        rain_occur[idx1]=1
        ro=data.frame(RO=rain_occur,PCS)
        
        #fit_stations_ro[[i]]=try(vglm(fmlaro ,binomialff(mv=TRUE), data=ro, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=F)
        fit_stations_ro[[i]]=try(vglm(fmlaro ,binomialff(multiple.responses=TRUE), data=ro, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=F)
        # mv=TRUE pour que vglm tienne compte des PCs (sinon, il calucle une proba constante).
        # x.arg= FALSE pour qu'il ne garde pas les elements servant au calcul des parametres.
        
        if(length(slotNames(fit_stations_ro[[i]]))==0){
          fit_stations_ro[[i]] = length(which(grids_total_precip[idx_dates,i] > threshocc))/length(grids_total_precip[idx_dates,i])
          statio_pix_ro[i]='statio'
        }
        
        if(length(slotNames(fit_stations_ro[[i]]))!=0){
          if(length(slotNames(fit_stations_ro[[i]]))>0 & is.finite(fit_stations_ro[[i]]@criterion$loglikelihood)!=TRUE){
            fit_stations_ro[[i]] = length(which(grids_total_precip[idx_dates,i] > threshocc))/length(grids_total_precip[idx_dates,i])
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
      
      idx1 = which(grids_total_precip[idx_dates,i] > threshfall) ### & grids_total_precip[idx_dates,i]<1000
      
      if(length(idx1)<2){
        fit_stations_rf[[i]]=NaN
      }
      
      if(length(idx1)>=2){    
        
        rain=grids_total_precip[idx_dates[idx1],i]
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
    filename=paste('~/Documents/LSCE/SWG/tp/SWG_ERAI_PR_',season[seas],'_fold_',ifold,'.RData',sep="")
    cat(filename,'\n')
    save(statio_pix_ro,fit_stations_ro,statio_pix_rf,fit_stations_rf, file = filename)
    
    
  } # end for seas
} # end for ifold








