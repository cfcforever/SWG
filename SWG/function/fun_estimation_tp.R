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
  
  load(paste0('~/Documents/LSCE/SWG/slp_SD_diagnosis/predictor/predictor_', SEAS[seas], '_1.RData'))
  colnames(PCS) = paste0("PC", 1:ncol(PCS))
  
  load(paste0("~/Documents/LSCE/SWG/data/tp_", city, "_1979_2017.RData"))
  dat = tp
  
  NB_STATION = 1
  
  idx_control_dates = which((DATE_ERAI['Y']>=year.begin) & (DATE_ERAI['Y']<=year.end))  
  idx_proj_dates    = which((DATE_ERAI['Y']>=1999) & (DATE_ERAI['Y']<=2017))  
  
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
  output = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SWG_ERAI_ESD_tp_", season[seas], "_cross-val_", year.begin, "_", year.end, "_1.RData")
  cat(output,'\n')
  save(statio_pix_ro,fit_stations_ro,statio_pix_rf,fit_stations_rf, file = output)
  
} # end for seas