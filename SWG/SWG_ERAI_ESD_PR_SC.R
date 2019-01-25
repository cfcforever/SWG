
######## Estimation des parametres de SWG pour PR pour CORDEX-ESD ################


#experiment = "cross-val"
experiment = "kfold"


library(stats4)
library(VGAM)
library(timeDate)
#datadir="/home/estimr1/pvait/RDATA/"

threshocc=0
threshfall=0
thresh_name=0
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')



############################################
############################################


#DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2008-12-31",by='day'))###limite à 2008 mem si donnees jusque 2012
DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2012-12-31",by='day'))

DATE_OBS= atoms(timeSequence(from="1979-01-01",to="2012-12-10",by='day'))

############################################
if(experiment == "cross-val"){
  
}
##################

############################################
if(experiment == "kfold"){
  
  for(ifold in 1:5){
    # ifold = 1
    a=date()
    cat("fold",ifold,a,'\n')
    
    for(seas in 1:4){
      # seas = 1
      cat(season[seas],'\n')
      
      #DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2008-12-31",by='day'))###limite à 2008 mem si donnees jusque 2012
      DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2012-12-31",by='day'))
      
      DATE_OBS= atoms(timeSequence(from="1979-01-01",to="2012-12-10",by='day'))
      
      
      idxdates=which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
      DATE_ERAI=DATE_ERAI[idxdates,]
      
      #datafile=paste(datadir,"ERAI_1.125/VALUE/",SEAS[seas],"_PCS_ERAI_79_12_EU_CORDEX.RData",sep="")
      #load(datafile)
      load(paste("~/Documents/LSCE/RDATA/",SEAS[seas],"_PCS_ERAI_79_12_SAM_CORDEX.RData",sep=""))
      
      
      
      load("~/Documents/LSCE/RDATA/Claris_stations_obs/Claris_stationobs.RData")
      #load(paste(datadir,"ECA_stations_obs/ECA_stationobs.RData",sep=""))
      
      NB_STATION = length(LAT)
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
        idx_control_dates=which( (DATE_ERAI['Y']>=1984) & (DATE_ERAI['Y']<=2003) ) # 1984 - 2003
        idx_proj_dates=which(DATE_ERAI['Y']<=1983) # 1979 - 1983
      }
      if(ifold==2){
        idx_control_dates=which((DATE_ERAI['Y']<=1983) | ((DATE_ERAI['Y']>=1989) & (DATE_ERAI['Y']<=2003) )) # 1979-1983 & 1989-2003
        idx_proj_dates=which((DATE_ERAI['Y']>=1984) & (DATE_ERAI['Y']<=1988)) # 1984-1988
      }
      if(ifold==3){
        idx_control_dates=which((DATE_ERAI['Y']<=1988) | ((DATE_ERAI['Y']>=1994) & (DATE_ERAI['Y']<=2003) )) # 1979-1988 & 1994-2003
        idx_proj_dates=which((DATE_ERAI['Y']>=1989) & (DATE_ERAI['Y']<=1993)) # 1989-1993
      }
      if(ifold==4){
        idx_control_dates=which((DATE_ERAI['Y']<=1993) | ((DATE_ERAI['Y']>=1999) & (DATE_ERAI['Y']<=2003)) ) # 1979-1993 & 1999-2003
        idx_proj_dates=which((DATE_ERAI['Y']>=1994) & (DATE_ERAI['Y']<=1998)) # 1994-1998
      }
      if(ifold==5){
        idx_control_dates=which(DATE_ERAI['Y']<=1998) # 1979 - 1998
        idx_proj_dates=which((DATE_ERAI['Y']>=1999) & (DATE_ERAI['Y']<=2003)) # 1999-2003
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
        
        idx1 = which(stations_claris_precip[idx_dates,i] > threshocc) ### & stations_eca_precip[idx_dates,i]<1000
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
            fit_stations_ro[[i]] = length(which(stations_claris_precip[idx_dates,i] > threshocc))/length(stations_claris_precip[idx_dates,i])
            statio_pix_ro[i]='statio'
          }
          
          if(length(slotNames(fit_stations_ro[[i]]))!=0){
            if(length(slotNames(fit_stations_ro[[i]]))>0 & is.finite(fit_stations_ro[[i]]@criterion$loglikelihood)!=TRUE){
              fit_stations_ro[[i]] = length(which(stations_claris_precip[idx_dates,i] > threshocc))/length(stations_claris_precip[idx_dates,i])
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
        
        idx1 = which(stations_claris_precip[idx_dates,i] > threshfall) ### & stations_claris_precip[idx_dates,i]<1000
        
        if(length(idx1)<2){
          fit_stations_rf[[i]]=NaN
        }
        
        if(length(idx1)>=2){    
          
          rain=stations_claris_precip[idx_dates[idx1],i]
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
      filename=paste('~/Documents/LSCE/SWG/SWG_ERAI_ESD_PR_',season[seas],'_fold_',ifold,'.RData',sep="")
      cat(filename,'\n')
      save(statio_pix_ro,fit_stations_ro,statio_pix_rf,fit_stations_rf, file = filename)
      
      
    } # end for seas
  } # end for ifold
  
  
} # end if xp kfold

############################################


q(save='no')


