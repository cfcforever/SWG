
######## Simulations par SWG 2004-2012 PR pour CORDEX-ESD ################



#experiment = "cross-val"
experiment = "kfold"


library(stats4)
library(VGAM)
library(timeDate)



###DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2008-12-31",by='day'))
#DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2012-12-31",by='day'))
DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2003-12-31",by='day'))
DATE_ERAI2= atoms(timeSequence(from="2004-01-01",to="2012-12-31",by='day'))

#datadir="/home/estimr1/pvait/RDATA/"

P1 = array(NaN,c(dim(DATE_ERAI)[1],81)) # matrice de probas d'occurrence
shape = array(NaN,c(dim(DATE_ERAI)[1],81))
rate = array(NaN,c(dim(DATE_ERAI)[1],81))

P12 = array(NaN,c(5,dim(DATE_ERAI2)[1],81)) # matrice de probas d'occurrence
shape2 = array(NaN,c(5,dim(DATE_ERAI2)[1],81))
rate2 = array(NaN,c(5,dim(DATE_ERAI2)[1],81))



MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

for(ifold in 1:5){
  a=date()
  cat("fold",ifold,a,'\n')
  
  ###DATE_ECA= atoms(timeSequence(from="1979-01-01",to="2012-12-31",by='day'))
  DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2012-12-10",by='day'))
  
  
  for(seas in 1:4){
    
    cat(season[seas],'\n')
    
    #DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2008-12-31",by='day'))###limite à 2008 mem si donnees jusque 2012
    DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2012-12-31",by='day'))
    ###DATE_ERAI= atoms(timeSequence(from="1979-01-01",to="2003-12-31",by='day'))
    idxdates=which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
    DATE_ERAI = DATE_ERAI[idxdates,]
    
    ##DATE_ERAI2= atoms(timeSequence(from="2004-01-01",to="2012-12-31",by='day'))
    #idxdates2=which(DATE_ERAI2['m']==MON[seas,1] | DATE_ERAI2['m']==MON[seas,2] | DATE_ERAI2['m']==MON[seas,3])
    #DATE_ERAI2 = DATE_ERAI2[idxdates,]
    
    
    #datafile=paste(datadir,"ERAI_1.125/VALUE/",SEAS[seas],"_PCS_ERAI_79_12_EU_CORDEX.RData",sep="")
    #load(datafile)
    load(paste("~/Documents/LSCE/RDATA/",SEAS[seas],"_PCS_ERAI_79_12_SAM_CORDEX.RData",sep=""))
    
    
    ############################################
    ###########proj##########################
    ############################################
    
    
    ############################################
    if(experiment == "kfold"){
      
      
      #    if(ifold==1){idx_proj_dates=which(DATE_ERAI['Y']<=1984)}
      #    if(ifold==2){idx_proj_dates=which((DATE_ERAI['Y']>=1985) & (DATE_ERAI['Y']<=1990))}
      #    if(ifold==3){idx_proj_dates=which((DATE_ERAI['Y']>=1991) & (DATE_ERAI['Y']<=1996))}
      #    if(ifold==4){idx_proj_dates=which((DATE_ERAI['Y']>=1997) & (DATE_ERAI['Y']<=2002))}
      #    if(ifold==5){idx_proj_dates=which(DATE_ERAI['Y']>=2003)}
      
      
      
      if(ifold==1){
        #      idx_control_dates=which( (DATE_ERAI['Y']>=1984) & (DATE_ERAI['Y']<=2003) ) # 1984 - 2003
        #      idx_proj_dates=which(DATE_ERAI['Y']<=1983) # 1979 - 1983
      }
      if(ifold==2){
        #      idx_control_dates=which((DATE_ERAI['Y']<=1983) | ((DATE_ERAI['Y']>=1989) & (DATE_ERAI['Y']<=2003) )) # 1979-1983 & 1989-2003
        #      idx_proj_dates=which((DATE_ERAI['Y']>=1984) & (DATE_ERAI['Y']<=1988)) # 1984-1988
      }
      if(ifold==3){
        #      idx_control_dates=which((DATE_ERAI['Y']<=1988) | ((DATE_ERAI['Y']>=1994) & (DATE_ERAI['Y']<=2003) )) # 1979-1988 & 1994-2003
        #      idx_proj_dates=which((DATE_ERAI['Y']>=1989) & (DATE_ERAI['Y']<=1993)) # 1989-1993
      }
      if(ifold==4){
        #      idx_control_dates=which((DATE_ERAI['Y']<=1993) | ((DATE_ERAI['Y']>=1999) & (DATE_ERAI['Y']<=2003)) ) # 1979-1993 & 1999-2003
        #      idx_proj_dates=which((DATE_ERAI['Y']>=1994) & (DATE_ERAI['Y']<=1998)) # 1994-1998
      }
      if(ifold==5){
        #      idx_control_dates=which(DATE_ERAI['Y']<=1998) # 1979 - 1998
        #      idx_proj_dates=which((DATE_ERAI['Y']>=1999) & (DATE_ERAI['Y']<=2003)) # 1999-2003
      }
      
      # Common to all ifolds
      idx_proj2_dates=which((DATE_ERAI['Y']>=2004) & (DATE_ERAI['Y']<=2012)) # 2004 - 2012
      
      
      #    PCS=PCS[idx_proj_dates,] # Pour PR
      PCS=PCS[idx_proj2_dates,] # Pour PR
      ##    PCS=PCS[idx_proj_dates, 3:12] # Pour T
      NB_pred = dim(PCS)[2]
      
      
      # Nom du fichier sauvegardé contenant les estimations
      #    filename=paste('/home/users/mvrac/R/ESD/Output/SWG/SWG_ERAI_ESD_PR_',season[seas],'_fold_',ifold,'.RData',sep="")
      #paste('SWG_ERAI_VALUE_ECA_PR_',season[seas],'_0mm_fold_',ifold,'.RData',sep="")
      filename=paste('~/Documents/LSCE/SWG/SWG_ERAI_ESD_PR_',season[seas],'_fold_',ifold,'.RData',sep="")
      load(filename)
      
      #    tempP1 = array(NaN,c(length(idx_proj_dates),length(fit_stations_ro)))
      #    tempshape = array(NaN,c(length(idx_proj_dates),length(fit_stations_ro)))
      #    temprate = array(NaN,c(length(idx_proj_dates),length(fit_stations_ro)))
      
      tempP12 = array(NaN,c(length(idx_proj2_dates),length(fit_stations_ro)))
      tempshape2 = array(NaN,c(length(idx_proj2_dates),length(fit_stations_ro)))
      temprate2 = array(NaN,c(length(idx_proj2_dates),length(fit_stations_ro)))
      
      
      slot=array(NA,length(fit_stations_rf))
      
      for (i in 1:length(fit_stations_rf)){
        slot[i]=length(slotNames(fit_stations_rf[[i]]))
      }
      idx=which(slot>0)
      
      for (k in idx){
        
        ## occurrence
        if(statio_pix_ro[[k]]=='statio'){
          par_ro=fit_stations_ro[[k]]
          temp=par_ro[1]
          #        tempP1[,k] = replicate(length(idx_proj_dates),temp)
          tempP12[,k] = replicate(length(idx_proj2_dates),temp)
        }
        
        if(statio_pix_ro[[k]]!='statio'){
          par_ro = coef(fit_stations_ro[[k]],matrix=T)
          temp = par_ro[1] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_ro[2:(NB_pred+1)])
          #        tempP1[,k] = exp(temp)/(1+exp(temp))
          tempP12[,k] = exp(temp)/(1+exp(temp))
        }
        
        
        #intensity
        par_rf=coef(fit_stations_rf[[k]],matrix=T)
        if(statio_pix_rf[[k]]=='statio'){
          temp2 =par_rf
          #        tempshape[,k] =replicate(length(idx_proj_dates),exp (temp2[2]))
          #        temprate[,k] = replicate(length(idx_proj_dates),exp (temp2[1]))
          tempshape2[,k] =replicate(length(idx_proj2_dates),exp (temp2[2]))
          temprate2[,k] = replicate(length(idx_proj2_dates),exp (temp2[1]))
        }
        
        if(statio_pix_rf[[k]]!='statio'){
          #        tempshape[,k] = exp(par_rf[1,2] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_rf[2:(NB_pred+1),2]))
          #        temprate[,k] = exp(par_rf[1,1] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_rf[2:(NB_pred+1),1]))
          tempshape2[,k] = exp(par_rf[1,2] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_rf[2:(NB_pred+1),2]))
          temprate2[,k] = exp(par_rf[1,1] + as.matrix(PCS[,1:(NB_pred)])%*%as.matrix(par_rf[2:(NB_pred+1),1]))
        }
        
      } # end for k
      
      #    P1[idxdates[idx_proj_dates],] = tempP1
      #    shape[idxdates[idx_proj_dates],] = tempshape
      #    rate[idxdates[idx_proj_dates],] = temprate
      P12[ifold, (idxdates[idx_proj2_dates] - 9131),] = tempP12
      shape2[ifold, (idxdates[idx_proj2_dates] - 9131),] = tempshape2
      rate2[ifold, (idxdates[idx_proj2_dates] - 9131),] = temprate2
      
      ############################################
    } # end if xp kfold
    
    
  } # end for seas
  
} # end for ifold

############################################


#save(P12,shape2,rate2,DATE_ERAI2, file = paste('/home/users/mvrac/R/ESD/Output/SWG/FieldParameters_SWG_ERAI_ESD_Claris_PR_',experiment,'_2004_2012.RData',sep=""))
save(P12,shape2,rate2,DATE_ERAI2, file = paste('~/Documents/LSCE/SWG/FieldParameters_SWG_ERAI_ESD_Claris_PR_',experiment,'_2004_2012.RData',sep=""))


##################
#####sampling########
##################

library(timeDate)
DATE = atoms(timeSequence(from="2004-01-01",to="2012-12-31",by='day'))

#datadir="/home/estimr1/pvait/RDATA/"

##################
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
#######################


for(i in 1:1){
  
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  cat('Run',NUM,'\n')
  
  for(ifold in 1:5){
    
    #idx=which(P1[1,]!='NaN')
    idx=which(P12[ifold,1,]!='NaN')
    
    cat(" ifold")
    Sample_temp = Sample_binomocc_gammaity(P12[ifold,,idx],shape2[ifold,,idx],rate2[ifold,,idx])
    Sample=array(NaN,dim=dim(P12))
    Sample[ifold,,idx]=Sample_temp
    
    cat("\n")
    
    Sample=round(Sample,digits=2)
    
    filesimu = paste('~/Documents/LSCE/SWG/SIMU_SWG_PR_ERAI_ESD_Claris_PR_',experiment,'_2004_2012_ifold_',ifold,'_run_',NUM,'.RData',sep="")
    save(Sample,DATE,file = filesimu)
    
    
  } # end for ifold
  
  
}

## q(save='no')

load("/Users/shengchen/Documents/LSCE/SWG/SIMU_SWG_PR_ERAI_ESD_Claris_PR_kfold_2004_2012_ifold_1_run_001.RData")
sum(is.na(Sample))/sum(!is.na(Sample))
Sample_1 = Sample[1,,] ;sum(is.na(Sample_1))
range(Sample_1)

load("/Users/shengchen/Documents/LSCE/SWG/SIMU_SWG_PR_ERAI_ESD_Claris_PR_kfold_2004_2012_ifold_2_run_001.RData")
sum(is.na(Sample))/sum(!is.na(Sample))
Sample_2 = Sample[2,,] ;sum(is.na(Sample_2))
range(Sample_2)

load("/Users/shengchen/Documents/LSCE/SWG/SIMU_SWG_PR_ERAI_ESD_Claris_PR_kfold_2004_2012_ifold_3_run_001.RData")
sum(is.na(Sample))/sum(!is.na(Sample))
Sample_3 = Sample[3,,] ;sum(is.na(Sample_3))
range(Sample_3)

load("/Users/shengchen/Documents/LSCE/SWG/SIMU_SWG_PR_ERAI_ESD_Claris_PR_kfold_2004_2012_ifold_4_run_001.RData")
sum(is.na(Sample))/sum(!is.na(Sample))
Sample_4 = Sample[4,,] ;sum(is.na(Sample_4))
range(Sample_4)

load("/Users/shengchen/Documents/LSCE/SWG/SIMU_SWG_PR_ERAI_ESD_Claris_PR_kfold_2004_2012_ifold_5_run_001.RData")
sum(is.na(Sample))/sum(!is.na(Sample))
Sample_5 = Sample[5,,] ;sum(is.na(Sample_5))
range(Sample_5)
