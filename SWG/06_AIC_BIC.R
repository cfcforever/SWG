##-----------------------------------------------
# change t2m  to tp
# change geom to slp
##-----------------------------------------------

library(ggplot2)
library(reshape2)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')
casename = c("PCS", "SysDyn", "PCS_SysDyn")

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

for (case in 3:3){
  for (seas in 1:4){
    
    load(paste0("~/Documents/LSCE/SWG/t2m/case_", case, "/geop_", SEAS[seas], "_", casename[case], ".RData"))
    
    n = ncol(PCS); n
    
    AIC = rep(NA, n)
    BIC = rep(NA, n)
    
    for (k in 1:n){
      
      # DATE_ERAI as the same as DATE_OBS
      DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
      
      # choose for the right months for season 
      idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
      
      # DATE_ERAI only in season
      DATE_ERAI = DATE_ERAI[idxdates,]
      
      # load tmean - this is main data
      load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
      dat = tmean
      
      # number of station - here only 1
      NB_STATION = 1
      
      # chosen period for estimation
      # 1979 - 2017: 39 years
      # 1979 - 2004: 26 years for estimation 
      # 2005 - 2017: 13 years for validation
      idx_control_dates = which((DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']<=2004))
      
      # PCS in estimation years
      load(paste0("~/Documents/LSCE/SWG/t2m/case_", case, "/geop_", SEAS[seas], "_", casename[case], ".RData"))
      PCS=as.data.frame(PCS[idx_control_dates,1:k])
      pc_name=paste0("PC", 1:ncol(PCS))[1:k]
      colnames(PCS) = pc_name
      
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
          fit_stations_tt[[(i)]] = try(vgam(fmlatt, uninormal(zero=NULL), data=tt,maxit=1000, x.arg= FALSE, y.arg= FALSE, qr.arg = FALSE),silent=T)
        }                          
        
      } # end for i NB_STATION
      
      AIC[k] = 2*(2+k*2) - 2*logLik(fit_stations_tt[[1]])
      BIC[k] = -2*logLik(fit_stations_tt[[1]]) + (2+k*2)*log(length(temp), base = exp(1))
      
    }
    
    data = cbind.data.frame(AIC, BIC)
    data$PC = c(1:n)
    data = melt(data, id.vars = "PC")
    
    p = ggplot(data)
    p = p + geom_point(aes(x=PC, y=value))
    p = p + geom_point(data = data[data$variable=="AIC",], aes(x=PC[which(value==min(value))], y=value[which(value==min(value))], col = "red"))
    p = p + geom_point(data = data[data$variable=="BIC",], aes(x=PC[which(value==min(value))], y=value[which(value==min(value))], col = "red"))
    p = p + facet_wrap(~variable, nrow = 2) + theme_bw() + theme(legend.position = "none")
    p = p + labs(title = paste0("AIC and BIC for t2m in ", season[seas], " (", casename[case], ")"), 
                 x = "Number of composants",
                 y = "")
    print(p)
    
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/t2m/Images/AIC_BIC_t2m_", season[seas], "_", casename[case], ".pdf"), width = 9, height = 7)
  
  }
}
