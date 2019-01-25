fun_estimation <- function(predictor, NUM, input.tmean, output.name, year.begin, year.end){
  
  library(timeDate)
  
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
    DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))
    
    # choose for the right months for season 
    idxdates = which(DATE_ERAI['m']==MON[seas,1] | DATE_ERAI['m']==MON[seas,2] | DATE_ERAI['m']==MON[seas,3])
    
    # DATE_ERAI only in season
    DATE_ERAI = DATE_ERAI[idxdates,]
    
    # /!\ load PCS of season - VERY IMPORTANT
    load(paste0(predictor, SEAS[seas], "_", NUM, ".RData"))
    colnames(PCS) = paste0("PC", 1:ncol(PCS))
    
    # this is main data
    load(input.tmean)
    dat = tmean
    
    # number of station - here only 1
    NB_STATION = 1
    
    idx_control_dates = which((DATE_ERAI['Y']>=year.begin) & (DATE_ERAI['Y']<=year.end))
    n = length(idx_control_dates)
    
    # PCS in estimation years
    PCS=PCS[idx_control_dates,]
    pc_name=names(PCS)
    npc = length(pc_name)
    
    # cherche les indices de jour dans DATE_ECA qui correspondent aux dates ERAI pour la période de calibration
    idx_dates=array(0,dim=length(idx_control_dates))
    for (i in 1:length(idx_control_dates)){
      idx_dates[i]=which(((DATE_OBS['Y'])[,1]==(DATE_ERAI['Y'])[idx_control_dates[i],1]) & 
                           ((DATE_OBS['m'])[,1]==(DATE_ERAI['m'])[idx_control_dates[i],1]) &
                           ((DATE_OBS['d'])[,1]==(DATE_ERAI['d'])[idx_control_dates[i],1]))
    } # end for i
    
    fmlatt=as.formula(paste("TT ~ ",paste(pc_name, collapse= "+")))
    fit_stations_tt <- vector("list", NB_STATION)
    
    i = 1
    for (i in 1:NB_STATION){
      
      temp=dat[idx_dates,i]
      
      if(length(which(is.na(temp)!=TRUE))==0){
        fit_stations_tt[[i]]=NaN
        output = paste0(output.name, season[seas], "_cross-val_", year.begin, "_", year.end, "_", NUM, ".RData")
        cat(output,'\n')
        cat("all NAs, no data", '\n')
        save(fit_stations_tt, file = output)
      }
      if(length(which(is.na(temp)!=TRUE))!=0){
        tt = data.frame(TT=temp,PCS)   # ncol = 1+20
        fit_stations_tt[[(i)]] = try(
          vgam(fmlatt,
               uninormal(zero=NULL),
               data=tt,
               maxit=1000,
               x.arg= FALSE,
               y.arg= FALSE,
               qr.arg = FALSE,
               trace = T
               # coefstart = as.numeric(coefstart)
               # coefstart = c(17, 1.1, 0.05, -0.003, 0.05, 0.004, -0.133, 0.02, -0.32, 0.02)
          ), silent=F)
        AIC = 2*(2+npc*2) - 2*logLik(fit_stations_tt[[1]])
        BIC = (2+npc*2)*log(length(temp), base = exp(1)) - 2*logLik(fit_stations_tt[[1]]) 
        
        output = paste0(output.name, season[seas], "_cross-val_", year.begin, "_", year.end, "_", NUM, ".RData")
        cat(output,'\n')
        save(fit_stations_tt, AIC, BIC, file = output)
      }                          
      
    } # end for i NB_STATION
  } # end for seas
}