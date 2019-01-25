library(ggplot2)
library(forecast)
library(tseries)

#### load t2m data - obs ####
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData") 


#### estimation ####
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
  
  idx_control_dates = which((DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']<=2004))
  
  # cherche les indices de jour dans DATE_ECA qui correspondent aux dates ERAI pour la période de calibration
  idx_dates=array(0,dim=length(idx_control_dates))
  for (i in 1:length(idx_control_dates)){
    idx_dates[i]=which(((DATE_OBS['Y'])[,1]==(DATE_ERAI['Y'])[idx_control_dates[i],1])
                       & ((DATE_OBS['m'])[,1]==(DATE_ERAI['m'])[idx_control_dates[i],1])
                       & ((DATE_OBS['d'])[,1]==(DATE_ERAI['d'])[idx_control_dates[i],1]))
  } # end for i
  
  temp=dat[idx_dates,1]
  
  model <- auto.arima(temp, start.p = 0, max.p = 10, d = 0, start.q = 0, max.q = 10)
  
  save(model, file = paste0("~/Documents/LSCE/SWG/t2m/case_0/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val.RData"))
} 


#### simulation ####
# 2005-2017
DATE_proj = atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day'))
LENGTH_proj = dim(atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day')))[1]

DATE = atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day'))

for(i in 1:100){
  
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  cat('Run', NUM, '\n')
  
  Sample = array(NaN,c(LENGTH_proj,1))
  for (seas in 1:4){
    
    idxdates = which(DATE_proj['m']==MON[seas,1] | DATE_proj['m']==MON[seas,2] | DATE_proj['m']==MON[seas,3])
    
    load(paste0("~/Documents/LSCE/SWG/t2m/case_0/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val.RData"))
    
    Sample[idxdates, ] = arima.sim(n=length(idxdates), list(ar = model$coef[1:model$arma[1]], ma = model$coef[(model$arma[1]+1):(model$arma[1]+model$arma[2])]),sd = sqrt(model$sigma2)) + tail(model$coef,1)
  }
  
  Sample=round(Sample,digits=2)
  
  filesimu = paste('~/Documents/LSCE/SWG/t2m/case_0/SIMU/SIMU_SWG_tmean', '_run_', NUM, '.RData', sep="")
  save(Sample,DATE,file = filesimu)
  
}


Obs = tmean[which(rownames(tmean)=="2005-01-01"):which(rownames(tmean)=="2017-12-31"),]
range(Obs)
# fun_qqplot(Obs, Sample, main = "tmean - 100 simulations - Stationary")

Sample_all = c()
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(paste0("~/Documents/LSCE/SWG/t2m/case_0/SIMU/SIMU_SWG_tmean_run_", NUM, ".RData"))
  Sample_all = c(Sample_all, Sample)
}
range(Sample_all)

fun_qqplot(Obs, Sample_all, main = "tmean - 100 simulations - Stationary")
dev.print(pdf, file = "~/Documents/LSCE/SWG/t2m/Images/qqplot_Stationary.pdf", width = 7, height = 7)


## ACF
# load tmean
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")   

for (seas in 4:4){
  # Obs as tmean - daily mean temperature from 2005 to 2017 in the zone [48N, 50N]x[01E, 04E]
  Obs = tmean[which(rownames(tmean)=="2005-01-01"):which(rownames(tmean)=="2017-12-31"),]
  DATE_Obs = atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day'))
  
  MON    = c(12, 1:11)
  MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
  SEAS   = c('DJF', 'MAM', 'JJA', 'SON')
  Obs = Obs[DATE_Obs$m==MON[seas,1] | DATE_Obs$m==MON[seas,2] | DATE_Obs$m==MON[seas,3]]
  
  for (k in 1:100){
    # NUM: number with 3 digits
    # k = 2
    NUM=cbind(formatC(k, digits = 0, width = 3, format = "f", flag = "0"))
    
    load(paste0("~/Documents/LSCE/SWG/t2m/case_", 0, "/SIMU/SIMU_tmean_", SEAS[seas], "_run_", NUM, ".RData"))
    Sample = Sample_seas
    
    # set number of lag
    nlag = 31
    
    # data_acf contains acf_obs and acf_sim
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
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/t2m/Images/AIC/acf_", SEAS[seas],"_case_0_run_", NUM, ".pdf"), width = 7, height = 7)
  }
}


#### draft ####
adf.test(temp, alternative = "stationary")

acf(temp)
pacf(temp)


model <- auto.arima(temp, start.p = 0, max.p = 10, d = 0, start.q = 0, max.q = 10)
model <- arima(temp, order = c(3,0,1))

y = arima.sim(n=1000, list(ar = model$coef[1:model$arma[1]], ma = model$coef[(model$arma[1]+1):(model$arma[1]+model$arma[2])]),sd = sqrt(model$sigma2)) + tail(model$coef,1)
range(y);range(temp)

auto.arima(temp, start.p = 0, max.p = 10, d = 0, start.q = 0, max.q = 10)
auto.arima(y, start.p = 0, max.p = 10, d = 0, start.q = 0, max.q = 10)
