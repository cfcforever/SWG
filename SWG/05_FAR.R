#### FAR for t2m from 1998 to 2017 ####
list_mean_sd_nat = vector("list", 3)
list_mean_sd_rea = vector("list", 3)
for (k in 1:3){
  load(paste0("~/Documents/LSCE/SWG/t2m/case_", k, "/tmean_mean_sd_1998_2017_nat.RData"))
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  load(paste0("~/Documents/LSCE/SWG/t2m/case_", k, "/tmean_mean_sd_1998_2017_rea.RData"))
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_rea[[k]] = m
}

# load("~/Documents/LSCE/SWG/t2m/case_2/tmean_mean_sd_1998_2017_nat.RData")
# mean_nat = mean
# sd_nat   = sd
# 
# load("~/Documents/LSCE/SWG/t2m/case_2/tmean_mean_sd_1998_2017_rea.RData")
# mean_rea = mean
# sd_rea   = sd


MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

mDATE = atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day'))
n = nrow(mDATE)

DATE = as.Date((0:(n-1)), origin = "1998-01-01")

case_name = c("PCS", "SysDyn", "PCS_SysDyn")

# thres = 10
for (year in 2003:2003){
  for (thres in 10:30){
    library(lubridate)
    
    list_prob_thres_nat = lapply(list_mean_sd_nat, FUN = function(m){
      pnorm(thres, mean = m[,"mean"], sd = m[,"sd"], lower.tail = F)
    })
    prob_thres_nat = do.call(rbind.data.frame, list_prob_thres_nat)
    colnames(prob_thres_nat) = "nat"
    
    list_prob_thres_rea = lapply(list_mean_sd_rea, FUN = function(m){
      pnorm(thres, mean = m[,"mean"], sd = m[,"sd"], lower.tail = F)
    })
    prob_thres_rea = do.call(rbind.data.frame, list_prob_thres_rea)
    colnames(prob_thres_rea) = "rea"
    
    prob_thres = cbind(prob_thres_nat, prob_thres_rea)
    prob_thres$date = rep(DATE, 3)
    # prob_thres$FAR = 1 - prob_thres_nat/prob_thres_rea
    # prob_thres$FAR = prob_thres$nat/prob_thres$rea
    prob_thres$FAR = (prob_thres$nat-prob_thres$rea)/(prob_thres$nat+prob_thres$rea)
    prob_thres$case = rep(case_name, each = length(DATE))
    
    # ## add case_0 stationary
    # load(paste0("~/Documents/LSCE/SWG/t2m/case_0/Sample_1998_2017_nat.RData"))
    # nat = apply(Sample.nat, MARGIN = 1, FUN = function(x){sum(x>=thres)/length(x)})
    # load(paste0("~/Documents/LSCE/SWG/t2m/case_0/Sample_1998_2017_rea.RData"))
    # rea = apply(Sample.rea, MARGIN = 1, FUN = function(x){sum(x>=thres)/length(x)})
    
    # statio = as.data.frame(cbind(nat, rea))
    # statio$date = DATE
    # statio$FAR = nat/rea
    # statio$case = "statio"
    
    # prob_thres = rbind.data.frame(prob_thres, statio)
    
    period_sel = intersect(which(year(prob_thres$date)==year) , which(month(prob_thres$date)==7|month(prob_thres$date)==6|month(prob_thres$date)==8))
    FAR_sel = prob_thres[period_sel,]
    FAR_sel$num = rep(1:(length(period_sel)/3), 3)
    
    FAR_sel$case <- factor(FAR_sel$case, levels = case_name)
    # FAR_sel$num = rep(1:(length(period_sel)/4), 4)
    
    # breaks = rep(NA,5)
    # labels = rep(NA,5)
    # for (k in 1:10){
    #   breaks[k] = which(year(FAR_sel$date)==unique(year(FAR_sel$date))[4*k-1])[1]
    #   labels[k] = as.character(FAR_sel$date[breaks[k]])
    # }
    
    breaks = rep(NA,3)
    labels = rep(NA,3)
    for (k in 1:3){
      breaks[k] = which(month(FAR_sel$date)==unique(month(FAR_sel$date))[k])[1]
      labels[k] = as.character(FAR_sel$date[breaks[k]])
    }
    
    p = ggplot(FAR_sel) 
    p = p + geom_line(aes(x=num, y=FAR, group=case, col=case)) + theme_bw()
    p = p + scale_x_continuous(breaks = breaks, labels = labels) + ylim(-1,1)
    p = p + labs(title = paste0("FAR_t2M_summer_over_", thres, "_degree"))
    print(p)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/t2m/Images/FAR/2PC/FAR_tmean_summer_", year, "_over_", thres, "_degree.pdf"), width = 7, height=7)
  }
}


case = 3
thres = 10
prob_threshold = array(NaN,c(n,2))
for(k in 1:n){
  prob_threshold[k,1] = pnorm(thres, mean = mean_nat[k], sd = sd_nat[k], lower.tail = F)
  prob_threshold[k,2] = pnorm(thres, mean = mean_rea[k], sd = sd_rea[k], lower.tail = F)
}
range(prob_threshold)

period_sel = which(mDATE$m==6 | mDATE$m==7 | mDATE$m==8)
prob_thres_sel = prob_threshold[period_sel, ]
range(prob_thres_sel)

FAR = 1 - prob_thres_sel[,1]/prob_thres_sel[,2]
range(FAR)

DATE_sel = DATE[period_sel]

dat_FAR = as.data.frame(matrix(NA, nrow = length(period_sel), ncol = 3))
colnames(dat_FAR) = c("date", "FAR", "num")
dat_FAR[,1] = DATE_sel
dat_FAR[,2] = FAR
dat_FAR[,3] = 1:length(period_sel)

breaks = rep(NA,10)
labels = rep(NA,10)
for (k in 1:10){
  breaks[k] = which(year(DATE_sel)==unique(year(DATE_sel))[2*k-1])[1]
  labels[k] = as.character(DATE_sel[breaks[k]])
}

p = ggplot(dat_FAR, aes(x=num, y=FAR)) + geom_point() + theme_bw()
p = p + scale_x_continuous(breaks = breaks, labels = labels)
p = p + labs(title = paste0("FAR_t2M_summer_over_", thres, "_degree_case_", name_case[case]))
p
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/t2m/Images/FAR/FAR_t2M_summer_over_", thres, "_degree_case_", case, ".pdf"), width = 9, height = 9)


#### precip : P1,shape,rate from 1979 to 2017 ####
DATE_proj = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
LENGTH_proj = dim(atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day')))[1]

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
  
  idx_proj_dates = which( (DATE_ERAI['Y']>=1979) & (DATE_ERAI['Y']<=2017) )    # 2005 - 2017
  
  load(paste("~/Documents/LSCE/SWG/tp/case_3/slp_",SEAS[seas],"_PCS_SysDyn.RData",sep=""))
  
  PCS=PCS[idx_proj_dates,]
  NB_pred=dim(PCS)[2]
  
  
  # Nom du fichier sauvegardé contenant les estimations
  filename=paste('~/Documents/LSCE/SWG/tp/case_3/SWG_ERAI_ESD_precip_',season[seas],'_cross-val.RData',sep="")
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
  
  td = which(DATE_OBS$Y==1979)[1]-1
  
  P1[(idxdates[idx_proj_dates] - td),] = tempP1       ; cat(range(tempP1),    "\n")
  shape[(idxdates[idx_proj_dates] - td),] = tempshape ; cat(range(tempshape), "\n")
  rate[(idxdates[idx_proj_dates] - td),] = temprate   ; cat(range(temprate),  "\n")
}
save(P1, shape, rate, file = "~/Documents/LSCE/SWG/tp/case_3/precip_P1_shape_rate_1979_2017.RData")


#### FAR for precip ####
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

mDATE = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
n = nrow(mDATE)

DATE = as.Date((0:(n-1)), origin = "1979-01-01")

case = 3
load(file = paste0("~/Documents/LSCE/SWG/tp/case_", case, "/precip_P1_shape_rate_1979_2017.RData"))
for (thres in 10:30){
  # thres = 30
  prob_threshold = array(NaN,c(n,1))
  for(k in 1:n){
    prob_threshold[k,1] = pgamma(thres, shape = shape[k], scale = 1/rate[k], lower.tail = F)*P1[k]
  }
  range(prob_threshold)
  # plot(x = DATE, y = prob_threshold, type = "l", main = "Probability of daily temperature exceeded to 10°C")
  
  period_sel = mDATE$m==MON[3,1] | mDATE$m==MON[3,2] | mDATE$m==MON[3,3]
  prob_thres_sel = prob_threshold[period_sel, ]
  range(prob_thres_sel)
  DATE_sel = DATE[period_sel]
  plot(x = DATE_sel, y = prob_thres_sel)
  
  period_sel_nat = mDATE>=1979 & mDATE <=1997
  prob_thres_sel_nat = prob_threshold[period_sel_nat]
  
  period_sel_rea = mDATE>=1998 & mDATE <=2016
  prob_thres_sel_rea = prob_threshold[period_sel_rea]
  
  FAR = 1 - prob_thres_sel_nat/prob_thres_sel_rea
  DATE_FAR = DATE[period_sel_rea]
  
  plot(x=DATE_FAR, y=FAR, "l", main = paste0("FAR of daily precipitation over ", thres, " mm"))
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/tp/Images/FAR/FAR_daily_precipitation_over_", thres, "_mm_case_", case, ".pdf"), width = 7, height = 7)
}
