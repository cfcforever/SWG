source("function/load_packages.R")
####

# common variables --------------------------------------------------------
city.names = c("Paris", "Madrid", "Stockholm")
city = city.names[3]
load(paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
DATE = as.Date(0:(nrow(DATE_OBS)-1), origin = "1979-01-01")
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

# Preliminary -------------------------------------------------------------
tam = as.data.frame(matrix(NA, nrow = length(unique(DATE_OBS$Y)), ncol = 1))
colnames(tam) = "tam"
date = rep(NA, length(unique(DATE_OBS$Y)))
for (y in 1979:2017){
  i = which(y==1979:2017)
  date[i] = as.character(DATE[DATE_OBS$Y==y][which(tmean[DATE_OBS$Y==y,] == max(tmean[DATE_OBS$Y==y,]))])
  tam[i,] = max(tmean[DATE_OBS$Y==y,])
}
tam$date = date
tam$year = 1979:2017

seas = 3
load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/predictor/predictor_", SEAS[seas], "_1.RData"))
idxdates = which(DATE_OBS['m']==MON[seas,1] | DATE_OBS['m']==MON[seas,2] | DATE_OBS['m']==MON[seas,3])
tam$dim = NA
tam$theta = NA
for (i in 1:39){
  tam[i,c("dim","theta")] = PCS[which(DATE[idxdates] == tam$date[i]),]
}

p1 <- ggplot(tam, aes(x = year, y = tam)) +
  geom_line() + geom_smooth(se = FALSE, method = "loess")
p2 <- ggplot(tam, aes(x = dim, y = tam)) +
  geom_point() + geom_smooth(se = FALSE, method = "loess")
p3 <- ggplot(tam, aes(x = theta, y = tam)) +
  geom_point() + geom_smooth(se = FALSE, method = "loess")
plot_grid(p1, p2, p3,  ncol = 1)   
# dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/gev/Image/temperature_annual_max_1979-2017_", city, ".pdf"), width = 7, height = 7)


# Estimation --------------------------------------------------------------
pc_name=names(PCS)
fmlatt=as.formula(paste("TT ~ ",paste(pc_name, collapse= "+")))
fit_stations_tt <- vector("list", 1)
PCS = tam[1:20,c("dim","theta")]
names(PCS) = c("PC1", "PC2")

temp = tam[1:20, "tam"]
tt = data.frame(TT=temp,PCS)
# fit_stations_tt[[(1)]] = try(
#   vglm(fmlatt,
#        gev,
#        data=tt,
#        maxit=1000,
#        # x.arg= FALSE,
#        # y.arg= FALSE,
#        # qr.arg = FALSE,
#        trace = T
#        # coefstart
#   ), silent=F)
# save(fit_stations_tt, file = "~/Documents/LSCE/SWG/slp_SD_diagnosis/gev/fit_stations_tt.RData")

fit_tt = fevd(TT, data = tt, location.fun = ~PC1+PC2, scale.fun = ~PC1+PC2)


# Simulation --------------------------------------------------------------
seas = 3
load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/predictor/predictor_", SEAS[seas], "_1.RData"))
SD = PCS; colnames(SD) = c("dim", "theta")
idxdates = which(DATE_OBS['m']==MON[seas,1] | DATE_OBS['m']==MON[seas,2] | DATE_OBS['m']==MON[seas,3])
DATE_summer = DATE_OBS[idxdates,]

idxselect = which(DATE_summer$Y>=1999 &  DATE_summer$Y<=2017)
SD = SD[idxselect,]
DATE_summer = DATE_summer[idxselect,]

SD_mean = as.data.frame(matrix(NA, nrow = 19, ncol = 2))
colnames(SD_mean) = colnames(SD)
for (y in 1999:2017){
  i = which(1999:2017==y)
  SD_mean[i,] = colMeans(SD[DATE_summer$Y==y,])
}

par = fit_tt$results$par
parnames = fit_tt$parnames

mu = array(NaN, c(19,1))
sigma = array(NaN, c(19,1))
shape = array(NaN, c(19,1))
for (k in 1:19){
  mu[k]    = par[1] + par[2]*SD_mean$dim[k] + par[3]*SD_mean$theta[k]
  sigma[k] = par[4] + par[5]*SD_mean$dim[k] + par[6]*SD_mean$theta[k]
  shape[k] = par[7]
}
# save(mu, sigma, shape, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/gev/mu_sigma_shape_", city, ".RData"))

load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/gev/mu_sigma_shape_", city, ".RData"))
## make 100 runs
for(i in 1:100){
  
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  cat('Run', NUM, '\n')
  
  Sample = array(NaN, c(19,1))
  for (k in 1:19){
    Sample[k,] = rgev(1, loc = mu[k], scale = sigma[k], shape = shape[k])
  }
  Sample=round(Sample,digits=2)
  
  filesimu = paste0('~/Documents/LSCE/SWG/slp_SD_diagnosis/gev/SIMU/', city, '/SIMU_SWG_tmax_1999_2017_', NUM, '.RData')
  save(Sample,file = filesimu)
}


# Boxplot - conditional model - stationary model - obs - GEV model --------
city.names = c("Paris", "Madrid", "Stockholm")
city = city.names[1]
load(file = paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
tmean = tmean[DATE_OBS$Y>=1999 &  DATE_OBS$Y<=2017,]

DATE_SIM = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))
DATE = as.Date(0:(nrow(DATE_SIM)-1), origin = "1999-01-01")

## tmax_obs
tmax_obs = as.data.frame(matrix(NA, nrow = length(unique(DATE_SIM$Y)), ncol = 2))
colnames(tmax_obs) = c("tmax", "year")
rownames(tmax_obs) = unique(DATE_SIM$Y)
for (y in 1999:2017){
  tmax_obs[paste(y),] = c(max(tmean[DATE_SIM$Y==y]), y)
}
tmax_obs$type = "obs"

## tmax_cond
tmax.list = vector("list", 100)
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_1_run_", NUM, ".RData"))
  tmax_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
  colnames(tmax_sim) = c("tmax", "year")
  rownames(tmax_sim) = unique(DATE$Y)
  for (y in unique(DATE$Y)){
    tmax_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
  }
  tmax.list[[i]] = tmax_sim
}
tmax_cond = do.call(rbind.data.frame, tmax.list)
tmax_cond$type = "cond"

## tmax_cond_2
tmax.list = vector("list", 100)
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(file = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_1_run_", NUM, ".RData"))
  tmax_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
  colnames(tmax_sim) = c("tmax", "year")
  rownames(tmax_sim) = unique(DATE$Y)
  for (y in unique(DATE$Y)){
    tmax_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
  }
  tmax.list[[i]] = tmax_sim
}
tmax_cond2 = do.call(rbind.data.frame, tmax.list)
tmax_cond2$type = "cond2"

## tmax_stat
tmax.list = vector("list", 100)
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_0_run_", NUM, ".RData"))
  tmax_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
  colnames(tmax_sim) = c("tmax", "year")
  rownames(tmax_sim) = unique(DATE$Y)
  for (y in unique(DATE$Y)){
    tmax_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
  }
  tmax.list[[i]] = tmax_sim
}
tmax_stat = do.call(rbind.data.frame, tmax.list)
tmax_stat$type = "stat"

## tmax_gev
tmax.list = vector("list", 100)
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/gev/SIMU/", city, "/SIMU_SWG_tmax_1999_2017_", NUM, ".RData"))
  tmax_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
  colnames(tmax_sim) = c("tmax", "year")
  rownames(tmax_sim) = unique(DATE$Y)
  tmax_sim$tmax = Sample[,1]
  tmax_sim$year = 1999:2017
  tmax.list[[i]] = tmax_sim
}
tmax_gev = do.call(rbind.data.frame, tmax.list)
tmax_gev$type = "gev"

## combine obs + stat + gev + cond(slp_SD) + cond2(slp_z500_SD)
tmax_data = rbind(tmax_obs, tmax_stat, tmax_gev, tmax_cond, tmax_cond2)

output = paste0("boxplot_tmax_annual_1999_2017_", city)
p <- ggplot(tmax_data) + theme_bw() +
  geom_boxplot(data = tmax_data[tmax_data$type!="obs",], mapping = aes(x = factor(year), y = tmax, col = type)) +
  scale_color_manual(values = matlab.like2(4)) +
  geom_point(data = tmax_data[tmax_data$type=="obs",], mapping = aes(x = factor(year), y = tmax), col = "black") + 
  # scale_color_manual(values = c(cond="blue", gev="red", stat="cyan", obs="black", cond2="green")) +
  labs(x = "year", col = "model", title = output)
print(p)
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/Image/", output, ".pdf"), width = 9, height = 5)


## CRPS
tCRPS = as.data.frame(array(NaN, c(19,4)))
colnames(tCRPS) = c("cond", "stat", "gev", "cond2")
ys = 1999:2017
for (t in 1:19){
  y = ys[t]
  tCRPS[t,1] = crps_sample(y = tmax_obs[t,1], dat = tmax_cond[tmax_cond$year==y, 1])
  tCRPS[t,2] = crps_sample(y = tmax_obs[t,1], dat = tmax_stat[tmax_stat$year==y, 1])
  tCRPS[t,3]  = crps_sample(y = tmax_obs[t,1], dat = tmax_gev[tmax_gev$year==y, 1])
  tCRPS[t,4] = crps_sample(y = tmax_obs[t,1], dat = tmax_cond2[tmax_cond2$year==y, 1])
}
tCRPS$year = ys

melt_CRPS = melt(tCRPS, id.vars = "year")
colnames(melt_CRPS) = c("year", "model", "CRPS")
melt_CRPS$model <- factor(melt_CRPS$model, levels = c("cond", "gev", "stat", "cond2"))

p <- ggplot(melt_CRPS) + theme_bw() + 
  geom_line(aes(x = year, y = CRPS, col = model)) +
  scale_color_manual(values = matlab.like2(4))
print(p)
save_pdf(filename = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/Image/tmax_CRPS_", city))
