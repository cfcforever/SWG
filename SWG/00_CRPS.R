# load package ------------------------------------------------------------
# Continuous Rank Probability Score (CRPS)
library(scoringRules)

SEAS   = c('DJF','MAM','JJA','SON')
season = c('Winter','Spring','Summer','Fall')
MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)


# CRPS  ----------------------------------------------------------
city.names = c("Paris", "Madrid", "Stockholm")
city = city.names[3]

load(file = paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
D = as.Date(0:(nrow(DATE_OBS)-1), origin = "1979-01-01")
idxdates = which(DATE_OBS["Y"]>=1999 & DATE_OBS["Y"]<=2017)
DATE = D[idxdates]

# obs: 1999 - 2017
obs = as.data.frame(matrix(NA, nrow = length(idxdates), ncol = 2))
colnames(obs) = c("tmean", "date")
obs[,1] = tmean[idxdates,1]
obs[,2] = D[idxdates]
obs$type = "obs"

load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/t2m_mean_sd_1999_2017_nat_1.RData"))
mean1 = mean; sd1 = sd
load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/t2m_mean_sd_1999_2017_nat_0.RData"))
mean0 = mean; sd0 = sd

tCRPS = as.data.frame(matrix(NA, nrow = length(idxdates), ncol = 2))
colnames(tCRPS) = c("stat", "cond")
for (i in 1:length(idxdates)){
  x = obs[i,1]
  tCRPS[i,] = c(crps(y = x, family = "normal", mean = mean0[i], sd = sd0[i]), 
                crps(y = x, family = "normal", mean = mean1[i], sd = sd1[i]))
}
tCRPS$date = D[idxdates]

s = rep(NaN, length(idxdates))
for (seas in 1:4){
  idx = which(month(tCRPS$date)==MON[seas,1] | month(tCRPS$date)==MON[seas,2] | month(tCRPS$date)==MON[seas,3])
  s[idx] = season[seas]
}
tCRPS$seas = s
tCRPS$seas <- factor(tCRPS$seas, levels = )

vmax = max(tCRPS[,1:2])
p <- ggplot(data = tCRPS, aes(x = stat, y = cond)) + theme_bw() +
  geom_point() + 
  geom_abline(intercept =0 , slope = 1, col = "red") +
  facet_wrap(~seas, nrow = 2) +
  coord_equal(ratio=1) + 
  xlim(0, vmax) + ylim(0,vmax) +
  labs(x = "non-conditional model", y = "conditional slp-SD model", 
       title = paste0("Comparison of daily CRPS for ", city))
print(p)
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/Image/CRPS_", city, ".pdf"), width = 7, height = 7)

# CRPS bad version --------------------------------------------------------------------
city.names = c("Paris", "Madrid", "Stockholm")
city = city.names[2]

list.city = vector("list", 3)
names(list.city) = city.names
for (city in city.names){
  load(file = paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))
  DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  D = as.Date(0:(nrow(DATE_OBS)-1), origin = "1979-01-01")
  idxdates = which(DATE_OBS["Y"]>=1999 & DATE_OBS["Y"]<=2017)
  
  # obs: 1999 - 2017
  obs = as.data.frame(matrix(NA, nrow = length(idxdates), ncol = 2))
  colnames(obs) = c("tmean", "date")
  obs[,1] = tmean[idxdates,1]
  obs[,2] = D[idxdates]
  obs$type = "obs"
  
  # sim_conditional: 1999 - 2017
  sim.cond.list = vector("list", 100)
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_1_run_", NUM, ".RData"))
    sim_cond = as.data.frame(matrix(NA, nrow = length(idxdates), ncol = 2))
    colnames(sim_cond) = c("tmean", "date")
    sim_cond[,1] = Sample[,1]
    sim_cond[,2] = D[idxdates]
    sim_cond$type = "contional-slp_SD"
    sim.cond.list[[i]] = sim_cond
  }
  sim_cond = do.call(rbind.data.frame, sim.cond.list)
  
  # sim_stationary: 1999 - 2017
  sim.stat.list = vector("list", 100)
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_0_run_", NUM, ".RData"))
    sim_stat = as.data.frame(matrix(NA, nrow = length(idxdates), ncol = 2))
    colnames(sim_stat) = c("tmean", "date")
    sim_stat[,1] = Sample[,1]
    sim_stat[,2] = D[idxdates]
    sim_stat$type = "non-conditional"
    sim.stat.list[[i]] = sim_stat
  }
  sim_stat = do.call(rbind.data.frame, sim.stat.list)
  
  dat = rbind(obs, sim_cond, sim_stat)
  # save(dat, file = "/Volumes/Data-ExFAT/CRPS_dat.RData")
  
  tCRPS = as.data.frame(matrix(NA, nrow = length(idxdates), ncol = 3))
  colnames(tCRPS) = c("non-conditional", "contional-slp_SD", "year")
  for (i in 1:length(idxdates)){
    dat0 <- dat[dat$date==D[idxdates[i]] & dat$type=="non-conditional", 1]
    dat1 <- dat[dat$date==D[idxdates[i]] & dat$type=="contional-slp_SD", 1]
    x = obs[i,1]
    tCRPS[i,] = c(crps_sample(y = x, dat = dat0),
                  crps_sample(y = x, dat = dat1),
                  D[idxdates[i]])
  }
}
plot(tCRPS$`non-conditional`, tCRPS$`contional-slp_SD`)
abline(0,1,col="red")

p <- ggplot(tCRPS, aes(x=year,y=value,col=variable)) + theme_bw() +
  geom_line() + 
  facet_wrap(~city, nrow = 3) +
  ylab("CRPS")
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/Image/CRPS_mean_1999-2017.pdf"), width = 7, height = 7)


# CRPS max ----------------------------------------------------------------
list.city = vector("list", 3)
names(list.city) = city.names
for (city in city.names){
  load(file = paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))
  DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  D = as.Date(0:(nrow(DATE_OBS)-1), origin = "1979-01-01")
  
  # tam_obs: 1999 - 2017
  tam = as.data.frame(matrix(NA, nrow = length(unique(DATE_OBS$Y)), ncol = 2))
  colnames(tam) = c("tam", "date")
  rownames(tam) = unique(DATE_OBS$Y)
  for (y in 1979:2017){
    Dy = D[DATE_OBS$Y==y][which(tmean[DATE_OBS$Y==y,] == max(tmean[DATE_OBS$Y==y,]))]
    tam[paste(y),] = c(max(tmean[DATE_OBS$Y==y,]), as.character(Dy))
  }
  
  tam_obs = as.data.frame(matrix(NA, nrow = length(unique(DATE_OBS$Y)), ncol = 2))
  colnames(tam_obs) = c("tam", "year")
  rownames(tam_obs) = unique(DATE_OBS$Y)
  tam_obs[,1] = as.numeric(tam[as.character(unique(DATE_OBS$Y)),1])
  tam_obs[,2] = unique(DATE_OBS$Y)
  tam_obs$type = "obs"
  
  tam_obs = tam_obs[tam_obs$year>=1999 & tam_obs$year<=2017,]
  
  # tam_sim_conditional: 1999 - 2017
  tam.list = vector("list", 100)
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_1_run_", NUM, ".RData"))
    tam_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
    colnames(tam_sim) = c("tam", "year")
    rownames(tam_sim) = unique(DATE$Y)
    for (y in unique(DATE$Y)){
      tam_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
    }
    tam.list[[i]] = tam_sim
  }
  tam_sim_1 = do.call(rbind.data.frame, tam.list)
  tam_sim_1$type = "sim"
  
  # tam_sim_stationary: 1999 - 2017
  tam.list = vector("list", 100)
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_0_run_", NUM, ".RData"))
    tam_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
    colnames(tam_sim) = c("tam", "year")
    rownames(tam_sim) = unique(DATE$Y)
    for (y in unique(DATE$Y)){
      tam_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
    }
    tam.list[[i]] = tam_sim
  }
  tam_sim_0 = do.call(rbind.data.frame, tam.list)
  tam_sim_0$type = "sim"
  
  tCRPS = as.data.frame(matrix(NA, nrow = length(1999:2017), ncol = 3))
  colnames(tCRPS) = c("non-conditional", "contional-slp_SD", "year")
  rownames(tCRPS) = c(1999:2017)
  for (y in 1999:2017){
    dat0 <- tam_sim_0[tam_sim_0$year==y,1]
    dat1 <- tam_sim_1[tam_sim_1$year==y,1]
    x <- tam_obs[tam_obs$year==y,1]
    tCRPS[rownames(tCRPS)==y,] = c(crps_sample(y = x, dat = dat0),
                                   crps_sample(y = x, dat = dat1),
                                   y)
    # tCRPS[rownames(tCRPS)==y,] = c(crps(obs = x, pred = dat0),
    #                                crps(obs = x, pred = dat1),
    #                                y)
  }
  tCRPS = melt(tCRPS, id.vars = "year")
  tCRPS$city = city
  
  list.city[[city]] = tCRPS
}
tCRPS = do.call(rbind.data.frame, list.city)

p <- ggplot(tCRPS, aes(x=year,y=value,col=variable)) + theme_bw() +
  geom_line() + 
  facet_wrap(~city, nrow = 3) +
  ylab("CRPS")
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/Image/CRPS_max_1999-2017.pdf"), width = 7, height = 7)


