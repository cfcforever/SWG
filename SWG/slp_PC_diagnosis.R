SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')
MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)

ens_pc = c(0,2:10,15,20)
npc = length(ens_pc)

library("factoextra")
for (seas in 1:4){
  load(file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[seas], "_pca.RData"))
  eig.val <- get_eigenvalue(pca)
  message(which(eig.val$cumulative.variance.percent>=90)[1])
  save(eig.val, file = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/eig_val_", SEAS[seas], ".RData"))
}

#### create PC with different numbers ####
ens_pc = c(0,2:10,15,20)
npc = length(ens_pc)
for (n in ens_pc){
  if (!dir.exists(paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", n, "PC"))){
    dir.create(paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", n, "PC"))
    message(paste0("The ", n, "PC", " dictionary created !"))
  }else{
    message(paste0("The ", n, "PC", " dictionary exists !"))
  }
  
  for (seas in 1:4){
    load(file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[seas], "_pca.RData"))
    PCS = as.data.frame(pca$x)[,1:n]
    
    # normalize PC 
    PCS = apply(PCS, MARGIN = 2, scale)
    colnames(PCS) = paste0("PC", 1:ncol(PCS))
    PCS = as.data.frame(PCS)
    save(PCS, file = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", n, "PC/predictor_",SEAS[seas],"_1.RData"))
  }
}

#### estimation for each case ####
## 1979 - 1998
for (n in ens_pc){
  for (k in 1:1)(
    fun_estimation(predictor = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/predictor_'),
                   NUM = k,
                   input.tmean = "~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData",
                   output.name = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/SWG_ERAI_ESD_tmean_'),
                   year.begin = 1979,
                   year.end = 1998)
  )
}


## 1999 - 2017
for (n in ens_pc){
  for (k in 1:1)(
    fun_estimation(predictor = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/predictor_'),
                   NUM = k,
                   input.tmean = "~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData",
                   output.name = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/SWG_ERAI_ESD_tmean_'),
                   year.begin = 1999,
                   year.end = 2017)
  )
}

# #### contribution ####
# library(ggplot2)
# library(reshape2)
# 
# MON=c(12,1:11)
# MON=matrix(MON,nrow=4,ncol=3,byrow=T)
# SEAS=c('DJF','MAM','JJA','SON')
# season=c('Winter','Spring','Summer','Fall')
# 
# casename = c(1:16)
# load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")
# casename = case_name
# ncase = length(casename)
# 
# load("~/Documents/LSCE/SWG/NA_slp_sst/cb_table.RData")         # cb_table
# load("~/Documents/LSCE/SWG/NA_slp_sst/predictor_name.RData")   # predictor_name
# 
# for (case in 1:ncase){
#   for (seas in 3){
#     
#     PC_name = unlist(predictor_name[which(!is.na(cb_table[case,]))])
#     
#     load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/1999_2017/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1999_2017_", case, ".RData"))
#     par_tt = coef(fit_stations_tt[[1]], matrix = T)
#     colnames(par_tt) = c("mean", "sd")
#     contrib = as.data.frame(apply(abs(par_tt[-1,]), 2, FUN = function(x){x/sum(x)}))
#     contrib$PC = PC_name
#     # contrib$type = "rea"
#     contrib$type = "1999-2017"
#     contrib = melt(contrib, id.vars = c("PC", "type"))
#     contrib$PC <- factor(contrib$PC, levels = PC_name)
#     contribution = contrib
#     
#     load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/1979_1998/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1979_1998_", case, ".RData"))
#     par_tt = coef(fit_stations_tt[[1]], matrix = T)
#     colnames(par_tt) = c("mean", "sd")
#     contrib = as.data.frame(apply(abs(par_tt[-1,]), 2, FUN = function(x){x/sum(x)}))
#     contrib$PC = PC_name
#     # contrib$type = "nat"
#     contrib$type = "1979-1998"
#     contrib = melt(contrib, id.vars = c("PC", "type"))
#     contrib$PC <- factor(contrib$PC, levels = PC_name)
#     contribution = rbind(contribution, contrib)
#     
#     p = ggplot(contribution)
#     p = p + geom_point(aes(x = PC, y = value, col = type)) 
#     p = p + facet_wrap(~variable, nrow = 2) + theme_bw()
#     p = p + labs(title = paste0("Contribution of principl components for t2m in ", season[seas], " (", casename[case], ")"), 
#                  x = "Composant",
#                  y = "Percent (X100%)")
#     print(p)
#     
#     dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Contribution/Contribution_t2m_1999_2017_", season[seas], "_", case, "_", casename[case], ".pdf"), width = 14, height = 7)
#   }
# }
# 

#### simulation ####
for (n in ens_pc){
  for (k in 1:1){
    cat(n, "\n")
    fun_simulation(predictor = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/predictor_'),
                   parameter = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/SWG_ERAI_ESD_tmean_'),
                   NUM = k, 
                   type = "nat",
                   output = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/'), 
                   y1 = 1999, y2 = 2017,
                   year.begin = 1979, year.end = 1998)
  }
}

for (n in ens_pc){
  for (k in 1:1){
    cat(n, "\n")
    fun_simulation(predictor = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/predictor_'),
                   parameter = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/SWG_ERAI_ESD_tmean_'),
                   NUM = k, 
                   type = "rea",
                   output = paste0('~/Documents/LSCE/SWG/slp_PC_diagnosis/', n, 'PC/'), 
                   y1 = 1999, y2 = 2017,
                   year.begin = 1999, year.end = 2017)
  }
}


#### SAMPLE 100 simulations ####
DATE = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))

## Gaussain distribution for temperature
Sample_gaussian_temp<- function(par_tt1,par_tt2){
  #a=date()
  #cat(a,'\n')
  Nb_stations=dim(par_tt1)[2]
  echant=array(NaN,dim=dim(par_tt1))
  for (station in 1:Nb_stations){
    for (j in 1:(dim(par_tt1)[1])){
      echant[j,station]=rnorm(1,mean=par_tt1[j,station],sd=par_tt2[j,station])
    }
  }
  #a=date()
  #cat(a,'\n')
  return(echant)
}

ens_pc = c(0,2:10,15,20)
npc = length(ens_pc)

for (k in 1:npc){
  i = ens_pc[k]
  dir.create(path = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", i, "PC/SIMU"))
  dir.create(path = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", i, "PC/SIMU/nat"))
  dir.create(path = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", i, "PC/SIMU/rea"))
  dir.create(path = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", i, "PC/SIMU/rea"))
}

for (k in 1:npc){
  pc = ens_pc[k]
  w = "rea"
  load(paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", pc, "PC/t2m_mean_sd_1999_2017_", w, "_1.RData"))
  range(mean)
  range(sd)
  
  ## make 100 runs
  for(i in 1:100){
    
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    cat('Run', NUM, '\n')
    
    Sample_temp = Sample_gaussian_temp(mean,sd)
    Sample=array(NaN,dim=dim(sd))
    Sample[,1]=Sample_temp
    
    Sample=round(Sample,digits=2)
    cat(range(Sample), "\n")
    
    filesimu = paste('~/Documents/LSCE/SWG/slp_PC_diagnosis/', pc, 'PC/SIMU/', w, '/SIMU_SWG_tmean_1999_2017_', w, '_', 1, '_run_', NUM, '.RData', sep="")
    save(Sample,DATE,file = filesimu)
  }
}

#### separator in season ####
MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
SEAS   = c('DJF', 'MAM', 'JJA', 'SON')

for (k in 1:npc){
  pc = ens_pc[k]
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", pc, "PC/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_", 1, "_run_", NUM, ".RData"))
    for (seas in 1:4){
      Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
      DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
      save(DATE_seas, Sample_seas, file = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", pc, "PC/SIMU/nat/SIMU_tmean_1999_2017_nat_1_", SEAS[seas], "_run_", NUM, ".RData"))
    }
  }
}


#### Distribution ####
source("function/CramerVonMisesTwoSamples.R")

load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
Obs = tmean[DATE_OBS$Y>=1999&DATE_OBS$Y<=2017,]

res = as.data.frame(array(NA, c(100,npc)))
colnames(res) = paste0(ens_pc, "PC")
pva = as.data.frame(array(NA, c(100,npc)))
colnames(pva) = paste0(ens_pc, "PC")
for (k in 1:npc){
  pc = ens_pc[k]
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", pc, "PC/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_1_run_", NUM, ".RData"))
    Sample = Sample[,1]
    cvm = CramerVonMisesTwoSamples(Sample, Obs)
    res[i, k] = cvm
    pva[i, k] = 1/6*exp(-cvm)
    # for (seas in 1:4){
    #   Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
    #   DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
    #   save(DATE_seas, Sample_seas, file = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", pc, "PC/SIMU/nat/SIMU_tmean_1999_2017_nat_1_", SEAS[seas], "_run_", NUM, ".RData"))
    # }
  }
}

res.melt = melt(res)
# res.melt$variable <- factor(res.melt$variable, labels = c("SD", levels(res.melt$variable)[-1]))
p <- ggplot(res.melt) + geom_boxplot(aes(x=variable, y=value)) + theme_bw() + 
  labs(x="slp - number of PC", y="Cramer von Mises", 
       title = paste0("Cramer von Mises - slp PC - 100 simulations (1999-2017)"))
print(p)
dev.print(pdf, file = "~/Documents/LSCE/SWG/slp_PC_diagnosis/Image/boxplot_Cramer-von_Mises_with_slp_PC.pdf", width = 7, height = 7)

pva.melt = melt(pva)
# pva.melt$variable <- factor(pva.melt$variable, labels = c("SD", levels(pva.melt$variable)[-1]))
p <- ggplot(pva.melt) + geom_boxplot(aes(x=variable, y=value)) + theme_bw() + 
  labs(x="slp - number of PC", y="p value", 
       title = paste0("p-value - slp PC - 100 simulations (1999-2017)"))
print(p)
dev.print(pdf, file = "~/Documents/LSCE/SWG/slp_PC_diagnosis/Image/boxplot_pvalue_with_slp_PC.pdf", width = 7, height = 7)


res = as.data.frame(array(NA, c(100,npc)))
colnames(res) = paste0(ens_pc, "PC")
pva = as.data.frame(array(NA, c(100,npc)))
colnames(pva) = paste0(ens_pc, "PC")

for (seas in 1:4){
  Obs_seas = Obs[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3]]
  for (k in 1:npc){
    pc = ens_pc[k]
    for (i in 1:100){
      NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
      load(paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", pc, "PC/SIMU/nat/SIMU_tmean_1999_2017_nat_1_", SEAS[seas], "_run_", NUM, ".RData"))
      cvm = CramerVonMisesTwoSamples(Sample_seas, Obs_seas)
      res[i, k] = cvm
      pva[i, k] = 1/6*exp(-cvm)
    }
  }
  res.melt = melt(res)
  res.melt$variable <- factor(res.melt$variable, labels = c("SD", levels(res.melt$variable)[-1]))
  p <- ggplot(res.melt) + geom_boxplot(aes(x=variable, y=value)) + theme_bw() + 
    labs(x="slp - number of PC", y="Cramer von Mises", 
         title = paste0("Cramer von Mises - slp PC - 100 simulations (", season[seas], " 1999-2017)"))
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/Image/boxplot_Cramer-von_Mises_with_slp_PC_", season[seas], ".pdf"), width = 7, height = 7)
  
  pva.melt = melt(pva)
  pva.melt$variable <- factor(pva.melt$variable, labels = c("SD", levels(pva.melt$variable)[-1]))
  p <- ggplot(pva.melt) + geom_boxplot(aes(x=variable, y=value)) + theme_bw() + 
    labs(x="slp - number of PC", y="p value", 
         title = paste0("p-value - slp PC - 100 simulations (", season[seas], " 1999-2017)"))
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/Image/boxplot_pvalue_with_slp_PC_", season[seas], ".pdf"), width = 7, height = 7)
}


#### Intensity ####
#### load tmean OBS
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

Obs = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
n = length(Obs)
DATE = as.Date((0:(n-1)), origin = "1999-01-01")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

ens_pc = c(0,2:10,15,20)
npc = length(ens_pc)
ncase = npc+1

list_mean_sd_nat = vector("list", ncase)
list_mean_sd_rea = vector("list", ncase)
for (k in 1:npc){
  i = ens_pc[k]
  
  load(paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", i, "PC/t2m_mean_sd_1999_2017_nat_1.RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  load(paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", i, "PC/t2m_mean_sd_1999_2017_rea_1.RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_rea[[k]] = m
}
for (k in ncase){
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/nat/t2m_mean_sd_1999_2017_nat_", 16, ".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/rea/t2m_mean_sd_1999_2017_rea_", 16, ".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_rea[[k]] = m
}

p1 = as.data.frame(matrix(NA, nrow = n, ncol = ncase))
for (l in 1:ncase){
  for (k in 1:n){
    p1[k,l] = pnorm(Obs[k], mean = list_mean_sd_rea[[l]][k,"mean"], sd = list_mean_sd_rea[[l]][k,"sd"], lower.tail = T)
  }
}

p0 = p1

i1 = as.data.frame(matrix(NA, nrow = n, ncol = ncase))
for (l in 1:ncase){
  for (k in 1:n){
    i1[k,l] = qnorm(p0[k,l], mean = list_mean_sd_nat[[l]][k,"mean"], sd = list_mean_sd_nat[[l]][k,"sd"], lower.tail = T)
  }
}

#### d_intensity: 1999-2017
d = Obs - i1; range(d)
colnames(d) = c(paste0(ens_pc, "PC"), "statio")
DATE = as.Date((0:(n-1)), origin = "1999-01-01")
d$date = DATE
save(d, file = "~/Documents/LSCE/SWG/slp_PC_diagnosis/d_intensity.RData")


#### daily `d_intensity` from 1999 to 2017 ####
ens_pc = c(0,2:10,15,20)
npc = length(ens_pc)
ncase = npc+1
load(file = "~/Documents/LSCE/SWG/slp_PC_diagnosis/d_intensity.RData")

melt_d = melt(d, id.vars = "date")
colnames(melt_d) = c("date", "case", "d")

ncase = ncol(d)-1

# # all_models
# {
#   p = ggplot(melt_d)
#   p = p + geom_line(aes(x=date, y=d, group=case, col=case)) + theme_bw()
#   cc <- colorRamps::matlab.like2(n = (ncol(d)-1))
#   p = p + scale_color_manual(values = cc)
#   output = "daily_d_intensity_1999_2017_all_models"
#   p = p + labs(title = output,
#                x = "date", y = "difference of intensity (°C)")
#   print(p)
#   dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/Image/", output, ".pdf"), width=11, height=7)
# 
# }

#### daily `d_intensity_tmean_summer` for each year from 1999 to 2017
# RMSE = function(m, o){
#   sqrt(mean((m - o)^2))
# }

for (year in 1999:2017){
  period_sel = intersect(which(year(d$date)==year) , which(month(d$date)==6|month(d$date)==7|month(d$date)==8))
  
  d1 = d[period_sel,]
  
  # rmse_table = as.data.frame(array(NA, c(1,npc)))
  # colnames(rmse_table) = paste0(ens_pc, "PC")
  # for (k in 1:npc){
  #   rmse_table[,k] = RMSE(d1[,k], d1[,(npc+1)])
  # }
  # rmse_table = sprintf("%.3f", rmse_table)
  
  d1 = melt(d1, id.vars = "date")
  colnames(d1) = c("date", "case", "d")
  d1$num = rep(1:length(period_sel), ncase)
  
  breaks = rep(NA,3)
  labels = rep(NA,3)
  for (k in 1:3){
    breaks[k] = which(month(d1$date)==unique(month(d1$date))[k])[1]
    labels[k] = as.character(d1$date[breaks[k]])
  }
  
  p = ggplot(d1)
  p = p + geom_line(aes(x=num, y=d, group=case, col=case)) + theme_bw()
  cc <- colorRamps::matlab.like2(n = ncase)
  # p = p + scale_color_manual(values = cc, labels = paste0(ens_pc, "PC (", rmse_table,")"))
  p = p + scale_color_manual(values = cc, labels = paste0(c("SysDyn", paste0(ens_pc[-1], "PC"), "statio")))
  p = p + scale_x_continuous(breaks = breaks, labels = labels) 
  p = p + labs(title = paste0("delta_Intensity_tmean_summer_in_", year),
               x = "date", y = expression(Delta*"I"))
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/Image/d_intensity_tmean_summer_", year, ".pdf"), width = 7, height=7)
  print(year)
}

#### Anomaly slp summer ####
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")
dI = d

load(file = "~/Documents/LSCE/SWG/slp_PC_diagnosis/d_intensity.RData")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

var = "slp"; 
# var = "sst"

cas = 3
if (cas == 1){
  load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1979_2017.RData"))
  cas = "1979_2017"
}else if(cas == 2){
  load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1979_1998.RData"))
  cas = "1979_1998"
}else if(cas == 3){
  load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1999_2017.RData"))
  cas = "1999_2017"
}

ens_pc = c(0,2:10,15,20)
npc = length(ens_pc)
mods = ens_pc

seas = 3

par(mfrow = c(6, 2))
## slp SD
{
  t_sup = which(dI[,2] - d[,"statio"]>0); length(t_sup)
  t_inf = which(dI[,2] - d[,"statio"]<0); length(t_inf)
  
  data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
  DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]
  
  for (seas in seas){
    idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
    t_sup_sea = intersect(t_sup, idxdates)
    data = apply(data_anomaly[,,t_sup_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_sup_1999-2017_", season[seas], "_slp_SD")
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    #dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Anomaly/", var, "/", cas, "/", output,".pdf"), width = 11, height = 5)
    
    
    t_inf_sea = intersect(t_inf, idxdates)
    data = apply(data_anomaly[,,t_inf_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_inf_1999-2017_", season[seas], "_slp_SD")
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    #dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Anomaly/", var, "/", cas, "/", output,".pdf"), width = 11, height = 5)
  }
}

## slp PC
for (model in c(1:5)){
  t_sup = which(d[,model] - d[,"statio"]>0); length(t_sup)
  t_inf = which(d[,model] - d[,"statio"]<0); length(t_inf)
  
  data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
  DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]
  
  for (seas in seas){
    idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
    t_sup_sea = intersect(t_sup, idxdates)
    data = apply(data_anomaly[,,t_sup_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_sup_1999-2017_", season[seas], "_slp_", ens_pc[model], "PC")
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    #dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Anomaly/", var, "/", cas, "/", output,".pdf"), width = 11, height = 5)
    
    
    t_inf_sea = intersect(t_inf, idxdates)
    data = apply(data_anomaly[,,t_inf_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_inf_1999-2017_", season[seas], "_slp_", ens_pc[model], "PC")
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    #dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Anomaly/", var, "/", cas, "/", output,".pdf"), width = 11, height = 5)
  }
}
dev.print(pdf,file = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/Image/Anomaly_slp_",season[seas],"_1.pdf"), width = 10, height = 10)
