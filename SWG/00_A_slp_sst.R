#### predictor ####
## combine all possibilities with slp and sst
nvar = 4
cb_table = matrix(NA, nrow = 2^nvar-1, ncol = nvar)
r = 1
for (k in 1:nvar){
  temp = combn(1:nvar, k)
  m = nrow(temp)
  for (i in 1:ncol(temp)){
    # cb_table[n, 1:m] = temp[,i]
    cb_table[r, temp[,i]] = "X"
    r = r + 1
  }
}
cb_table = as.data.frame(cb_table)
colnames(cb_table) = c("slp_PC", "slp_SD", "sst_PC", "sst_SD")
save(cb_table, file = "~/Documents/LSCE/SWG/NA_slp_sst/cb_table.RData")

case_name = rep(NA, nrow(cb_table))
for (k in 1:nrow(cb_table)){
  case_name[k] = paste(cb_table[k,!is.na(cb_table[k,])], collapse = "+")
}
save(case_name, file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")

library(xtable)
xtable(cb_table)

# predictor_name = vector("list", 4)
# predictor_name[[1]] = paste0("slp_PC_", 1:2)
# predictor_name[[2]] = paste0("slp_SD_", c("d", "t"))
# predictor_name[[3]] = paste0("sst_PC_", 1:2)
# predictor_name[[4]] = paste0("sst_PC_", c("d", "t"))

predictor_name = vector("list", 4)
predictor_name[[1]] = paste0("slp_PC_", 1:4)
predictor_name[[2]] = paste0("slp_SD_", c("d", "t"))
predictor_name[[3]] = paste0("sst_PC_", 1:4)
predictor_name[[4]] = paste0("sst_PC_", c("d", "t"))
save(predictor_name, file = "~/Documents/LSCE/SWG/NA_slp_sst/predictor_name.RData")

## It's OK for cb_table and case_name
load("~/Documents/LSCE/SWG/NA_slp_sst/cb_table.RData")         # cb_table
load("~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")        # case_name 
load("~/Documents/LSCE/SWG/NA_slp_sst/predictor_name.RData")   # predictor_name


#### normalize msl PCS : predictor_season_1
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  # load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/msl_",SEAS[seas],"_PCS.RData"))
  {
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/msl_", SEAS[seas], "_pca.RData"))
    eig.val <- get_eigenvalue(pca)
    print(eig.val[4,])
    PCS = as.data.frame(pca$x)[,1:4]
  }
  PCS = apply(PCS, MARGIN = 2, scale)
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_1.RData"))
}


#### normalize msl SysDyn : predictor_season_2
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/msl_",SEAS[seas],"_SysDyn.RData"))
  # print(c(colMeans(SysDyn), apply(SysDyn, 2, FUN = sd)))
  PCS = apply(SysDyn, MARGIN = 2, scale)
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_2.RData"))
}


#### normalize sst PCS : predictor_season_3
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  # load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/sst_",SEAS[seas],"_PCS.RData"))
  {
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/sst_", SEAS[seas], "_pca.RData"))
    eig.val <- get_eigenvalue(pca)
    print(eig.val[4,])
    PCS = as.data.frame(pca$x)[,1:4]
  }
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  PCS = apply(PCS, MARGIN = 2, scale)
  # PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_3.RData"))
}


#### normalize sst SysDyn : predictor_season_4
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/sst_",SEAS[seas],"_SysDyn.RData"))
  PCS = apply(SysDyn, MARGIN = 2, scale)
  # PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_4.RData"))
}


#### for the rest predictor_season_XXX: combn(1:4, 2)
cb = combn(1:4, 2)
SEAS=c('DJF','MAM','JJA','SON')
for (k in 1:ncol(cb)){
  for (seas in 1:4){
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[1,k], ".RData"))
    predictor = PCS
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[2,k], ".RData"))
    predictor = cbind(predictor, PCS)
    print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
    PCS = predictor
    colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
    PCS = as.data.frame(PCS)
    n = 4 + k
    save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", n, ".RData"))
  }
}


#### combn(1:4, 3)
cb = combn(1:4, 3)
SEAS=c('DJF','MAM','JJA','SON')
for (k in 1:ncol(cb)){
  for (seas in 1:4){
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[1,k], ".RData"))
    predictor = PCS
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[2,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[3,k], ".RData"))
    predictor = cbind(predictor, PCS)
    print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
    PCS = predictor
    colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
    PCS = as.data.frame(PCS)
    n = 10 + k
    save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", n, ".RData"))
  }
}

#### combn(1:4, 4)
cb = combn(1:4, 4)
SEAS=c('DJF','MAM','JJA','SON')
for (k in 1:ncol(cb)){
  for (seas in 1:4){
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[1,k], ".RData"))
    predictor = PCS
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[2,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[3,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", cb[4,k], ".RData"))
    predictor = cbind(predictor, PCS)
    print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
    PCS = predictor
    colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
    PCS = as.data.frame(PCS)
    n = 14 + k
    save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_",SEAS[seas],"_", n, ".RData"))
  }
}

#### estimation t2m ####
for (k in 1:15)(
  fun_estimation(predictor = '~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_',
                 NUM = k,
                 input.tmean = "~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData",
                 output.name = '~/Documents/LSCE/SWG/NA_slp_sst/tmean/1979_1998/SWG_ERAI_ESD_tmean_',
                 year.begin = 1979,
                 year.end = 1998)
)

for (k in 1:15)(
  fun_estimation(predictor = '~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_',
                 NUM = k,
                 input.tmean = "~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData",
                 output.name = '~/Documents/LSCE/SWG/NA_slp_sst/tmean/1999_2017/SWG_ERAI_ESD_tmean_',
                 year.begin = 1999,
                 year.end = 2017)
)

#### simulation t2m ####
for (k in 1:15){
  cat(k, "\n")
  fun_simulation(predictor = '~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_',
                 parameter = '~/Documents/LSCE/SWG/NA_slp_sst/tmean/1979_1998/SWG_ERAI_ESD_tmean_',
                 NUM = k, 
                 type = "nat",
                 output = '~/Documents/LSCE/SWG/NA_slp_sst/tmean/nat/', 
                 y1 = 1999, y2 = 2017,
                 year.begin = 1979, year.end = 1998)
}

for (k in 1:15){
  cat(k, "\n")
  fun_simulation(predictor = '~/Documents/LSCE/SWG/NA_slp_sst/pca/predictor_',
                 parameter = '~/Documents/LSCE/SWG/NA_slp_sst/tmean/1999_2017/SWG_ERAI_ESD_tmean_',
                 NUM = k, 
                 type = "rea",
                 output = '~/Documents/LSCE/SWG/NA_slp_sst/tmean/rea/', 
                 y1 = 1999, y2 = 2017,
                 year.begin = 1999, year.end = 2017)
}


#### Sampling for simulations ####
for (k in 1:15){
  dir.create(path = paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/SIMU/", k))
}

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

##
# mods = c(1,2,4,5,7,9,12,16)
for (k in 1:15){
  w = "nat"
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/", w, "/t2m_mean_sd_1999_2017_", w, "_", k, ".RData"))
  # sum(is.na(mean))
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
    
    filesimu = paste('~/Documents/LSCE/SWG/NA_slp_sst/tmean/SIMU/', k, '/SIMU_SWG_tmean_1999_2017_', w, '_', k, '_run_', NUM, '.RData', sep="")
    save(Sample,DATE,file = filesimu)
  }
}

#### separator in season ####
MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
SEAS   = c('DJF', 'MAM', 'JJA', 'SON')

# mods = c(1,2,4,5,7,9,12)
for (case in 1:15){
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/SIMU/", case, "/SIMU_SWG_tmean_1999_2017_nat_", case, "_run_", NUM, ".RData"))
    for (seas in 1:4){
      Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
      DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
      save(DATE_seas, Sample_seas, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/SIMU/", case, "/SIMU_tmean_", SEAS[seas], "_run_", NUM, ".RData"))
    }
  }
}

#### Distribution ####
MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
SEAS   = c('DJF', 'MAM', 'JJA', 'SON')

CramerVonMisesTwoSamples <- function(S1, S2){
  xS1 = sort(S1)
  M = length(xS1)
  xS2 = sort(S2)
  N = length(xS2)
  a = data.frame(val = xS1, rang = seq(M), ens = rep(1, M))
  b = data.frame(val = xS2, rang = seq(N), ens = rep(2, N))
  c = rbind(a, b)
  c = c[order(c$val), ]
  c = data.frame(c, rangTot = seq(M + N))
  dtfM = c[which(c$ens == 1), ]
  dtfN = c[which(c$ens == 2), ]
  somN = sum((dtfN$rang - dtfN$rangTot)^2)
  somM = sum((dtfM$rang - dtfM$rangTot)^2)
  U = N * somN + M * somM
  CvM = ((U/(as.numeric(N) * as.numeric(M)))/(N + M)) - ((4 * as.numeric(M) * as.numeric(N) - 1)/(6 * (M + N)))
  return(CvM)
}

load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
Obs = tmean[DATE_OBS$Y>=1999&DATE_OBS$Y<=2017,]

load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")
ncase = 15
res = as.data.frame(array(NA, c(100, ncase)))
colnames(res) = case_name
pva = as.data.frame(array(NA, c(100, ncase)))
colnames(pva) = case_name
for (seas in 1:4){
  Obs_seas = Obs[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3]]
  for (case in 1:ncase){
    for (i in 1:100){
      NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
      load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/SIMU/", case, "/SIMU_tmean_", SEAS[seas], "_run_", NUM, ".RData"))
      cvm = CramerVonMisesTwoSamples(Sample_seas, Obs_seas)
      res[i, case] = cvm
      pva[i, case] = 1/6*exp(-cvm)
    }
  }
  res.melt = res
  colnames(res.melt) = 1:ncase
  res.melt = melt(res.melt)
  p <- ggplot(res.melt) + geom_boxplot(aes(x=variable, y=value, col=variable)) + theme_bw() + 
    scale_color_discrete(label = paste0(1:ncase, "_", case_name)) + 
    labs(x="predictor", y="Cramer von Mises", 
         title = paste0("Cramer von Mises - 100 simulations (", season[seas], " 1999-2017)"))
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/boxplot_Cramer-von_Mises_", season[seas], ".pdf"), width = 7, height = 7)
  
  pva.melt = pva
  colnames(pva.melt) = 1:ncase
  pva.melt = melt(pva.melt)
  p <- ggplot(pva.melt) + geom_boxplot(aes(x=variable, y=value, col=variable)) + theme_bw() +
    scale_color_discrete(label = paste0(1:ncase, "_", case_name)) + 
    labs(x="predictor", y="p value",
         title = paste0("p-value - 100 simulations (", season[seas], " 1999-2017)"))
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/boxplot_p-value_", season[seas], ".pdf"), width = 7, height = 7)
}

res = as.data.frame(array(NA, c(100, ncase)))
colnames(res) = case_name
pva = as.data.frame(array(NA, c(100, ncase)))
colnames(pva) = case_name

for (case in 1:ncase){
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/SIMU/", case, "/SIMU_SWG_tmean_1999_2017_nat_", case, "_run_", NUM, ".RData"))
    cvm = CramerVonMisesTwoSamples(Sample, Obs)
    res[i, case] = cvm
    pva[i, case] = 1/6*exp(-cvm)
    # for (seas in 1:4){
    #   Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
    #   DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
    #   save(DATE_seas, Sample_seas, file = paste0("~/Documents/LSCE/SWG/slp_PC_diagnosis/", pc, "PC/SIMU/nat/SIMU_tmean_1999_2017_nat_1_", SEAS[seas], "_run_", NUM, ".RData"))
    # }
  }
}

res.melt = res
colnames(res.melt) = 1:ncase
res.melt = melt(res.melt)
p <- ggplot(res.melt) + geom_boxplot(aes(x=variable, y=value, col=variable)) + theme_bw() + 
  scale_color_discrete(label = paste0(1:ncase, "_", case_name)) + 
  labs(x="predictor", y="Cramer von Mises", 
       title = paste0("Cramer von Mises - 100 simulations (1999-2017)"))
print(p)
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/boxplot_Cramer-von_Mises.pdf"), width = 7, height = 7)

pva.melt = pva
colnames(pva.melt) = 1:ncase
pva.melt = melt(pva.melt)
p <- ggplot(pva.melt) + geom_boxplot(aes(x=variable, y=value, col=variable)) + theme_bw() +
  scale_color_discrete(label = paste0(1:ncase, "_", case_name)) + 
  labs(x="predictor", y="p value",
       title = paste0("p-value - 100 simulations (1999-2017)"))
print(p)
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/boxplot_p-value.pdf"), width = 7, height = 7)


#### Diagnostic - acf ####
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")

for (seas in 1:4){
  # Obs as tmean - daily mean temperature from 1999 to 2017 in the zone [48N, 50N]x[01E, 04E]
  Obs = tmean[which(rownames(tmean)=="1999-01-01"):which(rownames(tmean)=="2017-12-31"),]
  DATE_Obs = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))
  
  MON    = c(12, 1:11)
  MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
  SEAS   = c('DJF', 'MAM', 'JJA', 'SON')
  Obs = Obs[DATE_Obs$m==MON[seas,1] | DATE_Obs$m==MON[seas,2] | DATE_Obs$m==MON[seas,3]]
  
  # NUM: number with 3 digits
  k = 1 
  NUM=cbind(formatC(k, digits = 0, width = 3, format = "f", flag = "0"))
  
  # mods = c(1,2,4,5,7,9,12)
  for (ncase in 1:15){
    # load one Sample (see k and NUM)
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/SIMU/", ncase, "/SIMU_tmean_", SEAS[seas], "_run_001.RData"))
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
      p = ggplot(data = data_acf, mapping = aes(x = lag, y = acf, width=.75))
      p = p + geom_bar(stat = "identity", aes(fill = type), position = "dodge")
      p = p + scale_fill_manual(values = color)
      p = p + xlab("day")
      p = p + ggtitle(paste0("acf - Obs vs Sample_", k, " for ", case_name[ncase])) + theme_bw()
      print(p)
    }
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/acf/acf_", SEAS[seas],"_case_", ncase, ".pdf"), width = 9, height = 9)
    
  }
}


#### stationnary model ####
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

data1 = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
DATE1 = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))

data2 = tmean[DATE_OBS$Y>=1979 & DATE_OBS$Y<=1998,]
DATE2 = atoms(timeSequence(from="1979-01-01",to="1998-12-31",by='day'))

## nat
mean = array(NaN,c(length(data1),1)) # matrice de probas d'occurrence
sd = array(NaN,c(length(data1),1))

for (seas in 1:4){
  idxdates1 = which(DATE1['m']==MON[seas,1] | DATE1['m']==MON[seas,2] | DATE1['m']==MON[seas,3])
  idxdates2 = which(DATE2['m']==MON[seas,1] | DATE2['m']==MON[seas,2] | DATE2['m']==MON[seas,3])
  
  mean[idxdates1,] = mean(data2[idxdates2])   ; cat(mean(data2[idxdates2]), "\n"); cat(mean(data1[idxdates1]), "\n")
  sd[idxdates1,]   = sd(data2[idxdates2])     ; #cat(sd(data2[idxdates2])  , "\n"); cat(sd(data1[idxdates1]), "\n")
}
range(mean)
range(sd)
DATE = DATE1
save(mean, sd, DATE, file = "~/Documents/LSCE/SWG/NA_slp_sst/tmean/nat/t2m_mean_sd_1999_2017_nat_16.RData")

## rea
mean = array(NaN,c(length(data1),1)) # matrice de probas d'occurrence
sd = array(NaN,c(length(data1),1))

for (seas in 1:4){
  idxdates1 = which(DATE1['m']==MON[seas,1] | DATE1['m']==MON[seas,2] | DATE1['m']==MON[seas,3])
  
  mean[idxdates1,] = mean(data1[idxdates1])   ; cat(mean(data1[idxdates1]), "\n")
  sd[idxdates1,]   = sd(data1[idxdates1]) 
}
range(mean)
range(sd)
DATE = DATE1
save(mean, sd, DATE, file = "~/Documents/LSCE/SWG/NA_slp_sst/tmean/rea/t2m_mean_sd_1999_2017_rea_16.RData")


#### Intensity ####
#### load tmean OBS
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

Obs = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
n = length(Obs)
DATE = as.Date((0:(n-1)), origin = "1999-01-01")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

#### case 
case_name = c(1:16); case_sel = c(1:16)
case_name = case_name[case_sel]
ncase = length(case_name)

list_mean_sd_nat = vector("list", ncase)
list_mean_sd_rea = vector("list", ncase)
for (k in 1:ncase){
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/nat/t2m_mean_sd_1999_2017_nat_", case_sel[k], ".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/rea/t2m_mean_sd_1999_2017_rea_", case_sel[k], ".RData"))
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
colnames(d) = case_name
DATE = as.Date((0:(n-1)), origin = "1999-01-01")
d$date = DATE
save(d, file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")


#### daily `d_intensity` from 1999 to 2017
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")
load("~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")

case_name = c(case_name, "statio")
colnames(d) = c(case_name, "date")

melt_d = melt(d, id.vars = "date")
colnames(melt_d) = c("date", "case", "d")

## all_models
p = ggplot(melt_d)
p = p + geom_line(aes(x=date, y=d, group=case, col=case)) + theme_bw()
cc <- colorRamps::matlab.like2(n = 16)
p = p + scale_color_manual(values = cc)
output = "daily_d_intensity_1999_2017_all_models"
p = p + labs(title = output,
             x = "date", y = "difference of intensity (°C)")
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/", output, ".pdf"), width=11, height=7)

## different_models_in_one_figure
p = ggplot(melt_d)
p = p + geom_line(aes(x=date, y=d)) + theme_bw()
p = p + facet_wrap(.~case, nrow = 4)
output = "daily_d_intensity_1999_2017_different_models"
p = p + labs(title = output, x = "date", y = "difference of intensity (°C)")
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/", output, ".pdf"), width=11, height=11)

## different_model_in_different_figure
for (k in 1:16){
  p = ggplot(melt_d[melt_d$case==case_name[k],])
  p = p + geom_line(aes(x=date, y=d)) + theme_bw()
  cc <- colorRamps::matlab.like2(n = 16)
  p = p + scale_color_manual(values = cc)
  output = paste0("daily_d_intensity_1999_2017_of_model_", k, "_", case_name[k])
  p = p + labs(title = output,
               x = "date", y = "difference of intensity (°C)")
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/", output, ".pdf"), width=11, height=7)
}

## different_model_in_different_figure_in_different_year
DATE = d$date
for (year in 1999:2017){
  for (k in 2){
    p = ggplot(melt_d[melt_d$case==case_name[k] & year(DATE)==year,])
    p = p + geom_line(aes(x=date, y=d)) + theme_bw()
    output = paste0("daily_d_intensity_of_model_", year, "_", k, "_", case_name[k])
    p = p + labs(title = output,
                 x = "date", y = "difference of intensity (°C)")
    print(p)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/", output, ".pdf"), width=11, height=7)
  }
}


#### daily `d_intensity_tmean_summer` for each year from 1999 to 2017
for (year in 1999:2017){
  period_sel = intersect(which(year(d$date)==year) , which(month(d$date)==6|month(d$date)==7|month(d$date)==8))
  
  d1 = d[period_sel,]
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
  cc <- colorRamps::matlab.like2(n = 16)
  p = p + scale_color_manual(values = cc)
  p = p + scale_x_continuous(breaks = breaks, labels = labels) 
  p = p + labs(title = paste0("delta_Intensity_tmean_summer_in_", year),
               x = "date", y = "difference of intensity (°C)")
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/d_intensity_tmean_summer_", year, ".pdf"), width = 7, height=7)
  print(year)
}


#### d_mean
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")

case_name = c(1:16); case_sel = c(1:16)
case_name = case_name[case_sel]
ncase = length(case_name)

load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")
case_name = c(case_name, "statio")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

summer = which(month(d$date)==6|month(d$date)==7|month(d$date)==8)

d_mean = as.data.frame(colMeans(d[summer,1:ncase]))
colnames(d_mean) = "average"
d_mean$number = 1:16
d_mean$model = case_name
d_mean$model = factor(d_mean$model, levels = case_name)

p = ggplot(d_mean) 
p = p + geom_text(aes(x=number, y=average, label = number)) +  theme_bw()
p = p + geom_point(aes(x=number, y=average), col = "red", size = 0.5)
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/average_delta_of_intensity_summer_1999-2017.pdf"), width = 9, height=7)

mods = c(1,2,4,5,7,9,12,16)
# mods = c(5,7,9,12,16)

p = ggplot(d_mean[mods,]) 
p = p + geom_point(aes(x=number, y=average)) +  theme_bw()
p = p + geom_text(aes(x=number, y=average, label=model), hjust = 0, nudge_x = 0.5) 
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/average_delta_of_intensity_summer_1999-2017_best.pdf"), width = 9, height=7)


## d_year_1
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")
case_name = c(case_name, "statio")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
summer = which(month(d$date)==6|month(d$date)==7|month(d$date)==8)

ncase = ncol(d)-1

d = d[summer,]
d_year = matrix(NA, nrow = (2017-1999+1), ncol = ncase)
yy = unique(year(d$date))
for (y in 1:nrow(d_year)){
  d_year[y,] = colMeans(d[year(d$date)==yy[y],1:ncase])
}
d_year = as.data.frame(d_year)
colnames(d_year) = case_name
d_year$year = 1999:2017

d_year_1 = melt(d_year, id.vars = "year")
colnames(d_year_1) = c("year", "case", "intensity")

## annual delta of intensity from 1999 to 2017
mods = c(1,2,4,5,7,9,12,16)
d_year_1_mods = d_year_1[d_year_1$case%in%case_name[mods],]
p = ggplot(d_year_1_mods)
p = p + geom_line(data = d_year_1_mods[d_year_1_mods$case!="statio",], aes(x=year, y=intensity, col=case)) + theme_bw()
cc <- colorRamps::matlab.like2(n = length(mods))
p = p + scale_color_manual(values = cc)
p = p + geom_line(data = d_year_1_mods[d_year_1_mods$case=="statio",], aes(x = year, y=intensity), col="black")
output = paste0("annual_delta_of_intensity_summer_1999_2017_best")
p = p + labs(title = output,
             x = "date", y = "difference of intensity (°C)")
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/", output, ".pdf"), width = 9, height=7)

## annual delta of intensity from 1999 to 2016 for case XXX
for (k in mods){
  p = ggplot(d_year_1_mods[d_year_1_mods$case == case_name[k],])
  p = p + geom_line(aes(x=year, y=intensity), col="red") + theme_bw()
  p = p + geom_line(data = d_year_1_mods[d_year_1_mods$case=="statio",], aes(x = year, y=intensity), col="black")
  output = paste0("annual_delta_of_intensity_summer_1999_2017_", k, "_", case_name[k])
  p = p + labs(title = output,
               x = "date", y = "difference of intensity (°C)")
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Intensity/", output, ".pdf"), width = 7, height=7)
}


#### FAR ####
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")
case_name = c(case_name, "statio")

case_sel = 1:16
case_name = case_name[case_sel]
ncase = length(case_name)

list_mean_sd_nat = vector("list", ncase)
list_mean_sd_rea = vector("list", ncase)
for (k in 1:ncase){
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/nat/t2m_mean_sd_1999_2017_nat_", case_sel[k], ".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/rea/t2m_mean_sd_1999_2017_rea_", case_sel[k], ".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_rea[[k]] = m
}

## load tmean OBS
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

Obs = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
n = length(Obs)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

DATE = as.Date((0:(n-1)), origin = "1999-01-01")

##
FAR_case = c("FAR", "RR", "FARanda")
FAR_case = FAR_case[2]

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
    prob_thres$date = rep(DATE, ncase)
    
    if (FAR_case == "FAR"){
      prob_thres$FAR = 1 - prob_thres_nat/prob_thres_rea
    }else if (FAR_case == "RR"){
      prob_thres$FAR = prob_thres$rea/prob_thres$nat
    }else if (FAR_case == "FARanda"){
      prob_thres$FAR = (prob_thres$rea-prob_thres$nat)/(prob_thres$nat+prob_thres$rea)
    }
    
    prob_thres$case = rep(case_name, each = length(DATE))
    
    period_sel = intersect(which(year(prob_thres$date)==year) , which(month(prob_thres$date)==7|month(prob_thres$date)==6|month(prob_thres$date)==8))
    output = paste0(FAR_case, "_", year, "_summer_tmean_over_", thres, "_degree")
    
    FAR_sel = prob_thres[period_sel,]
    FAR_sel$num = rep(1:(length(period_sel)/ncase), ncase)
    
    FAR_sel$case <- factor(FAR_sel$case, levels = case_name)
    
    breaks = rep(NA,3)
    labels = rep(NA,3)
    for (k in 1:3){
      breaks[k] = which(month(FAR_sel$date)==unique(month(FAR_sel$date))[k])[1]
      labels[k] = as.character(FAR_sel$date[breaks[k]])
    }
    
    p = ggplot(FAR_sel) 
    p = p + geom_line(aes(x=num, y=FAR, group=case, col=case)) + theme_bw()
    p = p + scale_x_continuous(breaks = breaks, labels = labels) # + ylim(-1,1)
    cc <- colorRamps::matlab.like2(n = 16)
    p = p + scale_color_manual(values = cc)
    p = p + labs(title = output)
    print(p)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/FAR/", output, ".pdf"), width = 14, height=7)
  }
}


## average
for (year in 2003:2003){
  for (thres in 10:30){
    # thres = 20
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
    prob_thres$date = rep(DATE, ncase)
    
    if (FAR_case == "FAR"){
      prob_thres$FAR = 1 - prob_thres_nat[,1]/prob_thres_rea[,1]
    }else if (FAR_case == "RR"){
      prob_thres$FAR = prob_thres$rea/prob_thres$nat
    }else if (FAR_case == "FARanda"){
      prob_thres$FAR = (prob_thres$rea-prob_thres$nat)/(prob_thres$nat+prob_thres$rea)
    }
    
    prob_thres$case = rep(case_name, each = length(DATE))
    
    period_sel = intersect(which(year(prob_thres$date)==year) , which(month(prob_thres$date)==7|month(prob_thres$date)==6|month(prob_thres$date)==8))
    output = paste0("average_", FAR_case, "_", year, "_summer_tmean_over_", thres, "_degree")
    
    FAR_sel = prob_thres[period_sel,]
    FAR_sel$num = rep(1:(length(period_sel)/ncase), ncase)
    FAR_sel$case <- factor(FAR_sel$case, levels = case_name)
    
    FAR = as.data.frame(matrix(NA, nrow = 16, ncol = 2))
    colnames(FAR) = c("value", "model")
    for (k in 1:16){
      FAR[k, 1] = mean(FAR_sel$FAR[FAR_sel$case==case_name[k]])
      FAR[k, 2] = case_name[k]
    }
    FAR$number = c(1:16)
    
    p = ggplot(FAR) 
    p = p + geom_text(aes(x=number, y=value, label = number)) +  theme_bw()
    p = p + geom_point(aes(x=number, y=value), col = "red", size = 0.5)
    p = p + labs(title = output)
    print(p)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/FAR/", output, ".pdf"), width = 14, height=7)
  }
}
#

#### Criterion ####
DATE = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))
date = as.Date(0:(nrow(DATE)-1), origin = "1999-01-01")

## NAO
NAO = read.table(file = "~/Documents/LSCE/RDATA/Climate_indices/icpc_nao_daily.dat.txt")
colnames(NAO) = c("Y", "m", "d", "value")
NAO = NAO[NAO$Y>=1999 & NAO$Y<=2017, ]
data = NAO
for (k in 1:nrow(DATE)){
  if (!identical(as.numeric(data[k,1:3]), as.numeric(DATE[k,1:3]))){
    print(k)
    nl = matrix(NA, nrow = 1, ncol = 4)
    colnames(nl) = c("Y", "m", "d", "value")
    nl[1,1:3] = as.numeric(DATE[k,1:3])
    
    data = rbind(data[1:(k-1),], nl, data[k:nrow(data),])
  }
}
NAO_1 = data
# save(NAO_1, file = "~/Documents/LSCE/RDATA/Climate_indices/NAO_1.RData")

plot(date, NAO_1[,4], type = "l")

## AO
AO = read.table(file = "~/Documents/LSCE/RDATA/Climate_indices/icpc_ao_daily.dat.txt")
colnames(AO) = c("Y", "m", "d", "value")
AO = AO[AO$Y>=1999 & AO$Y<=2017, ]
data = AO
for (k in 1:nrow(DATE)){
  if (!identical(as.numeric(data[k,1:3]), as.numeric(DATE[k,1:3]))){
    print(k)
    nl = matrix(NA, nrow = 1, ncol = 4)
    colnames(nl) = c("Y", "m", "d", "value")
    nl[1,1:3] = as.numeric(DATE[k,1:3])
    
    data = rbind(data[1:(k-1),], nl, data[k:nrow(data),])
  }
}
AO_1 = data
# save(AO_1, file = "~/Documents/LSCE/RDATA/Climate_indices/AO_1.RData")

## ENSO
ENSO = read.table(file = "~/Documents/LSCE/RDATA/Climate_indices/inino34_daily.dat.txt")
head(ENSO,1); tail(ENSO,1)
DATE_ENSO = atoms(timeSequence(from="1981-09-01",to="2018-11-01",by='day'))
ENSO = cbind(DATE_ENSO[,1:3], ENSO[,2])
colnames(ENSO) = c("Y", "m", "d", "value")
ENSO = ENSO[ENSO$Y>=1999 & ENSO$Y<=2017, ]
# save(ENSO, file = "~/Documents/LSCE/RDATA/Climate_indices/ENSO.RData")


#### load d_intensity
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")

mods = c(1,2,4,5,7,9,12,16)

cor_tb = matrix(data = NA, nrow = 3, ncol = length(mods))
colnames(cor_tb) = c(case_name, "statio")[mods]
rownames(cor_tb) = c("NAO", "AO", "ENSO")

method = "spearman"
for (k in 1:length(mods)){
  cor_tb[1,k] = cor(NAO_1[,4], d[,mods[k]], use = "complete.obs", method = method)
  cor_tb[2,k] = cor( AO_1[,4], d[,mods[k]], use = "complete.obs", method = method)
  cor_tb[3,k] = cor( ENSO[,4], d[,mods[k]], use = "complete.obs", method = method)
}
save(cor_tb, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/cor_tb_", method, ".RData"))

# load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/cor_tb_", method, ".RData"))


cor_tb = as.data.frame(matrix(data = NA, nrow = length(mods), ncol = 3))
rownames(cor_tb) = c(case_name, "statio")[mods]
colnames(cor_tb) = c("NAO", "AO", "ENSO")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

seas = 4
method = "spearman"
for (k in 1:length(mods)){
  s = which(month(d$date)==MON[seas,1]|month(d$date)==MON[seas,2]|month(d$date)==MON[seas,3])
  cor_tb[k,1] = cor(NAO_1[s,4], d[s,mods[k]], use = "complete.obs", method = method)
  cor_tb[k,2] = cor( AO_1[s,4], d[s,mods[k]], use = "complete.obs", method = method)
  cor_tb[k,3] = cor( ENSO[s,4], d[s,mods[k]], use = "complete.obs", method = method)
}
cor_tb$"NAO+ENSO" = abs(cor_tb$NAO) + abs(cor_tb$ENSO)
cor_tb$"AO+ENSO" = abs(cor_tb$AO) + abs(cor_tb$ENSO)
# save(cor_tb, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/cor_tb_", season[seas], "_", method, ".RData"))

xtable(cor_tb)

#### Anomaly : based on "statio" model ####
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

var = "slp"; 
var = "sst"

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

# mods = c(1,2,4,5,7,9,12)
mods = 1:15
for (model in mods){
  t_sup = which(d[,model] - d[,16]>0); length(t_sup)
  t_inf = which(d[,model] - d[,16]<0); length(t_inf)
  
  data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
  DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]
  
  for (seas in c(1,3)){
    idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
    t_sup_sea = intersect(t_sup, idxdates)
    data = apply(data_anomaly[,,t_sup_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_sup_1999-2017_", season[seas], "_model_", model, "_", case_name[model])
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Anomaly/", var, "/", cas, "/", output,".pdf"), width = 11, height = 5)
    
    
    t_inf_sea = intersect(t_inf, idxdates)
    data = apply(data_anomaly[,,t_inf_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_inf_1999-2017_", season[seas], "_model_", model, "_", case_name[model])
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Anomaly/", var, "/", cas, "/", output,".pdf"), width = 11, height = 5)
  }
}


#### Anomaly : based on other model ####
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

var = "sst"

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

mods = 7
mods_ref = 1
for (model in mods){
  t_sup = which(d[,model] - d[,mods_ref]>0); length(t_sup)
  t_inf = which(d[,model] - d[,mods_ref]<0); length(t_inf)
  
  data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
  DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]
  
  for (seas in c(1,3)){
    idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
    t_sup_sea = intersect(t_sup, idxdates)
    data = apply(data_anomaly[,,t_sup_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_sup_1999-2017_", season[seas], "_model_", model, "_", case_name[model], "_based_on_", case_name[mods_ref], "_", cas)
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range, 
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Anomaly/", output,".pdf"), width = 11, height = 5)
    
    t_inf_sea = intersect(t_inf, idxdates)
    data = apply(data_anomaly[,,t_inf_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_inf_1999-2017_", season[seas], "_model_", model, "_", case_name[model], "_based_on_", case_name[mods_ref], "_", cas)
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range, 
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Anomaly/", output,".pdf"), width = 11, height = 5)
  }
}

#### dim and theta ####
var = "slp"
if (var == "slp"){
  load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_dim.RData")
  load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_theta.RData")
  plot(dim, theta, main = paste0("theta vs dim of ", var))
}else if (var == "z500"){
  load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_dim.RData")
  load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_theta.RData")
  plot(dim, theta, main = paste0("theta vs dim of ", var))
}else if (var == "sst"){
  load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_small_dim.RData")
  load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_small_theta.RData")
  plot(dim, theta, main = paste0("theta vs dim of ", var))
}

load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")

DATE = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))
dim   = dim[DATE$Y>=1999 & DATE$Y<=2017]
theta = theta[DATE$Y>=1999 & DATE$Y<=2017]
DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]
plot(dim, theta, main = paste0("theta vs dim of ", var))

data = as.data.frame(cbind(dim, theta))

load(file = "~/Documents/LSCE/SWG/NA_slp_sst/d_intensity.RData")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

mods = c(1,2,4,5,7,9,12)
# mods = 10:15
for (model in mods){
  t_sup = which(d[,model] - d[,16]>0); length(t_sup)
  t_inf = which(d[,model] - d[,16]<0); length(t_inf)
  data$type = NA
  for (seas in 1:4){
    
    idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
    
    t_sup_sea = intersect(t_sup, idxdates)
    t_inf_sea = intersect(t_inf, idxdates)
    
    data$type[t_sup_sea] = "sup"
    data$type[t_inf_sea] = "inf"
    
    data_sea = data[idxdates,]
    sum(is.na(data_sea))
    
    library(ggplot2)
    output = paste0("dim_theta_1999-2017_by_", var, "_in_", season[seas], "_for_model_", model, "_", case_name[model])
    p = ggplot(data_sea) + geom_point(aes(x = dim, y = theta, col = type)) + theme_bw()
    p = p + labs(title = output)
    print(p)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/SysDyn/", output, ".pdf"), width = 7, height = 7)
  }
}

model = 1; seas = 3
data$diff = d[,model] - d[,16]

idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
data_sea = data[idxdates,]

var = "slp"; 
output = paste0("dim_theta_1999-2017_by_", var, "_in_", season[seas], "_for_model_", model, "_", case_name[model], "_rainbow")
p = ggplot(data_sea) + geom_point(aes(x = dim, y = theta, col = diff)) + theme_bw()
p = p + scale_color_gradientn(colours = matlab.like2(10))
p
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/SysDyn/", output, ".pdf"), width = 7, height = 7)


## dim & theta selon NAO
load("~/Documents/LSCE/RDATA/Climate_indices/NAO_1.RData")

data_NAO = data
data_NAO$NAO = NAO_1[,4]
data_sea = data_NAO

## all seasons 
output = paste0("dim_theta_with_NAO_1999-2017_by_", var, "_rainbow")
p = ggplot(data_sea) + geom_point(aes(x = dim, y = theta, col = NAO)) + theme_bw()
# cc <- colorRamps::matlab.like2(n = 16)
# p = p + scale_color_manual(values = cc)
p = p + scale_color_gradientn(colours = matlab.like2(10))
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/SysDyn/", output, ".pdf"), width = 7, height = 7)

cc = scale_color_gradientn(colours = matlab.like2(100), limits = c(-3,3))
max(abs(range(data_sea$NAO[!is.na(data_sea$NAO)])))


p = ggplot(data_sea) + geom_point(aes(x = dim, y = theta, col = NAO)) + theme_bw()
p = p + cc
print(p)

## each season 
for (seas in 1:4){
  idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
  data_sea = data_NAO[idxdates,]
  
  output = paste0("dim_theta_with_NAO_1999-2017_", season[seas], "_by_", var, "_rainbow")
  p = ggplot(data_sea) + geom_point(aes(x = dim, y = theta, col = NAO)) + theme_bw()
  # cc <- colorRamps::matlab.like2(n = 16)
  # p = p + scale_color_manual(values = cc)
  p = p + scale_color_gradientn(colours = matlab.like2(10))
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/SysDyn/", output, ".pdf"), width = 7, height = 7)
}


#### Contribution ####
library(ggplot2)
library(reshape2)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

casename = c(1:16)
load(file = "~/Documents/LSCE/SWG/NA_slp_sst/case_name.RData")
casename = case_name
ncase = length(casename)

load("~/Documents/LSCE/SWG/NA_slp_sst/cb_table.RData")         # cb_table
load("~/Documents/LSCE/SWG/NA_slp_sst/predictor_name.RData")   # predictor_name

for (case in 1:ncase){
  for (seas in 3){
    
    PC_name = unlist(predictor_name[which(!is.na(cb_table[case,]))])
    
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/1999_2017/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1999_2017_", case, ".RData"))
    par_tt = coef(fit_stations_tt[[1]], matrix = T)
    colnames(par_tt) = c("mean", "sd")
    contrib = as.data.frame(apply(abs(par_tt[-1,]), 2, FUN = function(x){x/sum(x)}))
    contrib$PC = PC_name
    # contrib$type = "rea"
    contrib$type = "1999-2017"
    contrib = melt(contrib, id.vars = c("PC", "type"))
    contrib$PC <- factor(contrib$PC, levels = PC_name)
    contribution = contrib
    
    load(paste0("~/Documents/LSCE/SWG/NA_slp_sst/tmean/1979_1998/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1979_1998_", case, ".RData"))
    par_tt = coef(fit_stations_tt[[1]], matrix = T)
    colnames(par_tt) = c("mean", "sd")
    contrib = as.data.frame(apply(abs(par_tt[-1,]), 2, FUN = function(x){x/sum(x)}))
    contrib$PC = PC_name
    # contrib$type = "nat"
    contrib$type = "1979-1998"
    contrib = melt(contrib, id.vars = c("PC", "type"))
    contrib$PC <- factor(contrib$PC, levels = PC_name)
    contribution = rbind(contribution, contrib)
    
    p = ggplot(contribution)
    p = p + geom_point(aes(x = PC, y = value, col = type)) 
    p = p + facet_wrap(~variable, nrow = 2) + theme_bw()
    p = p + labs(title = paste0("Contribution of principl components for t2m in ", season[seas], " (", casename[case], ")"), 
                 x = "Composant",
                 y = "Percent (X100%)")
    print(p)
    
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Contribution/Contribution_t2m_1999_2017_", season[seas], "_", case, "_", casename[case], ".pdf"), width = 14, height = 7)
  }
}


#### contribution with heatmap : slp ####
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

seas = 3

load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/msl_", SEAS[seas], "_pca.RData"))

pca.loadings <- pca$rotation


for (PC in (1:4)){
  load("~/Documents/LSCE/SWG/coord/LON_LAT_NA_1.5x1.5.RData")
  
  contrib_PC = matrix(pca.loadings[,PC], nrow = length(LAT), ncol = length(LON))
  contrib_PC = t(contrib_PC)
  
  colnames(contrib_PC) = LAT; rownames(contrib_PC) = LON
  contrib_PC = contrib_PC[,rev(seq_along(LAT))]
  
  col.range = c(-max(abs(range(contrib_PC[!is.na(contrib_PC)]))), max(abs(range(contrib_PC[!is.na(contrib_PC)]))))
  cc <- colorRamps::matlab.like2(n = 100)
  output = paste0("contribution of PC_", PC, " for slp in ", season[seas])
  image.plot(contrib_PC, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
             main = output)
  world_map = map('world', plot = FALSE)
  lines(world_map, col = "white")
  contour(x = LON, y = rev(LAT), z = contrib_PC, add = TRUE)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Contribution/", output,".pdf"), width = 11, height = 5)
}


#### contribution with heatmap : sst ####
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

seas = 3

load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/pca/sst_", SEAS[seas], "_pca.RData"))

pca.loadings <- pca$rotation

for (PC in (1:4)){
  load("~/Documents/LSCE/SWG/coord/LON_LAT_NA_1.5x1.5_small.RData")
  load("~/Documents/LSCE/SWG/coord/available_pixels_for_sst_in_NA_1.5x1.5_small.RData")
  
  contrib_PC = matrix(NA, nrow = length(LAT), ncol = length(LON))
  contrib_PC[s] = pca.loadings[,PC]
  contrib_PC = t(contrib_PC)
  
  colnames(contrib_PC) = LAT; rownames(contrib_PC) = LON
  contrib_PC = contrib_PC[,rev(seq_along(LAT))]
  
  col.range = c(-max(abs(range(contrib_PC[!is.na(contrib_PC)]))), max(abs(range(contrib_PC[!is.na(contrib_PC)]))))
  cc <- colorRamps::matlab.like2(n = 100)
  output = paste0("contribution of PC_", PC, " for sst in ", season[seas])
  image.plot(contrib_PC, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
             main = output)
  world_map = map('world', plot = FALSE)
  lines(world_map, col = "white")
  contour(x = LON, y = rev(LAT), z = contrib_PC, add = TRUE)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_slp_sst/Image/Contribution/", output,".pdf"), width = 11, height = 5)
}
