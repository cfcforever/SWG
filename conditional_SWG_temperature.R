#### load packages -----------------------------------------------------------
source("function/load_packages.R")
####



#### common variables --------------------------------------------------------
SEAS   = c('DJF','MAM','JJA','SON')
season = c('Winter','Spring','Summer','Fall')
MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
####



#### load slp SysDyn ---------------------------------------------------------
if (!dir.exists("output")){
  dir.create("output")
  dir.create("output/predictor")
  for (seas in 1:4){
    load(paste0("DATA/msl_",SEAS[seas],"_SysDyn.RData"))
    PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){scale(x)})
    print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
    colnames(PCS) = paste0("PC", 1:ncol(PCS))
    PCS = as.data.frame(PCS)
    save(PCS, file = paste0("output/predictor/predictor_",SEAS[seas],"_1.RData"))
  }
}
####



#### Create directoires for output -------------------------------------------
if (!dir.exists(paste0("output/", city))){
  dir.create(paste0("output/", city))
}
if (!dir.exists(paste0("output/", city, "/SIMU"))){
  dir.create(path = paste0("output/", city, "/SIMU"))
  dir.create(path = paste0("output/", city, "/SIMU/nat"))
  dir.create(path = paste0("output/", city, "/SIMU/rea"))
}
if (!dir.exists(paste0("output/", city, "/Image"))){
  dir.create(path = paste0("output/", city, "/Image"))
}
####



#### Estimation for mean and sd ----------------------------------------------
source("function/fun_estimation_t2m.R")

## 1979 - 1998
fun_estimation(predictor = paste0('output/predictor/predictor_'),
               NUM = 1,
               input.tmean = paste0('DATA/tmean_', city, '_1979_2017.RData'),
               output.name = paste0('output/', city, '/SWG_ERAI_ESD_tmean_'),
               year.begin = 1979,
               year.end = 1998)

## 1999 - 2017
fun_estimation(predictor = paste0('output/predictor/predictor_'),
               NUM = 1,
               input.tmean = paste0('DATA/tmean_', city, '_1979_2017.RData'),
               output.name = paste0('output/', city, '/SWG_ERAI_ESD_tmean_'),
               year.begin = 1999,
               year.end = 2017)
####



#### Simulation for mean and sd ----------------------------------------------
source("function/fun_simulation_t2m.R")

## 1979 - 1998
fun_simulation(predictor = paste0('output/predictor/predictor_'),
               parameter = paste0('output/', city, '/SWG_ERAI_ESD_tmean_'),
               NUM = 1, 
               type = "nat",
               output = paste0('output/', city, '/'), 
               y1 = 1999, y2 = 2017,
               year.begin = 1979, year.end = 1998)

## 1999 - 2017
fun_simulation(predictor = paste0('output/predictor/predictor_'),
               parameter = paste0('output/', city, '/SWG_ERAI_ESD_tmean_'),
               NUM = 1, 
               type = "rea",
               output = paste0('output/', city, '/'), 
               y1 = 1999, y2 = 2017,
               year.begin = 1999, year.end = 2017)
####



#### Estimate stationary model -----------------------------------------------
load(paste0("DATA/tmean_", city, "_1979_2017.RData"))
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

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
save(mean, sd, DATE, file = paste0("output/", city, "/t2m_mean_sd_1999_2017_nat_0.RData"))

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
save(mean, sd, DATE, file = paste0("output/", city, "/t2m_mean_sd_1999_2017_rea_0.RData"))
####



#### SAMPLE 100 simulations --------------------------------------------------
DATE = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))

for (w in c("nat", "rea")){
  for (k in 0:1){ # 0: stationary // 1: conditional - slp_SD
    load(paste0("output/", city, "/t2m_mean_sd_1999_2017_", w, "_", k, ".RData"))
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
      
      filesimu = paste('output/', city, '/SIMU/', w, '/SIMU_SWG_tmean_1999_2017_', w, '_', k, '_run_', NUM, '.RData', sep="")
      save(Sample,DATE,file = filesimu)
    }
  }
}
####



#### separation in season ####
for (k in 0:1){ # 0: stationary // 1: conditional - slp_SD
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(paste0("output/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_", k, "_run_", NUM, ".RData"))
    for (seas in 1:4){
      Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
      DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
      save(DATE_seas, Sample_seas, file = paste0("output/", city, "/SIMU/nat/SIMU_tmean_1999_2017_nat_", k, "_", SEAS[seas], "_run_", NUM, ".RData"))
    }
  }
}
####

