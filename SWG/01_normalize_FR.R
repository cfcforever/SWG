library(lubridate)

#### dim and theta : msl for temperature from 1979 to 2018 in FR ####
load("~/msl_1979-2018_FR_dim.RData")
load("~/msl_1979-2018_FR_theta.RData")


## separate in 4 seasons
ndays = as.numeric(difftime(as.Date("2018-07-31"), as.Date("1979-01-01"), units = "days"))
DATE = as.Date((0:ndays), origin = "1979-01-01") ; head(DATE); tail(DATE)

MONTH = month(DATE)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

for (seas in 1:4){
  dim_sea   = dim[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  theta_sea = theta[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  
  SysDyn = cbind(dim_sea, theta_sea)
  save(SysDyn, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[seas], "_SysDyn.RData"))
}

#### dim and theta : z500 for temperature from 1979 to 2018 in FR ####
load("~/z500_1979-2018_FR_dim.RData")
load("~/z500_1979-2018_FR_theta.RData")
range(dim); range(theta)

## separate in 4 seasons
ndays = as.numeric(difftime(as.Date("2018-07-31"), as.Date("1979-01-01"), units = "days"))
DATE = as.Date((0:ndays), origin = "1979-01-01") ; head(DATE); tail(DATE)

MONTH = month(DATE)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

for (seas in 1:4){
  dim_sea   = dim[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  theta_sea = theta[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  
  SysDyn = cbind(dim_sea, theta_sea)
  save(SysDyn, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[seas], "_SysDyn.RData"))
}

#### predictor list ####
# predictor_season_1 : msl PCS only
# predictor_season_2 : msl SysDyn only
# predictor_season_3 : z500 PCS only
# predictor_season_4 : z500 SysDyn only
# predictor_season_5 : 1 + 2
# predictor_season_6 : 1 + 3
# predictor_season_7 : 1 + 4
# predictor_season_8 : 2 + 3
# predictor_season_9 : 2 + 4
# predictor_season_10 : 3 + 4
# predictor_season_11 : 1 + 2 + 3
# predictor_season_12 : 1 + 2 + 4
# predictor_season_13 : 1 + 3 + 4
# predictor_season_14 : 2 + 3 + 4
# predictor_season_15 : 1 + 2 + 3 + 4

#### normalize msl PCS : predictor_season_1 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/msl_",SEAS[seas],"_PCS.RData"))
  PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_1.RData"))
}

#### normalize msl SysDyn : predictor_season_2 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/msl_",SEAS[seas],"_SysDyn.RData"))
  PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_2.RData"))
}


#### normalize z500 PCS : predictor_season_3 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/z500_",SEAS[seas],"_PCS.RData"))
  PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_3.RData"))
}


#### normalize z500 SysDyn : predictor_season_4 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/z500_",SEAS[seas],"_SysDyn.RData"))
  PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_4.RData"))
}

#### predictor_season_5 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_1.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_2.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_5.RData"))
}

#### predictor_season_6 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_1.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_3.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_6.RData"))
}

#### predictor_season_7 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_1.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_4.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_7.RData"))
}

#### predictor_season_8 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_2.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_3.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_8.RData"))
}

#### predictor_season_9 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_2.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_4.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_9.RData"))
}

#### predictor_season_10 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_3.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_4.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_10.RData"))
}

#### predictor_season_11 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_1.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_2.RData"))
  predictor = cbind(predictor, PCS)
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_3.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_11.RData"))
}

#### predictor_season_12 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_1.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_2.RData"))
  predictor = cbind(predictor, PCS)
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_4.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_12.RData"))
}

#### predictor_season_13 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_1.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_3.RData"))
  predictor = cbind(predictor, PCS)
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_4.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_13.RData"))
}

#### predictor_season_14 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_2.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_3.RData"))
  predictor = cbind(predictor, PCS)
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_4.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_14.RData"))
}

#### predictor_season_15 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_1.RData"))
  predictor = PCS
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_2.RData"))
  predictor = cbind(predictor, PCS)
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_3.RData"))
  predictor = cbind(predictor, PCS)
  load(paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_4.RData"))
  predictor = cbind(predictor, PCS)
  print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
  PCS = predictor
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/predictor_",SEAS[seas],"_15.RData"))
}