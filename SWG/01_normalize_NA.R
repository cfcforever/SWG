library(lubridate)

#### dim and theta : msl for temperature from 1979 to 2018 in NA ####
load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_dim.RData")
load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_theta.RData")
plot(dim, theta)

load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_dim.RData")
load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_theta.RData")
plot(dim, theta)

# load("/Volumes/Data-ExFAT/nc/FR/msl_1979-2018_FR_dim.RData")
# load("/Volumes/Data-ExFAT/nc/FR/msl_1979-2018_FR_theta.RData")
# range(dim); range(theta)
# plot(dim, theta)

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
  save(SysDyn, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[seas], "_SysDyn.RData"))
}

#### dim and theta : z500 for temperature from 1979 to 2018 in NA ####
# load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_dim.RData")
# load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_theta.RData")
# range(dim); range(theta)

load("/Volumes/Data-ExFAT/nc/FR/z500_1979-2018_FR_dim.RData")
load("/Volumes/Data-ExFAT/nc/FR/z500_1979-2018_FR_theta.RData")
range(dim); range(theta)
plot(dim, theta)

load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_dim.RData")
load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_theta.RData")
plot(dim, theta)

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
  save(SysDyn, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[seas], "_SysDyn.RData"))
}


#### dim and theta : sst for temperature from 1979 to 2018 in NA ####
load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_dim.RData")
load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_theta.RData")
range(dim); range(theta)
plot(dim, theta)

load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_dim.RData")
load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_theta.RData")
plot(dim, theta)

load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_small_dim.RData")
load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_small_theta.RData")
plot(dim, theta)

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
  save(SysDyn, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[seas], "_SysDyn.RData"))
}


#### predictor list ####
# predictor_season_1 : msl PCS only
# predictor_season_2 : msl SysDyn only
# predictor_season_3 : z500 PCS only
# predictor_season_4 : z500 SysDyn only
# predictor_season_5 : sst PCS only
# predictor_season_6 : sst SysDyn only


#### normalize msl PCS : predictor_season_1 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_",SEAS[seas],"_PCS.RData"))
  PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_1.RData"))
}

# for (seas in 1:4){
#   load(paste0("~/Documents/LSCE/SWG/NA/pca/msl_",SEAS[seas],"_PCS.RData"))
#   # PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
#   # print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
#   # colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
#   # PCS = as.data.frame(PCS)
#   save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/predictor_",SEAS[seas],"_1.RData"))
# }

#### normalize msl SysDyn : predictor_season_2 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_",SEAS[seas],"_SysDyn.RData"))
  # print(c(colMeans(SysDyn), apply(SysDyn, 2, FUN = sd)))
  PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_2.RData"))
}


#### normalize z500 PCS : predictor_season_3 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_",SEAS[seas],"_PCS.RData"))
  PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_3.RData"))
}

# for (seas in 1:4){
#   load(paste0("~/Documents/LSCE/SWG/NA/pca/z500_",SEAS[seas],"_PCS.RData"))
#   # PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
#   # print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
#   # colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
#   # PCS = as.data.frame(PCS)
#   save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/predictor_",SEAS[seas],"_3.RData"))
# }


#### normalize z500 SysDyn : predictor_season_4 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_",SEAS[seas],"_SysDyn.RData"))
  PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_4.RData"))
}


#### normalize sst PCS : predictor_season_5 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_",SEAS[seas],"_PCS.RData"))
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_5.RData"))
}

# for (seas in 1:4){
#   load(paste0("~/Documents/LSCE/SWG/NA/pca/sst_",SEAS[seas],"_PCS.RData"))
#   # PCS = apply(PCS, MARGIN = 2, FUN = function(x){scale(x)})
#   # print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
#   # colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
#   # PCS = as.data.frame(PCS)
#   save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/predictor_",SEAS[seas],"_5.RData"))
# }


#### normalize sst SysDyn : predictor_season_6 ####
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_",SEAS[seas],"_SysDyn.RData"))
  PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_6.RData"))
}


#### for the rest predictor_season_XXX: combn(1:6, 2) ####
cb = combn(1:6, 2)
SEAS=c('DJF','MAM','JJA','SON')
for (k in 1:ncol(cb)){
  for (seas in 1:4){
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[1,k], ".RData"))
    predictor = PCS
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[2,k], ".RData"))
    predictor = cbind(predictor, PCS)
    print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
    PCS = predictor
    colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
    PCS = as.data.frame(PCS)
    n = 6 + k
    save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", n, ".RData"))
  }
}


#### combn(1:6, 3) ####
cb = combn(1:6, 3)
SEAS=c('DJF','MAM','JJA','SON')
for (k in 1:ncol(cb)){
  for (seas in 1:4){
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[1,k], ".RData"))
    predictor = PCS
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[2,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[3,k], ".RData"))
    predictor = cbind(predictor, PCS)
    print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
    PCS = predictor
    colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
    PCS = as.data.frame(PCS)
    n = 21 + k
    save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", n, ".RData"))
  }
}

#### combn(1:6, 4) ####
cb = combn(1:6, 4)
SEAS=c('DJF','MAM','JJA','SON')
for (k in 1:ncol(cb)){
  for (seas in 1:4){
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[1,k], ".RData"))
    predictor = PCS
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[2,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[3,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[4,k], ".RData"))
    predictor = cbind(predictor, PCS)
    print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
    PCS = predictor
    colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
    PCS = as.data.frame(PCS)
    n = 41 + k
    save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", n, ".RData"))
  }
}

#### combn(1:6, 5) ####
cb = combn(1:6, 5)
SEAS=c('DJF','MAM','JJA','SON')
for (k in 1:ncol(cb)){
  for (seas in 1:4){
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[1,k], ".RData"))
    predictor = PCS
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[2,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[3,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[4,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[5,k], ".RData"))
    predictor = cbind(predictor, PCS)
    print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
    PCS = predictor
    colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
    PCS = as.data.frame(PCS)
    n = 56 + k
    save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", n, ".RData"))
  }
}

#### combn(1:6, 6) ####
cb = combn(1:6, 6)
SEAS=c('DJF','MAM','JJA','SON')
for (k in 1:ncol(cb)){
  for (seas in 1:4){
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[1,k], ".RData"))
    predictor = PCS
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[2,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[3,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[4,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[5,k], ".RData"))
    predictor = cbind(predictor, PCS)
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", cb[6,k], ".RData"))
    predictor = cbind(predictor, PCS)
    print(c(colMeans(predictor), apply(predictor, 2, FUN = sd)))
    PCS = predictor
    colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
    PCS = as.data.frame(PCS)
    n = 62 + k
    save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_",SEAS[seas],"_", n, ".RData"))
  }
}

