library(lubridate)

#### dim and theta : geop_500 for temperature from 1979 to 2017 in NA ####
load("~/Documents/LSCE/Dynamic System/geop_500.1979-2017_NA_dim.RData")
load("~/Documents/LSCE/Dynamic System/geop_500.1979-2017_NA_theta.RData")


## separate in 4 seasons
ndays = as.numeric(difftime(as.Date("2017-12-31"), as.Date("1979-01-01"), units = "days"))
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
  save(SysDyn, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[seas], "_SysDyn.RData"))
}


#### 3 cases ####
# case_1: only PCS
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/t2m/geop_",SEAS[seas],"_PCS.RData"))
  PCS = apply(PCS, MARGIN = 2, FUN = function(x){(x-min(x))/(max(x)-min(x))})
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/t2m/case_1/geop_",SEAS[seas],"_PCS.RData"))
}


# case_2: only SysDyn
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/t2m/geop_",SEAS[seas],"_SysDyn.RData"))
  PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){(x-min(x))/(max(x)-min(x))})
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/t2m/case_2/geop_",SEAS[seas],"_SysDyn.RData"))
}


# case_3: PCS + SysDyn
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  
  dim_sea   = dim[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  theta_sea = theta[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  
  SysDyn = cbind(dim_sea, theta_sea)
  
  load(paste0("~/Documents/LSCE/SWG/t2m/geop_",SEAS[seas],"_PCS.RData"))
  
  PCS_SysDyn = cbind(PCS, SysDyn)
  
  PCS = apply(PCS_SysDyn, MARGIN = 2, FUN = function(x){(x-min(x))/(max(x)-min(x))})
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  
  save(PCS, file = paste0("~/Documents/LSCE/SWG/t2m/case_3/geop_", SEAS[seas], "_PCS_SysDyn.RData"))
}


#### dim and theta : slp for precipitation from 1979 to 2017 in NA ####
load("~/Documents/LSCE/Dynamic System/slp_1979-2017_NA_dim.RData")
load("~/Documents/LSCE/Dynamic System/slp_1979-2017_NA_theta.RData")


## separate in 4 seasons
ndays = as.numeric(difftime(as.Date("2017-12-31"), as.Date("1979-01-01"), units = "days"))
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
  save(SysDyn, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[seas], "_SysDyn.RData"))
}

#### 3 cases ####
# case_1: only PCS
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/tp/slp_",SEAS[seas],"_PCS.RData"))
  PCS = apply(PCS, MARGIN = 2, FUN = function(x){(x-min(x))/(max(x)-min(x))})
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/tp/case_1/slp_",SEAS[seas],"_PCS.RData"))
}


# case_2: only SysDyn
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/tp/slp_",SEAS[seas],"_SysDyn.RData"))
  PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){(x-min(x))/(max(x)-min(x))})
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/tp/case_2/slp_",SEAS[seas],"_SysDyn.RData"))
}


# case_3: PCS + SysDyn
SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  
  dim_sea   = dim[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  theta_sea = theta[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  
  SysDyn = cbind(dim_sea, theta_sea)
  
  load(paste0("~/Documents/LSCE/SWG/tp/slp_",SEAS[seas],"_PCS.RData"))
  
  PCS_SysDyn = cbind(PCS, SysDyn)
  
  PCS = apply(PCS_SysDyn, MARGIN = 2, FUN = function(x){(x-min(x))/(max(x)-min(x))})
  colnames(PCS) = paste0("PC", 1:ncol(PCS)) # for 90% variance
  PCS = as.data.frame(PCS)
  
  save(PCS, file = paste0("~/Documents/LSCE/SWG/tp/case_3/slp_", SEAS[seas], "_PCS_SysDyn.RData"))
}

