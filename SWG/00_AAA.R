#### usual definitions ####
MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))


#### check nc file ####
nc <- nc_open("~/t2m_1979-01-01to2018-07-31_mean_FR.nc")
# nc <- nc_open("/Volumes/Data-OSX/msl_1979-01-01to2018-07-31_NA.nc")
# nc <- nc_open("/Volumes/Data-OSX/LSCE/nc/slp.1948-2014_NA.nc")
# nc <- nc_open("~/z500_1979-01-01to2018-07-31_1.5_NA.nc")
print(nc)
data = ncvar_get(nc, "z") 
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)
range(data)

nc1 <- nc_open("/Volumes/Data-ExFAT/msl_1979-01-01to2018-07-31_NA.nc")
print(nc1)
data1 = ncvar_get(nc1, "msl") 
lon1 = ncvar_get(nc1, "longitude")
lat1 = ncvar_get(nc1, "latitude")
time1 = ncvar_get(nc1, "time")
nc_close(nc1)

DATE_nc  = atoms(timeSequence(from="1900-01-01",to="2018-12-31",by='hour'))
# head(DATE_OBS); tail(DATE_OBS)
head(DATE_nc[time+1,])
tail(DATE_nc[time+1,])

load("~/msl_1979-2018_FR_dim.RData")
load("~/msl_1979-2018_FR_theta.RData")


#### sst check ####
load(file = "~/Documents/LSCE/SWG/NA_1.5/pca/sst_DJF_pca.RData")
source("fun_search.PCAlevel.R")
search.PCAlevel(pca, 65)


#### create propre data for function of system dynamics ####
var = c("msl", "z500", "sst")
v = c("msl", "z", "sst")

k = 1

nc <- nc_open(paste0("~/", var[k], "_1979-01-01to2018-07-31_1.5_NA_mean.nc"))
print(nc)
data = ncvar_get(nc, v[k]) 
lon = ncvar_get(nc, "lon")
lat = ncvar_get(nc, "lat")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc)
range(data)

DATE_nc  = atoms(timeSequence(from="1900-01-01",to="2018-12-31",by='hour'))
head(DATE_nc[time,])
tail(DATE_nc[time,])


var = c("msl", "z500", "sst")
v = c("msl", "z", "sst")

k = 3

nc1 <- nc_open(paste0("~/", var[k], "_1979-01-01to2018-07-31_1.5_NA.nc"))
data1 = ncvar_get(nc1, v[k]); range(data1)
lon = ncvar_get(nc1, "lon")
lat = ncvar_get(nc1, "lat")
time = ncvar_get(nc1, "time"); length(time); length(time)/4
nt = length(time)
nc_close(nc1)

data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)
save(data, file = paste0("~/", var[k], "_1979-01-01to2018-07-31_mean_NA.RData"))

nx=dim(data1)[2];ny=dim(data1)[1]
nt=dim(data1)[3]
dat1 = data1*NA; dim(dat1) <- c(nt,ny,nx)
for (i in 1:nt) dat1[i,,] <- t(as.matrix(data1[,,i]))
# two dimentions
dim(dat1)=c(nt,nx*ny)
table(apply(dat1, 2, FUN = function(x){sum(!is.na(x))/nt}))

sel_1 = which(apply(dat1, 2, FUN = function(x){sum(!is.na(x))/nt}) == 1)
dat = dat1[,sel_1]
data = (dat[seq(from=1,to=nt,by = 4), ] + dat[seq(from=2,to=nt,by = 4), ] + dat[seq(from=3,to=nt,by = 4), ] + dat[seq(from=4,to=nt,by = 4), ])/4
dim(data)
range(data)
save(data, file = paste0("~/", var[k], "_1979-01-01to2018-07-31_mean_NA.RData"))



load("~/msl_1979-01-01to2018-07-31_mean_NA.RData"); dim(data); range(data)
load("~/z500_1979-01-01to2018-07-31_mean_NA.RData"); dim(data); range(data)
load("~/sst_1979-01-01to2018-07-31_mean_NA.RData"); dim(data); range(data)


#### rotation ####
load(paste0("~/Documents/LSCE/SWG/NA/pca/msl_DJF_pca.RData"))
rot = pca$rotation

r = array(rot, dim = c(175, 65, dim(rot)[2]))



SEAS=c('DJF','MAM','JJA','SON')
for (seas in 1:4){
  load(paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[seas], "_pca.RData"))
  
  rot = pca$rotation
  
  r = array(rot, dim = c(175, 65, dim(rot)[2]))
  
  #create dimension: time, longitude, latitude
  longitude <- ncdim_def( name = 'longitude', units = 'degrees_east', vals = LON)
  
  latitude <- ncdim_def( name = 'latitude', units = 'degrees_north', vals = LAT)
  
  t <- ncdim_def( name = 'time', units = 'day', vals = 1:dim(rot)[2] )
  
  # create variable: rotation
  rotation <- ncvar_def( name = 'rotation', units = '', dim = list(longitude,latitude,t), missval = NA, prec = 'double' )
  
  # create .nc - only one variable
  ncnew <- nc_create( filename = paste0('~/rotation_sst_', SEAS[seas], '.nc'), vars = rotation)
  
  # write values
  ncvar_put( nc = ncnew, varid = rotation, vals = r )
  
  # write attribution
  # ncatt_put( nc = ncnew, varid = 0, attname = 'description', attval = 'sst data in NINO3 area during 2009')
  
  # close .nc
  nc_close(ncnew)
}

# check .nc
nc <- nc_open('rotation.nc')
print(nc)
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
nc_close(nc)


#### dimension et theta ####
# NA - 1.5x1.5 - 1979-2018
load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_dim.RData")
load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_theta.RData")
plot(dim, theta)

load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_dim.RData")
load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_theta.RData")
plot(dim, theta)

load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_dim.RData")
load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_theta.RData")
plot(dim, theta)
load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_small_dim.RData")
load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_small_theta.RData")
plot(dim, theta)



load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_dim.RData")
load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_theta.RData")
plot(dim, theta)

load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_dim.RData")
load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_theta.RData")
plot(dim, theta)

load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_dim.RData")
load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_theta.RData")
plot(dim, theta)

load("~/z500_1979-2018_FR_dim.RData")
load("~/z500_1979-2018_FR_theta.RData")
plot(dim, theta)

load("~/msl_1979-2018_FR_dim.RData")
load("~/msl_1979-2018_FR_theta.RData")
plot(dim, theta)


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
save(mean, sd, DATE, file = "~/Documents/LSCE/SWG/NA/t2m_mean_sd_1999_2017_statio_nat.RData")


#### select the events ####
load(file = "~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

# winter events
seas = 1
idxdatas = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])

t_seas = tmean[idxdatas,]
t_seas = cbind(t_seas, DATE[idxdatas,])
winter_hot  = t_seas[t_seas[,1] %in% tail(sort(t_seas[,1], decreasing = F), 10),]
winter_cold = t_seas[t_seas[,1] %in% head(sort(t_seas[,1], decreasing = F), 10),]

# summer events
seas = 3
idxdatas = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])

t_seas = tmean[idxdatas,]
t_seas = cbind(t_seas, DATE[idxdatas,])
summer_hot  = t_seas[t_seas[,1] %in% tail(sort(t_seas[,1], decreasing = F), 10),]
summer_cold = t_seas[t_seas[,1] %in% head(sort(t_seas[,1], decreasing = F), 10),]


