library(ncdf4)
library(timeDate)

source("fun_search.PCAlevel.R")

#### calculate eigenvalue of pca
# library("factoextra")
# load(file = "~/Documents/LSCE/SWG/NA_slp_sst/pca/sst_JJA_pca.RData")
# eig.val <- get_eigenvalue(pca)
# eig.val

#### load msl_1979-01-01to2018-07-31_NA_1.5x1.5.nc ####
nc <- nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/msl_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
data1 = ncvar_get(nc, "msl"); range(data1)
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
# save(LON, LAT, file = "~/Documents/LSCE/SWG/coord/LON_LAT_NA_1.5x1.5.RData")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc)

nt = length(time)
data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)
# save(data, file = "/Volumes/Data-ExFAT/msl_1979-01-01to2018-07-31_mean_NA.RData")

nlon = length(LON)
nlat = length(LAT)
nt = dim(data)[3]

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(data1, data, dat)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/msl_", SEAS[4], "_PCS.RData"))



#### load z500_1979-01-01to2018-07-31_NA_1.5x1.5.nc ####
nc <- nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/z500_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
data1 = ncvar_get(nc, "z"); range(data1)
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc)

nt = length(time)
data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)
# save(data, file = "/Volumes/Data-ExFAT/z500_1979-01-01to2018-07-31_mean_NA.RData")

nlon = length(LON)
nlat = length(LAT)
nt = dim(data)[3]

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(data1, dat, data)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/z500_", SEAS[4], "_PCS.RData"))


#### load sst_1979-01-01to2018-07-31_NA_1.5x1.5.nc ####
# nc <- nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/sst_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
nc <- nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/sst_1979-01-01to2018-07-31_NA_1.5x1.5_small.nc")
data1 = ncvar_get(nc, "sst"); range(data1)
LON = ncvar_get(nc, "longitude") 
LAT = ncvar_get(nc, "latitude")
# save(LON, LAT, file = "~/Documents/LSCE/SWG/coord/LON_LAT_NA_1.5x1.5_small.RData")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc)

nt = length(time)
data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)
# save(data, file = "/Volumes/Data-ExFAT/sst_1979-01-01to2018-07-31_mean_NA.RData")
# save(data, file = "/Volumes/Data-ExFAT/sst_1979-01-01to2018-07-31_mean_NA_small.RData")

nlon = length(LON)
nlat = length(LAT)
nt = dim(data)[3]

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(data1, dat, data)

s = which(apply(w_dat, 2, FUN = function(x){sum(!is.na(x))/nt})==1)
# save(s, file = "~/Documents/LSCE/SWG/coord/available_pixels_for_sst_in_NA_1.5x1.5_small.RData")
w_dat = w_dat[,s]
range(w_dat)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA_1.5/pca/sst_", SEAS[4], "_PCS.RData"))




#### load msl_1979-01-01to2018-07-31_mean_NA.nc ####
nc <- nc_open("~/msl_1979-01-01to2018-07-31_NA.nc")
data1 = ncvar_get(nc, "msl"); range(data1)
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc)

# nc <- nc_open("/Volumes/Data-ExFAT/nc/NA/msl_1979-01-01to2018-07-31_mean_NA.nc")
# data1 = ncvar_get(nc, "msl"); range(data1)
# LON = ncvar_get(nc, "longitude")
# LAT = ncvar_get(nc, "latitude")
# nc_close(nc)

nt = length(time)
data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)
# save(data, file = "/Volumes/Data-ExFAT/msl_1979-01-01to2018-07-31_mean_NA.RData")

nlon = length(LON)
nlat = length(LAT)
nt = dim(data)[3]

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}

# rm(data, dat)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[4], "_PCS.RData"))



#### load z500_1979-01-01to2018-07-31_mean_NA.nc ####
nc <- nc_open("~/z500_1979-01-01to2018-07-31_NA.nc")
data1 = ncvar_get(nc, "z"); range(data1)
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc)

nt = length(time)
data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)
# save(data, file = "/Volumes/Data-ExFAT/z500_1979-01-01to2018-07-31_mean_NA.RData")

nlon = length(LON)
nlat = length(LAT)
nt = dim(data)[3]

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(dat, data)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[4], "_PCS.RData"))


#### load sst_1979-01-01to2018-07-31_mean_NA.nc ####
nc <- nc_open("~/sst_1979-01-01to2018-07-31_NA.nc")
data1 = ncvar_get(nc, "sst"); range(data1)
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc)

nt = length(time)
data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)
# save(data, file = "/Volumes/Data-ExFAT/sst_1979-01-01to2018-07-31_mean_NA.RData")

nlon = length(LON)
nlat = length(LAT)
nt = dim(data)[3]

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}

s = which(apply(w_dat, 2, FUN = function(x){sum(!is.na(x))/nt})==1)
w_dat = w_dat[,s]
range(w_dat)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = F)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/sst_", SEAS[4], "_PCS.RData"))



#### load msl_1979-01-01to2018-07-31_mean_NA.nc ####
nc <- nc_open("~/msl_1979-01-01to2018-07-31_mean_NA.nc")
print(nc)
data = ncvar_get(nc, "msl") 
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)

nlon = length(LON)
nlat = length(LAT)
nt = length(time)

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(dat, data)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/msl_", SEAS[4], "_PCS.RData"))



#### load z500_1979-01-01to2018-07-31_mean_NA.nc ####
nc <- nc_open("~/z500_1979-01-01to2018-07-31_mean_NA.nc")
print(nc)
data = ncvar_get(nc, "z") 
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)

nlon = length(LON)
nlat = length(LAT)
nt = length(time)

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(dat, data)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/NA/pca/z500_", SEAS[4], "_PCS.RData"))



#### load geop.nc ####
nc <- nc_open("~/Documents/LSCE/nc/geop_500.1979-2017_NA.nc")
print(nc)
data = ncvar_get(nc, "z") 
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)

# LON = LON[which(LON>=-5.25 & LON<=15.00)]
# LAT = LAT[which(LAT<=55.50 & LAT>=34.50)]

# data = data[which(LON>=-5.25 & LON<=15.00), which(LAT<=55.50 & LAT>=34.50),]

nlon = length(LON)
nlat = length(LAT)
nt = length(time)

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = T)
# pca3 = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = T, rank. = 3)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[1], "_pca.RData"))
# load(file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:20]
save(PCS, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[2], "_pca.RData"))
# load(file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:20]
save(PCS, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[3], "_pca.RData"))
# load(file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:27]
save(PCS, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[4], "_pca.RData"))
# load(file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:20]
save(PCS, file = paste0("~/Documents/LSCE/SWG/t2m/geop_", SEAS[4], "_PCS.RData"))


#### load slp.nc ####
nc <- nc_open("~/Documents/LSCE/nc/slp.1979-2017_NA.nc")
print(nc)
data = ncvar_get(nc, "msl") 
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)

# LON = LON[which(LON>=-5.25 & LON<=15.00)]
# LAT = LAT[which(LAT<=55.50 & LAT>=34.50)]

# data = data[which(LON>=-5.25 & LON<=15.00), which(LAT<=55.50 & LAT>=34.50),]

nlon = length(LON)
nlat = length(LAT)
nt = length(time)

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(dat, data)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[1], "_pca.RData"))
# load(file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:21]
save(PCS, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[2], "_pca.RData"))
# load(file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:24]
save(PCS, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[3], "_pca.RData"))
# load(file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:32]
save(PCS, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[4], "_pca.RData"))
# load(file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:24]
save(PCS, file = paste0("~/Documents/LSCE/SWG/tp/slp_", SEAS[4], "_PCS.RData"))

#### load msl_1979-01-01to2018-07-31_mean_FR.nc ####
nc <- nc_open("~/msl_1979-01-01to2018-07-31_mean_FR.nc")
print(nc)
data = ncvar_get(nc, "msl") 
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)

nlon = length(LON)
nlat = length(LAT)
nt = length(time)

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(dat, data)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/msl_", SEAS[4], "_PCS.RData"))



#### load z500_1979-01-01to2018-07-31_mean_FR.nc ####
nc <- nc_open("~/z500_1979-01-01to2018-07-31_mean_FR.nc")
print(nc)
data = ncvar_get(nc, "z") 
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)

nlon = length(LON)
nlat = length(LAT)
nt = length(time)

dat = data*NA; dim(dat) <- c(nt, nlat, nlon)
dim(dat)

for (i in 1:nt){ 
  dat[i,,] = t(data[,,i]) 
}
dim(dat) = c(nt, nlat*nlon)

pond = 1/sqrt(cos(LAT*pi/180))
scale = rep(pond, length(LON))

w_dat = dat*NA
for(pix in 1:ncol(w_dat)){
  w_dat[,pix] = dat[,pix]/scale[pix]
}
rm(dat, data)

SEAS=c('DJF','MAM','JJA','SON')
DATE_ERAI = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(12,1,2), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[1], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[1], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(3,4,5), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[2], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[2], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(6,7,8), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[3], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[3], "_PCS.RData"))

pca = prcomp(w_dat[DATE_ERAI$m%in%c(9,10,11), ], retx = T, scale = T)
search.PCAlevel(pca, 90)
save(pca, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[4], "_pca.RData"))
PCS = as.data.frame(pca$x)[,1:2]
save(PCS, file = paste0("~/Documents/LSCE/SWG/pca/z500_", SEAS[4], "_PCS.RData"))

