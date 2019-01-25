####/!\ ATTENTION: this code is not finished !!! ####

library(ncdf4)
library(timeDate)

#### load nc ####
var = "t2m"

if (var == "t2m"){
  nc <- nc_open("~/Documents/LSCE/nc/temp_2m.1979-2017_FR.nc")
  print(nc) # t2m[longitude,latitude,time]
  data = ncvar_get(nc, "t2m")  
  
  time = ncvar_get(nc, "time")
  nt = length(time)
  
  dat = data[,,seq(from = 1, by = 4, to = nt)] + data[,,seq(from = 2, by = 4, to = nt)] +
    data[,,seq(from = 3, by = 4, to = nt)] + data[,,seq(from = 4, by = 4, to = nt)]
  dat = dat/4 - 273.15
  
  nt = dim(dat)[3]
}

if (var == "tp"){
  nc <- nc_open("~/Documents/LSCE/nc/precip_total.1979-2017_FR.nc")
  print(nc)
  data = ncvar_get(nc, "tp")
  
  time = ncvar_get(nc, "time")
  nt = length(time)
  
  dat = data[,,seq(from = 1, by = 2, to = nt)] + data[,,seq(from = 2, by = 2, to = nt)]
  dat = dat * 10^3
  
  nt = dim(dat)[3]
}


#### create LAT and LON ####
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")

nlon = length(LON)
nlat = length(LAT)


#### create DATE ####
DATE = atoms(timeSequence(from="1979-01-01", to="2017-12-31", by="day"))


#### create grids ####
grids = as.data.frame(matrix(NA, nrow = nlat*nlon+1, ncol = 4))
grids[,1] = c("grid_id", 1:(nlat*nlon))
grids[,2] = c("location", paste0("grid_", 1:(nlat*nlon)))
grids[,3] = c("longitude", rep(LON, each = nlat))
grids[,4] = c("latitude", rep(LAT, nlon))
View(grids)


#### create grids_var ####
grids_var = dat*NA; dim(grids_var) <- c(nt, nlat, nlon)
for (i in 1:nt){ grids_var[i,,] = t(dat[,,i]) }
dim(grids_var) = c(nt, nlat*nlon)
colnames(grids_var) = paste0("grid_", 1:(nlat*nlon))
View(grids_var)

# give a proper name for grids_var
grids_tp = grids_var


#### save .RData ####
save.image(file = paste0("~/Documents/LSCE/SWG/", var, "/gridobs_", var, ".RData"))


#### clean up ####
nc_close(nc); rm(nc)
rm(nlat, nlon, nt, i, var)
rm(data, time, dat, grids_var)



