# ReadMe ------------------------------------------------------------------
#  msl_1979-01-01to2018-07-31_NA_1.5x1.5.nc
#: mean_sea_level_pressure extracted from ECMWF with 1.5° (space) and 6h (time) resolution
#
#  z500_1979-01-01to2018-07-31_NA_1.5x1.5.nc
#: 500 hPa geopotential extracted from ECMWF with 1.5° (space) and 6h (time) resolution
#
#  sst_1979-01-01to2018-07-31_NA_1.5x1.5.nc
#: sea_surface_temperature extracted from ECMWF with 1.5° (space) and 6h (time) resolution


# Load packages -----------------------------------------------------------
library(ncdf4)
library(timeDate)


# msl_1979-01-01to2018-07-31_NA_1.5x1.5.nc --------------------------------
nc <- nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/msl_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
data1 = ncvar_get(nc, "msl"); range(data1)
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc) 

nt = length(time)
data =   (data1[,,seq(from = 1, to = nt, by = 4)] + 
            data1[,,seq(from = 2, to = nt, by = 4)] + 
            data1[,,seq(from = 3, to = nt, by = 4)] + 
            data1[,,seq(from = 4, to = nt, by = 4)])/4
dim(data)
range(data)

DATE = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

# save(data, LON, LAT, DATE, file = "/Volumes/Data-ExFAT/DATA/RData/NA_1.5/msl_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")
# load("/Volumes/Data-ExFAT/DATA/RData/NA_1.5/msl_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")



# z500_1979-01-01to2018-07-31_NA_1.5x1.5.nc -------------------------------
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

DATE = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

# save(data, LON, LAT, DATE, file = "/Volumes/Data-ExFAT/DATA/RData/NA_1.5/z500_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")
# load(file = "/Volumes/Data-ExFAT/DATA/RData/NA_1.5/z500_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")



# sst_1979-01-01to2018-07-31_NA_1.5x1.5.nc --------------------------------
# nc <- nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/sst_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
nc <- nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/sst_1979-01-01to2018-07-31_NA_1.5x1.5_small.nc")
data1 = ncvar_get(nc, "sst"); range(data1)
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc) 

nt = length(time)
data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)

DATE = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

# save(data, LON, LAT, DATE, file = "/Volumes/Data-ExFAT/DATA/RData/NA_1.5/sst_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")
# load(file = "/Volumes/Data-ExFAT/DATA/RData/NA_1.5/sst_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")

# save(data, LON, LAT, DATE, file = "/Volumes/Data-ExFAT/DATA/RData/NA_1.5/sst_1979-01-01to2018-07-31_NA_1.5x1.5_daily_small.RData")
# load(file = "/Volumes/Data-ExFAT/DATA/RData/NA_1.5/sst_1979-01-01to2018-07-31_NA_1.5x1.5_daily_small.RData")
