library(ncdf4)
library(fields)
#

#### load msl ####
nc <- nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/msl_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
data1 = ncvar_get(nc, "msl"); range(data1)
LON = ncvar_get(nc, "longitude")
LAT = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time"); length(time); length(time)/4
nc_close(nc) 

nt = length(time)
data = (data1[,,seq(from=1,to=nt,by = 4)] + data1[,,seq(from=2,to=nt,by = 4)] + data1[,,seq(from=3,to=nt,by = 4)] + data1[,,seq(from=4,to=nt,by = 4)])/4
dim(data)
range(data)

DATE = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))

save(data, LON, LAT, DATE, file = "~/Documents/LSCE/SWG/NA_1.5/msl_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")

#### load z500 ####
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

save(data, LON, LAT, DATE, file = "~/Documents/LSCE/SWG/NA_1.5/z500_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")

#### load sst ####
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

# save(data, LON, LAT, DATE, file = "~/Documents/LSCE/SWG/NA_1.5/sst_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")
save(data, LON, LAT, DATE, file = "~/Documents/LSCE/SWG/NA_1.5/sst_1979-01-01to2018-07-31_NA_1.5x1.5_daily_small.RData")


#### Create anomaly.RData ####
var = "msl"; # var = "sst"

input.file = paste0("~/Documents/LSCE/SWG/NA_1.5/", var, "_1979-01-01to2018-07-31_NA_1.5x1.5_daily.RData")
# input.file = paste0("~/Documents/LSCE/SWG/NA_1.5/", var, "_1979-01-01to2018-07-31_NA_1.5x1.5_daily_small.RData")
load(file = input.file)

check_list = DATE[DATE$Y==2000,c("m","d")]

data_mean = array(NaN, c(dim(data)[1:2], 366))
data_sd   = array(NaN, c(dim(data)[1:2], 366))
for (t in 1:366){
  time_select = which(DATE$Y>=1999 & DATE$Y<=2017 & DATE$m==check_list[t,"m"] & DATE$d==check_list[t,"d"])
  data_mean[,,t] = apply(data[,,time_select], c(1,2), mean)
  data_sd[,,t] = apply(data[,,time_select], c(1,2), sd)
}
range(data_mean[!is.na(data_mean)]); range(data_sd[!is.na(data_sd)]); sum(!is.na(data_mean)); sum(!is.na(data_sd))

data_anomaly = array(NaN, dim(data))
for (k in 1:dim(data)[3]){
  t = which(check_list$m==DATE$m[k] & check_list$d==DATE$d[k])
  data_anomaly[,,k] = (data[,,k] - data_mean[,,t])
  # data_anomaly[,,k] = (data[,,k] - data_mean[,,t])/data_sd[,,t]
}
save(data_anomaly, LON, LAT, DATE, file = "~/Documents/LSCE/SWG/NA_slp_sst/slp_anomaly_1999_2017.RData")
# save(data_anomaly, LON, LAT, DATE, file = "~/Documents/LSCE/SWG/NA_1.5/sst_anomaly.RData")

#### RUN code ####
load(file = "~/Documents/LSCE/SWG/NA_1.5/sst_anomaly.RData")
load(file = "~/Documents/LSCE/SWG/NA_1.5/tmean/d_intensity.RData")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

for (model in c(25, 6, 15, 2, 11, 7, 1, 20, 36)){
  tp = which(d[,model] - d[,64]>0); length(tp)
  tn = which(d[,model] - d[,64]<0); length(tn)
  
  data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
  DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]
  
  for (seas in c(1:3)){
    idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
    tp_sea = intersect(tp, idxdates)
    data = apply(data_anomaly[,,tp_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    output = paste0("sst_positive_1999-2017_", season[seas], "_model_", model)
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", 
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/Image/Anomaly/", output,".pdf"), width = 11, height = 5)
    
    
    tn_sea = intersect(tn, idxdates)
    data = apply(data_anomaly[,,tn_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    
    output = paste0("sst_negative_1999-2017_", season[seas], "_model_", model)
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", 
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/Image/Anomaly/", output,".pdf"), width = 11, height = 5)
  }
}


#### draft ####
ap = as.data.frame(apply(data_anomaly[,,tp], 1:2, mean)); colnames(ap) = LAT; rownames(ap) = LON
an = as.data.frame(apply(data_anomaly[,,tn], 1:2, mean)); colnames(an) = LAT; rownames(an) = LON

plot_anomaly <- function(data, lon, lat, output){
  colnames(data) = lat; rownames(data) = lon
  data = data[,rev(seq_along(lat))]
  image.plot(data, x = lon, y = rev(lat), xlab = "Lon", ylab = "Lat", 
             main = "slp_positive_1999-2017_model_1")
  world_map = map('world', plot = FALSE)
  lines(world_map, col = "white")
  contour(x = lon, y = rev(lat), z = data, add = TRUE)
  
  dev.print(pdf, file=output, width = 11, height = 5)
}

plot_anomaly(data = apply(data_anomaly[,,tp], 1:2, mean), 
             lon = LON, lat = LAT,
             output = paste0("~/Documents/LSCE/SWG/NA_1.5/Image/Anomaly/slp_positive_1999-2017_model_1.pdf"))


#### sign test ####
load(file = "~/Documents/LSCE/SWG/NA_1.5/sst_anomaly.RData")
load(file = "~/Documents/LSCE/SWG/NA_1.5/d_intensity.RData")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

model = 7
tp = which(d[,model] - d[,64]>0); length(tp)
tn = which(d[,model] - d[,64]<0); length(tn)

data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]

for (seas in 1:4){
  idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
  
  # tp_sea
  tp_sea = intersect(tp, idxdates)
  
  nb_sign_positive = apply(data_anomaly[,,tp_sea], MARGIN = c(1:2), FUN = function(x){sum(x>0)})
  nb_sign_negative = apply(data_anomaly[,,tp_sea], MARGIN = c(1:2), FUN = function(x){sum(x<0)})
  
  data = nb_sign_positive/length(tp_sea)
  colnames(data) = LAT; rownames(data) = LON
  data = data[,rev(seq_along(LAT))]
  output = paste0("sst_positive_1999-2017_", season[seas], "_model_", model, "_avec_sign_positive")
  image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", 
             main = output)
  world_map = map('world', plot = FALSE)
  lines(world_map, col = "white")
  contour(x = LON, y = rev(LAT), z = data, levels = 2/3, add = TRUE)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/Image/Anomaly/", output,".pdf"), width = 11, height = 5)
  
  data = nb_sign_negative/length(tp_sea)
  colnames(data) = LAT; rownames(data) = LON
  data = data[,rev(seq_along(LAT))]
  output = paste0("sst_positive_1999-2017_", season[seas], "_model_", model, "_avec_sign_negative")
  image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", 
             main = output)
  world_map = map('world', plot = FALSE)
  lines(world_map, col = "white")
  contour(x = LON, y = rev(LAT), z = data, levels = 2/3, add = TRUE)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/Image/Anomaly/", output,".pdf"), width = 11, height = 5)
  
  # tn_sea
  tn_sea = intersect(tn, idxdates)
  nb_sign_positive = apply(data_anomaly[,,tn_sea], MARGIN = c(1:2), FUN = function(x){sum(x>0)})
  nb_sign_negative = apply(data_anomaly[,,tn_sea], MARGIN = c(1:2), FUN = function(x){sum(x<0)})
  
  data = nb_sign_positive/length(tn_sea)
  colnames(data) = LAT; rownames(data) = LON
  data = data[,rev(seq_along(LAT))]
  output = paste0("sst_negative_1999-2017_", season[seas], "_model_", model, "_avec_sign_positive")
  image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", 
             main = output)
  world_map = map('world', plot = FALSE)
  lines(world_map, col = "white")
  contour(x = LON, y = rev(LAT), z = data, levels = 2/3, add = TRUE)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/Image/Anomaly/", output,".pdf"), width = 11, height = 5)
  
  data = nb_sign_negative/length(tn_sea)
  colnames(data) = LAT; rownames(data) = LON
  data = data[,rev(seq_along(LAT))]
  output = paste0("sst_negative_1999-2017_", season[seas], "_model_", model, "_avec_sign_negative")
  image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", 
             main = output)
  world_map = map('world', plot = FALSE)
  lines(world_map, col = "white")
  contour(x = LON, y = rev(LAT), z = data, levels = 2/3, add = TRUE)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/Image/Anomaly/", output,".pdf"), width = 11, height = 5)
}


#### dim and theta ####
var = "slp"
if (var == "slp"){
  load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_dim.RData")
  load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_theta.RData")
  plot(dim, theta)
}else if (var == "z500"){
  load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_dim.RData")
  load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_theta.RData")
  plot(dim, theta)
}else if (var == "sst"){
  load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_dim.RData")
  load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_theta.RData")
  plot(dim, theta)
}

DATE = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))
dim   = dim[DATE$Y>=1999 & DATE$Y<=2017]
theta = theta[DATE$Y>=1999 & DATE$Y<=2017]
DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]

data = as.data.frame(cbind(dim, theta))

load(file = "~/Documents/LSCE/SWG/NA_1.5/d_intensity.RData")

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

model = 7
tp = which(d[,model] - d[,64]>0); length(tp)
tn = which(d[,model] - d[,64]<0); length(tn)
data$type = NA
for (seas in 1:4){
  
  idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
  
  tp_sea = intersect(tp, idxdates)
  tn_sea = intersect(tn, idxdates)
  
  data$type[tp_sea] = "p"
  data$type[tn_sea] = "n"
  
  data_sea = data[idxdates,]
  
  library(ggplot2)
  output = paste0("dim_theta_1999-2017_by_", var, "_in_", season[seas], "_for_model_", model)
  p = ggplot(data_sea) + geom_point(aes(x = dim, y = theta, col = type)) + theme_bw()
  p = p + labs(title = output)
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/Image/SysDyn/", output, ".pdf"), width = 7, height = 7)
}
  

model = 7
thres = 0.5
tp = which(d[,model] - d[,64]> thres)
tn = which(d[,model] - d[,64]<(-thres))

data$type = NA
for (seas in 1:4){
  
  idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
  
  tp_sea = intersect(tp, idxdates)
  tn_sea = intersect(tn, idxdates)
  
  data$type[tp_sea] = "p"
  data$type[tn_sea] = "n"
  
  data_sea = data[union(tp_sea, tn_sea),]
  
  library(ggplot2)
  output = paste0("dim_theta_1999-2017_by_", var, "_in_", season[seas], "_for_model_", model, "_thres_", thres)
  p = ggplot(data_sea) + geom_point(aes(x = dim, y = theta, col = type)) + theme_bw()
  p = p + labs(title = output)
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/Image/SysDyn/", output, ".pdf"), width = 7, height = 7)
}

