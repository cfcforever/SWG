# load package ------------------------------------------------------------
source("function/load_packages.R")


# Preliminary -------------------------------------------------------------

# nc <- nc_open("/Volumes/Data-ExFAT/DATA/EOBS_v17/tg_0.50deg_reg_v17.0.nc")
nc <- nc_open("/Volumes/Data-ExFAT/DATA/EOBS_v17/tg_1.50deg_reg_cdo_remapbli.nc")
print(nc)
data = ncvar_get(nc, "tg") 
lon = ncvar_get(nc, "longitude")
lat = ncvar_get(nc, "latitude")
time = ncvar_get(nc, "time")
nc_close(nc)
range(data)

# DATE from 1979-01-01 to 2017-12-31
date = as.Date(0:(length(time)-1), origin = "1950-01-01")
head(date); tail(date)
date_chosen = which(date=="1979-01-01"):which(date=="2017-12-31")
head(date[date_chosen]); tail(date[date_chosen])
DATE = date[date_chosen]

# EUROPE domain: [35N, 60N]X[10W,25E]
lon_chosen = which(lon>=-10 & lon<=25)
LON = lon[lon_chosen]
lat_chosen = which(lat>=30 & lat<=60)
LAT = lat[lat_chosen]
dat = data[lon_chosen, lat_chosen, date_chosen]

save(dat, LON, LAT, DATE, file = "/Volumes/Data-ExFAT/DATA/EOBS_v17/tg_EOBS-1.5_1979-2017_EU.RData")
####

# common variables --------------------------------------------------------
{
  SEAS   = c('DJF','MAM','JJA','SON')
  season = c('Winter','Spring','Summer','Fall')
  MON    = c(12, 1:11)
  MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
  
  case_name = c("statio", "slp_SD")
  ncase = length(case_name)
}
## Function
source("function/fun_limits.R")
####

# Start -------------------------------------------------------------------
load(file = "/Volumes/Data-ExFAT/DATA/EOBS_v17/tg_EOBS-1.5_1979-2017_EU.RData")
dim(dat)   # 23 x 21 x 14245

data = melt(dat)
colnames(data) = c("lon", "lat", "time", "tg")

# Estimation --------------------------------------------------------------
source("function/fun_estimation_t2m.R")

DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
DATE_nat = atoms(timeSequence(from="1979-01-01",to="1998-12-31",by='day'))
DATE_rea = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))

idxnat1 = which(DATE_nat['m']==MON[1,1] | DATE_nat['m']==MON[1,2] | DATE_nat['m']==MON[1,3])
idxnat2 = which(DATE_nat['m']==MON[2,1] | DATE_nat['m']==MON[2,2] | DATE_nat['m']==MON[2,3])
idxnat3 = which(DATE_nat['m']==MON[3,1] | DATE_nat['m']==MON[3,2] | DATE_nat['m']==MON[3,3])
idxnat4 = which(DATE_nat['m']==MON[4,1] | DATE_nat['m']==MON[4,2] | DATE_nat['m']==MON[4,3])

idxrea1 = which(DATE_rea['m']==MON[1,1] | DATE_rea['m']==MON[1,2] | DATE_rea['m']==MON[1,3])
idxrea2 = which(DATE_rea['m']==MON[2,1] | DATE_rea['m']==MON[2,2] | DATE_rea['m']==MON[2,3])
idxrea3 = which(DATE_rea['m']==MON[3,1] | DATE_rea['m']==MON[3,2] | DATE_rea['m']==MON[3,3])
idxrea4 = which(DATE_rea['m']==MON[4,1] | DATE_rea['m']==MON[4,2] | DATE_rea['m']==MON[4,3])

list_tg = vector("list", length(LON)*length(LAT))
names(list_tg) = paste0(rep(1:length(LON), each = length(LAT)), "-", rep(1:length(LAT), length(LON)))
for (lon in 1:length(LON)){
  for (lat in 1:length(LAT)){
    tg = data[data$lon==lon & data$lat==lat, "tg"]
    # if (sum(!is.na(tg[idxnat1]))==0 | sum(!is.na(tg[idxnat2]))==0 | sum(!is.na(tg[idxnat3]))==0 | sum(!is.na(tg[idxnat4]))==0 |
    #     sum(!is.na(tg[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017][idxrea1]))==0 | 
    #     sum(!is.na(tg[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017][idxrea1]))==0 |
    #     sum(!is.na(tg[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017][idxrea1]))==0 |
    #     sum(!is.na(tg[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017][idxrea1]))==0){
    #   list_tg[paste0(lon,"-",lat)] = NA
    # }
    if (sum(is.na(tg))!=0){
      list_tg[paste0(lon,"-",lat)] = NA
    }
    else{
      print(sum(!is.na(tg)))
      tg = as.data.frame(matrix(tg,ncol = 1))
      list_tg[[paste0(lon,"-",lat)]] = tg

      ## 1979 - 1998
      fun_estimation(predictor = paste0('~/Documents/LSCE/SWG/spatial_maps/predictor/predictor_'),
                     NUM = 1,
                     input.tmean = list_tg[[paste0(lon,"-",lat)]],
                     output.name = paste0('~/Documents/LSCE/SWG/spatial_maps/estimation/SWG_ERAI_ESD_', lon, '-', lat, '_tmean_'),
                     year.begin = 1979,
                     year.end = 1998)

      ## 1999 - 2017
      fun_estimation(predictor = paste0('~/Documents/LSCE/SWG/spatial_maps/predictor/predictor_'),
                     NUM = 1,
                     input.tmean = list_tg[[paste0(lon,"-",lat)]],
                     output.name = paste0('~/Documents/LSCE/SWG/spatial_maps/estimation/SWG_ERAI_ESD_', lon, '-', lat, '_tmean_'),
                     year.begin = 1999,
                     year.end = 2017)
    }
  }
}
save(list_tg, file = "~/Documents/LSCE/SWG/spatial_maps/list_tg.RData")

lon = 1
lat = 20
tmean = list_tg[[paste0(lon,"-",lat)]]; sum(!is.na(tmean))


## 1979 - 1998
fun_estimation(predictor = paste0('~/Documents/LSCE/SWG/spatial_maps/predictor/predictor_'),
               NUM = 1,
               input.tmean = list_tg[[paste0(lon,"-",lat)]],
               output.name = paste0('~/Documents/LSCE/SWG/spatial_maps/SWG_ERAI_ESD_tmean_'),
               year.begin = 1979,
               year.end = 1998)

## 1999 - 2017
fun_estimation(predictor = paste0('~/Documents/LSCE/SWG/spatial_maps/predictor/predictor_'),
               NUM = 1,
               input.tmean = list_tg[[paste0(lon,"-",lat)]],
               output.name = paste0('~/Documents/LSCE/SWG/spatial_maps/SWG_ERAI_ESD_tmean_'),
               year.begin = 1999,
               year.end = 2017)
####


# Simulation --------------------------------------------------------------
source("function/fun_simulation_t2m.R")

load(file = "~/Documents/LSCE/SWG/spatial_maps/list_tg.RData")
for (lon in 1:length(LON)){
  for (lat in 1:length(LAT)){
    tg = list_tg[[paste0(lon,"-",lat)]]
    if (sum(is.na(tg))==0){
      ## 1979 - 1998
      fun_simulation(predictor = paste0('~/Documents/LSCE/SWG/spatial_maps/predictor/predictor_'),
                     parameter = paste0('~/Documents/LSCE/SWG/spatial_maps/estimation/SWG_ERAI_ESD_', lon, '-', lat, '_tmean_'),
                     NUM = 1, 
                     type = "nat",
                     output = paste0('~/Documents/LSCE/SWG/spatial_maps/simulation/sim_', lon, '-', lat, '_'), 
                     y1 = 1999, y2 = 2017,
                     year.begin = 1979, year.end = 1998)
      
      ## 1999 - 2017
      fun_simulation(predictor = paste0('~/Documents/LSCE/SWG/spatial_maps/predictor/predictor_'),
                     parameter = paste0('~/Documents/LSCE/SWG/spatial_maps/estimation/SWG_ERAI_ESD_', lon, '-', lat, '_tmean_'),
                     NUM = 1, 
                     type = "rea",
                     output = paste0('~/Documents/LSCE/SWG/spatial_maps/simulation/sim_', lon, '-', lat, '_'), 
                     y1 = 1999, y2 = 2017,
                     year.begin = 1999, year.end = 2017)
    }
  }
}
####


# Intensity: SAVE d_intensity.RData ---------------------------------------
n = sum(DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017)
dI = array(NA, dim = c(length(LON), length(LAT), n))

DATE_OBS = atoms(timeSequence(from="1979-01-01", to="2017-12-31", by='day'))

for (lon in 1:length(LON)){
  for (lat in 1:length(LAT)){
    tg = list_tg[[paste0(lon,"-",lat)]]
    
    if (sum(is.na(tg))==0){
      Obs = tg[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
      DATE = as.Date((0:(n-1)), origin = "1999-01-01")
      
      load(paste0("~/Documents/LSCE/SWG/spatial_maps/simulation/sim_", lon, "-", lat, "_t2m_mean_sd_1999_2017_nat_1.RData"))
      cat(range(mean),'\n')
      cat(range(sd),'\n')
      mean_sd_nat = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
      colnames(mean_sd_nat) = c("mean", "sd")
      mean_sd_nat[,"mean"] = mean
      mean_sd_nat[,"sd"]   = sd
      
      load(paste0("~/Documents/LSCE/SWG/spatial_maps/simulation/sim_", lon, "-", lat, "_t2m_mean_sd_1999_2017_rea_1.RData"))
      cat(range(mean),'\n')
      cat(range(sd),'\n')
      mean_sd_rea = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
      colnames(mean_sd_rea) = c("mean", "sd")
      mean_sd_rea[,"mean"] = mean
      mean_sd_rea[,"sd"]   = sd
      
      p1 = as.data.frame(matrix(NA, nrow = n, ncol = 1))
      for (k in 1:n){
        p1[k,1] = pnorm(Obs[k], mean = mean_sd_rea[k,"mean"], sd = mean_sd_rea[k,"sd"], lower.tail = F)
      }
      
      p0 = p1
      
      i1 = as.data.frame(matrix(NA, nrow = n, ncol = 1))
      for (k in 1:n){
        i1[k,1] = qnorm(1-p0[k,1], mean = mean_sd_nat[k,"mean"], sd = mean_sd_nat[k,"sd"], lower.tail = T)
      }
      
      d = Obs - i1
      
      dI[lon,lat,] = d[,1]
    }
  }
}

save(dI, LON, LAT, DATE, file = paste0("~/Documents/LSCE/SWG/spatial_maps/dI.RData"))
####

# plot maps ---------------------------------------------------------------
library(ggplot2)

load("~/Documents/LSCE/SWG/spatial_maps/dI.RData")

## just for one specific spatial map
data = dI[,,1]
colnames(data) = LAT
rownames(data) = LON

dat = melt(data, varnames = c("long", "lat"))

world <- map_data("world")
col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
worldmap <- ggplot() + theme_bw() +
  geom_raster(data = dat, mapping = aes(x = long, y = lat, fill = value), interpolate = TRUE) +
  scale_fill_gradientn(colours = matlab.like2(10), limits = col.range) +
  geom_path(world, mapping = aes(x = long, y = lat, group = group)) +
  coord_cartesian(xlim = range(LON)+c(-0.5,0.5), ylim = range(LAT)+c(-0.5,0.5)) + 
  scale_y_continuous(breaks = (3:7) * 10, expand = c(0,0)) +
  scale_x_continuous(breaks = (-4:2) * 20, expand = c(0,0)) + 
  labs(x = "Lon", y = "Lat", fill = "°C")
print(worldmap)


## for multiple spatial maps
dd = melt(dI)
colnames(dd) = c("long", "lat", "time", "value")
dd$long = rep(LON, length(LAT)*length(DATE))
dd$lat  = rep(rep(LAT, each = length(LON)), length(DATE))
dd$time = rep(DATE, each = length(LON)*length(LAT))

t0 = which(DATE=="2006-01-01")
t1 = which(DATE=="2006-12-31")

dat = dd[dd$time%in%DATE[t0:t1],]
col.range = c(-max(abs(range(dat$value[!is.na(dat$value)]))), max(abs(range(dat$value[!is.na(dat$value)]))))

world <- map_data("world")
for (i in t0:t1){
  worldmap <- ggplot() + theme_bw() +
    geom_raster(data = dd[dd$time==DATE[i],], mapping = aes(x = long, y = lat, fill = value), interpolate = TRUE) +
    scale_fill_gradientn(colours = matlab.like2(100), limits = col.range, breaks = seq(-9,9,3)) +
    geom_path(world, mapping = aes(x = long, y = lat, group = group)) +
    coord_cartesian(xlim = range(LON)+c(-0.5,0.5), ylim = range(LAT)+c(-0.5,0.5)) + 
    scale_y_continuous(breaks = (3:7) * 10, expand = c(0,0)) +
    scale_x_continuous(breaks = (-4:2) * 20, expand = c(0,0)) + 
    labs(x = "Lon", y = "Lat", fill = "°C") +
    ggtitle(paste0("changes of intensities in Europe at ", DATE[i]))
  ggsave(plot = worldmap, filename = paste0("~/Documents/LSCE/SWG/spatial_maps/Image/", i, ".png"), units = "in", dpi = 300, width = 7, height = 7)
}

## make animation
library(animation)
saveGIF(
  for (i in 1:92){
    worldmap <- ggplot() + theme_bw() +
      geom_raster(data = dd[dd$time==DATE[i],], mapping = aes(x = long, y = lat, fill = value), interpolate = TRUE) +
      scale_fill_gradientn(colours = matlab.like2(100), limits = col.range, breaks = seq(-9,9,3)) +
      geom_path(world, mapping = aes(x = long, y = lat, group = group)) +
      coord_cartesian(xlim = range(LON)+c(-0.5,0.5), ylim = range(LAT)+c(-0.5,0.5)) + 
      scale_y_continuous(breaks = (3:7) * 10, expand = c(0,0)) +
      scale_x_continuous(breaks = (-4:2) * 20, expand = c(0,0)) + 
      labs(x = "Lon", y = "Lat", fill = "°C") +
      ggtitle(paste0("changes of intensities in Europe at ", DATE[i]))
    print(worldmap)
}, ani.width = 200*7, ani.height = 200*7, ani.res = 200, interval = 0.4, movie.name = "~/Documents/LSCE/SWG/spatial_maps/Image/animation.gif")



