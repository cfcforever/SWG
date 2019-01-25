#### load packages ####
source("function/load_packages.R")

#### common variables ####
{
  SEAS   = c('DJF','MAM','JJA','SON')
  season = c('Winter','Spring','Summer','Fall')
  MON    = c(12, 1:11)
  MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
  
  city.names = c("Paris", "Madrid", "Stockholm")
  city = city.names[1]
  
  case_name = c("statio", "slp_z500_SD")
  ncase = length(case_name)
}
####


#### load slp_z500 SysDyn ####
load(file = "~/Documents/LSCE/Dynamic System/msl_z500_19790101-20180731_NA_1.5x1.5_dim.RData")
load(file = "~/Documents/LSCE/Dynamic System/msl_z500_19790101-20180731_NA_1.5x1.5_theta.RData")
plot(dim, theta)

ndays = as.numeric(difftime(as.Date("2018-07-31"), as.Date("1979-01-01"), units = "days"))
DATE = as.Date((0:ndays), origin = "1979-01-01") ; head(DATE); tail(DATE)

MONTH = month(DATE)

for (seas in 1:4){
  dim_sea   = dim[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  theta_sea = theta[which(MONTH == MON[seas,1] | MONTH == MON[seas,2] | MONTH == MON[seas,3])]
  
  SysDyn = cbind(dim_sea, theta_sea)
  save(SysDyn, file = paste0("~/Documents/LSCE/SWG/slp_z500_SD/predictor/slp_z500_", SEAS[seas], "_SysDyn.RData"))
  
  PCS = apply(SysDyn, MARGIN = 2, FUN = function(x){scale(x)})
  print(c(colMeans(PCS), apply(PCS, 2, FUN = sd)))
  colnames(PCS) = paste0("PC", 1:ncol(PCS))
  PCS = as.data.frame(PCS)
  save(PCS, file = paste0("~/Documents/LSCE/SWG/slp_z500_SD/predictor/predictor_",SEAS[seas],"_1.RData"))
}
####


#### Estimation ####
city = city.names[3]

source("function/fun_estimation_t2m.R")

## 1979 - 1998
fun_estimation(predictor = paste0('~/Documents/LSCE/SWG/slp_z500_SD/predictor/predictor_'),
               NUM = 1,
               input.tmean = paste0('~/Documents/LSCE/SWG/data/tmean_', city, '_1979_2017.RData'),
               output.name = paste0('~/Documents/LSCE/SWG/slp_z500_SD/', city, '/SWG_ERAI_ESD_tmean_'),
               year.begin = 1979,
               year.end = 1998)

## 1999 - 2017
fun_estimation(predictor = paste0('~/Documents/LSCE/SWG/slp_z500_SD/predictor/predictor_'),
               NUM = 1,
               input.tmean = paste0('~/Documents/LSCE/SWG/data/tmean_', city, '_1979_2017.RData'),
               output.name = paste0('~/Documents/LSCE/SWG/slp_z500_SD/', city, '/SWG_ERAI_ESD_tmean_'),
               year.begin = 1999,
               year.end = 2017)
####


#### Simulation ####
source("function/fun_simulation_t2m.R")

## 1979 - 1998
for (k in 1:1){
  fun_simulation(predictor = paste0('~/Documents/LSCE/SWG/slp_z500_SD/predictor/predictor_'),
                 parameter = paste0('~/Documents/LSCE/SWG/slp_z500_SD/', city, '/SWG_ERAI_ESD_tmean_'),
                 NUM = k, 
                 type = "nat",
                 output = paste0('~/Documents/LSCE/SWG/slp_z500_SD/', city, '/'), 
                 y1 = 1999, y2 = 2017,
                 year.begin = 1979, year.end = 1998)
}

## 1999 - 2017
for (k in 1:1){
  fun_simulation(predictor = paste0('~/Documents/LSCE/SWG/slp_z500_SD/predictor/predictor_'),
                 parameter = paste0('~/Documents/LSCE/SWG/slp_z500_SD/', city, '/SWG_ERAI_ESD_tmean_'),
                 NUM = k, 
                 type = "rea",
                 output = paste0('~/Documents/LSCE/SWG/slp_z500_SD/', city, '/'), 
                 y1 = 1999, y2 = 2017,
                 year.begin = 1999, year.end = 2017)
}
####

#### create stationnary model: DON'T RUN every time !!! ####
for (city in city.names){
  load(paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))
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
  save(mean, sd, DATE, file = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/t2m_mean_sd_1999_2017_nat_0.RData"))
  
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
  save(mean, sd, DATE, file = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/t2m_mean_sd_1999_2017_rea_0.RData"))
}
####

#### SAMPLE 100 simulations ####
DATE = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))

## Gaussain distribution for temperature
Sample_gaussian_temp<- function(par_tt1,par_tt2){
  #a=date()
  #cat(a,'\n')
  Nb_stations=dim(par_tt1)[2]
  echant=array(NaN,dim=dim(par_tt1))
  for (station in 1:Nb_stations){
    for (j in 1:(dim(par_tt1)[1])){
      echant[j,station]=rnorm(1,mean=par_tt1[j,station],sd=par_tt2[j,station])
    }
  }
  #a=date()
  #cat(a,'\n')
  return(echant)
}

if (FALSE){
  for (city in city.names){
    dir.create(path = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/SIMU"))
    dir.create(path = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/SIMU/nat"))
    dir.create(path = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/SIMU/rea"))
    dir.create(path = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/Image"))
  }
}

for (city in city.names){
  for (w in c("nat", "rea")){
    for (k in 1){ # 0: stationary // 1: conditional - slp_SD
      load(paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/t2m_mean_sd_1999_2017_", w, "_", k, ".RData"))
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
        
        filesimu = paste('~/Documents/LSCE/SWG/slp_z500_SD/', city, '/SIMU/', w, '/SIMU_SWG_tmean_1999_2017_', w, '_', k, '_run_', NUM, '.RData', sep="")
        save(Sample,DATE,file = filesimu)
      }
    }
  }
}
####


#### separator in season ####
for (city in city.names){
  for (k in 0){ # 0: stationary // 1: conditional - slp_SD
    for (i in 1:100){
      NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
      load(paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_", k, "_run_", NUM, ".RData"))
      for (seas in 1:4){
        Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
        DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
        save(DATE_seas, Sample_seas, file = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/SIMU/nat/SIMU_tmean_1999_2017_nat_", k, "_", SEAS[seas], "_run_", NUM, ".RData"))
      }
    }
  }
}
####


#### FAR ####
city.names = c("Paris", "Madrid", "Stockholm")
ens_thres = 10:30
FAR_cases = c("FAR", "RR", "FARanda")
ens_year = 1999:2017

for (city in city.names){
  
  # city = city.names[1]
  load(paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))
  DATE_OBS = atoms(timeSequence(from="1979-01-01", to="2017-12-31", by='day'))
  
  list_mean_sd_nat = vector("list", 2)
  list_mean_sd_rea = vector("list", 2)
  for (k in 1:2){
    i = k-1
    load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/t2m_mean_sd_1999_2017_nat_", i,".RData"))
    cat(range(mean),'\n')
    cat(range(sd),'\n')
    m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
    colnames(m) = c("mean", "sd")
    m[,"mean"] = mean
    m[,"sd"]   = sd
    list_mean_sd_nat[[k]] = m
    
    load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/t2m_mean_sd_1999_2017_rea_", i,".RData"))
    cat(range(mean),'\n')
    cat(range(sd),'\n')
    m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
    colnames(m) = c("mean", "sd")
    m[,"mean"] = mean
    m[,"sd"]   = sd
    list_mean_sd_rea[[k]] = m
  }
  
  Obs = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
  n = length(Obs)
  DATE = as.Date((0:(n-1)), origin = "1999-01-01")
  
  # ens_thres = 10:30
  list_prob_thres = vector("list", length(ens_thres))
  names(list_prob_thres) = ens_thres
  for (thres in ens_thres){
    list_prob_thres_nat = lapply(list_mean_sd_nat, FUN = function(m){
      pnorm(thres, mean = m[,"mean"], sd = m[,"sd"], lower.tail = F)
    })
    prob_thres_nat = do.call(rbind.data.frame, list_prob_thres_nat)
    colnames(prob_thres_nat) = "nat"
    
    list_prob_thres_rea = lapply(list_mean_sd_rea, FUN = function(m){
      pnorm(thres, mean = m[,"mean"], sd = m[,"sd"], lower.tail = F)
    })
    prob_thres_rea = do.call(rbind.data.frame, list_prob_thres_rea)
    colnames(prob_thres_rea) = "rea"
    
    prob_thres = cbind(prob_thres_nat, prob_thres_rea)
    prob_thres$date = rep(DATE, ncase)
    
    prob_thres$FAR     = 1 - prob_thres$nat/prob_thres$rea
    prob_thres$RR      = prob_thres$rea/prob_thres$nat
    prob_thres$FARanda = (prob_thres$rea-prob_thres$nat)/(prob_thres$nat+prob_thres$rea)
    
    prob_thres$case = rep(case_name, each = length(DATE))
    prob_thres$thres = thres
    prob_thres$year = year(prob_thres$date)
    prob_thres$num = as.numeric(prob_thres$date)
    
    list_prob_thres[[as.character(thres)]] = prob_thres
  }
  
  prob_thres = do.call(rbind.data.frame, list_prob_thres)
  
  cc <- c("red", "blue")
  
  for (FAR_case in FAR_cases){
    for (year in ens_year){
      
      # FAR_case = FAR
      # year = 2003
      s = which(month(prob_thres$date)==7|month(prob_thres$date)==6|month(prob_thres$date)==8)
      period_sel = intersect(which(year(prob_thres$date)==year), s)
      FAR_sel = prob_thres[period_sel,]
      FAR_sel$case <- factor(FAR_sel$case, levels = case_name)
      
      breaks = rep(NA,3)
      labels = rep(NA,3)
      for (k in 1:3){
        breaks[k] = FAR_sel$num[which(month(FAR_sel$date)==unique(month(FAR_sel$date))[k])[1]]
        labels[k] = as.character(FAR_sel$date[which(month(FAR_sel$date)==unique(month(FAR_sel$date))[k])[1]])
      }
      
      saveGIF(
        for (thres in ens_thres){
          p <- ggplot(FAR_sel, aes(x=num, y=FARanda, group=case, col=case)) + 
            geom_blank() +
            geom_line(data = FAR_sel[FAR_sel$thres==thres,], mapping = aes(x=num, y=FARanda, group=case, col=case)) + theme_bw() +
            scale_color_manual(values = cc, name = "") + 
            scale_x_continuous(breaks = breaks, labels = labels) +
            labs(x = "date") +
            ggtitle(paste0(FAR_case, " for summer in ", year, " over ", thres, " degree"))
          print(p)
          # ggsave(plot = p, filename = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", i, ".png"), units = "in", dpi = 300, width = 7, height = 7)
          # dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", year, ".pdf"), width = 7, height = 7)
          # dev.print(png, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", i, ".png"), width = 700, height = 700)
        }, 
        ani.width = 200*7, ani.height = 200*7, ani.res = 200, interval = 0.5, autobrowse = F,
        movie.name = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/FAR/", FAR_case, "_summer_in_", year, ".gif")
      )
    }
  }
}
####

#### Intensity: SAVE d_intensity.RData ####
for (city in city.names){
  load(paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))
  DATE_OBS = atoms(timeSequence(from="1979-01-01", to="2017-12-31", by='day'))
  
  Obs = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
  n = length(Obs)
  DATE = as.Date((0:(n-1)), origin = "1999-01-01")
  
  list_mean_sd_nat = vector("list", 2)
  list_mean_sd_rea = vector("list", 2)
  for (k in 1:2){
    i = k-1
    load(paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/t2m_mean_sd_1999_2017_nat_", i,".RData"))
    cat(range(mean),'\n')
    cat(range(sd),'\n')
    m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
    colnames(m) = c("mean", "sd")
    m[,"mean"] = mean
    m[,"sd"]   = sd
    list_mean_sd_nat[[k]] = m
    
    load(paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/t2m_mean_sd_1999_2017_rea_", i,".RData"))
    cat(range(mean),'\n')
    cat(range(sd),'\n')
    m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
    colnames(m) = c("mean", "sd")
    m[,"mean"] = mean
    m[,"sd"]   = sd
    list_mean_sd_rea[[k]] = m
  }
  
  p1 = as.data.frame(matrix(NA, nrow = n, ncol = ncase))
  for (l in 1:ncase){
    for (k in 1:n){
      p1[k,l] = pnorm(Obs[k], mean = list_mean_sd_rea[[l]][k,"mean"], sd = list_mean_sd_rea[[l]][k,"sd"], lower.tail = F)
    }
  }
  
  p0 = p1
  
  i1 = as.data.frame(matrix(NA, nrow = n, ncol = ncase))
  for (l in 1:ncase){
    for (k in 1:n){
      i1[k,l] = qnorm(1-p0[k,l], mean = list_mean_sd_nat[[l]][k,"mean"], sd = list_mean_sd_nat[[l]][k,"sd"], lower.tail = T)
    }
  }
  ## d_intensity: 1999-2017
  d = Obs - i1; range(d)
  colnames(d) = c("statio", "slp_SD")
  DATE = as.Date((0:(n-1)), origin = "1999-01-01")
  d$date = DATE
  save(d, file = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/d_intensity.RData"))
}
####


#### Intensity: plot for each year from 1999 to 2017 ####
for (city in city.names){
  load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/d_intensity.RData"))
  
  melt_d = melt(d, id.vars = "date")
  colnames(melt_d) = c("date", "case", "d")
  
  for (year in 1999:2017){
    period_sel = intersect(which(year(d$date)==year) , which(month(d$date)==6|month(d$date)==7|month(d$date)==8))
    
    d1 = d[period_sel,]
    d1 = melt(d1, id.vars = "date")
    colnames(d1) = c("date", "case", "d")
    d1$num = rep(1:length(period_sel), ncase)
    
    breaks = rep(NA,3)
    labels = rep(NA,3)
    for (k in 1:3){
      breaks[k] = which(month(d1$date)==unique(month(d1$date))[k])[1]
      labels[k] = as.character(d1$date[breaks[k]])
    }
    
    p = ggplot(d1)
    p = p + geom_line(aes(x=num, y=d, group=case, col=case)) + theme_bw()
    cc <- c("red", "blue")
    p = p + scale_color_manual(values = cc, labels = paste0(c("statio", "slp_SD")))
    p = p + scale_x_continuous(breaks = breaks, labels = labels) 
    p = p + labs(title = paste0("delta_Intensity_tmean_summer_in_", year),
                 x = "date", y = expression(Delta*"I"))
    print(p)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/d_intensity_tmean_summer_", year, ".pdf"), width = 7, height=7)
    print(year)
  }
}
####


#### Anomaly slp version 0 ####
var = "slp"
cas = 3

if (cas == 1){
  load(file = paste0("/Volumes/Data-ExFAT/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1979_2017.RData"))
  cas = "1_1979_2017"
}else if(cas == 2){
  load(file = paste0("/Volumes/Data-ExFAT/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1979_1998.RData"))
  cas = "2_1979_1998"
}else if(cas == 3){
  load(file = paste0("/Volumes/Data-ExFAT/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1999_2017.RData"))
  cas = "3_1999_2017"
}

data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]

list.anomaly.maps = vector("list", 3)
for (city in city.names){
  load(file = paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/d_intensity.RData"))
  dI = d
  
  {
    t_sup = which(dI[,"slp_SD"] - dI[,"statio"]>0); length(t_sup)
    t_inf = which(dI[,"slp_SD"] - dI[,"statio"]<0); length(t_inf)
    
    # list.data = vector("list", 2)
    for (seas in c(1,3)){
      idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
      
      t_sup_sea = intersect(t_sup, idxdates)
      data = apply(data_anomaly[,,t_sup_sea], 1:2, mean)
      colnames(data) = LAT; rownames(data) = LON
      data = data[,rev(seq_along(LAT))]
      print(range(data))
      
      col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
      cc <- colorRamps::matlab.like2(n = 100)
      
      output = paste0(var, "_sup_1999-2017_", season[seas], "_slp_SD_", cas)
      {
        dat = melt(data, varnames = c("long", "lat"))
        
        world <- map_data("world")
        col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
        brks = seq(from = -180, to = 180, by = 30)
        label = seq(from = -180, to = 180, by = 30)
        col.contour = "white"
        worldmap <- ggplot() + theme_bw() +
          geom_raster(data = dat, mapping = aes(x = long, y = lat, fill = value), interpolate = TRUE) +
          scale_fill_gradientn(colours = matlab.like2(100), limits = col.range) +
          geom_contour(data = dat, mapping = aes(x = long, y = lat, z = value), 
                       col = col.contour, breaks = brks) +
          geom_dl(data = dat, aes(x = long, y = lat, z = value, label=..level..), 
                  method = "bottom.pieces", stat = "contour", 
                  col = col.contour, breaks = brks) +
          geom_path(world, mapping = aes(x = long, y = lat, group = group)) +
          coord_cartesian(xlim = range(LON)+c(-0.5,0.5), ylim = range(LAT)+c(-0.5,0.5)) + 
          scale_y_continuous(breaks = (3:7) * 10, expand = c(0,0)) +
          scale_x_continuous(breaks = (-4:2) * 20, expand = c(0,0)) + 
          labs(x = "Lon", y = "Lat", fill = "Pa")
        print(worldmap)
      }
      dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/Image/", output,".pdf"), width = 11, height = 5)
      
      
      t_inf_sea = intersect(t_inf, idxdates)
      data = apply(data_anomaly[,,t_inf_sea], 1:2, mean)
      colnames(data) = LAT; rownames(data) = LON
      data = data[,rev(seq_along(LAT))]
      print(range(data))
      
      col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
      cc <- colorRamps::matlab.like2(n = 100)
      
      output = paste0(var, "_inf_1999-2017_", season[seas], "_slp_SD_", cas)
      {
        dat = melt(data, varnames = c("long", "lat"))
        
        world <- map_data("world")
        col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
        brks = seq(from = -180, to = 180, by = 30)
        label = seq(from = -180, to = 180, by = 30)
        col.contour = "white"
        worldmap <- ggplot() + theme_bw() +
          geom_raster(data = dat, mapping = aes(x = long, y = lat, fill = value), interpolate = TRUE) +
          scale_fill_gradientn(colours = matlab.like2(100), limits = col.range) +
          geom_contour(data = dat, mapping = aes(x = long, y = lat, z = value), 
                       col = col.contour, breaks = brks) +
          geom_dl(data = dat, aes(x = long, y = lat, z = value, label=..level..), 
                  method = "bottom.pieces", stat = "contour", 
                  col = col.contour, breaks = brks) +
          geom_path(world, mapping = aes(x = long, y = lat, group = group)) +
          coord_cartesian(xlim = range(LON)+c(-0.5,0.5), ylim = range(LAT)+c(-0.5,0.5)) + 
          scale_y_continuous(breaks = (3:7) * 10, expand = c(0,0)) +
          scale_x_continuous(breaks = (-4:2) * 20, expand = c(0,0)) + 
          labs(x = "Lon", y = "Lat", fill = "Pa")
        print(worldmap)
      }
      dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_z500_SD/", city, "/Image/", output,".pdf"), width = 11, height = 5)
    }
  }
}
####


#### Anomaly slp version I ####
var = "slp"
cas = 3

if (cas == 1){
  load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1979_2017.RData"))
  cas = "1979_2017"
}else if(cas == 2){
  load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1979_1998.RData"))
  cas = "1979_1998"
}else if(cas == 3){
  load(file = paste0("~/Documents/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1999_2017.RData"))
  cas = "1999_2017"
}


for (city in city.names){
  load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/d_intensity.RData"))
  dI = d
  
  {
    t_sup = which(dI[,"slp_SD"] - dI[,"statio"]>0); length(t_sup)
    t_inf = which(dI[,"slp_SD"] - dI[,"statio"]<0); length(t_inf)
    
    data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
    DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]
    
    for (seas in c(1,3)){
      idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
      t_sup_sea = intersect(t_sup, idxdates)
      data = apply(data_anomaly[,,t_sup_sea], 1:2, mean)
      colnames(data) = LAT; rownames(data) = LON
      data = data[,rev(seq_along(LAT))]
      print(range(data))
      
      col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
      cc <- colorRamps::matlab.like2(n = 100)
      
      output = paste0(var, "_sup_1999-2017_", season[seas], "_slp_SD_2")
      {
        # image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
        #            main = output)
        image.plot(data, x = LON, y = rev(LAT), xlab = "", ylab = "", col = cc, zlim = col.range)
        world_map = map('world', plot = FALSE)
        lines(world_map, col = "black")
        contour(x = LON, y = rev(LAT), z = data, add = TRUE, col = "white") 
      }
      dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output,".pdf"), width = 11, height = 5)
      
      
      t_inf_sea = intersect(t_inf, idxdates)
      data = apply(data_anomaly[,,t_inf_sea], 1:2, mean)
      colnames(data) = LAT; rownames(data) = LON
      data = data[,rev(seq_along(LAT))]
      print(range(data))
      
      col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
      cc <- colorRamps::matlab.like2(n = 100)
      
      output = paste0(var, "_inf_1999-2017_", season[seas], "_slp_SD_2")
      {
        # image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
        #            main = output)
        image.plot(data, x = LON, y = rev(LAT), xlab = "", ylab = "", col = cc, zlim = col.range)
        world_map = map('world', plot = FALSE)
        lines(world_map, col = "black")
        contour(x = LON, y = rev(LAT), z = data, add = TRUE, col = "white")
      }
      dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output,".pdf"), width = 11, height = 5)
    }
  }
}
####


#### dim and theta depend on temperature ####
var = "slp"
if (var == "slp"){
  load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_dim.RData")
  load("~/Documents/LSCE/Dynamic System/msl_1979-2018_NA_1.5x1.5_theta.RData")
  plot(dim, theta, main = paste0("theta vs dim of ", var))
}else if (var == "z500"){
  load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_dim.RData")
  load("~/Documents/LSCE/Dynamic System/z500_1979-2018_NA_1.5x1.5_theta.RData")
  plot(dim, theta, main = paste0("theta vs dim of ", var))
}else if (var == "sst"){
  load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_small_dim.RData")
  load("~/Documents/LSCE/Dynamic System/sst_1979-2018_NA_1.5x1.5_small_theta.RData")
  plot(dim, theta, main = paste0("theta vs dim of ", var))
}

DATE  = atoms(timeSequence(from="1979-01-01",to="2018-07-31",by='day'))
dim   = dim[DATE$Y>=1979 & DATE$Y<=2017]
theta = theta[DATE$Y>=1979 & DATE$Y<=2017]
DATE  = DATE[DATE$Y>=1979 & DATE$Y<=2017,]
D = as.Date(0:(nrow(DATE)-1), origin = "1979-01-01")
plot(dim, theta, main = paste0("theta vs dim of ", var))

data = as.data.frame(cbind(dim, theta))

## choose city
city = city.names[3]
load(paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))
dat = cbind(data, tmean, D)

## dim vs temperature
p1 <- ggplot(dat) + geom_blank() + theme_bw() +
  geom_point(aes(x=tmean, y=dim)) + 
  xlab("Temperature") + 
  ggtitle(paste0("dim and theta vs temperature from 1979 to 2017 for ", city))
p2 <- ggplot(dat) + geom_blank() + theme_bw() +
  geom_point(aes(x=tmean, y=theta)) +
  xlab("Temperature")
plot_grid(p1, p2, ncol = 1)
save_pdf(filename = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/Image/SysDyn_vs_temperature_", city), 
         width = 7, height = 7)


## correlation between dim and T (theta and T)
cor_dim_T = cor(data$dim,   data$tmean, method = "spearman")
cor_the_T = cor(data$theta, data$tmean, method = "spearman")
for (seas in 1:4){
  idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
  data_sea = data[idxdates,]
  cor_dim_T_sea = cor(data_sea$dim, data_sea$tmean, method = "spearman")
  cor_the_T_sea = cor(data_sea$theta, data_sea$tmean, method = "spearman")
  
  cor_SD_T_sea_year = as.data.frame(matrix(NA, nrow = length(1979:2017), ncol = 3))
  colnames(cor_SD_T_sea_year) = c("dim_T", "the_T", "year")
  rownames(cor_SD_T_sea_year) = 1979:2017
  for (year in 1979:2017){
    idxdates = intersect(which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3]),
                         which(DATE['Y']==year))
    data_sea_year = data[idxdates,]
    
    cor_SD_T_sea_year[rownames(cor_SD_T_sea_year)==year,] = c(cor(data_sea_year$dim, data_sea_year$tmean, method = "spearman"),
                                                              cor(data_sea_year$theta, data_sea_year$tmean, method = "spearman"),
                                                              year)
  }
  
  p <- ggplot() + theme_bw() +
    geom_point(cor_SD_T_sea_year, mapping = aes(x=year,y=dim_T, col = "black")) +
    geom_line(mapping = aes(x=1979:2017,y=cor_dim_T_sea, col = "blue")) + 
    geom_line(mapping = aes(x=1979:2017,y=cor_dim_T, col = "red")) +
    scale_colour_manual(name = 'colour', guide = 'legend',
                        values =c('black'='black', 'blue'='blue', 'red'='red'), 
                        labels = c('by year', paste0("all ", season[seas]), 'all data')) +
    ylab("correlation dim vs T") + ggtitle(paste0(season[seas], " ", city))
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/dim_theta/cor_dim_T_", season[seas], ".pdf"), width=7, height=7)
  
  p <- ggplot() + theme_bw() +
    geom_point(cor_SD_T_sea_year, mapping = aes(x=year,y=the_T, col = "black")) +
    geom_line(mapping = aes(x=1979:2017,y=cor_the_T_sea, col = "blue")) + 
    geom_line(mapping = aes(x=1979:2017,y=cor_the_T, col = "red")) +
    scale_colour_manual(name = 'colour', guide = 'legend',
                        values =c('black'='black', 'blue'='blue', 'red'='red'), 
                        labels = c('by year', paste0("all ", season[seas]), 'all data')) +
    ylab("correlation the vs T") + ggtitle(paste0(season[seas], " ", city))
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/dim_theta/cor_the_T_", season[seas], ".pdf"), width=7, height=7)
}

## dim theta vs T (when T max)
for (year in 1979:2017){
  for (seas in 3){
    idxdates = intersect(which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3]),
                         which(DATE['Y']==year))
    data_sea = data[idxdates,]
    
    pt = which(data_sea$tmean%in%sort(data_sea$tmean, decreasing = T)[1:5])
    # pt = which(data_sea$tmean%in%sort(data_sea$tmean, decreasing = T)[1:10])
    
    p <- ggplot(data_sea) + theme_bw() +
      geom_line(aes(x = D, y = dim), col = "black") + 
      geom_point(data_sea[pt,], mapping = aes(x = D, y = dim), col = "red") +
      # geom_point(data_sea[pt,], mapping = aes(x = D, y = tmean), col = "red") +
      xlab("date") + ggtitle(paste0(season[seas], " ", year, " ", city))
    print(p)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/dim_theta/dim_vs_T_", season[seas], "_", year, ".pdf"), width = 11, height = 5)
    
    p <- ggplot(data_sea) + theme_bw() +
      geom_line(aes(x = D, y = theta), col = "black") + 
      geom_point(data_sea[pt,], mapping = aes(x = D, y = theta), col = "red") +
      # geom_point(data_sea[pt,], mapping = aes(x = D, y = tmean), col = "red") +
      xlab("date") + ggtitle(paste0(season[seas], " ", year, " ", city))
    print(p)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/dim_theta/theta_vs_T_", season[seas], "_", year, ".pdf"), width = 11, height = 5)
  }
}

## all seasons 
output = paste0("dim_theta_with_tmean_of_", city, "_1999-2017")
p = ggplot(data) + geom_point(aes(x = dim, y = theta, col = tmean)) + theme_bw()
p = p + scale_color_gradientn(colours = matlab.like2(10), limits = fun_limits(data[,1:3]))
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output, ".pdf"), width = 7, height = 7)

## all seasons with threshold
for (thres in 20:30){
  # thres = 27
  if (T){
    output = paste0("dim_theta_with_tmean_of_", city, "_1999-2017_over_", thres)
  }else{
    output = paste0("dim_theta_with_tmean_of_", city, "_1999-2017_under_", thres)
  }
  print(output)
  
  p <- ggplot(data) + geom_point(aes(x = dim, y = theta), col = "black", size = 0.5) +
    geom_point(data = data[data$tmean>=thres,], mapping = aes(x = dim, y = theta, col = tmean)) + 
    theme_bw() + scale_color_gradientn(colours = matlab.like2(10), limits = fun_limits(data))
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output, ".pdf"), width = 7, height = 7) 
}

# each season
for (seas in 1:4){
  idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
  data_sea = data[idxdates,]
  
  output = paste0("dim_theta_with_tmean_of_", city, "_1999-2017_", season[seas])
  p = ggplot(data, mapping = aes(x = dim, y = theta)) + geom_blank()
  p = p + geom_point(data_sea, mapping = aes(x = dim, y = theta, col = tmean)) + theme_bw()
  p = p + scale_color_gradientn(colours = matlab.like2(10), limits = fun_limits(data_sea))
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output, ".pdf"), width = 7, height = 7)
}

#### max annual temperature ####
for (city in city.names){
  
  ## annual_max_in_plot_theta_dim
  load(file = paste0("~/Documents/LSCE/SWG/data/tmean_", city, "_1979_2017.RData"))
  DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  D = as.Date(0:(nrow(DATE_OBS)-1), origin = "1979-01-01")
  
  tam = as.data.frame(matrix(NA, nrow = length(unique(DATE_OBS$Y)), ncol = 2))
  colnames(tam) = c("tam", "date")
  rownames(tam) = unique(DATE_OBS$Y)
  for (y in 1979:2017){
    Dy = D[DATE_OBS$Y==y][which(tmean[DATE_OBS$Y==y,] == max(tmean[DATE_OBS$Y==y,]))]
    tam[paste(y),] = c(max(tmean[DATE_OBS$Y==y,]), as.character(Dy))
  }
  
  p <- ggplot(data) + theme_bw() + 
    geom_point(aes(x = dim, y = theta)) + 
    geom_point(data = data[as.character(D)%in%tam$date, ], aes(x = dim, y = theta),  col = "red")
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/annual_max_in_plot_theta_dim.pdf"), width = 7, height = 7)
  
  ## tmean_max_annual_city_1999_2017_boxplot_rea
  tam.list = vector("list", 100)
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/SIMU/rea/SIMU_SWG_tmean_1999_2017_rea_1_run_", NUM, ".RData"))
    tam_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
    colnames(tam_sim) = c("tam", "year")
    rownames(tam_sim) = unique(DATE$Y)
    for (y in unique(DATE$Y)){
      tam_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
    }
    tam.list[[i]] = tam_sim
  }
  tam_sim = do.call(rbind.data.frame, tam.list)
  tam_sim$type = "sim"
  
  tam_obs = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
  colnames(tam_obs) = c("tam", "year")
  rownames(tam_obs) = unique(DATE$Y)
  tam_obs[,1] = as.numeric(tam[as.character(unique(DATE$Y)),1])
  tam_obs[,2] = unique(DATE$Y)
  tam_obs$type = "obs"
  
  tam_data = rbind(tam_obs, tam_sim)
  
  output = paste0("tmean_max_annual_", city, "_1999_2017_boxplot_rea")
  p <- ggplot(tam_data) + theme_bw() +
    geom_boxplot(data = tam_data[tam_data$type=="sim",], mapping = aes(x = year, y = tam, group = year)) +
    geom_point(data = tam_data[tam_data$type=="obs",], mapping = aes(x = year, y = tam), col = "red") + 
    ggtitle(output)
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output, ".pdf"), width = 7, height = 7)
  
  ## tmean_max_annual_city_1999_2017_boxplot_nat
  tam.list = vector("list", 100)
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/SIMU/nat/SIMU_SWG_tmean_1999_2017_nat_1_run_", NUM, ".RData"))
    tam_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
    colnames(tam_sim) = c("tam", "year")
    rownames(tam_sim) = unique(DATE$Y)
    for (y in unique(DATE$Y)){
      tam_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
    }
    tam.list[[i]] = tam_sim
  }
  tam_sim = do.call(rbind.data.frame, tam.list)
  tam_sim$type = "sim"
  
  tam_obs = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
  colnames(tam_obs) = c("tam", "year")
  rownames(tam_obs) = unique(DATE$Y)
  tam_obs[,1] = as.numeric(tam[as.character(unique(DATE$Y)),1])
  tam_obs[,2] = unique(DATE$Y)
  tam_obs$type = "obs"
  
  tam_data = rbind(tam_obs, tam_sim)
  
  output = paste0("tmean_max_annual_", city, "_1999_2017_boxplot_nat")
  p <- ggplot(tam_data) + theme_bw() +
    geom_boxplot(data = tam_data[tam_data$type=="sim",], mapping = aes(x = year, y = tam, group = year)) +
    geom_point(data = tam_data[tam_data$type=="obs",], mapping = aes(x = year, y = tam), col = "red") + 
    ggtitle(output)
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output, ".pdf"), width = 7, height = 7)
  
}


# for (k in 1:3){
#   city = city.names[k]
#   system2('open', args = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/annual_max_in_plot_theta_dim.pdf"), wait = F)
# }
