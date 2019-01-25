#### load packages ####
source("function/load_packages.R")
####

#### common variables ####
{
  SEAS   = c('DJF','MAM','JJA','SON')
  season = c('Winter','Spring','Summer','Fall')
  MON    = c(12, 1:11)
  MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
  
  city.names = c("Paris", "Madrid", "Stockholm")
  city = city.names[1]
  
  case_name = c("statio", "slp_SD")
  ncase = length(case_name)
}
## Function
source("function/fun_limits.R")
####


#### Estimation ####
source("function/fun_estimation_t2m.R")

## 1979 - 1998
year.begin = 1979
year.end = 1998
for (city in city.names){
  source("function/fun_estimation_tp.R")
}

## 1999 - 2017
year.begin = 1999
year.end = 2017
for (city in city.names){
  source("function/fun_estimation_tp.R")
}
####


#### Simulation ####
source("function/fun_simulation_t2m.R")

## 1979 - 1998
y1 = 1999; y2 = 2017; year.begin = 1979; year.end = 1998; type = "nat"
for (city in city.names){
  source("function/fun_simulation_tp.R")
}

## 1999 - 2017
y1 = 1999; y2 = 2017; year.begin = 1999; year.end = 2017; type = "rea"
for (city in city.names){
  source("function/fun_simulation_tp.R")
}
####


#### SAMPLE 100 simulations ####
if (FALSE){
  for (city in city.names){
    dir.create(path = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SIMU"))
    dir.create(path = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SIMU/nat"))
    dir.create(path = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SIMU/rea"))
  }
}

Sample_binomocc_gammaity<- function(par_ro,par_rf1,par_rf2,thresh=0){
  a=date()
  cat(a,'\n')
  
  Nb_stations=dim(par_ro)[2]
  echant=array(NaN,dim=dim(par_ro))
  for (station in 1:Nb_stations){
    simul_pluie=runif(dim(par_ro)[1],0,1)
    for (j in 1:(dim(par_ro)[1])){
      if (simul_pluie[j]<par_ro[j,station]){
        p=runif(1)
        C=pgamma(thresh,shape=par_rf1[j,station],scale=1/par_rf2[j,station])
        echant[j,station]=qgamma(p*(1-C)+C,shape=par_rf1[j,station],scale=1/par_rf2[j,station])
      }
      else {echant[j,station]=0}
    }
  }
  a=date()
  cat(a,'\n')
  return(echant)
}

for (city in city.names){
  for (w in c("nat", "rea")){
    for (k in 1:1){
      load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/precip_P1_shape_rate_1999_2017_", w, "_1.RData"))
      DATE = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))
      
      ## make 100 runs
      idx = 1
      for(i in 1:100){
        
        NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
        cat('Run',NUM,'\n')
        
        Sample_temp = Sample_binomocc_gammaity(par_ro = P1,par_rf1 = shape,par_rf2 = rate)
        Sample=array(NaN,dim=dim(P1))
        Sample[,idx]=Sample_temp
        
        Sample=round(Sample,digits=2)
        
        filesimu = paste('~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/', city, '/SIMU/', w, '/SIMU_SWG_precip_1999_2017_', w, '_', 1, '_run_', NUM, '.RData', sep="")
        save(Sample,DATE,file = filesimu)
      }
    }
  } 
}
####


#### separator in season ####
for (k in 1:1){
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SIMU/nat/SIMU_SWG_precip_1999_2017_nat_", 1, "_run_", NUM, ".RData"))
    for (seas in 1:4){
      Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
      DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
      save(DATE_seas, Sample_seas, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SIMU/nat/SIMU_precip_1999_2017_nat_1_", SEAS[seas], "_run_", NUM, ".RData"))
    }
  }
}
####


#### create stationnary model: DON'T RUN every time !!! ####
for (city in city.names){
  load(paste0("~/Documents/LSCE/SWG/data/tp_", city, "_1979_2017.RData"))
  DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  
  data1 = tp[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
  DATE1 = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))
  
  data2 = tp[DATE_OBS$Y>=1979 & DATE_OBS$Y<=1998,]
  DATE2 = atoms(timeSequence(from="1979-01-01",to="1998-12-31",by='day'))
  
  ## nat
  P1 = array(NaN,c(length(data1),1)) # matrice de probas d'occurrence
  shape = array(NaN,c(length(data1),1))
  rate = array(NaN,c(length(data1),1))
  
  for (seas in 1:4){
    idxdates1 = which(DATE1['m']==MON[seas,1] | DATE1['m']==MON[seas,2] | DATE1['m']==MON[seas,3])
    idxdates2 = which(DATE2['m']==MON[seas,1] | DATE2['m']==MON[seas,2] | DATE2['m']==MON[seas,3])
    
    P1[idxdates1, ]   = sum(data2[idxdates2]>0)/length(data2[idxdates2]);
    fit.gamma <- fitdist(data2[idxdates2][data2[idxdates2]>0], distr = "gamma", method = "mle")
    shape[idxdates1,] = fit.gamma$estimate[1];
    rate[idxdates1,]  = fit.gamma$estimate[2]
  }
  DATE = DATE1
  save(P1, shape, rate, DATE, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/precip_P1_shape_rate_1999_2017_nat_0.RData"))
  
  ## rea
  P1 = array(NaN,c(length(data1),1)) # matrice de probas d'occurrence
  shape = array(NaN,c(length(data1),1))
  rate = array(NaN,c(length(data1),1))
  
  for (seas in 1:4){
    idxdates1 = which(DATE1['m']==MON[seas,1] | DATE1['m']==MON[seas,2] | DATE1['m']==MON[seas,3])
    
    P1[idxdates1, ]   = sum(data1[idxdates1]>0)/length(data1[idxdates1]);
    fit.gamma <- fitdist(data1[idxdates1][data1[idxdates1]>0], distr = "gamma", method = "mle")
    shape[idxdates1,] = fit.gamma$estimate[1];
    rate[idxdates1,]  = fit.gamma$estimate[2]
  }
  DATE = DATE1
  save(P1, shape, rate, DATE, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/precip_P1_shape_rate_1999_2017_rea_0.RData")) 
}
####

#### FAR ####
city.names = c("Paris", "Madrid", "Stockholm")
ens_thres = 0:30
FAR_cases = c("FAR", "RR", "FARanda")
ens_year = 1999:2017

for (city in city.names){
  
  # city = city.names[1]
  load(paste0("~/Documents/LSCE/SWG/data/tp_", city, "_1979_2017.RData"))
  DATE_OBS = atoms(timeSequence(from="1979-01-01", to="2017-12-31", by='day'))
  
  list_P1_shape_rate_nat = vector("list", 2)
  list_P1_shape_rate_rea = vector("list", 2)
  for (k in 1:2){
    i = k-1
    load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/precip_P1_shape_rate_1999_2017_nat_", i,".RData"))
    cat("nat",'\n')
    cat(range(P1),'\n')
    cat(range(shape),'\n')
    cat(range(rate),'\n')
    m = as.data.frame(matrix(NA, nrow = nrow(P1), ncol = 3))
    colnames(m) = c("P1", "shape", "rate")
    m[,"P1"]    = P1
    m[,"shape"] = shape
    m[,"rate"]  = rate
    list_P1_shape_rate_nat[[k]] = m
    
    load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/precip_P1_shape_rate_1999_2017_rea_", i,".RData"))
    cat("rea",'\n')
    cat(range(P1),'\n')
    cat(range(shape),'\n')
    cat(range(rate),'\n')
    m = as.data.frame(matrix(NA, nrow = nrow(P1), ncol = 3))
    colnames(m) = c("P1", "shape", "rate")
    m[,"P1"]    = P1
    m[,"shape"] = shape
    m[,"rate"]  = rate
    list_P1_shape_rate_rea[[k]] = m
  }
  
  Obs = tp[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
  n = length(Obs)
  DATE = as.Date((0:(n-1)), origin = "1999-01-01")
  
  # ens_thres = 1:20
  list_prob_thres = vector("list", length(ens_thres))
  names(list_prob_thres) = ens_thres
  for (thres in ens_thres){
    
    # thres = 0
    list_prob_thres_nat = lapply(list_P1_shape_rate_nat, FUN = function(m){
      1 - (pgamma(thres, shape = m[, "shape"], scale = 1/m[,"rate"], lower.tail = T)*m[,"P1"] +
           pgamma(0, shape = m[, "shape"], scale = 1/m[,"rate"], lower.tail = T)*(1-2*m[,"P1"]))
    })
    prob_thres_nat = do.call(rbind.data.frame, list_prob_thres_nat)
    colnames(prob_thres_nat) = "nat"
    
    list_prob_thres_rea = lapply(list_P1_shape_rate_rea, FUN = function(m){
      1 - (pgamma(thres, shape = m[, "shape"], scale = 1/m[,"rate"], lower.tail = T)*m[,"P1"] +
             pgamma(0, shape = m[, "shape"], scale = 1/m[,"rate"], lower.tail = T)*(1-2*m[,"P1"]))
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
  
  for (FAR_case in FAR_cases[1]){
    for (year in ens_year){
      
      # FAR_case = "FAR"
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
      
      FAR_sel = FAR_sel[, c("date", FAR_case, "case", "thres", "year", "num")]
      colnames(FAR_sel) = c("date", "value", "case", "thres", "year", "num")
      
      saveGIF(
        for (thres in ens_thres){
          p <- ggplot(FAR_sel, aes(x=num, y=value, group=case, col=case)) + 
            geom_blank() +
            geom_line(data = FAR_sel[FAR_sel$thres==thres,], mapping = aes(x=num, y=value, group=case, col=case)) + theme_bw() +
            scale_color_manual(values = cc, name = "") + 
            scale_x_continuous(breaks = breaks, labels = labels) +
            labs(x = "date", y = FAR_case) +
            ggtitle(paste0(FAR_case, " for summer in ", year, " over ", thres, " mm"))
          print(p)
          # ggsave(plot = p, filename = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", i, ".png"), units = "in", dpi = 300, width = 7, height = 7)
          # dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", year, ".pdf"), width = 7, height = 7)
          # dev.print(png, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", i, ".png"), width = 700, height = 700)
        }, 
        ani.width = 200*7, ani.height = 200*7, ani.res = 200, interval = 0.5, autobrowse = F,
        movie.name = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/Image/FAR/", FAR_case, "_summer_in_", year, ".gif")
      )
    }
  }
}

#### Intensity: SAVE d_intensity.RData ####
## load tmean OBS
load(paste0("~/Documents/LSCE/SWG/data/tp_", city, "_1979_2017.RData"))
DATE_OBS = atoms(timeSequence(from="1979-01-01", to="2017-12-31", by='day'))

Obs = tp[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
n = length(Obs)
DATE = as.Date((0:(n-1)), origin = "1999-01-01")

list_P1_shape_rate_nat = vector("list", 2)
list_P1_shape_rate_rea = vector("list", 2)
for (k in 1:2){
  i = k-1
  load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/precip_P1_shape_rate_1999_2017_nat_", i,".RData"))
  cat("nat",'\n')
  cat(range(P1),'\n')
  cat(range(shape),'\n')
  cat(range(rate),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(P1), ncol = 3))
  colnames(m) = c("P1", "shape", "rate")
  m[,"P1"]    = P1
  m[,"shape"] = shape
  m[,"rate"]  = rate
  list_P1_shape_rate_nat[[k]] = m
  
  load(paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/precip_P1_shape_rate_1999_2017_rea_", i,".RData"))
  cat("rea",'\n')
  cat(range(P1),'\n')
  cat(range(shape),'\n')
  cat(range(rate),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(P1), ncol = 3))
  colnames(m) = c("P1", "shape", "rate")
  m[,"P1"]    = P1
  m[,"shape"] = shape
  m[,"rate"]  = rate
  list_P1_shape_rate_rea[[k]] = m
}

p1 = as.data.frame(matrix(NA, nrow = n, ncol = ncase))
for (l in 1:ncase){
  for (k in 1:n){
    if (Obs[k]==0){
      p1[k,l] = list_P1_shape_rate_rea[[l]][k,"P1"]
    }else{
      p1[k,l] = pgamma(Obs[k], shape = list_P1_shape_rate_rea[[l]][k,"shape"], rate = list_P1_shape_rate_rea[[l]][k,"rate"], lower.tail = T)
    }
  }
}

p0 = p1

i1 = as.data.frame(matrix(NA, nrow = n, ncol = ncase))
for (l in 1:ncase){
  for (k in 1:n){
    i1[k,l] = qnorm(p0[k,l], mean = list_mean_sd_nat[[l]][k,"mean"], sd = list_mean_sd_nat[[l]][k,"sd"], lower.tail = T)
  }
}

## d_intensity: 1999-2017
d = Obs - i1; range(d)
colnames(d) = c("statio", "slp_SD")
DATE = as.Date((0:(n-1)), origin = "1999-01-01")
d$date = DATE
save(d, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/d_intensity.RData"))


#### Intensity: plot for each year from 1999 to 2017 ####
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


#### Anomaly slp summer ####
load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/d_intensity.RData"))
dI = d

var = "slp"; 
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

col.range = c(-171, 171)
{
  t_sup = which(dI[,"slp_SD"] - d[,"statio"]>0); length(t_sup)
  t_inf = which(dI[,"slp_SD"] - d[,"statio"]<0); length(t_inf)
  
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
    
    output = paste0(var, "_sup_1999-2017_", season[seas], "_slp_SD")
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output,".pdf"), width = 11, height = 5)
    
    
    t_inf_sea = intersect(t_inf, idxdates)
    data = apply(data_anomaly[,,t_inf_sea], 1:2, mean)
    colnames(data) = LAT; rownames(data) = LON
    data = data[,rev(seq_along(LAT))]
    print(range(data))
    
    col.range = c(-max(abs(range(data[!is.na(data)]))), max(abs(range(data[!is.na(data)]))))
    cc <- colorRamps::matlab.like2(n = 100)
    
    output = paste0(var, "_inf_1999-2017_", season[seas], "_slp_SD")
    image.plot(data, x = LON, y = rev(LAT), xlab = "Lon", ylab = "Lat", col = cc, zlim = col.range,
               main = output)
    world_map = map('world', plot = FALSE)
    lines(world_map, col = "white")
    contour(x = LON, y = rev(LAT), z = data, add = TRUE)
    dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output,".pdf"), width = 11, height = 5)
  }
}


#### dim and theta depend on precipitation ####
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
plot(dim, theta, main = paste0("theta vs dim of ", var))

for (city in city.names){
  data = as.data.frame(cbind(dim, theta))
  load(paste0("~/Documents/LSCE/SWG/data/tp_", city, "_1979_2017.RData"))
  data = cbind(data, tp)
  
  ## all seasons 
  output = paste0("dim_theta_with_precip_of_", city, "_1999-2017")
  p = ggplot(data) + geom_point(aes(x = dim, y = theta, col = tp)) + theme_bw()
  p = p + scale_color_gradientn(colours = matlab.like2(10), limits = fun_limits(data,min0=T))
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/Image/", output, ".pdf"), width = 7, height = 7)
  
  ## all seasons with threshold
  thres = 15
  if (T){
    output = paste0("dim_theta_with_tmean_of_", city, "_1999-2017_over_", thres)
  }else{
    output = paste0("dim_theta_with_tmean_of_", city, "_1999-2017_under_", thres)
  }
  print(output)
  
  p <- ggplot(data, mapping = aes(x = dim, y = theta)) + geom_blank() +
    geom_point(data = data[data$tp>=thres,], mapping = aes(x = dim, y = theta, col = tp)) + 
    theme_bw() + scale_color_gradientn(colours = matlab.like2(10), limits = fun_limits(data, min0=T))
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/Image/", output, ".pdf"), width = 7, height = 7)
  
  ## each season 
  # for (seas in 1:4){
  #   idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
  #   data_sea = data[idxdates,]
  #   
  #   output = paste0("dim_theta_with_tmean_of_", city, "_1999-2017_", season[seas])
  #   p = ggplot(data_sea) + geom_point(aes(x = dim, y = theta, col = tmean)) + theme_bw()
  #   p = p + scale_color_gradientn(colours = matlab.like2(10), limits = fun_limits(data_sea))
  #   print(p)
  #   dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/", output, ".pdf"), width = 7, height = 7)
  # } 
}
####

#### max annual precip ####
for (city in city.names){
  
  ## annual_max_in_plot_theta_dim
  load(file = paste0("~/Documents/LSCE/SWG/data/tp_", city, "_1979_2017.RData"))
  DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
  D = as.Date(0:(nrow(DATE_OBS)-1), origin = "1979-01-01")
  
  pam = as.data.frame(matrix(NA, nrow = length(unique(DATE_OBS$Y)), ncol = 2))
  colnames(pam) = c("pam", "date")
  rownames(pam) = unique(DATE_OBS$Y)
  for (y in 1979:2017){
    Dy = D[DATE_OBS$Y==y][which(tp[DATE_OBS$Y==y,] == max(tp[DATE_OBS$Y==y,]))]
    pam[paste(y),] = c(max(tp[DATE_OBS$Y==y,]), as.character(Dy))
  }
  
  p <- ggplot(data) + theme_bw() + 
    geom_point(aes(x = dim, y = theta)) + 
    geom_point(data = data[as.character(D)%in%pam$date, ], aes(x = dim, y = theta),  col = "red")
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/Image/annual_max_in_plot_theta_dim.pdf"), width = 7, height = 7)
  
  ## tp_max_annual_city_1999_2017_boxplot_rea
  pam.list = vector("list", 100)
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SIMU/rea/SIMU_SWG_precip_1999_2017_rea_1_run_", NUM, ".RData"))
    pam_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
    colnames(pam_sim) = c("pam", "year")
    rownames(pam_sim) = unique(DATE$Y)
    for (y in unique(DATE$Y)){
      pam_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
    }
    pam.list[[i]] = pam_sim
  }
  pam_sim = do.call(rbind.data.frame, pam.list)
  pam_sim$type = "sim"
  
  pam_obs = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
  colnames(pam_obs) = c("pam", "year")
  rownames(pam_obs) = unique(DATE$Y)
  pam_obs[,1] = as.numeric(pam[as.character(unique(DATE$Y)),1])
  pam_obs[,2] = unique(DATE$Y)
  pam_obs$type = "obs"
  
  pam_data = rbind(pam_obs, pam_sim)
  
  output = paste0("tp_max_annual_", city, "_1999_2017_boxplot_rea")
  p <- ggplot(pam_data) + theme_bw() +
    geom_boxplot(data = pam_data[pam_data$type=="sim",], mapping = aes(x = year, y = pam, group = year)) +
    geom_point(data = pam_data[pam_data$type=="obs",], mapping = aes(x = year, y = pam), col = "red") + 
    ggtitle(output)
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/Image/", output, ".pdf"), width = 7, height = 7)
  
  ## tp_max_annual_city_1999_2017_boxplot_nat
  pam.list = vector("list", 100)
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    load(file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/SIMU/nat/SIMU_SWG_precip_1999_2017_nat_1_run_", NUM, ".RData"))
    pam_sim = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
    colnames(pam_sim) = c("pam", "year")
    rownames(pam_sim) = unique(DATE$Y)
    for (y in unique(DATE$Y)){
      pam_sim[paste(y),] = c(max(Sample[DATE$Y==y,]), y)
    }
    pam.list[[i]] = pam_sim
  }
  pam_sim = do.call(rbind.data.frame, pam.list)
  pam_sim$type = "sim"
  
  pam_obs = as.data.frame(matrix(NA, nrow = length(unique(DATE$Y)), ncol = 2))
  colnames(pam_obs) = c("pam", "year")
  rownames(pam_obs) = unique(DATE$Y)
  pam_obs[,1] = as.numeric(pam[as.character(unique(DATE$Y)),1])
  pam_obs[,2] = unique(DATE$Y)
  pam_obs$type = "obs"
  
  pam_data = rbind(pam_obs, pam_sim)
  
  output = paste0("tp_max_annual_", city, "_1999_2017_boxplot_nat")
  p <- ggplot(pam_data) + theme_bw() +
    geom_boxplot(data = pam_data[pam_data$type=="sim",], mapping = aes(x = year, y = pam, group = year)) +
    geom_point(data = pam_data[pam_data$type=="obs",], mapping = aes(x = year, y = pam), col = "red") + 
    ggtitle(output)
  print(p)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/tp/", city, "/Image/", output, ".pdf"), width = 7, height = 7)
}
####


# for (k in 1:3){
#   city = city.names[k]
#   system2('open', args = paste0("~/Documents/LSCE/SWG/slp_SD_diagnosis/", city, "/Image/annual_max_in_plot_theta_dim.pdf"), wait = F)
# }


