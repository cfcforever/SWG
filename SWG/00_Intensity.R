#### load tmean OBS ####
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

Obs = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
n = length(Obs)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

mDATE = atoms(timeSequence(from="1999-01-01",to="2017-12-31",by='day'))
n = nrow(mDATE)

DATE = as.Date((0:(n-1)), origin = "1999-01-01")


#### Intensity for t2m from 1998 to 2017 ####
case_name = c(1:64); case_sel = c(1:64)
case_sel = c(63, 61, 57, 48, 59, 50, 43, 27, 64)   # best BIC
case_sel = c(1, 2, 4, 7, 9, 13, 23, 64)            # best intensity
case_sel = c(25,15,11, 64)            # best intensity
case_name = case_name[case_sel]
ncase = length(case_name)


list_mean_sd_nat = vector("list", ncase)
list_mean_sd_rea = vector("list", ncase)
for (k in 1:ncase){
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/nat/t2m_mean_sd_1999_2017_nat_", case_sel[k], ".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/rea/t2m_mean_sd_1999_2017_rea_", case_sel[k], ".RData"))
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
    p1[k,l] = pnorm(Obs[k], mean = list_mean_sd_rea[[l]][k,"mean"], sd = list_mean_sd_rea[[l]][k,"sd"], lower.tail = T)
  }
}

p0 = p1

i1 = as.data.frame(matrix(NA, nrow = n, ncol = ncase))
for (l in 1:ncase){
  for (k in 1:n){
    i1[k,l] = qnorm(p0[k,l], mean = list_mean_sd_nat[[l]][k,"mean"], sd = list_mean_sd_nat[[l]][k,"sd"], lower.tail = T)
  }
}

d = Obs - i1; range(d)
colnames(d) = case_name
DATE = as.Date((0:(n-1)), origin = "1999-01-01")
d$date = DATE
# save(d, file = "~/Documents/LSCE/SWG/NA/Intensity_3.RData")
save(d, file = "~/Documents/LSCE/SWG/NA_1.5/tmean/d_intensity.RData")

load("~/Documents/LSCE/SWG/NA_1.5/tmean/d_intensity.RData")
library(lubridate)
## d_intensity: year - season - model
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
  p = p + scale_x_continuous(breaks = breaks, labels = labels) 
  p = p + labs(title = paste0("delta_Intensity_tmean_summer_in_", year),
               x = "date", y = "difference of intensity (°C)")
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/Image/Intensity/d_intensity_tmean_summer_", year, "_3.pdf"), width = 7, height=7)
  print(year)
}


## d_year_1
load(file = "~/Documents/LSCE/SWG/NA_1.5/tmean/d_intensity.RData")

d_year = matrix(NA, nrow = (2017-1999+1), ncol = ncase)
yy = unique(year(d$date))
for (y in 1:nrow(d_year)){
  d_year[y,] = colMeans(d[year(d$date)==yy[y],1:ncase])
}
d_year = as.data.frame(d_year)
load(file = "~/Documents/LSCE/SWG/NA/case_name.RData")
colnames(d_year) = c(case_name, "statio")[case_sel]
d_year$year = 1999:2017

d_year_1 = melt(d_year, id.vars = "year")
colnames(d_year_1) = c("year", "case", "intensity")

## annual delta of intensity from 1999 to 2017
p = ggplot(d_year_1)
p = p + geom_line(aes(x=year, y=intensity, col=case)) + theme_bw()
p = p + geom_line(data = d_year_1[d_year_1$case=="statio",], aes(x = year, y=intensity), col="black")
p = p + labs(title = "annual delta of intensity from 1999 to 2017",
             x = "date", y = "difference of intensity (°C)")
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/Image/Intensity/annual delta of intensity from 1999 to 2017.pdf"), width = 20, height=7)


## annual delta of intensity from 1999 to 2016 for case XXX
for (k in 32:63){
  p = ggplot(d_year_1[d_year_1$case == case_name[k],])
  p = p + geom_line(aes(x=year, y=intensity), col="red") + theme_bw()
  p = p + geom_line(data = d_year_1[d_year_1$case=="statio",], aes(x = year, y=intensity), col="black")
  p = p + labs(title = paste0("annual delta of intensity from 1999 to 2017 by ", case_name[k]),
               x = "date", y = "difference of intensity (°C)")
  print(p)
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/Image/Intensity/annual delta of intensity from 1999 to 2017 for ", k, ".pdf"), width = 7, height=7)
}

## average comparing to stationary model
d_year_mean = as.data.frame(colMeans(d_year[,1:64]))
colnames(d_year_mean) = "average"
d_year_mean$number = 1:64
d_year_mean$model = rownames(d_year_mean)

d_year_mean$model = factor(d_year_mean$model, levels = rownames(d_year_mean))
# d_year_mean$model = factor(d_year_mean$model, levels = 1:64)

d_year_mean_to_statio = d_year_mean
d_year_mean_to_statio$average = abs(d_year_mean$average - d_year_mean["statio", 1])
mods = which(d_year_mean_to_statio$average %in% head(sort(d_year_mean_to_statio$average, decreasing = F), 6))
# 1  2  7 13 49 64 

p = ggplot(d_year_mean) 
p = p + geom_text(aes(x=number, y=average, label = number)) +  theme_bw()
p = p + geom_point(aes(x=number, y=average), col = "red", size = 0.5)
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/Image/Intensity/average_annual_delta_of_intensity_1999-2017.pdf"), width = 9, height=7)

p = ggplot(d_year_mean[c(1,2,4,7,9,13,23,64),]) 
p = p + geom_point(aes(x=number, y=average)) +  theme_bw()
p = p + geom_text(aes(x=number, y=average, label=model), hjust = 0, nudge_x = 0.5) 
# p = p + geom_text(aes(x=number, y=average, label=number), hjust = 0, nudge_x = -2) 
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/Image/Intensity/average_annual_delta_of_intensity_1999-2017_2.pdf"), width = 9, height=7)

mods = c(64,63,61,57,48,59,50,43,27)
p = ggplot(d_year_mean[mods,]) 
p = p + geom_point(aes(x=number, y=average, col=model)) +  theme_bw()
# p = p + geom_text(aes(x=number, y=average, label=model), hjust = 0, nudge_x = 0.5) 
print(p)
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/NA_1.5/Image/Intensity/average_annual_delta_of_intensity_1999-2017_best_BIC.pdf"), width = 9, height=7)

