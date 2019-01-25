#### FAR for t2m from 1998 to 2017 ####
list_mean_sd_nat = vector("list", 3)
list_mean_sd_rea = vector("list", 3)
for (k in 1:3){
  load(paste0("~/Documents/LSCE/SWG/t2m/case_", k, "/2PC/tmean_mean_sd_1998_2017_nat.RData"))
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  load(paste0("~/Documents/LSCE/SWG/t2m/case_", k, "/2PC/tmean_mean_sd_1998_2017_rea.RData"))
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_rea[[k]] = m
}

load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
Obs = tmean[which(rownames(tmean)=="1998-01-01"):which(rownames(tmean)=="2017-12-31"),]
n = length(Obs)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

mDATE = atoms(timeSequence(from="1998-01-01",to="2017-12-31",by='day'))
n = nrow(mDATE)

DATE = as.Date((0:(n-1)), origin = "1998-01-01")

case_name = c("PCS", "SysDyn", "PCS_SysDyn")


p1 = as.data.frame(matrix(NA, nrow = n, ncol = 3))
for (l in 1:3){
  for (k in 1:n){
    p1[k,l] = pnorm(Obs[k], mean = list_mean_sd_rea[[l]][k,"mean"], sd = list_mean_sd_rea[[l]][k,"sd"], lower.tail = T)
  }
}

p0 = p1

i1 = as.data.frame(matrix(NA, nrow = n, ncol = 3))
for (l in 1:3){
  for (k in 1:n){
    i1[k,l] = qnorm(p0[k,l], mean = list_mean_sd_nat[[l]][k,"mean"], sd = list_mean_sd_nat[[l]][k,"sd"], lower.tail = T)
  }
}


d = Obs - i1
colnames(d) = case_name
d$date = DATE

# year = 2013
for (year in 2011:2017){
  # year = 2005
  period_sel = intersect(which(year(d$date)==year) , which(month(d$date)==7|month(d$date)==6|month(d$date)==8))
  
  d1 = d[period_sel,]
  d1 = melt(d1, id.vars = "date")
  colnames(d1) = c("date", "case", "d")
  d1$num = rep(1:length(period_sel), 3)
  
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
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/t2m/Images/Intensity/2PC/d_intensity_tmean_summer_", year, ".pdf"), width = 7, height=7)
  print(year)
}
