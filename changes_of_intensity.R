#### Create d_intensity.RData ------------------------------------------------
load(paste0("DATA/tmean_", city, "_1979_2017.RData"))
DATE_OBS = atoms(timeSequence(from="1979-01-01", to="2017-12-31", by='day'))

Obs = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
n = length(Obs)
DATE = as.Date((0:(n-1)), origin = "1999-01-01")

list_mean_sd_nat = vector("list", 2)
list_mean_sd_rea = vector("list", 2)
for (k in 1:2){
  i = k-1
  load(paste0("output/", city, "/t2m_mean_sd_1999_2017_nat_", i,".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  load(paste0("output/", city, "/t2m_mean_sd_1999_2017_rea_", i,".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_rea[[k]] = m
}

p1 = as.data.frame(matrix(NA, nrow = n, ncol = 2))
for (l in 1:2){
  for (k in 1:n){
    p1[k,l] = pnorm(Obs[k], mean = list_mean_sd_rea[[l]][k,"mean"], sd = list_mean_sd_rea[[l]][k,"sd"], lower.tail = F)
  }
}

p0 = p1

i1 = as.data.frame(matrix(NA, nrow = n, ncol = 2))
for (l in 1:2){
  for (k in 1:n){
    i1[k,l] = qnorm(1-p0[k,l], mean = list_mean_sd_nat[[l]][k,"mean"], sd = list_mean_sd_nat[[l]][k,"sd"], lower.tail = T)
  }
}

## d_intensity: 1999-2017
d = Obs - i1; range(d)
colnames(d) = c("stationary", "conditional")
DATE = as.Date((0:(n-1)), origin = "1999-01-01")
d$date = DATE
save(d, file = paste0("output/", city, "/d_intensity.RData"))
####

