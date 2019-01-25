MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

DATE = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
load(file = "~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
tmean = tmean[DATE$Y>=1979 & DATE$Y<=2017,]
DATE = DATE[DATE$Y>=1979 & DATE$Y<=2017,]

# casename = c("1_msl_PC", "2_msl_SD", "3_z500_PC", "4_z500_SD", 
#              "1+2", "1+3", "1+4", "2+3", "2+4", "3+4",
#              "1+2+3", "1+2+4", "1+3+4", "2+3+4", "1+2+3+4")

casename = c(1:64)

ncase = length(casename)

BIC_table = matrix(NA, nrow = 4, ncol = 64)
for (seas in 1:4){
  # seas = 1
  cat(season[seas],'\n')
  
  BIC_sea = rep(NaN, ncase)
  
  idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
  temp = tmean[idxdates]
  m = mean(temp)
  sd = sd(temp)
  cat(m, sd, '\n')
  # f1n <- fitdistr(temp,"normal")
  # cat(coef(f1n), '\n')
  loglik = sum(dnorm(temp, mean = m, sd = sd, log = TRUE))
  BIC_sea[64] = 2*log(length(idxdates)) - 2*f1n$loglik
  
  for (NUM in 1:63){
    # load(paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/1979_1998/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1979_1998_", NUM, ".RData"))
    load(paste0("~/Documents/LSCE/SWG/NA_1.5//tmean/1979_2017/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1979_2017_", NUM, ".RData"))
    load(paste0("~/Documents/LSCE/SWG/NA_1.5/pca/predictor_", SEAS[seas], "_", NUM, ".RData"))
    # load(paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/1999_2017/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1999_2017_", NUM, ".RData"))
    
    npc = ncol(PCS)
    BIC = (2+npc*2)*log(length(temp)) - 2*logLik(fit_stations_tt[[1]]) 
    
    BIC_sea[NUM] = BIC
  }
  
  for (k in 1:64){
    BIC_table[seas, k] = which(rank(BIC_sea)==k)
  }
  
  BIC = BIC_sea
  BIC = as.data.frame(BIC)
  BIC$type = casename[1:ncase]
  BIC$type <- factor(BIC$type, levels = casename)
  
  p = ggplot(BIC, aes(x=type, y=BIC)) + geom_point() + theme_bw()
  p = p + labs(title = paste0("BIC", season[seas]))
  print(p)
  # dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA/Image/BIC/BIC_", season[seas], "_1979_1998.pdf"), width = 14, height = 7)
  dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/Image/BIC/BIC_", season[seas], "_1979_2017.pdf"), width = 14, height = 7)
  # dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA/Image/BIC/BIC_", season[seas], "_1999_2017.pdf"), width = 14, height = 7)
}

BIC_rank = as.data.frame(matrix(NA, nrow = 64, ncol = 1))
colnames(BIC_rank) = "BIC"

load(file = "~/Documents/LSCE/SWG/NA/case_name.RData")
BIC_rank$model = case_name
# BIC_rank$model = c(case_name, "statio")

# all 4 seasons 
for (k in 1:64){
  BIC_rank[k,1] = which(BIC_table[1,]==k) + which(BIC_table[2,]==k) + which(BIC_table[3,]==k) + which(BIC_table[4,]==k)
}
save(BIC_rank, file = "~/Documents/LSCE/SWG/NA_1.5/BIC_rank.RData")
load("~/Documents/LSCE/SWG/NA_1.5/BIC_rank.RData")

# just summer 
for (k in 1:64){
  BIC_rank[k,1] = which(BIC_table[3,]==k)
}
save(BIC_rank, file = "~/Documents/LSCE/SWG/NA_1.5/tmean/BIC_rank_summer.RData")

BIC_rank = BIC_rank[order(BIC_rank$BIC, decreasing = F), ]
xtable(BIC_rank)
