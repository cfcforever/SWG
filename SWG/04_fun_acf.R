# load tmean
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")   

for (seas in 1:4){
  # Obs as tmean - daily mean temperature from 2005 to 2017 in the zone [48N, 50N]x[01E, 04E]
  Obs = tmean[which(rownames(tmean)=="2005-01-01"):which(rownames(tmean)=="2017-12-31"),]
  DATE_Obs = atoms(timeSequence(from="2005-01-01",to="2017-12-31",by='day'))
  
  MON    = c(12, 1:11)
  MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
  SEAS   = c('DJF', 'MAM', 'JJA', 'SON')
  Obs = Obs[DATE_Obs$m==MON[seas,1] | DATE_Obs$m==MON[seas,2] | DATE_Obs$m==MON[seas,3]]
  
  # NUM: number with 3 digits
  k = 1 
  NUM=cbind(formatC(k, digits = 0, width = 3, format = "f", flag = "0"))
  
  for (ncase in 1:3){
    # load one Sample (see k and NUM)
    load(paste0("~/Documents/LSCE/SWG/t2m/case_", ncase, "/SIMU/SIMU_tmean_", SEAS[seas], "_run_001.RData"))
    Sample = Sample_seas
    
    # set number of lag
    nlag = 31
    
    # data_acf contains acf_obs and acf_sim
    data_acf = as.data.frame(matrix(data = NA, nrow = 2*nlag, ncol = 3))
    
    acf_obs = acf(Obs, plot = F, lag.max = nlag-1)
    acf_obs = with(acf_obs, data.frame(lag, acf))
    acf_obs = acf_obs[-1,]
    acf_sim = acf(Sample, plot = F, lag.max = nlag-1)
    acf_sim = with(acf_sim, data.frame(lag, acf))
    acf_sim = acf_sim[-1,]
    data_acf = rbind(acf_obs, acf_sim)
    data_acf$type = c(rep("obs", nlag-1), rep("sim", nlag-1))
    
    color = c("black", "red")
    {
      p = ggplot(data = data_acf, mapping = aes(x = lag, y = acf))
      p = p + geom_bar(stat = "identity", aes(fill = type), position = "dodge")
      p = p + scale_fill_manual(values = color)
      p = p + xlab("day")
      p = p + ggtitle(paste0("acf - Obs vs Sample_", k))
      print(p)
    }
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/t2m/Images/acf_", SEAS[seas],"_case_", ncase, ".pdf"), width = 7, height = 7)
    
  }
}
