MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)
SEAS   = c('DJF', 'MAM', 'JJA', 'SON')

for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(paste0("~/Documents/LSCE/SWG/t2m/case_0/SIMU/SIMU_SWG_tmean_run_", NUM, ".RData"))
  for (seas in 1:4){
    Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
    DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
    save(DATE_seas, Sample_seas, file = paste0("~/Documents/LSCE/SWG/t2m/case_0/SIMU/SIMU_tmean_", SEAS[seas], "_run_", NUM, ".RData"))
  }
}

for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(paste0("~/Documents/LSCE/SWG/tp/case_3/SIMU/SIMU_SWG_precip_run_", NUM, ".RData"))
  for (seas in 1:4){
    Sample_seas = Sample[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
    DATE_seas   = DATE[DATE$m==MON[seas,1] | DATE$m==MON[seas,2] | DATE$m==MON[seas,3], ]
    save(DATE_seas, Sample_seas, file = paste0("~/Documents/LSCE/SWG/tp/case_3/SIMU/SIMU_precip_", SEAS[seas], "_run_", NUM, ".RData"))
  }
}


