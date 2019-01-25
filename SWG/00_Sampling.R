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

##
k = 1
load(paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/nat/t2m_mean_sd_1999_2017_nat_", k, ".RData"))
# sum(is.na(mean))
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
  
  filesimu = paste('~/Documents/LSCE/SWG/NA_1.5/tmean/SIMU/', k, '/SIMU_SWG_tmean_1999_2017_nat_', k, '_run_', NUM, '.RData', sep="")
  save(Sample,DATE,file = filesimu)
}

