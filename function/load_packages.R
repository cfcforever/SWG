# set study case and directory --------------------------------------------
if (!dir.exists("output")){
  dir.create("output")
}

# /!\ CHOOSE your case !!! 
case = 3

if (case == 1){
  SD = "slp_SD"
  output.dir = paste0("output/", SD, "/")   # where store output
  dim.dir   = "DATA/dim_theta/msl_1979-2018_NA_1.5x1.5_dim.RData"
  theta.dir = "DATA/dim_theta/msl_1979-2018_NA_1.5x1.5_theta.RData"
}else if(case == 2){
  SD = "slp_z500_SD"    
  output.dir = paste0("output/", SD, "/")
  dim.dir   = "DATA/dim_theta/msl_z500_19790101-20180731_NA_1.5x1.5_dim.RData"
  theta.dir = "DATA/dim_theta/msl_z500_19790101-20180731_NA_1.5x1.5_theta.RData"
}else if(case == 3){
  SD = "slp_z500_SD_scale"
  output.dir = paste0("output/", SD, "/")  
  dim.dir   = "DATA/dim_theta/msl_z500_19790101-20180731_NA_1.5x1.5_dim_scale.RData"
  theta.dir = "DATA/dim_theta/msl_z500_19790101-20180731_NA_1.5x1.5_theta_scale.RData"
}



# load packages -----------------------------------------------------------
{
  # maybe useless
  library(utils)
  
  # NetCDF
  library(ncdf4)
  
  # atoms
  library(timeDate)
  
  library(stats4)
  library(VGAM)
  
  # # CramerVonMisesTwoSamples
  # library(CDFt)
  
  library(ggplot2)
  library(reshape2)
  
  # library(evd)
  
  library(fields)
  
  library(lubridate)
  
  # library(fitdistrplus)
  
  # library(xtable)
  
  library(RColorBrewer)
  library(colorRamps)
  
  # library(factoextra)  
  
  # install.packages("devtools")  
  # library(jcolors)
  
  # library(animation)
  
  # anomaly_map
  library(maps)
  library(directlabels)
  
  # multiplot
  # library(cowplot)
  
  # conditional gev
  library(extRemes)
  
  # crps
  library(scoringRules)
  
  # rainfall model P.Naveau
  library(mev)
}

# load function -----------------------------------------------------------
source("function/fun_estimation_t2m.R")
source("function/fun_simulation_t2m.R")

source("function/fun_Sample_gaussian_temp.R")

source("function/plot_worldmap.R")


# load common variable ----------------------------------------------------
SEAS   = c('DJF','MAM','JJA','SON')
season = c('Winter','Spring','Summer','Fall')
MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)







































