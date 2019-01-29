# load packages -----------------------------------------------------------
{
  library(ncdf4)
  
  # atoms
  library(timeDate)
  
  library(stats4)
  library(VGAM)
  
  # CramerVonMisesTwoSamples
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
# source("function/fun_save_pdf.R")
# source("function/fun_limits.R")

source("function/fun_estimation_t2m.R")
source("function/fun_simulation_t2m.R")

source("function/fun_Sample_gaussian_temp.R")

source("function/plot_worldmap.R")


# load common variable ----------------------------------------------------
SEAS   = c('DJF','MAM','JJA','SON')
season = c('Winter','Spring','Summer','Fall')
MON    = c(12, 1:11)
MON    = matrix(MON, nrow=4, ncol=3, byrow=T)

