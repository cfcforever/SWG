# load packages -----------------------------------------------------------
{
  library(ncdf4)
  
  # atoms
  library(timeDate)
  
  library(stats4)
  library(VGAM)
  
  library(CDFt)
  
  library(ggplot2)
  library(reshape2)
  
  library(evd)
  
  library(fields)
  
  library(lubridate)
  
  library(fitdistrplus)
  
  library(xtable)
  
  library(RColorBrewer)
  library(colorRamps)
  
  library(factoextra)  
  
  # install.packages("devtools")  
  # library(jcolors)
  
  library(animation)
  
  # anomaly_map
  library(maps)
  library(directlabels)
  
  # multiplot
  library(cowplot)
  
  # conditional gev
  library(extRemes)
  
  # crps
  library(scoringRules)
}

# load function -----------------------------------------------------------
source("function/fun_save_pdf.R")
source("function/fun_limits.R")

source("function/fun_estimation_t2m.R")
source("function/fun_estimation_t2m.R")

source("function/plot_worldmap.R")
