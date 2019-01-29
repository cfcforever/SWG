#### load packages -----------------------------------------------------------
edit(file = "function/load_packages.R")
source("function/load_packages.R")
####



#### Anomaly slp -------------------------------------------------------------
var = "slp"
cas = 3

if (cas == 1){
  load(file = paste0("/Volumes/Data-ExFAT/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1979_2017.RData"))
  cas = "1_1979_2017"
}else if(cas == 2){
  load(file = paste0("/Volumes/Data-ExFAT/LSCE/SWG/NA_slp_sst/", var, "_anomaly_1979_1998.RData"))
  cas = "2_1979_1998"
}else if(cas == 3){
  load(file = paste0("DATA/", var, "_anomaly_1999_2017.RData"))
  cas = "1999_2017"
}

LON = data_anomaly$lon
LAT = data_anomaly$lat
DATE = data_anomaly$date
data_anomaly = data_anomaly$data

data_anomaly = data_anomaly[,,DATE$Y>=1999 & DATE$Y<=2017]
DATE = DATE[DATE$Y>=1999 & DATE$Y<=2017,]

## load d_intensity.RData
load(file = paste0(output.dir, city, "/d_intensity.RData"))
dI = d

## load dim and theta
load(dim.dir)
load(theta.dir)
dim_theta = cbind.data.frame(dim, theta)

{
  t_sup = which(dI[,"conditional"] - dI[,"stationary"]>0); length(t_sup)
  t_inf = which(dI[,"conditional"] - dI[,"stationary"]<0); length(t_inf)
  
  for (seas in c(1,3)){
    # seas = 1
    idxdates = which(DATE['m']==MON[seas,1] | DATE['m']==MON[seas,2] | DATE['m']==MON[seas,3])
    
    ## SUP
    {
      t_sup_sea = intersect(t_sup, idxdates)
      
      output = paste0(var, "_sup_1999-2017_", season[seas], "_", SD, "_dim_theta")
      p <- ggplot(dim_theta, aes(x = dim, y = theta)) + geom_blank() +
        geom_point(dim_theta[setdiff(1:nrow(dim_theta), t_sup_sea), ], mapping = aes(x=dim,y=theta), col="black") +
        geom_point(dim_theta[t_sup_sea, ], mapping = aes(x=dim,y=theta), col="red") +
        theme_bw()
      print(p)
      dev.print(pdf, file=paste0(output.dir, city, "/Image/", output,".pdf"), width = 7, height = 7)
      message(paste0("plot for ", output, " of ", city, " is saved!"))
        
      data = apply(data_anomaly[,,t_sup_sea], 1:2, mean)
      colnames(data) = LAT; rownames(data) = LON
      data = data[,rev(seq_along(LAT))]
      print(range(data))
      
      output = paste0(var, "_sup_1999-2017_", season[seas], "_", SD)
      dat = melt(data, varnames = c("long", "lat"))
      plot_worldmap(data = dat, val.limits = c(-100, 100))
      dev.print(pdf, file=paste0(output.dir, city, "/Image/", output,".pdf"), width = 11, height = 5)
      message(paste0("plot for ", output, " of ", city, " is saved!"))
    }
    
    ## INF
    {
      t_inf_sea = intersect(t_inf, idxdates)
      
      output = paste0(var, "_inf_1999-2017_", season[seas], "_", SD, "_dim_theta")
      p <- ggplot(dim_theta, aes(x = dim, y = theta)) + geom_blank() +
        geom_point(dim_theta[setdiff(1:nrow(dim_theta), t_inf_sea), ], mapping = aes(x=dim,y=theta), col="black") +
        geom_point(dim_theta[t_inf_sea, ], mapping = aes(x=dim,y=theta), col="red") +
        theme_bw()
      print(p)
      dev.print(pdf, file=paste0(output.dir, city, "/Image/", output,".pdf"), width = 7, height = 7)
      message(paste0("plot for ", output, " of ", city, " is saved!"))
      
      data = apply(data_anomaly[,,t_inf_sea], 1:2, mean)
      colnames(data) = LAT; rownames(data) = LON
      data = data[,rev(seq_along(LAT))]
      print(range(data))
      
      output = paste0(var, "_inf_1999-2017_", season[seas], "_", SD)
      dat = melt(data, varnames = c("long", "lat"))
      plot_worldmap(data = dat, val.limits = c(-100, 100))
      dev.print(pdf, file=paste0(output.dir, city, "/Image/", output,".pdf"), width = 11, height = 5)
      message(paste0("plot for ", output, " of ", city, " is saved!"))
    }
  }
}
####

