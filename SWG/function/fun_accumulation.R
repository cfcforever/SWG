## load data
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
load("~/Documents/LSCE/SWG/precip_47_51_N_10_15_E_1979_2017.RData")

## create DATE calendar
mDate = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))
DATE = as.Date((0:(nrow(mDate)-1)), origin = "1979-01-01")

## define number of accumulation days
naccu = 4

## fun_accumulation
fun_accumulation <- function(x, naccu, period, fun){
  
  # x: 1d-sequence
  # naccu: number of accumulation
  # period: 1d-sequence of date with format "yyyy-mm-dd"
  # fun: specific function to accumulate several time steps
  
  library(lubridate)
  y = year(period)
  ny = length(unique(y))
  
  ly = vector("list", ny); names(ly) = unique(y)
  ty = vector("list", ny); names(ty) = unique(y)
  dy = vector("list", ny); names(dy) = unique(y)
  for (i in unique(y)){
    xi = x[y==i]
    ly[[paste(i)]] = xi
    l = length(xi)
    istep = l-naccu+1
    
    m = matrix(NA, nrow = naccu, ncol = istep)
    for (k in 1:naccu){
      m[k,] = xi[c(k:(k+istep-1))]
    }
    
    ty[[paste(i)]] = apply(m, 2, FUN = fun)
    
    yi = period[y==i]
    dat = as.data.frame(matrix(NA, nrow = istep, ncol = 2))
    colnames(dat) = c("date", "value")
    dat[,1] = yi[c(1:istep)]
    dat[,2] = ty[[paste(i)]]
    dy[[paste(i)]] = dat

  }
  
  xaccu = do.call(rbind.data.frame, dy)
  rownames(xaccu) = NULL
  
  return(xaccu)
}


plot_annual_max = function(x, var){
  y = year(x[,1]); ny = length(unique(y))
  y_max = as.data.frame(matrix(NA, nrow = ny, ncol = 2))
  colnames(y_max) = c("year", "value")
  for (k in 1:ny){
    y_max[k,1] = unique(y)[k]
    y_max[k,2] = max(x[y==unique(y)[k],2]) 
  }
  
  library(ggplot2)
  p = ggplot(y_max, aes(x = year, y = value)) + geom_point() + theme_bw()
  p = p + labs(y=var)
  p
}

naccu = 4
tm = fun_accumulation(x = tmean[,1], naccu, period = DATE, fun = mean)
plot_annual_max(tm, var = paste0("annual max of ", naccu, "-day Temperature Mean"))
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/t2m_", naccu, "_day.pdf"))

tp = fun_accumulation(precip[,1], naccu, period = DATE, fun = sum)
plot_annual_max(tp, var = paste0("annual max of ", naccu, "-day Total Precipitation"))
dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/tp_", naccu, "_day.pdf"))

for(naccu in 4){
  tm = fun_accumulation(x = tmean[,1], naccu, period = DATE, fun = mean)
  plot_annual_max(tm, var = paste0("annual max of ", naccu, "-day Temperature Mean"))
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/t2m_", naccu, "_day.pdf"), width=7, height=7)
  
  tp = fun_accumulation(precip[,1], naccu, period = DATE, fun = sum)
  plot_annual_max(tp, var = paste0("annual max of ", naccu, "-day Total Precipitation"))
  dev.print(pdf, file=paste0("~/Documents/LSCE/SWG/tp_", naccu, "_day.pdf"), width=7, height=7)
}
