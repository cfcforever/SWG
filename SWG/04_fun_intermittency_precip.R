fun_intermittency = function(x, thres = 0){
  y = rep(NaN, length(x))
  for (i in 1:length(x)){
    if (x[i]>0){
      y[i] = 1
    }else{
      y[i] = 0
    }
  }
  return (y)
}

fun_prob_0_1 <- function(x){
  y1 = x[1:(length(x)-1)]
  y2 = x[2:length(x)]
  y = cbind(y1,y2)
  
  z = as.data.frame(matrix(0, 2, 2))
  colnames(z) = c("0", "1")
  rownames(z) = c("0", "1")
  for (k in 1:(length(x)-1)){
    if (y[k,1]==0 & y[k,2]==0){
      z["0","0"] = z["0","0"] + 1
    }else if(y[k,1]==0 & y[k,2]==1){
      z["0","1"] = z["0","1"] + 1
    }else if(y[k,1]==1 & y[k,2]==0){
      z["1","0"] = z["1","0"] + 1
    }else if(y[k,1]==1 & y[k,2]==1){
      z["1","1"] = z["1","1"] + 1
    }else{
      k = k + 1
    }
  }
  return(z/(length(x)-1))
}

#
time = as.Date(1:4748, origin = "2004-12-31")
range(time)

load("~/Documents/LSCE/SWG/precip_47_51_N_10_15_E_1979_2017.RData")
Obs = precip[which(rownames(precip)=="2005-01-01"):which(rownames(precip)=="2017-12-31"),]
plot(x = time, y = Obs)
range(Obs)   # 0.0000 36.0525

Obs_occ = fun_intermittency(Obs)
Obs_0_1 = fun_prob_0_1(Obs_occ)

load("~/Documents/LSCE/SWG/tp/case_1/SIMU/SIMU_SWG_precip_run_026.RData")   # 19,26,88
Sample = Sample[,1]
plot(x = time, y = Sample)
range(Sample)


## esay version 
tab_0_1 = matrix(as.matrix(Obs_0_1), 1, 4)
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  
  load(paste0("~/Documents/LSCE/SWG/tp/case_", 1, "/SIMU/SIMU_SWG_precip_run_", NUM, ".RData"))
  Sample_occ = fun_intermittency(Sample)
  Sample_0_1 = fun_prob_0_1(Sample_occ)
  
  tab_0_1 = rbind(tab_0_1, matrix(as.matrix(Sample_0_1), 1, 4))
} 

colnames(tab_0_1) = c("00", "01", "10", "11")
tab_0_1 = as.data.frame(tab_0_1)
melt_0_1 = melt(tab_0_1)
melt_0_1$type = rep(c("obs", rep("sim", 100)), 4) 


p = ggplot(melt_0_1)
p = p + geom_boxplot(data = melt_0_1[melt_0_1$type=="sim",], aes(x = variable, y = value))
p = p + geom_point(data = melt_0_1[melt_0_1$type=="obs",], aes(x = variable, y = value), col = "red")
p = p + theme_bw()
p = p + xlab("intermittency transition") + ylab("probability")
p = p + labs(title = "case_3")
p
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/tp/Images/Intermittency_transition_case_3.pdf"), width = 7, height = 7)


## better version 
tab_0_1 = matrix(as.matrix(Obs_0_1), 1, 4)
for (k in 1:3){
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    
    load(paste0("~/Documents/LSCE/SWG/tp/case_", k, "/SIMU/SIMU_SWG_precip_run_", NUM, ".RData"))
    Sample_occ = fun_intermittency(Sample)
    Sample_0_1 = fun_prob_0_1(Sample_occ)
    
    tab_0_1 = rbind(tab_0_1, matrix(as.matrix(Sample_0_1), 1, 4))
  } 
}

colnames(tab_0_1) = c("00", "01", "10", "11")
tab_0_1 = as.data.frame(tab_0_1)
melt_0_1 = melt(tab_0_1)
melt_0_1$type = rep(c("obs", rep("sim", 300)), 4) 
melt_0_1$case = rep(c("obs", rep("PCS", 100), rep("SysDyn", 100), rep("PCS+SysDyn", 100)), 4)

p = ggplot(melt_0_1)
p = p + geom_boxplot(data = melt_0_1[melt_0_1$type=="sim",], aes(x = variable, y = value, col = case))
p = p + geom_point(data = melt_0_1[melt_0_1$type=="obs",], aes(x = variable, y = value), col = "black")
p = p + theme_bw()
p = p + xlab("intermittency transition") + ylab("probability")
p
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/tp/Images/Intermittency_transition.pdf"), width = 7, height = 7)


