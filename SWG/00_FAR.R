#### load tmean OBS ####
load("~/Documents/LSCE/SWG/tmean_48_50_N_01_04_E_1979_2017.RData")
DATE_OBS = atoms(timeSequence(from="1979-01-01",to="2017-12-31",by='day'))

Obs = tmean[DATE_OBS$Y>=1999 & DATE_OBS$Y<=2017,]
n = length(Obs)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)

DATE = as.Date((0:(n-1)), origin = "1999-01-01")


#### FAR ####
load(file = "~/Documents/LSCE/SWG/NA/case_name.RData")
case_name = c(case_name, "statio")

# case_sel = c(1:64)
# case_sel = c(8,9,29,48,57,63,64)
case_sel = c(63, 61, 57, 48, 59, 50, 43, 27, 64)   # best BIC
case_sel = c(1, 2, 4, 7, 9, 13, 23, 64)            # best intensity
case_name = c(case_name[case_sel])
ncase = length(case_name)


list_mean_sd_nat = vector("list", ncase)
list_mean_sd_rea = vector("list", ncase)
for (k in 1:ncase){
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/nat/t2m_mean_sd_1999_2017_nat_", case_sel[k], ".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_nat[[k]] = m
  
  cat(k,'\n')
  load(paste0("~/Documents/LSCE/SWG/NA_1.5/tmean/rea/t2m_mean_sd_1999_2017_rea_", case_sel[k], ".RData"))
  cat(range(mean),'\n')
  cat(range(sd),'\n')
  m = as.data.frame(matrix(NA, nrow = nrow(mean), ncol = 2))
  colnames(m) = c("mean", "sd")
  m[,"mean"] = mean
  m[,"sd"]   = sd
  list_mean_sd_rea[[k]] = m
}

##
for (year in 2003:2003){
  for (thres in 10:30){
    library(lubridate)
    
    list_prob_thres_nat = lapply(list_mean_sd_nat, FUN = function(m){
      pnorm(thres, mean = m[,"mean"], sd = m[,"sd"], lower.tail = F)
    })
    prob_thres_nat = do.call(rbind.data.frame, list_prob_thres_nat)
    colnames(prob_thres_nat) = "nat"
    
    list_prob_thres_rea = lapply(list_mean_sd_rea, FUN = function(m){
      pnorm(thres, mean = m[,"mean"], sd = m[,"sd"], lower.tail = F)
    })
    prob_thres_rea = do.call(rbind.data.frame, list_prob_thres_rea)
    colnames(prob_thres_rea) = "rea"
    
    prob_thres = cbind(prob_thres_nat, prob_thres_rea)
    prob_thres$date = rep(DATE, ncase)
    
    # prob_thres$FAR = 1 - prob_thres_nat[,1]/prob_thres_rea[,1]
    # output = paste0("~/Documents/LSCE/SWG/NA_1.5/Image/FAR/average_FAR_tmean_over_", thres, "_degree.pdf")
    
    # prob_thres$FAR = prob_thres$nat/prob_thres$rea
    # output = paste0("~/Documents/LSCE/SWG/NA_1.5/Image/FAR/average_RR_tmean_over_", thres, "_degree.pdf")
    
    prob_thres$FAR = (prob_thres$nat-prob_thres$rea)/(prob_thres$nat+prob_thres$rea)
    output = paste0("~/Documents/LSCE/SWG/NA_1.5/Image/FAR/average_FARD_tmean_over_", thres, "_degree_2.pdf")
    
    prob_thres$case = rep(case_name, each = length(DATE))
    
    # ## add case_0 stationary
    # load(paste0("~/Documents/LSCE/SWG/t2m/case_0/Sample_1998_2017_nat.RData"))
    # nat = apply(Sample.nat, MARGIN = 1, FUN = function(x){sum(x>=thres)/length(x)})
    # load(paste0("~/Documents/LSCE/SWG/t2m/case_0/Sample_1998_2017_rea.RData"))
    # rea = apply(Sample.rea, MARGIN = 1, FUN = function(x){sum(x>=thres)/length(x)})
    
    # statio = as.data.frame(cbind(nat, rea))
    # statio$date = DATE
    # statio$FAR = nat/rea
    # statio$case = "statio"
    
    # prob_thres = rbind.data.frame(prob_thres, statio)
    
    period_sel = intersect(which(year(prob_thres$date)==year) , which(month(prob_thres$date)==7|month(prob_thres$date)==6|month(prob_thres$date)==8))
    FAR_sel = prob_thres[period_sel,]
    FAR_sel$num = rep(1:(length(period_sel)/ncase), ncase)
    
    FAR_sel$case <- factor(FAR_sel$case, levels = case_name)
    # FAR_sel$num = rep(1:(length(period_sel)/4), 4)
    
    # breaks = rep(NA,5)
    # labels = rep(NA,5)
    # for (k in 1:10){
    #   breaks[k] = which(year(FAR_sel$date)==unique(year(FAR_sel$date))[4*k-1])[1]
    #   labels[k] = as.character(FAR_sel$date[breaks[k]])
    # }
    
    breaks = rep(NA,3)
    labels = rep(NA,3)
    for (k in 1:3){
      breaks[k] = which(month(FAR_sel$date)==unique(month(FAR_sel$date))[k])[1]
      labels[k] = as.character(FAR_sel$date[breaks[k]])
    }
    
    p = ggplot(FAR_sel) 
    p = p + geom_line(aes(x=num, y=FAR, group=case, col=case)) + theme_bw()
    p = p + scale_x_continuous(breaks = breaks, labels = labels) + ylim(-1,1)
    p = p + labs(title = paste0("FAR_t2M_summer_over_", thres, "_degree"))
    print(p)
    dev.print(pdf, file=output, width = 14, height=7)
  }
}


## average
for (thres in 10:30){
  # thres = 10
  library(lubridate)
  
  list_prob_thres_nat = lapply(list_mean_sd_nat, FUN = function(m){
    pnorm(thres, mean = m[,"mean"], sd = m[,"sd"], lower.tail = F)
  })
  prob_thres_nat = do.call(rbind.data.frame, list_prob_thres_nat)
  colnames(prob_thres_nat) = "nat"
  
  list_prob_thres_rea = lapply(list_mean_sd_rea, FUN = function(m){
    pnorm(thres, mean = m[,"mean"], sd = m[,"sd"], lower.tail = F)
  })
  prob_thres_rea = do.call(rbind.data.frame, list_prob_thres_rea)
  colnames(prob_thres_rea) = "rea"
  
  prob_thres = cbind(prob_thres_nat, prob_thres_rea)
  prob_thres$date = rep(DATE, ncase)
  
  # prob_thres$FAR = 1 - prob_thres_nat[,1]/prob_thres_rea[,1]
  # output = paste0("~/Documents/LSCE/SWG/NA_1.5/Image/FAR/FAR_tmean_summer_", year, "_over_", thres, "_degree.pdf")
  
  # prob_thres$FAR = prob_thres$nat/prob_thres$rea
  # output = paste0("~/Documents/LSCE/SWG/NA_1.5/Image/FAR/FAR_tmean_summer_", year, "_over_", thres, "_degree.pdf")
  
  prob_thres$FAR = (prob_thres$nat-prob_thres$rea)/(prob_thres$nat+prob_thres$rea)
  output = paste0("~/Documents/LSCE/SWG/NA_1.5/Image/FAR/FARD_tmean_summer_", year, "_over_", thres, "_degree.pdf")
  
  prob_thres$case = rep(case_name, each = length(DATE))
  
  FAR = as.data.frame(matrix(NA, nrow = 64, ncol = 2))
  colnames(FAR) = c("value", "model")
  for (k in 1:64){
    FAR[k,] = c(mean(prob_thres$FAR[prob_thres$case==case_name[k]]), case_name[k])
  }
  FAR$number = c(1:64)
  FAR$value <- factor(FAR$value, levels = order(FAR$value))
  
  p = ggplot(FAR) 
  p = p + geom_text(aes(x=number, y=value, label = number)) +  theme_bw()
  p = p + geom_point(aes(x=number, y=value), col = "red", size = 0.5)
  p = p + labs(title = paste0("FAR_t2M_summer_over_", thres, "_degree"))
  print(p)
  dev.print(pdf, file=output, width = 14, height=7)
}

