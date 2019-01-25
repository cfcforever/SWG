source("~/Dropbox/Irstea/Code/Fonction/fun_nb_recurrence.R")
source("fun_intermittency_precip.R")

time = as.Date(1:4748, origin = "2004-12-31")
range(time)

load("~/Documents/LSCE/SWG/precip_47_51_N_10_15_E_1979_2017.RData")
Obs = precip[which(rownames(precip)=="2005-01-01"):which(rownames(precip)=="2017-12-31"),]
range(Obs)   # 0.0000 36.0525

Obs_occ = fun_intermittency(Obs)
Obs_duration = fun_nb_recurrence(Obs_occ)
Obs_dur_prob = apply(Obs_duration, MARGIN = 1, FUN = function(x){x/sum(x)})
colnames(Obs_dur_prob) = c("wet", "dry")
Obs_dur_prob = as.data.frame(Obs_dur_prob)
Obs_dur_prob$type = "obs"
Obs_dur_prob$case = "obs"

n=30
tab_dur_prob = Obs_dur_prob[1:n,]
tab_dur_prob$duration = c(1:n)
for (k in 1:3){
  for (i in 1:100){
    NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
    
    load(paste0("~/Documents/LSCE/SWG/tp/case_", k, "/SIMU/SIMU_SWG_precip_run_", NUM, ".RData"))
    Sample_occ = fun_intermittency(Sample)
    Sample_duration = fun_nb_recurrence(Sample_occ)
    Sample_dur_prob = apply(Sample_duration, MARGIN = 1, FUN = function(x){x/sum(x)})
    colnames(Sample_dur_prob) = c("wet", "dry")
    Sample_dur_prob = as.data.frame(Sample_dur_prob)
    Sample_dur_prob$type = paste0("sim_",i)
    Sample_dur_prob$case = paste0("case_",k)
    
    Sample_dur_prob = Sample_dur_prob[1:n,]
    Sample_dur_prob$duration = c(1:n)
    
    tab_dur_prob = rbind(tab_dur_prob, Sample_dur_prob)
  } 
}

p = ggplot(tab_dur_prob)
p = p + geom_boxplot(data = tab_dur_prob[tab_dur_prob$case=="case_3", ], aes(x = duration, y = dry, group = duration))
p = p + geom_line(data = tab_dur_prob[tab_dur_prob$type=="obs", ], aes(x = duration, y = dry, col = "red"))
p = p + theme_bw()
p = p + theme(legend.position="none")
p = p + labs(x="duration (day)", y = "probability (wet)", title = "case_3")
p
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/tp/Images/duration_distribution_dry_case_3.pdf"), width = 7, height = 7)
