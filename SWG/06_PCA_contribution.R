##-----------------------------------------------
# change t2m   to tp
# change tmean to precip 
##-----------------------------------------------

#### tmean ####
library(ggplot2)
library(reshape2)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')
casename = c("PCS", "SysDyn", "PCS_SysDyn")

for (case in 1:3){
  for (seas in 1:4){
    
    load(paste0("~/Documents/LSCE/SWG/t2m/case_", case, "/2PC/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val.RData"))
    
    par_tt = coef(fit_stations_tt[[1]], matrix = T)
    colnames(par_tt) = c("mean", "sd")
    contrib = as.data.frame(apply(abs(par_tt[-1,]), 2, FUN = function(x){x/sum(x)}))
    contrib$PC = rownames(contrib)
    contrib = melt(contrib, id.vars = "PC")
    contrib$PC <- factor(contrib$PC, levels = unique(contrib$PC))
    
    p = ggplot(contrib)
    p = p + geom_point(aes(x = PC, y = value)) 
    p = p + facet_wrap(~variable, nrow = 2) + theme_bw()
    p = p + labs(title = paste0("Contribution of principl components for t2m in ", season[seas], " (", casename[case], ")"), 
                 x = "Composant",
                 y = "Percent (X100%)")
    print(p)
    
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/t2m/Images/Contribution_composants_2PC_t2m_", season[seas], "_", casename[case], ".pdf"), width = 9, height = 7)
  }
}


#### precip ####
library(ggplot2)
library(reshape2)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')
casename = c("PCS", "SysDyn", "PCS_SysDyn")

for (case in 3:3){
  for (seas in 1:4){
    
    load(paste0("~/Documents/LSCE/SWG/tp/case_", case, "/SWG_ERAI_ESD_precip_", season[seas], "_cross-val.RData"))
    
    ## rf with operation
    par = coef(fit_stations_rf[[1]], matrix = T)
    par = exp(par)
    colnames(par) = c("rate", "shape")
    contrib = as.data.frame(apply(abs(par[-1,]), 2, FUN = function(x){x/sum(x)}))
    contrib$PC = rownames(contrib)
    contrib = melt(contrib, id.vars = "PC")
    contrib$PC <- factor(contrib$PC, levels = unique(contrib$PC))
    
    p = ggplot(contrib)
    p = p + geom_point(aes(x = PC, y = value)) 
    p = p + facet_wrap(~variable, nrow = 2) + theme_bw()
    p = p + labs(title = paste0("Contribution of principl components for precip in ", season[seas], " (", casename[case], ")"), 
                 x = "Composant",
                 y = "Percent (X100%)")
    print(p)
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/tp/Images/Contribution_composants_rf_", season[seas], "_", casename[case], ".pdf"), width = 9, height = 7)
    
    
    ## rf with NON operation
    par = coef(fit_stations_rf[[1]], matrix = T)
    colnames(par) = c("rate", "shape")
    contrib = as.data.frame(apply(abs(par[-1,]), 2, FUN = function(x){x/sum(x)}))
    contrib$PC = rownames(contrib)
    contrib = melt(contrib, id.vars = "PC")
    contrib$PC <- factor(contrib$PC, levels = unique(contrib$PC))
    
    p = ggplot(contrib)
    p = p + geom_point(aes(x = PC, y = value)) 
    p = p + facet_wrap(~variable, nrow = 2) + theme_bw()
    p = p + labs(title = paste0("Contribution of principl components for precip in ", season[seas], " (", casename[case], ")"), 
                 x = "Composant",
                 y = "Percent (X100%)")
    print(p)
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/tp/Images/Contribution_composants_rf_0_", season[seas], "_", casename[case], ".pdf"), width = 9, height = 7)
    
    
    ## ro with operation
    par = coef(fit_stations_ro[[1]], matrix = T)
    par = (exp(par))/(1+exp(par))
    contrib = as.data.frame(abs(par[-1,])/sum(abs(par[-1,])))
    colnames(contrib) = c("prob")
    contrib$PC = rownames(contrib)
    contrib$PC <- factor(contrib$PC, levels = unique(contrib$PC))
    
    p = ggplot(contrib)
    p = p + geom_point(aes(x = PC, y = prob)) + theme_bw()
    p = p + labs(title = paste0("Contribution of principl components for precip in ", season[seas], " (", casename[case], ")"), 
                 x = "Composant",
                 y = "Percent (X100%)")
    print(p)
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/tp/Images/Contribution_composants_ro_", season[seas], "_", casename[case], ".pdf"), width = 9, height = 7)
    
    
    ## ro with NON operation
    par = coef(fit_stations_ro[[1]], matrix = T)
    contrib = as.data.frame(abs(par[-1,])/sum(abs(par[-1,])))
    colnames(contrib) = c("prob")
    contrib$PC = rownames(contrib)
    contrib$PC <- factor(contrib$PC, levels = unique(contrib$PC))
    
    p = ggplot(contrib)
    p = p + geom_point(aes(x = PC, y = prob)) + theme_bw()
    p = p + labs(title = paste0("Contribution of principl components for precip in ", season[seas], " (", casename[case], ")"), 
                 x = "Composant",
                 y = "Percent (X100%)")
    print(p)
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/tp/Images/Contribution_composants_ro_0_", season[seas], "_", casename[case], ".pdf"), width = 9, height = 7)
    
  }
}

#### draft ####
seas = 1
case = 1
load(paste0("~/Documents/LSCE/SWG/t2m/case_", case, "/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val.RData"))

slotNames(fit_stations_tt[[1]])
# [1] "Bspline"          "nl.chisq"         "nl.df"            "spar"             "s.xargument"     
# [6] "var"              "extra"            "family"           "iter"             "predictors"      
# [11] "assign"           "callXm2"          "contrasts"        "df.residual"      "df.total"        
# [16] "dispersion"       "effects"          "offset"           "qr"               "R"               
# [21] "rank"             "ResSS"            "smart.prediction" "terms"            "Xm2"             
# [26] "Ym2"              "xlevels"          "call"             "coefficients"     "constraints"     
# [31] "control"          "criterion"        "fitted.values"    "misc"             "model"           
# [36] "na.action"        "post"             "preplot"          "prior.weights"    "residuals"       
# [41] "weights"          "x"                "y"   

par_tt = coef(fit_stations_tt[[1]], matrix = T)
class(par_tt)
contrib = as.data.frame(apply(abs(par_tt[-1,]), 2, FUN = function(x){x/sum(x)}))
contrib$PC = rownames(contrib)
contrib = melt(contrib, id.vars = "PC")
contrib$PC <- factor(contrib$PC, levels = unique(contrib$PC))

p = ggplot(contrib)
p = p + geom_point(aes(x = PC, y = value)) 
p = p + facet_wrap(~variable, nrow = 2) + theme_bw()
p = p + labs(title = paste0("Contribution of principl components for t2m in ", season[seas], " (", casename[case], ")"))
print(p)
