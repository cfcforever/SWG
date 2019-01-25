library(ggplot2)
library(reshape2)

MON=c(12,1:11)
MON=matrix(MON,nrow=4,ncol=3,byrow=T)
SEAS=c('DJF','MAM','JJA','SON')
season=c('Winter','Spring','Summer','Fall')

casename = c(1:63)
load(file = "~/Documents/LSCE/SWG/NA/case_name.RData")
casename = case_name
ncase = length(casename)

for (case in 52:ncase){
  for (seas in 1:4){
    
    load(paste0("~/Documents/LSCE/SWG/NA/tmean/1979_2004/SWG_ERAI_ESD_tmean_", season[seas], "_cross-val_1979_2004_", case, ".RData"))
    
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
    
    dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA/Image/Contribution/Contribution_t2m_1979_2004_", season[seas], "_", case, ".pdf"), width = 14, height = 7)
  }
}
