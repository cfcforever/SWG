## fun_qqplot
fun_qqplot <- function(Obs, Sample, main){
  qqplot(Obs, Sample, xlab = "obs", ylab = "sim", pch = 20,
         main = main, cex.main = 1.5, cex.lab = 1, 
         xlim = c(min(c(Obs, Sample)), max(c(Obs, Sample))),
         ylim = c(min(c(Obs, Sample)), max(c(Obs, Sample))))
  abline(0,1,col="red")
}


## load Obs data
load("~/Documents/LSCE/SWG/data/tmean_48_50_N_01_04_E_1979_2017.RData")
Obs = tmean[which(rownames(tmean)=="2005-01-01"):which(rownames(tmean)=="2017-12-31"),]

## validation period: 1999 - 2017
time = as.Date(1:length(Obs), origin = "2004-12-31")
range(time)

plot(x = time, y = Obs)
range(Obs) # -7.878333 27.243333

load(file = "~/Documents/LSCE/SWG/NA/case_name.RData")

## load one simulation (Sample)
k = 63
case = case_name[k]
load(paste0("~/Documents/LSCE/SWG/NA/tmean/SIMU/SIMU_SWG_tmean_2005_2017_nat_", k, "_run_001.RData"))
Sample = Sample[,1]
plot(x = time, y = Sample)
range(Sample)


## load all simulations
Sample_all = c()
for (i in 1:100){
  NUM=cbind(formatC(i, digits = 0, width = 3, format = "f", flag = "0"))
  load(paste0("~/Documents/LSCE/SWG/NA/tmean/SIMU/SIMU_SWG_tmean_2005_2017_nat_", k, "_run_", NUM, ".RData"))
  Sample_all = c(Sample_all, Sample)
}
range(Sample_all)

fun_qqplot(Obs, Sample_all, main = paste0("tmean - 100 simulations - ", case))
dev.print(pdf, file = paste0("~/Documents/LSCE/SWG/NA/Image/qqplot/qqplot_", case, ".pdf"), width = 9, height = 9)


## CramerVonMisesTwoSamples
CramerVonMisesTwoSamples <- function(S1, S2){
  xS1 = sort(S1)
  M = length(xS1)
  xS2 = sort(S2)
  N = length(xS2)
  a = data.frame(val = xS1, rang = seq(M), ens = rep(1, M))
  b = data.frame(val = xS2, rang = seq(N), ens = rep(2, N))
  c = rbind(a, b)
  c = c[order(c$val), ]
  c = data.frame(c, rangTot = seq(M + N))
  dtfM = c[which(c$ens == 1), ]
  dtfN = c[which(c$ens == 2), ]
  somN = sum((dtfN$rang - dtfN$rangTot)^2)
  somM = sum((dtfM$rang - dtfM$rangTot)^2)
  U = N * somN + M * somM
  CvM = ((U/(as.numeric(N) * as.numeric(M)))/(N + M)) - ((4 * as.numeric(M) * as.numeric(N) - 1)/(6 * (M + N)))
  return(CvM)
}

res <- CramerVonMisesTwoSamples(Sample_all, Obs)
pvalue = 1/6*exp(-res)
pvalue

# case_1:  0.1032542   (0.4788021)
# case_2:  0.006354681 (3.266804)
# case_3:  0.1005227   (0.5056121)


