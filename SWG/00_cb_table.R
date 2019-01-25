cb_table = matrix(NA, nrow = 63, ncol = 6)
n = 1
for (k in 1:6){
  temp = combn(1:6, k)
  m = nrow(temp)
  for (i in 1:ncol(temp)){
    cb_table[n, 1:m] = temp[,i]
    n = n + 1
  }
}

cb_table = as.data.frame(cb_table)
save(cb_table, file = "~/cb_table.RData")

cb_table[cb_table==1] = "slp_PC"
cb_table[cb_table==2] = "slp_SysDyn"
cb_table[cb_table==3] = "z500_PC"
cb_table[cb_table==4] = "z500_SysDyn"
cb_table[cb_table==5] = "sst_PC"
cb_table[cb_table==6] = "sst_SysDyn"
save(cb_table, file = "~/cb_table_2.RData")

case_name = rep(NA, nrow(cb_table))
for (k in 1:nrow(cb_table)){
  case_name[k] = paste(cb_table[k,!is.na(cb_table[k,])], collapse = "+")
}
save(case_name, file = "~/Documents/LSCE/SWG/NA/case_name.RData")

library(xtable)
xtable(cb_table)

library(gridExtra)
pdf("~/cb_table.pdf", height = 20, width = 8.5)
grid.table(cb_table)
dev.off()



