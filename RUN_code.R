city.names = c("Paris", "Madrid", "Stockholm")
for (city in city.names[2:3]){
  source("conditional_SWG_temperature.R") 
  rm(list=ls(all=TRUE))
}
