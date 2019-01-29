#### Conditional SWG for temperature in different cities ---------------------
city.names = c("Paris", "Madrid", "Stockholm")
for (city in city.names){
  source("conditional_SWG_temperature.R") 
  rm(list=ls(all=TRUE))
}
####



#### FAR ---------------------------------------------------------------------

####



#### Changes of intensity ----------------------------------------------------
city.names = c("Paris", "Madrid", "Stockholm")
for (city in city.names){
  source("changes_of_intensity.R") 
  rm(list=ls(all=TRUE))
}
####



#### Anomaly maps ------------------------------------------------------------
city.names = c("Paris", "Madrid", "Stockholm")
for (city in city.names){
  source("anomaly_maps.R") 
  rm(list=ls(all=TRUE))
}
####
