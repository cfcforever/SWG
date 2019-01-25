plot_worldmap <- function(data, val.limits = c(-200,200), col.contour = "white",
                          brks = seq(from = -180, to = 180, by = 30),
                          label = seq(from = -180, to = 180, by = 30),
                          brk_x = (3:7)*10, brk_y = (-4:2)*20, fill.name = "Pa"){
  #  data
  #= a data.frame with 3 columns which are c("LON", "LAT", "value").
  #
  #  var.limits
  #= the limits for "value" showing in the plot, out of limits will squish with min and max.
  
  library(ggplot2)
  library(scales)
  library(directlabels)
  
  colnames(data) = c("long", "lat", "value")
  
  world <- map_data("world")
  worldmap <- ggplot() + theme_bw() +
    geom_raster(data = data, mapping = aes(x = long, y = lat, fill = value), interpolate = TRUE) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                         limits = val.limits, oob = squish) +
    geom_contour(data = data, mapping = aes(x = long, y = lat, z = value),
                 col = col.contour, breaks = brks) +
    geom_dl(data = data, aes(x = long, y = lat, z = value, label=..level..),
            method = "bottom.pieces", stat = "contour", col = col.contour, breaks = brks) +
    geom_path(world, mapping = aes(x = long, y = lat, group = group)) +
    coord_cartesian(xlim = range(data$long)+c(-0.5,0.5), ylim = range(data$lat)+c(-0.5,0.5)) + 
    scale_y_continuous(breaks = brk_x, expand = c(0,0)) +
    scale_x_continuous(breaks = brk_y, expand = c(0,0)) + 
    labs(x = "Lon", y = "Lat", fill = fill.name)
  print(worldmap)
}

# save(data, LON, LAT, file = "/Volumes/Data-ExFAT/data_for_plot_worldmap.RData")
# load(file = "/Volumes/Data-ExFAT/data_for_plot_worldmap.RData")
# dat = melt(data, varnames = c("long", "lat"))
# plot_worldmap(data = dat, val.limits = c(-100,100))


