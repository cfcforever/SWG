save_pdf <- function(filename = "output", width = 7, height = 7){
  dev.print(pdf, file = paste0(filename, ".pdf"), width = width, height = height)
}
