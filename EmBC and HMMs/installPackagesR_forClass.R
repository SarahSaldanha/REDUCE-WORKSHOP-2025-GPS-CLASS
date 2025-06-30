rm(list = ls())

Sys.setenv(TZ='UTC')

options(repos = c(CRAN = "http://cran.fiocruz.br/"))

pkg.list = c("momentuHMM", "lubridate", "ggplot2", "marmap", "viridis", "raster", "ncdf4", "reshape2", "circular", "gridExtra", "viridisLite")
install.packages(pkg.list, dependencies = TRUE, repos = "http://cran.fiocruz.br/")



