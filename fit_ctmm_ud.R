# Code to analyze fish tracking data, to generate CTMM objects and UD using
# autocorrelated kernel density estimation. The telemetries objects and initial
# guesses for CTMM have been generated already.

# load libraries
library(foreach)
library(parallel)
library(doParallel)
library(future)
library(ctmm)
library(sp)

# Prepare parallel workers
cl <- makeCluster(length(availableWorkers()))
registerDoParallel(cl)

# read in data
telemetries <- readRDS("data/telemetries_w1_day.rds")
GUESS <- readRDS("data/GUESS_w1_day.rds")
pond <- sf::read_sf("data/polygon_pond.shp")
SP <- SpatialPolygons(list(Polygons(list(Polygon(sf::st_coordinates(pond$geometry)[,1:2])), 1)), proj4string = CRS("+proj=tmerc"))

# create CTMM objects in parallel
CTMMs <- foreach(i = seq_along(telemetries), .packages = "ctmm") %dopar% {
  # Fit the CTMM models. 1 took about 15 minutes on a single core.
  ctmm.select(telemetries[[i]], GUESS[[i]])
}

# Estimate the UDs
UD <- akde(telemetries, CTMMs, SP = SP, grid = list(dr = c(0.5,0.5)))

# Export the CTMM and UDs
saveRDS(UD, "w1_day_UDS.rds")
saveRDS(CTMMs, "w1_day_FITS.rds")

