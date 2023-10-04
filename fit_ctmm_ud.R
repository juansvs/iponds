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
length(availableWorkers())
cl <- makeCluster(length(availableWorkers()))



#registerDoParallel(cl)

# read in data
telemetries <- readRDS("data/telemetries_w1_day.rds")
GUESS <- readRDS("data/GUESS_w1_day.rds")
SP <- raster::shapefile("data/polygon_pond.shp")
#SP <- SpatialPolygons(list(Polygons(list(Polygon(sf::st_coordinates(pond$geometry)[,1:2])), 1)), proj4string = CRS("+proj=tmerc"))
# Define functions to run in parallel, with intermediate output
parctmm<-function(x) {
    out<-ctmm::ctmm.select(telemetries[[x]],GUESS[[x]])
    saveRDS(out,paste0("outputs/FIT_",x,".rds"))
    out
}

parakde <- function(x) {
    out<-ctmm::akde(telemetries[[x]],CTMMs[[x]],SP=SP, grid = list(dr = c(0.5,0.5), align.to.origin = T))
    saveRDS(out,paste0("outputs/UD_",x,".rds"))
    out
}

# Send the functions and variables to all workers
clusterExport(cl = cl, c("telemetries", "GUESS", "SP", "parctmm", "parakde"))

# Run in parallel
CTMMs <- clusterApplyLB(cl = cl, seq_along(telemetries), parctmm)
saveRDS(CTMMs, "w1_day_FITS.rds")
UD <- clusterApplyLB(cl = cl, seq_along(telemetries), parakde)
saveRDS(UD, "w1_day_UDS.rds")
# create CTMM objects in parallel
#CTMMs <- foreach(i = seq_along(telemetries), .packages = "ctmm") %dopar% {
  # Fit the CTMM models. 1 took about 15 minutes on a single core.
#  ctmm.select(telemetries[[i]], GUESS[[i]])
#}

# Estimate the UDs

#UD <- akde(telemetries, CTMMs, SP = SP, grid = list(dr = c(0.5,0.5), align.to.origin=TRUE))

# Export the CTMM and UDs

