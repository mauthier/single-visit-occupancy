##--------------------------------------------------------------------------------------------------------
## SCRIPT : Single-Visit Site Occupancy modelling
##
## Authors : Matthieu Authier & Ghislain Vieilledent
## Last update : 2016-04-21
## R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
##--------------------------------------------------------------------------------------------------------

rm (list = ls ())
lapply(c("coda", "mvtnorm", "binom", "hSDM", "rgdal", "shapefiles", "maptools", "gdata", "maps", "sp", "raster", "fields", "dplyr", "ggplot2", "marmap"), library, character.only=TRUE)

WorkDir <- "D:/Dossier mathieu/Desktop/hSDM_deldel"
DataDir <- "D:/Dossier mathieu/Desktop/hSDM_deldel"

setwd(WorkDir)

### load a bunch of custom functions for this analysis
## these functions are all built on hSDM main functions
source (paste (DataDir, "single_visit_occupancy_fct2source.r", sep = "/"))

### prepare Landscape for hSDM
lat <- seq (43.00, 48.40, 0.20)
lon <- seq (-6.00, -0.60, 0.20)

### projection
my_proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
xLand <- length (lon) - 1
yLand <- length (lat) - 1
Landscape <- raster (ncol = xLand, nrow = yLand, crs = my_proj)

Landscape[] <- 0
extent (Landscape) <- c (range (lon), range (lat))
coords <- coordinates (Landscape)
ncells <- ncell (Landscape)

### Neighbours
neighbors.mat <- adjacent (Landscape, cells = c (1:ncells), 
                           directions = 8, pairs = TRUE, sorted = TRUE
                           )
n.neighbors <- as.data.frame (table (as.factor (neighbors.mat[, 1])))[, 2]
summary (n.neighbors)
adj <- neighbors.mat[, 2]

### Generate symmetric adjacency matrix, A
A <- matrix (0, ncells, ncells)
index.start <- 1
for (i in 1:ncells) {
  index.end <- index.start + n.neighbors[i] - 1
  A[i, adj[c (index.start:index.end)]] <- 1
  index.start <- index.end + 1
}

### load 2012-2014 data: dd means "delphinus delphis"
dd1214 <- read.table (paste (DataDir, "CommonDolphinBayBiscayPELGAS.txt", sep = "/"), 
                      header = TRUE, colClasses = "numeric"
                      )

### extract cells that were visited in 2012, 2013 and 2014 during the pelgas survey
pelgas_box <- summarise (group_by (dd1214, Site), 
                         sum (ifelse (is.na (Detected), 0, 1)))$Site[which (summarise (group_by (dd1214, Site), sum (ifelse (is.na (Detected), 0, 1)))[2] == 3)]

### prepare data
dd1214 <- arrange (dd1214, Site, Year)
dd1214_suit <- summarise (group_by(dd1214, Site), unique(Bathy), unique(gBathy), unique(Iso200), unique(Lon), unique(Lat))
names(dd1214_suit) <- c ("Site", "Bathy", "gBathy", "Iso200", "Lon", "Lat")
dd1214_suit$cells <- extract (Landscape,
                              with (dd1214_suit,SpatialPoints (coords = cbind (Lon, Lat))),
                              cell = TRUE)[, 1]

### fit site-occupancy models
n_burnin <- 10000 ; n_mcmc <- 50000 ; n_thin <- 10 ; n_sim <- (n_mcmc - n_burnin)/n_thin

### without covariates
mod.hSDM.siteocc1 <- hSDM_siteocc (dd = dd1214, dd_suit = dd1214_suit)

### with covariates
mod.hSDM.siteocc2 <- hSDM_siteocc (dd = dd1214, dd_suit = dd1214_suit, cov = TRUE)

### check parameters convergence
gelman.diag (list (mod.hSDM.siteocc1[[1]]$mcmc,
                   mod.hSDM.siteocc1[[2]]$mcmc,
                   mod.hSDM.siteocc1[[3]]$mcmc), 
             autoburnin = FALSE)

gelman.diag (list (mod.hSDM.siteocc2[[1]]$mcmc,
                   mod.hSDM.siteocc2[[2]]$mcmc,
                   mod.hSDM.siteocc2[[3]]$mcmc), 
             autoburnin = FALSE)

### posterior estimates of occurrence at an average site
get_ci (inv_logit (mod.hSDM.siteocc1[[1]]$mcmc[,1]))
### posterior estimates of detection at an average site
get_ci (inv_logit (mod.hSDM.siteocc1[[1]]$mcmc[,2]))

### posterior estimates of occurrence at an average site
get_ci (inv_logit (mod.hSDM.siteocc2[[1]]$mcmc[,1]))
### posterior estimates of detection at an average site
get_ci (inv_logit (mod.hSDM.siteocc2[[1]]$mcmc[,2]))

### WAIC
# model without covariates
compute_waic(model = mod.hSDM.siteocc1[[1]], 
             dd = dd1214, 
             dd_suit = dd1214_suit, 
             cov = FALSE
             )
# model with covariates
compute_waic(model = mod.hSDM.siteocc2[[1]], 
             dd = dd1214, 
             dd_suit = dd1214_suit, 
             cov = TRUE
             )

### analyze 2014 data with single-visit model
dd2014 <- subset (dd1214, Year == 2014)
dd2014 <- arrange (dd2014, Site, Year)
dd2014$cells <- extract (Landscape,
                         with (dd2014, SpatialPoints (coords = cbind (Lon, Lat))),
                         cell = TRUE)[,1]

dd2014_suit <- summarise (group_by(dd2014, Site), unique (Bathy), unique (gBathy), unique (Iso200), unique (Lon), unique (Lat))
names (dd2014_suit) <- c("Site", "Bathy", "gBathy", "Iso200", "Lon", "Lat")
dd2014_suit$cells <- extract (Landscape,
                              with (dd2014_suit, SpatialPoints (coords = cbind (Lon,Lat))),
                              cell = TRUE)[,1]
### fit models
# single visit occupancy model
mod.hSDM.singlevisit1 <- hSDM_siteocc (dd = dd2014, dd_suit = dd2014_suit, cov = TRUE)
# single visit binomial model
mod.hSDM.singlevisit2 <- hSDM_binom (dd = dd2014, cov = TRUE)

### check parameters convergence
gelman.diag (list (mod.hSDM.singlevisit1[[1]]$mcmc,
                   mod.hSDM.singlevisit1[[2]]$mcmc,
                   mod.hSDM.singlevisit1[[3]]$mcmc), 
             autoburnin = FALSE)

gelman.diag (list (mod.hSDM.singlevisit2[[1]]$mcmc,
                   mod.hSDM.singlevisit2[[2]]$mcmc,
                   mod.hSDM.singlevisit2[[3]]$mcmc), 
             autoburnin = FALSE)

### posterior estimates of occurrence at an average site
get_ci (inv_logit (mod.hSDM.singlevisit1[[1]]$mcmc[,1]))
### posterior estimates of detection at an average site
get_ci (inv_logit (mod.hSDM.singlevisit1[[1]]$mcmc[,2]))

### posterior estimates of occurrence at an average site
get_ci (inv_logit (mod.hSDM.singlevisit2[[1]]$mcmc[,1]))

### WAIC
# single visit occupancy model
compute_waic(model = mod.hSDM.singlevisit1[[1]], 
             dd = dd2014, dd_suit=dd2014_suit, cov = TRUE)
# single visit binomial model
compute_waic(model = mod.hSDM.singlevisit2[[1]], 
             dd = dd2014, dd_suit=dd2014_suit, cov = TRUE)

### proportion of occupied area (POA) analysis
## use Agresti & Coull (1998)'s rule: +2 successes, +2 failures
# 2012-2014, site occupancy model w/o covariates
get_ci(
  (apply (psi_hat (model = mod.hSDM.siteocc1[[1]], 
                   dd = dd1214, 
                   dd_suit = dd1214_suit, cov = FALSE)[, which(dd1214_suit$Site %in% pelgas_box)], 1, sum) + 2)/(length(pelgas_box) + 4)
)

# 2012-2014, site occupancy model with covariates
get_ci(
  (apply (psi_hat (model = mod.hSDM.siteocc2[[1]], 
                   dd = dd1214, 
                   dd_suit = dd1214_suit, cov = TRUE)[, which(dd1214_suit$Site %in% pelgas_box)], 1, sum) + 2)/(length(pelgas_box) + 4)
)

# 2014, site occupancy model with covariates
get_ci(
  (apply (psi_hat (model = mod.hSDM.singlevisit1[[1]], 
                   dd = dd2014, 
                   dd_suit = dd2014_suit, cov = TRUE)[, which(dd2014_suit$Site %in% pelgas_box)], 1, sum) + 2)/(length(pelgas_box) + 4)
)

# 2014, binomial model with covariates
get_ci(
  (apply (psi_hat (model = mod.hSDM.singlevisit2[[1]], 
                   dd = dd2014, 
                   dd_suit = dd2014_suit, cov = TRUE)[, which(dd2014_suit$Site %in% pelgas_box)], 1, sum) + 2)/(length(pelgas_box) + 4)
)

save.image("hSDM_github.RData", safe = TRUE)

### simple plots
# dd4$psi <- apply(psi_hat(model=mod2[[1]], dd=dd4, dd_suit=dd4_suit, cov=TRUE),2,mean)
# dd4 <- dd4[, c("Lon", "Lat", "psi")]
# coordinates(dd4) <- ~Lon+Lat
# # coerce to SpatialPixelsDataFrame
# gridded(dd4) <- TRUE
# # coerce to raster
# dd4 <- raster(dd4)
# projection(dd4) <-  my_proj
# plot(dd4)