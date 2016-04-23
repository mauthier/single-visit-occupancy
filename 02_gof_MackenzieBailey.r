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

### load fitted models
load (paste (WorkDir, "hSDM_deldel.RData", sep = "/"))

### Goodness of Fit test of MacKenzie & Bailey
dd1214_gof_p <- array (NA, dim = c (n_sim, nrow (dd1214_suit), 3, 2))
dd1214_p_cov <- array (1, dim = c (nrow (dd1214_suit), 3, 4))
dd1214_gof_psi <- array (NA, dim = c (n_sim, nrow (dd1214_suit), 2))

## design matrix for psi
X_occ <- matrix (NA, nrow = nrow (dd1214_suit), ncol = 4)
X_occ[, 1] <- 1
X_occ[, 2] <- rescale (dd1214_suit$Bathy)
X_occ[, 3] <- rescale (dd1214_suit$Lat)
X_occ[, 4] <- rescale (dd1214_suit$Iso200)

for( i in 1:nrow (dd1214_suit) ){
  for( j in 1:3) {
    dd1214_p_cov[i, j, 1] <- ifelse (nrow(subset(subset(dd1214, Site == dd1214_suit$Site[i]), Year == 2011 + j)) == 1, subset (subset (dd1214, Site == dd1214_suit$Site[i]), Year == 2011 + j)$Detected, NA)
    dd1214_p_cov[i, j, 3] <- ifelse (nrow(subset(subset(dd1214, Site == dd1214_suit$Site[i]), Year == 2011 + j)) == 1, (subset (subset (dd1214, Site == dd1214_suit$Site[i]), Year == 2011 + j)$Effort_km - mean (dd1214$Effort_km)) / (2 * sd (dd1214$Effort_km)), NA)
    dd1214_p_cov[i, j, 4] <- ifelse (nrow(subset(subset(dd1214, Site == dd1214_suit$Site[i]), Year == 2011 + j)) == 1, subset (subset (dd1214, Site == dd1214_suit$Site[i]), Year == 2011 + j)$Subj/2, NA)
    
    dd1214_gof_p[, i, j, 1] <- inv_logit (mod.hSDM.siteocc2[[1]]$mcmc[- c (1:(n_burnin/n_thin)), 5:7] %*% dd1214_p_cov[i, j, 2:4])
    dd1214_gof_p[, i, j, 2] <- rbinom (n_sim, size = 1, prob = dd1214_gof_p[, i, j, 1])
  }
  dd1214_gof_psi[, i, 1] <- inv_logit (mod.hSDM.siteocc2[[1]]$mcmc[- c (1:(n_burnin/n_thin)), 1:4] %*% X_occ[i, ])
  dd1214_gof_psi[, i, 2] <- rbinom (n_sim, size = 1, prob = dd1214_gof_psi[, i, 1])
}

### posterior predictive check
h <- apply (dd1214_p_cov[, , 1], 1, function (x) { paste (x[1], x[2], x[3], sep = "") })
history <- as.matrix (expand.grid (c (0, 1, NA), c (0, 1, NA), c (0, 1, NA)))[-27, ] # remove the sequence "NA-NA-NA" which is not observed
possible_h <- apply(history, 1, function (x){ paste (x[1], x[2], x[3], sep = "") })

## store results in an array
ppc <- array (NA, dim = c (length (possible_h), n_sim, 3)) # 1: obs, 2: rep, 3: expected

for( k in 1:dim (ppc)[1]) {
  if( any (names(table (h)) == possible_h[k]) ) {
    ppc[k, , 1] <- as.numeric (table (h)[grep (possible_h[k], names (table (h)))])
  }
  for (l in 1:n_sim) {
    kk <- apply (dd1214_gof_p[l, , , 2] * dd1214_gof_psi[l, , 2], 1, function (x){ paste(x[1], x[2], x[3], sep = "") })
    if( any (names (table(kk)) == possible_h[k]) ) {
      ppc[k, l, 2] <- as.numeric (table (kk)[grep (possible_h[k], names (table (kk)))])
    }
    
    pr_h <- dd1214_gof_psi[l, which (h == possible_h[k]), 1] * dd1214_gof_p[l, which (h == possible_h[k]), , 1]^(history[k, ])*(1 - dd1214_gof_p[l, which (h == possible_h[k]), , 1])^(1 - history[k, ]) + ifelse(k == 1, 1, 0) * (1 - dd1214_gof_psi[l, which (h == possible_h[k]), 1])
    if (!is.null (nrow (pr_h))) {
      if (nrow (pr_h) == 1){ 
        if (length (which (is.na (pr_h))) == 2) { ppc[k, l, 3] <- pr_h[which (!is.na (pr_h))] }
        else { ppc[k, l, 3] <- sum (prod (pr_h, na.rm = TRUE)) }
      }
      else { ppc[k, l, 3] <- sum (apply (pr_h, 1, prod, na.rm = TRUE)) }
    }
  }
}
rm(i, j, k, l, kk)
ppc[, 1, 1:3]

X2_obs <- apply (((ppc[, , 1] - ppc[, , 3])^2) / (ppc[, , 3]), 2, sum, na.rm = TRUE)
X2_rep <- apply (((ppc[, , 2] - ppc[, , 3])^2) / (ppc[, , 3]), 2, sum, na.rm = TRUE)
mean (ifelse (X2_rep > X2_obs, 1, 0)) # 0.436
get_ci (X2_obs / X2_rep) # 0.381 1.310 1.815 (median 1.097)
median (X2_obs / X2_rep)
