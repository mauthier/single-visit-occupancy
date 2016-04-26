### useful functions
get_ci <- function(x, alpha=0.80){ return(round(c(HPDinterval(as.mcmc(x),prob=alpha)[1],mean(x),HPDinterval(as.mcmc(x),prob=alpha)[2]), 3)) }
which.cell <- function(y,x){sapply(1:length(x),function(i){which(y==x[i])})}

# rescale covariates
rescale <- function(x){(x-mean(x))/(2*sd(x))}
unscale <- function(x,y){2*sd(x, na.rm=TRUE)*y+mean(x, na.rm=TRUE)}

# link functions
inv_logit <- function(x){1/(1+exp(-x))}
logit <- function(x){log(x/(1-x))}

### Fit models
### imperfect detection, no iCAR model, no covariate at all
hSDM_siteocc <- function(dd, dd_suit, cov=FALSE){
  if(cov){
    det_formula <- '~ 1 + I(rescale(Effort_km)) + I(Subj/2)'
    occ_formula <- '~ 1 + I(rescale(Bathy)) + I(rescale(Lat)) + I(rescale(Iso200))'
    betaMu  <- rep(0,4); betaVar <-  c(5,rep(2.5,3))^2 ;
    gammaMu <- rep(0,3); gammaVar <- c(5,rep(2.5,2))^2 ;
  }
  else{
    det_formula <- occ_formula <- '~ 1'
    betaMu <-  rep(0,1); betaVar <-  c(5,rep(2.5,0))^2 ;
    gammaMu <- rep(0,1); gammaVar <- c(5,rep(2.5,0))^2 ;
  }
  return(replicate(3,
                   hSDM.siteocc(# Observations
                     presence=dd$Detected,
                     observability=as.formula(det_formula),
                     site=dd$Site,
                     data.observability=dd,
                     # Habitat
                     suitability=as.formula(occ_formula), 
                     data.suitability=dd_suit,
                     # Predictions
                     suitability.pred=NULL,
                     # Chains
                     burnin=n_burnin, mcmc=n_mcmc, thin=n_thin,
                     # Starting values
                     beta.start=0, gamma.start=0,
                     # Priors --> informative normal: exploit rescaling as in Gelman & al. 2008
                     mubeta=betaMu, Vbeta=betaVar, mugamma=gammaMu, Vgamma=gammaVar,
                     # Various
                     verbose=1, save.p=0, seed=runif(1,1,19811021)
                   ), simplify=FALSE
                   )
         )
}
hSDM_siteocc_iCAR <- function(dd, dd_suit, cov=FALSE){
  if(cov){
    det_formula <- '~ 1 + I(rescale(Effort_km)) + I(Subj/2)'
    occ_formula <- '~ 1 + I(rescale(Bathy)) + I(rescale(Lat)) + I(rescale(Iso200))'
    betaMu  <- rep(0,4); betaVar <-  c(5,rep(2.5,3))^2 ;
    gammaMu <- rep(0,3); gammaVar <- c(5,rep(2.5,2))^2 ;
    prior="1/Gamma";
  }
  else{
    det_formula <- occ_formula <- '~ 1'
    betaMu <-  0 ; betaVar  <- 5^2 ;
    gammaMu <- 0 ; gammaVar <- 5^2 ;
    prior="Uniform";
  }
  return(replicate(3,
                   hSDM.siteocc.iCAR(# Observations
                     presence=dd$Detected,
                     observability=as.formula(det_formula),
                     site=dd$Site,
                     data.observability=dd,
                     # Habitat
                     suitability=as.formula(occ_formula),
                     data.suitability=dd_suit,
                     # Spatial structure
                     spatial.entity=dd_suit$cells, n.neighbors=n.neighbors, neighbors=adj,
                     # Predictions
                     suitability.pred=NULL,
                     # Chains
                     burnin=n_burnin, mcmc=n_mcmc, thin=n_thin,
                     # Starting values
                     beta.start=0, gamma.start=0, Vrho.start=1,
                     # Priors --> informative normal: exploit rescaling as in Gelman & al. 2008
                     mubeta=betaMu, Vbeta=betaVar, mugamma=gammaMu, Vgamma=gammaVar,
                     priorVrho=prior, shape=2, rate=1, Vrho.max=(pi*pi)/6,
                     # Various
                     verbose=1, save.p=0, save.rho=0, seed=runif(1,1,19811021)
                   ), simplify = FALSE
                   )
         )
}

### perfect detection
hSDM_binom <- function(dd, cov=FALSE){
  if(cov){
    occ_formula <- '~ 1 + I(rescale(Effort_km)) + I(Subj/2) + I(rescale(Bathy)) + I(rescale(Lat)) + I(rescale(Iso200))'
    betaMu  <- rep(0,6); betaVar <-  c(5,rep(2.5,5))^2 ;
  }
  else{
    occ_formula <- '~ 1'
    betaMu <-  rep(0,1); betaVar <-  c(5,rep(2.5,0))^2 ;
  }
  return(replicate(3,
                   hSDM.binomial(# Observations
                     presence=dd$Detected,
                     trials=rep(1, nrow(dd)),
                     # Habitat
                     suitability=as.formula(occ_formula), 
                     data=dd,
                     # Predictions
                     suitability.pred=NULL,
                     # Chains
                     burnin=n_burnin, mcmc=n_mcmc, thin=n_thin,
                     # Starting values
                     beta.start=0,
                     # Priors --> informative normal: exploit rescaling as in Gelman & al. 2008
                     mubeta=betaMu, Vbeta=betaVar,
                     # Various
                     verbose=1, save.p=0, seed=runif(1,1,19811021)
                   ), simplify=FALSE
  )
  )
}

hSDM_binom_iCAR <- function(dd, cov=FALSE){
  if(cov){
    occ_formula <- '~ 1 + I(rescale(Effort_km)) + I(Subj/2) + I(rescale(Bathy)) + I(rescale(Lat)) + I(rescale(Iso200))'
    betaMu  <- rep(0,6); betaVar <-  c(5,rep(2.5,5))^2 ;
    prior="1/Gamma";
  }
  else{
    occ_formula <- '~ 1'
    betaMu <-  0 ; betaVar  <- 5^2 ;
    prior="Uniform";
  }
  return(replicate(3,
                   hSDM.binomial.iCAR(# Observations
                     presence=dd$Detected,
                     trials=rep(1, nrow(dd)),
                     # Habitat
                     suitability=as.formula(occ_formula), 
                     data=dd,
                     # Spatial structure
                     spatial.entity=dd$cells, n.neighbors=n.neighbors, neighbors=adj,
                     # Predictions
                     suitability.pred=NULL,
                     # Chains
                     burnin=n_burnin, mcmc=n_mcmc, thin=n_thin,
                     # Starting values
                     beta.start=0, Vrho.start=1,
                     # Priors --> informative normal: exploit rescaling as in Gelman & al. 2008
                     mubeta=betaMu, Vbeta=betaVar,
                     priorVrho=prior, shape=2, rate=1, Vrho.max=(pi*pi)/6,
                     # Various
                     verbose=1, save.p=0, save.rho=0, seed=runif(1,1,19811021)
                   ), simplify = FALSE
  )
  )
}

# function to predict
psi_hat <- function(model, dd, dd_suit=NULL, cov){
  ### spatial
  spatial <- any(names(as.data.frame(model$mcmc))=="Vrho")
  if(spatial) {
    rho <- matrix( rep( model$rho.pred[ifelse(is.null(dd_suit), dd$cells, dd_suit$cells)], 
                        each=nrow(model$mcmc[-c(1:(n_burnin/n_thin)),])), 
                   nrow=nrow(model$mcmc[-c(1:(n_burnin/n_thin)),]), 
                   ncol=ifelse(is.null(dd_suit), nrow(dd), nrow(dd_suit)), 
                   byrow=FALSE
    )
  }
  else{
    rho <- matrix( rep( 0, 
                        ifelse(is.null(dd_suit), nrow(dd), nrow(dd_suit))*nrow(model$mcmc[-c(1:(n_burnin/n_thin)),])), 
                   nrow=nrow(model$mcmc[-c(1:(n_burnin/n_thin)),]), 
                   ncol=ifelse(is.null(dd_suit), nrow(dd) ,nrow(dd_suit)), 
                   byrow=FALSE
    )
  }
  
  ### binomial likelihood
  if(is.null(dd_suit)){
    ### covariables
    if(cov){
      # occurrence
      psi <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),1] + 
                         model$mcmc[-c(1:(n_burnin/n_thin)),2]%*%t(rescale(dd$Effort_km)) + 
                         model$mcmc[-c(1:(n_burnin/n_thin)),3]%*%t(dd$Subj/2) + 
                         model$mcmc[-c(1:(n_burnin/n_thin)),4]%*%t(rescale(dd$Bathy)) + 
                         model$mcmc[-c(1:(n_burnin/n_thin)),5]%*%t(rescale(dd$Lat)) + 
                         model$mcmc[-c(1:(n_burnin/n_thin)),6]%*%t(rescale(dd$Iso200)) + 
                         rho
      )
    }
    else{
      # occurrence
      psi <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),1] + rho)
    }
  }
  else{
    ### covariables
    if(cov){
      # occurrence
      psi <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),1] + 
                         model$mcmc[-c(1:(n_burnin/n_thin)),2]%*%t(rescale(dd_suit$Bathy)) + 
                         model$mcmc[-c(1:(n_burnin/n_thin)),3]%*%t(rescale(dd_suit$Lat)) + 
                         model$mcmc[-c(1:(n_burnin/n_thin)),4]%*%t(rescale(dd_suit$Iso200)) + 
                         rho
      )
    }
    else{
      # occurrence
      psi <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),1] + rho)
    }
  }
  return(psi)
}

### WAIC
compute_waic <- function(model, dd, dd_suit=NULL, cov){
  WAIC <- function(data, theta){
    lppd <- 0
    pD <- 0
    for ( i in 1:ncol(theta) ) {
      # compute log of average likelihood
      ll <- dbinom( data[i], 1, theta[,i], log=FALSE )
      lppd <- lppd + log(mean(ll))
      pD <- pD + var(log(ll))
    }
    return(round(c(lppd, pD, -2*(lppd - pD)),2))
  }
  
  # function to predict
  theta_hat <- function(model, dd, dd_suit, cov){
    ### spatial
    spatial <- any(names(as.data.frame(model$mcmc))=="Vrho")
    if(spatial) {
      rho <- matrix( rep( model$rho.pred[ifelse(is.null(dd_suit), dd$cells, dd_suit$cells)], 
                          each=nrow(model$mcmc[-c(1:(n_burnin/n_thin)),])), 
                     nrow=nrow(model$mcmc[-c(1:(n_burnin/n_thin)),]), 
                     ncol=ifelse(is.null(dd_suit), nrow(dd), nrow(dd_suit)), 
                     byrow=FALSE
      )
    }
    else{
      rho <- matrix( rep( 0, 
                          ifelse(is.null(dd_suit), nrow(dd), nrow(dd_suit))*nrow(model$mcmc[-c(1:(n_burnin/n_thin)),])), 
                     nrow=nrow(model$mcmc[-c(1:(n_burnin/n_thin)),]), 
                     ncol=ifelse(is.null(dd_suit), nrow(dd) ,nrow(dd_suit)), 
                     byrow=FALSE
      )
    }
    proba <- array(NA, dim=c(nrow(model$mcmc[-c(1:(n_burnin/n_thin)),]), nrow(dd)))
    
    ### binomial likelihood
    if(is.null(dd_suit)){
      ### covariables
      if(cov){
        # occurrence
        proba <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),1] + 
                             model$mcmc[-c(1:(n_burnin/n_thin)),2]%*%t(rescale(dd$Effort_km)) + 
                             model$mcmc[-c(1:(n_burnin/n_thin)),3]%*%t(dd$Subj/2) + 
                             model$mcmc[-c(1:(n_burnin/n_thin)),4]%*%t(rescale(dd$Bathy)) + 
                             model$mcmc[-c(1:(n_burnin/n_thin)),5]%*%t(rescale(dd$Lat)) + 
                             model$mcmc[-c(1:(n_burnin/n_thin)),6]%*%t(rescale(dd$Iso200)) + 
                             rho
        )
      }
      else{
        # occurrence
        proba <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),1] + rho)
      }
    }
    ### imperfect detection
    else{
      if(cov){
        # detection
        p   <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),5] + 
                           model$mcmc[-c(1:(n_burnin/n_thin)),6]%*%t(rescale(dd$Effort_km)) + 
                           model$mcmc[-c(1:(n_burnin/n_thin)),7]%*%t(dd$Subj/2)
                         )
        # occurrence
        psi <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),1] + 
                           model$mcmc[-c(1:(n_burnin/n_thin)),2]%*%t(rescale(dd_suit$Bathy)) + 
                           model$mcmc[-c(1:(n_burnin/n_thin)),3]%*%t(rescale(dd_suit$Lat)) + 
                           model$mcmc[-c(1:(n_burnin/n_thin)),4]%*%t(rescale(dd_suit$Iso200)) + 
                           rho
                         )
        for(j in 1:ncol(proba)){
          proba[,j] <- p[,j]*psi[,which(dd_suit$Site==dd$Site[j])]
        }
      }
      else{
        # detection
        p   <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),2])
        # occurrence
        psi <- inv_logit(model$mcmc[-c(1:(n_burnin/n_thin)),1] + rho)
        for(j in 1:ncol(proba)){
          proba[,j] <- p*psi[,which(dd_suit$Site==dd$Site[j])]
        }
      }
    }
    return(proba)
  }
  
  return( WAIC(data=dd$Detected, theta=theta_hat(model=model, dd, dd_suit=dd_suit, cov)) )
}
