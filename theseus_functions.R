## Dave Osthus
## 3-14-23
## Functions for the minimum working example for the 
## manuscript "Towards Improved Heliosphere Sky Map Estimation with Theseus"

# © 2023. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
# National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
# Department of Energy/National Nuclear Security Administration. All rights in the program are.
# reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
# Security Administration. The Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
# derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
# others to do so.


####################################
####################################
## get bootstrap samples
parametric_bootstrap_data <- function(Xobs, Xpsf, Xpix, bootstrap){
  ## Xobs is the binned direct event data
  ## Xpsf is the point spread function
  ## Xpix is the pixel data
  ## bootstrap = T means bootstrap the data set.
  ## bootstrap = F means use the full data set and bck means
  
  ## order observations
  Xobs <- Xobs[order(Xobs$obs_id),]
  Xobs$row_id <- 1:nrow(Xobs)
  
  ## compute MOM
  Xobs$ena_rate <- (Xobs$counts/Xobs$time - Xobs$background)
  Xobs$lambda <- Xobs$time*(Xobs$ena_rate + Xobs$background)
  
  # stratified bootstrap (by orbit)
  if(bootstrap==T){
    unqgroupdf <- ddply(Xobs, .(orbit_number), summarise, size = length(lon))
    bootstrapids <- NULL
    for(jj in 1:nrow(unqgroupdf)){
      tempbootstrapids <- sort(sample(x=Xobs[Xobs$orbit_number == unqgroupdf$orbit_number[jj],]$row_id,
                                      size=unqgroupdf$size[jj],
                                      replace=T))
      bootstrapids <- c(bootstrapids, tempbootstrapids)
    }
  }else{
    bootstrapids <- sort(Xobs$obs_id)
  }

  ## get the new observed ids
  Xobsa <- Xobs[bootstrapids,]
  Xobsa <- Xobsa[order(Xobsa$obs_id),]
  if(bootstrap==T){
    Xobsa$new_counts <- rpois(nrow(Xobsa), lambda=Xobsa$lambda) 
  }else{
    Xobsa$new_counts <- Xobsa$counts
  }
    
  ## update the point spread function matrix, Xpsfa
  Xpsfa <- Xpsf[Xobsa$obs_id,]
  
  ## update Xpix
  Xpixa <- Xpix
  Xpixa$time <- as.numeric(t(Xpsfa)%*%Xobsa$time)
  if(min(Xpixa$time) == 0){
    Xpixa[which(Xpixa$time == 0),]$time <- min(Xpixa$time[which(Xpixa$time > 0)])
  }
  Xpixa$geographic_lon <- Xpixa$ecliptic_lon - 180
  Xpixa$geographic_lat <- Xpixa$ecliptic_lat
  
  ## change name for background
  Xobsa$bck <- Xobsa$background
  
  ## pack up and go
  return(list(Xobsa = Xobsa,
              Xpsfa = Xpsfa,
              Xpixa = Xpixa))
}
  


####################################
####################################
## make non-parametric regression computing function
getnonparreg <- function(myXobs, myXpix, myKK, mysm.method, mypprnterms, mypprweights, myfeaturename){
  ## myXobs is the binned direct event data
  ## myXpix is the pixel sky map data frame
  ## myKK is the pixel sky map blurring matrix
  ## mysm.method determines the candidate map to make
  ## mypprnterms is the nterms argument in ppr()
  ## mypprweights is the weight vector
  ## myfeaturename is the variable name to output

  
  ## myk is the number of basis functions to use in GAM
  myk <- 6 # a compromise between computational time and smoothness
  
  ## mygamma is as suggested on page 224 of Simon Wood's Generalized Additive Models: An Introduction with R to prevent overfitting when using GCV
  mygamma <- 1.4 

  ##fit ppr if mysm.method is supsmu or gcvspline
  if(mysm.method %in% c("supsmu","gcvspline")){
    candidatefit <- ppr(ena_rate_mom ~ x + y + z,
                  data = myXobs,
                  nterms = mypprnterms,
                  weights = mypprweights/mean(mypprweights),
                  sm.method = mysm.method)
  }
  
  ## fit GAM with cubic spline basis if mysm.method is gamtecr
  if(mysm.method == "gamtecr"){
    candidatefit <- suppressWarnings(bam(ena_rate_mom ~ te(x, bs="cr", k=myk) + te(y, bs="cr", k=myk) + te(z, bs="cr", k=myk) +
                                                  te(x,y, bs="cr", k=myk) + te(x,z, bs="cr", k=myk) + te(y,z, bs="cr", k=myk),
                                                  gamma = mygamma,
                                                  data = myXobs,
                                                  method = "fREML",
                                                  # discrete = T,
                                                  weights = mypprweights/mean(mypprweights)))
  }
  
  ## fit GAM with P-spline basis if mysm.method is gamteps
  if(mysm.method == "gamteps"){
    candidatefit <- suppressWarnings(bam(ena_rate_mom ~ te(x, bs="ps", k=myk) + te(y, bs="ps", k=myk) + te(z, bs="ps", k=myk) +
                                                  te(x,y, bs="ps", k=myk) + te(x,z, bs="ps", k=myk) + te(y,z, bs="ps", k=myk),
                                   gamma = mygamma,
                                   data = myXobs,
                                   method = "fREML",
                                   # discrete = T,
                                   weights = mypprweights/mean(mypprweights)))
  }



  ## add in fitted values
  myXobs$candidatefit <- pmax(0,as.numeric(predict(candidatefit, newdata=myXobs)))
  myXpix$candidatefit <- pmax(0,as.numeric(predict(candidatefit, newdata=myXpix)))


  ## rename
  names(myXpix)[which(names(myXpix) == "candidatefit")] <- paste0(myfeaturename,"_combined")
  names(myXobs)[which(names(myXobs) == "candidatefit")] <- paste0(myfeaturename,"_combined")



  return(list(myXpix = myXpix,
              myXobs = myXobs))
}



####################################
####################################
## estimate delta hat
getdeltahat <- function(x, hingebasis, blurred_map, twice_blurred_map, KK){
  ## x are the MARS regression coefficients
  ## hingebasis is the MARS design matrix
  ## blurred_map is the ppr() estimated blurred map from step 1
  ## twice_blurred_map is psi hat
  ## KK is the blurring matrix
  
  ## get initial delta
  delta <- as.numeric(hingebasis %*% matrix(x,ncol=1))
  
  ## make observation signal
  estimated_twice_blurred_map <- as.numeric(KK %*% delta)
  
  ## compute weighted mse
  wt <- length(twice_blurred_map)*(abs(twice_blurred_map)/sum(abs(twice_blurred_map)))
  wmse <- mean(wt*(twice_blurred_map - estimated_twice_blurred_map)^2)
  
  ## pack up and go
  return(wmse)
  
}
