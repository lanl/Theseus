## Dave Osthus
## 7-27-23
## Minimum working example for the manuscript "Towards Improved Heliosphere Sky Map Estimation with Theseus"
## This script will reproduce Figures 2, 3, 4, 6, and 7

# © 2023. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
# National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
# Department of Energy/National Nuclear Security Administration. All rights in the program are.
# reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
# Security Administration. The Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
# derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
# others to do so.

###########################
## install packages if not already installed
list.of.packages <- c("ggplot2","grid","gridExtra","mgcv","ggpubr","glmnet")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(new.packages)
if(length(new.packages)) install.packages(new.packages)

## load packages
library(ggplot2)
library(grid)
library(gridExtra)
library(mgcv)
library(ggpubr)
library(glmnet)
theme_set(theme_bw())


###########################
# setwd("/Users/dosthus/Documents/ibex/theseus/manuscript/Technometrics/MWE/")

###########################
## source in the functions 
source("theseus_functions.R")

## define terms for proper coordinate plotting
noselongitude <- 265
center = 180-(360 - noselongitude)
orig_360 = seq(0,300,60)
new_360 = orig_360-center+0.01
new_360[new_360<0.01]=new_360[new_360<0.01]+360

###########################
## load the data
data <- readRDS("data_illustration.RDS")

## binned direct event data
Xobs <- data$Xobs

## 2-degree sky map grid
Xpix <- data$Xpix

## matrix relating 2-degree pixels to binned direct event data
Xpsf <- data$Xpsf

## blurring matrix
KK   <- data$KK

## ibex color palette
ibex_palette <- data$ibex_palette

####################################################################################
####################################################################################
###  The below code will reproduce Figures 2, 3, 4, 6, and 7 of the manuscript.  ###
####################################################################################
####################################################################################

#############################
## bootstrap the data set
## setting bootstrap == F will reproduce the figure in the paper
## setting bootstrap == T will return a random bootstrap data set and the figures will look different
step0 <- parametric_bootstrap_data(Xobs, Xpsf, Xpix, bootstrap=F)
Xobsa <- step0$Xobsa
Xpsfa <- step0$Xpsfa
Xpixa <- step0$Xpixa


#######################################################################
#######################################################################
###  Theseus Step 1 of Stage 1: Estimate Candidate blurred sky map  ###
###  Make Figure 2  ###################################################
#######################################################################
#######################################################################

## compute the method of moments
Xobsa$ena_rate_mom <- (Xobsa$new_counts/Xobsa$time - Xobsa$bck)

## make the weights for the candidate fitting
Xobsa$time_wt <- Xobsa$time/sum(Xobsa$time)

############################################
## define the non-parametric regression scenarios
stage1componentsdf <- expand.grid(wt = c("nowt","wt"),
                                  smoother = c("supsmu","gcvspline","gamtecr","gamteps"))
stage1componentsdf$regmod <- "gam"
stage1componentsdf[stage1componentsdf$smoother %in% c("supsmu","gcvspline"),]$regmod <- "ppr"

## create scenario and quandrant name
stage1componentsdf$scenarioname <- paste0(stage1componentsdf$wt,"_",stage1componentsdf$smoother)
stage1componentsdf$quadrantname <- paste0(stage1componentsdf$regmod)


#### Stage 1; Step 1: Estimate candidate maps
## fit ppr/gam to each scenario
for(kk in 1:nrow(stage1componentsdf)){
  
  ## get names
  tempwt          <- as.character(stage1componentsdf$wt[kk])
  tempsmoother    <- as.character(stage1componentsdf$smoother[kk])
  tempfeaturename <- as.character(stage1componentsdf$scenarioname[kk])
  
  ## set weight vector
  if(tempwt == "nowt"){
    mypprweights = rep(1, nrow(Xobsa))
  }else{
    mypprweights = Xobsa$time_wt
  }
  
  ## get new feature
  temppprfit <- getnonparreg(myXobs = Xobsa,
                             myXpix = Xpixa, 
                             myKK = KK, 
                             mysm.method = tempsmoother,
                             mypprnterms = 100, 
                             mypprweights = mypprweights, 
                             myfeaturename = tempfeaturename)
  
  
  ## update data
  Xpixa <- data.frame(temppprfit$myXpix)
  Xobsa <- data.frame(temppprfit$myXobs)
  
  print(paste0("Finished with Candidate ",kk))
}

## get the ENA rate limits for plotting
step1mn <- min(c(Xpixa$nowt_supsmu_combined,  Xpixa$wt_supsmu_combined, Xpixa$nowt_gcvspline_combined, Xpixa$wt_gcvspline_combined,
                 Xpixa$nowt_gamtecr_combined, Xpixa$wt_gamtecr_combined, Xpixa$nowt_gamteps_combined, Xpixa$wt_gamteps_combined))
step1mx <- 1.1*max(c(Xpixa$nowt_supsmu_combined, Xpixa$wt_supsmu_combined,Xpixa$nowt_gcvspline_combined, Xpixa$wt_gcvspline_combined,
                 Xpixa$nowt_gamtecr_combined, Xpixa$wt_gamtecr_combined, Xpixa$nowt_gamteps_combined, Xpixa$wt_gamteps_combined))
step1breaksmin <- max(0,round(step1mn + .01,2))
step1breaksmid <- round(0.5*(step1mn + step1mx),2)
step1breaksmax <- round(step1mx - .01,2)


## Make Figure 2
fig2a <- ggpubr::ggarrange(
  ##
  ggplot(data=Xpixa)+
    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=nowt_supsmu_combined))+
    scale_fill_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
    ggtitle("Candidate 1 (PPR)")+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5))+
    scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
    ylab("Latitude")+
    xlab("Longitude"),
  ggplot(data=Xpixa)+
    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=wt_supsmu_combined))+
    scale_fill_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
    ggtitle("Candidate 2 (PPR)")+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5))+
    scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
    ylab("Latitude")+
    xlab("Longitude"),
  ##
  ggplot(data=Xpixa)+
    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=nowt_gcvspline_combined))+
    scale_fill_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
    ggtitle("Candidate 3 (PPR)")+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5))+
    scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
    ylab("Latitude")+
    xlab("Longitude"),
  ggplot(data=Xpixa)+
    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=wt_gcvspline_combined))+
    scale_fill_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
    ggtitle("Candidate 4 (PPR)")+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5))+
    scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
    ylab("Latitude")+
    xlab("Longitude"),
  ####
  ggplot(data=Xpixa)+
    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=nowt_gamtecr_combined))+
    scale_fill_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
    ggtitle("Candidate 5 (GAM)")+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5))+
    scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
    ylab("Latitude")+
    xlab("Longitude"),
  ggplot(data=Xpixa)+
    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=wt_gamtecr_combined))+
    scale_fill_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
    ggtitle("Candidate 6 (GAM)")+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5))+
    scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
    ylab("Latitude")+
    xlab("Longitude"),
  ####
  ggplot(data=Xpixa)+
    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=nowt_gamteps_combined))+
    scale_fill_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
    ggtitle("Candidate 7 (GAM)")+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5))+
    scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
    ylab("Latitude")+
    xlab("Longitude"),
  ggplot(data=Xpixa)+
    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=wt_gamteps_combined))+
    scale_fill_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
    ggtitle("Candidate 8 (GAM)")+
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5))+
    scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
    ylab("Latitude")+
    xlab("Longitude"),nrow=2, ncol=4, legend="bottom", common.legend=T)


## make the binned direct event data plot   
fig2b <- ggpubr::ggarrange(grid.rect(gp=gpar(col="white")),
                           ggplot(data=Xobsa)+
                             geom_point(aes(x=ecliptic_lon_center, y=ecliptic_lat, color=ena_rate_mom, size=I(1)))+
                             scale_color_gradientn(colors = ibex_palette$hex[-c(1,nrow(ibex_palette))], limits=c(step1mn, step1mx), breaks=c(step1breaksmin, step1breaksmid, step1breaksmax), name="ENAs/sec")+
                             ggtitle("Simulated Data")+
                             theme(legend.position="bottom",
                                   plot.title = element_text(hjust = 0.5))+
                             scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
                             scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
                             ylab("Latitude")+
                             xlab("Longitude"),
                           grid.rect(gp=gpar(col="white")),ncol=1,legend="bottom", common.legend=T, heights = c(1,4,1))


## Save Figure 2
fig2 <- arrangeGrob(fig2a, fig2b,  widths=c(4,1.25), nrow=1)
print("Save Fig 2")
ggsave(filename = "fig2.pdf", plot = fig2, width = 16, height = 6)           
             
             
             


############################################################################################################
############################################################################################################
###  Theseus Step 2 of Stage 1: Estimate a single, initial blurred sky map as an ensemble of candidates  ###
###  Make Figure 3  ########################################################################################                                               
############################################################################################################
############################################################################################################

###############################
#### Stage 1, Step 2: Combine the candidates with a gam
## do a linear adjustment
Xobsa$offset <- Xobsa$time*Xobsa$bck
mygamnames <- c(paste0(stage1componentsdf$scenarioname,"_combined"),"time","offset")
mybs <- "cr"

## Most of the time (but not always), the model of Equations 6-8 can be fit (due to the quasipoisson identity link)
fitgam <- try(suppressWarnings(bam(counts ~ time + 
                                            s(nowt_supsmu_combined,    bs=mybs, by=time) +
                                            s(wt_supsmu_combined,      bs=mybs, by=time) +
                                            s(nowt_gcvspline_combined, bs=mybs, by=time) +
                                            s(wt_gcvspline_combined,   bs=mybs, by=time) +
                                            s(nowt_gamtecr_combined,   bs=mybs, by=time) +
                                            s(wt_gamtecr_combined,     bs=mybs, by=time) +
                                            s(nowt_gamteps_combined,   bs=mybs, by=time) +
                                            s(wt_gamteps_combined,     bs=mybs, by=time) - 1,
                                   offset = offset, 
                                   data=Xobsa, 
                                   method = "fREML",
                                   family=quasipoisson(link="identity"),
                                   gamma=1.4)), silent=T)

## If class(fitgam) == "try-error", then there is non-finite deviance and some covariates need do be culled until it evaluates to true
if(class(fitgam)[1] == "try-error"){
  ## fit a gam to each component individually and make the initial map as a weighted average
  estXobsa <- rep(0, nrow(Xobsa))
  estXpixa <- rep(0, nrow(Xpixa))
  cnt <- 0
  for(jj in 1:nrow(stage1componentsdf)){
    tempgamdf <- subset(Xobsa, select=c("counts","time","offset",paste0(stage1componentsdf$scenarioname[jj],"_combined")))
    names(tempgamdf) <- c("counts","time","offset","cov")
    fitgamtemp <- try(suppressWarnings(bam(counts ~ time + 
                                                    s(cov, bs=mybs, by=time) - 1,
                                           offset = offset, 
                                           data = tempgamdf, 
                                           method = "fREML",
                                           family = quasipoisson(link="identity"),
                                           gamma = 1.4)), silent=T)
    if(class(fitgamtemp)[1] != "try-error"){
      cnt <- cnt + 1
      print(cnt)
      
      ## extract fits
      llmm <- {
        pdf(NULL)
        res <- plot(fitgam, pages=1)
        invisible(dev.off())
        res
      }
      
      estXobsa <- estXobsa + (as.numeric(fitgamtemp$coefficients[1]) + approx(x = llmm[[1]]$x, y=llmm[[1]]$fit, xout = tempgamdf$cov)$y)
      estXpixa <- estXpixa + (as.numeric(fitgamtemp$coefficients[1]) + approx(x = llmm[[1]]$x, y=llmm[[1]]$fit, xout = pmin(max(tempgamdf$cov),pmax(min(tempgamdf$cov),Xpixa[,paste0(stage1componentsdf$scenarioname[jj],"_combined")])), rule=2)$y)
    }
  }
  if(cnt > 2){
    ## if at least 6 of the candidate maps can be fit via gam, then take an average of the gam fits
    Xobsa$init_step1_estimated_blurred_map <- estXobsa/cnt
    Xpixa$init_step1_estimated_blurred_map <- estXpixa/cnt
  }else{
    ## if fewer than 3 individual gams fit, then set the map as the average of the 8 candidate maps
    Xobsa$init_step1_estimated_blurred_map <- pmax(0,rowMeans(subset(Xobsa, select=paste0(stage1componentsdf$scenarioname,"_combined"))))
    Xpixa$init_step1_estimated_blurred_map <- pmax(0,rowMeans(subset(Xpixa, select=paste0(stage1componentsdf$scenarioname,"_combined"))))
  }
}else{
  ## if fitgam has a finite deviance, proceed
  fitgam <- bam(counts ~ time + 
                         s(nowt_supsmu_combined,    bs=mybs, by=time) +
                         s(wt_supsmu_combined,      bs=mybs, by=time) +
                         s(nowt_gcvspline_combined, bs=mybs, by=time) +
                         s(wt_gcvspline_combined,   bs=mybs, by=time) +
                         s(nowt_gamtecr_combined, bs=mybs, by=time) +
                         s(wt_gamtecr_combined,   bs=mybs, by=time) +
                         s(nowt_gamteps_combined, bs=mybs, by=time) +
                         s(wt_gamteps_combined,   bs=mybs, by=time) - 1,
                offset = offset, 
                data=Xobsa, 
                method = "fREML",
                family=quasipoisson(link="identity"),
                gamma=1.4)
  
  ## extract fits
  llmm <- {
    pdf(NULL)
    res <- plot(fitgam, pages=1)
    invisible(dev.off())
    res
  }
  
  ## make blurred fit to observations
  Xobsa$init_step1_estimated_blurred_map <- pmax(0, as.numeric(fitgam$coefficients[1]) + 
                                                    approx(x = llmm[[1]]$x, y=llmm[[1]]$fit, xout = Xobsa$nowt_supsmu_combined)$y +
                                                    approx(x = llmm[[2]]$x, y=llmm[[2]]$fit, xout = Xobsa$wt_supsmu_combined)$y +
                                                    approx(x = llmm[[3]]$x, y=llmm[[3]]$fit, xout = Xobsa$nowt_gcvspline_combined)$y +
                                                    approx(x = llmm[[4]]$x, y=llmm[[4]]$fit, xout = Xobsa$wt_gcvspline_combined)$y +
                                                    approx(x = llmm[[5]]$x, y=llmm[[5]]$fit, xout = Xobsa$nowt_gamtecr_combined)$y +
                                                    approx(x = llmm[[6]]$x, y=llmm[[6]]$fit, xout = Xobsa$wt_gamtecr_combined)$y + 
                                                    approx(x = llmm[[7]]$x, y=llmm[[7]]$fit, xout = Xobsa$nowt_gamteps_combined)$y +
                                                    approx(x = llmm[[8]]$x, y=llmm[[8]]$fit, xout = Xobsa$wt_gamteps_combined)$y)
  
  
  ## make Step 2, Stage 1 blurred map
  Xpixa$init_step1_estimated_blurred_map <- pmax(0, as.numeric(fitgam$coefficients[1]) + 
                                                    approx(x = llmm[[1]]$x, y=llmm[[1]]$fit, xout = pmin(max(Xobsa$nowt_supsmu_combined),pmax(min(Xobsa$nowt_supsmu_combined),Xpixa$nowt_supsmu_combined)), rule=2)$y +
                                                    approx(x = llmm[[2]]$x, y=llmm[[2]]$fit, xout = pmin(max(Xobsa$wt_supsmu_combined),pmax(min(Xobsa$wt_supsmu_combined),Xpixa$wt_supsmu_combined)), rule=2)$y + 
                                                    approx(x = llmm[[3]]$x, y=llmm[[3]]$fit, xout = pmin(max(Xobsa$nowt_gcvspline_combined),pmax(min(Xobsa$nowt_gcvspline_combined),Xpixa$nowt_gcvspline_combined)), rule=2)$y +
                                                    approx(x = llmm[[4]]$x, y=llmm[[4]]$fit, xout = pmin(max(Xobsa$wt_gcvspline_combined),pmax(min(Xobsa$wt_gcvspline_combined),Xpixa$wt_gcvspline_combined)), rule=2)$y +
                                                    approx(x = llmm[[5]]$x, y=llmm[[5]]$fit, xout = pmin(max(Xobsa$nowt_gamtecr_combined),pmax(min(Xobsa$nowt_gamtecr_combined),Xpixa$nowt_gamtecr_combined)), rule=2)$y +
                                                    approx(x = llmm[[6]]$x, y=llmm[[6]]$fit, xout = pmin(max(Xobsa$wt_gamtecr_combined),pmax(min(Xobsa$wt_gamtecr_combined),Xpixa$wt_gamtecr_combined)), rule=2)$y + 
                                                    approx(x = llmm[[7]]$x, y=llmm[[7]]$fit, xout = pmin(max(Xobsa$nowt_gamteps_combined),pmax(min(Xobsa$nowt_gamteps_combined),Xpixa$nowt_gamteps_combined)), rule=2)$y +
                                                    approx(x = llmm[[8]]$x, y=llmm[[8]]$fit, xout = pmin(max(Xobsa$wt_gamteps_combined),pmax(min(Xobsa$wt_gamteps_combined),Xpixa$wt_gamteps_combined)), rule=2)$y)
  
  #######################################
  ## Plot the components of 
  Xpixa$stage1_step2_intercept                <- as.numeric(fitgam$coefficients[1]) 
  Xpixa$stage1_step2_nowt_supsmu_component    <- approx(x = llmm[[1]]$x, y=llmm[[1]]$fit, xout = pmin(max(Xobsa$nowt_supsmu_combined),pmax(min(Xobsa$nowt_supsmu_combined),Xpixa$nowt_supsmu_combined)), rule=2)$y
  Xpixa$stage1_step2_wt_supsmu_component      <- approx(x = llmm[[2]]$x, y=llmm[[2]]$fit, xout = pmin(max(Xobsa$wt_supsmu_combined),pmax(min(Xobsa$wt_supsmu_combined),Xpixa$wt_supsmu_combined)), rule=2)$y
  Xpixa$stage1_step2_nowt_gcvspline_component <- approx(x = llmm[[3]]$x, y=llmm[[3]]$fit, xout = pmin(max(Xobsa$nowt_gcvspline_combined),pmax(min(Xobsa$nowt_gcvspline_combined),Xpixa$nowt_gcvspline_combined)), rule=2)$y
  Xpixa$stage1_step2_wt_gcvspline_component   <- approx(x = llmm[[4]]$x, y=llmm[[4]]$fit, xout = pmin(max(Xobsa$wt_gcvspline_combined),pmax(min(Xobsa$wt_gcvspline_combined),Xpixa$wt_gcvspline_combined)), rule=2)$y
  Xpixa$stage1_step2_nowt_gamtecr_component   <- approx(x = llmm[[5]]$x, y=llmm[[5]]$fit, xout = pmin(max(Xobsa$nowt_gamtecr_combined),pmax(min(Xobsa$nowt_gamtecr_combined),Xpixa$nowt_gamtecr_combined)), rule=2)$y
  Xpixa$stage1_step2_wt_gamtecr_component     <- approx(x = llmm[[6]]$x, y=llmm[[6]]$fit, xout = pmin(max(Xobsa$wt_gamtecr_combined),pmax(min(Xobsa$wt_gamtecr_combined),Xpixa$wt_gamtecr_combined)), rule=2)$y
  Xpixa$stage1_step2_nowt_gamteps_component   <- approx(x = llmm[[7]]$x, y=llmm[[7]]$fit, xout = pmin(max(Xobsa$nowt_gamteps_combined),pmax(min(Xobsa$nowt_gamteps_combined),Xpixa$nowt_gamteps_combined)), rule=2)$y
  Xpixa$stage1_step2_wt_gamteps_component     <- approx(x = llmm[[8]]$x, y=llmm[[8]]$fit, xout = pmin(max(Xobsa$wt_gamteps_combined),pmax(min(Xobsa$wt_gamteps_combined),Xpixa$wt_gamteps_combined)), rule=2)$y
  
  ## get limits for plot
  step2min <- min(c(Xpixa$stage1_step2_nowt_supsmu_component, Xpixa$stage1_step2_wt_supsmu_component, Xpixa$stage1_step2_nowt_gcvspline_component, Xpixa$stage1_step2_wt_gcvspline_component,
                    Xpixa$stage1_step2_nowt_gamtecr_component, Xpixa$stage1_step2_wt_gamtecr_component, Xpixa$stage1_step2_nowt_gamteps_component, Xpixa$stage1_step2_wt_gamteps_component))
  step2max <- max(c(Xpixa$stage1_step2_nowt_supsmu_component, Xpixa$stage1_step2_wt_supsmu_component, Xpixa$stage1_step2_nowt_gcvspline_component, Xpixa$stage1_step2_wt_gcvspline_component,
                    Xpixa$stage1_step2_nowt_gamtecr_component, Xpixa$stage1_step2_wt_gamtecr_component, Xpixa$stage1_step2_nowt_gamteps_component, Xpixa$stage1_step2_wt_gamteps_component))
  step2breaksmin <- round(step2min +.005, 2)
  step2breaksmid <- round(0.5*(step2min + step2max))
  step2breaksmax <- round(step2max - .005, 2)
  

  ### Theseus Stage 1; Step 2: plot each map's contribution
  fig3a <- ggpubr::ggarrange(
    ggplot(data=Xpixa)+
      geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=stage1_step2_nowt_supsmu_component))+
      scale_fill_gradient2(limits=c(step2min, step2max), name="", breaks=c(step2breaksmin, step2breaksmid, step2breaksmax))+
      ggtitle("Candidate 1 (PPR)")+
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5))+
      scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
      scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
      ylab("Latitude")+
      xlab("Longitude"),
    ggplot(data=Xpixa)+
      geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=stage1_step2_wt_supsmu_component))+
      scale_fill_gradient2(limits=c(step2min, step2max), name="", breaks=c(step2breaksmin, step2breaksmid, step2breaksmax))+
      ggtitle("Candidate 2 (PPR)")+
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5))+
      scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
      scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
      ylab("Latitude")+
      xlab("Longitude"),
    ##
    ggplot(data=Xpixa)+
      geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=stage1_step2_nowt_gcvspline_component))+
      scale_fill_gradient2(limits=c(step2min, step2max), name="", breaks=c(step2breaksmin, step2breaksmid, step2breaksmax))+
      ggtitle("Candidate 3 (PPR)")+
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5))+
      scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
      scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
      ylab("Latitude")+
      xlab("Longitude"),
    ggplot(data=Xpixa)+
      geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=stage1_step2_wt_gcvspline_component))+
      scale_fill_gradient2(limits=c(step2min, step2max), name="", breaks=c(step2breaksmin, step2breaksmid, step2breaksmax))+
      ggtitle("Candidate 4 (PPR)")+
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5))+
      scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
      scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
      ylab("Latitude")+
      xlab("Longitude"),
    ####
    ggplot(data=Xpixa)+
      geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=stage1_step2_nowt_gamtecr_component))+
      scale_fill_gradient2(limits=c(step2min, step2max), name="", breaks=c(step2breaksmin, step2breaksmid, step2breaksmax))+
      ggtitle("Candidate 5 (GAM)")+
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5))+
      scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
      scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
      ylab("Latitude")+
      xlab("Longitude"),
    ggplot(data=Xpixa)+
      geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=stage1_step2_wt_gamtecr_component))+
      scale_fill_gradient2(limits=c(step2min, step2max), name="", breaks=c(step2breaksmin, step2breaksmid, step2breaksmax))+
      ggtitle("Candidate 6 (GAM)")+
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5))+
      scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
      scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
      ylab("Latitude")+
      xlab("Longitude"),
    ####
    ggplot(data=Xpixa)+
      geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=stage1_step2_nowt_gamteps_component))+
      scale_fill_gradient2(limits=c(step2min, step2max), name="", breaks=c(step2breaksmin, step2breaksmid, step2breaksmax))+
      ggtitle("Candidate 7 (GAM)")+
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5))+
      scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
      scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
      ylab("Latitude")+
      xlab("Longitude"),
    ggplot(data=Xpixa)+
      geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=stage1_step2_wt_gamteps_component))+
      scale_fill_gradient2(limits=c(step2min, step2max), name="", breaks=c(step2breaksmin, step2breaksmid, step2breaksmax))+
      ggtitle("Candidate 8 (GAM)")+
      theme(legend.position="bottom",
            plot.title = element_text(hjust = 0.5))+
      scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
      scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
      ylab("Latitude")+
      xlab("Longitude"),ncol=4,nrow=2,common.legend=T,legend="bottom")

  ## Figure 3b
  fig3b <- ggpubr::ggarrange(grid.rect(gp=gpar(col="white")),
                           ggplot(data=Xpixa)+
                             geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=init_step1_estimated_blurred_map))+
                             scale_fill_gradientn(colors=ibex_palette$hex[-c(1, nrow(ibex_palette))], name="ENAs/sec")+
                             ggtitle("Initial Blurred Sky Map")+
                             theme(legend.position="bottom",
                                   plot.title = element_text(hjust = 0.5))+
                             scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
                             scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
                             ylab("Latitude")+
                             xlab("Longitude"),
                           grid.rect(gp=gpar(col="white")),ncol=1, heights=c(1,4,1))
  
  ## Save Figure 3
  ## Note: This may not exactly match Fig 3.
  fig3 <- arrangeGrob(fig3a, fig3b, nrow=1, widths = c(4, 1.25))
  print("Save Fig 3")
  ggsave("fig3.pdf", plot=fig3, width=16, height=6)
}



###################################################################################
###################################################################################
###  Theseus Step 3 of Stage 1: Estimate the residual-adjusted blurred sky map  ###
###  Make Figure 4  ###############################################################
###################################################################################
###################################################################################

## is there any evidence for lack of fit?
Xobsa$step1_lambda <- Xobsa$time*(Xobsa$init_step1_estimated_blurred_map + Xobsa$bck)

## fit the residual map
Xobsa$resid <- (Xobsa$counts - Xobsa$step1_lambda)
residfit <- suppressWarnings(bam(resid ~ te(init_step1_estimated_blurred_map,     bs=mybs) +
                                   te(x,init_step1_estimated_blurred_map,   bs=mybs) +
                                   te(y,init_step1_estimated_blurred_map,   bs=mybs) +
                                   te(z,init_step1_estimated_blurred_map,   bs=mybs) +
                                   te(x, bs=mybs) +
                                   te(y, bs=mybs) +
                                   te(z, bs=mybs) +
                                   te(x,y, bs=mybs) +
                                   te(x,z, bs=mybs) +
                                   te(y,z, bs=mybs),
                                 data=Xobsa, 
                                 method = "fREML",
                                 gamma=1.4))

## add the fitted residual to Xobsa and Xpixa
Xobsa$fit_resid <- predict(residfit, data=Xobsa)
Xpixa$fit_resid <- predict(residfit, newdata=Xpixa)

## make the new blurred map, adjusting for the spatial trend in the residuals
fitgam2 <- try(suppressWarnings(bam(counts ~ time + 
                                             s(init_step1_estimated_blurred_map, bs=mybs, by=time) +
                                             s(fit_resid, bs=mybs, by=time)-1, 
                                    offset = offset, 
                                    data = Xobsa, 
                                    method = "fREML",
                                    family = quasipoisson(link="identity"),
                                    gamma = 1.4)), silent=T)

### update map based on whether gam() was able to fit
if(class(fitgam2)[1] == "try-error"){
  ## if fitgam2 wasn't able to fit, then just set the blurred map to the initial map (i.e., don't residual adjust)
  Xobsa$step1_estimated_blurred_map <- Xobsa$init_step1_estimated_blurred_map
  Xpixa$step1_estimated_blurred_map <- Xpixa$init_step1_estimated_blurred_map
}else{
  ## extract fits
  # llmm2 <- plot.gam(fitgam2,pages=1,plot=F)
  llmm2 <- {
    pdf(NULL)
    res <- plot(fitgam2)
    invisible(dev.off())
    res
  }
  
  ## make blurred fit to observations
  Xobsa$step1_estimated_blurred_map <- pmax(0, as.numeric(fitgam2$coefficients[1]) + 
                                              approx(x = llmm2[[1]]$x, y=llmm2[[1]]$fit, xout = Xobsa$init_step1_estimated_blurred_map)$y +
                                              approx(x = llmm2[[2]]$x, y=llmm2[[2]]$fit, xout = Xobsa$fit_resid)$y)
  
  ## make step 1 blurred map
  Xpixa$step1_estimated_blurred_map_int      <- as.numeric(fitgam2$coefficients[1])
  Xpixa$step1_estimated_blurred_map_mainmap  <- approx(x = llmm2[[1]]$x, y=llmm2[[1]]$fit, xout = pmin(max(Xobsa$init_step1_estimated_blurred_map),pmax(min(Xobsa$init_step1_estimated_blurred_map),Xpixa$init_step1_estimated_blurred_map)), rule=2)$y 
  Xpixa$step1_estimated_blurred_map_residadj <- approx(x = llmm2[[2]]$x, y=llmm2[[2]]$fit, xout = pmin(max(Xobsa$fit_resid),pmax(min(Xobsa$fit_resid),Xpixa$fit_resid)), rule=2)$y
  Xpixa$step1_estimated_blurred_map          <- pmax(0, Xpixa$step1_estimated_blurred_map_int + Xpixa$step1_estimated_blurred_map_mainmap + Xpixa$step1_estimated_blurred_map_residadj)                                                 
  
  
  step3min <- min(c(Xpixa$step1_estimated_blurred_map_int + Xpixa$step1_estimated_blurred_map_mainmap,
                    Xpixa$step1_estimated_blurred_map))
  step3max <- max(c(Xpixa$step1_estimated_blurred_map_int + Xpixa$step1_estimated_blurred_map_mainmap,
                    Xpixa$step1_estimated_blurred_map))
  
  step3breaksmin <- max(0, round(step3min + .01, 2))
  step3breaksmid <- round(0.5*(step3min + step3max),2)
  step3breaksmax <- round(step3max - .01, 2)
  
  
  ## Save Figure 4
  ## Note: This may not exactly match Fig 3.
  fig4 <- arrangeGrob(
            ggplot(data=Xpixa)+
              geom_raster(aes(x = ecliptic_lon_center, y = ecliptic_lat, fill  = step1_estimated_blurred_map_int + step1_estimated_blurred_map_mainmap))+
              scale_fill_gradientn(colors=ibex_palette$hex[-c(1,nrow(ibex_palette))], name="ENAs/sec", limits = c(step3min,step3max), breaks=c(step3breaksmin, step3breaksmid, step3breaksmax))+
              theme(legend.position="bottom",
                    plot.title = element_text(hjust = 0.5))+
              scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
              scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
              ggtitle("Initial Blurred Sky Map Component")+
              xlab("Longitude")+
              ylab("Latitude"),
            ggplot(data=Xpixa)+
              geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=step1_estimated_blurred_map_residadj))+
              scale_fill_gradient2(name="",
                                   breaks = c(round(.95*min(Xpixa$step1_estimated_blurred_map_residadj),3),0,
                                              round(.95*max(Xpixa$step1_estimated_blurred_map_residadj),3)),
                                   limits = c(1.1*min(Xpixa$step1_estimated_blurred_map_residadj),
                                              1.1*max(Xpixa$step1_estimated_blurred_map_residadj)))+
              theme(legend.position="bottom",
                    plot.title = element_text(hjust = 0.5))+
              scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
              scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
              ggtitle("Residual Adjustment Component")+
              xlab("Longitude")+
              ylab("Latitude"),
            ggplot(data=Xpixa)+
              geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=step1_estimated_blurred_map))+
              scale_fill_gradientn(colors=ibex_palette$hex[-c(1,nrow(ibex_palette))],
                                   name="ENAs/sec",
                                   limits = c(step3min,step3max),
                                   breaks=c(step3breaksmin, step3breaksmid, step3breaksmax))+
              theme(legend.position="bottom",
                    plot.title = element_text(hjust = 0.5))+
              scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
              scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
              ggtitle("Residual-Adjusted Blurred Sky Map")+
              xlab("Longitude")+
              ylab("Latitude"),nrow=1)
  print("Save Fig 4")
  ggsave("fig4.pdf", plot=fig4, width=12, height=4)
  
}



###########################################################
###########################################################
###  Theseus Stage 2: Deblur the point spread function  ###
###  Make Figure 6  #######################################
###########################################################
###########################################################

## prepare for ridge regression
Xpixa <- Xpixa[order(Xpixa$pix_id),]

## compute twice blurred map
Xpixa$estimated_blurred2_map <- Xpixa$step1_estimated_blurred_map - as.numeric(KK %*% Xpixa$step1_estimated_blurred_map)

## get penalty factors
## this regularizes smaller pixels more than larger pixels
## which is the same as regularizing pixels near the poles (latitude -90 and 90) more than near the equator
Xpixa$penaltyfactor <- 1/Xpixa$wt_pix

#### fit ridge to get the cross-validated lambda
ridge.cv.fit <- cv.glmnet(x = KK,
                          y = Xpixa$estimated_blurred2_map,
                          intercept = F,
                          alpha = 0, # ridge regression
                          penalty.factor = Xpixa$penaltyfactor, # regularizes smaller pixels more than larger ones
                          lower.limits = -pmax(0,Xpixa$step1_estimated_blurred_map), # ensure the blurred map plus the sharpening map is non-negative
                          nfolds = 5)

## fit ridge with best lambda
ridge.fit <- glmnet(x = KK,
                    y = Xpixa$estimated_blurred2_map,
                    intercept = F,
                    alpha = 0, # ridge regression
                    lambda = ridge.cv.fit$lambda.min,
                    penalty.factor = Xpixa$penaltyfactor, # regularizes smaller pixels more than larger ones
                    lower.limits = -pmax(0,Xpixa$step1_estimated_blurred_map)) # ensure the blurred map plus the sharpening map is non-negative


## add the initial ridge_delta to Xpixa
Xpixa$delta_dot <- as.numeric(ridge.fit$beta)
Xpixa$K_delta_dot <- as.numeric(KK %*% Xpixa$delta_dot)

## delta_dot may be biased as a result of the ridge regularization 
## 
hingebasis = data.frame(b0 = 1, 
                        b1 = Xpixa$delta_dot, 
                        b2 = Xpixa$delta_dot)
hingebasis[hingebasis$b1 >= 0,]$b1 <- 0
hingebasis[hingebasis$b2 <  0,]$b2 <- 0
hingebasis <- as.matrix(hingebasis, ncol=ncol(hingebasis))

## calculate delta_hat
get_delta_hat <- optim(par = c(0,rep(1,ncol(hingebasis)-1)),
                       fn = getdeltahat,
                       hingebasis = hingebasis,
                       blurred_map = Xpixa$step1_estimated_blurred_map,
                       KK = KK,
                       twice_blurred_map = Xpixa$estimated_blurred2_map,
                       control = list(maxit = 10000))

## make delta_hat
Xpixa$delta_hat <- as.numeric(hingebasis %*% as.matrix(get_delta_hat$par, ncol=1))
Xpixa$K_delta_hat <- as.numeric(KK %*% Xpixa$delta_hat)

## compute limits for plotting
stage2deltamin <- min(c(Xpixa$delta_dot, Xpixa$K_delta_dot, Xpixa$estimated_blurred2_map, Xpixa$delta_hat, Xpixa$K_delta_hat))
stage2deltamax <- max(c(Xpixa$delta_dot, Xpixa$K_delta_dot, Xpixa$estimated_blurred2_map, Xpixa$delta_hat, Xpixa$K_delta_hat))

## get limits for plots
stage2sharpeningminx <- min(c(Xpixa$K_delta_dot, Xpixa$K_delta_hat))
stage2sharpeningmaxx <- max(c(Xpixa$K_delta_dot, Xpixa$K_delta_hat))

## get psi breaks
psibreaks <- c(ceiling(min(Xpixa$estimated_blurred2_map)*1000)/1000,
                 0,
                 floor(max(Xpixa$estimated_blurred2_map)*1000)/1000)


## Save Figure 6
fig6 <- arrangeGrob(
          arrangeGrob(
            ggplot(data=Xpixa)+
              geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=estimated_blurred2_map)) +
              scale_fill_gradient2(name=expression(hat(psi)), breaks = psibreaks)+
              theme(legend.position="right")+
              scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
              scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
              xlab("Longitude") + ylab("Latitude")),
            ggplot(data=Xpixa)+
              geom_point(aes(x=K_delta_dot, y=estimated_blurred2_map), color=I("black"), fill=I("darkgrey"), shape=I(21)) +
              geom_abline(color=I("red"))+
              xlab(expression(K*dot(delta)))+
              ylab(expression(hat(psi)))+
              xlim(c(stage2sharpeningminx,stage2sharpeningmaxx)),
            ##
            ggplot(data=Xpixa)+
              geom_point(aes(x=K_delta_hat, y=estimated_blurred2_map), color=I("black"), fill=I("darkgrey"), shape=I(21)) +
              geom_abline(color=I("red"))+
              xlab(expression(K*hat(delta)))+
              ylab(expression(hat(psi)))+
              xlim(c(stage2sharpeningminx,stage2sharpeningmaxx)),
          nrow=1, widths = c(2.5,2,2))
print("Save Fig 6")
ggsave("fig6.pdf",plot=fig6, width=15, height=4)



###########################################################
###########################################################
###  Theseus Stage 2: Shwo the final result of Theseus  ###
###  Make Figure 7  #######################################
###########################################################
###########################################################

## make the unblurred map
Xpixa$estimated_unblurred_map <- pmax(0, Xpixa$step1_estimated_blurred_map + Xpixa$delta_hat)

## blur the unblurred map
Xpixa$estimated_blurred_map <- as.numeric(KK%*%Xpixa$estimated_unblurred_map)

## compute plotting limits
mnstep2 <- min(c(Xpixa$step1_estimated_blurred_map, Xpixa$estimated_unblurred_map, Xpixa$true), na.rm=T)
mxstep2 <- max(c(Xpixa$step1_estimated_blurred_map, Xpixa$estimated_unblurred_map, Xpixa$true), na.rm=T)

## compute more limies
stage2min <- min(c(Xpixa$init_step1_estimated_blurred_map, Xpixa$estimated_blurred_map, Xpixa$estimated_unblurred_map), na.rm=T)
stage2max <- max(c(Xpixa$init_step1_estimated_blurred_map, Xpixa$estimated_blurred_map, Xpixa$estimated_unblurred_map), na.rm=T)
stage2mid <- round((stage2max - stage2min)/2 + stage2min,2)

## Figure 7
fig7 <- arrangeGrob(
          ggplot(data=Xpixa)+
            geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=step1_estimated_blurred_map))+
            scale_fill_gradientn(colors=ibex_palette$hex[-c(1,nrow(ibex_palette))],
                                 name="ENAs/sec", 
                                 limits=c(stage2min, stage2max),
                                 breaks = c(ceiling(stage2min*100)/100, stage2mid, floor(stage2max*100)/100))+
            ylab("Latitude")+
            xlab("Longitude")+
            scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
            scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
            ggtitle("Residual-Adjusted Blurred Sky Map")+
            theme(legend.position = "bottom",
                  plot.title = element_text(hjust = 0.5)),
          ####
          ggplot(data=Xpixa)+
            geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=delta_hat))+
            scale_fill_gradient2(name="",
                                 breaks=c(round(min(Xpixa$delta_hat*.95),3),
                                                                       0,
                                                                       round(max(Xpixa$delta_hat*.95),3)))+
            ggtitle("Sharpening Map")+
            ylab("Latitude")+
            xlab("Longitude")+
            scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
            scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
            theme(legend.position = "bottom",
                  plot.title = element_text(hjust = 0.5)),
          ####
          ggplot(data=Xpixa)+
            geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=estimated_unblurred_map))+
            scale_fill_gradientn(colors=ibex_palette$hex[-c(1,nrow(ibex_palette))],
                                 name="ENAs/sec", 
                                 limits=c(stage2min, stage2max),
                                 breaks = c(ceiling(stage2min*100)/100, stage2mid, floor(stage2max*100)/100))+
            ylab("Latitude")+
            xlab("Longitude")+
            scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
            scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
            ggtitle("Unblurred Sky Map")+
            theme(legend.position = "bottom",
                  plot.title = element_text(hjust = 0.5)),
          ####
          ggplot(data=Xpixa)+
            geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=estimated_blurred_map))+
            scale_fill_gradientn(colors=ibex_palette$hex[-c(1,nrow(ibex_palette))],
                                 name="ENAs/sec", 
                                 limits=c(stage2min, stage2max),
                                 breaks = c(ceiling(stage2min*100)/100, stage2mid, floor(stage2max*100)/100))+
            ylab("Latitude")+
            xlab("Longitude")+
            scale_x_reverse(expand = c(.0,.0), breaks = new_360, labels = orig_360) +
            scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
            ggtitle("Blurred Sky Map")+
            theme(legend.position = "bottom",
                  plot.title = element_text(hjust = 0.5)),nrow=1)
print("Save Fig 7")
ggsave("fig7.pdf",plot=fig7, width=16, height=4)


