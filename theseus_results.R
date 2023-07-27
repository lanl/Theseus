## Dave Osthus
## 7-27-23
## Minimum working example for the manuscript "Towards Improved Heliosphere Sky Map Estimation with Theseus"
## This script plots the ESA 4, "A" maps from ISOC and Theseus, point estimates and 95% confidence interval widths

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
setwd("/Users/dosthus/Documents/ibex/theseus/manuscript/Technometrics/MWE/")

## define terms for proper coordinate plotting
noselongitude <- 265
center = 180-(360 - noselongitude)
orig_360 = seq(0,300,60)
new_360 = orig_360-center+0.01
new_360[new_360<0.01]=new_360[new_360<0.01]+360

###########################
## binned direct event data
df <- readRDS("data_results.RDS")$real_sky_maps

## ibex color palette
ibex_palette <- readRDS("data_illustration.RDS")$ibex_palette


###########################
## Plot Figure 13
mn13 <- min(subset(df, type == "mean")$value, na.rm=T)
mx13 <- max(subset(df, type == "mean")$value, na.rm=T)

fig13isoc<- ggplot(data=subset(df, method == "isoc" & type == "mean"))+
                geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=value))+
                facet_wrap(~map, ncol=2)+
                scale_fill_gradientn(colors = ibex_palette$hex, name="ENAs/sec", limits=c(mn13, mx13), na.value = 'darkgrey')+
                scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
                scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
                xlab("Longitude")+
                ylab("Latitude")+
                ggtitle(paste0("ISOC"))+
                theme(legend.position="bottom",
                      plot.title = element_text(hjust = 0.5))

fig13theseus <- ggplot(data=subset(df, method == "theseus" & type == "mean"))+
                    geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=value))+
                    facet_wrap(~map, ncol=2)+
                    scale_fill_gradientn(colors = ibex_palette$hex, name="ENAs/sec", limits=c(mn13, mx13), na.value = 'darkgrey')+
                    scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
                    scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
                    xlab("Longitude")+
                    ylab("Latitude")+
                    ggtitle(paste0("Theseus"))+
                    theme(legend.position="bottom",
                          plot.title = element_text(hjust = 0.5))

## Figure 13
grid.arrange(fig13isoc, fig13theseus, nrow=1)



###########################
## Plot Figure 14
mn14 <- min(subset(df, type == "ci")$value, na.rm=T)
mx14 <- quantile(subset(df, !is.na(value) & type == "ci")$value, probs=.999)
LEVELS <- seq(mn14, mx14,length.out=6)

fig14isoc<- ggplot(data=subset(df, method == "isoc" & type == "ci"))+
              geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=value))+
              facet_wrap(~map, ncol=2)+
              scale_fill_gradientn(name = "95% CI\nWidth",
                                   breaks   = round(c(0,mean(LEVELS),1.1*max(LEVELS)),2),
                                   colors   = c('white', 'gold',  'orange', 'red', 'darkred','black'),
                                   limits   = c(0,1.1*mx14),
                                   na.value = 'darkgrey') +                
              scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
              scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
              xlab("Longitude")+
              ylab("Latitude")+
              ggtitle(paste0("ISOC"))+
              theme(legend.position="bottom",
                    plot.title = element_text(hjust = 0.5))

fig14theseus <- ggplot(data=subset(df, method == "theseus" & type == "ci"))+
                  geom_raster(aes(x=ecliptic_lon_center, y=ecliptic_lat, fill=value))+
                  facet_wrap(~map, ncol=2)+
                  scale_fill_gradientn(name = "95% CI\nWidth",
                                       breaks   = round(c(0,mean(LEVELS),1.1*max(LEVELS)),2),
                                       colors   = c('white', 'gold',  'orange', 'red', 'darkred','black'),
                                       limits   = c(0,1.1*mx14),
                                       na.value = 'darkgrey') +                  
                  scale_x_reverse(expand = c(0,0), breaks = new_360, labels = orig_360)+
                  scale_y_continuous(expand = c(.0,.0), breaks = seq(-45,45,45)) +
                  xlab("Longitude")+
                  ylab("Latitude")+
                  ggtitle(paste0("Theseus"))+
                  theme(legend.position="bottom",
                        plot.title = element_text(hjust = 0.5))

## Figure 14
grid.arrange(fig14isoc, fig14theseus, nrow=1)


