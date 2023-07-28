# Theseus
Minimum working example for the sky map estimation method Theseus
---
## Introduction

Theseus is a two-stage procedure that estimates spatial maps of heliospheric energetic neutral atom (ENA) rates, referred to as sky map. In Stage 1, Theseus estimates a blurred sky map from noisy and irregularly spaced binned direct event data using an ensemble approach that leverages projection pursuit regression and generalized additive models. In Stage 2, Theseus deblurs the sky map by deconvolving the point spread function (PSF) with the blurred map using regularization. Unblurred sky map uncertainties are computed via bootstrapping.


## System requirements

The code is supported on all operating systems for which the requisite downloads (see below) are possible. The example code was tested on a MacBook Pro running macOS Monterey 12.6.3, using R version 4.2.2.

## Installation

To downloading and install software and packages:
 - R (>= 2.14.0) follow instructions at https://www.r-project.org/

Installation should take less than 15 minutes on a normal desktop computer.


## Demonstration

A reduced data set is provided, named **data_illustration.RDS**. This example code will reproduce Figures 2, 3, 4, 6, and 7 of the manuscript, "Towards Improved Heliosphere Sky Map Estimation with Theseus". The R code to produce the figures is **theseus_illustration.R** and will source in the functions found within **theseus_functions.R**. This should only take 1-2 minutes to run. 

The binned direct event data for the ESA 4 "A" maps along with the Theseus and ISOC sky map estimates are provided in **data_results.RDS**. The code to reproduce Figures 13 and 14 of the manuscript is **theseus_results.R** 

NOTE: Anyone wishing to use the binned direct event data for **space science** (not statistical methodological development) should contact the LANL IBEX team first (email dreisenfeld@lanl.gov).


## Instructions for use

After R is installed, run **theseus_illustration.R** to reproduce Figures 2, 3, 4, 6, and 7, or run **theseus_results.R** to reproduce Figures 13 and 14. 

Users may need to setwd('Path/to/Theseus/Directory/') in line 34 of both .R files.


## Data Details

### data_illustration.RDS

There are **5** data products in **data_illustration.RDS**.

**Xobs** is an 11 column data frame:
- obs_id: is a unique label for the binned direct event data, from 1 to the number of rows of Xobs
- ecliptic_lon: is the ecliptic longitude (between 0 and 360)
- ecliptic_lat: is the ecliptic latitude (between -90 and 90)
- ecliptic_lon_center: is the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- x, y, and z: are the sperical coordinates for the ecliptic longitude and latitude and are used for PPR and GAM regression
- orbit_number: is the orbit number from the data collection and is used in the bootstrap sampling
- counts: are the number of direct events
- time: is the exposure time (in seconds)
- background: is the background rate (background particles per second)

**Xpix** is an 8 column data frame:
- pix_id: is a unique label for each pixel, from 1 to the number of rows of Xpix (for a 2 degree map, that's 16,200)
- ecliptic_lon: is the ecliptic longitude (between 0 and 360)
- ecliptic_lat: is the ecliptic latitude (between -90 and 90)
- ecliptic_lon_center: is the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- x, y, and z: are the sperical coordinates for the ecliptic longitude and latitude and are used for PPR and GAM regression
- wt_pix: is proportional to the size of the pixel (area on a unit sphere) where larger pixels are near ecliptic_latitude 0 and smaller pixels are near ecliptic latitudes -90 and 90

**Xpsf** is a nrow(Xobs) by nrow(Xpix) sparse matrix where each entry is non-neagative and each row sums to 1. It is used to relate the unblurred sky map to the blurred binned direct events

**KK** is an nrow(Xpix) by nrow(Xpix) sparse matrix where each entry is non-negative and each row sums to 1. It is used to relate the unblurred sky map to the blurred sky map

**ibex_palette** is a 4 column data frame:
- red: red numeric value
- green: green numeric value
- blue: blue numeric value
- hex: hex value

### data_results.RDS 

There are **2** data products in **data_results.RDS**.

**ibex_data** is an 11 column data frame:
- esa: is the energy step for the IBEX data. Only ESA 4 is provided
- map: is the 6-month map id (e.g., "2013A"). Only "A" maps are provided
- ecliptic_lon: is the ecliptic longitude (between 0 and 360)
- ecliptic_lat: is the ecliptic latitude (between -90 and 90)
- ecliptic_lon_center: is the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- x, y, and z: are the sperical coordinates for the ecliptic longitude and latitude
- counts: are the number of direct events
- time: is the exposure time (in seconds)
- background: is the background rate (background particles per second)

**real_sky_maps** is an 8 column data frame:
- esa: is the energy step for the IBEX data. Only ESA 4 is provided
- map: is the 6-month map id (e.g., "2013A"). Only "A" maps are provided
- ecliptic_lon: is the ecliptic longitude (between 0 and 360)
- ecliptic_lat: is the ecliptic latitude (between -90 and 90)
- ecliptic_lon_center: is the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- method: the sky map estimation method. Either "theseus" of "isoc"
- quantity: the type of response, either "mean" or "ci" for 95% confidence interval width
- value: either the ENA rate (ENAs/sec when q

  
## Attribution

If you use the Theseus in your research work, please cite the following paper:

D Osthus, BP Weaver, LJ Beesley, KR Moran, MA Stricklin, EJ Zirnstein, PH Janzen, and DB Reisenfeld. Towards Improved Heliosphere Sky Map Estimation with Theseus. [Under consideration at Technometrics](https://arxiv.org/pdf/2210.12005.pdf).

---
Copyright 2023 for **CO4627**

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
