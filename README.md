# AIS Biomechanics Validation Project

This repository contains code to replicate the estimation of discrete and
functional limits of agreement as part of the AIS guide to validation of new
technologies.

To run the associated code you will require the following `R` packages installed:

- `tidyverse`
- `refund`
- `lme4`
- `fda`
- `mvtnorm`
- `viridis`
- `reshape`
- `future`
- `future.apply`
- `FunInf` (requires `remotes` to install as `remotes::install_github("arlohook/FunInf")`)

A brief guide to the files in this repository is as follows:

- **DGP.R** is the script to replicate the data generation process to produce the file *Demo Data.rds* 
- **LoA.R** is the raw R code to run the discrete and functional limits of agreement on the simulated data
- **Limits-of-Agreement.html** is the html file that details the thought process behind these validation methods
