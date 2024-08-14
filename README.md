# bbmixcp
Title: A Bayesian Beta-Binomial Piecewise Growth Mixture Model for Longitudinal Overdispersed Binomial Data

Journal: Statistical Methods in Medical Research (In press)

Author: Chun-Che Wen, Nathaniel Baker, Rajib Paul, Elizabeth Hill1, Kelly Hunt, Hong Li, Kevin Gray, and Brian Neelon

This repository includes the simulation  for the manuscript submitted in Statistical Methods in Medical Research. The beta-binomial piecewise growth mixture model (BB PGMM) is specifically
designed to address longitudinal overdispersed binomial responses. Within each class, we fit a piecewise linear beta-binomial mixed model with random changepoints for each study group to
detect critical windows of treatment efficacy. Using this model, we can cluster subjects who share similar characteristics, estimate the class-specific mean abstinence trends for each study group,
and quantify the treatment effect over time within each class.

The simulation (Simulation_bbmixcp.R) includes data generating, an MCMC algorithm, label-switching detection, figures (group- and subject-level trends), a triangle plot, and a traceplot.

## Files
Simulation_bbmixcp.R   =  R code simulation with sample size N=250

Simulation_bbmixcp.Rda =  MCMC samples for (generated from Simulation_bbmixcp.R script)

Simulation_bbmixcp.pdf =  R code simulation markdown  

Figures folder = All figures from the manuscript and supplement 

## Parameters

Beta1/Beta2/Beta3   = Fixed effects in each BB class

Rho1/Rho2/Rho3    = Correlation parameter in each BB class (0<rho<1)

Kappa = Changepoints (class 1: pl & tx, class 2: pl & tx, class 3: pl & tx)

Sigmab1/Sigmab2/Sigmab3 = 3 x 3 matrix of Random effect variance (random intercept/slope/slope after cp)
