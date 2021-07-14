# JULES-ES-1p0-common.R
# A script to underly JULES-ES-1p0 ensemble work. Loading data, helper functions etc.
# Loads ensemble and constrains ensemble to "level 0", which excludes inputs where the model doesn't run.
# Doug McNeall July 2021


## ----------------------------------------------------------------------
## Load packages
## ----------------------------------------------------------------------

library(RColorBrewer)
library(fields)
library(MASS)
library(DiceKriging)
library(DiceEval)
library(ncdf4)
library(ncdf4.helpers)
library(readxl)

library(foreach)

library(emtools)
library(imptools)
library(viztools)
library(julesR)

# move the right functions here
#source('~/brazilCSSP/code/brazil_cssp/per_pft.R') # eventually, move the relevant functions
#source('explore-JULES-ES-1p0_PPE_functions.R')


## ----------------------------------------------------------------------
## Data locations and constants
## ----------------------------------------------------------------------
ensloc <- '/project/carbon_ppe/JULES-ES-1p0_PPE/'

# Some pallete options
yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(11, "RdYlBu")
byr = rev(ryb)
rb = brewer.pal(11, "RdBu")
br = rev(rb)
blues = brewer.pal(9, 'Blues')
cbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ysec = 60*60*24*365
years <- 1850:2013

## ----------------------------------------------------------------------
## Helper functions
## ----------------------------------------------------------------------


reset <- function() {
  # Allows annotation of graphs, resets axes
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

makeTransparent<-function(someColor, alpha=100)
  # Transparent colours for plotting
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

anomalizeTSmatrix = function(x, ix){
  # Anomalise a timeseries matrix
  subx = x[ ,ix]
  sweepstats = apply(subx, 1, FUN=mean)
  anom = sweep(x, 1, sweepstats, FUN = '-')
  anom
}

extractTimeseries <- function(nc, variable){
  dat <- ncvar_get(nc, variable)
  out <- dat
  out
}

makeTimeseriesEnsemble <- function(variable, nens = 499, nts = 164, cn = 1850:2013){
  
  # nens is number of ensemble members
  # nts length of timeseries
  # cn is colnames()
  datmat <- matrix(NA, nrow = nens, ncol = nts)
  colnames(datmat) <- cn
  
  enslist <- paste("P", formatC(0:(nens-1), width=4, flag="0"), sep="")
  #floc <- paste0(ensloc,ensmember,subdir)
  
  for(i in 1:nens){
    
    vec <- rep(NA,nts)
    
    ensmember <- enslist[i] 
    
    fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    
    
    try(nc <- nc_open(paste0(fn)))
    try(dat <- extractTimeseries(nc, variable))
    
    datmat[i, ] <- dat
    nc_close(nc)
  }
  datmat
}

## Function to extract the "modern value" direct from the file (last 20 years of the timeseries)
# for each of the variables in the file, average the last 20 years as the "modern" value,
# and then place in a matrix

modernValue <- function(nc, variable, ix){
  # A basic function to read a variable and 
  # take the mean of the timeseries at locations ix
  dat <- ncvar_get(nc, variable)
  out <- mean(dat[ix])
  out
}

twoStep_glmnet <- function(X, y, nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL, seed=NULL, trace=FALSE, maxit=100,
                           REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, popsize=100){
  # Use lasso to reduce input dimension of emulator before
  # building.
  control_list = list(trace=trace, maxit=maxit, REPORT=REPORT, factr=factr, pgtol=pgtol, pop.size=popsize)
  xvars = colnames(X)
  data = data.frame(response=y, x=X)
  colnames(data) <- c("response", xvars)
  nval = length(y)
  
  # fit a lasso by cross validation
  library(glmnet)
  fit_glmnet_cv = cv.glmnet(x=X,y=y)
  
  # The labels of the retained coefficients are here
  # (retains intercept at index zero)
  coef_i = (coef(fit_glmnet_cv, s = "lambda.1se"))@i
  labs = labels(coef(fit_glmnet_cv, s = "lambda.1se"))[[1]]
  labs = labs[-1] # remove intercept
  glmnet_retained = labs[coef_i]
  
  start_form = as.formula(paste("~ ", paste(glmnet_retained , collapse= "+")))
  m = km(start_form, design=X, response=y, nugget=nugget, parinit=parinit,
         nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
  
  return(list(x=X, y=y, nugget=nugget, nuggetEstim=nuggetEstim,
              noiseVar=noiseVar, emulator=m, seed=seed, coefs=m@covariance@range.val,
              trends=m@trend.coef, meanTerms=all.vars(start_form), fit_glmnet_cv=fit_glmnet_cv))
}


## ----------------------------------------------------------------------
## Load ensemble
## ----------------------------------------------------------------------

# primary carbon cycle outputs
npp_ens <- makeTimeseriesEnsemble(variable = "npp_nlim_lnd_sum") / (1e12/ysec)
nbp_ens <-  makeTimeseriesEnsemble(variable = "nbp_lnd_sum") / (1e12/ysec)
cSoil_ens <-  makeTimeseriesEnsemble(variable = "cSoil_lnd_sum") / 1e12
cVeg_ens <-  makeTimeseriesEnsemble(variable = "cVeg_lnd_sum") / 1e12

total_land_carbon_ens <- cSoil_ens + cVeg_ens

lai_lnd_mean_ens <- makeTimeseriesEnsemble(variable = "lai_lnd_mean")

# fluxes
rh_lnd_sum_ens <- makeTimeseriesEnsemble(variable = "rh_lnd_sum") / (1e12/ysec)
fLuc_lnd_sum_ens <- makeTimeseriesEnsemble(variable = "fLuc_lnd_sum") / (1e12/ysec)
fHarvest_lnd_sum_ens <- makeTimeseriesEnsemble(variable = "fHarvest_lnd_sum") / (1e12/ysec)


# fractions
treeFrac_lnd_mean_ens <- makeTimeseriesEnsemble(variable = "treeFrac_lnd_mean")
shrubFrac_lnd_mean_ens <- makeTimeseriesEnsemble(variable = "shrubFrac_lnd_mean")
baresoilFrac_lnd_mean_ens <- makeTimeseriesEnsemble(variable = "baresoilFrac_lnd_mean")
#c3PftFrac_lnd_mean_ens <- makeTimeseriesEnsemble(variable = "c3PftFrac_lnd_mean_ens")
#c4PftFrac_lnd_mean_ens <- makeTimeseriesEnsemble(variable = "c4PftFrac_lnd_mean_ens")


## ----------------------------------------------------------------------
## Anomalize ensemble
## ----------------------------------------------------------------------

npp_ens_anom <- anomalizeTSmatrix(npp_ens, 1:20)
nbp_ens_anom <- anomalizeTSmatrix(nbp_ens, 1:20)
cSoil_ens_anom <- anomalizeTSmatrix(cSoil_ens, 1:20)
cVeg_ens_anom <- anomalizeTSmatrix(cVeg_ens, 1:20)
total_land_carbon_anom <- anomalizeTSmatrix(total_land_carbon_ens, 1:20)



#144:164 is the 1993:2013
#modernValue(nc = nc, variable = "npp_nlim_lnd_mean", ix = 144:164)

# apply to the test file to check it works
#vec <- sapply(varlist, FUN = modernValue, nc = nc, ix = 144:164)


## --------------------------------------------------------------------------------------
## Loop to extract the "modern value" of a number of model outputs
## Generate ensemble numbers, mean of the last 20 years of the timeseries (1994-2013)
##
## --------------------------------------------------------------------------------------


if (file.exists("ensemble.rdata")) {
  load("ensemble.rdata")
} else {
  
  nens = 499
  datmat <- matrix(nrow = nens, ncol = length(varlist))
  colnames(datmat) <- varlist
  
  enslist <- paste("P", formatC(0:(nens-1), width=4, flag="0"), sep="")
  floc <- paste0(ensloc,ensmember,subdir)
  
  for(i in 1:nens){
    
    vec <- rep(NA, length(varlist))
    
    ensmember <- enslist[i] 
    
    fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    #print(fn)
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(varlist, FUN = modernValue, nc = nc, ix = 144:164))
    datmat[i, ] <- vec
    nc_close(nc)
  }
  
  save(nens, datmat,enslist,floc, file ="ensemble.rdata")
}


## --------------------------------------------------------------------------------------
## Calculate an ensemble of anomalies for all variables
## For each ensemble member and each variable, calculate the change from the 20 years at the start of the run, 
## to the twenty years at the end of the run.
## ---------------------------------------------------------------------------------------


tsAnomaly <- function(nc, variable, startix = 1:20, endix = 144:164){
  
  # A basic function to read a variable and calculate the anomaly at the end of the run
  dat <- ncvar_get(nc, variable)
  endMean <- mean(dat[endix])
  startMean <- mean(dat[startix])
  out <- endMean - startMean
  out
}

#tsAnomaly(nc = nc, variable = "npp_nlim_lnd_mean")

# Generate ensemble  mean of the last 20 years of the timeseries (1994-2013)

if (file.exists("anomaly_ensemble.rdata")) {
  load("anomaly_ensemble.rdata")
} else {
  
  nens = 499
  datmatAnom <- matrix(nrow = nens, ncol = length(varlist))
  colnames(datmatAnom) <- varlist
  
  enslist <- paste("P", formatC(0:(nens-1), width=4, flag="0"), sep="")
  floc <- paste0(ensloc,ensmember,subdir)
  
  for(i in 1:nens){
    
    vec <- rep(NA, length(varlist))
    
    ensmember <- enslist[i] 
    
    fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    #print(fn)
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(varlist, FUN = tsAnomaly, nc = nc))
    datmatAnom[i, ] <- vec
    nc_close(nc)
  }
  
  save(nens, datmatAnom, enslist,floc, file ="anomaly_ensemble.rdata")
}


## ---------------------------------------------------------------------------------------------
## Clean data sets to "level 0"
## Initial clean of data set, removing variables that don't change, and removing NAs (models that didn't run).
## ---------------------------------------------------------------------------------------------


Y_nlevel0_ix <- which(is.na(datmat[,'year']))

YAnom_nlevel0_ix <- which(is.na(datmatAnom[,'year']))

# This should be true to proceed, or we'll have to start excluding the combined set.
identical(Y_nlevel0_ix, YAnom_nlevel0_ix)

# Y is the whole data set
Y <- datmat
# Y_level0 is the 'cleaned' data set, truncated to variables that change, and removing NAs
Y_level0 <- datmat[-Y_nlevel0_ix, -c(2,30,31)]
Y_nlevel0 <- datmat[Y_nlevel0_ix, -c(2, 30, 31)]


# Y is the whole data set
YAnom <- datmatAnom
# Y.level0 is the 'cleaned' data set, truncated to variables that change, and removing NAs
YAnom_level0 <- datmatAnom[-YAnom_nlevel0_ix, -c(2,30,31)]
YAnom_nlevel0 <- datmatAnom[YAnom_nlevel0_ix, -c(2,30,31)]



## ---------------------------------------------------------------------------------------------
## load the original design and input space, normalize to [0-1]
##
## ---------------------------------------------------------------------------------------------

# Load up the data
lhs_i = read.table('~/brazilCSSP/code/brazil_cssp/analyze_u-ao732/data/ES_PPE_i/lhs_u-ao732.txt', header = TRUE)
lhs_ii = read.table('~/brazilCSSP/code/brazil_cssp/analyze_u-ao732/data/ES_PPE_ii/lhs_u-ao732a.txt', header = TRUE)

toplevel_ix = 1:499

# The raw input data is a latin hypercube
lhs = rbind(lhs_i, lhs_ii)[toplevel_ix, ]
lhs_level0 <- lhs[-Y_nlevel0_ix,]

X = normalize(lhs)
colnames(X) = colnames(lhs)

X_level0 <- X[-Y_nlevel0_ix,]
X_nlevel0 <- X[Y_nlevel0_ix,]

d = ncol(X)
# lower and higher bound on the normalised matrix for visualisation
rx = rbind(rep(0,32), rep(1,32))


## ---------------------------------------------------------------------------------------------
## Find Level 1 constraints (F0 < 0.9)
## Defined from Level 0 basis
##
## ---------------------------------------------------------------------------------------------

level1_ix <- which(X_level0[, 'f0_io'] < 0.9)

X_level1 <- X_level0[level1_ix, ]
Y_level1 <- Y_level0[level1_ix, ]

YAnom_level1 <- YAnom_level0[level1_ix, ]


## ---------------------------------------------------------------------------------------------
## Find Level 1 constraints (F0 < 0.9) & b_wl_io > 0.15
## Defined from Level 0 basis
##
## ---------------------------------------------------------------------------------------------
level1a_ix <- which(X_level0[, 'f0_io'] < 0.9 & X_level0[, 'b_wl_io'] > 0.15 )

X_level1a <- X_level0[level1a_ix, ]
Y_level1a <- Y_level0[level1a_ix,]

YAnom_level1a <- YAnom_level0[level1a_ix, ]






