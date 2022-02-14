# JULES-ES-1p0-common.R
# A script to underly JULES-ES-1p0 ensemble work. Loading data, helper functions etc.
# Loads ensemble and constrains ensemble to "level 1a", which excludes inputs where the model doesn't run, and removes 
# ensemble members above/below (normalised) thresholds in f0_io and b_wl_io.
# Doug McNeall July 2021


# Note: at the moment, makeTimeSeriesEnsemble opens a file and extracts a single timesries, which is really inefficient.  
# Refactor to extract all the necessary timeseries.

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
library(wesanderson)

library(foreach)

library(emtools)
library(imptools)
library(viztools)
library(julesR)

# move the right functions here
source('~/brazilCSSP/code/brazil_cssp/per_pft.R') # eventually, move the relevant functions
source('explore-JULES-ES-1p0_PPE_functions.R')
source('~/myRpackages/julesR/vignettes/default_jules_parameter_perturbations.R')


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

zissou5 <- wes_palette('Zissou1', 5, type = c('discrete', 'continuous'))
zblue <- makeTransparent(as.character(zissou5)[1], 150)
zred <- makeTransparent(as.character(zissou5)[5], 150)

ysec = 60*60*24*365
years <- 1850:2013


# We're just interested in the "sum" (global totals) data, not the "mean" (global means) data
y_names_sum <- c('nbp_lnd_sum', 'fLuc_lnd_sum', 'npp_nlim_lnd_sum', 'cSoil_lnd_sum',
                 'cVeg_lnd_sum', 'landCoverFrac_lnd_sum', 'fHarvest_lnd_sum',
                 'lai_lnd_sum', 'rh_lnd_sum', 'treeFrac_lnd_sum', 'c3PftFrac_lnd_sum', 
                 'c4PftFrac_lnd_sum', 'shrubFrac_lnd_sum', 'baresoilFrac_lnd_sum')

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

makeTimeseriesEnsemble <- function(ensloc, variable, nens = 499, nts = 164, cn = 1850:2013){
  
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

getStandardMember <- function(ensloc, variable, nts = 164, cn = 1850:2013){
  
  datmat <- matrix(NA, nrow = 1, ncol = nts)
  colnames(datmat) <- cn
  
  ensmember <- 'S3'
  fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
  
  try(nc <- nc_open(paste0(fn)))
  try(dat <- extractTimeseries(nc, variable))
  
  datmat[1, ] <- dat
  nc_close(nc)
  datmat
  
}


mat2list <- function(X){
  
  # Turns the p columns of a matrix into a p length list,
  # with each column becoming an element of the list
  
  out <- vector(mode = 'list', length = ncol(X))
  for(i in 1:ncol(X)){
    out[[i]] <- X[ , i]
  }
  out
  
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


if (file.exists("ensemble_timeseries.rdata")) {
  load("ensemble_timeseries.rdata")
} else {
  
  # primary carbon cycle outputs
  npp_ens <- makeTimeseriesEnsemble(ensloc = ensloc, variable = "npp_nlim_lnd_sum") / (1e12/ysec)
  nbp_ens <-  makeTimeseriesEnsemble(ensloc = ensloc,variable = "nbp_lnd_sum") / (1e12/ysec)
  cSoil_ens <-  makeTimeseriesEnsemble(ensloc = ensloc,variable = "cSoil_lnd_sum") / 1e12
  cVeg_ens <-  makeTimeseriesEnsemble(ensloc = ensloc,variable = "cVeg_lnd_sum") / 1e12
  
  
  lai_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc,variable = "lai_lnd_mean")
  
  # fluxes
  rh_lnd_sum_ens <- makeTimeseriesEnsemble(ensloc = ensloc, variable = "rh_lnd_sum") / (1e12/ysec)
  fLuc_lnd_sum_ens <- makeTimeseriesEnsemble(ensloc = ensloc, variable = "fLuc_lnd_sum") / (1e12/ysec)
  fHarvest_lnd_sum_ens <- makeTimeseriesEnsemble(ensloc = ensloc, variable = "fHarvest_lnd_sum") / (1e12/ysec)
  
  
  # fractions
  treeFrac_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc, variable = "treeFrac_lnd_mean")
  shrubFrac_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc, variable = "shrubFrac_lnd_mean")
  baresoilFrac_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc, variable = "baresoilFrac_lnd_mean")
  #c3PftFrac_lnd_mean_ens <- makeTimeseriesEnsemble(variable = "c3PftFrac_lnd_mean_ens")
  #c4PftFrac_lnd_mean_ens <- makeTimeseriesEnsemble(variable = "c4PftFrac_lnd_mean_ens")
  
  save(npp_ens, nbp_ens, cSoil_ens, cVeg_ens, lai_lnd_mean_ens, rh_lnd_sum_ens, fLuc_lnd_sum_ens, fHarvest_lnd_sum_ens,
       treeFrac_lnd_mean_ens, shrubFrac_lnd_mean_ens, baresoilFrac_lnd_mean_ens, file = "ensemble_timeseries.rdata" )
  
}

total_land_carbon_ens <- cSoil_ens + cVeg_ens

# ------------------------------------------------------------------------------
# Get standard members
#
# ------------------------------------------------------------------------------

npp_stan <- getStandardMember(ensloc = ensloc, variable = "npp_nlim_lnd_sum") / (1e12/ysec)
nbp_stan <- getStandardMember(ensloc = ensloc, variable = "nbp_lnd_sum") / (1e12/ysec)
cSoil_stan <- getStandardMember(ensloc = ensloc, variable = "cSoil_lnd_sum") / 1e12
cVeg_stan <- getStandardMember(ensloc = ensloc, variable = "cVeg_lnd_sum") / 1e12
lai_lnd_mean_stan <- getStandardMember(ensloc = ensloc, variable = "lai_lnd_mean")


# fluxes
rh_lnd_sum_stan <- getStandardMember(ensloc = ensloc, variable = "rh_lnd_sum") / (1e12/ysec)
fLuc_lnd_sum_stan <- getStandardMember(ensloc = ensloc, variable = "fLuc_lnd_sum") / (1e12/ysec)
fHarvest_lnd_sum_stan <- getStandardMember(ensloc = ensloc, variable = "fHarvest_lnd_sum") / (1e12/ysec)


# fractions
treeFrac_lnd_mean_stan <- getStandardMember(ensloc = ensloc, variable = "treeFrac_lnd_mean")
shrubFrac_lnd_mean_stan <- getStandardMember(ensloc = ensloc, variable = "shrubFrac_lnd_mean")
baresoilFrac_lnd_mean_stan <- getStandardMember(ensloc = ensloc, variable = "baresoilFrac_lnd_mean")

## ----------------------------------------------------------------------
## Anomalize ensemble
## ----------------------------------------------------------------------

npp_ens_anom <- anomalizeTSmatrix(npp_ens, 1:20)
nbp_ens_anom <- anomalizeTSmatrix(nbp_ens, 1:20)
cSoil_ens_anom <- anomalizeTSmatrix(cSoil_ens, 1:20)
cVeg_ens_anom <- anomalizeTSmatrix(cVeg_ens, 1:20)
total_land_carbon_anom <- anomalizeTSmatrix(total_land_carbon_ens, 1:20)


## --------------------------------------------------------------------------------------
## Loop to extract the "modern value" of a number of model outputs
## Generate ensemble numbers, mean of the last 20 years of the timeseries (1994-2013)
##
## --------------------------------------------------------------------------------------


if (file.exists("ensemble.rdata")) {
  load("ensemble.rdata")
} else {
  
  nens = 499
  datmat <- matrix(nrow = nens, ncol = length(y_names_sum))
  colnames(datmat) <- y_names_sum
  
  enslist <- paste("P", formatC(0:(nens-1), width=4, flag="0"), sep="")
  floc <- paste0(ensloc,ensmember,subdir)
  
  for(i in 1:nens){
    
    vec <- rep(NA, length(y_names_sum))
    
    ensmember <- enslist[i] 
    
    fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    #print(fn)
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(y_names_sum, FUN = modernValue, nc = nc, ix = 144:164))
    datmat[i, ] <- vec
    nc_close(nc)
  }
  
  save(nens, datmat,enslist,floc, file ="ensemble.rdata")
}


# -------------------------------------------------------------------------------
# Standard Member modern value
  
ensmember <- 'S3'
fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
try(nc <- nc_open(paste0(fn)))
try(standard_modern_value <- sapply(y_names_sum, FUN = modernValue, nc = nc, ix = 144:164))
  
  

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
  datmatAnom <- matrix(nrow = nens, ncol = length(y_names_sum))
  colnames(datmatAnom) <- y_names_sum
  
  enslist <- paste("P", formatC(0:(nens-1), width=4, flag="0"), sep="")
  floc <- paste0(ensloc,ensmember,subdir)
  
  for(i in 1:nens){
    
    vec <- rep(NA, length(y_names_sum))
    
    ensmember <- enslist[i] 
    
    fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    #print(fn)
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(y_names_sum, FUN = tsAnomaly, nc = nc))
    datmatAnom[i, ] <- vec
    nc_close(nc)
  }
  
  save(nens, datmatAnom, enslist,floc, file ="anomaly_ensemble.rdata")
}



#npp_stan_anom <- anomalizeTSmatrix(npp_stan, 1:20)

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


p <- ncol(X_level0)

## ---------------------------------------------------------------------------------------------
# Outputs used for constraining the model
## ---------------------------------------------------------------------------------------------
ynames_const <- c('nbp_lnd_sum', 'npp_nlim_lnd_sum', 'cSoil_lnd_sum', 'cVeg_lnd_sum')
yunits_const <- c('GtC/year', 'GtC/year', 'GtC', 'GtC')
Y_const_level1a <- Y_level1a[, ynames_const]
Y_const_stan <- standard_modern_value[ynames_const]


scalevec <- c(1e12/ysec, 1e12/ysec, 1e12, 1e12)
Y_const_level1a_scaled <- sweep(Y_const_level1a, 2, STATS = scalevec, FUN = '/' )
Y_const_stan_scaled <- Y_const_stan / scalevec


# This is a "normalisation vector", for making the output numbers more manageable.
#cs_gb       cv    gpp_gb        nbp npp_n_gb    runoff
norm_vec = c(1e12, 1e12, 1e12/ysec , 1e12, 1e12, 1e9)


## ---------------------------------------------------------------------------------------------
# km emulator lists for Y and YAnom
# Standard emulator, the same for each of the 'sum' outputs, and the corresponding
# anomaly. Uses a linear prior, normalised inputs.
## ---------------------------------------------------------------------------------------------


Y_sum_level1a <- Y_level1a[ , y_names_sum]
YAnom_sum_level1a <- YAnom_level1a[ , y_names_sum]


Y_sum_level1a_list <- mat2list(Y_sum_level1a)
YAnom_sum_level1a_list <- mat2list(YAnom_sum_level1a)

if (file.exists("emlist_km_Y_level1a.rdata")) {
  load("emlist_km_Y_level1a.rdata")
} else {
  
  # Here, the list is a list version of the matrix Y_
  emlist_km_Y_level1a <- mclapply(X = Y_sum_level1a_list, FUN = km, formula = ~., design = X_level1a, mc.cores = 4) 
  
  save(emlist_km_Y_level1a, file = "emlist_km_Y_level1a.rdata")
  
}


if (file.exists("emlist_km_YAnom_level1a.rdata")) {
  load("emlist_km_YAnom_level1a.rdata")
} else {
  
  
  emlist_km_YAnom_level1a <- mclapply(X = YAnom_sum_level1a_list, FUN = km, formula = ~., design = X_level1a, mc.cores = 4) 
  
  save(emlist_km_YAnom_level1a, file = "emlist_km_YAnom_level1a.rdata")
  
}







