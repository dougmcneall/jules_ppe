#JULES-ES-1p0-common-data.R


# (C) Crown copyright, Met Office
# default jules parameters and perturbation limits in the ensemble (all PFTs).
source('default_jules_parameter_perturbations.R')

## ----------------------------------------------------------------------
## Data locations and constants
## ----------------------------------------------------------------------
#ensloc <- '/project/carbon_ppe/JULES-ES-1p0_PPE/'
ensloc_wave00 <- '/data/users/hadaw/JULES_ES_PPE/u-au932/'
ensloc_wave01 <- '/data/users/hadaw/JULES_ES_PPE/u-ck006/'

 
# Some pallete options
yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(11, "RdYlBu")
byr = rev(ryb)
rb = brewer.pal(11, "RdBu")
br = rev(rb)
blues = brewer.pal(9, 'Blues')
cbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

wave00col <- 'skyblue2'
wave01col <- 'tomato2'

#wave00col <- 'dodgerblue2'
#wave01col <- 'firebrick'
#rangecol <- 'grey'


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

y_names_all <-  c("nbp_lnd_sum", "year" ,"nbp_lnd_mean", "fLuc_lnd_sum", "fLuc_lnd_mean" , "npp_nlim_lnd_sum",  
                  "npp_nlim_lnd_mean" , "cSoil_lnd_sum" ,"cSoil_lnd_mean" ,
                  "cVeg_lnd_sum"  ,"cVeg_lnd_mean"  ,"landCoverFrac_lnd_sum" ,
                  "landCoverFrac_lnd_mean","fHarvest_lnd_sum","fHarvest_lnd_mean", 
                  "lai_lnd_sum" ,"lai_lnd_mean" ,"rh_lnd_sum" ,
                  "rh_lnd_mean" ,"treeFrac_lnd_sum", "treeFrac_lnd_mean" ,
                  "c3PftFrac_lnd_sum" ,"c3PftFrac_lnd_mean"  , "c4PftFrac_lnd_sum" ,
                  "c4PftFrac_lnd_mean" ,"shrubFrac_lnd_sum"  ,   "shrubFrac_lnd_mean",
                  "baresoilFrac_lnd_sum" ,"baresoilFrac_lnd_mean","residualFrac_lnd_sum" ,
                  "residualFrac_lnd_mean")

y_names_select <-  c("npp_nlim_lnd_sum", "nbp_lnd_sum", "cSoil_lnd_sum", "cVeg_lnd_sum",
                     "lai_lnd_mean",
                     "rh_lnd_sum" , "fLuc_lnd_sum", "fHarvest_lnd_sum",  
                     "landCoverFrac_lnd_mean", 
                     "treeFrac_lnd_mean" , "baresoilFrac_lnd_mean",
                     "shrubFrac_lnd_mean", "c3PftFrac_lnd_mean",
                     "c4PftFrac_lnd_mean"   
)

select_units <- c('GtC/year', 'GtC/year', 'GtC', 'GtC', 'index', 'GtC/year','GtC/year', 'GtC/year', '%', '%', '%', '%', '%', '%')
names(select_units) <- y_names_select

## --------------------------------------------------------------------------------------------------------
## Load Wave00 global means timeseries
## --------------------------------------------------------------------------------------------------------

global_mean_timeseries_wave00_file <- "data/global_mean_timeseries_wave00_2022-09-13.rdata"

if (file.exists(global_mean_timeseries_wave00_file)) {
  load(global_mean_timeseries_wave00_file)
} else {
  
  # primary carbon cycle outputs
  npp_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00,nstart = 0, nend = 498, variable = "npp_nlim_lnd_sum") / (1e12/ysec)
  nbp_ens_wave00 <-  makeTimeseriesEnsemble(ensloc = ensloc_wave00,nstart = 0, nend = 498,variable = "nbp_lnd_sum") / (1e12/ysec)
  cSoil_ens_wave00 <-  makeTimeseriesEnsemble(ensloc = ensloc_wave00,nstart = 0, nend = 498,variable = "cSoil_lnd_sum") / 1e12
  cVeg_ens_wave00 <-  makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "cVeg_lnd_sum") / 1e12
  
  
  lai_lnd_mean_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498, variable = "lai_lnd_mean")
  
  # fluxes
  rh_lnd_sum_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498, variable = "rh_lnd_sum") / (1e12/ysec)
  fLuc_lnd_sum_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "fLuc_lnd_sum") / (1e12/ysec)
  fHarvest_lnd_sum_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498, variable = "fHarvest_lnd_sum") / (1e12/ysec)
  
  
  # fractions
  treeFrac_lnd_mean_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498, variable = "treeFrac_lnd_mean")
  shrubFrac_lnd_mean_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00,nstart = 0, nend = 498,variable = "shrubFrac_lnd_mean")
  baresoilFrac_lnd_mean_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "baresoilFrac_lnd_mean")
  c3PftFrac_lnd_mean_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "c3PftFrac_lnd_mean")
  c4PftFrac_lnd_mean_ens_wave00 <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "c4PftFrac_lnd_mean")
  
  save(npp_ens_wave00, nbp_ens_wave00, cSoil_ens_wave00, 
       cVeg_ens_wave00, lai_lnd_mean_ens_wave00, rh_lnd_sum_ens_wave00,
       fLuc_lnd_sum_ens_wave00, fHarvest_lnd_sum_ens_wave00,
       treeFrac_lnd_mean_ens_wave00, shrubFrac_lnd_mean_ens_wave00, 
       baresoilFrac_lnd_mean_ens_wave00, c3PftFrac_lnd_mean_ens_wave00, c4PftFrac_lnd_mean_ens_wave00,
       file = global_mean_timeseries_wave00_file)
  
}

# ------------------------------------------------------------------------------
# Get standard members
#
# ------------------------------------------------------------------------------
global_mean_timeseries_stan_file <- "data/global_mean_timeseries_stan_2022-09-13.rdata"

if (file.exists(global_mean_timeseries_stan_file)) {
  load(global_mean_timeseries_stan_file)
} else {

npp_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "npp_nlim_lnd_sum") / (1e12/ysec)
nbp_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "nbp_lnd_sum") / (1e12/ysec)
cSoil_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "cSoil_lnd_sum") / 1e12
cVeg_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "cVeg_lnd_sum") / 1e12
lai_lnd_mean_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "lai_lnd_mean")


# fluxes
rh_lnd_sum_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "rh_lnd_sum") / (1e12/ysec)
fLuc_lnd_sum_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "fLuc_lnd_sum") / (1e12/ysec)
fHarvest_lnd_sum_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "fHarvest_lnd_sum") / (1e12/ysec)


# fractions
treeFrac_lnd_mean_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "treeFrac_lnd_mean")
shrubFrac_lnd_mean_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "shrubFrac_lnd_mean")
baresoilFrac_lnd_mean_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "baresoilFrac_lnd_mean")
c3PftFrac_lnd_mean_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "c3PftFrac_lnd_mean")
c4PftFrac_lnd_mean_stan <- getStandardMember(ensloc = ensloc_wave00, variable = "c4PftFrac_lnd_mean")

save(npp_stan, nbp_stan, cSoil_stan, cVeg_stan, lai_lnd_mean_stan,
     rh_lnd_sum_stan, fLuc_lnd_sum_stan, fHarvest_lnd_sum_stan, 
     treeFrac_lnd_mean_stan, shrubFrac_lnd_mean_stan, baresoilFrac_lnd_mean_stan,
     c3PftFrac_lnd_mean_stan,c4PftFrac_lnd_mean_stan, 
     file = global_mean_timeseries_stan_file)
}

## ----------------------------------------------------------------------
## Anomalize ensemble
## ----------------------------------------------------------------------

npp_ens_anom_wave00 <- anomalizeTSmatrix(npp_ens_wave00, 1:20)
nbp_ens_anom_wave00 <- anomalizeTSmatrix(nbp_ens_wave00, 1:20)
cSoil_ens_anom_wave00 <- anomalizeTSmatrix(cSoil_ens_wave00, 1:20)
cVeg_ens_anom_wave00 <- anomalizeTSmatrix(cVeg_ens_wave00, 1:20)

rh_lnd_sum_ens_anom_wave00 <- anomalizeTSmatrix(rh_lnd_sum_ens_wave00, 1:20)
fLuc_lnd_sum_ens_anom_wave00 <- anomalizeTSmatrix(fLuc_lnd_sum_ens_wave00, 1:20)
lai_lnd_mean_ens_anom_wave00 <- anomalizeTSmatrix(lai_lnd_mean_ens_wave00, 1:20) 

fHarvest_lnd_sum_ens_anom_wave00 <- anomalizeTSmatrix(fHarvest_lnd_sum_ens_wave00, 1:20)
treeFrac_lnd_mean_ens_anom_wave00 <- anomalizeTSmatrix(treeFrac_lnd_mean_ens_wave00, 1:20)
shrubFrac_lnd_mean_ens_anom_wave00 <- anomalizeTSmatrix(shrubFrac_lnd_mean_ens_wave00, 1:20)
baresoilFrac_lnd_mean_ens_anom_wave00 <- anomalizeTSmatrix(baresoilFrac_lnd_mean_ens_wave00, 1:20)
c3PftFrac_lnd_mean_ens_anom_wave00 <- anomalizeTSmatrix(c3PftFrac_lnd_mean_ens_wave00, 1:20)
c4PftFrac_lnd_mean_ens_anom_wave00 <- anomalizeTSmatrix(c4PftFrac_lnd_mean_ens_wave00, 1:20)


## --------------------------------------------------------------------------------------
## Loop to extract the "modern value" of a number of model outputs
## Generate ensemble numbers, mean of the last 20 years of the timeseries (1994-2013)
##
## --------------------------------------------------------------------------------------

global_mean_modern_value_wave00_file <- "data/global_mean_modern_value_wave00_2022-13-09.rdata"

if (file.exists(global_mean_modern_value_wave00_file)) {
  load(global_mean_modern_value_wave00_file)
} else {
  
  nens_wave00 = 499
  datmat_wave00 <- matrix(nrow = nens_wave00, ncol = length(y_names_all))
  colnames(datmat_wave00) <- y_names_all
  
  enslist_wave00 <- paste("P", formatC(0:(nens_wave00 - 1), width=4, flag="0"), sep="")
  
  for(i in 1:nens_wave00){
    
    vec <- rep(NA, length(y_names_all))
    
    ensmember <- enslist_wave00[i] 
    
    fn <- paste0(ensloc_wave00,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(y_names_all, FUN = modernValue, nc = nc, ix = 144:164))
    datmat_wave00[i, ] <- vec
    nc_close(nc)
  }
  
  save(nens_wave00, datmat_wave00, enslist_wave00, file = global_mean_modern_value_wave00_file)
}


# -------------------------------------------------------------------------------
# Standard Member modern value
modern_value_stan_file <- "data/modern_value_stan_2022-09-13.rdata"

if (file.exists(modern_value_stan_file)) {
  load(modern_value_stan_file)
} else {
  
  ensmember <- 'S3'
  fn <- paste0(ensloc_wave00,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
  try(nc <- nc_open(paste0(fn)))
  try(standard_modern_value <- sapply(y_names_all, FUN = modernValue, nc = nc, ix = 144:164))
  
  
  save(standard_modern_value, file = modern_value_stan_file)
  
}
# -------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------
## Calculate an ensemble of anomalies for all variables
## For each ensemble member and each variable, calculate the change from the 20 years at the start of the run, 
## to the twenty years at the end of the run.
## ---------------------------------------------------------------------------------------


# Generate ensemble  mean of the last 20 years of the timeseries (1994-2013)

global_mean_modern_value_anomaly_wave00_file <- "data/global_mean_modern_value_anomaly_wave00_2022-13-09.rdata"

if (file.exists(global_mean_modern_value_anomaly_wave00_file)) {
  load(global_mean_modern_value_anomaly_wave00_file)
} else {
  
  datmatAnom_wave00 <- matrix(nrow = nens_wave00, ncol = length(y_names_all))
  colnames(datmatAnom_wave00) <- y_names_all
  
  enslist <- paste("P", formatC(0:(nens_wave00 - 1), width=4, flag="0"), sep="")
  
  for(i in 1:nens_wave00){
    
    vec <- rep(NA, length(y_names_all))
    
    ensmember <- enslist_wave00[i] 
    
    fn <- paste0(ensloc_wave00,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(y_names_all, FUN = tsAnomaly, nc = nc))
    datmatAnom_wave00[i, ] <- vec
    nc_close(nc)
  }
  
  save(nens_wave00, datmatAnom_wave00, enslist_wave00, file = global_mean_modern_value_anomaly_wave00_file)
}



## ---------------------------------------------------------------------------------------------
## Clean data sets to "level 0"
## Initial clean of data set, removing variables that don't change, and removing NAs (models that didn't run).


## ---------------------------------------------------------------------------------------------

# this is wrong now, for some reason
Y_nlevel0_ix <- which(is.na(datmat_wave00[,'year']))

YAnom_nlevel0_ix <- which(is.na(datmatAnom_wave00[,'year']))

# This should be true to proceed, or we'll have to start excluding the combined set.
identical(Y_nlevel0_ix, YAnom_nlevel0_ix)

# Y is the whole data set
Y <- datmat_wave00
# Y_level0 is the 'cleaned' data set removing NAs
Y_level0 <- datmat_wave00[-Y_nlevel0_ix, ]
Y_nlevel0 <- datmat_wave00[Y_nlevel0_ix, ]


# Y is the whole data set
YAnom <- datmatAnom_wave00

# Y.level0 is the 'cleaned' data set  removing NAs
YAnom_level0 <- datmatAnom_wave00[-YAnom_nlevel0_ix, ]
YAnom_nlevel0 <- datmatAnom_wave00[YAnom_nlevel0_ix, ]

## ---------------------------------------------------------------------------------------------
## load the original design and input space, normalize to [0-1]
##
## ---------------------------------------------------------------------------------------------

# Load up the data
lhs_i = read.table('data/lhs_u-ao732.txt', header = TRUE)
lhs_ii = read.table('data/lhs_u-ao732a.txt', header = TRUE)

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


# ---------------------------------------------------------------------------------------------
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


emlist_km_Y_level1a_file <- "../emlist_km_Y_level1a_2022-04-08.rdata"

if (file.exists(emlist_km_Y_level1a_file )) {
  load(emlist_km_Y_level1a_file )
} else {
  
  # Here, the list is a list version of the matrix Y_
  emlist_km_Y_level1a <- mclapply(X = Y_sum_level1a_list, FUN = km, formula = ~., design = X_level1a, mc.cores = 4) 
  
  save(emlist_km_Y_level1a, file = emlist_km_Y_level1a_file )
  
}

emlist_km_YAnom_level1a_file <- "../emlist_km_YAnom_level1a_2022-04-08.rdata"


if (file.exists(emlist_km_YAnom_level1a_file )) {
  load(emlist_km_YAnom_level1a_file )
} else {
  
  
  emlist_km_YAnom_level1a <- mclapply(X = YAnom_sum_level1a_list, FUN = km, formula = ~., design = X_level1a, mc.cores = 4) 
  
  save(emlist_km_YAnom_level1a, file = emlist_km_YAnom_level1a_file )
  
}

## ------------------------------------------------------------------------------------
## Wave01 (second wave) specific stuff
## ------------------------------------------------------------------------------------ 

# Number of ensemble members (out of 500) to use for training in wave01
ntrain_wave01 <- 400


# Modern value JULES ensemble Wave01
nstart <- 499
nend <- (nstart + ntrain_wave01) - 1


ensemble_wave01_file <- "data/ensemble_wave01_2022-04-08.rdata"
  
if (file.exists(ensemble_wave01_file)) {
  load(ensemble_wave01_file)
} else {
  
  ens_wave01_mv <- makeJulesEnsembleModernValue(ensloc = ensloc_wave01, 
                                                varlist = y_names_sum,
                                                nstart = nstart,
                                                nend = nend, 
                                                ix = 144:164) 
  
  save(ens_wave01_mv, file = ensemble_wave01_file)
}

ensemble_wave01_anom_file <- "data/ensemble_wave01_2022-05-25.rdata"

if (file.exists(ensemble_wave01_anom_file )) {
  load(ensemble_wave01_anom_file )
} else {
  
  ens_wave01_anom <- makeJulesEnsembleAnomaly(ensloc = ensloc_wave01, 
                                              varlist = y_names_sum,
                                              nstart = nstart,
                                              nend = nend, 
                                              ix = 144:164) 
  
  save(ens_wave01_anom, file = ensemble_wave01_anom_file )
}


# the "select" wave01 ensemble, to tie in with constraints later
ens_select_wave01_mv_file <- "data/ens_select_wave01_mv_file_2022-09-26.rdata"

if (file.exists(ens_select_wave01_mv_file)) {
  load(ens_select_wave01_mv_file)
} else {
  
  ens_select_wave01_mv <- makeJulesEnsembleModernValue(ensloc = ensloc_wave01, 
                                                varlist = y_names_select,
                                                nstart = nstart,
                                                nend = nend, 
                                                ix = 144:164) 
  
  save(ens_select_wave01_mv, file = ens_select_wave01_mv_file)
}

ens_select_wave01_mv_anom_file <- "data/ens_select_wave01_mv_anom_file_2022-09-26.rdata"

if (file.exists(ens_select_wave01_mv_anom_file)) {
  load(ens_select_wave01_mv_anom_file )
} else {
  
  ens_select_wave01_mv_anom <- makeJulesEnsembleAnomaly(ensloc = ensloc_wave01, 
                                              varlist = y_names_select,
                                              nstart = nstart,
                                              nend = nend, 
                                              ix = 144:164) 
  
  save(ens_select_wave01_mv_anom, file = ens_select_wave01_mv_anom_file )
}



# Load input matrices and bind with wave00 inputs
lhs_wave01 <- read.table( 'data/lhs_example.txt', header = TRUE)

X_wave01 = normalize(lhs_wave01, wrt = rbind(lhs_i, lhs_ii, lhs_wave01))
colnames(X_wave01) = colnames(lhs_wave01)

# Match the 400 outputs we're using in the training data
X_wave01_train <- X_wave01[1:ntrain_wave01, ]


# Modern values that we use for constraints

Y_const_wave01 <- ens_wave01_mv$datmat[, ynames_const]
Y_const_wave01_scaled <- sweep(Y_const_wave01, 2, STATS = scalevec, FUN = '/' )

## -----------------------------------------------------------------------------
## Timeseries wave01
##
## 

## -----------------------------------------------------------------------------

ensemble_timeseries_wave01_file <- "data/ensemble_timeseries_wave01_2022-09-14.rdata"

if (file.exists(ensemble_timeseries_wave01_file)) {
  load(ensemble_timeseries_wave01_file)
} else {
  
  # primary carbon cycle outputs
  npp_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend, variable = "npp_nlim_lnd_sum") / (1e12/ysec)
  nbp_ens_wave01 <-  makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend,variable = "nbp_lnd_sum") / (1e12/ysec)
  cSoil_ens_wave01 <-  makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend,variable = "cSoil_lnd_sum") / 1e12
  cVeg_ens_wave01 <-  makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend,variable = "cVeg_lnd_sum") / 1e12
  # 
  # 
  lai_lnd_mean_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend,variable = "lai_lnd_mean")
  # 
  # # fluxes
  rh_lnd_sum_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend, variable = "rh_lnd_sum") / (1e12/ysec)
  fLuc_lnd_sum_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend, variable = "fLuc_lnd_sum") / (1e12/ysec)
  fHarvest_lnd_sum_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend, variable = "fHarvest_lnd_sum") / (1e12/ysec)
  # 
  # 
  # # fractions
  treeFrac_lnd_mean_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend,  variable = "treeFrac_lnd_mean")
  shrubFrac_lnd_mean_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend,  variable = "shrubFrac_lnd_mean")
  baresoilFrac_lnd_mean_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01,nstart = nstart, nend = nend,  variable = "baresoilFrac_lnd_mean")
  
  c3PftFrac_lnd_mean_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01, nstart = nstart, nend = nend,variable = "c3PftFrac_lnd_mean")
  c4PftFrac_lnd_mean_ens_wave01 <- makeTimeseriesEnsemble(ensloc = ensloc_wave01, nstart = nstart, nend = nend,variable = "c4PftFrac_lnd_mean")
  
  save(npp_ens_wave01,
       nbp_ens_wave01,
       cSoil_ens_wave01,
       cVeg_ens_wave01,
       lai_lnd_mean_ens_wave01,
       rh_lnd_sum_ens_wave01,
       fLuc_lnd_sum_ens_wave01,
       fHarvest_lnd_sum_ens_wave01,
       treeFrac_lnd_mean_ens_wave01,
       shrubFrac_lnd_mean_ens_wave01,
       baresoilFrac_lnd_mean_ens_wave01,
       c3PftFrac_lnd_mean_ens_wave01,
       c4PftFrac_lnd_mean_ens_wave01,
       file = ensemble_timeseries_wave01_file )
}

lhs_wave0_wave01_all <- rbind(lhs, lhs_wave01)



## Anomaly timeseries wave01

npp_ens_anom_wave01 <- anomalizeTSmatrix(npp_ens_wave01, 1:20)
nbp_ens_anom_wave01 <- anomalizeTSmatrix(nbp_ens_wave01, 1:20)
cSoil_ens_anom_wave01 <- anomalizeTSmatrix(cSoil_ens_wave01, 1:20)
cVeg_ens_anom_wave01 <- anomalizeTSmatrix(cVeg_ens_wave01, 1:20)

rh_lnd_sum_ens_anom_wave01 <- anomalizeTSmatrix(rh_lnd_sum_ens_wave01, 1:20)
fLuc_lnd_sum_ens_anom_wave01 <- anomalizeTSmatrix(fLuc_lnd_sum_ens_wave01, 1:20)
lai_lnd_mean_ens_anom_wave01 <- anomalizeTSmatrix(lai_lnd_mean_ens_wave01, 1:20) 

fHarvest_lnd_sum_ens_anom_wave01 <- anomalizeTSmatrix(fHarvest_lnd_sum_ens_wave01, 1:20)
treeFrac_lnd_mean_ens_anom_wave01 <- anomalizeTSmatrix(treeFrac_lnd_mean_ens_wave01, 1:20)
shrubFrac_lnd_mean_ens_anom_wave01 <- anomalizeTSmatrix(shrubFrac_lnd_mean_ens_wave01, 1:20)
baresoilFrac_lnd_mean_ens_anom_wave01 <- anomalizeTSmatrix(baresoilFrac_lnd_mean_ens_wave01, 1:20)
c3PftFrac_lnd_mean_ens_anom_wave01 <- anomalizeTSmatrix(c3PftFrac_lnd_mean_ens_wave01, 1:20)
c4PftFrac_lnd_mean_ens_anom_wave01 <- anomalizeTSmatrix(c4PftFrac_lnd_mean_ens_wave01, 1:20)



npp_stan_anom <- anomalizeTS(npp_stan)
nbp_stan_anom <- anomalizeTS(nbp_stan)
cSoil_stan_anom <-anomalizeTS(cSoil_stan)
cVeg_stan_anom <- anomalizeTS(cVeg_stan)
lai_lnd_mean_stan_anom <- anomalizeTS(lai_lnd_mean_stan)


# fluxes
rh_lnd_sum_stan_anom <- anomalizeTS(rh_lnd_sum_stan)
fLuc_lnd_sum_stan_anom <- anomalizeTS(fLuc_lnd_sum_stan)
fHarvest_lnd_sum_stan_anom <- anomalizeTS(fHarvest_lnd_sum_stan)


# fractions
treeFrac_lnd_mean_stan_anom <- anomalizeTS(treeFrac_lnd_mean_stan)
shrubFrac_lnd_mean_stan_anom <- anomalizeTS(shrubFrac_lnd_mean_stan)
baresoilFrac_lnd_mean_stan_anom <- anomalizeTS(baresoilFrac_lnd_mean_stan)
c3PftFrac_lnd_mean_stan_anom <- anomalizeTS(c3PftFrac_lnd_mean_stan)
c4PftFrac_lnd_mean_stan_anom <- anomalizeTS(c4PftFrac_lnd_mean_stan)


## -----------------------------------------------------------------------------------------------------
## Fix outliers in wave01
## find timeseries outliers
## 
## -----------------------------------------------------------------------------------------------------

# Timeseries that have problems. NBP, RH and cSoil seems to have large outliers

# These indices reference the separate ensembles
# cSoil over 6000
# rh over 200
# nbp less than -15
cSoil_outlier_ix_wave00 <- unique(which(cSoil_ens_wave00 > 6000, arr.ind = TRUE)[,'row'])
cSoil_outlier_ix_wave01 <- unique(which(cSoil_ens_wave01 > 6000, arr.ind = TRUE)[,'row'])

nbp_outlier_ix_wave00 <- unique(which(nbp_ens_wave00 < -15, arr.ind = TRUE)[,'row'])
nbp_outlier_ix_wave01 <- unique(which(nbp_ens_wave01 < -15, arr.ind = TRUE)[,'row'])

rh_lnd_sum_outlier_ix_wave00 <- unique(which(rh_lnd_sum_ens_wave00 > 200, arr.ind = TRUE)[,'row'])
rh_lnd_sum_outlier_ix_wave01 <- unique(which(rh_lnd_sum_ens_wave01 > 200, arr.ind = TRUE)[,'row'])


# are there additional excluded indices to those already excluded by the constraint

wave01_all_ix <- 1:ntrain_wave01


AW_constraints <- matrix(nrow = 2, ncol = length(ynames_const))

AW_constraints[1,] <- c(0, 35, 750, 300)
AW_constraints[2,] <- c(100, 80, 3000, 800)

colnames(AW_constraints) <- ynames_const
rownames(AW_constraints) <- c('min', 'max')


# conform to Andy's basic constraints
#level2_ix_wave01 <- which(apply(Y_const_wave01_scaled, 1, FUN = withinRange, maxes  = AW_constraints[2,], mins = AW_constraints[1,] ))

#nlevel2_ix_wave01 <- which(apply(Y_const_wave01_scaled, 1, FUN = withinRange, maxes  = AW_constraints[2,], mins = AW_constraints[1,] ) == FALSE)

level2_ix_wave01 <- which(Y_const_wave01_scaled[,'nbp_lnd_sum'] > 0 &
                            Y_const_wave01_scaled[,'npp_nlim_lnd_sum'] > 35 & Y_const_wave01_scaled[,'npp_nlim_lnd_sum'] < 80 &
                            Y_const_wave01_scaled[,'cSoil_lnd_sum'] > 750 & Y_const_wave01_scaled[,'cSoil_lnd_sum'] < 3000 &
                            Y_const_wave01_scaled[,'cVeg_lnd_sum'] > 300 & Y_const_wave01_scaled[,'cVeg_lnd_sum'] < 800
)


# Indices excluded in wave01 level2
level2_nix_wave01 <- setdiff(wave01_all_ix, level2_ix_wave01)

# would be interesting to see if these look normal in other ways
ts_outliers_ix_wave01 <- unique(c(cSoil_outlier_ix_wave01,nbp_outlier_ix_wave01, rh_lnd_sum_outlier_ix_wave01))
# are there any that are not excluded by level 2? (I assume so)
intersect(ts_outliers_ix_wave01, level2_nix_wave01)

without_outliers_ix_wave01 <- setdiff(wave01_all_ix,ts_outliers_ix_wave01)

# Remove these from the wave01 ensemble to remove outliers and excluded ensemble members
level2_and_ts_outliers_nix_wave01 <- union(level2_nix_wave01, ts_outliers_ix_wave01)

level2a_ix_wave01 <- setdiff(wave01_all_ix, level2_and_ts_outliers_nix_wave01)

wave00_all_ix <- 1:499
ts_outliers_ix_wave00 <- unique(c(cSoil_outlier_ix_wave00,nbp_outlier_ix_wave00, rh_lnd_sum_outlier_ix_wave00))
without_outliers_ix_wave00 <- setdiff(wave00_all_ix,ts_outliers_ix_wave00)

# Build "clean" complete dataset that conforms to level 1a (all training runs)
X_level1a_wave01 <- rbind(X_level1a, X_wave01_train[without_outliers_ix_wave01, ])
Y_const_level1a_wave01_scaled <- rbind(Y_const_level1a_scaled, Y_const_wave01_scaled[without_outliers_ix_wave01, ])

Y_sum_level1a_wave01 <- rbind(Y_sum_level1a, ens_wave01_mv$datmat[without_outliers_ix_wave01, ])

YAnom_sum_level1a_wave01 <- rbind(YAnom_sum_level1a, ens_wave01_anom$datmat[without_outliers_ix_wave01, ])


