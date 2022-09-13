#JULES-ES-1p0-common-data.R

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
## Load Wave00 timeseries
## --------------------------------------------------------------------------------------------------------

global_mean_timeseries_wave00_file <- "global_mean_timeseries_wave00_2022-09-13.rdata"

if (file.exists(global_mean_timeseries_wave00_file)) {
  load(global_mean_timeseries_wave00_file)
} else {
  
  # primary carbon cycle outputs
  npp_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00,nstart = 0, nend = 498, variable = "npp_nlim_lnd_sum") / (1e12/ysec)
  nbp_ens <-  makeTimeseriesEnsemble(ensloc = ensloc_wave00,nstart = 0, nend = 498,variable = "nbp_lnd_sum") / (1e12/ysec)
  cSoil_ens <-  makeTimeseriesEnsemble(ensloc = ensloc_wave00,nstart = 0, nend = 498,variable = "cSoil_lnd_sum") / 1e12
  cVeg_ens <-  makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "cVeg_lnd_sum") / 1e12
  
  
  lai_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498, variable = "lai_lnd_mean")
  
  # fluxes
  rh_lnd_sum_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498, variable = "rh_lnd_sum") / (1e12/ysec)
  fLuc_lnd_sum_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "fLuc_lnd_sum") / (1e12/ysec)
  fHarvest_lnd_sum_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498, variable = "fHarvest_lnd_sum") / (1e12/ysec)
  
  
  # fractions
  treeFrac_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498, variable = "treeFrac_lnd_mean")
  shrubFrac_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00,nstart = 0, nend = 498,variable = "shrubFrac_lnd_mean")
  baresoilFrac_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "baresoilFrac_lnd_mean")
  c3PftFrac_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "c3PftFrac_lnd_mean")
  c4PftFrac_lnd_mean_ens <- makeTimeseriesEnsemble(ensloc = ensloc_wave00, nstart = 0, nend = 498,variable = "c4PftFrac_lnd_mean")
  
  save(npp_ens, nbp_ens, cSoil_ens, cVeg_ens, lai_lnd_mean_ens, rh_lnd_sum_ens, fLuc_lnd_sum_ens, fHarvest_lnd_sum_ens,
       treeFrac_lnd_mean_ens, shrubFrac_lnd_mean_ens, baresoilFrac_lnd_mean_ens,c3PftFrac_lnd_mean_ens, c4PftFrac_lnd_mean_ens,
       file = global_mean_timeseries_wave00_file )
  
}

# ------------------------------------------------------------------------------
# Get standard members
#
# ------------------------------------------------------------------------------

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

## ----------------------------------------------------------------------
## Anomalize ensemble
## ----------------------------------------------------------------------

npp_ens_anom <- anomalizeTSmatrix(npp_ens, 1:20)
nbp_ens_anom <- anomalizeTSmatrix(nbp_ens, 1:20)
cSoil_ens_anom <- anomalizeTSmatrix(cSoil_ens, 1:20)
cVeg_ens_anom <- anomalizeTSmatrix(cVeg_ens, 1:20)

rh_lnd_sum_ens_anom <- anomalizeTSmatrix(rh_lnd_sum_ens, 1:20)
fLuc_lnd_sum_ens_anom <- anomalizeTSmatrix(fLuc_lnd_sum_ens, 1:20)
lai_lnd_mean_ens_anom <- anomalizeTSmatrix(lai_lnd_mean_ens, 1:20) 

fHarvest_lnd_sum_ens_anom <- anomalizeTSmatrix(fHarvest_lnd_sum_ens, 1:20)
treeFrac_lnd_mean_ens_anom <- anomalizeTSmatrix(treeFrac_lnd_mean_ens, 1:20)
shrubFrac_lnd_mean_ens_anom <- anomalizeTSmatrix(shrubFrac_lnd_mean_ens, 1:20)
baresoilFrac_lnd_mean_ens_anom <- anomalizeTSmatrix(baresoilFrac_lnd_mean_ens, 1:20)
c3PftFrac_lnd_mean_ens_anom <- anomalizeTSmatrix(c3PftFrac_lnd_mean_ens, 1:20)
c4PftFrac_lnd_mean_ens_anom <- anomalizeTSmatrix(c4PftFrac_lnd_mean_ens, 1:20)

## --------------------------------------------------------------------------------------
## Loop to extract the "modern value" of a number of model outputs
## Generate ensemble numbers, mean of the last 20 years of the timeseries (1994-2013)
##
## --------------------------------------------------------------------------------------


if (file.exists("ensemble_2022-04-08.rdata")) {
  load("ensemble_2022-04-08.rdata")
} else {
  
  nens = 499
  datmat <- matrix(nrow = nens, ncol = length(y_names_all))
  colnames(datmat) <- y_names_all
  
  enslist <- paste("P", formatC(0:(nens-1), width=4, flag="0"), sep="")
  
  for(i in 1:nens){
    
    vec <- rep(NA, length(y_names_all))
    
    ensmember <- enslist[i] 
    
    fn <- paste0(ensloc_wave00,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(y_names_all, FUN = modernValue, nc = nc, ix = 144:164))
    datmat[i, ] <- vec
    nc_close(nc)
  }
  
  save(nens, datmat,enslist, file ="ensemble_2023-04-08.rdata")
}


# -------------------------------------------------------------------------------
# Standard Member modern value

ensmember <- 'S3'
fn <- paste0(ensloc_wave00,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
try(nc <- nc_open(paste0(fn)))
try(standard_modern_value <- sapply(y_names_all, FUN = modernValue, nc = nc, ix = 144:164))
# -------------------------------------------------------------------------------


