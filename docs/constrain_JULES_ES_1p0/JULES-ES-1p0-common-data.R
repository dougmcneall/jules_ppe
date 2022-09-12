#JULES-ES-1p0-common-data.R

source('~/myRpackages/julesR/vignettes/default_jules_parameter_perturbations.R')

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
