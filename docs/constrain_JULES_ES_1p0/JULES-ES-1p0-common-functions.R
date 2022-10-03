#JULES-ES-1p0-common-functions.R

## (C) Crown copyright, Met Office

## ----------------------------------------------------------------------
## Helper functions
## ----------------------------------------------------------------------


anomalizeTS <- function(x, ix = 1:20){x - mean(x[ix]) } 

tsAnomaly <- function(nc, variable, startix = 1:20, endix = 144:164){
  
  # A basic function to read a variable and calculate the anomaly at the end of the run
  dat <- ncvar_get(nc, variable)
  endMean <- mean(dat[endix])
  startMean <- mean(dat[startix])
  out <- endMean - startMean
  out
}


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

# makeTimeseriesEnsemble <- function(ensloc, variable, nens = 499, nts = 164, cn = 1850:2013){
#   
#   # nens is number of ensemble members
#   # nts length of timeseries
#   # cn is colnames()
#   datmat <- matrix(NA, nrow = nens, ncol = nts)
#   colnames(datmat) <- cn
#   
#   enslist <- paste("P", formatC(0:(nens-1), width=4, flag="0"), sep="")
#   #floc <- paste0(ensloc,ensmember,subdir)
#   
#   for(i in 1:nens){
#     
#     vec <- rep(NA,nts)
#     
#     ensmember <- enslist[i] 
#     
#     #fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
#     fn <- paste0(ensloc,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
#     
#     
#     try(nc <- nc_open(paste0(fn)))
#     try(dat <- extractTimeseries(nc, variable))
#     
#     datmat[i, ] <- dat
#     nc_close(nc)
#   }
#   datmat
# }

makeTimeseriesEnsemble <- function(ensloc, variable, nstart, nend, cn = 1850:2013){
  
  ysec <- 31536000
  
  nens <- (nend - nstart) + 1
  # nens is number of ensemble members
  # nts length of timeseries
  # cn is colnames()
  datmat <- matrix(NA, nrow = nens, ncol = length(cn))
  colnames(datmat) <- cn
  
  enslist <- paste("P", formatC(nstart:nend, width=4, flag="0"), sep="")
  
  for(i in 1:nens){
    
    ensmember <- enslist[i]
    
    fn <- paste0(ensloc,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    
    try(nc <- nc_open(paste0(fn)))
    try(localtime <- ncvar_get(nc, 'time'))
    
    # This part compensates for the fact that sometimes years are missing
    try(localyear <- floor(2010 + (localtime / ysec)))
    try(ix <- which(cn%in%localyear))
    
    try(dat <- extractTimeseries(nc, variable))
    
    try(datmat[i, ix] <- dat)
    nc_close(nc)
  }
  datmat
}

getStandardMember <- function(ensloc, variable, nts = 164, cn = 1850:2013){
  
  datmat <- matrix(NA, nrow = 1, ncol = nts)
  colnames(datmat) <- cn
  
  ensmember <- 'S3'
  #fn <- paste0(ensloc,ensmember,'/stats/','JULES-ES-1p0_',ensmember,'_Annual_global.nc')
  fn <- paste0(ensloc,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
  
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

makeJulesEnsembleModernValue <- function(ensloc, varlist, nstart, nend, ix = 144:164){
  
  nens <- (nend - nstart) + 1
  datmat <- matrix(nrow = nens, ncol = length(varlist))
  colnames(datmat) <- varlist
  
  enslist <- paste("P", formatC(nstart:nend, width=4, flag="0"), sep="")
  
  for(i in 1:nens){
    
    vec <- rep(NA, length(varlist))
    
    ensmember <- enslist[i]
    
    fn <- paste0(ensloc,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(varlist, FUN = modernValue, nc = nc, ix = ix))
    datmat[i, ] <- vec
    nc_close(nc)
  }
  return(list(datmat = datmat, enslist = enslist))
}


anomalyValue <- function(nc, variable, ix){
  # A basic function to read a variable and 
  # take the mean of the timeseries at locations ix
  dat <- ncvar_get(nc, variable)
  out <- mean(dat[ix]) - mean(dat[1:20])
  out
}

makeJulesEnsembleAnomaly <- function(ensloc, varlist, nstart, nend, ix = 144:164){
  
  nens <- (nend - nstart) + 1
  datmat <- matrix(nrow = nens, ncol = length(varlist))
  colnames(datmat) <- varlist
  
  enslist <- paste("P", formatC(nstart:nend, width=4, flag="0"), sep="")
  
  for(i in 1:nens){
    
    vec <- rep(NA, length(varlist))
    
    ensmember <- enslist[i]
    
    fn <- paste0(ensloc,'JULES-ES-1p0_',ensmember,'_Annual_global.nc')
    
    try(nc <- nc_open(paste0(fn)))
    try(vec <- sapply(varlist, FUN = anomalyValue, nc = nc, ix = ix))
    datmat[i, ] <- vec
    nc_close(nc)
  }
  return(list(datmat = datmat, enslist = enslist))
  
  
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



multiPred <- function(Y, Xpred, fit_list ){
  
  pred_mean <- matrix(NA, ncol = ncol(Y), nrow = nrow(Xpred))
  pred_sd <- matrix(NA, ncol = ncol(Y), nrow = nrow(Xpred))
  
  colnames(pred_mean) <- colnames(Y)
  colnames(pred_sd) <- colnames(Y)
  
  for(i in 1:length(fit_list)){
    
    pred <- predict.km(object=fit_list[[i]], newdata = Xpred, type = 'UK')
    pred_mean[, i] <- pred$mean
    pred_sd[, i] <- pred$sd
    
  }
  
  return(list(pred_mean = pred_mean, pred_sd = pred_sd))
  
}

