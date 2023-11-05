#################
## library to store processing functions
####################
library(parallel)
library(forecast)
library(TSA)

## the function to implement IDAR algorithm
iterative_decorr_ts_autoAR = function(node_data, 
                                      TR=2, 
                                      ic = 'aicc', # information criteria to select AR order
                                      max.iter = 5, # max number of iterations to fit AR models
                                      max.cond.num = 1e8, # max estimated correlation matrix condition number
                                      save.cov = F # if save estimated covariance matrices for each whitening iterations
){
  
  continue = T
  detect_ar = 0 ## fitted AR order
  detect_ma = 0 ## fitted MA order
  cum_use_ar = 0 ## free ar parameters to estimate
  cum_use_ma = 0 ## free ma parameters to estimate
  processed_data = node_data
  steps = 0
  store_ar_ma = list()
  store_cov = list()
  whiten_cov_sqrt_inv = diag(1, length(node_data)) ## store the whitening matrix, to be directly multiplied to the data
  
  while(continue & steps < max.iter){
    # cat(steps)
    new_processed_data_res = decorr_ts_autoAR(node_data= processed_data, 
                                              TR = TR, 
                                              ic = ic, 
                                              max.cond.num = max.cond.num,
                                              current_iter = steps+1)
    
    if(steps >0  & all(c(new_processed_data_res$store$detect_ar, new_processed_data_res$store$detect_ma)==0)){
      continue = F
    }else{
      cum_use_ar = cum_use_ar + new_processed_data_res$store$use_ar
      cum_use_ma = cum_use_ma + new_processed_data_res$store$use_ma
      processed_data = new_processed_data_res$decorr_resid
      steps = steps + 1
      
      whiten_cov_sqrt_inv = new_processed_data_res$corr_matrix_sqrt_inv%*% whiten_cov_sqrt_inv
      
      store_ar_ma[[steps]] = c('detect_ar' = new_processed_data_res$store$detect_ar, 
                               'detect_ma' = new_processed_data_res$store$detect_ma, 
                               'use_ar' = new_processed_data_res$store$use_ar, 
                               'use_ma' = new_processed_data_res$store$use_ma)
      # acf(processed_data, lag.max = 100, main=paste0('iter', steps))
      if(save.cov){
        store_cov[[steps]] = new_processed_data_res$corr_matrix_toeplitz_row
      }
    }
    
  }
  
  
  
  return(list(whitened_data = processed_data,
              whiten_cov_sqrt_inv = whiten_cov_sqrt_inv, ## this whitening matrix not necessarily symmetric, so need to save the whole matrix
              store = list(steps = steps, cum_use_ar = cum_use_ar, cum_use_ma = cum_use_ma, store_ar_ma = store_ar_ma),
              store_cov = store_cov ## covariance matrices can be saved if needed.  
  ))
}

## the function to decorrelate data with one AR model
decorr_ts_autoAR = function(node_data, 
                            TR = 2, 
                            ic = 'aicc', 
                            max.cond.num = 1e8, ## the maximum condition number allowed for the serial correlation matrix, if larger than this threshold we will set tail correlations to zero
                            current_iter =1 ## for later iterations, may only have a few high lag serial correlations, may try to set lower lag coefficients at zero
){
  ## set the start.p according to TR 
  start.p = round(10/TR)
  
  ## use stepwise selection and compute at several start.p values
  
  ## first iteration: full high order AR model
  
  fit = list()
  
  ## always fit this moderate order procedure
  fit[[1]] = auto.arima(node_data, 
                        max.p = start.p*2, 
                        max.q = 0, ic = ic, 
                        trace= F, max.order = start.p*2, 
                        start.p= start.p, 
                        start.q=0,
                        d = 0,
                        stepwise= T,
                        approximation = T,
                        nmodels = 10, ## restrict search models to improve speed
                        allowmean=F, ## since we have de-mean data, set to false to speed up search
                        allowdrift = F
  )
  
  ## this may be very time consuming for later iterations, so skip for iteration>2
  if(current_iter > 1){
    fit[[2]] = list() 
  }else{
    fit[[2]] = auto.arima(node_data, 
                          max.p = start.p*2, 
                          max.q = 0, ic = ic, 
                          trace= F, max.order = start.p*2, 
                          start.p = start.p*2,
                          start.q=0,
                          d = 0,
                          stepwise= T,
                          approximation = T,
                          nmodels = 5,
                          allowmean=F,
                          allowdrift = F
    )
  }
  
  ## always fit this low order procedure
  fit[[3]] = auto.arima(node_data, 
                        max.p = start.p*2, 
                        max.q = 0, ic = ic, 
                        trace= F, max.order = start.p*2, 
                        start.p = round(start.p/2),
                        start.q=0,
                        d = 0,
                        stepwise= T,
                        approximation = T,
                        nmodels = 10,
                        allowmean=F,
                        allowdrift = F
  )
  
  
  fit_ic = sapply(fit, function(tmp){
    if(length(tmp)>0){
      compute_ic(tmp$loglik, p=tmp$arma[1], q=0, Time = length(node_data), k = 0)[[ic]]
    }else{
      Inf
    }})
  
  
  fit = fit[[which(fit_ic == min(fit_ic))[1]]]    
  use_p = fit$arma[1]
  use_q = fit$arma[2]
  
  ## then try option 2: roughly look at remaining high serial correlations, and fit with a few high order terms, and check if those provide better fit
  ## reason: may not need all the lags to explain the correlations and can reduce number of parameters to estimate (though is post-hoc selection of order, not rigorous)
  ## reason: include high order terms because we may skip high orders previously due to computation time
  acfs = acf(node_data, lag.max = round(2*10/TR), plot = F)
  signif_lag  = which( p.adjust(2*(1-pnorm(abs(acfs$acf) * sqrt(length(node_data)))), method='none') <0.05) ## do not adjust, just rough selection
  
  if(length(signif_lag)>0){
    ## this model only have lags corresponding to signif_lag
    tryCatch(
      {
        tmp = rep(0, max(signif_lag))
        tmp[signif_lag]<-NA
        ## this may fail due to non-stationary model estimates
        fit_sparse_lag = arima(node_data, order = c(max(signif_lag), 0, 0), include.mean = F, fixed = tmp)
        fit_ic_sparse_lag = compute_ic(fit_sparse_lag$loglik, p=sum(is.na(tmp)), q=0, k=0, Time=length(node_data))[[ic]]
        
        ## check if the new model gives better fit
        if(fit_ic_sparse_lag<min(fit_ic)){
          fit_ic = fit_ic_sparse_lag
          fit = fit_sparse_lag
          use_p = sum(is.na(tmp))
          use_q = 0
        }
      }, error = function(e){
        print('sparse sig lag fail: \n')
        print(e)
      }
    )
    
    
    ## this model have high order lags and signif_lag;
    tryCatch({
      fit_my_arima = function(p){
        p.zero = round(start.p*1.5)
        set_zero = rep(0, p.zero)
        set_zero[signif_lag[signif_lag<=p.zero]] <- NA
        return(list(res = stats::arima(node_data, order = c(p, 0, 0), fixed = c(set_zero, rep(NA, p-p.zero )), include.mean = F),
                    use_p = p-p.zero  + sum(is.na(set_zero)),
                    use_q = 0)
        )
      }
      ## this may fail due to non-stationary model fit
      res = mclapply((start.p*2-0):(start.p*2) , fit_my_arima) ## just fit one model of high order
      
      ic_tmp = sapply(res, function(x){
        if(class(x) == 'list'){
          compute_ic(x$res$loglik, 
                     p=x$use_p, 
                     q=x$use_q, k=0, Time = length(node_data))[[ic]]
        }else{
          NA
        }
      }
      )
      
      if(sum(!is.na(ic_tmp))>0){
        best_tmp = which(sapply(ic_tmp, function(x) as.numeric(x)==min(ic_tmp, na.rm=T)))
        
        fit_high_order_only = res[[best_tmp]]$res
        fit_ic_high_order_only = ic_tmp[[best_tmp]]
        
        if(fit_ic_high_order_only < min(fit_ic)){
          fit=fit_high_order_only 
          fit_ic = fit_ic_high_order_only
          use_p = res[[best_tmp]]$use_p
          use_q = res[[best_tmp]]$use_q
        }
      }
    }, error = function(e){
      print('high order + sparse sig lag fail \n')
      print(e)
    }
    )
    
    
  }else{
    ## if no significant lags, we just fit the AR0
    fit_sparse_lag = arima(node_data, order = c(0, 0, 0), include.mean = F)
    fit_ic_sparse_lag = compute_ic(fit_sparse_lag$loglik, p=0, q=0, k=0, Time=length(node_data))[[ic]]
    
    ## check if the new model gives better fit
    if(fit_ic_sparse_lag<min(fit_ic)){
      fit_ic = fit_ic_sparse_lag
      fit = fit_sparse_lag
      use_p = 0
      use_q = 0
    }
  }
  
  
  
  {  ## compute theoretical ACF based on estimated ARMA model
    arp = fit$arma[1]
    maq = fit$arma[2]
    if(arp > 0){
      arpara = fit$coef[paste0('ar', 1:arp)]
    }else{
      arpara = 0 
    }
    if(maq >0){
      mapara = fit$coef[paste0('ma', 1:maq)]
    }else{
      mapara = 0
    }
    
    acfs = ARMAacf(ar = arpara, 
                   ma = mapara, 
                   lag.max=length(node_data), pacf=F)
    
    ################
    ## account for all correlations by ARMA theoretical ACF
    
    ## new problem: COV may have bad condition number. We may still want to force some small correlations to 0.
    ## new problem: COV may not be positive definite??
    COV0 = toeplitz(acfs[1:length(node_data)])
    k=1
    while(kappa(COV0) > max.cond.num & k < floor(length(node_data)/10)){
      # print(k)
      COV0 = toeplitz(c(COV0[1,1:(ncol(COV0)-k*10)], rep(0, k*10)))
      k=k+1
    }
    
    corr_matrix_sqrt_inv = sqrt_mat_decorr(COV0, diag(1, length(node_data)))
    decorr_resid = corr_matrix_sqrt_inv %*% node_data
    
    ## problem with outlier removal: this makes the de-correlating matrix not really matching the data; this may make the decor_mat * regressor invalid 
    
    # tmp = tsoutliers(decorr_resid) ## in case there are very extreme outliers
    # decorr_resid[tmp$index] <- tmp$replacements
  }
  
  # plot(decorr_resid, type='l')
  # acf(decorr_resid, lag.max = 100)
  
  
  return(list(decorr_resid = decorr_resid, ## the whitened data after one iteration of auto.arima
              ## !!! the sqrt-inv of COV0 is not toeplitz any more! to save space, we can only store the COV0 first row for all iterations
              corr_matrix_sqrt_inv = corr_matrix_sqrt_inv, 
              corr_matrix_toeplitz_row = COV0[1,],
              store = list(detect_ar = arp,
                           detect_ma = maq,
                           use_ar  = use_p,
                           use_ma = use_q)
  )
  )
}



###########################################################
## other utility functions
###########################################################

## the function to de-correlation data vector v based on covariance matrix A
sqrt_mat_decorr = function(A, v){
    A.eigen = eigen(A)
    A.eigen$vectors %*% diag(1/sqrt(A.eigen$values)) %*% t(A.eigen$vectors) %*% v
}


## function to compute AIC, BIC, AICc by hand in order to unify these numbers across different packages
## for auto.arima only use when no drift and no mean is fitted
## note that the default included aic of arima seems to be wrong
## https://otexts.com/fpp2/arima-estimation.html
compute_ic = function(loglik, p, q, k, Time){ ## enter free p and q values; k=1 if include mean and k=0 if not
  aic = -2*loglik + 2*(p+q+k+1)
  bic = aic + (p+q+k+1)*(log(Time)-2)
  aicc = aic + 2*(p+q+k+1)*(p+q+k+2)/(Time-p-q-k-2)
  return(list(aic = aic, bic=bic, aicc=aicc))
}





