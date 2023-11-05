#########################
## function to whiten the data for a node of a provided subject's observation using IDAR, the IDAR-iter1 method, the AR1 and ARMA11 models

whiten_one_subj = function(node_data, 
                           node, 
                           max.iter){
  
  node_data = scale(node_data)
  cat('node', node, '\n')
  
  
  tmp = tsoutliers(node_data) ## in case there are very extreme outliers: e.g. subj1 node1
  node_data[tmp$index] <- tmp$replacements
  
  ##############################################
  ## implement whitening procedures
  
  tmp_whiten <- tmp_iter1 <- tmp_AR1 <- tmp_ARMA11 <- 
    matrix(0, nrow=nrow(node_data), ncol = 1, dimnames = list(NULL, paste0('node', node)))
  tmp_cov_whiten <- tmp_cov_iter1 <- tmp_cov_AR1 <- tmp_cov_ARMA11 <- tmp_cov_whiten_all_steps <- list()
  tmp_store_arma = list()
  
  tryCatch({
    
    
    
    
    ### use the proposed procedure for 5 iterations
    fit_res = iterative_decorr_ts_autoAR(node_data = node_data, TR = TR, ic = 'aicc', max.iter=max.iter)
    tmp_whiten[,1] <- fit_res$whitened_data
    tmp_cov_whiten[[1]] <- fit_res$whiten_cov_sqrt_inv
    tmp_store_arma[[1]] = fit_res$store
    tmp_cov_whiten_all_steps[[1]] = fit_res$store_cov
    
    
    
    
    ### use the proposed procedure but with only 1 iteration
    fit_res_iter1 = iterative_decorr_ts_autoAR(node_data = node_data, TR = TR, ic = 'aicc', max.iter=1)
    tmp_iter1[,1] <- fit_res_iter1$whitened_data
    tmp_cov_iter1[[1]] <- fit_res_iter1$whiten_cov_sqrt_inv
    
    
    ### use the AR1 model
    tryCatch({
      fit_AR1 = arima(node_data, order = c(1, 0, 0))
      tmp_cov_AR1[[1]]<- sqrt_mat_decorr(
        toeplitz(c(ARMAacf(ar = fit_AR1$coef, lag.max = length(node_data)-1))),
        diag(1, length(node_data))
      )
      tmp_AR1[,1] = tmp_cov_AR1[[1]] %*% node_data
    }, error = function(e)print('AR1 fail \n'))
    
    
    ### use the ARMA(1,1) model
    tryCatch({
      fit_ARMA11 = arima(node_data, order = c(1, 0, 1))
      tmp_cov_ARMA11[[1]]<- sqrt_mat_decorr(
        toeplitz(c(ARMAacf(ar = fit_ARMA11$coef['ar1'], ma=fit_ARMA11$coef['ma1'], lag.max = length(node_data)-1))),
        diag(1, length(node_data))
      )
      tmp_ARMA11[,1] = tmp_cov_ARMA11[[1]] %*% node_data
      
    }, error = function(e)print('ARMA11 fail \n'))
    
    
    
  }, error = function(e)print(e))
  
  
  names(tmp_store_arma) <- paste0('node', node)
  
  
  return(list(
    
    subj_time_series_whiten =  tmp_whiten,
    subj_time_series_iter1  = tmp_iter1,
    subj_time_series_AR1  = tmp_AR1,
    subj_time_series_ARMA11  = tmp_ARMA11,
    subj_time_series_nopw = node_data,
    
    subj_cov_whiten = tmp_cov_whiten[[1]],
    subj_cov_iter1 = tmp_cov_iter1[[1]],
    subj_cov_AR1=  tmp_cov_AR1[[1]],
    subj_cov_ARMA11= tmp_cov_ARMA11[[1]],
    
    subj_store_ar_ma = tmp_store_arma,
    subj_store_cov_whiten_all_steps = tmp_cov_whiten_all_steps[[1]]
    
  ))  
  
}



## run an example

whiten_one_subj(
  matrix(rnorm(100*8), ncol=8),
  2,
  5
)