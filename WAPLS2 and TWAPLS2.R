#' WA-PLS training function 2
#' 
#' WA-PLS training function, which can perform \code{fx} correction. 1/fx 
#' correction will be applied at step 2 and step 7.
#' 
#' @importFrom stats lm
#' 
#' @param modern_taxa The modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate The modern climate value at each sampling site.
#' @param nPLS The number of components to be extracted.
#' @param usefx Boolean flag on whether or not use \code{fx} correction.
#' @param fx_method Binned or p-spline smoothed \code{fx} correction: if 
#'     \code{usefx = FALSE}, this should be \code{NA}; otherwise, 
#'     \code{\link{fx}} function will be used when choosing "bin";
#'     \code{\link{fx_pspline}} function will be used when choosing "pspline".
#' @param bin Binwidth to get fx, needed for both binned and p-splined method.
#'     if \code{usefx = FALSE}, this should be \code{NA};
#' @return A list of the training results, which will be used by the predict 
#'     function. Each element in the list is described below:
#'     \describe{
#'     \item{\code{fit}}{the fitted values using each number of components.}
#'     \item{\code{x}}{the observed modern climate values.}
#'     \item{\code{taxon_name}}{the name of each taxon.}
#'     \item{\code{optimum}}{the updated taxon optimum (u* in the WA-PLS 
#'     paper).}
#'     \item{\code{comp}}{each component extracted (will be used in step 7 
#'     regression).}
#'     \item{\code{u}}{taxon optimum for each component (step 2).}
#'     \item{\code{z}}{a parameter used in standardization for each component 
#'     (step 5).}
#'     \item{\code{s}}{a parameter used in standardization for each component 
#'     (step 5).}
#'     \item{\code{orth}}{a list that stores orthogonalization parameters 
#'     (step 4).}
#'     \item{\code{alpha}}{a list that stores regression coefficients (step 7).}
#'     \item{\code{meanx}}{mean value of the observed modern climate values.}
#'     \item{\code{nPLS}}{the total number of components extracted.}
#'     }
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' 
#' # MTCO
#' fit_Tmin2 <- fxTWAPLS::WAPLS.w2(taxa, modern_pollen$Tmin, nPLS = 5)
#' fit_f_Tmin2 <- fxTWAPLS::WAPLS.w2(taxa, 
#'                                 modern_pollen$Tmin, 
#'                                 nPLS = 5, 
#'                                 usefx = TRUE, 
#'                                 fx_method = "bin",
#'                                 bin = 0.02)
#' }
#' 
#' @seealso \code{\link{fx}}, \code{\link{TWAPLS.w}}, and
#'     \code{\link{WAPLS.predict.w}}
WAPLS.w2 <- function(modern_taxa, 
                    modern_climate, 
                    nPLS = 5, 
                    usefx = FALSE, 
                    fx_method = "bin",
                    bin = NA) {
  # Step 0. Centre the environmental variable by subtracting the weighted mean
  x <- modern_climate
  y <- modern_taxa
  y <- as.matrix(y)
  y <- y/rowSums(y)
  
  nc <- ncol(modern_taxa)
  nr <- nrow(modern_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  # Define some matrix to store the values
  u <- matrix(NA, nc, nPLS) # u of each component
  # u of each component, standardized the same way as r
  u_sd <- matrix(NA, nc, nPLS) 
  optimum <- matrix(NA, nc, nPLS) # u updated
  r <- matrix(NA, nr, nPLS) # site score
  z <- matrix(NA, 1, nPLS) # standardize
  s <- matrix(NA, 1, nPLS) # standardize
  orth <- list() # store orthogonalization parameters
  alpha <- list() # store regression coefficients
  comp <- matrix(NA, nr, nPLS) # each component
  fit <- matrix(NA, nr, nPLS) # current estimate
  fr <- matrix(NA, nr, nPLS) # frequency of site score

  
  pls <- 1
  # Step 1. Take the centred environmental variable(xi) as initial site 
  # scores (ri). 
  r[, pls] <- x - mean(x)
  
  # Step 2. Calculate new species scores (uk* by weighted averaging of the 
  # site scores) 
  if(usefx==FALSE){
    # uk = sumi_yik*xi/sumi_yik; 
    u[, pls] <- t(y) %*% r[, pls] / sumi_yik 
    
  }else{
    if(fx_method=="bin"){fr[,pls]<-fxTWAPLS::fx(r[,pls],bin=bin)}else{
      fr[,pls]<-fxTWAPLS::fx_pspline(r[,pls],bin=bin)}
    #uk=sumi_(yik*xi/fxi)/sumi_(yik/fxi); 
    u[, pls] <- t(y/fr[,pls]) %*% r[, pls] / colSums(y/fr[,pls]) 
  }
  
  # Step 3. Calculate new site scores (ri) by weighted averaging of the species 
  # scores
  r[, pls] <- y %*% u[, pls] / rowSums(y) # xi=sumk_yik*uk/sumk_yik; 1*nsite
  
  # Step 4. For the first axis go to Step 5.
  
  # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
  z[, pls] <- mean(r[, pls], na.rm = TRUE)
  s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  
  # Step 6. Take the standardized score as the new component
  comp[, pls] <- r[, pls]
  
  # Step 7. Regress the environmental variable (xJ on the components obtained 
  # so far using weights and take the fitted values as current estimates
  if(usefx == FALSE) {
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                    weights = sumk_yik / Ytottot)
  } else{
    
    if(fx_method=="bin"){
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                      weights = 1 / fxTWAPLS::fx(x,bin) )
    }else{
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                      weights = 1 / fxTWAPLS::fx_pspline(x,bin) )
    }
  }
  
  fit[, pls] <- lm[["fitted.values"]]
  alpha[[pls]] <- lm[["coefficients"]]
  u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
  optimum[, pls] <- alpha[[pls]][1] + u_sd[, pls] * alpha[[pls]][2]
  
  for(pls in 2:nPLS) {
    # Go to Step 2 with the residuals of the regression as the new site 
    # scores (rJ.
    r[, pls] <- lm[["residuals"]]
    
    # Step 2. Calculate new species scores (uk* by weighted averaging of the 
    # site scores) uk=sumi_(yik*xi/fxi)/sumi_(yik/fxi); 
    if(usefx==FALSE){
      # uk = sumi_yik*xi/sumi_yik; 
      u[, pls] <- t(y) %*% r[, pls] / sumi_yik 
      
    }else{
      if(fx_method=="bin"){fr[,pls]<-fxTWAPLS::fx(r[,pls],bin=bin)}else{
        fr[,pls]<-fxTWAPLS::fx_pspline(r[,pls],bin=bin)}
      #uk=sumi_(yik*xi/fxi)/sumi_(yik/fxi); 
      u[, pls] <- t(y/fr[,pls]) %*% r[, pls] / colSums(y/fr[,pls]) 
    }
    
    # Step 3. Calculate new site scores (ri) by weighted averaging of the 
    # species scores
    # xi=sumk_yik*uk/sumk_yik; 1*nsite
    r[, pls] <- y %*% u[, pls] / rowSums(y) 
    
    # Step 4. For second and higher components, make the new site scores (ri) 
    # uncorrelated with the previous components by orthogonalization 
    # (Ter Braak, 1987 : Table 5 .2b)
    v <- rep(NA, pls - 1)
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      v[pls - j] <- sum(sumk_yik * fi * xi) / Ytottot
      xinew <- xi - v[pls - j] * fi
    }
    orth[[pls]] <- v
    r[, pls] <- xinew
    
    # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
    z[, pls] <- mean(r[, pls], na.rm = TRUE)
    s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    
    # Step 6. Take the standardized score as the new component
    comp[, pls] <- r[, pls]
    
    # Step 7. Regress the environmental variable on the components obtained so 
    # far using weights and take the fitted values as current estimates 
    if(usefx == FALSE) {
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                      weights = sumk_yik / Ytottot)
    } else{
      
      if(fx_method=="bin"){
        lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                        weights = 1 / fxTWAPLS::fx(x,bin) )
      }else{
        lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                        weights = 1 / fxTWAPLS::fx_pspline(x,bin) )
      }
    }
    
    fit[, pls] <- lm[["fitted.values"]]
    alpha[[pls]] <- lm[["coefficients"]]
    
    u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
    optimum[, pls] <-
      alpha[[pls]][1] + u_sd[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, 
               modern_climate, 
               colnames(modern_taxa), 
               optimum,
               comp, 
               u, 
               z, 
               s, 
               orth, 
               alpha, 
               mean(modern_climate), 
               nPLS)
  names(list) <- c("fit", 
                   "x", 
                   "taxon_name", 
                   "optimum", 
                   "comp", 
                   "u", 
                   "z", 
                   "s", 
                   "orth", 
                   "alpha", 
                   "meanx", 
                   "nPLS")
  return(list)
}

#' TWA-PLS training function 2
#' 
#' TWA-PLS training function, which can perform \code{fx} correction. 1/fx
#' correction will be applied at step 2 and step 7.
#' 
#' @importFrom stats lm
#' 
#' @inheritParams WAPLS.w
#'
#' @return A list of the training results, which will be used by the predict 
#'     function. Each element in the list is described below:
#'     \describe{
#'     \item{\code{fit}}{the fitted values using each number of components.}
#'     \item{\code{x}}{the observed modern climate values.}
#'     \item{\code{taxon_name}}{the name of each taxon.}
#'     \item{\code{optimum}}{the updated taxon optimum}
#'     \item{\code{comp}}{each component extracted (will be used in step 7 
#'     regression).}
#'     \item{\code{u}}{taxon optimum for each component (step 2).}
#'     \item{\code{t}}{taxon tolerance for each component (step 2).}
#'     \item{\code{z}}{a parameter used in standardization for each component 
#'     (step 5).}
#'     \item{\code{s}}{a parameter used in standardization for each component 
#'     (step 5).}
#'     \item{\code{orth}}{a list that stores orthogonalization parameters 
#'     (step 4).}
#'     \item{\code{alpha}}{a list that stores regression coefficients (step 7).}
#'     \item{\code{meanx}}{mean value of the observed modern climate values.}
#'     \item{\code{nPLS}}{the total number of components extracted.}
#'     }
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' 
#' # MTCO
#' fit_t_Tmin2 <- fxTWAPLS::TWAPLS.w2(taxa, modern_pollen$Tmin, nPLS = 5)
#' fit_tf_Tmin2 <- fxTWAPLS::TWAPLS.w2(taxa, 
#'                                   modern_pollen$Tmin, 
#'                                   nPLS = 5, 
#'                                   usefx = TRUE, 
#'                                   fx_method = "bin",
#'                                   bin = 0.02)
#' }
#' 
#' @seealso \code{\link{fx}}, \code{\link{TWAPLS.predict.w}}, and
#'     \code{\link{WAPLS.w}}
TWAPLS.w2 <- function(modern_taxa,
                     modern_climate,
                     nPLS = 5,
                     usefx = FALSE,
                     fx_method = "bin",
                     bin = NA){
  # Step 0. Centre the environmental variable by subtracting the weighted mean
  x <- modern_climate
  y <- modern_taxa
  y <- as.matrix(y)
  y <- y/rowSums(y)
  
  nc <- ncol(modern_taxa)
  nr <- nrow(modern_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  #Define some matrix to store the values
  u <- matrix(NA, nc, nPLS) # u of each component
  # u of each component, standardized the same way as r
  u_sd <- matrix(NA, nc, nPLS)
  optimum <- matrix(NA, nc, nPLS) # u updated
  t <- matrix(NA, nc, nPLS) # tolerance
  r <- matrix(NA, nr, nPLS) # site score
  z <- matrix(NA, 1, nPLS) # standardize
  s <- matrix(NA, 1, nPLS) # standardize
  orth <- list() # store orthogonalization parameters
  alpha <- list() # store regression coefficients
  comp <- matrix(NA, nr, nPLS) # each component
  fit <- matrix(NA, nr, nPLS) # current estimate
  fr <- matrix(NA, nr, nPLS) # frequency of site score
  
  pls <- 1
  # Step 1. Take the centred environmental variable (xi) as initial site 
  # scores (ri). 
  r[, pls] <- x - mean(x)
  
  # Step 2. Calculate uk and tk
  if(usefx==FALSE){
    u[, pls] <- t(y) %*% r[, pls] / sumi_yik # uk=sumi_yik*xi/sumi_yik; 
    n2 <- matrix(NA, nc, 1)
    for (k in 1:nc) {
      t[k, pls] <- sqrt(sum(y[, k] * (r[, pls] - u[k, pls]) ^ 2) / sumi_yik[k])
      n2[k] <- 1 / sum((y[, k] / sum(y[, k])) ^ 2)
      t[k, pls] <- t[k, pls] / sqrt(1 - 1 / n2[k])
    }
    
  }else{
    if(fx_method=="bin"){fr[,pls]<-fxTWAPLS::fx(r[,pls],bin=bin)}else{
      fr[,pls]<-fxTWAPLS::fx_pspline(r[,pls],bin=bin)}
    #uk=sumi_(yik*xi/fxi)/sumi_(yik/fxi); 
    u[, pls] <- t(y/fr[,pls]) %*% r[, pls] / colSums(y/fr[,pls]) 
    n2 <- matrix(NA, nc, 1)
    for (k in 1:nc) {
      t[k, pls] <- sqrt(sum(y[, k] /fr[,pls]* (r[, pls] - u[k, pls]) ^ 2) / colSums(y/fr[,pls])[k])
      n2[k] <- 1 / sum((y[, k]/fr[,pls] / sum(y[, k]/fr[,pls])) ^ 2)
      t[k, pls] <- t[k, pls] / sqrt(1 - 1 / n2[k])
    }
  }
  
  
  
  # Step 3. Calculate new site scores (ri)
  #xi; 1*nsite
  r[, pls] <- (y %*% (u[, pls] / t[, pls] ^ 2)) / (y %*% (1 / t[, pls] ^ 2))
  
  # Step 4. For the first axis go to Step 5.
  
  # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
  z[, pls] <- mean(r[, pls], na.rm = TRUE)
  s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  
  # Step 6. Take the standardized score as the new component
  comp[, pls] <- r[, pls]
  
  # Step 7. Regress the environmental variable on the components obtained so far
  # using weights and take the fitted values as current estimates 
  if (usefx == FALSE) {
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                    weights = sumk_yik / Ytottot)
  } else{
    
    if(fx_method=="bin"){
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                      weights = 1 / fxTWAPLS::fx(x,bin) )
    }else{
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                      weights = 1 / fxTWAPLS::fx_pspline(x,bin) )
    }
  }
  
  fit[, pls] <- lm[["fitted.values"]]
  alpha[[pls]] <- lm[["coefficients"]]
  u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
  optimum[, pls] <- alpha[[pls]][1] + u_sd[, pls] * alpha[[pls]][2]
  
  for (pls in 2:nPLS) {
    # Go to Step 2 with the residuals of the regression as the new site 
    # scores (ri).
    r[, pls] <- lm[["residuals"]]
    
    # Step 2. Calculate new uk and tk
    if(usefx==FALSE){
      u[, pls] <- t(y) %*% r[, pls] / sumi_yik # uk=sumi_yik*xi/sumi_yik; 
      n2 <- matrix(NA, nc, 1)
      for (k in 1:nc) {
        t[k, pls] <- sqrt(sum(y[, k] * (r[, pls] - u[k, pls]) ^ 2) / sumi_yik[k])
        n2[k] <- 1 / sum((y[, k] / sum(y[, k])) ^ 2)
        t[k, pls] <- t[k, pls] / sqrt(1 - 1 / n2[k])
      }
      
    }else{
      if(fx_method=="bin"){fr[,pls]<-fxTWAPLS::fx(r[,pls],bin=bin)}else{
        fr[,pls]<-fxTWAPLS::fx_pspline(r[,pls],bin=bin)}
      #uk=sumi_(yik*xi/fxi)/sumi_(yik/fxi); 
      u[, pls] <- t(y/fr[,pls]) %*% r[, pls] / colSums(y/fr[,pls]) 
      n2 <- matrix(NA, nc, 1)
      for (k in 1:nc) {
        t[k, pls] <- sqrt(sum(y[, k] /fr[,pls]* (r[, pls] - u[k, pls]) ^ 2) / colSums(y/fr[,pls])[k])
        n2[k] <- 1 / sum((y[, k]/fr[,pls] / sum(y[, k]/fr[,pls])) ^ 2)
        t[k, pls] <- t[k, pls] / sqrt(1 - 1 / n2[k])
      }
    }
    
    # Step 3. Calculate new site scores (r;) by weighted averaging of the 
    # species scores, i.e. new
    r[,pls] <- (y%*%(u[,pls]/t[,pls]^2))/(y%*%(1/t[,pls]^2)); #xi; 1*nsite
    
    # Step 4. For second and higher components, make the new site scores (r;) 
    # uncorrelated with the previous components by orthogonalization 
    # (Ter Braak, 1987 : Table 5 .2b)
    v <- rep(NA, pls - 1)
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      v[pls - j] <- sum(sumk_yik * fi * xi) / Ytottot
      xinew <- xi - v[pls - j] * fi
    }
    orth[[pls]] <- v
    # plot(xinew~r[,pls]);abline(0,1)
    r[, pls] <- xinew
    
    # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
    z[, pls] <- mean(r[, pls], na.rm = TRUE)
    s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    
    # Step 6. Take the standardized scores as the new component
    comp[, pls] <- r[, pls]
    
    # Step 7. Regress the environmental variable (xJ on the components obtained 
    # so far using weights and take the fitted values as current estimates 
    if (usefx == FALSE) {
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                      weights = sumk_yik / Ytottot)
    } else{
      
      if(fx_method=="bin"){
        lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                        weights = 1 / fxTWAPLS::fx(x,bin) )
      }else{
        lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                        weights = 1 / fxTWAPLS::fx_pspline(x,bin) )
      }
    }
    
    fit[, pls] <- lm[["fitted.values"]]
    alpha[[pls]] <- lm[["coefficients"]]
    
    u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
    optimum[, pls] <-
      alpha[[pls]][1] + u_sd[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, 
               modern_climate, 
               colnames(modern_taxa), 
               optimum, 
               comp, 
               u, 
               t, 
               z, 
               s, 
               orth, 
               alpha, 
               mean(modern_climate), 
               nPLS)
  names(list) <- c("fit", 
                   "x", 
                   "taxon_name", 
                   "optimum", 
                   "comp", 
                   "u", 
                   "t", 
                   "z", 
                   "s", 
                   "orth", 
                   "alpha", 
                   "meanx", 
                   "nPLS")
  return(list)
}