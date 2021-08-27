# Functions Associated with SNP-wise Regression

# log10BF
# To calculate log10 Bayes Factor for any additive model, given a prior
# variance (often set to 0.5)
#
# Inputs: matrix of predictors g, outcome vector y, prior variance g
#
# Returns: log10 Bayes Factor of model with all predictors included
log10BF = function(g,y,sigmaa) {
  p=dim(g)[2]
  n=dim(g)[1]
  if (is.null(dim(g)[2])){
    g = g - mean(g)
    y = y - mean(y)
    n=length(g)
    X = g
    invnu = 1/sigmaa^2
    invOmega = invnu + t(X) %*% X
    B = (t(X) %*% cbind(y))/invOmega
    invOmega0 = n
    return(-0.5*log10(det(invOmega)) + 0.5*log10(invOmega0) - log10(sigmaa) - (n/2)*(log10(t(y- X %*% B) %*% y) - log10(t(y) %*% y)))
  }
  else {
    g = scale(g, scale = FALSE)
    y = scale(y, scale = FALSE)
    X = g
    invnu = diag(rep(1/sigmaa^2, p))
    invOmega = invnu + t(X) %*% X
    B = solve(invOmega, t(X) %*% cbind(y))
    invOmega0 = n
    return(-0.5*log10(det(invOmega)) + 0.5*log10(invOmega0) - p*log10(sigmaa) - (n/2)*(log10(t(y- X %*% B) %*% y) - log10(t(y) %*% y - n*mean(y)^2)))
  }
}


# snpwise_weights
# To calculate weights for models in the SNP-wise regression
#
# Inputs: predictor matrix X, outcome vector y, prior variance, optional
# prior value (defaults to 1/p)
#
# Returns: model by variables included, log10 posterior scores, posterior
# model probabilities (PMPs), log10 normalizing constant, log10 Bayes
# Factor
snpwise_weights = function(X,y,sigmaa,priorpi) {
  models <- c()
  log10post <- c()
  log10bf <- c()
  
  #get number of parameters
  if (is.null(dim(X)[2])) {
    par <- 1
  }
  else {
    par <- dim(X)[2]
  }
  
  #if no prior input, set to 1/number of parameters
  if (missing(priorpi)){
    p <- par
  }
  else {
    p <- 1/priorpi
  }
  
  #get all possible parameter combinations
  l <- rep(list(0:1),par)
  combs <- expand.grid(l)
  m <- dim(combs)[1]
  
  #loop through each possible model
  for (i in c(1:m)){
    #matrix of values of all variables in model
    incl <- c()
    #number of variables in model
    numincl <- 0
    
    #check each variable to see if in model
    for (j in c(1:par)){
      if (combs[i,j]==1){
        incl <- cbind(incl, X[,j])
        numincl <- numincl + 1
      }
    }
    
    #if only one variable in model
    if (numincl==1) {
      for (k in c(1:par)){
        if (combs[i,k]==1){
          varin <- k
        }
      }
      newmod <- paste0("x", varin)
      
      #get weight
      prior <- (1/p)*(1-1/p)^(par-1)
      lbf <- log10BF(X[,varin], y, sigmaa)
      log10bf <- c(log10bf, lbf)
      post <- log10(prior) + lbf
    }
    
    #if multiple variables in model
    else if (numincl != 0){
      newmod <- ""
      for (j in c(1:par)){
        
        if (combs[i,j]==1){
          var <- paste0("x", j)
          newmod <- paste(newmod, sep =",", var)
        }
      }
      newmod <- sub('.', '', newmod)
      #get weight
      prior <- (1/p)^numincl * (1-1/p)^(numincl-1)
      lbf <- log10BF(incl, y, sigmaa)
      log10bf <- c(log10bf, lbf)
      post <- log10(prior) + lbf
    }
    
    #if null model
    else {
      log10bf <- c(log10bf, 0)
      post <- log10((1-1/p)^(p))
      newmod <- "NULL"
    }
    
    models <- c(models, newmod)
    #add weight to list
    log10post <- c(log10post, post)
  }
  
  #find log10 normalizing constant
  maxlog10post <- max(log10post)
  log10norm <- maxlog10post + log10(sum(10^(log10post-maxlog10post)))
  wts <- 10^((log10post)-log10norm)
  log10weights <- log10(wts)
  
  #return dataframe of each model and its posterior weight
  toreturn <- data.frame(models, log10post, wts, log10norm, log10bf)
  names(toreturn)[1] <- "models"
  names(toreturn)[2] <- "log10postscores"
  names(toreturn)[3] <- "pmps"
  names(toreturn)[4] <- "log10norm"
  names(toreturn)[5] <- "log10bayesfactor"
  
  return(toreturn)
}

# snpwisereg
# Runs full SNP-wise regression algorithm
#
# Inputs: predictor matrix X, outcome vector y, prior variance, optional
# prior value (defaults to 1/p)
#
# Returns: posterior inclusion probabilities (PIPs) for each predictor

snpwisereg = function(X,y,sigmaa, priorpi){
  normalize <- 0
  #if no prior input, set to 1/number of parameters
  if (missing(priorpi)){
    p <- par
  }
  else {
    p <- 1/priorpi
  }
  
  df <- data.frame(X,y)
  #get number of variables
  snps <- dim(X)[2]
  #vector of variable names
  snpnames <- c()
  #vector of variable PIPs
  pips <- rep(0,snps)
  
  for (i in c(1:snps)){
    colnames(df)[i] <- paste0("x",i)
    snpnames <- c(snpnames, paste0("x",i))
  }
  genweights <- snpwise_weights(X,y,sigmaa,priorpi)
  
  #loop through each variable and stepwise regress while forcing to 
  #include the variable
  for (i in c(1:snps)){
    reg <- stepwise(df, y = colnames(df)[snps+1], include = colnames(df)[i], selection = "forward")
    selected <- reg$variate
    #filter out occasional duplicates in selected variables
    selected <- unique(selected)
    if (selected[1]=="intercept"){
      selected <- selected[-1]
    }
    incl <- ""
    numincl <- 0
    
    #make list of variables selected
    for (j in c(1:snps)) {
      if (paste0("x",j) %in% selected) {
        numincl <- numincl + 1
        if (incl == "") {
          incl <- paste0("x",j)
        }
        else {
          toincl <- paste0("x",j)
          incl <- paste(incl, toincl, sep=",")
        }
      }
    }
    
    #filter data to only include selected variables
    ourrow <- filter(genweights, models == incl)
    #update normalizing constant
    normalize <- normalize + ourrow$pmps
    #add model PMP to included variables' PIPs
    for (k in c(1:snps)){
      if (colnames(df)[k] %in% selected){
        pips[k] <- pips[k] + ourrow$pmps}
    }
  }
  
  pips <- pips/normalize
  
  #return snp names and corresponding PIPs
  toreturn <- data.frame(snpnames, pips)
  names(toreturn)[1] <- "snp"
  names(toreturn)[2] <- "pip"
  return(toreturn)
}