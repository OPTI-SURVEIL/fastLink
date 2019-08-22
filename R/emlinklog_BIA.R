#' emlinklog_BIA
#'
#' Expectation-Maximization with Bayesian initial value selection (https://arxiv.org/pdf/1504.06870.pdf) for Record Linkage 
#' allowing for dependencies across linkage fields
#'
#' @usage emlinklog_BIA(patterns, nobs.a, nobs.b, iter.max, tol, varnames, starts, st.iters)
#'
#' @param patterns table that holds the counts for each unique agreement
#' pattern. This object is produced by the function: tableCounts.
#' @param nobs.a Number of observations in dataset A
#' @param nobs.b Number of observations in dataset B
#' @param iter.max Max number of iterations. Default is 5000
#' @param tol Convergence tolerance. Default is 1e-05
#' @param varnames The vector of variable names used for matching. Automatically provided if using \code{fastLink()} wrapper. Used for
#' clean visualization of EM results in summary functions.
#' @param starts The number of initial parameter sets to be explored. Defaults to 150.
#' @param st.iters The number of EM iterations for each initial parameter set prior to evaluation. Defaults to 150.
#'
#' @return \code{emlinklog_cv} returns a list with the following components:
#' \item{zeta.j}{The posterior match probabilities for each unique pattern.}
#' \item{p.m}{The probability of finding a match.}
#' \item{p.u}{The probability of finding a non-match.}
#' \item{p.gamma.j.m}{The probability of observing a particular agreement pattern conditional on being in the set of matches.}
#' \item{p.gamma.j.u}{The probability of observing a particular agreement pattern conditional on being in the set of non-matches.}
#' \item{patterns.w}{Counts of the agreement patterns observed, along with the Felligi-Sunter Weights.}
#' \item{iter.converge}{The number of iterations it took the EM algorithm to converge.}
#' \item{nobs.a}{The number of observations in dataset A.}
#' \item{nobs.b}{The number of observations in dataset B.}
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com> and Benjamin Fifield
#'
#' @examples
#' \dontrun{
#' ## Calculate gammas
#' g1 <- gammaCKpar(dfA$firstname, dfB$firstname)
#' g2 <- gammaCKpar(dfA$middlename, dfB$middlename)
#' g3 <- gammaCKpar(dfA$lastname, dfB$lastname)
#' g4 <- gammaKpar(dfA$birthyear, dfB$birthyear)
#'
#' ## Run tableCounts
#' tc <- tableCounts(list(g1, g2, g3, g4), nobs.a = nrow(dfA), nobs.b = nrow(dfB))
#'
#' ## Run EM
#' em.log <- emlinklog_BIA(tc, nobs.a = nrow(dfA), nobs.b = nrow(dfB))
#' }
#'
#' @export
#' @importFrom gtools rdirichlet
#' @importFrom stats glm model.matrix
emlinklog_BIA <- function(patterns, nobs.a, nobs.b, iter.max = 5000, tol = 1e-5, varnames = NULL, starts = 150, st.iters = 150) {
  dmultinom_robust = function (x, size = NULL, prob, log = FALSE){
    K <- length(prob)
    if (length(x) != K) 
      stop("x[] and prob[] must be equal length vectors.")
    if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 
        0) 
      stop("probabilities must be finite, non-negative and not all 0")
    prob <- prob/s
    x <- round(x)
    if (any(x < 0)) 
      stop("'x' must be non-negative")
    N <- sum(x)
    if (is.null(size)) 
      size <- N
    else if (size != N) 
      stop("size != sum(x), i.e. one is wrong")
    i0 <- prob == 0
    if (any(i0)) {
      if (any(x[i0] != 0)) 
        return(if (log) -Inf else 0)
      if (all(i0)) 
        return(if (log) 0 else 1)
      x <- x[!i0]
      prob <- prob[!i0]
    }
    r <- lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
    if (log) 
      r
    else exp(r)
  }
  if(is.null(varnames)) varnames = 1:(ncol(patterns)-1)
    ## OPTIONS  
    ## draw p.m from a reasonable beta distribution
    p.m.mode = min(c(nobs.a,nobs.b)/sum(patterns[,'counts']))
    if(p.m.mode > 1) p.m.mode = 0.5
    alpha = 1.05; beta = pmax(.05,3*p.m.mode)/p.m.mode - 1.95
    
    options(digits=16)
    
    #cluster initialization
     n.cores = parallel::detectCores() - 1
     cl = makeCluster(n.cores)
     registerDoParallel(cl)
    
    ## Edge case
    if(is.null(nrow(patterns))){
        patterns <- as.data.frame(t(as.matrix(patterns)))
    }
    nvar = length(varnames)
    keep = sapply(1:nvar,function(j) length(unique(patterns[,j])) > 1)
    
    ## Number of fields
    nfeatures <- sum(keep)
    
    ## Patterns:
    gamma.j.k <- as.matrix(patterns[, which(keep)])

    ## Patterns counts:
    n.j <- patterns[,'counts']  # Counts
    
    ## Number of unique patterns:
    N <- nrow(gamma.j.k)
    
    ## parallelize initial parameter sampling and EM
    p.m = rbeta(starts,alpha, beta)
    p.u <- 1 - p.m
    
    p.gamma.k.m <- p.gamma.k.u <- NULL
    ## Field specific probability of observing gamma.k conditional on M
    ## Mild assumption that gamma != 0 more likely
    if (is.null(p.gamma.k.m)) {
      p.gamma.k.m <- list()
      for (i in 1:nfeatures) {
        l.m <- length(unique(na.omit(gamma.j.k[, i])))
        c.m <- seq(from = 1, to = 4*l.m, by = 4)
        p.gamma.k.m[[i]] <- rdirichlet(starts, c.m)
      }
    }
    
    ## Field specific probability of observing gamma.k conditional on U
    ## Assume smaller than among M
    if (is.null(p.gamma.k.u)) {
      p.gamma.k.u <- list()
      for (i in 1:nfeatures) {
        l.u <- length(unique(na.omit(gamma.j.k[, i])))
        tempmat = p.gamma.k.m[[i]]
        
        for(l in 2:l.u) tempmat[,l] = runif(starts, 0, tempmat [,l])
        tempmat[,1] = 1 - rowSums(as.matrix(tempmat[,2:l.u]))
        p.gamma.k.u[[i]] <- tempmat
      }
    }
    
    p.gamma.k.j.m <- matrix(rep(NA, N * nfeatures), nrow = nfeatures, ncol = N)
    p.gamma.k.j.u <- matrix(rep(NA, N * nfeatures), nrow = nfeatures, ncol = N)
    
    p.gamma.j.m <- matrix(rep(NA, N), nrow = N, ncol = 1)
    p.gamma.j.u <- matrix(rep(NA, N), nrow = N, ncol = 1)
    
    logxpy <- function(lx,ly) {
      temp <- cbind(lx, ly)
      apply(temp, 1, max) + log1p(exp(-abs(lx-ly)))
    }
    
    pat <- data.frame(gamma.j.k)
    pat[is.na(pat)] <- -1
    pat <- replace(pat, TRUE, lapply(pat, factor))
    factors <- model.matrix(~ ., pat)[,-1] #remove intercept column
    c <- 1e-06
    
    logsumexp = function(x1,x2){
      mx = max(x1,x2)
      mx + log(exp(x1-mx)+exp(x2-mx))
    }
    
    g <- function(fit) {
      logwt = predict(fit)
      logwt = logwt - max(logwt)
      wt = exp(logwt)
      wt/sum(wt)
    }
    
    sumlog <- function(x) { sum(log(x), na.rm = T) }
    
    inits = foreach(s = 1:starts) %dopar% {
      for (i in 1:nfeatures) {
        temp.01 <- temp.02 <- gamma.j.k[, i]
        temp.1 <- unique(na.omit(temp.01))
        temp.2 <- p.gamma.k.m[[i]][s,]
        temp.3 <- p.gamma.k.u[[i]][s,]
        for (j in 1:length(temp.1)) {
          temp.01[temp.01 == temp.1[j]] <- temp.2[j]
          temp.02[temp.02 == temp.1[j]] <- temp.3[j]
        }
        p.gamma.k.j.m[i, ] <- temp.01
        p.gamma.k.j.u[i, ] <- temp.02
      }
      
      p.gamma.j.m <- as.matrix((apply(p.gamma.k.j.m, 2, sumlog)))
      p.gamma.j.m <- exp(p.gamma.j.m)
      
      p.gamma.j.u <- as.matrix((apply(p.gamma.k.j.u, 2, sumlog)))
      p.gamma.j.u <- exp(p.gamma.j.u)
      
      #First E step
      ## ------
      ## E-Step:
      ## ------
      
      log.prod <- log(p.gamma.j.m) + log(p.m[s])
      max.log.prod <- max(log.prod)
      
      log.sum <- logxpy(log(p.gamma.j.m) + log(p.m[s]), log(p.gamma.j.u) + log(p.u[s]))
      zeta.j <- exp(log.prod - max.log.prod)/(exp(log.sum - max.log.prod))
      
      for (iter in 1:st.iters) {
        ## --------
        ## M-step :
        ## --------
        ## get theta.m and theta.u
        num.prod <- n.j * zeta.j
        p.m. <- sum(num.prod)/sum(n.j)
        p.u. <- 1 - p.m.
   
        matches <- glm(count ~ ., data = data.frame(count = ((zeta.j * n.j) + c), factors), family = "quasipoisson")
        
        non.matches <- glm(count ~ .*., data = data.frame(count = ((1 - zeta.j) * n.j + c), factors), family = "quasipoisson")
        
        ## Predict & renormalization fn as in Murray 2017
        
        p.gamma.j.m = as.matrix(g(matches))
        p.gamma.j.u = as.matrix(g(non.matches))
        
        ## --------
        ## E-step :
        ## --------
        
        log.prod <- log(p.gamma.j.m) + log(p.m.)
        max.log.prod <- max(log.prod)
        
        log.sum <- logxpy(log(p.gamma.j.m) + log(p.m.), log(p.gamma.j.u) + log(p.u.))
        zeta.j <- exp(log.prod - max.log.prod)/(exp(log.sum - max.log.prod))
        
      }
      
      LL = dmultinom_robust(n.j * zeta.j, prob = p.gamma.j.m, log = T) + dmultinom_robust(n.j * (1-zeta.j), prob = p.gamma.j.u ,log = T)
      
      list(zeta.j = zeta.j, p.m = p.m., p.u = p.u., 
           BIC = -2 * LL + (1 + length(coef(matches)-1) + length(coef(non.matches)) - 1) * log(sum(n.j)))
    }
    stopCluster(cl)
    
    #averaging weights
    w = sapply(inits,'[[','BIC')
    w = w - min(w)
    w = exp(-0.5 * w) / sum(exp(-0.5 *w))
    
    count <- 1
    warn.once <- 1
    delta <- 1
    
    zeta.j = do.call(cbind,list(x = sapply(inits,'[[','zeta.j')))
    
    zeta.j = apply(zeta.j,1,weighted.mean,w = w)
    
    matches <- glm(count ~ ., data = data.frame(count = ((zeta.j * n.j) + c), factors),
                   family = "quasipoisson")
    
    non.matches <- glm(count ~ .*., data = data.frame(count = ((1 - zeta.j) * n.j + c), factors),
                       family = "quasipoisson")
    
    p.gamma.j.m = as.matrix(g(matches))
    p.gamma.j.u = as.matrix(g(non.matches))
    num.prod <- n.j * zeta.j
    p.m <- sum(num.prod)/sum(n.j)
    p.u <- 1 - p.m
    ## The EM Algorithm presented in the paper starts here:
    #now stop if average training likelihood fails to improve after M step
  
    
    while (abs(delta) >= tol) {
      
      if((count %% 100) == 0) {
        cat("Iteration number", count, "\n")
        cat("Maximum difference in log-likelihood =", round(delta, 4), "\n")
      }
      
      ## Old Paramters
      p.old <- c(p.m, p.u, unlist(p.gamma.j.m), unlist(p.gamma.j.u))
      
      ## ------
      ## E-Step:
      ## ------
      
      log.prod <- log(p.gamma.j.m) + log(p.m)
      max.log.prod <- max(log.prod)
      
      log.sum <- logxpy(log(p.gamma.j.m) + log(p.m), log(p.gamma.j.u) + log(p.u))
      zeta.j <- exp(log.prod - max.log.prod)/(exp(log.sum - max.log.prod))
      
      ## --------
      ## M-step :
      ## --------
      num.prod <- n.j * zeta.j
      p.m <- sum(num.prod)/sum(n.j)
      p.u <- 1 - p.m
      
      ## get theta.m and theta.u
      matches <- glm(count ~ ., data = data.frame(count = ((zeta.j * n.j) + c), factors),
                     family = "quasipoisson")
      
      non.matches <- glm(count ~ .*., data = data.frame(count = ((1 - zeta.j) * n.j + c), factors),
                         family = "quasipoisson")
      
      ## Predict & renormalization fn as in Murray 2017
      
      p.gamma.j.m = as.matrix(g(matches))
      p.gamma.j.u = as.matrix(g(non.matches))
      
      ## Updated parameters:
      p.new <- c(p.m, p.u, unlist(p.gamma.j.m), unlist(p.gamma.j.u))
      
      if(p.m < 1e-13 & warn.once == 1) {
        warning("The overall probability of finding a match is too small. Increasing the amount of overlap between the datasets might help, see e.g., clusterMatch()")
        warn.once <- 0
      }
      
      ## Max difference between the updated and old parameters:
      delta <- max(abs(p.new - p.old))
      count <- count + 1
      
      if(count > iter.max) {
        warning("The EM algorithm has run for the specified number of iterations but has not converged yet.")
        break
      }
    }
   
    weights <- log(p.gamma.j.m) - log(p.gamma.j.u)
    
    data.w <- cbind(patterns, weights, p.gamma.j.m, p.gamma.j.u)
    nc <- ncol(data.w)
    colnames(data.w)[nc-2] <- "weights"
    colnames(data.w)[nc-1] <- "p.gamma.j.m"
    colnames(data.w)[nc] <- "p.gamma.j.u"
    
    inf <- which(data.w == Inf, arr.ind = T)
    ninf <- which(data.w == -Inf, arr.ind = T)
    
    data.w[inf[, 1], unique(inf[, 2])] <- 150
    data.w[ninf[, 1], unique(ninf[, 2])] <- -150

    if(!is.null(varnames)){
        output <- list("zeta.j"= zeta.j,"p.m"= p.m, "p.u" = p.u, 
                       "p.gamma.j.m" = p.gamma.j.m, "p.gamma.j.u" = p.gamma.j.u, "patterns.w" = data.w, "iter.converge" = count,
                       "nobs.a" = nobs.a, "nobs.b" = nobs.b, "varnames" = varnames,match_mod = matches,
                       nonmatch_mod = non.matches)
    }else{
        output <- list("zeta.j"= zeta.j,"p.m"= p.m, "p.u" = p.u, 
                       "p.gamma.j.m" = p.gamma.j.m, "p.gamma.j.u" = p.gamma.j.u, "patterns.w" = data.w, "iter.converge" = count,
                       "nobs.a" = nobs.a, "nobs.b" = nobs.b, "varnames" = paste0("gamma.", 1:nfeatures),match_mod = matches,
                       nonmatch_mod = non.matches)
    }

    
    class(output) <- c("fastLink", "fastLink.EM")
    
    return(output)
}

dmultinom_robust = function (x, size = NULL, prob, log = FALSE){
  K <- length(prob)
  if (length(x) != K) 
    stop("x[] and prob[] must be equal length vectors.")
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 
      0) 
    stop("probabilities must be finite, non-negative and not all 0")
  prob <- prob/s
  x <- round(x)
  if (any(x < 0)) 
    stop("'x' must be non-negative")
  N <- sum(x)
  if (is.null(size)) 
    size <- N
  else if (size != N) 
    stop("size != sum(x), i.e. one is wrong")
  i0 <- prob == 0
  if (any(i0)) {
    if (any(x[i0] != 0)) 
      return(if (log) -Inf else 0)
    if (all(i0)) 
      return(if (log) 0 else 1)
    x <- x[!i0]
    prob <- prob[!i0]
  }
  r <- lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
  if (log) 
    r
  else exp(r)
}
