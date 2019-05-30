#' gammaCKpar
#'
#' Field comparisons for string variables. Three possible agreement patterns are considered:
#' 0 total disagreement, 1 partial agreement, 2 agreement.
#' The distance between strings is calculated using a Jaro-Winkler distance.
#'
#' @usage gammaCKpar(matAp, matBp, n.cores, cut.a, cut.p, method, w)
#'
#' @param matAp vector storing the comparison field in data set 1
#' @param matBp vector storing the comparison field in data set 2
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#' @param cut.a Lower bound for full match, ranging between 0 and 1. Default is 0.92
#' @param cut.p Lower bound for partial match, ranging between 0 and 1. Default is 0.88
#' @param dedupe Logical. Whether internal linkage or linkage between two frames is being performed.
#' @param transform An optional function returning transformed strings or a list of transformed strings. 
#' First argument should be the strings to be transformed
#' @param transform.args An optional list of arguments to the transform function
#' @param method String distance method, options are: "jw" Jaro-Winkler (Default), "jaro" Jaro, and "lv" Edit. 
#' May also be passed as a custom function with argument list stringdist.args. The first two arguments of 
#' any custom function provided should be the strings to be compared
#' @param method.args List of arguments for custom function provided to method
#' @param w Parameter that describes the importance of the first characters of a string (only needed if method = "jw"). Default is .10
#'
#' @return \code{gammaCKpar} returns a list with the indices corresponding to each
#' matching pattern, which can be fed directly into \code{tableCounts} and \code{matchesLink}.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield <benfifield@gmail.com>, and Kosuke Imai
#'
#' @examples
#' \dontrun{
#' g1 <- gammaCKpar(dfA$firstname, dfB$lastname)
#' }
#'
#' @export
#' @importFrom stringdist stringdistmatrix
## ------------------------
## gammaCKpar: Now it takes values 0, 1, 2
## This function applies gamma.k
## in parallel
## ------------------------

gammaCKpar <- function(matAp, matBp, n.cores = NULL, cut.a = 0.92, cut.p = 0.88, dedupe = F, transform = NULL,
                       transform.args = NULL, method = "jw",method.args=NULL,
                       w = .10) {

    if(any(class(matAp) %in% c("tbl_df", "data.table"))){
        matAp <- as.data.frame(matAp)[,1]
    }
    if(any(class(matBp) %in% c("tbl_df", "data.table"))){
        matBp <- as.data.frame(matBp)[,1]
    }
    
    matAp[matAp == ""] <- NA
    matBp[matBp == ""] <- NA

    if(sum(is.na(matAp)) == length(matAp) | length(unique(matAp)) == 1){
        cat("WARNING: You have no variation in this variable, or all observations are missing in dataset A.\n")
    }
    if(sum(is.na(matBp)) == length(matBp) | length(unique(matBp)) == 1){
        cat("WARNING: You have no variation in this variable, or all observations are missing in dataset B.\n")
    }
    
    if( !is.function(method)){
      if(!(method %in% c("jw", "jaro", "lv")))
        stop("Invalid string distance method. Method should be one of 'jw', 'jaro', or 'lv'.")
      if(method == "jw" & !is.null(w)){
        if(w < 0 | w > 0.25){
          stop("Invalid value provided for w. Remember, w in [0, 0.25].")
        }
      }
    } 
        
    if(is.function(method) && is.null(method.args))
      stop("You must provide a list of arguments if using a custom string comparison function")

    if(is.null(n.cores)) {
        n.cores <- detectCores() - 1
    }

    matrix.1 <- as.matrix(as.character(matAp))
    matrix.2 <- as.matrix(as.character(matBp))

    matrix.1[is.na(matrix.1)] <- "9999"
    matrix.2[is.na(matrix.2)] <- "9998"

    u.values.1 <- unique(matrix.1)
    u.values.2 <- unique(matrix.2)
    
    if(!is.null(transform)){
      u.trans.1 = do.call(transform, c(list(u.values.1), transform.args))
      u.trans.2 = do.call(transform, c(list(u.values.2), transform.args))
    }
    
    exports = pkgs = character(0)
    
    if(is.function(method)){
      depsm = find_dependencies(method)
      depst = find_dependencies(transform)
      
      exports = unique(c(depsm$depends$calls[sapply(depsm$depends$pkgs,function(x) '.GlobalEnv' %in% x)],
                         depst$depends$calls[sapply(depst$depends$pkgs,function(x) '.GlobalEnv' %in% x)]))
      
      pkgs = unique(c(unlist(depsm$depends$pkgs),unlist(depst$depends$pkgs)))
      pkgs = pkgs[!(pkgs %in% c('base','.GlobalEnv'))]
    }
    
    if(dedupe){
      n.slices <- max(round(length(u.values.1)/(300), 0), 1) 
      limit = round(seq(0,length(u.values.1),length(u.values.1)/n.slices),0)
      
      temp.1 <- temp.2 <- list()
      
      n.cores2 <- min(n.cores, n.slices*(n.slices + 1)/2)
      temp = list()
      
      for(i in 1:n.slices) {
        if(is.null(transform)){
          temp[[i]] <- list(u.values.1[(limit[i]+1):limit[i+1]], limit[i])
        }else{
          temp[[i]] <- list(lapply(u.trans.1, '[',(limit[i]+1):limit[i+1]), limit[i])
        }
      }
      
      stringvec <- function(m, y, cut, strdist = method, p1 = w, identical) {
        if(!is.list(m[[1]])){
          x <- as.matrix(m[[1]])
          e <- as.matrix(y[[1]])
          nr = nrow(e)
        }else{
          x = m[[1]]
          e = y[[1]]
          nr = length(e[[1]])
        }     
        
        if(is.function(strdist)){
          if(identical){
            t <- do.call(strdist,c(list(e), method.args))
            # diag(t) = sapply(1:length(e[[1]]), function(i){
            #   inp = lapply(e,'[',i)
            #   do.call(strdist,c(list(inp),list(inp), method.args))
            #})
            diag(t) = cut[[2]]+1
          }else{
            t <- do.call(strdist,c(list(e, x), method.args))
          }
          t[ t < cut[[2]] ] <- 0
          t <- Matrix(t, sparse = T)
        }else{
          if(strdist == "jw") {
            if(identical){
              t <- 1 - stringdistmatrix(a = e, method = "jw", p = p1, nthread = 1)
              t = as.matrix(t) * lower.tri(t)
              
              #diag(t) = 1 - stringdist(e , e,  method = "jw", p = p1, nthread = 1)
              diag(t) = cut[[2]]+1
            }else{
              t <- 1 - stringdistmatrix(e, x, method = "jw", p = p1, nthread = 1)
            }
            t[ t < cut[[2]] ] <- 0
            t <- Matrix(t, sparse = T)
          }
          
          if(strdist == "jaro") {
            if(identical){
              t <- 1 - stringdistmatrix(a=e, method = "jw", nthread = 1)  
              t = as.matrix(t) * lower.tri(t)
              
              #diag(t) = 1 - stringdist(e, e, method = 'jw', nthread = 1)
              diag(t) = cut[[2]]+1
            }else{
              t <- 1 - stringdistmatrix(e, x, method = "jw", nthread = 1)  
            }
            t[ t < cut[[2]] ] <- 0
            t <- Matrix(t, sparse = T)
          }
          
          if(strdist == "lv") {
            if(identical){
              t <- stringdistmatrix(e, method = method, nthread = 1)
              t = as.matrix(t) * lower.tri(t)
              
              #diag(t) = stringdist(e, e, method = 'lv', nthread = 1)
              diag(t) = 0
            }else{
              t <- stringdistmatrix(e, x, method = method, nthread = 1)
            }
            
            t.1 <- nchar(as.matrix(e))
            t.2 <- nchar(as.matrix(x))
            o <- t(apply(t.1, 1, function(w){ ifelse(w >= t.2, w, t.2)}))
            t <- 1 - t * (1/o)
            t[ t < cut[[2]] ] <- 0
            t <- Matrix(t, sparse = T)
          }
        }
        
        t@x[t@x >= cut[1]] <- 2
        t@x[t@x >= cut[2] & t@x < cut[1]] <- 1; gc()
        slice.1 <- m[[2]]
        slice.2 <- y[[2]]
        indexes.2 <- which(t == 2, arr.ind = T)
        indexes.2[, 1] <- indexes.2[, 1] + slice.2
        indexes.2[, 2] <- indexes.2[, 2] + slice.1
        indexes.1 <- which(t == 1, arr.ind = T)
        indexes.1[, 1] <- indexes.1[, 1] + slice.2
        indexes.1[, 2] <- indexes.1[, 2] + slice.1
        list(indexes.2, indexes.1)
      }
      
      if (n.cores2 == 1) '%oper%' <- foreach::'%do%'
      else { 
        '%oper%' <- foreach::'%dopar%'
        cl <- snow::makeCluster(n.cores2,type = 'SOCK')
        registerDoSNOW(cl)
        on.exit(stopCluster(cl))
      }
      
      do <- expand.grid.jc(1:n.slices,1:n.slices)
      drop = do[,2] < do[,1]
      do = matrix(do[!drop,],ncol=2)
      
      pb = txtProgressBar(0,nrow(do),style = 1)
      
      progress <- function(n) setTxtProgressBar(pb, n)
      
      opts <- list(progress = progress)
      
      temp.f <- foreach(i = 1:nrow(do), .packages = c("stringdist", "Matrix",pkgs),.export = exports,
                        .options.snow=opts) %oper% { 
                          r1 <- do[i, 1]
                          r2 <- do[i, 2]
                          #if(i %% 100 == 0) setTxtProgressBar(pb,i)
                          stringvec(temp[[r1]], temp[[r2]], c(cut.a, cut.p),
                                    identical = r1==r2)
                        }
      close(pb)
      
    }else{
      n.slices1 <- max(round(length(u.values.1)/(300), 0), 1) 
      n.slices2 <- max(round(length(u.values.2)/(300), 0), 1) 
      
      limit.1 <- round(quantile((0:nrow(u.values.2)), p = seq(0, 1, 1/n.slices2)), 0)
      limit.2 <- round(quantile((0:nrow(u.values.1)), p = seq(0, 1, 1/n.slices1)), 0)
      
      n.cores <- min(n.cores, n.slices1 * n.slices2)
      
      temp.1 <- temp.2 <- list()
      
      for(i in 1:n.slices2) {
        if(is.null(transform)){
          temp.1[[i]] <- list(u.values.2[(limit.1[i]+1):limit.1[i+1]], limit.1[i])
        }else{
          temp.1[[i]] <- list(lapply(u.trans.2, '[',(limit.1[i]+1):limit.1[i+1]), limit.1[i])
        }
      }
      
      for(i in 1:n.slices1) {
        if(is.null(transform)){
          temp.2[[i]] <- list(u.values.1[(limit.2[i]+1):limit.2[i+1]], limit.2[i])
        }else{
          temp.2[[i]] <- list(lapply(u.trans.1,'[',(limit.2[i]+1):limit.2[i+1]), limit.2[i])
        }
      }
      
      stringvec <- function(m, y, cut, strdist = method, p1 = w) {
        if(!is.list(m[[1]])){
          x <- as.matrix(m[[1]])
          e <- as.matrix(y[[1]])
          nr = nrow(e)
        }else{
          x = m[[1]]
          e = y[[1]]
          nr = length(e[[1]])
        } 
        
        if(is.function(strdist)){
          t <- do.call(strdist,c(list(e, x), method.args))
          t[ t < cut[[2]] ] <- 0
          t <- Matrix(t, sparse = T, nrow = nr)
        }else{
          if(strdist == "jw") {
            t <- 1 - stringdistmatrix(e, x, method = "jw", p = p1, nthread = 1)
            t[ t < cut[[2]] ] <- 0
            t <- Matrix(t, sparse = T)
          }
          
          if(strdist == "jaro") {
            t <- 1 - stringdistmatrix(e, x, method = "jw", nthread = 1)
            t[ t < cut[[2]] ] <- 0
            t <- Matrix(t, sparse = T)
          }
          
          if(strdist == "lv") {
            t <- stringdistmatrix(e, x, method = method, nthread = 1)
            t.1 <- nchar(as.matrix(e))
            t.2 <- nchar(as.matrix(x))
            o <- t(apply(t.1, 1, function(w){ ifelse(w >= t.2, w, t.2)}))
            t <- 1 - t * (1/o)
            t[ t < cut[[2]] ] <- 0
            t <- Matrix(t, sparse = T)
          }
        }
        
        t@x[t@x >= cut[1]] <- 2
        t@x[t@x >= cut[2] & t@x < cut[1]] <- 1; gc()
        slice.1 <- m[[2]]
        slice.2 <- y[[2]]
        indexes.2 <- which(t == 2, arr.ind = T)
        indexes.2[, 1] <- indexes.2[, 1] + slice.2
        indexes.2[, 2] <- indexes.2[, 2] + slice.1
        indexes.1 <- which(t == 1, arr.ind = T)
        indexes.1[, 1] <- indexes.1[, 1] + slice.2
        indexes.1[, 2] <- indexes.1[, 2] + slice.1
        list(indexes.2, indexes.1)
      }
      
      do <- expand.grid(1:n.slices2, 1:n.slices1)
      
      if (n.cores == 1) '%oper%' <- foreach::'%do%'
      else { 
        '%oper%' <- foreach::'%dopar%'
        cl <- snow::makeCluster(n.cores,type = 'SOCK')
        registerDoSNOW(cl)
        on.exit(stopCluster(cl))
      }
      
      pb = txtProgressBar(0,nrow(do),style = 1)
      
      progress <- function(n) setTxtProgressBar(pb, n)
      
      opts <- list(progress = progress)
      
      temp.f <- foreach(i = 1:nrow(do), .packages = c("stringdist", "Matrix",pkgs),.export = exports,
                        .options.snow=opts) %oper% { 
                          r1 <- do[i, 1]
                          r2 <- do[i, 2]
                          #if(i %% 100 == 0) setTxtProgressBar(pb,value = i)
                          stringvec(temp.1[[r1]], temp.2[[r2]], c(cut.a, cut.p))
                        }
      close(pb)
    }
    
    
    gc()

    reshape2 <- function(s) { s[[1]] }
    reshape1 <- function(s) { s[[2]] }
    temp.2 <- lapply(temp.f, reshape2)
    temp.1 <- lapply(temp.f, reshape1)

    indexes.2 <- do.call('rbind', temp.2)
    indexes.1 <- do.call('rbind', temp.1)
    
    .identical = indexes.2[,1] == indexes.2[,2]
    
    ht1 <- new.env(hash=TRUE)
    ht2 <- new.env(hash=TRUE)

    n.values.2 <- as.matrix(cbind(u.values.1[indexes.2[, 1]], u.values.2[indexes.2[, 2]]))
    n.values.1 <- as.matrix(cbind(u.values.1[indexes.1[, 1]], u.values.2[indexes.1[, 2]]))

    matches.2 <- lapply(seq_len(nrow(n.values.2)), function(i) n.values.2[i, ])
    matches.1 <- lapply(seq_len(nrow(n.values.1)), function(i) n.values.1[i, ])

    if(Sys.info()[['sysname']] == 'Windows') {
        if (n.cores == 1) '%oper%' <- foreach::'%do%'
        else { 
            '%oper%' <- foreach::'%dopar%'
            cl <- makeCluster(n.cores)
            registerDoParallel(cl)
        }
        if(length(matches.2) > 0) {
            final.list2 <- foreach(i = 1:length(matches.2)) %oper% {
            ht1 <- which(matrix.1 == matches.2[[i]][[1]]); ht2 <- which(matrix.2 == matches.2[[i]][[2]])
            list(ht1, ht2)
            }
        }
        if(length(matches.1) > 0) {
                final.list1 <- foreach(i = 1:length(matches.1)) %oper% {
                ht1 <- which(matrix.1 == matches.1[[i]][[1]]); ht2 <- which(matrix.2 == matches.1[[i]][[2]])
                list(ht1, ht2)
            }
        }
        if(n.cores > 1){
            stopCluster(cl)
        }
    } else {
      no_cores <- n.cores
    	final.list2 <- mclapply(matches.2, function(s){
            ht1 <- which(matrix.1 == s[1]); ht2 <- which(matrix.2 == s[2]);
            list(ht1, ht2) }, mc.cores = getOption("mc.cores", no_cores))

    	final.list1 <- mclapply(matches.1, function(s){
            ht1 <- which(matrix.1 == s[1]); ht2 <- which(matrix.2 == s[2]);
            list(ht1, ht2) }, mc.cores = getOption("mc.cores", no_cores))
    }

    if(length(matches.2) == 0){ 
      final.list2 <- list()
      warning("There are no identical (or nearly identical) matches. We suggest either changing the value of cut.p") 
    }
    
    if(length(matches.1) == 0){ 
    	final.list1 <- list()
    	warning("There are no partial matches. We suggest either changing the value of cut.p or using gammaCK2par() instead") 
    }
    
    na.list <- list()
    na.list[[1]] <- which(matrix.1 == "9999")
    na.list[[2]] <- which(matrix.2 == "9998")

    out <- list()
    out[["matches2"]] <- final.list2
    out[["matches1"]] <- final.list1
    out[["nas"]] <- na.list
    out[['.identical']] = .identical
    class(out) <- c("fastLink", "gammaCKpar")

    return(out)
}


## ------------------------
## End of gammaCKpar
## ------------------------

