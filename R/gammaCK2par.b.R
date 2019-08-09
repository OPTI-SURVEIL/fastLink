#' gammaCK2par.b
#'
#' Field comparisons for string variables for blocked data. Two possible agreement patterns are considered:
#' 0 total disagreement, 2 agreement.
#' Distance between strings may be calculated using various metrics, including custom functions
#'
#' @usage gammaCK2par.b(matAp, matBp, blocklist, dedupe.blocks, n.cores, cut.a, transform, transform.args, method, method.args, w)
#'
#' @param matAp vector storing the comparison field in data set 1
#' @param matBp vector storing the comparison field in data set 2
#' @param blocklist list of blocking indices
#' @param dedupe.blocks Logical indicators for whether internal linkage and duplication is taking place
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#' @param cut.a Lower bound for full match. Default is 0.92
#' @param transform An optional function returning transformed strings or a list of transformed strings. 
#' First argument should be the strings to be transformed
#' @param transform.args An optional list of arguments to the transform function
#' @param method String distance method, options are: "jw" Jaro-Winkler (Default), "jaro" Jaro, and "lv" Edit. 
#' May also be passed as a custom function with argument list stringdist.args. The first two arguments of 
#' any custom function provided should be the strings to be compared
#' @param method.args List of arguments for custom function provided to method
#' @param w Parameter that describes the importance of the first characters of a string (only needed if method = "jw"). Default is .10
#' @return \code{gammaCK2par.b} returns a list of lists (for each block) with the indices corresponding to each
#' matching pattern, which can be fed directly into \code{tableCounts.b} and \code{matchesLink.b}.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield <benfifield@gmail.com>, and Kosuke Imai
#' 
#' @export

gammaCK2par.b <- function(matAp, matBp, blocklist, n.cores = NULL, cut.a = 0.92, transform = NULL, 
                          transform.args, method = 'jw', w = 0.1, method.args, dedupe.blocks) {
  varnm = names(matAp)
  if(any(class(matAp) %in% c("tbl_df", "data.table",'data.frame'))){
    matAp <- as.data.frame(matAp)[,1]
  }
  if(any(class(matBp) %in% c("tbl_df", "data.table",'data.frame'))){
    matBp <- as.data.frame(matBp)[,1]
  }
  
  matAp[matAp == ""] <- NA
  matBp[matBp == ""] <- NA
  
  if(sum(is.na(matAp)) == length(matAp) | length(unique(matAp)) == 1){
    cat("WARNING: You have no variation in variable '", varnm,"', or all observations are missing in dataset A.\n",sep = '')
  }
  if(sum(is.na(matBp)) == length(matBp) | length(unique(matBp)) == 1){
    cat("WARNING: You have no variation in variable '", varnm,"', or all observations are missing in dataset B.\n",sep = '')
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
  
  matrix.1[is.na(matrix.1)] <- "1234MF"
  matrix.2[is.na(matrix.2)] <- "9876ES"
  
  #progress bar for blocklist level operations
  pb1 = txtProgressBar(0, length(blocklist),style = 1,char = '+')
  progress1 = function(n) setTxtProgressBar(pb1, n)
  opts1 = list(progress = progress1)
  
  cl1 = makeCluster(n.cores)
  registerDoSNOW(cl1)
  
  if(!is.null(transform)){
    u.values.1 <- unique(matrix.1)
    u.trans.1 = as.data.frame(do.call(transform, c(list(u.values.1), transform.args)),stringsAsFactors = F)
    tr_names = colnames(u.trans.1)
    u.trans.1 = lapply(1:nrow(u.trans.1), function(i) u.trans.1[i,])
    names(u.trans.1) = u.values.1; u.trans.1 = list2env(u.trans.1)
    if(identical(matAp,matBp)){
      u.trans.2 = u.trans.1
    }else{
      u.values.2 <- unique(matrix.2)
      u.trans.2 = as.data.frame(do.call(transform, c(list(u.values.2), transform.args)),stringsAsFactors = F)
      u.trans.2 = lapply(1:nrow(u.trans.2), function(i) u.trans.2[i,])
      names(u.trans.2) = u.values.2; u.trans.2 = list2env(u.trans.2)
    } 
  }
  
  clusterExport(cl1, c('matrix.1','matrix.2','blocklist','dedupe.blocks','cut.a','u.trans.1','u.trans.2','method','method.args'))
  
  #find necessary packages and objects to export
  exports = pkgs = character(0)
  if(is.function(method)){
    depsm = fastLink:::find_dependencies(method)
    #depst = find_dependencies(transform)
    
    exports = unique(c(depsm$depends$calls[sapply(depsm$depends$pkgs,function(x) '.GlobalEnv' %in% x)]))#,
    #depst$depends$calls[sapply(depst$depends$pkgs,function(x) '.GlobalEnv' %in% x)]))
    
    pkgs = unique(c(unlist(depsm$depends$pkgs)))#,unlist(depst$depends$pkgs)))
    pkgs = pkgs[!(pkgs %in% c('base','.GlobalEnv'))]
  }
  
  clusterExport(cl1, exports)
  
  #find number of comparisons in each block
  cat('finding number of unique pairs of values of variable *', varnm, '* across all blocks \n', sep = '')
  
  compvals = foreach(i = 1:length(blocklist), .options.snow = opts1, .combine = c) %dopar% {
    b = blocklist[[i]]
    u.1 = unique(matrix.1[b[[1]]]); u.2 = unique(matrix.2[b[[2]]])
    
    as.numeric(length(u.1)) * as.numeric(length(u.2)) / 300^2
  }
  
  inpar = compvals <= 1
  
  #define comparison routine
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
        diag(t) = cut + 1
      }else{
        t <- do.call(strdist,c(list(e, x), method.args))
      }
      t[ t < cut ] <- 0
      t <- Matrix(t, sparse = T)
    }else{
      if(strdist == "jw") {
        if(identical){
          t <- 1 - stringdistmatrix(a = e, method = "jw", p = p1, nthread = 1)
          t = as.matrix(t) * lower.tri(t)
          diag(t) = cut + 1
        }else{
          t <- 1 - stringdistmatrix(e, x, method = "jw", p = p1, nthread = 1)
        }
        t[ t < cut ] <- 0
        t <- Matrix(t, sparse = T)
      }
      
      if(strdist == "jaro") {
        if(identical){
          t <- 1 - stringdistmatrix(a=e, method = "jw", nthread = 1)  
          t = as.matrix(t) * lower.tri(t)
          diag(t) = cut + 1
        }else{
          t <- 1 - stringdistmatrix(e, x, method = "jw", nthread = 1)  
        }
        t[ t < cut ] <- 0
        t <- Matrix(t, sparse = T)
      }
      
      if(strdist == "lv") {
        if(identical){
          t <- stringdistmatrix(e, method = 'lv', nthread = 1)
          t = as.matrix(t) * lower.tri(t)
          diag(t) = 0
        }else{
          t <- stringdistmatrix(e, x, method = 'lv', nthread = 1)
        }
        
        t.1 <- nchar(as.matrix(e))
        t.2 <- nchar(as.matrix(x))
        o <- t(apply(t.1, 1, function(w){ ifelse(w >= t.2, w, t.2)}))
        t <- 1 - t * (1/o)
        t[ t < cut ] <- 0
        t <- Matrix(t, sparse = T)
      }
    }
    
    t@x[t@x >= cut] <- 2; gc()       	
    slice.1 <- m[[2]]
    slice.2 <- y[[2]]
    indexes.2 <- which(t == 2, arr.ind = T)
    indexes.2[, 1] <- indexes.2[, 1] + slice.2
    indexes.2[, 2] <- indexes.2[, 2] + slice.1
    list(indexes.2)
  }
  clusterExport(cl1, c('stringvec','inpar'))
  out = vector('list',length(blocklist))
  #perform comparisons on all small blocks
  
  if(sum(inpar)>0){
    cat('\nparallelizing over small blocks\n', sep = '')
  
    out[inpar] = foreach(i = 1:sum(inpar), 
                         .options.snow = opts1,.packages = c('Matrix','fastLink','RcppAlgos',pkgs)) %dopar% {
      b = blocklist[inpar][[i]]
      ddp = dedupe.blocks[inpar][i]
      m.1 = matrix.1[b[[1]],]; m.2 = matrix.2[b[[2]],]
      
      u.values.1 <- unique(m.1)
      u.values.2 <- unique(m.2)
      
      if(!is.null(transform)){
        temp.1 = mget(u.values.1,u.trans.1)
        temp.1 = as.list(data.frame(do.call(rbind,temp.1),stringsAsFactors = F)); names(temp.1) = tr_names
        temp.2 = mget(u.values.2,u.trans.2)
        temp.2 = as.list(data.frame(do.call(rbind,temp.2),stringsAsFactors = F)); names(temp.2) = tr_names
      }else{
        temp.1 = u.values.1
        temp.2 = u.values.2
      }
      
      if(ddp & length(u.values.1) == 1){
        final.list2 <- list(1:length(m.1), 1:length(m.2))
        final.list1 <- list()
        .identical = T
      }else{
        temp.f <- stringvec(list(temp.1,0), list(temp.2,0), cut.a, identical = ddp)
        
        indexes.2 <- temp.f[[1]]
        .identical = indexes.2[,1] == indexes.2[,2]
        
        n.values.2 <- as.matrix(cbind(u.values.1[indexes.2[, 1]], u.values.2[indexes.2[, 2]]))
        
        if(sum(n.values.2 == "1234MF") > 0) {
          t1 <- which(n.values.2 == "1234MF", arr.ind = T)[1]
          n.values.2 <- n.values.2[-t1, ]; rm(t1)
        }
        
        if(sum(n.values.2 == "9876ES") > 0) {
          t1 <- which(n.values.2 == "9876ES", arr.ind = T)[1]
          n.values.2 <- n.values.2[-t1, ]; rm(t1)
        }
        
        matches.2 <- lapply(seq_len(nrow(n.values.2)), function(i) n.values.2[i, ])
        
        values2store.1 = unique(n.values.2[,1])
        values2store.2 = unique(n.values.2[,2])
        
        ht.1 = lapply(values2store.1, function(x) which(m.1 == x))
        names(ht.1) = values2store.1
        ht.1 = list2env(ht.1)
        
        ht.2 = lapply(values2store.2, function(x) which(m.2 == x))
        names(ht.2) = values2store.2
        ht.2 = list2env(ht.2)
        
        final.list2 <- lapply(matches.2, function(x){
          ht1 <- ht.1[[as.character(x[[1]])]]; ht2 <- ht.2[[as.character(x[[2]])]]
          list(ht1, ht2)
        })
      }
      
      na.list <- list()
      na.list[[1]] <- which(m.1 == "1234MF")
      na.list[[2]] <- which(m.2 == "9876ES")  
      
      out = list(matches2 = final.list2, nas = na.list, .identical = .identical)
      class(out) <- c("fastLink", "gammaCK2par")
      out
    }
  }
  if(sum(!inpar) > 0){
    cat('\niterating over larger blocks\n', sep = '')
    
    for(j in (sum(inpar)+1):length(blocklist)){
      b = blocklist[!inpar][[j-sum(inpar)]]
      ddp = dedupe.blocks[!inpar][j-sum(inpar)]
      m.1 = matrix.1[b[[1]],]; m.2 = matrix.2[b[[2]],]
      
      u.values.1 <- unique(m.1)
      u.values.2 <- unique(m.2)
      
      if(!is.null(transform)){
        tmp.1 = mget(u.values.1,u.trans.1)
        tmp.1 = as.list(data.frame(do.call(rbind,tmp.1),stringsAsFactors = F)); names(tmp.1) = tr_names
        tmp.2 = mget(u.values.2,u.trans.2)
        tmp.2 = as.list(data.frame(do.call(rbind,tmp.2),stringsAsFactors = F)); names(tmp.2) = tr_names
      }else{
        tmp.1 = u.values.1
        tmp.2 = u.values.2
      }
      
      if(ddp){
        n.slices <- max(round(length(u.values.1)/(300), 0), 1) 
        limit = round(seq(0,length(u.values.1),length(u.values.1)/n.slices),0)
        
        temp = list()
        
        for(i in 1:n.slices) {
          if(is.null(transform)){
            temp[[i]] <- list(u.values.1[(limit[i]+1):limit[i+1]], limit[i])
          }else{
            temp[[i]] <- list(lapply(tmp.1, '[',(limit[i]+1):limit[i+1]), limit[i])
          }
        }
        
        do <- expand.grid.jc(1:n.slices,1:n.slices)
        drop = do[,2] < do[,1]
        do = matrix(do[!drop,],ncol=2)
      }else{
        n.slices1 <- max(round(length(u.values.1)/(300), 0), 1) 
        n.slices2 <- max(round(length(u.values.2)/(300), 0), 1) 
        
        limit.1 <- round(quantile((0:length(u.values.2)), p = seq(0, 1, 1/n.slices2)), 0)
        limit.2 <- round(quantile((0:length(u.values.1)), p = seq(0, 1, 1/n.slices1)), 0)
        
        temp.1 <- temp.2 <- list()
        
        for(i in 1:n.slices2) {
          if(is.null(transform)){
            temp.1[[i]] <- list(u.values.2[(limit.1[i]+1):limit.1[i+1]], limit.1[i])
          }else{
            temp.1[[i]] <- list(lapply(tmp.2, '[',(limit.1[i]+1):limit.1[i+1]), limit.1[i])
          }
        }
        
        for(i in 1:n.slices1) {
          if(is.null(transform)){
            temp.2[[i]] <- list(u.values.1[(limit.2[i]+1):limit.2[i+1]], limit.2[i])
          }else{
            temp.2[[i]] <- list(lapply(tmp.1,'[',(limit.2[i]+1):limit.2[i+1]), limit.2[i])
          }
        }
        
        do <- expand.grid(1:n.slices2, 1:n.slices1)
      }
      
      pb2 = txtProgressBar(0,nrow(do))
      
      progress2 = function(n) setTxtProgressBar(pb2,n)
      opts2 = list(progress = progress2)
      
      temp.f <- foreach(i = 1:nrow(do), .packages = c("Matrix",'fastLink','RcppAlgos',pkgs), .options.snow = opts2) %dopar% {
        r1 <- do[i, 1]
        r2 <- do[i, 2]
        
        stringvec(temp[[r1]], temp[[r2]], cut.a, identical = r1==r2 & ddp)
      }
      close(pb2)
      
      indexes.2 <- do.call('rbind', lapply(temp.f,'[[',1))
      .identical = indexes.2[,1] == indexes.2[,2]
      
      n.values.2 <- as.matrix(cbind(u.values.1[indexes.2[, 1]], u.values.2[indexes.2[, 2]]))
      
      if(sum(n.values.2 == "1234MF") > 0) {
        t1 <- which(n.values.2 == "1234MF", arr.ind = T)[1]
        n.values.2 <- n.values.2[-t1, ]; rm(t1)
      }
      
      if(sum(n.values.2 == "9876ES") > 0) {
        t1 <- which(n.values.2 == "9876ES", arr.ind = T)[1]
        n.values.2 <- n.values.2[-t1, ]; rm(t1)
      }
      
      matches.2 <- lapply(seq_len(nrow(n.values.2)), function(i) n.values.2[i, ])
      
      values2store.1 = unique(n.values.2[,1])
      values2store.2 = unique(n.values.2[,2])
      
      ht.1 = foreach(i = 1:length(values2store.1)) %dopar% {
        which(m.1 == values2store.1[i])
      }; names(ht.1) = values2store.1
      ht.1 = list2env(ht.1)
      
      ht.2 = foreach(i = 1:length(values2store.2)) %dopar% {
        which(m.2 == values2store.2[i])
      }; names(ht.2) = values2store.2
      ht.2 = list2env(ht.2)
      
      final.list2 <- foreach(i = 1:length(matches.2)) %dopar% {
        ht1 <- ht.1[[as.character(matches.2[[i]][[1]])]]; ht2 <- ht.2[[as.character(matches.2[[i]][[2]])]]
        list(ht1, ht2)
      }
      
      na.list <- list()
      na.list[[1]] <- which(matrix.1 == '1234MF')
      na.list[[2]] <- which(matrix.2 == '9876ES')
      out[!inpar][[j-sum(inpar)]] = list()
      out[!inpar][[j-sum(inpar)]][["matches2"]] <- final.list2
      out[!inpar][[j-sum(inpar)]][["nas"]] <- na.list
      out[!inpar][[j-sum(inpar)]][[".identical"]] <- .identical
      class(out[!inpar][[j-sum(inpar)]]) <- c("fastLink", "gammaCK2par")
      
      capture.output(setTxtProgressBar(pb1,1))
      setTxtProgressBar(pb1,j)
      gc()
    }
  }
      
  close(pb1)
  stopCluster(cl1)
  
  return(out)
}


## ------------------------
## End of gammaNUMCKpar
## ------------------------

