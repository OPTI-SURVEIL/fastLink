#' CKpar.b
#'
#' posterior match probability calculation based on string similarity and Fellegi-Sunter priors
#'
#' @usage CKpar.b(matAp, matBp, blocklist, dedupe.blocks, n.cores, transform, transform.args, 
#' method, method.args, isoreg, prior, resolution)
#'
#' @param matAp vector storing the comparison field in data set 1
#' @param matBp vector storing the comparison field in data set 2
#' @param blocklist list of blocking indices
#' @param dedupe.blocks Logical indicators for whether internal linkage and duplication is taking place
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#' @param transform An optional function returning transformed strings or a list of transformed strings. 
#' First argument should be the strings to be transformed
#' @param transform.args An optional list of arguments to the transform function
#' @param method String distance method, passed as a custom function with argument list stringdist.args. The first two arguments of 
#' any custom function provided should be the strings to be compared
#' @param method.args List of arguments for custom function provided to method
#' @param isoreg Isotonic regression object that will be scaled by prior weights to predict conditional probability of matching based on similarity scores
#' @param priors prior probability of match within the stratum examined for each block
#' @param resolution Number of probability bins over which to aggregate results
#' @param drop.lowest Whether to drop matched indices from the lowest probability bin - this can save a lot of memory and processing wasted on nonmatches
#' @return \code{CKpar.b} returns a list of lists (for each block) with the indices corresponding to posterior match probabilities within a number of probability bins given by the resolution parameter
#'
#' @export

CKpar.b <- function(matAp, matBp, blocklist, n.cores = NULL, transform = NULL, 
                          transform.args, method = NULL, method.args, dedupe.blocks, isoreg = NULL, priors = NULL,
                    resolution = 1e-2, drop.lowest = T){
  varnm = names(matAp)
  if(any(class(matAp) %in% c("tbl_df", "data.table",'data.frame'))){
    matAp <- as.data.frame(matAp)[,1]
  }
  if(any(class(matBp) %in% c("tbl_df", "data.table",'data.frame'))){
    matBp <- as.data.frame(matBp)[,1]
  }
  
  matAp[matAp == ""] <- NA
  matBp[matBp == ""] <- NA
  
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
  
  clusterExport(cl1, c('matrix.1','matrix.2','blocklist','dedupe.blocks','u.trans.1','u.trans.2','method','method.args'),
                envir = environment())
  
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
  
  clusterExport(cl1, exports,
                envir = environment())
  
  #find number of comparisons in each block
  cat('finding number of unique pairs of values of variable *', varnm, '* across all blocks \n', sep = '')
  
  compvals = foreach(i = 1:length(blocklist), .options.snow = opts1, .combine = c) %dopar% {
    b = blocklist[[i]]
    u.1 = unique(matrix.1[b[[1]]]); u.2 = unique(matrix.2[b[[2]]])
    
    as.numeric(length(u.1)) * as.numeric(length(u.2)) / 300^2
  }
  
  inpar = compvals <= 1
  
  #define comparison routine
  stringvec <- function(m, y, strdist = method, identical, binres = resolution, 
                        isoreg_ = isoreg, prior, drop.lowest = drop.lowest) {
    
    if(!is.list(m[[1]])){
      x <- as.matrix(m[[1]])
      e <- as.matrix(y[[1]])
      nr = nrow(e)
    }else{
      x = m[[1]]
      e = y[[1]]
      nr = length(e[[1]])
    }     
    
    newzeta_ = seq(0,1,binres)
    
    if(identical){
      t <- do.call(strdist,c(list(e), method.args))
    }else{
      t <- do.call(strdist,c(list(e, x), method.args))
    }
    
    isoreg_ = isoreg_scale(isoreg_, prior)
    
    t = fit.isored(isoreg_,t)
    
    bins = matrix(.bincode(t,newzeta_, include.lowest = T),ncol = ncol(t))
    
    newzeta = newzeta_[2:length(newzeta_)]
    counts = 0 * newzeta
    
    
    for(b in unique(bins)){
      counts[b] = sum(bins == b)
      newzeta[b] = mean(t[bins == b])
    }
    
    if(length(t) > 1e4) gc()       	
    slice.1 <- m[[2]]
    slice.2 <- y[[2]]
    res = list(meanprob = newzeta[1], count = counts[1], matches = matrix(numeric(0),ncol = 2))
    if(!drop.lowest) res$matches = which(bins == i, arr.ind = T) + 
      cbind(rep(slice.2,sum(bins==i)),rep(slice.1,sum(bins==i)))
    res = c(list(res), lapply(2:length(newzeta), function(i){
      list(meanprob = newzeta[i], count = counts[i], 
           matches = which(bins == i, arr.ind = T) + 
             cbind(rep(slice.2,sum(bins==i)),rep(slice.1,sum(bins==i))))
    }))
    
    names(res) = newzeta_[2:length(newzeta_)] - binres/2
    res
  }
  
  clusterExport(cl1, c('stringvec','inpar'),
                envir = environment())
  out = vector('list',length(blocklist))
  #perform comparisons on all small blocks
  
  if(sum(inpar)>0){
    cat('\nparallelizing over small blocks\n', sep = '')
  
    out[inpar] = foreach(i = 1:sum(inpar), 
                         .options.snow = opts1,
                   .packages = c('Matrix','fastLink','RcppAlgos',pkgs)) %dopar% {
      b = blocklist[inpar][[i]]
      ddp = dedupe.blocks[inpar][i]
      m.1 = matrix.1[b[[1]],]; m.2 = matrix.2[b[[2]],]
      
      u.values.1 <- unique(m.1)
      u.values.2 <- unique(m.2)
      
      u.values.1 <- u.values.1[u.values.1 != '1234MF']
      u.values.2 <- u.values.2[u.values.2 != '9876ES']
      
      ht.1 = lapply(u.values.1, function(x) which(m.1 == x))
      names(ht.1) = u.values.1
      ht.1 = list2env(ht.1)
      
      ht.2 = lapply(u.values.2, function(x) which(m.2 == x))
      names(ht.2) = u.values.2
      ht.2 = list2env(ht.2)
      
      if(!is.null(transform)){
        temp.1 = mget(u.values.1,u.trans.1)
        temp.1 = as.list(data.frame(do.call(rbind,temp.1),stringsAsFactors = F)); names(temp.1) = tr_names
        temp.2 = mget(u.values.2,u.trans.2)
        temp.2 = as.list(data.frame(do.call(rbind,temp.2),stringsAsFactors = F)); names(temp.2) = tr_names
      }else{
        temp.1 = u.values.1
        temp.2 = u.values.2
      }
                         
        temp.f <- stringvec(list(temp.1,0), list(temp.2,0), 
                            identical = ddp, prior = priors[[i]])
        #now turn u.value indices into lists of row indices
                   
        #first, aggregate by each unique value in column with least unique values, then substitute values and indices
        
        temp.f = lapply(temp.f, function(x){
          if(length(x$matches)==0) return(x)
          ms = lapply(1:2, function(i) sort(unique(x$matches[,i])))
          major = which.min(sapply(ms,length))
          ord = c(major,3-major)
          
          ms = ms[[major]]
          mms = tapply(x$matches[,3-major],x$matches[,major], function(x) x, simplify = F)
          mmids = as.integer(factor(sapply(mms,paste0,collapse = '.')))
          
          ms = lapply(unique(mmids), function(i) ms[mmids == i])
          mms = mms[!duplicated(mmids)]
          
          if(ddp){
            .identical = c(rep(F, length(ms)),rep(T, length(ms)))
            temp = c(lapply(seq_along(ms),function(i) list(ms[[i]], setdiff(mms[[i]],ms[[i]]))[ord]),
                     lapply(seq_along(ms),function(i) list(ms[[i]], intersect(mms[[i]],ms[[i]]))[ord]))
            .identical = .identical[sapply(temp, function(x) length(x[[2]]) > 0)]
            temp = temp[sapply(temp, function(x) length(x[[2]]) > 0)]
                     
          }else{
            temp = lapply(seq_along(ms),function(i) list(ms[[i]], mms[[i]])[ord])
            .identical = rep(F, length(temp))
          }
          temp = lapply(temp, function(x) list(u.values.2[x[[1]]],u.values.1[x[[2]]]))
          temp = lapply(temp, function(x){ res1 = b[[1]][unlist(mget(x[[2]],ht.1))]
                                           res2 = b[[2]][unlist(mget(x[[1]],ht.2))]
                                           #names(res1) <- NULL; names(res2) <- NULL
                                           list(res1,res2)})

          x$matches = temp
          x$.identical = .identical
          x})
        temp.f
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
      ht.1 = lapply(u.values.1, function(x) which(m.1 == x))
      names(ht.1) = u.values.1
      ht.1 = list2env(ht.1)
      
      ht.2 = lapply(u.values.2, function(x) which(m.2 == x))
      names(ht.2) = u.values.2
      ht.2 = list2env(ht.2)
      
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
        if(ddp) return(stringvec(temp[[r1]], temp[[r2]], identical = r1==r2 & ddp, prior = priors[[j]]))
        stringvec(temp.2[[r1]], temp.1[[r2]], identical = r1==r2 & ddp, prior = priors[[j]])
        
      }
      aggr.fun = function(...){
        meanprobs = sapply(..., '[[','meanprob')
        counts = sapply(..., '[[','count')
        matches = lapply(..., '[[','matches')
        list(meanprob = ifelse(sum(counts) ==0, meanprobs[1],sum(meanprobs * counts)/sum(counts)),
             count = sum(counts), matches = do.call(rbind,matches))
      }
      
      tmp = lapply(seq_along(temp.f[[1]]), function(i) lapply(temp.f,'[[', i))
      tmp = lapply(tmp,aggr.fun); names(tmp) = names(temp.f[[1]])
      
      close(pb2)
      
      temp.f = lapply(tmp, function(x){
        if(length(x$matches)==0) return(x)
        ms = lapply(1:2, function(i) sort(unique(x$matches[,i])))
        major = which.min(sapply(ms,length))
        ord = c(major,3-major)
        
        ms = ms[[major]]
        mms = tapply(x$matches[,3-major],x$matches[,major], function(x) x,simplify = F)
        mmids = as.integer(factor(sapply(mms,paste0,collapse = '.')))
        
        ms = lapply(unique(mmids), function(i) ms[mmids == i])
        mms = mms[!duplicated(mmids)]
        
        if(ddp){
          .identical = c(rep(F, length(ms)),rep(T, length(ms)))
          temp = c(lapply(seq_along(ms),function(i) list(ms[[i]], setdiff(mms[[i]],ms[[i]]))[ord]),
                   lapply(seq_along(ms),function(i) list(ms[[i]], intersect(mms[[i]],ms[[i]]))[ord]))
          .identical = .identical[sapply(temp, function(x) length(x[[2]]) > 0)]
          temp = temp[sapply(temp, function(x) length(x[[2]]) > 0)]
          
        }else{
          temp = lapply(seq_along(ms),function(i) list(ms[[i]], mms[[i]])[ord])
          .identical = rep(F, length(temp))
        }
        temp = lapply(temp, function(x) list(u.values.2[x[[1]]],u.values.1[x[[2]]]))
        temp = lapply(temp, function(x){ res1 = b[[1]][unlist(mget(x[[2]],ht.1))]
        res2 = b[[2]][unlist(mget(x[[1]],ht.2))]
        #names(res1) <- NULL; names(res2) <- NULL
        list(res1,res2)})
        
        x$matches = temp
        x$.identical = .identical
        x})
      
      out[!inpar][[j-sum(inpar)]] = temp.f
      capture.output(setTxtProgressBar(pb1,1))
      setTxtProgressBar(pb1,j)
      gc()
    }
  }
  #consolidate results from all blocks together
  aggr.fun = function(...){
    meanprobs = sapply(..., '[[','meanprob')
    counts = sapply(..., '[[','count')
    matches = lapply(..., '[[','matches')
    res = list(meanprob = ifelse(sum(counts) ==0, meanprobs[1],sum(meanprobs * counts)/sum(counts)),
         count = sum(counts), matches = do.call(c,matches))
  }
    
  tmp = lapply(seq_along(out[[1]]), function(i) lapply(out,'[[',i))
  tmp = lapply(tmp, aggr.fun)
  tmp = lapply(seq_along(tmp), function(i){
    names(tmp[[i]]) = paste0(names(tmp[[i]]),names(out[[1]])[i])
    tmp[[i]]
  })
  tmp = unlist(tmp,recursive = F)
  return(tmp)
}

