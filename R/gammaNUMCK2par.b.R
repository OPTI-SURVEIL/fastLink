#' gammaNUMCK2par.b
#'
#' Field comparisons for numeric variables for blocked data. Three possible agreement patterns are considered:
#' 0 total disagreement, 2 agreement.
#' The distance between numbers is calculated using their absolute distance.
#'
#' @usage gammaNUMCK2par.b(matAp, matBp, blocklist, dedupe.blocks, n.cores, cut.a)
#'
#' @param matAp vector storing the comparison field in data set 1
#' @param matBp vector storing the comparison field in data set 2
#' @param blocklist list of blocking indices
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#' @param cut.a Lower bound for full match. Default is 1
#' @param dedupe.blocks Logical indicators for whether internal linkage and duplication is taking place
#' @return \code{gammaNUMCK2par.b} returns a list of lists (for each block) with the indices corresponding to each
#' matching pattern, which can be fed directly into \code{tableCounts.b} and \code{matchesLink.b}.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield <benfifield@gmail.com>, and Kosuke Imai
#' 
#' @export
#' @importFrom RcppAlgos comboGeneral

gammaNUMCK2par.b <- function(matAp, matBp, blocklist, n.cores = NULL, cut.a = 1, dedupe.blocks) {
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
  
  if(is.null(n.cores)) {
    n.cores <- detectCores() - 1
  }
  
  matrix.1 <- as.matrix(as.numeric(matAp))
  matrix.2 <- as.matrix(as.numeric(matBp))
  
  max <- max(max(matrix.1, na.rm = T), max(matrix.2, na.rm = T))   
  end.points <- c((round((max), 0) + 1), (round(max + cut.p, 0) + 3))
  matrix.1[is.na(matrix.1)] <- end.points[2]
  matrix.2[is.na(matrix.2)] <- end.points[1]
  
  #progress bar for blocklist level operations
  pb1 = txtProgressBar(0, length(blocklist),style = 1,char = '+')
  progress1 = function(n) setTxtProgressBar(pb1, n)
  opts1 = list(progress = progress1)
  
  cl1 = makeCluster(n.cores)
  registerDoSNOW(cl1)
  clusterExport(cl1, c('matrix.1','matrix.2','blocklist','dedupe.blocks','cut.a'))
  
  #find number of comparisons in each block
  cat('finding number of unique pairs of values of variable *', varnm, '* across all blocks \n', sep = '')
  
  compvals = foreach(i = 1:length(blocklist), .options.snow = opts1, .combine = c) %dopar% {
    b = blocklist[[i]]
    u.1 = unique(matrix.1[b[[1]]]); u.2 = unique(matrix.2[b[[2]]])
    
    as.numeric(length(u.1)) * as.numeric(length(u.2)) / 4500^2
  }
  
  inpar = compvals <= 1
  
  #perform comparisons on all small blocks
  
  difference <- function(m, y, cut,identical) {
    x <- as.matrix(m[[1]])
    e <- as.matrix(y[[1]])        
    
    if(identical){
      inds = comboGeneral(1:nrow(x),2)
      res <- abs(x[inds[,1]] - x[inds[,2]])
      t <- matrix(cut + 1,ncol = nrow(x),nrow = nrow(x))
      t[lower.tri(t)] = res
      diag(t) = 0
      diag(t)[x == end.points[2]] <- cut + 1 #keep out NAs
    }else{
      t <- fastLink:::calcPWDcpp(as.matrix(x), as.matrix(e))
    }
    t[ t == 0 ] <- cut
    t[ t > cut ] <- 0
    t <- Matrix(t, sparse = T)
    
    t@x[t@x <= cut] <- 2; gc()       	                
    
    slice.1 <- m[[2]]
    slice.2 <- y[[2]]
    indexes.2 <- which(t == 2, arr.ind = T)
    indexes.2[, 1] <- indexes.2[, 1] + slice.2
    indexes.2[, 2] <- indexes.2[, 2] + slice.1
    list(indexes.2)
  }
  
  clusterExport(cl1, c('difference','inpar'))
  out = vector('list',length(blocklist))
  
  if(sum(inpar)>0){
    cat('\nparallelizing over small blocks\n', sep = '')
  
    out[inpar] = foreach(i = 1:sum(inpar), 
                         .options.snow = opts1,.packages = c('Matrix','fastLink','RcppAlgos')) %dopar% {
      b = blocklist[inpar][[i]]
      ddp = dedupe.blocks[inpar][i]
      m.1 = matrix.1[b[[1]],]; m.2 = matrix.2[b[[2]],]
      
      u.values.1 <- unique(m.1)
      u.values.2 <- unique(m.2)
      
      if(ddp & length(u.values.1) == 1){
        final.list2 <- list(1:length(m.1), 1:length(m.2))
        final.list1 <- list()
        .identical = T
      }else{
        temp.f <- difference(list(u.values.1,0), list(u.values.2,0), cut.a, identical = ddp)
        
        indexes.2 <- temp.f[[1]]
        .identical = indexes.2[,1] == indexes.2[,2]
        
        n.values.2 <- as.matrix(cbind(u.values.1[indexes.2[, 2]], u.values.2[indexes.2[, 1]]))
        
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
      na.list[[1]] <- which(m.1 == end.points[2])
      na.list[[2]] <- which(m.2 == end.points[1])  
      
      out = list(matches2 = final.list2, nas = na.list, .identical = .identical)
      class(out) <- c("fastLink", "gammaNUMCK2par")
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
      
      if(ddp){
        n.slices <- max(round(length(u.values.1)/(4500), 0), 1) 
        limit = round(seq(0,length(u.values.1),length(u.values.1)/n.slices),0)
        
        temp = list()
        
        for(i in 1:n.slices) {
          temp[[i]] <- list(u.values.1[(limit[i]+1):limit[i+1]], limit[i])
        }
        
        do <- expand.grid.jc(1:n.slices,1:n.slices)
        drop = do[,2] < do[,1]
        do = matrix(do[!drop,],ncol=2)
      }else{
        n.slices1 <- max(round(length(u.values.1)/(4500), 0), 1) 
        n.slices2 <- max(round(length(u.values.2)/(4500), 0), 1) 
        
        limit.1 <- round(quantile((0:length(u.values.2)), p = seq(0, 1, 1/n.slices2)), 0)
        limit.2 <- round(quantile((0:length(u.values.1)), p = seq(0, 1, 1/n.slices1)), 0)
        
        temp.1 <- temp.2 <- list()
        
        for(i in 1:n.slices2) {
          temp.1[[i]] <- list(u.values.2[(limit.1[i]+1):limit.1[i+1]], limit.1[i])
        }
        
        for(i in 1:n.slices1) {
          temp.2[[i]] <- list(u.values.1[(limit.2[i]+1):limit.2[i+1]], limit.2[i])
        }
        
        do <- expand.grid(1:n.slices2, 1:n.slices1)
      }
      
      temp.f <- foreach(i = 1:nrow(do), .packages = c("Matrix",'fastLink','RcppAlgos')) %dopar% {
        r1 <- do[i, 1]
        r2 <- do[i, 2]
        
        difference(temp[[r1]], temp[[r2]], cut.a, identical = r1==r2 & ddp)
      }
      
      indexes.2 <- do.call('rbind', lapply(temp.f,'[[',1))
      .identical = indexes.2[,1] == indexes.2[,2]
      
      n.values.2 <- as.matrix(cbind(u.values.1[indexes.2[, 2]], u.values.2[indexes.2[, 1]]))
      
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
      na.list[[1]] <- which(matrix.1 == end.points[2])
      na.list[[2]] <- which(matrix.2 == end.points[1])
      out[!inpar][[j-sum(inpar)]] = list()
      out[!inpar][[j-sum(inpar)]][["matches2"]] <- final.list2
      out[!inpar][[j-sum(inpar)]][["nas"]] <- na.list
      out[!inpar][[j-sum(inpar)]][[".identical"]] <- .identical
      class(out[!inpar][[j-sum(inpar)]]) <- c("fastLink", "gammaNUMCK2par")
      
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

