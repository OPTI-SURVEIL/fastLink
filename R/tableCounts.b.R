#' tableCounts.b
#'
#' Count pairs with the same pattern in the cross product between two datasets. Version for blocked data.
#'
#' @usage tableCounts.b(gammalist, blocklist, dedupe.blocks, n.cores)
#'
#' @param gammalist A list (across blocks) of lists of objects produced by gammaKpar, gammaCK2par, or
#' gammaCKpar. 
#' @param blocklist A list of indices for the blocks
#' @param dedupe.blocks vector of logical indicators for whether a block represents self-comparison of a set of data
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield <benfifield@gmail.com>, and Kosuke Imai
#'
#' @return \code{tableCounts.b} returns counts of all unique mathching patterns, which can be
#' fed directly into \code{emlinkMAR} to get posterior matching probabilities for each unique pattern.
#'
#' @export
#' @importFrom parallel detectCores makeCluster stopCluster mclapply
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach "%dopar%" "%do%" foreach
## ------------------------
## To count unique patterns:
## tableCounts is the
## functions that does the trick
## ------------------------

tableCounts.b <- function(gammalist, blocklist, dedupe.blocks, n.cores = NULL) {
  
  nobs.a = sapply(blocklist, function(l) length(l[[1]]))
  nobs.b = sapply(blocklist, function(l) length(l[[2]]))
  if(is.null(names(gammalist))) names(gammalist) = 1:length(gammalist)
  
  if(is.null(n.cores)) {
    n.cores <- detectCores() - 1
  }
  
  #find slice dimensions for each block
  max_slice = 4500
  slice_a_adj = sapply(nobs.a, function(n) n/max_slice)
  slice_b_adj = sapply(nobs.b, function(n) n/max_slice)
  
  adj_ = cbind(slice_a_adj,slice_b_adj)
  adj_ = 1 / pmin(adj_[,1],adj_[,2],1)
  
  #slice each block
  n.slices1 <- pmax(round(nobs.a/(max_slice * adj_), 0), 1) 
  n.slices2 <- pmax(round(nobs.b/(max_slice * adj_), 0), 1) 
  
  ## Prep objects for m_func_par
  limit.1 <- lapply(1:length(blocklist), function(i) round(quantile((0:nobs.a[[i]]), p = seq(0, 1, 1/n.slices1[[i]])), 0))
  limit.2 <- lapply(1:length(blocklist), function(i) round(quantile((0:nobs.b[[i]]), p = seq(0, 1, 1/n.slices2[[i]])), 0))
  
  last1 <- sapply(limit.1, length)
  last2 <- sapply(limit.2, length)
  
  n.lim.1 <- lapply(1:length(blocklist), function(i) limit.1[[i]][-1] - limit.1[[i]][-last1[i]])
  n.lim.2 <- lapply(1:length(blocklist), function(i) limit.2[[i]][-1] - limit.2[[i]][-last2[i]])
  
  ind.i <- lapply(n.slices1, function(x) 1:x)
  ind.j <- lapply(n.slices2, function(x) 1:x)
  ind <- lapply(seq_along(blocklist), function(i) cbind(expand.grid.jc(ind.i[[i]], ind.j[[i]]),i))
  
  ind[dedupe.blocks] = lapply(seq_len(sum(dedupe.blocks)), function(i){
    keep_rows = which(ind[dedupe.blocks][[i]][,1]<=ind[dedupe.blocks][[i]][,2])
    matrix(ind[dedupe.blocks][[i]][keep_rows,],ncol = 3)
  })

  ez_rows = lapply(ind[dedupe.blocks], function(x) which(x[,1] == x[,2]))
  excess_zeroes = rep(0,length(blocklist))
  excess_zeroes[dedupe.blocks] = sapply(seq_len(sum(dedupe.blocks)), function(i){
    nlm = n.lim.1[dedupe.blocks][[i]]; ind. = ind[dedupe.blocks][[i]]; ezr = ez_rows[[i]]
    sum(nlm[ind.[ezr,1]] * (nlm[ind.[ezr,1]] + 1)/2)
  }) 
  
  ind = do.call(rbind,ind)
  
  ## Lists of indices:
  ##     temp - exact
  ##     ptemp - partial
  ##     natemp - NAs
  
  temp = lapply(seq_along(blocklist), function(i) lapply(seq_along(gammalist),function(j) gammalist[[j]][[i]]$matches2))
  ptemp = lapply(seq_along(blocklist), function(i) lapply(seq_along(gammalist),function(j) gammalist[[j]][[i]]$matches1))
  natemp = lapply(seq_along(blocklist), function(i) lapply(seq_along(gammalist),function(j) gammalist[[j]][[i]]$nas))
  identical = lapply(seq_along(blocklist), function(i) lapply(seq_along(gammalist),function(j) gammalist[[j]][[i]]$.identical))
  
  ## Take care of small blocks in parallel
  gammas = m_func_par_b(temp = temp, ptemp = ptemp, natemp = natemp,
                limit1 = limit.1, limit2 = limit.2,
                nlim1 = n.lim.1, nlim2 = n.lim.2,
                ind = ind, listid = rep(1, 2),
                identical = identical,
                dedupe = dedupe.blocks,
                matchesLink = FALSE, threads = n.cores)
  gc()
      
  gammas_mat = do.call(rbind,lapply(gammas,function(x){
    cbind(x[[1]],x[[2]])
  }))
  
  rm(gammas); gc()
  
  counts.f <- as.matrix(tapply(as.numeric(gammas_mat[, 2]), gammas_mat[, 1], sum))
  counts.d <- cbind( as.numeric(row.names(counts.f)), counts.f)
  
  if(any(dedupe.df)){
    counts.d[counts.d[,1]==0,2] = counts.d[counts.d[,1]==0,2] - sum(excess_zeroes) #remove excess zero counts
  }
  colnames(counts.d) <- c("pattern.id", "count")
  
  ## Merge Counts
  seq <- 1:(length(gammalist)*3)
  b <- 2^(seq)
  patterns.vec <- matrix(NA, 4, length(gammalist))
  for(i in 1:length(gammalist)){
    patterns.vec[,i] <- c(b[1:3 + (i-1)*3], 0)
  }
  patterns <- expand.grid(as.data.frame(patterns.vec))
  pattern.id <- rowSums(patterns)
  patterns <- cbind(patterns, pattern.id)
  data.new.0 <- merge(patterns, counts.d, by = "pattern.id")
  data.new.0 <- data.new.0[,-1]
  
  b<-2
  patterns.2vec <- c()
  for(i in 1:length(gammalist)){
    patterns.2vec <- c(patterns.2vec, 1/b^(1 + (i-1)*3))
  }
  patterns.2 <- t((patterns.2vec) * t(data.new.0[,1:length(gammalist)]))
  data.new.1 <- cbind(patterns.2, data.new.0[,length(gammalist)+1])
  names <- c(paste0("gamma.", names(gammalist)), "counts")
  colnames(data.new.1) <- names
  sub.nc <- which(colSums(data.new.1) == 0)
  sub.nc <- sub.nc[sub.nc > length(gammalist)]
  if(length(sub.nc) > 0){
    data.new.1 <- data.new.1[, -sub.nc]
  }
  nc <- ncol(data.new.1)
  na.data.new <- data.new.1[, -c(nc), drop = FALSE]
  na.data.new[na.data.new == 4] <- NA
  data.new <- cbind(na.data.new, data.new.1[, nc])
  colnames(data.new)[nc] <- "counts"
  data.new <- data.new[data.new[, nc] > 0, ]
  class(data.new) <- c("fastLink", "tableCounts")
  
  return(data.new)
  
}
