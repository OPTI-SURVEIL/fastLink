#' matchesLink.b
#'
#' matchesLink.b, operating on blocked data, produces two dataframes that store
#' all the pairs that share a pattern that conforms
#' to the an interval of the Fellegi-Sunter
#' weights
#'
#' @usage matchesLink.b(gammalist, blocklist, em, thresh, dedupe.blocks, n.cores = NULL)
#'
#' @param gammalist A list of objects produced by either gammaKpar.b or
#' gammaCKpar.b. 
#' @param blocklist A list containing indices of blocks for record linkage.
#' @param em parameters obtained from the Expectation-Maximization algorithm under the MAR assumption. These estimates are
#' produced by emlinkMARmov
#' @param thresh is the interval of posterior zeta values for the agreements that we want to examine closer. Ranges between 0 and 1.
#' Can be a vector of length 1 (from specified value to 1) or 2 (from first specified value to second specified value).
#' @param dedupe.blocks vector of logical indicators for whether a block represents a self-comparison set of data
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#'
#' @return \code{matchesLink.b} returns an nmatches X 2 matrix with the indices of the
#' matches rows in dataset A and dataset B.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield <benfifield@gmail.com>, and Kosuke Imai
#'
#' @export

## ------------------------
## To recover the matches (their indices)
## we use matchesLink
## ------------------------

matchesLink.b <- function(gammalist, blocklist, em, thresh, dedupe.blocks, n.cores = NULL) {
  
  nobs.a = sapply(blocklist, function(l) length(l[[1]]))
  nobs.b = sapply(blocklist, function(l) length(l[[2]]))
  
  if(is.null(n.cores)) {
    n.cores <- detectCores() - 1
  }
  
  if(min(thresh) < 0 | max(thresh) > 1){
    stop("The specified threshold values are not valid posterior probabilities. These must range between 0 and 1.")
  }
  
  ## Slicing the data:
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
  
  ## Get the correct cuts
  em.obj <- data.frame(em$patterns.w)
  em.obj$zeta.j <- em$zeta.j
  em.obj <- em.obj[order(em.obj[, "weights"]), ]
  l.t <- thresh[1]
  u.t <- thresh[2]
  
  l.b <- suppressWarnings(min(em.obj$weights[em.obj$zeta.j >= l.t]))
  if(is.na(u.t)){
    u.b <- 1e10
  }else{
    u.b <- max(em.obj$weights[em.obj$zeta.j < u.t])
  }
  
  tablem <- em$patterns.w[em$patterns.w[, "weights"] >= l.b & em$patterns.w[, "weights"] < u.b, ]
  list <- tablem
  list[is.na(list)] <- 4
  
  if(is.null(dim(list))) {
    list <- t(as.matrix(list))
  }
  
  list <- list[, !colnames(list) %in% c("counts", "weights", "p.gamma.j.m", "p.gamma.j.u")]
  
  if(is.null(dim(list))) {
    list <- t(as.matrix(list))
  }
  
  ncol <- ncol(list)
  power <- rep(NA, length(gammalist))
  for(i in 1:length(gammalist)){
    power[i] <- 1 + (i-1)*3
  }
  power.s <- power[1:ncol]
  base <- 2^(power.s)
  list <- t(base * t(list))
  list.id <- rowSums(list)
  
  ## Lists of indices:
  ##     temp - exact
  ##     ptemp - partial
  ##     natemp - NAs
  temp = lapply(seq_along(blocklist), function(i) lapply(seq_along(gammalist),function(j) gammalist[[j]][[i]]$matches2))
  ptemp = lapply(seq_along(blocklist), function(i) lapply(seq_along(gammalist),function(j) gammalist[[j]][[i]]$matches1))
  natemp = lapply(seq_along(blocklist), function(i) lapply(seq_along(gammalist),function(j) gammalist[[j]][[i]]$nas))
  identical = lapply(seq_along(blocklist), function(i) lapply(seq_along(gammalist),function(j) gammalist[[j]][[i]]$.identical))
 
  ind.i <- lapply(n.slices1, function(x) 1:x)
  ind.j <- lapply(n.slices2, function(x) 1:x)
  ind <- lapply(seq_along(blocklist), function(i) cbind(expand.grid.jc(ind.i[[i]], ind.j[[i]]),i))
  
  ind[dedupe.blocks] = lapply(seq_len(sum(dedupe.blocks)), function(i){
    keep_rows = which(ind[dedupe.blocks][[i]][,1]<=ind[dedupe.blocks][[i]][,2])
    matrix(ind[dedupe.blocks][[i]][keep_rows,],ncol = 3)
  })
  
  ind = do.call(rbind,ind)
  
  ## Run main function
  gammas = m_func_par_b(temp = temp, ptemp = ptemp, natemp = natemp,
                                   limit1 = limit.1, limit2 = limit.2,
                                   nlim1 = n.lim.1, nlim2 = n.lim.2,
                                   ind = ind, listid = list.id,
                                   identical = identical,
                                   dedupe = dedupe.blocks,
                                   matchesLink = T, threads = n.cores)
    
  gammas_mat <- lapply(1:length(gammas), function(i){
    x = gammas[[i]]
    blk = ind[i,3]
    res = cbind(blocklist[[blk]]$dfA.inds[x[[1]] + 1],
                blocklist[[blk]]$dfB.inds[x[[2]] + 1])
  })
  
  temp <- do.call('rbind', gammas_mat)
  
  rm(gammas, gammas_mat); gc()
  
  temp <- data.frame(inds.a = temp[,1], inds.b = temp[,2])
  
  class(temp) <- c("fastLink", "matchesLink")
  
  return(temp)
}

## ------------------------
## End of matcheLink
## ------------------------

