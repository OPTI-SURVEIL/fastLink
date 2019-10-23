#' post_proc_gamma
#'
#' posterior match probability calculation based on string similarity and Fellegi-Sunter priors
#'
#' @usage post_proc_gamma(dfA, dfB, varname, fastlinkres, gammalist, isoreg, method, method.args, 
#' transform, transform.args, n.cores, chunksize, resolution, min.posterior, drop.lowest)
#'
#' @param dfA First dataframe to be linked
#' @param dfB Second dataframe to be linked
#' @param varname Name of variable that posterior processing is based on. Defaults to 'Name'. Must be present in both dfA and dfB
#' @param fastlinkres Prior linkage result using exact matching on varname field, which gives fine-grained prior probabilities that any record pair is a match 
#' @param gammalist Gammalist object containing sparser representation of field agreement from the prior linkage result
#' @param isoreg A compact format isotonic regression object, used to convert classifier scores to match probabilities under prior class imbalance
#' @param method String distance method, passed as a custom function with argument list stringdist.args. The first two arguments of 
#' any custom function provided should be the strings to be compared
#' @param method.args List of arguments for custom function provided to method
#' @param transform An optional function returning transformed strings or a list of string transformations. 
#' First argument should be the strings to be transformed
#' @param transform.args An optional list of arguments to the transform function
#' @param n.cores Number of cores to be used for parallel operations. Defaults to the number of available cores minus 1
#' @param chunksize Target number of record pairs that will be recovered and compared at a time. 
#' This parameter defaults to 90000, but may need to be adjusted for best performance, 
#' depending on the number of record pairs evaluated and the complexity of the fuzzy matching method
#' @param resolution 1/Number of probability bins over which to aggregate results. Defaults to 1e-2 (100 bins)
#' @param unblocked.size Number of comparisons considered before resorting to parallel operations - can be tweaked to speed up performance in some scenarios
#' @param min.posterior Sets a cap to decide which agreement patterns are evaluated for fuzzy matching. Blocks for which the maximum attainable probability 
#' of match given the classifier score falls below this threshold are not evaluated. Defaults to 0.5.
#' @param drop.lowest Whether to drop matched indices from the lowest probability bin - this can save a lot of memory and processing wasted on nonmatches
#' @param keyvar optional variable name used to distinguish true match identity
#' @param cluster optionally provide a pre-made cluster object. can help with repeated cluster startup times.
#' @param diceroll.gamma if set to TRUE (the default), a gammalist-like object is constructed (binmatches) by sampling inclusion of each record pair 
#' according to the posterior probability of name matching, which provides an easy way to re-integrate the results into the fastLink protocol. If set to FALSE, 
#' binmatches is returned for ranges of posterior match scores binned according to the \code{resolution} parameter
#' @return \code{post_proc_gamma} returns a list with the binned posterior counts, probabilities, and matching indices at each level but the lowest
#'
#' @export

post_proc_gamma = function(dfA,dfB,varname = 'Name', fastlinkres, gammalist, isoreg, ecdf, method, method.args,
                           transform = NULL, transform.args = NULL, n.cores = NULL, chunksize = 300^2, resolution = 1e-2, unblocked.size = 2e6,
                           min.posterior = 0.1, drop.lowest = T, keyvar = NULL, cluster = NULL, diceroll.gamma = T){
  
  if(!is.null(cluster)){
    parallel::clusterEvalQ(cluster, rm(list = ls()))
    parallel::clusterEvalQ(cluster, gc())
  }
    
    
  dedupe = identical(dfA,dfB)
  if(is.null(n.cores)) n.cores = parallel::detectCores() - 1
  
  glist = gammalist[intersect(names(gammalist), fastlinkres$varnames)]
  nv = length(fastlinkres$varnames)
  varnames = fastlinkres$varnames
  
  predseq = ecdf$scores
  isopreds = fit.isored(isoreg,predseq)
  pos = fastlinkres$zeta.j * fastlinkres$patterns.w[,'counts']
  p.score.m = ecdf$m
  recov.pos = pos * 0
  
  for(i in 1:length(pos)){
    z = fastlinkres$zeta.j[i]
    ips = 1 / ((1 / isopreds - 1) / ((1-isoreg$p.m)/isoreg$p.m) * ((1-z)/z) + 1)
    if(!any(ips >= min.posterior)) next
    ind = which(ips >= min.posterior)[sum(ips >= min.posterior)]
    p.s.m = p.score.m[ind]
    recov.pos[i] = pos[i] * p.s.m
  }
  
  #min.prior = 1/((1 / min.posterior - 1)/((1 / max(isoreg$y) - 1) / ((1-isoreg$p.m)/isoreg$p.m)) + 1)
  
  getinds = recov.pos >= 1 & fastlinkres$patterns.w[,grep(varname,varnames)] == 0
  if(sum(getinds, na.rm=T) == 0){
    cat('No agreement patterns meet criteria for post-processing, returning NULL result \n')
    return(NULL)
  } 
  if(is.data.frame(fastlinkres$patterns.w)) fastlinkres$patterns.w = as.matrix(fastlinkres$patterns.w)
  matchpats = matrix(fastlinkres$patterns.w[which(getinds),1:nv], ncol = nv)
  zeta_old = fastlinkres$zeta.j[which(getinds)]
  counts_old = fastlinkres$patterns.w[which(getinds),nv+1]
  
  matchkeys = do.call(paste0,as.data.frame(matchpats))
  targkeys = do.call(paste0,as.data.frame(fastlinkres$patterns.w[,1:nv]))
  
  dummyEM = fastlinkres
  dummyEM$zeta.j[!(targkeys %in% matchkeys)] = -1
  dummyEM$patterns.w[!(targkeys %in% matchkeys),'weights'] = -1e24
  
  chunksize = floor(sqrt(chunksize))
  n.chunks = round(sum(fastlinkres$patterns.w[,'counts'])/sum(counts_old))
  
  n.slices1 <- max(round(as.numeric(fastlinkres$nobs.a)/(chunksize), 0), 1) 
  n.slices2 <- max(round(as.numeric(fastlinkres$nobs.b)/(chunksize), 0), 1) 
  nc <- min(n.cores, n.slices1 * n.slices2)
  
  limit.1 <- round(quantile((0:fastlinkres$nobs.a), p = seq(0, 1, 1/n.slices1)), 0)
  limit.2 <- round(quantile((0:fastlinkres$nobs.b), p = seq(0, 1, 1/n.slices2)), 0)
  
  last1 <- length(limit.1)
  last2 <- length(limit.2)
  
  n.lim.1 <- limit.1[-1] - limit.1[-last1]
  n.lim.2 <- limit.2[-1] - limit.2[-last2]
  
  power <- rep(NA, length(glist))
  for(i in 1:length(glist)){
    power[i] <- 1 + (i-1)*3
  }
  power.s <- power[1:nv]
  base <- 2^(power.s)
  list = matchpats; list[is.na(list)] <- 4
  list <- t(base * t(list))
  list.id <- rowSums(list)
  
  temp <- vector(mode = "list", length = length(glist))
  ptemp <- vector(mode = "list", length = length(glist))
  natemp <- vector(mode = "list", length = length(glist))
  identical <- vector(mode = "list", length = length(glist))
  
  for(i in 1:length(glist)){
    temp[[i]] <- glist[[i]]$matches2
    if(!is.null(glist[[i]]$matches1)) {
      ptemp[[i]] <- glist[[i]]$matches1
    }
    natemp[[i]] <- glist[[i]]$nas
    identical[[i]] <- glist[[i]]$.identical
  }
  
  ind.i <- 1:n.slices1
  ind.j <- 1:n.slices2
  ind <- as.matrix(expand.grid(ind.i, ind.j))
  
  if(dedupe){
    keep_rows = which(ind[,1]<=ind[,2])
    ind = matrix(ind[keep_rows,],ncol = 2) #removing empty blocks for internal linkage represented by upper right tri matrix
  }
  
  chunkseq = lapply(seq(1,nrow(ind),n.chunks), function(i) seq(i,min(i+n.chunks-1,nrow(ind)),1))
  
  
  
  if(sum(counts_old)>unblocked.size){
    if(is.null(cluster)){cl = makeCluster(min(n.cores, length(chunkseq)))} else{ cl = cluster} 
    registerDoSNOW(cl)
    pb = txtProgressBar(0,length(chunkseq))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    '%oper%' <- foreach::'%dopar%'
    depsm = find_dependencies(method)
    exports = unique(c(depsm$depends$calls[sapply(depsm$depends$pkgs,function(x) '.GlobalEnv' %in% x)]))
    pkgs = unique(c(unlist(depsm$depends$pkgs)))#,unlist(depst$depends$pkgs)))
    pkgs = pkgs[!(pkgs %in% c('base','.GlobalEnv'))]
    n.thread = 1

  }else{
    '%oper%' <- foreach::'%do%'
    chunkseq = list(1:nrow(ind))
    pkgs = NULL; exports = NULL
    opts = list()
    n.thread = n.cores
  }
  if(class(method.args$model) == 'xgb.Booster'){
    method.args$model = xgb.Booster.complete(method.args$model)
    if(sum(counts_old)>unblocked.size) {xgb.parameters(method.args$model) <- list(nthread = 1)}else{
      xgb.parameters(method.args$model) <- list(nthread = n.cores)
    }
  }
  
  uvals1 = unique(dfA[[varname]])
  uvals2 = unique(dfB[[varname]])
  
  if(!is.null(transform)){
    utrans1 = do.call(transform,c(list(uvals1), transform.args))
    utrans2 = do.call(transform,c(list(uvals2), transform.args))
  }else{
    utrans1 = uvals1
    utrans2 = uvals2
  }
  
  # clusterExport(cl, c('ind','uvals1','uvals2','utrans1',
  #                     'utrans2','dfA','dfB','temp','ptemp','natemp','limit.1','limit.2','n.lim.1','n.lim.2','list.id',
  #                     'identical','dedupe','targkeys','fastlinkres','keyvar','varname','method','method.args',
  #                     'isoreg',exports),envir = environment())

  res = foreach(i = chunkseq, .packages = unique(c('Matrix',pkgs)), .options.snow = opts, .export = exports) %oper% {
    gc()
    newzeta = seq(0,1,resolution)
    posteriors = newzeta[1:(1/resolution)] + resolution/2
    binmatches = lapply(posteriors,function(x)list())
    names(binmatches) = paste0('matches', posteriors)
    true.counts = counts = 0 * posteriors
    
    matches = m_func_par(temp = temp, ptemp = ptemp, natemp = natemp,
                         limit1 = limit.1, limit2 = limit.2,
                         nlim1 = n.lim.1, nlim2 = n.lim.2,
                         ind = matrix(ind[i, ],ncol=2), listid = list.id,
                         identical = identical,
                         dedupe = dedupe,
                         matchesLink = TRUE, threads = n.thread)
    
    gammas_mat <- lapply(matches, function(x){
      as.matrix(data.frame(x[[1]], x[[2]]))
    })
    
    matches <- as.data.frame(do.call('rbind', gammas_mat) + 1)
    
    pats = do.call(paste0,lapply(varnames, function(v) (dfA[[v]][matches[[1]]] == dfB[[v]][matches[[2]]]) * 2))
    
    priors = fastlinkres$zeta.j[match(pats,targkeys)]
    tlist1 = lapply(utrans1,'[', match(dfA[[varname]][matches[[1]]], uvals1))
    tlist2 = lapply(utrans2,'[', match(dfB[[varname]][matches[[2]]], uvals2))
    
    if(!is.null(keyvar)){
      keys1 = dfA[[keyvar]][matches[[1]]]
      keys2 = dfB[[keyvar]][matches[[2]]]
      truths = keys1 == keys2
    }
    scores = do.call(method, c(list(tlist1, tlist2), method.args,nthread = n.thread))
    scores = fit.isored(isoreg,scores)
    
    for(z in unique(priors)){
      inds = priors==z
      scores[inds] = 1 / ((1 / scores[inds]- 1) / ((1-isoreg$p.m)/isoreg$p.m) * ((1-z)/z) + 1)
    }
    bins = .bincode(scores, newzeta, include.lowest = T)
    if(diceroll.gamma){
      bin.incl = rbinom(length(scores),1,scores) == 1
      ms = lapply(matches[1:2], '[', bin.incl)
      bincl.matches = list()
      swtch2 = T
      while(swtch2){
        ums = lapply(ms, function(x) sort(unique(x)))
        major = which.min(sapply(ums,length))
        ord = c(major,3-major)
        
        mms = tapply(ms[[3-major]],ms[[major]],function(x) x)
        temp1 = sparseMatrix(i = ms[[major]], j = ms[[3-major]])
        temp2 = match(which(rowSums(temp1)/length(ums[[3-major]]) >= 0.5),ums[[major]])
        if(length(temp2) == 0){swtch2 = F; break}
        xs = ums[[major]][temp2]; ys = Reduce(intersect,mms[temp2])
        bincl.matches = c(bincl.matches, list(list(xs,ys)[ord]))
        dropinds = ms[[major]] %in% xs & ms[[3-major]] %in% ys
        if(sum(dropinds)==0 | sum(dropinds) == length(ms[[1]])){swtch2 = F; break}
        ms = lapply(ms,'[',!dropinds)
      }
      mmsl = sapply(mms,paste0,collapse = '.')
      mmsl = as.integer(factor(mmsl,levels = unique(mmsl)))
      
      for(l in unique(mmsl)){
        xs = ums[[major]][mmsl == l]
        ys = mms[[which(mmsl==l)[1]]]
        bincl.matches = c(bincl.matches, list(list(xs,ys)[ord]))
      }
      
    }
    for(b in sort(unique(bins))){
        ms = lapply(matches[1:2], '[', bins==b)
        if(!is.null(keyvar)) true.counts[b] = sum(truths[bins==b])
        posteriors[b] = mean(scores[bins==b])
        counts[b] = length(ms[[1]])
        
      if(!diceroll.gamma){
        if(b>1 | !drop.lowest){
          swtch2 = T
          bmatches = list()
          while(swtch2){
            ums = lapply(ms, function(x) sort(unique(x)))
            major = which.min(sapply(ums,length))
            ord = c(major,3-major)
            
            mms = tapply(ms[[3-major]],ms[[major]],function(x) x)
            temp1 = sparseMatrix(i = ms[[major]], j = ms[[3-major]])
            temp2 = match(which(rowSums(temp1)/length(ums[[3-major]]) >= 0.5),ums[[major]])
            if(length(temp2) == 0){swtch2 = F; break}
            xs = ums[[major]][temp2]; ys = Reduce(intersect,mms[temp2])
            bmatches = c(bmatches, list(list(xs,ys)[ord]))
            dropinds = ms[[major]] %in% xs & ms[[3-major]] %in% ys
            if(sum(dropinds)==0 | sum(dropinds) == length(ms[[1]])){swtch2 = F; break}
            ms = lapply(ms,'[',!dropinds)
          }
          
          mmsl = sapply(mms,paste0,collapse = '.')
          mmsl = as.integer(factor(mmsl,levels = unique(mmsl)))
          
          for(l in unique(mmsl)){
            xs = ums[[major]][mmsl == l]
            ys = mms[[which(mmsl==l)[1]]]
            bmatches = c(bmatches, list(list(xs,ys)[ord]))
          }
          binmatches[[b]] = bmatches
      }
      }
    }
    if(diceroll.gamma) binmatches = bincl.matches
    res = list(posteriors = posteriors, counts = counts, binmatches = binmatches)
    if(!is.null(keyvar)) res$true.counts = true.counts
    #ret = counts > 0
    #res = lapply(res,'[',ret)
    res
  }
  if(sum(counts_old)>unblocked.size) stopCluster(cl)
  
  counts = lapply(res,'[[','counts')
  posteriors = lapply(res,'[[','posteriors')
  binmatches = lapply(res,'[[','binmatches')
  if(!is.null(keyvar)) true.counts = lapply(res,'[[','true.counts')
  
  posteriors = rowSums(do.call(cbind,posteriors) * do.call(cbind,counts)) / rowSums(do.call(cbind,counts))
  counts = Reduce('+', counts)
  if(diceroll.gamma){
    binmatches = do.call(c,binmatches)
  }else{
    binmatches = lapply(seq_along(binmatches[[1]]), function(i) do.call(c,lapply(binmatches,'[[',i)))
    names(binmatches) = paste0('matches',seq(0,1-resolution,resolution) + resolution/2)
  }
  
  
  if(!is.null(keyvar)) true.pmatch = Reduce('+',true.counts) / counts  
  
  incl = !is.na(posteriors)
  res = list(posteriors = posteriors[incl], counts = counts[incl], binmatches = binmatches[incl], matchkeys = matchkeys)
  if(!is.null(keyvar)) res$true.pmatch = true.pmatch[incl]
  res
}

# old version using blocking for large indices
# post_proc_gamma = function(dfA,dfB,varname = 'Name', fastlinkres, gammalist, isoreg, method, method.args,
#                            transform = NULL, transform.args = NULL, n.cores = NULL, unblocked.size = 2e6, resolution = 1e-2,
#                            min.posterior = 0.5, drop.lowest = T){
#   
#   dedupe = identical(dfA,dfB)
#   if(is.null(n.cores)) n.cores = parallel::detectCores() - 1
#   
#   glist = gammalist[intersect(names(gammalist), fastlinkres$varnames)]
#   nv = length(fastlinkres$varnames)
#   varnames = fastlinkres$varnames
#   
#   min.prior = 1/((1 / min.posterior - 1)/((1 / max(isoreg$y) - 1) / ((1-isoreg$p.m)/isoreg$p.m)) + 1)
#   
#   getinds = fastlinkres$zeta.j >= min.prior & fastlinkres$patterns.w[,grep(varname,varnames)] == 0
#   if(sum(getinds, na.rm=T) == 0){
#     cat('No agreement patterns meet criteria for post-processing, returning original gammalist \n')
#     return(glist)
#   } 
#   
#   matchpats = matrix(fastlinkres$patterns.w[which(getinds),1:nv], ncol = nv)
#   counts_old = fastlinkres$patterns.w[which(getinds),nv+1]
#   zeta_old = fastlinkres$zeta.j[which(getinds)]
#   
#   ord = order(counts_old)
#   matchpats = matchpats[ord,]; counts_old = counts_old[ord]; zeta_old = zeta_old[ord]
#   
#   matchkeys = do.call(paste0,as.data.frame(matchpats))
#   targkeys = do.call(paste0,as.data.frame(fastlinkres$patterns.w[,1:nv]))
#   
#   np = which.min(abs(cumsum(counts_old)-unblocked.size))
#   smallpats = matrix(matchpats[1:np,], ncol = nv)
#   
#   dummyEM = fastlinkres
#   dummyEM$zeta.j[!(targkeys %in% matchkeys[1:np])] = -1
#   dummyEM$patterns.w[!(targkeys %in% matchkeys[1:np]),'weights'] = -1e24
#   
#   matches = do.call(cbind,matchesLink(gammalist,fastlinkres$nobs.a, fastlinkres$nobs.b,dummyEM,0,dedupe = dedupe))
#   
#   pats = do.call(paste0,lapply(varnames, function(v) (dfA[[v]][matches[,1]] == dfB[[v]][matches[,2]]) * 2))
#   priors = fastlinkres$zeta.j[match(pats,targkeys)]
#   
#   uvals1 = unique(dfA[[varname]][matches[,1]])
#   uvals2 = unique(dfB[[varname]][matches[,2]])
#   
#   if(!is.null(transform)){
#     utrans1 = do.call(transform,c(list(uvals1), transform.args))
#     utrans2 = do.call(transform,c(list(uvals2), transform.args))
#   }else{
#     utrans1 = uvals1
#     utrans2 = uvals2
#   }
#   
#   tlist1 = lapply(utrans1,'[', match(dfA[[varname]][matches[,1]], uvals1))
#   tlist2 = lapply(utrans2,'[', match(dfB[[varname]][matches[,2]], uvals2))
#   
#   #may need to come up with a better solution to this
#   method.args$combo = F
#   
#   scores = do.call(method, c(list(tlist1, tlist2), method.args,nthread = n.cores))
#   
#   for(z in unique(priors)){
#     scores[priors == z] = fit.isored(isoreg_scale(wbe_isofill_noem, z), scores[priors==z])
#   }
#   
#   newzeta = seq(0,1,resolution)
#   
#   bins = .bincode(scores, newzeta, include.lowest = T)
#   
#   for(b in sort(unique(bins))){
#     ms = matrix(matches[bins == b,],ncol = 2)
#     ums = sort(unique(ms[,1]))
#     breaks = tapply(ms[,2],ms[,1],function(x) x)
#     if(b==1 & drop.lowest){
#       glist[[varname]][[paste0('matches',newzeta[b]+resolution/2)]] = numeric(0)
#     }else{
#       glist[[varname]][[paste0('matches',newzeta[b]+resolution/2)]] = lapply(seq_along(ums),function(i)
#         list(ums[i],breaks[[i]]))
#     }
#     
#     glist[[varname]][[paste0('count',newzeta[b]+1e-2/2)]] = nrow(ms)
#     glist[[varname]][[paste0('meanprob',newzeta[b]+1e-2/2)]] = mean(scores[bins==b])
#     
#   }
#   
#   if(np < nrow(matchpats)){
#     #may need to come up with a better solution to this
#     method.args$combo = T
#     
#     largepats = matrix(matchpats[(np+1):nrow(matchpats),],ncol = nv)
#     blocklist = list()
#     priors = vector()
#     
#     nonmatch_vars = apply(largepats,2,function(x) 0 %in% x)
#     NA_vars = apply(largepats,2,function(x) any(is.na(x)))
#     NAblocks = lapply(glist[NA_vars], function(x) list(matches2 = list(list(dfA.inds = x$nas[[1]], dfB.inds = 1:nrow(dfB)),
#                                                                        list(dfA.inds = setdiff(1:nrow(dfA),x$nas[[1]]), 
#                                                                             dfB.inds = x$nas[[2]]))))
#     #need to add in handling of NAs as well
#     
#     unblocks = lapply(fastlinkres$varnames[nonmatch_vars],
#                       function(v) fastLink:::combineBlocks(gammalist[v],c(nrow(dfA),nrow(dfB))))
#     names(unblocks) = fastlinkres$varnames[nonmatch_vars]
#     unblocks_ = lapply(unblocks, function(x) lapply(1:ncol(x[[1]]), function(j) list(dfA.inds = which(x[[1]][,j]),
#                                                                                      dfB.inds = which(x[[2]][,j]))))
#     for(i in 1:nrow(largepats)){
#       p = largepats[i,]
#       # may be problems here if 0 or 2 are not in the pattern
#       
#       blocks = fastLink:::combineBlocks(c(glist[which(p==2)],NAblocks[varnames[is.na(p)]]),c(nrow(dfA),nrow(dfB)))
#       blocks = lapply(1:ncol(blocks[[1]]), function(j) list(dfA.inds = which(blocks[[1]][,j]),
#                                                             dfB.inds = which(blocks[[2]][,j])))
#       
#       blocks = fastLink:::thin_blocks(blocks, do.call(c,unblocks_[varnames[which(p==0)]]), dims = c(nrow(dfA),nrow(dfB)),
#                                       dedupe = dedupe)
#       
#       blocklist = c(blocklist, blocks)
#       priors = c(priors, rep(zeta_old[np + i], length(blocks)))
#       print(i/nrow(largepats))
#     }
#     if(dedupe){
#       dedupe.blocks = sapply(blocklist,function(x) identical (x[[1]],x[[2]]))
#     }else{
#       dedupe.blocks = rep(F, length(blocklist))
#     }
#     glist2 = fastLink::CKpar.b(dfA[varname],dfB[varname],blocklist,priors = priors,drop.lowest = drop.lowest,
#                                dedupe.blocks = dedupe.blocks,isoreg = isoreg,
#                                method = method, method.args = method.args,
#                                transform = transform, transform.args = transform.args)
#     
#     common = intersect(names(glist2), names(glist[varname]))
#     countvars = grep('count',common,val=T)
#     meanvars = grep('mean',common,val=T)
#     indvars = grep('match',common,val=T)
#     
#     glist2[meanvars] = as.list((unlist(glist2[meanvars]) * unlist(glist2[countvars]) + 
#                                   unlist(glist[varname][meanvars]) * unlist(glist[varname][countvars]))/ 
#                                  (unlist(glist[varname][countvars]) + unlist(glist2[countvars])))
#     
#     glist2[countvars] = as.list((unlist(glist[varname][countvars]) + unlist(glist2[countvars])))
#     
#     glist2[indvars] = lapply(indvars, function(v) c(glist2[[v]],glist[varname][[v]]))
#     glist[varname][common] = NULL
#     glist[varname] = c(glist[varname], glist2)
#   }
#   glist
# }
