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
#' of match given the classifier score falls below this threshold are not evaluated. Defaults to 0.5. If set to 0, all blocks are evaluated.
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
                           min.posterior = 0.1, drop.lowest = T, keyvar = NULL, cluster = NULL, diceroll.gamma = F){
  
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
  
  getinds = (recov.pos >= 1 | min.posterior==0) & fastlinkres$patterns.w[,grep(varname,varnames)] == 0
  if(sum(getinds, na.rm=T) == 0){
    cat('No agreement patterns meet criteria for post-processing, returning NULL result \n')
    return(NULL)
  } 
  if(is.data.frame(fastlinkres$patterns.w)) fastlinkres$patterns.w = as.matrix(fastlinkres$patterns.w)
  matchpats = matrix(fastlinkres$patterns.w[which(getinds),1:nv], ncol = nv)
  colnames(matchpats) = colnames(fastlinkres$patterns.w)[1:nv]
  counts_old = fastlinkres$patterns.w[which(getinds),nv+1]
  
  matchkeys = do.call(paste,as.data.frame(matchpats))
  targkeys = do.call(paste,as.data.frame(fastlinkres$patterns.w[,1:nv]))
  
  interpats = matchpats; interpats[,grep(varname,varnames)] = 2
  interkeys = do.call(paste,data.frame(interpats))
  
  zeta_old = fastlinkres$zeta.j[which(getinds)]
  zeta_high = fastlinkres$zeta.j[match(interkeys,targkeys,nomatch = length(fastlinkres$zeta.j)+100)]
  levs = lapply(lapply(data.frame(fastlinkres$patterns.w[,1:nv]),factor),levels)
  for(ind in which(is.na(zeta_high))){
    newdat = data.frame(interpats)[ind,]
    newdat[is.na(newdat)] = -1
    newdat[] = lapply(1:ncol(newdat),function(i) factor(newdat[,i],levels = levs[[i]]))
    terms.m = names(coef(fastlinkres$match_mod))[-1]
    terms.m = gsub('2$|1$|0$','',terms.m); terms.m = gsub('2\\.|1\\.|0\\.',':',terms.m)
    newdat.m = model.matrix(as.formula(paste('~',paste(terms.m,collapse ='+'))), newdat)
    
    terms.u = names(coef(fastlinkres$nonmatch_mod))[-1]
    terms.u = gsub('2$|1$|0$','',terms.u); terms.u = gsub('2\\.|1\\.|0\\.',':',terms.u)
    newdat.u = model.matrix(as.formula(paste('~',paste(terms.u,collapse ='+'))), newdat)
    
    pos = predict(fastlinkres$match_mod,newdata = data.frame(newdat.m),type = 'response')
    neg = predict(fastlinkres$nonmatch_mod,newdata = data.frame(newdat.u),type = 'response')
    zeta_high[ind] = pos/(pos + neg)
  }
  
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
    
    pats = do.call(paste,lapply(varnames, function(v) (dfA[[v]][matches[[1]]] == dfB[[v]][matches[[2]]]) * 2))
    
    priors = fastlinkres$zeta.j[match(pats,targkeys)]
    uppers = zeta_high[match(pats,matchkeys)]
    tlist1 = lapply(utrans1,'[', match(dfA[[varname]][matches[[1]]], uvals1))
    tlist2 = lapply(utrans2,'[', match(dfB[[varname]][matches[[2]]], uvals2))
    
    if(!is.null(keyvar)){
      keys1 = dfA[[keyvar]][matches[[1]]]
      keys2 = dfB[[keyvar]][matches[[2]]]
      truths = keys1 == keys2
    }
    scores = do.call(method, c(list(tlist1, tlist2), method.args,nthread = n.thread))
    scores = fit.isored(isoreg,scores)
    #scores = uppers / ((1 / scores- 1) / ((1-isoreg$p.m)/isoreg$p.m) * ((1-priors)/priors) + 1) #if isoreg is based on name matching
    scores = 1 / ((1 / scores- 1) / ((1-isoreg$p.m)/isoreg$p.m) * ((1-priors)/priors) + 1) #if isoreg is based on record pairs
    
    bins = .bincode(scores, newzeta, include.lowest = T)
    if(diceroll.gamma){
      bin.incl = rbinom(length(scores),1,scores) == 1
      bincl.matches = list()
      if(sum(bin.incl > 0)){
        
        ms = lapply(matches[1:2], '[', bin.incl)
        
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
            temp2 = match(which(rowSums(temp1)/length(ums[[3-major]]) >= 0.5),ums[[major]])[1]
            if(length(temp2) == 0){swtch2 = F; break}
            xs = ums[[major]][temp2]; ys = Reduce(intersect,mms[temp2])
            dropinds = ms[[major]] %in% xs & ms[[3-major]] %in% ys
            if(sum(dropinds)==0 | sum(dropinds) == length(ms[[1]])){swtch2 = F; break}
            bmatches = c(bmatches, list(list(xs,ys)[ord]))
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
