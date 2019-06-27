#' stringSubset
#'
#' Removes as candidate matches any observations with no close matches on
#' string-distance measures.
#'
#' @usage stringSubset(vecA, vecB, similarity.threshold, stringdist.method,
#' jw.weight, n.cores)
#' @param vecA A character or factor vector from dataset A
#' @param vecB A character or factor vector from dataset B
#' @param similarity.threshold Lower bound on string-distance measure for being considered a possible match.
#' If an observation has no possible matches above this threshold, it is discarded from the match. Default is 0.8.
#' @param stringdist.method The method to use for calculating string-distance similarity. Possible values are
#' 'jaro' (Jaro Distance), 'jw' (Jaro-Winkler), and 'lv' (Levenshtein). Default is 'jw'.
#' @param jw.weight Parameter that describes the importance of the first characters of a string (only needed if stringdist.method = "jw"). Default is .10.
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#'
#' @return A list of length two, where the both entries are a vector of indices to be included in the match from dataset A (entry 1) and dataset B (entry 2). 
#'
#' @examples
#' \dontrun{
#' subset_out <- stringSubset(dfA$firstname, dfB$lastname, n.cores = 1)
#' fl_out <- fastLink(dfA[subset_out$dfA.block == 1,], dfB[subset_out$dfB.block == 1,],
#' varnames = c("firstname", "lastname", "streetname", "birthyear"), n.cores = 1)
#' }
#' @export
stringSubset <- function(vecA, vecB,
                         similarity.threshold = .8, stringdist.method = "jw",
                         jw.weight = .10, n.cores = NULL){

    if(class(vecA) == "factor"){
        vecA <- as.character(vecA)
    }
    if(class(vecB) == "factor"){
        vecB <- as.character(vecB)
    }
    if(class(vecA) != "character" | class(vecB) != "character"){
        stop("vecA and vecB must be of class factor or character.")
    }
    if(!(stringdist.method %in% c("jw", "jaro", "lv"))){
        stop("Invalid string distance method. Method should be one of 'jw', 'jaro', or 'lv'.")
    }
    if(similarity.threshold < 0 | similarity.threshold > 1){
        stop("similarity.threshold must be between 0 and 1.")
    }
    if(stringdist.method == "jw" & !is.null(jw.weight)){
        if(jw.weight < 0 | jw.weight > 0.25){
            stop("Invalid value provided for jw.weight. Remember, jw.weight in [0, 0.25].")
        }
    }
    
    ## Remove any very unlikely matches by first name
    gamma.out <- gammaCK2par(vecA, vecB, cut.a = similarity.threshold, method = stringdist.method, n.cores = n.cores)
    gamma.sub <- do.call(Map, c(c, gamma.out[[1]]))

    ## Get the voter file ids
    ids.A <- unique(gamma.sub[[1]])
    ids.B <- unique(gamma.sub[[2]])
    out <- list(dfA.inds = ids.A[order(ids.A)], dfB.inds = ids.B[order(ids.B)])
    class(out) <- "fastLink.block"
    
    return(out)
    
}

#' blockData
#'
#' Contains functionalities for blocking two data sets on one or more variables prior to
#' conducting a merge.
#'
#' @usage blockData(dfA, dfB, varnames, window.block, window.size,
#' kmeans.block, nclusters, iter.max, n.cores)
#' @param dfA Dataset A - to be matched to Dataset B
#' @param dfB Dataset B - to be matched to Dataset A
#' @param varnames A vector of variable names to use for blocking.
#' Must be present in both dfA and dfB. Blocks for each variable will be combined conjunctively, i.e. EACH blocking condition has to be true
#' @param window.block A vector of variable names indicating that the variable should be
#' blocked using windowing blocking. Must be present in varnames. If dfA and dfB are identical, the window blocks are deduplicated
#' @param window.size The size of the window for window blocking. Default is 1
#' (observations +/- 1 on the specified variable will be blocked together).
#' @param kmeans.block A vector of variable names indicating that the variable should be
#' blocked using k-means blocking. Must be present in varnames.
#' @param nclusters Number of clusters to create with k-means. Default value is the
#' number of clusters where the average cluster size is 100,000 observations.
#' @param iter.max Maximum number of iterations for the k-means algorithm to run. Default is 5000
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#'
#' @return A list with an entry for each block. Each list entry contains two vectors --- one with the indices indicating the block members in dataset A,
#' and another containing the indices indicating the block members in dataset B.
#'
#' @usage
#' \dontrun{
#' block_out <- blockData(dfA, dfB, varnames = c("city", "birthyear"))
#' }
#'
#' @export
blockData <- function(dfA, dfB, varnames, window.block = NULL,
                      window.size = 1,
                      kmeans.block = NULL,
                      nclusters = max(round(min(nrow(dfA), nrow(dfB)) / 100000, 0), 1),                      
                      iter.max = 5000, 
                      n.cores = NULL){

    cat("\n")
    cat(c(paste(rep("=", 20), sep = "", collapse = ""), "\n"))
    cat("blockData(): Blocking Methods for Record Linkage\n")
    cat(c(paste(rep("=", 20), sep = "", collapse = ""), "\n\n"))
    
    ## ---------------------------
    ## Clean data and check inputs
    ## ---------------------------
    if(any(class(dfA) %in% c("tbl_df", "data.table"))){
        dfA <- as.data.frame(dfA)
    }
    if(any(class(dfB) %in% c("tbl_df", "data.table"))){
        dfB <- as.data.frame(dfB)
    }
    
    dedupe = identical(dfA,dfB)
    
    if(any(!(varnames %in% names(dfA)))){
        stop("Some variables in varnames are not present in dfA.")
    }
    if(any(!(varnames %in% names(dfB)))){
        stop("Some variables in varnames are not present in dfB.")
    }
    if(any(!(window.block %in% varnames))){
        stop("You have provided a variable name for window.block that is not in 'varnames'.")
    }
    if(any(!(kmeans.block %in% varnames))){
        stop("You have provided a variable name for kmeans.block that is not in 'varnames'.")
    }
    classA <- lapply(dfA[,varnames], class)
    classB <- lapply(dfB[,varnames], class)
    if(any(unlist(classA)[names(classA) %in% window.block] != "numeric") |
       any(unlist(classB)[names(classB) %in% window.block] != "numeric")){
        stop("You have specified that a variable be blocked using window blocking, but that variable is not of class 'numeric'. Please check your variable classes.")
    }
    if(is.null(n.cores)){
        n.cores <- parallel::detectCores() - 1
    }
    
    
    ## ----------
    ## Block data
    ## ----------
    cat("Blocking variables.\n")
    gammalist <- vector(mode = "list", length = length(varnames))
    for(i in 1:length(varnames)){
        blocktype <- ifelse(varnames[i] %in% window.block, "window", ifelse(varnames[i] %in% kmeans.block, "k-means", "exact"))
        cat("    Blocking variable", varnames[i], "using", blocktype, "blocking.\n")
        if(varnames[i] %in% window.block){
            gammalist[[i]] <- gammaNUMCK2par(dfA[,varnames[i]], dfB[,varnames[i]], n.cores = n.cores, cut.a = window.size,dedupe = dedupe)
        }else if(varnames[i] %in% kmeans.block){
            gammalist[[i]] <- kmeansBlock(dfA[,varnames[i]], dfB[,varnames[i]],
                                          nclusters = nclusters,
                                          iter.max = iter.max,
                                          n.cores = n.cores)
        }else{
            gammalist[[i]] <- gammaKpar(dfA[,varnames[i]], dfB[,varnames[i]], gender = FALSE, n.cores = n.cores)
        }
    }
    #add comparisons of all missing indices to all nonmissing, and to each other
    na.lists = lapply(gammalist, '[[','nas')
    allinds = list(1:nrow(dfA), 1:nrow(dfB))
    
    for(i in 1:length(varnames)){
      sublist = na.lists[[i]]
      
      nas.present = sapply(sublist, length) > 0
      if(dedupe){
        if(!nas.present[[1]]) next
        na.inds = sublist[[1]]
        comp.inds = allinds[[1]][-na.inds]
        
        templist = list(
          list(comp.inds, na.inds),
          list(na.inds, na.inds)
        )
        
      }else{
        templist = vector(mode = 'list', length = sum(nas.present))
        for(j in 1:2){
          if(!nas.present[j]) next
          na.inds = sublist[[j]]
          comp.inds = allinds[[3-j]]
          templist[[j]] = vector(mode = 'list', length = 2)
          templist[[j]][[j]] = na.inds
          templist[[j]][[3-j]] = comp.inds
          }
      }
      gammalist[[i]]$matches2 = c(gammalist[[i]]$matches2, templist)
    }
    
    
    cat("\n")
    # if(combine.method == 'AND') pat_template = c(NA,2)
    # if(combine.method == 'OR') pat_template = c(NA,0,2)
    #   
    # str = paste0('expand.grid(',paste(rep('pat_template',length(varnames)),collapse = ','),')')
    # patts = eval(parse(text = str))
    # patts = patts[apply(patts,1,sum,na.rm=T)>0,]
    # patts = cbind(patts,100)
    # zeta.j = rep(0,nrow(patts))
    # weights = rep(0,nrow(patts))
    # 
    # block_groups = vector(mode = 'list', length = nrow(patts))
    # for(i in 1:nrow(patts)){
    #   zeta.j_ = zeta.j
    #   zeta.j_[i] = 1
    #   weights_ = weights
    #   weights_[i] = 100
    #   
    #   em = list(patterns.w = cbind(patts,weights),zeta.j = zeta.j_)
    #   block_groups[[i]] = matchesLink(gammalist, nobs.a = nrow(dfA), nobs.b = nrow(dfB),
    #                                   em = em, thresh = 1,
    #                                   n.cores = n.cores, dedupe = dedupe)
    # } Problem - too many record pairs!
    
    
    
    ## --------------
    ## Combine blocks
    ## --------------
    cat("Combining blocked variables for final blocking assignments.\n\n")
    combineblocks_out <- combineBlocks(gammalist)
    indlist_a <- which(combineblocks_out$dfA.block, arr.ind = T)
    indlist_a <- tapply(indlist_a[,1],indlist_a[,2], function(x) x)
    indlist_b <- which(combineblocks_out$dfB.block, arr.ind = T)
    indlist_b <- tapply(indlist_b[,1],indlist_b[,2], function(x) x)
    
    ## Clean up
    blocklist_out <- vector(mode = "list", length = length(indlist_a))
    for(i in 1:length(blocklist_out)){
        blocklist_out[[i]] <- list(dfA.inds = indlist_a[[i]], dfB.inds = indlist_b[[i]])
    }
    if(dedupe){
      lengths = sapply(blocklist_out,function(l) length(l$dfA.inds) * length(l$dfB.inds))
      blocklist_out = blocklist_out[lengths>1]
    }
    names(blocklist_out) <- paste0("block.", 1:length(blocklist_out))
    class(blocklist_out) <- "fastLink.block"
    return(blocklist_out)
    
}

is.prime <- function(x)
  vapply(x, function(y) sum(y / 1:y == y %/% 1:y), integer(1L)) == 2L


combineBlocks <- function(blocklist){
  if(length(blocklist) > 1){
    primes = (2:500)[is.prime(2:500)]
    
    allindsA = vector(mode = 'list', length = length(blocklist))
    allindsB = vector(mode = 'list', length = length(blocklist))
    
    for(i in 1:length(blocklist)){
      allindsA[[i]] = do.call(rbind,lapply(1:length(blocklist[[i]][[1]]), function(j) cbind(blocklist[[i]][[1]][[j]][[1]],j)))
      allindsB[[i]] = do.call(rbind,lapply(1:length(blocklist[[i]][[1]]), function(j) cbind(blocklist[[i]][[1]][[j]][[2]],j)))
    }
    
    levs = length(blocklist)
    
    block_combos = lapply(levs,function(l){
      domats = RcppAlgos::comboGeneral(1:length(blocklist),l)
      
      blockmats = vector(mode = 'list', length = nrow(domats))
      
      for(i in 1:nrow(domats)){
        i_ = domats[i,]
        overlapsA = as.data.frame(allindsA[[i_[1]]])
        for(n in 2:length(i_))
          overlapsA = merge(overlapsA, as.data.frame(allindsA[[i_[n]]]), by = 'V1')
        
        overlapsA = overlapsA[,-1]
        
        idA = as.matrix(overlapsA) %*% log(primes[1:length(i_)])
        
        overlapsA = overlapsA[!duplicated(idA),]
        idA = idA[!duplicated(idA)]
        
        overlapsB = as.data.frame(allindsB[[i_[1]]])
        for(n in 2:length(i_))
          overlapsB = merge(overlapsB, as.data.frame(allindsB[[i_[n]]]), by = 'V1')
        
        overlapsB = overlapsB[,-1]
        
        idB = as.matrix(overlapsB) %*% log(primes[1:length(i_)])
        overlapsB = overlapsA[!duplicated(idB),]
        idB = idB[!duplicated(idB)]
        
        blockmats[[i]] = overlapsA[which(idA %in% idB),]
        
      }
      list(domat = domats, blockmats = blockmats)
    })
    
    
    blkgrps = block_combos[[1]]$blockmats[[1]]
  }else blkgrps = matrix(1:length(blocklist[[1]][[1]]),ncol = 1)
  
  indsA_out <- vector(mode = "list", length = nrow(blkgrps))
  indsB_out <- vector(mode = "list", length = nrow(blkgrps))
  indsA <- vector(mode = "list", length = ncol(blkgrps))
  indsB <- vector(mode = "list", length = ncol(blkgrps))
    
  for(i in 1:nrow(blkgrps)){
      
    for(j in 1:ncol(blkgrps)){
      indsA[[j]] <- blocklist[[j]][[1]][[blkgrps[i,j]]][[1]]
      indsB[[j]] <- blocklist[[j]][[1]][[blkgrps[i,j]]][[2]]
    }
    indsA_intersect <- Reduce(intersect, indsA)
    if(length(indsA_intersect) == 0) next
    indsB_intersect <- Reduce(intersect, indsB)
    if(length(indsB_intersect) == 0) next
    
    indsA_out[[i]] <- cbind(indsA_intersect, i)
    indsB_out[[i]] <- cbind(indsB_intersect, i)
    
  }
    
  indsA_out <- do.call(rbind, indsA_out)
  indsB_out <- do.call(rbind, indsB_out)
  matA_out <- sparseMatrix(i = indsA_out[,1], j = indsA_out[,2])
  matB_out <- sparseMatrix(i = indsB_out[,1], j = indsB_out[,2])
    
  out <- list(dfA.block = matA_out, dfB.block = matB_out)
  class(out) <- "fastLink.block"
  return(out)
  }
  
  
  


kmeansBlock <- function(vecA, vecB, nclusters, iter.max, n.cores){
    
    if(class(vecA) == "factor"){
        vecA <- as.character(vecA)
    }
    if(class(vecB) == "factor"){
        vecB <- as.character(vecB)
    }
    
    ## Clean and combine
    vec <- c(vecA, vecB)
    setid <- c(rep("A", length(vecA)), rep("B", length(vecB)))
    dims <- as.numeric(as.factor(vec))
    
    ## Run kmeans
    km.out <- kmeans(na.omit(dims), centers = nclusters, iter.max = iter.max)
    cluster <- rep(NA, length(vec))
    cluster[which(!is.na(vec))] <- km.out$cluster
    
    ## Run gammaKpar
    out <- gammaKpar(cluster[setid == "A"], cluster[setid == "B"], gender = FALSE, n.cores = n.cores)
    return(out)
    
}

#' consolidate_blocks
#'
#' Allows for disjunctive blocking using multiple criteria, i.e. ONE OR MORE blocking conditions must be met
#'
#' @usage consolidate_blocks(blocklist1, blocklist2, ..., dedupe = F)
#' @param ... Any number of blocklists output by \code{blockData()}, separated by commas
#' @param dedupe Boolean indicator as to whether the blocks represent an internal linkage problem (\code{TRUE}) or linkage between two datasets (\code{FALSE})
#'
#' @return A list of non-overlapping blocks covering the same record pairs as the original inputs
#'
#' @usage
#' \dontrun{
#' blocks_sex_OR_yob <- consolidate_blocks(blocks_sex, blocks_yob, dedupe = T)
#' }
#'
#' @export

consolidate_blocks = function(...,dedupe = F){
  blocklist = unique(c(...))
  
  #step 1: combine blocks that fully overlap in dimension A
  {
    Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
      cbind(blocklist[[i]]$dfA.inds,i)
    }))
    
    Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
    
    pwrs = seq(0,10,length.out = nrow(Ablockmat))
    
    blockkey = matrix(2^pwrs,nrow = 1)
    
    Akeys = (blockkey %*% Ablockmat)[1,]
    
    dup.keys = unique(Akeys[duplicated(Akeys)])
    
    inds = which(Akeys %in% dup.keys)
    if(length(inds) > 0){
      consolidate_a = tapply(inds,Akeys[inds],function(x)x)
      
      cat('Step 1: Combining fully overlapped dfA blocks \n')
      pb = txtProgressBar(0,length(consolidate_a))
      i = 1
      
      for(v in consolidate_a){
        templist = list(dfA.inds = blocklist[[v[1]]]$dfA.inds,
                        dfB.inds = unique(do.call(c,lapply(v, function(i) blocklist[[i]]$dfB.inds))))
        
        blocklist[[v[1]]] = templist
        
        setTxtProgressBar(pb,i)
        i = i+1
      }
      close(pb)
      dumpinds = unlist(lapply(consolidate_a, '[',-1))
      blocklist = blocklist[-dumpinds]
    }
  }
  #step 2: combine blocks that fully overlap in dimension B
  {
    Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
      cbind(blocklist[[i]]$dfB.inds,i)
    }))
    
    Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
    
    pwrs = seq(0,10,length.out = nrow(Bblockmat))
    
    blockkey = matrix(2^pwrs,nrow = 1)
    
    Bkeys = (blockkey %*% Bblockmat)[1,]
    
    dup.keys = unique(Bkeys[duplicated(Bkeys)])
    
    inds = which(Bkeys %in% dup.keys)
    if(length(inds)>0){
      consolidate_b = tapply(inds,Bkeys[inds],function(x)x)
      
      cat('Step 2: Combining fully overlapped dfB blocks \n')
      pb = txtProgressBar(0,length(consolidate_b))
      i = 1
      for(v in consolidate_b){
        templist = list(dfA.inds = unique(do.call(c,lapply(v, function(i) blocklist[[i]]$dfA.inds))),
                        dfB.inds = blocklist[[v[1]]]$dfB.inds)
        
        blocklist[[v[1]]] = templist
        
        setTxtProgressBar(pb,i)
        i = i+1
      }
      close(pb)
      dumpinds = unlist(lapply(consolidate_b, '[',-1))
      blocklist = blocklist[-dumpinds]
    }
    
  }
  ##Step 3: reduce blocks with dimension A fully contained in another block
  swtch1 = 1
  swtch2 = 1
  while(swtch1 + swtch2 > 0){
    Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
      cbind(blocklist[[i]]$dfA.inds,i)
    }))
    
    Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
      cbind(blocklist[[i]]$dfB.inds,i)
    }))  
    
    Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
    Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
    nA = sapply(blocklist,function(l) length(l$dfA.inds))
    
    Acomps = triu(t(Ablockmat) %*% Ablockmat, k = 1)
    Bcomps = triu(t(Bblockmat) %*% Bblockmat, k = 1)
    
    nAcomp = (Acomps * (Bcomps>0))@x
    
    overlaps = which(Acomps * Bcomps>0,arr.ind = T)
    nAs = overlaps
    nAs[] = nA[nAs]
    
    todo = matrix(overlaps[pmin(nAs[,1],nAs[,2])==nAcomp,],ncol = 2)
    if(nrow(todo)>0){
      swtch1 = 1
      
      donors = apply(nAs[pmin(nAs[,1],nAs[,2])==nAcomp,], 1, which.min)
      todo = cbind(todo[cbind(1:nrow(todo), donors)],todo[cbind(1:nrow(todo), 3-donors)])
      
      cat('Step 3: Reducing partially contained dfA blocks \n')
      pb = txtProgressBar(0,nrow(todo))
      del_list = NULL
      
      for(i in 1:nrow(todo)){
        donor = blocklist[[todo[i,1]]]; receiver = blocklist[[todo[i,2]]]
        
        b_olap = intersect(donor$dfB.inds, receiver$dfB.inds)
        b_resid = setdiff(donor$dfB.inds, b_olap)
        
        if(length(b_resid) == 0){del_list = c(del_list, todo[i,1]); next}
        
        blocklist[[todo[i,1]]]$dfB.inds = b_resid
        setTxtProgressBar(pb, i)
      }
      close(pb)
      del_list = unique(del_list)
      if(!is.null(del_list)) blocklist = blocklist[-del_list]
      
    }else{
      swtch1 = 0
    }
    
    ###Step 4: reduce blocks with dimension B fully contained in another block
    
    Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
      cbind(blocklist[[i]]$dfA.inds,i)
    }))
    
    Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
      cbind(blocklist[[i]]$dfB.inds,i)
    }))  
    
    Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
    Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
    
    nB = sapply(blocklist,function(l) length(l$dfB.inds))
    
    Acomps = triu(t(Ablockmat) %*% Ablockmat, k = 1)
    Bcomps = triu(t(Bblockmat) %*% Bblockmat, k = 1)
    
    nBcomp = (Bcomps * (Acomps > 0))@x
    
    overlaps = which(Acomps * Bcomps>0,arr.ind = T)
    nBs = overlaps
    nBs[] = nB[nBs]
    
    todo = matrix(overlaps[pmin(nBs[,1],nBs[,2])==nBcomp,],ncol = 2)
    
    if(nrow(todo) > 0){
      swtch2 = 1
      donors = apply(matrix(nBs[pmin(nBs[,1],nBs[,2])==nBcomp,],ncol = 2), 1, which.min)
      todo = cbind(todo[cbind(1:nrow(todo), donors)],todo[cbind(1:nrow(todo), 3-donors)])
      
      cat('Step 4: Reducing partially contained dfB blocks \n')
      pb = txtProgressBar(0,nrow(todo))
      del_list = NULL
      
      for(i in 1:nrow(todo)){
        donor = blocklist[[todo[i,1]]]; receiver = blocklist[[todo[i,2]]]
        
        a_olap = intersect(donor$dfA.inds, receiver$dfA.inds)
        a_resid = setdiff(donor$dfA.inds, a_olap)
        
        if(length(a_resid) == 0){del_list = c(del_list, todo[i,1]); next}
        
        blocklist[[todo[i,1]]]$dfA.inds = a_resid
        setTxtProgressBar(pb, i)
      }
      close(pb)
      del_list = unique(del_list)
      if(!is.null(del_list)) blocklist = blocklist[-del_list]
    }else{
      swtch2 = 0
    }
    
  }
  
  
  #Step 5a: if dedupe is true, go ahead and separate out blocks on the diagonal
  
  if(dedupe){
    resolve_diag = sapply(blocklist,function(l){
      !identical(l$dfA.inds,l$dfB.inds) & any(l$dfA.inds %in% l$dfB.inds)
    })
    if(sum(resolve_diag) > 0){
      templists = vector(mode = 'list', length = sum(resolve_diag))
      
      cat('Step 5a: separating blocks along diagonal from rectangular blocks \n')
      pb = txtProgressBar(0,sum(resolve_diag))
      
      i = 1
      for(lst in blocklist[resolve_diag]){
        diag_inds = intersect(lst$dfA.inds,lst$dfB.inds)
        
        rem_A = setdiff(lst$dfA.inds,diag_inds)
        rem_B = setdiff(lst$dfB.inds,diag_inds)
        
        indicator = (length(rem_A) > 0) + 2*(length(rem_B) > 0)
        
        templists[[i]] = list(list(dfA.inds = diag_inds, dfB.inds = diag_inds))
        
        if(indicator == 1)
          templists[[i]] = list(list(dfA.inds = diag_inds, dfB.inds = diag_inds), 
                                list(dfA.inds = rem_A, dfB.inds = diag_inds))
        if(indicator == 2)
          templists[[i]] = list(list(dfA.inds = diag_inds, dfB.inds = diag_inds), 
                                list(dfA.inds = diag_inds, dfB.inds = rem_B))
        
        if(indicator == 3){
          templists[[i]] = list(list(dfA.inds = diag_inds, dfB.inds = diag_inds), 
                                list(dfA.inds = rem_A, dfB.inds = diag_inds),
                                list(dfA.inds = lst$dfA.inds, dfB.inds = rem_B))
        }
        
        setTxtProgressBar(pb,i)
        i = i+1
      }
      close(pb)
      
      templists = unlist(templists,recursive = F)
      
      blocklist = blocklist[!resolve_diag]
      blocklist = c(blocklist,templists)
    }
  }
  
  #Step 5: Resolve overlapping blocks
  
  Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
    cbind(blocklist[[i]]$dfA.inds,i)
  }))
  
  Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
  
  Acomps = triu(t(Ablockmat) %*% Ablockmat, k = 1)
  
  Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
    cbind(blocklist[[i]]$dfB.inds,i)
  }))
  
  Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
  
  Bcomps = triu(t(Bblockmat) %*% Bblockmat, k = 1)
  
  overlaps = which(Bcomps * Acomps > 0, arr.ind = T)
  pass = 1
  while(nrow(overlaps)>0){
    nrp = sapply(blocklist,function(l) as.numeric(length(l$dfB.inds)) * as.numeric(length(l$dfA.inds)))
    olap_nrp = (Bcomps * Acomps)@x
    
    nrps = overlaps
    nrps[] = nrp[nrps]
    
    donors = apply(nrps,1,which.min)
    
    
    
    overlaps = cbind(overlaps[cbind(1:nrow(overlaps),donors)],overlaps[cbind(1:nrow(overlaps),3-donors)])
    
    donors = unique(overlaps[,1])
    
    donors = donors[order(nrp[donors], decreasing = T)]
    
    cat('Step 5: Resolving remaining overlap between blocks, pass', pass,'\n')
    
    i = 1
    pb = txtProgressBar(0,length(donors))
    for(d in donors){
      templist = blocklist[[d]]
      
      todo = matrix(overlaps[overlaps[,1] == d,],ncol = 2)
      todo = matrix(todo[order(olap_nrp[overlaps[,1] == d], decreasing = T),],ncol =2)
      
      olaps_a = lapply(1:nrow(todo), function(i){
        l2 = blocklist[[todo[i,2]]]
        intersect(templist$dfA.inds,l2$dfA.inds)
      })
      
      olaps_b = lapply(1:nrow(todo), function(i){
        l2 = blocklist[[todo[i,2]]]
        intersect(templist$dfB.inds,l2$dfB.inds)
      })
      
      templist = list(templist)
      
      for(ii in 1:nrow(todo)){
        do = sapply(templist, function(l) any(olaps_a[[ii]] %in% l$dfA.inds) && any(olaps_b[[ii]] %in% l$dfB.inds))
        if(!any(do)) next
        
        rem_A = lapply(templist[do], function(l) setdiff(l$dfA.inds, olaps_a[[ii]]))
        ol_A = lapply(templist[do], function(l) intersect(l$dfA.inds, olaps_a[[ii]]))
        rem_B = lapply(templist[do], function(l) setdiff(l$dfB.inds, olaps_b[[ii]]))
        ol_B = lapply(templist[do], function(l) intersect(l$dfB.inds, olaps_b[[ii]]))
        
        
        templist[do] = lapply(1:sum(do), function(i){
          l = templist[do][[i]]
          dupe_i = identical(l$dfA.inds,l$dfB.inds)
          dupe_j = identical(ol_A[[i]],ol_B[[i]])
          
          indicator = (length(rem_A[[i]])>0) + 2*(length(rem_B[[i]])>0)
          
          if(indicator == 0) return(NULL)
          if(indicator == 1) return(list(list(
            dfA.inds = rem_A[[i]], dfB.inds = l$dfB.inds
          )))
          if(indicator == 2) return(list(list(
            dfA.inds = l$dfA.inds, dfB.inds = rem_B[[i]]
          )))
          if(indicator == 3){
            indicator2 = dupe_i + 2 * dupe_j
            
            if(indicator2 %in% c(0,2)){
              res1 = list(dfA.inds = ol_A[[i]], dfB.inds = rem_B[[i]])
              res2 = list(dfA.inds = rem_A[[i]], dfB.inds = l$dfB.inds)
              return(list(res1,res2)[c(length(rem_B[[i]])>0, length(rem_A[[i]]) > 0)])
            }
            
            if(indicator2 == 1){
              res1 = list(dfA.inds = ol_A[[i]], dfB.inds = ol_A[[i]])
              res2 = list(dfA.inds = rem_A[[i]], dfB.inds = rem_A[[i]])
              res3 = list(dfA.inds = ol_A[[i]], dfB.inds = setdiff(l$dfB.inds,c(rem_B[[i]],rem_A[[i]])))
              return(list(res1,res2,res3)[c(length(ol_A[[i]])>1, length(rem_A[[i]])>1,length(res3$dfB.inds)>0)])
            }
            
            if(indicator2 == 3){
              res1 = list(dfA.inds = rem_A[[i]], dfB.inds = rem_A[[i]])
              res2 = list(dfA.inds = ol_A[[i]], dfB.inds = rem_B[[i]])
              return(list(res1,res2)[c(length(rem_A[[i]])>1, length(rem_B[[i]])>0)])
            }
          }})
        templist = c(unlist(templist[do], recursive = F),templist[!do])
        templist = templist[sapply(templist,length)>0]
        
      }
      blocklist = c(blocklist, templist)
      
      setTxtProgressBar(pb,i)
      i = i+1
    }
    close(pb)
    
    blocklist = blocklist[-donors]
    
    # Make sure all overlap has been eliminated
    Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
      cbind(blocklist[[i]]$dfA.inds,i)
    }))
    
    Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
    
    Acomps = triu(t(Ablockmat) %*% Ablockmat,k = 1)
    
    Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
      cbind(blocklist[[i]]$dfB.inds,i)
    }))
    
    Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
    
    Bcomps = triu(t(Bblockmat) %*% Bblockmat, k = 1)
    
    overlaps = which(Bcomps * Acomps > 0, arr.ind = T)
    pass = pass + 1
    cat(sum(Bcomps * Acomps),' overlapping record pairs remain \n')
  }
  #Step 6, remove redundant record pairs that got scrambled A-B if dedupe is true
  if(dedupe){
    
    ABcomps = triu(t(Ablockmat) %*% Bblockmat, k = 1)
    BAcomps = triu(t(Bblockmat) %*% Ablockmat, k = 1)
    
    overlaps = which(ABcomps * BAcomps > 0, arr.ind = T)
    
    pass = 1
    while(nrow(overlaps)>0){
      nrp = sapply(blocklist,function(l) as.numeric(length(l$dfB.inds)) * as.numeric(length(l$dfA.inds)))
      olap_nrp = (ABcomps * BAcomps)@x
      
      nrps = overlaps
      nrps[] = nrp[nrps]
      
      donors = apply(nrps,1,which.min)
      
      overlaps = cbind(overlaps[cbind(1:nrow(overlaps),donors)],overlaps[cbind(1:nrow(overlaps),3-donors)])
      
      donors = unique(overlaps[,1])
      
      donors = donors[order(nrp[donors], decreasing = T)]
      
      cat('Step 6: Resolving redundancies where indices are switched between dfA and dfB, pass', pass,'\n')
      
      i = 1
      pb = txtProgressBar(0,length(donors))
      for(d in donors){
        templist = blocklist[[d]]
        
        todo = matrix(overlaps[overlaps[,1] == d,],ncol = 2)
        todo = matrix(todo[order(olap_nrp[overlaps[,1] == d], decreasing = T),],ncol =2)
        
        olaps_ab = lapply(1:nrow(todo), function(i){
          l2 = blocklist[[todo[i,2]]]
          intersect(templist$dfA.inds,l2$dfB.inds)
        })
        
        olaps_ba = lapply(1:nrow(todo), function(i){
          l2 = blocklist[[todo[i,2]]]
          intersect(templist$dfB.inds,l2$dfA.inds)
        })
        
        templist = list(templist)
        
        for(ii in 1:nrow(todo)){
          do = sapply(templist, function(l) any(olaps_ab[[ii]] %in% l$dfA.inds) && any(olaps_ba[[ii]] %in% l$dfB.inds))
          if(!any(do)) next
          
          rem_A = lapply(templist[do], function(l) setdiff(l$dfA.inds, olaps_ab[[ii]]))
          ol_A = lapply(templist[do], function(l) intersect(l$dfA.inds, olaps_ab[[ii]]))
          rem_B = lapply(templist[do], function(l) setdiff(l$dfB.inds, olaps_ba[[ii]]))
          ol_B = lapply(templist[do], function(l) intersect(l$dfB.inds, olaps_ba[[ii]]))
          
          templist[do] = lapply(1:sum(do), function(i){
            l = templist[do][[i]]
            dupe_i = identical(l$dfA.inds,l$dfB.inds)
            dupe_j = identical(ol_A[[i]],ol_B[[i]])
            
            indicator = (length(rem_A[[i]])>0) + 2*(length(rem_B[[i]])>0)
            
            if(indicator == 0) return(NULL)
            if(indicator == 1) return(list(list(
              dfA.inds = rem_A[[i]], dfB.inds = l$dfB.inds
            )))
            if(indicator == 2) return(list(list(
              dfA.inds = l$dfA.inds, dfB.inds = rem_B[[i]]
            )))
            if(indicator == 3){
              indicator2 = dupe_i + 2 * dupe_j
              
              if(indicator2 %in% c(0,2)){
                res1 = list(dfA.inds = ol_A[[i]], dfB.inds = rem_B[[i]])
                res2 = list(dfA.inds = rem_A[[i]], dfB.inds = l$dfB.inds)
                return(list(res1,res2)[c(length(rem_B[[i]])>0, length(rem_A[[i]]) > 0)])
              }
              
              if(indicator2 == 1){
                res1 = list(dfA.inds = ol_A[[i]], dfB.inds = ol_A[[i]])
                res2 = list(dfA.inds = rem_A[[i]], dfB.inds = rem_A[[i]])
                res3 = list(dfA.inds = ol_A[[i]], dfB.inds = setdiff(l$dfB.inds,c(rem_B[[i]],rem_A[[i]])))
                return(list(res1,res2,res3)[c(length(ol_A[[i]])>1, length(rem_A[[i]])>1,length(res3$dfB.inds)>0)])
              }
              
              if(indicator2 == 3){
                res1 = list(dfA.inds = rem_A[[i]], dfB.inds = rem_A[[i]])
                res2 = list(dfA.inds = ol_A[[i]], dfB.inds = rem_B[[i]])
                return(list(res1,res2)[c(length(rem_A[[i]])>1, length(rem_B[[i]])>0)])
              }
            }})
          templist = c(unlist(templist[do], recursive = F),templist[!do])
          templist = templist[sapply(templist,length)>0]
          
        }
        blocklist = c(blocklist, templist)
        
        setTxtProgressBar(pb,i)
        i = i+1
      }
      close(pb)
      
      blocklist = blocklist[-donors]
      
      # Make sure all overlap has been eliminated
      Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
        cbind(blocklist[[i]]$dfA.inds,i)
      }))
      
      Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
      
      Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
        cbind(blocklist[[i]]$dfB.inds,i)
      }))
      
      Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
      
      ABcomps = triu(t(Ablockmat) %*% Bblockmat, k = 1)
      BAcomps = triu(t(Bblockmat) %*% Ablockmat, k = 1)
      
      overlaps = which(ABcomps * BAcomps > 0, arr.ind = T)
      
      pass = pass + 1
      cat(sum(Bcomps * Acomps),' redundant record pairs remain \n')
    }
  }
  return(blocklist)
}

