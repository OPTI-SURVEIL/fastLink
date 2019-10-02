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
#' @param inverse.block Vector of variable names for which blocked records should NOT agree. Used during posterior match recalculation with fuzzy name matching, for instance
#' @param na.block Vector of variable names for blocked records should have missing comparisons. Used during posterior match recalculation with fuzzy name matching, for instance
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
                      n.cores = NULL, inverse.block = NULL, na.block = NULL){
  
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
        blocktype <- ifelse(varnames[i] %in% window.block, "window", 
                            ifelse(varnames[i] %in% kmeans.block, "k-means", 
                                   ifelse(varnames[i] %in% inverse.block, "inverse", "exact")))
        cat("    Blocking variable", varnames[i], "using", blocktype, "blocking.\n")
        if(varnames[i] %in% window.block){
            gammalist[[i]] <- gammaNUMCK2par(dfA[,varnames[i]], dfB[,varnames[i]], n.cores = n.cores, cut.a = window.size,dedupe = dedupe)
        }else if(varnames[i] %in% kmeans.block){
            gammalist[[i]] <- kmeansBlock(dfA[,varnames[i]], dfB[,varnames[i]],
                                          nclusters = nclusters,
                                          iter.max = iter.max,
                                          n.cores = n.cores)
        }else{
            gammalist[[i]] <- gammaKpar(dfA[,varnames[i]], dfB[,varnames[i]], gender = FALSE, n.cores = n.cores,
                                        inverse = blocktype == 'inverse')
        }
    }
    #add comparisons of all missing indices to all nonmissing, and to each other, only if inverse block is null
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

    ## --------------
    ## Combine blocks
    ## --------------
    cat("Combining blocked variables for final blocking assignments.\n\n")
    combineblocks_out <- combineBlocks(gammalist, c(nrow(dfA), nrow(dfB)))
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
      length1 = sapply(blocklist_out,function(l) length(l$dfA.inds) == 1 && length(l$dfB.inds) == 1 && identical(l$dfA.inds,l$dfB.inds))
      blocklist_out = blocklist_out[!length1]
    }
    names(blocklist_out) <- paste0("block.", 1:length(blocklist_out))
    class(blocklist_out) <- "fastLink.block"
    return(blocklist_out)
    
}

is.prime <- function(x)
  vapply(x, function(y) sum(y / 1:y == y %/% 1:y), integer(1L)) == 2L


combineBlocks <- function(blocklist,df.dims){
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
  matA_out <- sparseMatrix(i = indsA_out[,1], j = indsA_out[,2],dims = c(df.dims[1],nrow(blkgrps)))
  matB_out <- sparseMatrix(i = indsB_out[,1], j = indsB_out[,2],dims = c(df.dims[2],nrow(blkgrps)))
    
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
  
  Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
    cbind(blocklist[[i]]$dfA.inds,i)
  }))
  
  Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
  
  Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
    cbind(blocklist[[i]]$dfB.inds,i)
  }))
  
  Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
  
  olap = triu(t(Ablockmat) %*% Ablockmat * (t(Bblockmat) %*% Bblockmat),k = 1)
  
  swtch = T
  
  pass = 1
  
  recalc = F
  
  mode = 'same_dim'
  
  while(swtch){
    
    #step 1: combine blocks that fully overlap in dimension A
    swtch1 = T
    while(swtch1){
      while(swtch1){
        if(recalc){
          Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfA.inds,i)
          }))
          
          Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
          
          Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfB.inds,i)
          }))
          
          Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
        }
        
        pwrs = seq(0,10,length.out = nrow(Ablockmat))
        
        blockkey = matrix(2^pwrs,nrow = 1)
        
        Akeys = (blockkey %*% Ablockmat)[1,]
        
        dup.keys = unique(Akeys[duplicated(Akeys)])
        
        inds = which(Akeys %in% dup.keys)
        if(length(inds) > 0){
          cat('Pass', pass,'Step 1: Combining fully overlapped dfA blocks \n')
          blocklist = full_block_merge(inds,Akeys[inds],blocklist,dim = 1)
          #if(!checkres(blocklist,keys)) return(blocklist)
          recalc = T
          swtch1 = T
        }else{
          recalc = F
          swtch1 = F
        }
        
        #step 1: combine blocks that fully overlap in dimension B
        
        if(recalc){
          Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfA.inds,i)
          }))
          
          Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
          
          Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfB.inds,i)
          }))
          
          Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
        }
        
        pwrs = seq(0,10,length.out = nrow(Bblockmat))
        
        blockkey = matrix(2^pwrs,nrow = 1)
        
        Bkeys = (blockkey %*% Bblockmat)[1,]
        
        dup.keys = unique(Bkeys[duplicated(Bkeys)])
        
        inds = which(Bkeys %in% dup.keys)
        if(length(inds)>0){
          cat('Pass', pass, 'Step 1: Combining fully overlapped dfB blocks \n')
          blocklist = full_block_merge(inds,Bkeys[inds],blocklist,dim = 2)
          #if(!checkres(blocklist,keys)) return(blocklist)
          recalc = T
          swtch1 = T
        }else{
          recalc = F
          swtch1 = F
        }
      }
      
      #Step 1a: If dedupe is true, combine blocks that fully overlap in one dimension, but that have been switched A-B
      if(dedupe){
        dupe_blocks = sapply(blocklist,function(x) identical(x[[1]],x[[2]]))
        
        if(recalc){
          Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfA.inds,i)
          }))
          
          Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
          
          Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfB.inds,i)
          }))
          
          Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
        }
        pwrs = seq(0,10,length.out = nrow(Ablockmat))
        
        blockkey = matrix(2^pwrs,nrow = 1)
        
        Akeys = (blockkey %*% Ablockmat)[1,!dupe_blocks]
        Bkeys = (blockkey %*% Bblockmat)[1,!dupe_blocks]
        
        dup.keys = intersect(Bkeys,Akeys)
        
        indsA = which(Akeys %in% dup.keys)
        indsB = which(Bkeys %in% dup.keys)
        
        blocksA = which(!dupe_blocks)[indsA]; blocksB = which(!dupe_blocks)[indsB]
        if(length(indsA)>0 & length(indsB) > 0){
          cat('Pass', pass,'Step 1: Combining blocks that fully overlap A to B \n')
          #debug(full_block_merge)
          #if(pass > 3) debug(full_block_merge)
          blocklist = full_block_merge(c(blocksA,blocksB),c(Akeys[indsA],Bkeys[indsB]),blocklist,dim = 1, AtoB = T)
          #if(!checkres(blocklist,keys)) return(blocklist)
          recalc = T
          swtch1 = T
        }else{
          recalc = F
          swtch1 = F
        }
      }
    }
    #Step 2
    
    if(recalc){
      Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
        cbind(blocklist[[i]]$dfA.inds,i)
      }))
      
      Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
        cbind(blocklist[[i]]$dfB.inds,i)
      }))  
      
      Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
      Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
    }
    
    nA = sapply(blocklist,function(l) length(l$dfA.inds))
    
    Acomps = triu(t(Ablockmat) %*% Ablockmat, k = 1)
    Bcomps = triu(t(Bblockmat) %*% Bblockmat, k = 1)
    
    nAcomp = (Acomps * (Bcomps>0))@x
    
    overlaps = which(Acomps * Bcomps>0,arr.ind = T)
    nAs = overlaps
    nAs[] = nA[nAs]
    
    do = pmin(nAs[,1],nAs[,2])==nAcomp
    
    todo = matrix(overlaps[do,],ncol = 2)
    nAs = matrix(nAs[do,], ncol = 2)
    
    if(nrow(todo)>0){
      do = sapply(todo[,1], function(i) length(blocklist[[i]][[1]]) > length(blocklist[[i]][[2]]))
      todo = matrix(todo[do,],ncol = 2)
      nAs = matrix(nAs[do,], ncol = 2)
    } 
    
    if(nrow(todo)>0){
      cat('Pass', pass,'Step 2: Reducing partially contained dfA blocks \n')
      blocklist = block_reduce(todo,nAs,blocklist,1)
      recalc = T
    }else{
      recalc = F
    }
    
    ###Step 2: reduce blocks with dimension B fully contained in another block ONLY IF the donor block will become more square
    
    if(recalc){
      Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
        cbind(blocklist[[i]]$dfA.inds,i)
      }))
      
      Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
        cbind(blocklist[[i]]$dfB.inds,i)
      }))  
      
      Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
      Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
      
      Acomps = triu(t(Ablockmat) %*% Ablockmat, k = 1)
      Bcomps = triu(t(Bblockmat) %*% Bblockmat, k = 1)
    }
    
    nB = sapply(blocklist,function(l) length(l$dfB.inds))
    
    nBcomp = (Bcomps * (Acomps > 0))@x
    
    overlaps = which(Acomps * Bcomps>0,arr.ind = T)
    nBs = overlaps
    nBs[] = nB[nBs]
    
    do = pmin(nBs[,1],nBs[,2])==nBcomp
    
    todo = matrix(overlaps[do,],ncol = 2)
    nBs = matrix(nBs[do,],ncol = 2)
    if(nrow(todo) > 0){
      do = sapply(todo[,1], function(i) length(blocklist[[i]][[2]]) > length(blocklist[[i]][[1]]))
      todo = matrix(todo[do,],ncol = 2)
      nBs = matrix(nBs[do,], ncol = 2)
    }
    
    if(nrow(todo) > 0){
      cat('Pass', pass,'Step 2: Reducing partially contained dfB blocks \n')
      
      blocklist = block_reduce(todo,nBs,blocklist,2)
      #if(!checkres(blocklist,keys)) return(blocklist)
      recalc = T
    }else{
      recalc = F
    }
    
    if(dedupe){
      #Step 2: Reduce blocks with dimension A fully contained within dimension B of some other block
      
      if(recalc){
        Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
          cbind(blocklist[[i]]$dfA.inds,i)
        }))
        
        Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
          cbind(blocklist[[i]]$dfB.inds,i)
        }))  
        
        Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
        Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
      }
      
      nA = sapply(blocklist,function(l) length(l$dfA.inds))
      
      ABcomps = triu(t(Ablockmat) %*% Bblockmat, k = 1)
      BAcomps = triu(t(Bblockmat) %*% Ablockmat, k = 1)
      
      nABcomp = (ABcomps * (BAcomps > 0))@x
      
      overlaps = which(ABcomps * BAcomps>0,arr.ind = T)
      nAs = nA[overlaps[,1]]
      
      todo = matrix(overlaps[nAs==nABcomp,] ,ncol = 2)
      nAs = nAs[nAs==nABcomp]
      
      if(nrow(todo)>0){
        do = sapply(todo[,1], function(i) length(blocklist[[i]][[1]]) > length(blocklist[[i]][[2]]))
        todo = matrix(todo[do,],ncol = 2)
        nAs = nAs[do]
      }
      
      if(nrow(todo) > 0){
        cat('Pass', pass,'Step 2: Reducing blocks with dfA fully contained in another block\'s dfB \n')
        
        blocklist = block_reduce(todo,nAs,blocklist,1,crossdim = T)
        #if(!checkres(blocklist,keys)) return(blocklist)
        recalc = T
      }else{
        recalc = F
      }
      
      #Step 2: Reduce blocks with dimension B fully contained within dimension A of some other block
      
      if(recalc){
        Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
          cbind(blocklist[[i]]$dfA.inds,i)
        }))
        
        Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
          cbind(blocklist[[i]]$dfB.inds,i)
        }))  
        
        Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
        Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
        
        ABcomps = triu(t(Ablockmat) %*% Bblockmat, k = 1)
        BAcomps = triu(t(Bblockmat) %*% Ablockmat, k = 1)
      }
      
      nB = sapply(blocklist,function(l) length(l$dfB.inds))
      
      nBAcomp = (BAcomps * (ABcomps > 0))@x
      
      overlaps = which(BAcomps * ABcomps>0,arr.ind = T)
      nBs = nB[overlaps[,1]]
      
      todo = matrix(overlaps[nBs==nBAcomp,] ,ncol = 2)
      nBs = nBs[nBs==nBAcomp]
      
      if(nrow(todo)>0){
        do = sapply(todo[,1], function(i) length(blocklist[[i]][[2]]) > length(blocklist[[i]][[1]]))
        todo = matrix(todo[do,],ncol = 2)
        nBs = nBs[do]
      }
      
      if(nrow(todo) > 0){
        cat('Pass', pass,'Step 2: Reducing blocks with dfB fully contained in another block\'s dfA \n')
        blocklist = block_reduce(todo,nBs,blocklist,2,crossdim = T)
        #if(!checkres(blocklist,keys)) return(blocklist)
        recalc = T
      }else{
        recalc = F
      }
    }
    
    #Step 4: Resolve overlapping blocks
    if(mode == 'same_dim'){
      if(recalc){
        Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
          cbind(blocklist[[i]]$dfA.inds,i)
        }))
        
        Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
        
        Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
          cbind(blocklist[[i]]$dfB.inds,i)
        }))
        
        Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
      }
      
      Acomps = triu(t(Ablockmat) %*% Ablockmat, k = 1)
      
      Bcomps = triu(t(Bblockmat) %*% Bblockmat, k = 1)
      
      overlaps = which(Bcomps * Acomps > 0, arr.ind = T)
      
      if(nrow(overlaps)>0){
        cat('Pass', pass,'Step 3: Resolving partial overlap between blocks \n')
        
        blocklist = de_overlap_blocks(blocklist,overlaps,A_olaps = Acomps[overlaps],B_olaps = Bcomps[overlaps])
        
        #if(!checkres(blocklist,keys)) return(blocklist)
        recalc = T
      }else{
        recalc = F
      }
    }
    
    if(mode == 'diff_dim'){
      #Step 3, remove redundant record pairs that got scrambled A-B if dedupe is true
      if(dedupe){
        
        if(recalc){
          Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfA.inds,i)
          }))
          
          Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
          
          Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfB.inds,i)
          }))
          
          Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
        }
        
        ABcomps = triu(t(Ablockmat) %*% Bblockmat, k = 1)
        BAcomps = triu(t(Bblockmat) %*% Ablockmat, k = 1)
        
        overlaps = which(ABcomps * BAcomps > 0, arr.ind = T)
        
        if(nrow(overlaps)>0){
          cat('Pass', pass,'Step 4: Resolving redundancies where indices are switched between dfA and dfB \n')
          
          blocklist = de_overlap_blocks(blocklist,overlaps,ABcomps[overlaps],BAcomps[overlaps],cross_ind = T)
          #if(!checkres(blocklist,keys)) return(blocklist)
          recalc = T
        }else{
          recalc = F
        }
      }
    }
    
    #Check for residual overlap
    if(recalc){
      Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
        cbind(blocklist[[i]]$dfA.inds,i)
      }))
      
      Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2], x = 1)
      
      Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
        cbind(blocklist[[i]]$dfB.inds,i)
      }))
      
      Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2], x = 1)
    }
    
    if(mode == 'same_dim'){
      olap1 = triu(t(Ablockmat) %*% Ablockmat,k = 1) 
      olap2 = triu(t(Bblockmat) %*% Bblockmat,k = 1) 
      olap = olap1 * olap2
      rm(olap1, olap2); gc()
      if(length(olap@x) == 0){
        if(!dedupe){
          swtch = F
        }else{
          mode = 'diff_dim' 
        } 
      }
      cat('After pass ',pass,', ',sum(olap),' overlapping record pairs remain \n',sep = '')
    }
    
    if(mode == 'diff_dim'){
      olap1 = triu(t(Ablockmat) %*% Bblockmat,k = 1) 
      olap2 = triu(t(Bblockmat) %*% Ablockmat,k = 1) 
      olap_swtch = olap1 * olap2
      rm(olap1, olap2); gc()
      
      if(length(olap_swtch@x) == 0){
        
        olap1 = triu(t(Ablockmat) %*% Ablockmat,k = 1) 
        olap2 = triu(t(Bblockmat) %*% Bblockmat,k = 1) 
        olap = olap1 * olap2
        rm(olap1, olap2); gc()
        
        if(length(olap@x) > 0){
          mode = 'same_dim'
        }else swtch = F
      } 
      
      cat('After pass ',pass,', ',sum(olap_swtch),' redundant record pairs remain \n',sep = '')
    }
    
    recalc = F
    
    pass = pass + 1
  }
  
  if(dedupe){
    cat('Pass', pass,'Step 4a: separating blocks along diagonal from rectangular blocks \n')
    blocklist = resolve_diags(blocklist)
    #checkres(blocklist,keys)
    recalc = T
  }
  
  ident = sapply(blocklist, function(l) identical(l$dfA.inds,l$dfB.inds))
  nrp1 = sapply(blocklist, function(l) length(l$dfA.inds) == 1 &&  length(l$dfB.inds) == 1)
  drop = ident & nrp1 
  if(dedupe) blocklist = blocklist[!drop]
  #checkres(blocklist,keys)
  return(blocklist)
}

full_block_merge = function(inds,keys,blocklist,dim, AtoB = F){
  
  merge_list = tapply(inds,keys,function(x)x)
  
  duped = duplicated(lapply(merge_list,sort))
  if(AtoB){
    bigtab = t(do.call(cbind,merge_list))
    duped = unlist(c(F, sapply(seq_along(merge_list)[-1],
                       function(i) any(bigtab[i,] %in% bigtab[1:(i-1),]))))
  } 
    
  merge_list = merge_list[!duped]
  keepblocks = sapply(merge_list,'[',1)
  dropped_blocks = sapply(merge_list,function(v) v[-1])
  
  #kbs = intersect(keepblocks, unlist(dropped_blocks))
  
  if(length(merge_list) == 0) return(blocklist)
  
  pb = txtProgressBar(0,length(merge_list))
  
  i = 1
  
  for(v in merge_list){
    templist = vector(mode = 'list', length = 2)
    names(templist) = c('dfA.inds','dfB.inds')
    templist[[dim]] = blocklist[[v[1]]][[dim]]
    templist[[3-dim]] = blocklist[[v[1]]][[3-dim]]
    templist[[3-dim]] = union(templist[[3-dim]],
                              do.call(c,lapply(v[-1], function(i) blocklist[[i]][[3 - dim + (-AtoB)^dim]])))
    
    blocklist[[v[1]]] = templist
    
    setTxtProgressBar(pb,i)
    i = i+1
  }
  close(pb)
  dumpinds = na.omit(unlist(lapply(merge_list, '[',-1)))
  #dumpinds = na.omit(setdiff(dumpinds, keepblocks))
  if(length(dumpinds) > 0) blocklist = blocklist[-dumpinds]
  #checkres(blocklist,keys1)
  blocklist
}

block_reduce = function(domat, nmat, blocklist, dim, crossdim = F){
  
  if(!crossdim){
    donors = apply(nmat, 1, which.min)
    domat = cbind(domat[cbind(1:nrow(domat), donors)],domat[cbind(1:nrow(domat), 3-donors)])
  }
  
  # do = sapply(domat[,1], function(i) length(blocklist[[i]][[dim]]) > length(blocklist[[i]][[3-dim]]))
  # if(sum(do)==0) return(blocklist)
  # domat = matrix(domat[do,],ncol = 2)
  # 
  pb = txtProgressBar(0,nrow(domat))
  del_list = NULL
  
  for(i in 1:nrow(domat)){
    donor = blocklist[[domat[i,1]]]; receiver = blocklist[[domat[i,2]]]
    
    olap = intersect(donor[[3-dim]], receiver[[3-dim + (-crossdim)^dim]])
    resid = setdiff(donor[[3-dim]], olap)
    
    if(length(resid) == 0){del_list = c(del_list, domat[i,1]); next}
    
    blocklist[[domat[i,1]]][[3-dim]] = resid
    setTxtProgressBar(pb, i)
  }
  close(pb)
  del_list = unique(del_list)
  if(!is.null(del_list)) blocklist = blocklist[-del_list]
  #checkres(blocklist,keys1)
  blocklist
}

resolve_diags = function(blocklist){
  resolve_diag = sapply(blocklist,function(l){
    !identical(l$dfA.inds,l$dfB.inds) & any(l$dfA.inds %in% l$dfB.inds)
  })
  if(sum(resolve_diag) > 0){
    templists = vector(mode = 'list', length = sum(resolve_diag))
    
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
  blocklist
}

de_overlap_blocks = function(blocklist,olap_mat,A_olaps, B_olaps,cross_ind = F, nround = 25){
  
  ssize = sapply(blocklist, function(l) length(unlist(l)))
  allA = matrix(sapply(olap_mat,function(i) as.numeric(length(blocklist[[i]]$dfA.inds))),ncol = 2)
  allB = matrix(sapply(olap_mat,function(i) as.numeric(length(blocklist[[i]]$dfB.inds))),ncol = 2)
  if(cross_ind){
    temp = cbind(allA[,1],allB[,2]); allB = cbind(allB[,1],allA[,2]); allA = temp
  }
  
  nrp = allA * allB
  #the donor will be the block for which resulting number of indices is smallest
  #min(dfB + olA  vs dfA + olB)
  
  ssizes = new_ssizes = olap_mat
  ssizes[] = ssize[ssizes]
  new_ssizes = matrix(pmin(allA + B_olaps, allB + A_olaps), ncol = 2)
  
  donors = apply(nrp,1,which.min)
  
  olap_mat = cbind(olap_mat[cbind(1:nrow(olap_mat),donors)],olap_mat[cbind(1:nrow(olap_mat),3-donors)])
  
  donors = setdiff(olap_mat[,1],olap_mat[,2])
  delta = 1
  tempmat = matrix(olap_mat[!(olap_mat[,1] %in% donors),],ncol = 2)
  nround = nround - 1
  while(delta > 0 & nround > 0){
    donors = c(donors, setdiff(tempmat[,1],tempmat[,2]))
    delta = length(setdiff(tempmat[,1],tempmat[,2]))
    tempmat = matrix(olap_mat[!(olap_mat[,1] %in% donors),],ncol = 2)
    nround = nround - 1
  } #order such that we can erase donor after finishing its processing
  
  i = 1
  pb = txtProgressBar(0,length(donors))
  for(d in donors){
    templist = blocklist[[d]]
    
    todo = olap_mat[olap_mat[,1] == d,2]
    
    olaps_a = lapply(todo, function(i){
      l2 = blocklist[[i]]
      if(cross_ind){
        intersect(templist$dfA.inds,l2$dfB.inds)
      }else{
        intersect(templist$dfA.inds,l2$dfA.inds)
      } 
    })
    
    olaps_b = lapply(todo, function(i){
      l2 = blocklist[[i]]
      if(cross_ind){
        intersect(templist$dfB.inds,l2$dfA.inds)
      }else{
        intersect(templist$dfB.inds,l2$dfB.inds)
      } 
    })
    
    templist = list(templist)
    
    for(ii in 1:length(todo)){
      do = sapply(templist, function(l) any(olaps_a[[ii]] %in% l$dfA.inds) && any(olaps_b[[ii]] %in% l$dfB.inds))
      if(!any(do)) next
      
      rem_A = lapply(templist[do], function(l) setdiff(l$dfA.inds, olaps_a[[ii]]))
      ol_A = lapply(templist[do], function(l) intersect(l$dfA.inds, olaps_a[[ii]]))
      rem_B = lapply(templist[do], function(l) setdiff(l$dfB.inds, olaps_b[[ii]]))
      ol_B = lapply(templist[do], function(l) intersect(l$dfB.inds, olaps_b[[ii]]))
      
      templist[do] = lapply(1:sum(do), function(i){
        l = templist[do][[i]]
        
        indicator = (length(rem_A[[i]])>0) + 2*(length(rem_B[[i]])>0)
        
        if(indicator == 0) return(NULL)
        if(indicator == 1) return(list(list(
          dfA.inds = rem_A[[i]], dfB.inds = l$dfB.inds
        )))
        if(indicator == 2) return(list(list(
          dfA.inds = l$dfA.inds, dfB.inds = rem_B[[i]]
        )))
        if(indicator == 3){
          cost1 = length(ol_A[[i]]) + length(rem_B[[i]]) + length(rem_A[[i]]) + length(l$dfB.inds)
          
          cost2 = length(rem_A[[i]]) + length(ol_B[[i]]) + length(l$dfA.inds) + length(rem_B[[i]])
          
          if(cost1 < cost2){
            return(list(
              list(
                dfA.inds = ol_A[[i]],dfB.inds = rem_B[[i]]
              ),list(
                dfA.inds = rem_A[[i]], dfB.inds = l$dfB.inds
              )))
          }else{
            return(list(
              list(
                dfA.inds = rem_A[[i]],dfB.inds = ol_B[[i]]
              ),list(
                dfA.inds = l$dfA.inds, dfB.inds = rem_B[[i]]
              )))
          }
        }})
      templist = c(unlist(templist[do], recursive = F),templist[!do])
      templist = templist[sapply(templist,length)>0]
    }
    blocklist = c(blocklist, templist)
    blocklist[[d]] = list()
    setTxtProgressBar(pb,i)
    i = i+1
  }
  close(pb)
  
  blocklist = blocklist[-donors]
  blocklist
}

checkres = function(blocklist, keys){
  
  ident = sapply(blocklist,function(l) identical(l$dfA.inds,l$dfB.inds))
  
  compares2 = lapply(blocklist, function(l) expand.grid.jc(l[[1]],l[[2]]))
  compares2[ident] = lapply(blocklist[ident], function(l) RcppAlgos::comboGeneral(l[[1]],2))
  compares2 = do.call(rbind,compares2)
  
  lefts = compares2[,1]
  for(i in 1:length(lefts)){
    lefts[i] = which.min(compares2[i,])
  }
  compares2 = cbind(compares2[cbind(1:nrow(compares2),lefts)],compares2[cbind(1:nrow(compares2),3-lefts)])
  compares2 = compares2[compares2[,1] != compares2[,2],]
  
  compares2 = as_tibble(compares2)
  compares2 = compares2 %>% group_by_all() %>% summarize(count = length(V1))
  
  keys2 = log(2) * compares2$V1 + log(3) * compares2$V2
  check = all(keys %in% keys2)
  if(!check) cat("SOME INDICES MISSING AFTER PREVIOUS STEP!!!!!!!! \n")
  return(check)
}


thin_blocks = function(blocklist,antiblocklist,dedupe = F, dims){
  
  ab = unique(antiblocklist) #antiblocklist[order(sapply(antiblocklist,length))]
  
  #for(ab in antiblocklist){
    
    Ainds_anti = do.call(rbind,lapply(1:length(ab),function(i){
      cbind(ab[[i]]$dfA.inds,i)
    }))
    Ablockmat_anti = sparseMatrix(i = Ainds_anti[,1], j = Ainds_anti[,2], x=1, dims = c(dims[1],length(ab)))
    
    #Ablockmat_anti = sparseMatrix(i = Ainds_anti[,1], j = Ainds_anti[,2])
    
    Binds_anti = do.call(rbind,lapply(1:length(ab),function(i){
      cbind(ab[[i]]$dfB.inds,i)
    }))
    Bblockmat_anti = sparseMatrix(i = Binds_anti[,1], j = Binds_anti[,2], x=1, dims = c(dims[2],length(ab)))
    
    #step 1: delete parts of blocks that are fully contained
    swtch = T
    recalc = T
    while(swtch){
      swtch = F
      if(recalc){
        Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
          cbind(blocklist[[i]]$dfA.inds,i)
        }))
        Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
          cbind(blocklist[[i]]$dfB.inds,i)
        }))
        Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2],x=1, dims = c(dims[1],length(blocklist)))
        Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2],x=1, dims = c(dims[2],length(blocklist)))
        recalc = F
      }
      
      overs = t(Ablockmat) %*% Ablockmat_anti * (t(Bblockmat) %*% Bblockmat_anti)
      overlaps = matrix(which(overs>0, arr.ind = T),ncol = 2)
      overs = overs@x
      
      do = unique(overlaps[,1])
      
      
      if(length(do)>0){
        doblock = sapply(do, function(i){
          inds = overlaps[,1] == i
          overlaps[inds,2][which.max(overs[inds])]
        })
        
        newblocks = do.call(c,lapply(seq_along(do), function(i){
          i_ = do[i]
          b = blocklist[[i_]]; u = ab[[doblock[i]]]
          
          Acomm = intersect(b$dfA.inds, u$dfA.inds)
          Adiff = setdiff(b$dfA.inds, u$dfA.inds)
          Bdiff = setdiff(b$dfB.inds, u$dfB.inds)
          
          incl = c(length(Bdiff) > 0, length(Adiff) > 0)
            
          list(list(dfA.inds = Acomm, dfB.inds = Bdiff), list(dfA.inds = Adiff, dfB.inds = b$dfB.inds))[incl]
        }))
        blocklist = c(blocklist[-do], newblocks)
        recalc = T
        swtch = T
      }
    }
    #if dedupe, dim A to B
    if(dedupe){
      swtch = T
      while(swtch){
        swtch = F
        if(recalc){
          Ainds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfA.inds,i)
          }))
          Binds = do.call(rbind,lapply(1:length(blocklist),function(i){
            cbind(blocklist[[i]]$dfB.inds,i)
          }))
          Ablockmat = sparseMatrix(i = Ainds[,1], j = Ainds[,2],x=1, dims = c(dims[1],length(blocklist)))
          Bblockmat = sparseMatrix(i = Binds[,1], j = Binds[,2],x=1, dims = c(dims[2],length(blocklist)))
          recalc = F
        }
        
        overs = t(Ablockmat) %*% Bblockmat_anti * (t(Bblockmat) %*% Ablockmat_anti)
        overlaps = matrix(which(overs>0, arr.ind = T),ncol = 2)
        overs = overs@x
        
        do = unique(overlaps[,1])

        if(length(do)>0){
          doblock = sapply(do, function(i){
            inds = overlaps[,1] == i
            overlaps[inds,2][which.max(overs[inds])]
          })
          
          newblocks = do.call(c,lapply(seq_along(do), function(i){
            i_ = do[i]
            b = blocklist[[i_]]; u = ab[[doblock[i]]]
            
            Acomm = intersect(b$dfA.inds, u$dfB.inds)
            Adiff = setdiff(b$dfA.inds, u$dfB.inds)
            Bdiff = setdiff(b$dfB.inds, u$dfA.inds)
            
            incl = c(length(Bdiff) > 0, length(Adiff) > 0)
            
            list(list(dfA.inds = Acomm, dfB.inds = Bdiff), list(dfA.inds = Adiff, dfB.inds = b$dfB.inds))[incl]
          }))
          blocklist = c(blocklist[-do], newblocks)
          recalc = T
          swtch = T
        }
      }
      #and finally, subset out duplicate indices as separate blocks
      duped_inds = sapply(blocklist, function(x) any(x[[1]] %in% x[[2]]) & !identical(x[[1]], x[[2]]))
      if(any(duped_inds)){
        newblocks = do.call(c,lapply(blocklist[duped_inds], function(x){
          comm = intersect(x[[1]],x[[2]])
          diff1 = setdiff(x[[1]],x[[2]])
          diff2 = setdiff(x[[2]],x[[1]])
          
          incl = c(T, length(diff1)>0, length(diff2)>0)
          list(list(dfA.inds = comm, dfB.inds = comm),
               list(dfA.inds = diff1, dfB.inds = x[[2]]),
               list(dfA.inds = x[[1]], dfB.inds = diff2))[incl]
        }))
        blocklist = c(blocklist[!duped_inds], newblocks)
      }
    }
  #}
  blocklist
}
  