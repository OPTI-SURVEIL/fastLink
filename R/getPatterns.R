#' getPatterns
#'
#' Get the full matching patterns for all matched pairs in dataset A and dataset B
#'
#' @param matchesA Indices of matches from dataframe A.
#' @param matchesB Indices of matches from dataframe B.
#' @param varnames A vector of variable names to use for matching.
#' @param partial.match A vector of booleans, indicating whether to include
#' a partial matching category for each column. Must be same length
#' as varnames. Default is FALSE for all variables.
#' @param gammalist A list of agreement pointers returned by gammaCKpar or similar functions for each column
#'
#' @return \code{getPatterns()} returns a dataframe with a row for each matched pair,
#' where each column indicates the matching pattern for each matching variable.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com> and Ben Fifield <benfifield@gmail.com>
#' @export
getPatterns <- function(matchesA, matchesB, varnames,
                        partial.match, gammalist){
    
  
    ## --------------
    ## Start function
    ## --------------
    
    ## ----------
    ## Get gammas
    ## ----------
    gammalist2 <- vector(mode = "list", length = length(varnames))
    namevec <- rep(NA, length(varnames))
    
    for(i in 1:length(gammalist)){
      mtab1 = do.call(rbind,
                lapply(1:length(gammalist[[i]]$matches2), function(j) expand.grid.jc(gammalist[[i]]$matches2[[j]][[1]],j)))
      mtab2 = do.call(rbind,
                lapply(1:length(gammalist[[i]]$matches2), function(j) expand.grid.jc(gammalist[[i]]$matches2[[j]][[2]],j)))
      
      i1s = mtab1[match(matchesA,mtab1[,1]),2]; i2s = mtab2[match(matchesB,mtab2[,1]),2]
      
      gammalist2[[i]] = (i1s == i2s) * 2
      gammalist2[[i]][matchesA %in% gammalist[[i]]$nas[[1]] | matchesB %in% gammalist[[i]]$nas[[2]]] = NA
      
      if(partial.match[i]){
        pmtab1 = do.call(rbind,
                  lapply(1:length(gammalist[[i]]$matches1), function(j) expand.grid.jc(gammalist[[i]]$matches1[[j]][[1]],j)))
        pmtab2 = do.call(rbind,
                  lapply(1:length(gammalist[[i]]$matches1), function(j) expand.grid.jc(gammalist[[i]]$matches1[[j]][[2]],j)))
        
        i1s = pmtab1[match(matchesA,pmtab1[,1]),2]; i2s = pmtab2[match(matchesB,pmtab2[,1]),2]
        
        doinds = which(i1s == i2s)
        gammalist2[[i]][doinds] = 1
      }
        
      namevec[i] <- paste0("gamma.", varnames[i])
        
    }
    gammalist2 <- data.frame(do.call(cbind, gammalist2))
    names(gammalist2) <- namevec

    return(gammalist2)

}

