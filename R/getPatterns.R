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
      hits = rep(0,length(matchesA))
      if(length(gammalist[[i]]$matches2) > 0){
        mtab1 = do.call(rbind,
                        lapply(1:length(gammalist[[i]]$matches2), function(j) expand.grid.jc(gammalist[[i]]$matches2[[j]][[1]],j)))
        mtab1 = matrix(mtab1[mtab1[,1] %in% matchesA,],ncol = 2)
        
        mtab2 = do.call(rbind,
                        lapply(1:length(gammalist[[i]]$matches2), function(j) expand.grid.jc(gammalist[[i]]$matches2[[j]][[2]],j)))
        mtab2 = matrix(mtab2[mtab2[,1] %in% matchesB,],ncol = 2)
        
        hits = sapply(1:length(matchesA), function(i) any(mtab1[mtab1[,1] == matchesA[i],2] %in% mtab2[mtab2[,1] == matchesB[i],2]))
      }
      
      gammalist2[[i]] = 2 * hits
      gammalist2[[i]][matchesA %in% gammalist[[i]]$nas[[1]] | matchesB %in% gammalist[[i]]$nas[[2]]] = NA
      
      if(partial.match[i]){
        doinds = integer(0)
        if(length(gammalist[[i]]$matches1) > 0){
          pmtab1 = do.call(rbind,
                           lapply(1:length(gammalist[[i]]$matches1), function(j) expand.grid.jc(gammalist[[i]]$matches1[[j]][[1]],j)))
          pmtab1 = pmtab1[pmtab1[,1] %in% matchesA,]
          
          pmtab2 = do.call(rbind,
                           lapply(1:length(gammalist[[i]]$matches1), function(j) expand.grid.jc(gammalist[[i]]$matches1[[j]][[2]],j)))
          pmtab2 = pmtab2[pmtab2[,1] %in% matchesA,]
          
          doinds = sapply(1:length(matchesA), function(i) any(pmtab1[pmtab1[,1] == matchesA[i],2] %in% pmtab2[pmtab2[,1] == matchesB[i],2]))
        }
        
        gammalist2[[i]][doinds] = 1
      }
        
      namevec[i] <- paste0("gamma.", varnames[i])
        
    }
    gammalist2 <- data.frame(do.call(cbind, gammalist2))
    names(gammalist2) <- namevec

    return(gammalist2)

}

