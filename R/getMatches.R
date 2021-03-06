globalVariables(c('V2', 'V3', '.'))

#' getMatches
#'
#' Subset two data frames to the matches returned by \code{fastLink()}
#' or \code{matchesLink()}. Can also return a single deduped data frame
#' if dfA and dfB are identical and fl.out is of class 'fastLink.dedupe'. 
#' The deduped data frame contains an extra column, \code{prior.indices}, 
#' which contains NA values for all rows not identified as duplicates, and otherwise 
#' lists the row indices of all prior rows linked to the current row. 
#'
#' @usage getMatches(dfA, dfB, fl.out, threshold.match, combine.dfs)
#' @param dfA Dataset A - matched to Dataset B by \code{fastLink()}.
#' @param dfB Dataset B - matches to Dataset A by \code{fastLink()}.
#' @param fl.out Either the output from \code{fastLink()} or \code{matchesLink()}.
#' @param threshold.match A number between 0 and 1 indicating the lower bound that the
#' user wants to declare a match. For instance, threshold.match = .85 will return all pairs with posterior probability greater than .85 as matches.
#' Default is 0.85.
#' @param combine.dfs Whether to combine the two data frames being merged into a single data frame. If FALSE, two data frames are returned in a list. Does nothing in the case of deduping a single dataframe. Default is TRUE.
#' @param inspect.matches Boolean. Whether to  return a single dataframe with each matched pair presented in subsequent rows, with spaces separating the pairs. Overrides \code{combine.dfs}
#'
#' @return \code{getMatches()} returns a list of two data frames:
#' \item{dfA.match}{A subset of \code{dfA} subsetted down to the successful matches.}
#' \item{dfB.match}{A subset of \code{dfB} subsetted down to the successful matches.}
#'
#' @author Ben Fifield  <benfifield@gmail.com>
#'
#' @examples
#' \dontrun{
#' fl.out <- fastLink(dfA, dfB,
#' varnames = c("firstname", "lastname", "streetname", "birthyear"),
#' n.cores = 1)
#' ret <- getMatches(dfA, dfB, fl.out)
#' }
#' @export
getMatches <- function(dfA, dfB, fl.out, threshold.match = 0.85, combine.dfs = TRUE, inspect.matches = FALSE){

    ## Convert data frames
    if(any(class(dfA) %in% c("tbl_df", "data.table"))){
        dfA <- as.data.frame(dfA)
    }
    if(any(class(dfB) %in% c("tbl_df", "data.table"))){
        dfB <- as.data.frame(dfB)
    }

    if(inherits(fl.out, "fastLink.dedupe") & !identical(dfA, dfB)){
        stop("You have provided a fastLink object from deduping a single data frame, but dfA and dfB are not identical. Please check your inputs.")
    }
    if(identical(dfA, dfB) & !inherits(fl.out, "fastLink.dedupe")){
        stop("dfA and dfB are identical, but fl.out is not of class 'fastLink.dedupe.' Please check your inputs.")
    }

    ## Depending on class
    if(inherits(fl.out, "matchesLink")){
        dfA.match <- dfA[fl.out$inds.a,]
        dfB.match <- dfB[fl.out$inds.b,]
        if(combine.dfs){
            names.dfB <- names(dfB.match)[!(names(dfB.match) %in% names(dfA.match))]
            if(length(names.dfB) > 0){
                df.match <- cbind(dfA.match, dfB.match[,names.dfB])
            }else{
                df.match <- dfA.match
            }
            out <- df.match
        }else{
            out <- list(dfA.match = dfA.match, dfB.match = dfB.match)
        }
    }else if(inherits(fl.out,'fastLink') && inspect.matches){
      ord = order(fl.out$posterior,decreasing = T)
      varnames = fl.out$EM$varnames
      dfA.match <- dfA[fl.out$matches$inds.a[ord],varnames]
      dfB.match <- dfB[fl.out$matches$inds.b[ord],varnames]
      
      combinedframe = matrix('',nrow = nrow(dfA.match) * 4, ncol = ncol(dfA.match) + 1 + 'posterior' %in% names(fl.out))
      namevec = c('row.index',names(dfA.match),'p_match')
      colnames(combinedframe) = namevec[1:ncol(combinedframe)]
      
      prefixes = c('dfA','dfB')
      if(inherits(fl.out, "fastLink.dedupe")) prefixes[2] = 'dfA'
      
      combinedframe[seq(1,nrow(combinedframe)-3,4),1:(ncol(dfA.match)+1)] = as.matrix(cbind(paste0(prefixes[1],'.',fl.out$matches$inds.a[ord]),dfA.match))
      combinedframe[seq(2,nrow(combinedframe)-2,4),1:(ncol(dfA.match)+1)] = as.matrix(cbind(paste0(prefixes[2],'.',fl.out$matches$inds.b[ord]),dfB.match))
      
      patmat = as.matrix(cbind('agreement pattern:',fl.out$patterns[ord,]))
      if('posterior' %in% names(fl.out)){
        patmat = cbind(patmat, round(fl.out$posterior[ord],4))
      }
      
      combinedframe[seq(3,nrow(combinedframe)-1,4),] = patmat
      out <- as.data.frame(combinedframe)
    }else if(inherits(fl.out, "fastLink.dedupe")){
        #idea - mark each duplicated row (i.e. with larger index), with the list of earlier row indices linked to it
        ## Get ID
        id_tmp <- data.frame(id = 1:nrow(dfA))

        ## Get matches
        matches <- data.table(
            cbind(
                fl.out$matches$inds.b[fl.out$posterior >= threshold.match],
                fl.out$matches$inds.a[fl.out$posterior >= threshold.match]
            )
        )

        pasteT <- function(x) {
            x = sort(x) # 1,2,3 is the same as 3,2,1
            x = paste(x, collapse = ",") 
            x
        }

        ## Dedupe
        matches[, V3 := pasteT(V2), by = "V1"]    
        ans <- matches[, .(id_2 = unique(V3)), by = "V1"]
        #ans$id_2 <- as.numeric(as.factor(ans$id_2))
        colnames(ans) <- c("id", "id_2")

        ## Merge and output new df
        out_df <- merge(id_tmp, ans, by = "id",all = T)
        dfA$prior.indices <- out_df$id_2
        out <- dfA
        
        
    }else{
      dfA.match <- dfA[fl.out$matches$inds.a,]
      dfB.match <- dfB[fl.out$matches$inds.b,]
        if(combine.dfs){
            names.dfB <- names(dfB.match)[!(names(dfB.match) %in% names(dfA.match))]
            if(length(names.dfB) > 0){
                df.match <- cbind(dfA.match, dfB.match[,names.dfB])
            }else{
                df.match <- dfA.match
            }
            df.match <- cbind(df.match, fl.out$patterns)
            if("posterior" %in% names(fl.out)){
                df.match$posterior <- fl.out$posterior
                df.match <- df.match[df.match$posterior >= threshold.match,]
            }
            out <- df.match
        }else{
          dfA.match <- cbind(dfA.match, fl.out$patterns)
          dfB.match <- cbind(dfB.match, fl.out$patterns)
          if("posterior" %in% names(fl.out)){
            dfA.match$posterior <- fl.out$posterior
            dfB.match$posterior <- fl.out$posterior
            dfA.match <- dfA.match[dfA.match$posterior >= threshold.match,]
            dfB.match <- dfB.match[dfB.match$posterior >= threshold.match,]
          }
          out <- list(dfA.match = dfA.match, dfB.match = dfB.match)
        }
    }

    return(out)
    
}

