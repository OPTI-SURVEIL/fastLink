#' gammaKpar
#'
#' Field comparisons: 0 disagreement, 2 total agreement.
#'
#' @usage gammaKpar(matAp, matBp, gender, n.cores)
#' 
#' @param matAp vector storing the comparison field in data set 1
#' @param matBp vector storing the comparison field in data set 2
#' @param gender Whether the matching variable is gender. Will override
#' standard warnings of missingness/nonvariability. Default is FALSE.
#' @param dedupe Logical indicator for whether internal linkage and duplication is taking place
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#'
#' @return \code{gammaKpar} returns a list with the indices corresponding to each
#' matching pattern, which can be fed directly into \code{tableCounts} and \code{matchesLink}.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield <benfifield@gmail.com>, and Kosuke Imai
#'
#' @examples
#' \dontrun{
#' g1 <- gammaKpar(dfA$birthyear, dfB$birthyear)
#' }
#' @export

## ------------------------
## gamma.k.par
## This function applies gamma.k
## in parallel
## ------------------------

gammaKpar <- function(matAp, matBp, gender = FALSE, dedupe = F,n.cores = NULL) {

    ## For visible bindings
    i <- NULL
    varnm = names(matAp)

    if(any(class(matAp) %in% c("tbl_df", "data.table",'data.frame'))){
        matAp <- as.data.frame(matAp)[,1]
    }
    if(any(class(matBp) %in% c("tbl_df", "data.table",'data.frame'))){
        matBp <- as.data.frame(matBp)[,1]
    }

    if(is.null(n.cores)) {
        n.cores <- detectCores() - 1
    }

    if(is.character(matAp) | is.factor(matAp)) matAp[matAp == ""] <- NA
    if(is.character(matBp) | is.factor(matBp)) matBp[matBp == ""] <- NA

    if(!gender){
        if(sum(is.na(matAp)) == length(matAp) | length(unique(matAp)) == 1){
            cat("WARNING: You have no variation in variable '", varnm,"', or all observations are missing in dataset A.\n",sep = '')
        }
        if(sum(is.na(matBp)) == length(matBp) | length(unique(matBp)) == 1){
            cat("WARNING: You have no variation in variable '", varnm, "', or all observations are missing in dataset B.\n",sep = '')
        }
    }else{
        if(sum(is.na(matAp)) == length(matAp)){
            cat("WARNING: You have no variation in variable '", varnm, "', or all observations are missing in dataset A.\n",sep = '')
        }
        if(sum(is.na(matBp)) == length(matBp)){
            cat("WARNING: You have no variation in variable '", varnm, "', or all observations are missing in dataset B.\n",sep = '')
        }
    }

    matrix.1 <- as.matrix(as.character(matAp))
    matrix.2 <- as.matrix(as.character(matBp))

    matrix.1[is.na(matrix.1)] <- "1234MF"
    matrix.2[is.na(matrix.2)] <- "9876ES"

    
    u.values.1 <- unique(matrix.1)
    u.values.2 <- unique(matrix.2)

    matches <- u.values.1[u.values.1 %in% u.values.2]
    
    .identical = rep(T,length(matches))
    
    ht1 <- new.env(hash=TRUE)
    ht2 <- new.env(hash=TRUE)
    matches.l <- as.list(matches)
    out <- list()
    if(length(matches) < 100) n.cores = 1
    if(length(matches.l) > 0){
      if(Sys.info()[['sysname']] == "Windows") {
        if (n.cores == 1) '%oper%' <- foreach::'%do%'
        else { 
          '%oper%' <- foreach::'%dopar%'
          cl <- makeCluster(n.cores)
          registerDoSNOW(cl)
          on.exit(stopCluster(cl))
        }
        
        pb = txtProgressBar(0, length(matches.l),style = 1)
        progress = function(n) setTxtProgressBar(pb, n)
        opts = list(progress = progress)
        
        final.list <- foreach(i = 1:length(matches.l), .options.snow = opts) %oper% {
          ht1 <- which(matrix.1 == matches.l[[i]]); 
          if(dedupe){ht2 = ht1} else{ht2 <- which(matrix.2 == matches.l[[i]])}
          list(ht1, ht2)
        }
        
        close(pb)
        
      } else {
        final.list <- mclapply(matches.l, function(s){
          ht1[[s]] <- which(matrix.1 == s) 
          if(dedupe){ ht2[[s]] <- ht1[[s]]}else{ht2[[s]] = which(matrix.2 == s)}
          list(ht1[[s]], ht2[[s]]) }, mc.cores = getOption("mc.cores", n.cores))
      }
      
      
      out[["matches2"]] <- final.list
    }else{
      #warning(paste0('There are no exact matches for variable ', varnm,', and it will not be considered during linkage \nThis may indicate that there are no true matches in your data, and estimated match probabilities may be unreliable'))
      out[["matches2"]] <- vector(length = 0,mode = 'list') 
    }
    
    na.list <- list()
    na.list[[1]] <- which(matrix.1 == "1234MF")
    na.list[[2]] <- which(matrix.2 == "9876ES")

    out[["nas"]] <- na.list
    out[['.identical']] <- .identical
    class(out) <- c("fastLink", "gammaKpar")

    return(out)
}

## ------------------------
## End of gamma.k.par
## ------------------------
