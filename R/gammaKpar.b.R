#' gammaKpar.b
#'
#' Field comparisons for blocked data: 0 disagreement, 2 total agreement.
#'
#' @usage gammaKpar.b(matAp, matBp, blocklist, gender, dedupe.blocks, n.cores)
#' 
#' @param matAp vector storing the comparison field in data set 1
#' @param matBp vector storing the comparison field in data set 2
#' @param blocklist list of blocking indices
#' @param gender Whether the matching variable is gender. Will override
#' standard warnings of missingness/nonvariability. Default is FALSE.
#' @param dedupe.blocks Logical indicators for whether internal linkage and duplication is taking place
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#'
#' @return \code{gammaKpar.b} returns a list (over blocks) of lists with the indices corresponding to each
#' matching pattern, which can be fed directly into \code{tableCounts.b} and \code{matchesLink.b}.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield <benfifield@gmail.com>, and Kosuke Imai
#'
#' @export

## ------------------------
## gamma.k.par
## This function applies gamma.k
## in parallel
## ------------------------

gammaKpar.b <- function(matAp, matBp, blocklist, gender = FALSE, dedupe.blocks, n.cores = NULL) {
  
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
  
  pb1 = txtProgressBar(0, length(blocklist),style = 1,char = '+')
  progress1 = function(n) setTxtProgressBar(pb1, n)
  opts1 = list(progress = progress1)
  
  #get unique values coocurring across all blocks
  
  cat('getting matching unique values of variable *', varnm, '* across all blocks \n', sep = '')
  
  cl1 = makeCluster(n.cores)
  registerDoSNOW(cl1)
  clusterExport(cl1, c('matrix.1','matrix.2','blocklist','dedupe.blocks'))
  
  allmatches <- foreach(i = 1:length(blocklist),.options.snow = opts1) %dopar% {
    b = blocklist[[i]]
    m.1 = matrix.1[b[[1]],]; m.2 = matrix.2[b[[2]],]
    matches = intersect(unique(m.1),unique(m.2))
  }
    
  #get matching indices per block
  cat('\ngetting matching indices for variable *', varnm, '* across all blocks \n', sep = '')
  
  out <- foreach(i = 1:length(blocklist), .options.snow = opts1) %dopar% {
    b = blocklist[[i]]
    ddp = dedupe.blocks[i]
    m.1 = matrix.1[b[[1]],]; m.2 = matrix.2[b[[2]],]
    m.b = allmatches[[i]]
    
    if(length(m.b) > 0){
      final.list = lapply(m.b,function(x){
        ht2 = ht1 = which(m.1 == x) 
        if(!ddp) ht2 = which(m.2 == x)
        list(ht1,ht2)
      })
    }else{
      final.list = vector('list',length = 0)
    }
    
    na.list = list(which(m.1 == '1234MF'), which(m.2 == '9876ES'))
    
    .identical = rep(T,length(m.b))
    
    out = list(matches2 = final.list, nas = na.list, .identical = .identical)
    class(out) <- c("fastLink", "gammaKpar")
    out
  }
  
  stopCluster(cl1)
  close(pb1)
  
  return(out)
}

## ------------------------
## End of gamma.k.par.b
## ------------------------
