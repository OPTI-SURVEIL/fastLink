#' fastLink
#'
#' Run the fastLink algorithm to probabilistically match
#' two datasets.
#'

#' @usage fastLink(dfA, dfB, varnames, stringdist.match, stringdist.method,
#'   numeric.match, partial.match, cut.a, cut.p, jw.weight, cut.a.num,
#'   cut.p.num, priors.obj, w.lambda, w.pi, address.field, gender.field,
#'   estimate.only, em.obj, dedupe.matches, linprog.dedupe, reweight.names,
#'   firstname.field, cond.indep, n.cores, tol.em, threshold.match, return.all,
#'   return.df, verbose)
#'
#' @param dfA Dataset A - to be matched to Dataset B
#' @param dfB Dataset B - to be matched to Dataset A
#' @param blocklist An optional list of blocks of record pairs produced by blockData. If included, agreement patterns are counted within each block,
#' then aggregated together prior to EM.
#' @param varnames A vector of variable names to use for matching. Must be
#'   present in both dfA and dfB
#' @param stringdist.match A vector of variable names indicating which variables
#'   should use string distance matching. Must be a subset of 'varnames' and
#'   must not be present in 'numeric.match'.
#' @param string.transform An optional function returning a transformation or
#'   list of transformations of the input strings
#'   
#' @param string.transform.args A list of arguments to the provided string transformation function
#' 
#' @param stringdist.method String distance method for calculating similarity,
#'   options are: "jw" Jaro-Winkler (Default), "jaro" Jaro, and "lv" Edit. May
#'   also be a function for custom methods. The custom function should return a N1 by N2 matrix of similarity scores based on 
#'   the strings to be compared, given in the first two arguments, as well as any additional arguments, which must be passed 
#'   as a list to stringdist.args. If internal linkage and deduplication is desired, the custom function should be return
#'   a lower triangular matrix (N1 by N1) of scores for permutations of a vector of strings when the second argument is missing.
#'   
#' @param stringdist.args List of arguments to be provided to stringdist.method
#'   if a custom function is used.
#'   
#' @param numeric.match A vector of variable names indicating which variables
#'   should use numeric matching. Must be a subset of 'varnames' and must not be
#'   present in 'stringdist.match'.
#' @param partial.match A vector of variable names indicating whether to include
#'   a partial matching category for the string distances. Must be a subset of
#'   'varnames' and 'stringdist.match'.
#' @param cut.a Lower bound for full string-distance match, ranging between 0
#'   and 1. Default is 0.94
#' @param cut.p Lower bound for partial string-distance match, ranging between 0
#'   and 1. Default is 0.88
#' @param jw.weight Parameter that describes the importance of the first
#'   characters of a string (only needed if stringdist.method = "jw"). Default
#'   is .10
#' @param cut.a.num Lower bound for full numeric match. Default is 1
#' @param cut.p.num Lower bound for partial numeric match. Default is 2.5
#' @param priors.obj A list containing priors for auxiliary movers information,
#'   as output from calcMoversPriors(). Default is NULL
#' @param w.lambda How much weight to give the prior on lambda versus the data.
#'   Must range between 0 (no weight on prior) and 1 (weight fully on prior).
#'   Default is NULL (no prior information provided).
#' @param w.pi How much weight to give the prior on pi versus the data. Must
#'   range between 0 (no weight on prior) and 1 (weight fully on prior). Default
#'   is NULL (no prior information provided).
#' @param address.field The name of the address field. To be used when
#'   'pi.prior' is included in 'priors.obj'. Default is NULL (no matching
#'   variables should have address prior applied). Must be present in
#'   'varnames'.
#' @param gender.field The name of the field indicating gender. If provided, the
#'   exact-matching gender prior is used in the EM algorithm. Default is NULL
#'   (do not implement exact matching on gender). Must be present in 'varnames'.
#' @param estimate.only Whether to stop running the algorithm after the EM step
#'   (omitting getting the matched indices of dataset A and dataset B). Only the
#'   EM object will be returned. Can be used when running the match on a random
#'   sample and applying to a larger dataset, or for out-of-sample prediction of
#'   matches. Default is FALSE.
#' @param counts.only Whether to stop running the algorithm after counting agreement patterns. 
#' Will return the frequency table used in the EM aglorithm, as well as the list of matching indices for each field, 
#' which may be used to retrieve matches at later points. Default is FALSE.
#' @param em.obj An EM object from a prior run of 'fastLink' or 'emlinkMARmov'.
#'   Parameter estimates will be applied to the matching patterns in 'dfA' and
#'   'dfB'. If provided. 'estimate.only' is set to FALSE. Often provided when
#'   parameters have been estimated on a smaller sample, and the user wants to
#'   apply them to the full dataset. Default is NULL (EM will be estimated from
#'   matching patterns in 'dfA' and 'dfB').
#' @param dedupe.matches Whether to dedupe the set of matches returned by the
#'   algorithm. Default is TRUE.
#' @param linprog.dedupe If deduping matches, whether to use Winkler's linear
#'   programming solution to dedupe. Default is FALSE.
#' @param reweight.names Whether to reweight the posterior match probabilities
#'   by the frequency of individual first names. Default is FALSE.
#' @param firstname.field The name of the field indicating first name. Must be
#'   provided if reweight.names = TRUE.
#' @param cond.indep Estimates for the parameters of interest are obtained from
#'   the Fellegi-Sunter model under conditional independence. Default is TRUE.
#'   If set to FALSE parameters estimates are obtained from a model that allows
#'   for dependencies across linkage fields.
#' @param n.cores Number of cores to parallelize over. Default is NULL.
#' @param tol.em Convergence tolerance for the EM Algorithm. Default is 1e-04.
#' @param threshold.match A number between 0 and 1 indicating either the lower
#'   bound (if only one number provided) or the range of certainty that the user
#'   wants to declare a match. For instance, threshold.match = .85 will return
#'   all pairs with posterior probability greater than .85 as matches, while
#'   threshold.match = c(.85, .95) will return all pairs with posterior
#'   probability between .85 and .95 as matches.
#' @param auto.threshold Boolean. If set to \code{TRUE}, the default, the threshold for accepting matches is automatically calculated 
#' to return the closest number of matches to the overall estimated proportion. Will be set to FALSE if threshold.match is provided
#' or return.all is TRUE
#' @param return.all Whether to return the most likely match for each
#'   observation in dfA and dfB. Overrides user setting of
#'   \code{threshold.match} by setting \code{threshold.match} to 0.0001, and
#'   automatically dedupes all matches. Default is TRUE.
#' @param return.df Whether to return the entire dataframe of dfA and dfB
#'   instead of just the indices. Default is FALSE.
#' @param verbose Whether to print elapsed time for each step. Default is FALSE.
#'
#' @return \code{fastLink} returns a list of class 'fastLink' containing the
#'   following components if calculating matches: \item{matches}{An nmatches X 2
#'   matrix containing the indices of the successful matches in \code{dfA} in
#'   the first column, and the indices of the corresponding successful matches
#'   in \code{dfB} in the second column.} \item{EM}{A list with the output of
#'   the EM algorithm, which contains the exact matching patterns and the
#'   associated posterior probabilities of a match for each matching pattern.}
#'   \item{patterns}{A matrix with the observed matching patterns for each
#'   successfully matched pair.} \item{nobs.a}{The number of observations in
#'   dataset A.} \item{nobs.b}{The number of observations in dataset B.}
#'   \item{zeta.name}{If reweighting by name, the posterior probability of a
#'   match for each match in dataset A and B.}
#'
#'   If only running the EM and not returning the matched indices,
#'   \code{fastLink} only returns the EM object.
#'
#' @author Ted Enamorado <ted.enamorado@gmail.com>, Ben Fifield
#'   <benfifield@gmail.com>, and Kosuke Imai
#'
#' @examples
#' \dontrun{
#' fl.out <- fastLink(dfA, dfB,
#' varnames = c("firstname", "lastname", "streetname", "birthyear"),
#' n.cores = 1)
#' }
#' @export
fastLink <- function(dfA, dfB, blocklist = NULL, varnames,
                     stringdist.match = NULL, 
                     string.transform = NULL,
                     string.transform.args = NULL,
                     stringdist.method = "jw",
                     stringdist.args = NULL,
                     numeric.match = NULL, 
                     partial.match = NULL,
                     cut.a = 0.94, cut.p = 0.88,
                     jw.weight = .10,
                     cut.a.num = 1, cut.p.num = 2.5,
                     priors.obj = NULL,
                     w.lambda = NULL, w.pi = NULL, address.field = NULL,
                     gender.field = NULL, estimate.only = FALSE, counts.only = FALSE, em.obj = NULL,
                     dedupe.matches = TRUE, linprog.dedupe = FALSE,
                     reweight.names = FALSE, firstname.field = NULL, cond.indep = TRUE,
                     n.cores = NULL, tol.em = 1e-04, threshold.match = NULL,auto.threshold = T,
                     return.all = FALSE, return.df = FALSE, verbose = FALSE){

    cat("\n")
    cat(c(paste(rep("=", 20), sep = "", collapse = ""), "\n"))
    cat("fastLink(): Fast Probabilistic Record Linkage\n")
    cat(c(paste(rep("=", 20), sep = "", collapse = ""), "\n\n"))

    if(!is.null(blocklist)){
      blocklist = lapply(blocklist,function(b){
        b$dfA = dfA[b$dfA.inds,]
        b$dfB = dfB[b$dfB.inds,]
        b
      })
    }
    ## --------------------------------------
    ## Process inputs and stop if not correct
    ## --------------------------------------
    
    if(any(class(dfA) %in% c("tbl_df", "data.table"))){
        dfA <- as.data.frame(dfA)
    }
    if(any(class(dfB) %in% c("tbl_df", "data.table"))){
        dfB <- as.data.frame(dfB)
    }
    if(any(!(varnames %in% names(dfA)))){
        stop("Some variables in varnames are not present in dfA.")
    }
    if(any(!(varnames %in% names(dfB)))){
        stop("Some variables in varnames are not present in dfB.")
    }
    if(any(!(stringdist.match %in% varnames))){
        stop("You have provided a variable name for stringdist.match that is not in 'varnames'.")
    }
    if(any(!(numeric.match %in% varnames))){
        stop("You have provided a variable name for numeric.match that is not in 'varnames'.")
    }
    if(length(intersect(numeric.match, stringdist.match)) > 0){
        stop("There is a variable present in both 'numeric.match' and 'stringdist.match'. Please select only one matching metric for each variable.")
    }
    if(is.null(numeric.match)) {
      if (any(!(partial.match %in% varnames)) | any(!(partial.match %in% 
                                                      stringdist.match))) {
        stop("You have provided a variable name for 'partial.match' that is not present in either 'varnames', 'numeric.match', or 'stringdist.match'.")
      }
    } else {
      if (any(!(partial.match %in% varnames)) | any(!(partial.match %in% unique(c(stringdist.match, numeric.match))))) {
        stop("You have provided a variable name for 'partial.match' that is not present in either 'varnames', 'numeric.match', or 'stringdist.match'.")
      }
    }    
    if(!is.null(address.field)){
        if(length(address.field) > 1 | length(gender.field) > 1){
            stop("'address.field' must have at most only one variable name.")
        }
        if(!(address.field %in% varnames)){
            stop("You have provided a variable name for 'address.field' that is not in 'varnames'.")
        }
    }
    if(!is.null(gender.field)){
        if(length(gender.field) > 1){
            stop("'gender.field' must have at most one variable name.")
        }
        if(!(gender.field %in% varnames)){
            stop("You have provided a variable name for 'gender.field' that is not in 'varnames'.")
        }
    }
    if(reweight.names == TRUE & is.null(firstname.field)){
        stop("If reweighting the match probability by first name, you must provide the name of the field representing first name.")
    }
    if(!is.null(firstname.field)){
        if(length(firstname.field) > 1){
            stop("'firstname.field' must have at most one variable name.")
        }
        if(!(firstname.field %in% varnames)){
            stop("You have provided a variable name for 'firstname.field' that is not in 'varnames'.")
        }
    }
    if(!is.null(em.obj)){
        if(!("fastLink.EM" %in% class(em.obj))){
            stop("If providing an EM object, it must be of class 'fastLink.EM'.")
        }
    }
    if(!is.null(em.obj) & estimate.only){
        estimate.only <- FALSE
        cat("You have provided an EM object but have set 'estimate.only' to TRUE. Setting 'estimate.only' to FALSE so that matched indices are returned.\n")
    }
    if(!is.function(stringdist.method)){
      if(!(stringdist.method %in% c("jw", "jaro", "lv")))
        stop("Invalid string distance method. Method should be one of 'jw', 'jaro', or 'lv'.")
      if(stringdist.method == "jw" & !is.null(jw.weight)){
        if(jw.weight < 0 | jw.weight > 0.25){
          stop("Invalid value provided for jw.weight. Remember, jw.weight in [0, 0.25].")
        }
      }
    }
    
    if(is.function(stringdist.method) & is.null(stringdist.args))
      stop("You must provide a list of arguments if using a custom string comparison function")

    if(return.all){
        threshold.match <- 0.001
        if(!dedupe.matches){
            cat("You have specified that all matches be returned but are not deduping the matches. The resulting object may be very large.\n")
        }
    }else{
        cat("If you set return.all to FALSE, you will not be able to calculate a confusion table as a summary statistic.\n")
    }
    if(!is.null(priors.obj) & cond.indep == FALSE){
        cat("The current implementation of fastLink can only incorporate prior information under the conditionally independent model. Ignoring prior information in estimation.")
        priors.obj <- NULL
        w.lambda <- NULL
        w.pi <- NULL
        address.field <- NULL
        gender.field <- NULL
    }

    ## Check class of numeric indicators
    classA <- lapply(dfA[,varnames], class)
    classB <- lapply(dfB[,varnames], class)
    if(any(unlist(classA)[names(classA) %in% numeric.match] != "numeric") |
       any(unlist(classB)[names(classB) %in% numeric.match] != "numeric")){
        stop("You have specified that a variable be compared using numeric matching, but that variable is not of class 'numeric'. Please check your variable classes.")
    }

    ## Check if data frames are identical
    dedupe.df <- FALSE
    if(identical(dfA, dfB)){
        cat("dfA and dfB are identical, assuming deduplication of a single data set.\nSetting return.all to FALSE.\n\n")
        dedupe.matches <- FALSE
        return.all <- FALSE
        dedupe.df <- TRUE
    }

    if(!is.null(threshold.match) || return.all == T) auto.threshold = F
    
    ## Create boolean indicators
    sm.bool <- which(varnames %in% stringdist.match)
    stringdist.match <- rep(FALSE, length(varnames))
    if(length(sm.bool) > 0){
        stringdist.match[sm.bool] <- TRUE
    }

    nm.bool <- which(varnames %in% numeric.match)
    numeric.match <- rep(FALSE, length(varnames))
    if(length(nm.bool) > 0){
        numeric.match[nm.bool] <- TRUE
    }

    pm.bool <- which(varnames %in% partial.match)
    partial.match <- rep(FALSE, length(varnames))
    if(length(pm.bool) > 0){
        partial.match[pm.bool] <- TRUE
    }

    af.bool <- which(varnames %in% address.field)
    address.field <- rep(FALSE, length(varnames))
    if(length(af.bool) > 0){
        address.field[af.bool] <- TRUE
    }

    gf.bool <- which(varnames %in% gender.field)
    gender.field <- rep(FALSE, length(varnames))
    if(length(gf.bool) > 0){
        gender.field[gf.bool] <- TRUE
    }

    fn.bool <- which(varnames %in% firstname.field)
    firstname.field <- rep(FALSE, length(varnames))
    if(length(fn.bool) > 0){
        firstname.field[fn.bool] <- TRUE
    }

    ## ----------------------------
    ## Calculate agreement patterns
    ## ----------------------------
    cat("Calculating matches for each variable.\n")
    start <- Sys.time()
    
    gammalist <- vector(mode = "list", length = length(varnames))
    if(!is.null(blocklist)) gammalist.recover = gammalist
    
    for(i in 1:length(varnames)){
      if(verbose){
        matchtype <- ifelse(stringdist.match[i], "string-distance", ifelse(numeric.match[i], "numeric", "exact"))
        cat("    Matching variable", varnames[i], "using", matchtype, "matching.\n")
      }
      ## Convert to character
      if(is.factor(dfA[,varnames[i]]) | is.factor(dfB[,varnames[i]])){
        dfA[,varnames[i]] <- as.character(dfA[,varnames[i]])
        dfB[,varnames[i]] <- as.character(dfB[,varnames[i]])
      }
      ## Warn if no variation (except for gender blocking)
      if(!gender.field[i]){
        if(sum(is.na(dfA[,varnames[i]])) == nrow(dfA) | length(unique(dfA[,varnames[i]])) == 1){
          cat(paste("WARNING: You have no variation in dataset A for", varnames[i], "or all observations are missing."))
        }
        if(sum(is.na(dfB[,varnames[i]])) == nrow(dfB) | length(unique(dfB[,varnames[i]])) == 1){
          cat(paste("WARNING: You have no variation in dataset B for", varnames[i], "or all observations are missing."))
        }
      }
      if(sum(dfA[,varnames[i]] %in% dfB[,varnames[i]]) == 0){
        cat(paste0("WARNING: You have no exact matches for ", varnames[i], "."))
      }

      ## Get patterns
      if(stringdist.match[i]){
        if(partial.match[i]){
          if(!is.null(blocklist)){
            templists = lapply(1:length(blocklist), function(b){
              if(verbose) cat(paste('Processing variable', varnames[i],'for block', b,'\n'))
              b = blocklist[[b]]
              capture.output(res <- gammaCKpar(
                b$dfA[,varnames[i]], b$dfB[,varnames[i]], cut.a = cut.a, cut.p = cut.p, method = stringdist.method, transform = string.transform,
                transform.args = string.transform.args,method.args = stringdist.args, w = jw.weight, n.cores = n.cores,
                dedupe = dedupe.df))
              res
            })
            gammalist[[i]] = templists
            if(!estimate.only){
              gammalist.recover[[i]] = lapply(1:length(gammalist[[i]]), function(j){
                glist = gammalist[[i]][[j]]; b = blocklist[[j]]
                glist$matches2 = lapply(glist$matches2, function(l){
                  l[[1]] = b$dfA.inds[l[[1]]]; l[[2]] = b$dfB.inds[l[[2]]]
                  l
                })
                glist$matches1 = lapply(glist$matches1, function(l){
                  l[[1]] = b$dfA.inds[l[[1]]]; l[[2]] = b$dfB.inds[l[[2]]]
                  l
                })
                glist$nas = lapply(1:2, function(i){
                  b[[i]][glist$nas[[i]]]
                })
                glist
              })
              gammalist.recover[[i]]$matches2 = do.call('c',lapply(gammalist.recover[[i]],'[[','matches2'))
              gammalist.recover[[i]]$matches1 = do.call('c',lapply(gammalist.recover[[i]],'[[','matches1'))
              gammalist.recover[[i]]$nas = list(unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',1))),
                                                unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',2))))
              gammalist.recover[[i]]$.identical = unlist(lapply(gammalist.recover[[i]],'[[','.identical'))
            }
          }else{
            gammalist[[i]] <- gammaCKpar(
              dfA[,varnames[i]], dfB[,varnames[i]], cut.a = cut.a, cut.p = cut.p, method = stringdist.method, transform = string.transform,
              transform.args = string.transform.args,method.args = stringdist.args, w = jw.weight, n.cores = n.cores,
              dedupe = dedupe.df)
          }
          
        }else{
          if(!is.null(blocklist)){
            templists = lapply(1:length(blocklist), function(b){
              if(verbose) cat(paste('Processing variable', varnames[i],'for block', b,'\n'))
              b = blocklist[[b]]
              capture.output(res<-gammaCK2par(
                b$dfA[,varnames[i]], b$dfB[,varnames[i]], cut.a = cut.a, method = stringdist.method, transform = string.transform,
                transform.args = string.transform.args,method.args = stringdist.args, w = jw.weight, n.cores = n.cores,
                dedupe = dedupe.df))
              res
            })
            gammalist[[i]] = templists
            if(!estimate.only){
              gammalist.recover[[i]] = lapply(1:length(gammalist[[i]]), function(j){
                glist = gammalist[[i]][[j]]; b = blocklist[[j]]
                glist$matches2 = lapply(glist$matches2, function(l){
                  l[[1]] = b$dfA.inds[l[[1]]]; l[[2]] = b$dfB.inds[l[[2]]]
                  l
                })
                glist$nas = lapply(1:2, function(i){
                  b[[i]][glist$nas[[i]]]
                })
                glist
              })
              gammalist.recover[[i]]$matches2 = do.call('c',lapply(gammalist.recover[[i]],'[[','matches2'))
              gammalist.recover[[i]]$nas = list(unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',1))),
                                                unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',2))))
              gammalist.recover[[i]]$.identical = unlist(lapply(gammalist.recover[[i]],'[[','.identical'))
            }
          }else{
            gammalist[[i]] <- gammaCK2par(dfA[,varnames[i]], dfB[,varnames[i]], cut.a = cut.a, 
                                          method = stringdist.method, transform = string.transform,
                                          transform.args = string.transform.args,
                                          method.args = stringdist.args,
                                          w = jw.weight, n.cores = n.cores,
                                          dedupe = dedupe.df)
          }
        }
      }else if(numeric.match[i]){
        if(partial.match[i]){
          if(!is.null(blocklist)){
            templists = lapply(1:length(blocklist), function(b){
              if(verbose) cat(paste('Processing variable', varnames[i],'for block', b,'\n'))
              b = blocklist[[b]]
              capture.output(res<-gammaNUMCKpar(
                b$dfA[,varnames[i]], b$dfB[,varnames[i]], cut.a = cut.a.num, cut.p = cut.p.num, n.cores = n.cores,
                dedupe = dedupe.df))
              res
            })
            gammalist[[i]] = templists
            if(!estimate.only){
              gammalist.recover[[i]] = lapply(1:length(gammalist[[i]]), function(j){
                glist = gammalist[[i]][[j]]; b = blocklist[[j]]
                glist$matches2 = lapply(glist$matches2, function(l){
                     l[[1]] = b$dfA.inds[l[[1]]]; l[[2]] = b$dfB.inds[l[[2]]]
                     l
                   })
                glist$matches1 = lapply(glist$matches1, function(l){
                     l[[1]] = b$dfA.inds[l[[1]]]; l[[2]] = b$dfB.inds[l[[2]]]
                     l
                  })
                  glist$nas = lapply(1:2, function(i){
                     b[[i]][glist$nas[[i]]]
                  })
                  glist
              })
              gammalist.recover[[i]]$matches2 = do.call('c',lapply(gammalist.recover[[i]],'[[','matches2'))
              gammalist.recover[[i]]$matches1 = do.call('c',lapply(gammalist.recover[[i]],'[[','matches1'))
              gammalist.recover[[i]]$nas = list(unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',1))),
                                                unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',2))))
              gammalist.recover[[i]]$.identical = unlist(lapply(gammalist.recover[[i]],'[[','.identical'))
            }
            
          }else{
            gammalist[[i]] <- gammaNUMCKpar(
              dfA[,varnames[i]], dfB[,varnames[i]], cut.a = cut.a.num, cut.p = cut.p.num, n.cores = n.cores,
              dedupe = dedupe.df
            )
          }
        }else{
          if(!is.null(blocklist)){
            templists = lapply(1:length(blocklist), function(b){
              if(verbose) cat(paste('Processing variable', varnames[i],'for block', b,'\n'))
              b = blocklist[[b]]
              capture.output(res<-gammaNUMCK2par(
                b$dfA[,varnames[i]], b$dfB[,varnames[i]], cut.a = cut.a.num, n.cores = n.cores,
                dedupe = dedupe.df))
              res
            })
            gammalist[[i]] = templists
            if(!estimate.only){
              gammalist.recover[[i]] = lapply(1:length(gammalist[[i]]), function(j){
                glist = gammalist[[i]][[j]]; b = blocklist[[j]]
                glist$matches2 = lapply(glist$matches2, function(l){
                  l[[1]] = b$dfA.inds[l[[1]]]; l[[2]] = b$dfB.inds[l[[2]]]
                  l
                })
                glist$nas = lapply(1:2, function(i){
                  b[[i]][glist$nas[[i]]]
                })
                glist
              })
              gammalist.recover[[i]]$matches2 = do.call('c',lapply(gammalist.recover[[i]],'[[','matches2'))
              gammalist.recover[[i]]$nas = list(unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',1))),
                                                unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',2))))
              gammalist.recover[[i]]$.identical = unlist(lapply(gammalist.recover[[i]],'[[','.identical'))
            }
          }else{
            gammalist[[i]] <- gammaNUMCK2par(
              dfA[,varnames[i]], dfB[,varnames[i]], cut.a = cut.a.num, n.cores = n.cores,
              dedupe = dedupe.df
            )
          }
        }
      }else{
        if(!is.null(blocklist)){
          templists = lapply(1:length(blocklist), function(b){
            if(verbose) cat(paste('Processing variable', varnames[i],'for block', b,'\n'))
            b = blocklist[[b]]
            capture.output(res<-gammaKpar(
              b$dfA[,varnames[i]], b$dfB[,varnames[i]], gender = gender.field[i], n.cores = n.cores,
              dedupe = dedupe.df))
            res
          })
          gammalist[[i]] = templists
          if(!estimate.only){
            gammalist.recover[[i]] = lapply(1:length(gammalist[[i]]), function(j){
              glist = gammalist[[i]][[j]]
              b = blocklist[[j]]
              glist$matches2 = lapply(glist$matches2, function(l){
                l[[1]] = b$dfA.inds[l[[1]]]; l[[2]] = b$dfB.inds[l[[2]]]
                l
              })
              glist$nas = lapply(1:2, function(i){
                b[[i]][glist$nas[[i]]]
              })
              glist
            })
            gammalist.recover[[i]]$matches2 = do.call('c',lapply(gammalist.recover[[i]],'[[','matches2'))
            gammalist.recover[[i]]$nas = list(unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',1))),
                                              unique(unlist(lapply(lapply(gammalist.recover[[i]],'[[','nas'),'[',2))))
            gammalist.recover[[i]]$.identical = unlist(lapply(gammalist.recover[[i]],'[[','.identical'))
          }
        }else{
          gammalist[[i]] <- gammaKpar(dfA[,varnames[i]], dfB[,varnames[i]], gender = gender.field[i], n.cores = n.cores,
                                      dedupe = dedupe.df)
        }
      }
    }
      
    end <- Sys.time()
    if(verbose){
        cat("Calculating matches for each variable took", round(difftime(end, start, units = "mins"), 2), "minutes.\n\n")
    }

    ## Get row numbers
    nr_a <- nrow(dfA)
    nr_b <- nrow(dfB)

    ## ------------------------------
    ## Get counts for zeta parameters
    ## ------------------------------
    cat("Getting counts for parameter estimation.\n")
    start <- Sys.time()
    if(!is.null(blocklist)){
      counts = lapply(1:length(blocklist),function(i){
        if(verbose) cat('Counting patterns for block ', i)
        glist = lapply(gammalist,'[[',i); nobs.a = length(blocklist[[i]]$dfA.inds); nobs.b = length(blocklist[[i]]$dfB.inds)
        capture.output(res <- tableCounts(glist, nobs.a, nobs.b, n.cores =n.cores,# min(c(nobs.a * nobs.b %/% 1e6 + 1, parallel::detectCores()-1)), 
                                          dedupe = dedupe.df))
        res
      })
      counts = tibble::as_tibble(do.call(rbind, counts))
      counts = dplyr::group_by_at(counts,dplyr::vars(-counts))
      counts = dplyr::summarize(counts,counts = sum(counts))
      counts = as.matrix(counts)
    }else{
      counts <- tableCounts(gammalist, nobs.a = nr_a, nobs.b = nr_b, n.cores = n.cores,
                            dedupe = dedupe.df)
    }
    
    colnames(counts)[seq_along(varnames)] = paste0('gamma.',varnames)
  
    end <- Sys.time()
    if(verbose){
        cat("Getting counts for parameter estimation took", round(difftime(end, start, units = "mins"), 2), "minutes.\n\n")
    }
    if(counts.only){
      out = list(patterns = counts,
                 gammalist = gammalist)
      return(out)
    }
    
    ## ------------------------------
    ## Run or impute the EM algorithm
    ## ------------------------------
    if(is.null(em.obj)){
        ## Run EM algorithm
        cat("Running the EM algorithm.\n")
        start <- Sys.time()
        if(is.null(priors.obj)){
            lambda.prior <- NULL
            pi.prior <- NULL
        }else{
            if("lambda.prior" %in% names(priors.obj)){
                lambda.prior <- priors.obj$lambda.prior
            }
            if("pi.prior" %in% names(priors.obj)){
                if(!("lambda.prior" %in% names(priors.obj))){
                    stop("Must specify a prior for lambda if providing a prior for pi.")
                }
                pi.prior <- priors.obj$pi.prior
            }else{
                pi.prior <- NULL
            }
        }
        if(cond.indep == FALSE){
            resultsEM <- emlinklog(patterns = counts, nobs.a = nr_a, nobs.b = nr_b,
                                   tol = tol.em, varnames = varnames)  
        }else{
            resultsEM <- emlinkMARmov(patterns = counts, nobs.a = nr_a, nobs.b = nr_b,
                                      tol = tol.em,
                                      prior.lambda = lambda.prior, w.lambda = w.lambda,
                                      prior.pi = pi.prior, w.pi = w.pi,
                                      address.field = address.field, 
                                      gender.field = gender.field,
                                      varnames = varnames)
        }
        end <- Sys.time()
        if(verbose){
            cat("Running the EM algorithm took", round(difftime(end, start, units = "secs"), 2), "seconds.\n\n")
        }
    }else{
        cat("Imputing matching probabilities using provided EM object.\n")
        resultsEM <- emlinkRS(counts, em.obj, nr_a, nr_b)
    }
    if(auto.threshold){
      trgt_count = resultsEM$p.m * sum(resultsEM$patterns.w[,'counts'])
      ordered_counts = resultsEM$patterns.w[order(resultsEM$zeta.j,decreasing = T),'counts']
      thresh_ind = which.min(abs(cumsum(ordered_counts) - trgt_count))
      threshold.match = sort(resultsEM$zeta.j,decreasing = T)[thresh_ind]
      if(threshold.match < 0.5 & thresh_ind == 1){
        cat('No patterns have match probabilities above 0.5. Setting threshold to 0.85 and returning no matches. \n')
        threshold.match = 0.85
      }else{
        cat('Selected match probability threshold is: ', threshold.match,'\n')
      }
     }
    resultsEM$threshold.match = threshold.match
    if(max(resultsEM$zeta.j) < threshold.match) {
        warning(paste0("No matches found for the threshold value used. We recommend trying a lower threshold.match value. Note that you currently have threshold.match set to ", threshold.match, "."))
    }

    ## -----------------------------------------------
    ## Get the estimated matches, dedupe, and reweight
    ## -----------------------------------------------
    if(!estimate.only){
        
        ## Get matches
        cat("Getting the indices of estimated matches.\n")
        start <- Sys.time()
        if(!is.null(blocklist)){
          matches <- matchesLink(gammalist.recover, nobs.a = nr_a, nobs.b = nr_b,
                                 em = resultsEM, thresh = threshold.match,
                                 n.cores = n.cores, dedupe = dedupe.df)
        }else{
          matches <- matchesLink(gammalist, nobs.a = nr_a, nobs.b = nr_b,
                                 em = resultsEM, thresh = threshold.match,
                                 n.cores = n.cores, dedupe = dedupe.df)
        } 
        
        
        end <- Sys.time()
        if(verbose){
            cat("Getting the indices of estimated matches took", round(difftime(end, start, units = "mins"), 2), "minutes.\n\n")
        }

        ## Get the patterns
        if(!is.null(blocklist)){
          patterns <- getPatterns(matchesA = matches$inds.a, matchesB = matches$inds.b,
                                  varnames = varnames, partial.match = partial.match,
                                  gammalist = gammalist.recover)
        }else{
          patterns <- getPatterns(matchesA = matches$inds.a, matchesB = matches$inds.b,
                                  varnames = varnames, partial.match = partial.match,
                                  gammalist = gammalist)
        } 
        
        ## Run deduplication
        if(dedupe.matches & length(matches$inds.a) > 0){
            cat("Deduping the estimated matches.\n")
            start <- Sys.time()
            ddm.out <- dedupeMatches(matchesA = dfA[matches$inds.a,], matchesB = dfB[matches$inds.b,],
                                     EM = resultsEM, matchesLink = matches, patterns = patterns,
                                     linprog = linprog.dedupe)
            matches <- ddm.out$matchesLink
            resultsEM <- ddm.out$EM
            end <- Sys.time()
            if(verbose){
                cat("Deduping the estimated matches took", round(difftime(end, start, units = "mins"), 2), "minutes.\n\n")
            }
        }else if(length(matches$inds.a) > 0){
            cat("Calculating the posterior for each pair of matched observations.\n")
            start <- Sys.time()
            zeta <- getPosterior(dfA[matches$inds.a,], dfB[matches$inds.b,], EM = resultsEM,
                                 patterns = patterns)
            end <- Sys.time()
            if(verbose){
                cat("Calculating the posterior for each matched pair took", round(difftime(end, start, units = "mins"), 2), "minutes.\n\n")
            }
        }

        ## Get the patterns
        cat("Getting the match patterns for each estimated match.\n")
        start <- Sys.time()
        if(!is.null(blocklist)){
          patterns <- getPatterns(matchesA = matches$inds.a, matchesB = matches$inds.b,
                                  varnames = varnames, partial.match = partial.match,
                                  gammalist = gammalist.recover)
        }else{
          patterns <- getPatterns(matchesA = matches$inds.a, matchesB = matches$inds.b,
                                  varnames = varnames, partial.match = partial.match,
                                  gammalist = gammalist)
        } 
        end <- Sys.time()
        if(verbose){
            cat("Getting the match patterns for each estimated match took", round(difftime(end, start, units = "mins"), 2), "minutes.\n\n")
        }

        ## Reweight first names or get zeta
        if(reweight.names & length(matches$inds.a) > 0){
            cat("Reweighting match probabilities by frequency of occurrence.\n")
            start <- Sys.time()
            rwn.out <- nameReweight(dfA, dfB, EM = resultsEM, gammalist = gammalist, matchesLink = matches,
                                    varnames = varnames, firstname.field = firstname.field,
                                    patterns = patterns, threshold.match = threshold.match, n.cores = n.cores)
            end <- Sys.time()
            if(verbose){
                cat("Reweighting by first name took", round(difftime(end, start, units = "mins"), 2), "minutes.\n\n")
            }
        }

        ## Return object
        out <- list()
        if(return.df){
            out[["dfA.match"]] <- dfA[matches$inds.a,]
            out[["dfB.match"]] <- dfB[matches$inds.b,]
        }
        out[["matches"]] <- matches
        out[["EM"]] <- resultsEM
        out[["patterns"]] <- patterns
        if(dedupe.matches & length(matches$inds.a) > 0){
            out[["posterior"]] <- ddm.out$max.zeta
        }else if(length(matches$inds.a) > 0){
            out[["posterior"]] <- zeta
        }
        if(reweight.names & length(matches$inds.a) > 0){
            out[["posterior"]] <- rwn.out
        }
        out[["nobs.a"]] <- nr_a
        out[["nobs.b"]] <- nr_b
        if(return.all){
            class(out) <- c("fastLink", "confusionTable")
        }else{
            class(out) <- "fastLink"
        }
        if(dedupe.df){
            class(out) <- c(class(out), "fastLink.dedupe")
        }
    }else{
        out <- resultsEM
    }

    return(out)

}

