#' find_dependencies
#'
#' get functions and packages required by a function. Used for parallelization
#' of custom string comparison functions
#'
#' @usage find_dependencies(f)
#'
#' @param f the function to be analyzed
#' @param drops function calls to ignore. Used to limit redundancy during
#'   recursion
#' @param lyr the current level of recursion, used to control warning messages
#'
#' @return \code{find_dependencies} returns a list with members 'depends', a
#'   data table containing packages and functions called by the target, and
#'   'missing', a vector of functions referred to within the call stack that are
#'   not currently defined. Note that missing functions are present in the call
#'   stack, but may not actually be called depending on other function
#'   parameters 
#'

find_dependencies = function(f,drops = NULL,lyr = 1){
  calls = codetools::findGlobals(f)
  
  calls = calls[!(calls %in% drops)]
  drops = c(drops,calls)
  
  fcalls = sapply(calls,function(x) try(get(x),silent = T))
  missing = sapply(fcalls,class) == 'character'
  
  missingcalls = calls[missing]; calls = calls[!missing];fcalls = fcalls[!missing]
  prims = sapply(fcalls,is.primitive)
  
  calls = calls[!prims]
  
  pkgs = sapply(sapply(calls,find),function(c) gsub('package:','',c))
  
  defs = pkgs %in% options()$defaultPackages | pkgs == 'base'
  
  calls = calls[!defs]; pkgs = pkgs[!defs]
  
  res = list(depends = data.table::data.table(calls = calls, pkgs = pkgs,stringsAsFactors = F),
             missing = missingcalls)
  
  if(length(calls) == 0) return(res)
  
  reslyr = lapply(res$depends$calls,function(c) find_dependencies(get(c),drops = drops,lyr = lyr+1))
  
  depends = data.table::rbindlist(lapply(c(list(res),reslyr), '[[','depends'))
  depends = depends[!duplicated(depends$calls)]
  missing = unique(do.call('c',lapply(c(list(res),reslyr), '[[','missing')))
  
  if(lyr == 1 && length(missing)>0) warning(paste0('Missing calls in function tree: ',paste0(missing,collapse=', ')))
  
  list(depends = depends, missing = missing)
}
