#utility functions

#fast expand.grid
expand.grid.jc <- function(seq1,seq2) {
  cbind(Var1 = rep.int(seq1, length(seq2)), 
        Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2))))
}

#get indices from  self-comparison inclusive combination

combo_deindexer = function(inds,N){
  a = -0.5; b_u = N + 1.5
  
  i = floor((-b_u + sqrt(b_u^2 - 4 * a * -(N+inds)))/(2*a))
  
  j = inds + 0.5*i^2 - (N+1.5)*i + N +i
  cbind(i,j)
}

#scale isotonic regression by prior mixing proportion of matches to nonmatches
isoreg_scale = function(isoreg, p.m){
  ys = 1 / ((1 / isoreg$y - 1) / ((1-isoreg$p.m)/isoreg$p.m) * ((1-p.m)/p.m) + 1)
  isoreg$y = ys
  isoreg$p.m = p.m
  isoreg
}

#fit isotonic regression to similarity scores
fit.isored <- function(isored, x0){
  x = isored$x
  y = isored$y
  ind = cut(x0, breaks = x, labels = FALSE, include.lowest = TRUE)
  min.x <- min(x)
  max.x <- max(x)
  ind[x0 > max.x] = length(x)
  ind[x0 < min.x] = 1
  
  slope = c(diff(y)/diff(x),0)
  val = y[ind] + pmax(x0 - x[ind],0) * slope[ind]
  val
}



#unpack blocked gammalist
unpack_blocked_gammalist = function(gammalist, blocklist){
  gammalist.recover = gammalist

  for(i in 1:length(blocklist)){
    indsa = blocklist[[i]]$dfA.inds; indsb = blocklist[[i]]$dfB.inds
    for(j in 1:length(gammalist.recover)){
      gammalist.recover[[j]][[i]]$matches2 = 
        lapply(gammalist.recover[[j]][[i]]$matches2, function(l) list(indsa[l[[1]]],indsb[l[[2]]]))
      
      if(!is.null(gammalist.recover[[j]][[i]]$matches1)){
        gammalist.recover[[j]][[i]]$matches1 = 
          lapply(gammalist.recover[[j]][[i]]$matches1, function(l) list(indsa[l[[1]]],indsb[l[[2]]]))
      }
      gammalist.recover[[j]][[i]]$nas = list(indsa[gammalist.recover[[j]][[i]]$nas[[1]]],
                                             indsb[gammalist.recover[[j]][[i]]$nas[[2]]])
    }
  }
  for(j in 1:length(gammalist.recover)){
    gammalist.recover[[j]]$matches2 = do.call(c,lapply(gammalist.recover[[j]], '[[', 'matches2'))
    if(!is.null(gammalist.recover[[j]][[i]]$matches1)){
      gammalist.recover[[j]]$matches1 = do.call(c,lapply(gammalist.recover[[j]], '[[', 'matches1'))
    }
    gammalist.recover[[j]]$nas = 
      list(unique(do.call(c,lapply(gammalist.recover[[j]], function(l) l$nas[[1]]))),
           unique(do.call(c,lapply(gammalist.recover[[j]], function(l) l$nas[[2]]))))
    gammalist.recover[[j]] = 
      gammalist.recover[[j]][names(gammalist.recover[[j]])[nchar(names(gammalist.recover[[j]]))>0]]
  }
  gammalist.recover
}
