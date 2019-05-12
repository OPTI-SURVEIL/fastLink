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

