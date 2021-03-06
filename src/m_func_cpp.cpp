// [[Rcpp::depends(RcppArmadillo,RcppEigen)]]
#define EIGEN_NO_DEBUG

#include <RcppArmadillo.h>
#include "progress.hpp"
#include "eta_progress_bar.hpp"
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

typedef Eigen::Triplet<double> Trip;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<double>::InnerIterator InIt;

arma::mat indexing(const std::vector<arma::vec> s, const int l1, const int l2,
		   const int l3, const int l4, const bool identical, const bool dedupe){
  
  // Get elements, declare matrix object
  arma::vec s0 = s[0]; arma::vec s1 = s[1];
  arma::mat index_out;

  // Subset down to right elements
  arma::uvec s0_l1 = s0 > l1;
  arma::uvec s0_l2 = s0 <= l2;
  arma::uvec s1_l3 = s1 > l3;
  arma::uvec s1_l4 = s1 <= l4;
  
  arma::uvec s0_bool = s0_l1 % s0_l2;
  arma::uvec s1_bool = s1_l3 % s1_l4;
  
  if(sum(s0_bool) >= 1){

    if(sum(s1_bool) >= 1){
      // Subset vectors to indices that fall in range
      arma::vec temp0 = s0.elem(find(s0_bool == true)); 
      arma::vec temp1 = s1.elem(find(s1_bool == true));
      
      if(dedupe){
        // Rcout << "dedupe is true" << std::endl;
        if(identical){
          // Rcout << "identical is true" << std::endl;
        //maximum row index has to be less than maximum column index
        arma::uvec s0_bool2 = temp0 < temp1.max();
        //minimum column index has to be greater than minimum row index
        arma::uvec s1_bool2 = temp1 > temp0.min();
        if(sum(s0_bool2) >= 1)
          if(sum(s1_bool2) >= 1){
            arma::vec temp0_ = temp0.elem(find(s0_bool2 == true));
            arma::vec temp1_ = temp1.elem(find(s1_bool2 == true));
            
            unsigned i; unsigned j; // Expand grid, declare size of matrix
            int jlow = 0;
            
            int rowcount = 0;
            int nrow = 0;
            
            for(i = 0; i < temp0_.n_elem; i++){
              nrow += sum(temp1_ > temp0_(i));
            }
            index_out.set_size(nrow, 2);
            
            for(i = 0; i < temp0_.n_elem; i++){
              while(temp1_(jlow) <= temp0_(i))
                jlow++;
            for(j = jlow; j < temp1_.n_elem; j ++){
                index_out(rowcount,0) = temp0_[i];
                index_out(rowcount,1) = temp1_[j];
                rowcount++;
              }
            }
          }
        }else{
          // Rcout << "identical is false" << std::endl;
          unsigned i; unsigned j;// Expand grid, declare size of matrix
          int rowcount = 0;
          
          index_out.set_size(temp0.n_elem * temp1.n_elem, 2);
          for(i = 0; i < temp0.n_elem; i++)
            for(j = 0; j < temp1.n_elem; j++){
              index_out(rowcount,0) = std::min(temp0[i],temp1[j]);
              index_out(rowcount,1) = std::max(temp0[i],temp1[j]);
              rowcount++;
            }
            
        }
      }else{
        // Rcout << "dedupe is false" << std::endl;
        unsigned i; unsigned j;// Expand grid, declare size of matrix
        int rowcount = 0;
        
        index_out.set_size(temp0.n_elem * temp1.n_elem, 2);
        for(i = 0; i < temp0.n_elem; i++)
          for(j = 0; j < temp1.n_elem; j++){
            index_out(rowcount,0) = temp0[i];
            index_out(rowcount,1) = temp1[j];
            rowcount++;
          }
        }
      if(!index_out.is_empty()){
        // Rcout << "Original lists:" << std::endl;
        // Rcout << s0 << std::endl;
        
        index_out.col(0) -= l1;
        index_out.col(1) -= l3;//converted to index within submatrix
      }
    }
  }
  return index_out;
}

std::vector<arma::vec> indexing_na(const std::vector<arma::vec> s,
				   const int l1, const int l2,
				   const int l3, const int l4,
				   const bool dedupe){

  // Unpack
  const arma::vec s0 = s[0]; const arma::vec s1 = s[1];
  
  // Subset
  arma::uvec s0_l1 = s0 > l1;
  arma::uvec s0_l2 = s0 <= l2;
  arma::uvec s1_l3 = s1 > l3;
  arma::uvec s1_l4 = s1 <= l4;
  
  // if(dedupe && (l1 == l3) && (l2 == l4)){
  //   arma::uvec s0_l2 = s0 <= (l2-1);
  //   arma::uvec s1_l3 = s1 > (l3+1);
  // }
  
  arma::uvec s0_bool = s0_l1 % s0_l2;
  arma::uvec s1_bool = s1_l3 % s1_l4;
  arma::vec temp0 = s0.elem(find(s0_bool == true)) - l1; //converted to index within submatrix
  arma::vec temp1 = s1.elem(find(s1_bool == true)) - l3;

  // Output
  std::vector<arma::vec> out(2);
  out[0] = temp0;
  out[1] = temp1;
  return out;
  
}

std::vector<SpMat> unpack_matches(const std::vector< std::vector<arma::mat> > x,
				  const arma::vec dims, const bool match){

  // Declare objects
  int len_x = x.size(); int i; unsigned j; unsigned k; std::vector<arma::mat> feature_adj;
  int matrix_length; 
  std::vector<SpMat> list_out(len_x); arma::mat adj_store;
  double val; arma::mat feature_adj_j;
  
  // Loop over features
  for(i = 0; i < len_x; i++){

    // Unpack that feature
    feature_adj = x[i]; //list of matrices of matching indices for unique values in the feature
    matrix_length = 0;
    arma::vec is_not_null(feature_adj.size()); //check whether each list is empty

    // Loop over entries in feature, delete if null
    for(j = 0; j < feature_adj.size(); j++){
      feature_adj_j = feature_adj[j];
      if(!feature_adj_j.is_empty()){
	is_not_null[j] = 1;
	matrix_length += feature_adj_j.n_rows; //keep track of total number of matches (by matrix rows)
      } 
    }

    // Get correct entry for the sparse list
    if(match){
      val = std::pow(2.0, 2 + (i * 3));
    }else{
      val = std::pow(2.0, 1 + (i * 3));
    }

    // Create tripletList out of feature_adj
    std::vector<Trip> tripletList;
    tripletList.reserve(matrix_length);
    for(j = 0; j < feature_adj.size(); j++){
      if(is_not_null[j] == 1){
	adj_store = feature_adj[j];
	for(k = 0; k < adj_store.n_rows; k++){
	  tripletList.push_back(Trip(adj_store(k,0)-1, adj_store(k,1)-1, val));
	}
      }
    }

    // Convert to sparse matrix
    SpMat sp(dims(0), dims(1));
    sp.setFromTriplets(tripletList.begin(), tripletList.end());

    // Store in list_out
    list_out[i] = sp;
    
  }

  return list_out;
  
}

arma::vec getNotIn(const arma::vec vec1, const arma::vec vec2){

  unsigned i; unsigned j; bool match; arma::vec store_notin(vec1.n_elem); int counter = 0;
  for(i = 0; i < vec1.n_elem; i++){
    match = false;
    for(j = 0; j < vec2.n_elem; j++){
      if(vec1(i) == vec2(j)){
	match = true;
	break;
      }
    }
    if(!match){
      store_notin(counter) = vec1(i);
      counter++;
    }
  }

  return store_notin.subvec(0, counter-1);
  
}

std::vector<SpMat> create_sparse_na(const std::vector< std::vector<arma::vec> > nas,
				    const arma::vec dims, const arma::vec lowerlims, const bool dedupe){
//lowerlims is the lower index of each dimension

  int i; unsigned j; //unsigned k; 
  double val; const int nobs_a = dims[0]; const int nobs_b = dims[1];
  arma::vec nas_a; arma::vec nas_b; std::vector<SpMat> list_out(nas.size());
  arma::vec nobs_a_notnull_inb; std::vector<arma::vec> nas_extract(2);
  
  // Create comparison vector of indices of nobs_a
  arma::vec nobs_a_vec(nobs_a);
  for(i = 0; i < nobs_a; i++){
    nobs_a_vec[i] = i+1;
  }
  //iterating over 
  int nassize = nas.size();
  for(i = 0; i < nassize; i++){

    // Get exponent value
    val = std::pow(2.0, 3 + (i * 3));

    // Extract indices of NAs
    nas_extract = nas[i];
    nas_a = nas_extract[0];
    nas_b = nas_extract[1];
    // Create triplet
    std::vector<Trip> tripletList;
    tripletList.reserve(nas_a.size() * nobs_b + nas_b.size() * nobs_a);
    
    if(dedupe && (lowerlims(0) == lowerlims(1))){
      //Rcout << "going into diagonal subroutine" << std::endl;
      for(j = 0; j < nas_a.size(); j++){
        for(int k = std::max(0.0,nas_a[j]); k < nobs_b; k++){
          tripletList.push_back(Trip(nas_a[j]-1, k, val));
        }
      }
      
      bool same_vec = approx_equal(nobs_a_vec, nas_a, "absdiff", 0.002);
      
      if(!same_vec){
        arma::uvec temp;
        nobs_a_notnull_inb = getNotIn(nobs_a_vec, nas_a); 
        for(j = 0; j < nas_b.size(); j++){
          temp = nobs_a_notnull_inb < nas_b[j];
          arma::vec noannib = nobs_a_notnull_inb.elem(find(temp == true));
          
          for(unsigned k = 0; k < noannib.size(); k++){
            tripletList.push_back(Trip(noannib[k]-1, nas_b[j]-1, val));
          }
        }
      }
      
    }else{
      // Create triplet
      //Rcout << "going into nondiagonal subroutine" << std::endl;
      tripletList.reserve(nas_a.size() * nobs_b + nas_b.size() * nobs_a);
      for(j = 0; j < nas_a.size(); j++){
        for(int k = 0; k < nobs_b; k++){
          tripletList.push_back(Trip(nas_a[j]-1, k, val));
        }
      }
      
      bool same_vec = approx_equal(nobs_a_vec, nas_a, "absdiff", 0.002);
      
      if(!same_vec){
        nobs_a_notnull_inb = getNotIn(nobs_a_vec, nas_a); 
        for(j = 0; j < nas_b.size(); j++){
          for(unsigned k = 0; k < nobs_a_notnull_inb.size(); k++){
            tripletList.push_back(Trip(nobs_a_notnull_inb[k]-1, nas_b[j]-1, val));
          }
        }
      }
    }
    
    // Convert to sparse matrix
    SpMat sp(dims(0), dims(1));
    sp.setFromTriplets(tripletList.begin(), tripletList.end());
    
    // Store in list.out
    list_out[i] = sp;

  }

  return list_out;

}

std::vector<arma::vec> m_func(const std::vector< std::vector<arma::mat> > matches,
			      const std::vector< std::vector<arma::mat> > pmatches,
			      const std::vector< std::vector<arma::vec> > nas,
			      const arma::vec lims,
			      const arma::vec lims_2,
			      const arma::vec listid,
			      const bool dedupe,
			      const bool matchesLink
			      ){
  
  // Create sparse matches, pmatches object
  const std::vector<SpMat> matches_up  = unpack_matches(matches,  lims, true); //one sparse matrix for matches in each field
  const std::vector<SpMat> pmatches_up = unpack_matches(pmatches, lims, false);
  
  //Rcout << "Completed unpacking successfully in mfunc " << std::endl;
  // Create sparse NA matrix
  
  const std::vector<SpMat> nas_sp = create_sparse_na(nas, lims, lims_2, dedupe); //one sparse matrix for nas in each field
  
  //Rcout << "Completed sparse Nas successfully in mfunc " << std::endl;
  
  // Add up everything
  SpMat sp(lims(0), lims(1));
  SpMat match_pmatch(lims(0), lims(1));
  SpMat match_pmatch_na(lims(0), lims(1)); //creating sparse matrices of approp dimension
  unsigned i;
  for(i = 0; i < matches_up.size(); i++){ //iterating over fields
    match_pmatch = matches_up[i] + pmatches_up[i];
    match_pmatch_na = match_pmatch + nas_sp[i];
    sp = sp + match_pmatch_na;
  } //resultant matrix contains information on matches, partial matches, and nas hashed for each field

  std::vector<arma::vec> nz_out(2); //list of 2 vectors .. number of zeros?
  if(matchesLink == false){

    // Create table by iterating through
    arma::vec nz(sp.nonZeros()); //vector for all nonzero entries, i.e. missing or matched at any level in any field
    int counter = 0;
    for(i = 0; i < sp.outerSize(); i++){ //iterate over columns
      for(InIt it(sp,i); it; ++it){ //iterate over nonzero elements in each column
	nz[counter] = it.value(); //stores all nonzero elements in nz
	counter++;
      }
    }
    
    // Get unique values and create table
    arma::vec nz_unique = unique(nz); //unique hash values of nonzero elements
    arma::vec nz_unique_counts(nz_unique.n_elem);
    arma::uvec nz_unique_i;
    for(i = 0; i < nz_unique_counts.n_elem; i++){
      nz_unique_i = find(nz == nz_unique(i));
      nz_unique_counts(i) = nz_unique_i.n_elem;
    }
    
    int num_zeros = lims(0)*lims(1) - sum(nz_unique_counts);
    
    // Add zeros, create nz_out
    nz_unique.resize(nz_unique.n_elem + 1);
    nz_unique_counts.resize(nz_unique_counts.n_elem + 1);
    nz_unique_counts(nz_unique_counts.n_elem-1) = num_zeros;
    nz_out[0] = nz_unique;
    nz_out[1] = nz_unique_counts;
    
  }else{

    // Find elements of sp that are in listid, replace with 9999
    arma::vec nz_out_rows(sp.nonZeros());
    arma::vec nz_out_cols(sp.nonZeros());
    int counter = 0;
    for(i = 0; i < sp.outerSize(); i++){
      for(InIt it(sp,i); it; ++it){
	if(any(listid == it.value())){
	  nz_out_rows(counter) = it.row() + lims_2(0);
	  nz_out_cols(counter) = it.col() + lims_2(1);
	  counter++;
	}
      }
    }

    // Store in nz_out, only take 0:(counter-1)
    if(counter > 0){
      nz_out[0] = nz_out_rows.subvec(0, counter-1);
      nz_out[1] = nz_out_cols.subvec(0, counter-1);
    }
    
  }

  return nz_out;
  
}

// [[Rcpp::export]]
std::vector< std::vector<arma::vec> > m_func_par(const std::vector< std::vector< std::vector<arma::vec> > > temp,
						 const std::vector< std::vector< std::vector<arma::vec> > > ptemp,
						 const std::vector< std::vector<arma::vec> > natemp,
						 const arma::vec limit1, const arma::vec limit2,
						 const arma::vec nlim1, const arma::vec nlim2,
						 const arma::mat ind,
						 const arma::vec listid,
						 const std::vector< std::vector <bool> > identical,
						 const bool dedupe = false,
						 const bool matchesLink = false,
						 const int threads = 1){

  // Declare objects (shared)
  std::vector< std::vector<arma::vec> > ind_out(ind.n_rows);

  // Declare objects (private)
  int n; int m;
  std::vector< std::vector<arma::vec> > temp_feature;
  std::vector< std::vector<arma::vec> > ptemp_feature;
  std::vector <bool> ident_feature;

  // Declare objects (firstprivate)
  std::vector< std::vector<arma::mat> > templist(temp.size());
  std::vector< std::vector<arma::mat> > ptemplist(ptemp.size());
  std::vector< std::vector<arma::vec> > natemplist(natemp.size());
  std::vector<arma::vec> mf_out(2);
  arma::vec lims(2);
  arma::vec lims_2(2);
  ETAProgressBar pb;
  Progress p(ind.n_rows, true, pb);
  // Declare pragma environment
#ifdef _OPENMP
  omp_set_num_threads(threads);
  int threadsused = omp_get_max_threads();
  Rcout << "    Parallelizing calculation using OpenMP. "
	<< threadsused << " threads out of "
	<< omp_get_num_procs() << " are used."
	<< std::endl;
#pragma omp parallel for private(n, m, temp_feature, ptemp_feature, ident_feature) firstprivate(lims, lims_2, templist, ptemplist, natemplist, mf_out) 
#endif
  

  for(unsigned i = 0; i < ind.n_rows; i++){
      // Get indices of the rows
      if(p.increment()){
        n = ind(i,0)-1; m = ind(i, 1)-1;
        lims(0) = nlim1(n); lims(1) = nlim2(m); //size
        lims_2(0) = limit1(n), lims_2(1) = limit2(m); //start
        
        // Loop over the number of features
        for(unsigned j = 0; j < temp.size(); j++){
          
          // Within this, loop over the list of each feature
          temp_feature = temp[j];
          ptemp_feature = ptemp[j];
          ident_feature = identical[j];
          std::vector<arma::mat> indlist(temp_feature.size());
          std::vector<arma::mat> pindlist(ptemp_feature.size());
          unsigned k;
          for(k = 0; k < temp_feature.size(); k++){
            if(temp_feature.size() > 0){
              indlist[k] = indexing(temp_feature[k], limit1[n], limit1[n+1],
                                    limit2[m], limit2[m+1], ident_feature[k], dedupe); //returns 2 column matrix of matching indices for field
            }
            
          }
          for(k = 0; k < ptemp_feature.size(); k++){
            if(ptemp_feature.size() > 0){
              pindlist[k] = indexing(ptemp_feature[k], limit1[n], limit1[n+1],
                                     limit2[m], limit2[m+1], false, dedupe); //returns 2 column matrix of partially matching indices for field
            }
          }
          templist[j] = indlist;
          ptemplist[j] = pindlist;
          natemplist[j] = indexing_na(natemp[j], limit1[n], limit1[n+1],
                                      limit2[m], limit2[m+1],dedupe); //returns 2 member list of na indices
        }
        
        // Run m_func, initial arguments are a list of lists of 2 column matrices of matches, partial matches,
        // and a 2 member list of missing indices for each field
        mf_out = m_func(templist, ptemplist, natemplist, lims, lims_2, listid, 
                        dedupe, matchesLink);
        ind_out[i] = mf_out;
      }
      
    }
    

  return ind_out;
  
}

// [[Rcpp::export]]
std::vector< std::vector<arma::vec> > m_func_par_b(const std::vector< std::vector< std::vector< std::vector<arma::vec> > > > temp,
                                                 const std::vector< std::vector< std::vector< std::vector<arma::vec> > > > ptemp,
                                                 const std::vector< std::vector< std::vector<arma::vec> > > natemp,
                                                 const std::vector< arma::vec > limit1, const std::vector < arma::vec > limit2,
                                                 const std::vector< arma::vec > nlim1, const std::vector < arma::vec > nlim2,
                                                 const arma::mat ind,
                                                 const arma::vec listid,
                                                 const std::vector< std::vector< std::vector <bool> > > identical,
                                                 const std::vector<bool> dedupe,
                                                 const bool matchesLink = false,
                                                 const int threads = 1){
  
  //use 3rd column in ind to index list positions for blocks
  //Rcout << "started m_func_par_b " << std::endl;
  
  // Declare objects (shared)
  std::vector< std::vector<arma::vec> > ind_out(ind.n_rows);
  //Rcout << "declared shared objects " << std::endl;
  // Declare objects (private)
  int n; int m; int blk;
  std::vector< std::vector<arma::vec> > temp_feature;
  std::vector< std::vector<arma::vec> > ptemp_feature;
  std::vector <bool> ident_feature;
  //Rcout << "declared private objects " << std::endl;
  // Declare objects (firstprivate)
  std::vector< std::vector<arma::mat> > templist((temp[1]).size()); //stores matches for each feature
  std::vector< std::vector<arma::mat> > ptemplist((ptemp[1]).size()); //stores partial matches for each feature
  std::vector< std::vector<arma::vec> > natemplist((natemp[1]).size()); //stores missing indices for each feature
//Rcout << "declared firstprivate objects " << std::endl;
  std::vector<arma::vec> mf_out(2);
  
  arma::vec lims(2);
  arma::vec lims_2(2);
  //Rcout << "declared all firstprivate objects " << std::endl;
  ETAProgressBar pb;
  Progress p(ind.n_rows, true, pb);
  //Rcout << "initialized progress bar " << std::endl;
  // Declare pragma environment
#ifdef _OPENMP
  omp_set_num_threads(threads);
  int threadsused = omp_get_max_threads();
  Rcout << "    Parallelizing calculation using OpenMP. "
        << threadsused << " threads out of "
        << omp_get_num_procs() << " are used."
        << std::endl;
#pragma omp parallel for private(n, m, blk, temp_feature, ptemp_feature, ident_feature) firstprivate(lims, lims_2, templist, ptemplist, natemplist, mf_out) 
#endif
  
  
  for(unsigned i = 0; i < ind.n_rows; i++){
    // Get indices of the rows
    if(p.increment()){
      n = ind(i,0)-1; m = ind(i, 1)-1; blk = ind(i,2)-1;
      lims(0) = (nlim1[blk])(n); lims(1) = (nlim2[blk])(m); //size
      lims_2(0) = (limit1[blk])(n), lims_2(1) = (limit2[blk])(m); //start
      
      //Rcout << "Checkpoint 1 - starting loop over features" << std::endl;
      
      // Loop over the number of features
      for(unsigned j = 0; j < (temp[1]).size(); j++){
        
        // Within this, loop over the list of each feature
        temp_feature = (temp[blk])[j];
        ptemp_feature = (ptemp[blk])[j];
        ident_feature = (identical[blk])[j];
        std::vector<arma::mat> indlist(temp_feature.size());
        std::vector<arma::mat> pindlist(ptemp_feature.size());
//        Rcout << "feature" << std::endl;
//        Rcout << j << std::endl;
        
        
        unsigned k;
        for(k = 0; k < temp_feature.size(); k++){
          if(temp_feature.size() > 0){
            indlist[k] = indexing(temp_feature[k], (limit1[blk])[n], (limit1[blk])[n+1],
                                  (limit2[blk])[m], (limit2[blk])[m+1], ident_feature[k], dedupe[blk]); //returns 2 column matrix of matching indices for field
          }
//          Rcout << "Completed temp for feature" << j << std::endl;
          
        }
        for(k = 0; k < ptemp_feature.size(); k++){
          if(ptemp_feature.size() > 0){
            pindlist[k] = indexing(ptemp_feature[k], (limit1[blk])[n], (limit1[blk])[n+1],
                                   (limit2[blk])[m], (limit2[blk])[m+1], false, dedupe[blk]); //returns 2 column matrix of partially matching indices for field
          }
        }
        templist[j] = indlist;
        ptemplist[j] = pindlist;
        natemplist[j] = indexing_na((natemp[blk])[j], (limit1[blk])[n], (limit1[blk])[n+1],
                                    (limit2[blk])[m], (limit2[blk])[m+1],dedupe[blk]); //returns 2 member list of na indices
//        Rcout << "Completed natemp for feature" << j << "in round " << i << std::endl;
      }
      
      // Run m_func, initial arguments are a list of lists of 2 column matrices of matches, partial matches,
      // and a 2 member list of missing indices for each field
      
//      Rcout << "Checkpoint 2 - starting m_func in block "<< i << std::endl;
      
      mf_out = m_func(templist, ptemplist, natemplist, lims, lims_2, listid, 
                      dedupe[blk], matchesLink);
      ind_out[i] = mf_out;
    }
    
  }
  
  return ind_out;
  
}
