#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]




List assign_mos_haps(NumericVector moi,
                      NumericVector freq,
                      NumericVector haps){
  int n=moi.size();
  List infec_m(n);
  NumericVector infection;
  for(int i=0; i<n; i++){
     infec_m[i] = sample(haps,moi[i],FALSE,freq);
   
  }
  
  return(infec_m);
}


// [[Rcpp::export]]

NumericVector get_min_age_haps(NumericMatrix age_haps){
  int n=age_haps.nrow();
  NumericVector min_age_haps_nonzero;
  for(int i=0; i<n; i++){
    NumericVector age_haps_nonzero;
    NumericVector age_haps_mos;
    
    age_haps_mos = age_haps.row(i);
    age_haps_nonzero = age_haps_mos[age_haps_mos!=0];
    if(age_haps_nonzero.size()!=0){
    min_age_haps_nonzero[i] = min(age_haps_nonzero);
    }else{
      min_age_haps_nonzero[i]=0 ;
    }
    
  }
  
  return(min_age_haps_nonzero);
}





// [[Rcpp::export]]
NumericMatrix increment_non_zero_age(NumericMatrix x) {
  for(int i=0; i<x.rows(); i++){
    for(int j=0; j<x.cols(); j++){
      if(x(i, j) != 0){
        x(i, j)++;
      }
    }
  }
  return(x);
}



// [[Rcpp::export]]
NumericVector get_biting_status(NumericVector x) {
  return x[x > 0];
}

// [[Rcpp::export]]
List first_three(List x) {
  IntegerVector idx = IntegerVector::create(0, 1, 2);
  return x[idx];
}

// [[Rcpp::export]]
List with_names(List x, CharacterVector y) {
  return x[y];
}



// [[Rcpp::export]]
NumericVector between(NumericVector x, double s, double e) {
  return x[(x >= s)*(x <= e)];
}



// [[Rcpp::export]]
NumericVector count_infs(NumericVector x, NumericVector s_vec, 
                         NumericVector e_vec){
  NumericVector contain_vec(x.size());

  for(int i=0; i<s_vec.size(); i++){
    for(int j=0; j<x.size(); j++){
      if(x(j) <= e_vec(i) && x(j) >= s_vec(i)){
        contain_vec(j) += 1;
      }
    }
  }
  
  return contain_vec;
}



// [[Rcpp::export]]
NumericMatrix count_infs_all_people(NumericVector x, List s_vec_all, 
                         List e_vec_all,
                         int n_people){
  NumericMatrix contain_vec(x.size(), n_people);
  
  for(int person_i=0; person_i<n_people; person_i++){
    NumericVector s_vec = s_vec_all(person_i);
    NumericVector e_vec = e_vec_all(person_i);
    for(int i=0; i<s_vec.size(); i++){
      for(int j=0; j<x.size(); j++){
        if(x(j) <= e_vec(i) && x(j) >= s_vec(i)){
          contain_vec(j, person_i) += 1;
        }
      }
    }
  }

  
  return contain_vec;
}



// [[Rcpp::export]]
NumericVector get_infs_including_positive(NumericVector x, 
                                          List s_vec_all, 
                                          List e_vec_all,
                                          List positive_times,
                                          int n_people){
  NumericMatrix contain_vec(x.size(), n_people);
  int i_x=0;
  for(int person_i=0; person_i<n_people; person_i++){
    NumericVector s_vec = s_vec_all(person_i);
    NumericVector e_vec = e_vec_all(person_i);
    NumericVector positive_times_i = positive_times(person_i);
    
    for(int i=0; i<s_vec.size(); i++){
      for(int j=0; j<positive_times_i.size(); j++){
        if(positive_times_i(j) <= e_vec(i) && positive_times_i(j) >= s_vec(i)){
          x(i_x) = e_vec(i) - s_vec(i);
          i_x++;
          break;
        }
      }
    }
  }
  
  
  return x;
}


