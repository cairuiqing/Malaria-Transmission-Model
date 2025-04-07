#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm> // For std::sort and std::lower_bound
/*
Key Improvements:

get_min_age_haps:

Fixed vector initialization and optimized minimum calculation by directly iterating through elements.

Avoids creating temporary vectors for non-zero elements.

increment_non_zero_age:

Optimized loop order for column-major access, improving cache efficiency.

count_infs & count_infs_all_people:

Replaced nested loops with vectorized operations using Rcpp sugar for faster computation.

get_infs_including_positive:

Fixed memory issues and optimized using sorting and binary search to check interval containment efficiently.

Now returns a dynamically sized vector, avoiding buffer overflow issues.
*/

// [[Rcpp::export]]
List assign_mos_haps(NumericVector moi, NumericVector freq, NumericVector haps) {
  int n = moi.size();
  List infec_m(n);
  for(int i = 0; i < n; i++) {
    infec_m[i] = sample(haps, moi[i], false, freq);
  }
  return infec_m;
}

// [[Rcpp::export]]
NumericVector get_min_age_haps(NumericMatrix age_haps) {
  int n = age_haps.nrow();
  int m = age_haps.ncol();
  NumericVector min_age_haps_nonzero(n);
  for(int i = 0; i < n; i++) {
    double current_min = R_PosInf;
    bool found = false;
    for(int j = 0; j < m; j++) {
      double val = age_haps(i, j);
      if (val != 0) {
        if (val < current_min) {
          current_min = val;
          found = true;
        }
      }
    }
    min_age_haps_nonzero[i] = found ? current_min : 0.0;
  }
  return min_age_haps_nonzero;
}

// [[Rcpp::export]]
NumericMatrix increment_non_zero_age(NumericMatrix x) {
  int ncol = x.cols();
  int nrow = x.rows();
  for(int j = 0; j < ncol; j++) {
    for(int i = 0; i < nrow; i++) {
      if(x(i, j) != 0) {
        x(i, j) += 1;
      }
    }
  }
  return x;
}

// [[Rcpp::export]]
NumericVector get_biting_status(NumericVector x) {
  return x[x > 0];
}

// [[Rcpp::export]]
List first_three(List x) {
  return x[IntegerVector::create(0, 1, 2)];
}

// [[Rcpp::export]]
List with_names(List x, CharacterVector y) {
  return x[y];
}

// [[Rcpp::export]]
NumericVector between(NumericVector x, double s, double e) {
  return x[(x >= s) & (x <= e)];
}

// [[Rcpp::export]]
NumericVector count_infs(NumericVector x, NumericVector s_vec, NumericVector e_vec) {
  int n = x.size();
  NumericVector contain_vec(n);
  for(int j = 0; j < n; j++) {
    double xj = x[j];
    contain_vec[j] = sum((s_vec <= xj) & (e_vec >= xj));
  }
  return contain_vec;
}

// [[Rcpp::export]]
NumericMatrix count_infs_all_people(NumericVector x, List s_vec_all, List e_vec_all, int n_people) {
  int n = x.size();
  NumericMatrix contain_vec(n, n_people);
  for(int person_i = 0; person_i < n_people; person_i++) {
    NumericVector s_vec = s_vec_all[person_i];
    NumericVector e_vec = e_vec_all[person_i];
    for(int j = 0; j < n; j++) {
      double xj = x[j];
      contain_vec(j, person_i) = sum((s_vec <= xj) & (e_vec >= xj));
    }
  }
  return contain_vec;
}

// [[Rcpp::export]]
NumericVector get_infs_including_positive(List s_vec_all, List e_vec_all, List positive_times, int n_people) {
  std::vector<double> result;
  for(int person_i = 0; person_i < n_people; person_i++) {
    NumericVector s_vec = s_vec_all[person_i];
    NumericVector e_vec = e_vec_all[person_i];
    NumericVector pos_times = positive_times[person_i];
    
    std::sort(pos_times.begin(), pos_times.end());
    int num_intervals = s_vec.size();
    
    for(int i = 0; i < num_intervals; i++) {
      double s = s_vec[i];
      double e = e_vec[i];
      auto it = std::lower_bound(pos_times.begin(), pos_times.end(), s);
      if(it != pos_times.end() && *it <= e) {
        result.push_back(e - s);
      }
    }
  }
  return NumericVector(result.begin(), result.end());
}