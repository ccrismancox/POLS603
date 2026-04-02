#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector fi(int r, int d, NumericVector eXB_i) {
  // Create storage for (r, d) pairs satisfying r >= d-1
  std::vector<int> r_vals;
  std::vector<int> d_vals;
  
  for (int ri = 0; ri <= r; ++ri) {
    for (int di = 0; di <= d; ++di) {
      if (ri >= di - 1) {
        r_vals.push_back(ri);
        d_vals.push_back(di);
      }
    }
  }
  
  int n = r_vals.size();
  NumericVector gout(n, 0.0);
  
  // Set gout where d == 0 to 1
  for (int i = 0; i < n; ++i) {
    if (d_vals[i] == 0) {
      gout[i] = 1.0;
    }
  }
  
  
  // Main recursion
  for (int r0 = 1; r0 <= r; ++r0) {
    for (int d0 = 1; d0 <= std::min(r0,d); ++d0) {
      for (int i = 0; i < n; ++i) {
        if (r_vals[i] == r0 && d_vals[i] == d0) {
          double val1 = 0.0, val2 = 0.0;
          for (int j = 0; j < n; ++j) {
            if (r_vals[j] == r0 - 1 && d_vals[j] == d0)
              val1 = gout[j];
            if (r_vals[j] == r0 - 1 && d_vals[j] == d0 - 1)
              val2 = gout[j]* eXB_i[r0-1];
          }
          
          gout[i] =  val1 + val2 ; 
        }
      }
    }
  }
  
  return gout;
  
  
}


// [[Rcpp::export]]
NumericVector f(const NumericVector T, const NumericVector k, List& eXB) {
  int N = eXB.length();
  NumericVector out (N);
  
  for(int i=0; i < N; i ++){
    NumericVector fout = fi(T[i], k[i], eXB[i]);
    int n = fout.length();
    out[i] = log(fout[n-1]);
  }
  return out;
}



// [[Rcpp::export]]
NumericVector D_fi(int r, int d, NumericVector eXB_i, NumericMatrix Xi) {
  // Build (r, d) grid where r >= d - 1
  std::vector<int> r_vals;
  std::vector<int> d_vals;
  
  for (int ri = 0; ri <= r; ++ri) {
    for (int di = 0; di <= d; ++di) {
      if (ri >= di - 1) {
        r_vals.push_back(ri);
        d_vals.push_back(di);
      }
    }
  }
  
  int nrow_out = r_vals.size();
  int ncol_out = Xi.ncol();
  
  // Initialize output matrix
  NumericMatrix Dgout(nrow_out, ncol_out);
  NumericVector fout = fi(r, d, eXB_i);
  
  // Recursive update
  for (int r0 = 1; r0 <= r; ++r0) {
    for (int d0 = 1; d0 <= std::min(r0, d); ++d0) {
      // fi0(r0-1, d0-1, eXB.j)
      
      
      for (int i = 0; i < nrow_out; ++i) {
        if (r_vals[i] == r0 && d_vals[i] == d0) {
          NumericVector term1(ncol_out), term2(ncol_out), term3(ncol_out);
          double exp_val = eXB_i[r0-1];
          
          for(int j =0; j < nrow_out; j++){
            if (r_vals[j] == r0 - 1 && d_vals[j] == d0) {
              
              term2 = Dgout(j, _);
            }
            if (r_vals[j] == r0 - 1 && d_vals[j] == d0 - 1) {
              term1 = fout[j] * exp_val * Xi(r0-1, _ );
              term3 = Dgout(j, _) * exp_val;
              
            }
          }
          Dgout(i, _) = term1 + term2 + term3;
        }
      }
    }
  }
  // Return last row
  return Dgout(nrow_out - 1, _);
}



// [[Rcpp::export]]
NumericMatrix D_f(const NumericVector& T, const NumericVector& k, 
                  const List & eXB, List X, const int p) {
  int N = eXB.length();
  int i = 0;
  NumericMatrix out (p, N);
  NumericVector numer (p);
  double denom;
  for(i=0; i < N; i ++){
    NumericVector eXBi = eXB[i];
    NumericVector fout = fi(T[i],k[i],eXBi);
    int n = fout.length();
    denom = fout[n-1];
    NumericMatrix Xi = X[i];
    numer  = D_fi(T[i],k[i],eXBi, Xi);
    out(_, i) = numer/denom;
  }
  return  out;
}
