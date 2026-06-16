#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

static double quantile_type7(std::vector<double> values, double p) {
  if (values.empty()) {
    Rcpp::stop("Cannot compute a quantile of an empty vector.");
  }

  std::sort(values.begin(), values.end());

  double h = (values.size() - 1) * p;
  size_t lo = static_cast<size_t>(std::floor(h));
  size_t hi = static_cast<size_t>(std::ceil(h));
  double frac = h - lo;

  return values[lo] + frac * (values[hi] - values[lo]);
}

static std::vector<double> column_to_vector(const arma::mat& X, arma::uword j) {
  std::vector<double> values(X.n_rows);

  for (arma::uword i = 0; i < X.n_rows; ++i) {
    values[i] = X(i, j);
  }

  return values;
}

// [[Rcpp::export]]
arma::mat matrixInverse_cpp(const arma::mat& A) {
  return inv(A);
}

// [[Rcpp::export]]
arma::mat matrixMultiply_cpp(const arma::mat& A,
                             const arma::mat& B,
                             bool transA = false,
                             bool transB = false) {
  if (transA && transB) {
    return A.t() * B.t();
  }
  if (transA) {
    return A.t() * B;
  }
  if (transB) {
    return A * B.t();
  }

  return A * B;
}

// [[Rcpp::export]]
arma::vec matrixVectorMultiply_cpp(const arma::mat& A, const arma::vec& b) {
  return A * b;
}

// [[Rcpp::export]]
arma::mat matrixCor_cpp(const arma::mat& A) {
  int n = A.n_rows;
  arma::mat B = A.each_row() - mean(A, 0);
  arma::mat S = trans(B) * B / (n - 1);
  arma::vec D = 1.0 / sqrt(S.diag());
  D.replace(arma::datum::nan, 0);
  arma::mat R = S.each_col() % D;
  R = R.each_row() % trans(D);
  return R;
}

// [[Rcpp::export]]
Rcpp::List matrixEigen_cpp(const arma::mat& A) {
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A);
  arma::uvec indices = sort_index(eigval, "descend");
  eigval = eigval(indices);
  eigvec = eigvec.cols(indices);
  return Rcpp::List::create(Rcpp::Named("values") = eigval,
                            Rcpp::Named("vectors") = eigvec);
}

// [[Rcpp::export]]
arma::mat matrixKronecker_cpp(const arma::mat& A, const arma::mat& B) {
  return arma::kron(A, B);
}

// [[Rcpp::export]]
arma::mat matrixListProduct_cpp(const Rcpp::List& matrixList) {
  int n = matrixList.size();
  arma::mat result = Rcpp::as<arma::mat>(matrixList[0]);

  for (int i = 1; i < n; i++) {
    result *= Rcpp::as<arma::mat>(matrixList[i]);
  }

  return result;
}

// [[Rcpp::export]]
arma::mat matrixGeneralizedInverse_cpp(const arma::mat& A, double tol = 5e-16) {
  return arma::pinv(A, tol);
}

// [[Rcpp::export]]
arma::mat matrixSylvesterEigen_cpp(const arma::mat& UA,
                         const arma::vec& DA,
                         const arma::mat& UB,
                         const arma::vec& DB,
                         const arma::mat& C)
{
arma::mat tC = UA.t() * C * UB;

arma::mat tX(tC.n_rows, tC.n_cols, arma::fill::none);
for (size_t i = 0; i < tC.n_rows; ++i) {
  for (size_t j = 0; j < tC.n_cols; ++j) {
    tX(i, j) = tC(i, j) / (DA(i) + DB(j));
  }
}

arma::mat X = UA * tX * UB.t();
return X;
}

// [[Rcpp::export]]
arma::mat matrixSylvester_cpp(const arma::mat& A,
                        const arma::mat& B,
                        const arma::mat& C)
{
  arma::mat X;
  bool success = arma::sylvester(X, A, B, -C);

  if(!success) {
    Rcpp::stop("Sylvester equation solver did not converge or was not successful.");
  }

  return X;
}

// [[Rcpp::export]]
Rcpp::List matrixSVD_cpp(const arma::mat& A) {
  arma::mat U, V;
  arma::vec s;

  bool success = arma::svd_econ(U, s, V, A, "both");

  if (!success) {
    Rcpp::stop("SVD computation failed.");
  }

  return Rcpp::List::create(
    Rcpp::Named("d") = s,   //
    Rcpp::Named("u") = U,   // U: m × min(m,p)
    Rcpp::Named("v") = V    // V: p × min(m,p)
  );
}

// [[Rcpp::export]]
arma::mat matrixSolveMat_cpp(const arma::mat& A, const arma::mat& B) {
  return arma::solve(A, B);
}

// [[Rcpp::export]]
Rcpp::List matrixScale_cpp(const arma::mat& X,
                           bool center = true,
                           bool standardized = true,
                           bool robust = false) {
  const double lower = 0.015;
  const double upper = 0.975;
  const double winsor_inflation = 1.036992766476;

  arma::mat result = X;
  arma::rowvec centers(X.n_cols, arma::fill::zeros);
  arma::rowvec scales(X.n_cols, arma::fill::ones);

  if (center || standardized) {
    for (arma::uword j = 0; j < X.n_cols; ++j) {
      if (center) {
        if (robust) {
          centers(j) = quantile_type7(column_to_vector(X, j), 0.5);
        } else {
          centers(j) = arma::mean(X.col(j));
        }
      }

      if (standardized) {
        if (robust) {
          std::vector<double> values = column_to_vector(X, j);
          double lo = quantile_type7(values, lower);
          double hi = quantile_type7(values, upper);
          arma::vec winsorized(X.n_rows);

          for (arma::uword i = 0; i < X.n_rows; ++i) {
            winsorized(i) = std::min(std::max(X(i, j), lo), hi);
          }

          scales(j) = arma::stddev(winsorized, 0) * winsor_inflation;
        } else {
          scales(j) = arma::stddev(X.col(j), 0);
        }
      }
    }
  }

  if (center) {
    result.each_row() -= centers;
  }

  if (standardized) {
    result.each_row() /= scales;
  }

  return Rcpp::List::create(
    Rcpp::Named("x") = result,
    Rcpp::Named("center") = centers,
    Rcpp::Named("scale") = scales
  );
}
