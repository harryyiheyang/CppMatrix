#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat matrixInverse(const arma::mat& A) {
  return inv(A);
}

// [[Rcpp::export]]
arma::mat matrixMultiply(const arma::mat& A, const arma::mat& B) {
  return A * B;
}

// [[Rcpp::export]]
arma::vec matrixVectorMultiply(const arma::mat& A, const arma::vec& b) {
  return A * b;
}

// [[Rcpp::export]]
arma::mat matrixCor(const arma::mat& A) {
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
Rcpp::List matrixEigen(const arma::mat& A) {
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
arma::mat matrixKronecker(const arma::mat& A, const arma::mat& B) {
  return arma::kron(A, B);
}

// [[Rcpp::export]]
arma::mat matrixListProduct(const Rcpp::List& matrixList) {
  int n = matrixList.size();
  arma::mat result = Rcpp::as<arma::mat>(matrixList[0]);

  for (int i = 1; i < n; i++) {
    result *= Rcpp::as<arma::mat>(matrixList[i]);
  }

  return result;
}

// [[Rcpp::export]]
arma::mat matrixGeneralizedInverse(const arma::mat& A, double tol = 5e-16) {
  return arma::pinv(A, tol);
}

// [[Rcpp::export]]
arma::mat matrixSylvesterEigen(const arma::mat& UA,
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
arma::mat matrixSylvester(const arma::mat& A,
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
Rcpp::List matrixSVD(const arma::mat& A) {
  arma::mat U, V;
  arma::vec s;

  bool use_econ = (A.n_rows * A.n_cols > 100000000);

  bool success;
  if (use_econ) {
    success = arma::svd_econ(U, s, V, A, "both");
  } else {
    success = arma::svd(U, s, V, A);
  }

  if (!success) {
    Rcpp::stop("SVD failed.");
  }

  return Rcpp::List::create(
    Rcpp::Named("d") = s,
    Rcpp::Named("u") = U,
    Rcpp::Named("v") = V
  );
}
