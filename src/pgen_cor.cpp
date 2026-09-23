#include <Rcpp.h>
#include <cmath>
#include <cstddef>

// [[Rcpp::export]]
Rcpp::NumericMatrix pgen_pair_cor_cpp(const Rcpp::NumericMatrix& A,
                                      const Rcpp::NumericMatrix& B,
                                      bool self) {
  if (A.nrow() != B.nrow()) Rcpp::stop("PGEN blocks have different sample counts.");
  if (self && A.ncol() != B.ncol()) Rcpp::stop("A self-correlation block must be square.");

  const std::size_t n = A.nrow();
  Rcpp::NumericMatrix out(A.ncol(), B.ncol());
  for (int j = 0; j < B.ncol(); ++j) {
    const double* y = B.begin() + n * static_cast<std::size_t>(j);
    const int end = self ? j + 1 : A.ncol();
    for (int i = 0; i < end; ++i) {
      const double* x = A.begin() + n * static_cast<std::size_t>(i);
      std::size_t observed = 0;
      double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
      for (std::size_t k = 0; k < n; ++k) {
        if (!std::isfinite(x[k]) || !std::isfinite(y[k])) continue;
        ++observed;
        sx += x[k];
        sy += y[k];
        sxx += x[k] * x[k];
        syy += y[k] * y[k];
        sxy += x[k] * y[k];
      }
      double r = NA_REAL;
      if (observed >= 2) {
        const double count = static_cast<double>(observed);
        const double vx = sxx - sx * sx / count;
        const double vy = syy - sy * sy / count;
        if (vx > 0.0 && vy > 0.0) {
          r = (sxy - sx * sy / count) / std::sqrt(vx * vy);
          if (r > 1.0) r = 1.0;
          if (r < -1.0) r = -1.0;
        }
      }
      out(i, j) = r;
      if (self) out(j, i) = r;
    }
    if ((j & 15) == 0) Rcpp::checkUserInterrupt();
  }
  return out;
}
