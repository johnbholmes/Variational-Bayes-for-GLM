#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// [[Rcpp::export]]
static inline double normal_cdf(const double x) {
  // Standard normal CDF using erfc for numerical stability
  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// [[Rcpp::export]]
std::vector<double> EPdiff01_logitnormalkpartcplus(double mu, double sigma) {
  // In part one we calculate the component functions.
  const double sigma2 = sigma * sigma;

  double sum_EP_terms = 0.0;
  double sum_EdiffP = 0.0;
  double L = 1.507306;
  std::size_t k = 11;
  std::vector<double> coef_nearzero = {0,0.999991820, -0.499903021,  0.332646462, -0.246745120,  0.189036359, -0.139389342,  0.091239940, -0.048514467,  0.018953113, -0.004725329, 0.000556454};
 
  for (std::size_t i = 1; i <= k; ++i) {
    const double N = static_cast<double>(i);

    double p0a = std::exp(N * mu + 0.5 * N * N * sigma2);
    double p0b = std::exp(-N * mu + 0.5 * N * N * sigma2);

    // if element is infinity, the cdf to be multiplied
    // is zero. Adding this line stops NaN.
    if (!std::isfinite(p0a)) p0a = 0.0;
    if (!std::isfinite(p0b)) p0b = 0.0;

    const double p0c = coef_nearzero[i] * N;

    const double p1 = ((i % 2 == 1) ? 1.0 : -1.0) - p0c; // (-1)^(N-1) - p0c

    const double p2 = p0a * normal_cdf((-L - mu - N * sigma2) / sigma);
    const double p3 = p0b * normal_cdf((mu - L - N * sigma2) / sigma);
    const double p4 = p0c * (p0a * normal_cdf((-mu - N * sigma2) / sigma));
    const double p5 = p0c * (p0b * normal_cdf((mu - N * sigma2) / sigma));

    // In part two, we use the results of part one to find E(P) and the first three derivatives.
    sum_EP_terms += p1 * (p2 - p3) + p4 - p5;
    sum_EdiffP += N * (p1 * (p2 + p3) + p4 + p5);
  }

  // In part three, calculate the moments needed for EP;
  const double EP = normal_cdf(mu / sigma) + sum_EP_terms;


  std::vector<double> results;
  results.reserve(2);
  results.push_back(EP);
  results.push_back(sum_EdiffP);

  return results;
}