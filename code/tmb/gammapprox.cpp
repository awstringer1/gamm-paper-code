// This file fits a GAMM

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // DATA
  DATA_VECTOR(y); // Full vector of data
  int n = y.size();
  DATA_MATRIX(X); // Spline design matrix, single covariate for now
  int p = X.cols();
  if (X.rows() != n) throw std::runtime_error("Inappropriate number of rows in X");
  DATA_SPARSE_MATRIX(Z); // Random intercept design matrix
  int m = Z.cols();
  if (Z.rows() != n) throw std::runtime_error("Inappropriate number of rows in Z");
  DATA_MATRIX(S); // Penalty matrix. Store as dense matrix, it's small enough that this should be
                  // efficient, and its dimension doesn't scale with n
  
  // PARAMETERS
  PARAMETER(logprec); // -2*log(random effects standard deviation)
  Type sigma = exp(-logprec / 2.);
  PARAMETER(loglambda); // log(smoothing penalty parameter), single function only for now
  Type lambda = exp(loglambda);
  PARAMETER(alpha); // Intercept
  PARAMETER_VECTOR(beta); // Spline weights
  if (beta.size() != p) throw std::runtime_error("Inappropriate number of columns in X");
  PARAMETER_VECTOR(u); // Random effects
  if (u.size() != m) throw std::runtime_error("Inappropriate number of columns in Z");

  // Bernoulli log-likelihood
  vector<Type> eta = X * beta + Z * u; // Linear predictor
  for (int i = 0; i < n; i++) eta[i] += alpha;
  vector<Type> mu(n);
  for (int i = 0; i < n; i++) mu[i] = 1. / (1. + exp(-eta[i]));
  Type nll = 0;
  for (int i = 0; i < n; i++) nll -= dbinom(y[i], Type(1.), mu[i], 1); // 1 == log  

  // Normal prior
  Type nlrep = 0.;
  for (int i = 0; i < m; i++) nlrep -= dnorm(u[i], Type(0.), sigma, 1); // 1 == log

  // Spline penalty
  Type nlsp = -((p - 1) / 2.) * loglambda + ((p - 1) / 2.) * log(2. * 3.141593); // Initialize it with the Gaussian normalizing constant part that I forget to add every time I do this.
  nlsp += 0.5 * lambda * ((S * beta) * beta).sum();

  return nll + nlrep + nlsp;
}
