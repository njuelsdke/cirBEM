/// @file beta_glmm.h
/// Mixed-effects beta-regression cosinor model with a random subject intercept.
#ifndef beta_glmm_h
#define beta_glmm_h

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type beta_glmm(objective_function<Type>* obj) {
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_IVECTOR(id);
  DATA_SCALAR(prior_mean);
  DATA_SCALAR(prior_sd);
  DATA_INTEGER(use_prior);

  PARAMETER_VECTOR(beta);
  PARAMETER(log_phi);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(u);

  Type phi = exp(log_phi);
  Type tau = exp(log_tau);
  vector<Type> eta = X * beta;
  int n = y.size();

  Type nll = 0.0;
  for (int i = 0; i < u.size(); i++) {
    nll -= dnorm(u(i), Type(0), tau, true);
  }
  for (int j = 0; j < n; j++) {
    Type mu = Type(1) / (Type(1) + exp(-(eta(j) + u(id(j)))));
    Type a = mu * phi;
    Type b = (Type(1) - mu) * phi;
    Type ld = lgamma(a + b) - lgamma(a) - lgamma(b)
              + (a - Type(1)) * log(y(j))
              + (b - Type(1)) * log(Type(1) - y(j));
    nll -= ld;
  }
  if (use_prior) nll -= dnorm(log_phi, prior_mean, prior_sd, true);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
