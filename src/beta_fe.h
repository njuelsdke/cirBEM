/// @file beta_fe.h
/// Fixed-effects beta-regression cosinor model.
#ifndef beta_fe_h
#define beta_fe_h

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type beta_fe(objective_function<Type>* obj) {
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_SCALAR(prior_mean);
  DATA_SCALAR(prior_sd);
  DATA_INTEGER(use_prior);

  PARAMETER_VECTOR(beta);
  PARAMETER(log_phi);

  Type phi = exp(log_phi);
  vector<Type> eta = X * beta;
  int n = y.size();

  Type nll = 0.0;
  for (int j = 0; j < n; j++) {
    Type mu = Type(1) / (Type(1) + exp(-eta(j)));
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
