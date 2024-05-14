// Stan version of state-space spawner-recruitment model with time-varying alpha

data{
  int T;           // time steps, for pinks this is generations, i.e. rows of data
  vector[T] S_obs; // observed spawners
  vector[T] H_obs; // observed harvest
  vector[T] S_cv;  // spawner observation error CV
  vector[T] H_cv;  // harvest observation error CV
  real pSmax_mean; // mean on prior for Smax
  real pSmax_sig;  // SD on prior for Smax
}

transformed data{
real ln_beta_pr;
real ln_beta_pr_sig;

ln_beta_pr_sig=sqrt(log(1+((1/pSmax_sig)*(1/pSmax_sig))/((1/pSmax_mean)*(1/pSmax_mean)))); //this converts sigma on the untransformed scale to a log scale
ln_beta_pr=log(1/pSmax_mean)-0.5*ln_beta_pr_sig*ln_beta_pr_sig; //convert smax prior to per capita slope - transform to log scale with bias correction
}

parameters{
  //latent states
  vector<lower=0>[T] lnR;             // log recruitment states
  vector<lower=0.01,upper=0.99>[T] U; // harvest rate

  //variance components
  real<lower = 0> sigma;
  real<lower = 0> sigma_alpha;

  //parms
  real ln_alpha0;         // initial productivity (on log scale)
  real<upper=0> ln_beta;  // Ricker b - force it negative
  vector[T-1] alpha_dev;  // time varying year-to-year deviations in a
}

transformed parameters{
  // latent states in SS model
  vector<lower=0>[T] S;  // spawner states
  vector[T] lnS;         // log spawner states
  vector<lower=0>[T] C;  // catch states
  vector[T] lnC;         // log catch states
  vector<lower=0>[T] R;  // recruitment states
  // Ricker model parms
  real beta;             // beta transformed to normal
  vector[T] ln_alpha;    // ln_alpha in each year

  beta = exp(ln_beta);
  R = exp(lnR);

  // Calculate latent states of returns, spawners and catch by return year
  for(t in 1:T){
    S[t] = R[t] * (1 - U[t]);
    lnS[t] = log(S[t]);
    C[t] = R[t] * U[t];
    lnC[t] = log(C[t]);
  }

  // random walk on time varying alpha
  ln_alpha[1] = ln_alpha0;                                    // initial value
  for(t in 2:T){
    ln_alpha[t] = ln_alpha[t-1] + alpha_dev[t-1]*sigma_alpha; // random walk of log_a
  }
}

model{
  // Priors
  ln_alpha0 ~ normal(0,5);
  ln_beta ~ normal(ln_beta_pr,ln_beta_pr_sig); //per capita capacity parameter - wide prior
  alpha_dev ~ std_normal();                    //standardized (z-scales) deviances

  //variance terms
  sigma ~ normal(0,1);       //half normal on variance (lower limit of zero)
  sigma_alpha ~ normal(0,1); //half normal on variance (lower limit of zero)

  // Likelihoods
  // Process model
  for(t in 1:T) R[t] ~ normal(lnS[t] + ln_alpha[t] - beta * S[t], sigma);

  // Observation model
  for(t in 1:T){
    U[t] ~ beta(1,1);
    H_obs[t] ~ lognormal(lnC[t], sqrt(log((H_cv[t]^2) + 1)));
    S_obs[t] ~ lognormal(lnS[t], sqrt(log((S_cv[t]^2) + 1)));
  }
}

// generate some posterior predictive checks
