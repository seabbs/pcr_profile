data {
  int<lower=0> N; // number of data points
  int<lower=0> P; // number of patients
  int patient_ID[N];
  vector[N] day_of_test; // day of test (integer)
  int test_result[N];
  int symptomatic[N];
  vector[P] time_first_symptom;
  real time_last_asym[P];
  real<lower = 0> incub_alpha;
  real<lower = 0> incub_beta;
  vector[P] time_first_test;
  // int T_e[P];
}

parameters {
  vector[P] noise;
  // real <lower = 0> alpha;
  // real <lower = 0> beta;
  real mu;
}

transformed parameters {
    vector[P] T_e; 
    for(i in 1:P) {
      T_e[i] = time_last_asym[i] + mu + noise[i];
    }
}


model {
  real upp[P];
  real low[P];

  mu ~ cauchy(0, 2.5);
  noise ~ normal(0, 1);
  
  // loop over all people
    for(j in 1:P) {
    // Symptom likelihood
    upp[j] = time_first_symptom[j] - T_e[j] <= 0 ? 0.00001 : time_first_symptom[j] - T_e[j];
    low[j] = time_last_asym[j] - T_e[j] <= 0 ? 0.00001 : time_last_asym[j] - T_e[j];
    target += gamma_lcdf(upp[j] | incub_alpha, incub_beta) - gamma_lcdf(low[j] | incub_alpha, incub_beta);
  }
  
   // alpha ~ cauchy(0, 2.5) T[0,];
   // beta  ~ cauchy(0, 2.5) T[0,];

}

  // loop over all data points
  // for(i in 1:N) {
  //   // PCR positivity likelihood
  //   diff = day_of_test[i] - T_e[patient_ID[i]];
  //   target += gamma_lpdf(diff | alpha, beta);
  // }
  
