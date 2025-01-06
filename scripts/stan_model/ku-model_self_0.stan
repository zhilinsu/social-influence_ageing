// Implementation of the preference-uncertainty (KU) model 
// for the social discounting task.
// KU Model #0: Hierarchical structure with uninformative priors. 
// 
// Ref: Moutoussis et al., 2016, *PLoS Comput. Biol.*
// 
// Zhilin Su 
// zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
//
// Project: Zhilin's first PhD project - Social discounting across the adult lifespan. 
// Last revision: 20 Nov 2023, Zhilin Su. 

data { // define our data 
    int<lower=1> n_trials; 
    int<lower=1> n_subjects; 
    int<lower=1, upper=2> s_dec[n_subjects, n_trials]; // participants' decisions 
    real<lower=0> s_rs[n_subjects, n_trials]; // smaller and sooner reward for oneself
    real<lower=0> s_rl[n_subjects, n_trials]; // larger and later reward for oneself 
    real<lower=0> s_dl[n_subjects, n_trials]; // delay period in days for oneself 
}

parameters { // define our free parameters 
    // Hyper(group)-parameters 
    real km_mean_raw; 
    real ku_mean_raw; 
    real<lower=0> km_sd_raw; 
    real<lower=0> ku_sd_raw;
    
    // Subject-level parameters 
    vector[n_subjects] km_raw; // the mean of normal distribution for k (in log10 space)
    vector[n_subjects] ku_raw; // the standard deviation of normal distribution for k (in log10 space)
}

transformed parameters { // transform our parameters
    // Transform subject-level raw parameters 
    vector<upper=0>[n_subjects] km; 
    vector<lower=0>[n_subjects] ku;
    
    for (s in 1:n_subjects) {
        km[s] = -log(1 + exp(km_mean_raw + km_sd_raw * km_raw[s]));
    		ku[s] = log(1 + exp(ku_mean_raw + ku_sd_raw * ku_raw[s]));  
    }
}

model { // define the model/likelihood function 
    // Hyper-parameters priors 
    km_mean_raw ~ normal(0, 3); 
    km_sd_raw ~ cauchy(0, 2); 
    ku_mean_raw ~ normal(0, 3); 
    ku_sd_raw ~ cauchy(0, 2); 
    
    // Individual parameters priors 
    km_raw ~ normal(0, 1); 
    ku_raw ~ normal(0, 1);
    
    for (s in 1:n_subjects) {
        // Define values 
        vector[2] p; // choice probability 
        
        for (i in 1:n_trials) {
            p[2] = normal_cdf(log((s_rl[s, i] / s_rs[s, i] - 1)/ s_dl[s, i]), km[s], ku[s]); 
            p[1] = 1 - p[2]; 
            
            s_dec[s, i] ~ categorical(p); 
        }
    }
} 

generated quantities {
    real<upper=0> km_mean; 
    real<lower=0> ku_mean; 
    
    real log_lik[n_subjects];
    
    km_mean = -log(1 + exp(km_mean_raw)); 
    ku_mean = log(1 + exp(ku_mean_raw));  
    
    {
        for (s in 1:n_subjects) {
            vector[2] p; 
            
            log_lik[s] = 0; 
            
            for (i in 1:n_trials) {
                p[2] = normal_cdf(log((s_rl[s, i] / s_rs[s, i] - 1)/ s_dl[s, i]), km[s], ku[s]); 
                p[1] = 1 - p[2]; 
                
                log_lik[s] = log_lik[s] + categorical_lpmf(s_dec[s, i] | p); 
            }
        }
    }
}
