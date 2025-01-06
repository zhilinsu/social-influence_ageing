// Implementation of the preference-uncertainty (KU) model 
// for the social discounting task.
// The script is for Other blocks.
// KU Model #2: Hierarchical structure with other-noise parameters (tau_o), 
//              and with uninformative priors 
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
    real<lower=0> o_rs[n_subjects, n_trials]; // smaller and sooner reward for others
    real<lower=0> o_rl[n_subjects, n_trials]; // larger and later reward for others
    real<lower=0> o_dl[n_subjects, n_trials]; // delay period in days for others
}

parameters { // define our free parameters 
    // Hyper(group)-parameters 
    real km_mean_raw; 
    real ku_mean_raw; 
    real tau_o_mean_raw; 
    real<lower=0> km_sd_raw; 
    real<lower=0> ku_sd_raw;
    real<lower=0> tau_o_sd_raw;
    
    // Subject-level parameters 
    vector[n_subjects] km_raw; // the mean of normal distribution for k (in log10 space)
    vector[n_subjects] ku_raw; // the standard deviation of normal distribution for k (in log10 space)
    vector[n_subjects] tau_o_raw; // other-noise parameters 
}

transformed parameters {
    // Transform subject-level raw parameters 
    vector<upper=0>[n_subjects] km; 
    vector<lower=0>[n_subjects] ku;
    vector<lower=0, upper=10>[n_subjects] tau_o; 
    
    for (s in 1:n_subjects) {
        km[s] = -log(1 + exp(km_mean_raw + km_sd_raw * km_raw[s])); 
        ku[s] = log(1 + exp(ku_mean_raw + ku_sd_raw * ku_raw[s])); 
        tau_o[s] = Phi_approx(tau_o_mean_raw + tau_o_sd_raw * tau_o_raw[s]) * 10; 
    }
}

model { // define the model/likelihood function
    // Hyper-parameters priors
    km_mean_raw ~ normal(0, 3); 
    km_sd_raw ~ cauchy(0, 2); 
    ku_mean_raw ~ normal(0, 3); 
    ku_sd_raw ~ cauchy(0, 2); 
    tau_o_mean_raw ~ normal(0, 7.5); 
    tau_o_sd_raw ~ cauchy(0, 5);
    
    // Individual parameters priors 
    km_raw ~ normal(0, 1); 
    ku_raw ~ normal(0, 1); 
    tau_o_raw ~ normal(0, 1); 
    
    for (s in 1:n_subjects) {
        // Define values 
        vector[2] p; // choice probability 
        
        real inverse_tau_o; 
    	  inverse_tau_o = 1 / tau_o[s]; 
        
        for (i in 1:n_trials) {
        	
            p[2] = normal_cdf(log((o_rl[s, i] / o_rs[s, i] - 1)/ o_dl[s, i]), km[s], ku[s]);
            p[2] = pow(p[2], inverse_tau_o) / (pow(p[2], inverse_tau_o) + pow(1 - p[2], inverse_tau_o)); 
            p[1] = 1 - p[2]; 
            
            s_dec[s, i] ~ categorical(p); 
        }
    }
} 

generated quantities {
    real<upper=0> km_mean; 
    real<lower=0> ku_mean;
    real<lower=0, upper=10> tau_o_mean; 
    
    real log_lik[n_subjects];
    
    km_mean = -log(1 + exp(km_mean_raw));
    ku_mean = log(1 + exp(ku_mean_raw));
    tau_o_mean = Phi_approx(tau_o_mean_raw) * 10 ;
    
    {
        for (s in 1:n_subjects) {
            vector[2] p; 
            
            log_lik[s] = 0; 
            
            for (i in 1:n_trials) {
                p[2] = normal_cdf(log((o_rl[s, i] / o_rs[s, i] - 1)/ o_dl[s, i]), km[s], ku[s]); 
                p[2] = pow(p[2], (1 / tau_o[s])) / (pow(p[2], (1 / tau_o[s])) + pow(1 - p[2], (1 / tau_o[s])));  
                p[1] = 1 - p[2]; 
                
                log_lik[s] = log_lik[s] + categorical_lpmf(s_dec[s, i] | p); 
            }
        }
    }
}