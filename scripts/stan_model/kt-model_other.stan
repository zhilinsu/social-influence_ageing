// Implementation of the preference-temperature (KT) model 
// for the social discounting task.
// The script is for Other blocks.
// KT model #0: Hierarchical structure with uninformative priors.  
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
    real k_mean_raw; 
    real t_mean_raw; 
    real<lower=0> k_sd_raw; 
    real<lower=0> t_sd_raw; 
    
    // Subject-level raw parameters
    vector[n_subjects] k_raw; // discount rates in log10 space  
    vector[n_subjects] t_raw; // softmax temperature in log10 space
}

transformed parameters { // transform our free parameters 
    // Transform subject-level raw parameters 
    vector<upper=0>[n_subjects] k; 
    vector<lower=-1, upper=1>[n_subjects] t;
    vector[n_subjects] K; // discount rates
    vector[n_subjects] T; // softmax temperature
    
    for (s in 1:n_subjects) {
        k[s] = -log(1 + exp(k_mean_raw + k_sd_raw * k_raw[s])); 
        t[s] = Phi_approx(t_mean_raw + t_sd_raw * t_raw[s]) * 2 - 1; 
        
        K[s] = pow(10, k[s]); 
        T[s] = pow(10, t[s]); 
    }
}

model { // define the model/likelihood function 
    // Hyper-parameters priors 
    k_mean_raw ~ normal(0, 3); 
    k_sd_raw ~ cauchy(0, 2);  
    t_mean_raw ~ normal(0, 3); 
    t_sd_raw ~ cauchy(0, 2); 
    
    // Individual parameters priors 
    k_raw ~ normal(0, 1); 
    t_raw ~ normal(0, 1);
    
    for (s in 1:n_subjects) {
        // Define values 
        vector[2] v; // subjective values; v[1] is SS, v[2] is LL
        vector[2] p; // choice probability 
        
        for (i in 1:n_trials) {
            // Subjective values
            v[1] = o_rs[s, i]; 
            v[2] = o_rl[s, i] / (1 + K[s] * o_dl[s, i]); 
            
            // Softmax action prob 
            p = softmax(T[s] * v); 
            
            // choice likelihood 
            s_dec[s, i] ~ categorical(p); 
        }
    }
}

generated quantities {
	real<upper=0> k_mean; 
    real<lower=-1, upper=1> t_mean; 
    
    real log_lik[n_subjects];
    
    k_mean = -log(1 + exp(k_mean_raw));
    t_mean = Phi_approx(t_mean_raw) * 2 - 1; 
    
    {
        for (s in 1:n_subjects) {
            vector[2] v; 
            vector[2] p; 
            
            log_lik[s] = 0; 
            
            for (i in 1:n_trials) {
                v[1] = o_rs[s, i]; 
                v[2] = o_rl[s, i] / (1 + K[s] * o_dl[s, i]); 
                
                p = softmax(T[s] * v); 
                
                log_lik[s] = log_lik[s] + categorical_lpmf(s_dec[s, i] | p); 
            }
        } 
    }
}
