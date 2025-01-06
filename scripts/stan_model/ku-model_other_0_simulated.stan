// Implementation of the preference-uncertainty (KU) model 
// for the social discounting task.
// KU Model #0: Hierarchical structure with uninformative priors. 
// This is for simulating data in Other blocks, using individual-level posteriors. 
// 
// Ref: Moutoussis et al., 2016, *PLoS Comput. Biol.*
// 
// Zhilin Su 
// zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
//
// Project: Zhilin's first PhD project - Social discounting across the adult lifespan. 
// Last revision: 21 August 2024, Zhilin Su.

// Define our data
data { 
  int<lower=1> n_trials; 
  int<lower=1> n_subjects; 
  real<lower=0> o_rs[n_subjects, n_trials]; // smaller and sooner reward for others
  real<lower=0> o_rl[n_subjects, n_trials]; // larger and later reward for others
  real<lower=0> o_dl[n_subjects, n_trials]; // delay period in days for others 
  vector<upper=0>[n_subjects] km; 
  vector<lower=0>[n_subjects] ku; 
}

// Generate data according to the hierarchical model 
generated quantities {
	// Decisions to be generated
	int<lower=1, upper=2> s_dec_sim[n_subjects, n_trials]; 
  
  for (s in 1:n_subjects) {
  	vector[2] p; 
    
    for (i in 1:n_trials) {
    	p[2] = normal_cdf(log((o_rl[s, i] / o_rs[s, i] - 1)/ o_dl[s, i]), km[s], ku[s]); 
    	p[1] = 1 - p[2]; 
    	
    	s_dec_sim[s, i] = categorical_rng(p);  
    } // trial loop
  } // subject loop 
}
