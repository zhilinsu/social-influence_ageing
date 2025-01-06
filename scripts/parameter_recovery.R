# ============================================================================.
# Info ####
# ============================================================================.
# Perform parameter recovery for the KU model #0 (without any noise parameters).  
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Project: Zhilin's first PhD project - Social discounting across the adult lifespan 
# Last revision: 12 Dec 2023, Zhilin Su 

rm(list = ls())

library(rstan)
library(Hmisc)
library(ggplot2)
library(ggtext)
library(dplyr)

set.seed(1213)

# ============================================================================.
# Function definitions ####
# ============================================================================.
generate_options <- function(min_reward, max_reward) {
	temporal_options <- c(0, 1, 7, 14, 28, 42, 63, 90)
	
	# Generate all combinations
	option_pairs <- expand.grid(
		s_rs = min_reward:max_reward, # self reward sooner 
		s_ds = temporal_options, # self delay sooner
		s_rl = min_reward:max_reward, # self reward later
		s_dl = temporal_options # self delay later
	)
	
	option_pairs <- option_pairs[
		option_pairs$s_rl > 5 &
		option_pairs$s_rs < option_pairs$s_rl &
		option_pairs$s_ds < option_pairs$s_dl &
		option_pairs$s_ds == 0, ]
	
	return(option_pairs)
}

logspace <- function(a, b, n = 50) {
	
	return(10^(seq(a, b, length.out = n)))
}

compute_values <- function(k, option_pairs) {
	
	v_sooner <- option_pairs[, 1] / (1 + k * option_pairs[, 2])
	v_later <- option_pairs[, 3] / (1 + k * option_pairs[, 4])
	
	return(list(v_sooner = v_sooner, v_later = v_later))
}

update_k_values <- function(choices_count, n_trials, min_k, max_k) {
	
	if (max(choices_count) < n_trials) max_k <- max_k - 0.1
	if (min(choices_count) > 0) min_k <- min_k + 0.1
	
	return(list(min_k = min_k, max_k = max_k))
}

compute_choices <- function(min_k, max_k, n_trials, option_pairs) {
	
	choices_count <- 1
	
	while (min(choices_count) > 0 || max(choices_count) < n_trials) {
		v_sooner <- matrix(0, nrow = nrow(option_pairs), ncol = n_trials)
		v_later <- matrix(0, nrow = nrow(option_pairs), ncol = n_trials)
		for (k in logspace(min_k, max_k, n_trials)) {
			values <- compute_values(k, option_pairs)
			v_sooner <- cbind(v_sooner, values$v_sooner)
			v_later <- cbind(v_later, values$v_later)
		}
		v_chosen <- (v_later > v_sooner)
		v_chosen_info <- rowSums(v_chosen)
		sorted_choices_count <- sort(v_chosen_info, index.return = TRUE)
		choices_count <- sorted_choices_count$x
		i_choices_count <- sorted_choices_count$ix
		
		k_values <- update_k_values(choices_count, n_trials, min_k, max_k)
		min_k <- k_values$min_k
		max_k <- k_values$max_k
	}
	
	return(list(choices_count = choices_count, i_choices_count = i_choices_count))
}

generate_final_options <- function(results, option_pairs, n_trials) {
	
	final_option_pairs <- matrix(0, nrow = n_trials, ncol = ncol(option_pairs))
	colnames(final_option_pairs) <- c("s_rs", "s_ds", "s_rl", "s_dl")
	final_option_pairs <- as.data.frame(final_option_pairs)
	for (i in 1:n_trials) {
		j <- i
		any_choices_count <- which(results$choices_count == i)
		
		while(length(any_choices_count) == 0) {
			j <- j + 1
			any_choices_count <- which(results$choices_count == j)
		}
		random_row <- sample(any_choices_count, 1)
		final_option_pairs[i, ] <- option_pairs[random_row, ]
	}
	final_option_pairs <- final_option_pairs[sample(nrow(final_option_pairs)), ]
	
	return(final_option_pairs)
}

# ============================================================================.
# Generate option pairs ####
# ============================================================================.
option_pairs <- generate_options(1, 20)
results <- compute_choices(-4, 0, 50, option_pairs)

# ============================================================================.
# Initialisation ####
# ============================================================================.
parm_recovery_list <- list()

# ============================================================================.
# Set Stan parameters ####
# ============================================================================.
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

n_iter     <- 4000
n_chains   <- 4 
n_warmup   <- floor(n_iter / 2)
n_thin     <- 1
n_seed     <- 1213 # set the seed

n_trials <- 50
n_pt <- 160
n_sim <- 20

for (s in 1:n_sim) { # repeat 20 times 
	# ============================================================================.
	# Simulate data ####
	# ============================================================================.
	self_decision <- array(0, dim = c(n_pt, n_trials)) 
	reward_sooner <- array(0, dim = c(n_pt, n_trials))
	reward_later <- array(0, dim = c(n_pt, n_trials))
	delay_later <- array(0, dim = c(n_pt, n_trials))
	true_km <- array(0, n_pt)
	true_ku <- array(0, n_pt)
	
	for (i in 1:n_pt) {
		final_options <- generate_final_options(results, option_pairs, 50)
		
		self_decision[i, ] <- sample(1:2, n_trials, replace = T)
		reward_sooner[i, ] <- final_options[, "s_rs"]
		reward_later[i, ] <- final_options[, "s_rl"]
		delay_later[i, ] <- final_options[, "s_dl"]
	}
	
	simulated_model_file <- "scripts/stan_model/ku-model_self_0_simulated.stan"
	
	data_list <-  list(n_trials   = n_trials,
										 n_subjects = n_pt,
										 s_dec      = self_decision,
										 s_rs       = reward_sooner,
										 s_rl       = reward_later,
										 s_dl       = delay_later)
	
	sim_out <- stan(simulated_model_file,
									data = data_list,
									chains = n_chains, 
									iter   = n_iter, 
									warmup = n_warmup, 
									thin   = n_thin,
									init   = "random", 
									seed   = n_seed)
	
	for (i in 1:n_pt) {
		while(TRUE) {
			draw <- sample(1:(n_chains * n_warmup), 1)
			d_dec_sim <- extract(sim_out, pars = "s_dec_sim")[[1]][draw, i, ]
			
			# Check if any values in d_dec_sim are NaN
			if (!any(is.nan(d_dec_sim))) {
				true_km[i] <- extract(sim_out, pars = "km")[[1]][draw, i]
				true_ku[i] <- extract(sim_out, pars = "ku")[[1]][draw, i]
				self_decision[i, ] <- d_dec_sim
				break # exit the loop if no NaN values
			}
		}
	}
	
	# ============================================================================.
	# Recover data ####
	# ============================================================================.
	model_file <- "scripts/stan_model/ku_0-3/ku-model_self_0.stan"
	
	data_list <-  list(n_trials   = n_trials,
										 n_subjects = n_pt,
										 s_dec      = self_decision,
										 s_rs       = reward_sooner,
										 s_rl       = reward_later,
										 s_dl       = delay_later)
	
	fit <- stan(model_file,
							data = data_list,
							chains = n_chains, 
							iter   = n_iter, 
							warmup = n_warmup, 
							thin   = n_thin,
							init   = "random", 
							seed   = n_seed)
	
	#saveRDS(fit, paste0("data/recovered/recovered_data_", s, ".RDS"))  
	
	estimated_km <- summary(fit, pars = "km")$summary[, "mean"]
	estimated_ku <- summary(fit, pars = "ku")$summary[, "mean"]
	
	corr_km_km <- rcorr(true_km, estimated_km, type = "spearman")
	corr_ku_ku <- rcorr(true_ku, estimated_ku, type = "spearman")
	corr_km_ku <- rcorr(true_km, estimated_ku, type = "spearman")
	corr_ku_km <- rcorr(true_ku, estimated_km, type = "spearman")
	
	parm_recovery <- tibble(
		simulated = c("km", "km", "ku", "ku"), 
		recovered = c("km", "ku", "km", "ku"), 
		corr = c(corr_km_km$r[1, 2], 
						 corr_km_ku$r[1, 2], 
						 corr_ku_km$r[1, 2], 
						 corr_ku_ku$r[1, 2])
	)
	
	parm_recovery_list[[s]] <- parm_recovery
	
	#save(true_km, true_ku, estimated_km, estimated_ku, 
	#		 file = paste0("behaviour/processed-data/simulation_", s, ".RData"))
}

#saveRDS(parm_recovery_list, "data/result_parameter_recovery.RDS")

# ============================================================================.
# Summarise the simulations ####
# ============================================================================.
parm_recovery_list <- readRDS("data/result_parameter_recovery.RDS")

vec_corr_km_km <- numeric(n_sim)
vec_corr_km_ku <- numeric(n_sim)
vec_corr_ku_km <- numeric(n_sim)
vec_corr_ku_ku <- numeric(n_sim)

vec_corr_km_km <- vapply(parm_recovery_list, function(tb) tb[[1, "corr"]], numeric(1))
vec_corr_km_ku <- vapply(parm_recovery_list, function(tb) tb[[2, "corr"]], numeric(1))
vec_corr_ku_km <- vapply(parm_recovery_list, function(tb) tb[[3, "corr"]], numeric(1))
vec_corr_ku_ku <- vapply(parm_recovery_list, function(tb) tb[[4, "corr"]], numeric(1))

parm_recovery_final <- tibble(
	simulated = c("km", "km", "ku", "ku"), 
	recovered = c("km", "ku", "km", "ku"), 
	corr = c(mean(vec_corr_km_km), 
					 mean(vec_corr_km_ku), 
					 mean(vec_corr_ku_km), 
					 mean(vec_corr_ku_ku))
)

# ============================================================================.
# Visualisation ####
# ============================================================================.
# Create colour palette
pal <- c("#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30")

plt <- ggplot(data = parm_recovery_final, 
							aes(x = simulated, y = recovered, fill = corr)) + 
	geom_tile() + 
	scale_fill_gradientn(colours = pal,
											 limit = c(-0.25, 1.05), 
											 name = "Spearman's r") + 
	guides(fill = guide_colorbar(title.position = "right", 
															 title.hjust = 0.5, 
															 frame.colour = NULL, 
															 barheight = 16, 
															 ticks.colour = "black", 
															 ticks.linewidth = 0.5,  
															 draw.ulim = FALSE, 
															 draw.llim = FALSE, 
															 override.aes = list(limits = c(0, 1)))) + 
	labs(x = "Simulated", y = "Recovered") + 
	theme_classic() + 
	theme(text = element_text(family = "Arial"), 
				axis.text.x = element_text(size = 16, face = "italic"),
				axis.text.y = element_text(size = 16, face = "italic"), 
				axis.title = element_text(size = 16), 
				legend.title = element_text(size = 12, angle = 270), 
				legend.text = element_text(size = 12))

ggsave("data/plots/recovered.png", plt, 
			 height = 4, width = 4 * 1.3, dpi = 1200)
