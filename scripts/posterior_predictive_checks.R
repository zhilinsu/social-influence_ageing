# ============================================================================.
# Info ####
# ============================================================================.
# The script is to perform posterior predictive checks for the KU model #0. 
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Social influence and ageing project. 
# Last revision: 21 August 2024, Zhilin Su 

rm(list = ls())

library(rstan)
library(tidyverse)
library(lme4)
library(kableExtra)
library(ggsignif)
library(BayesFactor)
source("scripts/helper_functions/helper_functions.R")

set.seed(1213)

# 1 for young impulsive; 2 for older impulsive; 
# 3 for young patient; 4 for older patient. 
acc_result_list <- vector(mode = "list", length = 4)

# 1 for healthy young adults, and 2 for healthy older adults
pop <- 2

# ============================================================================.
# Preparation #### 
# ============================================================================.
load("data/sd_data.RData")

if (pop == 1) {
	pt <- sd_ya
	other_k <- k_ya 
	excl_id <- c(110, 111, 118, 158)
} else if (pop == 2) {
	pt <- sd_oa
	other_k <- k_oa 
	excl_id <- c(220, 261, 273)
}

# Pre-process other_order 
other_k <- other_k %>% 
	as_tibble() %>%
	mutate(id = as.numeric(id), 
				 other1_k = as.numeric(other1_k), 
				 other2_k = as.numeric(other2_k), 
				 other1_pref = as.factor(other1_pref), 
				 other2_pref = as.factor(other2_pref)) %>%
	mutate(other_i_k = if_else(other1_pref == "i", other1_k, other2_k),
				 other_p_k = if_else(other1_pref == "p", other1_k, other2_k))

# Exclude participants 
excl_condition <- !pt[1, "id", ] %in% excl_id 
pt <- pt[,, excl_condition]
other_k <- other_k %>% 
	filter(!id %in% excl_id)

# Extract metadata
sz <- dim(pt) # sz = [# of trials, # of experiment parameters, # of participants]
n_trials <- sz[1] # number of trials 
n_pt <- sz[3] # number of participants

# Reshape data 
id_vec <- vector("numeric", length = n_pt)
id_vec <- pt[1, "id", 1:n_pt]
names(id_vec) <- NULL

if (pop == 1) {
	id_vec_ya <- id_vec
} else if (pop == 2) {
	id_vec_oa <- id_vec
}

pt_new <- array(0, dim = sz)
col_names <- c("id", "s_dec_1", "s_rs_1", "s_rl_1", "s_dl_1", 
							 "so_dec_i", "o_rs_i", "o_rl_i", "o_dl_i", "o_dec_i", 
							 "s_dec_i", "s_rs_i", "s_rl_i", "s_dl_i", 
							 "so_dec_p", "o_rs_p", "o_rl_p", "o_dl_p", "o_dec_p", 
							 "s_dec_p", "s_rs_p", "s_rl_p", "s_dl_p")
colnames(pt_new) <- col_names

# Reshape data 
for (i in 1:n_pt) {
	
	other1_pref <- other_k$other1_pref[i]
	
	if (other1_pref == "i") {
		pt_new[, 1:5, i] <- pt[, 1:5, i] # from id to s_dl_1
		pt_new[, 6:14, i] <- pt[, 6:14, i] # from so_dec_i to s_dl_i 
		pt_new[, 15:23, i] <- pt[, 15:23, i] # from so_dec_p to s_dl_p 
	} else if (other1_pref == "p") {
		pt_new[, 1:5, i] <- pt[, 1:5, i] # from id to s_dl_1 
		pt_new[, 6:14, i] <- pt[, 15:23, i] # from so_dec_i to s_dl_i 
		pt_new[, 15:23, i] <- pt[, 6:14, i] # from so_dec_p to s_dl_p 
	}
} 

# ============================================================================.
# Set Stan parameters ####
# ============================================================================.
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_iter     <- 4000
n_chains   <- 4 
n_warmup   <- floor(n_iter / 2)
n_thin     <- 1
n_seed     <- 1213

# ============================================================================.
# Data simulation #### 
# ============================================================================.
## Prepare experimental parameters ####
# Impulsive block 
i_reward_sooner_sim <- array(0, dim = c(n_pt, n_trials))
i_reward_later_sim <- array(0, dim = c(n_pt, n_trials))
i_delay_later_sim <- array(0, dim = c(n_pt, n_trials))

# Patient block 
p_reward_sooner_sim <- array(0, dim = c(n_pt, n_trials))
p_reward_later_sim <- array(0, dim = c(n_pt, n_trials))
p_delay_later_sim <- array(0, dim = c(n_pt, n_trials))

for (i in 1:n_pt) {
	i_reward_sooner_sim[i, ] <- pt_new[, "o_rs_i", i]
	i_reward_later_sim[i, ] <- pt_new[, "o_rl_i", i]
	i_delay_later_sim[i, ] <- pt_new[, "o_dl_i", i]
	
	p_reward_sooner_sim[i, ] <- pt_new[, "o_rs_p", i]
	p_reward_later_sim[i, ] <- pt_new[, "o_rl_p", i]
	p_delay_later_sim[i, ] <- pt_new[, "o_dl_p", i]
}

## Get empirical posterior values 
km_i <- array(0, n_pt)
ku_i <- array(0, n_pt)

km_p <- array(0, n_pt)
ku_p <- array(0, n_pt)

if (pop == 1) {
	ku_model_other_i_ya <- 
		readRDS("data/stanfit/stanfit_ku0_ya_other_i.rds")
	ku_model_other_p_ya <- 
		readRDS("data/stanfit/stanfit_ku0_ya_other_p.rds")
	
	km_i <- colMeans(rstan::extract(ku_model_other_i_ya)$km)
	ku_i <- colMeans(rstan::extract(ku_model_other_i_ya)$ku)
	
	km_p <- colMeans(rstan::extract(ku_model_other_p_ya)$km)
	ku_p <- colMeans(rstan::extract(ku_model_other_p_ya)$ku)
	
} else if (pop == 2) {
	ku_model_other_i_oa <- 
		readRDS("data/stanfit/stanfit_ku0_oa_other_i.rds")
	ku_model_other_p_oa <- 
		readRDS("data/stanfit/stanfit_ku0_oa_other_p.rds")
	
	km_i <- colMeans(rstan::extract(ku_model_other_i_oa)$km)
	ku_i <- colMeans(rstan::extract(ku_model_other_i_oa)$ku)
	
	km_p <- colMeans(rstan::extract(ku_model_other_p_oa)$km)
	ku_p <- colMeans(rstan::extract(ku_model_other_p_oa)$ku)
}

simulated_model_file <- 
	"scripts/stan_model/ku-model_other_0_simulated.stan"

i_sim_data_list <-  list(n_trials   = n_trials,
												 n_subjects = n_pt,
												 o_rs       = i_reward_sooner_sim,
												 o_rl       = i_reward_later_sim,
												 o_dl       = i_delay_later_sim,
												 km         = km_i, 
												 ku         = ku_i)

p_sim_data_list <-  list(n_trials   = n_trials,
												 n_subjects = n_pt,
												 o_rs       = p_reward_sooner_sim,
												 o_rl       = p_reward_later_sim,
												 o_dl       = p_delay_later_sim,
												 km         = km_p, 
												 ku         = ku_p)

if (pop == 1) {
	saveRDS(i_sim_data_list, "data/i_sim_data_list_ya.RDS")
	saveRDS(p_sim_data_list, "data/p_sim_data_list_ya.RDS")
} else if (pop == 2) {
	saveRDS(i_sim_data_list, "data/i_sim_data_list_oa.RDS")
	saveRDS(p_sim_data_list, "data/p_sim_data_list_oa.RDS")
}

i_sim_data <- stan(simulated_model_file,
									 data = i_sim_data_list,
									 chains = n_chains,
									 iter   = n_iter,
									 warmup = n_warmup,
									 thin   = n_thin,
									 init   = "random", 
									 seed   = n_seed, 
									 algorithm = "Fixed_param")

p_sim_data <- stan(simulated_model_file,
									 data = p_sim_data_list,
									 chains = n_chains,
									 iter   = n_iter,
									 warmup = n_warmup,
									 thin   = n_thin,
									 init   = "random",
									 seed   = n_seed,
									 algorithm = "Fixed_param")

if (pop == 1) {
	saveRDS(i_sim_data, "data/i_sim_data_ya.RDS")
	saveRDS(p_sim_data, "data/p_sim_data_ya.RDS")
} else if (pop == 2) {
	saveRDS(i_sim_data, "data/i_sim_data_oa.RDS")
	saveRDS(p_sim_data, "data/p_sim_data_oa.RDS")
}

# ============================================================================.
# Accuracy calculation #### 
# ============================================================================.
# rstan::extract(i_sim_data)$s_dec_sim[a, b, c] 
# a: iteration; b: participant; c: trial

value_later <- array(0, dim = c(n_pt, n_trials))
model_choice <- array(0, dim = c(n_pt, n_trials)) 
acc_result <- array(0, n_pt)

## For impulsive #### 
for (j in 1:n_pt) {
	k_value <- 10^(as.numeric(other_k[j, "other_i_k"]))
	value_later[j, ] <- i_reward_later_sim[j, ] / (1 + k_value * i_delay_later_sim[j, ])
	model_choice[j, ] <- (value_later[j, ] > i_reward_sooner_sim[j, ]) + 1
	
	acc_vec <- array(0, 8000)
	s_dec_sim_matrix <- rstan::extract(i_sim_data)$s_dec_sim[, j, ]
	
	for (i in 1:8000) {
		sim_dec <- s_dec_sim_matrix[i, ]
		acc_vec[i] <- mean(sim_dec == model_choice[j, ])
		
	} # iteration loop 
	acc_result[j] = mean(acc_vec)
} # participant loop 

if (pop == 1) {
	acc_result_list[[1]] <- acc_result 
} else if (pop == 2) {
	acc_result_list[[2]] <- acc_result 
}

## For patient ####
for (j in 1:n_pt) {
	k_value <- 10^(as.numeric(other_k[j, "other_p_k"]))
	value_later[j, ] <- p_reward_later_sim[j, ] / (1 + k_value * p_delay_later_sim[j, ])
	model_choice[j, ] <- (value_later[j, ] > p_reward_sooner_sim[j, ]) + 1
	
	acc_vec <- array(0, 8000)
	s_dec_sim_matrix <- rstan::extract(p_sim_data)$s_dec_sim[, j, ]
	
	for (i in 1:8000) {
		sim_dec <- s_dec_sim_matrix[i, ]
		acc_vec[i] <- mean(sim_dec == model_choice[j, ])
		
	} # iteration loop 
	acc_result[j] = mean(acc_vec)
} # participant loop 

if (pop == 1) {
	acc_result_list[[3]] <- acc_result 
} else if (pop == 2) {
	acc_result_list[[4]] <- acc_result 
}

saveRDS(acc_result_list, "data/acc_result_list.RDS")

# ============================================================================.
# Analyses #### 
# ============================================================================.
# Convert list to tibble 
sim_acc_result <- bind_rows(
	tibble(group = "1", preference = "i", accuracy = acc_result_list[[1]], 
				 id = id_vec_ya),
	tibble(group = "1", preference = "p", accuracy = acc_result_list[[3]],
				 id = id_vec_ya),
	tibble(group = "2", preference = "i", accuracy = acc_result_list[[2]],
				 id = id_vec_oa),
	tibble(group = "2", preference = "p", accuracy = acc_result_list[[4]],
				 id = id_vec_oa)
) %>% 
	mutate(group = factor(group), 
				 preference = factor(preference), 
				 accuracy = as.numeric(accuracy)) %>%
	mutate(id = factor(id)) %>%
	mutate(accuracy = accuracy * 100) %>%
	mutate(condition = interaction(group, preference, sep = ":"))

saveRDS(sim_acc_result, "data/sim_acc_result.RDS")

# LMM 
sim_acc_result <- sim_acc_result %>%
	mutate(accuracy = accuracy / 100.0)

lmer_model <- lmer(accuracy ~ 1 + group * preference + (1|id), 
									 data = sim_acc_result)
lmer_model_stats <- stats.rlmerMod.f(lmer_model)
row.names(lmer_model_stats) = c("Intercept",
																"Group (older vs young)",
																"Others (patient vs impulsive)",
																"Group * Others")

kable(lmer_model_stats, align = rep("c", 6),
			caption = "Linear mixed effects model for learning accuracy") %>%
	kable_styling()

## Bayes factor
### BF_01 for the main effect of age group
set.seed(1213)
bf_null <- lmBF(accuracy ~ id, data = sim_acc_result, whichRandom = "id")
bf_group <- lmBF(accuracy ~ group + id, data = sim_acc_result, 
								 whichRandom = "id")
bf01_group <- bf_null / bf_group 

### BF_01 for the interaction 
bf_interaction_null <- lmBF(accuracy ~ group + preference + id, 
														data = sim_acc_result, whichRandom = "id")
bf_interaction <- lmBF(accuracy ~ group + preference + group:preference + id,
											 data = sim_acc_result, whichRandom = "id")
bf01_interaction <- bf_interaction_null / bf_interaction 

# Plot  
vars <- "accuracy"

interaction_lims <- c("1:i", "1:p", "2:i", "2:p")
interaction_labs <- c("", " Young", "", " Older")
interaction_mean_cols <- c("#0057b7", "#ffd700", "#dbe0ff", "#f1e1bc")
interaction_data_cols <- c("#314498", "#fcb34b", "#9e9dd5", "#fecd91")

set.seed(1213)
plt_sim_acc_result <- sim_acc_result |> 
	mutate(accuracy = accuracy * 100)

plt_data <- plt_sim_acc_result %>% 
	group_by(group, preference, condition) %>% 
	filter(!is.na(get(vars))) %>%
	summarise(n = n(), 
						mean = mean(get(vars), na.rm = TRUE), 
						sd = sd(get(vars), na.rm = TRUE), 
						se = sd / sqrt(n)) %>% 
	ungroup()

# discontinuity marker 
diagonal_line <- grid::linesGrob(
	x = c(-0.03, 0.03),   
	y = c(0.5, 0),    
	gp = grid::gpar(col = "black", lwd = 1.5)
)

# Main function ---- 
plt <- ggplot(plt_data, aes(x = preference, y = mean, group = condition)) + 
	geom_jitter(data = plt_sim_acc_result,
							aes(x = preference, y = get(vars), 
									colour = condition), 
							position = position_jitterdodge(dodge.width = 0.75, 
																							jitter.height = 0, 
																							jitter.width = 0.3), 
							size = 1.6, alpha = 0.2, shape = 19, stroke = 0,
							show.legend = FALSE) + 
	geom_errorbar(aes(ymin = mean - se, 
										ymax = mean + se), 
								position = position_dodge(0.75), 
								colour = "black", width = 0, size = 0.25) + 
	geom_jitter(aes(y = mean, 
									fill = condition), 
							position = position_dodge(0.75), 
							colour = "black", shape = 21, stroke = 0.25, size = 2) + 
	scale_x_discrete(labels = c("Impulsive", "Patient")) +
	scale_y_continuous(limits = c(40, 105),
										 breaks = seq(40, 100, by = 10), 
										 labels = c(0, seq(50, 100, by = 10)),
										 expand = c(0, 0)) + 
	labs(x = NULL, y = "Simulated learning accuracy (%)") +
	scale_fill_manual(name = "condition", 
										limits = interaction_lims, 
										values = interaction_mean_cols, 
										labels = interaction_labs) +
	scale_colour_manual(name = "condition", 
											limits = interaction_lims, 
											values = interaction_data_cols, 
											labels = interaction_labs) +
	theme_classic() + 
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(2, 2), colour = NA),
				axis.text.x = element_text(size = 12, colour = "black"),
				axis.text.y = element_text(size = 10), 
				axis.title.x = element_text(size = 10),
				axis.title.y = element_text(size = 12), 
				legend.text = element_text(size = 8, margin = margin(l = 3)),
				legend.key.spacing.x = unit(-0.1, "cm"),
				legend.key.spacing.y = unit(-0.1, "cm"),
				legend.key.height = unit(0.4, "cm"),
				legend.key.width = unit(0.4, "cm"), 
				legend.box.spacing = unit(0, "cm"), 
				legend.position = c(0.5, 0.25)) +
	guides(fill = guide_legend(ncol = 2, byrow = TRUE, 
														 override.aes = list(size = 2)), 
				 colour = guide_legend(ncol = 2, byrow = TRUE)) +
	geom_signif(comparisons = list(c("i", "p")),
							y_position = 95, tip_length = 0, vjust = 0.5, size = 0.25,
							textsize = 3.5, family = "Arial",
							annotation = "***") +
	geom_signif(xmin = c(0.85, 1.85), xmax = c(1.15, 2.15),  
							y_position = c(89, 89), tip_length = 0, vjust = 0, size = 0.25,
							textsize = 3, family = "Arial",  
							annotation = "ns") + 
	annotate("point", x = 0.9, y = 83.35, colour = "red", alpha = 0.5, size = 1) +
	annotate("point", x = 1.27, y = 82.03, colour = "red", alpha = 0.5, size = 1) +
	annotate("point", x = 1.9, y = 86.39, colour = "red", alpha = 0.5, size = 1) +
	annotate("point", x = 2.27, y = 84.79, colour = "red", alpha = 0.5, size = 1) +
	geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.3) + 
	annotation_custom(diagonal_line, xmin = -Inf, xmax = Inf, ymin = 43, ymax = 45) +
	annotation_custom(diagonal_line, xmin = -Inf, xmax = Inf, ymin = 45, ymax = 47) +
	coord_cartesian(clip = "off")

ggsave("data/plots/acc_sim.png", plt,
			 height = 4, width = 4 * 0.65, dpi = 1200) 