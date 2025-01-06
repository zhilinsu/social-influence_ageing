# ============================================================================.
# Info ####
# ============================================================================.
# The main analysis script generates plots and performs statistical analyses.  
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Social influence and ageing project.  
# Last revision: 21 August 2024, Zhilin Su 
#
# ============================================================================.
# 0-1. Preparation ####
# ============================================================================.
# Clear the workspace
rm(list = ls())

# Load libraries
## Statistics and modelling 
library(rstatix) # pipe-friendly statistics 
library(lme4) # linear mixed-effects modelling 
library(Hmisc) # rcorr correlations
library(BayesFactor) # Bayes' factor 
library(psych) # corr.test
library(cocor) # correlations comparisons 
library(BFpack) # Bayes' factor 

## Visualisation 
library(ggplot2)
library(ggpubr)
library(ggtext) 
library(ggpattern)
library(ggsignif)
library(grid) 
library(gridExtra)
library(patchwork)

## Wrangling 
library(dplyr)
library(readr)
library(forcats)
library(abind)
library(stringr)

## Others
library(finalfit) 
library(kableExtra)
#library(mixedpower) # simulation-based power analysis 
source("scripts/helper_functions/helper_functions.R")
source("scripts/helper_functions/rankBasedCommonFunctions.R") # Bayes' factor
source("scripts/helper_functions/spearmanSampler.R") # Bayes' factor 

# General settings of plotting 
resolution <- 1200
plot_h <- 4
ax_text <- 10
ax_title <- 12

# ============================================================================.
# 0-2. Exclusion ####
# ============================================================================.
excl_id <- c(110, 111, 118, 158, 220, 261, 273)
## 110, 111: previous study of psychology 
## 158: diagnosis of ADHD
## 118, 220, 273: missing task data 
## 261: potential diagnosis of dementia 

excl_id_qnr <- c(106, 266)
## 106, 266: missing data

# ============================================================================.
# 0-3. Data pre-processing ####
# ============================================================================.
load("data/sd_data.RData")

pt_1 <- sd_ya
pt_2 <- sd_oa
other_k_1 <- k_ya 
other_k_2 <- k_oa 

pt <- abind(pt_1, pt_2, along = 3)
other_k <- rbind(other_k_1, other_k_2)

rm(sd_ya, sd_oa, pt_1, pt_2, k_ya, k_oa, other_k_1, other_k_2)

# Convert data type
other_k <- other_k |> 
	as_tibble() |>
	mutate(id = as.numeric(id), 
				 other1_k = as.numeric(other1_k), 
				 other2_k = as.numeric(other2_k), 
				 other1_pref = as.factor(other1_pref), 
				 other2_pref = as.factor(other2_pref))

# Pre-process other_order 
other_k <- other_k |>
	mutate(other_i_k = if_else(other1_pref == "i", other1_k, other2_k),
				 other_p_k = if_else(other1_pref == "p", other1_k, other2_k))

# Exclude participants 
excl_condition <- !pt[1, "id", ] %in% excl_id 
pt <- pt[,, excl_condition]
other_k <- other_k |> 
	filter(!id %in% excl_id)

# Extract metadata
sz <- dim(pt) # sz = [# of trials, # of experiment parameters, # of participants]
n_trials <- sz[1] # number of trials 
n_pt <- sz[3] # number of participants

# Reshape data 
id_vec <- vector("numeric", length = n_pt)
id_vec <- pt[1, "id", 1:n_pt]
names(id_vec) <- NULL

pt_new <- array(0, dim = sz)
col_names <- c("id", "s_dec_1", "s_rs_1", "s_rl_1", "s_dl_1", 
							 "so_dec_i", "o_rs_i", "o_rl_i", "o_dl_i", "o_dec_i", 
							 "s_dec_i", "s_rs_i", "s_rl_i", "s_dl_i", 
							 "so_dec_p", "o_rs_p", "o_rl_p", "o_dl_p", "o_dec_p", 
							 "s_dec_p", "s_rs_p", "s_rl_p", "s_dl_p")
colnames(pt_new) <- col_names

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

# Calculate accuracy  
col_names <- c("id", "acc_i", "acc_p")
acc_mat <- matrix(nrow = n_pt, ncol = length(col_names))
colnames(acc_mat) <- col_names
acc_tb <- as_tibble(acc_mat)
acc_tb$id <- id_vec

rm(acc_mat)

self_decision <- array(0, dim = c(n_pt, n_trials)) 
reward_sooner <- array(0, dim = c(n_pt, n_trials))
reward_later <- array(0, dim = c(n_pt, n_trials))
delay_later <- array(0, dim = c(n_pt, n_trials))

value_later <- array(0, dim = c(n_pt, n_trials)) # subjective value of later offer
model_choice <- array(0, dim = c(n_pt, n_trials)) 

for (o_block in c("i", "p")) {
	other_k_vec <- other_k[, paste0("other_", o_block, "_k")]
	
	for (i_id in 1:n_pt) {
		self_decision[i_id, ] <- pt_new[, paste0("so_dec_", o_block), i_id] 
		reward_sooner[i_id, ] <- pt_new[, paste0("o_rs_", o_block), i_id]
		reward_later[i_id, ] <- pt_new[, paste0("o_rl_", o_block), i_id]
		delay_later[i_id, ] <- pt_new[, paste0("o_dl_", o_block), i_id] 
		
		k_value <- as.numeric(10^(other_k_vec[i_id, 1]))
		
		value_later[i_id, ] <- reward_later[i_id, ] / (1 + k_value * delay_later[i_id, ])
		
		model_choice[i_id, ] <- (value_later[i_id, ] > reward_sooner[i_id, ]) + 1
		
		acc <- mean(self_decision[i_id, ] == model_choice[i_id, ])
		acc_tb[i_id, paste0("acc_", o_block)] <- acc
	}
}

# Questionnaires
qnr <- read_csv("data/questionnaires.csv") |>
	rename(id = ID_Code, group = Group, age = Age, gender = Gender, 
				 ami_ba = AMI_BA, # behavioural activation 
				 ami_sm = AMI_SM, # social motivation 
				 ami_es = AMI_ES, # emotional sensitivity 
				 srp_interpers = SRP_interpersonal, # interpersonal manipulation 
				 srp_affect = SRP_affective, # callous affect 
				 srp_life = SRP_lifestyle, # erratic lifestyle 
				 srp_antisoc = SRP_antisocial, # anti-social behaviour 
				 qcae_affect = QCAE_AE, # affective empathy 
				 qcae_cogn = QCAE_CE, # cognitive empathy (mentalising)
				 tas_ddf = TAS_ddf, # difficulty describing feelings 
				 tas_dif = TAS_dif, # difficulty identifying feelings 
				 tas_eot = TAS_eot, # externally-oriented thinking 
				 aq_sk = AQ_SK, # social skills 
				 aq_as = AQ_AS, # attention switching 
				 aq_ad = AQ_AD, # attention to detail 
				 aq_cm = AQ_CM, # communication 
				 aq_im = AQ_IM # imagination 
				 ) |> 
	filter(!id %in% excl_id) 

# signed KL divergence 
load("data/ku0_kld_ya.RData")
kld_ya <- result_tb
load("data/ku0_kld_oa.RData")
kld_oa <- result_tb 

kld <- rbind(kld_ya, kld_oa) |> 
  filter(!id %in% excl_id)

rm(result_tb, kld_ya, kld_oa)

# Merging 
data_all <- qnr |>
  left_join(acc_tb, by = "id") |>
  left_join(kld, by = "id") |>
  left_join(other_k, by = "id") |>
  mutate(group = as.factor(group))

# Handle exceptions: when subjects had no Others of different preferences 
## Data were analysed for Others with the greatest distance between
## the subject’s own discount rate, and the Other's discount rate estimated by subject.

for (i_id in 1:n_pt) {
  diff_i <- data_all[i_id, "self_1_km"] - data_all[i_id, "other_i_km"]
  diff_p <- data_all[i_id, "self_1_km"] - data_all[i_id, "other_p_km"]
  
  # no Impulsive other
  if (check_sign(diff_i, diff_p) && diff_i > 0) { 
    if (abs(diff_i) > abs(diff_p)) {
      data_all[i_id, "acc_p"] <- data_all[i_id, "acc_i"]
      data_all[i_id, "acc_i"] <- NA
    } else if (abs(diff_i) < abs(diff_p)) {
      data_all[i_id, "acc_i"] <- NA
    }
  }
  
  # no Patient other 
  if (check_sign(diff_i, diff_p) && diff_i < 0) { 
    if (abs(diff_i) > abs(diff_p)) {
      data_all[i_id, "acc_p"] <- NA
    } else if (abs(diff_i) < abs(diff_p)) {
      data_all[i_id, "acc_i"] <- data_all[i_id, "acc_p"]
      data_all[i_id, "acc_p"] <- NA
    }
  }
}

# Save results 
saveRDS(data_all, 
				file = "data/data_all.RData")
write_csv(data_all, 
					file = "data/data_all.csv")

# ============================================================================.
# 1. Demographics ####
# ============================================================================.
# Sample size
demogr_sample_size <- calculate_sample_size(data_all, "id", "group")

# Age
demogr_age <- get_summary(data_all, "age", "group")

# Gender
demogr_female <- calculate_females(data_all, "group")
gender_test <- chisq_test(data_all$gender, data_all$group)

# Education years 
education_summary <- get_summary(data_all, "Education_years", "group")
education_norm <- check_normality(data_all, "Education_years", "group")

## Independent two-sample Wilcoxon t-test 
wilcox_education <- data_all |> 
	wilcox_test(Education_years ~ group, paired = FALSE)

wilcox_education_z <- qnorm(wilcox_education$p / 2)

wilcox_education_es <- data_all |>
	wilcox_effsize(Education_years ~ group, paired = FALSE, ci = TRUE)

## Output for JASP analysis (Bayes' factor)
jasp_education <- data_all |>
	select(id, group, Education_years)
write_csv(jasp_education, 
					file = "data/jasp_education.csv")

# Standardised IQ 
iq_summary <- get_summary(data_all, "WTAW_Standard", "group")
iq_norm <- check_normality(data_all, "WTAW_Standard", "group")

## Independent two-sample Wilcoxon t-test 
wilcox_iq <- data_all |> 
	wilcox_test(WTAW_Standard ~ group, paired = FALSE)

wilcox_iq_z <- qnorm(wilcox_iq$p / 2)

wilcox_iq_es <- data_all |>
	wilcox_effsize(WTAW_Standard ~ group, paired = FALSE, ci = TRUE)

## Output for JASP analysis (Bayes' factor)
jasp_iq <- data_all |>
	select(id, group, WTAW_Standard)
write_csv(jasp_iq, 
					file = "data/jasp_iq.csv")

# ============================================================================.
# 2-1. Accuracy - descriptive ####
# ============================================================================.
# Preparation 
data_acc_1 <- data_all |> 
	select(id, group, acc_i) |> 
	rename(acc = acc_i) |> 
	mutate(other_pref = "i", 
				 other_pref = as.factor(other_pref)) |>
	mutate(condition1 = interaction(group, other_pref, sep = ":"))

data_acc_2 <- data_all |> 
	select(id, group, acc_p) |> 
	rename(acc = acc_p) |> 
	mutate(other_pref = "p", 
				 other_pref = as.factor(other_pref)) |>
	mutate(condition1 = interaction(group, other_pref, sep = ":"))

data_acc <- rbind(data_acc_1, data_acc_2) |>
	mutate(acc = acc * 100) |>
	arrange(id, other_pref)

rm(data_acc_1, data_acc_2)

# Save results for JASP 
write_csv(data_acc, 
					file = "data/jasp_acc.csv")

# Summary
acc_summary <- get_summary(data_acc, "acc", "condition1")

# Binomial test 
## Young, impulsive
acc_ya_impulsive_all <- data_acc |>
	filter(group == 1 & other_pref == "i") |> 
	filter(!is.na(acc)) |>
	count() |> 
	pull(n)
acc_ya_impulsive_above50 <- data_acc |>
	filter(group == 1 & other_pref == "i") |> 
	filter(!is.na(acc)) |> 
	filter(acc > 50) |> 
	count() |> 
	pull(n)
acc_ya_impulsive_binom <- binom.test(x = acc_ya_impulsive_above50, 
																		 n = acc_ya_impulsive_all, 
																		 p = 0.5, alternative = "greater")

## Young, patient
acc_ya_patient_all <- data_acc |>
	filter(group == 1 & other_pref == "p") |> 
	filter(!is.na(acc)) |>
	count() |> 
	pull(n)
acc_ya_patient_above50 <- data_acc |>
	filter(group == 1 & other_pref == "p") |> 
	filter(!is.na(acc)) |> 
	filter(acc > 50) |> 
	count() |> 
	pull(n)
acc_ya_patient_binom <- binom.test(x = acc_ya_patient_above50, 
																	 n = acc_ya_patient_all, 
																	 p = 0.5, alternative = "greater")

## Older, impulsive
acc_oa_impulsive_all <- data_acc |>
	filter(group == 2 & other_pref == "i") |> 
	filter(!is.na(acc)) |>
	count() |> 
	pull(n)
acc_oa_impulsive_above50 <- data_acc |>
	filter(group == 2 & other_pref == "i") |> 
	filter(!is.na(acc)) |> 
	filter(acc > 50) |> 
	count() |> 
	pull(n)
acc_oa_impulsive_binom <- binom.test(x = acc_oa_impulsive_above50, 
																		 n = acc_oa_impulsive_all, 
																		 p = 0.5, alternative = "greater")

## Older, patient 
acc_oa_patient_all <- data_acc |>
	filter(group == 2 & other_pref == "p") |> 
	filter(!is.na(acc)) |>
	count() |> 
	pull(n)
acc_oa_patient_above50 <- data_acc |>
	filter(group == 2 & other_pref == "p") |> 
	filter(!is.na(acc)) |> 
	filter(acc > 50) |> 
	count() |> 
	pull(n)
acc_oa_patient_binom <- binom.test(x = acc_oa_patient_above50, 
																	 n = acc_oa_patient_all, 
																	 p = 0.5, alternative = "greater")

# ============================================================================.
# 2-2. Accuracy - plot (Figure 1B) ####
# ============================================================================.
vars <- "acc"

interaction_lims <- c("1:i", "1:p", "2:i", "2:p")
interaction_labs <- c("", " Young", "", " Older")
interaction_alphas <- c(0.3, 0.3, 0.3, 0.3)
interaction_mean_cols <- c("#0057b7", "#ffd700", "#dbe0ff", "#f1e1bc")
interaction_data_cols <- c("#314498", "#fcb34b", "#9e9dd5", "#fecd91")

# For 95% confidence interval 
critical_value <- qnorm(0.975)

plt_data <- data_acc |> 
	group_by(group, other_pref, condition1) |> 
	filter(!is.na(get(vars))) |>
	summarise(n = n(), 
						mean = mean(get(vars), na.rm = TRUE), 
						sd = sd(get(vars), na.rm = TRUE), 
						se = sd / sqrt(n)) |> 
	mutate(ci = critical_value * se) |> 
	ungroup()

limits <- data_acc |> 
	summarise(max = max(get(vars), na.rm = TRUE), 
						min = min(get(vars), na.rm = TRUE))

set.seed(1213)

# Main function  
plt <- ggplot(plt_data, aes(x = other_pref, y = mean, group = condition1)) + 
	geom_jitter(data = data_acc,
							aes(x = other_pref, y = get(vars), 
									colour = condition1), 
							position = position_jitterdodge(dodge.width = 0.75, 
																							jitter.height = 0, 
																							jitter.width = 0.3), 
							size = 1.6, alpha = 0.2, shape = 19, stroke = 0,
							show.legend = FALSE) + 
	geom_errorbar(aes(ymin = mean - se, 
										ymax = mean + se), 
								position = position_dodge(0.75), 
								colour = "black", width = 0, linewidth = 0.25) + 
	geom_jitter(aes(y = mean, 
									fill = condition1), 
							position = position_dodge(0.75), 
							colour = "black", shape = 21, stroke = 0.25, size = 2) + 
	scale_x_discrete(labels = c("Impulsive", "Patient")) +
	scale_y_continuous(limits = c(50, 105), expand = c(0, 0)) + 
	labs(x = NULL, y = "Learning accuracy (%)") +
	scale_fill_manual(name = "condition1", 
										limits = interaction_lims, 
										values = interaction_mean_cols, 
										labels = interaction_labs) +
	scale_colour_manual(name = "condition1", 
											limits = interaction_lims, 
											values = interaction_data_cols, 
											labels = interaction_labs) +
	theme_classic() + 
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(2, 2), colour = NA),
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_text), 
				axis.title.x = element_text(size = ax_text),
				axis.title.y = element_text(size = ax_title), 
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
							annotation = "**") +
	geom_signif(xmin = c(0.85, 1.85), xmax = c(1.15, 2.15),  
							y_position = c(89, 89), tip_length = 0, vjust = 0, size = 0.25,
							textsize = 3, family = "Arial",  
							annotation = "ns")

## Discontinuity marker 
diagonal_line <- grid::linesGrob(
	x = c(-0.03, 0.03),   
	y = c(0.5, 0),    
	gp = grid::gpar(col = "black", lwd = 1.5)
)

plt <- plt + 
	scale_y_continuous(limits = c(40, 105),
										 breaks = seq(40, 100, by = 10), 
										 labels = c(0, seq(50, 100, by = 10)),
										 expand = c(0, 0)) +
	geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.3) + 
	annotation_custom(diagonal_line, xmin = -Inf, xmax = Inf, ymin = 43, ymax = 45) +
	annotation_custom(diagonal_line, xmin = -Inf, xmax = Inf, ymin = 45, ymax = 47) +
	coord_cartesian(clip = "off")

ggsave("data/plots/acc.png", plt, 
			 height = plot_h, width = plot_h * 0.75, dpi = resolution)

# ============================================================================.
# 2-3. Accuracy - LMM ####
# ============================================================================.
data_acc <- data_acc |>
	mutate(acc = acc / 100.0)

lmer_model <- lmer(acc ~ 1 + group * other_pref + (1|id), 
									 data = data_acc)

lmer_model_stats <- stats.rlmerMod.f(lmer_model)
row.names(lmer_model_stats) = c("Intercept",
																"Group (older vs young)",
																"Others (patient vs impulsive)",
																"Group * Others")

kable(lmer_model_stats, align = rep("c", 6),
			caption = "Linear mixed effects model for learning accuracy") |>
	kable_styling()

# ============================================================================.
# 2-4. Accuracy - Bayes factor ####
# ============================================================================.
data_acc <- data_acc |>
	mutate(id = as.factor(id)) |> 
	filter(!is.na(acc))

set.seed(1213)
# BF_01 for the main effect of group 
bf_null <- lmBF(acc ~ id, data = data_acc, whichRandom = "id")
bf_group <- lmBF(acc ~ group + id, data = data_acc, 
								 whichRandom = "id")
bf01_group <- bf_null / bf_group

# BF_01 for the interaction 
bf_interaction_null <- lmBF(acc ~ group + other_pref + id, data = data_acc, 
														whichRandom = "id")
bf_interaction <- lmBF(acc ~ group + other_pref + group:other_pref + id,
											 data = data_acc, whichRandom = "id")
bf01_interaction <- bf_interaction_null / bf_interaction

# ============================================================================.
# 3-1. Self-report confidence ####
# ============================================================================.
# Preparation 
data_confidence_1 <- data_all |> 
	select(id, group, SD_learn_p1, other1_pref) |> 
	rename(sd_confidence = SD_learn_p1, 
				 other_pref = other1_pref) |> 
	mutate(condition1 = interaction(group, other_pref, sep = ":"))

data_confidence_2 <- data_all |> 
	select(id, group, SD_learn_p2, other2_pref) |> 
	rename(sd_confidence = SD_learn_p2, 
				 other_pref = other2_pref) |> 
	mutate(condition1 = interaction(group, other_pref, sep = ":"))

data_confidence <- rbind(data_confidence_1, data_confidence_2) |>
	arrange(id, other_pref)

rm(data_confidence_1, data_confidence_2)

# Descriptive 
confidence_summary <- get_summary(data_confidence, "sd_confidence", "condition1")

# LMM 
lmer_model <- lmer(sd_confidence ~ 1 + group * other_pref + (1|id), 
									 data = data_confidence)

lmer_model_stats <- stats.rlmerMod.f(lmer_model)
row.names(lmer_model_stats) = c("Intercept",
																"Group (older vs young)",
																"Others (patient vs impulsive)",
																"Group * Others")

kable(lmer_model_stats, align = rep("c", 6),
			caption = "Linear mixed effects model for self-report confidence.") |>
	kable_styling()

# Bayes factor 
data_confidence <- data_confidence |>
	mutate(id = as.factor(id)) |> 
	filter(!is.na(sd_confidence))

set.seed(1213)
## BF_01 for main effect of other's preference 
bf_null <- lmBF(sd_confidence ~ id, data = data_confidence, 
								whichRandom = "id")
bf_others <- lmBF(sd_confidence ~ other_pref + id, data = data_confidence, 
									whichRandom = "id")

bf01_others <- bf_null / bf_others

## BF_01 for the interaction 
bf_interaction_null <- lmBF(sd_confidence ~ group + other_pref + id, 
														data = data_confidence, whichRandom = "id")
bf_interaction <- lmBF(sd_confidence ~ group + other_pref + group:other_pref + id,
											 data = data_confidence, whichRandom = "id")

bf01_interaction <- bf_interaction_null / bf_interaction

# ============================================================================.
# 3-2. Self-report similarity ####
# ============================================================================.
# Preparation
data_similarity_1 <- data_all |> 
	select(id, group, SD_similar_p1, other1_pref) |> 
	rename(sd_similarity = SD_similar_p1, 
				 other_pref = other1_pref) |> 
	mutate(condition1 = interaction(group, other_pref, sep = ":"))

data_similarity_2 <- data_all |> 
	select(id, group, SD_similar_p2, other2_pref) |> 
	rename(sd_similarity = SD_similar_p2, 
				 other_pref = other2_pref) |> 
	mutate(condition1 = interaction(group, other_pref, sep = ":"))

data_similarity <- rbind(data_similarity_1, data_similarity_2) |>
	arrange(id, other_pref)

rm(data_similarity_1, data_similarity_2)

# Descriptive 
similarity_summary <- get_summary(data_similarity, "sd_similarity", "condition1")

# LMM 
lmer_model <- lmer(sd_similarity ~ 1 + group * other_pref + (1|id), 
									 data = data_similarity)

lmer_model_stats <- stats.rlmerMod.f(lmer_model)
row.names(lmer_model_stats) = c("Intercept",
																"Group (older vs young)",
																"Others (patient vs impulsive)",
																"Group * Others")

kable(lmer_model_stats, align = rep("c", 6),
			caption = "Linear mixed effects model for self-report similarity.") |>
	kable_styling()

# Bayes factor
data_similarity <- data_similarity |>
	mutate(id = as.factor(id)) |> 
	filter(!is.na(sd_similarity))

set.seed(1213)
## BF_01 for main effect of group
bf_null <- lmBF(sd_similarity ~ id, data = data_similarity, 
								whichRandom = "id")
bf_group <- lmBF(sd_similarity ~ group + id, data = data_similarity, 
								 whichRandom = "id")

bf01_group <- bf_null / bf_group

## BF_01 for the interaction 
bf_interaction_null <- lmBF(sd_confidence ~ group + other_pref + id, 
														data = data_confidence, whichRandom = "id")
bf_interaction <- lmBF(sd_confidence ~ group + other_pref + group:other_pref + id,
											 data = data_confidence, whichRandom = "id")

bf01_interaction <- bf_interaction_null / bf_interaction

# ============================================================================.
# 4-1. Mean of discounting distribution (km) ####
# ============================================================================.
# Preparation 
data_km <- data_all |> 
	select(id, group, self_1_km)

# Save results for JASP analysis (Bayes factor)
write_csv(data_km,
					file = "data/jasp_km.csv")

# Descriptive 
km_summary <- get_summary(data_km, "self_1_km", "group") 
km_norm <- check_normality(data_km, "self_1_km", "group")

# Independent two-sample Wilcoxon t-test
ind_wilcox_test <- data_km |>
	wilcox_test(self_1_km ~ group, paired = FALSE)

ind_wilcox_test_z <- qnorm(ind_wilcox_test$p / 2)

ind_wilcox_test_es <- data_km |>
	wilcox_effsize(self_1_km ~ group, paired = FALSE, 
								 ci = TRUE)

## Plot (Figure S3) ####
vars <- "self_1_km"

set.seed(1213)
critical_value <- qnorm(0.975)

plt_data <- data_km %>% 
	group_by(group) %>% 
	filter(!is.na(get(vars))) %>%
	summarise(n = n(), 
						mean = mean(get(vars), na.rm = TRUE), 
						sd = sd(get(vars), na.rm = TRUE), 
						se = sd / sqrt(n)) %>% 
	mutate(ci = critical_value * se) %>% 
	ungroup()

plt <- ggplot(plt_data, aes(x = group, y = mean)) + 
	geom_jitter(data = data_km,
							aes(x = group, y = get(vars), 
									colour = group, alpha = group), 
							position = position_jitterdodge(dodge.width = 0, 
																							jitter.height = 0, 
																							jitter.width = 0.3), 
							size = 1.3, shape = 19, stroke = 0) + 
	geom_errorbar(aes(ymin = mean - se, 
										ymax = mean + se), 
								position = position_dodge(0), 
								colour = "black", width = 0, size = 0.25) +
	geom_jitter(aes(y = mean, 
									fill = group), 
							position = position_dodge(0), 
							colour = "black", shape = 21, stroke = 0.25, size = 2) + 
	scale_x_discrete(labels = c("Young", "Older")) +
	scale_y_continuous(limits = c(-12.5, 0.5), 
										 breaks = c(-12, -10, -8, -6, -4, -2, 0), 
										 expan = expansion(mult = c(0, 0))) + 
	labs(x = NULL, y = "Mean of discounting distribution<br>(self baseline)") +
	scale_fill_manual(values = c("#00a087", "#9cddcd"), 
										labels = c("Young", "Older")) +
	scale_colour_manual(values = c("#00a087", "#9cddcd"),
											labels = c("Young", "Older")) +
	scale_alpha_manual(values = c(0.3, 0.5)) + 
	theme_classic() + 
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(2, 2), colour = NA),
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_text), 
				axis.title.x = element_text(size = ax_text), 
				axis.title.y = element_markdown(size = ax_title, lineheight = 1.2), 
				legend.position = "none") +
	geom_hline(yintercept = 0, linetype = "dashed") + 
	geom_signif(comparisons = list(c("1", "2")),
							y_position = -2.2, tip_length = 0, vjust = 0, size = 0.25, 
							textsize = 3, family = "Arial", 
							annotation = "ns")

ggsave("data/plots/km.png", plt, 
			 height = plot_h, width = plot_h * 0.6, dpi = resolution)

# ============================================================================.
# 4-2. Standard deviation of discounting distribution (ku) ####
# ============================================================================.
# Preparation 
data_ku <- data_all |> 
	select(id, group, self_1_ku)

# Save results for JASP analysis (Bayes factor)
write_csv(data_ku,
					file = "data/jasp_ku.csv")

# Descriptive 
ku_summary <- get_summary(data_ku, "self_1_ku", "group")
ku_norm <- check_normality(data_ku, "self_1_ku", "group")

# Independent two-sample Wilcoxon t-test
ind_wilcox_test <- data_ku |>
	wilcox_test(self_1_ku ~ group, paired = FALSE)

ind_wilcox_test_z <- qnorm(ind_wilcox_test$p/2)

ind_wilcox_test_es <- data_ku |>
	wilcox_effsize(self_1_ku ~ group, paired = FALSE, 
								 ci = TRUE)

## Plot (Figure S3) ####
vars <- "self_1_ku"

set.seed(1213)
critical_value <- qnorm(0.975)

plt_data <- data_ku %>% 
	group_by(group) %>% 
	filter(!is.na(get(vars))) %>%
	summarise(n = n(), 
						mean = mean(get(vars), na.rm = TRUE), 
						sd = sd(get(vars), na.rm = TRUE), 
						se = sd / sqrt(n)) %>%
	mutate(ci = critical_value * se) %>% 
	ungroup()

limits <- data_ku %>% 
	summarise(max = max(get(vars)), 
						min = min(get(vars)))

plt <- ggplot(plt_data, aes(x = group, y = mean)) + 
	geom_jitter(data = data_ku,
							aes(x = group, y = get(vars), colour = group, alpha = group), 
							position = position_jitterdodge(dodge.width = 0, 
																							jitter.height = 0, 
																							jitter.width = 0.3), 
							size = 1.3, shape = 19, stroke = 0) + 
	geom_errorbar(aes(ymin = mean - se, 
										ymax = mean + se), 
								position = position_dodge(0), 
								colour = "black", width = 0, size = 0.25) + 
	geom_jitter(aes(y = mean,
									fill = group), 
							position = position_dodge(0), 
							colour = "black", shape = 21, stroke = 0.25, size = 2) + 
	scale_x_discrete(labels = c("Young", "Older")) +
	scale_y_continuous(limits = c(0, 4.5), expan = expansion(mult = c(0, 0))) + 
	labs(x = NULL, y = "Standard deviation of discounting distribution<br>(self baseline)") +
	scale_fill_manual(values = c("#00a087", "#9cddcd"),
										labels = c("Young", "Older")) +
	scale_colour_manual(values = c("#00a087", "#9cddcd"),
											labels = c("Young", "Older")) +
	scale_alpha_manual(values = c(0.3, 0.5)) + 
	theme_classic() + 
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(2, 2), colour = NA),
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_text), 
				axis.title.x = element_text(size = ax_text), 
				axis.title.y = element_markdown(size = ax_title, lineheight = 1.2), 
				legend.position = "none") +
	geom_hline(yintercept = 0, linetype = "dashed") + 
	geom_signif(comparisons = list(c("1", "2")),
							y_position = 3, tip_length = 0, vjust = 0, size = 0.25, 
							textsize = 3, family = "Arial", 
							annotation = "ns")

ggsave("data/plots/ku.png", plt, 
			 height = plot_h, width = plot_h * 0.6, dpi = resolution)

# ============================================================================.
# 5-1. Signed KL divergence (D_KL) - descriptive ####
# ============================================================================.
# Preparation 
data_kl_i <- data_all |> 
	dplyr::select(id, group, self_1_km, self_1_ku, norm_kl_i, WTAW_Standard) |> 
	mutate(other_pref = "i") |> 
	dplyr::rename(norm_kl = norm_kl_i)

data_kl_p <- data_all |> 
	dplyr::select(id, group, self_1_km, self_1_ku, norm_kl_p, WTAW_Standard) |> 
	mutate(other_pref = "p") |> 
	dplyr::rename(norm_kl = norm_kl_p)

data_kl <- bind_rows(data_kl_i, data_kl_p) |> 
	arrange(id, other_pref) |> 
	mutate(other_pref = as.factor(other_pref)) |> 
	mutate(condition1 = interaction(group, other_pref, sep = ":"))

rm(data_kl_i, data_kl_p)

# Self_1_km centring on the grand mean
self_1_km_grand_mean <- data_kl |> 
	summarise(mean = mean(self_1_km)) |> 
	as.numeric() 

data_kl <- data_kl |>
	mutate(centred_self_1_km = self_1_km - self_1_km_grand_mean)

# Save results for JASP 
data_kl_p <- data_kl |>
	filter(other_pref == "p")

write_csv(data_kl_p, 
					file = "data/jasp_kl_p.csv")

data_kl_older <- data_all |>
	select(id, group, norm_kl_i, norm_kl_p) |> 
	filter(group == 2) 

write_csv(data_kl_older, 
					file = "data/jasp_kl_older.csv")

# Summary
kl_summary <- get_summary(data_kl, "norm_kl", "condition1")
kl_norm <- check_normality(data_kl, "norm_kl", "condition1")

# One-sample nonparametric t-tests 
## Impulsive 
ind_wilcox_test_i <- data_kl |>
	filter(other_pref == "i") |> 
	wilcox_test(norm_kl ~ 1, mu = 0)

ind_wilcox_test_i_z <- qnorm(ind_wilcox_test_i$p/2)

ind_wilcox_test_i_es <- data_kl |>
	filter(other_pref == "i") |>
	wilcox_effsize(norm_kl ~ 1, paired = FALSE, 
								 ci = TRUE)

## Patient 
ind_wilcox_test_p <- data_kl |>
	filter(other_pref == "p") |> 
	wilcox_test(norm_kl ~ 1, mu = 0)

ind_wilcox_test_p_z <- qnorm(ind_wilcox_test_p$p/2)

ind_wilcox_test_p_es <- data_kl |>
	filter(other_pref == "p") |>
	wilcox_effsize(norm_kl ~ 1, paired = FALSE, 
								 ci = TRUE)

# Paired nonparametric t-test 
ind_wilcox_test <- data_kl |> 
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_test(norm_kl ~ other_pref, paired = TRUE)

ind_wilcox_test_z <- qnorm(ind_wilcox_test$p / 2)

ind_wilcox_test_es <- data_kl |> 
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_effsize(norm_kl ~ other_pref, paired = TRUE, ci = TRUE)

# ============================================================================.
# 5-2. Signed D_KL - plot (Figure 2A) ####
# ============================================================================.
vars <- "norm_kl"

interaction_lims <- c("1:i", "1:p", "2:i", "2:p")
interaction_labs <- c("", " Young", "", " Older")
interaction_bar_cols <- c("#0057b7", "#ffd700", "#dbe0ff", "#f1e1bc")
interaction_data_cols <- c("#314498", "#fcb34b", "#8d8cbe", "#fecd91")

# For 95% confidence interval
critical_value <- qnorm(0.975)

plt_data <- data_kl |> 
	group_by(group, other_pref, condition1) |> 
	filter(!is.na(get(vars))) |>
	dplyr::summarise(n = n(), 
									 mean = mean(get(vars), na.rm = TRUE), 
									 sd = sd(get(vars), na.rm = TRUE), 
									 se = sd / sqrt(n)) |> 
	mutate(ci = critical_value * se) |> 
	ungroup()

limits <- data_kl |> 
	summarise(max = max(get(vars), na.rm = TRUE), 
						min = min(get(vars), na.rm = TRUE))

set.seed(1213)
# Main function
plt <- ggplot(plt_data, aes(x = other_pref, y = mean, group = condition1)) + 
	geom_bar(aes(fill = condition1), 
					 stat = "identity", 
					 position = position_dodge(0.75), 
					 width = 0.75, linewidth = 0.3, 
					 colour = "black") + 
	geom_jitter(data = data_kl,
							aes(x = other_pref, y = get(vars), colour = condition1), 
							position = position_jitterdodge(dodge.width = 0.75, 
																							jitter.height = 0, 
																							jitter.width = 0.3), 
							size = 1.3, alpha = 0.3, shape = 19, stroke = 0,   
							show.legend = FALSE) + 
	geom_errorbar(aes(ymin = mean - se, 
										ymax = mean + se), 
								position = position_dodge(0.75), 
								colour = "black", width = 0, alpha = 1, size = 0.35) + 
	scale_x_discrete(labels = c("Impulsive", "Patient")) +
	scale_y_continuous(limits = c(floor(limits$min) * 1.1, ceiling(limits$max) * 1.1), 
										 expand = c(0, 0)) + 
	labs(x = NULL, y = "Signed KL divergence") +
	scale_fill_manual(name = "condition1", 
										limits = interaction_lims, 
										values = interaction_bar_cols, 
										labels = interaction_labs) +
	scale_colour_manual(name = "condition1", 
											limits = interaction_lims, 
											values = interaction_data_cols, 
											labels = interaction_labs) +
	theme_classic() + 
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(0.5, 0.5), colour = NA),
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_title), 
				axis.title.x = element_text(size = ax_title),
				axis.title.y = element_text(size = ax_title), 
				legend.text = element_text(size = 8, margin = margin(l = 3)),
				legend.key.spacing.x = unit(-0.13, "cm"),
				legend.key.spacing.y = unit(0.2, "cm"),
				legend.key.height = unit(0.35, "cm"),
				legend.key.width = unit(0.35, "cm"), 
				legend.position = c(0.5, 0.2)) +
	guides(fill = guide_legend(ncol = 2, byrow = TRUE), 
				 colour = guide_legend(ncol = 2, byrow = TRUE)) +
	geom_hline(yintercept = 0, linetype = "dashed") + 
	geom_signif(xmin = c(0.85), xmax = c(1.15),  
							y_position = c(1.5), tip_length = 0, vjust = 0.5, size = 0.25,
							textsize = 3.5, family = "Arial",  
							annotation = c("**")) +
	geom_signif(xmin = c(1.85), xmax = c(2.15),  
							y_position = c(1.5), tip_length = 0, vjust = 0, size = 0.25,
							textsize = 3, family = "Arial",  
							annotation = c("ns")) +
	geom_signif(xmin = c(0.85), xmax = c(1.85),  
							y_position = c(3), tip_length = 0, vjust = 0.5, size = 0.25,
							textsize = 3.5, family = "Arial",  
							annotation = c("***")) +
	geom_signif(xmin = c(1.15), xmax = c(2.15),  
							y_position = c(2.25), tip_length = 0, vjust = 0, size = 0.25,
							textsize = 3, family = "Arial",  
							annotation = c("ns")) +
	annotation_custom(grob = grid::textGrob("\u2190 shift away", x = -0.23, y = 0.20, rot = 90,
																					gp = grid::gpar(fontsize = 9, fontface = "italic", fontfamily = "Arial"))) +
	annotation_custom(grob = grid::textGrob("shift toward \u2192", x = -0.23, y = 0.80, rot = 90,
																					gp = grid::gpar(fontsize = 9, fontface = "italic", fontfamily = "Arial"))) +
	coord_cartesian(expand = TRUE, clip = "off") + 
	theme(plot.margin = unit(c(0.2, 0.2, 0, 1), "cm"))

# Save the plot 
ggsave("data/plots/norm_kl.png", plt, 
			 height = plot_h, width = plot_h * 0.8, dpi = resolution)

# ============================================================================.
# 5-3. Signed D_KL - LMM ####
# ============================================================================.
# LMM 1 - without self km as covariates
lmer_model_1 <- lmer(norm_kl ~ 1 + group * other_pref + (1|id),
										 data = data_kl)

lmer_model_stats_1 <- stats.rlmerMod.f(lmer_model_1)
row.names(lmer_model_stats_1) = c("Intercept",
																	"Group (older vs young)",
																	"Others (patient vs impulsive)",
																	"Group * Others")

kable(lmer_model_stats_1, align = rep("c", 6),
			caption = "Linear mixed effects model for signed KL divergence (without self km as covariates)") |>
	kable_styling()

## Control for IQ
lmer_model_1_iq <- lmer(norm_kl ~ 1 + WTAW_Standard + group * other_pref + (1|id),
												data = data_kl)

lmer_model_stats_1_iq <- stats.rlmerMod.f(lmer_model_1_iq)
row.names(lmer_model_stats_1_iq) = c("Intercept",
																		 "IQ",
																		 "Group (older vs young)",
																		 "Others (patient vs impulsive)",
																		 "Group * Others")

kable(lmer_model_stats_1_iq, align = rep("c", 6),
			caption = "Linear mixed effects model for signed KL divergence (without self km as covariates) controlling for IQ") |>
	kable_styling()

# LMM 2
lmer_model_2 <- lmer(norm_kl ~ 1 + group * other_pref * centred_self_1_km + (1|id),
										 data = data_kl)
lmer_model_stats_2 <- stats.rlmerMod.f(lmer_model_2)
row.names(lmer_model_stats_2) = c("Intercept",
																	"Group (older vs young)", 
																	"Others (patient vs impulsive)",
																	"Centred self km",
																	"Group * Others",
																	"Group * Centred self km",
																	"Others * Centred self km",
																	"Group * Others * Centred self km")

kable(lmer_model_stats_2, align = rep("c", 6),
			caption = "Linear mixed effects model for signed KL divergence (with self km as covariates)") |>
	kable_styling()

## Control for IQ
lmer_model_2_iq <- lmer(norm_kl ~ 1 + WTAW_Standard + group * other_pref * centred_self_1_km + (1|id),
												data = data_kl)
lmer_model_stats_2_iq <- stats.rlmerMod.f(lmer_model_2_iq)
row.names(lmer_model_stats_2_iq) = c("Intercept",
																		 "IQ", 
																		 "Group (older vs young)",
																		 "Others (patient vs impulsive)",
																		 "Centred self km",
																		 "Group * Others",
																		 "Group * Centred self km",
																		 "Others * Centred self km",
																		 "Group * Others * Centred self km")

kable(lmer_model_stats_2_iq, align = rep("c", 6),
			caption = "Linear mixed effects model for signed KL divergence (with self km as covariates) controlling for IQ") |>
	kable_styling()

# ============================================================================.
# 5-4. Signed D_KL - Bayes factor ####
# ============================================================================.
data_kl <- data_kl |>
	mutate(id = as.factor(id)) |> 
	filter(!is.na(norm_kl))

set.seed(1213)
# BF_01 for the main effect of self_km 
bf_null <- lmBF(norm_kl ~ id, data = data_kl, whichRandom = "id")
bf_km <- lmBF(norm_kl ~ centred_self_1_km + id, data = data_kl, 
							whichRandom = "id")

bf01_km <- bf_null / bf_km

# BF_01 for the interaction between age and self_km
bf_interaction_full <- lmBF(norm_kl ~ group + other_pref + centred_self_1_km + 
															group:other_pref + group:centred_self_1_km + 
															other_pref:centred_self_1_km + id, 
														data = data_kl, whichRandom = "id")

bf_interaction_group_km_null <- lmBF(norm_kl ~ group + other_pref + centred_self_1_km +
																		 	group:other_pref + 
																		 	other_pref:centred_self_1_km + id,
																		 data = data_kl, whichRandom = "id")

bf01_interaction_group_km <- bf_interaction_group_km_null / bf_interaction_full

# BF_01 for the interaction between others and self_km
bf_interaction_others_km_null <- lmBF(norm_kl ~ group + other_pref + centred_self_1_km +
																				group:other_pref + group:centred_self_1_km +
																				id,
																			data = data_kl, whichRandom = "id")

bf01_interaction_age_others_km <- bf_interaction_others_km_null / bf_interaction_full

# BF_01 for the interaction between age, others ,and self_km
bf_full <- lmBF(norm_kl ~ group + other_pref + centred_self_1_km + 
									group:other_pref + group:centred_self_1_km + other_pref:centred_self_1_km + 
									group:other_pref:centred_self_1_km + id, 
								data = data_kl, whichRandom = "id")

bf01_three_way_interaction <- bf_interaction_full / bf_full

# BF_01 for the main effect of IQ 
bf_iq <- lmBF(norm_kl ~ WTAW_Standard + id, data = data_kl, whichRandom = "id")
bf_01_iq <- bf_null / bf_iq

# ============================================================================.
# 5-5. Signed D_KL - post hoc ####
# ============================================================================.
# Impulsive, young vs older 
wil_kl_impulsive_yo <- data_kl |> 
	filter(other_pref == "i") |> 
	wilcox_test(norm_kl ~ group, paired = FALSE)

z_kl_impulsive_yo <- qnorm(wil_kl_impulsive_yo$p / 2)

wil_kl_impulsive_yo_es <- data_kl |>
	filter(other_pref == "i") |>
	wilcox_effsize(norm_kl ~ group, paired = FALSE, ci = TRUE)

# Patient, young vs older 
wil_kl_patient_yo <- data_kl |> 
	filter(other_pref == "p") |> 
	wilcox_test(norm_kl ~ group, paired = FALSE)

z_kl_patient_yo <- qnorm(wil_kl_patient_yo$p / 2)

wil_kl_patient_yo_es <- data_kl |>
	filter(other_pref == "p") |>
	wilcox_effsize(norm_kl ~ group, paired = FALSE, ci = TRUE)

# Young, impulsive vs patient
wil_kl_young_ip <- data_kl |> 
	filter(group == 1) |> 
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_test(norm_kl ~ other_pref, paired = TRUE)

z_kl_young_ip <- qnorm(wil_kl_young_ip$p / 2)

wil_kl_young_ip_es <- data_kl |>
	filter(group == 1) |> 
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_effsize(norm_kl ~ other_pref, paired = TRUE, ci = TRUE)

# Older, impulsive vs patient
wil_kl_older_ip <- data_kl |> 
	filter(group == 2) |> 
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_test(norm_kl ~ other_pref, paired = TRUE)

z_kl_older_ip <- qnorm(wil_kl_older_ip$p / 2)

wil_kl_older_ip_es <- data_kl |>
	filter(group == 2) |>
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_effsize(norm_kl ~ other_pref, paired = TRUE, ci = TRUE)

# ============================================================================.
# 5-6. Simulation-based power analysis ####
# ============================================================================.
# set.seed(1213)
# 
# # Information about model used for simulation 
# model_1 <- lmer_model_1 # model used for simulation 
# model_2 <- lmer_model_2 
# 
# data <- data_kl # data used to fit the model 
# 
# fixed_effects_1 <- c("group", "other_pref") # all fixed effects specified in the model
# fixed_effects_2 <- c("group", "other_pref", "centred_self_1_km")
# 
# simvar <- "id" # which random effect do we want to vary in the simulation?
# 
# # Simulation parameters 
# steps <- 160 # sample size 
# critical_value <- 2 # which t/z value do we want to use to test for significance?
# n_sim <- 1000 # how many single simulations should be used to estimate power?
# 
# # Run simulation 
# power_kl_1 <- mixedpower(model = model_1, data = data, 
# 											   fixed_effects = fixed_effects_1, 
# 											   simvar = simvar, steps = steps,
# 											   critical_value = critical_value, 
# 											   n_sim = n_sim)
# 
# power_kl_2 <- mixedpower(model = model_2, data = data, 
# 											   fixed_effects = fixed_effects_2, 
# 											   simvar = simvar, steps = steps,
# 											   critical_value = critical_value, 
# 											   n_sim = n_sim)
# 											   
# save(power_kl_1, power_kl_2, 
# 		 file = "data/power.RData")

# ============================================================================.
# 6.1 Factor analysis ####
# ============================================================================.
# Preparation 
qnr <- qnr |> 
	select(id, ami_ba, ami_sm, ami_es, srp_interpers, srp_affect, srp_life, srp_antisoc, 
				 qcae_affect, qcae_cogn, tas_ddf, tas_dif, tas_eot,
				 aq_sk, aq_as, aq_ad, aq_cm, aq_im) %>%
	filter(complete.cases(.))
	
# Convert the data type 
qnr_wo_id <- qnr[, 2:length(qnr)]
qnr_wo_id <- as.data.frame(qnr_wo_id)

# Determine the number of factors 
pa_result <- fa.parallel(qnr_wo_id)
vss_result <- vss(qnr_wo_id)

# Factor analysis 
fa_result_3 <- fa(r = qnr_wo_id, nfactors = 3, rotate = "oblimin", scores = "Thurstone",
									fm = "ml")

# Factor scores 
fa_scores <- fa_result_3$scores
colnames(fa_scores) <- c("autism_apathy_alexithymia",
												 "psychopathic",
												 "empathy_emotional_motivation")

fa_scores <- as_tibble(fa_scores)
fa_scores$id <- qnr$id

write_csv(fa_scores, file = "data/efa_3_factor_r.csv")

# Arrange the subscales based on their loadings 
v_order_subscales <- c("aq_cm", "aq_sk", "aq_as", "ami_sm", "tas_ddf", "tas_dif", 
											 "qcae_cogn", "tas_eot", "aq_im", "ami_ba", 
											 "srp_interpers", "srp_life", "srp_affect", "srp_antisoc", 
											 "aq_ad", "qcae_affect", "ami_es")

label_subscales <- c("(AQ) Communication", 
										 "(AQ) Social skills", 
										 "(AQ) Attention switching", 
										 "(AMI) Social motivation", 
										 "(TAS) Difficulty describing feelings", 
										 "(TAS) Difficulty identifying feelings", 
										 "(QCAE) Cognitive empathy", 
										 "(TAS) Externally-oriented thinking", 
										 "(AQ) Imagination", 
										 "(AMI) Behavioural motivation", 
										 "(SRP) Interpersonal manipulation", 
										 "(SRP) Erratic lifestyle", 
										 "(SRP) Callous affect", 
										 "(SRP) Anti-social behaviour", 
										 "(AQ) Attention to detail", 
										 "(QCAE) Affective empathy", 
										 "(AMI) Emotional sensitivity")

# ============================================================================.
# 6.2 Factor analysis - plot (Figure S6A) ####
# ============================================================================.
loadings <- unclass(fa_result_3$loadings)
row_names <- rownames(loadings)

df_loadings <- as_tibble(loadings) |> 
	mutate(subscales = row_names) |>
	mutate(which_factor = case_when(
		abs(ML1) < 0.4 & abs(ML2) < 0.4 & abs(ML3) < 0.4 ~ "ML0",
		abs(ML1) >= abs(ML2) & abs(ML1) >= abs(ML3) ~ "ML1",
		abs(ML2) >= abs(ML1) & abs(ML2) >= abs(ML3) ~ "ML2",
		abs(ML3) >= abs(ML1) & abs(ML3) >= abs(ML2) ~ "ML3"
	))

df_loadings_long <- df_loadings |>
	tidyr::pivot_longer(cols = -c(subscales, which_factor),  
							 names_to = "factor",    
							 values_to = "loadings") |> 
	mutate(abs_loadings = abs(loadings), 
				 subscales = factor(subscales),
				 which_factor = factor(which_factor), 
				 factor = factor(factor)) |>
	mutate(subscales = factor(subscales,
														levels = rev(v_order_subscales))) |> 
	mutate(condition = interaction(factor, which_factor, sep = ":"))

# Factor 1 - autism_apathy_alexithymia 
df_loadings_1 <- df_loadings_long |> 
	filter(factor == "ML1") 

plt_1 <- df_loadings_1 |>
	ggplot(aes(x = subscales, y = abs_loadings, fill = which_factor)) +
	geom_bar(stat = "identity", colour = "#fcb34b", linewidth = 0.75) + 
	scale_fill_manual(
		values = c("#fcb34b", "NA", "NA", "NA"), 
		breaks = c("ML1", "ML2", "ML3", "ML0")
	) +
	coord_flip() + 
	scale_x_discrete(breaks = rev(v_order_subscales),
									 labels = rev(label_subscales)) + 
	scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), 
										 expand = c(0, 0)) +
	labs(x = "", 
			 y = "Factor loadings") + 
	theme_classic() + 
	theme(text = element_text(family = "Arial"), 
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_text, colour = "black"),
				axis.title = element_text(size = 16),
				axis.line = element_line(linewidth = 0.7, colour = "black")) + 
	guides(fill = "none")

# Factor 2 - psychopathic traits 
df_loadings_2 <- df_loadings_long |> 
	filter(factor == "ML2") 

plt_2 <- df_loadings_2 |>
	ggplot(aes(x = subscales, y = abs_loadings, fill = which_factor)) +
	geom_bar(stat = "identity", colour = "#00a087", linewidth = 0.75) + 
	scale_fill_manual(
		values = c("NA", "#00a087", "NA", "NA"), 
		breaks = c("ML1", "ML2", "ML3", "ML0")
	) +
	coord_flip() + 
	scale_x_discrete(breaks = rev(v_order_subscales),
									 labels = rev(label_subscales)) + 
	scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), 
										 expand = c(0, 0)) +
	labs(x = "", 
			 y = "Factor loadings") + 
	theme_classic() + 
	theme(text = element_text(family = "Arial"), 
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_text, colour = "black"),
				axis.title = element_text(size = 16),
				axis.line = element_line(linewidth = 0.7, colour = "black")) + 
	guides(fill = "none")

# Factor 3 - empathy_emotional_motivation
df_loadings_3 <- df_loadings_long |> 
	filter(factor == "ML3") 

plt_3 <- df_loadings_3 |> 
	ggplot(aes(x = subscales, y = abs_loadings, fill = which_factor)) +
	geom_bar(stat = "identity", colour = "#0057b7", linewidth = 0.75) + 
	scale_fill_manual(
		values = c("NA", "NA", "#0057b7", "NA"), 
		breaks = c("ML1", "ML2", "ML3", "ML0")
	) +
	coord_flip() + 
	scale_x_discrete(breaks = rev(v_order_subscales),
									 labels = rev(label_subscales)) + 
	scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), 
										 expand = c(0, 0)) +
	labs(x = "", 
			 y = "Factor loadings") + 
	theme_classic() + 
	theme(text = element_text(family = "Arial"), 
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_text, colour = "black"),
				axis.title = element_text(size = 16),
				axis.line = element_line(linewidth = 0.7, colour = "black")) + 
	guides(fill = "none")

# Combine plots
combined_plot <- df_loadings_long |> 
	ggplot(aes(x = subscales, y = abs_loadings, fill = condition, colour = factor)) +
	geom_bar(stat = "identity", linewidth = 0.75) +
	scale_fill_manual(
		values = c("NA", "#fcb34b", "NA", "NA",
							 "NA", "NA", "#00a087", "NA",
							 "NA", "NA", "NA", "#0057b7"), 
		breaks = c("ML1:ML0", "ML1:ML1", "ML1:ML2", "ML1:ML3",
							 "ML2:ML0", "ML2:ML1", "ML2:ML2", "ML2:ML3",
							 "ML3:ML0", "ML3:ML1", "ML3:ML2", "ML3:ML3")
	) + 
	scale_colour_manual(
		values = c("#fcb34b", "#00a087", "#0057b7"), 
		breaks = c("ML1", "ML2", "ML3")
	) + 
	coord_flip() +
	scale_x_discrete(breaks = rev(v_order_subscales), labels = rev(label_subscales)) +
	scale_y_continuous(limits = c(0, 1), 
										 breaks = c(0, 0.25, 0.5, 0.75, 1), 
										 labels = c(0, 0.25, 0.5, 0.75, 1), 
										 expand = c(0, 0)) +
	labs(x = "", y = "Factor loadings") +
	theme_classic() +
	theme() +
	guides(fill = "none", colour = "none") +
	facet_grid(~factor) +
	theme(
		text = element_text(family = "Arial"), 
		axis.text.x = element_text(size = 12, colour = "black"),
		axis.text.y = element_text(size = 12, colour = "black"),
		axis.title = element_text(size = 16, face = "bold"),
		axis.line = element_line(linewidth = 0.7, colour = "black"), 
		strip.background = element_blank(),
		strip.text.x = element_blank(), 
		panel.spacing = unit(1, "lines")
	) 

ggsave("data/plots/fa_factor.png", combined_plot, 
			 height = 6, width = 6 * 1.5, dpi = resolution)

# ============================================================================.
# 6.3 Factor scores comparison ####
# ============================================================================.
fa_scores <- fa_scores |>
	mutate(group = case_when(
		id %/% 100 == 1 ~ 1, 
		id %/% 100 == 2 ~ 2
	)) |>  
	mutate(group = factor(group))

fa_scores_longer <- fa_scores |>
	tidyr::pivot_longer(cols = -c(id, group), 
							 names_to = "factor", 
							 values_to = "scores") |>
	mutate(condition = interaction(group, factor, sep = ":")) |> 
	mutate(factor = factor(factor,
												 levels = c("autism_apathy_alexithymia", 
												 					 "psychopathic", 
												 					 "empathy_emotional_motivation"))) 

# Between-group comparisons 
## Descriptive 
fc_sample_size <- calculate_sample_size(fa_scores_longer, "scores", "condition")
fc_summary <- get_summary(fa_scores_longer, "scores", "condition")
fc_outliers <- detect_outliers(fa_scores_longer, "scores", "condition")
fc_norm <- check_normality(fa_scores_longer, "scores", "condition")

## Wilcox tests
### Factor 1 - Autism & alexithymia 
wil_fc1 <- fa_scores_longer |>
	filter(factor == "autism_apathy_alexithymia") |>
	wilcox_test(scores ~ group)

z_fc1 <- qnorm(wil_fc1$p / 2)

wil_fc1_es <- fa_scores_longer |> 
	filter(factor == "autism_apathy_alexithymia") |>
	wilcox_effsize(scores ~ group, ci = TRUE)

fa_scores_longer_fc1 <- fa_scores_longer |> 
	filter(factor == "autism_apathy_alexithymia")
write_csv(fa_scores_longer_fc1, 
					file = "data/jasp_fc1.csv") # for Bayes' factor 

### Factor 2 - Psychopathic traits 
wil_fc2 <- fa_scores_longer |>
	filter(factor == "psychopathic") |>
	wilcox_test(scores ~ group)

z_fc2 <- qnorm(wil_fc2$p / 2)

wil_fc2_es <- fa_scores_longer |> 
	filter(factor == "psychopathic") |>
	wilcox_effsize(scores ~ group, ci = TRUE)

### Factor 3 - Affective empathy & motivation 
wil_fc3 <- fa_scores_longer |>
	filter(factor == "empathy_emotional_motivation") |>
	wilcox_test(scores ~ group)

z_fc3 <- qnorm(wil_fc3$p / 2)

wil_fc3_es <- fa_scores_longer |> 
	filter(factor == "empathy_emotional_motivation") |>
	wilcox_effsize(scores ~ group, ci = TRUE)

# ============================================================================.
# 6.4 Factor scores - correlations ####
# ============================================================================.
data_qnr_factor <- fa_scores |> 
	mutate(id = as.factor(id))

data_qnr_factor_kl <- data_qnr_factor |> 
	left_join(data_kl, by = c("id", "group")) |>  
	select(autism_apathy_alexithymia, psychopathic, empathy_emotional_motivation, 
				 group, norm_kl, other_pref, condition1)

data_factor_kl_i_ya <- data_qnr_factor_kl |>
	filter(group == "1", other_pref == "i")
data_factor_kl_i_oa <- data_qnr_factor_kl |>
	filter(group == "2", other_pref == "i")
data_factor_kl_p_ya <- data_qnr_factor_kl |>
	filter(group == "1", other_pref == "p")
data_factor_kl_p_oa <- data_qnr_factor_kl |>
	filter(group == "2", other_pref == "p")

# Factor 1 - autistic & alexithymic traits
## Young, impulsive 
corr_factor1_kl_i_ya <- rcorr(data_factor_kl_i_ya$norm_kl, 
															data_factor_kl_i_ya$autism_apathy_alexithymia, 
															type = "spearman")

corr_factor1_kl_i_ya_ci <- corr.test(data_factor_kl_i_ya$norm_kl,
																		 data_factor_kl_i_ya$autism_apathy_alexithymia,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_i_ya$norm_kl,
																		yVals = data_factor_kl_i_ya$autism_apathy_alexithymia)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Older, impulsive 
corr_factor1_kl_i_oa <- rcorr(data_factor_kl_i_oa$norm_kl, 
															data_factor_kl_i_oa$autism_apathy_alexithymia, 
															type = "spearman")

corr_factor1_kl_i_oa_ci <- corr.test(data_factor_kl_i_oa$norm_kl,
																		 data_factor_kl_i_oa$autism_apathy_alexithymia,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_i_oa$norm_kl,
																		yVals = data_factor_kl_i_oa$autism_apathy_alexithymia)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Young, patient 
corr_factor1_kl_p_ya <- rcorr(data_factor_kl_p_ya$norm_kl, 
															data_factor_kl_p_ya$autism_apathy_alexithymia, 
															type = "spearman")

corr_factor1_kl_p_ya_ci <- corr.test(data_factor_kl_p_ya$norm_kl,
																		 data_factor_kl_p_ya$autism_apathy_alexithymia,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_p_ya$norm_kl,
																		yVals = data_factor_kl_p_ya$autism_apathy_alexithymia)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Older, impulsive 
corr_factor1_kl_p_oa <- rcorr(data_factor_kl_p_oa$norm_kl, 
															data_factor_kl_p_oa$autism_apathy_alexithymia, 
															type = "spearman")

corr_factor1_kl_p_oa_ci <- corr.test(data_factor_kl_p_oa$norm_kl,
																		 data_factor_kl_p_oa$autism_apathy_alexithymia,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_p_oa$norm_kl,
																		yVals = data_factor_kl_p_oa$autism_apathy_alexithymia)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

# Factor 2 - psychopathic traits 
## Young, impulsive 
corr_factor2_kl_i_ya <- rcorr(data_factor_kl_i_ya$norm_kl, 
															data_factor_kl_i_ya$psychopathic, 
															type = "spearman")

corr_factor2_kl_i_ya_ci <- corr.test(data_factor_kl_i_ya$norm_kl,
																		 data_factor_kl_i_ya$psychopathic,
																		 method = "spearman")
 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_i_ya$norm_kl,
																		yVals = data_factor_kl_i_ya$psychopathic)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Older, impulsive 
corr_factor2_kl_i_oa <- rcorr(data_factor_kl_i_oa$norm_kl, 
															data_factor_kl_i_oa$psychopathic, 
															type = "spearman")

corr_factor2_kl_i_oa_ci <- corr.test(data_factor_kl_i_oa$norm_kl,
																		 data_factor_kl_i_oa$psychopathic,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_i_oa$norm_kl,
																		yVals = data_factor_kl_i_oa$psychopathic)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Young, patient 
corr_factor2_kl_p_ya <- rcorr(data_factor_kl_p_ya$norm_kl, 
															data_factor_kl_p_ya$psychopathic, 
															type = "spearman")

corr_factor2_kl_p_ya_ci <- corr.test(data_factor_kl_p_ya$norm_kl,
																		 data_factor_kl_p_ya$psychopathic,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_p_ya$norm_kl,
																		yVals = data_factor_kl_p_ya$psychopathic)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Older, patient 
corr_factor2_kl_p_oa <- rcorr(data_factor_kl_p_oa$norm_kl, 
															data_factor_kl_p_oa$psychopathic, 
															type = "spearman")

corr_factor2_kl_p_oa_ci <- corr.test(data_factor_kl_p_oa$norm_kl,
																		 data_factor_kl_p_oa$psychopathic,
																		 method = "spearman")
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_p_oa$norm_kl,
																		yVals = data_factor_kl_p_oa$psychopathic)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

# Factor 3 - emotional sensitivity, affective empathy
## Young, impulsive 
corr_factor3_kl_i_ya <- rcorr(data_factor_kl_i_ya$norm_kl, 
															data_factor_kl_i_ya$empathy_emotional_motivation, 
															type = "spearman")

corr_factor3_kl_i_ya_ci <- corr.test(data_factor_kl_i_ya$norm_kl,
																		 data_factor_kl_i_ya$empathy_emotional_motivation,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_i_ya$norm_kl,
																		yVals = data_factor_kl_i_ya$empathy_emotional_motivation)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Older, impulsive
corr_factor3_kl_i_oa <- rcorr(data_factor_kl_i_oa$norm_kl, 
															data_factor_kl_i_oa$empathy_emotional_motivation, 
															type = "spearman")

corr_factor3_kl_i_oa_ci <- corr.test(data_factor_kl_i_oa$norm_kl,
																		 data_factor_kl_i_oa$empathy_emotional_motivation,
																		 method = "spearman")
 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_i_oa$norm_kl,
																		yVals = data_factor_kl_i_oa$empathy_emotional_motivation)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Young, patient 
corr_factor3_kl_p_ya <- rcorr(data_factor_kl_p_ya$norm_kl, 
															data_factor_kl_p_ya$empathy_emotional_motivation, 
															type = "spearman")

corr_factor3_kl_p_ya_ci <- corr.test(data_factor_kl_p_ya$norm_kl,
																		 data_factor_kl_p_ya$empathy_emotional_motivation,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_p_ya$norm_kl,
																		yVals = data_factor_kl_p_ya$empathy_emotional_motivation)
bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

## Older, patient 
corr_factor3_kl_p_oa <- rcorr(data_factor_kl_p_oa$norm_kl, 
															data_factor_kl_p_oa$empathy_emotional_motivation, 
															type = "spearman")

corr_factor3_kl_p_oa_ci <- corr.test(data_factor_kl_p_oa$norm_kl,
																		 data_factor_kl_p_oa$empathy_emotional_motivation,
																		 method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_factor_kl_p_oa$norm_kl,
																		yVals = data_factor_kl_p_oa$empathy_emotional_motivation)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

# Correction for multiple comparisons 
## Older, impulsive 
factor_kl_i_oa <- c(corr_factor1_kl_i_oa$P[1, 2],
										corr_factor2_kl_i_oa$P[1, 2],
										corr_factor3_kl_i_oa$P[1, 2])

factor_kl_i_oa_adjusted <- p.adjust(factor_kl_i_oa, method = "fdr")

## Older, patient 
factor_kl_p_oa <- c(corr_factor1_kl_p_oa$P[1, 2],
										corr_factor2_kl_p_oa$P[1, 2],
										corr_factor3_kl_p_oa$P[1, 2])

factor_kl_p_oa_adjusted <- p.adjust(factor_kl_p_oa, method = "fdr")

# ============================================================================.
# 6.5 Factor scores - correlational plot (Figures 2B, S6B) ####
# ============================================================================.
# Factor 3 & older impulsive 
plt <- ggplot(data_factor_kl_i_oa, aes(x = empathy_emotional_motivation, y = norm_kl)) + 
	geom_point(colour = "#98a6e2", size = 2.5, alpha = 0.5, shape = 19, stroke = 0, 
						 show.legend = FALSE) +
	geom_smooth(colour = "#98a6e2", fill = "#98a6e2", alpha = 0.15, size = 0.75,
							method = "lm", formula = y ~ x, show.legend = FALSE) + 
	scale_x_continuous(name = paste0("Affective empathy & emotional motivation")) +
	scale_y_continuous(name = paste0("Impulsive", " <i>D</i>", "<i><sub>KL</sub></i>", " (Older)")) +
	theme_classic() +
	theme(text = element_text(family = "Arial"),
				axis.text = element_text(size = ax_text),
				axis.title.x = element_text(size = ax_title),
				axis.title.y = element_markdown(size = ax_title)) +
	annotate(geom = "richtext", x = -1, y = -2,
					 label = paste0("_r_", "<sub>s</sub>(71)=",
					 							 round(corr_factor3_kl_i_oa$r[1,2], 3), 
					 							 ", _p_=",
					 							 round(corr_factor3_kl_i_oa$P[1,2], 3), "*"),
					 fill = NA, label.color = NA, size = 4, 
					 family = "Arial")

ggsave("data/plots/corr_fc3_kld_i_oa.png",
			 height = plot_h, width = plot_h * 1.0, dpi = resolution)

# Factor 1 & older patient 
plt <- ggplot(data_factor_kl_p_oa, aes(x = autism_apathy_alexithymia, y = norm_kl)) + 
	geom_point(colour = "#ffbf85", size = 2.5, alpha = 0.5, shape = 19, stroke = 0, 
						 show.legend = FALSE) +
	geom_smooth(colour = "#ffbf85", fill = "#ffbf85", alpha = 0.15, size = 0.75,
							method = "lm", formula = y ~ x, show.legend = FALSE) + 
	scale_x_continuous(name = paste0("Autistic & alexithymic traits")) +
	scale_y_continuous(name = paste0("Patient", " <i>D</i>", "<i><sub>KL</sub></i>", " (Older)")) +
	theme_classic() +
	theme(text = element_text(family = "Arial"),
				axis.text = element_text(size = ax_text),
				axis.title.x = element_markdown(size = ax_title),
				axis.title.y = element_markdown(size = ax_title)) +
	annotate(geom = "richtext", x = 0, y = -2,
					 label = paste0("_r_", "<sub>s</sub>(66)=",
					 							 round(corr_factor1_kl_p_oa$r[1,2], 3), 
					 							 ", _p_=",
					 							 round(corr_factor1_kl_p_oa$P[1,2], 3), "*"),
					 fill = NA, label.color = NA, size = 4, 
					 family = "Arial")

ggsave("data/plots/corr_fc1_kld_p_oa.png",
			 height = plot_h, width = plot_h * 1.0, dpi = resolution)

# ============================================================================.
# 6.6 Factor scores - between-group comparisons of correlations (Table S9) ####
# ============================================================================.
# Preparation 
df_factor_kl_i_ya <- as.data.frame(data_factor_kl_i_ya)
df_factor_kl_i_oa <- as.data.frame(data_factor_kl_i_oa)
df_factor_kl_p_ya <- as.data.frame(data_factor_kl_p_ya)
df_factor_kl_p_oa <- as.data.frame(data_factor_kl_p_oa)

# Factor 1, impulsive 
cocor_fct1_i <- cocor.indep.groups(corr_factor1_kl_i_ya$r[1, 2],
																	 corr_factor1_kl_i_oa$r[1, 2],
																	 n1 = 68, n2 = 73)

## Bayes factor 
set.seed(1213)
factor1_kl_i_ya <- df_factor_kl_i_ya[, c("autism_apathy_alexithymia", "norm_kl")]
factor1_kl_i_oa <- df_factor_kl_i_oa[, c("autism_apathy_alexithymia", "norm_kl")]

fit <- cor_test(factor1_kl_i_ya, factor1_kl_i_oa)
BF_result <- BF(fit, hypothesis = "norm_kl_with_autism_apathy_alexithymia_in_g1 = norm_kl_with_autism_apathy_alexithymia_in_g2")

# Factor 1, patient 
cocor_fct1_p <- cocor.indep.groups(corr_factor1_kl_p_ya$r[1, 2],
																	 corr_factor1_kl_p_oa$r[1, 2],
																	 n1 = 71, n2 = 68)

## Bayes factor 
set.seed(1213)
factor1_kl_p_ya <- df_factor_kl_p_ya[, c("autism_apathy_alexithymia", "norm_kl")]
factor1_kl_p_oa <- df_factor_kl_p_oa[, c("autism_apathy_alexithymia", "norm_kl")]

fit <- cor_test(factor1_kl_p_ya, factor1_kl_p_oa)
BF_result <- BF(fit, hypothesis = "norm_kl_with_autism_apathy_alexithymia_in_g1 = norm_kl_with_autism_apathy_alexithymia_in_g2")

# Factor 2, impulsive 
cocor_fct2_i <- cocor.indep.groups(corr_factor2_kl_i_ya$r[1, 2],
																	 corr_factor2_kl_i_oa$r[1, 2],
																	 n1 = 68, n2 = 73)
## Bayes factor 
set.seed(1213)
factor2_kl_i_ya <- df_factor_kl_i_ya[, c("psychopathic", "norm_kl")]
factor2_kl_i_oa <- df_factor_kl_i_oa[, c("psychopathic", "norm_kl")]

fit <- cor_test(factor2_kl_i_ya, factor2_kl_i_oa)
BF_result <- BF(fit, hypothesis = "norm_kl_with_psychopathic_in_g1 = norm_kl_with_psychopathic_in_g2")

# Factor 2, patient 
cocor_fct2_p <- cocor.indep.groups(corr_factor2_kl_p_ya$r[1, 2],
																	 corr_factor2_kl_p_oa$r[1, 2],
																	 n1 = 71, n2 = 68)
## Bayes factor 
set.seed(1213)
factor2_kl_p_ya <- df_factor_kl_p_ya[, c("psychopathic", "norm_kl")]
factor2_kl_p_oa <- df_factor_kl_p_oa[, c("psychopathic", "norm_kl")]

fit <- cor_test(factor2_kl_p_ya, factor2_kl_p_oa)
BF_result <- BF(fit, hypothesis = "norm_kl_with_psychopathic_in_g1 = norm_kl_with_psychopathic_in_g2")

# Factor 3, impulsive 
cocor_fct3_i <- cocor.indep.groups(corr_factor3_kl_i_ya$r[1, 2],
																	 corr_factor3_kl_i_oa$r[1, 2],
																	 n1 = 68, n2 = 73)

## Bayes factor 
set.seed(1213)
factor3_kl_i_ya <- df_factor_kl_i_ya[, c("empathy_emotional_motivation", "norm_kl")]
factor3_kl_i_oa <- df_factor_kl_i_oa[, c("empathy_emotional_motivation", "norm_kl")]

fit <- cor_test(factor3_kl_i_ya, factor3_kl_i_oa)
BF_result <- BF(fit, hypothesis = "norm_kl_with_empathy_emotional_motivation_in_g1 = norm_kl_with_empathy_emotional_motivation_in_g2")

# Factor 3, patient 
cocor_fct3_p <- cocor.indep.groups(corr_factor3_kl_p_ya$r[1, 2],
																	 corr_factor3_kl_p_oa$r[1, 2],
																	 n1 = 71, n2 = 68)

## Bayes factor 
set.seed(1213)
factor3_kl_p_ya <- df_factor_kl_p_ya[, c("empathy_emotional_motivation", "norm_kl")]
factor3_kl_p_oa <- df_factor_kl_p_oa[, c("empathy_emotional_motivation", "norm_kl")]

fit <- cor_test(factor3_kl_p_ya, factor3_kl_p_oa)
BF_result <- BF(fit, hypothesis = "norm_kl_with_empathy_emotional_motivation_in_g1 = norm_kl_with_empathy_emotional_motivation_in_g2")

# ============================================================================.
# S.1 Correlations between simulation and estimated parameters (Figure S2) ####
# ============================================================================.
## km estimated from the KU model vs simulation k, for impulsive others 
corr_km_k_i <- rcorr(data_all$other_i_km,
										 data_all$other_i_k,
										 type = "spearman")

corr_km_k_i_ci <- corr.test(data_all$other_i_km,
														data_all$other_i_k,
														method = "spearman")

plt <- ggplot(data_all, aes(x = other_i_km, y = other_i_k)) + 
	geom_point(colour = "#0057b7", size = 2.5, alpha = 0.5, shape = 19, stroke = 0, 
						 show.legend = FALSE) + 
	geom_smooth(colour = "#0057b7", fill = "#0057b7", method = "lm", alpha = 0.15, 
							formula = y ~ x, show.legend = FALSE) + 
	scale_x_continuous(name = paste0("<i>km</i>", " estimated from the KU model (impulsive)")) +  
	scale_y_continuous(name = paste0("<i>k</i>", " for simulation (impulsive)")) + 
	scale_linetype_manual(values = "solid") +
	theme_classic() + 
	theme(legend.position = "none") + 
	theme(text = element_text(family = "Arial"),
				axis.text = element_text(size = ax_text),
				axis.title.x = element_markdown(size = ax_title),
				axis.title.y = element_markdown(size = ax_title)) +
	guides(shape = guide_legend(override.aes = list(size = 2))) +  
	annotate(geom = "richtext", x = -2.5, y = -3,
					 label = paste0("_r_", "<sub>s</sub>(152)=",
					 							 round(corr_km_k_i$r[1,2], 3), 
					 							 ", _p_<0.001***"), 
					 fill = NA, label.color = NA, size = 4, 
					 family = "Arial") 

ggsave("data/plots/corr_km_k_i.png", plt, 
			 height = plot_h, width = plot_h * 1, dpi = resolution)

## km estimated from the KU model vs simulation k, for patient others 
corr_km_k_p <- rcorr(data_all$other_p_km,
										 data_all$other_p_k,
										 type = "spearman")

corr_km_k_p_ci <- corr.test(data_all$other_p_km,
														data_all$other_p_k,
														method = "spearman")

plt <- ggplot(data_all, aes(x = other_p_km, y = other_p_k)) + 
	geom_point(colour = "#fcb34b", size = 2.5, alpha = 0.5, shape = 19, stroke = 0, 
						 show.legend = FALSE) + 
	geom_smooth(colour = "#fcb34b", fill = "#fcb34b", method = "lm", alpha = 0.15, 
							formula = y ~ x, show.legend = FALSE) + 
	scale_x_continuous(name = paste0("<i>km</i>", " estimated from the KU model (patient)")) +  
	scale_y_continuous(name = paste0("<i>k</i>", " for simulation (patient)")) + 
	scale_linetype_manual(values = "solid") +
	theme_classic() + 
	theme(legend.position = "none") + 
	theme(text = element_text(family = "Arial"),
				axis.text = element_text(size = ax_text),
				axis.title.x = element_markdown(size = ax_title),
				axis.title.y = element_markdown(size = ax_title)) +
	guides(shape = guide_legend(override.aes = list(size = 2))) +  
	annotate(geom = "richtext", x = -3.5, y = -5,
					 label = paste0("_r_", "<sub>s</sub>(152)=",
					 							 round(corr_km_k_p$r[1,2], 3), 
					 							 ", _p_<0.001***"), 
					 fill = NA, label.color = NA, size = 4, 
					 family = "Arial") 

ggsave("data/plots/corr_km_k_p.png", plt, 
			 height = plot_h, width = plot_h * 1, dpi = resolution)

# ============================================================================.
# S.2 Mean & SD of km and ku in the KU model for all blocks (Table S3) ####
# ============================================================================.
# Self baseline 
data_self_1 <- data_all |> 
	select(id, group, self_1_km, self_1_ku)

km_self_1_summary <- get_summary(data_self_1, "self_1_km", "group")
ku_self_1_summary <- get_summary(data_self_1, "self_1_ku", "group")

# Impulsive Other
data_other_i <- data_all |>
	select(id, group, other_i_km, other_i_ku)

km_other_i_summary <- get_summary(data_other_i, "other_i_km", "group")
ku_other_i_summary <- get_summary(data_other_i, "other_i_ku", "group")

# Self after Impulsive Other
data_self_i <- data_all |>
	select(id, group, self_i_km, self_i_ku)

km_self_i_summary <- get_summary(data_self_i, "self_i_km", "group")
ku_self_i_summary <- get_summary(data_self_i, "self_i_ku", "group")

# Patient Other 
data_other_p <- data_all |>
	select(id, group, other_p_km, other_p_ku)

km_other_p_summary <- get_summary(data_other_p, "other_p_km", "group")
ku_other_p_summary <- get_summary(data_other_p, "other_p_ku", "group")

# Self after Patient Other 
data_self_p <- data_all |>
	select(id, group, self_p_km, self_p_ku)

km_self_p_summary <- get_summary(data_self_p, "self_p_km", "group")
ku_self_p_summary <- get_summary(data_self_p, "self_p_ku", "group")

# ============================================================================.
# S.3 Effect of the order of other's preferences ####
# ============================================================================.
# Preparation 
other_order <- data_all |>
	select(id, group, other1_pref) |> 
	mutate(id = as.factor(id))

data_kl_with_order <- data_kl |>
	left_join(other_order, by = c("id", "group"))

# LMM 1 - without self km as covariates
lmer_model_1 <- lmer(norm_kl ~ 1 + group * other_pref + other1_pref + (1|id),
										 data = data_kl_with_order)

lmer_model_stats_1 <- stats.rlmerMod.f(lmer_model_1)
row.names(lmer_model_stats_1) = c("Intercept",
																	"Group (older vs young)",
																	"Others (patient vs impulsive)",
																	"Order (patient first vs impulsive first)", 
																	"Group * Others")

kable(lmer_model_stats_1, align = rep("c", 6),
			caption = "Linear mixed effects model for signed KL divergence (with the order of others and without self km as covariates)") |>
	kable_styling()

## Bayes factor 
### BF_01 for the main effect of other's order 
set.seed(1213)
bf_null <- lmBF(norm_kl ~ group * other_pref + id, 
								data = data_kl_with_order, 
								whichRandom = "id")

bf_order <- lmBF(norm_kl ~ group * other_pref + other1_pref + id, 
								 data = data_kl_with_order,
								 whichRandom = "id")

bf01_order <- bf_null / bf_order

# LMM 2 - with self km as covariates
lmer_model_2 <- lmer(norm_kl ~ 1 + group * other_pref * centred_self_1_km 
										 + other1_pref + (1|id), data = data_kl_with_order)

lmer_model_stats_2 <- stats.rlmerMod.f(lmer_model_2)
row.names(lmer_model_stats_2) = c("Intercept",
																	"Group (older vs young)", 
																	"Others (patient vs impulsive)",
																	"Centred self km",
																	"Order (patient first vs impulsive first)", 
																	"Group * Others",
																	"Group * Centred self km",
																	"Others * Centred self km",
																	"Group * Others * Centred self km")

kable(lmer_model_stats_2, align = rep("c", 6),
			caption = "Linear mixed effects model for signed KL divergence (with the order of others and with self km as covariates)") |>
	kable_styling()

## Bayes factor 
### BF_01 for the main effect of other's order 
set.seed(1213)
bf_null <- lmBF(norm_kl ~ group * other_pref * centred_self_1_km + id, 
								data = data_kl_with_order, 
								whichRandom = "id")

bf_order <- lmBF(norm_kl ~ group * other_pref * centred_self_1_km + other1_pref + id, 
								 data = data_kl_with_order,
								 whichRandom = "id")

bf01_order_2 <- bf_null / bf_order

# ============================================================================.
# S.4 Partial correlations (Tables S1, S2) ####
# ============================================================================.
# Factor 3 & older impulsive 
data_pcorr_i_oa_factor <- data_all |>
	dplyr::select(id, group, norm_kl_i,
								WTAW_Standard, `ACE_(Memory/26)'`, `ACE_(Attention/18)'`) |>
	dplyr::filter(group == 2) |>
	dplyr::filter(!is.na(norm_kl_i)) |>
	mutate(id = factor(id)) |>
	left_join(data_qnr_factor, by = "id") |>
	dplyr::filter(!is.na(empathy_emotional_motivation))

## IQ 
fit1 <- lm(norm_kl_i ~ WTAW_Standard, data = data_pcorr_i_oa_factor)
fit2 <- lm(empathy_emotional_motivation ~ WTAW_Standard, data = data_pcorr_i_oa_factor)
pcorr_iq_i <- corr.test(resid(fit1), resid(fit2), method = "spearman")

## Attention 
fit1 <- lm(norm_kl_i ~ `ACE_(Attention/18)'`, data = data_pcorr_i_oa_factor)
fit2 <- lm(empathy_emotional_motivation ~ `ACE_(Attention/18)'`, data = data_pcorr_i_oa_factor)
pcorr_attention_i <- corr.test(resid(fit1), resid(fit2), method = "spearman")

## Memory 
fit1 <- lm(norm_kl_i ~ `ACE_(Memory/26)'`, data = data_pcorr_i_oa_factor)
fit2 <- lm(empathy_emotional_motivation ~ `ACE_(Memory/26)'`, data = data_pcorr_i_oa_factor)
pcorr_memory_i <- corr.test(resid(fit1), resid(fit2), method = "spearman")

# Factor 1 & older patient 
data_pcorr_p_oa_factor <- data_all |>
	dplyr::select(id, group, norm_kl_p,
								WTAW_Standard, `ACE_(Memory/26)'`, `ACE_(Attention/18)'`) |>
	dplyr::filter(group == 2) |>
	dplyr::filter(!is.na(norm_kl_p)) |>
	mutate(id = factor(id)) |>
	left_join(data_qnr_factor, by = "id") |>
	dplyr::filter(!is.na(autism_apathy_alexithymia))

## IQ 
fit1 <- lm(norm_kl_p ~ WTAW_Standard, data = data_pcorr_p_oa_factor)
fit2 <- lm(autism_apathy_alexithymia ~ WTAW_Standard, data = data_pcorr_p_oa_factor)
pcorr_iq_p <- corr.test(resid(fit1), resid(fit2), method = "spearman")

## Attention 
fit1 <- lm(norm_kl_p ~ `ACE_(Attention/18)'`, data = data_pcorr_p_oa_factor)
fit2 <- lm(autism_apathy_alexithymia ~ `ACE_(Attention/18)'`, data = data_pcorr_p_oa_factor)
pcorr_attention_p <- corr.test(resid(fit1), resid(fit2), method = "spearman")

## Memory 
fit1 <- lm(norm_kl_p ~ `ACE_(Memory/26)'`, data = data_pcorr_p_oa_factor)
fit2 <- lm(autism_apathy_alexithymia ~ `ACE_(Memory/26)'`, data = data_pcorr_p_oa_factor)
pcorr_memory_p <- corr.test(resid(fit1), resid(fit2), method = "spearman")

# ============================================================================.
# S.5 Model-free analysis ####
# ============================================================================.
# Preparation 
col_names <- c("id", "self_prop_ll_1", "self_prop_ll_i", "self_prop_ll_p")
self_prop_ll_mat <- matrix(nrow = n_pt, ncol = length(col_names))
colnames(self_prop_ll_mat) <- col_names
self_prop_ll_tb <- as_tibble(self_prop_ll_mat)
self_prop_ll_tb$id <- id_vec

rm(self_prop_ll_mat)

self_decision <- array(0, dim = c(n_pt, n_trials)) 

for (s_block in c("1", "i", "p")) {
	for (i_id in 1:n_pt) {
		self_decision[i_id, ] <- pt_new[, paste0("s_dec_", s_block), i_id] 
		
		self_ll <- self_decision[i_id, ] == 2 
		self_prop_ll <- sum(self_ll) / n_trials * 100 
		self_prop_ll_tb[i_id, paste0("self_prop_ll_",s_block)] <- self_prop_ll
	}
}

# Differences in proportions of LL choices 
new_cols <- as.data.frame(
	matrix(NA, nrow = nrow(self_prop_ll_tb), ncol = 2, 
				 dimnames = list(NULL, c("diff_prop_ll_i", "diff_prop_ll_p")))
)

self_prop_ll_tb <- cbind(self_prop_ll_tb, new_cols)

for (i_id in 1:n_pt) {
	self_prop_ll_tb[i_id, "diff_prop_ll_i"] <- 
		self_prop_ll_tb[i_id, "self_prop_ll_i"] - self_prop_ll_tb[i_id, "self_prop_ll_1"]
	
	self_prop_ll_tb[i_id, "diff_prop_ll_p"] <- 
		self_prop_ll_tb[i_id, "self_prop_ll_p"] - self_prop_ll_tb[i_id, "self_prop_ll_1"]
}

# Handle exceptions: when subjects had no Others of different preferences 
## Data were analysed for Others with the greatest distance between
## the subject’s own discount rate, and the Other's discount rate estimated by subject.
for (i_id in 1:n_pt) {
	diff_i <- data_all[i_id, "self_1_km"] - data_all[i_id, "other_i_km"]
	diff_p <- data_all[i_id, "self_1_km"] - data_all[i_id, "other_p_km"]
	
	# no Impulsive other
	if (check_sign(diff_i, diff_p) && diff_i > 0) { 
		if (abs(diff_i) > abs(diff_p)) {
			self_prop_ll_tb[i_id, "diff_prop_ll_p"] <- self_prop_ll_tb[i_id, "diff_prop_ll_i"]
			self_prop_ll_tb[i_id, "diff_prop_ll_i"] <- NA
		} else if (abs(diff_i) < abs(diff_p)) {
			self_prop_ll_tb[i_id, "diff_prop_ll_i"] <- NA
		}
	}
	
	# no Patient other 
	if (check_sign(diff_i, diff_p) && diff_i < 0) { 
		if (abs(diff_i) > abs(diff_p)) {
			self_prop_ll_tb[i_id, "diff_prop_ll_p"] <- NA
		} else if (abs(diff_i) < abs(diff_p)) {
			self_prop_ll_tb[i_id, "diff_prop_ll_i"] <- self_prop_ll_tb[i_id, "diff_prop_ll_p"] 
			self_prop_ll_tb[i_id, "diff_prop_ll_p"] <- NA
		}
	}
}

self_prop_ll_tb <- self_prop_ll_tb |>
	select(id, diff_prop_ll_i, diff_prop_ll_p)

df_self_prop_ll	<- tidyr::pivot_longer(self_prop_ll_tb,
																			 cols = starts_with("diff_prop_ll"),
																			 names_to = "other_pref",
																			 values_to = "diff_prop_ll") |> 
	mutate(other_pref = str_replace(other_pref, "diff_prop_ll_", "")) |>
	mutate(group = id %/% 100) |>
	mutate(other_pref = factor(other_pref), 
				 group = factor(group)) |> 
	mutate(condition = interaction(group, other_pref, sep = ":"))

df_self_prop_ll_unsigned <- df_self_prop_ll |>
	mutate(diff_prop_ll = case_when(other_pref == "i" ~ -diff_prop_ll, 
																	other_pref == "p" ~ diff_prop_ll))

## LMM (Table S5) ####  
lmer_model <- lmer(diff_prop_ll ~ 1 + other_pref * group + (1|id),
									 data = df_self_prop_ll_unsigned)

lmer_model_stats <- stats.rlmerMod.f(lmer_model)

## Post-hoc 
### Impulsive, young vs older 
wil_prop_impulsive_yo <- df_self_prop_ll_unsigned |> 
	filter(other_pref == "i") |> 
	wilcox_test(diff_prop_ll ~ group, paired = FALSE)

z_prop_impulsive_yo <- qnorm(wil_prop_impulsive_yo$p / 2)

wil_prop_impulsive_yo_es <- df_self_prop_ll_unsigned |>
	filter(other_pref == "i") |>
	wilcox_effsize(diff_prop_ll ~ group, paired = FALSE, ci = TRUE)

### Patient, young vs older 
wil_prop_patient_yo <- df_self_prop_ll_unsigned |> 
	filter(other_pref == "p") |> 
	wilcox_test(diff_prop_ll ~ group, paired = FALSE)

z_prop_patient_yo <- qnorm(wil_prop_patient_yo$p / 2)

wil_prop_patient_yo_es <- df_self_prop_ll_unsigned |>
	filter(other_pref == "p") |>
	wilcox_effsize(diff_prop_ll ~ group, paired = FALSE, ci = TRUE)

#### Save results for JASP 
data_prop_p <- df_self_prop_ll_unsigned |>
	filter(other_pref == "p")

write_csv(data_prop_p, 
					file = "data/jasp_prop_p.csv")

### Young, impulsive vs patient
wil_prop_young_ip <- df_self_prop_ll_unsigned |> 
	filter(group == 1) |> 
	group_by(id) |> 
	filter(!any(is.na(diff_prop_ll))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_test(diff_prop_ll ~ other_pref, paired = TRUE)

z_prop_young_ip <- qnorm(wil_prop_young_ip$p / 2)

wil_prop_young_ip_es <- df_self_prop_ll_unsigned |>
	filter(group == 1) |> 
	group_by(id) |> 
	filter(!any(is.na(diff_prop_ll))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_effsize(diff_prop_ll ~ other_pref, paired = TRUE, ci = TRUE)

### Older, impulsive vs patient
wil_prop_older_ip <- df_self_prop_ll_unsigned |> 
	filter(group == 2) |> 
	group_by(id) |> 
	filter(!any(is.na(diff_prop_ll))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_test(diff_prop_ll ~ other_pref, paired = TRUE)

z_prop_older_ip <- qnorm(wil_prop_older_ip$p / 2)

wil_prop_older_ip_es <- df_self_prop_ll_unsigned |>
	filter(group == 2) |>
	group_by(id) |> 
	filter(!any(is.na(diff_prop_ll))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_effsize(diff_prop_ll ~ other_pref, paired = TRUE, ci = TRUE)

#### Save results for JASP 
data_prop_older <- self_prop_ll_tb |>
	as_tibble() |>
	filter(id > 200) |>
	mutate(diff_prop_ll_i = -diff_prop_ll_i)

write_csv(data_prop_older, 
					file = "data/jasp_prop_older.csv")

## Plot (Figure S4) #### 
vars <- "diff_prop_ll"

interaction_lims <- c("1:i", "1:p", "2:i", "2:p")
interaction_labs <- c("", " Young", "", " Older")
interaction_bar_cols <- c("#0057b7", "#ffd700", "#dbe0ff", "#f1e1bc")
interaction_data_cols <- c("#314498", "#fcb34b", "#8d8cbe", "#fecd91")

# Unsigned 
plt_data <- df_self_prop_ll_unsigned |> 
	group_by(group, other_pref, condition) |> 
	filter(!is.na(get(vars))) |>
	dplyr::summarise(n = n(), 
									 mean = mean(get(vars), na.rm = TRUE), 
									 sd = sd(get(vars), na.rm = TRUE), 
									 se = sd / sqrt(n)) |> 
	ungroup()

plt <- ggplot(plt_data, aes(x = other_pref, y = mean, group = condition)) + 
	geom_bar(aes(fill = condition), 
					 stat = "identity", 
					 position = position_dodge(0.75), 
					 width = 0.75, linewidth = 0.3, 
					 colour = "black") + 
	geom_jitter(data = df_self_prop_ll_unsigned,
							aes(x = other_pref, y = get(vars), colour = condition), 
							position = position_jitterdodge(dodge.width = 0.75, 
																							jitter.height = 0, 
																							jitter.width = 0.3), 
							size = 1.3, alpha = 0.3, shape = 19, stroke = 0,   
							show.legend = FALSE) + 
	geom_errorbar(aes(ymin = mean - se, 
										ymax = mean + se), 
								position = position_dodge(0.75), 
								colour = "black", width = 0, alpha = 1, size = 0.35) + 
	scale_x_discrete(labels = c("Impulsive", "Patient")) +
	scale_y_continuous(limits = c(-50, 55), expand = c(0, 0)) + 
	labs(x = NULL, y = "Differences in choice proportions of LL (%)") +  
	scale_fill_manual(name = "condition", 
										limits = interaction_lims, 
										values = interaction_bar_cols, 
										labels = interaction_labs) + 
	scale_colour_manual(name = "condition", 
											limits = interaction_lims, 
											values = interaction_data_cols, 
											labels = interaction_labs) +
	theme_classic() +
	geom_hline(yintercept = 0, linetype = "dashed") + 
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(0.5, 0.5), colour = NA),
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_title), 
				axis.title.x = element_text(size = ax_title),
				axis.title.y = element_text(size = ax_title), 
				legend.text = element_text(size = 8, margin = margin(l = 3)),
				legend.key.spacing.x = unit(-0.13, "cm"),
				legend.key.spacing.y = unit(0.2, "cm"),
				legend.key.height = unit(0.35, "cm"),
				legend.key.width = unit(0.35, "cm"), 
				legend.position = c(0.5, 0.1)) +
	guides(fill = guide_legend(ncol = 2, byrow = TRUE), 
				 colour = guide_legend(ncol = 2, byrow = TRUE)) + 
	geom_signif(xmin = c(0.85), xmax = c(1.85),  
							y_position = c(40), tip_length = 0, vjust = 0.5, size = 0.25,
							textsize = 3.5, family = "Arial",  
							annotation = c("*")) +
	geom_signif(xmin = c(1.15), xmax = c(2.15),  
							y_position = c(30), tip_length = 0, vjust = 0, size = 0.25,
							textsize = 3, family = "Arial",  
							annotation = c("ns")) + 
	geom_signif(xmin = c(0.85), xmax = c(1.15),  
							y_position = c(20), tip_length = 0, vjust = 0.5, size = 0.25,
							textsize = 3.5, family = "Arial",  
							annotation = c("**")) +
	geom_signif(xmin = c(1.85), xmax = c(2.15),  
							y_position = c(20), tip_length = 0, vjust = 0, size = 0.25,
							textsize = 3, family = "Arial",  
							annotation = c("ns")) 

ggsave("data/plots/prop_ll_unsigned_self1.png", plt, 
			 height = plot_h, width = plot_h * 0.8, dpi = resolution)

# ============================================================================.
# S.6 Outliers ####
# ============================================================================.
# Set outlier criteria: above or below three standard deviations
## For young impulsive 
lower_ya_i <- kl_summary[["1:i"]]$mean - 3 * kl_summary[["1:i"]]$sd
upper_ya_i <- kl_summary[["1:i"]]$mean + 3 * kl_summary[["1:i"]]$sd

out_id_ya_i <- data_all |>
	filter(group == 1) |>
	filter(norm_kl_i < lower_ya_i | norm_kl_i > upper_ya_i) |> 
	pull(id)

## For older impulsive 
lower_oa_i <- kl_summary[["2:i"]]$mean - 3 * kl_summary[["2:i"]]$sd
upper_oa_i <- kl_summary[["2:i"]]$mean + 3 * kl_summary[["2:i"]]$sd

out_id_oa_i <- data_all |>
	filter(group == 2) |>
	filter(norm_kl_i < lower_oa_i | norm_kl_i > upper_oa_i) |>
	pull(id)

## For young patient 
lower_ya_p <- kl_summary[["1:p"]]$mean - 3 * kl_summary[["1:p"]]$sd
upper_ya_p <- kl_summary[["1:p"]]$mean + 3 * kl_summary[["1:p"]]$sd

out_id_ya_p <- data_all |>
	filter(group == 1) |>
	filter(norm_kl_p < lower_ya_p | norm_kl_p > upper_ya_p) |>
	pull(id)

## For older patient 
lower_oa_p <- kl_summary[["2:p"]]$mean - 3 * kl_summary[["2:p"]]$sd
upper_oa_p <- kl_summary[["2:p"]]$mean + 3 * kl_summary[["2:p"]]$sd

out_id_oa_p <- data_all |>
	filter(group == 2) |>
	filter(norm_kl_p < lower_oa_p | norm_kl_p > upper_oa_p) |> 
	pull(id)

data_kl_wo_outliers <- data_kl |> 
	mutate(norm_kl = replace(norm_kl, 
													 id %in% out_id_ya_i & other_pref == "i", 
													 NA)) |>
	mutate(norm_kl = replace(norm_kl, 
													 id %in% out_id_oa_i & other_pref == "i", 
													 NA)) |>
	mutate(norm_kl = replace(norm_kl, 
													 id %in% out_id_ya_p & other_pref == "p", 
													 NA)) |>
	mutate(norm_kl = replace(norm_kl, 
													 id %in% out_id_oa_p & other_pref == "p", 
													 NA)) |>
	na.omit(norm_kl)

## Plot (Figure S5) #### 
vars <- "norm_kl"

interaction_lims <- c("1:i", "1:p", "2:i", "2:p")
interaction_labs <- c("", " Young", "", " Older")
interaction_bar_cols <- c("#0057b7", "#ffd700", "#dbe0ff", "#f1e1bc")
interaction_data_cols <- c("#314498", "#fcb34b", "#8d8cbe", "#fecd91")

plt_data <- data_kl_wo_outliers |> 
	group_by(group, other_pref, condition1) |> 
	filter(!is.na(get(vars))) |>
	dplyr::summarise(n = n(), 
									 mean = mean(get(vars), na.rm = TRUE), 
									 sd = sd(get(vars), na.rm = TRUE), 
									 se = sd / sqrt(n)) |> 
	ungroup()

limits <- data_kl_wo_outliers |> 
	summarise(max = max(get(vars), na.rm = TRUE), 
						min = min(get(vars), na.rm = TRUE))

set.seed(1213)
# Main function
plt <- ggplot(plt_data, aes(x = other_pref, y = mean, group = condition1)) + 
	geom_bar(aes(fill = condition1), 
					 stat = "identity", 
					 position = position_dodge(0.75), 
					 width = 0.75, linewidth = 0.3, 
					 colour = "black") + 
	geom_jitter(data = data_kl_wo_outliers,
							aes(x = other_pref, y = get(vars), colour = condition1), 
							position = position_jitterdodge(dodge.width = 0.75, 
																							jitter.height = 0, 
																							jitter.width = 0.3), 
							size = 1.3, alpha = 0.3, shape = 19, stroke = 0,   
							show.legend = FALSE) + 
	geom_errorbar(aes(ymin = mean - se, 
										ymax = mean + se), 
								position = position_dodge(0.75), 
								colour = "black", width = 0, alpha = 1, size = 0.35) + 
	scale_x_discrete(labels = c("Impulsive", "Patient")) +
	scale_y_continuous(limits = c(floor(limits$min) * 1.1, ceiling(limits$max) * 1.1), 
										 expand = c(0, 0)) + 
	labs(x = NULL, y = "Signed KL divergence") +
	scale_fill_manual(name = "condition1", 
										limits = interaction_lims, 
										values = interaction_bar_cols, 
										labels = interaction_labs) +
	scale_colour_manual(name = "condition1", 
											limits = interaction_lims, 
											values = interaction_data_cols, 
											labels = interaction_labs) +
	theme_classic() + 
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(0.5, 0.5), colour = NA),
				axis.text.x = element_text(size = ax_title, colour = "black"),
				axis.text.y = element_text(size = ax_title), 
				axis.title.x = element_text(size = ax_title),
				axis.title.y = element_text(size = ax_title), 
				legend.text = element_text(size = 8, margin = margin(l = 3)),
				legend.key.spacing.x = unit(-0.13, "cm"),
				legend.key.spacing.y = unit(0.2, "cm"),
				legend.key.height = unit(0.35, "cm"),
				legend.key.width = unit(0.35, "cm"), 
				legend.position = c(0.5, 0.2)) +
	guides(fill = guide_legend(ncol = 2, byrow = TRUE), 
				 colour = guide_legend(ncol = 2, byrow = TRUE)) +
	geom_hline(yintercept = 0, linetype = "dashed") + 
	geom_signif(xmin = c(0.85), xmax = c(1.15),  
							y_position = c(1), tip_length = 0, vjust = 0.5, size = 0.25,
							textsize = 3.5, family = "Arial",  
							annotation = c("*")) +
	geom_signif(xmin = c(1.85), xmax = c(2.15),  
							y_position = c(1), tip_length = 0, vjust = 0, size = 0.25,
							textsize = 3, family = "Arial",  
							annotation = c("ns")) +
	geom_signif(xmin = c(0.85), xmax = c(1.85),  
							y_position = c(2.5), tip_length = 0, vjust = 0.5, size = 0.25,
							textsize = 3.5, family = "Arial",  
							annotation = c("***")) +
	geom_signif(xmin = c(1.15), xmax = c(2.15),  
							y_position = c(1.75), tip_length = 0, vjust = 0, size = 0.25,
							textsize = 3, family = "Arial",  
							annotation = c("ns")) +
	annotation_custom(grob = grid::textGrob("\u2190 shift away", x = -0.23, y = 0.20, rot = 90,
																					gp = grid::gpar(fontsize = 9, fontface = "italic", fontfamily = "Arial"))) +
	annotation_custom(grob = grid::textGrob("shift toward \u2192", x = -0.23, y = 0.80, rot = 90,
																					gp = grid::gpar(fontsize = 9, fontface = "italic", fontfamily = "Arial"))) +
	coord_cartesian(expand = TRUE, clip = "off") + 
	theme(plot.margin = unit(c(0.2, 0.2, 0, 1), "cm"))

# Save the plot 
ggsave("data/plots/norm_kl_wo_outliers.png", plt, 
			 height = plot_h, width = plot_h * 0.8, dpi = resolution)

# LMM
## LMM 1 - without self km as covariates (Table S6) ####
lmer_model_1 <- lmer(norm_kl ~ 1 + group * other_pref + (1|id),
										 data = data_kl_wo_outliers)

lmer_model_stats_1 <- stats.rlmerMod.f(lmer_model_1)
row.names(lmer_model_stats_1) = c("Intercept",
																	"Group (older vs young)",
																	"Others (patient vs impulsive)",
																	"Group * Others")

kable(lmer_model_stats_1, align = rep("c", 6),
			caption = "Linear mixed effects model for signed KL divergence 
			(without outliers and without self km as covariates)") |>
	kable_styling()

## LMM 2 - with self km as covariates (Table S7) #### 
lmer_model_2 <- lmer(norm_kl ~ 1 + group * other_pref * centred_self_1_km + (1|id),
										 data = data_kl_wo_outliers)

lmer_model_stats_2 <- stats.rlmerMod.f(lmer_model_2)
row.names(lmer_model_stats_2) = c("Intercept",
																	"Group (older vs young)",
																	"Others (patient vs impulsive)",
																	"Self baseline km",
																	"Group * Others", 
																	"Group * self baseline km", 
																	"Others * self baseline km", 
																	"Group * Others * self baseline km")

kable(lmer_model_stats_2, align = rep("c", 6),
			caption = "Linear mixed effects model for signed KL divergence 
			(without outliers and with self km as covariates)") |>
	kable_styling()

## Post-hoc 
### Impulsive, young vs older 
wil_kl_impulsive_yo_wo_outliers <- data_kl_wo_outliers |> 
	filter(other_pref == "i") |> 
	wilcox_test(norm_kl ~ group, paired = FALSE)

z_kl_impulsive_yo_wo_outliers <- qnorm(wil_kl_impulsive_yo_wo_outliers$p / 2)

wil_kl_impulsive_yo_es_wo_outliers <- data_kl_wo_outliers |>
	filter(other_pref == "i") |>
	wilcox_effsize(norm_kl ~ group, paired = FALSE, ci = TRUE)

### Patient, young vs older 
wil_kl_patient_yo_wo_outliers <- data_kl_wo_outliers |> 
	filter(other_pref == "p") |> 
	wilcox_test(norm_kl ~ group, paired = FALSE)

z_kl_patient_yo_wo_outliers <- qnorm(wil_kl_patient_yo_wo_outliers$p / 2)

wil_kl_patient_yo_es_wo_outliers <- data_kl_wo_outliers |>
	filter(other_pref == "p") |>
	wilcox_effsize(norm_kl ~ group, paired = FALSE, ci = TRUE)

### Young, impulsive vs patient
wil_kl_young_ip_wo_outliers <- data_kl_wo_outliers |> 
	filter(group == 1) |> 
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_test(norm_kl ~ other_pref, paired = TRUE)

z_kl_young_ip_wo_outliers <- qnorm(wil_kl_young_ip_wo_outliers$p / 2)

wil_kl_young_ip_es_wo_outliers <- data_kl_wo_outliers |>
	filter(group == 1) |> 
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_effsize(norm_kl ~ other_pref, paired = TRUE, ci = TRUE)

### Older, impulsive vs patient
wil_kl_older_ip_wo_outliers <- data_kl_wo_outliers |> 
	filter(group == 2) |> 
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_test(norm_kl ~ other_pref, paired = TRUE)

z_kl_older_ip_wo_outliers <- qnorm(wil_kl_older_ip_wo_outliers$p / 2)

wil_kl_older_ip_es_wo_outliers <- data_kl_wo_outliers |>
	filter(group == 2) |>
	group_by(id) |> 
	filter(!any(is.na(norm_kl))) |>
	filter(n() == 2) |>
	ungroup() |> 
	wilcox_effsize(norm_kl ~ other_pref, paired = TRUE, ci = TRUE)

# Save results for JASP 
data_kl_p_wo_outliers <- data_kl_wo_outliers |>
	filter(other_pref == "p")

write_csv(data_kl_p_wo_outliers, 
					file = "data/jasp_kl_p_wo_outliers.csv")

data_kl_older_wo_outliers <- data_all |>
	select(id, group, norm_kl_i, norm_kl_p) |> 
	filter(group == 2) |> 
	filter(!(id %in% out_id_oa_i | id %in% out_id_oa_p))

write_csv(data_kl_older_wo_outliers, 
					file = "data/jasp_kl_older_wo_outliers.csv")

## Correlations (Table S10) #### 
data_qnr_factor_kl_wo_outliers <- data_qnr_factor |>
	left_join(data_kl_wo_outliers, by = c("id", "group")) |>
	select(autism_apathy_alexithymia, psychopathic, empathy_emotional_motivation, 
				 group, norm_kl, other_pref, condition1)

### Young impulsive
#### Factor 1
data_qnr_factor_kl_wo_outliers_i_ya <- data_qnr_factor_kl_wo_outliers |>
	filter(group == 1, other_pref == "i")

corr_fct1_kl_wo_outliers_i_ya <- rcorr(data_qnr_factor_kl_wo_outliers_i_ya$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_i_ya$autism_apathy_alexithymia, 
																			 type = "spearman")

corr_fct1_kl_wo_outliers_i_ya_ci <- corr.test(data_qnr_factor_kl_wo_outliers_i_ya$norm_kl,
																							data_qnr_factor_kl_wo_outliers_i_ya$autism_apathy_alexithymia,
																							method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_i_ya$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_i_ya$autism_apathy_alexithymia)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

#### Factor 2 
corr_fct2_kl_wo_outliers_i_ya <- rcorr(data_qnr_factor_kl_wo_outliers_i_ya$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_i_ya$psychopathic, 
																			 type = "spearman")

corr_fct2_kl_wo_outliers_i_ya_ci <- corr.test(data_qnr_factor_kl_wo_outliers_i_ya$norm_kl,
																							data_qnr_factor_kl_wo_outliers_i_ya$psychopathic,
																							method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_i_ya$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_i_ya$psychopathic)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

#### Factor 3
corr_fct3_kl_wo_outliers_i_ya <- rcorr(data_qnr_factor_kl_wo_outliers_i_ya$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_i_ya$empathy_emotional_motivation, 
																			 type = "spearman")

corr_fct3_kl_wo_outliers_i_ya_ci <- corr.test(data_qnr_factor_kl_wo_outliers_i_ya$norm_kl,
																							data_qnr_factor_kl_wo_outliers_i_ya$empathy_emotional_motivation,
																							method = "spearman")

set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_i_ya$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_i_ya$empathy_emotional_motivation)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)
bf01 <- 1 / bf10

# Older impulsive 
data_qnr_factor_kl_wo_outliers_i_oa <- data_qnr_factor_kl_wo_outliers |>
	filter(group == 2, other_pref == "i")

## fct 1: r = -0.08, p = 0.4898
corr_fct1_kl_wo_outliers_i_oa <- rcorr(data_qnr_factor_kl_wo_outliers_i_oa$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_i_oa$autism_apathy_alexithymia, 
																			 type = "spearman")

corr_fct1_kl_wo_outliers_i_oa_ci <- corr.test(data_qnr_factor_kl_wo_outliers_i_oa$norm_kl,
																							data_qnr_factor_kl_wo_outliers_i_oa$autism_apathy_alexithymia,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_i_oa$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_i_oa$autism_apathy_alexithymia)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

## fct 2: r = 0.05, p = 0.6713
corr_fct2_kl_wo_outliers_i_oa <- rcorr(data_qnr_factor_kl_wo_outliers_i_oa$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_i_oa$psychopathic, 
																			 type = "spearman")

corr_fct2_kl_wo_outliers_i_oa_ci <- corr.test(data_qnr_factor_kl_wo_outliers_i_oa$norm_kl,
																							data_qnr_factor_kl_wo_outliers_i_oa$psychopathic,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_i_oa$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_i_oa$psychopathic)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

## fct 3: r = 0.28, p = 0.018
corr_fct3_kl_wo_outliers_i_oa <- rcorr(data_qnr_factor_kl_wo_outliers_i_oa$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_i_oa$empathy_emotional_motivation, 
																			 type = "spearman")

corr_fct3_kl_wo_outliers_i_oa_ci <- corr.test(data_qnr_factor_kl_wo_outliers_i_oa$norm_kl,
																							data_qnr_factor_kl_wo_outliers_i_oa$empathy_emotional_motivation,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_i_oa$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_i_oa$empathy_emotional_motivation)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

# Young patient
data_qnr_factor_kl_wo_outliers_p_ya <- data_qnr_factor_kl_wo_outliers |>
	filter(group == 1, other_pref == "p")

## fct 1: r = -0.03, p = 0.8368
corr_fct1_kl_wo_outliers_p_ya <- rcorr(data_qnr_factor_kl_wo_outliers_p_ya$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_p_ya$autism_apathy_alexithymia, 
																			 type = "spearman")

corr_fct1_kl_wo_outliers_p_ya_ci <- corr.test(data_qnr_factor_kl_wo_outliers_p_ya$norm_kl,
																							data_qnr_factor_kl_wo_outliers_p_ya$autism_apathy_alexithymia,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_p_ya$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_p_ya$autism_apathy_alexithymia)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

## fct 2: r = -0.07, p = 0.5673
corr_fct2_kl_wo_outliers_p_ya <- rcorr(data_qnr_factor_kl_wo_outliers_p_ya$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_p_ya$psychopathic, 
																			 type = "spearman")

corr_fct2_kl_wo_outliers_p_ya_ci <- corr.test(data_qnr_factor_kl_wo_outliers_p_ya$norm_kl,
																							data_qnr_factor_kl_wo_outliers_p_ya$psychopathic,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_p_ya$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_p_ya$psychopathic)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

## fct 3: r = -0.09, p = 0.4834
corr_fct3_kl_wo_outliers_p_ya <- rcorr(data_qnr_factor_kl_wo_outliers_p_ya$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_p_ya$empathy_emotional_motivation, 
																			 type = "spearman")

corr_fct3_kl_wo_outliers_p_ya_ci <- corr.test(data_qnr_factor_kl_wo_outliers_p_ya$norm_kl,
																							data_qnr_factor_kl_wo_outliers_p_ya$empathy_emotional_motivation,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_p_ya$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_p_ya$empathy_emotional_motivation)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

# Older patient 
data_qnr_factor_kl_wo_outliers_p_oa <- data_qnr_factor_kl_wo_outliers |>
	filter(group == 2, other_pref == "p")

## fct 1: r = 0.35, p = 0.0033
corr_fct1_kl_wo_outliers_p_oa <- rcorr(data_qnr_factor_kl_wo_outliers_p_oa$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_p_oa$autism_apathy_alexithymia, 
																			 type = "spearman")

corr_fct1_kl_wo_outliers_p_oa_ci <- corr.test(data_qnr_factor_kl_wo_outliers_p_oa$norm_kl,
																							data_qnr_factor_kl_wo_outliers_p_oa$autism_apathy_alexithymia,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_p_oa$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_p_oa$autism_apathy_alexithymia)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

## fct 2: r = 0.16, p = 0.2016
corr_fct2_kl_wo_outliers_p_oa <- rcorr(data_qnr_factor_kl_wo_outliers_p_oa$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_p_oa$psychopathic, 
																			 type = "spearman")

corr_fct2_kl_wo_outliers_p_oa_ci <- corr.test(data_qnr_factor_kl_wo_outliers_p_oa$norm_kl,
																							data_qnr_factor_kl_wo_outliers_p_oa$psychopathic,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_p_oa$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_p_oa$psychopathic)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10

## fct 3: r = -0.11, p = 0.3788
corr_fct3_kl_wo_outliers_p_oa <- rcorr(data_qnr_factor_kl_wo_outliers_p_oa$norm_kl, 
																			 data_qnr_factor_kl_wo_outliers_p_oa$empathy_emotional_motivation, 
																			 type = "spearman")

corr_fct3_kl_wo_outliers_p_oa_ci <- corr.test(data_qnr_factor_kl_wo_outliers_p_oa$norm_kl,
																							data_qnr_factor_kl_wo_outliers_p_oa$empathy_emotional_motivation,
																							method = "spearman")

### Bayes factor 
set.seed(1213)
rho_samples <- spearmanGibbsSampler(xVals = data_qnr_factor_kl_wo_outliers_p_oa$norm_kl,
																		yVals = data_qnr_factor_kl_wo_outliers_p_oa$empathy_emotional_motivation)

bf10 <- computeBayesFactorOneZero(rho_samples$rhoSamples,
																	whichTest = "Spearman",
																	priorParameter = 1)

bf01 <- 1 / bf10



