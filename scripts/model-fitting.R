# ============================================================================.
# Info ####
# ============================================================================.
# Estimate and visualise the fitting of Stan models using WAIC and LOO-IC.
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Project: Zhilin's first PhD project - Social discounting across the adult lifespan 
# Last revision: 20 Nov 2023, Zhilin Su 

# ============================================================================.
# Preparation ####
# ============================================================================.
rm(list = ls())

library(rstan)
library(loo)
library(ggplot2)
library(ggpattern)
library(ggtext)
library(dplyr)
library(grid)

# Index vectors 
model_file <- c("kt", "ku0", "ku1", "ku2")
group <- c("ya", "oa")
self_block <- c("1", "i", "p")
other_block <- c("i", "p")

# Import files
folder_path <- "data/stanfit/"

stanfit_list <- list()

for (i_model in 1:length(model_file)) {
	file <- model_file[i_model] 
	
	for (i_group in group) {
		for (i_self_block in self_block) {
			file_name <- paste("stanfit_", file, "_", i_group, "_self_", i_self_block, ".rds", sep = "")
			
			full_file_name <- paste(folder_path, "/", file_name, sep = "")
			object_name <- paste(file, "_", i_group, "_self_", i_self_block, sep = "")
			stanfit_list[[object_name]] <- readRDS(full_file_name)
		}
		
		for (i_other_block in other_block) {
			file_name <- paste("stanfit_", file, "_", i_group, "_other_", i_other_block, ".rds", sep = "")
			
			full_file_name <- paste(folder_path, "/", file_name, sep = "")
			object_name <- paste(file, "_", i_group, "_other_", i_other_block, sep = "")
			stanfit_list[[object_name]] <- readRDS(full_file_name)
		}
	}
}

# ============================================================================.
# Model fitting ####
# ============================================================================.
# Initialise lists 
fit_tb <- tibble(model = character(), waic = numeric(), loo = numeric())

for (i_model in 1:length(stanfit_list)) {
	stanfit <- stanfit_list[[i_model]]
	model_name <- names(stanfit_list)[[i_model]]
	
	log_lik <- extract_log_lik(stanfit, "log_lik")
	
	waic <- waic(log_lik)
	waic_estimate <- waic$estimates["waic", "Estimate"]
	
	loo <- loo(log_lik)
	looic_estimate <- loo$estimates["looic", "Estimate"]
	
	new_row <- tibble(model = model_name, 
										waic = waic_estimate, 
										loo = looic_estimate)
	
	fit_tb <- rbind(fit_tb, new_row)
}

saveRDS(fit_tb, file = "data/fit.rds")

# ============================================================================.
# Plotting ####
# ============================================================================.
# Load files ---- 
fit_tb <- readRDS("data/fit_tb.rds")

# General settings of plotting 
resolution <- 1200
plot_h <- 4
ax_text <- 12
ax_title <- 14

model_vec <- fit_tb$model

witht_grouping_index <- grep("_ya_|_oa_", model_vec)
fit_with_grouping <- fit_tb[witht_grouping_index, ]

with_kt_index <- grep("kt_", fit_with_grouping$model)
with_ku0_index <- grep("ku0_", fit_with_grouping$model)
with_ku1_index <- grep("ku1_", fit_with_grouping$model)
with_ku2_index <- grep("ku2_", fit_with_grouping$model)

with_ya_index <- grep("_ya_", fit_with_grouping$model)
with_oa_index <- grep("_oa_", fit_with_grouping$model)

with_self_index  <- grep("_self_", fit_with_grouping$model)
with_other_index <- grep("_other_", fit_with_grouping$model)

with_1_index <- grep("_1", fit_with_grouping$model)
with_i_index <- grep("_i", fit_with_grouping$model)
with_p_index <- grep("_p", fit_with_grouping$model)

# Create helper columns
fit_with_grouping <- fit_with_grouping %>% 
	mutate(model_type = case_when(
		row_number() %in% with_kt_index ~ "kt",
		row_number() %in% with_ku0_index ~ "ku0",
		row_number() %in% with_ku1_index ~ "ku1",
		row_number() %in% with_ku2_index ~ "ku2",
		TRUE ~ NA_character_
	)) %>%
	mutate(group = case_when(
		row_number() %in% with_ya_index ~ "young", 
		row_number() %in% with_oa_index ~ "older", 
		TRUE ~ NA_character_
	)) %>% 
	mutate(agent = case_when(
		row_number() %in% with_self_index  ~ "self",
		row_number() %in% with_other_index ~ "other",
		TRUE ~ NA_character_
	)) %>%
	mutate(block = case_when(
		row_number() %in% with_1_index ~ "1", 
		row_number() %in% with_i_index ~ "i", 
		row_number() %in% with_p_index ~ "p", 
		TRUE ~ NA_character_
 	)) %>%
	mutate(model_type = factor(model_type, levels = c("kt", "ku0", "ku1", "ku2")),
				 group = factor(group, levels = c("young", "older")), 
				 agent = factor(agent, levels = c("self", "other")), 
				 block = factor(block, levels = c("1", "i", "p"))) %>%
	mutate(condition1 = interaction(model_type, group, sep = ":")) %>% 
	mutate(condition2 = interaction(model_type, group, agent, sep = ":")) %>%
	mutate(condition3 = interaction(model_type, group, agent, block, sep = ":"))

# Plot (fig. 1C)
## Preparation 
plt_data <- fit_with_grouping %>%
	group_by(model_type) %>%
	summarise(waic = sum(waic), loo = sum(loo)) %>% 
	ungroup()

limits <- plt_data %>% 
	summarise(min_loo = min(loo), 
						max_loo = max(loo))

plt_data <- plt_data %>%
	mutate(delta_loo = loo - limits$min_loo)

## Main function - delta_loo 
plt_delta_loo <- ggplot(plt_data, aes(x = model_type, y = delta_loo)) + 
	geom_bar(stat = "identity",
					 width = 0.5, linewidth = 0.3, 
					 colour = "black", fill = "#00a087") +
	geom_text(aes(label = round(delta_loo)), 
						vjust = -0.5, size = 6, colour = "black", alpha = 1, 
						family = "Arial") + 
	scale_x_discrete(labels = c("KT", 
															"KU<br>
															<span style='font-size:10pt'>w/o noise</span>", 
															"KU<br>
															<span style='font-size:10pt'>w/ self noise</span>", 
															"KU <br>
															<span style='font-size:10pt'>w/ other noise</span>")) +
	scale_y_continuous(limits = c(0, 5050), expand = c(0, 0)) + 
	labs(x = NULL, y = "\u0394LOO-IC") +
	theme_classic() + 
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(2, 2), color = NA),
				axis.text.x = element_markdown(size = ax_text, 
																	 colour = "black"),
				axis.text.y = element_markdown(size = ax_text), 
				axis.title.x = element_markdown(size = ax_title), 
				axis.title.y = element_markdown(size = ax_title),
				legend.text = element_text(size = ax_text, margin = margin(l = 3)),
				legend.spacing.x = unit(-0.07, "cm"),
				legend.spacing.y = unit(-0.05, "cm"),
				legend.key.height = unit(1, "cm"),
				legend.key.width = unit(0.5, "cm")) +
	annotation_custom(grob = textGrob("← better model fit", x = -0.23, y = 0.2, rot = 90,
										gp = gpar(fontsize = 11, fontface = "italic", fontfamily = "Arial"))) +
	coord_cartesian(expand = TRUE, clip = "off") + 
	theme(plot.margin = unit(c(0.2, 0.2, 0, 1), "cm"))

ggsave("data/plots/fit.png", plt_delta_loo, 
			 height = plot_h, width = plot_h * 1.25, dpi = resolution)
