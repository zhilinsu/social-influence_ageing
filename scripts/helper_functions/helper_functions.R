# ============================================================================.
# Info ####
# ============================================================================.
# The helper functions for the analysis of social discounting task.   
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Project: Zhilin's first PhD project - Social discounting across the adult lifespan 
# Last revision: 07 Aug 2023, Zhilin Su 

# ============================================================================.
# Linear mixed effects models ####
# ============================================================================.
# Function to produce reporting stats for rlmer objects
#
# Wald confidence interval code adapted from Ben Bolker at 
# https://gist.github.com/kamermanpr/aaa598485b6e990017375359ff5f4533)
#
# Then z-value and p-value generation added by Anthony Gabay 12/08/19
#
# The script was fetched from Jo Culer's OSF at https://osf.io/xgw7h/ 

stats.rlmerMod <- function(object, level = 0.95) {
	
	# Extract beta coefficients
	beta <- fixef(object)
	
	# Extract names of coefficients
	parm <- names(beta)
	
	# Extract standard errors for the coefficients
	se <- sqrt(diag(vcov(object)))
	
	# Set level of confidence interval
	conf.level <- qnorm((1 + level) / 2)
	
	# Calculate z value
	z = beta/se
	
	# Calculate CI and create table
	ctab <- cbind(beta,
								beta - (conf.level * se), 
								beta + (conf.level * se),
								se,
								z,
								2*pnorm(-abs(z)))
	
	
	# label column names
	colnames(ctab) <- c('beta',
											paste(100 * ((1 - level) / 2), '%'),
											paste(100 * ((1 + level) / 2), '%'),
											'SE',
											'z',
											'p')
	
	# Output
	return(ctab)
	
}

## Function to format p values taken from 
## https://stackoverflow.com/questions/23018256/printing-p-values-with-0-001

pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {
	
	roundr <- function(x, digits = 1) {
		res <- sprintf(paste0('%.', digits, 'f'), x)
		zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
		res[res == paste0('-', zzz)] <- zzz
		res
	}
	
	sapply(pvals, function(x, sig.limit) {
		if (x < sig.limit)
			if (html)
				return(sprintf('&lt; %s', format(sig.limit))) else
					return(sprintf('<%s', format(sig.limit)))
		if (x > .1)
			return(roundr(x, digits = 2)) else
				return(roundr(x, digits = digits))
	}, sig.limit = sig.limit)
}

## Functions to format tables adapted from stats.rlmerMod above by Jo Cutler 25/02/20

stats.rlmerMod.f <- function(object, level = 0.95) {
	
	# Extract beta coefficients
	beta <- fixef(object)
	
	# Extract names of coefficients
	parm <- names(beta)
	
	# Extract standard errors for the coefficients
	se <- sqrt(diag(vcov(object)))
	
	# Set level of confidence interval
	conf.level <- qnorm((1 + level) / 2)
	
	# Calculate z value
	z = beta/se
	
	# Calculate v value
	p <- 2*pnorm(-abs(z))
	
	# Calculate CI and create table
	ctab <- cbind(sprintf("%.3f", signif(beta,2)),
								sprintf("%.3f", signif(beta - (conf.level * se))), 
								sprintf("%.3f", signif(beta + (conf.level * se))),
								sprintf("%.3f", signif(se,2)),
								round(z,2),
								pvalr(p))
	
	
	# label column names
	colnames(ctab) <- c('beta',
											paste(100 * ((1 - level) / 2), '%'),
											paste(100 * ((1 + level) / 2), '%'),
											'SE',
											'z',
											'p')
	
	# Output
	return(ctab)
	
}

# ============================================================================.
# Other functions ####
# ============================================================================.
# The purpose of the functions below should be self-explanatory.
 
calculate_sample_size <- function(data, var, filter_cond) {
	# Get the unique conditions
	conditions <- unique(data[[filter_cond]])
	
	results <- list()
	
	for (condition in conditions) {
		filtered_data <- data[data[[filter_cond]] == condition, ] # subset 
		filtered_data <- filtered_data[[var]]
		
		sample_size <- length(filtered_data[!is.na(filtered_data)])
		
		results[[condition]] <- list("sample_size" = sample_size)
	}
	return(results)
}

calculate_females <- function(data, filter_cond) {
	# Get the unique conditions 
	conditions <- unique(data[[filter_cond]])
	
	results <- list()
	
	for (condition in conditions) {
		filtered_data <- data[data[[filter_cond]] == condition, ] # subset
		filtered_data <- filter(filtered_data, gender == 2)
		female_num <- nrow(filtered_data)
		
		results[[condition]] <- list("numbers" = female_num)
	}
	return(results)
}

get_summary <- function(data, var, filter_cond) {
	# Get the unique conditions
	conditions <- unique(data[[filter_cond]])
	
	results <- list()
	
	for (condition in conditions) {
		filtered_data <- data[data[[filter_cond]] == condition, ] # subset 
		filtered_data <- filtered_data[[var]]
		
		result_mean <- mean(filtered_data, na.rm = TRUE)
		result_median <- median(filtered_data, na.rm = TRUE)
		result_sd <- sd(filtered_data, na.rm = TRUE)
		
		sample_size <- length(filtered_data[!is.nan(filtered_data)])
		result_se <- result_sd / sqrt(sample_size)
		
		result_max <- max(filtered_data, na.rm = TRUE)
		result_min <- min(filtered_data, na.rm = TRUE)
		
		results[[condition]] <- list("mean" = result_mean, 
																 "median" = result_median,
																 "sd" = result_sd,
																 "se" = result_se, 
																 "max" = result_max,
																 "min" = result_min)
	}
	return(results)
}

detect_outliers <- function(data, var, filter_cond) {
	# Get the unique conditions
	conditions <- unique(data[[filter_cond]])
	
	results <- list()
	
	for (condition in conditions) {
		filtered_data <- data[data[[filter_cond]] == condition, ] # subset
		filtered_data <- filtered_data[, c("id", var, filter_cond)]
		outliers <- identify_outliers(filtered_data, all_of(var))
		
		result_outliers <- outliers[outliers$is.outlier, ]$id
		result_extremes <- outliers[outliers$is.extreme, ]$id
		
		results[[condition]] <- list("outliers" = result_outliers, 
																 "extreme_outliers" = result_extremes)
	}
	return(results)
}

check_normality <- function(data, var, filter_cond) {
	# Shapiro-Wilk Normality Test
	# Get the unique conditions 
	conditions <- unique(data[[filter_cond]])
	
	results <- list()
	
	for (condition in conditions) {
		filtered_data <- data[data[[filter_cond]] == condition, ] # subset
		result_norm <- shapiro_test(filtered_data, !!sym(var))
		
		# Get the statistic and p-value 
		result_statistic <- result_norm$statistic
		result_p <- result_norm$p
		result_p_format <- p_tidy(result_p, digits = 3)
		
		results[[condition]] <- list("statistic" = result_statistic,
																 "p" = result_p, 
																 "p_format" = result_p_format)
	}
	return(results)
}

perform_t_test <- function(data, var, filter_cond, mu = 0, alternative = "greater") {
	# Get the unique conditions
	conditions <- unique(data[[filter_cond]])
	
	results <- list()
	
	for (condition in conditions) {
		filtered_data <- data[data[[filter_cond]] == condition, ] # subset 
		result <- t_test(filtered_data, 
										 formula = as.formula(paste0(var, "~ 1")), 
										 mu = mu, alternative = alternative, detailed = TRUE)
		
		# Get the statistic and p-value
		result_statistic <- result$statistic
		result_p <- result$p
		result_p_format <- p_tidy(result_p, digits = 3)
		
		results[[condition]] <- list("result" = result,
																 "statistic" = result_statistic, 
																 "p" = result_p, 
																 "p_format" = result_p_format)
	}
	return(results)
}

# Function to test whether two numbers are of the same sign 
check_sign <- function(num1, num2) {
	if ((num1 > 0 && num2 > 0) || (num1 < 0 && num2 < 0)) {
		return(TRUE) 
	} else {
		return(FALSE)
	}
}

specify_decimal <- function(x, k) {
	ans <- format(round(x, k), nsmall = k)
	return(ans)
}