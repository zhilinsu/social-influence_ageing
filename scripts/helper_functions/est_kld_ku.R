# ============================================================================.
# Info #### 
# ============================================================================.
#' Estimate the Kullback-Leibler divergence (KLD) between two probability distributions (p, q), 
#' represented by two sets of mean and sd, using kernel density estimation and interpolation. 
#'
#' @param mean_p A numeric value of mean for the first distribution (p) 
#' @param sd_p A numeric value of sd for the first distribution (p)
#' @param mean_q A numeric value of mean for the second distribution (q)
#' @param sd_q A numeric value of sd for the second distribution (q)
#' @param min_range The minimum range of interpolation. Defaults to -4.
#' @param max_range The maximum range of interpolation. Defaults to 0.
#' @param n_points The number of points to use in the interpolation. Defaults to 1000.
#'
#' @return The estimated Kullback-Leibler divergence of the two distributions.
#'
#' @examples
#'  # Get samples from stanfit objects, for the parameter k  
#'  stanfit_summary_1 <- summary(stanfit1, c("km", "ku"))$summary
#'  stanfit_summary_2 <- summary(stanfit2, c("km", "ku"))$summary
#'  
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Project: Zhilin's first PhD project - Social discounting across the adult lifespan. 
# Last revision: 03 Aug 2023, Zhilin Su. 

est_kld_ku <- function(mean_p, sd_p, mean_q, sd_q, min_range = -4, max_range = 0, n_points = 10000) {
	# Set seed for reproducibility 
	set.seed(1213)
	
	# Generate random data points from a normal distribution with the specified mean 
	# and standard deviation
	data_points_p <- rnorm(n_points, mean = mean_p, sd = sd_p)
	data_points_q <- rnorm(n_points, mean = mean_q, sd = sd_q)
	
	density_p <- density(data_points_p)
	density_q <- density(data_points_q)
	
	# Create a sequence of points covering the specified range
	x_common <- seq(min_range, max_range, length.out = n_points)
	
	# Create functions that interpolate the y values of the densities,
	# rule = 2 to avoid returning NAs 
	interpolator_p <- approxfun(density_p$x, density_p$y, method = "linear", rule = 2)
	interpolator_q <- approxfun(density_q$x, density_q$y, method = "linear", rule = 2)

	# Interpolate the densities over the common range 
	prob_interpolated_p <- interpolator_p(x_common)
	prob_interpolated_q <- interpolator_q(x_common)
	
	# Replace value below 1e-10
	replace_below_threshold <- function (vec) {
		if_else(vec < 1e-10, 1e-10, vec)
	}
	
	prob_interpolated_p <- sapply(prob_interpolated_p, replace_below_threshold)
	prob_interpolated_q <- sapply(prob_interpolated_q, replace_below_threshold)
	
	# Normalise probability distributions 
	prob_interpolated_p <- prob_interpolated_p / sum(prob_interpolated_p)
	prob_interpolated_q <- prob_interpolated_q / sum(prob_interpolated_q)
	
	# Compute KL divergence 
	kl_divergence <- sum(prob_interpolated_p * log(prob_interpolated_p / prob_interpolated_q),
											 na.rm = TRUE)
	
	return(kl_divergence)
}