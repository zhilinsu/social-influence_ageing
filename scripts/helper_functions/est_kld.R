# ============================================================================.
#### Info #### 
# ============================================================================.
#' Estimate the Kullback-Leibler divergence (KLD) between two probability distributions (p, q), 
#' represented by two sets of samples, using kernel density estimation and interpolation. 
#'
#' @param sample_p A numeric vector of samples from the first distribution (p).
#' @param sample_q A numeric vector of samples from the second distribution (q).
#' @param min_range The minimum range of interpolation. Defaults to -4.
#' @param max_range The maximum range of interpolation. Defaults to 0.
#' @param n_points The number of points to use in the interpolation. Defaults to 1000.
#'
#' @return The estimated Kullback-Leibler divergence of the two distributions.
#'
#' @examples
#'  # Get samples from stanfit objects, for the parameter k  
#'  sample1 <- extract(stanfit1)$k
#'  sample2 <- extract(stanfit2)$k
#'  kld <- est_kld(sample1, sample2)
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Project: Zhilin's first PhD project - Social discounting across the adult lifespan. 
# Last revision: 12 Jun. 2023, Zhilin Su. 

est_kld <- function(sample_p, sample_q, min_range = -4, max_range = 0, n_points = 10000) {
	# Estimate kernel density estimations 
	density_p <- density(sample_p)
	density_q <- density(sample_q)
	
	# Create a sequence of points covering the specified range
	x_common <- seq(min_range, max_range, length.out = n_points)
	
	# Create functions that interpolate the y values of the densities,
	# rule = 2 to avoid returning NAs 
	interpolator_p <- approxfun(density_p$x, density_p$y, method = "linear", rule = 2)
	interpolator_q <- approxfun(density_q$x, density_q$y, method = "linear", rule = 2)
	
	# Interpolate the densities over the common range 
	prob_interpolated_p <- interpolator_p(x_common)
	prob_interpolated_q <- interpolator_q(x_common)
	
	# (a) Add small constant to densities to avoid division by zero 
	# epsilon <- 1e-10
	# prob_interpolated_p <- prob_interpolated_p + epsilon 
	# prob_interpolated_q <- prob_interpolated_q + epsilon 
	
	# (b) Replace value below 1e-10
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