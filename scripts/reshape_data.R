# ============================================================================.
#### Info #### 
# ============================================================================.
#' Convert MATLAB arrays from .mat files into R 3D arrays, for the subsequent 
#' preference uncertainty modelling using R and Stan.
#'
#' @param pop A number to specify the age group population.
#'
#' @return A 3D R arrays.
#' @export
#' @examples
#' reshape_data(1) # for young adults
#' reshape_data(2) # for older adults
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
# 
# Project: Zhilin's first PhD project - Social discounting across the adult lifespan.
# Last revision: 20 Nov 2023, Zhilin Su.

# ============================================================================.
#### Preparation #### 
# ============================================================================.
rm(list = ls())

# Load libraries
library(R.matlab)

# ============================================================================.
#### Function definition 1 #### 
# ============================================================================.
# This function is for experiment parameters and participants' choices.
reshape_data <- function(pop) {
  
  if (pop == 1) { # 1 for young adults
    file_name <- "data/sd_ya_preference-uncertainty.mat"
  } else if (pop == 2) { # 2 for older adults
    file_name <- "data/sd_oa_preference-uncertainty.mat"
  }
  data <- readMat(file_name)

  # Extract the dimensions 
  dimension <- dim(data$data.all.array) # the third dim is sample size

  # Convert MATLAB arrays into R arrays
  data <- array(unlist(data), dim = dimension)

  # col_names 
  col_names <- readMat("data/data-variable-names.mat")
  col_names <- unlist(col_names) # convert from list into char arrays

  # row_names 
  row_names <- list() # initialise an empty list
  for (i in 1:50) {
    row_names[i] <- paste("trial", i, sep = "")
  }

  # dim_names 
  dim_names <- list() # initialise an empty list
  for (i in 1:dimension[3]) {
    dim_names[i] <- paste("pt", i, sep = "")
  }

  # Rename dimension names 
  rownames(data) <- row_names
  colnames(data) <- col_names
  dimnames(data)[[3]] <- dim_names

  return(data)
}

# ============================================================================.
#### Function definition 2 #### 
# ============================================================================.
# This function is for Others' discount rates (k = log10(K)).
reshape_k <- function(pop) {
  
  if (pop == 1) { # 1 for young adults
    file_name <- "data/sd_ya_k-values.mat"
  } else if (pop == 2) { # 2 for older adults 
    file_name <- "data/sd_oa_k-values.mat"
  }
  data <- readMat(file_name)
  data <- data$k.all.array
  
  # Extract the dimensions 
  dimension <- dim(data)
  
  # col_names 
  col_names <- readMat("data/k-variable-names.mat")
  col_names <- unlist(col_names) # convert from list into char arrays
  
  # Rename column names 
  colnames(data) <- col_names
  
  # Initialise two arrays to store the results
  other1_pref <- vector("character", length = dimension[1])
  other2_pref <- vector("character", length = dimension[1])
  
  # Specify the preferences of Others 
  ## larger k --> more impulsive; smaller k --> more patient 
  for (i in 1:dimension[1]) {
    if (data[i, "other1_k"] > data[i, "other2_k"]) { 
      other1_pref[i] <- "i"
      other2_pref[i] <- "p"
    } else if (data[i, "other1_k"] < data[i, "other2_k"]) { 
      other1_pref[i] <- "p"
      other2_pref[i] <- "i"
    }
  }
  
  # Combine results 
  data <- cbind(data, other1_pref, other2_pref)
}

# ============================================================================.
#### Run the function #### 
# ============================================================================.
sd_ya <- reshape_data(1)
sd_oa <- reshape_data(2)
k_ya <- reshape_k(1)
k_oa <- reshape_k(2)

# Save the data
save(sd_ya, sd_oa, k_ya, k_oa, file = "data/sd_data.RData")
