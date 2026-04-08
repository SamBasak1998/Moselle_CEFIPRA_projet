
## set directory where RDS files are stored
rds_dir <- "C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Mathematical HMSC\\HMSC_Outputs_new\\"

## load Gaussian
m_fit_normal  <- readRDS(paste0(rds_dir, "m_fit_normal.rds"))
mpost_normal  <- readRDS(paste0(rds_dir, "mpost_normal.rds"))

## load Poisson
m_fit_poisson <- readRDS(paste0(rds_dir, "m_fit_poisson.rds"))
mpost_poisson <- readRDS(paste0(rds_dir, "mpost_poisson.rds"))

#load log nomral poisson
#HMSC_MODEL_LP <- readRDS("HMSC_")
load(file.path(rds_dir, "HMSC_model.RData"))

m_fit_lognormal <- HMSC_model
mpost_lognormal <- convertToCodaObject(m_fit_lognormal)

# Path to your RData file
rds_dir <- "C:/Fish Data - Dr. Camara/Moselle/GTN_results_First/NB LCMM/Mathematical HMSC/HMSC_outputs/"

# Load the RData file and check what objects it contains
loaded_objs <- load(file.path(rds_dir, "HMSC_model.RData"))
print(loaded_objs)   # Shows the names of objects loaded

# Suppose the object is called "m_fit_lognormal" inside the file:
m_fit_lognormal <- get(loaded_objs)  # Assign it to a known name

# Convert to coda object for diagnostics
mpost_lognormal <- convertToCodaObject(m_fit_lognormal)


# Rename for clarity
m_fit_lognormal <- m_fit
mpost_lognormal <- mpost

library(coda)
library(dplyr)

#
# Extract MCMC chains for Beta
beta_chains <- mpost_normal$Beta

# Effective Sample Size
ess_values <- effectiveSize(beta_chains)

# Gelman-Rubin PSRF
psrf_values <- gelman.diag(beta_chains, multivariate = FALSE)$psrf

# Build convergence table
convergence_table <- data.frame(
  Parameter = names(ess_values),
  ESS = as.numeric(ess_values)
)

head(convergence_table)

# Now you can remove the old temporary names to keep environment clean
rm(m_fit, mpost)

## sanity check
class(m_fit_normal)
class(mpost_normal)
class(m_fit_poisson)
class(mpost_poisson)
#Import libraries
library(Hmsc)
library(coda)
library(ggplot2)
library(gridExtra)

# Species and covariates to inspect
species_sel <- c("ABL", "BOU", "BRB")          # adapt if names differ
cov_sel     <- c("helophytes_m_2_", "hydrophytes_m_2_", "GTN")

out_dir <- file.path(rds_dir, "convergence_plots/")
dir.create(out_dir, showWarnings = FALSE)
# -----------------------------
# Convergence Check: Beta Plots
# -----------------------------

# Required libraries
library(Hmsc)
library(ggplot2)
library(tidyr)
library(dplyr)

# Output directory for PNGs
out_dir <- "C:/HMSC_outputs/convergence/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.exists(out_dir)
# Close any open graphics devices
graphics.off()

# Load your saved Hmsc objects
# m_fit_normal, m_fit_poisson, HMSC_model (lognormal poisson)
# mpost_normal, mpost_poisson, HMSC_model$mpost (if saved)

# Species and covariates to plot
species_to_check <- c("ABL", "BOU", "BRB")        # Example species
covariates_to_check <- c("helophytes_m_2_", "hydrophytes_m_2_", "GTN")  # Example covariates

# -----------------------------
# Helper function to extract beta array
# -----------------------------
# Helper function to extract beta array from mpost object
get_beta_array <- function(mpost_obj) {
  if(!is.null(mpost_obj$Beta)) {
    return(mpost_obj$Beta)  # standard structure
  } else {
    stop("mpost object does not contain Beta array")
  }
}

get_beta_matrix <- function(mpost_obj) {
  # Unwrap the list of mcmc objects to a matrix
  if(is.list(mpost_obj$Beta)) {
    beta_mcmc <- mpost_obj$Beta[[1]]           # first (and usually only) chain
    beta_mat <- as.matrix(beta_mcmc)           # convert to plain matrix
    colnames(beta_mat) <- colnames(beta_mcmc) # keep column names
    return(beta_mat)
  } else {
    stop("mpost$Beta is not a list as expected")
  }
}
# Loop over distributions
distributions <- list(
  Gaussian = list(m_fit = m_fit_normal, mpost = mpost_normal),
  Poisson  = list(m_fit = m_fit_poisson, mpost = mpost_poisson),
  LogNormalPoisson = list(m_fit = m_fit_lognormal, mpost = mpost_lognormal)
)

str(mpost_normal)

# -----------------------------
# Close all open graphics
# -----------------------------
graphics.off()

# -----------------------------
# Loop over distributions and create trace plots
# -----------------------------
for(dist_name in names(distributions)) {
  
  cat("============================================\n")
  cat("Processing distribution:", dist_name, "\n")
  cat("============================================\n")
  
  mpost_obj <- distributions[[dist_name]]$mpost
  beta_mat  <- get_beta_matrix(mpost_obj)
  
  coef_names <- colnames(beta_mat)
  
  for(sp in species_to_check) {
    for(cv in covariates_to_check) {
      
      # Match column names exactly as in HMSC
      pattern <- paste0("\\[", cv, ".*\\(C[0-9]+\\), ", sp, ".*\\(S[0-9]+\\)\\]")
      coef_idx <- grep(pattern, coef_names)
      
      if(length(coef_idx) == 0) {
        warning(paste("Coefficient not found for species", sp, "and covariate", cv))
        next
      }
      
      beta_samples <- beta_mat[, coef_idx]
      
      # Safe filename
      safe_name <- gsub("[^A-Za-z0-9_]", "_", paste(dist_name, sp, cv, sep="_"))
      png_file <- file.path(out_dir, paste0("beta_trace_", safe_name, ".png"))
      
      # Save PNG
      png(png_file, width = 1600, height = 1200, res = 150)
      plot(beta_samples, type = "l",
           main = paste("Trace plot:", dist_name, "-", cv, "-", sp),
           xlab = "Iteration",
           ylab = expression(beta),
           col = "steelblue")
      dev.off()
      
      cat("Saved trace plot:", png_file, "\n")
    }
  }
}

cat("All beta trace plots saved as PNGs in:", out_dir, "\n")

# Libraries
library(Hmsc)
library(coda)
library(ggplot2)

# Paths
rds_dir <- "C:/HMSC_outputs/"
out_dir <- file.path(rds_dir, "convergence")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Subfolders for organization
gamma_dir     <- file.path(out_dir, "Gamma_trace")
ess_dir       <- file.path(out_dir, "ESS")
lastcloud_dir <- file.path(out_dir, "Last_cloud")
dir.create(gamma_dir, showWarnings = FALSE)
dir.create(ess_dir, showWarnings = FALSE)
dir.create(lastcloud_dir, showWarnings = FALSE)

# -----------------------------
# Define distributions
# -----------------------------
distributions <- list(
  Gaussian = list(mpost = mpost_normal),
  Poisson  = list(mpost = mpost_poisson),
  LogNormalPoisson = list(mpost = mpost_lognormal)
)

# -----------------------------
# Helper Functions
# -----------------------------

# Extract Gamma matrix from mpost object
get_gamma_matrix <- function(mpost_obj) {
  # mpost_obj$Gamma is a list of mcmc objects
  gamma_mcmc <- mpost_obj$Gamma[[1]]           # assuming only 1 chain per mpost element
  gamma_mat <- as.matrix(gamma_mcmc)
  return(gamma_mat)
}

# Calculate ESS for a matrix of samples
get_ess <- function(mat) {
  apply(mat, 2, effectiveSize)
}

# Save last cloud plot for first n_cols parameters
plot_last_cloud <- function(mat, n_cols = 5, dist_name, param_type = "Gamma", out_dir) {
  n_cols <- min(n_cols, ncol(mat))
  last_samples <- tail(mat, 5000)  # last 5000 iterations
  png_file <- file.path(out_dir, paste0("LastCloud_", param_type, "_", dist_name, ".png"))
  png(png_file, width = 1600, height = 1200, res = 150)
  pairs(last_samples[, 1:n_cols],
        main = paste("Last cloud -", param_type, "-", dist_name),
        pch = 20, cex = 0.5, col = "steelblue")
  dev.off()
  cat("Saved last cloud plot:", png_file, "\n")
}

# -----------------------------
# Loop over distributions
# -----------------------------
for(dist_name in names(distributions)) {
  cat("\n============================================\n")
  cat("Processing distribution:", dist_name, "\n")
  cat("============================================\n")
  
  mpost_obj <- distributions[[dist_name]]$mpost
  
  # --- Gamma Trace Plots ---
  gamma_mat <- get_gamma_matrix(mpost_obj)
  coef_names <- colnames(gamma_mat)
  
  for(i in seq_len(ncol(gamma_mat))) {
    beta_samples <- gamma_mat[, i]
    
    # Safe file name
    safe_name <- gsub("[^A-Za-z0-9_]", "_", paste(dist_name, coef_names[i], sep="_"))
    png_file <- file.path(gamma_dir, paste0("Gamma_trace_", safe_name, ".png"))
    
    # Trace plot
    png(png_file, width = 1600, height = 1200, res = 150)
    plot(beta_samples, type = "l",
         main = paste("Trace plot: Gamma -", coef_names[i], "-", dist_name),
         xlab = "Iteration",
         ylab = expression(Gamma),
         col = "darkorange")
    dev.off()
  }
  cat("Saved all Gamma trace plots for", dist_name, "\n")
  
  # --- ESS Boxplot ---
  gamma_ess <- get_ess(gamma_mat)
  ess_file <- file.path(ess_dir, paste0("ESS_Gamma_", dist_name, ".png"))
  png(ess_file, width = 1600, height = 1200, res = 150)
  boxplot(gamma_ess,
          main = paste("ESS for Gamma -", dist_name),
          ylab = "Effective Sample Size",
          col = "lightblue")
  dev.off()
  cat("Saved ESS boxplot for Gamma -", dist_name, "\n")
  
  # --- Last Cloud ---
  plot_last_cloud(gamma_mat, n_cols = 5, dist_name = dist_name,
                  param_type = "Gamma", out_dir = lastcloud_dir)
}

cat("\nAll Gamma trace plots, ESS boxplots, and last cloud plots saved.\n")
cat("Check folders:\n", gamma_dir, "\n", ess_dir, "\n", lastcloud_dir, "\n")

# -----------------------------
# Last Cloud for Beta
# -----------------------------
beta_dir <- file.path(out_dir, "Last_cloud_Beta")
dir.create(beta_dir, showWarnings = FALSE)

# Loop over distributions
for(dist_name in names(distributions)) {
  cat("Creating last cloud for Beta -", dist_name, "\n")
  
  mpost_obj <- distributions[[dist_name]]$mpost
  beta_mat <- as.matrix(mpost_obj$Beta[[1]])  # Convert mcmc.list to matrix
  coef_names <- colnames(beta_mat)
  
  # Plot last 5000 iterations for first 5 coefficients
  n_cols <- min(5, ncol(beta_mat))
  last_samples <- tail(beta_mat, 5000)
  
  png_file <- file.path(beta_dir, paste0("LastCloud_Beta_", dist_name, ".png"))
  png(png_file, width = 1600, height = 1200, res = 150)
  pairs(last_samples[, 1:n_cols],
        main = paste("Last cloud - Beta -", dist_name),
        pch = 20, cex = 0.5, col = "forestgreen")
  dev.off()
  
  cat("Saved last cloud plot for Beta -", dist_name, ":", png_file, "\n")
}

# -----------------------------
# ESS for Beta coefficients - all priors
# -----------------------------
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)

# Output folder
ess_dir <- file.path(out_dir, "ESS/")
dir.create(ess_dir, showWarnings = FALSE, recursive = TRUE)

# List of priors and their corresponding mpost objects
priors <- list(
  Gaussian = mpost_normal,
  Poisson = mpost_poisson,
  LogNormalPoisson = mpost_lognormal
)

for(prior_name in names(priors)) {
  
  mpost_obj <- priors[[prior_name]]
  beta_mcmc <- mpost_obj$Beta[[1]]  # mcmc object
  
  # Compute ESS for each coefficient
  ess_values <- effectiveSize(beta_mcmc)
  
  # Convert to data.frame
  ess_df <- data.frame(
    Coefficient = names(ess_values),
    ESS = as.numeric(ess_values)
  )
  
  # Extract species and covariate from coefficient names
  # Format: "B[helophytes_m_2_ (C2), ABL (S1)]"
  ess_df <- ess_df %>%
    mutate(
      Covariate = gsub("B\\[|\\s\\(C[0-9]+\\),.*\\]", "", Coefficient),
      Species   = gsub(".*,(.*)\\(S[0-9]+\\)\\]", "\\1", Coefficient),
      Species   = trimws(Species)
    )
  
  # Boxplot of ESS by Species
  p <- ggplot(ess_df, aes(x = Species, y = ESS)) +
    geom_boxplot(fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 400, linetype = "dashed", color = "red") +
    labs(
      title = paste0("Effective Sample Size (ESS) for Beta - ", prior_name),
      y = "ESS", x = "Species"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot
  ess_file <- file.path(ess_dir, paste0("ess_beta_", prior_name, ".png"))
  ggsave(ess_file, plot = p, width = 10, height = 6, dpi = 150)
  
  cat("ESS plot for Beta coefficients saved for", prior_name, "at:", ess_file, "\n")
}
