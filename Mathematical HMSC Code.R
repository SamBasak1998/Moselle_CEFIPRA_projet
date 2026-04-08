# ========================================================
# SPATIO-TEMPORAL HMSC MODEL - MOSSELLE STUDY
# ========================================================

# Load packages
library(Hmsc)    # Hierarchical Modelling of Species Communities
library(ape)     # Phylogenetic tree
library(coda)    # MCMC diagnostics
library(dplyr)   # Data manipulation
library(tibble)  # Clean tibbles

# ========================================================
# 1. PREPARE RESPONSE Matrices

#Y <- as.matrix(abundance_data[, species_cols])
#dim(Y)  # check dimensions
Y <- read_excel(file.choose())     # response
X <- read_excel(file.choose())   # environmental covariates
#studyDesign <- read_csv(file.choose(), show_col_types = FALSE)
traits_all <- read_excel(file.choose())
View(Y)

# ========================================================
# 2. PREPARE ENVIRONMENTAL COVARIATES X
# ========================================================

# Raw environmental covariates
X_raw <- X
# Scale ONLY between_sector_sim
X_raw$between_sector_sim <- scale(X_raw$between_sector_sim)

View(X_raw)
colnames(X_raw)
# Add centered time
# Get unique years in order
unique_years <- sort(unique(X_raw$Year))  # e.g., 2014:2019

# Create a lookup table for centered time
year_lookup <- setNames(unique_years - mean(unique_years), unique_years)

# Assign centered time per row
X_raw$time_centered <- year_lookup[as.character(X_raw$Year)]
View(X_raw)
X_raw
# Continuous variables already scaled except between_sector_sim
X_scaled <- X_raw
X_scaled
# ========================================================
# 3. PREPARE TRAITS
# ========================================================
View(traits_all)
# Exclude scientific_name, code, family
trait_cols <- setdiff(colnames(traits_all), c("Scientific_name", "SP_code", "FAMILY", "Scientific_Code"))
Traits <- traits_all[, trait_cols]
row.names(Traits) <- traits_all$SP_Code # rownames = species names

# ========================================================
# 4. LOAD PHYLOGENY
# ========================================================

phylo <- read.tree("C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Mathematical HMSC\\Files\\phylogeny_tree_nwk\\moselle_fish_phylogeny.nwk")
phylo <- reorder(phylo, "postorder")  # optional, for HMSC

#Match Phylo name with Sp code:

library(ape)
library(dplyr)
library(readxl)

# -----------------------------
# 1. Load phylogeny
# -----------------------------
#phylo <- read.tree("phylogeny_19species.nwk")  # adjust path
phylo <- read.tree("C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Mathematical HMSC\\Files\\phylogeny_tree_nwk\\moselle_fish_phylogeny.nwk")
phylo
# -----------------------------
# 2. Create mapping from scientific names to SP_Code
# -----------------------------
# Fix: convert underscores to spaces in tip labels
phylo_names_with_spaces <- gsub("_", " ", phylo$tip.label)

# Map to species codes
species_mapping <- c(
  "Alburnus alburnus" = "ABL",
  "Rhodeus amarus" = "BOU",
  "Blicca bjoerkna" = "BRB",
  "Abramis brama" = "BRC",
  "Squalius cephalus" = "CHE",
  "Rutilus rutilus" = "GAR",
  "Proterorhinus semilunaris" = "GDL",
  "Neogobius fluviatilis" = "GFL",
  "Ponticola kessleri" = "GKR",
  "Gobio gobio" = "GOU",
  "Cobitis taenia" = "LOR",
  "Perca fluviatilis" = "PER",
  "Lepomis gibbosus" = "PES",
  "Pseudorasbora parva" = "PSR",
  "Scardinius erythrophthalmus" = "ROT",
  "Sander lucioperca" = "SAN",
  "Silurus glanis" = "SIL",
  "Tinca tinca" = "TAN",
  "Leuciscus leuciscus" = "VAN"
)

# Rename phylo tips to species codes
phylo$tip.label <- species_mapping[phylo_names_with_spaces]

# Check
phylo$tip.label
# Drop any tips not in Y (safety)
# Y should be loaded already with SP_Code columns
phylo <- drop.tip(phylo, setdiff(phylo$tip.label, colnames(Y)))

# -----------------------------
# 4. Prepare Traits matrix
# -----------------------------
traits_all <- read_excel(file.choose())  # choose your traits Excel file

# Keep only SP_Code as rownames, drop scientific names & code columns
Traits <- traits_all %>%
  select(-Scientific_Name, -Scientific_Code, -SP_Code, -FAMILY) %>%
  as.data.frame()

rownames(Traits) <- traits_all$SP_Code

# Optional: reorder Traits to match Y column order
Traits <- Traits[colnames(Y), ]

# -----------------------------
# 5. Check alignment
# -----------------------------
setdiff(colnames(Y), phylo$tip.label)
setdiff(phylo$tip.label, colnames(Y))

# Confirm the order of species in Y and Traits
species_order <- colnames(Y)  # should match Traits rownames

# Check tree tip labels
phylo$tip.label
# likely something like "Alburnus_alburnus", "Rhodeus_amarus", etc.

# Reassign tip labels to match species codes
phylo$tip.label <- species_order

# Now the alignment check should pass
stopifnot(all(colnames(Y) == rownames(Traits)))       # TRUE
stopifnot(all(colnames(Y) %in% phylo$tip.label))     # TRUE


stopifnot(all(colnames(Y) == rownames(Traits)))       # must be TRUE
stopifnot(all(colnames(Y) %in% phylo$tip.label))     # must be TRUE
stopifnot(all(phylo$tip.label %in% colnames(Y)))     # must be TRUE

cat("All species align across Y, Traits, and Phylogeny ✅\n")

#Optional Do a correlation matrix with phylotree

# ========================================================
# 5. STUDY DESIGN
# ========================================================
# Number of rows (observations)
n_obs <- nrow(X_raw)  # 70 in your case

# Get unique stations in the order you want
stations <- unique(X_raw$Station)

# Get unique years in order
unique_years <- sort(unique(X_raw$Year))

# Spatial factor: repeated per row in X_raw
spatial_factor <- factor(X_raw$Station)  # Each row has the correct station

# Temporal factor: centered time
year_lookup <- setNames(unique_years - mean(unique_years), unique_years)
temporal_factor <- factor(year_lookup[as.character(X_raw$Year)])

# Spatio-temporal interaction
spatiotemporal_factor <- factor(paste(X_raw$Station, X_raw$Year, sep = "_"))

# Assemble studyDesign
studyDesign <- data.frame(
  sample = 1:n_obs,
  spatial = spatial_factor,
  temporal = temporal_factor,
  spatiotemporal = spatiotemporal_factor
)

# Check
head(studyDesign)
dim(studyDesign)  # Should be 70 x 4
write.csv(studyDesign, "studyDesign_updated.csv")

studyDesign <- read.csv(file.choose())
str(studyDesign)
# ========================================================
# 6. RANDOM LEVELS
# ========================================================

# --- Spatial (Gower similarity)
S_total <- as.matrix(read.csv(file.choose(), row.names = 1))
View(S_total)

rL_spatial <- HmscRandomLevel(
  units = levels(studyDesign$spatial),
  sKnot = S_total
)

View(rL_spatial)
# Ensure row/column names of S_total match levels(studyDesign$spatial)
all(rownames(S_total) == levels(studyDesign$spatial))
all(colnames(S_total) == levels(studyDesign$spatial))

# --- Temporal (AR1)
rho_t <- 0.8
K_temporal_AR1 <- matrix(0, 6, 6)
for (k in 1:6) for (k_prime in 1:6) K_temporal_AR1[k, k_prime] <- rho_t^abs(k - k_prime)

rL_temporal <- HmscRandomLevel(
  units = levels(studyDesign$temporal),
  sKnot = K_temporal_AR1,
  sType = "Full"
)


# --- Temporal (AR1)
rho_t <- 0.8
n_time <- length(unique_years)
K_temporal_AR1 <- matrix(0, n_time, n_time)
for(k in 1:n_time) for(kp in 1:n_time) K_temporal_AR1[k,kp] <- rho_t^abs(k - kp)

rL_temporal <- HmscRandomLevel(
  units = levels(studyDesign$temporal),
  sKnot = K_temporal_AR1,
)

# --- Spatio-temporal interaction (IID)
rL_spacetime <- HmscRandomLevel(
  units = levels(studyDesign$spatiotemporal),
)

# ========================================================
# 7. DEFINE HMSC MODEL
# ========================================================
View(X_raw)
colnames(X_raw)
vars -> colnames(X_raw)

X_raw$between_sector_sim <- as.numeric(scale(X_raw$between_sector_sim))


X_scaled <- X_raw %>%
  select(-Station, -Sector, -Year) %>%
  as.data.frame()

View(X_scaled)
trait_cols <- c(
  "FFG",
  "LS",
  "HAB_RHEO",
  "HAB_TROPH",
  "HAB_LIV",
  "MOR_LMAX",
  "MOR_MIGR",
  "REP_MODE",
  "REP_GUILD",
  "REP_INCUB_PER",
  "TOL_HAB",
  "TOL_OXY"
)

Traits[] <- lapply(Traits, factor)

TrFormula <- as.formula(
  paste("~", paste(trait_cols, collapse = " + "))
)
TrFormula
stopifnot(
  all(colnames(Y) == rownames(Traits)),
  all(sapply(Traits, is.factor)),
  !anyNA(Traits)
)

studyDesign$sample <- factor(studyDesign$sample)
studyDesign$temporal <- factor(studyDesign$temporal)
studyDesign$spatial <- factor(studyDesign$spatial)
studyDesign$spatiotemporal <- factor(studyDesign$spatiotemporal)
str(studyDesign)

m <- Hmsc(
  Y = Y,
  XData = X_scaled,
  XFormula = ~ helophytes_m_2_ + hydrophytes_m_2_ + substrate_m_2_ + trees_m_2_ + GTN + between_sector_sim + time_centered,
  TrData = Traits,
  TrFormula = TrFormula,
  phyloTree = phylo,
  studyDesign = studyDesign,
  ranLevels = list(
    spatial = rL_spatial,
    temporal = rL_temporal,
    spatiotemporal = rL_spacetime
  ),
  distr = "lognormal poisson"  # default, can adjust later
)

# Check model structure
print(m)

# Check structure
str(X_scaled)

# Convert all columns to numeric if needed
X_scaled_clean <- X_scaled %>%
  mutate(across(everything(), ~ as.numeric(.)))

# Check again
str(X_scaled_clean)

library(Hmsc)
library(coda)

# --- Paths ---
out_dir <- "C:/Fish Data - Dr. Camara/Moselle/GTN_results_First/NB LCMM/Mathematical HMSC/Files/HMSC_Outputs_new/"
dir.create(out_dir, showWarnings = FALSE)

rds_dir <- "C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Mathematical HMSC\\HMSC_outputs\\"

m_fit_normal <- readRDS()
m_fit_normal  <- readRDS(paste0(rds_dir, "m_fit_normal.rds"))
# --- Model setup (update your objects) ---
# Y, X_scaled, Traits, TrFormula, phylo, studyDesign, rL_spatial, rL_temporal, rL_spacetime

distributions <- c("normal") #check parameters -> mean and variance values -> check interval of variance 

for(dist in distributions){
  
  cat("============================================\n")
  cat("Starting HMSC with distribution:", dist, "\n")
  cat("============================================\n")
  
  # Create HMSC model object
  m <- Hmsc(
    Y = Y,
    XData = X_scaled_clean,
    XFormula = ~ helophytes_m_2_ + hydrophytes_m_2_ + substrate_m_2_ + trees_m_2_ + GTN + between_sector_sim + time_centered,
    TrData = Traits,
    TrFormula = TrFormula,
    phyloTree = phylo,
    studyDesign = studyDesign,
    ranLevels = list(
      spatial = rL_spatial,
      temporal = rL_temporal,
      spatiotemporal = rL_spacetime
    ),
    distr = dist
  )
  
  cat("Sampling MCMC... \n")
  
  m_fit <- sampleMcmc(
    m,
    thin = 10,
    samples = 30000,
    transient = 5000,
    nChains = 2,
    nParallel = 1,
    verbose = 1000, 
    updater = list(GammaEta = FALSE)
  )
  
  cat("Finished MCMC for", dist, "\n")
  cat("Saving results...\n")
  
  # Save results
  saveRDS(m_fit, file = paste0(out_dir, "m_fit_", gsub(" ", "_", dist), ".rds"))
  
  # Convert to coda object for diagnostics
  mpost <- convertToCodaObject(m_fit)
  saveRDS(mpost, file = paste0(out_dir, "mpost_", gsub(" ", "_", dist), ".rds"))
  
  cat("Results saved for distribution:", dist, "\n\n")
}

cat("All distributions finished. Check convergence diagnostics on saved mpost objects.\n")

#RUn with existing m_fit#
m_fit <- sampleMcmc(
  m_fit_normal,
  samples = 30000,
  transient = 5000,
  thin = 10,
  nChains = 2,
  verbose = 1000,
  updater = list(GammaEta = FALSE)
)

##Mpost converted to Coda##
mpost_normal <- convertToCodaObject(m_fit)

ess_beta <- effectiveSize(mpost_normal$Beta)
summary(ess_beta)

psrf_beta <- gelman.diag(mpost_normal$Beta, multivariate = FALSE)$psrf
summary(psrf_beta[,1])

ess_gamma <- effectiveSize(mpost_normal$Gamma)
summary(ess_gamma)

psrf_gamma <- gelman.diag(mpost_normal$Gamma, multivariate = FALSE)$psrf
summary(psrf_gamma[,1])

mpost_normal$Omega

# Compute predicted values from the posterior
VP_model <- computePredictedValues(m_fit, expected = FALSE)

# Evaluate model fit (RMSE and R²)
VP_model <- data.frame(
  Species = colnames(VP_model),
  RMSE = evaluateModelFit(hM = m_fit, predY = VP_model)$RMSE,
  R2   = evaluateModelFit(hM = m_fit, predY = VP_model)$R2
)

# Quick summary
summary(VP_model)
VP_model   # full table


saveRDS(mpost_normal, file = paste0(out_dir, "mpost_", gsub(" ", "_", dist), ".rds"))


# ========================================================
# 8. MCMC SAMPLING PARAMETERS (run overnight)
# ========================================================

nChains <- 4
samples <- 20000
thin <- 10
transient <- 5000
nParallel <- 4

set.seed(123)

# Run MCMC (commented out for overnight run)
m_fit <- sampleMcmc(
   m,
   samples = samples,
   thin = thin,
   transient = transient,
   nChains = nChains,
   nParallel = nParallel,
   verbose = 100,
   updater = list(GammaEta = FALSE)
 )


# Check if posterior samples were actually drawn
summary(m_fit)

# Check the number of iterations per chain
m_fit$MCMC$nSamples
m_fit$MCMC$burnin
m_fit$MCMC$thin

# Check which components have posterior samples
names(m_fit$Post)

m_fit <- sampleMcmc(
  m,
  samples = samples,        # e.g., 20000
  thin = thin,              # e.g., 10
  transient = transient,    # e.g., 5000
  nChains = nChains,        # e.g., 4
  nParallel = nParallel,    # e.g., 4
  verbose = 100
  # Do NOT disable GammaEta unless you really want to skip trait-env updates
)
# Save model
# saveRDS(m_fit, "hmsc_spatiotemporal_fit.rds")

# Compute variance partitioning from HMSC model
VP <- computeVariancePartitioning(m_fit)  # or m, depending on your object
str(VP)
# Plot VP per species
plotVariancePartitioning(m_fit, VP)  # this gives the stacked bar per species

# Check posterior slots
names(m_fit$postList)

# Check Beta samples
lapply(m_fit$postList, function(x) dim(x$Beta))

# Or more directly
str(m_fit$postList, max.level = 1)

beta_post <- getPostEstimate(m_fit, "Beta")
str(beta_post)

# save in current session
save(m_fit, mpost, file = "HMSC_model.RData")

# in new script
load("HMSC_model.RData")

# ========================================================
# 9. CONVERGENCE DIAGNOSTICS
# ========================================================
library(coda)
library(lattice)

# Convert to coda object
mpost <- convertToCodaObject(m_fit)
names(mpost)
# Effective sample size for Beta
# ess_beta <- effectiveSize(mpost$Beta)
# cat("ESS Beta - Min:", min(ess_beta), "Max:", max(ess_beta), "\n")
# ESS for Beta (environmental responses)
ess_beta <- effectiveSize(mpost$Beta)

cat("ESS Beta\n")
cat("  Min:", min(ess_beta), "\n")
cat("  Median:", median(ess_beta), "\n")
cat("  Max:", max(ess_beta), "\n")

gelman_beta <- gelman.diag(
  mpost$Beta,
  multivariate = FALSE
)

cat("Max R-hat (Beta):",
    max(gelman_beta$psrf[, 1]), "\n")

##Some more plots@
# Evaluate model fit
EF <- evaluateModelFit(m_fit)
R2_df <- data.frame(Species = names(EF$R2), R2 = EF$R2)

ggplot(R2_df, aes(y = R2)) +
  geom_boxplot(fill = "steelblue") +
  theme_minimal() +
  ggtitle("Explanatory power (R²) per species") +
  ylab("R²") +
  xlab("")
ggsave(paste0(out_dir, "/R2_boxplot.png"), width = 6, height = 4)

OmegaCor <- computeAssociations(m)
supportLevel <- 0.95

for (r in 1:m$nr){
  plotOrder <- corrMatOrder(OmegaCor[[r]]$mean, order="AOE")
  toPlot <- ((OmegaCor[[r]]$support > supportLevel) +
               (OmegaCor[[r]]$support < (1-supportLevel)) > 0) * OmegaCor[[r]]$mean
  colnames(toPlot) <- rownames(toPlot) <- gsub("_"," ",x=colnames(toPlot))
  png(paste0(out_dir, "/Omega_association_run_", r, ".png"), width=1000, height=1000)
  corrplot(toPlot[plotOrder, plotOrder],
           method = "color",
           col = colorRampPalette(c("blue","white","red"))(200),
           title = "",
           type = "lower",
           tl.col = "black",
           tl.cex = 0.7,
           mar = c(0,0,6,0))
  dev.off()
}

# R-hat for Beta
# gelman_beta <- gelman.diag(mpost$Beta, multivariate = FALSE)
# cat("R-hat Beta - Max:", max(gelman_beta$psrf[,1]), "\n")

# ========================================================
# 10. EXTRACT POSTERIORS & VARIANCE PARTITIONING
# ========================================================

library(lattice)

# Beta[ covariate 1 , species 1 ]
xyplot(
  mpost$Beta[, 1],
  main = "Trace plot: Beta (covariate 1, species 1)",
  xlab = "Iteration",
  ylab = "Beta value"
)

# =========================================================
# 0. SETUP
# =========================================================

library(Hmsc)
library(coda)
library(lattice)
library(pheatmap)
library(ggplot2)

# Output directory
outdir <- "HMSC_diagnostics_MCMC1"
dir.create(outdir, showWarnings = FALSE)

# =========================================================
# 1. CONVERT TO CODA OBJECT
# =========================================================

mpost <- convertToCodaObject(m_fit)

# =========================================================
# 2. ESS & R-HAT DIAGNOSTICS
# =========================================================

ess_beta <- effectiveSize(mpost$Beta)

gelman_beta <- gelman.diag(
  mpost$Beta,
  multivariate = FALSE
)

# Save diagnostics summary
diag_summary <- data.frame(
  ESS_min    = min(ess_beta),
  ESS_median = median(ess_beta),
  ESS_max    = max(ess_beta),
  Rhat_max   = max(gelman_beta$psrf[, 1])
)

write.csv(
  diag_summary,
  file = file.path(outdir, "diagnostics_summary.csv"),
  row.names = FALSE
)

print(diag_summary)

# =========================================================
# 3. TRACE PLOTS (example Betas)
# =========================================================

pdf(file.path(outdir, "traceplots_beta.pdf"), width = 8, height = 6)

# First few beta parameters
for (i in 1:min(6, ncol(mpost$Beta))) {
  xyplot(
    mpost$Beta[, i],
    main = paste("Trace plot: Beta parameter", i),
    xlab = "Iteration",
    ylab = "Beta value"
  )
}

dev.off()

#Trace Plots
library(ggplot2)
library(dplyr)


# Directory to save individual traceplots
library(coda)
#library(ggplot2)

# Convert to coda object if not done
mpost <- convertToCodaObject(m_fit)
str(mpost)

# Extract Beta chains: list of matrices [iterations x covariates]
beta_chains <- mpost$Beta  # each element is a chain
nChains <- length(beta_chains)
cov_names <- colnames(beta_chains[[1]])
species_names <- rownames(beta_chains[[1]])

n_iter <- nrow(beta_chains[[1]])
n_cov <- ncol(beta_chains[[1]])

library(coda)

beta_post <- getPostEstimate(m_fit, "Beta")
# Set output directory
trace_dir <- "traceplots_species"
dir.create(trace_dir, showWarnings = FALSE)

# Convert to coda
mpost <- convertToCodaObject(m_fit)

# For each covariate
cov_names <- colnames(m_fit$XData)
sp_names <- m_fit$spNames

summary(mpost)
str(mpost)

summary(m_fit)
str(m_fit$postList$Beta)
# Loop over species
# Extract posterior chains for Beta
# This returns a list: each element = chain, each chain is array [covariates x species x iterations]
beta_chains_raw <- getPostEstimate(m_fit, "Beta", what="chains")

# Loop over species
for (sp_index in seq_along(sp_names)) {
  sp <- sp_names[sp_index]
  
  pdf(file.path(trace_dir, paste0("trace_", sp, ".pdf")), width=10, height=6)
  
  n_cov <- length(cov_names)
  par(mfrow=c(2, ceiling(n_cov/2)))
  
  for (cov_index in seq_along(cov_names)) {
    cov <- cov_names[cov_index]
    
    # Combine iterations from all chains
    all_iters <- do.call(c, lapply(beta_chains_raw, function(ch) {
      # ch is [covariates x species x iterations]
      # extract the iteration vector for this covariate & species
      as.vector(ch[cov_index, sp_index, ])
    }))
    
    plot(all_iters, type="l",
         main=paste(sp, "-", cov),
         xlab="Iteration", ylab="Beta")
  }
  
  dev.off()
}

cat("Trace plots saved to:", trace_dir, "\n")
# =========================================================
# 4. EXTRACT POSTERIORS
# =========================================================

beta_post  <- getPostEstimate(m_fit, "Beta")
gamma_post <- getPostEstimate(m_fit, "Gamma")

beta_mean    <- beta_post$mean
beta_support <- beta_post$support
beta_q025    <- beta_post$q025
beta_q975    <- beta_post$q975

gamma_mean <- gamma_post$mean

# Save numeric outputs
saveRDS(
  list(
    beta_mean    = beta_mean,
    beta_support = beta_support,
    beta_q025    = beta_q025,
    beta_q975    = beta_q975,
    gamma_mean   = gamma_mean
  ),
  file = file.path(outdir, "posterior_estimates.rds")
)
# -----------------------------
# HMSC Posterior Plots Script
# -----------------------------
# Input: mpost (coda object from convertToCodaObject(m_fit))
# Outputs: Trace plots, posterior densities, credible intervals, Omega heatmap, variance partitioning
# -----------------------------

# Load packages
library(coda)
library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)

outdir <- "HMSC_diagnostics_MCMC1_1"
dir.create(outdir, showWarnings = FALSE)

# -----------------------------
# 1. Posterior summaries
# -----------------------------
posterior_summary <- function(mpost_component, comp_name){
  stats <- summary(mpost_component)$statistics
  quant <- summary(mpost_component)$quantiles
  df <- data.frame(
    Parameter = rownames(stats),
    Mean = stats[,"Mean"],
    SD = stats[,"SD"],
    `2.5%` = quant[,"2.5%"],
    `97.5%` = quant[,"97.5%"]
  )
  df$Component <- comp_name
  return(df)
}

# Components to summarize
components <- c("Beta", "Gamma", "V", "Sigma")
all_post <- do.call(rbind, lapply(components, function(comp){
  posterior_summary(mpost[[comp]], comp)
}))

# Save summaries
write.csv(all_post, "HMSC_posterior_summaries.csv", row.names=FALSE)

# -----------------------------
# 2. Trace plots & posterior densities (ggplot)
# -----------------------------
# Folder to save plots
output_dir <- "HMSC_diagnostics_MCMC1_1"

# -----------------------------
# Trace plots & posterior densities (ggplot)
# -----------------------------
plot_mcmc_density <- function(mpost_component, comp_name){
  # Melt all chains
  df_list <- lapply(seq_along(mpost_component), function(i){
    chain <- as.data.frame(mpost_component[[i]])
    chain$Iteration <- 1:nrow(chain)
    chain$Chain <- paste0("Chain", i)
    melt(chain, id.vars=c("Iteration","Chain"))
  })
  df_all <- do.call(rbind, df_list)
  
  # Density plot
  p <- ggplot(df_all, aes(x=value, fill=Chain)) +
    geom_density(alpha=0.4) +
    facet_wrap(~variable, scales="free") +
    theme_bw() +
    labs(title=paste("Posterior densities for", comp_name))
  
  ggsave(filename=file.path(output_dir, paste0("Density_", comp_name, ".png")),
         plot=p, width=12, height=8)
  
  # Trace plots for first 6 parameters
  first_vars <- unique(df_all$variable)[1:min(6,length(unique(df_all$variable)))]
  df_trace <- df_all %>% filter(variable %in% first_vars)
  
  p_trace <- ggplot(df_trace, aes(x=Iteration, y=value, color=Chain)) +
    geom_line() +
    facet_wrap(~variable, scales="free_y") +
    theme_bw() +
    labs(title=paste("Trace plots for", comp_name))
  
  ggsave(filename=file.path(output_dir, paste0("Trace_", comp_name, ".png")),
         plot=p_trace, width=12, height=8)
}

# -----------------------------
# Credible interval plot
# -----------------------------
plot_credible_intervals <- function(df, comp_name){
  p <- ggplot(df, aes(x=Parameter, y=Mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=0.2) +
    coord_flip() +
    theme_bw() +
    labs(title=paste("Posterior mean & 95% CI for", comp_name))
  
  ggsave(filename=file.path(output_dir, paste0("CI_", comp_name, ".png")),
         plot=p, width=8, height=10)
}

# -----------------------------
# Omega heatmap
# -----------------------------
omega_mean <- matrix(colMeans(as.matrix(mpost$Omega[[1]])), nrow=sqrt(ncol(mpost$Omega[[1]])))
rownames(omega_mean) <- colnames(omega_mean) <- paste0("Sp", 1:ncol(omega_mean))

pheatmap(omega_mean, cluster_rows=FALSE, cluster_cols=FALSE, 
         main="Residual correlation Omega", 
         filename=file.path(output_dir, "Omega_heatmap.png"))

# -----------------------------
# Variance partitioning V plot
# -----------------------------
v_mean <- colMeans(as.matrix(mpost$V[[1]]))
v_df <- data.frame(Parameter=names(v_mean), Mean=v_mean)

p_v <- ggplot(v_df, aes(x=Parameter, y=Mean)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  theme_bw() +
  labs(title="Variance partitioning (V) per species/covariate")

ggsave(filename=file.path(output_dir, "VariancePartitioning_V.png"),
       plot=p_v, width=8, height=10)
# beta_post <- getPostEstimate(m_fit, "Beta")
# gamma_post <- getPostEstimate(m_fit, "Gamma")
# alpha_j <- beta_post$mean[, "time_centered"]
# VP <- computeVariancePartitioning(m_fit)
# print(VP$vals)
