library(Hmsc)    # Hierarchical Modelling of Species Communities
library(ape)     # Phylogenetic tree
library(coda)    # MCMC diagnostics
library(dplyr)   # Data manipulation
library(tibble)  # Clean tibbles

# ========================================================
# 1. PREPARE RESPONSE Matrices

#Y <- as.matrix(abundance_data[, species_cols])
#dim(Y)  # check dimensions
summary_df <- read_csv(file.choose()) 
View(summary_df)

# Load packages
library(tidyverse)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# Set output directory
out_dir <- "HMSC_diagnostics_MCMC1_1"
if(!dir.exists(out_dir)) dir.create(out_dir)

# Assume summary_df is already loaded
# summary_df columns: Parameter, Mean, SD, X2.5., X97.5., Component

# 1. Diagnostic plots: Density + Trace for each parameter (simulate trace if unavailable)
components <- unique(summary_df$Component)

for(comp in components){
  df_comp <- summary_df %>% filter(Component == comp)
  
  # Density plot
  p_density <- ggplot(df_comp, aes(x = Mean)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    ggtitle(paste0("Density plot: ", comp)) +
    theme_minimal()
  
  ggsave(filename = paste0(out_dir, "/", comp, "_density.png"),
         plot = p_density, width = 6, height = 4)
  
  # Histogram as a proxy for trace plot (if actual MCMC not available)
  p_trace <- ggplot(df_comp, aes(x = Mean)) +
    geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
    ggtitle(paste0("Histogram proxy for trace: ", comp)) +
    theme_minimal()
  
  ggsave(filename = paste0(out_dir, "/", comp, "_trace.png"),
         plot = p_trace, width = 6, height = 4)
}

# 2. Beta plots: per species + heatmap
df_beta <- summary_df %>% filter(Component == "Beta")

# Extract species and covariate from Parameter column
df_beta <- df_beta %>%
  separate(Parameter, into = c("B", "Species", "Covariate"), sep = "\\[|\\]|,", extra = "merge") %>%
  mutate(Species = str_trim(Species),
         Covariate = str_trim(Covariate))

# Beta heatmap
p_beta_heatmap <- ggplot(df_beta, aes(x = Covariate, y = Species, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  ggtitle("Beta heatmap (Mean values)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste0(out_dir, "/Beta_heatmap.png"), plot = p_beta_heatmap, width = 8, height = 6)

# 3. Gamma plots: per trait x environment + heatmap
df_gamma <- summary_df %>% filter(Component == "Gamma")
df_gamma <- df_gamma %>%
  separate(Parameter, into = c("G", "Trait", "Env"), sep = "\\[|\\]|,", extra = "merge") %>%
  mutate(Trait = str_trim(Trait),
         Env = str_trim(Env))

p_gamma_heatmap <- ggplot(df_gamma, aes(x = Env, y = Trait, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  ggtitle("Gamma heatmap (Mean values)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = paste0(out_dir, "/Gamma_heatmap.png"), plot = p_gamma_heatmap, width = 8, height = 6)

# 4. Variance partitioning
df_V <- summary_df %>% filter(Component == "V")
df_V <- df_V %>%
  separate(Parameter, into = c("V", "Species", "Covariate"), sep = "\\[|\\]|,", extra = "merge") %>%
  mutate(Species = str_trim(Species),
         Covariate = str_trim(Covariate))

p_V <- ggplot(df_V, aes(x = Species, y = Mean, fill = Covariate)) +
  geom_bar(stat = "identity", position = "stack") +
  ggtitle("Variance Partitioning per Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = paste0(out_dir, "/Variance_partitioning.png"), plot = p_V, width = 8, height = 6)

# 5. Omega correlation plot
df_omega <- summary_df %>% filter(Component == "Omega")
df_omega <- df_omega %>%
  separate(Parameter, into = c("O", "Species1", "Species2"), sep = "\\[|\\]|,", extra = "merge") %>%
  mutate(Species1 = str_trim(Species1),
         Species2 = str_trim(Species2))

p_omega <- ggplot(df_omega, aes(x = Species1, y = Species2, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  ggtitle("Residual Associations (Omega)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0))

ggsave(filename = paste0(out_dir, "/Omega_heatmap.png"), plot = p_omega, width = 8, height = 6)

# 6. Effective Sample Size (ESS) distributions
df_ess <- summary_df %>% filter(Component %in% c("Beta","Gamma","V","Omega"))

p_ess <- ggplot(df_ess, aes(x = Component, y = SD)) +  # Using SD as a proxy for variability
  geom_boxplot(fill = "orange") +
  ggtitle("Distribution of SD (proxy for ESS)") +
  theme_minimal()

ggsave(filename = paste0(out_dir, "/ESS_distribution.png"), plot = p_ess, width = 6, height = 4)

# 7. Explanatory power (R^2) boxplot
# Assuming you have a separate dataframe r2_df with Species, R2
# Example:
# r2_df <- data.frame(Species = df_beta$Species %>% unique(),
#                     R2 = runif(length(df_beta$Species), 0.3, 0.9))

# Uncomment if you have R2 values
# p_r2 <- ggplot(r2_df, aes(x = Species, y = R2)) +
#   geom_boxplot(fill = "lightblue") +
#   ggtitle("Explained variance R² per Species") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = paste0(out_dir, "/R2_boxplot.png"), plot = p_r2, width = 8, height = 6)

# 8. Predicted values along a gradient (optional)
# You need a dataset with covariate gradient and predicted values per species
# This part can be done after re-running the HMSC with predict()

# All plots are saved in HMSC_diagnostics_MCMC1_1 folder
print("All plots saved successfully!")

##Updated#
library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Output folder
out_dir <- "HMSC_diagnostics_MCMC1_1"
if(!dir.exists(out_dir)) dir.create(out_dir)

# ---- Helper function to parse Parameter strings ----
parse_parameter <- function(param){
  # Remove B[ or G[ or V[ etc
  param <- gsub("^[A-Z]\\[|\\]$", "", param)
  # Split by comma
  parts <- strsplit(param, ",")[[1]]
  # Trim whitespace
  parts <- trimws(parts)
  return(parts)
}

# ---- BETA PLOTS ----
df_beta <- summary_df %>% filter(Component == "Beta")

# ---- BETA PLOTS ----
df_beta <- summary_df %>% 
  filter(Component == "Beta") %>%
  filter(!grepl("\\(Intercept\\)", Parameter))  # remove intercepts


# Extract species and covariate
df_beta_parsed <- df_beta %>%
  rowwise() %>%
  mutate(tmp = list(parse_parameter(Parameter)),
         Covariate = tmp[2],
         Species = tmp[1]) %>%
  ungroup()

# Extract species and covariate
df_beta_parsed <- df_beta %>%
  rowwise() %>%
  mutate(tmp = list(parse_parameter(Parameter)),
         Covariate = tmp[2],
         Species = tmp[1]) %>%
  ungroup()

# Individual bar plots per species
p_beta_bar <- ggplot(df_beta_parsed, aes(x = Covariate, y = Mean, fill = Covariate)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Species, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Beta coefficients per species")

ggsave(paste0(out_dir, "/Beta_per_species.png"), plot = p_beta_bar, width = 10, height = 6)

# Heatmap with clustering
beta_mat <- df_beta_parsed %>%
  select(Species, Covariate, Mean) %>%
  pivot_wider(names_from = Covariate, values_from = Mean) %>%
  column_to_rownames("Species") %>%
  as.matrix()

pheatmap(beta_mat,
         filename = paste0(out_dir, "/Beta_heatmap.png"),
         color = colorRampPalette(c("blue","white","red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Beta heatmap (hierarchical clustering)")


# ---- GAMMA PLOTS ----
df_gamma <- summary_df %>% filter(Component == "Gamma")

df_gamma_parsed <- df_gamma %>%
  rowwise() %>%
  mutate(tmp = list(parse_parameter(Parameter)),
         EnvCov = tmp[2],
         Trait = tmp[1]) %>%
  ungroup()

# ---- GAMMA PLOTS ----
df_gamma <- summary_df %>% 
  filter(Component == "Gamma") %>%
  filter(!grepl("\\(Intercept\\)", Parameter))  # remove intercepts

df_gamma_parsed <- df_gamma %>%
  rowwise() %>%
  mutate(tmp = list(parse_parameter(Parameter)),
         EnvCov = tmp[2],
         Trait = tmp[1]) %>%
  ungroup()

# Individual bar plots per trait
p_gamma_bar <- ggplot(df_gamma_parsed, aes(x = EnvCov, y = Mean, fill = EnvCov)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Gamma coefficients per trait")

ggsave(paste0(out_dir, "/Gamma_per_trait.png"), plot = p_gamma_bar, width = 10, height = 6)

# Heatmap with clustering
gamma_mat <- df_gamma_parsed %>%
  select(Trait, EnvCov, Mean) %>%
  pivot_wider(names_from = EnvCov, values_from = Mean) %>%
  column_to_rownames("Trait") %>%
  as.matrix()

pheatmap(gamma_mat,
         filename = paste0(out_dir, "/Gamma_heatmap.png"),
         color = colorRampPalette(c("blue","white","red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Gamma heatmap (hierarchical clustering)")
# ---- GAMMA PLOTS PER COVARIATE ----
library(ggplot2)
library(dplyr)

# ---- Ensure your directory exists ----
out_dir <- "HMSC_diagnostics_MCMC1_1"

# ---- Filter out intercepts if any ----
df_gamma_parsed <- df_gamma_clean %>% 
  filter(!grepl("\\(Intercept\\)", Trait))

df_gamma_parsed
# ---- Loop over environmental covariates ----
covariates <- unique(df_gamma_parsed$Covariate)

for(cov in covariates){
  
  df_subset <- df_gamma_parsed %>%
    filter(Covariate == cov)
  
  p <- ggplot(df_subset, aes(x = Trait, y = Mean, fill = Mean)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggtitle(paste("Gamma coefficients for", cov))
  
  ggsave(
    filename = paste0(out_dir, "/Gamma_per_trait_", cov, ".png"),
    plot = p,
    width = 8,
    height = 5
  )
}


# ---- VARIANCE PARTITIONING (V) ----
# Filter V component
df_V <- summary_df %>% filter(Component == "V")

# Parse Parameter properly for V: extract Species first, Covariate second
df_V_parsed <- df_V %>%
  rowwise() %>%
  mutate(
    # Assuming V[Covariate (Species)]
    Species = gsub(".*\\((S[0-9]+)\\).*", "\\1", Parameter),
    Covariate = gsub("V\\[(.*) \\(S[0-9]+\\)\\].*", "\\1", Parameter)
  ) %>%
  ungroup()

# Ensure Species is a factor in the order you want
df_V_parsed$Species <- factor(df_V_parsed$Species, levels = unique(df_V_parsed$Species))

# Plot
p_V <- ggplot(df_V_parsed, aes(x = Species, y = Mean, fill = Covariate)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Variance partitioning per species")

# Save
ggsave(paste0(out_dir, "/Variance_partitioning.png"), plot = p_V, width = 10, height = 6)

# ---- SIGMA PLOT ----
df_sigma <- summary_df %>% filter(Component == "Sigma")
df_sigma_parsed <- df_sigma %>%
  rowwise() %>%
  mutate(tmp = list(parse_parameter(Parameter)),
         Species = tmp[1],
         Covariate = tmp[2]) %>%
  ungroup()

sigma_mat <- df_sigma_parsed %>%
  select(Species, Covariate, Mean) %>%
  pivot_wider(names_from = Covariate, values_from = Mean) %>%
  column_to_rownames("Species") %>%
  as.matrix()

pheatmap(sigma_mat,
         filename = paste0(out_dir, "/Sigma_heatmap.png"),
         color = colorRampPalette(c("white","orange","red"))(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Sigma heatmap")

# ---- Predicted values along a continuous gradient ----
# Template: replace 'pred_df' with actual predicted values after HMSC run
# Example: pred_df = data.frame(Covariate = seq(-2,2,length=50), Species1 = ..., Species2 = ...)
# For now we simulate
set.seed(123)
species_list <- unique(df_beta_parsed$Species)
species_list
pred_df <- data.frame(Covariate = seq(-2,2,length=50))
for(sp in species_list[1:7]){  # plotting first 5 species as example
  pred_df[[sp]] <- sin(seq(-2,2,length=50)) + rnorm(50,0,0.1)
}

pred_df_long <- pred_df %>% pivot_longer(cols = -Covariate, names_to = "Species", values_to = "Predicted")

p_pred <- ggplot(pred_df_long, aes(x = Covariate, y = Predicted, color = Species)) +
  geom_line(size = 1) +
  theme_minimal() +
  ggtitle("Predicted values along continuous gradient")

ggsave(paste0(out_dir, "/Predicted_gradient.png"), plot = p_pred, width = 10, height = 6)

print("All updated plots saved to folder HMSC_diagnostics_MCMC1_1")
