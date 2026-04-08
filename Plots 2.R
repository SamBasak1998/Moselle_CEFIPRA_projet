# -------------------------------
# HMSC plotting pipeline
# -------------------------------

library(Hmsc)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(corrplot)
library(RColorBrewer)

# ---- OUTPUT DIRECTORY ----
out_dir <- "HMSC_diagnostics_MCMC1_1"
dir.create(out_dir, showWarnings = FALSE)

# ---- Load model if new script ----
load("C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Mathematical HMSC\\HMSC_diagnostics_MCMC1_1\\HMSC_model.RData")  

# ---- FUNCTIONS ----
parse_parameter <- function(param){
  # Remove intercepts if needed
  param <- gsub("\\(Intercept\\)", "", param)
  # Extract species / covariate / trait from HMSC parameter string
  res <- str_match(param, "\\[(.*) \\(C\\d+\\), (.*) \\(S\\d+\\)\\]")[,2:3]
  return(res)
}

# ---- 1. BETA PLOTS ----
df_beta <- summary_df %>% filter(Component == "Beta") %>%
  filter(!str_detect(Parameter, "\\(Intercept\\)")) %>%
  rowwise() %>%
  mutate(tmp = list(parse_parameter(Parameter)),
         Species = tmp[2],
         Covariate = tmp[1]) %>%
  ungroup()

# 1a. Bar plot per species
p_beta_bar <- ggplot(df_beta, aes(x = Species, y = Mean, fill = Covariate)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Beta coefficients per species")
ggsave(paste0(out_dir, "/Beta_per_species.png"), p_beta_bar, width = 10, height = 6)

# 1b. Heatmap with dendrograms
beta_mat <- df_beta %>%
  select(Species, Covariate, Mean) %>%
  pivot_wider(names_from = Covariate, values_from = Mean) %>%
  column_to_rownames("Species") %>%
  as.matrix()

pheatmap(beta_mat, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Beta heatmap",
         filename = paste0(out_dir, "/Beta_heatmap.png"),
         width = 8, height = 6)

# ---- 2. GAMMA PLOTS ----
df_gamma <- summary_df %>% filter(Component == "Gamma") %>%
  rowwise() %>%
  mutate(tmp = str_match(Parameter, "\\[(.*) \\(C\\d+\\), (.*) \\(T\\d+\\)\\]")[,2:3],
         EnvCov = tmp[1],
         Trait = tmp[2]) %>%
  ungroup()

# 2a. Main heatmap (with dendrogram)
gamma_mat <- df_gamma %>%
  select(Trait, EnvCov, Mean) %>%
  pivot_wider(names_from = EnvCov, values_from = Mean, values_fn = mean) %>%
  column_to_rownames("Trait") %>%
  as.matrix()

pheatmap(gamma_mat, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Gamma heatmap",
         fontsize_col = 8,
         filename = paste0(out_dir, "/Gamma_heatmap.png"),
         width = 10, height = 8)

# 2b. Individual bar plots per environmental covariate
env_list <- unique(df_gamma$EnvCov)
for(env in env_list){
  p <- ggplot(df_gamma %>% filter(EnvCov == env), 
              aes(x = Trait, y = Mean, fill = Trait)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0("Gamma coefficients for ", env))
  ggsave(paste0(out_dir, "/Gamma_", env, ".png"), p, width = 10, height = 6)
}

# ---- 3. VARIANCE PARTITIONING ----
VP <- computeVariancePartitioning(m_fit)
VP_df <- as.data.frame(VP$vals) %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Covariate", values_to = "Mean")

p_V <- ggplot(VP_df, aes(x = Species, y = Mean, fill = Covariate)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Variance Partitioning per species")
ggsave(paste0(out_dir, "/VariancePartitioning.png"), p_V, width = 10, height = 6)

# ---- 4. SIGMA PLOTS ----
sigma_df <- as.data.frame(mpost$sigma) %>%
  pivot_longer(everything(), names_to = "Species", values_to = "Value")
p_sigma <- ggplot(sigma_df, aes(x = Species, y = Value)) +
  geom_violin(fill = "skyblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Sigma distributions per species")
ggsave(paste0(out_dir, "/Sigma_violin.png"), p_sigma, width = 10, height = 6)

# ---- 5. R² / explanatory power ----
# Compute predictions first
predY <- predict(m_fit, type = "response")

str(predY)
predY_use <- predY[[1]]
str(predY_use)
EF <- evaluateModelFit(m_fit, predY = predY_use)

EF <- evaluateModelFit(
  hM = m_fit,
  predY = predict(m_fit, expected = TRUE)[[1]]
)
R2_df <- data.frame(Species = names(EF$R2), R2 = EF$R2)
p_R2 <- ggplot(R2_df, aes(x = Species, y = R2)) +
  geom_boxplot(fill = "orange") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Explanatory power (R²) per species")
ggsave(paste0(out_dir, "/ExplanatoryPower_R2.png"), p_R2, width = 10, height = 6)

# ---- 6. OMEGA heatmap ----
OmegaCor <- computeAssociations(m_fit)
supportLevel <- 0.95

for(r in 1:m_fit$nr){
  
  plotOrder <- corrMatOrder(OmegaCor[[r]]$mean, order = "AOE")
  toPlot <- ((OmegaCor[[r]]$support > supportLevel) +
               (OmegaCor[[r]]$support < (1 - supportLevel)) > 0) *
    OmegaCor[[r]]$mean
  
  png(
    filename = paste0(out_dir, "/Omega_heatmap_rep", r, ".png"),
    width = 1600,
    height = 1400,
    res = 150
  )
  
  corrplot(
    toPlot[plotOrder, plotOrder],
    method = "color",
    col = colorRampPalette(c("blue","white","red"))(200),
    title = paste0("Omega heatmap, replicate ", r),
    tl.col = "black",
    tl.cex = 0.7
  )
  
  dev.off()
}
# ---- 7. RHO heatmap (optional) ----
rho_mat <- mpost$Rho[[1]]  # assuming 1 replicate
pheatmap(rho_mat, main="Latent factor correlations (rho)",cluster_rows = FALSE,
         cluster_cols = FALSE, filename=paste0(out_dir,"/Rho_heatmap.png"))

K <- 5  # or 4, something readable
rho_mat <- mpost$Rho[[1]][1:K, 1:K]

pheatmap(
  rho_mat,
  main = "Latent factor correlations (rho)",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  filename = paste0(out_dir, "/Rho_heatmap_top", K, ".png"),
  width = 6,
  height = 5
)

# ---- 8. ETA distributions per latent variable ----
eta_df <- as.data.frame(mpost$eta[[1]]) %>%
  pivot_longer(everything(), names_to="LatentFactor", values_to="Value")
p_eta <- ggplot(eta_df, aes(x=LatentFactor, y=Value)) +
  geom_violin(fill="lightgreen") +
  theme_minimal() +
  ggtitle("Eta distributions per latent factor")
ggsave(paste0(out_dir,"/Eta_violin.png"), p_eta, width=10, height=6)

eta_mat <- mpost$eta[[1]]

# Make sure it's matrix
eta_df <- as.data.frame(eta_mat)
colnames(eta_df) <- paste0("LF", seq_len(ncol(eta_df)))

eta_df$Site <- seq_len(nrow(eta_df))

eta_long <- eta_df %>%
  pivot_longer(
    cols = starts_with("LF"),
    names_to = "LatentFactor",
    values_to = "Value"
  )

# ---- 9. PREDICTIONS ALONG CONTINUOUS GRADIENT (example) ----
# ---- Prediction along environmental gradient (fixed) ----
library(Hmsc)
#library(HMSC)
library(tidyverse)
library(ggplot2)

# --- Parameters ---
species_sel <- "helophytes_m_2_"
covariate <- "ABL"
out_dir <- "HMSC_diagnostics_MCMC1_1"  # your folder

# --- Find species index in Y ---
species_idx <- which(colnames(m_fit$Y) == species_sel)

# --- Predict using training X ---
predY <- predict(m_fit, X = m_fit$X)

# --- Extract mean and 95% credible intervals ---
Y_mean  <- apply(predY[, , species_idx], 2, mean)
Y_lower <- apply(predY[, , species_idx], 2, quantile, probs = 0.025)
Y_upper <- apply(predY[, , species_idx], 2, quantile, probs = 0.975)

# --- Create data frame for plotting ---
df_plot <- data.frame(
  CovariateValue = m_fit$X[, covariate],
  PredMean = Y_mean,
  PredLower = Y_lower,
  PredUpper = Y_upper
)

# --- Plot ---
p <- ggplot(df_plot, aes(x = CovariateValue, y = PredMean)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper), alpha = 0.2, fill = "blue") +
  xlab(covariate) + ylab(species_sel) +
  theme_minimal() +
  ggtitle(paste("Predicted abundance for", species_sel, "vs", covariate))

# --- Save the plot ---
ggsave(filename = paste0(out_dir, "/", species_sel, "_prediction.png"),
       plot = p, width = 7, height = 5, dpi = 300)
