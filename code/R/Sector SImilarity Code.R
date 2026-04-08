library(dplyr)
library(FactoMineR)
library(cluster)
library(vegan)  # For vegdist


long_data <- read_excel(file.choose())
View(long_data)
# ---------------------------
# 1. Prepare station-level data
# ---------------------------
station_data <- long_data %>%
  select(Station, Sector,
         # Block A continuous
         AverageDepth_m_, pente_, galet, gravier, sable, Bloc, Litere, vase,
         # Block A categorical
         presence_nav_canal, pr_nav_ch,
         # Block B categorical
         fr_low, location_sp_sec, median_height, nearby_open_space, total_height) %>%
  distinct(Station, .keep_all = TRUE)
View(station_data)
# ---------------------------
# 2. Block A - Continuous
# ---------------------------
cont_A <- c("AverageDepth_m_", "pente_", "galet", "gravier", "sable", "Bloc", "Litere", "vase")
X_A_cont <- scale(station_data[, cont_A])

# PCA for weights
pca_A <- prcomp(X_A_cont, center = TRUE, scale. = TRUE)
var_explained <- summary(pca_A)$importance[3,]
npc <- which(cumsum(var_explained) >= 0.85)[1]

loadings_A <- pca_A$rotation[, 1:npc]
eigenvalues_A <- pca_A$sdev[1:npc]^2
weights_A_cont <- rowSums(sweep(loadings_A^2, 2, eigenvalues_A, "*"))
total_weight_A_cont <- sum(weights_A_cont)

# Pairwise weighted Gower distances for continuous block
D_A_cont_raw <- as.matrix(vegdist(X_A_cont, method = "gower", na.rm = TRUE))
D_A_cont <- D_A_cont_raw * total_weight_A_cont / total_weight_A_cont  # Normalize (essentially identity)

# ---------------------------
# 3. Block A - Categorical
# ---------------------------
cat_A <- c("presence_nav_canal", "pr_nav_ch")
X_A_cat <- station_data[, cat_A] %>%
  mutate(across(everything(), factor))

if(ncol(X_A_cat) > 0){
  # Check if more than one unique level exists per column
  valid_cols <- sapply(X_A_cat, function(x) nlevels(x) > 1)
  if(sum(valid_cols) > 0){
    X_A_cat <- X_A_cat[, valid_cols, drop=FALSE]
    D_A_cat_raw <- as.matrix(daisy(X_A_cat, metric = "gower"))
    weights_A_cat <- rep(1, ncol(X_A_cat))   # Equal weight if MCA not meaningful
    total_weight_A_cat <- sum(weights_A_cat)
    D_A_cat <- D_A_cat_raw * total_weight_A_cat / total_weight_A_cat  # Normalize
  } else {
    D_A_cat <- matrix(0, nrow = nrow(X_A_cont), ncol = nrow(X_A_cont))
  }
} else {
  D_A_cat <- matrix(0, nrow = nrow(X_A_cont), ncol = nrow(X_A_cont))
}

# Block A similarity
S_A <- 1 - D_A_cont
S_A
S_A[S_A < 0] <- 0  # Force [0,1]
S_A
# ---------------------------
# 4. Block B - Categorical
# ---------------------------
cat_B <- c("fr_low", "location_sp_sec", "median_height", "nearby_open_space", "total_height")
X_B_cat <- station_data[, cat_B] %>%
  mutate(across(everything(), factor))

# Only keep columns with >1 level
valid_cols_B <- sapply(X_B_cat, function(x) nlevels(x) > 1)
X_B_cat <- X_B_cat[, valid_cols_B, drop=FALSE]

if(ncol(X_B_cat) > 0){
  D_B_raw <- as.matrix(daisy(X_B_cat, metric = "gower"))
  # MCA for weights (optional)
  mca_B <- MCA(X_B_cat, graph = FALSE)
  weights_B <- rep(1, ncol(X_B_cat))  # use equal weight or mean(mca_B$var$contrib)
  total_weight_B <- sum(weights_B)
  D_B <- D_B_raw * total_weight_B / total_weight_B  # Normalize
} else {
  D_B <- matrix(0, nrow = nrow(X_A_cont), ncol = nrow(X_A_cont))
}

# Block B similarity
S_B <- 1 - D_B
S_B[S_B < 0] <- 0  # Force [0,1]

S_B
# ---------------------------
# 5. Compute block weights via MAD
# ---------------------------
V_A <- mad(as.vector(D_A_cont))
V_B <- mad(as.vector(D_B))
weight_A <- V_A / (V_A + V_B)
weight_B <- V_B / (V_A + V_B)

# ---------------------------
# 6. Compute total similarity
# ---------------------------
S_total <- weight_A * S_A + weight_B * S_B
S_total
# Force diagonals = 1
diag(S_total) <- 1

# Ensure all values in [0,1]
S_total[S_total < 0] <- 0
S_total[S_total > 1] <- 1
S_total
# ---------------------------
# 7. Optional: create vector for HMSC
# ---------------------------
between_sector_sim <- rowMeans(S_total)

View(station_data)
station_data
#Create the name order
stations <- station_data$Station
sectors  <- station_data$Sector
N <- length(stations)
N

stations <- rownames(S_total)

sector_lookup <- station_data %>%
  distinct(Station, Sector) %>%
  tibble::deframe()

#Create the similarity scores
within_sector_sim <- sapply(seq_len(N), function(i) {
  same_sector_idx <- which(sectors == sectors[i] & seq_len(N) != i)
  
  if (length(same_sector_idx) == 0) return(NA_real_)
  
  mean(S_total[i, same_sector_idx])
})

between_sector_sim <- sapply(seq_len(N), function(i) {
  other_sector_idx <- which(sectors != sectors[i] & seq_len(N) != i)
  
  if (length(other_sector_idx) == 0) return(NA_real_)
  
  mean(S_total[i, other_sector_idx])
})

within_sector_sim
between_sector_sim

#Export
similarity_df <- data.frame(
  Station = station_data$Station,
  Sector  = station_data$Sector,
  within_sector_sim  = within_sector_sim,
  between_sector_sim = between_sector_sim
)


saveRDS(
  list(
    block_weights = c(
      A = weight_A,
      B = weight_B
    ),
    weights_A_cont = weights_A_cont,
    weights_B_cat  = weights_B
  ),
  file = "block_and_variable_weights.rds"
)

library(ggplot2)
library(reshape2)

S_melt <- melt(S_total)

ggplot(S_melt, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(
    title = "Block-weighted station similarity (S_total)",
    x = "Station index",
    y = "Station index",
    fill = "Similarity"
  ) +
  theme_minimal()

sim_long <- data.frame(
  Station = station_data$Station,
  Type = rep(c("Within-sector", "Between-sector"), each = nrow(station_data)),
  Similarity = c(within_sector_sim, between_sector_sim)
)

ggplot(sim_long, aes(Type, Similarity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Within- vs Between-sector similarity",
    y = "Mean similarity"
  )

contrast_df <- data.frame(
  Station = station_data$Station,
  contrast = within_sector_sim - between_sector_sim
)

ggplot(contrast_df, aes(x = Station, y = contrast)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Within minus between sector similarity",
    y = "Similarity contrast"
  )

#Updated Heatmap
station_names <- c("AFE5", "AFE7", "ZE1", "ZE2", "ZE3")

rownames(S_total) <- station_names
colnames(S_total) <- station_names

library(ggplot2)
library(reshape2)

S_melt <- melt(S_total, varnames = c("Station_i", "Station_j"))

ggplot(S_melt, aes(Station_i, Station_j, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(limits = c(0, 1)) +
  coord_equal() +
  labs(
    title = "Block-weighted station similarity (S_total)",
    x = "Station",
    y = "Station",
    fill = "Similarity"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

S_total

#Dendogram
D_total <- as.dist(1 - S_total)
hc <- hclust(S_total, method = "average")
S_total

plot(
  hc,
  labels = station_names,
  main = "Hierarchical clustering of stations\n(Block-weighted similarity)",
  xlab = "Station",
  ylab = "Dissimilarity (1 − similarity)",
  sub = ""
)

abline(h = 0.6, col = "red", lty = 2)
