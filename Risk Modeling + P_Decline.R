# 0. Load required libraries
library(dplyr)
library(readr)
library(tidyr)
library(Matrix)
library(Hmsc)    # Hierarchical Modelling of Species Communities
library(ape)     # Phylogenetic tree
library(coda)    # MCMC diagnostics
library(dplyr)   # Data manipulation
library(tibble)  # Clean tibbles

# Optional: for logistic regression later
library(lme4)
library(dplyr)
library(tidyr)
library(cluster)
library(lme4)
library(ggplot2)
library(reshape2)
library(readxl)
library(readr)
library(RColorBrewer)
# ------------------------------------------------------------
# 1. Import data
# ------------------------------------------------------------
#rds_dir <- "C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Mathematical HMSC\\Files for vulnerability\\HMSC_outputs\\"

## load Gaussian
#m_fit_normal  <- readRDS(paste0(rds_dir, "m_fit_normal.rds"))
#mpost_normal  <- readRDS(paste0(rds_dir, "mpost_normal.rds"))
X = read_excel(file.choose())
# Adjust directories as needed
V_total <- read.csv(file.choose())              # station × species × year, long format
GTN_density <- read_excel(file.choose())  # station × year × species
#Species_density <- read.csv("Species_density.csv")  # same format
Traits <- read_excel(file.choose())     # species × traits
pGNN <- read_csv(file.choose())                       # station × GTN spread probability
time_since_arrival <- read_excel(file.choose())  # station × year
#Y <- read.csv("Y.csv")                             # species × station × year, optional
X
Traits
# ------------------------------
# 2️⃣ Compute Trait Overlap (f_trait)
# ------------------------------
species_names <- Traits$SP_Code

species_names

Traits_matrix <- as.data.frame(Traits[,-1])
rownames(Traits_matrix) <- species_names

gower_mat <- as.matrix(daisy(Traits_matrix, metric="gower"))

GTN_species <- colnames(GTN_density)[grepl("GTN", colnames(GTN_density))]  # adjust pattern if needed
library(cluster)
#Create overlapping
# Traits: species × traits (all categorical)
Traits_matrix <- Traits
rownames(Traits_matrix) <- Traits_matrix$SP_Code
#Traits_matrix$SP_Code <- NULL   # remove species column

library(cluster)

# Convert to standard data.frame if tibble
Traits_df <- as.data.frame(Traits_matrix)

# Set rownames to SP_Code
rownames(Traits_df) <- Traits_df$SP_Code

# Remove SP_Code column (keep only traits)
Traits_df$SP_Code <- NULL

# Convert all columns to factor (Gower handles factors)
Traits_df[] <- lapply(Traits_df, as.factor)

# Compute Gower distance
gower_dist <- daisy(Traits_df, metric = "gower")

# Convert to similarity matrix if you want (1 - distance)
gower_mat <- 1 - as.matrix(gower_dist)

# Check if GTN is now a row
"GTN" %in% rownames(gower_mat)
# Should return TRUE

# Example: compute trait-based overlap with GTN species
# Compute trait-based overlap with GTN
GTN_species <- c("GTN")  # your single GTN species

f_trait <- sapply(rownames(Traits_df), function(sp) {
  max(gower_mat[sp, GTN_species], na.rm = TRUE)
})

# Normalize to 0–1
f_trait <- f_trait / max(f_trait)

# ------------------------------
# 3️⃣ Compute Spatial Overlap (f_spatial)
# ------------------------------
species_cols <- species_names

f_spatial <- "GTN"
# species columns in GTN_density
species_cols <- species_names  # your 20 species including GTN

# initialize list to store results
f_spatial <- list()

for(sp in species_cols){
  # for GTN itself, spatial overlap = 1
  if(sp == "GTN"){
    f_spatial[[sp]] <- rep(1, nrow(GTN_density))
  } else {
    # proportional abundance relative to GTN
    f_spatial[[sp]] <- GTN_density[["GTN"]] / (GTN_density[[sp]] + 1e-6)
  }
}

# convert list to data.frame
f_spatial <- as.data.frame(f_spatial)

# normalize each column to [0,1]
f_spatial <- as.data.frame(lapply(f_spatial, function(x) x / max(x, na.rm = TRUE)))

# ------------------------------
# 4️⃣ Compute I_overlap
# ------------------------------
I_overlap <- f_spatial
for(sp in species_cols){
  I_overlap[[sp]] <- I_overlap[[sp]] * f_trait[sp]
}

# Optional threshold
threshold_overlap <- 0.2
for(sp in species_cols){
  I_overlap[[sp]] <- ifelse(I_overlap[[sp]] < threshold_overlap, 0, I_overlap[[sp]])
}

# ------------------------------
# 5️⃣ Compute T_cumul
# ------------------------------
lambda <- 0.2
library(dplyr)

lambda <- 0.2

T_cumul <- time_since_arrival %>%
  mutate(
    T_cumul = exp(lambda * pmax(0, Year - first_detection))
  )


T_cumul
summary(T_cumul$T_cumul)

# ------------------------------
# 6️⃣ Compute Structured Risk Index R_ijk
# ------------------------------
str(V_total[[sp]])
str(V_total)

str(I_overlap)
View(I_overlap)

#Prepare the data frames for Risk#

library(dplyr)
library(tidyr)

I_long <- I_overlap %>%
  pivot_longer(
    cols = species_cols,
    names_to = "Species",
    values_to = "I_overlap"
  )

str(I_long)

I_overlap_fixed <- I_overlap
I_overlap_fixed$Station  <- pGNN$Station
I_overlap_fixed$Year     <- pGNN$Year
I_overlap_fixed$Campaign <- pGNN$Campaign
stopifnot(nrow(I_overlap_fixed) == nrow(pGNN))

I_overlap_long <- I_overlap_fixed %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "Species",
    values_to = "I_overlap"
  )

str(I_overlap_long)
View(I_overlap_long)
#Prepare Vulnerability#

V_long <- V_total %>%
  mutate(
    V_total = abs(V_total)
  )

str(V_long)

names(pGNN)
str(pGNN)

##################################
library(dplyr)
library(tidyr)

# ------------------------------
# Build base grid
# ------------------------------
R_ijk <- expand.grid(
  Station = unique(pGNN$Station),
  Year    = unique(pGNN$Year),
  Species = species_cols,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------
# Join all components
# ------------------------------
R_ijk <- R_ijk %>%
  left_join(
    pGNN %>% 
      select(Station, Year, GTN_spread_prob),
    by = c("Station", "Year")
  ) %>%
  left_join(
    V_total %>%
      mutate(V_total = abs(V_total)),
    by = c("Station", "Species")
  ) %>%
  left_join(
    T_cumul %>%
      select(Station, Year, T_cumul),
    by = c("Station", "Year")
  ) %>%
  left_join(
    I_overlap_long,
    by = c("Station", "Year", "Species")
  )

# ------------------------------
# Compute risk
# ------------------------------
R_ijk <- R_ijk %>%
  mutate(
    R_ijk = GTN_spread_prob * V_total * T_cumul * I_overlap
  )

str(R_ijk)
View(R_ijk)
summary(R_ijk)
table(is.na(R_ijk$R_ijk))

#Clean risk table#
R_ijk_clean <- R_ijk %>%
  filter(Species != "GTN")

table(is.na(R_ijk_clean$R_ijk))

write.csv(R_ijk, "R_ijk.csv", row.names = FALSE)

#Prepare for plotting#

# Make sure columns are numeric
R_ijk_clean$V_total <- as.numeric(R_ijk_clean$V_total)
R_ijk_clean$T_cumul <- as.numeric(R_ijk_clean$T_cumul)
R_ijk_clean$R_ijk   <- as.numeric(R_ijk_clean$R_ijk)

str(R_ijk_clean)
# =========================================
# Compute predicted probability of decline
# =========================================

library(dplyr)

# ------------------------------
# 1️⃣ Start with R_ijk_clean
# ------------------------------
data_logit <- R_ijk_clean %>%
  filter(!is.na(R_ijk)) %>%  # skip rows with NA risk (e.g., GTN)
  select(Station, Species, Year, R = R_ijk, V = V_total, T_cumul)

# ------------------------------
# 2️⃣ Scale predictors
# ------------------------------
data_logit <- data_logit %>%
  mutate(
    R_s = as.numeric(scale(R)),
    V_s = as.numeric(scale(V)),
    T_s = as.numeric(scale(T_cumul))
  )

str(data_logit)
# ------------------------------
# 3️⃣ Compute P_decline using proxy coefficients
# ------------------------------
alpha_0 <- -1
alpha_R <- -1.5
alpha_V <- -0.8
alpha_T <- 0.5

data_logit <- data_logit %>%
  mutate(
    P_decline = plogis(alpha_0 + alpha_R * R_s +
                         alpha_V * V_s +
                         alpha_T * T_s)
  )

# ------------------------------
# 4️⃣ Quick check
# ------------------------------
summary(data_logit$P_decline)
hist(
  data_logit$P_decline,
  breaks = 50,
  main = "Predicted Probability of Decline",
  xlab = "P_decline",
  col = "skyblue",
  border = "white"
)

summary(data_logit$P_decline)
# ------------------------------
# 5️⃣ Save results
# ------------------------------
write.csv(data_logit, "P_decline_per_species_site_year.csv", row.names = FALSE)

# ------------------------------
# 4️⃣ Check results
# ------------------------------

summary(data_logit$P_decline)

hist(
  data_logit$P_decline,
  breaks = 50,
  main   = "Predicted Probability of Decline",
  xlab   = "P_decline",
  col    = "skyblue",
  border = "white"
)

# ------------------------------
# 5️⃣ Save results
# ------------------------------

write.csv(data_logit, "P_decline_scaled.csv", row.names = FALSE)
library(dplyr)
library(ggplot2)

library(dplyr)
library(ggplot2)
library(fs)  # for dir_create()
library(dplyr)
library(ggplot2)
library(fs)

# ------------------------------
# 0️⃣ Load P_decline CSV
# ------------------------------
#data_logit <- read.csv("P_decline_scaled.csv", stringsAsFactors = FALSE)
data_logit <- read.csv(file.choose())
# Ensure Station is a factor with correct downstream → upstream order
station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")
station_order <- intersect(station_order, unique(data_logit$Station))
data_logit$Station <- factor(data_logit$Station, levels = station_order)

# ------------------------------
# 1️⃣ Aggregate mean P_decline per Species × Year × Station
# ------------------------------
agg_data <- data_logit %>%
  group_by(Species, Year, Station) %>%
  summarize(P_decline_mean = mean(P_decline, na.rm = TRUE), .groups = "drop")

agg_data
# ------------------------------
# 2️⃣ Output directory for PNGs
# ------------------------------
output_dir <- "species_plots_1"
fs::dir_create(output_dir)

# ------------------------------
# 3️⃣ Generate species-wise plots
# ------------------------------
species_list <- unique(agg_data$Species)

for(sp in species_list){
  sp_data <- agg_data %>% filter(Species == sp)
  
  p <- ggplot(sp_data, aes(x = Station, y = P_decline_mean, color = factor(Year), group = Year)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Set1", name = "Year") +
    labs(title = paste0("Predicted Decline Probability: ", sp),
         x = "Station (Downstream → Upstream)",
         y = "Mean P_decline") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save as PNG
  ggsave(filename = fs::path(output_dir, paste0(sp, "_P_decline.png")),
         plot = p,
         width = 7,
         height = 5,
         dpi = 300)
}

# ------------------------------
# ✅ Done
# ------------------------------
cat("All species plots saved in:", output_dir, "\n")

# ====================================================
# P_decline Functional Group Visualization Script
# ====================================================

library(dplyr)
library(ggplot2)

# ------------------------------
# 0️⃣ Prepare output folder
# ------------------------------
output_dir <- "output"
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyr)

plot_functional_group <- function(data, group_var, filename_prefix, output_dir = "output") {
  
  # Station order: downstream → upstream
  station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")
  
  # Aggregate mean P_decline
  agg_data <- data %>%
    group_by(.data[[group_var]], Station, Year) %>%
    summarize(
      P_decline_mean = mean(P_decline, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create full combination grid
  all_combinations <- expand_grid(
    !!sym(group_var) := unique(agg_data[[group_var]]),
    Station = station_order,
    Year = sort(unique(agg_data$Year))
  )
  
  # Join & RE-FACTOR Station
  agg_data <- all_combinations %>%
    left_join(agg_data, by = c(group_var, "Station", "Year")) %>%
    mutate(
      Station = factor(Station, levels = station_order)
    )
  
  # Safety check
  if (nrow(agg_data) == 0) {
    stop("Aggregation resulted in 0 rows – check inputs.")
  }
  
  # Plot
  p_line <- ggplot(
    agg_data,
    aes(
      x = Station,
      y = P_decline_mean,
      color = factor(Year),
      group = Year
    )
  ) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    facet_wrap(as.formula(paste("~", group_var))) +
    scale_x_discrete(limits = station_order) +   # 🔥 THIS FORCES ORDER
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Mean P_decline by", group_var),
      x = "Station (Downstream → Upstream)",
      y = "Mean P_decline",
      color = "Year"
    )
  
  # Output folder
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Save
  ggsave(
    filename = file.path(output_dir, paste0(filename_prefix, "_line.png")),
    plot = p_line,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("Saved: ", file.path(output_dir, paste0(filename_prefix, "_line.png")))
}

plot_functional_group(data_logit, "Habitat_Living", "P_decline_HabitatLiving")
plot_functional_group(data_logit, "Habitat_Trophic", "P_decline_HabitatTrophic")
plot_functional_group(data_logit, "Cluster_Group", "P_decline_ClusterGroup")

library(dplyr)
library(ggplot2)
library(tidyr)

plot_functional_group_ribbon <- function(data, group_var, filename_prefix,
                                         output_dir = "output") {
  
  station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")
  
  # Aggregate mean & SD
  agg_data <- data %>%
    group_by(.data[[group_var]], Station, Year) %>%
    summarize(
      mean_P = mean(P_decline, na.rm = TRUE),
      sd_P   = sd(P_decline, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Station = factor(Station, levels = station_order),
      ymin = pmax(mean_P - sd_P, 0),
      ymax = pmin(mean_P + sd_P, 1)
    )
  
  p <- ggplot(
    agg_data,
    aes(x = Station, y = mean_P, color = factor(Year), group = Year)
  ) +
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax, fill = factor(Year)),
      alpha = 0.25,
      color = NA
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    facet_wrap(as.formula(paste("~", group_var))) +
    scale_x_discrete(limits = station_order) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Mean P_decline ± SD by", group_var),
      x = "Station (Downstream → Upstream)",
      y = "P_decline",
      color = "Year",
      fill  = "Year"
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(
    file.path(output_dir, paste0(filename_prefix, "_ribbon.png")),
    p, width = 11, height = 6, dpi = 300
  )
}

plot_functional_group_ribbon(data_logit, "Habitat_Living",  "P_decline_HabitatLiving")
plot_functional_group_ribbon(data_logit, "Habitat_Trophic", "P_decline_HabitatTrophic")
plot_functional_group_ribbon(data_logit, "Cluster_Group",   "P_decline_ClusterGroup")

#Violin Plots#
plot_violin_functional <- function(data, group_var, filename_prefix,
                                   output_dir = "output") {
  
  station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")
  
  data <- data %>%
    mutate(Station = factor(Station, levels = station_order))
  
  p <- ggplot(
    data,
    aes(x = Station, y = P_decline)
  ) +
    geom_violin(fill = "grey80", color = "grey30", trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
    facet_wrap(as.formula(paste("~", group_var))) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Distribution of P_decline by", group_var),
      x = "Station (Downstream → Upstream)",
      y = "P_decline"
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(
    file.path(output_dir, paste0(filename_prefix, "_violin.png")),
    p, width = 11, height = 6, dpi = 300
  )
}

plot_violin_functional(data_logit, "Habitat_Living",  "P_decline_HabitatLiving")
plot_violin_functional(data_logit, "Habitat_Trophic", "P_decline_HabitatTrophic")
plot_violin_functional(data_logit, "Cluster_Group",   "P_decline_ClusterGroup")

#Species Specific plots#
plot_species_violin <- function(data, filename_prefix, output_dir = "output") {
  
  station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")
  
  data <- data %>%
    mutate(Station = factor(Station, levels = station_order))
  
  p <- ggplot(
    data,
    aes(x = Station, y = P_decline, fill = Species)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    facet_wrap(~ Species, scales = "free_y") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = "Species-specific distribution of P_decline",
      x = "Station",
      y = "P_decline"
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(
    file.path(output_dir, paste0(filename_prefix, "_species_violin.png")),
    p, width = 14, height = 10, dpi = 300
  )
}

plot_species_violin(data_logit, "P_decline_Species")

#Violin plot by year#
library(dplyr)
library(ggplot2)

plot_species_violin_by_year <- function(data,
                                        filename_prefix,
                                        output_dir = "output") {
  
  station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")
  
  plot_data <- data %>%
    mutate(
      Station = factor(Station, levels = station_order),
      Year    = factor(Year)
    )
  
  p <- ggplot(
    plot_data,
    aes(x = Station, y = P_decline)
  ) +
    geom_violin(
      fill = "grey80",
      color = "grey30",
      trim = TRUE
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      alpha = 0.6
    ) +
    facet_grid(
      Species ~ Year,
      scales = "free_y"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text.y = element_text(face = "bold"),
      strip.text.x = element_text(face = "bold"),
      axis.text.x  = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Temporal evolution of species-specific P_decline distributions",
      x = "Station (Downstream → Upstream)",
      y = "P_decline"
    )
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(
    file.path(output_dir, paste0(filename_prefix, "_byYear.png")),
    p,
    width = 16,
    height = 12,
    dpi = 300
  )
}

plot_species_violin_by_year(data_logit, "P_decline_Species")

plot_species_violin_per_year <- function(data,
                                         output_dir = "output") {
  
  station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")
  
  years <- sort(unique(data$Year))
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (yr in years) {
    
    plot_data <- data %>%
      filter(Year == yr) %>%
      mutate(Station = factor(Station, levels = station_order))
    
    p <- ggplot(
      plot_data,
      aes(x = Station, y = P_decline)
    ) +
      geom_violin(fill = "grey80", color = "grey30", trim = TRUE) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
      facet_wrap(~ Species, scales = "free_y") +
      scale_y_continuous(limits = c(0, 1)) +
      theme_minimal(base_size = 11) +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = paste("Species-specific P_decline distribution – Year", yr),
        x = "Station (Downstream → Upstream)",
        y = "P_decline"
      )
    
    ggsave(
      file.path(output_dir, paste0("P_decline_Species_", yr, ".png")),
      p,
      width = 14,
      height = 10,
      dpi = 300
    )
  }
}

plot_species_violin_per_year(data_logit)

#Violin Plot per year - split by species groups#
station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")

data_logit <- data_logit %>%
  mutate(
    Station = factor(Station, levels = station_order),
    Year    = factor(Year)
  )

#Define palettes#
palette_habitat_living <- c(
  "Benthic"       = "#1f4e79",
  "Semi-Benthic"  = "#2a9d8f",
  "Non-Benthic"   = "#e76f51"
)

palette_habitat_trophic <- c(
  "Benthic"     = "#2a6f3e",
  "Non-Benthic" = "#6a4c93"
)

palette_cluster <- c(
  "Cluster1" = "#457b9d",
  "Cluster2" = "#c1121f",
  "Cluster3" = "#f4a261"
)

library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)

# ==============================
# Libraries
# ==============================
library(dplyr)
library(ggplot2)
library(ggridges)
library(forcats)
install.packages("ggridges")
# ==============================
# Station order (Downstream → Upstream)
# ==============================
station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")

# ==============================
# 1) READABLE SPECIES VIOLIN PLOTS (per year)
# ==============================
plot_species_violin_by_year <- function(data,
                                        group_var,
                                        palette,
                                        filename_prefix,
                                        output_dir = "output") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  data <- data %>%
    mutate(Station = factor(Station, levels = station_order))
  
  years <- sort(unique(data$Year))
  
  for (yr in years) {
    
    plot_data <- data %>%
      filter(Year == yr, !is.na(.data[[group_var]]))
    
    p <- ggplot(
      plot_data,
      aes(
        x = Station,
        y = P_decline,
        fill = .data[[group_var]]
      )
    ) +
      geom_violin(
        alpha = 0.85,
        color = "grey30",
        trim = TRUE
      ) +
      geom_boxplot(
        width = 0.10,
        outlier.shape = NA,
        alpha = 0.6,
        color = "black"
      ) +
      stat_summary(
        fun = median,
        geom = "point",
        size = 1.8,
        color = "black"
      ) +
      facet_wrap(~ Species, scales = "free_y") +
      scale_fill_manual(values = palette) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste("Species-level P_decline distributions – Year", yr),
        subtitle = paste("Median (●), distribution per station | Group:", group_var),
        x = "Station (Downstream → Upstream)",
        y = "P_decline",
        fill = group_var
      )
    
    ggsave(
      file.path(output_dir, paste0(filename_prefix, "_", yr, ".png")),
      plot = p,
      width = 14,
      height = 8,
      dpi = 300
    )
  }
}

# ==============================
# 2) RIDGELINE DISTRIBUTION PLOT
# ==============================
plot_species_ridgeline <- function(data,
                                   group_var,
                                   palette,
                                   filename,
                                   output_dir = "output") {
  
  plot_data <- data %>%
    filter(!is.na(P_decline)) %>%
    mutate(
      Station = factor(Station, levels = station_order),
      Species = fct_reorder(Species, P_decline, .fun = median)
    )
  
  p <- ggplot(
    plot_data,
    aes(
      x = P_decline,
      y = Species,
      fill = .data[[group_var]]
    )
  ) +
    geom_density_ridges(
      alpha = 0.85,
      scale = 1.1,
      color = "white",
      rel_min_height = 0.01
    ) +
    facet_wrap(~ Station, nrow = 1) +
    scale_fill_manual(values = palette) +
    theme_minimal(base_size = 12) +
    labs(
      title = "Species-specific P_decline distributions across stations",
      subtitle = "Ridgeline densities emphasize variability and uncertainty",
      x = "P_decline",
      y = "Species",
      fill = group_var
    )
  
  ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    p,
    width = 14,
    height = 8,
    dpi = 300
  )
}

# ==============================
# 3) MEDIAN SPATIAL TREND (± IQR)
# ==============================
plot_species_median_trend <- function(data,
                                      group_var,
                                      palette,
                                      filename,
                                      output_dir = "output") {
  
  summary_data <- data %>%
    filter(!is.na(P_decline)) %>%
    mutate(Station = factor(Station, levels = station_order)) %>%
    group_by(Species, Station, .data[[group_var]]) %>%
    summarize(
      median_P = median(P_decline),
      q25 = quantile(P_decline, 0.25),
      q75 = quantile(P_decline, 0.75),
      .groups = "drop"
    )
  
  p <- ggplot(
    summary_data,
    aes(
      x = Station,
      y = median_P,
      group = Species,
      color = .data[[group_var]]
    )
  ) +
    geom_line(alpha = 0.75, linewidth = 0.9) +
    geom_point(size = 2) +
    geom_ribbon(
      aes(ymin = q25, ymax = q75, fill = .data[[group_var]]),
      alpha = 0.2,
      color = NA
    ) +
    facet_wrap(~ Species, scales = "free_y") +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    theme_minimal(base_size = 11) +
    labs(
      title = "Spatial trends in species vulnerability",
      subtitle = "Median ± IQR across stations",
      x = "Station (Downstream → Upstream)",
      y = "Median P_decline",
      color = group_var,
      fill = group_var
    )
  
  ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    p,
    width = 14,
    height = 10,
    dpi = 300
  )
}

data_logit
plot_species_violin_by_year(
  data_logit,
  group_var = "Cluster_Group",
  filename_prefix = "P_decline_Cluster",
  palette = palette_cluster
)

plot_species_ridgeline(
  data_logit,
  "Habitat_Living",
  palette_habitat_living,
  "Ridgeline_P_decline_HabitatLiving"
)

plot_species_median_trend(
  data_logit,
  "Habitat_Living",
  palette_habitat_living,
  "MedianTrend_P_decline_HabitatLiving"
)

plot_species_ridgeline_by_year <- function(data,
                                           group_var,
                                           palette,
                                           filename_prefix,
                                           output_dir = "output") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  data <- data %>%
    mutate(Station = factor(Station, levels = c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")))
  
  years <- sort(unique(data$Year))
  
  for (yr in years) {
    
    plot_data <- data %>%
      filter(Year == yr, !is.na(.data[[group_var]])) %>%
      mutate(
        Species = forcats::fct_reorder(Species, P_decline, .fun = median)
      )
    
    p <- ggplot(
      plot_data,
      aes(
        x = P_decline,
        y = Species,
        fill = .data[[group_var]]
      )
    ) +
      geom_density_ridges(
        alpha = 0.85,
        scale = 1.2,
        color = "white",
        rel_min_height = 0.01
      ) +
      facet_wrap(~ Station, nrow = 1) +
      scale_fill_manual(values = palette) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste("Species P_decline distributions – Year", yr),
        subtitle = paste("Ridgeline by station | Group:", group_var),
        x = "P_decline",
        y = "Species",
        fill = group_var
      )
    
    ggsave(
      file.path(output_dir, paste0(filename_prefix, "_", yr, ".png")),
      p,
      width = 14,
      height = 8,
      dpi = 300
    )
  }
}

plot_species_median_trend_by_year <- function(data,
                                              group_var,
                                              palette,
                                              filename_prefix,
                                              output_dir = "output") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  data <- data %>%
    mutate(Station = factor(Station, levels = c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")))
  
  years <- sort(unique(data$Year))
  
  for (yr in years) {
    
    summary_data <- data %>%
      filter(Year == yr, !is.na(.data[[group_var]])) %>%
      group_by(Species, Station, .data[[group_var]]) %>%
      summarize(
        median_P = median(P_decline),
        q25 = quantile(P_decline, 0.25),
        q75 = quantile(P_decline, 0.75),
        .groups = "drop"
      )
    
    p <- ggplot(
      summary_data,
      aes(
        x = Station,
        y = median_P,
        group = Species,
        color = .data[[group_var]]
      )
    ) +
      geom_line(alpha = 0.8, linewidth = 0.9) +
      geom_point(size = 2) +
      geom_ribbon(
        aes(ymin = q25, ymax = q75, fill = .data[[group_var]]),
        alpha = 0.25,
        color = NA
      ) +
      facet_wrap(~ Species, scales = "free_y") +
      scale_color_manual(values = palette) +
      scale_fill_manual(values = palette) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = paste("Spatial trends in P_decline – Year", yr),
        subtitle = paste("Median ± IQR | Group:", group_var),
        x = "Station (Downstream → Upstream)",
        y = "Median P_decline",
        color = group_var,
        fill = group_var
      )
    
    ggsave(
      file.path(output_dir, paste0(filename_prefix, "_", yr, ".png")),
      p,
      width = 14,
      height = 9,
      dpi = 300
    )
  }
}

# Habitat Living
plot_species_ridgeline_by_year(
  data_logit,
  "Cluster_Group",
  palette_cluster,
  "Ridgeline_P_decline_Cluster"
)

#Median Trend by Year#
plot_species_median_trend_by_year(
  data_logit,
  "Cluster_Group",
  palette_cluster,
  "MedianTrend_P_decline_Cluster"
)

#Summary Plots#
# ------------------------------
# Condensed Summary Plots: Ridgeline + Median per Functional Group
# ------------------------------
library(dplyr)
library(ggplot2)
library(ggridges)
library(patchwork)  # for combining plots

plot_summary_functional_group <- function(data, group_var, palette, filename_prefix, output_dir="output") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Ensure Station order
  station_order <- c("ZE1", "ZE2", "ZE3", "AFE5", "AFE7")
  data <- data %>%
    mutate(Station = factor(Station, levels = station_order))
  
  years <- sort(unique(data$Year))
  
  # Aggregate median + IQR per species × station × year
  summary_data <- data %>%
    group_by(Species, Station, Year, .data[[group_var]]) %>%
    summarise(
      median_P = median(P_decline, na.rm = TRUE),
      Q25 = quantile(P_decline, 0.25, na.rm = TRUE),
      Q75 = quantile(P_decline, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Median trend per station across species (average median)
  median_trend <- summary_data %>%
    group_by(Station, Year, .data[[group_var]]) %>%
    summarise(
      median_P_station = median(median_P, na.rm = TRUE),
      Q25_station = quantile(median_P, 0.25, na.rm = TRUE),
      Q75_station = quantile(median_P, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Loop through years to create one figure per year
  for (yr in years) {
    
    plot_data <- summary_data %>% filter(Year == yr)
    trend_data <- median_trend %>% filter(Year == yr)
    
    # Ridgeline plot: species distributions
    p_ridge <- ggplot(plot_data, aes(x = median_P, y = Species, fill = .data[[group_var]])) +
      geom_density_ridges(alpha = 0.85, scale = 1) +
      scale_fill_manual(values = palette) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "P_decline",
        title = paste("Ridgeline: Species distributions – Year", yr),
        fill = group_var
      )
    
    # Median trend plot: station-wise
    p_median <- ggplot(trend_data, aes(x = Station, y = median_P_station, group = .data[[group_var]], color = .data[[group_var]])) +
      geom_line(size = 1.2) +
      geom_point(size = 2) +
      geom_ribbon(aes(ymin = Q25_station, ymax = Q75_station, fill = .data[[group_var]]), alpha = 0.25, color = NA) +
      scale_color_manual(values = palette) +
      scale_fill_manual(values = palette) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Station (Downstream → Upstream)",
        y = "Median P_decline",
        title = paste("Median ± IQR trend – Year", yr),
        color = group_var,
        fill = group_var
      )
    
    # Combine plots side by side
    combined <- p_ridge + p_median + plot_layout(widths = c(1, 1.2))
    
    # Save figure
    ggsave(
      filename = file.path(output_dir, paste0(filename_prefix, "_Year_", yr, ".png")),
      plot = combined,
      width = 16,
      height = 8,
      dpi = 300
    )
  }
  
  message("Summary figures saved in folder: ", output_dir)
}

# ------------------------------
# Example usage
# ------------------------------
palette_habitat_living <- c(
  "Benthic"       = "#1f4e79",
  "Semi-Benthic"  = "#2a9d8f",
  "Non-Benthic"   = "#e76f51"
)

palette_habitat_trophic <- c(
  "Benthic"     = "#2a6f3e",
  "Non-Benthic" = "#6a4c93"
)

palette_cluster <- c(
  "Cluster1" = "#457b9d",
  "Cluster2" = "#c1121f",
  "Cluster3" = "#f4a261"
)

# Run for all three functional groupings
plot_summary_functional_group(data_logit, "Habitat_Living", palette_habitat_living, "Summary_HabitatLiving")
plot_summary_functional_group(data_logit, "Habitat_Trophic", palette_habitat_trophic, "Summary_HabitatTrophic")
plot_summary_functional_group(data_logit, "Cluster_Group", palette_cluster, "Summary_ClusterGroup")

