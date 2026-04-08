library ( ape )

# 1. Load phylogenetic tree ( Newick format ) 
phylo <- read.tree ("C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Mathematical HMSC\\Files\\phylogeny_tree_nwk\\moselle_fish_phylogeny.nwk")
# 2. Compute patristic distances ( sum of branch lengths ) 

D_phylo <- cophenetic(phylo)
rho <- 0.1	# Decay parameter ( to calibrate ) 
C <- exp (-rho*D_phylo )
# 4. Normalize diagonal 
diag(C) <- 1

plot(phylo, type = "fan", main = "Moselle Fish Phylogeny (Fan Style)", cex = 0.7)

# ========================================
# Phylogenetic plots for Moselle fish
# ========================================

# Load libraries
library(ape)
library(pheatmap)
library(gplots)
library(ggtree)
library(ggplot2)
# ========================================
# Phylogenetic plots for Moselle fish
# ========================================

# Load libraries
library(ape)
library(pheatmap)
library(gplots)

# Output directory
out_dir <- "C:/Fish Data - Dr. Camara/Moselle/GTN_results_First/NB LCMM/Mathematical HMSC/HMSC_diagnostics_MCMC1_1/"

# ---- 1. Load phylogenetic tree ----
phylo <- read.tree("C:/Fish Data - Dr. Camara/Moselle/GTN_results_First/NB LCMM/Mathematical HMSC/Files/phylogeny_tree_nwk/moselle_fish_phylogeny.nwk")

# ---- 2. Compute phylogenetic correlation matrix C ----
D_phylo <- cophenetic(phylo)
rho <- 0.1
C <- exp(-rho * D_phylo)
diag(C) <- 1

write.csv(D_phylo, "Phylo Matrix.csv")
# ---- 3. Basic tree plot (rectangular) ----
png(paste0(out_dir, "phylo_tree_rectangular.png"), width = 1200, height = 800)
plot(phylo, main = "Moselle Fish Phylogeny (Rectangular)", cex = 0.8)
dev.off()

# ---- 4. Fan-style tree ----
png(paste0(out_dir, "phylo_tree_fan.png"), width = 1200, height = 800)
plot(phylo, type = "fan", main = "Moselle Fish Phylogeny (Fan Style)", cex = 0.7)
dev.off()

# ---- 5. Example tree colored by functional group / risk ----
set.seed(123)
groups <- sample(1:3, length(phylo$tip.label), replace = TRUE) # replace with real risk/vuln scores
tip_colors <- c("red", "blue", "green")[groups]

png(paste0(out_dir, "phylo_tree_colored.png"), width = 1200, height = 800)
plot(phylo, tip.color = tip_colors, main = "Phylogeny Colored by Group", cex = 0.8)
legend("topright", legend = paste("Group", 1:3), col = c("red", "blue", "green"), pch = 19)
dev.off()

# ---- 6. Heatmap of phylogenetic correlation matrix ----
pheatmap(C,
         clustering_method = "complete",
         main = "Phylogenetic Correlation Matrix (C)",
         fontsize_row = 7,
         fontsize_col = 7,
         filename = paste0(out_dir, "phylo_correlation_heatmap.png"))

# ---- 7. Heatmap with dendrogram from tree ----
png(paste0(out_dir, "phylo_heatmap_dendro.png"), width = 1200, height = 1000)
heatmap.2(C,
          trace = "none",
          dendrogram = "both",
          Rowv = as.dendrogram(as.hclust(phylo)),
          Colv = as.dendrogram(as.hclust(phylo)),
          col = colorRampPalette(c("white", "blue"))(100),
          main = "Phylogenetic Correlations with Tree Dendrogram",
          key = TRUE,
          cexRow = 0.7,
          cexCol = 0.7)
dev.off()

png(paste0(out_dir, "phylo_heatmap_dendro.png"), width = 1400, height = 1200)
heatmap.2(C,
          trace = "none",
          dendrogram = "both",
          Rowv = as.dendrogram(as.hclust(phylo)),
          Colv = as.dendrogram(as.hclust(phylo)),
          col = colorRampPalette(c("white", "blue"))(100),
          main = "Phylogenetic Correlations with Tree Dendrogram",
          key = TRUE,
          cexRow = 1.0,  # increase row label size
          cexCol = 1.0,  # increase column label size
          srtCol = 45,   # rotate column labels for readability
          adjCol = c(1,1) # adjust rotation alignment
)
dev.off()

C
View(C)
