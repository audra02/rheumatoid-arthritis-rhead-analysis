library(annmatrix)
library(ggplot2)

obj <- readRDS("rhead.rds")
anyNA(obj)

plots_dir <- file.path("..", "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# QC1

# CpG anotacijos
region <- rowanns(obj)$relation_to_island

# Regionų grupės
region2 <- rep(NA_character_, length(region))
region2[region == "Island"] <- "island"
region2[region %in% c("N_Shore", "S_Shore")] <- "shore"
region2[region %in% c("N_Shelf", "S_Shelf")] <- "shelf"
region2[region == "OpenSea"] <- "sea"

keep <- !is.na(region2)

# Vidutinis beta pagal CpG
df <- data.frame(
  beta = rowMeans(obj[keep, ]),
  region = factor(region2[keep], levels = c("island", "shore", "shelf", "sea"))
)

print(tapply(df$beta, df$region, summary))
print(tapply(df$beta, df$region, mean))

d_island <- density(df$beta[df$region == "island"])
d_shore  <- density(df$beta[df$region == "shore"])
d_shelf  <- density(df$beta[df$region == "shelf"])
d_sea    <- density(df$beta[df$region == "sea"])

ymax <- max(d_island$y, d_shore$y, d_shelf$y, d_sea$y) * 1.05

png(file.path(plots_dir, "qc1_cpg_regions.png"), width = 1200, height = 900, res = 150)

par(mar = c(5, 5, 4, 8), xpd = TRUE)

plot(
  d_island,
  col = "green3",
  lwd = 2,
  main = "DNR metilinimo lygis skirtinguose CpG regionuose",
  xlab = "Vidutinis DNR metilinimo lygis",
  ylab = "CpG pozicijų dažnis",
  xlim = c(0, 1),
  ylim = c(0, ymax)
)

lines(d_shore, col = "cornflowerblue", lwd = 2)
lines(d_shelf, col = "orange", lwd = 2)
lines(d_sea, col = "red", lwd = 2)

legend(
  "topright",
  inset = c(-0.18, 0),
  legend = c("island", "shore", "shelf", "sea"),
  col = c("green3", "cornflowerblue", "orange", "red"),
  lwd = 2,
  bty = "n"
)

dev.off()

# QC2

# Mėginių anotacijos
meta <- colanns(obj)

annot <- data.frame(
  sample_id = colnames(obj),
  cell_type = as.character(meta[["celltype"]]),
  diagnosis = as.character(meta[["diagnosis"]]),
  stringsAsFactors = FALSE
)

cat("Duomenų matricos dydis:", dim(obj), "\n")

# Koreliacija tarp mėginių
cor_mat <- cor(obj)

# Porų lentelė
make_pair_df <- function(cor_mat, annot, group_var, label_name) {
  idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
  
  df <- data.frame(
    sample1 = colnames(cor_mat)[idx[, 1]],
    sample2 = colnames(cor_mat)[idx[, 2]],
    correlation = cor_mat[idx],
    stringsAsFactors = FALSE
  )
  
  group_map <- annot[[group_var]]
  names(group_map) <- annot$sample_id
  
  df$group1 <- as.character(group_map[df$sample1])
  df$group2 <- as.character(group_map[df$sample2])
  
  df <- df[!is.na(df$group1) & !is.na(df$group2), , drop = FALSE]
  
  df$pair_type <- ifelse(df$group1 == df$group2, "same", "different")
  df$comparison <- factor(
    ifelse(df$group1 == df$group2,
           paste("Same", label_name),
           paste("Different", label_name)),
    levels = c(paste("Same", label_name), paste("Different", label_name))
  )
  
  df
}

pairs_celltype <- make_pair_df(cor_mat, annot, "cell_type", "cell type")
pairs_diagnosis <- make_pair_df(cor_mat, annot, "diagnosis", "diagnosis")

# Porinių koreliacijų grafikas
plot_pair_boxplot <- function(df, title_text, file_name) {
  p <- ggplot(df, aes(x = comparison, y = correlation)) +
    geom_boxplot(outlier.alpha = 0.25) +
    geom_jitter(width = 0.15, alpha = 0.08, size = 0.7) +
    labs(
      title = title_text,
      x = NULL,
      y = "Porinė mėginių koreliacija"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
  
  ggsave(
    file.path(plots_dir, file_name),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  print(p)
}

plot_pair_boxplot(
  pairs_celltype,
  "Mėginių panašumas pagal ląstelės tipą",
  "pairwise_similarity_celltype.png"
)

plot_pair_boxplot(
  pairs_diagnosis,
  "Mėginių panašumas pagal diagnozę",
  "pairwise_similarity_diagnosis.png"
)

# Koreliacijų santrauka
summarize_pairs <- function(df, label) {
  same_vals <- df$correlation[df$pair_type == "same"]
  diff_vals <- df$correlation[df$pair_type == "different"]
  
  cat("\n", label, "\n", sep = "")
  cat("N vienodų porų:", length(same_vals), "\n")
  cat("N skirtingų porų:", length(diff_vals), "\n")
  cat("Vidutinė koreliacija vienodose porose:", round(mean(same_vals), 4), "\n")
  cat("Vidutinė koreliacija skirtingose porose:", round(mean(diff_vals), 4), "\n")
  cat("Mediana vienodose porose:", round(median(same_vals), 4), "\n")
  cat("Mediana skirtingose porose:", round(median(diff_vals), 4), "\n")
  
  if (length(same_vals) > 1 && length(diff_vals) > 1) {
    wt <- wilcox.test(same_vals, diff_vals)
    cat("Wilcoxon p reikšmė:", format(wt$p.value, scientific = TRUE), "\n")
  }
}

summarize_pairs(pairs_celltype, "ląstelės tipas")
summarize_pairs(pairs_diagnosis, "diagnozė")

# PCA

# Pagrindinių komponenčių analizė
pca <- prcomp(t(obj), center = TRUE)
var_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

pca_df$cell_type <- annot$cell_type[match(pca_df$sample_id, annot$sample_id)]
pca_df$diagnosis <- annot$diagnosis[match(pca_df$sample_id, annot$sample_id)]

# PCA grafikas
plot_pca <- function(df, color_var, title_text, file_name) {
  p <- ggplot(df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      title = title_text,
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(
    file.path(plots_dir, file_name),
    p,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  print(p)
}

plot_pca(
  pca_df,
  "cell_type",
  "PCA pagal ląstelės tipą",
  "pca_celltype.png"
)

plot_pca(
  pca_df,
  "diagnosis",
  "PCA pagal diagnozę",
  "pca_diagnosis.png"
)