library(annmatrix)

# Užkraunami duomenys
obj <- readRDS("rhead.rds")

# =========================
# KOKYBĖS KONTROLĖ #1
# Ar duomenyse matoma tikėtina DNR metilinimo tendencija
# skirtinguose CpG regionuose?
# =========================

ann <- obj@''

region <- ann$relation_to_island
region2 <- rep(NA_character_, length(region))

region2[region == "Island"] <- "island"
region2[region %in% c("N_Shore", "S_Shore")] <- "shore"
region2[region %in% c("N_Shelf", "S_Shelf")] <- "shelf"
region2[region == "OpenSea"] <- "sea"

keep <- !is.na(region2)

mean_beta <- rowMeans(obj[keep, ], na.rm = TRUE)

df <- data.frame(
  beta = mean_beta,
  region = factor(region2[keep], levels = c("island", "shore", "shelf", "sea"))
)

print(tapply(df$beta, df$region, summary))
print(tapply(df$beta, df$region, mean))

d_island <- density(df$beta[df$region == "island"], na.rm = TRUE)
d_shore  <- density(df$beta[df$region == "shore"],  na.rm = TRUE)
d_shelf  <- density(df$beta[df$region == "shelf"],  na.rm = TRUE)
d_sea    <- density(df$beta[df$region == "sea"],    na.rm = TRUE)

ymax <- max(d_island$y, d_shore$y, d_shelf$y, d_sea$y) * 1.05

png("qc1_cpg_regions.png", width = 1200, height = 900, res = 150)

par(
  mar = c(5, 5, 4, 8),
  xpd = TRUE
)

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
lines(d_sea,   col = "red", lwd = 2)

legend(
  "topright",
  inset = c(-0.18, 0),
  legend = c("island", "shore", "shelf", "sea"),
  col = c("green3", "cornflowerblue", "orange", "red"),
  lwd = 2,
  bty = "n"
)

dev.off()

# Papildomai pavaizduoja ir R lange
par(
  mar = c(5, 5, 4, 8),
  xpd = TRUE
)

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
lines(d_sea,   col = "red", lwd = 2)

legend(
  "topright",
  inset = c(-0.18, 0),
  legend = c("island", "shore", "shelf", "sea"),
  col = c("green3", "cornflowerblue", "orange", "red"),
  lwd = 2,
  bty = "n"
)

# =========================
# KOKYBĖS KONTROLĖ #2
# Ar mėginiai, kurie turėtų būti panašūs,
# iš tiesų yra panašesni tarpusavyje?
# Vertinama naudojant koreliacijų matricą
# pagal ląstelės tipą.
# =========================

cann <- attr(obj, ".annmatrix.cann")

cell_col <- "celltype"

# Parodomas ląstelės tipo pasiskirstymas
cat("Ląstelės tipo pasiskirstymas:\n")
print(table(cann[[cell_col]], useNA = "ifany"))

# Kad skaičiavimai būtų greitesni,
# pasirenkamos labiausiai kintančios CpG pozicijos
vars <- apply(obj, 1, var, na.rm = TRUE)
top_n <- 10000
top_idx <- order(vars, decreasing = TRUE)[1:top_n]

mat <- obj[top_idx, ]

cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")

# Mėginiai surikiuojami pagal ląstelės tipą
ord_cell <- order(as.character(cann[[cell_col]]), colnames(obj))

cor_ord <- cor_mat[ord_cell, ord_cell]

png("qc2_all_samples_heatmap_by_celltype.png", width = 1400, height = 1200, res = 150)

par(mar = c(6, 6, 4, 2))
image(
  1:ncol(cor_ord), 1:ncol(cor_ord),
  t(cor_ord[nrow(cor_ord):1, ]),
  col = colorRampPalette(c("white", "steelblue", "navy"))(100),
  zlim = c(-1, 1),
  axes = FALSE,
  xlab = "Mėginiai",
  ylab = "Mėginiai",
  main = "Mėginių koreliacijų matrica pagal ląstelės tipą"
)
box()

dev.off()

# Papildomai pavaizduoja ir R lange
par(mar = c(6, 6, 4, 2))
image(
  1:ncol(cor_ord), 1:ncol(cor_ord),
  t(cor_ord[nrow(cor_ord):1, ]),
  col = colorRampPalette(c("white", "steelblue", "navy"))(100),
  zlim = c(-1, 1),
  axes = FALSE,
  xlab = "Mėginiai",
  ylab = "Mėginiai",
  main = "Mėginių koreliacijų matrica pagal ląstelės tipą"
)
box()