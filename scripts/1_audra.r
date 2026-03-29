library(annmatrix)

# Užkraunami duomenys
obj <- readRDS("rhead.rds")

# Sukuriamas kelias į plots aplanką
plots_dir <- file.path("..", "plots")

# Jei plots aplanko nėra, jis sukuriamas
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# =========================
# KOKYBĖS KONTROLĖ #1
# Ar duomenyse matoma tikėtina DNR metilinimo tendencija
# skirtinguose CpG regionuose?
# =========================

# Išsitraukiamos CpG anotacijos
ann <- obj@''

# Sukuriamos regionų grupės
region <- ann$relation_to_island
region2 <- rep(NA_character_, length(region))

region2[region == "Island"] <- "island"
region2[region %in% c("N_Shore", "S_Shore")] <- "shore"
region2[region %in% c("N_Shelf", "S_Shelf")] <- "shelf"
region2[region == "OpenSea"] <- "sea"

# Paliekamos tik reikalingos CpG pozicijos
keep <- !is.na(region2)

# Apskaičiuojamas vidutinis beta kiekvienai CpG pozicijai
mean_beta <- rowMeans(obj[keep, ], na.rm = TRUE)

# Sukuriamas duomenų rėmelis grafikui
df <- data.frame(
  beta = mean_beta,
  region = factor(region2[keep], levels = c("island", "shore", "shelf", "sea"))
)

# Santrauka
print(tapply(df$beta, df$region, summary))
print(tapply(df$beta, df$region, mean))

# Tankio kreivės
d_island <- density(df$beta[df$region == "island"], na.rm = TRUE)
d_shore  <- density(df$beta[df$region == "shore"],  na.rm = TRUE)
d_shelf  <- density(df$beta[df$region == "shelf"],  na.rm = TRUE)
d_sea    <- density(df$beta[df$region == "sea"],    na.rm = TRUE)

# Automatinė y ašies riba
ymax <- max(d_island$y, d_shore$y, d_shelf$y, d_sea$y) * 1.05

# Išsaugo į PNG
png(file.path(plots_dir, "qc1_cpg_regions.png"), width = 1200, height = 900, res = 150)

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

# Išsitraukiamos mėginių anotacijos
cann <- attr(obj, ".annmatrix.cann")

# Tikslus ląstelės tipo stulpelio pavadinimas
cell_col <- "celltype"

# Parodomas ląstelės tipo pasiskirstymas
cat("Ląstelės tipo pasiskirstymas:\n")
print(table(cann[[cell_col]], useNA = "ifany"))

# Kad skaičiavimai būtų greitesni,
# pasirenkamos labiausiai kintančios CpG pozicijos
vars <- apply(obj, 1, var, na.rm = TRUE)
top_n <- 10000
top_idx <- order(vars, decreasing = TRUE)[1:top_n]

# Sukuriama mažesnė matrica koreliacijos skaičiavimui
mat <- obj[top_idx, ]

# Apskaičiuojama mėginių koreliacijų matrica
cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")

# Mėginiai surikiuojami pagal ląstelės tipą
ord_cell <- order(as.character(cann[[cell_col]]), colnames(obj))

# Pertvarkoma koreliacijų matrica
cor_ord <- cor_mat[ord_cell, ord_cell]

# Išsaugoma į PNG
png(file.path(plots_dir, "qc2_all_samples_heatmap_by_celltype.png"), width = 1400, height = 1200, res = 150)

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