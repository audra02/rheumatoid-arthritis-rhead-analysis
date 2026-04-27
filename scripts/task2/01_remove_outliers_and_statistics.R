# ------------------------------------------------
# 0. DUOMENŲ ĮKĖLIMAS
# ------------------------------------------------

obj <- readRDS("rhead.rds")

# Pagrindinė beta reikšmių matrica
beta <- obj

# CpG ir mėginių anotacijos
probe_anno <- attr(obj, ".annmatrix.rann")
sample_anno <- attr(obj, ".annmatrix.cann")


# ------------------------------------------------
# 1. OUTLIER PAŠALINIMAS
# ------------------------------------------------

outliers <- c(
  "GSM3833612_9704031135_R03C01",
  "GSM3833638_9704031135_R02C01",
  "GSM3833716_9259684070_R01C02",
  "GSM3833615_9704031135_R05C01",
  "GSM3833773_9274651074_R06C01"
)

beta <- beta[, !colnames(beta) %in% outliers]
sample_anno <- sample_anno[!sample_anno$id %in% outliers, ]

# sulygiuojame sample_anno eiliškumą su beta stulpeliais
sample_anno <- sample_anno[match(colnames(beta), sample_anno$id), ]


# ------------------------------------------------
# 2. GRUPIŲ APIBRĖŽIMAS
# ------------------------------------------------

group <- tolower(sample_anno$diagnosis)

ra_idx <- which(group == "ra")
control_idx <- which(group == "control")


# ------------------------------------------------
# 3. HIPOTEZIŲ TESTAVIMAS
# ------------------------------------------------

p_values <- apply(beta, 1, function(row) {
  wilcox.test(row[ra_idx], row[control_idx])$p.value
})

effect_size <- apply(beta, 1, function(row) {
  mean(row[ra_idx]) - mean(row[control_idx])
})


# ------------------------------------------------
# 4. REZULTATŲ LENTELĖ
# ------------------------------------------------

results <- data.frame(
  CpG = rownames(beta),
  p_value = p_values,
  effect_size = effect_size
)


# ------------------------------------------------
# 5. IŠSAUGOJIMAS
# ------------------------------------------------

attr(beta, ".annmatrix.cann") <- sample_anno
attr(beta, ".annmatrix.rann") <- probe_anno

saveRDS(beta, "rhead_filtered.rds")
saveRDS(results, "results.rds")