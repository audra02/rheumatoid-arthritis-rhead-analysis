# ------------------------------------------------
# REZULTATŲ APRAŠYMAS
# ------------------------------------------------
# Šis skriptas pašalina iš ankstesnių QC žingsnių nustatytus outlier mėginius,
# kiekvienai CpG pozicijai palygina RA ir kontrolės grupes Wilcoxon testu
# ir sukuria rezultatų lentelę su p reikšmėmis, efekto dydžiais,
# koreguotomis p reikšmėmis bei reikšmingumo žyma.
#
# Sukuriami failai:
# - rhead_filtered.rds – filtruotas annmatrix objektas be outlier mėginių
# - results.rds        – visų CpG statistinės analizės rezultatai
# ------------------------------------------------

library(annmatrix)

# ------------------------------------------------
# 1. DUOMENŲ ĮKĖLIMAS
# ------------------------------------------------

obj <- readRDS("rhead.rds")


# ------------------------------------------------
# 2. OUTLIER PAŠALINIMAS
# ------------------------------------------------

outliers <- c(
  "GSM3833612_9704031135_R03C01",
  "GSM3833638_9704031135_R02C01",
  "GSM3833716_9259684070_R01C02",
  "GSM3833615_9704031135_R05C01",
  "GSM3833773_9274651074_R06C01",
  "GSM3833617_9704031129_R05C01",
  "GSM3833483_9406922147_R06C02",
  "GSM3833453_9406922031_R05C02",
  "GSM3833478_9407201070_R06C02",
  "GSM3833640_9704031140_R01C01"
)

obj <- obj[, !colnames(obj) %in% outliers]


# ------------------------------------------------
# 3. GRUPIŲ APIBRĖŽIMAS
# ------------------------------------------------

ra_idx <- which(obj$diagnosis == "ra")
control_idx <- which(obj$diagnosis == "control")


# ------------------------------------------------
# 4. HIPOTEZIŲ TESTAVIMAS
# ------------------------------------------------

p_values <- apply(obj, 1, function(x)
  wilcox.test(x[ra_idx], x[control_idx])$p.value
)

effect_size <- rowMeans(obj[, ra_idx]) - rowMeans(obj[, control_idx])


# ------------------------------------------------
# 5. REZULTATŲ LENTELĖ
# ------------------------------------------------

results <- data.frame(
  CpG = rownames(obj),
  p_value = p_values,
  effect_size = effect_size
)

results$p_adj <- p.adjust(results$p_value, method = "fdr")
results$significant <- results$p_adj < 0.05


# ------------------------------------------------
# 6. IŠSAUGOJIMAS
# ------------------------------------------------

saveRDS(obj, "rhead_filtered.rds")
saveRDS(results, "results.rds")