# =========================================================
# UŽDUOTIS 6. IŠSKIRČIŲ PAIEŠKA SU INTER-ARRAY CORRELATION
# Variantas pagal skaidres:
# outlier = mėginys, kurio vidutinė koreliacija yra
# daugiau nei 3 SD mažesnė už bendrą vidurkį
# =========================================================

# -----------------------------
# 0. DUOMENŲ ĮKĖLIMAS
# -----------------------------

obj <- readRDS("data/rhead.rds")

# Pagrindinė beta reikšmių matrica
beta <- obj

# CpG ir mėginių anotacijos
probe_anno <- attr(obj, ".annmatrix.rann")
sample_anno <- attr(obj, ".annmatrix.cann")

# -----------------------------
# 1. FUNKCIJA:
#    apskaičiuoti vidutinę koreliaciją
# -----------------------------

# Ši funkcija:
# 1) sudaro mėginių koreliacijų matricą
# 2) kiekvienam mėginiui apskaičiuoja vidutinę koreliaciją
#    su visais kitais mėginiais
get_mean_correlation <- function(beta_mat) {
  
  # Koreliacijos tarp mėginių
  cor_mat <- cor(beta_mat, use = "pairwise.complete.obs")
  
  # Savikoreliacijos nevertiname
  diag(cor_mat) <- NA
  
  # Kiekvieno mėginio vidutinė koreliacija su kitais
  mean_cor <- colMeans(cor_mat, na.rm = TRUE)
  
  return(mean_cor)
}

# -----------------------------
# 2. FUNKCIJA:
#    išskirčių žymėjimas pagal paskaitos taisyklę
# -----------------------------

# Pagal skaidres:
# išskirtis = mėginys, kurio vidutinė koreliacija
# yra daugiau nei 3 standartiniais nuokrypiais mažesnė už vidurkį
flag_outliers_sd <- function(x, n_sd = 3) {
  
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  
  lower_bound <- mean_x - n_sd * sd_x
  outlier_flag <- x < lower_bound
  
  return(list(
    flag = outlier_flag,
    mean_x = mean_x,
    sd_x = sd_x,
    lower_bound = lower_bound
  ))
}

# -----------------------------
# 3. FUNKCIJA:
#    viena IAC iteracija
# -----------------------------

run_iac_iteration <- function(beta_mat, anno_df, n_sd = 3) {
  
  # Apskaičiuojame mean correlation
  mean_cor <- get_mean_correlation(beta_mat)
  
  # Pažymime išskirtis pagal mean - 3*SD taisyklę
  out_res <- flag_outliers_sd(mean_cor, n_sd = n_sd)
  
  # Sudarome rezultatų lentelę
  result_df <- data.frame(
    id = colnames(beta_mat),
    mean_cor = mean_cor,
    outlier_iac = out_res$flag,
    overall_mean = as.numeric(out_res$mean_x),
    overall_sd = as.numeric(out_res$sd_x),
    lower_bound = as.numeric(out_res$lower_bound),
    stringsAsFactors = FALSE
  )
  
  # Pridedame mėginių anotacijas
  result_df <- merge(result_df, anno_df, by = "id")
  
  # Surikiuojame nuo mažiausios mean_cor
  result_df <- result_df[order(result_df$mean_cor), ]
  
  return(result_df)
}

# -----------------------------
# 4. 1 ITERACIJA
# -----------------------------

iter1_res <- run_iac_iteration(beta_mat = beta, anno_df = sample_anno, n_sd = 3)

# Išsitraukiame 1 iteracijos išskirtis
iter1_outliers <- iter1_res[iter1_res$outlier_iac,
                            c("id", "celltype", "diagnosis", "age", "smoking",
                              "mean_cor", "lower_bound",
                              "qc_iac", "qc_detection", "qc_purity")]

# -----------------------------
# 5. 2 ITERACIJA
# -----------------------------

# Laikinai pašaliname 1 iteracijoje rastas išskirtis
keep_iter2 <- !(sample_anno$id %in% iter1_outliers$id)

beta_iter2 <- beta[, keep_iter2, drop = FALSE]
sample_anno_iter2 <- sample_anno[keep_iter2, , drop = FALSE]

iter2_res <- run_iac_iteration(beta_mat = beta_iter2, anno_df = sample_anno_iter2, n_sd = 3)

iter2_outliers <- iter2_res[iter2_res$outlier_iac,
                            c("id", "celltype", "diagnosis", "age", "smoking",
                              "mean_cor", "lower_bound",
                              "qc_iac", "qc_detection", "qc_purity")]

# -----------------------------
# 6. 3 ITERACIJA
# -----------------------------

# Laikinai pašaliname 1 ir 2 iteracijose rastas išskirtis
ids_remove_iter3 <- c(iter1_outliers$id, iter2_outliers$id)

keep_iter3 <- !(sample_anno$id %in% ids_remove_iter3)

beta_iter3 <- beta[, keep_iter3, drop = FALSE]
sample_anno_iter3 <- sample_anno[keep_iter3, , drop = FALSE]

iter3_res <- run_iac_iteration(beta_mat = beta_iter3, anno_df = sample_anno_iter3, n_sd = 3)

iter3_outliers <- iter3_res[iter3_res$outlier_iac,
                            c("id", "celltype", "diagnosis", "age", "smoking",
                              "mean_cor", "lower_bound",
                              "qc_iac", "qc_detection", "qc_purity")]

# -----------------------------
# 7. VISŲ ITERACIJŲ SANTRAUKA
# -----------------------------

iter1_outliers$iteration <- 1
iter2_outliers$iteration <- 2

all_outliers <- rbind(iter1_outliers, iter2_outliers, iter3_outliers)

all_outliers

# Kiek unikalių mėginių buvo pažymėta per 3 iteracijas
length(unique(all_outliers$id))

# Kurie mėginiai kartojasi per iteracijas (sanity check - turi but 0)
sort(table(all_outliers$id), decreasing = TRUE)

# -----------------------------
# 8. GALUTINĖ RANKINĖ PERŽIŪRA
# -----------------------------

# Šią lentelę naudosime vertinimui:
# ar išskirtis labiau techninė, ar gali būti biologinės kilmės
all_outliers[order(all_outliers$iteration, all_outliers$mean_cor),
             c("iteration", "id", "celltype", "diagnosis", "age", "smoking",
               "mean_cor", "lower_bound",
               "qc_iac", "qc_detection", "qc_purity")]

# Papildomai: blogiausi mėginiai pagal originalų qc_iac
sample_anno[order(sample_anno$qc_iac),
            c("id", "celltype", "diagnosis", "age", "smoking",
              "qc_iac", "qc_detection", "qc_purity")][1:20, ]

# ============================================
# GRAFIKAS: IAC (mean_cor) ir outlieriai
# ============================================

library(ggplot2)

# -----------------------------
# 1. Paimame 1 iteracijos rezultatus
# (nes jie skaičiuoti visam rinkiniui)
# -----------------------------

plot_df <- iter1_res

# -----------------------------
# 2. Sukuriame outlier žymėjimą
# -----------------------------

# pažymime visus outlier iš visų iteracijų
all_outlier_ids <- unique(all_outliers$id)

plot_df$outlier_all_iter <- plot_df$id %in% all_outlier_ids

# -----------------------------
# 3. Surikiuojame pagal mean_cor
# -----------------------------

plot_df <- plot_df[order(plot_df$mean_cor), ]

# sukuriame indeksą x ašiai
plot_df$index <- 1:nrow(plot_df)

# -----------------------------
# 4. Grafikas
# -----------------------------

grafikas <- ggplot(plot_df, aes(x = index, y = mean_cor)) +
  
  # visi mėginiai
  geom_point(aes(color = outlier_all_iter), size = 2) +
  
  # spalvos
  scale_color_manual(values = c("black", "red")) +
  
  # riba (mean - 3SD)
  geom_hline(yintercept = unique(plot_df$lower_bound),
             linetype = "dashed", color = "blue") +
  
  theme_bw() +
  
  labs(
    title = "Inter-Array Correlation (IAC) analizė",
    x = "Mėginiai (surikiuoti pagal mean_cor)",
    y = "Vidutinė koreliacija (mean_cor)",
    color = "Outlier"
  )

# Išsaugome grafiką kaip PNG failą
ggsave(
  filename = "IAC_plot.png",
  plot = grafikas,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================
# 1. Tik bcell outlier (iš 1 iteracijos)
# ============================================

bcell_outliers <- c(
  "GSM3833716_9259684070_R01C02",
  "GSM3833638_9704031135_R02C01",
  "GSM3833615_9704031135_R05C01"
)

beta_bcell_out <- beta[, bcell_outliers]

# Koreliacijų matrica
cor_bcell_out <- cor(beta_bcell_out, use = "pairwise.complete.obs")

round(cor_bcell_out, 3)

# ============================================
# 2. Paimam kelis normalius bcell mėginius
# ============================================

# visi bcell
all_bcell <- sample_anno$id[sample_anno$celltype == "bcell"]

# išmetam outlier
normal_bcell <- setdiff(all_bcell, bcell_outliers)

# pasiimam pvz 5 atsitiktinius
set.seed(1)
normal_sample <- sample(normal_bcell, 5)

beta_bcell_normal <- beta[, normal_sample]

cor_bcell_normal <- cor(beta_bcell_normal, use = "pairwise.complete.obs")

round(cor_bcell_normal, 3)

bcell_outliers <- c(
  "GSM3833716_9259684070_R01C02",
  "GSM3833638_9704031135_R02C01",
  "GSM3833615_9704031135_R05C01"
)

sample_anno[sample_anno$id %in% bcell_outliers,
            c("id", "celltype", "diagnosis", "age", "smoking",
              "plate", "sentrix_id", "well_row", "well_col",
              "qc_iac", "qc_detection", "qc_purity")]

# ============================================
# 1. Tik tcd4m outlier (iš 2 iteracijos)
# ============================================

tcd4m_ids <- c(
  "GSM3833483_9406922147_R06C02",
  "GSM3833453_9406922031_R05C02",
  "GSM3833617_9704031129_R05C01"
)

beta_tcd4m <- beta[, tcd4m_ids]

round(cor(beta_tcd4m), 3)

# ============================================
# 2. Paimam kelis normalius tcd4m mėginius
# ============================================

# visi tcd4m
all_tcd4m <- sample_anno$id[sample_anno$celltype == "tcd4m"]

# išmetam šituos 3
normal_tcd4m <- setdiff(all_tcd4m, tcd4m_ids)

# paimam kelis normalius
set.seed(1)
normal_sample <- sample(normal_tcd4m, 5)

beta_norm <- beta[, normal_sample]

round(cor(beta_norm), 3)

tcd4m_outliers <- c(
  "GSM3833483_9406922147_R06C02",
  "GSM3833453_9406922031_R05C02",
  "GSM3833617_9704031129_R05C01"
)

sample_anno[sample_anno$id %in% tcd4m_outliers,
            c("id", "celltype", "diagnosis", "age", "smoking",
              "plate", "sentrix_id", "well_row", "well_col",
              "qc_iac", "qc_detection", "qc_purity")]


# ============================================
# 1. Nurodome mėginius, kuriuos šaliname
# ============================================

remove_ids <- c(
  "GSM3833612_9704031135_R03C01",
  "GSM3833638_9704031135_R02C01",
  "GSM3833716_9259684070_R01C02",
  "GSM3833615_9704031135_R05C01"
)

# ============================================
# 2. Patikriname, ar visi ID egzistuoja
# ============================================

all(remove_ids %in% colnames(beta))

# ============================================
# 3. Pašaliname mėginius iš beta matricos
# ============================================

beta_filtered <- beta[, !colnames(beta) %in% remove_ids]

# ============================================
# 4. Pašaliname iš sample anotacijų
# ============================================

sample_anno_filtered <- sample_anno[!sample_anno$id %in% remove_ids, ]

# ============================================
# 5. Patikriname, ar viskas sutampa
# ============================================

stopifnot(all(colnames(beta_filtered) == sample_anno_filtered$id))

# ============================================
# 6. Sukuriame naują annmatrix objektą
# ============================================

# nukopijuojame seną objektą
obj_filtered <- beta_filtered

# grąžiname anotacijas
attr(obj_filtered, ".annmatrix.rann") <- probe_anno
attr(obj_filtered, ".annmatrix.cann") <- sample_anno_filtered

# ============================================
# 7. Patikriname matmenis
# ============================================

dim(obj_filtered)        # turėtų būti 485512 x 367
dim(sample_anno_filtered)

# ============================================
# 8. Išsaugome naują .rds failą
# ============================================

saveRDS(obj_filtered, file = "rhead_filtered.rds")


# =========================================================
# UŽDUOTIS 7. Koreliacija
# =========================================================


# 0. Clean start
rm(list = ls())
cat("\014")
gc()

# 1. Load data
rhead <- readRDS("rhead_filtered.rds")

beta <- rhead
sample_anno <- attr(rhead, ".annmatrix.cann")
probe_anno  <- attr(rhead, ".annmatrix.rann")


# 2. Clustering (PAGRINDINIS)

# Koreliacija tarp mėginių
cor_mat <- cor(beta)

# Koreliacijos atstumas
dist_mat <- as.dist(1 - cor_mat)

# Hierarchinis klasterizavimas
hc <- hclust(dist_mat, method = "average")

# Dendrograma
png("dendrogram_full.png", width = 1200, height = 800)

plot(hc,
     labels = FALSE,
     main = "Hierarchinis klasterizavimas (1 - cor)",
     xlab = "Mėginiai",
     ylab = "Atstumas")

rect.hclust(hc, k = 5, border = "red")

dev.off()


# 3. Cluster priskyrimas

clusters <- cutree(hc, k = 5)

# Priskirti cluster prie sample_anno
sample_anno$cluster <- clusters[sample_anno$id]

# Patikrinti biologinę reikšmę
table(sample_anno$cluster, sample_anno$celltype)
table(sample_anno$cluster, sample_anno$diagnosis)

### =========================
### 4. SENSITIVITY ANALYSIS
### (be vieno bcell mėginio)
### =========================

outlier_id <- "GSM3833773_9274651074_R06C01"

# Pašaliname tik šiai analizei
beta_no1 <- beta[, colnames(beta) != outlier_id]

# Clustering iš naujo
cor_mat2 <- cor(beta_no1)
dist_mat2 <- as.dist(1 - cor_mat2)
hc2 <- hclust(dist_mat2, method = "average")

# Dendrograma be outlier
png("dendrogram_no_outlier.png", width = 1200, height = 800)

plot(hc2,
     labels = FALSE,
     main = "Clustering be vieno bcell mėginio",
     xlab = "Mėginiai",
     ylab = "Atstumas")

rect.hclust(hc2, k = 4, border = "blue")

dev.off()

# 5. Cluster priskyrimas (be outlier)

clusters2 <- cutree(hc2, k = 4)

# atitinkamai sumažinti sample_anno
sample_anno_no1 <- sample_anno[sample_anno$id != outlier_id, ]

sample_anno_no1$cluster <- clusters2[sample_anno_no1$id]

# Patikrinti struktūrą
table(sample_anno_no1$cluster, sample_anno_no1$celltype)
table(sample_anno_no1$cluster, sample_anno_no1$diagnosis)
