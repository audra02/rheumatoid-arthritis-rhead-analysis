### =========================
### 0. Aplinkos išvalymas
### =========================
rm(list = ls())
gc()

### =========================
### 1. Duomenų įkėlimas
### =========================
rhead <- readRDS("rhead_filtered.rds")

sample_anno <- attr(rhead, ".annmatrix.cann")

### =========================
### 2. Plokštelė vs diagnozė
### =========================

# Kryžminė lentelė
tab_plate_diag <- table(sample_anno$plate, sample_anno$diagnosis)
tab_plate_diag

# Proporcijos kiekvienoje plokštelėje
prop.table(tab_plate_diag, margin = 1)

### =========================
### 3. Grafiko išsaugojimas
### =========================
png("diagnozes_pasiskirstymas_per_ploksteles.png",
    width = 1400, height = 900, res = 150)

par(mar = c(10, 4, 4, 2))

barplot(t(tab_plate_diag),
        beside = TRUE,
        col = c("steelblue", "tomato"),
        legend.text = colnames(tab_plate_diag),
        args.legend = list(x = "topright"),
        las = 2,
        main = "Diagnozės pasiskirstymas per plokšteles",
        xlab = "Plokštelė",
        ylab = "Mėginių skaičius")

dev.off()
