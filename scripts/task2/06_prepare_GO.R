# ------------------------------------------------
# 06_prepare_GO.R
# ------------------------------------------------
# Paruošiami BED failai GREAT / GO analizei.
# Background: visos tirtos CpG pozicijos.
# Foreground: patikimos CpG pozicijos atskirai pagal metilinimo kryptį.
# ------------------------------------------------

library(annmatrix)

# ------------------------------------------------
# 1. DUOMENŲ ĮKĖLIMAS
# ------------------------------------------------

obj <- readRDS("rhead_filtered.rds")
results <- readRDS("results.rds")

stopifnot(all(c("CpG", "p_adj", "effect_size") %in% colnames(results)))

# ------------------------------------------------
# 2. ANOTACIJŲ PARUOŠIMAS
# ------------------------------------------------

ann <- rowanns(obj)

if (!"id" %in% colnames(ann)) {
  ann$id <- rownames(obj)
}

chr_col <- c("chr", "CHR", "chromosome", "Chromosome")[
  c("chr", "CHR", "chromosome", "Chromosome") %in% colnames(ann)
][1]

pos_col <- c("pos", "POS", "MAPINFO", "mapinfo", "position", "Position")[
  c("pos", "POS", "MAPINFO", "mapinfo", "position", "Position") %in% colnames(ann)
][1]

if (is.na(chr_col) || is.na(pos_col)) {
  stop("Anotacijose nerasti chromosomos arba pozicijos stulpeliai.")
}

ann <- ann[match(results$CpG, ann$id), ]
stopifnot(all(ann$id == results$CpG))

bed_dat <- data.frame(
  chr = as.character(ann[[chr_col]]),
  pos = as.numeric(ann[[pos_col]]),
  CpG = results$CpG,
  p_adj = results$p_adj,
  effect_size = results$effect_size,
  stringsAsFactors = FALSE
)

bed_dat <- bed_dat[!is.na(bed_dat$chr) & !is.na(bed_dat$pos), ]

bed_dat$chr <- ifelse(
  grepl("^chr", bed_dat$chr),
  bed_dat$chr,
  paste0("chr", bed_dat$chr)
)

# BED formatas yra 0-based: start = pos - 1, end = pos.
bed_dat$start <- as.integer(bed_dat$pos - 1)
bed_dat$end <- as.integer(bed_dat$pos)

# ------------------------------------------------
# 3. FUNKCIJA BED FAILUI IŠSAUGOTI
# ------------------------------------------------
# Naudojamas format(..., scientific = FALSE), kad koordinatės nebūtų išsaugotos scientific notation formatu

write_bed <- function(dat, file) {
  out <- data.frame(
    chr = dat$chr,
    start = format(dat$start, scientific = FALSE, trim = TRUE),
    end = format(dat$end, scientific = FALSE, trim = TRUE),
    CpG = dat$CpG,
    stringsAsFactors = FALSE
  )
  
  write.table(
    out,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# ------------------------------------------------
# 4. BACKGROUND IR FOREGROUND BED FAILAI
# ------------------------------------------------

fdr_cutoff <- 0.05

background <- bed_dat
ra_hyper <- bed_dat[bed_dat$p_adj < fdr_cutoff & bed_dat$effect_size > 0, ]
control_hyper <- bed_dat[bed_dat$p_adj < fdr_cutoff & bed_dat$effect_size < 0, ]

write_bed(
  background,
  "../../plots/task2/background_all_tested_cytosines.bed"
)

write_bed(
  ra_hyper,
  "../../plots/task2/foreground_RA_hypermethylated_cytosines.bed"
)

write_bed(
  control_hyper,
  "../../plots/task2/foreground_control_hypermethylated_cytosines.bed"
)

cat("Background CpG:", nrow(background), "\n")
cat("RA labiau metilintos CpG:", nrow(ra_hyper), "\n")
cat("Kontrolė labiau metilintos CpG:", nrow(control_hyper), "\n")
