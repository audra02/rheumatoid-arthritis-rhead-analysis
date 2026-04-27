### Remove outliers

# Load data
rhead <- readRDS("rhead.rds")

# Outlier sample IDs
outliers <- c(
  "GSM3833612_9704031135_R03C01",
  "GSM3833638_9704031135_R02C01",
  "GSM3833716_9259684070_R01C02",
  "GSM3833615_9704031135_R05C01",
  "GSM3833773_9274651074_R06C01"  # bcell outlier from clustering
)

# Remove outliers
rhead_filtered <- rhead[, !colnames(rhead) %in% outliers]

# Save result
saveRDS(rhead_filtered, "rhead_filtered.rds")

#################################################

#1. Palyginkite grupes naudojant statistinį testą

#################################################
