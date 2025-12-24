###############################################################
# GXE Stability Analysis Pipeline
# Mixed Model (BLUP) + WAASB + GGE
# Author: Dr. Kavyashree N M
###############################################################

## ---- Load required packages ----
## NOTE: Please install required packages before running this script
## tidyverse, readxl, metan, openxlsx

library(readxl)
library(metan)
library(tidyverse)
library(openxlsx)

##------Helper function to create workbook---
write_with_title <- function(wb, sheet, title, data) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, title, startRow = 1, startCol = 1)
  writeData(wb, sheet, data, startRow = 3, startCol = 1)
}

## ---- Define project directories ----
data_dir    <- "STAB_DATA"
results_dir <- "STAB_RESULTS"
figures_dir <- file.path(results_dir, "FIGURES")
qc_dir <- file.path(results_dir, "QC")

## Create directories if they do not exist
dir.create(data_dir, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)
dir.create(qc_dir, showWarnings = FALSE)

## Output workbook path
workbook <- file.path(results_dir, "MET_Summary.xlsx")

## ---- Load dataset ----
data_file <- file.path(data_dir, "DATA.xlsx")
if (!file.exists(data_file)) {
  stop("❌ DATA.xlsx not found in STAB_DATA folder.")
}
dataset <- read_excel(data_file)

## ---- Validate required columns ----
required_cols <- c("GENO", "ENV", "REP", "HT", "YLD")

missing_cols <- setdiff(required_cols, colnames(dataset))
if (length(missing_cols) > 0) {
  stop(paste("❌ Missing required columns:", 
             paste(missing_cols, collapse = ", ")))
}

## ---- Convert design variables to factors ----
dataset<- dataset %>%
  mutate(
    GENO = factor(GENO),
    ENV  = factor(ENV),
    REP  = factor(REP),
    HT   = as.numeric(HT),
    YLD  = as.numeric(YLD)
  )
## ---- Dataset structure ----
str(dataset)
summary(dataset)
## ---- Handle missing and zero values ----
dataset_clean <- dataset %>%
  remove_rows_na() %>%        # removes rows with NA
  replace_rows_zero()         # replaces 0 with NA
## ---- Outlier detection ----
png(file.path(qc_dir, "Outliers_HT.png"), width = 1600, height = 1200, res = 150)
find_outliers(dataset_clean, var = HT, plots = TRUE)
dev.off()

png(file.path(qc_dir, "Outliers_YLD.png"), width = 1600, height = 1200, res = 150)
find_outliers(dataset_clean, var = YLD, plots = TRUE)
dev.off()
#Now onwards Use cleaned data for analysis 
data_use <- dataset_clean #If data is already clean, no need to rename

## ---- Overall descriptive statistics ----
desc_all <- metan::desc_stat(
  data_use,
  stats = c("mean", "sd", "cv", "min", "max")
)

mean_geno <- metan::means_by(data_use, GENO) #Genotype-wise means
mean_env <- metan::means_by(data_use, ENV) #Environment-wise means
mean_ge <- data_use %>%
  group_by(ENV, GENO) %>%
  metan::desc_stat(HT, YLD, stats = "mean") #GxE Means

## ---- Create descriptive workbook ----
wb_desc <- createWorkbook()

write_with_title(wb_desc, "Overall_Stats",
                 "Descriptive Statistics (Cleaned Data)", desc_all)

write_with_title(wb_desc, "Genotype_Means",
                 "Mean Performance of Genotypes", mean_geno)

write_with_title(wb_desc, "Environment_Means",
                 "Mean Performance of Environments", mean_env)

write_with_title(wb_desc, "GE_Means",
                 "Genotype × Environment Means", mean_ge)

saveWorkbook(
  wb_desc,
  file.path(results_dir, "Descriptive_Statistics.xlsx"),
  overwrite = TRUE
)

##---Mixed Model (BLUP) using gamem_met()---
#BLUP model for HT
model_HT <- gamem_met(
  .data = data_use,
  env  = ENV,
  gen  = GENO,
  rep  = REP,
  resp = HT
)

plot(model_HT, type = "re") #Random effects distribution
vcomp_HT <- get_model_data(model_HT, "vcomp") #Variance components (HT)
print(vcomp_HT)

blupg_HT <- get_model_data(model_HT, "blupg") #Genotype BLUPs (HT)
blupge_HT <- get_model_data(model_HT, "blupge") #G×E BLUPs (HT)
bluege_HT <- get_model_data(model_HT, "bluege") #Predicted values (HT) (BLUPg + BLUPge)

#BLUP model for YLD
model_YLD <- gamem_met(
  .data = data_use,
  env  = ENV,
  gen  = GENO,
  rep  = REP,
  resp = YLD
)

plot(model_YLD, type = "re")
vcomp_YLD <- get_model_data(model_YLD, "vcomp")
blupg_YLD  <- get_model_data(model_YLD, "blupg")
blupge_YLD <- get_model_data(model_YLD, "blupge")
bluege_YLD <- get_model_data(model_YLD, "bluege")

##---WAASB (Weighted Average of Absolute Scores)---
#WAASB model (HT + YLD)
model_waasb <- waasb(
  .data = data_use,
  env  = ENV,
  gen  = GENO,
  rep  = REP,
  resp = c(HT, YLD)
)

class(model_waasb)
names(model_waasb)

#WAASB indices (Lower WAASB = more stable genotype)
waasb_indices <- get_model_data(model_waasb, "WAASB") #for all traits

ipca_weights <- get_model_data(model_waasb, "PctWAASB") #IPCA contribution (weights), This replaces old AMMI F-tests
lrt_results <- get_model_data(model_waasb, "lrt") #Likelihood Ratio Test (G×E significance)
genetic_params <- get_model_data(model_waasb, "genpar") #Genetic parameters

## ---- Mean and stability statistics ----
stab_stats <- ge_stats(data_use, ENV, GENO, REP, YLD) #gives all stab parameters for selected trait YIELD

yld_stats <- stab_stats$YLD #extract YLD table
colnames(yld_stats)

## ---- Mean vs Stability table ----
mean_stability_table <- yld_stats %>%
  select(GEN, Y, WAASB) %>%
  rename(Mean_YLD = Y) %>%
  mutate(
    Rank_YLD   = rank(-Mean_YLD),
    Rank_WAASB = rank(WAASB),
    Total_Score = Rank_YLD + Rank_WAASB
  ) %>%
  arrange(Total_Score)

#Identify top candidates
top_genotypes <- mean_stability_table %>%
  slice(1:10)

print(top_genotypes)

## ---- Save selection table ----
write.csv(
  mean_stability_table,
  file.path(results_dir, "Mean_vs_WAASB_YLD.csv"),
  row.names = FALSE
)

##---GGE Biplots---
## ---- GGE model for Yield ----
model_gge <- gge(
  .data = data_use,
  env  = ENV,
  gen  = GENO,
  resp = YLD
)

## ---- Mean vs Stability plot ----
plot(
  model_gge,
  type = 2,
  col.gen = "black",
  col.env = "blue",
  repel = TRUE,
  title = "GGE Biplot: Mean vs Stability (YLD)"
)

## ---- Who-Won-Where (Polygon) ----
plot(
  model_gge,
  type = 3,          # polygon view
  scaling = 0,
  col.gen = "blue",
  col.env = "darkgreen",
  repel = TRUE,
  title = "GGE Biplot: Who-Won-Where (YLD)"
)

## ---- Save GGE plots ----
ggsave(
  filename = file.path(figures_dir, "GGE_Mean_Stability_YLD.png"),
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(figures_dir, "GGE_Who_Won_Where_YLD.png"),
  width = 8,
  height = 6,
  dpi = 300
)




