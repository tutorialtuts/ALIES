# =============================================================================
# ALIES: Arm-Level Immune Exclusion Score
# Pan-Cancer Analysis + IMvigor210 External Validation
#
# Title:   Arm-Level Aneuploidy Encodes Immune Exclusion Across Human Cancers
# Author:  Kishore Kumar S
# Contact: kishorekumars.med@gmail.com
#
# Description:
#   Complete reproducible analysis pipeline for the ALIES study.
#   Covers: data loading, immune deconvolution, arm-immune correlation,
#   LASSO score construction, survival analysis, ICI subgroup analysis,
#   IMvigor210 external validation, sensitivity analysis, and all figures
#   and tables reported in the manuscript.
#
# Input files required (place in working directory):
#   PANCAN_Aneuploidy.xlsx       - Taylor et al. 2018 arm-level CNA calls
#   TCGA_clinical_data.tsv.gz    - TCGA Pan-Cancer clinical data (GDC)
#   TCGA_expression_matrix.tsv.gz- RNA-seq FPKM-UQ (GDC)
#   TCGA_mutation.maf.gz         - Somatic mutations MAF (GDC)
#
# IMvigor210 (Section 9):
#   Requires IMvigor210CoreBiologies R package
#   URL: http://research-pub.gene.com/IMvigor210CoreBiologies/
#
# Outputs:
#   outputs/tables/   - TSV tables (main + supplementary)
#   outputs/figures/  - PDF + PNG figures
#   outputs/rdata/    - Intermediate R objects
#   outputs/session_info.txt
#
# R version: >= 4.3
# Bioconductor version: >= 3.17
#
# Citation:
#   Kishore Kumar S. Arm-Level Aneuploidy Encodes Immune Exclusion Across
#   Human Cancers. [Journal]. [Year]. doi:[doi]
#
# Data availability:
#   TCGA: https://portal.gdc.cancer.gov/
#   Taylor 2018 aneuploidy dataset: doi:10.1016/j.ccell.2018.03.007
#   IMvigor210: doi:10.1038/nature25501
#
# Code licence: MIT
# =============================================================================


# =============================================================================
# SECTION 0: SETUP
# =============================================================================

# Set working directory — change this to your project folder
setwd(".")   # <-- UPDATE THIS PATH

# Create output directories
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/rdata",   recursive = TRUE, showWarnings = FALSE)

# Set global seed for reproducibility
set.seed(2024)


# =============================================================================
# SECTION 1: PACKAGE INSTALLATION AND LOADING
# =============================================================================

cran_pkgs <- c(
  "tidyverse", "data.table", "readxl", "ggplot2", "pheatmap",
  "RColorBrewer", "survival", "survminer", "scales", "ggpubr",
  "patchwork", "broom", "viridis", "glmnet", "pROC", "BiocManager"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

bioc_pkgs <- c("GSVA", "Biobase", "SummarizedExperiment")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("All packages loaded. R version:", R.version$version.string, "\n")


# =============================================================================
# SECTION 2: HELPER FUNCTIONS
# =============================================================================

# Standardise TCGA barcodes
trim_barcode    <- function(x) substr(as.character(x), 1, 15)
trim_barcode_12 <- function(x) substr(as.character(x), 1, 12)

# Retain only primary tumour samples (sample code 01–09)
is_primary_tumor <- function(bc15) {
  substr(bc15, 14, 15) %in% formatC(1:9, width = 2, flag = "0")
}

# Cancer-type-adjusted partial Spearman correlation
# Residualises y on cancer type then computes Spearman rho with x
partial_spearman <- function(data, x_col, y_col, covariate = "cancer_type") {
  df <- data %>%
    select(x = all_of(x_col), y = all_of(y_col), cov = all_of(covariate)) %>%
    filter(!is.na(x), !is.na(y), !is.na(cov))
  if (nrow(df) < 30)
    return(data.frame(rho = NA_real_, p_value = NA_real_, n = nrow(df)))
  y_resid <- tryCatch(
    residuals(lm(y ~ factor(cov), data = df)),
    error = function(e) df$y
  )
  ct <- suppressWarnings(cor.test(as.numeric(df$x), y_resid, method = "spearman"))
  data.frame(rho = as.numeric(ct$estimate), p_value = ct$p.value, n = nrow(df))
}

# Bootstrap 95% CI for cancer-type-adjusted Spearman rho (1,000 replicates)
bootstrap_rho_adjusted <- function(data, arm_col, cd8_col,
                                   covariate = "cancer_type",
                                   n_boot = 1000, seed = 42) {
  set.seed(seed)
  df <- data %>%
    select(arm_v = all_of(arm_col), cd8_v = all_of(cd8_col),
           cov = all_of(covariate)) %>%
    filter(!is.na(arm_v), !is.na(cd8_v), !is.na(cov)) %>%
    mutate(arm_v = suppressWarnings(as.numeric(arm_v)))
  cd8_resid <- tryCatch(
    residuals(lm(cd8_v ~ factor(cov), data = df)),
    error = function(e) df$cd8_v
  )
  obs_rho <- suppressWarnings(cor(df$arm_v, cd8_resid, method = "spearman"))
  boot_rhos <- replicate(n_boot, {
    idx   <- sample(nrow(df), replace = TRUE)
    d_b   <- df[idx, ]
    r_b   <- tryCatch(residuals(lm(cd8_v ~ factor(cov), data = d_b)),
                      error = function(e) d_b$cd8_v)
    suppressWarnings(cor(d_b$arm_v, r_b, method = "spearman"))
  })
  data.frame(
    arm              = arm_col,
    obs_rho          = round(obs_rho, 5),
    CI_2.5           = round(quantile(boot_rhos, 0.025, na.rm = TRUE), 5),
    CI_97.5          = round(quantile(boot_rhos, 0.975, na.rm = TRUE), 5),
    CI_excludes_zero = sign(quantile(boot_rhos, 0.025, na.rm = TRUE)) ==
                       sign(quantile(boot_rhos, 0.975, na.rm = TRUE)),
    n                = nrow(df)
  )
}


# =============================================================================
# SECTION 3: DATA LOADING
# =============================================================================

cat("\n--- Section 3: Loading data ---\n")

# ── 3a. Arm-level CNA (Taylor et al. 2018) ──────────────────────────────────
stopifnot("PANCAN_Aneuploidy.xlsx not found" = file.exists("PANCAN_Aneuploidy.xlsx"))
arm_sheet <- excel_sheets("PANCAN_Aneuploidy.xlsx")[1]
arm_raw   <- as.data.frame(read_excel("PANCAN_Aneuploidy.xlsx", sheet = arm_sheet))

# Auto-detect arm columns (e.g. 1p, 1q, ... 22q)
arm_cols <- names(arm_raw)[grepl("^([0-9]{1,2}|X|Y)[pq]$", names(arm_raw))]
if (length(arm_cols) == 0) {
  # Fallback: content-based detection (values in {-1, 0, 1})
  arm_cols <- names(arm_raw)[sapply(names(arm_raw), function(cn) {
    v <- suppressWarnings(as.numeric(as.character(arm_raw[[cn]])))
    v <- v[!is.na(v)]
    length(v) > 10 && all(v %in% c(-1, 0, 1))
  })]
}
stopifnot("No arm columns detected" = length(arm_cols) > 0)

# Auto-detect sample ID column
arm_id_col <- names(arm_raw)[grepl("sample|barcode|ID", names(arm_raw),
                                   ignore.case = TRUE)][1]
if (is.na(arm_id_col)) arm_id_col <- names(arm_raw)[1]

cat("Arm data:", nrow(arm_raw), "rows |", length(arm_cols), "arm columns\n")

# ── 3b. Clinical data ────────────────────────────────────────────────────────
clinical_raw <- fread("TCGA_clinical_data.tsv.gz", sep = "\t", header = TRUE,
                      fill = TRUE, quote = "")

# ── 3c. Expression data (RNA-seq FPKM-UQ) ───────────────────────────────────
expr_raw <- fread("TCGA_expression_matrix.tsv.gz", sep = "\t", header = TRUE,
                  nThread = 4, fill = TRUE)

# ── 3d. Mutation data (MAF) ──────────────────────────────────────────────────
maf_raw <- fread("TCGA_mutation.maf.gz", sep = "\t", header = TRUE, fill = TRUE)

cat("Data loading complete.\n")


# =============================================================================
# SECTION 4: SAMPLE STANDARDISATION AND FILTERING
# =============================================================================

cat("\n--- Section 4: Standardising barcodes and filtering primary tumours ---\n")

# ── 4a. Arm data ─────────────────────────────────────────────────────────────
arm_data <- arm_raw %>%
  as_tibble() %>%
  rename(sample_id = all_of(arm_id_col)) %>%
  mutate(
    sample_id_15 = trim_barcode(sample_id),
    patient_id   = trim_barcode_12(sample_id)
  ) %>%
  filter(is_primary_tumor(sample_id_15)) %>%
  mutate(across(all_of(arm_cols), ~ suppressWarnings(as.numeric(as.character(.))))) %>%
  distinct(sample_id_15, .keep_all = TRUE) %>%
  select(sample_id_15, patient_id, all_of(arm_cols)) %>%
  mutate(
    aneuploidy_score = rowSums(
      across(all_of(arm_cols), ~ as.integer(!is.na(.) & . != 0)),
      na.rm = TRUE
    )
  )

cat("Arm data (primary tumours):", nrow(arm_data), "samples\n")

# ── 4b. Clinical data ─────────────────────────────────────────────────────────
# Auto-detect key columns
detect_col <- function(data, candidates)
  names(data)[tolower(names(data)) %in% tolower(candidates)][1]

cli_id_col <- detect_col(clinical_raw,
  c("sampleID","sample","bcr_patient_barcode","submitter_id","Tumor_Sample_Barcode"))
if (is.na(cli_id_col))
  cli_id_col <- names(clinical_raw)[which(sapply(clinical_raw, function(v)
    any(grepl("^TCGA-", as.character(v)[1:20]))))[1]]

ct_col <- detect_col(clinical_raw,
  c("type","acronym","cancer_type","project_id","Study","_primary_disease"))

os_time_col  <- detect_col(clinical_raw,
  c("OS.time","os_time","days_to_death","OS_MONTHS","overall_survival"))
os_stat_col  <- detect_col(clinical_raw,
  c("OS","os_status","vital_status","OS_STATUS","Overall_Survival_Status"))

clinical <- clinical_raw %>%
  as_tibble() %>%
  rename(sample_id_raw = all_of(cli_id_col)) %>%
  mutate(
    sample_id_15 = trim_barcode(sample_id_raw),
    patient_id   = trim_barcode_12(sample_id_raw)
  ) %>%
  filter(is_primary_tumor(sample_id_15)) %>%
  distinct(sample_id_15, .keep_all = TRUE)

if (!is.na(ct_col)) {
  clinical <- clinical %>%
    mutate(cancer_type = toupper(trimws(gsub("^TCGA-", "",
                                             as.character(.data[[ct_col]])))))
}

# Add age column if available
age_col <- detect_col(clinical_raw,
  c("age_at_initial_pathologic_diagnosis","age_at_diagnosis","age_at_index","age"))
if (!is.na(age_col))
  clinical$age_at_diagnosis <- suppressWarnings(as.numeric(clinical[[age_col]]))

# Standardise survival columns
if (!is.na(os_time_col))
  clinical$OS_time <- suppressWarnings(as.numeric(clinical[[os_time_col]]))
if (!is.na(os_stat_col))
  clinical$OS_status_raw <- as.character(clinical[[os_stat_col]])

cat("Clinical data (primary tumours):", nrow(clinical), "samples |",
    length(unique(clinical$cancer_type)), "cancer types\n")


# =============================================================================
# SECTION 5: TUMOUR MUTATIONAL BURDEN (TMB)
# =============================================================================

cat("\n--- Section 5: Computing TMB ---\n")

# Non-synonymous mutation types (standard definition)
nonsynon_types <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
  "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
  "Nonstop_Mutation", "Translation_Start_Site"
)

maf_sample_col <- detect_col(maf_raw,
  c("Tumor_Sample_Barcode","sample_id","tumor_barcode"))
maf_type_col   <- detect_col(maf_raw,
  c("Variant_Classification","variant_classification","mutation_type"))

tmb_data <- maf_raw %>%
  as_tibble() %>%
  rename(sample_id_raw = all_of(maf_sample_col),
         variant_class = all_of(maf_type_col)) %>%
  filter(variant_class %in% nonsynon_types) %>%
  mutate(sample_id_15 = trim_barcode(sample_id_raw)) %>%
  group_by(sample_id_15) %>%
  summarise(
    mutation_count = n(),
    TMB            = mutation_count / 38,   # exome size ~38 Mb
    log2_TMB       = log2(TMB + 0.01),
    .groups = "drop"
  )

cat("TMB computed for", nrow(tmb_data), "samples. Median:",
    round(median(tmb_data$TMB), 2), "mutations/Mb\n")


# =============================================================================
# SECTION 6: IMMUNE CELL DECONVOLUTION (ssGSEA)
# =============================================================================

cat("\n--- Section 6: Immune cell deconvolution (ssGSEA) ---\n")

# Published marker gene sets (Bindea et al. 2013; Angelova et al. 2015)
IMMUNE_MARKERS <- list(
  "T_cell_CD8"    = c("CD8A","CD8B","GZMA","GZMB","PRF1","IFNG","TBX21","EOMES","CXCR3"),
  "T_cell_CD4"    = c("CD4","IL2","IL21","CXCR5","BCL6","ICOS","SH2D1A"),
  "Treg"          = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"),
  "B_cell"        = c("CD19","MS4A1","CD79A","CD79B","PAX5","BLK"),
  "NK_cell"       = c("KLRD1","NKG7","GNLY","NCAM1","FCGR3A","KLRB1"),
  "Macrophage_M1" = c("NOS2","TNF","IL1B","IL6","CXCL10","CD68"),
  "Macrophage_M2" = c("MRC1","CD163","MSR1","ARG1","IL10","TGFB1"),
  "Monocyte"      = c("CD14","LYZ","S100A8","S100A9","FCN1","VCAN"),
  "Neutrophil"    = c("ELANE","MPO","PRTN3","AZU1","CEACAM8","FCGR3B"),
  "DC_myeloid"    = c("ITGAX","CLEC9A","FLT3","THBD","CD1C","BATF3")
)

# Prepare expression matrix (genes x samples)
expr_df <- as.data.frame(expr_raw)
gene_col <- names(expr_df)[1]
gene_names <- as.character(expr_df[[gene_col]])
expr_df    <- expr_df[, -1, drop = FALSE]

# Deduplicate genes and samples; keep primary tumours only
gene_names <- gene_names[!duplicated(gene_names)]
expr_df    <- expr_df[!duplicated(as.character(expr_raw[[gene_col]])), ]
col_ids    <- trim_barcode(colnames(expr_df))
keep       <- !duplicated(col_ids) & is_primary_tumor(col_ids) &
              grepl("^TCGA-", col_ids)
expr_df    <- expr_df[, keep, drop = FALSE]
colnames(expr_df) <- col_ids[keep]

expr_mat   <- suppressWarnings(apply(expr_df, 2, as.numeric))
rownames(expr_mat) <- gene_names
expr_mat   <- expr_mat[rowSums(!is.na(expr_mat)) > 0, ]
expr_mat   <- expr_mat[rowSums(expr_mat, na.rm = TRUE) > 0, ]
expr_mat[is.na(expr_mat)] <- 0

# Back-transform if log2-scaled
if (quantile(expr_mat[expr_mat > 0], 0.999, na.rm = TRUE) < 50)
  expr_mat <- pmax(2^expr_mat - 1, 0)

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# Filter marker sets to genes present in matrix
markers_filt <- lapply(IMMUNE_MARKERS, function(g) g[g %in% rownames(expr_mat)])
markers_filt <- markers_filt[sapply(markers_filt, length) >= 2]
cat("Immune cell types with >= 2 markers:", length(markers_filt), "\n")

# Run ssGSEA (GSVA >= 1.40 API)
gsva_result <- tryCatch({
  param <- ssgseaParam(expr_mat, markers_filt, minSize = 2)
  gsva(param, verbose = FALSE)
}, error = function(e) {
  gsva(expr_mat, markers_filt, method = "ssgsea", min.sz = 2, verbose = FALSE)
})

# Normalise scores to [0, 1] per cell type
gsva_norm <- t(apply(gsva_result, 1, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}))

immune_fractions <- as.data.frame(t(gsva_norm)) %>%
  tibble::rownames_to_column("sample_id_15") %>%
  as_tibble() %>%
  mutate(sample_id_15 = trim_barcode(sample_id_15))

immune_cols <- setdiff(names(immune_fractions), "sample_id_15")
cd8_col     <- "T_cell_CD8"   # primary immune readout
cat("Immune fractions computed for", nrow(immune_fractions), "samples\n")

saveRDS(immune_fractions, "outputs/rdata/immune_fractions.rds")


# =============================================================================
# SECTION 7: MASTER DATASET ASSEMBLY
# =============================================================================

cat("\n--- Section 7: Assembling master dataset ---\n")

master <- arm_data %>%
  left_join(
    clinical %>% select(sample_id_15, cancer_type, patient_id,
                        any_of(c("OS_time","OS_status_raw","age_at_diagnosis"))),
    by = "sample_id_15"
  ) %>%
  left_join(
    tmb_data %>% select(sample_id_15, TMB, log2_TMB, mutation_count),
    by = "sample_id_15"
  ) %>%
  left_join(immune_fractions, by = "sample_id_15") %>%
  filter(!is.na(cancer_type), !is.na(aneuploidy_score))

has_tmb <- !all(is.na(master$TMB))
has_age <- "age_at_diagnosis" %in% names(master) &&
           !all(is.na(master$age_at_diagnosis))

cat("Master dataset:", nrow(master), "samples |",
    length(unique(master$cancer_type)), "cancer types\n")

# Retain cancer types with n >= 50 for correlation analyses
valid_cancers <- master %>% count(cancer_type) %>%
  filter(n >= 50) %>% pull(cancer_type)
master_filt <- master %>% filter(cancer_type %in% valid_cancers)
cat("After n >= 50 filter:", nrow(master_filt), "samples |",
    length(valid_cancers), "cancer types\n")

saveRDS(master_filt, "outputs/rdata/master_dataset.rds")
write_tsv(master_filt, "outputs/tables/master_dataset.tsv")


# =============================================================================
# SECTION 8: ARM-IMMUNE CORRELATION ANALYSIS (Aim 1)
# =============================================================================

cat("\n--- Section 8: Arm x immune correlation analysis ---\n")

# Pan-cancer: partial Spearman, cancer-type adjusted, BH-FDR correction
cat("Computing pan-cancer correlations (400 arm x immune pairs)...\n")
pancancer_df <- expand.grid(arm = arm_cols, immune_cell = immune_cols,
                             stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(res = list(partial_spearman(master_filt, arm, immune_cell))) %>%
  unnest(res) %>%
  filter(!is.na(rho)) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p_value, method = "BH"), significant = FDR < 0.05) %>%
  arrange(FDR)

cat("Significant pairs (FDR < 0.05):", sum(pancancer_df$significant), "of",
    nrow(pancancer_df), "\n")

# Per-cancer Spearman correlations
cat("Computing per-cancer correlations...\n")
per_cancer_df <- master_filt %>%
  group_by(cancer_type) %>%
  group_modify(function(d, k) {
    expand.grid(arm = arm_cols, immune_cell = immune_cols,
                stringsAsFactors = FALSE) %>%
      rowwise() %>%
      mutate(
        df  = list(d %>% select(arm_v = all_of(arm),
                                immune_v = all_of(immune_cell)) %>%
                     filter(!is.na(arm_v), !is.na(immune_v))),
        n   = nrow(df[[1]]),
        ct  = if (n >= 15) list(suppressWarnings(
              cor.test(as.numeric(df[[1]]$arm_v), df[[1]]$immune_v,
                       method = "spearman"))) else list(NULL),
        rho     = if (!is.null(ct[[1]])) as.numeric(ct[[1]]$estimate) else NA,
        p_value = if (!is.null(ct[[1]])) ct[[1]]$p.value else NA
      ) %>%
      select(arm, immune_cell, rho, p_value, n) %>%
      filter(!is.na(rho))
  }) %>%
  ungroup() %>%
  group_by(cancer_type) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"), significant = FDR < 0.05) %>%
  ungroup()

write_tsv(pancancer_df,  "outputs/tables/SupplTable_S1_pancancer_correlations.tsv")
write_tsv(per_cancer_df, "outputs/tables/SupplTable_S2_percancer_correlations.tsv")

# Bootstrap 95% CIs for key arms
cat("Bootstrap CIs (cancer-type adjusted; 1,000 replicates)...\n")
key_arms <- c("20q","1q","17q","16p","13q","7p","9p","9q","1p","4q","5q")
key_arms <- key_arms[key_arms %in% arm_cols]
boot_ci  <- map_dfr(key_arms,
                    ~ bootstrap_rho_adjusted(master_filt, .x, cd8_col))
write_tsv(boot_ci, "outputs/tables/SupplTable_S3_bootstrap_CI.tsv")
cat("Bootstrap CI computed for", nrow(boot_ci), "arms\n")

# Permutation test (1,000 replicates) to confirm enrichment over null
cat("Permutation test...\n")
set.seed(2024)
perm_counts <- replicate(1000, {
  perm_data <- master_filt %>%
    mutate(cancer_type_perm = sample(cancer_type))
  n_sig <- 0
  for (arm in arm_cols[1:min(5, length(arm_cols))]) {
    for (ic in immune_cols[1:2]) {
      res <- partial_spearman(perm_data %>%
               rename(cancer_type = cancer_type_perm), arm, ic)
      if (!is.na(res$rho) && res$p_value < 0.05) n_sig <- n_sig + 1
    }
  }
  n_sig
})
cat("Permutation null (mean ± SD significant pairs):",
    round(mean(perm_counts), 1), "±", round(sd(perm_counts), 1), "\n")
cat("Observed significant pairs:", sum(pancancer_df$significant), "\n")
cat("Empirical p:", mean(perm_counts >= sum(pancancer_df$significant)), "\n")

saveRDS(list(pancancer = pancancer_df, per_cancer = per_cancer_df,
             boot_ci = boot_ci),
        "outputs/rdata/correlation_results.rds")


# =============================================================================
# SECTION 8B: CANCER-TYPE-STRATIFIED EFFECT SIZES
# =============================================================================

cat("\n--- Section 8B: Cancer-type-stratified effect sizes ---\n")

# Weighted mean difference (gain vs neutral) within each cancer type,
# pooled by neutral-group n — removes cancer composition bias
stratified_effect <- function(data, arm_col, cd8_col) {
  df <- data.frame(
    arm_v  = suppressWarnings(as.numeric(data[[arm_col]])),
    cd8_v  = data[[cd8_col]],
    cancer = data$cancer_type
  ) %>% filter(!is.na(arm_v), !is.na(cd8_v))
  per_ct <- df %>%
    group_by(cancer, arm_status = case_when(
      arm_v == -1 ~ "Loss", arm_v == 0 ~ "Neutral", arm_v == 1 ~ "Gain")) %>%
    summarise(mean_cd8 = mean(cd8_v), n = n(), .groups = "drop") %>%
    pivot_wider(names_from = arm_status, values_from = c(mean_cd8, n),
                names_sep = "_")
  wtd_gain <- if (all(c("mean_cd8_Neutral","mean_cd8_Gain") %in% names(per_ct))) {
    d <- per_ct %>% filter(!is.na(mean_cd8_Neutral), !is.na(mean_cd8_Gain)) %>%
         mutate(diff = mean_cd8_Gain - mean_cd8_Neutral)
    weighted.mean(d$diff, d$n_Neutral, na.rm = TRUE)
  } else NA_real_
  data.frame(arm = arm_col, wtd_diff_gain = round(wtd_gain, 5))
}

strat_effects <- map_dfr(
  c("20q","1q","17q","16p","13q","7p","9p"),
  ~ stratified_effect(master_filt, .x, cd8_col)
) %>%
  mutate(
    cd8_grand_mean = mean(master_filt[[cd8_col]], na.rm = TRUE),
    pct_diff = round(wtd_diff_gain / cd8_grand_mean * 100, 1)
  )

print(strat_effects)
write_tsv(strat_effects, "outputs/tables/Table_S_effect_sizes.tsv")


# =============================================================================
# SECTION 9: MULTIVARIABLE REGRESSION FOR 1q
# =============================================================================

cat("\n--- Section 9: Multivariable regression ---\n")

if ("1q" %in% arm_cols) {
  mod_data <- master_filt %>%
    select(immune_val  = all_of(cd8_col),
           arm_val     = `1q`,
           aneuploidy_score, cancer_type,
           any_of("log2_TMB")) %>%
    filter(!is.na(immune_val), !is.na(arm_val), !is.na(aneuploidy_score)) %>%
    mutate(arm_val = factor(arm_val, levels = c(-1, 0, 1)))

  pred_str <- if (has_tmb) "arm_val + aneuploidy_score + log2_TMB + cancer_type" else
                           "arm_val + aneuploidy_score + cancer_type"
  fit_1q   <- lm(as.formula(paste("immune_val ~", pred_str)), data = mod_data)
  reg_1q   <- tidy(fit_1q) %>%
    filter(grepl("arm_val", term)) %>%
    mutate(FDR = p.adjust(p.value, method = "BH"),
           arm = "1q")

  cat("1q multivariable regression:\n")
  print(reg_1q)
  write_tsv(reg_1q, "outputs/tables/Table2B_1q_multivariable.tsv")
}


# =============================================================================
# SECTION 10: ALIES SCORE CONSTRUCTION (LASSO-REGULARISED)
# =============================================================================

cat("\n--- Section 10: ALIES construction (LASSO) ---\n")

# Build binary feature matrix: one feature per arm event (gain / loss)
arm_df    <- master_filt %>% select(all_of(arm_cols)) %>%
             mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))
X_gain    <- arm_df %>% mutate(across(everything(), ~ as.integer(!is.na(.) & . == 1)))
X_loss    <- arm_df %>% mutate(across(everything(), ~ as.integer(!is.na(.) & . == -1)))
names(X_gain) <- paste0(names(X_gain), "_gain")
names(X_loss) <- paste0(names(X_loss), "_loss")
X_lasso   <- as.matrix(cbind(X_gain, X_loss))

y_lasso   <- master_filt[[cd8_col]]
valid     <- !is.na(y_lasso)
X_lasso   <- X_lasso[valid, ]
y_resid   <- residuals(lm(y_lasso[valid] ~ factor(master_filt$cancer_type[valid])))

# LASSO (alpha = 1, 10-fold CV, lambda.1se)
set.seed(2024)
lasso_cv  <- cv.glmnet(X_lasso, y_resid, alpha = 1, nfolds = 10,
                       standardize = TRUE, type.measure = "mse")
lasso_coefs <- coef(lasso_cv, s = "lambda.1se")
lasso_df <- data.frame(
  feature = rownames(lasso_coefs),
  coef    = as.numeric(lasso_coefs)
) %>% filter(feature != "(Intercept)", coef != 0) %>%
  arrange(coef)

cat("LASSO selected", nrow(lasso_df), "features at lambda.1se =",
    round(lasso_cv$lambda.1se, 6), "\n")
cat("Immune-cold events (coef < 0):\n")
print(lasso_df %>% filter(coef < 0))

write_tsv(lasso_df, "outputs/tables/SupplTable_S4_LASSO_coefficients.tsv")

# Build ALIES from LASSO immune-cold features
cold_feats   <- lasso_df %>% filter(coef < 0)
alies_raw    <- numeric(nrow(master_filt))
for (i in seq_len(nrow(cold_feats))) {
  feat <- cold_feats$feature[i]
  w    <- abs(cold_feats$coef[i])
  if (grepl("_gain$", feat)) {
    a <- sub("_gain$", "", feat)
    if (!a %in% names(master_filt)) next
    v <- suppressWarnings(as.numeric(master_filt[[a]]))
    alies_raw <- alies_raw + w * as.integer(!is.na(v) & v == 1)
  } else if (grepl("_loss$", feat)) {
    a <- sub("_loss$", "", feat)
    if (!a %in% names(master_filt)) next
    v <- suppressWarnings(as.numeric(master_filt[[a]]))
    alies_raw <- alies_raw + w * as.integer(!is.na(v) & v == -1)
  }
}

master_filt$ALIES <- scales::rescale(alies_raw, to = c(0, 100))
master_filt$ALIES_tertile <- factor(
  case_when(
    master_filt$ALIES <= quantile(master_filt$ALIES, 1/3, na.rm = TRUE) ~ "Low",
    master_filt$ALIES >= quantile(master_filt$ALIES, 2/3, na.rm = TRUE) ~ "High",
    TRUE ~ "Intermediate"
  ),
  levels = c("Low", "Intermediate", "High")
)

# Validate: partial Spearman rho vs CD8
alies_partial <- partial_spearman(master_filt, "ALIES", cd8_col, "cancer_type")
aneu_partial  <- partial_spearman(master_filt, "aneuploidy_score", cd8_col, "cancer_type")
tmb_partial   <- if (has_tmb) partial_spearman(master_filt, "log2_TMB", cd8_col, "cancer_type")
                 else data.frame(rho = NA, p_value = NA, n = NA)

cat(sprintf("ALIES partial rho = %.4f (p = %.2e)\n",
            alies_partial$rho, alies_partial$p_value))
cat(sprintf("Aneuploidy partial rho = %.4f (p = %.2e)\n",
            aneu_partial$rho, aneu_partial$p_value))

biomarker_comp <- bind_rows(
  data.frame(biomarker = "ALIES",           rho = alies_partial$rho,
             p_value = alies_partial$p_value, n = alies_partial$n),
  data.frame(biomarker = "Aneuploidy_score", rho = aneu_partial$rho,
             p_value = aneu_partial$p_value, n = aneu_partial$n),
  if (has_tmb) data.frame(biomarker = "log2_TMB", rho = tmb_partial$rho,
                           p_value = tmb_partial$p_value, n = tmb_partial$n)
)
write_tsv(biomarker_comp, "outputs/tables/SupplTable_S5_biomarker_comparison.tsv")

# Per-cancer ALIES validation
per_cancer_alies <- map_dfr(valid_cancers, function(ct) {
  d <- master_filt %>% filter(cancer_type == ct)
  if (nrow(d) < 30 || var(d$ALIES, na.rm = TRUE) == 0) return(NULL)
  ct_res <- suppressWarnings(cor.test(d$ALIES, d[[cd8_col]], method = "spearman"))
  data.frame(cancer_type = ct, n = nrow(d),
             rho = round(ct_res$estimate, 4), p_value = ct_res$p.value)
}) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"), sig = FDR < 0.05)

cat("Per-cancer ALIES sig neg correlations:",
    sum(per_cancer_alies$sig & per_cancer_alies$rho < 0, na.rm = TRUE), "\n")
write_tsv(per_cancer_alies, "outputs/tables/SupplTable_S6_percancer_ALIES.tsv")

saveRDS(master_filt, "outputs/rdata/master_with_ALIES.rds")


# =============================================================================
# SECTION 11: SENSITIVITY ANALYSIS — EXCLUDE HYPERMUTATED TUMOURS
# =============================================================================

cat("\n--- Section 11: Sensitivity analysis (exclude top 5% TMB) ---\n")

if (has_tmb) {
  tmb_thr      <- quantile(master_filt$TMB, 0.95, na.rm = TRUE)
  master_sens  <- master_filt %>% filter(is.na(TMB) | TMB <= tmb_thr)
  cat("Samples retained after excluding top 5% TMB (>",
      round(tmb_thr, 1), "):", nrow(master_sens), "\n")

  sens_df <- map_dfr(c("20q","1q","17q","16p","13q"), function(a) {
    if (!a %in% names(master_sens)) return(NULL)
    res <- partial_spearman(master_sens, a, cd8_col, "cancer_type")
    data.frame(arm = a, rho_sens = round(res$rho, 5),
               FDR_sens = p.adjust(res$p_value, method = "BH"), n = res$n)
  })

  orig_rhos <- pancancer_df %>%
    filter(immune_cell == cd8_col, arm %in% c("20q","1q","17q","16p","13q")) %>%
    select(arm, rho_full = rho)

  sens_compare <- left_join(sens_df, orig_rhos, by = "arm") %>%
    mutate(pct_attenuation = round((rho_sens - rho_full) / abs(rho_full) * 100, 1))

  cat("Sensitivity results:\n")
  print(sens_compare)
  write_tsv(sens_compare, "outputs/tables/SupplTable_S7_sensitivity.tsv")
}


# =============================================================================
# SECTION 12: SURVIVAL ANALYSIS
# =============================================================================

cat("\n--- Section 12: Survival analysis ---\n")

# Standardise survival data
prep_survival <- function(data) {
  os_t <- names(data)[names(data) %in%
    c("OS_time","OS.time","os_time","days_to_death","overall_survival")][1]
  os_s <- names(data)[names(data) %in%
    c("OS_status_raw","OS","os_status","vital_status","OS_STATUS")][1]
  if (is.na(os_t) || is.na(os_s)) {
    cat("WARNING: survival columns not found\n"); return(NULL)
  }
  data %>%
    mutate(
      OS_time_d   = suppressWarnings(as.numeric(.data[[os_t]])),
      OS_event    = case_when(
        as.character(.data[[os_s]]) %in%
          c("1","Dead","DECEASED","dead","1:DECEASED") ~ 1L,
        as.character(.data[[os_s]]) %in%
          c("0","Alive","LIVING","alive","0:LIVING") ~ 0L,
        TRUE ~ NA_integer_
      ),
      OS_years = OS_time_d / 365.25
    ) %>%
    filter(!is.na(OS_time_d), !is.na(OS_event), OS_time_d > 0,
           !is.na(ALIES_tertile))
}

surv_df <- prep_survival(master_filt)

if (!is.null(surv_df) && nrow(surv_df) > 100) {

  # Multivariable Cox (stratified by cancer type; with age if available)
  cox_vars <- c("ALIES",
                if (has_tmb) "log2_TMB",
                "aneuploidy_score",
                if (has_age) "age_at_diagnosis")
  cox_vars <- cox_vars[cox_vars %in% names(surv_df) &
                       sapply(cox_vars, function(v) !all(is.na(surv_df[[v]])))]

  cox_fit  <- tryCatch(
    coxph(as.formula(paste("Surv(OS_years, OS_event) ~",
                           paste(cox_vars, collapse = " + "),
                           "+ strata(cancer_type)")),
          data = surv_df),
    error = function(e) { cat("Cox failed:", e$message, "\n"); NULL }
  )

  if (!is.null(cox_fit)) {
    cox_tidy <- tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
    cat("Pan-cancer Cox results:\n")
    print(cox_tidy)
    write_tsv(cox_tidy, "outputs/tables/Table3_Cox_pancancer.tsv")
  }

  # C-index: base model vs base + ALIES
  cox_base <- tryCatch(
    coxph(as.formula(paste("Surv(OS_years, OS_event) ~",
                           paste(setdiff(cox_vars, "ALIES"), collapse = " + "),
                           "+ strata(cancer_type)")),
          data = surv_df),
    error = function(e) NULL
  )
  if (!is.null(cox_base) && !is.null(cox_fit)) {
    cat(sprintf("C-index — base: %.3f  |  base + ALIES: %.3f\n",
                summary(cox_base)$concordance[1],
                summary(cox_fit)$concordance[1]))
  }

  # Kaplan–Meier by ALIES tertile
  km_fit  <- survfit(Surv(OS_years, OS_event) ~ ALIES_tertile, data = surv_df)
  km_plot <- ggsurvplot(
    km_fit, data = surv_df,
    risk.table = TRUE, pval = TRUE, conf.int = FALSE,
    palette = c("#2196F3","#FF9800","#F44336"),
    title = "Overall Survival by ALIES Tertile — Pan-Cancer",
    xlab = "Time (years)", ylab = "Survival Probability",
    legend.title = "ALIES",
    ggtheme = theme_classic(base_size = 12),
    risk.table.height = 0.25
  )
  ggsave("outputs/figures/Fig3A_KM_PanCancer_ALIES.pdf",
         print(km_plot), width = 8, height = 7)
  ggsave("outputs/figures/Fig3A_KM_PanCancer_ALIES.png",
         print(km_plot), width = 8, height = 7, dpi = 300)
  cat("Kaplan-Meier figure saved.\n")
}


# =============================================================================
# SECTION 13: ICI SUBGROUP ANALYSIS (TMB-high + 1q gain)
# =============================================================================

cat("\n--- Section 13: ICI subgroup analysis ---\n")

if (has_tmb && "1q" %in% names(master_filt)) {
  ici_data <- master_filt %>%
    filter(!is.na(TMB), !is.na(.data[[cd8_col]])) %>%
    mutate(
      TMB_high    = TMB >= quantile(TMB, 0.80, na.rm = TRUE),
      arm_1q      = suppressWarnings(as.numeric(`1q`)),
      arm_1q_gain = arm_1q == 1
    )

  # Fisher: 1q gain prevalence in TMB-high vs TMB-non-high
  fisher_1q <- fisher.test(table(ici_data$arm_1q_gain, ici_data$TMB_high))
  cat(sprintf("1q gain vs TMB-high: OR = %.3f, p = %.4f\n",
              fisher_1q$estimate, fisher_1q$p.value))

  # Within TMB-high: CD8 fraction by 1q status
  tmb_high_data <- ici_data %>% filter(TMB_high)
  wt_1q <- wilcox.test(
    tmb_high_data[[cd8_col]][tmb_high_data$arm_1q_gain],
    tmb_high_data[[cd8_col]][!tmb_high_data$arm_1q_gain],
    alternative = "less"
  )
  cat(sprintf("1q gain vs non-gain in TMB-high (n=%d): Wilcoxon p = %.4f\n",
              nrow(tmb_high_data), wt_1q$p.value))
  cat(sprintf("  Mean CD8 1q-gain: %.4f | non-gain: %.4f\n",
              mean(tmb_high_data[[cd8_col]][tmb_high_data$arm_1q_gain], na.rm = TRUE),
              mean(tmb_high_data[[cd8_col]][!tmb_high_data$arm_1q_gain], na.rm = TRUE)))

  write_tsv(
    data.frame(
      analysis = c("1q_gain_vs_TMBhigh_Fisher","1q_CD8_in_TMBhigh_Wilcoxon"),
      OR_p     = c(fisher_1q$estimate, NA),
      p_value  = c(fisher_1q$p.value, wt_1q$p.value),
      n        = c(nrow(ici_data), nrow(tmb_high_data))
    ),
    "outputs/tables/Table_ICI_subgroup.tsv"
  )
}


# =============================================================================
# SECTION 14: IMvigor210 EXTERNAL VALIDATION
# =============================================================================

cat("\n--- Section 14: IMvigor210 external validation ---\n")
cat("NOTE: Requires IMvigor210CoreBiologies package.\n")
cat("Install: install.packages('IMvigor210CoreBiologies_1.0.0.tar.gz',\n")
cat("          repos=NULL, type='source')\n")
cat("Download: http://research-pub.gene.com/IMvigor210CoreBiologies/\n\n")

if (requireNamespace("IMvigor210CoreBiologies", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(IMvigor210CoreBiologies)
    library(Biobase)
  })

  # Load data safely (avoids DESeq v1 conflict)
  data("cds", package = "IMvigor210CoreBiologies")
  pkg_data_dir <- system.file("data", package = "IMvigor210CoreBiologies")
  tmp_env <- new.env(parent = emptyenv())
  load(file.path(pkg_data_dir, "fmone.rda"), envir = tmp_env)
  fmone <- tmp_env[[ls(tmp_env)[1]]]; rm(tmp_env)

  # Extract expression (NChannelSet: channel "exprs")
  imv_expr <- tryCatch(assayData(fmone)[["exprs"]],
               error = function(e) exprs(fmone))

  # Clinical and response
  imv_clin <- pData(cds) %>%
    as_tibble(rownames = "sample_id") %>%
    mutate(
      response = case_when(
        grepl("CR|PR", `Best Confirmed Overall Response`) ~ "Responder",
        TRUE ~ "Non-responder"
      ),
      OS_days  = suppressWarnings(as.numeric(`Survival`)),
      OS_event = as.integer(`censOS` == 1)
    )

  cat("IMvigor210:", nrow(imv_clin), "samples |",
      mean(imv_clin$response == "Responder", na.rm = TRUE) * 100, "% ORR\n")

  # Immune deconvolution via ssGSEA
  markers_imv <- lapply(IMMUNE_MARKERS, function(g) g[g %in% rownames(imv_expr)])
  markers_imv <- markers_imv[sapply(markers_imv, length) >= 2]

  gsva_imv <- tryCatch({
    param <- ssgseaParam(imv_expr, markers_imv, minSize = 2)
    gsva(param, verbose = FALSE)
  }, error = function(e) {
    gsva(imv_expr, markers_imv, method = "ssgsea", min.sz = 2, verbose = FALSE)
  })

  imv_immune <- as.data.frame(t(gsva_imv)) %>%
    tibble::rownames_to_column("sample_id") %>%
    mutate(across(-sample_id, ~ {
      x <- as.numeric(.); rng <- range(x, na.rm = TRUE)
      if (diff(rng) == 0) rep(0, length(x)) else (x - rng[1]) / diff(rng)
    }))

  # Compute immune signatures used in Mariathasan et al. 2018
  # TIS (Tumour Inflammation Signature): proxy from marker genes
  tis_genes  <- c("CD8A","IFNG","GZMB","CXCL10","CXCL9","PDCD1","LAG3","TIGIT",
                  "CD27","NKG7","PSMB10","IDO1","CCL5","STAT1","HLA-DQA1")
  tis_genes_present <- tis_genes[tis_genes %in% rownames(imv_expr)]
  tis_scores <- colMeans(imv_expr[tis_genes_present, , drop = FALSE], na.rm = TRUE)

  # TGFb exclusion score: TGFB1 + TGFB2 + TGFB3 + ACTA2 + TAGLN
  tgfb_genes <- c("TGFB1","TGFB2","TGFB3","ACTA2","TAGLN")
  tgfb_present <- tgfb_genes[tgfb_genes %in% rownames(imv_expr)]
  tgfb_scores <- colMeans(imv_expr[tgfb_present, , drop = FALSE], na.rm = TRUE)

  z_score <- function(x) {
    x <- as.numeric(x)
    if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }

  # Non-circular ALIES proxy: TGFb_exclusion - TIS (no CD8 in formula)
  proxy_raw    <- z_score(tgfb_scores) - z_score(tis_scores)
  alies_proxy  <- scales::rescale(proxy_raw, to = c(0, 100))

  master_imv <- imv_clin %>%
    left_join(imv_immune, by = "sample_id") %>%
    mutate(ALIES_proxy = alies_proxy[match(sample_id, names(alies_proxy))],
           TGFb_score  = tgfb_scores[match(sample_id, names(tgfb_scores))],
           TIS_score   = tis_scores[match(sample_id, names(tis_scores))]) %>%
    filter(!is.na(ALIES_proxy))

  cat("IMvigor210 master:", nrow(master_imv), "samples\n")

  # Proxy validation: correlate with immune cell types
  imv_immune_cols <- intersect(names(imv_immune), names(IMMUNE_MARKERS))
  proxy_corr_df <- map_dfr(imv_immune_cols, function(ic) {
    df <- master_imv %>% filter(!is.na(.data[[ic]]))
    ct <- suppressWarnings(cor.test(df$ALIES_proxy, df[[ic]], method = "spearman"))
    data.frame(immune_cell = ic, rho = round(ct$estimate, 4),
               p_value = ct$p.value)
  }) %>% mutate(FDR = p.adjust(p_value, method = "BH"))
  write_tsv(proxy_corr_df, "outputs/tables/SupplTable_S8_IMvigor_proxy_validation.tsv")

  # ICI response: Wilcoxon + ROC
  resp_wilcox <- wilcox.test(
    master_imv$ALIES_proxy[master_imv$response == "Non-responder"],
    master_imv$ALIES_proxy[master_imv$response == "Responder"]
  )
  cat(sprintf("ALIES proxy responders vs non-responders: Wilcoxon p = %.2e\n",
              resp_wilcox$p.value))
  cat(sprintf("  Median — non-responders: %.2f  |  responders: %.2f\n",
              median(master_imv$ALIES_proxy[master_imv$response == "Non-responder"], na.rm = TRUE),
              median(master_imv$ALIES_proxy[master_imv$response == "Responder"], na.rm = TRUE)))

  roc_alies <- roc(master_imv$response, master_imv$ALIES_proxy,
                   levels = c("Responder","Non-responder"), quiet = TRUE)
  cat(sprintf("AUC (ALIES proxy): %.3f (95%% CI %.3f-%.3f)\n",
              auc(roc_alies),
              ci.auc(roc_alies)[1], ci.auc(roc_alies)[3]))

  if (has_tmb) {
    # TMB not directly available in IMvigor210; use FMOne mutation load proxy
    # (available in fmone pData)
    fmone_clin <- pData(fmone) %>%
      as_tibble(rownames = "sample_id") %>%
      select(sample_id, any_of(c("FMOne mutation burden per MB","TMB")))
    if (ncol(fmone_clin) > 1) {
      tmb_col_imv <- names(fmone_clin)[2]
      master_imv <- master_imv %>%
        left_join(fmone_clin, by = "sample_id") %>%
        mutate(TMB_imv = suppressWarnings(as.numeric(.data[[tmb_col_imv]])),
               log2_TMB_imv = log2(TMB_imv + 0.01))
      roc_tmb <- roc(master_imv$response, master_imv$log2_TMB_imv,
                     levels = c("Responder","Non-responder"), quiet = TRUE)
      roc_comb <- roc(master_imv$response,
                      predict(glm(I(response == "Non-responder") ~
                                    ALIES_proxy + log2_TMB_imv,
                                  data = master_imv, family = binomial)),
                      quiet = TRUE)
      cat(sprintf("AUC — TMB: %.3f  |  ALIES: %.3f  |  Combined: %.3f\n",
                  auc(roc_tmb), auc(roc_alies), auc(roc_comb)))
      write_tsv(
        data.frame(model = c("ALIES","TMB","ALIES+TMB"),
                   AUC = round(c(auc(roc_alies), auc(roc_tmb), auc(roc_comb)), 3)),
        "outputs/tables/Table4_IMvigor_AUC.tsv"
      )
    }
  }

  # IMvigor210 Cox regression
  imv_surv <- master_imv %>%
    filter(!is.na(OS_days), !is.na(OS_event), OS_days > 0, !is.na(ALIES_proxy))
  imv_cox_vars <- c("ALIES_proxy",
                    if ("log2_TMB_imv" %in% names(imv_surv)) "log2_TMB_imv",
                    if ("T_cell_CD8" %in% names(imv_surv)) "T_cell_CD8")
  imv_cox_vars <- imv_cox_vars[imv_cox_vars %in% names(imv_surv)]

  if (nrow(imv_surv) > 50 && length(imv_cox_vars) >= 1) {
    imv_cox <- tryCatch(
      coxph(as.formula(paste("Surv(OS_days/365.25, OS_event) ~",
                             paste(imv_cox_vars, collapse = " + "))),
            data = imv_surv),
      error = function(e) NULL
    )
    if (!is.null(imv_cox)) {
      imv_cox_tidy <- tidy(imv_cox, exponentiate = TRUE, conf.int = TRUE) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4)))
      cat("IMvigor210 Cox results:\n")
      print(imv_cox_tidy)
      write_tsv(imv_cox_tidy, "outputs/tables/Table4_IMvigor_Cox.tsv")
    }
  }

  saveRDS(master_imv, "outputs/rdata/master_imv.rds")

} else {
  cat("IMvigor210CoreBiologies not installed. Skipping external validation.\n")
  cat("To install: http://research-pub.gene.com/IMvigor210CoreBiologies/\n")
}


# =============================================================================
# SECTION 15: PUBLICATION FIGURES
# =============================================================================

cat("\n--- Section 15: Generating figures ---\n")

# Chromosomal arm display order
arm_order <- c(paste0(rep(1:22, each = 2), c("p","q")), "Xp","Xq")
arm_order <- arm_order[arm_order %in% arm_cols]

col_fun <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# ── Figure 1: Pan-cancer arm x CD8 heatmap ───────────────────────────────────
cat("Figure 1...\n")
hm_data <- pancancer_df %>%
  filter(immune_cell == cd8_col) %>%
  select(arm, cancer_type = immune_cell, rho) %>%
  { per_cancer_df %>% filter(immune_cell == cd8_col) %>%
    select(arm, cancer_type, rho) %>%
    pivot_wider(names_from = cancer_type, values_from = rho, values_fill = NA) }

if (nrow(hm_data) > 0) {
  hm_mat <- as.matrix(hm_data[, -1])
  rownames(hm_mat) <- hm_data$arm
  arm_ord <- arm_order[arm_order %in% rownames(hm_mat)]
  if (length(arm_ord)) hm_mat <- hm_mat[arm_ord, , drop = FALSE]
  rho_max <- min(0.6, max(abs(hm_mat), na.rm = TRUE))
  for (ext in c("pdf","png")) {
    f <- paste0("outputs/figures/Fig1_Heatmap_Arm_CD8.", ext)
    if (ext == "pdf") pdf(f, width = 14, height = 10) else
      png(f, width = 14, height = 10, units = "in", res = 300)
    pheatmap(hm_mat, color = col_fun,
             breaks = seq(-rho_max, rho_max, length.out = 101),
             cluster_rows = FALSE, cluster_cols = TRUE, na_col = "grey90",
             fontsize_row = 8, fontsize_col = 8, angle_col = 45,
             border_color = NA,
             main = "Cancer-type-adjusted Spearman rho: Arm vs CD8+ T cells")
    dev.off()
  }
}

# ── Figure 2: LASSO coefficient plot ─────────────────────────────────────────
cat("Figure 2...\n")
if (nrow(lasso_df) > 0) {
  p_lasso <- lasso_df %>%
    mutate(feature = reorder(feature, coef),
           direction = ifelse(coef < 0, "Immune-cold", "Immune-hot")) %>%
    ggplot(aes(x = feature, y = coef, fill = direction)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c("Immune-cold" = "#D32F2F",
                                 "Immune-hot"  = "#1565C0")) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(title = "ALIES LASSO coefficients",
         subtitle = paste0("lambda.1se; n = ", nrow(lasso_df), " features selected"),
         x = "Arm event", y = "LASSO coefficient", fill = NULL) +
    theme_classic(base_size = 11)
  ggsave("outputs/figures/Fig2A_LASSO_coefficients.pdf", p_lasso, width = 8, height = 5)
  ggsave("outputs/figures/Fig2A_LASSO_coefficients.png", p_lasso, width = 8, height = 5, dpi = 300)
}

# ALIES vs CD8 scatter
fig2b_data <- master_filt %>%
  select(ALIES, cd8 = all_of(cd8_col), cancer_type) %>%
  filter(!is.na(ALIES), !is.na(cd8))
rho_lbl <- sprintf("rho = %.3f, p = %.2e",
                   alies_partial$rho, alies_partial$p_value)
p_scatter <- ggplot(fig2b_data, aes(x = ALIES, y = cd8)) +
  geom_point(aes(color = cancer_type), alpha = 0.25, size = 0.6) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
  annotate("text", x = max(fig2b_data$ALIES) * 0.6,
           y = max(fig2b_data$cd8) * 0.95, label = rho_lbl, size = 3.5) +
  labs(title = "ALIES vs CD8+ T cell fraction",
       x = "ALIES (0 = immune-hot; 100 = immune-cold)",
       y = "CD8+ T cell fraction") +
  theme_classic(base_size = 11) + theme(legend.position = "none") +
  scale_color_viridis_d(option = "turbo")
ggsave("outputs/figures/Fig2B_ALIES_scatter.pdf", p_scatter, width = 7, height = 5)
ggsave("outputs/figures/Fig2B_ALIES_scatter.png", p_scatter, width = 7, height = 5, dpi = 300)

# Biomarker comparison bar
p_biocomp <- biomarker_comp %>%
  mutate(biomarker = reorder(biomarker, abs(rho))) %>%
  ggplot(aes(x = biomarker, y = rho, fill = rho < 0)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = sprintf("rho=%.3f", rho)),
            vjust = ifelse(biomarker_comp$rho < 0, 1.5, -0.5), size = 3.5) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = "#D32F2F", "FALSE" = "#1565C0"),
                    guide = "none") +
  labs(title = "Biomarker comparison: partial Spearman rho with CD8+",
       subtitle = "Cancer-type adjusted",
       x = NULL, y = "Partial Spearman rho") +
  theme_classic(base_size = 11)
ggsave("outputs/figures/Fig2C_biomarker_comparison.pdf", p_biocomp, width = 6, height = 4)
ggsave("outputs/figures/Fig2C_biomarker_comparison.png", p_biocomp, width = 6, height = 4, dpi = 300)

# ── Bootstrap CI forest plot ──────────────────────────────────────────────────
p_forest <- boot_ci %>%
  arrange(obs_rho) %>%
  mutate(arm = factor(arm, levels = arm)) %>%
  ggplot(aes(x = arm, y = obs_rho, ymin = CI_2.5, ymax = CI_97.5,
             color = obs_rho < 0)) +
  geom_pointrange(linewidth = 0.6, fatten = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  coord_flip() +
  scale_color_manual(values = c("TRUE" = "#D32F2F", "FALSE" = "#1565C0"),
                     guide = "none") +
  labs(title = "Bootstrap 95% CI: arm vs CD8+ T cell fraction",
       subtitle = "1,000 replicates; cancer-type adjusted",
       x = "Chromosome arm", y = "Spearman rho (95% CI)") +
  theme_classic(base_size = 11)
ggsave("outputs/figures/FigS_Bootstrap_CI.pdf", p_forest, width = 7, height = 5)
ggsave("outputs/figures/FigS_Bootstrap_CI.png", p_forest, width = 7, height = 5, dpi = 300)

# ── Per-cancer ALIES bar ──────────────────────────────────────────────────────
p_pc_alies <- per_cancer_alies %>%
  arrange(rho) %>%
  mutate(cancer_type = factor(cancer_type, levels = cancer_type),
         sig = ifelse(sig & rho < 0, "*", "")) %>%
  ggplot(aes(x = cancer_type, y = rho, fill = rho)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  geom_text(aes(label = sig), vjust = 0, size = 5) +
  coord_flip() +
  scale_fill_gradient2(low = "#D32F2F", mid = "white", high = "#1565C0",
                       midpoint = 0, name = "Spearman rho") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(title = "ALIES correlation with CD8+ by cancer type",
       subtitle = "* = FDR < 0.05",
       x = NULL, y = "Spearman rho") +
  theme_classic(base_size = 11)
ggsave("outputs/figures/FigS_PerCancer_ALIES.pdf", p_pc_alies, width = 7, height = 8)
ggsave("outputs/figures/FigS_PerCancer_ALIES.png", p_pc_alies, width = 7, height = 8, dpi = 300)

cat("All figures saved.\n")


# =============================================================================
# SECTION 16: MAIN TABLES
# =============================================================================

cat("\n--- Section 16: Generating tables ---\n")

# Table 1: Cohort characteristics
table1 <- master_filt %>%
  group_by(cancer_type) %>%
  summarise(
    N                     = n(),
    Median_age            = if (has_age) round(median(age_at_diagnosis, na.rm=TRUE), 1) else NA,
    Median_TMB            = round(median(TMB, na.rm = TRUE), 2),
    Median_aneuploidy     = round(median(aneuploidy_score, na.rm = TRUE), 1),
    High_ALIES_pct        = round(mean(ALIES_tertile == "High", na.rm = TRUE) * 100, 1),
    Median_CD8_fraction   = round(median(.data[[cd8_col]], na.rm = TRUE), 4),
    .groups = "drop"
  ) %>% arrange(cancer_type)
write_tsv(table1, "outputs/tables/Table1_cohort_characteristics.tsv")

# Table 2A: Effect sizes for leading immune-depleting arms
write_tsv(strat_effects, "outputs/tables/Table2A_arm_effect_sizes.tsv")

cat("Tables saved.\n")


# =============================================================================
# SECTION 17: SESSION INFO
# =============================================================================

sink("outputs/session_info.txt")
cat("ALIES Analysis — Session Information\n")
cat("Date:", format(Sys.time()), "\n\n")
sessionInfo()
sink()

cat("\n")
cat("=============================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("=============================================================\n")
cat("Samples analysed:  ", nrow(master_filt), "\n")
cat("Cancer types:      ", length(unique(master_filt$cancer_type)), "\n")
cat(sprintf("ALIES partial rho: %.4f (p = %.2e)\n",
            alies_partial$rho, alies_partial$p_value))
cat("Outputs:           outputs/tables/  |  outputs/figures/\n")
cat("=============================================================\n")
