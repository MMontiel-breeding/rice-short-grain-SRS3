################################################################################
# GRAIN MORPHOLOGY ANALYSIS
#
# Author: Montiel et al.
#
# Script 01 | Grain Morphology Statistical Analysis
# Manuscript:
# Genetic mapping of a putative spontaneous short-grain mutation identifies a
# predicted loss-of-function allele of SRS3 in rice
#
################################################################################

################################################################################
# IMPORT
################################################################################

rm(list = ls())

# Run from the repository root. The Excel file may also be beside the script.
# Required packages: dplyr, tidyr, tibble, ggplot2, writexl and openxlsx.
packages <- c("dplyr", "tidyr", "tibble", "ggplot2", "writexl", "openxlsx")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if(length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(writexl)
  library(openxlsx)
})

input_file <- file.path("S1_Data_Grain_Morphology_Measurements.xlsx") # <- Input file # Input: raw individual-grain measurements
if(!file.exists(input_file) && file.exists(basename(input_file))) input_file <- basename(input_file)

output_dir <- "supplementary_information"
exclude_outliers <- TRUE                 # FALSE retains valid IQR-flagged records

if(!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- openxlsx::read.xlsx(input_file, sheet = 1)

traits <- c("Length", "Width", "L_W_Ratio")
required <- c("ImageID", "Classification", "Genotype", traits)

if(!all(required %in% names(dat))) stop("Missing required columns.")

if(!nrow(dat) || !all(vapply(dat[traits], is.numeric, logical(1)))) {
  stop("Measurements must be numeric and contain records.")
}

# Standardize phenotype label used in the Excel workbook.
dat$Genotype <- ifelse(dat$Genotype == "Wild", "WT", dat$Genotype)

if(anyNA(dat[c("Classification", "Genotype")]) ||
   any(!dat$Classification %in% c("Rough", "Brown", "Milled")) ||
   any(!dat$Genotype %in% c("WT", "Mutant"))) {
  stop("Check phenotype and processing labels.")
}

dat$Record_ID <- seq_len(nrow(dat))         # Links each record to its original input row

cat("\nRows:", nrow(dat), "\n")

################################################################################
# HISTOGRAMS BEFORE QC
################################################################################

safe_hist <- function(values, ...) {
  values <- values[is.finite(values)]
  if(length(values)) hist(values, ...) else {
    plot.new()
    text(0.5, 0.5, "No finite observations")
  }
}

par(mfrow = c(3,6))

for(cls in unique(dat$Classification)){
  
  for(gen in unique(dat$Genotype)){
    
    x <- dat %>%
      filter(
        Classification == cls,
        Genotype == gen
      )
    
    safe_hist(
      x$Length,
      main = paste(cls, gen, "Length"),
      xlab = "Length",
      col = "lightgray"
    )
    
    safe_hist(
      x$Width,
      main = paste(cls, gen, "Width"),
      xlab = "Width",
      col = "lightgray"
    )
    
    safe_hist(
      x$L_W_Ratio,
      main = paste(cls, gen, "L/W"),
      xlab = "L/W",
      col = "lightgray"
    )
  }
}

################################################################################
# BOXPLOTS BEFORE QC
################################################################################

grain_long_before <- dat %>%
  pivot_longer(
    c(Length, Width, L_W_Ratio),
    names_to = "Trait",
    values_to = "Value"
  )

print(ggplot(
  grain_long_before,
  aes(Genotype, Value, fill = Genotype)
) +
  geom_boxplot() +
  facet_grid(
    Trait ~ Classification,
    scales = "free_y"
  ) +
  theme_bw())

################################################################################
# IQR OUTLIER FUNCTION
################################################################################

iqr_outlier <- function(x){
  
  valid <- is.finite(x) & x > 0
  flags <- rep(FALSE, length(x))
  
  if(!any(valid)) return(flags)
  
  q1 <- quantile(x[valid], .25, names = FALSE, type = 7)
  q3 <- quantile(x[valid], .75, names = FALSE, type = 7)
  
  iqr <- q3 - q1
  
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  
  flags[valid] <- x[valid] < lower | x[valid] > upper
  flags
}

################################################################################
# FLAG OUTLIERS
################################################################################

dat_qc <- dat %>%
  group_by(
    Classification,
    Genotype
  ) %>%
  mutate(
    
    Length_outlier =
      iqr_outlier(Length),
    
    Width_outlier =
      iqr_outlier(Width),
    
    Ratio_outlier =
      iqr_outlier(L_W_Ratio),
    
    Invalid_measurement =
      !is.finite(Length) | !is.finite(Width) | !is.finite(L_W_Ratio) |
      Length <= 0 | Width <= 0 | L_W_Ratio <= 0,
    
    IQR_outlier = Length_outlier | Width_outlier | Ratio_outlier,
    
    Remove = Invalid_measurement | (exclude_outliers & IQR_outlier)
    
  ) %>%
  ungroup()

################################################################################
# QC SUMMARY
################################################################################

cat("\n=================================\n")
cat("QC SUMMARY\n")
cat("=================================\n")

print(table(dat_qc$Remove))

cat(
  "\nRemoved:",
  sum(dat_qc$Remove),
  "\n"
)

cat(
  "Retained:",
  sum(!dat_qc$Remove),
  "\n"
)

################################################################################
# EXCLUDED OBSERVATIONS
################################################################################

extreme_points <- dat_qc %>%
  filter(Remove)

print(
  extreme_points %>%
    arrange(desc(Length)) %>%
    dplyr::select(
      ImageID,
      Classification,
      Genotype,
      Length,
      Width,
      L_W_Ratio
    ),
  n = 20
)

################################################################################
# COUNTS BY GROUP
################################################################################

qc_counts <- dat_qc %>%
  count(Classification, Genotype, Remove)

print(qc_counts, n = Inf)

################################################################################
# CLEAN DATASET
################################################################################

master <- dat_qc %>%
  filter(!Remove)

################################################################################
# HISTOGRAMS AFTER QC
################################################################################

par(mfrow = c(3,6))

for(cls in unique(master$Classification)){
  
  for(gen in unique(master$Genotype)){
    
    x <- master %>%
      filter(
        Classification == cls,
        Genotype == gen
      )
    
    safe_hist(
      x$Length,
      main = paste("QC", cls, gen, "Length"),
      xlab = "Length",
      col = "lightblue"
    )
    
    safe_hist(
      x$Width,
      main = paste("QC", cls, gen, "Width"),
      xlab = "Width",
      col = "lightblue"
    )
    
    safe_hist(
      x$L_W_Ratio,
      main = paste("QC", cls, gen, "L/W"),
      xlab = "L/W",
      col = "lightblue"
    )
  }
}

################################################################################
# BOXPLOTS AFTER QC
################################################################################

grain_long_qc <- master %>%
  pivot_longer(
    c(Length, Width, L_W_Ratio),
    names_to = "Trait",
    values_to = "Value"
  )

print(ggplot(
  grain_long_qc,
  aes(Genotype, Value, fill = Genotype)
) +
  geom_boxplot() +
  facet_grid(
    Trait ~ Classification,
    scales = "free_y"
  ) +
  theme_bw())

################################################################################
# SUMMARY STATISTICS + WELCH T-TEST
################################################################################

summary_table <- grain_long_qc %>%
  group_by(Classification, Trait) %>%
  group_modify(~{
    
    wt <- .x$Value[.x$Genotype == "WT"]
    mutant <- .x$Value[.x$Genotype == "Mutant"]
    
    wt <- wt[is.finite(wt)]
    mutant <- mutant[is.finite(mutant)]
    
    if(length(wt) < 2 || length(mutant) < 2) {
      stop("At least two retained grains per phenotype are required.")
    }
    
    wt_mean <- mean(wt)
    mutant_mean <- mean(mutant)
    
    test <- t.test(
      wt,
      mutant,
      var.equal = FALSE
    )
    
    tibble(
      
      WT_N = length(wt),
      WT_Mean = wt_mean,
      WT_SD = sd(wt, na.rm = TRUE),
      
      Mutant_N = length(mutant),
      Mutant_Mean = mutant_mean,
      Mutant_SD = sd(mutant, na.rm = TRUE),
      
      Difference = mutant_mean - wt_mean,
      Percent_Difference = 100 * (mutant_mean - wt_mean) / wt_mean,
      SG_WT_Ratio = mutant_mean / wt_mean,
      Difference_CI95_Lower = -unname(test$conf.int[2]),
      Difference_CI95_Upper = -unname(test$conf.int[1]),
      
      Overall_N = nrow(.x),
      Overall_Median = median(.x$Value, na.rm = TRUE),
      Overall_Min = min(.x$Value, na.rm = TRUE),
      Overall_Max = max(.x$Value, na.rm = TRUE),
      
      t_value = as.numeric(test$statistic),
      df = as.numeric(test$parameter),
      p_value = test$p.value
    )
  }) %>%
  ungroup()

print(summary_table, n = Inf)

# Differences and CIs: SG minus WT; t statistic: WT minus SG (original ordering).
ttest_table <- summary_table %>%
  mutate(
    FDR_BH = p.adjust(p_value, method = "BH")
  )

ttest_table %>%
  mutate(
    p_value = format(p_value, scientific = TRUE, digits = 10),
    FDR_BH = format(FDR_BH, scientific = TRUE, digits = 10)
  ) %>%
  print(n = Inf)

################################################################################
# SUPPLEMENTAL TABLE S4 - DESCRIPTIVE STATISTICS
################################################################################

summary_tableS4 <- grain_long_qc %>%
  group_by(
    Genotype,
    Classification,
    Trait
  ) %>%
  summarise(
    N = n(),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    SE = SD/sqrt(N),
    Median = median(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tableS4, n = Inf)

################################################################################
# SENSITIVITY TO IQR FILTERING
################################################################################

mean_comparison <- function(long_data) {
  
  long_data %>%
    group_by(Classification, Trait, Genotype) %>%
    summarise(
      N = n(),
      Mean = mean(Value),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Genotype,
      values_from = c(N, Mean),
      names_glue = "{Genotype}_{.value}"
    ) %>%
    mutate(
      Difference = Mutant_Mean - WT_Mean,
      Percent_Difference = 100 * Difference / WT_Mean,
      SG_WT_Ratio = Mutant_Mean / WT_Mean
    )
}

invalid_ids <- dat_qc$Record_ID[dat_qc$Invalid_measurement]

iqr_removed_ids <- dat_qc$Record_ID[
  dat_qc$Invalid_measurement |
    dat_qc$IQR_outlier
]

sensitivity_table <- mean_comparison(
  grain_long_before %>%
    filter(!Record_ID %in% invalid_ids)
) %>%
  rename_with(
    ~ paste0(.x, "_Before"),
    -c(Classification, Trait)
  ) %>%
  full_join(
    mean_comparison(
      grain_long_before %>%
        filter(!Record_ID %in% iqr_removed_ids)
    ) %>%
      rename_with(
        ~ paste0(.x, "_After"),
        -c(Classification, Trait)
      ),
    by = c("Classification", "Trait")
  ) %>%
  mutate(
    Change_percentage_points =
      Percent_Difference_After -
      Percent_Difference_Before
  )

print(sensitivity_table, n = Inf)

################################################################################
# SIMPLE BOXPLOTS
################################################################################

my_colors <- c(
  "WT" = "#4E79A7",
  "Mutant" = "#E15759"
)

fig_1_boxplot <- grain_long_qc %>%
  filter(
    Trait == "Length",
    Classification == "Milled"
  ) %>%
  ggplot(
    aes(x = Genotype, y = Value, fill = Genotype)
  ) +
  geom_boxplot(
    width = 0.55,
    linewidth = 1,
    staplewidth = 0.5,
    outlier.size = 1,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = my_colors,
    breaks = c("WT", "Mutant"),
    labels = c("WT", "SG"),
    name = "Phenotype"
  ) +
  scale_x_discrete(
    limits = c("WT", "Mutant"),
    labels = c("WT", "SG")
  ) +
  labs(
    title = "Grain length in WT and SG RU2002174",
    x = NULL,
    y = "Grain length (mm)"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman"),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    panel.grid.minor = element_blank()
  )

print(fig_1_boxplot)

################################################################################
# EXPORT ONE STATISTICAL TABLE S1 AND SEPARATE QC FILES
################################################################################

# One row per trait and processing type; WT and SG statistics are side by side.
# Publication P-values below 0.001 are displayed as <0.001, including underflowed zeros.
# Unformatted numerical test results remain available in ttest_table in R.
table_S4 <- summary_tableS4 %>%
  mutate(
    Genotype = ifelse(
      Genotype == "Mutant",
      "SG",
      "WT"
    )
  ) %>%
  pivot_wider(
    names_from = Genotype,
    values_from = c(
      N,
      Mean,
      SD,
      SE,
      Median,
      Min,
      Max
    ),
    names_glue = "{Genotype}_{.value}"
  ) %>%
  left_join(
    ttest_table %>%
      select(
        Classification,
        Trait,
        Difference,
        Percent_Difference,
        SG_WT_Ratio,
        Difference_CI95_Lower,
        Difference_CI95_Upper,
        t_value,
        df,
        p_value,
        FDR_BH
      ),
    by = c("Classification", "Trait")
  ) %>%
  mutate(
    Unit = ifelse(
      Trait == "L_W_Ratio",
      "ratio",
      "mm"
    )
  ) %>%
  arrange(
    match(
      Classification,
      c("Rough", "Brown", "Milled")
    ),
    match(
      Trait,
      c("Length", "Width", "L_W_Ratio")
    )
  ) %>%
  transmute(
    Processing_type = Classification,
    Trait = ifelse(
      Trait == "L_W_Ratio",
      "Length-to-width ratio",
      Trait
    ),
    Unit,
    WT_N,
    WT_Mean,
    WT_SD,
    WT_SE,
    WT_Median,
    WT_Min,
    WT_Max,
    SG_N,
    SG_Mean,
    SG_SD,
    SG_SE,
    SG_Median,
    SG_Min,
    SG_Max,
    Difference_SG_minus_WT = Difference,
    Percent_difference = Percent_Difference,
    SG_WT_ratio = SG_WT_Ratio,
    CI95_lower = Difference_CI95_Lower,
    CI95_upper = Difference_CI95_Upper,
    t_WT_minus_SG = t_value,
    df,
    P_value = ifelse(
      p_value < 0.001,
      "<0.001",
      format(
        p_value,
        scientific = TRUE,
        digits = 5
      )
    ),
    P_BH = ifelse(
      FDR_BH < 0.001,
      "<0.001",
      format(
        FDR_BH,
        scientific = TRUE,
        digits = 5
      )
    )
  )

# Full phenotype names and plain-language column labels in the exported table.
names(table_S4) <- c(
  "Rice processing type",
  "Grain trait",
  "Measurement unit",
  "Wild type (WT): number of grains",
  "Wild type (WT): mean",
  "Wild type (WT): standard deviation",
  "Wild type (WT): standard error",
  "Wild type (WT): median",
  "Wild type (WT): minimum",
  "Wild type (WT): maximum",
  "Short-grain mutant (SG): number of grains",
  "Short-grain mutant (SG): mean",
  "Short-grain mutant (SG): standard deviation",
  "Short-grain mutant (SG): standard error",
  "Short-grain mutant (SG): median",
  "Short-grain mutant (SG): minimum",
  "Short-grain mutant (SG): maximum",
  "Mean difference (SG minus WT)",
  "Change relative to WT (%)",
  "Mean ratio (SG divided by WT)",
  "95% CI for difference: lower limit",
  "95% CI for difference: upper limit",
  "Welch t statistic (WT minus SG)",
  "Welch degrees of freedom",
  "Welch test P-value",
  "BH-adjusted P-value"
)

print(
  table_S4,
  n = Inf,
  width = Inf
)

# One worksheet, one table, with a title and explanatory notes.
wb <- createWorkbook()

addWorksheet(
  wb,
  "Table S4",
  gridLines = FALSE
)

writeData(
  wb,
  "Table S4",
  "S4 Table. Grain morphology of wild-type RU2002174 and the RU2002174 short-grain mutant.",
  startRow = 2,
  colNames = FALSE
)

mergeCells(
  wb,
  "Table S4",
  cols = 1:10,
  rows = 2
)

writeData(
  wb,
  "Table S4",
  "WT = wild type; SG = short-grain mutant (coded as Mutant in the input).",
  startRow = 3,
  colNames = FALSE
)

mergeCells(
  wb,
  "Table S4",
  cols = 1:10,
  rows = 3
)

writeData(
  wb,
  "Table S4",
  table_S4,
  startRow = 5
)

header_style <- createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#294C60",
  textDecoration = "bold",
  wrapText = TRUE,
  halign = "center",
  valign = "center"
)

addStyle(
  wb,
  "Table S4",
  header_style,
  rows = 5,
  cols = 1:26,
  gridExpand = TRUE
)

setRowHeights(
  wb,
  "Table S4",
  rows = 5,
  heights = 65
)

setColWidths(
  wb,
  "Table S4",
  cols = 1:26,
  widths = 27
)

setColWidths(
  wb,
  "Table S4",
  cols = 2,
  widths = 25
)

addStyle(
  wb,
  "Table S4",
  createStyle(numFmt = "0.0000"),
  rows = 6:14,
  cols = 4:24,
  gridExpand = TRUE,
  stack = TRUE
)

addStyle(
  wb,
  "Table S4",
  createStyle(numFmt = "#,##0"),
  rows = 6:14,
  cols = c(4, 11),
  gridExpand = TRUE,
  stack = TRUE
)

for(rows in list(6:8, 12:14)) {
  addStyle(
    wb,
    "Table S4",
    createStyle(fgFill = "#F1F5F7"),
    rows = rows,
    cols = 1:26,
    gridExpand = TRUE,
    stack = TRUE
  )
}

notes <- c(
  "Length and width are in mm; length-to-width ratio is dimensionless. All values are after QC with the selected filtering setting.",
  "Mean difference = SG mean minus WT mean; change (%) = 100 x (SG mean - WT mean) / WT mean. Negative values indicate a lower SG mean.",
  "The 95% confidence interval describes the SG-minus-WT mean difference. The Welch t statistic retains the WT-minus-SG direction.",
  "BH-adjusted P-values account for nine comparisons. Values below 0.001 are displayed as <0.001.",
  "Individual grains are subsamples from bulk seed lots, not independent biological replicates."
)

for(i in seq_along(notes)) {
  writeData(
    wb,
    "Table S4",
    notes[i],
    startRow = 16 + i,
    colNames = FALSE
  )
  
  mergeCells(
    wb,
    "Table S4",
    cols = 1:10,
    rows = 16 + i
  )
}

freezePane(
  wb,
  "Table S4",
  firstActiveRow = 6,
  firstActiveCol = 4
)

saveWorkbook(
  wb,
  file.path(
    output_dir,
    "S4_Table_Grain_Morphology_Statistics.xlsx"
  ),
  overwrite = TRUE
)


################################################################################
# END OF SCRIPT
################################################################################
