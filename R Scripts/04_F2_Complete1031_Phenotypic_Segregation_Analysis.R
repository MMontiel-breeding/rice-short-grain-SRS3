################################################################################
# 
# Manuscript: Genetic mapping of a putative spontaneous short-grain mutation
# identifies a predicted loss-of-function allele of SRS3 in rice
#
# Author: Maria Guadalupe Montiel
# 
# Script 04 | F2 Complete 1031 Population Phenotype Segregation
#
################################################################################

# 1 | START CLEAN AND LOAD PACKAGES --------------------------------------------
rm(list = ls())

required_packages <- c("readxl", "dplyr", "tidyr", "stringr", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE,
          FUN.VALUE = logical(1))]

if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# 2 | ANALYSIS SETTINGS --------------------------------------------------------
input_file <- "S2_Data_Complete_1031_F2_Phenotypes.xlsx"
input_sheet <- "Sheet1"
expected_total_records <- 1031L  # -> Use the finalized file without the four unclassified records.

expected_wt_proportion <- 0.75
expected_sg_proportion <- 0.25
alpha <- 0.05
confidence_level <- 0.95

# -> All plants are analyzed together; no pedigree-specific analyses.
# -> No setwd(), automatic exports or changes to the original workbook.

# 3 | STANDARDIZE LABELS -------------------------------------------------------
clean_text <- function(x) {
  x <- str_squish(as.character(x))
  x[x == ""] <- NA_character_
  x
}

standardize_phenotype <- function(x) {
  y <- str_to_lower(clean_text(x))
  case_when(
    y %in% c("wild", "wt", "wild type", "wild-type",
             "wildtype", "long", "normal") ~ "WT",
    y %in% c("mutant", "sg", "short", "short grain", "short-grain") ~ "SG",
    TRUE ~ NA_character_
  )
}

format_p_for_paper <- function(p) {
  if (p < 0.001) {
    format(p, scientific = TRUE, digits = 3)
  } else {
    formatC(p, format = "f", digits = 3)
  }
}

# 4 | IMPORT THE PHENOTYPE DATA ------------------------------------------------
if (!file.exists(input_file)) stop("Input file not found; update input_file.")
if (!input_sheet %in% excel_sheets(input_file)) {
  stop("Input sheet not found; update input_sheet.")
}

raw_data <- read_excel(input_file, sheet = input_sheet)
required_columns <- c("Line ID", "Phenotype")
missing_columns <- setdiff(required_columns, names(raw_data))

if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

analysis_data <- raw_data %>% mutate(
  Record_ID = row_number(),
  Line_ID = clean_text(`Line ID`),
  Phenotype_original = clean_text(Phenotype),
  Phenotype_standard = standardize_phenotype(Phenotype)
)
# -> Original columns are preserved; only plant ID and phenotype are used for testing.

# 5 | VALIDATE BEFORE ANALYSIS -------------------------------------------------
validation_table <- tibble(
  Validation_item = c("Expected records", "Observed records",
                      "Valid WT/SG phenotypes", "Unclassified phenotypes",
                      "Missing plant IDs", "Duplicated plant IDs"),
  Value = c(expected_total_records, nrow(analysis_data),
            sum(!is.na(analysis_data$Phenotype_standard)),
            sum(is.na(analysis_data$Phenotype_standard)),
            sum(is.na(analysis_data$Line_ID)),
            sum(duplicated(analysis_data$Line_ID)))
)
print(validation_table, n = Inf)

problematic_records <- analysis_data %>% filter(
  is.na(Line_ID) | is.na(Phenotype_standard) |
    duplicated(Line_ID) | duplicated(Line_ID, fromLast = TRUE)
)

if (nrow(problematic_records)) {
  print(problematic_records, n = Inf, width = Inf)
  stop("Review these records. No plants were automatically excluded or assigned to WT.")
}

if (nrow(analysis_data) != expected_total_records) {
  stop("Expected ", expected_total_records, " records; found ", nrow(analysis_data),
       ". Confirm the finalized input dataset.")
}

expected_proportions <- c(WT = expected_wt_proportion, SG = expected_sg_proportion)
if (any(!is.finite(expected_proportions)) ||
    any(expected_proportions <= 0) ||
    abs(sum(expected_proportions) - 1) > 1e-10) {
  stop("Expected proportions must be finite, positive and sum to one.")
}

cat("\nValidation passed:", nrow(analysis_data),
    "unique plants with valid WT/SG classifications.\n")
# -> The initial 147 plants are already included and must not be added again.

# 6 | COUNTS AND FREQUENCIES ---------------------------------------------------
overall_summary <- analysis_data %>%
  count(Phenotype_standard, name = "Observed_count") %>%
  complete(Phenotype_standard = c("WT", "SG"), fill = list(Observed_count = 0)) %>%
  mutate(
    Total = sum(Observed_count),
    Frequency = Observed_count / Total,
    Percentage = 100 * Frequency,
    Expected_proportion = if_else(Phenotype_standard == "WT",
                                  expected_wt_proportion, expected_sg_proportion),
    Expected_count = Total * Expected_proportion,
    Observed_minus_expected = Observed_count - Expected_count
  ) %>%
  arrange(factor(Phenotype_standard, levels = c("WT", "SG")))

print(overall_summary, n = Inf, width = Inf)

observed_counts <- setNames(
  overall_summary$Observed_count[
    match(c("WT", "SG"), overall_summary$Phenotype_standard)],
  c("WT", "SG")
)

# 7 | PRIMARY TEST: PEARSON CHI-SQUARE -----------------------------------------
overall_chi <- chisq.test(
  observed_counts, p = expected_proportions, correct = FALSE
)
# -> Goodness-of-fit test against the expected 3 WT:1 SG segregation ratio.

overall_chi_results <- tibble(
  Total_F2 = sum(observed_counts),
  Observed_WT = unname(observed_counts["WT"]),
  Observed_SG = unname(observed_counts["SG"]),
  Expected_WT = sum(observed_counts) * expected_wt_proportion,
  Expected_SG = sum(observed_counts) * expected_sg_proportion,
  Chi_square = unname(overall_chi$statistic),
  Degrees_of_freedom = unname(overall_chi$parameter),
  P_value = overall_chi$p.value,
  Alpha = alpha,
  Statistical_conclusion = if_else(
    overall_chi$p.value < alpha,
    "Significant deviation from 3:1",
    "No significant deviation from 3:1"
  )
)
print(overall_chi_results, width = Inf)

# 8 | COMPLEMENTARY TEST: EXACT BINOMIAL ---------------------------------------
exact_binomial <- binom.test(
  unname(observed_counts["SG"]), sum(observed_counts),
  p = expected_sg_proportion, alternative = "two.sided",
  conf.level = confidence_level
)

exact_results <- tibble(
  Total_F2 = sum(observed_counts),
  SG_count = unname(observed_counts["SG"]),
  Observed_SG_frequency = unname(observed_counts["SG"] / sum(observed_counts)),
  Expected_SG_frequency = expected_sg_proportion,
  Exact_binomial_P_value = exact_binomial$p.value,
  Confidence_level = confidence_level,
  Lower_CI = unname(exact_binomial$conf.int[1]),
  Upper_CI = unname(exact_binomial$conf.int[2])
)
print(exact_results, width = Inf)
# -> An exact check of the same segregation expectation, not an independent experiment.
# -> The confidence interval describes the SG proportion in the complete population.

# 9 | DISPLAY PHENOTYPE COUNTS -------------------------------------------------
plot_data <- overall_summary %>% mutate(
  Phenotype_standard = factor(Phenotype_standard, levels = c("WT", "SG")),
  Count_label = paste0(Observed_count, "\n(",
                       sprintf("%.1f", Percentage), "%)")
)

segregation_plot <- ggplot(
  plot_data, aes(Phenotype_standard, Observed_count, fill = Phenotype_standard)
) +
  geom_col(width = 0.62) +
  geom_text(aes(label = Count_label), vjust = -0.35, size = 4.5) +
  scale_fill_manual(values = c(WT = "#0072B2", SG = "#D55E00")) +
  labs(
    title = "Phenotypic segregation in the complete IC206-derived F2 population",
    subtitle = paste0(
      "Expected 3 WT : 1 SG; chi-square = ",
      formatC(overall_chi_results$Chi_square, format = "f", digits = 2),
      ", df = ", overall_chi_results$Degrees_of_freedom,
      ", P = ", format_p_for_paper(overall_chi_results$P_value)
    ),
    x = "Phenotype", y = "Number of F2 plants"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold")) +
  expand_limits(y = max(plot_data$Observed_count) * 1.12)

print(segregation_plot)  # -> Display only; no automatic saving.

# 10 | RESULTS WORDING AND PROCESSING SUMMARY ----------------------------------
wt_count <- unname(observed_counts["WT"])
sg_count <- unname(observed_counts["SG"])
total_count <- sum(observed_counts)

results_verb <- ifelse(
  overall_chi$p.value < alpha,
  "differed significantly from", "did not differ significantly from"
)

results_text <- sprintf(
  paste0(
    "Among %d phenotyped F2 plants, %d were WT (%.2f%%) and %d were SG (%.2f%%). ",
    "Observed segregation %s the expected 3 WT:1 SG ratio ",
    "(expected counts: %.2f WT and %.2f SG; chi-square = %.2f, df = %d, P = %s)."
  ),
  total_count, wt_count, 100 * wt_count / total_count,
  sg_count, 100 * sg_count / total_count, results_verb,
  total_count * expected_wt_proportion, total_count * expected_sg_proportion,
  unname(overall_chi$statistic), unname(overall_chi$parameter),
  format_p_for_paper(overall_chi$p.value)
)
cat("\n", results_text, "\n", sep = "")

################################################################################
# END OF SCRIPT
################################################################################