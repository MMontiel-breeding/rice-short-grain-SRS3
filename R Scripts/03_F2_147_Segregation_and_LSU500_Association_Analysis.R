################################################################################
# Paper: Genetic mapping of a putative spontaneous short-grain mutation identifies
# a predicted loss-of-function allele of SRS3 in rice
#
# Author: Maria Guadalupe Montiel
# Script 02 | F2 LSU500 single-marker association analysis
#
# Workflow: Import -> Validate -> Filter -> Associate -> Adjust -> Summarize
# Methods: 3:1 segregation test, categorical ANOVA, binomial logistic LRT,
#          Benjamini-Hochberg FDR and Bonferroni correction.
# Results are displayed in R. No output files are saved automatically.
################################################################################

rm(list = ls())                                      # -> Start with a clean workspace.

# ── 1. PREPARE THE ANALYSIS ─────────────────────────────────────────────────────

required_packages <- c("readxl", "dplyr", "stringr", "ggplot2", "tibble", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Install: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr)
  library(ggplot2); library(tibble); library(tidyr)
})

# ── 2. SELECT INPUTS AND FILTERS ────────────────────────────────────────────────

input_file <- file.choose()                          # -> Select the revised F2 workbook.
marker_file <- file.choose()                         # -> Select the revised marker dictionary.
sheet_name <- 1

# For the repository, replace file.choose() with the corresponding relative paths.
minimum_call_rate <- 0.95                            # -> Valid calls / all evaluated plants.
apply_maf_filter <- TRUE                             # -> FALSE removes the MAF requirement.
minimum_maf <- 0.20                                  # -> Retain markers with MAF > 0.20.
fdr_alpha <- 0.05 ##--> we can do 0.001
expected_wt_proportion <- 0.75
expected_sg_proportion <- 0.25                        # -> Expected recessive segregation.

if (!all(file.exists(c(input_file, marker_file)))) stop("Check the input file paths.")

# ── 3. DEFINE VALIDATION AND ANALYSIS FUNCTIONS ─────────────────────────────────

clean_text <- function(x) {
  x <- str_trim(as.character(x)); x[x == ""] <- NA_character_
  x
}

canonical_phenotype <- function(x) {
  y <- str_to_lower(clean_text(x))
  case_when(
    y %in% c("mutant", "short", "sg", "short grain", "short-grain") ~ "SG",
    y %in% c("wild", "normal", "long", "wt", "wild type", "wild-type", "wildtype") ~ "WT",
    TRUE ~ NA_character_
  )
}

normalize_genotype <- function(x, allele1, allele2) {
  x <- str_replace_all(toupper(clean_text(x)), "\\s+", "")
  x[x %in% c("FAIL", "FAILED", "NOCALL", "NO_CALL", "NA", "N/A", ".", "-", "--")] <- NA_character_
  
  vapply(x, function(g) {
    if (is.na(g)) return(NA_character_)
    z <- if (str_detect(g, fixed("/"))) strsplit(g, "/", fixed = TRUE)[[1]] else rep(g, 2)
    if (length(z) != 2 || any(!z %in% c(allele1, allele2))) return("INVALID")
    paste(sort(z), collapse = "/")                   # -> A/G = G/A; A becomes A/A.
  }, character(1))
}

maf_from_calls <- function(g, allele1, allele2) {
  g <- g[!is.na(g)]
  if (!length(g)) return(NA_real_)
  alleles <- unlist(strsplit(g, "/", fixed = TRUE))
  min(sum(alleles == allele1), sum(alleles == allele2)) / length(alleles)
}

safe_anova <- function(pheno_binary, genotype) {
  d <- tibble(y = pheno_binary, genotype = factor(genotype)) %>%
    filter(!is.na(y), !is.na(genotype)) %>% droplevels()
  
  if (nrow(d) < 10 || nlevels(d$genotype) < 2) {
    return(tibble(anova_n = nrow(d), anova_df = NA_real_, anova_F = NA_real_,
                  anova_p = NA_real_, anova_R2 = NA_real_))
  }
  
  fit <- lm(y ~ genotype, data = d)                  # -> Historical categorical ANOVA.
  a <- anova(fit)
  ss_marker <- a$`Sum Sq`[1]; ss_total <- sum(a$`Sum Sq`)
  
  tibble(anova_n = nrow(d), anova_df = a$Df[1], anova_F = a$`F value`[1],
         anova_p = a$`Pr(>F)`[1],
         anova_R2 = ifelse(ss_total > 0, ss_marker / ss_total, NA_real_))
}

safe_logistic <- function(pheno_binary, genotype) {
  d <- tibble(y = pheno_binary, genotype = factor(genotype)) %>%
    filter(!is.na(y), !is.na(genotype)) %>% droplevels()
  
  if (nrow(d) < 10 || nlevels(d$genotype) < 2 || length(unique(d$y)) < 2) {
    return(tibble(logistic_n = nrow(d), logistic_df = NA_real_, logistic_LR_chisq = NA_real_,
                  logistic_p = NA_real_, mcfadden_R2 = NA_real_, boundary_class = NA))
  }
  
  counts <- d %>% group_by(genotype) %>% summarise(n = n(), sg = sum(y), .groups = "drop")
  
  # -> Class-count likelihood equals the categorical logistic-model likelihood.
  # This also handles genotype classes containing only WT or only SG plants.
  binomial_loglik <- function(sg, n, p) {
    wt <- n - sg; value <- numeric(length(n))
    i <- sg > 0; value[i] <- value[i] + sg[i] * log(p[i])
    i <- wt > 0; value[i] <- value[i] + wt[i] * log1p(-p[i])
    sum(value)
  }
  
  ll_full <- binomial_loglik(counts$sg, counts$n, counts$sg / counts$n)
  ll_null <- binomial_loglik(sum(d$y), nrow(d), mean(d$y))
  lr <- max(0, 2 * (ll_full - ll_null)); df_test <- nrow(counts) - 1
  
  tibble(logistic_n = nrow(d), logistic_df = df_test, logistic_LR_chisq = lr,
         logistic_p = pchisq(lr, df = df_test, lower.tail = FALSE),
         mcfadden_R2 = 1 - ll_full / ll_null,
         boundary_class = any(counts$sg == 0 | counts$sg == counts$n))
}

# ── 4. IMPORT AND VALIDATE THE DATA ─────────────────────────────────────────────

# -> The supplied workbook has one header row followed by the F2 plants.
raw_data <- read_excel(input_file, sheet = sheet_name, .name_repair = "minimal")
if (anyDuplicated(names(raw_data))) stop("Duplicate column names in the genotype workbook.")
if (ncol(raw_data) < 3 || names(raw_data)[1] != "Name" ||
    !str_detect(names(raw_data)[2], regex("phenotype", ignore_case = TRUE))) {
  stop("Expected one header row: Name, Phenotype, then marker IDs.")
}

marker_ids <- names(raw_data)[-(1:2)]
if (any(is.na(marker_ids) | marker_ids == "")) stop("Missing marker IDs.")

analysis_data <- raw_data %>% mutate(
  Sample_ID = clean_text(.data[[names(raw_data)[1]]]),
  Phenotype_original = .data[[names(raw_data)[2]]],
  Phenotype = canonical_phenotype(Phenotype_original),
  Phenotype_binary = case_when(Phenotype == "SG" ~ 1, Phenotype == "WT" ~ 0, TRUE ~ NA_real_)
)

if (anyNA(analysis_data$Sample_ID) || anyDuplicated(analysis_data$Sample_ID)) {
  stop("Missing or duplicate sample IDs.")
}
if (anyNA(analysis_data$Phenotype)) stop("Unrecognized phenotype labels.")
analysis_f2 <- analysis_data
if (nrow(analysis_f2) != 147) warning("Expected 147 plants; observed ", nrow(analysis_f2), ".")

marker_reference <- read_excel(marker_file, sheet = 1) %>%
  dplyr::rename(
    Marker_ID = Agriplex_ID_LSU500,
    IRGSPv_Position = `IRGSPv_Pos (bp)`
  )

required_fields <- c("Marker_ID", "Chr", "IRGSPv_Position", "Allele1", "Allele2")

if (!all(required_fields %in% names(marker_reference))) {
  stop("Missing marker-dictionary fields.")
}

marker_reference <- marker_reference %>% transmute(
  marker_id = clean_text(Marker_ID), chromosome = as.integer(Chr),
  position_raw = as.character(IRGSPv_Position),
  allele1 = toupper(clean_text(Allele1)), allele2 = toupper(clean_text(Allele2))
)
if (anyNA(marker_reference$marker_id) || anyDuplicated(marker_reference$marker_id)) {
  stop("Missing or duplicate dictionary IDs.")
}

marker_dictionary <- tibble(marker_column = marker_ids, marker_id = marker_ids) %>%
  left_join(marker_reference, by = "marker_id")
if (anyNA(marker_dictionary$allele1) || anyNA(marker_dictionary$allele2)) {
  stop("Some input markers lack documented alleles.")
}

# -> Preserve original coordinates. Unverified decimal values are not assumed to be bp.
# Missing physical positions do not exclude markers from association testing.
marker_dictionary <- marker_dictionary %>% mutate(
  position_value = suppressWarnings(as.numeric(position_raw)),
  position_bp = ifelse(!is.na(position_value) & position_value > 100 &
                         position_value == floor(position_value), position_value, NA_real_)
) %>% select(marker_column, marker_id, chromosome, position_raw, position_bp, allele1, allele2)

invalid_calls <- list()
for (i in seq_len(nrow(marker_dictionary))) {
  m <- marker_dictionary$marker_column[i]; original <- analysis_f2[[m]]
  calls <- normalize_genotype(original, marker_dictionary$allele1[i], marker_dictionary$allele2[i])
  bad <- which(!is.na(calls) & calls == "INVALID")
  if (length(bad)) invalid_calls[[m]] <- tibble(
    Sample_ID = analysis_f2$Sample_ID[bad], marker_id = m, original_call = as.character(original[bad])
  )
  analysis_f2[[m]] <- calls
}
invalid_calls <- bind_rows(invalid_calls)
if (nrow(invalid_calls)) {
  print(invalid_calls, n = 20); stop("Undocumented genotype calls; inspect before proceeding.")
}

cat("\n--- PHENOTYPE VALIDATION ---\n")
print(table(analysis_f2$Phenotype_original)); print(table(analysis_f2$Phenotype))

# ── 5. CHECK THE EXPECTED 3:1 SEGREGATION ────────────────────────────────────────

phenotype_counts <- tibble(
  Phenotype = c("WT", "SG"),
  Observed = c(sum(analysis_f2$Phenotype == "WT"), sum(analysis_f2$Phenotype == "SG")),
  Expected_proportion = c(expected_wt_proportion, expected_sg_proportion)
) %>% mutate(Expected = sum(Observed) * Expected_proportion)

chi_test <- chisq.test(phenotype_counts$Observed, p = phenotype_counts$Expected_proportion)
chi_results <- tibble(
  total_F2 = nrow(analysis_f2), observed_WT = phenotype_counts$Observed[1],
  observed_SG = phenotype_counts$Observed[2], chi_square = unname(chi_test$statistic),
  df = unname(chi_test$parameter), p_value = chi_test$p.value
)

cat("\n--- SEGREGATION TEST ---\n")
print(phenotype_counts); print(chi_results, width = Inf)

plot_data <- phenotype_counts %>% select(Phenotype, Observed, Expected) %>%
  pivot_longer(c(Observed, Expected), names_to = "Count_type", values_to = "Count")

segregation_plot <- ggplot(plot_data, aes(Phenotype, Count, fill = Count_type)) +
  geom_col(position = "dodge", width = 0.7) +
  theme_classic() +
  labs(
    title = "F2 phenotype segregation",
    subtitle = sprintf(
      "Expected ratio: 3 WT : 1 SG\nχ² = %.2f; df = %d; P = %.3f",
      chi_results$chi_square, chi_results$df, chi_results$p_value
    ),
    x = "Phenotype",
    y = "Number of plants",
    fill = NULL
  )

print(segregation_plot)

# ── 6. FILTER MARKERS AND SUMMARIZE COMPLETENESS ────────────────────────────────

# -> Failed and missing calls remain in the call-rate denominator.
# Every marker is evaluated with the same filtering criteria.
marker_qc <- lapply(seq_len(nrow(marker_dictionary)), function(i) {
  m <- marker_dictionary$marker_column[i]; g <- analysis_f2[[m]]; called <- !is.na(g)
  tibble(marker_column = m, n_total = length(g), n_called = sum(called),
         call_rate = mean(called), n_genotype_classes = length(unique(g[called])),
         maf = maf_from_calls(g, marker_dictionary$allele1[i], marker_dictionary$allele2[i]))
}) %>% bind_rows() %>% left_join(marker_dictionary, by = "marker_column") %>% mutate(
  passes_call_rate = call_rate >= minimum_call_rate,
  passes_maf = if (apply_maf_filter) !is.na(maf) & maf > minimum_maf else TRUE,
  passes_all_filters = passes_call_rate & passes_maf & n_genotype_classes >= 2,
  filter_reason = case_when(
    !passes_call_rate ~ "Low call rate", n_genotype_classes < 2 ~ "Monomorphic/no valid calls",
    !passes_maf ~ "MAF filter", TRUE ~ "Tested"
  )
)

cat("\n--- MARKER FILTERING ---\n")
print(marker_qc %>% count(filter_reason), n = Inf)
markers_for_analysis <- marker_qc %>% filter(passes_all_filters) %>% pull(marker_column)
if (!length(markers_for_analysis)) stop("No markers passed the selected filters.")

sample_qc <- tibble(
  Sample_ID = analysis_f2$Sample_ID, Phenotype = analysis_f2$Phenotype,
  completeness_QC_markers = rowMeans(!is.na(as.data.frame(
    analysis_f2[marker_qc$marker_column[marker_qc$passes_call_rate]]
  )))
)
cat("\nLowest sample completeness on QC-retained markers; no plants removed:\n")
print(sample_qc %>% arrange(completeness_QC_markers) %>% slice_head(n = 5))

# ── 7. TEST GENOTYPE–PHENOTYPE ASSOCIATION ──────────────────────────────────────

# -> Genotype is categorical; no additive, dominant or recessive model is imposed.
all_results <- lapply(markers_for_analysis, function(m) {
  marker_dictionary %>% filter(marker_column == m) %>% bind_cols(
    safe_anova(analysis_f2$Phenotype_binary, analysis_f2[[m]]),
    safe_logistic(analysis_f2$Phenotype_binary, analysis_f2[[m]])
  )
}) %>% bind_rows()

if (anyNA(all_results$anova_p) || anyNA(all_results$logistic_p)) {
  stop("A marker test failed; inspect all_results before multiple-testing correction.")
}

# -> Adjust each test family separately using all tested markers.
all_results <- all_results %>% mutate(
  anova_FDR = p.adjust(anova_p, "BH"), anova_Bonferroni = p.adjust(anova_p, "bonferroni"),
  logistic_FDR = p.adjust(logistic_p, "BH"), logistic_Bonferroni = p.adjust(logistic_p, "bonferroni"),
  logistic_neglog10P = -log10(pmax(logistic_p, .Machine$double.xmin))
) %>% left_join(
  marker_qc %>% select(marker_column, call_rate, maf, n_called, n_genotype_classes),
  by = "marker_column"
)

# -> SG proportions by genotype provide interpretable effect estimates.
# Exact 95% binomial intervals are descriptive and not adjusted simultaneous intervals.
genotype_effects <- lapply(markers_for_analysis, function(m) {
  tibble(genotype = analysis_f2[[m]], SG = analysis_f2$Phenotype_binary) %>%
    filter(!is.na(genotype)) %>% group_by(genotype) %>%
    summarise(n = n(), SG_count = sum(SG), .groups = "drop") %>% mutate(
      marker_id = m, WT_count = n - SG_count, SG_proportion = SG_count / n,
      SG_CI_lower = mapply(function(s, n) binom.test(s, n)$conf.int[1], SG_count, n),
      SG_CI_upper = mapply(function(s, n) binom.test(s, n)$conf.int[2], SG_count, n)
    )
}) %>% bind_rows()

effect_ranges <- genotype_effects %>% group_by(marker_id) %>%
  summarise(SG_proportion_range = max(SG_proportion) - min(SG_proportion), .groups = "drop")
all_results <- all_results %>% left_join(effect_ranges, by = "marker_id")

anova_ranked <- all_results %>% arrange(anova_p)
logistic_ranked <- all_results %>% arrange(logistic_p)
r2_ranked <- all_results %>% arrange(desc(anova_R2))
significant_markers <- logistic_ranked %>%
  filter(logistic_FDR < fdr_alpha & logistic_Bonferroni < fdr_alpha)

cat("\n--- TOP MARKERS: ANOVA ---\n")
print(anova_ranked %>% select(
  marker_id, chromosome, position_bp, anova_p, anova_R2, anova_FDR, anova_Bonferroni
) %>% slice_head(n = 10), width = Inf)

cat("\n--- TOP MARKERS: LOGISTIC LRT ---\n")
print(logistic_ranked %>% select(
  marker_id, chromosome, position_bp, logistic_p, logistic_FDR,
  logistic_Bonferroni, mcfadden_R2, boundary_class
) %>% slice_head(n = 10), width = Inf)

cat("\n--- SIGNIFICANT MARKERS AND GENOTYPE EFFECTS ---\n")
print(significant_markers, n = Inf, width = Inf)
print(genotype_effects %>% filter(marker_id %in% significant_markers$marker_id),
      n = Inf, width = Inf)

# ── 8. DISPLAY THE GENOME-WIDE ASSOCIATION SCAN ─────────────────────────────────

# -> Unresolved physical positions affect plotting only, not association testing.
plot_markers <- logistic_ranked %>% filter(!is.na(chromosome), !is.na(position_bp)) %>%
  arrange(chromosome, position_bp)
bonferroni_threshold <- fdr_alpha / nrow(all_results) # -> All tested markers, not only plotted ones.

if (nrow(plot_markers)) {
  chr_info <- plot_markers %>% group_by(chromosome) %>%
    summarise(chr_length = max(position_bp), .groups = "drop") %>%
    arrange(chromosome) %>% mutate(offset = lag(cumsum(chr_length), default = 0))
  plot_markers <- plot_markers %>% left_join(chr_info, by = "chromosome") %>%
    mutate(cumulative_position = position_bp + offset)
  axis_positions <- plot_markers %>% group_by(chromosome) %>%
    summarise(center = (min(cumulative_position) + max(cumulative_position)) / 2, .groups = "drop")
  
  marker_plot <- ggplot(plot_markers,
                        aes(cumulative_position, logistic_neglog10P, color = factor(chromosome %% 2))) +
    geom_point(size = 1.6) +
    geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("#4472C4", "#70AD47"), guide = "none") +
    scale_x_continuous(breaks = axis_positions$center, labels = axis_positions$chromosome) +
    theme_classic() +
    labs(title = "F2 single-marker logistic association",
         subtitle = paste("Bonferroni threshold based on", nrow(all_results), "tested markers"),
         x = "Chromosome", y = expression(-log[10](italic(P))))
  print(marker_plot)
}

# ── 9. SUMMARIZE FILTERING AND RESULTS ──────────────────────────────────────────

# -> ANOVA R2 concerns the observed 0/1 phenotype; McFadden R2 is likelihood-based.
# Neither measures heritability or quantitative grain-length variation explained.
analysis_summary <- tibble(
  Samples = nrow(analysis_f2), WT = sum(analysis_f2$Phenotype == "WT"),
  SG = sum(analysis_f2$Phenotype == "SG"), Original_markers = nrow(marker_qc),
  Call_rate_retained = sum(marker_qc$passes_call_rate),
  Variable_after_call_QC = sum(marker_qc$passes_call_rate & marker_qc$n_genotype_classes >= 2),
  MAF_filter_used = apply_maf_filter, Tested_markers = nrow(all_results),
  ANOVA_FDR_significant = sum(all_results$anova_FDR < fdr_alpha),
  ANOVA_Bonferroni_significant = sum(all_results$anova_Bonferroni < fdr_alpha),
  Logistic_FDR_significant = sum(all_results$logistic_FDR < fdr_alpha),
  Logistic_Bonferroni_significant = sum(all_results$logistic_Bonferroni < fdr_alpha)
)

cat("\n--- ANALYSIS SUMMARY ---\n")
print(analysis_summary, width = Inf)

cat(sprintf(
  "\nFiltering: %d plants and %d original markers. No markers excluded by identifier. ",
  analysis_summary$Samples, analysis_summary$Original_markers
))
cat(sprintf(
  "Valid calls / all plants defined call rate. %d markers passed call rate >= %.2f; %d were variable. ",
  analysis_summary$Call_rate_retained, minimum_call_rate, analysis_summary$Variable_after_call_QC
))
cat(if (apply_maf_filter) paste0("MAF > ", minimum_maf, " was additionally required. ")
    else "No MAF filter was applied. ")
cat(sprintf("%d markers were tested. No plants were removed.\n", nrow(all_results)))

cat(sprintf(
  "Results: %d WT and %d SG; 3:1 segregation chi-square = %.4f, df = %d, P = %.6g. ",
  chi_results$observed_WT, chi_results$observed_SG, chi_results$chi_square,
  chi_results$df, chi_results$p_value
))
cat(sprintf(
  "Logistic LRT detected %d BH-FDR-significant and %d Bonferroni-significant markers at alpha %.2f.\n",
  analysis_summary$Logistic_FDR_significant,
  analysis_summary$Logistic_Bonferroni_significant, fdr_alpha
))
cat("Full tables remain in marker_qc, all_results, genotype_effects and the ranked tables.\n")
cat("No output files were saved automatically.\n")

# -> Distinguish variable markers before MAF filtering from those ultimately tested.
markers_by_chromosome <- marker_qc %>% group_by(chromosome) %>% summarise(
  Original_markers = n(), Call_rate_retained = sum(passes_call_rate),
  Variable_after_call_QC = sum(passes_call_rate & n_genotype_classes >= 2),
  Tested_after_all_filters = sum(passes_all_filters), .groups = "drop"
) %>% arrange(chromosome)

cat("\n--- MARKER DISTRIBUTION BY CHROMOSOME ---\n")
print(markers_by_chromosome, n = Inf, width = Inf)

chromosome_summary <- markers_by_chromosome %>% filter(!is.na(chromosome)) %>% summarise(
  Chromosomes = n(), Variable_min = min(Variable_after_call_QC),
  Variable_max = max(Variable_after_call_QC), Variable_mean = mean(Variable_after_call_QC),
  Tested_min = min(Tested_after_all_filters), Tested_max = max(Tested_after_all_filters),
  Tested_mean = mean(Tested_after_all_filters)
)
print(chromosome_summary, width = Inf)

# Save Images

ggsave("Supplemental_Figure_S2_F2_Phenotypic_Segregation.png", plot = segregation_plot, width = 7, height = 5, dpi = 600, bg = "white")
ggsave("Supplemental_Figure_S3_LSU500_Association.png", plot = marker_plot, width = 9, height = 6, dpi = 600, bg = "white")
