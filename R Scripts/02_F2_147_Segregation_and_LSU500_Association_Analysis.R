################################################################################
# Paper: Genetic mapping of a putative spontaneous short-grain mutation identifies
# a predicted loss-of-function allele of SRS3 in rice
#
# Author: Montiel et al.
# Script 02 | F2 LSU500 single-marker association analysis
#
# ANALYSES
# 1. Pearson chi-square goodness-of-fit test of the expected 3 WT:1 SG ratio.
# 2. Marker QC using call rate, genotype variability and MAF > 0.20.
# 3. Preliminary categorical ANOVA of the binary WT/SG phenotype.
# 4. Primary categorical binomial logistic likelihood-ratio test.
# 5. Benjamini-Hochberg FDR and Bonferroni multiple-testing correction.
# 6. Descriptive WT/SG distributions within genotype classes as effect estimates.
#
#################################################################################

rm(list = ls())

# ── 1. PREPARE THE ANALYSIS ───────────────────────────────────────────────────

required_packages <- c("readxl", "dplyr", "stringr", "ggplot2", "tibble", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Install: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr)
  library(ggplot2); library(tibble); library(tidyr)
})

# ── 2. SELECT INPUTS AND FILTERS ──────────────────────────────────────────────

input_file <- file.choose()
marker_file <- file.choose()

sheet_name <- 1

minimum_call_rate <- 0.8
minimum_maf <- 0.20
fdr_alpha <- 0.01
expected_wt_proportion <- 0.75
expected_sg_proportion <- 0.25

if (!all(file.exists(c(input_file, marker_file)))) stop("Check the input file paths.")

# ── 3. DEFINE VALIDATION AND ANALYSIS FUNCTIONS ───────────────────────────────

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
    paste(sort(z), collapse = "/")
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
    return(tibble(
      anova_n = nrow(d),
      anova_df = NA_real_,
      anova_F = NA_real_,
      anova_p = NA_real_
    ))
  }
  
  fit <- lm(y ~ genotype, data = d)
  a <- anova(fit)
  
  tibble(
    anova_n = nrow(d),
    anova_df = a$Df[1],
    anova_F = a$`F value`[1],
    anova_p = a$`Pr(>F)`[1]
  )
}

safe_logistic <- function(pheno_binary, genotype) {
  d <- tibble(y = pheno_binary, genotype = factor(genotype)) %>%
    filter(!is.na(y), !is.na(genotype)) %>% droplevels()
  
  if (nrow(d) < 10 || nlevels(d$genotype) < 2 || length(unique(d$y)) < 2) {
    return(tibble(
      logistic_n = nrow(d),
      logistic_df = NA_real_,
      logistic_LR_chisq = NA_real_,
      logistic_p = NA_real_,
      boundary_class = NA
    ))
  }
  
  counts <- d %>%
    group_by(genotype) %>%
    summarise(n = n(), sg = sum(y), .groups = "drop")
  
  binomial_loglik <- function(sg, n, p) {
    wt <- n - sg
    value <- numeric(length(n))
    
    i <- sg > 0
    value[i] <- value[i] + sg[i] * log(p[i])
    
    i <- wt > 0
    value[i] <- value[i] + wt[i] * log1p(-p[i])
    
    sum(value)
  }
  
  ll_full <- binomial_loglik(counts$sg, counts$n, counts$sg / counts$n)
  ll_null <- binomial_loglik(sum(d$y), nrow(d), mean(d$y))
  lr <- max(0, 2 * (ll_full - ll_null))
  df_test <- nrow(counts) - 1
  
  tibble(
    logistic_n = nrow(d),
    logistic_df = df_test,
    logistic_LR_chisq = lr,
    logistic_p = pchisq(lr, df = df_test, lower.tail = FALSE),
    boundary_class = any(counts$sg == 0 | counts$sg == counts$n)
  )
}

# ── 4. IMPORT AND VALIDATE THE DATA ───────────────────────────────────────────

raw_data <- read_excel(input_file, sheet = sheet_name, .name_repair = "minimal")

if (anyDuplicated(names(raw_data))) stop("Duplicate column names in the genotype workbook.")

if (ncol(raw_data) < 3 || names(raw_data)[1] != "Name" ||
    !str_detect(names(raw_data)[2], regex("phenotype", ignore_case = TRUE))) {
  stop("Expected one header row: Name, Phenotype, then marker IDs.")
}

marker_ids <- names(raw_data)[-(1:2)]
if (any(is.na(marker_ids) | marker_ids == "")) stop("Missing marker IDs.")

analysis_data <- raw_data %>%
  mutate(
    Sample_ID = clean_text(.data[[names(raw_data)[1]]]),
    Phenotype_original = .data[[names(raw_data)[2]]],
    Phenotype = canonical_phenotype(Phenotype_original),
    Phenotype_binary = case_when(
      Phenotype == "SG" ~ 1,
      Phenotype == "WT" ~ 0,
      TRUE ~ NA_real_
    )
  )

if (anyNA(analysis_data$Sample_ID) || anyDuplicated(analysis_data$Sample_ID)) {
  stop("Missing or duplicate sample IDs.")
}

if (anyNA(analysis_data$Phenotype)) stop("Unrecognized phenotype labels.")

analysis_f2 <- analysis_data

if (nrow(analysis_f2) != 147) {
  warning("Expected 147 plants; observed ", nrow(analysis_f2), ".")
}

marker_reference <- read_excel(marker_file, sheet = 1) %>%
  dplyr::rename(
    Marker_ID = Agriplex_ID_LSU500,
    IRGSPv_Position = `IRGSPv_Pos (bp)`
  )

required_fields <- c("Marker_ID", "Chr", "IRGSPv_Position", "Allele1", "Allele2")
if (!all(required_fields %in% names(marker_reference))) {
  stop("Missing marker-dictionary fields.")
}

marker_reference <- marker_reference %>%
  transmute(
    marker_id = clean_text(Marker_ID),
    chromosome = as.integer(Chr),
    position_raw = as.character(IRGSPv_Position),
    allele1 = toupper(clean_text(Allele1)),
    allele2 = toupper(clean_text(Allele2))
  )

if (anyNA(marker_reference$marker_id) || anyDuplicated(marker_reference$marker_id)) {
  stop("Missing or duplicate dictionary IDs.")
}

marker_dictionary <- tibble(
  marker_column = marker_ids,
  marker_id = marker_ids
) %>%
  left_join(marker_reference, by = "marker_id")

if (anyNA(marker_dictionary$allele1) || anyNA(marker_dictionary$allele2)) {
  stop("Some input markers lack documented alleles.")
}

marker_dictionary <- marker_dictionary %>%
  mutate(
    position_value = suppressWarnings(as.numeric(position_raw)),
    position_bp = ifelse(
      !is.na(position_value) &
        position_value > 100 &
        position_value == floor(position_value),
      position_value,
      NA_real_
    )
  ) %>%
  select(
    marker_column,
    marker_id,
    chromosome,
    position_raw,
    position_bp,
    allele1,
    allele2
  )

invalid_calls <- list()

for (i in seq_len(nrow(marker_dictionary))) {
  m <- marker_dictionary$marker_column[i]
  original <- analysis_f2[[m]]
  
  calls <- normalize_genotype(
    original,
    marker_dictionary$allele1[i],
    marker_dictionary$allele2[i]
  )
  
  bad <- which(!is.na(calls) & calls == "INVALID")
  
  if (length(bad)) {
    invalid_calls[[m]] <- tibble(
      Sample_ID = analysis_f2$Sample_ID[bad],
      marker_id = m,
      original_call = as.character(original[bad])
    )
  }
  
  analysis_f2[[m]] <- calls
}

invalid_calls <- bind_rows(invalid_calls)

if (nrow(invalid_calls)) {
  print(invalid_calls, n = 20)
  stop("Undocumented genotype calls; inspect before proceeding.")
}

cat("\n--- PHENOTYPE VALIDATION ---\n")
print(table(analysis_f2$Phenotype_original))
print(table(analysis_f2$Phenotype))

# ── 5. CHECK THE EXPECTED 3:1 SEGREGATION ─────────────────────────────────────

phenotype_counts <- tibble(
  Phenotype = c("WT", "SG"),
  Observed = c(
    sum(analysis_f2$Phenotype == "WT"),
    sum(analysis_f2$Phenotype == "SG")
  ),
  Expected_proportion = c(
    expected_wt_proportion,
    expected_sg_proportion
  )
) %>%
  mutate(Expected = sum(Observed) * Expected_proportion)

chi_test <- chisq.test(
  phenotype_counts$Observed,
  p = phenotype_counts$Expected_proportion
)

chi_results <- tibble(
  total_F2 = nrow(analysis_f2),
  observed_WT = phenotype_counts$Observed[1],
  observed_SG = phenotype_counts$Observed[2],
  chi_square = unname(chi_test$statistic),
  df = unname(chi_test$parameter),
  p_value = chi_test$p.value
)

cat("\n--- SEGREGATION TEST ---\n")
print(phenotype_counts)
print(chi_results, width = Inf)

plot_data <- phenotype_counts %>%
  select(Phenotype, Observed, Expected) %>%
  pivot_longer(
    c(Observed, Expected),
    names_to = "Count_type",
    values_to = "Count"
  )

segregation_plot <- ggplot(
  plot_data,
  aes(Phenotype, Count, fill = Count_type)
) +
  geom_col(position = "dodge", width = 0.7) +
  theme_classic() +
  labs(
    title = "F2 phenotype segregation",
    subtitle = sprintf(
      "Expected ratio: 3 WT : 1 SG\nχ² = %.2f; df = %d; P = %.3f",
      chi_results$chi_square,
      chi_results$df,
      chi_results$p_value
    ),
    x = "Phenotype",
    y = "Number of plants",
    fill = NULL
  )

print(segregation_plot)

# ── 6. FILTER MARKERS AND SUMMARIZE COMPLETENESS ──────────────────────────────

marker_qc <- lapply(
  seq_len(nrow(marker_dictionary)),
  function(i) {
    m <- marker_dictionary$marker_column[i]
    g <- analysis_f2[[m]]
    called <- !is.na(g)
    
    tibble(
      marker_column = m,
      n_total = length(g),
      n_called = sum(called),
      call_rate = mean(called),
      n_genotype_classes = length(unique(g[called])),
      maf = maf_from_calls(
        g,
        marker_dictionary$allele1[i],
        marker_dictionary$allele2[i]
      )
    )
  }
) %>%
  bind_rows() %>%
  left_join(marker_dictionary, by = "marker_column") %>%
  mutate(
    passes_call_rate = call_rate >= minimum_call_rate,
    passes_maf = !is.na(maf) & maf > minimum_maf,
    passes_all_filters =
      passes_call_rate &
      passes_maf &
      n_genotype_classes >= 2,
    filter_reason = case_when(
      !passes_call_rate ~ "Low call rate",
      n_genotype_classes < 2 ~ "Monomorphic/no valid calls",
      !passes_maf ~ "MAF filter",
      TRUE ~ "Tested"
    )
  )

cat("\n--- MARKER FILTERING ---\n")
print(marker_qc %>% count(filter_reason), n = Inf)

markers_for_analysis <- marker_qc %>%
  filter(passes_all_filters) %>%
  pull(marker_column)

if (!length(markers_for_analysis)) {
  stop("No markers passed the selected filters.")
}

sample_qc <- tibble(
  Sample_ID = analysis_f2$Sample_ID,
  Phenotype = analysis_f2$Phenotype,
  completeness_QC_markers = rowMeans(
    !is.na(
      as.data.frame(
        analysis_f2[
          marker_qc$marker_column[
            marker_qc$passes_call_rate
          ]
        ]
      )
    )
  )
)

cat("\nLowest sample completeness on QC-retained markers; no plants removed:\n")
print(
  sample_qc %>%
    arrange(completeness_QC_markers) %>%
    slice_head(n = 5)
)

# ── 7. TEST GENOTYPE-PHENOTYPE ASSOCIATION ────────────────────────────────────

all_results <- lapply(
  markers_for_analysis,
  function(m) {
    marker_dictionary %>%
      filter(marker_column == m) %>%
      bind_cols(
        safe_anova(
          analysis_f2$Phenotype_binary,
          analysis_f2[[m]]
        ),
        safe_logistic(
          analysis_f2$Phenotype_binary,
          analysis_f2[[m]]
        )
      )
  }
) %>%
  bind_rows()

if (anyNA(all_results$anova_p) || anyNA(all_results$logistic_p)) {
  stop(
    "A marker test failed; inspect all_results before multiple-testing correction."
  )
}

all_results <- all_results %>%
  mutate(
    anova_FDR = p.adjust(anova_p, "BH"),
    anova_Bonferroni = p.adjust(anova_p, "bonferroni"),
    logistic_FDR = p.adjust(logistic_p, "BH"),
    logistic_Bonferroni = p.adjust(logistic_p, "bonferroni"),
    logistic_neglog10P = -log10(
      pmax(logistic_p, .Machine$double.xmin)
    )
  ) %>%
  left_join(
    marker_qc %>%
      select(
        marker_column,
        call_rate,
        maf,
        n_called,
        n_genotype_classes
      ),
    by = "marker_column"
  )


anova_ranked <- all_results %>%
  arrange(anova_p)

logistic_ranked <- all_results %>%
  arrange(logistic_p)

significant_markers <- logistic_ranked %>%
  filter(
    logistic_FDR < fdr_alpha &
      logistic_Bonferroni < fdr_alpha
  )

cat("\n--- TOP MARKERS: ANOVA ---\n")
print(
  anova_ranked %>%
    select(
      marker_id,
      chromosome,
      position_bp,
      anova_F,
      anova_df,
      anova_p,
      anova_FDR,
      anova_Bonferroni
    ) %>%
    slice_head(n = 10),
  width = Inf
)

cat("\n--- TOP MARKERS: LOGISTIC LRT ---\n")
print(
  logistic_ranked %>%
    select(
      marker_id,
      chromosome,
      position_bp,
      logistic_LR_chisq,
      logistic_df,
      logistic_p,
      logistic_FDR,
      logistic_Bonferroni,
      boundary_class
    ) %>%
    slice_head(n = 10),
  width = Inf
)

cat("\n--- SIGNIFICANT SNPs: MARKER, CHROMOSOME AND POSITION ---\n")

if (nrow(significant_markers)) {
  print(
    significant_markers %>%
      select(
        marker_id,
        chromosome,
        position_bp,
        logistic_LR_chisq,
        logistic_df,
        logistic_p,
        logistic_FDR,
        logistic_Bonferroni
      ),
    n = Inf,
    width = Inf
  )
} else {
  cat(
    "No markers were significant under both logistic corrections at alpha = ",
    fdr_alpha,
    ".\n",
    sep = ""
  )
}

# ── 8. DISPLAY THE GENOME-WIDE ASSOCIATION SCAN ───────────────────────────────

plot_markers <- logistic_ranked %>%
  filter(
    !is.na(chromosome),
    !is.na(position_bp)
  ) %>%
  arrange(
    chromosome,
    position_bp
  )

bonferroni_threshold <- fdr_alpha / nrow(all_results)

if (nrow(plot_markers)) {
  chr_info <- plot_markers %>%
    group_by(chromosome) %>%
    summarise(
      chr_length = max(position_bp),
      .groups = "drop"
    ) %>%
    arrange(chromosome) %>%
    mutate(
      offset = lag(
        cumsum(chr_length),
        default = 0
      )
    )
  
  plot_markers <- plot_markers %>%
    left_join(chr_info, by = "chromosome") %>%
    mutate(
      cumulative_position = position_bp + offset
    )
  
  axis_positions <- plot_markers %>%
    group_by(chromosome) %>%
    summarise(
      center = (
        min(cumulative_position) +
          max(cumulative_position)
      ) / 2,
      .groups = "drop"
    )
  
  marker_plot <- ggplot(
    plot_markers,
    aes(
      cumulative_position,
      logistic_neglog10P,
      color = factor(chromosome %% 2)
    )
  ) +
    geom_point(size = 1.6) +
    geom_hline(
      yintercept = -log10(bonferroni_threshold),
      linetype = "dashed",
      color = "red"
    ) +
    scale_color_manual(
      values = c("#4472C4", "#70AD47"),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = axis_positions$center,
      labels = axis_positions$chromosome
    ) +
    theme_classic() +
    labs(
      title = "F2 single-marker logistic association",
      subtitle = paste(
        "Bonferroni threshold based on",
        nrow(all_results),
        "tested markers"
      ),
      x = "Chromosome",
      y = expression(-logitalic(P))
    )
  
  print(marker_plot)
}

# ── 9. SUMMARIZE FILTERING AND RESULTS ────────────────────────────────────────

analysis_summary <- tibble(
  Samples = nrow(analysis_f2),
  WT = sum(analysis_f2$Phenotype == "WT"),
  SG = sum(analysis_f2$Phenotype == "SG"),
  Original_markers = nrow(marker_qc),
  Call_rate_retained = sum(marker_qc$passes_call_rate),
  Variable_after_call_QC = sum(
    marker_qc$passes_call_rate &
      marker_qc$n_genotype_classes >= 2
  ),
  Tested_markers = nrow(all_results),
  ANOVA_FDR_significant = sum(
    all_results$anova_FDR < fdr_alpha
  ),
  ANOVA_Bonferroni_significant = sum(
    all_results$anova_Bonferroni < fdr_alpha
  ),
  Logistic_FDR_significant = sum(
    all_results$logistic_FDR < fdr_alpha
  ),
  Logistic_Bonferroni_significant = sum(
    all_results$logistic_Bonferroni < fdr_alpha
  )
)

cat("\n--- ANALYSIS SUMMARY ---\n")
print(analysis_summary, width = Inf)

cat(sprintf(
  "\nFiltering: %d plants and %d original markers. No markers excluded by identifier. ",
  analysis_summary$Samples,
  analysis_summary$Original_markers
))

cat(sprintf(
  "Valid calls / all plants defined call rate. %d markers passed call rate >= %.2f; %d were variable. ",
  analysis_summary$Call_rate_retained,
  minimum_call_rate,
  analysis_summary$Variable_after_call_QC
))

cat(
  paste0(
    "MAF > ",
    minimum_maf,
    " was additionally required. "
  )
)

cat(sprintf(
  "%d markers were tested. No plants were removed.\n",
  nrow(all_results)
))

cat(sprintf(
  paste0(
    "Results: %d WT and %d SG; 3:1 segregation chi-square = %.4f, ",
    "df = %d, P = %.6g. "
  ),
  chi_results$observed_WT,
  chi_results$observed_SG,
  chi_results$chi_square,
  chi_results$df,
  chi_results$p_value
))

cat(sprintf(
  paste0(
    "Logistic LRT detected %d BH-FDR-significant and %d ",
    "Bonferroni-significant markers at alpha %.2f.\n"
  ),
  analysis_summary$Logistic_FDR_significant,
  analysis_summary$Logistic_Bonferroni_significant,
  fdr_alpha
))

# ── 10. COMPILE FINAL TABLE ───────────────────────────────────────────────────

results_complete <- marker_qc %>%
  left_join(
    all_results %>%
      select(
        marker_id,
        logistic_n,
        logistic_LR_chisq,
        logistic_df,
        logistic_p,
        logistic_FDR,
        logistic_Bonferroni,
        anova_p
      ),
    by = "marker_id"
  ) %>%
    arrange(
    is.na(logistic_p),
    logistic_p,
    chromosome,
    position_bp
  ) %>%
  transmute(
    `Marker ID` = marker_id,
    Chromosome = chromosome,
    `Position (bp)` = position_bp,
    `Allele 1` = allele1,
    `Allele 2` = allele2,
    `Plants assessed (N)` = n_total,
    `Valid genotype calls (N)` = n_called,
    `Call rate` = call_rate,
    `Minor allele frequency (MAF)` = maf,
    `Observed genotype classes (N)` = n_genotype_classes,
    `Analysis status` = filter_reason,
    `N plants analyzed` = logistic_n,
    `Logistic LRT chi-square` = logistic_LR_chisq,
    `Logistic LRT df` = logistic_df,
    `Logistic LRT P` = logistic_p,
    `Logistic BH-adjusted P` = logistic_FDR,
    `Logistic Bonferroni-adjusted P` = logistic_Bonferroni,
    `Preliminary ANOVA P` = anova_p,
    `Significant under both logistic corrections` = case_when(
      filter_reason != "Tested" ~ "Not tested",
      logistic_FDR < fdr_alpha &
        logistic_Bonferroni < fdr_alpha ~ "Yes",
      TRUE ~ "No"
    )
  )

writexl::write_xlsx(
  results_complete,
  path = "F2_147_LSU500_Association_Results.xlsx"
)

View(results_complete)


# Save images
ggsave(
  "Supplemental_Figure_S2_F2_Phenotypic_Segregation.png",
  plot = segregation_plot,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

if (exists("marker_plot")) {
  ggsave(
    "Supplemental_Figure_S3_LSU500_Association.png",
    plot = marker_plot,
    width = 9,
    height = 6,
    dpi = 600,
    bg = "white"
  )
}

cat(
  "\nOne result table and two supplemental figures were saved automatically.\n"
)

################################################################################
# END OF SCRIPT
################################################################################