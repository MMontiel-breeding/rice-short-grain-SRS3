
#################################################################################
# Paper: Genetic mapping of a putative spontaneous short-grain mutation
#        identifies a predicted loss-of-function allele of SRS3 in rice
#
# Script author: Maria G. Montiel
# Analysis: Puerto Rico WT–SG marker-genotype comparison
#
# Purpose: compare marker-genotype backgrounds of WT and SG plants selected
# from a single segregating panicle row.
#
# Outputs: results and filtering description displayed in the console.
# No files are saved automatically.
################################################################################

rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# -> Choose whether to exclude the two samples flagged for lower completeness.
exclude_low_completeness <- FALSE

flagged_samples <- c("SM-Long-41", "SM-Short-45")
call_rate_threshold <- 0.95

input_file <- "Table_S4_Puerto_Rico_98_Plants_Phenotypes_and_LSU500_Genotypes.xlsx"
marker_file <- "Table_S3_LSU500_SNP_Marker_Information.xlsx"

if (!all(file.exists(c(input_file, marker_file)))) {
  stop("Input files not found. Check the working directory and input paths.")
}

#--------------------------------------------------
# Read data
#--------------------------------------------------

# -> Validate sample identifiers and the original phenotype labels.
dat <- read_excel(input_file)
meta_cols <- c("Plate name", "Well", "Sample_ID", "Class")

if (!all(meta_cols %in% names(dat))) {
  stop("Required sample metadata columns are missing.")
}

dat <- dat %>%
  mutate(
    Sample_ID = str_trim(as.character(Sample_ID)),
    Class = str_trim(as.character(Class))
  )

if (any(is.na(dat$Sample_ID) | dat$Sample_ID == "")) {
  stop("Missing sample IDs.")
}

if (anyDuplicated(dat$Sample_ID)) {
  stop("Duplicate sample IDs.")
}

if (any(is.na(dat$Class) | !dat$Class %in% c("Long", "Short"))) {
  stop("Unrecognized phenotype class. Expected Long or Short.")
}

snp_cols <- setdiff(names(dat), meta_cols)
n_samples_original <- nrow(dat)
n_markers_original <- length(snp_cols)

if (n_samples_original != 98L || n_markers_original != 550L) {
  warning("Input differs from the expected 98 samples and 550 markers.")
}

cat("\nINPUT DATA\n")
print(table(dat$Class))
cat("Samples:", n_samples_original, "| Markers:", n_markers_original, "\n")

#--------------------------------------------------
# Marker information
#--------------------------------------------------

# -> Match markers by ID and obtain their documented allele definitions.
agriplex <- read_excel(
  marker_file,
  sheet = "final_agriplex_snps_list"
)

agriplex <- agriplex %>%
  rename(
    Marker_ID = Agriplex_ID_LSU500,
    IRGSPv_Position = `IRGSPv_Pos (bp)`
  )

required_cols <- c(
  "Marker_ID", "Chr", "IRGSPv_Position", "Allele1", "Allele2"
)

if (!all(required_cols %in% names(agriplex))) {
  stop("Required marker dictionary columns are missing.")
}

agriplex <- agriplex %>%
  mutate(
    Marker_ID = str_trim(as.character(Marker_ID)),
    Allele1 = toupper(str_trim(as.character(Allele1))),
    Allele2 = toupper(str_trim(as.character(Allele2)))
  )

if (anyNA(agriplex$Marker_ID) || anyDuplicated(agriplex$Marker_ID)) {
  stop("Missing or duplicate marker IDs in the dictionary.")
}

marker_info <- tibble(
  SNP = snp_cols,
  Marker_ID = str_remove(snp_cols, "^LSU_SNP[0-9]+_")
) %>%
  left_join(agriplex, by = "Marker_ID")

if (any(is.na(marker_info$Allele1) | marker_info$Allele1 == "") ||
    any(is.na(marker_info$Allele2) | marker_info$Allele2 == "")) {
  stop("Some markers lack documented allele definitions.")
}

# Source coordinates are preserved, not converted.
# Decimal coordinate entries require verification before reporting as bp.

#--------------------------------------------------
# Standardize genotype calls
#--------------------------------------------------

# -> Treat reversed heterozygotes as equivalent and preserve multi-base alleles.
# Examples: A / G = G / A; GC = GC/GC; GC / TT = GC/TT.
missing_labels <- c(
  "", "FAIL", "FAILED", "NA", "N/A",
  "NOCALL", "NO_CALL", ".", "-", "--"
)

normalize_gt <- function(x, allele1, allele2) {
  x <- toupper(str_remove_all(str_trim(as.character(x)), "\\s+"))
  x[is.na(x) | x %in% missing_labels] <- NA_character_
  
  vapply(x, function(g) {
    if (is.na(g)) return(NA_character_)
    
    alleles <- strsplit(g, "/", fixed = TRUE)[[1]]
    
    if (length(alleles) == 1L &&
        alleles %in% c(allele1, allele2)) {
      alleles <- rep(alleles, 2L)
    }
    
    if (length(alleles) != 2L ||
        any(!alleles %in% c(allele1, allele2))) {
      return(NA_character_)
    }
    
    paste(sort(alleles), collapse = "/")
  }, character(1), USE.NAMES = FALSE)
}

invalid_calls <- tibble(
  Sample_ID = character(),
  SNP = character(),
  Original_call = character()
)

for (j in seq_along(snp_cols)) {
  s <- snp_cols[j]
  original <- as.character(dat[[s]])
  
  clean_original <- toupper(
    str_remove_all(str_trim(original), "\\s+")
  )
  
  dat[[s]] <- normalize_gt(
    original,
    marker_info$Allele1[j],
    marker_info$Allele2[j]
  )
  
  originally_missing <- is.na(clean_original) |
    clean_original %in% missing_labels
  
  bad <- which(!originally_missing & is.na(dat[[s]]))
  
  if (length(bad)) {
    invalid_calls <- bind_rows(
      invalid_calls,
      tibble(
        Sample_ID = dat$Sample_ID[bad],
        SNP = s,
        Original_call = original[bad]
      )
    )
  }
}

if (nrow(invalid_calls)) {
  print(invalid_calls, n = Inf)
  stop("Unexpected allele calls detected. Resolve them before continuing.")
}

#--------------------------------------------------
# Sample completeness
#--------------------------------------------------

# -> Describe completeness before optionally excluding the two flagged samples.
sample_qc <- dat %>%
  select(Sample_ID, Class) %>%
  mutate(
    N_Called = rowSums(!is.na(as.matrix(dat[snp_cols]))),
    Call_Rate_All_Markers = N_Called / n_markers_original,
    Flagged = Sample_ID %in% flagged_samples,
    Excluded = exclude_low_completeness & Flagged
  )

cat("\nFLAGGED SAMPLE COMPLETENESS BEFORE MARKER QC\n")
print(filter(sample_qc, Flagged), width = Inf)

if (exclude_low_completeness) {
  if (!all(flagged_samples %in% dat$Sample_ID)) {
    warning("One or more flagged sample IDs were not present.")
  }
  
  dat <- dat %>%
    filter(!Sample_ID %in% flagged_samples)
}

if (!all(c("Long", "Short") %in% dat$Class)) {
  stop("Both phenotype classes are required.")
}

n_excluded <- n_samples_original - nrow(dat)

################################################################################
# QUALITY CONTROL
################################################################################

# -> Retain markers with >=95% valid calls among the selected samples.
# Failed and missing calls remain in the denominator.
qc_table <- tibble(
  SNP = snp_cols,
  N_Total = nrow(dat),
  N_Called = vapply(
    dat[snp_cols],
    function(x) sum(!is.na(x)),
    integer(1)
  )
) %>%
  mutate(
    N_Missing = N_Total - N_Called,
    CallRate = N_Called / N_Total,
    Pass_Call_Rate = CallRate >= call_rate_threshold
  ) %>%
  left_join(marker_info, by = "SNP")

snp_cols_qc <- qc_table %>%
  filter(Pass_Call_Rate) %>%
  pull(SNP)

if (!length(snp_cols_qc)) {
  stop("No markers passed call-rate QC.")
}

snp_cols <- snp_cols_qc

cat("\nMARKER FILTERING\n")
cat("Samples analyzed:", nrow(dat), "\n")
cat("Original markers:", n_markers_original, "\n")
cat("Retained markers:", length(snp_cols), "\n")
cat("Removed markers:", n_markers_original - length(snp_cols), "\n")
cat("No MAF or polymorphism filter was applied.\n")

# -> Reassess sample completeness using only the retained markers.
sample_qc_after <- dat %>%
  select(Sample_ID, Class) %>%
  mutate(
    N_Called = rowSums(!is.na(as.matrix(dat[snp_cols]))),
    Call_Rate_QC_Markers = N_Called / length(snp_cols)
  )

cat("\nLOWEST SAMPLE CALL RATES AFTER MARKER QC\n")
print(
  sample_qc_after %>% arrange(Call_Rate_QC_Markers) %>% slice_head(n = 5),
  width = Inf
)

################################################################################
# POLYMORPHIC MARKERS AFTER QC
################################################################################

# -> Count markers showing genotype variation; do not restrict the comparison.
polymorphic_snps <- vapply(
  dat[snp_cols],
  function(x) length(unique(x[!is.na(x)])) > 1L,
  logical(1)
)

n_poly <- sum(polymorphic_snps)

poly_chr <- tibble(SNP = snp_cols[polymorphic_snps]) %>%
  left_join(marker_info, by = "SNP")

poly_chr_summary <- poly_chr %>%
  count(Chr, name = "Variable_Markers") %>%
  arrange(Chr)

cat("\nVARIABLE MARKERS AFTER QC\n")
cat("Total:", n_poly, "\n")
print(poly_chr_summary, n = Inf)

################################################################################
# CONSENSUS COMPARISON
################################################################################

# -> Compare the unique modal genotype in WT and SG at each retained marker.
# A mode is the most frequent genotype; it need not exceed 50%.
consensus_gt <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  
  tab <- table(x)
  if (sum(tab == max(tab)) > 1L) return(NA_character_)
  
  names(tab)[which.max(tab)]
}

modal_status <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("No valid calls")
  
  tab <- table(x)
  if (sum(tab == max(tab)) > 1L) return("Tied mode")
  
  "Evaluable"
}

long_df <- dat %>% filter(Class == "Long")
short_df <- dat %>% filter(Class == "Short")

consensus_results <- lapply(snp_cols, function(s) {
  lg <- long_df[[s]]
  sg <- short_df[[s]]
  
  long_cons <- consensus_gt(lg)
  short_cons <- consensus_gt(sg)
  
  different <- if (is.na(long_cons) || is.na(short_cons)) {
    NA
  } else {
    long_cons != short_cons
  }
  
  tibble(
    SNP = s,
    Long_consensus = long_cons,
    Short_consensus = short_cons,
    Long_N_Called = sum(!is.na(lg)),
    Short_N_Called = sum(!is.na(sg)),
    Long_Status = modal_status(lg),
    Short_Status = modal_status(sg),
    Different = different
  )
}) %>%
  bind_rows() %>%
  left_join(marker_info, by = "SNP")

# -> Exclude tied or unavailable modes from the concordance denominator.
total_scored <- sum(!is.na(consensus_results$Different))

if (!total_scored) {
  stop("No evaluable modal-genotype comparisons.")
}

different_snps <- sum(consensus_results$Different, na.rm = TRUE)
matching_snps <- total_scored - different_snps
percent_concordance <- 100 * matching_snps / total_scored

cat("\nMODAL-GENOTYPE COMPARISON\n")
cat("Evaluable markers:", total_scored, "\n")
cat("Matching modes:", matching_snps, "\n")
cat("Different modes:", different_snps, "\n")
cat("Concordance:", round(percent_concordance, 2), "%\n")

cat("\nMARKERS WITH DIFFERENT MODAL GENOTYPES\n")
print(
  consensus_results %>%
    filter(Different %in% TRUE) %>%
    select(Marker_ID, Chr, Long_consensus, Short_consensus),
  n = Inf
)

cat("\nMARKERS EXCLUDED FROM MODAL COMPARISON\n")
print(
  consensus_results %>%
    filter(is.na(Different)) %>%
    select(Marker_ID, Long_Status, Short_Status),
  n = Inf
)

################################################################################
# GENOTYPE FREQUENCIES
################################################################################

# -> Summarize each genotype's frequency among valid calls within each class.
frequency_results <- lapply(snp_cols, function(s) {
  dat %>%
    transmute(Class, Genotype = .data[[s]]) %>%
    filter(!is.na(Genotype)) %>%
    count(Class, Genotype, name = "N") %>%
    group_by(Class) %>%
    mutate(N_Called = sum(N), Frequency = N / N_Called) %>%
    ungroup() %>%
    mutate(SNP = s)
}) %>%
  bind_rows()

################################################################################
# MOST DIVERGENT MARKERS
################################################################################

# -> Compare frequencies of the same genotype between WT and SG.
# Total variation distance ranges from 0 (identical) to 1 (no overlap).
# This is an additional descriptive summary, not a significance test.
wide_freq <- frequency_results %>%
  select(SNP, Class, Genotype, Frequency) %>%
  pivot_wider(
    names_from = Class,
    values_from = Frequency,
    values_fill = 0
  ) %>%
  left_join(
    select(consensus_results, SNP, Long_N_Called, Short_N_Called),
    by = "SNP"
  ) %>%
  mutate(
    Frequency_Difference = ifelse(
      Long_N_Called > 0 & Short_N_Called > 0,
      abs(Long - Short),
      NA_real_
    )
  )

divergence_summary <- wide_freq %>%
  group_by(SNP) %>%
  summarise(
    Total_Variation_Distance = if (all(is.na(Frequency_Difference))) {
      NA_real_
    } else {
      sum(Frequency_Difference) / 2
    },
    .groups = "drop"
  ) %>%
  left_join(marker_info, by = "SNP") %>%
  arrange(desc(Total_Variation_Distance))

cat("\nTEN MARKERS WITH THE LARGEST GENOTYPE-FREQUENCY DIFFERENCES\n")
print(
  divergence_summary %>%
    select(Marker_ID, Chr, Total_Variation_Distance) %>%
    slice_head(n = 10),
  width = Inf
)

################################################################################
# WITHIN-CLASS VARIABILITY
################################################################################

# -> Distinguish observed uniformity, variation, and absence of valid calls.
within_class_variation <- function(df, snp_cols) {
  lapply(snp_cols, function(s) {
    called <- df[[s]][!is.na(df[[s]])]
    n_genotypes <- length(unique(called))
    
    status <- if (n_genotypes == 0L) {
      "No valid calls"
    } else if (n_genotypes == 1L) {
      "One observed genotype"
    } else {
      "Variable"
    }
    
    tibble(
      SNP = s,
      N_Called = length(called),
      N_unique_genotypes = n_genotypes,
      Status = status
    )
  }) %>%
    bind_rows()
}

long_var <- within_class_variation(long_df, snp_cols)
short_var <- within_class_variation(short_df, snp_cols)

cat("\nWT WITHIN-CLASS VARIABILITY\n")
print(count(long_var, Status))

cat("\nSG WITHIN-CLASS VARIABILITY\n")
print(count(short_var, Status))

################################################################################
# FINAL SUMMARY AND RESULTS DESCRIPTION
################################################################################

# -> Assemble the counts needed to reconcile Methods and Results.
analysis_summary <- tibble(
  Samples_Analyzed = nrow(dat),
  WT_Samples = nrow(long_df),
  SG_Samples = nrow(short_df),
  Samples_Excluded = n_excluded,
  Original_Markers = n_markers_original,
  QC_Retained_Markers = length(snp_cols),
  Variable_Markers = n_poly,
  Evaluable_Markers = total_scored,
  Matching_Modes = matching_snps,
  Different_Modes = different_snps,
  Unevaluable_Markers = length(snp_cols) - total_scored,
  Concordance_Percent = percent_concordance
)

cat("\nFINAL SUMMARY\n")
print(analysis_summary, width = Inf)

filtering_description <- sprintf(
  paste0(
    "Genotype calls were standardized using documented allele definitions. ",
    "Failed and missing calls were counted against marker call rate. ",
    "Of %d markers, %d passed the >=%.0f%% call-rate threshold across ",
    "%d analyzed samples; %d markers were removed. ",
    "No MAF or polymorphism filter was applied."
  ),
  n_markers_original, length(snp_cols), 100 * call_rate_threshold,
  nrow(dat), n_markers_original - length(snp_cols)
)

results_description <- sprintf(
  paste0(
    "The analysis included %d WT and %d SG plants. ",
    "Among the retained markers, %d showed genotype variation. ",
    "Modal genotypes were evaluable in both classes at %d markers: ",
    "%d matched and %d differed, giving %.2f%% modal-genotype concordance. ",
    "%d markers were excluded from this comparison because a mode was ",
    "tied or unavailable in at least one class."
  ),
  nrow(long_df), nrow(short_df), n_poly, total_scored,
  matching_snps, different_snps, percent_concordance,
  length(snp_cols) - total_scored
)

cat("\nFILTERING DESCRIPTION\n", filtering_description, "\n", sep = "")
cat("\nRESULTS DESCRIPTION\n", results_description, "\n", sep = "")

cat(
  "\nSamples excluded:", n_excluded,
  "\nConcordance is descriptive class-level marker agreement, ",
  "not individual-level or whole-genome identity.",
  "\nSample counts do not represent a segregation ratio.",
  "\nMarker QC is recalculated when flagged samples are excluded.",
  "\nNo files were saved.\n",
  sep = ""
)

# Optional inspection of complete results:
# View(qc_table)
# View(consensus_results)
# View(frequency_results)
# sessionInfo()