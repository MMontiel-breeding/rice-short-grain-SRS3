################################################################################
# AVANT-DERIVED F2 POPULATION
# CHI-SQUARE SEGREGATION ANALYSIS
#
# Author: Maria Montiel
#
# Description:
# A subset of 100 F2 individuals from the Avant-derived population was
# selected and used for segregation analysis.
#
# The SG phenotype was evaluated against the expected 3:1 Mendelian
# segregation ratio for a single recessive locus using a chi-square
# goodness-of-fit test.
################################################################################

rm(list = ls())

################################################################################
# IMPORT DATA
################################################################################

dat <- read.csv(
  "Table_S7_Avant_F2_Phenotypes.csv",
  stringsAsFactors = FALSE
)

################################################################################
# PHENOTYPE COUNTS
################################################################################

counts <- table(dat$phenotype)

WT <- as.numeric(counts["Wild"])
SG <- as.numeric(counts["Mutant"])

total <- WT + SG

################################################################################
# EXPECTED COUNTS
################################################################################

expected_WT <- total * 0.75
expected_SG <- total * 0.25

################################################################################
# CHI-SQUARE TEST
################################################################################

chi <- chisq.test(
  x = c(WT, SG),
  p = c(0.75, 0.25),
  correct = FALSE
)

################################################################################
# RESULTS TABLE
################################################################################

results <- data.frame(
  Population = "Avant F2",
  WT = WT,
  SG = SG,
  Total = total,
  Expected_WT = expected_WT,
  Expected_SG = expected_SG,
  Chi_Square = as.numeric(chi$statistic),
  DF = as.numeric(chi$parameter),
  P_Value = chi$p.value
)

print(results)

