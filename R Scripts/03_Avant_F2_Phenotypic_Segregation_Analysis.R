################################################################################
# AVANT-DERIVED F2 POPULATION
# CHI-SQUARE SEGREGATION ANALYSIS
#
# Author: Montiel et al.
# Script 03 | Avant F2  Phenotypic Segregation
#
# Imput data: S6 Table. Grain phenotypes and segregation analysis of 100 F₂ 
# plants from the Avant × RU2002174 SG cross.
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
  "S6_Table",
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

