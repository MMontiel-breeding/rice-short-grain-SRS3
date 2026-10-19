################################################################################
# IC206 × SG F2 POPULATION
# SEGREGATION ANALYSIS
################################################################################

rm(list = ls())

library(readxl)
library(dplyr)

################################################################################
# IMPORT
################################################################################

dat <- read_excel("Data/05 F2 147 plants - segregation analysis.xlsx")

################################################################################
# PHENOTYPE COUNTS
################################################################################

obs <- table(dat$`Pheno/position`)

obs

wt <- as.numeric(obs["Wild"])
sg <- as.numeric(obs["Mutant"])

total <- wt + sg

################################################################################
# EXPECTED 3:1 RATIO
################################################################################

exp_wt <- total * 0.75
exp_sg <- total * 0.25

################################################################################
# CHI-SQUARE TEST
################################################################################

chi <- chisq.test(
  x = c(wt, sg),
  p = c(0.75, 0.25),
  correct = FALSE
)

################################################################################
# RESULTS TABLE
################################################################################

segregation_results <- data.frame(
  Total_Plants = total,
  Observed_WT = wt,
  Observed_SG = sg,
  Expected_WT = exp_wt,
  Expected_SG = exp_sg,
  Chi_Square = as.numeric(chi$statistic),
  DF = as.numeric(chi$parameter),
  P_Value = chi$p.value
)

segregation_results

################################################################################
# EXPORT
################################################################################

write.csv(
  segregation_results,
  "IC206_F2_Segregation_Analysis.csv",
  row.names = FALSE
)

################################################################################
# CONSOLE REPORT
################################################################################

cat("\n==============================\n")
cat("SEGREGATION ANALYSIS\n")
cat("==============================\n")
cat("WT:", wt, "\n")
cat("SG:", sg, "\n")
cat("Total:", total, "\n")
cat("Expected WT:", round(exp_wt,2), "\n")
cat("Expected SG:", round(exp_sg,2), "\n")
cat("Chi-square:", round(as.numeric(chi$statistic),4), "\n")
cat("P-value:", signif(chi$p.value,4), "\n")