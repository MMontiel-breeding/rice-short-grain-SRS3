# Genetic Mapping of a Putative Spontaneous Short-Grain Mutation Identifies a Predicted Loss-of-Function Allele of SRS3 in Rice

## Repository Overview

This repository contains datasets, analysis scripts, supplementary materials, and supporting information associated with the manuscript:

**Montiel M., Angira B., Richards J.K., Famoso A.**
*Genetic mapping of a putative spontaneous short-grain mutation identifies a predicted loss-of-function allele of SRS3 in rice.*

The study describes the discovery, inheritance, fine mapping, and candidate variant identification of a short-grain (SG) mutant detected within the LSU AgCenter rice breeding program.

WT refers to the wild-type, long-grain phenotype, and SG refers to the short-grain phenotype. In recoded marker-genotype tables, WT and SG identify homozygotes matching the corresponding parental controls, while HET identifies heterozygotes. Marker genotypes and observed grain phenotypes are recorded separately.

## Repository Status

This repository is currently under development and will continue to be updated as manuscript revisions, supplementary materials, and supporting datasets are finalized.

The file names below provide the intended naming convention for the repository. Some files are being renamed or updated locally and may not yet be available under these names. Supplementary numbering remains provisional until all materials are finalized.

The initial IC206 mapping analysis uses 147 F₂ plants. The expanded IC206 dataset is being reconciled with the original records; the current combined working files contain 1,051 plants and should not yet be treated as the final manuscript analysis dataset.

Additional files, sequencing information, supplementary tables, and updated scripts will be added as they are verified.

## Repository Contents

### Data

Phenotypic measurements, genotype calls, and supporting metadata are organized by analysis. Supplementary tables containing individual observations also serve as analysis input files.

| Dataset | File name | Contents |
| --- | --- | --- |
| Individual grain measurements — S9 Table | `S9_Table_Grain_Phenotype_Measurements.csv` | Original grain length, width, and length-to-width ratio measurements, image identifiers, processing classifications, and WT/mutant labels, before script-based quality control. |
| PACE marker and primer information — S2 Table | `Table_S2_SNP_Markers_and_Primer_Sequences.xlsx` | Marker identifiers, chromosome positions, alleles, and primer sequences, including the mutation-targeted assay. |
| LSU500 marker information — S3 Table | `Table_S3_LSU500_SNP_Marker_Information.xlsx` | Dictionary of 550 LSU500 markers, including marker identifiers, chromosomes, supplied physical positions, and allele definitions. |
| Puerto Rico population — S4 Table | `Table_S4_Puerto_Rico_98_Plants_Phenotypes_and_LSU500_Genotypes.xlsx` | Original phenotype classifications and LSU500 calls for 98 individual plants, with sample, plate, and well identifiers. Long corresponds to WT and Short to SG. |
| Initial IC206-derived F₂ population — S5 Table | `Table_S5_F2_147_Plants_Phenotypes_and_LSU500_Genotypes.xlsx` | Plant identifiers, observed phenotypes, and original LSU500 genotype calls for the initial 147-plant mapping population. The final presentation format is under review. |
| Avant-derived F₂ population — S7 Table | `Table_S7_Avant_F2_Phenotypes.csv` | Identifiers, population labels, pedigrees, and phenotypes of 100 F₂ plants used for segregation analysis. Wild corresponds to WT and Mutant to SG. |
| Expanded IC206 population: original calls — working file | `IC206_F2_1051_Phenotypes_and_Genotypes_Original_Calls.xlsx` | Combined 1,051-plant working dataset, including all 147 initial plants and the 513 plants in the later genotyping summary. Original first-round and later-summary calls are retained in separate columns. |
| Expanded IC206 population: recoded calls — working file | `IC206_F2_1051_Phenotypes_and_Genotypes_WT_SG_HET.xlsx` | The same plants with interpretable calls recoded as WT, SG, or HET. Ambiguous calls are labeled UNRESOLVED and amplification failures NO_CALL. A marker key documents the conversion. |

The working expanded-population files retain both first-round and later phenotype annotations where they differ. Blank later-summary fields indicate that no record was available in that source. These files require reconciliation before final supplementary numbering and analysis.

Additional datasets will cover F₃ progeny tests, fine-mapping recombinants, and candidate-variant surveys. Their file names will be added after review.

### Scripts

Analysis scripts are provided in R. Each script addresses a defined part of the manuscript and identifies its required input files.

**Grain morphology — `01_Grain_Morphology_QC_and_Statistical_Analysis.R`**

- **Input:** `S9_Table_Grain_Phenotype_Measurements.csv`.
- Examines grain length, width, and length-to-width ratio before and after quality control, separately by processing type and genotype.
- Flags outliers using the 1.5 × IQR rule within each processing-type × genotype group. When exclusion is enabled, an observation flagged for any of the three traits is removed from the filtered dataset; the original measurements remain unchanged.
- Calculates group summaries, Welch's two-sample t-tests, mean and percentage differences, and Benjamini–Hochberg-adjusted P-values across the nine trait-by-processing comparisons.
- Generates the grain morphology summary in S1 Table and the boxplots in S1 Fig, together with quality-control records. These comparisons describe individual measured grains and do not establish independent plant-level biological replication.

**Puerto Rico genotyping — `02_Puerto_Rico_Original_Population_Genotyping_Analysis.R`**

- **Inputs:** S4 Table and the LSU500 marker dictionary in S3 Table.
- Standardizes genotype calls, treats failed calls as missing, and recognizes reversed heterozygote calls as equivalent.
- Filters markers at a minimum call rate of 0.95, summarizes marker variability, and compares the modal genotype of WT and SG groups at each evaluable marker.
- Reports evaluable markers, matching and differing modal genotypes, and concordance. This comparison describes group-level marker similarity, not pairwise whole-genome identity.
- The manuscript analysis retains all 98 plants. The optional exclusion of two lower-completeness samples must remain disabled for that analysis. Results are displayed in R; this script does not perform a phenotype-segregation χ² test.

**Initial IC206 mapping and segregation — `03_F2_147_Segregation_and_LSU500_Association_Analysis.R`**

- **Inputs:** the combined phenotype/genotype file for the 147 plants in S5 Table and the marker dictionary in S3 Table. A separate phenotype-only input is unnecessary.
- Tests observed WT and SG counts against the expected 3 WT:1 SG ratio using Pearson's χ² goodness-of-fit test.
- Evaluates all 550 markers, retaining those with call rate ≥0.95, at least two observed genotype classes, and minor allele frequency >0.20. MAF is calculated from allele counts among valid genotype calls.
- Tests the 265 retained markers using categorical logistic likelihood-ratio tests and complementary ANOVA on the binary phenotype. Genotype is treated as categorical without imposing an additive inheritance model.
- Applies BH and Bonferroni corrections separately to each test family and reports genotype-specific WT/SG counts, SG proportions, and exact 95% binomial confidence intervals.
- Produces results supporting S6 Table, the segregation plot in S2 Fig, and the logistic association scan in S3 Fig. The dashed association threshold is based on 0.05 divided by the number of tested markers.
- The MAF switch permits a separate run without MAF filtering. Automatic sensitivity analysis and export of the complete S6 workbook are being incorporated into the final script; current versions may retain results only as R objects.

**Avant inheritance — `04_Avant_F2_Phenotypic_Segregation_Analysis.R`**

- **Input:** `Table_S7_Avant_F2_Phenotypes.csv`.
- Counts WT and SG phenotypes among 100 F₂ plants and calculates expected counts under a 3:1 ratio.
- Performs Pearson's χ² goodness-of-fit test and displays observed counts, expected counts, the χ² statistic, degrees of freedom, and P-value.
- Supports the second-background inheritance analysis reported in the manuscript. This script uses base R and does not require additional packages.

**Expanded IC206 segregation — `05_IC206_F2_Complete_Phenotypic_Segregation_Analysis.R`**

- **Input:** the expanded-population phenotype file, pending reconciliation. The earlier file was provisionally named `Table_S8_IC206_F2_Complete_Phenotypes.xlsx`.
- Checks plant identifiers, duplicate records, and phenotype labels before calculating observed and expected counts.
- Performs Pearson's χ² goodness-of-fit test against 3 WT:1 SG, with a complementary exact binomial test and confidence interval for the SG proportion.
- Displays a phenotype-count plot and a statistical summary. It does not infer recombination breakpoints or analyze PACE associations.
- The current script expects 1,031 records. Its expected sample size and input-column mapping must be updated to the verified final dataset before rerunning it; the new 1,051-plant working files are not yet direct replacements.

Scripts for PACE recombinant classification, increased-resolution mapping, F₃ progeny testing, and candidate-variant surveys will be documented when their review is complete.

### Supplementary Information

Supplementary materials associated with the manuscript are provided in this repository and/or as Supporting Information submitted with the manuscript.

| Item | File name | Description |
| --- | --- | --- |
| S1 Table | `Table_S1_Grain_Morphology_Statistical_Summary.xlsx` | A single statistical table comparing WT and SG grain morphology after quality control. |
| S2–S5 Tables | See the Data section | Marker/primer information, LSU500 marker metadata, and individual Puerto Rico and initial IC206 records. |
| S6 Table | `Table_S6_F2_147_LSU500_Marker_QC_and_Association_Results.xlsx` | Quality-control information for all 550 markers and association results for retained markers, with multiple-testing adjustments, genotype-specific phenotype summaries, and no-MAF sensitivity results. |
| S7 Table | `Table_S7_Avant_F2_Phenotypes.csv` | Individual phenotype records for the 100-plant Avant-derived population. |
| S8 Table | Pending final dataset verification | Expanded IC206 population records; the final sample list and presentation remain under review. |
| S9 Table | `S9_Table_Grain_Phenotype_Measurements.csv` | Original individual grain measurements used for morphology analysis. |
| S1 Fig | `Supplemental_Figure_S1_Grain_Morphology_Boxplots.png` | Grain morphology boxplots by phenotype and processing type after quality control. |
| S2 Fig | `Supplemental_Figure_S2_F2_Phenotypic_Segregation.png` | Observed and expected WT/SG counts for the initial 147-plant IC206 population. |
| S3 Fig | `Supplemental_Figure_S3_LSU500_Association.png` | Logistic single-marker association scan with the Bonferroni threshold. |

Further materials will include recombinant classifications, progeny-test records, and candidate-variant evidence. Final table/figure numbering and citations will be reconciled with the revised manuscript.

## Whole-Genome Sequencing Information

Whole-genome sequencing datasets, accession links, and supporting variant files will be added as their availability and annotations are verified. The documentation will identify the sequenced samples, reference assembly, variant-processing methods, and files supporting candidate-variant identification.

## Reproducibility

This repository is intended to support independent reproduction of the analyses reported in the manuscript.

- Use the input files specified for each script. Script numbers organize the analyses; they do not imply that all scripts must be run in one R session.
- Run each script from the beginning. Scripts clear the R environment with `rm(list = ls())` and load their required packages where applicable.
- Update input paths or use the file-selection prompts provided by the script. Folder names and capitalization must match the local repository layout.
- Check phenotype labels: Long, Wild, and WT refer to the long-grain phenotype; Short, Mutant, and SG refer to the short-grain phenotype. Preserve the original labels in source datasets.
- Missing and ambiguous genotype calls must not be assigned to WT, SG, or HET without supporting evidence. Refer to the marker dictionary and any file-specific coding key.
- The initial 147 IC206 plants are part of the expanded population and must not be counted a second time or treated as an independent population replicate.
- Working files and analyses awaiting reconciliation are identified above. Published results should be reproduced using the final verified datasets and script settings.

Software/package versions and session information will be included with the finalized analysis materials. Some current scripts display results without automatically saving them; export instructions are documented within the corresponding scripts.

## Contact

Maria Montiel  
Rice Breeding Program  
H. Rouse Caffey Rice Research Station  

GitHub Repository:  
https://github.com/MMontiel-breeding/rice-short-grain-SRS3

---
