# Genetic Mapping of a Putative Spontaneous Short-Grain Mutation Identifies a Predicted Loss-of-Function Allele of SRS3 in Rice

## Overview

This repository provides R scripts and the input datasets needed to reproduce the grain morphology, F₂ phenotypic segregation, and LSU500 single-marker association analyses associated with the manuscript by Montiel M., Angira B., Richards J.K., and Famoso A.

Other supporting tables and figures are provided with the manuscript. Tables used as script inputs are also included here so the analyses can be run using the repository files.

WT denotes the wild-type, long-grain phenotype; SG denotes the short-grain phenotype.

## R scripts and input data

| Script | Analysis | Required input datasets |
| --- | --- | --- |
| `01_Grain_Morphology_QC_and_Statistical_Analysis.R` | Quality control, descriptive statistics, WT–SG comparisons, and plots of grain length, width, and length-to-width ratio by processing type. | Raw grain morphology measurements. |
| `02_F2_147_Segregation_and_LSU500_Association_Analysis.R` | Phenotypic segregation against a 3 WT:1 SG ratio, LSU500 marker quality control, and single-marker association analysis in the initial 147-plant IC206-derived F₂ population. | Phenotypes and LSU500 genotype calls for the initial 147 plants; LSU500 marker information. |
| `03_Avant_F2_Phenotypic_Segregation_Analysis.R` | Phenotypic segregation against a 3 WT:1 SG ratio in the 100-plant Avant-derived F₂ population. | Grain phenotypes for the 100 Avant-derived F₂ plants. |
| `04_F2_Complete1031_Phenotypic_Segregation_Analysis.R` | Phenotypic segregation against a 3 WT:1 SG ratio in the expanded 1,031-plant IC206-derived F₂ population. | Grain phenotypes from the expanded 1,031-plant population. |

The five input datasets are:

1. **Raw grain morphology measurements:** individual grain measurements, processing classifications, and phenotype labels used in script 01.
2. **Marker information for the LSU500 SNP panel:** marker identifiers, chromosomes, reference-genome positions, and allele definitions used in script 02.
3. **Phenotypes and LSU500 genotype calls for the initial 147 IC206-derived F₂ plants:** plant identifiers, observed grain phenotypes, and marker genotype calls used in script 02.
4. **Grain phenotypes for the 100-plant Avant-derived F₂ population:** individual plant records used in script 03.
5. **Grain phenotypes and chromosome 5 marker genotypes from the expanded 1,031-plant IC206-derived F₂ population:** script 04 uses the plant identifiers and grain phenotypes; chromosome 5 marker genotypes are not analyzed by that script.

## Running the analyses

1. Download the scripts and their required input datasets.
2. Open the desired script in R or RStudio and install the packages it requires.
3. Set the input file paths or select the files when prompted. Use the worksheet and column names expected by the script.
4. Run the script from the beginning. Each analysis can be run separately.

The initial 147 IC206-derived F₂ plants are included in the expanded 1,031-plant population and must not be added to that total again.

## Repository status

The repository is being prepared for manuscript submission. Input filenames, script references, and supplementary table labels will be aligned with the final manuscript numbering. Datasets are identified above by description to avoid confusion with earlier table numbers.

## Contact

Maria Montiel
