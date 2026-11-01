# Genetic Mapping of a Putative Spontaneous Short-Grain Mutation Identifies a Predicted Loss-of-Function Allele of *SRS3* in Rice

## Overview

This repository contains the data, metadata, supporting tables, and R scripts associated with the manuscript **“Genetic mapping of a putative spontaneous short-grain mutation identifies a predicted loss-of-function allele of *SRS3* in rice”** by Maria Montiel, Brijesh Angira, Jonathan K. Richards, and Adam Famoso.

The study characterizes a recessive short-grain (SG) phenotype identified in the long-grain wild-type (WT) rice breeding line RU2002174. The analyses include grain-morphology comparisons, phenotypic segregation, LSU500 marker association, recombinant fine-mapping, F3 progeny testing, whole-genome sequence analysis, candidate-variant identification, and evaluation of the candidate *SRS3* variant in additional rice germplasm.

## Repository contents

### R scripts

- `01_Grain_Morphology_QC_and_Statistical_Analysis_Corrected.R`  
  Performs grain-measurement quality control, removes statistical outliers using Tukey's 1.5 x interquartile-range criterion, summarizes grain length, width, and length-to-width ratio, and compares WT and SG samples using Welch's t-tests with Benjamini-Hochberg adjustment.

- `02_F2_147_Segregation_and_LSU500_Association_Analysis.R`  
  Evaluates 3:1 phenotypic segregation in the initial 147-plant IC206 x RU2002174 SG F2 population and conducts LSU500 single-marker association analysis using logistic-regression likelihood-ratio tests, Bonferroni and Benjamini-Hochberg corrections, McFadden's pseudo-R2, and phenotype proportions by genotype class.

- `03_Avant_F2_Phenotypic_Segregation_Analysis.R`  
  Evaluates phenotypic segregation in the 100-plant Avant x RU2002174 SG F2 population using a chi-square goodness-of-fit test for the expected 3:1 WT:SG ratio.

- `04_F2_Complete1031_Phenotypic_Segregation_Analysis.R`  
  Evaluates phenotypic segregation in the complete 1,031-plant IC206 x RU2002174 SG F2 population using a chi-square goodness-of-fit test for the expected 3:1 WT:SG ratio.

### Data and supporting information

The repository and article Supporting Information provide the materials underlying the reported analyses, including:

- Grain-measurement data and quality-control records for rough, brown, and milled rice.
- F2 and F3 phenotypes and marker genotypes.
- Recombinant classifications and F3 progeny-test results used for fine-mapping.
- LSU500 and PACE marker metadata, genomic positions, alleles, and primer sequences, including the LSU1761 candidate-variant assay.
- Whole-genome sequencing metrics and candidate-site quality metrics.
- Filtered SNP and short-indel data from the chromosome 5 fine-mapped interval.
- The annotated interval-variant table used to identify the candidate G-to-T substitution at chromosome 5 position 3,204,936.
- Pairwise SNP identity results for the seven sequenced samples.
- Candidate-variant genotypes and phenotype-concordance results for F2 recombinants.
- Native Trait Panel survey results and the IRRI SNP-Seek 3K regional query.
- Statistical-analysis scripts used for the grain-morphology, segregation, and LSU500 association analyses.

The merged Supporting Information workbook contains Tables S1-S13. Supporting figures are provided with the article.

## Supporting tables

- **S1 Table:** LSU500 marker information.
- **S2 Table:** PACE marker information and primer sequences.
- **S3 Table:** LSU500 genotypes for RU2002174 SG, WT, and historical samples.
- **S4 Table:** Grain-morphology measurements, quality control, summary statistics, and exploratory WT-SG comparisons.
- **S5 Table:** Initial 147-plant IC206-derived F2 population and segregation analysis.
- **S6 Table:** Avant-derived F2 population and segregation analysis.
- **S7 Table:** LSU500 marker filtering and single-marker association results.
- **S8a-S8c Tables:** Expanded F2 mapping, recombinant classifications, high-resolution mapping, and F3 progeny tests.
- **S9 Table:** Whole-genome sequencing and candidate-site metrics.
- **S10 Table:** Pairwise SNP identity among the seven sequenced samples.
- **S11 Table:** Filtered and annotated variants in the chromosome 5 fine-mapped interval.
- **S12a-S12b Tables:** LSU1761 genotype-phenotype concordance and Native Trait Panel survey.
- **S13 Table:** IRRI SNP-Seek 3K regional query results.

## Reproducibility

The R scripts are numbered in the recommended order of use. Each script should be run with its corresponding input dataset available in the working directory or with the file path updated in the script. The scripts generate the statistical summaries and results reported in the manuscript and Supporting Information.

Analyses were conducted in R. Additional bioinformatic software and parameters used for sequence processing, alignment, variant calling, filtering, and annotation are documented in the manuscript Methods section.

## Data availability

Data and analysis scripts are available in this GitHub repository:

https://github.com/MMontiel-breeding/rice-short-grain-SRS3

The repository will be archived through Zenodo after the final GitHub release.

**Zenodo DOI:** To be added after repository archiving.

Raw sequencing reads are associated with **NCBI BioProject PRJNA1472089**. Public accessibility of the BioProject should be confirmed before final submission.

## Citation

Please cite the associated article when using these data or scripts. The complete journal citation and DOI will be added after publication.

## License

License information will be added before the final repository release.

## Contact

For questions about the repository or associated analyses, contact:

**Adam Famoso**  
H. Rouse Caffey Rice Research Station  
Louisiana State University Agricultural Center  
Email: afamoso@agcenter.lsu.edu
