# Genetic Mapping of a Putative Spontaneous Short-Grain Mutation Identifies a Predicted Loss-of-Function Allele of *SRS3* in Rice

## Overview

This repository contains the data, metadata, supporting tables, and R scripts associated with the manuscript **“Genetic mapping of a putative spontaneous short-grain mutation identifies a predicted loss-of-function allele of *SRS3* in rice”** by Maria Montiel, Brijesh Angira, Jonathan K. Richards, and Adam Famoso.

The study characterizes a recessive short-grain (SG) phenotype identified in the long-grain wild-type (WT) rice breeding line RU2002174. The analyses include grain-morphology comparisons, phenotypic segregation, LSU500 marker association, recombinant fine-mapping, F3 progeny testing, whole-genome sequence analysis, candidate-variant identification, and evaluation of the candidate *SRS3* variant in additional rice germplasm.

Supplementary tables retain the same numbering and nomenclature used in the manuscript. In the repository, each Supporting Table is provided as an individual Excel file named according to its manuscript designation (e.g., `S1_Table.xlsx`, `S2_Table.xlsx`, `S8a_Table.xlsx`). The full title and description of each table are provided below and correspond directly to the Supporting Information cited in the manuscript. Additional datasets required to reproduce the reported analyses but not included as Supporting Tables in the manuscript are designated separately as S1 Data, S2 Data, etc.

Required R packages are specified within each script. Scripts check for the packages required for the corresponding analysis before execution. Package installation is not performed automatically.

The repository and article Supporting Information provide the data and metadata underlying the reported phenotypic, genetic-mapping, fine-mapping, candidate-variant, and sequencing analyses. These materials include grain-measurement data and associated quality-control records; F2 and F3 phenotypes and marker genotypes; recombinant classifications and F3 progeny-test results; LSU500 and PACE marker metadata, genomic positions, alleles, and primer sequences, including the LSU1761 candidate-variant assay; whole-genome sequencing and candidate-site metrics; filtered and annotated SNP and short-indel data from the fine-mapped chromosome 5 interval; pairwise SNP identity results for the seven sequenced samples; candidate-variant genotypes and phenotype-concordance results for F2 recombinants; Native Trait Panel survey results; IRRI SNP-Seek 3K regional query results; and the R scripts used for the reported grain-morphology, segregation, and LSU500 marker-association analyses.

## Repository contents

### R scripts

- `01_Grain_Morphology_QC_and_Statistical_Analysis_Corrected.R`  
  Performs grain-measurement quality control, removes statistical outliers using Tukey's 1.5 × interquartile-range criterion, summarizes grain length, width, and length-to-width ratio, and compares WT and SG samples using two-sided Welch's t-tests with Benjamini-Hochberg adjustment.  
  **Input:** S1 Data (raw grain morphology measurements).  
  **Output:** S4 Table (grain morphology statistics and quality-control results).

- `02_F2_147_Segregation_and_LSU500_Association_Analysis.R`  
  Evaluates 3:1 phenotypic segregation in the initial 147-plant IC206 × RU2002174 SG F₂ population and conducts LSU500 single-marker association analysis using logistic-regression likelihood-ratio tests, Bonferroni and Benjamini-Hochberg corrections, McFadden's pseudo-R², and phenotype proportions by genotype class.  
  **Inputs:** S1 Table (LSU500 marker information) and S5 Table (phenotypes and LSU500 genotypes of the 147 F₂ plants).  
  **Outputs:** Segregation results reported in S5 Table; LSU500 marker quality-control and association results reported in S7 Table; and the association results represented in S1 Fig.

- `03_Avant_F2_Phenotypic_Segregation_Analysis.R`  
  Evaluates phenotypic segregation in the 100-plant Avant × RU2002174 SG F₂ population using a chi-square goodness-of-fit test for the expected 3:1 WT:SG ratio.  
  **Input:** S6 Table (Avant × RU2002174 SG F₂ phenotypes).  
  **Output:** Chi-square segregation statistics displayed in the R console and reported in S6 Table and the manuscript.

- `04_F2_Complete1031_Phenotypic_Segregation_Analysis.R`  
  Evaluates phenotypic segregation in the complete 1,031-plant IC206 × RU2002174 SG F₂ population using a chi-square goodness-of-fit test for the expected 3:1 WT:SG ratio.  
  **Input:** S2 Data (phenotypes of the complete 1,031-plant F₂ population).  
  **Output:** Chi-square segregation statistics displayed in the R console and reported in the manuscript.

### Data and supporting information

- **S1 Data. Grain morphology measurements.** Raw individual-grain measurements of length, width, and length-to-width ratio for rough, brown, and milled rice from wild-type (WT) and short-grain (SG) samples. These data are used as input for the grain morphology quality-control and statistical analyses reported in S4 Table.

- **S2 Data. Phenotypes of the complete 1,031-plant IC206 × RU2002174 SG F₂ population.** Plant identifiers and WT or SG grain-phenotype classifications used as input for the phenotypic segregation analysis of the complete F₂ population. The 1,031 plants include the initial 147-plant subset used for coarse mapping. Corresponding chromosome 5 marker genotypes and recombinant classifications are provided in S8a Table.

## Supporting tables

- **S1 Table. Marker information for the LSU500 AmpSeq panel.** Marker identifiers, chromosomes, physical positions in the Nipponbare IRGSP-1.0 reference genome, and designated alleles for the 550 SNP markers included in the LSU500 panel.

- **S2 Table. PACE markers and primer sequences used for chromosome 5 fine-mapping and candidate-variant genotyping.** The table provides marker numbers, LSU assay identifiers, chromosome positions in the Nipponbare IRGSP-1.0 reference genome, parental alleles, and allele-specific and common primer sequences. LSU1761 was developed in this study to genotype the candidate G→T substitution in SRS3 at chromosome 5 position 3,204,936 bp.

- **S3 Table. LSU500 genotypes used to compare short-grain, wild-type, and historical RU2002174 samples.** Genotypes are reported for three SG and nine WT plants selected from one segregating RU2002174 panicle row in the 2019–2020 Puerto Rico winter nursery and for six RU2002174 samples genotyped in breeding trials from 2018 through 2020. Of the 550 markers assayed, 18 were excluded because they had missing calls in at least five samples, and five were excluded because at least nine samples had heterozygous calls. Of the remaining 527 markers, 130 were classified as heterogeneous among the six historical RU2002174 samples. The remaining markers were used to compare the SG and WT plants with the RU2002174 consensus genotype.

- **S4 Table. Grain morphology of wild-type RU2002174 and the RU2002174 short-grain mutant.** Grain length, width, and length-to-width ratio are summarized for rough, brown, and milled rice. Comparisons used two-sided Welch’s t-tests with Benjamini–Hochberg adjustment across the nine trait-by-processing combinations. Statistical outliers were removed using the 1.5 × interquartile-range criterion. Because grains were subsampled from bulk seed lots rather than independently grown plants, the comparisons are exploratory and do not represent inference among biological replicates.

- **S5 Table. Phenotypes and LSU500 genotypes of 147 F₂ plants from the IC206 × RU2002174 SG cross.** The population included 105 WT and 42 SG plants and was used to evaluate segregation and conduct the initial genome-wide marker-association analysis. The table reports grain phenotype and genotype calls for the 550 markers in the LSU500 panel. Marker quality-control and association results are provided in S7 Table.

- **S6 Table. Grain phenotypes and segregation analysis of 100 F₂ plants from the Avant × RU2002174 SG cross.** Observed WT and SG phenotype counts used to compare with the 3:1 ratio expected for a single recessive locus using a chi-square goodness-of-fit test (expected: 75 WT and 25 SG; χ² = 0.21, df = 1, P = 0.644).

- **S7 Table. LSU500 marker quality control and genome-wide association results for the short-grain phenotype.** Results are reported for the 265 markers retained from the LSU500 panel after filtering for call rate ≥0.80 and minor allele frequency >0.20 in the 147-plant IC206 × RU2002174 SG F₂ population. For each marker, the table provides marker-quality metrics, categorical logistic-regression results, McFadden’s pseudo-R², unadjusted P values, Bonferroni-adjusted P values, Benjamini–Hochberg-adjusted P values, and phenotype proportions by genotype class.

- **S8a Table. Phenotypes and chromosome 5 marker genotypes of 1,031 F₂ plants used for initial recombinant mapping.** The table reports grain phenotype and marker genotypes for the expanded IC206 × RU2002174 SG F₂ population. Recombinant plants were identified from changes in parental marker haplotypes across chromosome 5 and were classified as single or double recombinants. These data supported initial localization of the SG locus between SNP5 and SNP7.

- **S8b Table. Additional marker genotypes used to refine the short-grain locus within the SNP5–SNP7 interval.** Eight additional PACE markers were genotyped among informative F₂ recombinants within the chromosome 5 interval from 3.00 to 3.95 Mb. Recombinant breakpoints refined the SG locus to the approximately 233-kb interval between SNP5B and SNP6B.

- **S8c Table. Recombinant F₂ genotypes and F₃ progeny-test results used for high-resolution mapping of the short-grain locus.** The table provides marker information, observed or inferred marker genotypes of 14 recombinant F₂ parents, F₂ phenotypes, F₃ family sizes, and WT and SG F₃ phenotype counts. F₂ marker genotypes were directly determined for nine parents and inferred from eight genotyped F₃ progeny for five parents. The combined recombinant and progeny-test results localized the SG locus to the 41.6-kb interval between SNP5B_6 and SNP5C_1 at chromosome 5 positions 3,189,947–3,231,570 bp.

- **S9 Table. Sequencing summary statistics and candidate-site genotype metrics for the seven sequenced samples.** The table reports total reads, mapped-read percentage, mean coverage based on mapped bases, and genotype metrics at the candidate SRS3 variant at chromosome 5 position 3,204,936. Total reads represent individual reads rather than read pairs. Candidate-site metrics include genotype, variant-call depth, reference and alternate allele depths, and genotype quality.

- **S10 Table. Pairwise SNP identity among the seven sequenced samples.** Pairwise identity was calculated from 2,985,920 genome-wide sites that were polymorphic among the seven samples after filtering for DP ≥10, GQ ≥30, and genotype calls in all samples. Genotype pairs sharing both alleles were scored 1.0, those sharing one allele 0.5, and those sharing no alleles 0.0. Values therefore represent relative SNP identity among samples rather than whole-genome sequence identity.

- **S11 Table. Filtered SNPs and short indels within the fine-mapped chromosome 5 interval.** The table reports 213 variants retained between chromosome 5 positions 3,190,345 and 3,229,775 bp within the marker-defined interval of 3,189,947–3,231,570 bp after requiring depth ≥10, genotype quality ≥30, biallelic status, and genotype calls in all seven sequenced samples. Effects were predicted using SnpEff v4.3t and a custom database based on the Nipponbare IRGSP-1.0 reference genome and MSU Release 7 gene models. REF indicates homozygosity for the reference allele, ALT indicates homozygosity for the alternate allele, and HET indicates heterozygosity. The G→T substitution at position 3,204,936 was the only retained variant consistently distinguishing the WT- and SG-class samples.

- **S12a Table. Candidate-variant genotypes and grain phenotypes of F₂ recombinants.** LSU1761 PACE genotype calls for the candidate G→T variant at chromosome 5 position 3,204,936 are reported with the grain phenotypes and mapping classifications of the recombinant F₂ plants. G:G indicates homozygosity for the WT allele, G:T or T:G indicates heterozygosity, and T:T indicates homozygosity for the SG-associated allele. Failed or ambiguous calls are identified.

- **S12b Table. Candidate-variant survey of the Native Trait Panel.** The table provides accession identifiers, germplasm metadata, grain class, and LSU1761 PACE genotype calls for the candidate SRS3 variant at chromosome 5 position 3,204,936. G:G indicates homozygosity for the WT allele, G:T or T:G indicates heterozygosity, and T:T indicates homozygosity for the SG-associated allele. “No Amp” indicates failure to obtain a genotype call. Grain class describes the recognized grain type of each accession and does not indicate the RU2002174 SG phenotype.

- **S13 Table. IRRI SNP-Seek regional query surrounding the candidate SRS3 variant on chromosome 5.** Regional SNP data were retrieved from the IRRI SNP-Seek 3K dataset for 3,024 accessions across chromosome 5 positions 3,204,800–3,205,000 bp. The query returned 12 variant positions, but the candidate position at 3,204,936 was not returned as a variant site. The table provides the returned SNP positions and accession-level allele calls. Blank cells indicate missing calls.

## Reproducibility

The R scripts are numbered according to the analyses performed in the manuscript. Each script identifies its required input file and should be run from the repository root unless otherwise indicated. Required R packages are specified within each script, and package installation is not performed automatically.

The scripts reproduce the statistical analyses and summaries reported in the manuscript and Supporting Information. Some scripts generate Supporting Tables as output files, whereas analyses for which no separate output file is required report the relevant statistics directly in the R console.

Statistical analyses were conducted in R. Bioinformatic software, versions, filtering criteria, and parameters used for sequence processing, alignment, variant calling, filtering, and annotation are documented in the Materials and Methods section of the manuscript.

## License

The author-generated R analysis scripts in this repository are licensed under the MIT License. See the `LICENSE` file for details.

## Citation

When using data or analysis scripts from this repository, please cite the associated manuscript:

Montiel, M., Angira, B., Richards, J. K., & Famoso, A. Genetic mapping of a putative spontaneous short-grain mutation identifies a predicted loss-of-function allele of *SRS3* in rice. Manuscript submitted to PLOS ONE.

## Contact

For questions about the repository or associated analyses, contact:

**Adam Famoso**  
H. Rouse Caffey Rice Research Station  
Louisiana State University Agricultural Center  
Email: afamoso@agcenter.lsu.edu
