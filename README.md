# Genome-wide association study (GWAS) and construction of polygenic risk scores (PRS) for height and weight

### Authors:
- Mark Zorin
- Dmitrii Iliushchenko
- Alexander Rakitko (Supervisor)

## Environment
- MacOS Monterey 12.4(M1 support) or Ubuntu 20.4 
- [R-4.2.0](https://www.r-project.org)
- [plink v1.9](https://www.cog-genomics.org/plink/)
- [PRice-2](https://www.prsice.info)
- [PRIMUS v1.9](https://primus.gs.washington.edu/primusweb/index.html)

## Content
`here description of scripts if there any of them`

## Description
In this project, the main task is to get acquainted with the genome wide association study (GWAS) and the computing a polygenic risk scores (PRS). GWAS analysis includes several steps of quality control (QC), population stratification, as well as the conduct of the GWAS analysis itself. At the end, the calculation of the PRS based on the statistics of the GWAS will be analyzed. All steps of QC and steps visualization were taken from [Marees et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)

## Data
The data is precented by 39041 people from the Russian cohort, gentified using the Illumina Infinium Global Screening Array (GSA) provided by Genotek. Since the Genotek data cannot be provided to third parties, we recommend using [simulated data]() for all steps.

## Structure of Data
`Describe Bam .fam, .bim, .bed or binary variant and how to extract to one folder. Input .png for better visualization`

## Quality Control
Before we start make sure you have [plink](https://www.cog-genomics.org/plink/)  and [R](https://www.r-project.org) installed.

1. Filtering samlpes and positions (SNP) by call rate
  
   First step is filtering positions and samples that are missing in a large proportion of the subjects
   
   a.
   ```
   ```
   b.
   ```
   ```
   c.
   ```
   ```
   d.
   ```
   ```
   
2. Sex discrepancy

   Filtering samples with a discrepancy between the genetic sex and the sex indicated in the personal data
   ```
   ```
3. Heterozygosity
   
   Excluding samples in which the observed heterozygosity deviated by more than 3 standard deviations from the sample mean. Heterozygosity was assessed for the cohort after filtering genetic variants in linkage disequilibrium (search window - 50 SNPs, number of SNPs for window shift at the end of the step - 5, r2 between SNPs < 0.2)
   ```
   ```
4. MAF filtering

   Filtering out of positions based on minor allele frequency (MAF), while SNP with low MAF often associated with genotyping errors. Threshold depends on number of positions contained in data, in our case we will use **threshold 0.01**
   ```
   ```
 5. Rremove positions on sex chromosomes and mitochondria
 
   ```
   ```
 
 6. Hardy–Weinberg equilibrium
 
    Filtering of positions with inclination of the Hardy-Weinberg equilibrium: positions with a significant difference between the observed genotype frequencies and those expected according to the exact Hardy-Weinberg test (p-value < 1e-05) were removed.
    ```
    ```
 7. Relatedness
    
    Identification of close relatives individuals was carried out using the [PRIMUS](https://primus.gs.washington.edu/primusweb/res/documentation.html) program (Rapid Reconstruction of Pedigrees from Genome-wide Estimates of Identity by Descent). Pairs with PI_HAT > 0.15 were considered relatives. With the help of PRIMUS, we can obtain dataset, consisting of pairwise unrelated individuals.
    
  ```
  ```
## Popilation Stratification


## GWAS

## Results 

## Literature
