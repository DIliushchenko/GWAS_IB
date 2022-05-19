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
The data is precented by 39041 people from the Russian cohort, gentified using the Illumina Infinium Global Screening Array (GSA) v1.0 / v2.0 / v3.0 provided by Genotek. Scan images processing and genotypes calling were performed using GenomeStudio v2.0. Since the Genotek data cannot be provided to third parties, we recommend using [simulated data](https://github.com/MareesAT/GWA_tutorial) for all steps.

## Structure of Data
`Describe Bam .fam, .bim, .bed or binary variant and how to extract to one folder. Input .png for better visualization`

## Quality Control
Before we start make sure you have [plink](https://www.cog-genomics.org/plink/)  and [R](https://www.r-project.org) installed.

1. Filtering samlpes and positions (SNP) by call rate
  
   First step is filtering positions and samples that are missing in a large proportion of the subjects
   
   a. Removal of SNPs with missingness of genotype in more than 20\% cases
   ```
   plink --bfile raw_file --geno 0.2 --make-bed --out file_1
   ```
   b. Removal of individuals with missingness of genotype in more than 20\%  cases
   ```
   plink --bfile file_1 --mind 0.2 --make-bed --out file_2
   ```
   c. Removal of SNPs with missingness of genotype in more than 2\% cases
   ```
   plink --bfile file_2 --geno 0.02 --make-bed --out file_3
   ```
   d. Removal of individuals with missingness of genotype in more than 2\%  cases
   ```
   plink --bfile file_3 --mind 0.02 --make-bed --out file_4
   ```
   
2. Filtering samples with a discrepancy between the genetic sex and the sex indicated in the personal data

   Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate. 
   
   a. Mark subjects who do not fulfil these requirements as "PROBLEM"
   ```
   plink --bfile file_4 --check-sex 
   ```
   b. Filter samples
   ```
   plink --bfile file_4 --impute-sex --make-bed --out file_5
   ```
   
3. Heterozygosity check
   
   Excluding samples in which the observed heterozygosity deviated by more than 3 standard deviations from the sample mean. Heterozygosity was assessed for the cohort after filtering genetic variants in linkage disequilibrium (search window - 50 SNPs, number of SNPs for window shift at the end of the step - 5, r2 between SNPs < 0.2)
   a. Generate a list of independent SNP
   ```
   plink --bfile file_5 --exclude high_inversion_regions.txt --range --indep-pairwise 50 5 0.2 --out independent_SNP
   ```
   b. Generate file contains your pruned data set
   ```
   plink --bfile file_5 --extract independent_SNP.prune.in --het --out R_check
   ```
   c. Generate a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean    
   ```
   Rscript --no-save heterozygosity_outliers_list.R
   ```
   d. Adapt output file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns    
   ```
   sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
   ```
   e. Remove heterozygosity rate outliers
   ```
   plink --bfile file_5 --remove het_fail_ind.txt --make-bed --out file_6
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
