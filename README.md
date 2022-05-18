# Genome-wide association study (GWAS) and construction of polygenic risk scores (PRS) for height and weight

### Authors:
- Mark Zorin
- Dmitrii Iliushchenko
- Alexander Rakitko (Supervisor)

## Environment

- MacOS Monterey 12.4 M1 Pro and Ubuntu 20.4 
- [R-4.2.0](https://www.r-project.org)
- plink [v1.9](https://www.cog-genomics.org/plink/)
- [PRice-2](https://www.prsice.info)

## Content
`here description of scripts`

## Description
In this project, the main task is to get acquainted with the genome wide association study (GWAS) and the computing a polygenic risk scores (PRS). GWAS analysis includes several steps of quality control (QC), population stratification, as well as the conduct of the GWAS analysis itself. At the end, the calculation of the PRS based on the statistics of the GWAS will be analyzed. All steps of QC and steps visualization were taken from [Marees et al. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)

## Data
The data is precented by 39041 people from the Russian cohort, gentified using the Illumina Infinium Global Screening Array (GSA) provided by Genotek. Since the Genotek data cannot be provided to third parties, we recommend using [simulated data]() for all steps.

### Structure of Data
`Describe Bam .fam, .bim, .bed or binary variant and how to extract to one folder`

### Quality Control
Before we start make sure you have [plink](https://www.cog-genomics.org/plink/)  and [R](https://www.r-project.org) installed.

1. Filtering samlpes and positions (SNP) by call rate
  
   First step is filtering positions and smaples that are missing in a large proportion of the subjects
   
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

3. 
   
  
