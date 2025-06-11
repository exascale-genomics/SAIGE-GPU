# SAIGE-GPU

## Introduction
SAIGE-GPU is a modified version of [SAIGE](https://saigegit.github.io/SAIGE-doc/) (1) where step 1 takes advantage of GPU infrastructures to fit the null model by using the full GRM approach.
To take full advantage of SAIGE-GPU you must have access to a local HPC with GPU infrastructure or have access to a cloud platform offering GPU machines.

Genome-wide association studies (GWAS) analysis enables the identification of broad phenotypic associations with all polymorphic DNA variants in the genome among participants.
A genome-wide Phenome-Wide Association Study (gwPheWAS) is even more challenging in that it essentially involves running a series of GWAS analyses on multiple traits.
One of the major computational hurdles in conducting GWAS analyses on thousands of phenotypes derived from electronic health records can be the repetitive step of fitting a null model.

To address the computational demands of running thousands of GWAS analyses, we made enhancements to the SAIGE step 1 method. 
While SAIGE was initially designed for Central Processing Unit (CPU) based computations, we adapted it to leverage both CPUs and Graphic Processing Units (GPUs) on the DOE Oak Ridge Leadership Computing Facility (OLCF) Summit High-Performance Computer (HPC). 
This adaptation resulted in a substantial acceleration of the analysis process, primarily due to the extensive matrix operations required for performing linear regressions in a GWAS analysis. 

An essential component of the GWAS analysis involves constructing a Genetic Relationship Matrix (GRM), which quantifies the genetic relatedness or similarity among individuals in the study cohort. 
SAIGE offers users the flexibility to create either a sparse GRM or a full GRM. While a sparse GRM is advantageous for its speed of execution and lower memory requirements, opting for a full GRM takes into account pairwise relatedness between all individuals, resulting in a more precise depiction of genetic relatedness. 
This can be highly beneficial for numerous downstream analyses, such as estimating heritability, assessing genetic correlations, and gaining a more comprehensive understanding of the genetic architecture of the trait under investigation (2). 
Furthermore, the full GRM can be applied in various applications of SAIGE-based tools.

The SAIGE method consists of two primary steps: 
- Step 1: involves fitting a null model to establish the relationship between a set of SNPs and the phenotypes; and
- Step 2: tests a larger set of SNPs for associations with the phenotypes.

In our modifications, we optimized step 1 by using GPUs for matrix-vector operations in the calculation of the GRM and employed MPI to distribute the data across multiple GPUs. 
While the standard SAIGE method on CPU-based machines was suitable for relatively small cohorts, it became impractical for larger cohort ancestries due to the substantial matrix-vector operations involved. 
Our modifications allowed us to distribute the genotype matrix required for the GRM calculations across multiple GPUs, resulting in reduced computation time through parallel processing. 

Step 2 in SAIGE can present a distinct challenge due to the need for millions of association tests for each trait that need to be managed. We are making significant efforts to reduce this challenge by implementing methods to use MPI and GPUs as well.

We have used this implementation across thousands of phenotypes derived from electronic health records of over [650,000 participants from the Veterans Affairs (VA) Million Veteran Program (MVP)](https://www.medrxiv.org/content/10.1101/2023.06.28.23291975v1) (3).

## Inputs
### Step 1
For the SAIGE-GPU step 1 implementation, we have used the inputs required by the SAIGE methods in the same formats:

- PLINK BED/BIM/FAM files for genotyping (required)
- Phenotype file (required)
- Covariates file (optional)
- SNP Include file (optioanl)

We have not tested step 1 with input VCF files.

### Step 2
#### Manifest File Usage

##### Overview

The manifest file is an optional tab-delimited file that allows you to specify multiple traits for analysis in a single file, rather than providing comma-separated lists via command line arguments.

##### File Format

The manifest file must be a **tab-delimited** text file with exactly **3 columns** and **no header row**.

###### Column Structure

| Column | Content | Description |
|--------|---------|-------------|
| 1 | GMMATmodelFiles | Path to GMMAT model files for each trait |
| 2 | varianceRatioFiles | Path to variance ratio files for each trait |
| 3 | SAIGEOutputFiles | Path to SAIGE output files for each trait |

###### Example File Content

```
/path/to/trait1_gmmat.rda	/path/to/trait1_variance.txt	/path/to/trait1_output.txt
/path/to/trait2_gmmat.rda	/path/to/trait2_variance.txt	/path/to/trait2_output.txt
/path/to/trait3_gmmat.rda	/path/to/trait3_variance.txt	/path/to/trait3_output.txt
```

##### Usage

###### Option 1: Command Line Arguments (without manifest file)
```bash
Rscript your_script.R \
  --GMMATmodelFile "/path/to/trait1.rda,/path/to/trait2.rda" \
  --varianceRatioFile "/path/to/var1.txt,/path/to/var2.txt" \
  --SAIGEOutputFile "/path/to/out1.txt,/path/to/out2.txt"
```

###### Option 2: Manifest File
```bash
Rscript your_script.R --manifestFile /path/to/manifest.txt
```

##### File Requirements

- **Format**: Tab-delimited text file
- **Columns**: Exactly 3 columns
- **Rows**: At least 1 row (no empty files)
- **Content**: No empty cells or NA values
- **Encoding**: Plain text (UTF-8 recommended)

##### Validation

The script automatically validates:
- ✅ File existence
- ✅ Correct number of columns (3)
- ✅ Non-empty file
- ✅ No missing or empty values
- ⚠️ File path existence (warnings for missing files)

##### Creating a Manifest File

###### Method 1: Text Editor
1. Open any text editor
2. Enter file paths separated by tabs
3. Save as `.txt` file

###### Method 2: Excel/Spreadsheet
1. Create a 3-column spreadsheet
2. Fill in the file paths
3. Save as "Tab-delimited text" (.txt)

###### Method 3: R Script
```r
# Create manifest programmatically
manifest_data <- data.frame(
  gmmat = c("/path/to/trait1.rda", "/path/to/trait2.rda"),
  variance = c("/path/to/var1.txt", "/path/to/var2.txt"),
  output = c("/path/to/out1.txt", "/path/to/out2.txt")
)

write.table(manifest_data, "manifest.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

##### Common Issues

| Issue | Solution |
|-------|----------|
| "Found X columns" error | Ensure exactly 3 tab-separated columns |
| "Empty file" error | Add at least one row of data |
| "Empty or NA values" error | Check for missing values in any column |
| File warnings | Verify all file paths exist and are accessible |

##### Tips

- Use absolute paths to avoid path resolution issues
- Check file permissions if you get access errors
- Keep file paths consistent (all relative or all absolute)
- Test with a small manifest file first

## Outputs
The outputs generated in step 1 SAIGE-GPU are exactly the same as the ones generated by the original SAIGE methods and can be used directly into step 2.
The outputs for step 2 are also the same as previous generated versions.

## References
1. W. Zhou, J. B. Nielsen, L. G. Fritsche, R. Dey, M. E. Gabrielsen, B. N. Wolford, J. LeFaive, P. VandeHaar, S. A. Gagliano, A. Gifford, L. A. Bastarache, W.-Q. Wei, J. C. Denny, M. Lin, K. Hveem, H. M. Kang, G. R. Abecasis, C. J. Willer, S. Lee, Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. Nat. Genet. 50, 1335–1341 (2018).
2. Lee, S. H., Yang, J., Goddard, M. E., Visscher, P. M., & Wray, N. R. (2012). Estimation of pleiotropy between complex diseases using single-nucleotide polymorphism-derived genomic relationships and restricted maximum likelihood. Bioinformatics, 28(19), 2540-2542.
3. Anurag Verma, Jennifer E Huffman, Alex Rodriguez, Mitchell Conery, Molei Liu, Yuk-Lam Ho, Youngdae Kim, David A Heise, Lindsay Guare, Vidul Ayakulangara Panickan, Helene Garcon, Franciel Linares, Lauren Costa, Ian Goethert, Ryan Tipton, Jacqueline Honerlaw, Laura Davies, Stacey Whitbourne, Jeremy Cohen, Daniel C Posner, Rahul Sangar, Michael Murray, Xuan Wang, Daniel R Dochtermann, Poornima Devineni, Yunling Shi, Tarak Nath Nandi, Themistocles L Assimes, Charles A Brunette, Robert J Carroll, Royce Clifford, Scott Duvall, Joel Gelernter, Adriana Hung, Sudha K Iyengar, Jacob Joseph, Rachel Kember, Henry Kranzler, Daniel Levey, Shiuh-Wen Luoh, Victoria C Merritt, Cassie Overstreet, Joseph D Deak, Struan F A Grant, Renato Polimanti, Panos Roussos, Yan V Sun, Sanan Venkatesh, Georgios Voloudakis, Amy Justice, Edmon Begoli, Rachel Ramoni, Georgia Tourassi, Saiju Pyarajan, Philip S Tsao, Christopher J O’Donnell, Sumitra Muralidhar, Jennifer Moser, Juan P Casas, Alexander G Bick, Wei Zhou, Tianxi Cai, Benjamin F Voight, Kelly Cho, Michael J Gaziano, Ravi K Madduri, Scott M Damrauer, Katherine P Liao. Diversity and Scale: Genetic Architecture of 2,068 Traits in the VA Million Veteran Program. medRxiv 2023.06.28.23291975; doi: https://doi.org/10.1101/2023.06.28.23291975
