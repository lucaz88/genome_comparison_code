# Collection of R scripts used in the genome comparison analysis

- _ATLAS_annotation.R_ is a wrapping workflow for all the functional annotation performed in the manuscript: prokka, KEGG Orthology (including phytohormones production, DHPS and taurine utilization), BioVx (for membrane transporters), AntiSMASH (for secondary metabolites), Vibrioferrin biosynthesis and transport (blastp against UniProt), DMSP degradation pathways (blastp against UniProt). It also loads and run the _KM_reconstruction.R_ script. Please note that, as specifide in the code, the wrapper depends on several packages/tools that need to be already installed in your system. Also, the related executables have to be directly available in the terminal by adding their folders to the system $PATH.

- _KM_reconstruction.R_ includes two functions that allow to recover the diagrams of all KEGG modules (KMs) through the KEGG API package and to reconstruct KMs completness from a table of annotated KEGG Orthology (KOs)

  - _KMdiagram_fetcher_ allows to retrieve the diagrams and store them as a list in a RData object. It uses parallelization to speed up the fetching process but don't exceed 6-7 instances as you might get errors due to a bottleneck presumibily in the KEGG API.
  - _KMreco_ uses the retrieved list of diagrams as masks to evaluate the completness of each KMs givens a table of annotated KOs. The input table is a presence/absence table (i.e. 0 and 1) and has samples as rows and annotated KOs as columns. Further arguments are "len_breaks" and "allowed_gaps" that can be used to specify the ammount of allowed gaps allowed in a KM to still call such KM complete depending on its lengths: len_breaks=c(3, 10), allowed_gaps=c(0,1,2) will allows 0 gaps for KM with &lt; 3 reactions, 1 gap for KM with 3 up to 9 reactions, 2 gaps for KM with ≥ 10 reactions. By defauls the function will not allow any gap regardless of the KM length, however, it's recomendable to set these values as issues in the annotation process that could lead to "artifical" gaps are always expected.
  
  USAGE EXAMPLE (type commands directly in R):
```  
# Download Repository from GitHub, copy the link from `Code:` -> `Download ZIP`
repo_url <- "https://github.com/lucaz88/genome_comparison_code/archive/refs/heads/main.zip" 
filename <- "genome_compar_git.zip"
download.file(url=repo_url, destfile=filename)
unzip(zipfile=filename, exdir=".")

# load and run the code
source("./genome_comparison_code-main/KM_reconstruction.R")
KM_str <- KMdiagram_fetcher(ncore = 7, create_RData = T, path = "~/Downloads/genome_comparison_code-main") # it takes a few minutes
myannotation <- read.csv("example_KO_table.csv", header = T, row.names = 1)
KMreco <- KMreco(indata = myannotation, KM_str = KM_str, len_breaks = c(3), allowed_gaps = c(0,1))

# example of the presence/absence table of KMs
KMreco$KM_pa[1:5, 1:5]
```

- _Statistical_analyses.R_ contains the statistical analysis for the creation of Genome Functional Clusters (GFCs) and the Linked-Trait Clusters (LTCs)


# Other files:
- _example_KO_table.csv_ is a test dataset of KOs presence/absence for 5 genomes to test the _KM_reconstruction_ code
- _KM_str_2020-01-09_plus.rds_ contains the KM structures updated to Jannuary 2020 plus the structures of transporters and two-component systems that were removed with the Release 92.0  of KEGG (October 1, 2019) as they were redundant with the BRITE hierarchies (not used in our analysis)
- _howto_BioVx_ my workaround to get a running version of BioVx suite on Ubuntu 20.04 LTS


# Information About the R Session:
```
> sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
 [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] future.apply_1.8.1  future_1.21.0       KEGGREST_1.32.0     xml2_1.3.2          stringr_1.4.0      
[6] BiocManager_1.30.16

loaded via a namespace (and not attached):
 [1] parallelly_1.26.0      rstudioapi_0.13        XVector_0.32.0         magrittr_2.0.1        
 [5] BiocGenerics_0.38.0    zlibbioc_1.38.0        IRanges_2.26.0         R6_2.5.0              
 [9] httr_1.4.2             GenomeInfoDb_1.28.0    globals_0.14.0         tools_4.1.1           
[13] parallel_4.1.1         png_0.1-7              digest_0.6.27          crayon_1.4.1          
[17] GenomeInfoDbData_1.2.6 codetools_0.2-18       S4Vectors_0.30.0       bitops_1.0-7          
[21] curl_4.3.2             RCurl_1.98-1.3         stringi_1.6.2          compiler_4.1.1        
[25] Biostrings_2.60.1      stats4_4.1.1           listenv_0.8.0   
      
> Sys.time()
[1] "2021-09-02 14:32:12 CEST"

> Sys.Date()
[1] "2021-09-02"
```
