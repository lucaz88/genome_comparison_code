# Collection of R scripts used in the genome comparison analysis

- _ATLAS_annotation.R_ is a wrapping workflow for all the functional annotation performed in the manuscript: prokka, KEGG Orthology (including phytohormones production, DHPS and taurine utilization), BioVx (for membrane transporters), AntiSMASH (for secondary metabolites), Vibrioferrin biosynthesis and transport (blastp against UniProt), DMSP degradation pathways (blastp against UniProt). It also loads and run the _KM_reconstruction.R_ script. Please note that, as specifide in the code, the wrapper depends on several packages/tools that need to be already installed in your system. Also, the related executables have to be directly available in the terminal by adding their folders to the system $PATH.

- _KM_reconstruction.R_ includes two functions that allow to recover the diagrams of all KEGG modules (KMs) through the KEGG API package and to reconstruct KMs completness from a table of annotated KEGG Orthology (KOs)

  - _KMdiagram_fetcher_ allows to retrieve the diagrams and store them as a list in a RData object. It uses parallelization to speed up the fetching process but don't exceed 6-7 instances as you might get errors due to a bottleneck presumibily in the KEGG API.
  - _KMreco_ uses the retrieved list of diagrams as masks to evaluate the completness of each KMs givens a table of annotated KOs. The input table is a presence/absence table (i.e. 0 and 1) and has samples as rows and annotated KOs as columns. Further arguments are "len_breaks" and "allowed_gaps" that can be used to specify the ammount of allowed gaps allowed in a KM to still call such KM complete depending on its lengths: len_breaks=c(3, 10), allowed_gaps=c(0,1,2) will allows 0 gaps for KM with &lt; 3 reactions, 1 gap for KM with 3 up to 9 reactions, 2 gaps for KM with ≥ 10 reactions. By defauls the function will not allow any gap regardless of the KM length, however, it's recomendable to set these values as issues in the annotation process that could lead to "artifical" gaps are always expected.
  
  USAGE EXAMPLE (type commands directly in R)
```  
# download and unpack the whole repository in a folder, e.g. ~/Downloads/genome_comparison_code-main
setwd("~/Downloads/genome_comparison_code-main")
source("KM_reconstruction.R")

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
