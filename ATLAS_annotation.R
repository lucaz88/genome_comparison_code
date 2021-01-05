# set working directory
setwd("/media/lucaz/Bioinfo/HB_complete_genomes_v4/")


# global setting
options(stringsAsFactors = F)
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8") # to get month names in english


#!!! following tools must be correctly installed (check required dependencies) in the system
#a prokka (https://github.com/tseemann/prokka)
#b kofamscan (https://github.com/takaram/kofam_scan)
#c gblast (http://www.tcdb.org/files/gblast.zip)
#d Antismash (https://antismash.secondarymetabolites.org/#!/download)
#e blast+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)


#!!! required files
# Supplementary tables' file
suppl_tabs = "/media/lucaz/DATA/Google_drive/Laboratorio/AAA_Progetti_NUOVI/2016_HFSP/4_paper_genome_comparison/6_REV_Comm_biology/Suppl_tables.xlsx" #! adjust the path
# Supplementary file 3 
VFds = "/media/lucaz/DATA/Google_drive/Laboratorio/AAA_Progetti_NUOVI/2016_HFSP/4_paper_genome_comparison/6_REV_Comm_biology/Suppl_file_3-uniprot-vibrioferrin_01mai20.fasta" #! adjust the path
# Supplementary file 4 
DMSPds = "/media/lucaz/DATA/Google_drive/Laboratorio/AAA_Progetti_NUOVI/2016_HFSP/4_paper_genome_comparison/6_REV_Comm_biology/Suppl_file_4-uniprot-DMSP_21dec20.fasta" #! adjust the path
# KM_reconstruction.R code
KM_code = "/home/lucaz/myscript/GitHub/gnm_compar/KM_reconstruction.R" 
# file with KM structures (optional, it can be downloaded anew) 
KM_str = "/home/lucaz/myscript/GitHub/gnm_compar/KM_str_2020-01-09_plus.rds"



### I annotation ----
dir.create("3_Ann")



### prokka
dir.create("3_Ann/Prokka")
list.gnm = list.files("1_gnm", pattern = ".fna$", full.names = T)

sapply(list.gnm, function(j) system2(command = "prokka", args = c(paste0("--outdir 3_Ann/Prokka/", gsub(".fna", "", basename(j))),
                                                                  j, "--cpus 7", "--compliant", "--rnammer"), stderr = F)) # "--norrna" "--notrna"

# clean up
system("gio trash 3_Ann/Prokka/*/*.fna")
system("gio trash 3_Ann/Prokka/*/*.fsa")
system("gio trash 3_Ann/Prokka/*/*.sqn")

# merge all genomes and add filename in header
#! speeds up other annotations (e.g. KOfam that iterates per KO)
system('ls 3_Ann/Prokka/*/*.faa | while read file; do cat $file | sed "s/^>/>$(basename $file | cut -f 1 -d ".")./g" >> 3_Ann/all_genes.faa; done')





### KEGG
dir.create("3_Ann/KEGG")

#!!! local KOfam search 
#! run on concatenated seqs of predicted peptides from prokka
#! ~12h for 370,000 genes on 7 cores
system2(command = "/home/lucaz/kofamscan-1.2.0/exec_annotation", 
        args = c("-o 3_Ann/KEGG/kofam_result.txt", "3_Ann/all_genes.faa"))





### BioV-transporters
dir.create("3_Ann/Transporter")
#! 80 genomes/day on 5 cores (out of 8 cores; leave ~40% free CPU for overhead usage of some steps)

# #!!! loop genome-wise - in R 
# all_gnm = list.files("3_Ann/Prokka", recursive = T, pattern = ".faa$", full.names = T)
# all_gnm = all_gnm[!sapply(strsplit(all_gnm, "/"), "[", 3) %in% 
#                     list.files("3_Ann/Transporter/", recursive = F, full.names = F)] # resume run
# for (i in all_gnm) {
#   system2(command = "python2.7",
#           args = c("/home/lucaz/BioVx/scripts/gblast3.py",
#                    paste("-i", i),
#                    file.path("-o 3_Ann/Transporter", strsplit(i, "/")[[1]][3])),
#           stdout = F, stderr = F
#   )
# }

#!!! loop genome-wise - in BASH (requires GNU parallel, https://www.gnu.org/software/parallel/)
#! run following command in shell (cd to same working dir); adjust '-j' to you available cores
conda activate BioVx; ls 3_Ann/Prokka/*/*.faa | parallel -j 4 'gblast3.py -i {} -o 3_Ann/Transporter/$(cut -d'/' -f3 <<< {})'

#! if you wanna keep only the essential output:
# tar -czf 3_Ann/Transporter/transp.tar.gz 3_Ann/Transporter/*/results.tsv





### Secondary Metabolites

dir.create("3_Ann/Antismash")

#! run in bash from the same working directory as the one used in R
#! in addition to AntiSMASH, it requires the GNU parallel tool (http://ftp.gnu.org/gnu/parallel/)
#! now R-code wrapping as it would be too messy, just change the '-j' parameter to meet your available cores

# ls 3_Ann/prokka_ann/*.gbk | parallel -j 30 --tmpdir _tmp 'antismash -c 2 --taxon bacteria --fullhmmer --cf-create-clusters --cb-general --cb-knownclusters --cb-subclusters --tta-threshold 0.65 --asf --pfam2go --smcog-trees --genefinding-tool none --output-dir 3_Ann/Antismash/$(basename {} | sed "s/.gbk//g") {}

#! if you wanna keep only the essential output:
# tar -czf 3_Ann/Antismash/antismash_html.tar.gz 3_Ann/Antismash/*/index.html





### Vibrioferrin

dir.create("3_Ann/Vibrioferrin")

#! the assembled reference dataset is provided as 'Suppl_file_3-uniprot-vibrioferrin_01mai20.faa'
#! UniProt query: "name:vibrioferrin length:[100 TO *]"
# parse ref db 
system2(command = "makeblastdb",
        args = c(paste("-in", VFds), 
                 "-dbtype prot"),
        stdout = T, stderr = T
)

# run blastp
#! adjust 'num_threads' to your available CPUs
system2(command = "blastp",
        args = c("-query 3_Ann/all_genes.faa",
                 paste("-db", VFds),
                 "-out 3_Ann/Vibrioferrin/Vibrioferrin_blastp.tsv",
                 "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'", 
                 "-evalue 1e-5", "-max_target_seqs 20", "-num_threads 7"),
        stdout = T, stderr = T
)





### DMSP degradation pathways

dir.create("3_Ann/DMSP_ddd")

#! the assembled reference dataset is provided as 'Suppl_file_4-uniprot-DMSP_19dec20.faa'
#! UniProt query: "gene:dmda OR gene:dmdb OR gene:dmdc OR gene:dmdd OR gene:dddl OR gene:dddy OR gene:dddd OR gene:dddq OR gene:dddw OR gene:dddp OR gene:dddk OR gene:alma1" # OR gene:ddda OR gene:dddc OR gene:acuk OR gene:acun AND reviewed:yes
# parse ref db 
system2(command = "makeblastdb",
        args = c(paste("-in", DMSPds), 
                 "-dbtype prot"),
        stdout = T, stderr = T
)

# run blastp
#! adjust 'num_threads' to your available CPUs
system2(command = "blastp",
        args = c("-query 3_Ann/all_genes.faa",
                 paste("-db", DMSPds),
                 "-out 3_Ann/DMSP_ddd/DMSP_blastp.tsv",
                 "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'", 
                 "-evalue 1e-5", "-max_target_seqs 20", "-num_threads 7"),
        stdout = T, stderr = T
)







### II parse annotation ----

### prokka
library(rtracklayer)
prok_list = list.files("3_Ann/Prokka", pattern = ".gff$", full.names = T, recursive = T)
prok_ann = lapply(prok_list, readGFF)
# to keep also tRNA and tmRNA use i$type != "gene"
prok_ann = lapply(prok_ann, function(i) as.data.frame(i[i$type == "CDS", c("ID","start","end","strand","gene","product","db_xref","eC_number")]))
names(prok_ann) = sapply(strsplit(prok_list, "/"), "[[", 3)
prok_ann = do.call(rbind, prok_ann)
prok_ann = data.frame(filename=gsub("\\..*", "", row.names(prok_ann)), prok_ann)

ann_master_TAB = prok_ann





### KEGG

## parse KOfam output
# KO_ann1 = read.delim("3_Ann/KEGG/kofam_result.txt", h=T, quote="", fill=F, check.names = F)
# KO_ann1 = strsplit(KO_ann1[2:nrow(KO_ann1), ], "[ ]+")
# KO_ann1 = lapply(KO_ann1, function(z) c(z[1:6], paste(z[7:length(z)], collapse = " ")))
# KO_ann1 = do.call(rbind, KO_ann1)
# colnames(KO_ann1) = c("sign","gene_name","KO","thrshld","score","E.value","KO.definition")
# KO_ann1 = as.data.frame(apply(KO_ann1, 2, unlist))
# write.table(KO_ann1, "3_Ann/KEGG/kofam_result.tsv", col.names = T, row.names = F, quote=F, sep="\t")
# KO_ann1 = read.delim("3_Ann/KEGG/kofam_result.tsv", h=T, quote="", fill=F)
# KO_ann2 = KO_ann1[grepl("\\*", KO_ann1$sign), ] # keep only signif
# write.table(KO_ann2, "3_Ann/KEGG/kofam_result_signif.tsv", col.names = T, row.names = F, quote=F, sep="\t")
KO_ann2 = read.delim("3_Ann/KEGG/kofam_result_signif.tsv", h=T, quote="", fill=F)
KO_ann2$gene_name = gsub(".*\\.", "", KO_ann2$gene_name)

## resolve promiscuity
library(dplyr)
KO_ann3 = KO_ann2 %>% group_by(gene_name) %>%
  # keep best match that is not an 'uncharacterized protein' if there is any
  filter( if (any(KO.definition != "uncharacterized protein")) { 
    KO.definition != "uncharacterized protein"
  } else KO.definition == KO.definition ) %>%
  filter(score == max(score)) %>%
  # keep match with highest threshold
  filter( if (length(thrshld) > 1) { 
    thrshld == max(thrshld)
  } else thrshld == thrshld )
#---

ann_master_TAB$KO = KO_ann3$KO[match(ann_master_TAB$ID, KO_ann3$gene_name)]





### BioV-transporters
library(stringr)
trans_list = list.files(path = Sys.glob("3_Ann/Transporter/*", dirmark = F), full.names = T, pattern = "results.tsv")
trans_ann1 = lapply(trans_list, function(i) read.delim(i, h=T, sep = "\t"))
trans_ann2 = do.call(rbind, trans_ann1)
e_val = 1e-6 # specify e-value threshold
TMP_over = 1 # min number of trans-membrane alpha-helical counts
trans_ann2 = trans_ann2[trans_ann2$TM_Overlap_Score >= TMP_over & trans_ann2$e.value <= e_val, 
               c("X.Query_id","Hit_tcid","Family_Abrv","Hit_desc","Predicted_Substrate")]
trans_ann3 = data.frame(trans_ann2[, 1:3],
                        accession=sapply(strsplit(trans_ann2$Hit_desc, " "), function(x) paste(x[c(1:2)], collapse = " ")), 
                        Hit_desc=sapply(strsplit(trans_ann2$Hit_desc, " "), function(x) paste(x[-c(1:2)], collapse = " ")),
                        Predicted_Substrate=gsub("CHEBI:[0-9]+;", "", trans_ann2[, 5]),
                        CHEBI=sapply(str_extract_all(trans_ann2[, 5], "CHEBI:[0-9]+"), paste0, collapse=";"))

## Select for relevant ones
trans_ann4 = trans_ann3[sapply(trans_ann3$CHEBI, function(i) {
  any(strsplit(i, ";")[[1]] %in% 
        c("CHEBI:15956", "CHEBI:41236", "CHEBI:13905", "CHEBI:22882", "CHEBI:22884", "CHEBI:3108" #B7
          , "CHEBI:30411" #B12
          , "CHEBI:17439", "CHEBI:60496", "CHEBI:48820", "CHEBI:3979", "CHEBI:14041", "CHEBI:23435" #CYANO-B12
          , "CHEBI:17154", "CHEBI:44258", "CHEBI:7556", "CHEBI:14645", "CHEBI:25521" #B3-nicotinamide
          , "CHEBI:15940", "CHEBI:44319", "CHEBI:7559", "CHEBI:25538" #B3-nicotinic acid
          , "CHEBI:9532" #B1-2P
          # , "CHEBI:9533" #B1-1P
          # , "CHEBI:18385", "CHEBI:46393", "CHEBI:9530", "CHEBI:15227", "CHEBI:26941" #B1-noP
        ))}), ]

## drop drug efflux pumps
bad_transp = unique(trans_ann4$Predicted_Substrate)[c(1,4,10,13,14,15,16)] #!!! manual inspect for suspicious transporters and edit
trans_ann4 = trans_ann4[!trans_ann4$Predicted_Substrate %in% bad_transp, ]

# ## if you wanna keep all transporters
# trans_ann4 = trans_ann3

ann_master_TAB$TCdb = trans_ann4$Hit_tcid[match(ann_master_TAB$ID, trans_ann4$X.Query_id)]





### Secondary Metabolites
#!!! code for AntiSMASH v5
library(rvest)
sm_list = list.files(path = Sys.glob("3_Ann/Antismash/*", dirmark = F), full.names = T, pattern = "index.html")
sec_metab = lapply(sm_list, function(i) tryCatch(read_html(i), error = function(x) cat("No secondary metabolites ", i,"\n")) )
sec_metab2 = lapply(sec_metab, function(i) html_table(html_nodes(i, "table"), fill = T))
sec_metab3 = lapply(seq_along(sec_metab2), function(i) {
  i2 = sec_metab2[[i]][[length(sec_metab2[[i]])]] # select the last table of the html page (is of all other tables that corresponds to different genomic regions)
  i3 = data.frame(filename=strsplit(sm_list[i], "/")[[1]][length(strsplit(sm_list[i], "/")[[1]])-1], i2) 
  return(i3) }) 
sec_metab4 = do.call(rbind, sec_metab3)
sec_metab4 = sec_metab4[!grepl("unknown|,", sec_metab4$Type), ]
sec_metab4$From = as.numeric(gsub(",", "", sec_metab4$From))
sec_metab4$To = as.numeric(gsub(",", "", sec_metab4$To))

## select relevant ones
sm_list = as.data.frame(readxl::read_xlsx(suppl_tabs, col_names = T, sheet = "S_tab_6", skip = 3)) 
sm_list = sm_list[!is.na(sm_list$`Interaction traits`), ] # remove extra lines in Excel table
sec_metab4 = sec_metab4[sec_metab4$Type %in% sm_list$ID[sm_list$`Interaction traits` != "-"], ]

ann_master_TAB$SecMetab = NA
invisible(apply(sec_metab4, 1, function(j) {
  ann_master_TAB$SecMetab[ann_master_TAB$filename == j[1] & 
                            (ann_master_TAB$start <= as.integer(j[4]) | ann_master_TAB$start < as.integer(j[5])) &
                            (ann_master_TAB$end > as.integer(j[4]) | ann_master_TAB$end >= as.integer(j[5]))
  ] <<- j[3]
}))





### Vibrioferrin
library(stringr)
VF_ann1 = read.delim("3_Ann/Vibrioferrin/Vibrioferrin_blastp.tsv", h=F)
VF_ann1$gene_ID = gsub(".*\\.","",VF_ann1$V1)
VF_ann2 = VF_ann1
VF_ann2 = VF_ann2[VF_ann2$V11 < 1e-6, ] # filter based on min e-value - Gärdes et al., 2013
# VF_ann2 = VF_ann2[VF_ann2$V12 > 60, ] # filter based on min bit score

## get coherent gene names
VF_ann2$gene_name = sapply(VF_ann2$V13, function(i) {
  ifelse(grepl("pvu", i, ignore.case = T), paste0("pvu", substr(gsub(".*pvu", "", i, ignore.case = T), 1, 1)),
         ifelse(grepl("pvs", i, ignore.case = T), paste0("pvs", substr(gsub(".*pvs", "", i, ignore.case = T), 1, 1)), NA
                ))
}) 

## parse blast output
VF_ann3 = do.call(rbind, lapply(split(VF_ann2, f = VF_ann2$V1), function(i) {
  #! find best hit
  i2 = i[i[,11] == min(i[,11]) & i[,12] == max(i[,12]), , drop=F]
  #! solve promiscuity (i.e. label as unclear if ties have different gene names)
  if (nrow(i2) > 1) {
    i2 = i2[1, , drop=F]
    if (length(unique(i2[, 15])) == 1) {
      i2[1, 15] = unique(i2[, 15])
    } else { i2[1, 15] = "unclear_ann" }
  }
  return(i2)
}))

## check for completeness
library(reshape2)
VF_ann3$gnm = gsub("_.*", "", VF_ann3$gene_ID)
VF_ann4 = dcast(VF_ann3, gnm~gene_name)
#!!! completeness: require at least 4 of 5 genes for both operons
VF_ann5 = data.frame(gnm=VF_ann4$gnm,
                     pvsABCDE=apply(VF_ann4[, grepl("pvs", colnames(VF_ann4))], 1, function(x) ifelse(sum(x>0) >= 4, 1, 0)),
                     pvuABCDE=apply(VF_ann4[, grepl("pvu", colnames(VF_ann4))], 1, function(x) ifelse(sum(x>0) >= 4, 1, 0)))
VF_ann6 = VF_ann3[(VF_ann3$gnm %in% VF_ann5$gnm[VF_ann5$pvsABCDE > 0] &
                     grepl("pvs", VF_ann3$gene_name)) |
                    (VF_ann3$gnm %in% VF_ann5$gnm[VF_ann5$pvuABCDE > 0] 
                     & grepl("pvu", VF_ann3$gene_name)), ]

ann_master_TAB$Vibrioferrin = VF_ann6$gene_name[match(ann_master_TAB$ID, VF_ann6$gene_ID)]
  




### Phytohormones
#!!! based on the KEGG annotation 

ph_list = as.data.frame(readxl::read_xlsx(suppl_tabs, col_names = T, sheet = "S_tab_7", skip = 3)) 
ph_list = ph_list[!is.na(ph_list$`Interaction traits`), ] # remove extra lines in Excel table

## preare a list of KO-path for each hormone
ph_list = split(ph_list$`KEGG Orthology`, f = ph_list$`Interaction traits`)
ph_list = lapply(ph_list, function(i) strsplit(gsub("\\n| ", "", i), "\\|"))
ph_list = lapply(seq_along(ph_list), function(i) data.frame(ph=names(ph_list)[i], KOs=do.call(cbind, ph_list[[i]])))
ph_list = do.call(rbind, ph_list)
ph_list$length = sapply(ph_list$KOs, function(i) length(strsplit(i, ",")[[1]]))

## check for any path in genomes using the annotated KOs
ph_ann1 = lapply(ph_list$KOs, function(i) {
  i2 = strsplit(i, ",")[[1]] 
  with(ann_master_TAB[ann_master_TAB$KO %in% i2, ], table(filename, KO))
})
ph_ann2 = lapply(seq_along(ph_ann1), function(i) {
  i2 = apply(ph_ann1[[i]], 1, function(j) sum(j>0))
  if (length(i2) > 0) {
    data.frame(ph_list[i,], filename=names(i2), hits=as.integer(i2))
  } # else {   data.frame(ph_list[i,], filename=NA, hits=NA) }
})
ph_ann3 = do.call(rbind, ph_ann2)
  
## check for completeness
len_breaks = NULL # set to NULL if no gap is allowed (i.e. allowed_gaps = 0)
allowed_gaps = c(0)
#! estimate gaps allowed for each path
if (!is.null(len_breaks)) {
  len_bin <- .bincode(ph_ann3$length, breaks = c(0, len_breaks, Inf), right = F, include.lowest = T)
} else {
  len_bin <- .bincode(ph_ann3$length, breaks = c(0, Inf), right = F, include.lowest = T)
}
allowed_gaps <- allowed_gaps[len_bin]
#! add 'allowed_gaps' as bonus
ph_ann4 = ph_ann3
ph_ann4$hits = ph_ann3$hits + allowed_gaps
#! filter for completeness
ph_ann5 = ph_ann4
ph_ann5$hits = ph_ann5$hits/ph_ann5$length
ph_ann5$hits[ph_ann5$hits >=1] = 1
ph_ann5$hits[ph_ann5$hits <1] = 0

## keep only complete ph path in each gnm, and the unique list of KOs
ph_ann5 = ph_ann5[ph_ann5$hits > 0, ]
ph_ann6 = ph_ann5 %>%
  group_by(ph,filename) %>%
  summarise(KOs=unique(unlist(strsplit(KOs, ","))))

ann_master_TAB$Phytohormones = NA
invisible(apply(ph_ann6, 1, function(j) {
  ann_master_TAB$Phytohormones[ann_master_TAB$filename %in% j[2] &
    ann_master_TAB$KO %in% j[3]
  ] <<- j[1]
}))





### DMSP
library(stringr)
DMSP_ann1 = read.delim("3_Ann/DMSP_ddd/DMSP_blastp.tsv", h=F)
DMSP_ann1$gene_ID = gsub(".*\\.","",DMSP_ann1$V1)
DMSP_ann1$filename = gsub("_.*", "", DMSP_ann1$gene_ID)
DMSP_ann2 = DMSP_ann1
DMSP_ann2 = DMSP_ann2[DMSP_ann2$V11 < 1e-70, ] # filter based on min e-value - Gärdes et al., 2013
# DMSP_ann2 = DMSP_ann2[DMSP_ann2$V12 > 60, ] # filter based on min bit score
  
## get coherent gene names
DMSP_ann2$gene_name = sapply(DMSP_ann2$V13, function(i) {
  ifelse(grepl("ddd", i, ignore.case = T), paste0("ddd", toupper(substr(gsub(".* ddd|.*=ddd", "", i, ignore.case = T), 1, 1))),
         ifelse(grepl("dmd", i, ignore.case = T), paste0("dmd", toupper(substr(gsub(".* dmd|.*=dmd", "", i, ignore.case = T), 1, 1))),
                ifelse(grepl("acu", i, ignore.case = T), paste0("acu", toupper(substr(gsub(".* acu|.*=acu", "", i, ignore.case = T), 1, 1))), "NA"
                )))
})

## parse blast output
DMSP_ann3 = do.call(rbind, lapply(split(DMSP_ann2, f = DMSP_ann2$V1), function(i) {
  #! find best hit
  i2 = i[i[,11] == min(i[,11]) & i[,12] == max(i[,12]), , drop=F]
  #! solve promiscuity (i.e. label as unclear if ties have different gene names)
  if (nrow(i2) > 1) {
    i2 = i2[1, , drop=F]
    if (length(unique(i2[, 15])) == 1) {
      i2[1, 15] = unique(i2[, 15])
    } else { i2[1, 15] = "unclear_ann" }
  }
  return(i2)
}))

## check for completeness
DMSP_ann4 = table(DMSP_ann3$filename, DMSP_ann3$gene_name)

DMSP_ann4a = DMSP_ann4[, grepl("ddd|ALMA", colnames(DMSP_ann4))]
DMSP_ann4a = apply(DMSP_ann4a, 1, function(j) sum(j>0))
DMSP_ann4a = names(DMSP_ann4a)[DMSP_ann4a > 0] #!!! completeness: require at least 1 gene (are all redundant)

DMSP_ann4b = DMSP_ann4[, grepl("dmd", colnames(DMSP_ann4))]
DMSP_ann4b = apply(DMSP_ann4b, 1, function(j) sum(j>0))
DMSP_ann4b = names(DMSP_ann4b)[DMSP_ann4b >= 3] #!!! completeness: require at least 3 of 4 genes

ann_master_TAB$DMSP = NA
ann_master_TAB$DMSP[match(DMSP_ann3$gene_ID[DMSP_ann3$filename %in% DMSP_ann4a], ann_master_TAB$ID)] = "DMSP_cleavage"
ann_master_TAB$DMSP[match(DMSP_ann3$gene_ID[DMSP_ann3$filename %in% DMSP_ann4b], ann_master_TAB$ID)] = "DMSP_demethylation"





### DHPS 
#!!! based on the KEGG annotation 

## seek for the key KO:
## "K15509" (sulfopropanediol 3-dehydrogenase, gene hpsN)
dhps_ann1 = ann_master_TAB[ann_master_TAB$KO %in% "K15509", ]

ann_master_TAB$DHPS = NA
invisible(apply(dhps_ann1, 1, function(j) {
  ann_master_TAB$DHPS[ann_master_TAB$filename %in% j[1] &
                        ann_master_TAB$KO %in% j[10]
  ] <<- "DHPS_catabolism"
}))





### Taurine 
#!!! based on the KEGG annotation 

## seek for the key KOs:
## 1) "K07255" + "K07256" (taurine dehydrogenase, tauXY genes) & K03852 (sulfoacetaldehyde acetyltransferase, xsc gene)
tau_ann1a = with(ann_master_TAB[ann_master_TAB$KO %in% c("K03851","K03852"), ], table(filename, KO))
tau_ann2a = apply(tau_ann1a, 1, function(j) sum(j>0))
tau_ann3a = names(tau_ann2a)[tau_ann2a > 1] #!!! completeness: require all genes (2)

ann_master_TAB$Taurine[ann_master_TAB$filename %in% tau_ann3a &
                         ann_master_TAB$KO %in% c("K03851","K03852")] = "Taurine_in_TCA"


## 2) "K03851" (taurine-pyruvate aminotransferase, tpa gene) & K03852 (sulfoacetaldehyde acetyltransferase, xsc gene)
tau_ann1b = with(ann_master_TAB[ann_master_TAB$KO %in% c("K07255","K07256","K03852"), ], table(filename, KO))
tau_ann2b = apply(tau_ann1b, 1, function(j) sum(j>0))
tau_ann3b = names(tau_ann2b)[tau_ann2b > 2] #!!! completeness: require all genes (3)

ann_master_TAB$Taurine[ann_master_TAB$filename %in% tau_ann3b &
                         ann_master_TAB$KO %in% c("K07255","K07256","K03852")] = "Taurine_in_TCA"







### III Format datasets & Interaction traits --------------------------------------------------

meta = as.data.frame(readxl::read_xlsx(suppl_tabs, col_names = T, sheet = "S_tab_2", skip = 3)) 
meta = meta[!is.na(meta$Filename), ] # remove extra lines in Excel table


### reconstruct KM
library(reshape2)
source(KM_code)

# # download new KM structures
# KM_str = KMdiagram_fetcher(ncore = 7, create_RData = T, path = "~/Downloads") # will save the fetched file in the current wd
# load KM structures used in the paper
KM_str = readRDS(KM_str)

ann_master_TAB$cnt = 1
KO_ann = dcast(ann_master_TAB, filename~KO, value.var = "cnt")
row.names(KO_ann) = KO_ann$filename
KO_ann = KO_ann[, -1]
KO_ann[KO_ann > 1] = 1
KM_ann = KMreco(indata = KO_ann, KM_str = KM_str, 
                len_breaks = c(3), allowed_gaps = c(0,1))


### create trait-table
fun_profi = cbind(KM_ann$KM_pa
                  , 
                  )






ITs_mapping = rbind(ITs_mapping, data.frame(col_lab=colnames(tr_VF), trait=c("Siderophore","Siderophore-trans"),
                                             text=c("Vibrioferrin biosynthesis","Vibrioferrin transporter"), db="VF"))





##### prepare input datasets
nonITs_mapping = data.frame(col_lab=colnames(KM_ann)[!colnames(KM_ann) %in% ITs_mapping$col_lab],
                             trait=NA,
                             text=gsub(" \\[.*", "", KM_hier$X4[match(colnames(KM_ann)[!colnames(KM_ann) %in% ITs_mapping$col_lab], gsub(" .*", "", KM_hier$X4))]),
                             db="KEGG_modules")
ITs_mapping = rbind(ITs_mapping, nonITs_mapping)
ITs_mapping = cbind(ITs_mapping, KM_hier[match(ITs_mapping$col_lab, gsub(" .*", "", KM_hier$X4)), c("X2","X3")])
# write.table(ITs_mapping, "trait_list.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
fun_profi = cbind(KM_ann, tr_PH, tr_TonB, tr_SM, tr_VF)
# write.table(data.frame(ID=rownames(fun_profi), fun_profi), "fun_profi_all.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
ITs_mapping = ITs_mapping[match(colnames(fun_profi), ITs_mapping$col_lab), ]
ITs_mapping$RA = colSums(fun_profi)/nrow(fun_profi)







