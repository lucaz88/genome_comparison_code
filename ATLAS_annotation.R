## set working directory
# setwd("/media/lucaz/Bioinfo/HB_complete_genomes_v4/")
setwd("/media/lucaz/DATA/HB_complete_genomes_v4/")


## global setting
options(stringsAsFactors = F)
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8") # to get month names in english


##!!! following tools must be correctly installed (check required dependencies) in the system
#a prokka (https://github.com/tseemann/prokka)
#b kofamscan (https://github.com/takaram/kofam_scan)
#c gblast (http://www.tcdb.org/files/gblast.zip)
#d Antismash (https://antismash.secondarymetabolites.org/#!/download)
#e blast+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)


##!!! required files
# Supplementary tables' file
suppl_tabs = "/media/lucaz/DATA/Google_drive/Laboratorio/AAA_Progetti_NUOVI/2016_HFSP/4_paper_genome_comparison/6_REV_Comm_biology/Suppl_tables.xlsx" #! adjust the path
# Supplementary file 3
PHds = "/media/lucaz/DATA/Google_drive/Laboratorio/AAA_Progetti_NUOVI/2016_HFSP/4_paper_genome_comparison/6_REV_Comm_biology/Suppl_file_3-phytohormone-KOs.fasta" #! adjust the path
# Supplementary file 4 
#! UniProt query: "name:vibrioferrin length:[100 TO *]"
VFds = "/media/lucaz/DATA/Google_drive/Laboratorio/AAA_Progetti_NUOVI/2016_HFSP/4_paper_genome_comparison/6_REV_Comm_biology/Suppl_file_4-uniprot-vibrioferrin_01mai20.fasta" #! adjust the path
# Supplementary file 5 
#! UniProt query: "gene:dmda OR gene:dmdb OR gene:dmdc OR gene:dmdd OR gene:dddl OR gene:dddy OR gene:dddd OR gene:dddq OR gene:dddw OR gene:dddp OR gene:dddk OR gene:alma1" # OR gene:ddda OR gene:dddc OR gene:acuk OR gene:acun AND reviewed:yes
DMSPds = "/media/lucaz/DATA/Google_drive/Laboratorio/AAA_Progetti_NUOVI/2016_HFSP/4_paper_genome_comparison/6_REV_Comm_biology/Suppl_file_5-uniprot-DMSP_21dec20.fasta" #! adjust the path
# KM_reconstruction.R code
KM_code = "/home/lucaz/myscript/GitHub/gnm_compar/KM_reconstruction.R" 
# file with KM structures (optional, it can be downloaded anew) 
KMds = "/home/lucaz/myscript/GitHub/gnm_compar/KM_str_2020-01-09_plus.rds"


load("Atlas_ann_v1.RData")
#



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
# conda activate BioVx; ls 3_Ann/Prokka/*/*.faa | parallel -j 4 'gblast3.py -i {} -o 3_Ann/Transporter/$(cut -d'/' -f3 <<< {})'

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





### Phytohormones

dir.create("3_Ann/Phytohormones")

# prepare ref db
system2(command = "makeblastdb",
        args = c(paste("-in", PHds), 
                 "-dbtype prot"),
        stdout = T, stderr = T
)

# run blastp
#! adjust 'num_threads' to your available CPUs
system2(command = "blastp",
        args = c("-query 3_Ann/all_genes.faa",
                 paste("-db", PHds),
                 "-out 3_Ann/Phytohormones/Phytohormones_blastp.tsv",
                 "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'", 
                 "-evalue 1e-5", "-max_target_seqs 20", "-num_threads 7"),
        stdout = T, stderr = T
)





### Vibrioferrin

dir.create("3_Ann/Vibrioferrin")

# prepare ref db
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

# prepare ref db
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

#! load packages
library(rtracklayer)
#! set variables
gff_list = list.files("3_Ann/Prokka", pattern = ".gff$", full.names = T, recursive = T)

prok_ann = lapply(gff_list, readGFF)
# to keep also tRNA and tmRNA use i$type != "gene"
prok_ann = lapply(prok_ann, function(i) as.data.frame(i[i$type == "CDS", c("ID","start","end","strand","gene","product","db_xref","eC_number")]))
names(prok_ann) = sapply(strsplit(gff_list, "/"), "[[", 3)
prok_ann = do.call(rbind, prok_ann)
prok_ann = data.frame(filename=gsub("\\..*", "", row.names(prok_ann)), prok_ann)

ann_master_TAB = prok_ann





### KEGG

#! load packages
library(dplyr)
#! set variables
KO_ann1 = read.delim("3_Ann/KEGG/kofam_result.txt", h=T, quote="", fill=F, check.names = F)

KO_ann1 = strsplit(KO_ann1[2:nrow(KO_ann1), ], "[ ]+")
KO_ann1 = lapply(KO_ann1, function(z) c(z[1:6], paste(z[7:length(z)], collapse = " ")))
KO_ann1 = do.call(rbind, KO_ann1)
colnames(KO_ann1) = c("sign","gene_name","KO","thrshld","score","E.value","KO.definition")
KO_ann1 = as.data.frame(apply(KO_ann1, 2, unlist))
write.table(KO_ann1, "3_Ann/KEGG/kofam_result.tsv", col.names = T, row.names = F, quote=F, sep="\t")
KO_ann1 = read.delim("3_Ann/KEGG/kofam_result.tsv", h=T, quote="", fill=F)
KO_ann2 = KO_ann1[grepl("\\*", KO_ann1$sign), ] # keep only signif
write.table(KO_ann2, "3_Ann/KEGG/kofam_result_signif.tsv", col.names = T, row.names = F, quote=F, sep="\t")
# KO_ann2 = read.delim("3_Ann/KEGG/kofam_result_signif.tsv", h=T, quote="", fill=F)
KO_ann2$gene_name = gsub(".*\\.", "", KO_ann2$gene_name)

#! resolve promiscuity
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

KO_ann = KO_ann3
ann_master_TAB$KO = KO_ann$KO[match(ann_master_TAB$ID, KO_ann$gene_name)]
remove(list = ls(pattern='^KO_ann[0-9]+'))




### BioV-transporters

#! load packages
library(stringr)
#! set variables
trans_list = list.files(path = Sys.glob("3_Ann/Transporter/*", dirmark = F), full.names = T, pattern = "results.tsv")

transp_ann1 = lapply(trans_list, function(i) read.delim(i, h=T, sep = "\t"))
transp_ann1 = do.call(rbind, transp_ann1)
e_val = 1e-6 # specify e-value threshold
TMP_over = 1 # min number of trans-membrane alpha-helical counts
transp_ann2 = transp_ann1[transp_ann1$TM_Overlap_Score >= TMP_over &
                            transp_ann1$TM_Overlap_Score != "None" &
                            transp_ann1$e.value <= e_val, ]
transp_ann3 = data.frame(transp_ann2[, c("X.Query_id","Hit_tcid","Family_Abrv")],
                        accession=sapply(strsplit(transp_ann2$Hit_desc, " "), function(x) paste(x[c(1:2)], collapse = " ")), 
                        Hit_desc=sapply(strsplit(transp_ann2$Hit_desc, " "), function(x) paste(x[-c(1:2)], collapse = " ")),
                        Predicted_Substrate=gsub("CHEBI:[0-9]+;", "", transp_ann2[, "Predicted_Substrate"]),
                        CHEBI=sapply(str_extract_all(transp_ann2[, "Predicted_Substrate"], "CHEBI:[0-9]+"), paste0, collapse=";"))

# #! Select B-vitamin transporters
# transp_ann4 = transp_ann3[sapply(transp_ann3$CHEBI, function(i) {
#   any(strsplit(i, ";")[[1]] %in% 
#         c("CHEBI:15956", "CHEBI:41236", "CHEBI:13905", "CHEBI:22882", "CHEBI:22884", "CHEBI:3108" #B7
#           , "CHEBI:30411" #B12
#           , "CHEBI:17439", "CHEBI:60496", "CHEBI:48820", "CHEBI:3979", "CHEBI:14041", "CHEBI:23435" #CYANO-B12
#           , "CHEBI:17154", "CHEBI:44258", "CHEBI:7556", "CHEBI:14645", "CHEBI:25521" #B3-nicotinamide
#           , "CHEBI:15940", "CHEBI:44319", "CHEBI:7559", "CHEBI:25538" #B3-nicotinic acid
#           , "CHEBI:9532" #B1-2P
#           # , "CHEBI:9533" #B1-1P
#           # , "CHEBI:18385", "CHEBI:46393", "CHEBI:9530", "CHEBI:15227", "CHEBI:26941" #B1-noP
#         ))}), ]
# 
# #! drop drug efflux pumps
# bad_transp = unique(transp_ann4$Predicted_Substrate)[c(1,4,10,13,14,15,16)] #!!! manual inspect for suspicious transporters and edit
# transp_ann4 = transp_ann4[!transp_ann4$Predicted_Substrate %in% bad_transp, ]

#! remove generic transporter (e.g. for ions, metals and without a predicted substrate)
bad_cpd = c("None","hydron","ion","cation","anion","proton","sodium(1+)",
            "potassium(1+)","calcium(2+)","lithium(1+)","chloride","electron",
            "fluoride","barium(2+)","strontium(2+)","cadmium(2+)","lead(2+)",
            "cobalt(2+)","mercury(2+)","nickel(2+)","zinc(2+)","dioxygen",
            "carbon dioxide","rubidium(1+)","iron(3+)","iron(2+)","manganese(2+)",
            "magnesium(2+)","metal cation","molecule","copper(1+)","zinc ion",
            "water","metabolite","dicarboxylic acid dianion","copper cation",
            "surfactant","vanadium oxoanion","monocarboxylic acid anion",
            "calcium ion","inorganic anion","base","inorganic cation",
            "monoatomic monocation","sodium atom","lithium atom","copper(2+)",
            "chromate(2-)","silver(1+)","selenite(2-)","detergent",
            "hydrogencarbonate","arsenite(3-)","arsenate(3-)",
            "arsenic molecular entity","benzalkonium chloride",
            "sodium dodecyl sulfate","CCCP","antimonite","aluminium(3+)",
            "sodium tungstate","molybdate","bromide","sodium iodide",
            "chlorate","bromate","periodate","thiocyanate","tetrafluoroborate(1-)",
            "nitric acid","silicic acid","sodium chloride","potassium chloride",
            "tungstate","silicate(4-)","silicon atom","caesium(1+)","mercury(0)")
transp_ann4 = transp_ann3[!sapply(transp_ann3$Predicted_Substrate, function(z) all(strsplit(z, ", ")[[1]] %in% bad_cpd)), ]

# #! if you wanna keep all transporters
# transp_ann4 = transp_ann3

transp_ann = transp_ann4
ann_master_TAB$TCdb = transp_ann$Hit_tcid[match(ann_master_TAB$ID, transp_ann$X.Query_id)]
remove(list = ls(pattern='^transp_ann[0-9]+'))




### Secondary Metabolites
#!!! code for AntiSMASH v5

#! load packages
library(rvest)
library(rtracklayer)
#! set variables
sm_list = list.files(path = Sys.glob("3_Ann/Antismash/*", dirmark = F), full.names = T, pattern = "index.html")
gff_list = list.files("3_Ann/Prokka", pattern = ".gff$", full.names = T, recursive = T)
gff_list = gff_list[match(sapply(strsplit(sm_list, "/"), function(i) rev(i)[2]), 
                          sapply(strsplit(gff_list, "/"), function(i) rev(i)[2]))]

sm_ann1 = lapply(sm_list, function(i) tryCatch(read_html(i), error = function(x) cat("No secondary metabolites ", i,"\n")) )
sm_ann2 = lapply(sm_ann1, function(i) {
  gkb_ann = html_table(html_nodes(i, "table"), fill = T)
  gkb_ann = gkb_ann[[length(gkb_ann)]] # last table include all the others
  gff_locus = html_text(html_nodes(i, "strong"))
  gff_locus = gff_locus[1:nrow(gkb_ann)] # keep only region with annotations
  cbind(gff_locus, gkb_ann)
})
sm_ann3 = lapply(seq_along(sm_ann2), function(i) { # add filename
  data.frame(filename=strsplit(sm_list[i], "/")[[1]][length(strsplit(sm_list[i], "/")[[1]])-1], sm_ann2[[i]]) 
})
sm_ann4 = lapply(seq_along(sm_ann3), function(i) { # map annotations to CDS
  i_sm = sm_ann3[[i]]
  i_sm$From = as.integer(gsub(",", "", i_sm$From))
  i_sm$To = as.integer(gsub(",", "", i_sm$To))
  i_gff = as.data.frame(readGFF(gff_list[i]))
  i_gff = i_gff[i_gff$type == "CDS", ]
  i_gff$seqid = sapply(strsplit(as.character(i_gff$seqid), "\\|"), function(j) rev(j)[[1]])
  i_sm = do.call(rbind, lapply(seq_len(nrow(i_sm)), function(j) {
    j = i_sm[j, , drop=T]
    i_gff_j = i_gff[i_gff$seqid == j["gff_locus"] &
                      (i_gff$start <= j["From"] | i_gff$start < j["To"]) &
                      (i_gff$end > j["From"] | i_gff$end >= j["To"]), c("seqid","ID")]
    merge(i_gff_j, j, by.x="seqid", by.y="gff_locus", all=T)
  }))
})
sm_ann4 = do.call(rbind, sm_ann4)
sm_ann4 = sm_ann4[!grepl("unknown|,", sm_ann4$Type), ]

# #! select only relevant ones
# sm_list = as.data.frame(readxl::read_xlsx(suppl_tabs, col_names = T, sheet = "S_tab_6", skip = 3)) 
# sm_list = sm_list[!is.na(sm_list$`Interaction traits`), ] # remove extra lines in Excel table
# sm_ann4 = sm_ann4[sm_ann4$Type %in% sm_list$ID[sm_list$`Interaction traits` != "-"], ]

sm_ann = sm_ann4
ann_master_TAB$SecMetab = sm_ann$Type[match(ann_master_TAB$ID, sm_ann$ID)]
remove(list = ls(pattern='^sm_ann[0-9]+|sm_list'))





### Phytohormones

#! load packages
library(dplyr)
#! set variables
ph_ann1 = read.delim("3_Ann/Phytohormones/Phytohormones_blastp.tsv", h=F)

ph_ann1$gene_ID = gsub(".*\\.","",ph_ann1$V1)
ph_ann1$gnm = gsub("_.*", "", ph_ann1$gene_ID)
ph_ann1$KO = gsub("\\|.*", "", ph_ann1$V2)
ph_ann1$KO = gsub("K01502,", "", ph_ann1$KO) # K01502 has same gene sequence as K01501, but K01502 is not in KEGG map00380

## parse blast output
ph_ann2 = do.call(rbind, lapply(split(ph_ann1, f = ph_ann1$V1), function(i) {
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
ph_ann2 = ph_ann2[ph_ann2$V11 < 5, ] # Zhang et al., 2019
ph_ann2 = ph_ann2[ph_ann2$V12 > 60, ] # Zhang et al., 2019

## check for any path in genomes using the annotated KOs
ph_list = as.data.frame(readxl::read_xlsx(suppl_tabs, col_names = T, sheet = "S_tab_7", skip = 3)) 
ph_list = ph_list[!is.na(ph_list$`Interaction traits`), ] # remove extra lines in Excel table
ph_list = split(ph_list$`KEGG Orthology`, f = ph_list$`Interaction traits`)
ph_list = lapply(ph_list, function(i) strsplit(gsub("\\n| ", "", i), "\\|"))
ph_list = lapply(seq_along(ph_list), function(i) data.frame(ph=names(ph_list)[i], KO=do.call(cbind, ph_list[[i]])))
ph_list = do.call(rbind, ph_list)
ph_list$len = sapply(ph_list$KO, function(i) length(strsplit(i, ",")[[1]]))

ph_ann3 = lapply(ph_list$KO, function(i) {
  i2 = strsplit(i, ",")[[1]] 
  with(ph_ann2[ph_ann2$KO %in% i2, ], table(gnm, KO))
})
ph_ann3 = lapply(seq_along(ph_ann3), function(i) {
  i2 = apply(ph_ann3[[i]], 1, function(j) sum(j>0))
  if (length(i2) > 0) {
    data.frame(ph_list[i,], gnm=names(i2), ann_len=as.integer(i2))
  } # else {   data.frame(ph_list[i,], gnm=NA, ann_len=NA) }
})
ph_ann3 = do.call(rbind, ph_ann3)

## check for completeness (use same rules as for KMs)
len_breaks = c(3) # set to NULL if no gap is allowed (i.e. allowed_gaps = 0)
allowed_gaps = c(0,1)
#! estimate gaps allowed for each path
if (!is.null(len_breaks)) {
  len_bin <- .bincode(ph_ann3$ann_len, breaks = c(0, len_breaks, Inf), right = F, include.lowest = T)
} else {
  len_bin <- .bincode(ph_ann3$ann_len, breaks = c(0, Inf), right = F, include.lowest = T)
}
allowed_gaps <- allowed_gaps[len_bin]
#! add 'allowed_gaps' as bonus
ph_ann4 = ph_ann3
ph_ann4$ann_len = ph_ann3$ann_len + allowed_gaps
#! filter for ann_leneteness
ph_ann5 = ph_ann4
ph_ann5$compl = ph_ann5$ann_len/ph_ann5$len
ph_ann5$compl[ph_ann5$compl >=1] = 1
ph_ann5$compl[ph_ann5$compl <1] = 0
ph_ann5 = ph_ann5[ph_ann5$compl > 0, ]

## reformat table
ph_ann6 = ph_ann5 %>%
  group_by(gnm, ph) %>%
  summarise(KO=unique(unlist(strsplit(KO, ",")))) %>%
  group_by(gnm, KO) %>%
  summarise(ph=paste0(unique(ph), collapse = ";"))
ph_ann7 = merge(ph_ann2, ph_ann6, by=c("gnm","KO"), all.x=T)
  
ph_ann = ph_ann7
ann_master_TAB$Phytohormones = ph_ann$ph[match(ann_master_TAB$ID, ph_ann$gene_ID)]
remove(list = ls(pattern='^ph_ann[0-9]+|ph_list'))





### Vibrioferrin

#! load packages
library(stringr)
library(reshape2)
#! set variables
vf_ann1 = read.delim("3_Ann/Vibrioferrin/Vibrioferrin_blastp.tsv", h=F)

vf_ann1$gene_ID = gsub(".*\\.","",vf_ann1$V1)
vf_ann2 = vf_ann1
vf_ann2 = vf_ann2[vf_ann2$V11 < 1e-5, ] # filter based on min e-value
# vf_ann2 = vf_ann2[vf_ann2$V12 > 60, ] # filter based on min bit score

#! get coherent gene names
vf_ann2$gene_name = sapply(vf_ann2$V13, function(i) {
  ifelse(grepl("pvu", i, ignore.case = T), paste0("pvu", substr(gsub(".*pvu", "", i, ignore.case = T), 1, 1)),
         ifelse(grepl("pvs", i, ignore.case = T), paste0("pvs", substr(gsub(".*pvs", "", i, ignore.case = T), 1, 1)), NA
                ))
}) 

#! parse blast output
vf_ann3 = do.call(rbind, lapply(split(vf_ann2, f = vf_ann2$V1), function(i) {
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

#! check for completeness
vf_ann3$gnm = gsub("_.*", "", vf_ann3$gene_ID)
vf_ann4 = dcast(vf_ann3, gnm~gene_name)
#! completeness: require at least 4 of 5 genes for both operons
vf_ann5 = data.frame(gnm=vf_ann4$gnm,
                     pvsABCDE=apply(vf_ann4[, grepl("pvs", colnames(vf_ann4))], 1, function(x) ifelse(sum(x>0) >= 4, 1, 0)),
                     pvuABCDE=apply(vf_ann4[, grepl("pvu", colnames(vf_ann4))], 1, function(x) ifelse(sum(x>0) >= 4, 1, 0)))
vf_ann6 = vf_ann3[(vf_ann3$gnm %in% vf_ann5$gnm[vf_ann5$pvsABCDE > 0] &
                     grepl("pvs", vf_ann3$gene_name)) |
                    (vf_ann3$gnm %in% vf_ann5$gnm[vf_ann5$pvuABCDE > 0] 
                     & grepl("pvu", vf_ann3$gene_name)), ]
vf_ann6$gene_name = ifelse(grepl("pvs", vf_ann6$gene_name), "pvsABCDE",
                           ifelse(grepl("pvu", vf_ann6$gene_name), "pvuABCDE", NA))

vf_ann = vf_ann6
ann_master_TAB$Vibrioferrin = vf_ann$gene_name[match(ann_master_TAB$ID, vf_ann$gene_ID)]
remove(list = ls(pattern='^vf_ann[0-9]+'))





### DMSP

#! load packages
library(stringr)
#! set variables
DMSP_ann1 = read.delim("3_Ann/DMSP_ddd/DMSP_blastp.tsv", h=F)

DMSP_ann1$gene_ID = gsub(".*\\.","",DMSP_ann1$V1)
DMSP_ann1$filename = gsub("_.*", "", DMSP_ann1$gene_ID)
DMSP_ann2 = DMSP_ann1
DMSP_ann2 = DMSP_ann2[DMSP_ann2$V11 < 1e-70, ] # filter based on min e-value - Gärdes et al., 2013
# DMSP_ann2 = DMSP_ann2[DMSP_ann2$V12 > 60, ] # filter based on min bit score
  
#! get coherent gene names
DMSP_ann2$gene_name = sapply(DMSP_ann2$V13, function(i) {
  ifelse(grepl("ddd", i, ignore.case = T), paste0("ddd", toupper(substr(gsub(".* ddd|.*=ddd", "", i, ignore.case = T), 1, 1))),
         ifelse(grepl("dmd", i, ignore.case = T), paste0("dmd", toupper(substr(gsub(".* dmd|.*=dmd", "", i, ignore.case = T), 1, 1))),
                ifelse(grepl("acu", i, ignore.case = T), paste0("acu", toupper(substr(gsub(".* acu|.*=acu", "", i, ignore.case = T), 1, 1))), "NA"
                )))
})

#! parse blast output
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

#! check for completeness
DMSP_ann4 = table(DMSP_ann3$filename, DMSP_ann3$gene_name)

DMSP_ann4a = DMSP_ann4[, grepl("ddd|ALMA", colnames(DMSP_ann4))]
DMSP_ann4a = apply(DMSP_ann4a, 1, function(j) sum(j>0))
DMSP_ann4a = names(DMSP_ann4a)[DMSP_ann4a > 0] #! completeness: require at least 1 gene (are all redundant)

DMSP_ann4b = DMSP_ann4[, grepl("dmd", colnames(DMSP_ann4))]
DMSP_ann4b = apply(DMSP_ann4b, 1, function(j) sum(j>0))
DMSP_ann4b = names(DMSP_ann4b)[DMSP_ann4b >= 3] #! completeness: require at least 3 of 4 genes

ann_master_TAB$DMSP = NA
ann_master_TAB$DMSP[match(DMSP_ann3$gene_ID[DMSP_ann3$filename %in% DMSP_ann4a], ann_master_TAB$ID)] = "DMSP_cleavage"
ann_master_TAB$DMSP[match(DMSP_ann3$gene_ID[DMSP_ann3$filename %in% DMSP_ann4b], ann_master_TAB$ID)] = "DMSP_demethylation"





### DHPS 
#!!! based on the KEGG annotation 

#! seek for the key KO:
#! "K15509" (sulfopropanediol 3-dehydrogenase, gene hpsN)
dhps_ann1 = ann_master_TAB[ann_master_TAB$KO %in% "K15509", ]

ann_master_TAB$DHPS = NA
invisible(apply(dhps_ann1, 1, function(j) {
  ann_master_TAB$DHPS[ann_master_TAB$filename %in% j[1] &
                        ann_master_TAB$KO %in% j[10]
  ] <<- "DHPS_catabolism"
}))





### Taurine 
#!!! based on the KEGG annotation 

#! seek for the key KOs:
#! 1) "K07255" + "K07256" (taurine dehydrogenase, tauXY genes) & K03852 (sulfoacetaldehyde acetyltransferase, xsc gene)
tau_ann1a = with(ann_master_TAB[ann_master_TAB$KO %in% c("K03851","K03852"), ], table(filename, KO))
tau_ann2a = apply(tau_ann1a, 1, function(j) sum(j>0))
tau_ann3a = names(tau_ann2a)[tau_ann2a > 1] #! completeness: require all genes (2)

ann_master_TAB$Taurine[ann_master_TAB$filename %in% tau_ann3a &
                         ann_master_TAB$KO %in% c("K03851","K03852")] = "Taurine_in_TCA"


#! 2) "K03851" (taurine-pyruvate aminotransferase, tpa gene) & K03852 (sulfoacetaldehyde acetyltransferase, xsc gene)
tau_ann1b = with(ann_master_TAB[ann_master_TAB$KO %in% c("K07255","K07256","K03852"), ], table(filename, KO))
tau_ann2b = apply(tau_ann1b, 1, function(j) sum(j>0))
tau_ann3b = names(tau_ann2b)[tau_ann2b > 2] #! completeness: require all genes (3)

ann_master_TAB$Taurine[ann_master_TAB$filename %in% tau_ann3b &
                         ann_master_TAB$KO %in% c("K07255","K07256","K03852")] = "Taurine_in_TCA"







### pre-III benchmark KMreco (optional) --------------------------------------------------

#!!! perform benchmarking of KM reconstruction 
#!!! investigate the best setting for KM completness



### run KMreco with different parameters for gap filling
library(reshape2)
source(KM_code)

# ## download new KM structures
# KM_str = KMdiagram_fetcher(ncore = 7, create_RData = T, path = "~/Downloads") # will save the fetched file in the current wd
# load KM structures used in the paper
KM_str = readRDS(KMds)

ann_master_TAB$cnt = 1
KO_mat = dcast(ann_master_TAB, filename~KO, value.var = "cnt")
row.names(KO_mat) = KO_mat$filename
KO_mat = KO_mat[, -1]
KO_mat[KO_mat > 1] = 1

# # !!! next command takes time to run, better to save the output
# KM_bm = list(KMreco(indata = KO_mat, KM_str = KM_str, len_breaks = c(2), allowed_gaps = c(0,1))$KM_pa,
#              KMreco(indata = KO_mat, KM_str = KM_str, len_breaks = c(3), allowed_gaps = c(0,1))$KM_pa,
#              KMreco(indata = KO_mat, KM_str = KM_str, len_breaks = c(4), allowed_gaps = c(0,1))$KM_pa,
#              KMreco(indata = KO_mat, KM_str = KM_str, len_breaks = c(5), allowed_gaps = c(0,1))$KM_pa,
#              KMreco(indata = KO_mat, KM_str = KM_str, len_breaks = c(6), allowed_gaps = c(0,1))$KM_pa
# )
# saveRDS(KM_bm, "KMreco_benchmark.RDS")
KM_bm = readRDS("KMreco_benchmark.RDS")


### investigate KMs that change in completeness frequency
KM_bm2 = lapply(KM_bm, function(i) apply(i, 2, function(j) sum(j)/nrow(i)))
KM_bm2 = lapply(seq_along(KM_bm2), function(i) data.frame(type=i,
                                                          KM=names(KM_bm2[[i]]), 
                                                          freq=as.numeric(KM_bm2[[i]])))
KM_bm3 = do.call(rbind, KM_bm2)
KM_bm3 = dcast(KM_bm3, KM~type, fill = 0)
names(KM_bm3)[-1] = c("1gap2rxn","1gap3rxn","1gap4rxn","1gap5rxn","1gap6rxn")
KM_bm4 = KM_bm3[apply(KM_bm3[, -1], 1, function(z) length(unique(z)) > 1), ] #remove unchanged
KM_bm4 = KM_bm4[apply(KM_bm4[, -1], 1, function(z) any(z > 0.03)), ] #remove KM below min frequency threshold (3% of genomes)

## make table
# ??? missing text names for KM of length 2
library(readxl)
KM_meta = as.data.frame(read_xlsx(suppl_tabs, col_names = T, sheet = 5, skip = 3))
KM_meta = KM_meta[!is.na(KM_meta$Description), ]
# KM_meta = read.delim("/media/lucaz/DATA/DBs_repository/KEGG/KEGG_M_hier_9jan20plus.txt")
# KM_meta$ID = gsub(" .*", "", KM_meta[, 4]) 
# names(KM_meta)[4]="Description"
# KM_meta=KM_meta[, c(5,5,5,4,2,3)]
tKM_bm = data.frame(KM=KM_meta$Description[match(KM_bm4[,1], KM_meta$ID)], 
                    KM_len=sapply(KM_bm4[,1], function(i) KM_str[[i]]$length),
                    KM_lenmacro=sapply(KM_bm4[,1], function(i) KM_str[[i]]$length_macro),
                    KM_bm4[, -1], 
                    CV= apply(KM_bm4[, -1], 1, function(i) (sd(i)/mean(i))*100),
                    KM_meta[match(KM_bm4[, 1], KM_meta$ID), c(5,6)])
tKM_bm$KM_lenmacro[tKM_bm$KM_lenmacro == ""] = "-"
tKM_bm$KM = gsub(" \\[.*", "", tKM_bm$KM)
write.table(tKM_bm, "KMreco_benchmark.res", col.names = T, row.names = F, quote = F, sep = "\t", na = "")

## plot
library(ggplot2)
library(ggthemes)
pKM_bm = melt(KM_bm4)
pKM_bm$KM_len = sapply(pKM_bm$KM, function(i) 
  ifelse(KM_str[[i]]$just_macro, KM_str[[i]]$length_macro, KM_str[[i]]$length))
ggplot(pKM_bm, aes(x = value, y = after_stat(count), colour = variable)) +
  geom_freqpoly(binwidth = 0.05) + labs(x="Frequency", y="#KMs") +
  xlim(c(0,1)) + theme_minimal()







### III Format datasets & Interaction traits --------------------------------------------------

##### reconstruct KM
library(reshape2)
source(KM_code)

# ## download new KM structures
# KM_str = KMdiagram_fetcher(ncore = 7, create_RData = T, path = "~/Downloads") # will save the fetched file in the current wd

## load KM structures used in the paper
KM_str = readRDS(KMds)

ann_master_TAB$cnt = 1
KO_mat = dcast(ann_master_TAB, filename~KO, value.var = "cnt")
row.names(KO_mat) = KO_mat$filename
KO_mat = KO_mat[, -1]
KO_mat[KO_mat > 1] = 1
KM_ann = KMreco(indata = KO_mat, KM_str = KM_str, 
                len_breaks = c(3), allowed_gaps = c(0,1))



### fix B1 completeness - M00127 is coded in a weird way in KEGG
b1_KOs = data.frame(KO=c("K03147","K00878","K03149","K00788","K00941","K14153","K14154","K21219","K21220","K00946"),
                     Gene=c("thiC","thiM","thiG","thiE","thiD","thiDE","THI6","thiDN","thiN","thiL"),
                     Step=c("Step_1","Step_1","Step_1","Step_2","Step_2","Step_2","Step_2","Step_2","Step_2","Step_3"))
b1_KOs = b1_KOs[b1_KOs$KO %in% colnames(KO_mat), ]
b1_gnm = KO_mat[, as.character(b1_KOs$KO)]

new_b1_compl = apply(b1_gnm, 1, function(j) {
  j2 = split(j, f = b1_KOs$Step) # split KOs list into unique levels
  j3 = sum(unlist(lapply(j2, function(k) any(k > 0)))) # levels with at least 1 ann KO
  ifelse(length(unique(b1_KOs$Step)) - j3 <= 0, 1, 0) #! completeness: require at least 2 of 3
})
new_b1_compl = new_b1_compl[match(row.names(KM_ann$KM_pa), names(new_b1_compl))]
KM_ann$KM_pa[, "M00127"] = new_b1_compl



### fix B12 completeness - M00122 is only last half
b12_KOs = data.frame(KO=c("K13786","K02303","K00798","K19221","K02232","K02225","K02227","K02231","K19712","K02233"),
                     Gene=c("cobR","cobA","MMAB,pdu0","cobO","cobQ,CbiP ","cobC","cobD,cbiB","cobP,cobU","cobY","cobS,cobV"),
                     Step=c("Step_1","Step_2","Step_2","Step_2","Step_3","Step_4","Step_4","Step_5","Step_5","Step_6"))
b12_KOs = b12_KOs[b12_KOs$KO %in% colnames(KO_mat), ]
b12_gnm = KO_mat[, as.character(b12_KOs$KO)]

new_b12_compl = apply(b12_gnm, 1, function(j) {
  j2 = split(j, f = b12_KOs$Step) # split KOs list into unique levels
  j3 = sum(unlist(lapply(j2, function(k) any(k > 0)))) # levels with at least 1 ann KO
  ifelse(length(unique(b12_KOs$Step)) - j3 <= 1, 1, 0) #! completeness: require at least 5 of 6
})
new_b12_compl = new_b12_compl[match(row.names(KM_ann$KM_pa), names(new_b12_compl))]
KM_ann$KM_pa[, "M00122"] = new_b12_compl





##### create trait-table
library(tidyr)
library(dplyr)
fun_profi = cbind(KM_ann$KM_pa
                  , ann_master_TAB %>% 
                    pivot_wider(id_cols = filename, names_from = TCdb, values_from = cnt, 
                                values_fn = sum, values_fill = 0) %>% 
                    select(!c("filename", "NA")) %>% 
                    replace(., . > 1, 1) 
                  , ann_master_TAB %>% 
                    pivot_wider(id_cols = filename, names_from = SecMetab, values_from = cnt, 
                                values_fn = sum, values_fill = 0) %>% 
                    select(!c("filename", "NA")) %>% 
                    replace(., . > 1, 1) 
                  , ann_master_TAB %>% 
                    pivot_wider(id_cols = filename, names_from = Vibrioferrin, values_from = cnt, 
                                values_fn = sum, values_fill = 0) %>% 
                    select(!c("filename", "NA")) %>% 
                    replace(., . > 1, 1) 
                  , ann_master_TAB %>% 
                    pivot_wider(id_cols = filename, names_from = Phytohormones, values_from = cnt, 
                                values_fn = sum, values_fill = 0) %>% 
                    select(!c("filename", "NA") & !contains(";")) %>% 
                    replace(., . > 1, 1) 
                  , ann_master_TAB %>% 
                    pivot_wider(id_cols = filename, names_from = DMSP, values_from = cnt, 
                                values_fn = sum, values_fill = 0) %>% 
                    select(!c("filename", "NA")) %>% 
                    replace(., . > 1, 1) 
                  , ann_master_TAB %>% 
                    pivot_wider(id_cols = filename, names_from = DHPS, values_from = cnt, 
                                values_fn = sum, values_fill = 0) %>% 
                    select(!c("filename", "NA")) %>% 
                    replace(., . > 1, 1) 
                  , ann_master_TAB %>% 
                    pivot_wider(id_cols = filename, names_from = Taurine, values_from = cnt, 
                                values_fn = sum, values_fill = 0) %>% 
                    select(!c("filename", "NA")) %>% 
                    replace(., . > 1, 1) 
)


### filter out noise
mtf = 0.03 # Minimum trait frequency
good_traits = colnames(fun_profi)[colSums(fun_profi)/nrow(fun_profi) > mtf]
# good_traits = good_traits[!good_traits %in% colnames(trait_cor)[colSums(is.na(trait_cor)) > 0]]
fun_profi_filt = fun_profi[, good_traits] 







### IV add trait metadata --------------------------------------------------

### create metadata table for traits
trait_meta = data.frame(ID = colnames(fun_profi),
                        RA = colSums(fun_profi)/nrow(fun_profi),
                        "Annotation" = c(rep("KEGG modules", ncol(KM_ann$KM_pa))
                                         , rep("TCdb", length(na.omit(unique(ann_master_TAB$TCdb))))
                                         , rep("AntiSMASH", length(na.omit(unique(ann_master_TAB$SecMetab))))
                                         , rep("Manual", length(na.omit(unique(ann_master_TAB$Phytohormones)[!grepl(";", unique(ann_master_TAB$Phytohormones))])))
                                         , rep("Manual", length(na.omit(unique(ann_master_TAB$Vibrioferrin))))
                                         , rep("Manual", length(na.omit(unique(ann_master_TAB$DMSP))))
                                         , rep("Manual", length(na.omit(unique(ann_master_TAB$DHPS))))
                                         , rep("Manual", length(na.omit(unique(ann_master_TAB$Taurine))))
                        ))
library(readxl)
suppl_info = lapply(4:7, function(i) {
  meta = as.data.frame(read_xlsx(suppl_tabs, col_names = T, sheet = i, skip = 3)) 
  meta = meta[!is.na(meta$Description), ] # remove extra lines in Excel table
})
suppl_info = Reduce(function(x, y) 
  merge(x, y, by=c("ID","Description","Frequency","Interaction traits","Notes","Reference"), 
        all=T), suppl_info)
trait_meta = merge(trait_meta, suppl_info, by="ID", all=T)
trait_meta = trait_meta[, c(1,2,3,4,6,10,11,7)]

## add details on transporter
trait_meta$Description[trait_meta$Annotation == "TCdb" & is.na(trait_meta$Description)] = 
  paste(transp_ann$Predicted_Substrate[match(trait_meta$ID[trait_meta$Annotation == "TCdb" & is.na(trait_meta$Description)], transp_ann$Hit_tcid)], "transporter")
# write.table(trait_meta, "trait_meta.tsv", 
#             col.names = T, row.names = F, quote = F, sep = "\t", na = "")

trait_meta_filt = trait_meta[match(colnames(fun_profi_filt), trait_meta$ID), ]







### save ----

save(ann_master_TAB,
     KO_ann, transp_ann, sm_ann, ph_ann, vf_ann,
     KM_ann, fun_profi, fun_profi_filt,
     trait_meta, trait_meta_filt,
     file = "Atlas_ann_v1.RData")
