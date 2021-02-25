## set working directory
setwd("/media/lucaz/DATA/HB_complete_genomes_v4/")


## global setting
options(stringsAsFactors = F)
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8") # to get month names in english


## libraries
# install.packages("BiocManager")
# BiocManager::install(c("apcluster"))
library(apcluster)
library(reshape2)
library(vegan)
library(qualpalr); library(DescTools) # creates colored and gray palettes


##!!! required files
# Genome metadata (use Supplementary table 2 to replicate the manuscript anaysis or any metadata file of the genomes you are investigating; it needs a column containing the filename of the genome fasta for the mapping)
gnm_meta = "/media/lucaz/DATA/Google_drive/Laboratorio/AAA_Progetti_NUOVI/2016_HFSP/4_paper_genome_comparison/6_REV_Comm_biology/Suppl_tables.xlsx" #! adjust the path
# ann_master_TAB, gnm_profi, gnm_profi_filt and trait_meta_filt objects from the code ATLAS_annotation.R
load("Atlas_ann_v1.RData")


#



##### Create Genome Functional Clusters (GFCs) ------------------------------------------------------

### calculate genome pairwise correlation (i.e. Pearson - r)
gnm_cor = cor(t(gnm_profi))
gnm_cor[gnm_cor < 0] = 0 #! there are no neg correlation (as expected)



### run affinity propagation

# #!!! code to check fitting of parameter "q"
# x = seq(from=0, to=0.95, length=100)
# y = sapply(x, function(i) {
#   length(apcluster(s = gnm_cor, details=F, q=i,
#                    lam=0.5, seed=1234, maxits=1000, convits=200)@clusters)
# })
# plot(x,y, xlab="Q-vaules", ylab="# GFCs")
# #---
apcl_gnm = apcluster(s = gnm_cor, details=T, q=0.5, lam=0.5, seed=1234, maxits=1000, convits=500)

# heatmap(apcl_gnm, gnm_cor); # plot(apcl_gnm, gnm_profi); 
gnm_GFC_tmp = do.call(rbind, lapply(1:length(apcl_gnm@clusters), function(i) data.frame(i, apcl_gnm@clusters[[i]])))
gnm_GFC = gnm_GFC_tmp$i[order(gnm_GFC_tmp$apcl_gnm.clusters..i.., decreasing = F)]
table(gnm_GFC)


## extract hierarchical cluster from apcluster results
gnm_dist = as.matrix(cophenetic(as.dendrogram(aggExCluster(s = gnm_cor, x = apcl_gnm))))
gnm_dist = gnm_dist[match(row.names(gnm_cor), row.names(gnm_dist)), 
                    match(colnames(gnm_cor), colnames(gnm_dist))]
gnm_hc = hclust(as.dist(gnm_dist), "complete")
gnm_hc = reorder(gnm_hc, wts = colSums(gnm_cor), agglo.FUN = "mean") # improve dendro sorting
gnm_hc$order = as.integer(gnm_hc$order) # otherwise it rises an issue when plotting with iheatmapr



### parse results
GFC_table = data.frame(gnm=row.names(gnm_profi), GFC=gnm_GFC)

## add metadata
gnm_meta2 = as.data.frame(readxl::read_xlsx(gnm_meta, col_names = T, sheet = "S_tab_2", skip = 3)) 
gnm_meta2 = gnm_meta2[!is.na(gnm_meta2$Filename), ] # remove extra lines in Excel table
#! the GFC info included in Supplementary table 2 will be overwritten by this code
GFC_table = cbind(GFC_table[, -c(1), drop=F], 
                  gnm_meta2[match(GFC_table$gnm, gnm_meta2$Filename), -1])

## add annotation info
GFC_table$`#genes` = sapply(GFC_table$Filename, function(i) sum(ann_master_TAB$filename == i))
GFC_table$`Gene annotated` = sapply(GFC_table$Filename, function(i) { 
  sum(apply(ann_master_TAB[ann_master_TAB$filename == i, c(6,8:17)], #! adjust the column range to include all you annotations in the master file
            1, function(j) any(!is.na(j))))
})
GFC_table$`Gene annotated (%)` = round(GFC_table$`Gene annotated` / GFC_table$`#genes`, 2)







##### Create Linked Trait CLusters (LTCs) ------------------------------------------------------

### calculate trait pairwise correlation (i.e. Pearson - r)

f_r = function(x, y) {
  cont_tab = table(factor(x, levels = c(0, 1)), factor(y, levels = c(0, 1)))
  
  ## add 1 fake absence to traits without O values (M00005, M00052)
  if (sum(cont_tab[2, ]) == sum(cont_tab) | (sum(cont_tab[, 2]) == sum(cont_tab))) { # no zeros for x
    cont_tab[1, 1] = 1
  }
  ##---
  
  joint_p = cont_tab/sum(cont_tab)
  Pab = joint_p[4]
  Pa = sum(joint_p[2, ])
  Pb = sum(joint_p[, 2])
  
  D = Pab - (Pa * Pb)
  r = D / sqrt(Pa * (1 - Pa) * Pb * (1 - Pb))
  return(r)
}

pairFUN = function(mydata, fun.xy, ncore) {
  # if function ('fun.xy') generates a vector with multiple output values, specify the index ('val_index') of the desired one 
  
  smpl_comb = as.matrix(combn(x = 1:ncol(mydata), m = 2))
  
  library(future.apply)
  plan(multiprocess, workers = ncore)
  comb_out = future_apply(X = smpl_comb, MARGIN = 2, FUN = function(z) 
    fun.xy(mydata[, z[1]], mydata[, z[2]]))
  
  comb_mat = matrix(NA, nrow = ncol(mydata), ncol = ncol(mydata))
  comb_mat[lower.tri(comb_mat, diag = F)] = comb_out
  comb_mat = as.matrix(as.dist(comb_mat))
  diag(comb_mat) = 1
  dimnames(comb_mat) = list(colnames(mydata), colnames(mydata))
  
  return(comb_mat)
}

trait_cor = pairFUN(gnm_profi_filt, fun.xy = f_r, ncore = 7)
trait_cor[trait_cor < 0] = 0


## test for significant correlation (r)
pair_r = trait_cor
pair_r[upper.tri(trait_cor, diag = T)] = NA
pair_r = melt(pair_r, na.rm = T)
pair_r$x2 = sapply(pair_r$value, function(i) i^2 * nrow(gnm_profi_filt)) # x2 test
pair_r$p.val = sapply(pair_r$x2, function(i) pchisq(i, df=2, lower.tail=F))
pair_r$p.val.adj = p.adjust(pair_r$p.val, method = "fdr")
min_signif_r = min(pair_r$value[pair_r$p.val.adj <= 0.05])
trait_cor[trait_cor < min_signif_r] = 0



### run affinity propagation

# #!!! code to check fitting of parameter "q"
# x = seq(from=0, to=0.95, length=100)
# y = sapply(x, function(i) {
#   length(apcluster(s = trait_cor, details=F, q=i,
#                    lam=0.5, seed=1234, maxits=1000, convits=200)@clusters)
# })
# plot(x,y, xlab="Q-vaules", ylab="# LTCs")
# #---
apcl_trait = apcluster(s = trait_cor, details=T, q=0.5, lam=0.5, seed=1234, maxits=1000, convits=500)

# heatmap(apcl_trait, trait_cor); # plot(apcl_trait, gnm_profi); 
trait_LTC_tmp = do.call(rbind, lapply(1:length(apcl_trait@clusters), function(i) data.frame(i, apcl_trait@clusters[[i]])))
trait_LTC = trait_LTC_tmp$i[order(trait_LTC_tmp$apcl_trait.clusters..i.., decreasing = F)]
table(trait_LTC)


## extract hierarchical cluster from apcluster results
trait_dist = as.matrix(cophenetic(as.dendrogram(aggExCluster(s = trait_cor, x = apcl_trait))))
trait_dist = trait_dist[match(row.names(trait_cor), row.names(trait_dist)), 
                        match(colnames(trait_cor), colnames(trait_dist))]
trait_hc = hclust(as.dist(trait_dist), "complete")
trait_hc = reorder(trait_hc, wts = colSums(trait_cor), agglo.FUN = "mean") # improve dendro sorting
trait_hc$order = as.integer(trait_hc$order) # otherwise it rises an issue when plotting with iheatmapr



### parse results

#! mark LTCs with ITs
LTC_with_IT = as.matrix(table(trait_LTC, !is.na(trait_meta_filt$`Interaction traits`)))
LTC_with_IT = LTC_with_IT[!grepl("uncl.", row.names(LTC_with_IT)), ] # avoid marking 'uncl.'
trait_LTC[trait_LTC %in% row.names(LTC_with_IT)[LTC_with_IT[, "TRUE"] > 0]] = paste0(trait_LTC[trait_LTC %in% row.names(LTC_with_IT)[LTC_with_IT[, "TRUE"] > 0]], "*")

LTC_table = data.frame(LTC=trait_LTC, 
                       LTC_size=as.integer(table(trait_LTC))[match(trait_LTC, names(table(trait_LTC)))])


## find bearing GFCs
#! LTC completeness in gnm (> 60% genetic traits are complete in a genome)
LTC_in_gnm = lapply(split(as.data.frame(t(gnm_profi_filt)), f=LTC_table$LTC), 
                    function(i) ifelse(colMeans(i) > 0.6, 1, 0))
LTC_in_gnm = do.call(rbind, LTC_in_gnm)
#! LTC completeness in GFC (> 60% genomes have that LTC complete)
LTC_in_GFC = lapply(split(as.data.frame(t(LTC_in_gnm)), f=GFC_table$GFC), 
                    function(i) ifelse(colMeans(i) > 0.6, 1, 0))
LTC_in_GFC = do.call(rbind, LTC_in_GFC)
LTC_in_GFC = melt(LTC_in_GFC, varnames = c("GFC", "LTC"))
LTC_in_GFC = LTC_in_GFC[LTC_in_GFC$value > 0, ]
LTC_in_GFC = LTC_in_GFC[order(LTC_in_GFC$GFC, decreasing = F), ]
LTC_in_GFC = aggregate(GFC~LTC, LTC_in_GFC, FUN = function(i) paste(i, collapse = ","))

LTC_table$bearing_GFC = LTC_in_GFC$GFC[match(LTC_table$LTC, LTC_in_GFC$LTC)]
LTC_table$LTC[is.na(LTC_table$bearing_GFC)] = "uncl."
LTC_table$LTC_size[is.na(LTC_table$bearing_GFC)] = sum(is.na(LTC_table$bearing_GFC))


## mean correlation (r) in each LTC
LTC_mean_r = as.data.frame(t(sapply(unique(LTC_table$LTC), function(i) {
  i2 = trait_cor[LTC_table$LTC == i, LTC_table$LTC == i]
  c(LTC=i, mean_r=mean(i2[upper.tri(i2, diag = F)]))
})))
LTC_table$mean_r = round(as.numeric(LTC_mean_r$mean_r[match(LTC_table$LTC, LTC_mean_r$LTC)]), 2)


## add metadata
LTC_table = cbind(LTC_table, trait_meta_filt)







### save ----

save(ann_master_TAB,
     KO_ann, transp_ann, sm_ann, ph_ann, vf_ann,
     KM_ann, gnm_profi, gnm_profi_filt,
     trait_meta, trait_meta_filt,
     gnm_cor, gnm_hc, GFC_table, trait_cor, trait_hc, LTC_table,
     file = "Atlas_ann_v1_stat.RData")
