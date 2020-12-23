KMdiagram_fetcher <- function(ncore, create_RData=T, path=getwd(), new_date=Sys.Date()) {
  #    "ncore": degree of parallelization, can be > of real number of threads but don't exagerate otherwise the fetching script will rise an error
  #    "create_RData": whether to save or not the fetched list of KM diagrams as RData object.
  #                    It saves also a list of KOs which are marked as optional in KM diagrams 
  #    "path": to the output folder for the RData object
  #    "new_date": present date used in the output filname e.g. KM_str_'2019-04-12'.RData
  
  
  ## (get &)load packages
  if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")
  pkg_list <- c("stringr", "xml2", "KEGGREST", "future.apply")
  new_pkg <- pkg_list[!(pkg_list %in% installed.packages()[,"Package"])]
  if (length(new_pkg)) BiocManager::install(new_pkg)
  sapply(pkg_list, library, character.only = T)
  
  
  ## get list of KM names
  x <- readLines("https://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir=")
  x <- x[grep("^D", x)]
  KM_list <- str_extract(x, "M[0-9]{5}")
  
  
  ## fetch diagrams
  plan(multiprocess, workers = ncore)
  fetched_KM <- future_lapply(KM_list, function(imod) { 
    txt <- read_html(paste0("http://www.genome.jp/kegg-bin/show_module?", imod))
    txt2 <- gsub(".*<map id=\"module\" name=\"module\">|</map>.*", "", txt)
    txt3 <- strsplit(txt2, "<area shape=\"rect\" ")[[1]][-1]
    txt4 <- t(as.data.frame(strsplit(txt3, '\\" ')))[, c(1,3,5), drop=F]
    txt5 <- data.frame(row.names = NULL, KM_lvl=gsub("_.*", "", gsub('id=\"', "", txt4[,1])), KO_id=gsub(".*_", "", gsub('id=\"', "", txt4[,1])),
                       KO=gsub('title="| .*', "", txt4[,2]), t(as.data.frame(strsplit(gsub('coords="|">.*', "", txt4[,3]), ","))))
    names(txt5)[4:7] <- c("X1","Y1","X2","Y2")
    ## remove optional terms '-K00000'
    txt_str <- keggGet(paste0("md:", imod))[[1]]$DEFINITION
    optKOs <- unlist(c(gsub("\\-", "", str_extract_all(txt_str, "\\-K[0-9]{5}")[[1]]),
                       str_extract_all(str_extract_all(txt_str, "\\-\\([K0-9,]*[^)]\\)")[[1]], "K[0-9]{5}")))
    txt_str_clean <- gsub(paste0("-", optKOs, collapse = "|"), "", txt_str)
    txt6 <- txt5[!txt5$KO %in% optKOs, ]
    #---
    
    KM_len <- length(unique(txt6$KM_lvl))
    
    ## if single-level KM check whether it just contain a single macromolecular complex (simplify completness/specificity analysis)
    if (KM_len == 1 & str_count(paste(txt_str, collapse = " "), "\\(|\\)|,| ") == 0) { just_macro = TRUE
    } else { just_macro = FALSE }
    #---
    
    KM_len_macro <- ""; if (just_macro == TRUE) { KM_len_macro <- str_count(txt_str_clean, "K[0-9]{5}") }
    iKM <- list(just_macro=just_macro, length=KM_len, length_macro=KM_len_macro, str=txt6)
    
    return(list(optKOs, iKM))
    } )
  
  KM_str <- lapply(fetched_KM, function(i) i[[2]])
  names(KM_str) <- KM_list
  KM_str <- KM_str[order(as.integer(gsub("M", "", KM_list)))]
  opt_KOs <- lapply(fetched_KM, function(i) i[[1]]) # list of otional KOs in modules
  names(opt_KOs) <- KM_list
  opt_KOs <- opt_KOs[order(as.integer(gsub("M", "", KM_list)))]
  
  if (create_RData) {
    dir.create(path, recursive = T)
    save(KM_str, opt_KOs, file = file.path(path, paste0("KM_str_", new_date,".RData")))
  }
  return(KM_str)
}

KMreco <- function(indata, KM_str, len_breaks=NULL, allowed_gaps=c(0)) {
  #!!!    "indata": file must be a presence/absence (i.e. 1/0) matrix/data.frame 
  #       with genomes as rows and annotated KOs as columns
  #!!!    "KM_str": contains KM diagrams, obtained with 'KMdiagram_fetcher' function
  #!!!    "len_breaks": vector of integer of open-right limits for binning of KM length; e.g. c(3, 10)
  #!!!    "allowed_gaps": vector of integer specfing nomber of gaps allowed in each bin; it has length "breaks"+1; e.g. c(0,1,2) 
  indata <- as.data.frame(indata)
  options(stringsAsFactors = F)
  
  
  
  ##### Parse genomes to get KM completness
  KMgnm_compl <- data.frame(KM=numeric(0))
  for (gnm in row.names(indata)) { # loop thru genomes
    i_KM_compl <- c()
    
    # get KO list
    i_gnm <- indata[gnm, , drop=F]
    i_gnmKO <- colnames(i_gnm)[i_gnm > 0]
    
    # loop through all KM
    for (KM in names(KM_str)) { 
      i_KM_str <- KM_str[[KM]]
      
      #1 look which KOs are present in the genome
      i_KM_str$str$exist <- i_KM_str$str$KO %in% i_gnmKO
      if (sum(i_KM_str$str$exist) == 0) { #!!! when no KO of the KM is present in the genome
        i_KM_compl <- rbind.data.frame(i_KM_compl, c(KM, 0))
      } else { #!!! when at least some KOs of the KM are present in the genome
        
        #2 per each KM level get completness
        if (i_KM_str$just_macro == T) { # when single-level macromolecular complex
          i_KM_existKO <- i_KM_str$str$exist

        } else { # when is not a single-level macromolecular complex
          i_KM_existKO <- c()
          for (i_lvl in unique(i_KM_str$str$KM_lvl)) { # loop through levels of that KM
            i_lvl_KM_str <- i_KM_str$str[i_KM_str$str$KM_lvl %in% i_lvl, , drop=F]
            i_lvl_KM_str[, c(4:7)] <- apply(i_lvl_KM_str[, c(4:7)], 2, function(x) as.numeric(x)) # fix
            
            if (nrow(i_lvl_KM_str) == 1) { #!!! if level has single-KO
              i_KM_existKO <- c(i_KM_existKO, i_lvl_KM_str$exist)
              
            } else { #!!! if level has multiple-KOs
              #2a     find first row of KOs
              sol_starts <- which(i_lvl_KM_str$Y1 == min(i_lvl_KM_str$Y1))
              #2b    looks for presence of gaps between gene1 X2 - gene2 X1
              sol_x2_lim <- c()
              if (length(sol_starts) > 1) { # multiple hypotetic solution
                for (l in 1:(length(sol_starts)-1)) {
                  if (i_lvl_KM_str[sol_starts[l], "X2"] != i_lvl_KM_str[sol_starts[l+1], "X1"] & # when there is a gap
                      i_lvl_KM_str[sol_starts[l+1], "X1"] - i_lvl_KM_str[sol_starts[l], "X2"] == 4) { #!! fix fake gap that can be with nested redundant level of different length (e.g. M00003)
                    sol_x2_lim <- c(sol_x2_lim, i_lvl_KM_str[sol_starts[l], "X2"])
                    if (l == length(sol_starts)-1) { sol_x2_lim <- c(sol_x2_lim, i_lvl_KM_str[sol_starts[l+1], "X2"]) } # add last level solution
                  } 
                }
              }
              if (length(sol_x2_lim) == 0) { sol_x2_lim <- c(sol_x2_lim, i_lvl_KM_str[sol_starts[length(sol_starts)], "X2"]) } # when there is only one solution or none of hypotetic split was good 
              #2c    split level in all possible solutions - X-wise
              i_lvl_KM_str2 <- list()
              for (l in 1:length(sol_x2_lim)) {
                i_lvl_KM_str2[[l]] <- i_lvl_KM_str[i_lvl_KM_str$X2 <= sol_x2_lim[l], ]
                i_lvl_KM_str <- i_lvl_KM_str[!i_lvl_KM_str$X2 <= sol_x2_lim[l], ] # remove parsed element from level
              }
              #2d    split each solution Y-wise
              for (l in 1:length(i_lvl_KM_str2)) {
                i_lvl_KM_str2[[l]] <- split(i_lvl_KM_str2[[l]], f = as.factor(i_lvl_KM_str2[[l]]$Y1))
              }
              #2e    pull out all possible KO-combinations from each solution
              sol_paths_exist <- list()
              for (l in 1:length(i_lvl_KM_str2)) { # loop thru solutions - for completness
                sol_paths_exist[[l]] <- expand.grid(lapply(i_lvl_KM_str2[[l]], function(z) z$exist)) # use z$KO for debug
              }
              #2f    discard each solution with missing KOs
              sol_paths_exist2 <- list()
              for (l in 1:length(sol_paths_exist)) {
                tmp_exist <- sol_paths_exist[[l]][apply(sol_paths_exist[[l]], 1, function(y) sum(y) == length(y)), ]
                if (length(tmp_exist) > 0) { # keep only sol without 0
                  sol_paths_exist2 <- c(sol_paths_exist2, list(tmp_exist))
                }
              }
              #2g    look whether there is a solution for this KM-level
              if (length(unlist(sol_paths_exist2)) > 0) {
                sol_exist <- TRUE
              } else { sol_exist <- FALSE}
              
              i_KM_existKO <- c(i_KM_existKO, sol_exist)
            }
          }
        }
        
        #3 estimate KM completness
        i_KM_compl <- rbind.data.frame(i_KM_compl, c(KM, sum(i_KM_existKO)))
      }
    }
    
    names(i_KM_compl) <- c("KM", gnm)
    KMgnm_compl <- merge(KMgnm_compl, i_KM_compl, by.x="KM", by.y="KM", all = T)
  }
  
  KMgnm_compl2 <- t(apply(KMgnm_compl[, -1], 2, as.numeric))
  colnames(KMgnm_compl2) <- KMgnm_compl[, 1]

  
  
  ### apply completness rules
  
  #! get KM lengths
  KM_len <- unlist(lapply(KM_str, "[[", "length"))
  
  
  #! parse macromolecular complex KMs with same completness rules as other KMs
  macro_complx <- unlist(lapply(KM_str, "[[", "just_macro"))
  KM_len[macro_complx] <- as.integer(unlist(lapply(KM_str, "[[", "length_macro")))[macro_complx]
  
  
  #! estimate gaps allowed for each KM
  if (!is.null(len_breaks)) {
    KM_len_bin <- .bincode(KM_len, breaks = c(0, len_breaks, Inf), right = F, include.lowest = T)
  } else {
    KM_len_bin <- .bincode(KM_len, breaks = c(0, Inf), right = F, include.lowest = T)
  }
  KM_allowed_gaps <- allowed_gaps[KM_len_bin]

  
  #! sort vectors
  KM_allowed_gaps = KM_allowed_gaps[match(colnames(KMgnm_compl2), names(KM_len))]
  KM_len = KM_len[match(colnames(KMgnm_compl2), names(KM_len))]
  
  
  #! add 'allowed_gaps' as bonus
  KMgnm_compl2_filt = t(apply(KMgnm_compl2, 1, function(i) i + KM_allowed_gaps))
  
  
  #! check for completness
  KMgnm_compl2_filt = t(t(KMgnm_compl2_filt)/KM_len)
  KMgnm_compl2_filt[KMgnm_compl2_filt >=1] = 1
  KMgnm_compl2_filt[KMgnm_compl2_filt <1] = 0
  
  
  #! prepare outputs
  KMgnm_compl3 <- KMgnm_compl2_filt[, colSums(KMgnm_compl2_filt) > 0]
  
  KMgnm_compl2_RA <- t(t(KMgnm_compl2)/KM_len)
  KMgnm_compl2_RA <- KMgnm_compl2_RA[, colSums(KMgnm_compl2_RA) > 0]

  # object in list are data.frame with genomes as rows and KMs as columns
  tmp <- list(KMgnm_compl2,                   # number of KOs in each KM
              KM_len,                         # length of KMs (min number of reactions)
              KMgnm_compl2_RA,                # percentage of KM completness
              KMgnm_compl3)                   # presence-absence of complete KMs
  names(tmp) <- c("KM_#KOs", "KM_lenght", "KM_compl", "KM_pa")
  return(tmp) 
  ##---
}


# ### USAGE EXAMPLE ###
# KM_str <- KMdiagram_fetcher(ncore = 7, create_RData = T, path = "~")
# KMreco <- KMreco(indata = myannotation, KM_str = KM_str, len_breaks = c(3), allowed_gaps = c(0,1))
# ### ------------- ###


KMredund <- function(indata, KM_str, len_breaks, allowed_gaps) { #!!! STILL TO CHECK AND DEBUG!!!
  #!!!    "indata": output oject of KMreco, i.e. 'KM_pa'
  #!!!    "KM_str": contains KM diagrams, obtained with 'KMdiagram_fetcher' function
  #!!!    "len_breaks": vector of integer of open-right limits for binning of KM length; e.g. c(3, 10)
  #!!!    "allowed_gaps": vector of integer specfing nomber of gaps allowed in each bin; it has length "breaks"+1; e.g. c(0,1,2) 
  options(stringsAsFactors = F)
  
  
  ##### Parse genomes to get KM specificity & redundanc
  #NBB!!! this script parse only KM which have been declared complete -> any gap has to be filled
  library(IRanges)
  KMgnm_spec <- data.frame(KM=numeric(0))
  KMgnm_redund <- data.frame(KM=numeric(0))
  KMgnm_redund_bis <- data.frame(KM=numeric(0))
  KMgnm_gapfillKOs <- list()
  
  for (gnm in 1:nrow(KMgnm_compl3)) { # loop thru genomes
    i_KM_adj <- c()
    i_KM_redund <- c()
    i_KM_redund_bis <- c()
    i_KMgnm_gapfillKOs <- c()
    
    # get KO and KM list
    i_gnm <- indata[gnm, , drop=F]
    i_gnmKO <- colnames(i_gnm)[i_gnm > 0]
    
    i_gnmKM <- KMgnm_compl3[gnm, , drop=F]
    i_gnmKM <- colnames(i_gnmKM)[i_gnmKM > 0]
    
    # get genome KO specificity
    # NB!!! scores depends only on KEGG modules that are complete in the i-genome
    # i.e. score KOs based on the number of times that each KO appears in a KM: 1 in only 1 KM, 0.2 in 5 KM
    # NB2!! KO-spec are adjusted by multipling the value per KO occurence in the i-genome
    # 0 will be assigned when KO is absent in a specific genome
    # NB3!! missing KOs from KM-solutions that are considered complete (see criteria L532-537) will be stored
    
    for (KM in 1:length(KM_str)) { # loop thru KM
      #1 look if KM i present in genome
      if (!names(KM_str[KM]) %in% i_gnmKM) { # when KM is missing
        i_KM_adj <- rbind.data.frame(i_KM_adj, c(names(KM_str[KM]), 0))
        i_KM_redund <- rbind.data.frame(i_KM_redund, c(names(KM_str[KM]), 0))
        i_KM_redund_bis <- rbind.data.frame(i_KM_redund_bis, c(names(KM_str[KM]), 0))
        
      } else { # when KM is present
        i_KM_str <- KM_str[[KM]]
        
        #0 set bonus for gap-filling
        if (i_KM_str$length < breaks[1]) { GF_bonus <- allowed_gaps[1]
        } else if (i_KM_str$length < breaks[2]) { GF_bonus <- allowed_gaps[2]
        } else if (i_KM_str$length < breaks[3]) { GF_bonus <- allowed_gaps[3]
        } else if (i_KM_str$length < breaks[4]) { GF_bonus <- allowed_gaps[4]
        } else if (i_KM_str$length >= breaks[4]) { GF_bonus <- allowed_gaps[5] }
        #---
        
        #1 look which KOs are present in the genome
        i_KM_str$str$exist <- i_KM_str$str$KO %in% i_gnmKO
        #1b estimate KO-spec
        i_gnm_KOspec <- table(unlist(lapply(KM_str[names(KM_str) %in% i_gnmKM], function(i) as.character(unique(i$str$KO)))))
        i_gnm_KOspec <- structure(as.integer(i_gnm_KOspec), names=names(i_gnm_KOspec))
        i_gnm_KOspec2 <- rep(1, length(i_gnm_KOspec))/i_gnm_KOspec
        i_KM_str$str$KO_spec <- i_gnm_KOspec2[match(i_KM_str$str$KO, names(i_gnm_KOspec2))]
        #1c get KO-spec-adj
        i_KM_str$str$KO_spec_adj <- i_KM_str$str$KO_spec # when KO is missing leave global value of KO-spec
        i_KM_str$str$KO_spec_adj[i_KM_str$str$exist == T] <- as.numeric(i_KM_str$str$KO_spec_adj[i_KM_str$str$exist == T] *
                                                                                  i_gnm[, match(i_KM_str$str$KO[i_KM_str$str$exist == T], colnames(i_gnm))])
        #1d trim adj-KO-spec at 1
        i_KM_str$str$KO_spec_adj[i_KM_str$str$KO_spec_adj > 1] <- 1
        #1e fix
        i_KM_str$str[, c(4:7,9,10)] <- apply(i_KM_str$str[, c(4:7,9,10)], 2, function(x) as.numeric(x)) # fix
        
        #2 per each KM level get completness
        if (i_KM_str$just_macro == T) { # when single-level macromolecular complex
          #!! get immediately KM specificity
          i_KM_adj <- rbind.data.frame(i_KM_adj, c(names(KM_str[KM]), sum(i_KM_str$str$KO_spec_adj)/i_KM_str$length_macro))
          i_KM_redund <- rbind.data.frame(i_KM_redund, c(names(KM_str[KM]), 1))
          i_KM_redund_bis <- rbind.data.frame(i_KM_redund_bis, c(names(KM_str[KM]), 1))
          # gap-fill missing KOs if there are
          if (any(i_KM_str$str$exist == F)) { i_KMgnm_gapfillKOs <- append(i_KMgnm_gapfillKOs, as.character(i_KM_str$str$KO[i_KM_str$str$exist == F])) }
          
        } else { # when it is not a single-level macromolecular complex
          # gap-fill single-KO levels
          for (i_lvl in unique(i_KM_str$str$KM_lvl)) {
            if (nrow(i_KM_str$str[i_KM_str$str$KM_lvl == i_lvl, , drop=F]) == 1) {
              if (i_KM_str$str$exist[i_KM_str$str$KM_lvl == i_lvl] == F) {
                i_KM_str$str$exist[i_KM_str$str$KM_lvl == i_lvl] <- T
                i_KMgnm_gapfillKOs <- append(i_KMgnm_gapfillKOs, as.character(i_KM_str$str$KO[i_KM_str$str$KM_lvl == i_lvl]))
                GF_bonus <- GF_bonus-1
              }}}
          if (GF_bonus < 0) print(paste0("Error in ", names(KM_str[KM]), " genomes ", row.names(KMgnm_compl3)[gnm],
                                         " - KM was called complete with too many gaps")) # control line
          
          i_KM_adjKO <- c()
          i_KM_redundKO <- c()
          i_KM_redundKO_bis <- c()
          for (i_lvl in unique(i_KM_str$str$KM_lvl)) { # loop through levels of that KM
            i_lvl_KM_str <- i_KM_str$str[i_KM_str$str$KM_lvl %in% i_lvl, , drop=F]
            
            if (nrow(i_lvl_KM_str) == 1) { #!!! if level has single-KO
              if (i_lvl_KM_str$exist == T) {
                i_KM_adjKO <- c(i_KM_adjKO, i_lvl_KM_str$KO_spec_adj)
                i_KM_redundKO <- c(i_KM_redundKO, 1)
                if (length(i_gnm[, colnames(i_gnm) == i_lvl_KM_str$KO]) > 0) {
                  i_KM_redundKO_bis <- c(i_KM_redundKO_bis, 1*i_gnm[, colnames(i_gnm) == i_lvl_KM_str$KO])
                } else { i_KM_redundKO_bis <- c(i_KM_redundKO_bis, 1) } # 1 which has been gap-filled
              } else {
                i_KM_adjKO <- c(i_KM_adjKO, 0)
                i_KM_redundKO <- c(i_KM_redundKO, 0)
                i_KM_redundKO_bis <- c(i_KM_redundKO_bis, 0)
              }
              
            } else { #!!! if level has multiple-KOs
              #2a    find all possible X-wise solutions (looks for non-overlapping stacks of KOs)
              sol_starts <- as.data.frame(reduce(IRanges(i_lvl_KM_str$X1, i_lvl_KM_str$X2)))
              sol_x2_lim <- sol_starts$end
              #2b    check for correct width of the gap
              if (length(sol_x2_lim) > 1) {
                wrong_gap <- unlist(sapply(2:nrow(sol_starts), function(x) if(sol_starts$start[x] - sol_starts$end[x-1] != 4) return(x-1))) #!! fix fake gap that can be with nested redundant level of different length (e.g. M00003)
                if (length(wrong_gap) > 0) { sol_x2_lim <- sol_x2_lim[-wrong_gap] }
              }
              #2c    split level in all possible solutions - X-wise
              i_lvl_KM_str2 <- list()
              for (l in 1:length(sol_x2_lim)) {
                i_lvl_KM_str2[[l]] <- i_lvl_KM_str[i_lvl_KM_str$X2 <= sol_x2_lim[l], ]
                i_lvl_KM_str <- i_lvl_KM_str[!i_lvl_KM_str$X2 <= sol_x2_lim[l], ] # remove parsed element from level
              }
              #2d    split each solution Y-wise
              for (l in 1:length(i_lvl_KM_str2)) {
                i_lvl_KM_str2[[l]] <- split(i_lvl_KM_str2[[l]], f = as.factor(i_lvl_KM_str2[[l]]$Y1))
              }
              #2e    pull out all possible KO-combinations from each solution
              sol_paths_adj <- list()
              sol_paths_exist <- list()
              sol_paths_KOs <- list()
              for (l in 1:length(i_lvl_KM_str2)) { # loop thru solutions - for specificity
                sol_paths_adj[[l]] <- expand.grid(lapply(i_lvl_KM_str2[[l]], function(z) z$KO_spec_adj)) # use z$KO for debug
                sol_paths_exist[[l]] <- expand.grid(lapply(i_lvl_KM_str2[[l]], function(z) z$exist)) # use z$KO for debug
                sol_paths_KOs[[l]] <- expand.grid(lapply(i_lvl_KM_str2[[l]], function(z) z$KO))
              }
              #2f    gap-fill multi-KO levels
              # if no solution is complete, find solution with more retrieved KOs and gap-fill it
              best_solset <- which.max(lapply(sol_paths_exist, function(x) max(apply(x, 1, sum))))
              if (!any(apply(sol_paths_exist[[best_solset]], 1, sum) == ncol(sol_paths_exist[[best_solset]]))) { # when none is complete
                best_sol <- which.max(apply(sol_paths_exist[[best_solset]], 1, sum)) # find the best complete solution
                if (sum(!sol_paths_exist[[best_solset]][best_sol, ]) <= GF_bonus) { # if is still possible to gap-fill
                  sol_paths_exist[[best_solset]][best_sol, which(sol_paths_exist[[best_solset]][best_sol, ] == F)] <- T # gap-fill the best
                  i_KMgnm_gapfillKOs <- append(i_KMgnm_gapfillKOs, as.character(sol_paths_KOs[[best_solset]][best_sol, which(sol_paths_exist[[best_solset]][best_sol, ] == F)]))
                  GF_bonus <- GF_bonus-sum(!sol_paths_exist[[best_solset]][best_sol, ])
                }
              }
              #2g    discard solutions with still missing KOs
              sol_paths_adj2 <- list()
              sol_paths_KOs2 <- list()
              for (l in 1:length(sol_paths_exist)) {
                tmp_exist <- sol_paths_exist[[l]][apply(sol_paths_exist[[l]], 1, function(y) sum(as.logical(y)) == length(y)), , drop=F]
                if (length(row.names(tmp_exist)) > 0) { # keep only sol without 0
                  tmp_adj <- sol_paths_adj[[l]][apply(sol_paths_exist[[l]], 1, function(y) sum(as.logical(y)) == length(y)), , drop=F]
                  tmp_KOs <- sol_paths_KOs[[l]][apply(sol_paths_exist[[l]], 1, function(y) sum(as.logical(y)) == length(y)), , drop=F]
                  sol_paths_adj2 <- c(sol_paths_adj2, list(tmp_adj))
                  sol_paths_KOs2 <- c(sol_paths_KOs2, list(tmp_KOs))
                }
              }
              
              if (length(unlist(sol_paths_adj2)) > 0) {
                #2h    calculate the adj-KO-spec average for each solution
                sol_paths_adj3 <- lapply(sol_paths_adj2, function(x) rowMeans(as.matrix(x)))
                #2i    take the solution with the highest average adj-KO-spec
                sol_best_adj <- max(unlist(sol_paths_adj3))
                sol_best_redund <- length((unlist(sol_paths_adj3)))
                sol_best_redund_bis <- sum(unlist(lapply(sol_paths_KOs2, function(z)
                  apply(z, 1, function(zz) { if(any(colnames(i_gnm) %in% zz)) { prod(i_gnm[, colnames(i_gnm) %in% zz])
                  } else { 1 } })))) # 1 which has been gap-filled
              } else {
                sol_best_adj <- 0
                sol_best_redund <- 0
                sol_best_redund_bis <- 0
              }
              
              i_KM_adjKO <- c(i_KM_adjKO, sol_best_adj)
              i_KM_redundKO <- c(i_KM_redundKO, sol_best_redund)
              i_KM_redundKO_bis <- c(i_KM_redundKO_bis, sol_best_redund_bis)
            }
          }
          
          #3 sum adj-KO-spec of all KM levels and standardyze by KM length
          i_KM_adj <- rbind.data.frame(i_KM_adj, c(names(KM_str[KM]), sum(i_KM_adjKO)/i_KM_str$length))
          
          #3_bis fix redundancy based on completness rules for missing gene
          if (length(i_KM_redundKO) < breaks[1]) { i_KM_redundKO <- i_KM_redundKO; i_KM_redundKO_bis <- i_KM_redundKO_bis
          } else if (length(i_KM_redundKO) < breaks[2]) { i_KM_redundKO[which(i_KM_redundKO == 0)[allowed_gaps[2]]] <- 1; i_KM_redundKO_bis[which(i_KM_redundKO_bis == 0)[allowed_gaps[2]]] <- 1
          } else if (length(i_KM_redundKO) < breaks[3]) { i_KM_redundKO[which(i_KM_redundKO == 0)[1:allowed_gaps[3]]] <- 1; i_KM_redundKO_bis[which(i_KM_redundKO_bis == 0)[1:allowed_gaps[3]]] <- 1
          } else if (length(i_KM_redundKO) < breaks[4]) { i_KM_redundKO[which(i_KM_redundKO == 0)[1:allowed_gaps[4]]] <- 1; i_KM_redundKO_bis[which(i_KM_redundKO_bis == 0)[1:allowed_gaps[4]]] <- 1
          } else if (length(i_KM_redundKO) >= breaks[4]) { i_KM_redundKO[which(i_KM_redundKO == 0)[1:allowed_gaps[5]]] <- 1; i_KM_redundKO_bis[which(i_KM_redundKO_bis == 0)[1:allowed_gaps[5]]] <- 1 }
          i_KM_redund <- rbind.data.frame(i_KM_redund, c(names(KM_str[KM]), prod(i_KM_redundKO)))
          i_KM_redund_bis <- rbind.data.frame(i_KM_redund_bis, c(names(KM_str[KM]), prod(i_KM_redundKO_bis)))
        }
      }
    }
    
    names(i_KM_adj) <- c("KM", row.names(indata)[gnm])
    KMgnm_spec <- merge(KMgnm_spec, i_KM_adj, by.x="KM", by.y="KM", all = T)
    names(i_KM_redund) <- c("KM", row.names(indata)[gnm])
    KMgnm_redund <- merge(KMgnm_redund, i_KM_redund, by.x="KM", by.y="KM", all = T)
    names(i_KM_redund_bis) <- c("KM", row.names(indata)[gnm])
    KMgnm_redund_bis <- merge(KMgnm_redund_bis, i_KM_redund_bis, by.x="KM", by.y="KM", all = T)
    KMgnm_gapfillKOs[[row.names(KMgnm_compl3)[gnm]]] <- i_KMgnm_gapfillKOs
  }
  
  KMgnm_spec2 <- t(apply(KMgnm_spec[, -1], 2, as.numeric)); colnames(KMgnm_spec2) <- KMgnm_spec[, 1]
  KMgnm_spec2 <- KMgnm_spec2[, match(colnames(KMgnm_compl3), colnames(KMgnm_spec2))]
  KMgnm_redund2 <- t(apply(KMgnm_redund[, -1], 2, as.numeric)); colnames(KMgnm_redund2) <- KMgnm_redund[, 1]
  KMgnm_redund2 <- KMgnm_redund2[, match(colnames(KMgnm_compl3), colnames(KMgnm_redund2))]
  KMgnm_redund_bis2 <- t(apply(KMgnm_redund_bis[, -1], 2, as.numeric)); colnames(KMgnm_redund_bis2) <- KMgnm_redund_bis[, 1]
  KMgnm_redund_bis2 <- KMgnm_redund_bis2[, match(colnames(KMgnm_compl3), colnames(KMgnm_redund_bis2))]
  
  ## prepare output
  # object in list are data.frame with genomes as rows and KMs as columns
  tmp <- list(KMgnm_spec2,                    # specificity of KMs (1: each KO of a KM are present only in that KM; <1: increase level of shareness of the KOs of that KM
              KMgnm_redund2,                  # how many solution does a KM have in each genome
              KMgnm_redund_bis2,              # how many solution does a KM have in each genome (multiple copy of a KO will be counted as extra KO for that KM level)
              KMgnm_gapfillKOs)               # KOs added to gapfill KMs
  names(tmp) <- c("KM_specificity", "KM_#sol", "KM_#sol+multicopy", "gapfilled_KOs")
  return(tmp)              
}


