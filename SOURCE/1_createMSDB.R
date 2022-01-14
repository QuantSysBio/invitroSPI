### invitroSPI ###
# description:  parse search result files and create MSDB
# input:        sample list
# output:       MSDB
# authors:      GS, JL, HPR

library(plyr)
library(dplyr)
library(seqinr)
library(stringr)

source("SOURCE/loadFunctions.r")

print("--------------------------------------------")
print("1) PARSE SEARCH RESULT FILES AND CREATE MSDB")
print("--------------------------------------------")

project_name = snakemake@params[["project_name"]]
include_scanNum = snakemake@params[["include_scanNum"]]

SEARCHRESULTS_PATH = "INPUT/search_results/"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       sep = ",", header = T, stringsAsFactors = F)

# sample_list = read.table(file = "INPUT/sample_list.csv",
#                          sep = ",", header = T, stringsAsFactors = F)

# filter sample list
sample_list = sample_list[sample_list$project_name == project_name, ]

# add meta-information
if (!any(is.na(sample_list$metainformation)) & length(unique(sample_list$metainformation)) == 1) {
  metainfo = read.csv(sample_list$metainformation[1], stringsAsFactors = F)
  meta = T
} else {
  meta = F
}

### MAIN PART ###

if (!dir.exists(paste0("OUTPUT/", project_name))) {
  dir.create(paste0("OUTPUT/", project_name))
}

columns = c(            "runID" ,
                        "substrateID" ,
                        "substrateSeq" ,
                        "digestTime" ,
                        "pepSeq" ,
                        "productType",
                        "spliceType",
                        "positions",
                        "scanNum",
                        "rank",
                        "ionScore",
                        "qValue",
                        "charge",
                        "PTM")

MSDB = list()
pb = txtProgressBar(min = 0, max = nrow(sample_list), style = 3)

for (i in 1:nrow(sample_list)) {
  # setTxtProgressBar(pb, i)
  
  # load search result file
  con = file(paste0(SEARCHRESULTS_PATH, sample_list$filename[i]),"r")
  i1 = readLines(con)
  close(con)
  
  # remove header from search result file
  # and remove not assigned queries if there are any
  start = grep("prot_hit_num",i1)
  
  currentSearchFile = read.table(paste(SEARCHRESULTS_PATH, sample_list$filename[i],sep="/"),
                                 skip = start-1,
                                 sep = ",", header = TRUE, fill = TRUE)
  
  kk = which(currentSearchFile$prot_hit_num == "")
  if (length(kk) > 0) {
    kk = min(kk)
    currentSearchFile = currentSearchFile[c(1:(kk-1)),]
  }
  
  substrateID = sample_list$substrateID[i]
  digestTime = sample_list$digestTime[i]
  
  
  if(nrow(currentSearchFile) > 0){
    
    # empty DB for current search result file
    currentDB = data.frame(matrix(ncol = length(columns), nrow = nrow(currentSearchFile)))
    colnames(currentDB) = columns
    
    # fill data frame with info straight outta .csv
    currentDB$charge <- currentSearchFile$pep_exp_z
    currentDB$rank <- currentSearchFile$pep_rank
    currentDB$pepSeq <- currentSearchFile$pep_seq
    currentDB$ionScore <- currentSearchFile$pep_score
    currentDB$qValue <- currentSearchFile$pep_expect
    
    
    currentDB$PTM <- currentSearchFile$pep_var_mod
    currentDB$scanNum <- currentSearchFile$pep_query
    
    # currentDB[which(duplicated(currentDB) | duplicated(currentDB, fromLast = T)),] %>% nrow()
    # cntS = currentSearchFile[, c("pep_exp_z","pep_rank","pep_seq", "pep_score", "pep_expect", "pep_var_mod", "pep_query")]
    # cntS[which(duplicated(cntS) | duplicated(cntS, fromLast = T)),] %>% nrow()
    
    #get actual scan numbers, not query numbers
    # varies between Mascot and Mascot Distiller
    if(include_scanNum != "no") {
      
      warn = F
      
      scans = rep(NA,dim(currentSearchFile)[1])
      for(ii in 1:dim(currentSearchFile)[1]){
        tit = currentSearchFile$pep_scan_title[ii]
        
        if (str_detect(tit, pattern = "scan=")) {
          
          scans[ii] = as.numeric(strsplit(strsplit(tit,split="scan=")[[1]][2],split="~")[[1]][1])
          
        } else if (str_detect(tit, pattern = "Scan ")) {
          
          scans[ii] = as.numeric(unlist(strsplit(tit, " "))[3])
          
        } else if (str_detect(tit, pattern = coll(".")) | str_detect(tit, pattern = "TSN")) {
          scans[ii] = as.numeric(strsplit(tit,split="\\.")[[1]][2])
          
        } else {
          warn = T
        }
        
      }
      currentDB$scanNum <- scans
      
      if (warn) {
        print("could not extract scan numbers - check search result formatting!")
        currentDB$scanNum = currentSearchFile$pep_query
      }
      
    }
    
    
    #extract product type and position information
    pepName <- as.character(currentSearchFile$prot_acc)
    pepType <- c()
    pepPos <- c()
    for (z in 1:length(pepName)){
      splitPepName <- unlist(strsplit(pepName[z], "_"))
      pepType <- c(pepType, splitPepName[1])
      splitPepName <- splitPepName[-1]
      pepPos <- c(pepPos, paste(splitPepName, collapse = "_"))
      
    }
    
    currentDB$productType <- pepType
    currentDB$positions <- pepPos
    
    # substrate sequence
    if (grepl(pattern = ".fasta", sample_list$substrateSeq[i])) {
      
      cntSeq = read.fasta(file = paste0("INPUT/sequences/", sample_list$substrateSeq[i]),
                          seqtype = "AA", strip.desc = T)[[1]] %>%
        unlist() %>%
        paste(collapse = "")
      
    } else {
      cntSeq = sample_list$substrateSeq[i]
    }
    
    # add metainfo
    currentDB$substrateID = substrateID
    currentDB$substrateSeq <- cntSeq
    currentDB$digestTime <- digestTime
    currentDB$runID = paste(substrateID, digestTime, sample_list$replicate[i], sep = "-")
    
    # add further meta-information if provided
    if (meta) {
      
      cntMeta = metainfo[metainfo$filename == sample_list$filename[i], ]
      currentDB = cbind(currentDB, cntMeta[rep(1,nrow(currentDB)),]) %>%
        as.data.frame()
      
      dp = which(duplicated(names(currentDB)))
      if (length(dp) > 0) { currentDB = currentDB[,-dp]}
      
    }
    
    # assign substrate hits as PCP
    substratehits = which(currentDB$substrateSeq == currentDB$pepSeq)
    if (length(substratehits) > 0) {
      currentDB$productType[substratehits] = "PCP"
      currentDB$positions[substratehits] = paste(1,
                                                 nchar(currentDB$substrateSeq[substratehits]),
                                                 sep = "_")
    }
    
    print("")
    paste0("obtain correct annotations for positions and splice types for run: ",
           currentDB$runID[1]) %>%
      print()
    
    currentDB = currentDB %>% mapping()
    currentDB$spliceType[currentDB$productType == "PCP"] = NA
    
  } else {
    currentDB = data.frame(matrix(ncol = length(columns), nrow = 0))
    colnames(currentDB) = columns
    print("No peptides listed in this file! Skipping...")
  }
  
  # append data frame for current search file to master data frame
  MSDB[[i]] = currentDB
  
}

MSDB = plyr::ldply(MSDB)

### OUTPUT ###
save(MSDB, file = unlist(snakemake@output[["MSDB"]]))

