### invitroSPI ###
# description:  parse search result files and create MSDB
# input:        sample list
# output:       MSDB
# authors:      GS, JL, modified by HPR

library(plyr)
library(dplyr)
library(seqinr)
library(stringr)

print("--------------------------------------------")
print("1) PARSE SEARCH RESULT FILES AND CREATE MSDB")
print("--------------------------------------------")

project_name = snakemake@params[["project_name"]]
SEARCHRESULTS_PATH = "INPUT/search_results"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       sep = ",", header = T, stringsAsFactors = F)

# sample_list = read.table(file = "INPUT/sample_list.csv",
#                          sep = ";", header = T, stringsAsFactors = F)

# filter sample list
sample_list = sample_list[sample_list$project_name == project_name, ]

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
  setTxtProgressBar(pb, i)
  
  # empty DB for current search result file
  currentDB = data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(currentDB) = columns
  
  # extract header
  con = file(paste(SEARCHRESULTS_PATH, sample_list$filename[i],sep="/"), "r")
  header = character()
  
  while(TRUE){
    line = readLines(con, n = 1)
    
    # append lines
    header = c(header,line)
    
    # check when to stop reading
    x = grep("prot_hit",line)
    if(length(x)>0) {break}
  }
  
  close(con)
  header = header[1:length(header)-1]
  
  substrateID = sample_list$substrateID[i]
  digestTime = sample_list$digestTime[i]
  
  # read in the current search file but skip header
  try(currentSearchFile <- read.table(paste(SEARCHRESULTS_PATH, sample_list$filename[i],sep="/"),
                                      skip = length(header), sep = ",", header = TRUE, fill = TRUE))
  
  if(nrow(currentSearchFile) != 0){
    
    currentDB[nrow(currentDB)+ nrow(currentSearchFile),] <- NA  
    
    # fill data frame with info straight outta .csv
    currentDB$charge <- currentSearchFile$pep_exp_z
    currentDB$rank <- currentSearchFile$pep_rank
    currentDB$pepSeq <- currentSearchFile$pep_seq
    currentDB$ionScore <- currentSearchFile$pep_score
    currentDB$qValue <- currentSearchFile$pep_expect
    currentDB$PTM <- currentSearchFile$pep_var_mod
    currentDB$scanNum <- currentSearchFile$pep_query
    
    #get actual scan numbers, not query numbers
    # varies between Mascot and Mascot Distiller
    scans = rep(NA,dim(currentSearchFile)[1])
    for(ii in 1:dim(currentSearchFile)[1]){
      tit = currentSearchFile$pep_scan_title[ii]
      
      if (str_detect(tit, pattern = "scan=")) {
        
        scans[ii] = as.numeric(strsplit(strsplit(tit,split="scan=")[[1]][2],split="~")[[1]][1])
        
      } else if (str_detect(tit, pattern = "Scan ")) {
        
        scans[[i]] = as.numeric(unlist(strsplit(tit, " "))[3])
      }
      
      
    }
    currentDB$scanNum <- scans
    
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
    
    # Determine spliceTypes of PSPs
    pepPos <- strsplit(pepPos, "_")
    spliceType <- c()
    for (n in 1:length(pepPos)){
      
      currentPos <- as.numeric(pepPos[[n]])
      
      if(length(currentPos) == 2){
        
        spliceType <- c(spliceType, NA)
        
      } else if(length(currentPos) == 4){
        
        # Check if splice reactants 1 and 2 overlap
        overlap <- any(currentPos[3]:currentPos[4] %in% currentPos[1]:currentPos[2])
        
        if(overlap == TRUE){
          spliceType <- c(spliceType, "trans")
          
        }else if(currentPos[2] > currentPos[3]){
          spliceType <- c(spliceType, "revCis")
          
        }else{
          spliceType <- c(spliceType, "cis")
          
        }
      }else{
        spliceType <- c(spliceType, NA)
      }
    }
    
    currentDB$spliceType <- spliceType
    
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
    
  } else {
    print("No peptides listed in this file! Skipping...")
  }
  
  # append data frame for current search file to master data frame
  MSDB[[i]] = currentDB
  
}

MSDB = plyr::ldply(MSDB)

### OUTPUT ###
save(MSDB, file = unlist(snakemake@output[["MSDB"]]))
