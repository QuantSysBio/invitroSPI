### IN VITRO SPI ###
# description:  functions for data preprocessing
# input:        -
# output:       -
# author:       HR


library(dplyr)
library(stringr)


########## plotting colours ########## 
plottingCols = c(
  PCP = "darkorange",
  PSP = "dodgerblue",
  cis = "darkslateblue",
  revCis = "lightskyblue",
  trans = "darkorchid",
  allcis = "darkslateblue",
  
  PaesM = "gray",
  invitroSPI = "red",
  
  SpechtDB = "olivedrab1",
  wholeDB = "olivedrab4",
  randomDB = "lightsteelblue"
)



########## disentange rows containing multimappers ########## 
disentangleMultimappers.Type = function(DB) {
  
  print("DISENTANGLE MULTI-MAPPERS FOR PRODUCT TYPE")
  
  k = which(str_detect(DB$positions, coll(";")))
  
  if (length(k) > 0) {
    
    DB_mm = DB[k, ]
    DB_nomm = DB[-k, ]
    
    mm = list()
    
    pb = txtProgressBar(min = 0, max = nrow(DB_mm), style = 3)
    for (r in 1:nrow(DB_mm)) {
      
      setTxtProgressBar(pb, r)
      
      cnt_types = str_split(DB_mm$spliceType[r], coll(";"), simplify = T) %>%
        paste()
      
      cnt_pos = str_split(DB_mm$positions[r], coll(";"), simplify = T) %>%
        paste()
      
      cntDB = DB_mm[r, ]
      
      # all multi-mappers correspond to the same product type
      if (cnt_types %>% unique() %>% length() == 1) {
        cntDB$spliceType = cnt_types %>% unique()
        
        # in case of cis and trans --> keep cis
      } else if ( any("trans" %in% cnt_types) ) {
        
        x = which(cnt_types == "trans")
        cnt_types = cnt_types[-x]
        cnt_pos = cnt_pos[-x]
        
        if (cnt_types %>% unique() %>% length() == 1) {
          cntDB$spliceType = cnt_types %>% unique()
          
        } else {
          cntDB$spliceType = "type_multi-mapper"
        }
        
        cntDB$positions = cnt_pos %>% paste(collapse = ";")
        
      } else {
        cntDB$spliceType = "type_multi-mapper"
        cntDB$positions = cnt_pos %>% paste(collapse = ";")
      }
      
      mm[[r]] = cntDB
    }
    
    mm = plyr::ldply(mm)
    DB = rbind(DB_nomm, mm) %>%
      as.data.frame()
  }
  
  return(DB)
}


disentangleMultimappers.AA = function(DB, retSinglePos = T) {
  
  print("DISENTANGLE MULTI-MAPPERS FOR AA AT sP1 AND sP1'")
  
  DB$AAmultimapper = "no"
  
  k = which(str_detect(DB$positions, coll(";")))
  
  if (length(k) > 0) {
    DB_mm = DB[k, ]
    DB_nomm = DB[-k, ]
    
    pb = txtProgressBar(min = 0, max = nrow(DB_mm), style = 3)
    for (r in 1:nrow(DB_mm)) {
      
      setTxtProgressBar(pb, r)
      
      cnt_pos = strsplit(DB_mm$positions[r], ";") %>%
        unlist() %>%
        str_split_fixed(pattern = coll("_"), n = Inf)
      
      sP1 = str_sub(DB_mm$substrateSeq[r], start = cnt_pos[, 2], end = cnt_pos[, 2])
      sP1d = str_sub(DB_mm$substrateSeq[r], start = cnt_pos[, 3], end = cnt_pos[, 3])
      
      if ((length(unique(sP1)) > 1) | (length(unique(sP1d)) > 1)) {
        DB_mm$AAmultimapper[r] = "yes"
        
      } else if (retSinglePos) {
        DB_mm$positions[r] = paste(cnt_pos[1, c(1:4)], collapse = "_")
      }
      
    }
    
    DB = rbind(DB_nomm, DB_mm) %>%
      as.data.frame()
  }
  
  
  return(DB)
}


disentangleMultimappers.SRlen = function(DB, retSinglePos = T) {
  
  print("DISENTANGLE MULTI-MAPPERS FOR SR length'")
  
  DB$SRmultimapper = "no"
  
  k = which(str_detect(DB$positions, coll(";")))
  
  if (length(k) > 0) {
    
    DB_mm = DB[k, ]
    DB_nomm = DB[-k, ]
    
    pb = txtProgressBar(min = 0, max = nrow(DB_mm), style = 3)
    for (r in 1:nrow(DB_mm)) {
      
      setTxtProgressBar(pb, r)
      
      cnt_pos = strsplit(DB_mm$positions[r], ";") %>%
        unlist() %>%
        str_split_fixed(pattern = coll("_"), n = Inf)
      
      SR1len = as.numeric(cnt_pos[, 2]) - as.numeric(cnt_pos[, 1]) + 1
      SR2len = as.numeric(cnt_pos[, 4]) - as.numeric(cnt_pos[, 3]) + 1
      
      if ((length(unique(SR1len)) > 1) | (length(unique(SR2len)) > 1)) {
        DB_mm$SRmultimapper[r] = "yes"
        
      }  else if (retSinglePos) {
        DB_mm$positions[r] = paste(cnt_pos[1, c(1:4)], collapse = "_")
      }
      
    }
    
    DB = rbind(DB_nomm, DB_mm) %>%
      as.data.frame()
    
    
  }
  
  
  return(DB)
}


removeMultimappers.Type = function(DB) {
  print("REMOVE PEPTIDES THAT CAN SOLELY BE MULTI-MAPPERS FOR PRODUCT TYPE")
  
  k = which(DB$spliceType == "type_multi-mapper")
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


removeMultimappers.AA = function(DB) {
  print("REMOVE PEPTIDES THAT ARE MULTI-MAPPERS IN TERMS OF AA AT sP1 OR sP1'")
  
  k = which(DB$AAmultimapper == "yes")
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


removeMultimappers.SRlen = function(DB) {
  print("REMOVE PEPTIDES THAT ARE MULTI-MAPPERS IN TERMS OF SR LENGTH")
  
  k = which(DB$SRmultimapper == "yes")
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


########## remove faulty substrate ########## 
removeTSN112 = function(DB) {
  print("REMOVE FAULTY SUBSTRATE TSN112")
  
  if ("TSN112" %in% DB$substrateID) {
    DB = DB[-which(DB$substrateID == "TSN112"), ]
  }
  
  return(DB)
}


########## remove substrate sequence from PCP hits########## 
removeSubstrateFromPCPs = function(DB) {
  print("REMOVE SUBSTRATE SEQUENCE FROM PCPs")
  
  k = which(DB$substrateSeq == DB$pepSeq & str_detect(DB$productType, "PCP"))
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}

########## filter for peptide length ########## 
filterPepLength = function(DB, cutoff=6) {
  print(paste0("REMOVE PEPTIDES THAT ARE SHORTER THAN ", cutoff, " AA"))
  
  k = which(nchar(DB$pepSeq) < cutoff)
  
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}

########## synthesis errors ########## 

remSynthErrors = function(DB) {
  print("REMOVE SYNTHESIS ERRORS")
  
  k = which(str_detect(DB$productType, "synError"))
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


keepSynthErrors = function(DB) {
  print("FILTER FOR SYNTHESIS ERRORS")
  
  k = which(str_detect(DB$productType, "synError"))
  if (length(k) > 0) {
    DB = DB[k, ]
  }
  
  return(DB)
}


########## early time points only ########## 
filterEarlyTimepoints = function(DB) {
  print("FILTER FOR EARLY TIME POINTS")
  
  DB = DB[which(DB$digestTime %in% c(2, 4)), ]
  
  return(DB)
}


########## 20S standard proteasome only ########## 
filter20Sstandard = function(DB) {
  print("FILTER FOR 20S STANDARD PROTEASOME")
  
  if ("protIsotype" %in% names(DB)) {
    DB = DB[DB$protIsotype %in% c("20S standard", "20S K562"), ]
  }
  
  return(DB)
}


########## remove peptides containing PTMS ########## 
remPTMs = function(DB) {
  print("REMOVE PTM CONTAINING PRODUCTS")
  
  k = which(! DB$PTM %in% c("", NA))
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}


########## no products containing substrate N or C ########## 
remSubstrateNorCterm = function(DB) {
  
  print("REMOVE PRODCUCTS THAT CONTAIN THE SUBSTRATE'S N OR C TERM")
  
  rem = c()
  
  pb = txtProgressBar(min = 0, max = nrow(DB), style = 3)
  for (i in 1:nrow(DB)) {
    
    setTxtProgressBar(pb, i)
    
    pos = strsplit(DB$positions[i], ";") %>%
      unlist() %>%
      str_split_fixed(pattern = "_", Inf)
    
    substrateLen = nchar(DB$substrateSeq[i])
    
    cntNCterm = rep(F, nrow(pos))
    
    for (j in 1:nrow(pos)) {
      
      if (str_detect(DB$productType[i], "PCP")) {
        if ((pos[j, 1] == 1) | (pos[j, 2] == substrateLen)) {
          cntNCterm[j] = T
        }
      } else if (str_detect(DB$productType[i], "PSP")) {
        if ((pos[j, 1] == 1) | (pos[j, 3] == 1) | (pos[j, 4] == substrateLen) | (pos[j, 2] == substrateLen)) {
          cntNCterm[j] = T
        }
      }
      
    }
    
    if (all(cntNCterm)) {
      rem = c(rem, i)
    }
    
  }
  
  if (length(rem) > 0) {
    DB = DB[-rem, ]
  }
  
  return(DB)
}


########## no products containing SR of 1 aa length ########## 
remSR1aa = function(DB) {
  
  print("REMOVE PRODCUCTS THAT CONTAIN A SR OF 1 AA LENGTH")
  
  rem = c()
  
  # PSPs with 1 aa SR
  for (i in which(str_detect(DB$productType, "PSP"))) {
    
    cnt_pos = strsplit(DB$positions[i], ";") %>%
      unlist() %>%
      str_split_fixed(pattern = "_", Inf)
    
    cnt_sr1 = as.numeric(cnt_pos[, 2]) - as.numeric(cnt_pos[, 1]) + 1
    cnt_sr2 = as.numeric(cnt_pos[, 4]) - as.numeric(cnt_pos[, 3]) + 1
    
    if ( all(cnt_sr1 == 1) | all(cnt_sr2 == 1) ) {
      rem = c(rem, i)
    }
  }
  
  if (length(rem) > 0) {
    DB = DB[-rem, ]
  }
  
  return(DB)
}

########## label products containing substrate's N/C term or a SR of 1 aa length ##########

filterTypes = function(DB) {
  
  DB.spliceTypes = DB$spliceType
  
  pb = txtProgressBar(min = 0, max = nrow(DB), style = 3)
  
  for(i in 1:nrow(DB)) {
    
    setTxtProgressBar(pb, i)
    
    pos = strsplit(DB$positions[i], ";") %>%
      unlist() %>%
      str_split_fixed(pattern = "_", Inf)
    
    if (all(pos == "")) {
      pos = str_locate(DB$substrateSeq[i], DB$pepSeq[i])
    }
    
    substrateLen = nchar(DB$substrateSeq[i])
    
    # N/C term containing PCPs and PSPs
    cntNCterm = rep(F, nrow(pos))
    
    for (j in 1:nrow(pos)) {
      
      if (str_detect(DB$productType[i], "PCP")) {
        
        if ((pos[j, 1] == 1) | (pos[j, 2] == substrateLen)) {
          cntNCterm[j] = T
        }
        
      } else if (str_detect(DB$productType[i], "PSP")) {
        
        if ((pos[j, 1] == 1) | (pos[j, 3] == 1) | (pos[j, 4] == substrateLen) | (pos[j, 2] == substrateLen)) {
          cntNCterm[j] = T
        }
        
      }
      
    }
    
    if (all(cntNCterm)) {
      DB.spliceTypes[i] = paste0(DB.spliceTypes[i], "_NCterm")
    } else {
      DB.spliceTypes[i] = DB.spliceTypes[i]
    }
    
    # PSPs with 1 aa SR
    if (str_detect(DB$productType[i], "PSP")) { 
      
      cnt_sr1 = as.numeric(pos[, 2]) - as.numeric(pos[, 1]) + 1
      cnt_sr2 = as.numeric(pos[, 4]) - as.numeric(pos[, 3]) + 1
      
      if ( all(cnt_sr1 == 1) | all(cnt_sr2 == 1) ) {
        DB.spliceTypes[i] = paste0(DB.spliceTypes[i], "_1aaSR")
      } else {
        DB.spliceTypes[i] = DB.spliceTypes[i]
      }
      
    } else {
      DB.spliceTypes[i] = DB.spliceTypes[i]
    }
    
  }
  
  DB$spliceType = DB.spliceTypes %>% unlist()
  
  return(DB)
}



########## I/L redundancy in product sequences #########
ILredundancy = function(DB) {
  
  print("REPLACE ALL ISOLEUCINS BY LEUCINS")
  
  DB$pepSeq = str_replace_all(DB$pepSeq, "I", "L")
  DB$substrateSeq = str_replace_all(DB$substrateSeq, "I", "L")
  
  return(DB)
}

ILredundancy.pepsOnly = function(DB) {
  
  print("REPLACE ALL ISOLEUCINS BY LEUCINS (only in products)")
  
  DB$pepSeq = str_replace_all(DB$pepSeq, "I", "L")
  
  return(DB)
}

########## unique product sequences ########## 
uniquePeptides = function(DB) {
  print("UNIQUE PEPTIDES")
  
  DB = DB[, c("substrateID", "substrateSeq", "pepSeq",
              "productType", "spliceType", "positions")] %>%
    unique()
  
  
  return(DB)
}


########## cosmetics ########## 
DBcosmetics = function(DB) {
  
  print("COSMETICS AND STATS")
  
  if ("X" %in% names(DB)) {
    DB$X = NULL
  }
  
  # replace NA in splice type column of PCPs
  DB$spliceType[str_detect(DB$productType, "PCP")] = "PCP"
  
  # statistics
  print(paste0("number of identified products: ", nrow(DB)))
  
  print(table(DB$spliceType))
  print(table(DB$spliceType) / nrow(DB))
  
  print("substrates")
  print(paste0("number of substrates: ", DB$substrateSeq %>% unique() %>% length()))
  DB$substrateID %>%
    unique() %>% 
    paste(collapse = ", ") %>%
    print()
  
  return(DB)
}


removeRowNames = function(DB) {
  if ("X" %in% names(DB)) {
    DB$X = NULL
  }
  
  return(DB)
}


########## load random DBs #########
loadRandom = function(fpath) {
  
  print("LOAD RANDOM DATABASES")
  
  # read cis, revCis, trans and PCP
  # concatenate and return list of data frames
  
  cis = read.csv(fpath[str_detect(fpath, "_cis")], stringsAsFactors = F)
  revCis = read.csv(fpath[str_detect(fpath, "_revCis")], stringsAsFactors = F)
  trans = read.csv(fpath[str_detect(fpath, "_trans")], stringsAsFactors = F)
  
  if (any(str_detect(fpath, "PCP_"))) {
    PCP = read.csv(fpath[str_detect(fpath, "_PCP")], stringsAsFactors = F)
  } else {
    PCP = NA
  }
  
  random = list(cis = cis,
                revCis = revCis,
                trans = trans,
                PCP = PCP)
  
  return(random)
}


