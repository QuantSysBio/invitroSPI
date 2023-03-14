### in vitro PCPS ###
# description:  calculate all possible spliced peptides from substrate
# input:        protein sequences
# output:       substrateID_allPSP.fasta (DB for Mascot search)
# authors:      HPR

library(data.table)
library(parallel)
library(dplyr)
library(seqinr)
library(stringr)
library(beepr)


# numCPU = if (detectCores() > 8) detectCores()-1 else if (detectCores() >= 4) 4 else 1
# cl <- makeForkCluster(numCPU)


# ----- user parameters -----
nmers = c(5, 40)
generateTrans = T  # if false, generate only cis peptides
outdir = "../DATA/polypeptideDBs/"
suppressWarnings(dir.create(outdir, recursive = T))

# change here (you can read in a dataframe with the name and the sequence)
sample_list = read.csv("~/Documents/SCIENCE/_tools/aSPIre+invitroSPI/data/sample_list_aSPIre.csv", stringsAsFactors = F)
substrates = sample_list %>%
  select(substrateID, substrateSeq) %>%
  unique()


# ----- functions -----

getCis = function(subSeq, N, Lext,L) {
  
  if (L-1 >= N) {
    
    allCis = sapply(seq(1,L-N), function(i){
      
      sapply(seq(i+Lext-1,N-Lext+i-1), function(j){
        
        sr2 = N-j+i-1
        sapply(seq(j+2, L-sr2+1), function(k){
          
          n = k+sr2-1
          pepSeq = paste(str_sub(subSeq, start=i, end=j),
                         str_sub(subSeq, start=k, end=n), sep = "")
          positions = paste(i,j,k,n, sep = "_")
          
          return(c(pepSeq, positions))
          
        })
      })
    })
    CIS = unlist(allCis)
    cisPep = data.frame(pepSeq = CIS[seq(1,length(CIS),2)],
                        PositionSubstrate = CIS[seq(2,length(CIS),2)])
  } else {
    cisPep=NA
  }
  
  return(cisPep)
}



getRevCis = function(subSeq, L, N, Lext=1) {
  
  
  if (L >= N) {
    
    allRevCis = sapply(seq(1,L-N+1), function(k){
      
      sapply(seq(k+Lext-1,N-Lext+k-1), function(n){
        
        sr1 = N-n+k-1
        sapply(seq(n+1, L-sr1+1), function(i){
          
          j = i+sr1-1
          
          pepSeq = paste(str_sub(subSeq, start=i, end=j),
                         str_sub(subSeq, start=k, end=n), sep = "")
          positions = paste(i,j,k,n, sep = "_")
          return(c(pepSeq, positions))
        })
      })
    })
    
    REVCIS = unlist(allRevCis)
    revcisPep = data.frame(pepSeq = REVCIS[seq(1,length(REVCIS),2)],
                           PositionSubstrate = REVCIS[seq(2,length(REVCIS),2)])
    
  } else {
    revcisPep = NA
  }
  
  
  return(revcisPep)
}




getAllTrans = function(L, N, subSeq, Lext=1) {
  
  allSRs = lapply(seq(Lext,N-Lext), function(sr1){
    
    # SR1
    Is = seq(1, L-sr1+1)
    Js = Is+sr1-1
    
    # SR2
    Ks = seq(1,L-N+sr1+1)
    Ns = Ks+N-sr1-1
    
    SR1_pos = paste(Is,Js,sep = "_")
    SR2_pos = paste(Ks,Ns,sep = "_")
    
    SR1 = str_sub(subSeq,Is,Js)
    SR2 = str_sub(subSeq,Ks,Ns)
    
    PEP = outer(SR1, SR2, paste, sep="") %>% as.vector()
    positions = outer(SR1_pos,SR2_pos,paste, sep = "_") %>% as.vector()
    
    return(data.frame(pepSeq = PEP,
                      PositionSubstrate = positions
    ))
    
  })
  
  transPep = plyr::ldply(allSRs)
  
  
  return(transPep)
}


getPCP = function(L,subSeq,N) {
 
    if (L > N) {
    
    allPCP = sapply(seq(1,L-N+1), function(i){
      j = N+i-1
      positions = paste(i,j, sep = "_")
      pepSeq=str_sub(subSeq,i,j)
      return(c(pepSeq,positions))
    })
    
    PCP = unlist(allPCP)
    PCP_Pep = data.frame(pepSeq = PCP[1,],
                         PositionSubstrate = PCP[2,])
  } else {
    PCP = NA
  }
  
  return(PCP_Pep)
}


# ----- compute DB -----

for (i in 1:nrow(substrates)) {
  
  
  Name = substrates$substrateID[i]
  print("...............")
  print(Name)
  print("...............")
  
  # get substrate sequence and length
  subSeq = substrates$substrateSeq[i]
  if (grepl(".fasta", subSeq)) {
    subSeq = read.fasta(file = subSeq, seqtype = "AA", strip.desc = T, as.string = T) %>%
      unlist() %>%
      suppressWarnings()
  }
  
  L = nchar(subSeq)
  
  # initiate empty fasta
  write.fasta(sequences = subSeq, names = Name,
              file.out=paste0(outdir, Name, "_allPSP_0.fasta"),
              open="w")
  
  start_time = Sys.time()
  # generate peptides
  if (generateTrans) {
    print("generating cis and trans-spliced peptides.....")
    
    # all cis and trans
    allPSP <- lapply(seq(nmers[1], nmers[2]), function(N){
      paste0("generating ",N,"-mers") %>% print()
      
      DB = getAllTrans(subSeq = subSeq, L = L, N = N, Lext=1)
      write.fasta(as.list(DB$pepSeq), names = paste0("PSP_", DB$PositionSubstrate),
                  file.out = paste0(outdir, Name, "_allPSP_0.fasta"),
                  open="a")
    })
    
    
    print("generating non-spliced peptides.....")
    Nmax = if(nmers[2] >= L-1) L-1 else nmers[2]
    
    allPCP <- lapply(seq(nmers[1], Nmax), function(N){
      paste0("generating ",N,"-mers") %>% print()
      
      DB = getPCP(subSeq = subSeq, L = L, N = N)
      write.fasta(as.list(DB$pepSeq), names = paste0("PSP_", DB$PositionSubstrate),
                  file.out = paste0(outdir, Name, "_allPSP_0.fasta"),
                  open="a")
    })
    
  } else {
    print("generating cis-spliced peptides.....")
    
    # forward cis
    allCis <- lapply(seq(nmers[1], nmers[2]), function(N){
      paste0("generating fwd cis ",N,"-mers") %>% print()
      
      DB = getCis(subSeq = subSeq, L = L, N = N, Lext=1)
      write.fasta(as.list(DB$pepSeq), names = paste0("PSP_", DB$PositionSubstrate),
                  file.out = paste0(outdir, Name, "_allPSP_0.fasta"),
                  open="a")
    })
    
    # reverse cis
    allRevCis <- lapply(seq(nmers[1], nmers[2]), function(N){
      paste0("generating rev cis ",N,"-mers") %>% print()
      
      DB = getRevCis(subSeq = subSeq, L = L, N = N, Lext=1)
      write.fasta(as.list(DB$pepSeq), names = paste0("PSP_", DB$PositionSubstrate),
                  file.out = paste0(outdir, Name, "_allPSP_0.fasta"),
                  open="a")
    })
    
    
  }
  
  end_time = Sys.time()
  paste0("..... finished after: ", difftime(end_time, start_time, units = "s")) %>% print()
  
}


print("FINISHED!")
beep("mario")

