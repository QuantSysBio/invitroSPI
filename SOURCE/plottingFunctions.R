### invitroSPI ###
# description:  plotting functions for database diagnostics
# input:        -
# output:       -
# authors:      HPR


# ----- general stats ------
generalStats = function(DB, tp) {
  cntTBL = table(DB$spliceType) %>% as.data.frame()
  cntTBL = cbind(cntTBL, 
                 (table(DB$spliceType) / nrow(DB)) %>% as.data.frame())
  names(cntTBL) = c(paste0("time point: ", tp), "# peptides", "type", "frequency")
  cntTBL$frequency = round(cntTBL$frequency, 4)
  
  return(cntTBL)
}


# ----- number of peptides over time ----
getNumberofPeptides = function(ProteasomeDB) {
  
  DB_filtered = ProteasomeDB %>%
    filterPepLength() %>%
    disentangleMultimappers.Type() %>%
    DBcosmetics()
  
  # fraction of PSMs
  psms = DB_filtered %>%
    group_by(substrateID,spliceType, digestTime) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  psms_allTP = DB_filtered %>%
    group_by(substrateID,spliceType) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    mutate(digestTime = 0)
  
  # number + frequency of unique peptides
  peps = DB_filtered %>%
    uniquePeptides() %>%
    group_by(substrateID,spliceType, digestTime) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  peps_allTP = DB_filtered %>%
    uniquePeptides() %>%
    group_by(substrateID,spliceType) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    mutate(digestTime = 0)
  
  return(list(psms = rbind(psms, psms_allTP[, names(psms)]),
              peps = rbind(peps, peps_allTP[, names(peps)])))
  
}


plotNumberofPeptides = function(ProteasomeDB, outname) {
  
  res = getNumberofPeptides(ProteasomeDB)
  
  psms = res$psms
  peps = res$peps
  
  if (length(unique(ProteasomeDB$digestTime)) == 1) {
    psms = psms[psms$digestTime == 0, ]
    peps = peps[peps$digestTime == 0, ]
  }
  
  # number of PSMs over time
  psms$digestTime = factor(psms$digestTime, levels = sort(unique(as.numeric(psms$digestTime))))
  psms$spliceType = factor(psms$spliceType, levels = c("PCP", "cis", "revCis", "trans", "type_multi-mapper"))
  
  theme_set(theme_classic())
  
  a = ggplot(psms, aes(x=digestTime, y=n, col=spliceType)) +
    geom_boxplot(position=position_dodge(width=0.8)) +
    geom_point(alpha = .3, position=position_dodge(width=0.8)) +
    scale_y_continuous(trans='log10') +
    scale_x_discrete(labels = gsub("^0$","all",sort(unique(psms$digestTime)))) +
    scale_color_manual("product type",
                       values = c(plottingCols["PCP"], plottingCols["cis"], plottingCols["revCis"],
                                  plottingCols["trans"], "gray"),
                       labels = c("cleavage", "forward cis", "reverse cis", "trans", "multi-mapper")) +
    ggtitle("number of PSMs") +
    ylab("# PSMs") +
    xlab("digestion time")
  
  a2 = ggplot(psms, aes(x=digestTime, y=freq, col=spliceType)) +
    geom_boxplot(position=position_dodge(width=0.8)) +
    scale_x_discrete(labels = gsub("^0$","all",sort(unique(psms$digestTime)))) +
    scale_color_manual("product type",
                       values = c(plottingCols["PCP"], plottingCols["cis"], plottingCols["revCis"],
                                  plottingCols["trans"], "gray"),
                       labels = c("cleavage", "forward cis", "reverse cis", "trans", "multi-mapper")) +
    ylim(c(0,1)) +
    ggtitle("frequency of PSMs") +
    ylab("frequency") +
    xlab("digestion time")
  
  
  peps$digestTime = factor(peps$digestTime, levels = sort(unique(as.numeric(peps$digestTime))))
  peps$spliceType = factor(peps$spliceType, levels = c("PCP", "cis", "revCis", "trans", "type_multi-mapper"))
  
  b = ggplot(peps, aes(x=digestTime, y=n, col=spliceType)) +
    geom_boxplot(position=position_dodge(width=0.8)) +
    geom_point(alpha = .3, position=position_dodge(width=0.8)) +
    scale_y_continuous(trans='log10') +
    scale_x_discrete(labels = gsub("^0$","all",sort(unique(peps$digestTime)))) +
    scale_color_manual("product type",
                       values = c(plottingCols["PCP"], plottingCols["cis"], plottingCols["revCis"],
                                  plottingCols["trans"], "gray"),
                       labels = c("cleavage", "forward cis", "reverse cis", "trans", "multi-mapper")) +
    ggtitle("number of unique peptides") +
    ylab("# unique peptides") +
    xlab("digestion time")
  
  b2 = ggplot(peps, aes(x=digestTime, y=freq, col=spliceType)) +
    geom_boxplot(position=position_dodge(width=0.8)) +
    scale_x_discrete(labels = gsub("^0$","all",sort(unique(peps$digestTime)))) +
    scale_color_manual("product type",
                       values = c(plottingCols["PCP"], plottingCols["cis"], plottingCols["revCis"],
                                  plottingCols["trans"], "gray"),
                       labels = c("cleavage", "forward cis", "reverse cis", "trans", "multi-mapper")) +
    ylim(c(0,1)) +
    ggtitle("frequency of unique peptides") +
    ylab("frequency") +
    xlab("digestion time")
  
  p = grid.arrange(a,b, a2,b2, nrow=2, ncol=2)
  ggsave(filename = outname,
         plot = p, device = cairo_pdf(), height = 10, width = 14, dpi = "retina")
}


# ----- coverage maps ------

plotCoverage = function(ProtesomeDB, tp, name) {
  
  k = which(str_detect(ProteasomeDB$positions, ";"))
  if (length(k) > 0) {
    print("removing all multi-mappers")
    ProteasomeDB = ProteasomeDB[-k,]
  }
  
  pcp_pos = ProteasomeDB$positions[ProteasomeDB$productType == "PCP"]
  psp_pos = ProteasomeDB$positions[ProteasomeDB$productType == "PSP"]
  S = ProteasomeDB$substrateSeq[1]
  
  # --- PCPs ---
  pcp_pos = str_split(pcp_pos, "_", simplify = T)
  pcp_pos = apply(pcp_pos,2,as.numeric) %>%
    as.data.frame()
  pcp_pos = pcp_pos[order(pcp_pos$V1, decreasing = F), ]
  
  par(mfrow = c(2,1))
  y0 =0.1
  plot(NULL, ylim=c(0,nrow(pcp_pos)*0.05),
       ylab = "", xlab="substrate",
       xlim = c(0, nchar(S)),
       main = paste0(name, ": ", tp, " hrs"),
       axes=F)
  
  for (i in 1:nrow(pcp_pos)) {
    
    if (length(which(pcp_pos$V1[1:i] <= pcp_pos$V1[i] & 
                     pcp_pos$V2[1:i] >= pcp_pos$V1[i])) > 1) {
      
      cnt.intercept = length(which(pcp_pos$V1[1:i] <= pcp_pos$V1[i] & 
                                     pcp_pos$V2[1:i] >= pcp_pos$V1[i]))*0.1
      
    } else {cnt.intercept = y0}
    
    segments(x0 = pcp_pos$V1[i], y0 = cnt.intercept,
             x1 = pcp_pos$V2[i], y1 = cnt.intercept,
             col = plottingCols["PCP"], lwd = .5)
  }
  
  xlabel = c("0",
             paste(S %>% strsplit("") %>% unlist(), seq(1, nchar(S)), sep = ""))
  text(x = seq(0, nchar(S)), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=.5)
  axis(2, labels = NA)
  
  
  # --- PSPs ---
  psp_pos = str_split(psp_pos, "_", simplify = T)
  psp_pos = apply(psp_pos,2,as.numeric) %>%
    as.data.frame()
  psppos = psp_pos[order(psp_pos$V1, decreasing = F), ]
  
  psp_pos = rbind(as.matrix(psppos[,c(1:2)]), as.matrix(psppos[,c(3:4)])) %>%
    as.data.frame()
  
  y0 =0.1
  plot(NULL, ylim=c(0,nrow(psp_pos)*0.05),
       ylab = "", xlab="substrate",
       xlim = c(0, nchar(S)),
       main = paste0(name, ": ", tp, " hrs"),
       axes=F)
  
  for (i in 1:nrow(psp_pos)) {
    
    if (length(which(psp_pos$V1[1:i] <= psp_pos$V1[i] & 
                     psp_pos$V2[1:i] >= psp_pos$V1[i])) > 1) {
      
      cnt.intercept = length(which(psp_pos$V1[1:i] <= psp_pos$V1[i] & 
                                     psp_pos$V2[1:i] >= psp_pos$V1[i]))*0.1
      
    } else {cnt.intercept = y0}
    
    segments(x0 = psp_pos$V1[i], y0 = cnt.intercept,
             x1 = psp_pos$V2[i], y1 = cnt.intercept,
             col = plottingCols["PSP"], lwd = .5)
  }
  
  xlabel = c("0",
             paste(S %>% strsplit("") %>% unlist(), seq(1, nchar(S)), sep = ""))
  text(x = seq(0, nchar(S)), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=.5)
  axis(2, labels = NA)
  
  out = recordPlot()
  return(out)
}


# ----- length distributions -----
# function for violin plots
ViolinPlot = function(data) {
  
  # print stats
  s = data %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(n = dplyr::n(),
                     mean = mean(value),
                     median = median(value),
                     std = sd(value))
  print.data.frame(s)
  
  theme_set(theme_classic())
  
  # plotting
  data$type = factor(data$type, levels = c("PCP", "cis", "revCis", "trans"))
  
  sc = ggplot(data = data, aes(x = type, y = value, fill = type)) +
    geom_violin(size = .1) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = .2,
                 position = position_dodge(width = 4)) +
    xlab("product type") +
    scale_fill_manual("product type",
                      values = c(plottingCols["PCP"], plottingCols["cis"], plottingCols["revCis"],
                                 plottingCols["trans"]),
                      labels = c("cleavage", "forward cis", "reverse cis", "trans")) +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none")
  
  # add labs and title downstream
  
  return(sc)
}


# peptide length
PepLength = function(DB, tp) {
  
  print("PEPTIDE LENGTH")
  DB = removeMultimappers.Type(DB)
  
  pl = data.frame(value = nchar(DB$pepSeq),
                  type = DB$spliceType)
  
  pepLen = ViolinPlot(data = pl) +
    ylab("peptide product length (aa residues)") +
    scale_x_discrete(labels = c("cleavage", "forward cis", "reverse cis", "trans")) +
    # scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
    ggtitle("peptide length distribution",
            subtitle = paste0(tp, " hrs"))
  
  return(pepLen)
}


# SR length
SRLength = function(DB, tp) {
  
  print("SPLICE-REACTANT LENGTH")
  
  DB = DB[str_detect(DB$productType, "PSP"), ] %>%
    removeMultimappers.Type() %>%
    removeMultimappers.SRlen()
  
  srlen = function(db) {
    pos = str_split_fixed(db$positions, pattern = "_", n = Inf)
    
    sr = data.frame(value_sr1 = as.numeric(pos[, 2]) - as.numeric(pos[, 1]) + 1,
                    value_sr2 = as.numeric(pos[, 4]) - as.numeric(pos[, 3]) + 1,
                    type = db$spliceType)
    sr$type = factor(sr$type, levels = c("cis", "revCis", "trans"))
    return(sr)
  }
  
  sr = srlen(db = DB)
  
  sr1 = sr[, c("value_sr1", "type")]
  names(sr1)[1] = "value"
  sr2 = sr[, c("value_sr2", "type")]
  names(sr2)[1] = "value"
  
  sr1.plot = ViolinPlot(data = sr1) +
    ylab("splice reactant length (aa residues)") +
    scale_x_discrete(labels = c("forward cis", "reverse cis", "trans")) +
    # scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
    ggtitle("SR1 length distribution",
            subtitle = paste0(tp, " hrs"))
  
  sr2.plot = ViolinPlot(data = sr2) +
    ylab("splice reactant length (aa residues)") +
    scale_x_discrete(labels = c("forward cis", "reverse cis", "trans")) +
    # scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
    ggtitle("SR2 length distribution",
            subtitle = paste0(tp, " hrs"))
  
  sr.plot = arrangeGrob(grobs = list(sr1.plot, sr2.plot), ncol=2)
  
  return(sr.plot)
}


# intervening sequence length
IVSeqLength = function(DB, tp) {
  
  print("INTERVENING SEQUENCE LENGTH")
  
  DB = DB[str_detect(DB$productType, "PSP"), ] %>%
    removeMultimappers.Type() %>%
    removeMultimappers.IVSeqLen()
  
  ivlen = function(db) {
    # intervening sequence length
    # calculated as in Specht et al., Paes et al. and SPI-ART
    pos = str_split_fixed(db$positions, pattern = "_", n = Inf)
    iv = rep(NA, nrow(db))
    
    cistrans.idx = which(db$spliceType %in% c("cis", "trans"))
    iv[cistrans.idx] = (abs(as.numeric(pos[cistrans.idx, 3]) - as.numeric(pos[cistrans.idx, 2])) - 1) %>%
      as.numeric()
    
    revcis.idx = which(db$spliceType == "revCis")
    iv[revcis.idx] = (abs(as.numeric(pos[revcis.idx, 1]) - as.numeric(pos[revcis.idx, 4])) - 1) %>%
      as.numeric()
    
    ivl = data.frame(value = iv,
                     type = db$spliceType)
    return(ivl)
  }
  
  iv = ivlen(db = DB)
  iv = iv[iv$type != "trans", ]
  
  iv.plot = ViolinPlot(data = iv) +
    ylab("intervening sequence length (aa residues)") +
    scale_x_discrete(labels = c("forward cis", "reverse cis")) +
    ggtitle("intervening sequence length distribution",
            subtitle = paste0(tp, " hrs"))
  
  return(iv.plot)
}


# ------ type frequencies -----

TypeFrequency = function(DB) {
  
  print("FREQUENCY OF PRODUCT TYPES")
  
  Freq = DB %>% 
    removeMultimappers.Type() %>%
    group_by(spliceType) %>% 
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    arrange(desc(spliceType)) %>%
    mutate(ypos = cumsum(freq)- 0.5*freq )
  
  print(Freq)
  
  Freq$spliceType = factor(Freq$spliceType, levels = c("PCP", "cis", "revCis", "trans"))
  
  freqs = ggplot(Freq, aes(x="", y=freq, fill=spliceType)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    scale_fill_manual("product type",
                      values = c(plottingCols["PCP"], plottingCols["cis"], plottingCols["revCis"],
                                 plottingCols["trans"]),
                      labels = c("cleavage", "forward cis", "reverse cis", "trans")) +
    theme_void() +
    geom_text(aes(y = ypos, label=paste0("n = ", n)), size = 4)
  
  return(freqs)
}


