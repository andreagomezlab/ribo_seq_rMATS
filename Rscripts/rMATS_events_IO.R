library(svglite)
source("src/rMATS_processing_events.R")


get_metatable <- function(cell_type='PV', timepoint='48hrs', treat="DOI", all=F) {
  if (treat == "DOI") {
    meta_table <- read.delim("files/90_samples_meta_data_DOI_fileID.csv", sep = ",")
  }else if (treat == "Psilocybin"){
    meta_table <- read.delim("files/28_samples_meta_data_Psilo_fileID.csv", sep = ",")
  }
  
  meta_table$Sample.ID = gsub("\\-", "", meta_table$Sample.ID)
  rownames(meta_table) <- meta_table$Sample.ID
  meta_table = meta_table[order(meta_table$CellType, meta_table$Treatment, meta_table$Timepoint, decreasing = T),]
  meta_table$Treatment[meta_table$Treatment == "Saline"] = "Ctrl"
  meta_table$Batch = paste0("B", meta_table$Batch)
  meta_table$Group = factor(paste(meta_table$Treatment,meta_table$Timepoint,sep="."))
  
  if (all) {
    return(meta_table)
  } else {
    meta_table = meta_table[meta_table$CellType %in% cell_type & meta_table$Timepoint %in% timepoint,]
    meta_table = droplevels(meta_table)
  }
  
  return(meta_table)
}


get_summary_pie_values <- function(loc_res, celltype="PV", timepoint="48hrs", treat="DOI") {
  loc_res_ = loc_res[loc_res$CellType == celltype & loc_res$Treatment %in% treat & loc_res$Timepoint %in% timepoint,]
  
  nall <- nFDR <- nPadj <- nPSI <- 0
  aux = NULL
  
  for (i in 1:nrow(loc_res_)) {
    events = read.delim(file = loc_res_$TotalEvents[i], sep = "\t")
    nall = nrow(events)
    nGene = length(unique(events$GeneID))
    aux = rbind(aux, data.frame(celltype = loc_res_$CellType[i], event_type = loc_res_$EventType[i], group = "Total", value=nall, nGene=nGene, e_down=NA, e_up=NA, g_down=NA, g_up=NA))
    events = calculate_median_SJC(events)
    events = calculate_median_IJC(events)
    events = mark_filter_median_SJC(events, 0)
    nGene = length(unique(events$GeneID[events$Median_passed == "yes"]))
    aux = rbind(aux, data.frame(celltype = loc_res_$CellType[i], event_type = loc_res_$EventType[i], group = "events_passed_lc", value=table(events$Median_passed)[2], nGene=nGene,e_down=NA, e_up=NA, g_down=NA, g_up=NA))
    events = read.delim(file = loc_res_$Filtered[i], sep = "\t")
    nrow(events)
    nPadj =nrow(events[events$pval.adjust <= 0.05,])
    nGene = length(unique(events$geneSymbol[events$pval.adjust <= 0.05]))
    down = table(events[events$pval.adjust <= 0.05,]$status)[1]
    up = table(events[events$pval.adjust <= 0.05,]$status)[2]
    aux = rbind(aux, data.frame(celltype = loc_res_$CellType[i], event_type = loc_res_$EventType[i], group = "Pval.adjust", value=nPadj, nGene=nGene, e_down=down, e_up=up, g_down=NA, g_up=NA))
    nPSI = nrow(events[events$psi_passed == "yes",])
    nGene = length(unique(events$geneSymbol[events$psi_passed == "yes"]))
    down = table(events[events$psi_passed == "yes",]$status)[1]
    up = table(events[events$psi_passed == "yes",]$status)[2]
    g_down = length(unique(events$GeneID[events$psi_passed == "yes" & events$status == "CTRL"]))
    g_up = length(unique(events$GeneID[events$psi_passed == "yes" & events$status == "TREAT"]))
    if (celltype == c("Rbp4") & loc_res_$Timepoint[i] %in% c("1week", "1month") || celltype == "Cux2ERT" & loc_res_$Timepoint[i] == "1month"){
      nPSI = nrow(events[events$psi_passed == "yes" & events$non_corrected_overlap == "no",])
      nGene = length(unique(events$geneSymbol[events$psi_passed == "yes" & events$non_corrected_overlap == "no"]))
      down = table(events[events$psi_passed == "yes" & events$non_corrected_overlap == "no",]$status)[1]
      up = table(events[events$psi_passed == "yes"& events$non_corrected_overlap == "no",]$status)[2]
      g_down = length(unique(events$GeneID[events$psi_passed == "yes" & events$non_corrected_overlap == "no" & events$status == "CTRL"]))
      g_up = length(unique(events$GeneID[events$psi_passed == "yes" & events$non_corrected_overlap == "no" & events$status == "TREAT"]))
    }
    aux = rbind(aux, data.frame(celltype = loc_res_$CellType[i], event_type = loc_res_$EventType[i], group = "PSI.filter", value=nPSI, nGene=nGene, e_down=down, e_up=up, g_down=g_down, g_up=g_up))
  }
  
  return(aux)
}

get_combine_gene_events <- function(loc_res, celltype="PV", treat="DOI", time="48hrs", batch=FALSE) {
  loc_res_ = loc_res[loc_res$CellType == celltype & loc_res$Treatment == treat & loc_res$Timepoint == time,]
  genes_events = NULL
  for (i in 1:nrow(loc_res_)) {
    events = read.delim(file = loc_res_$Filtered[i], sep = "\t")
    events_ = events[events$psi_passed == "yes",]
    if (batch) {
      events_ = events[events$psi_passed == "yes" & events$non_corrected_overlap == "no",]
    }
    print(paste0(length(unique(events_$geneSymbol)), " unique genes"))
    genes_events = unique(c(genes_events, events_$geneSymbol))
  }
  print(paste0(length(genes_events), " total unique genes"))
  return(genes_events)
}

get_rMATS_results <- function(rmats_output_path, count_type="JC", event_type="SE", fdr=0.05) {
  GRCm39_table = read.delim("files/sequence_report.tsv", sep = "\t")
  mats = NULL
  
  if (count_type == "JC" | count_type == "JCEC") {
    if (event_type == "SE" | event_type == "A5SS" | event_type == "A3SS" | event_type == "MXE" | event_type == "RI") {
      mats = read.delim(paste0(rmats_output_path, paste0(event_type, ".MATS.", count_type,".txt")))
      #mats = mats[mats$FDR <= fdr,]
      mats$avgPSI_1 <- sapply(strsplit(mats$IncLevel1, ","), function(z) mean(as.numeric(z)))
      mats$avgPSI_2 <- sapply(strsplit(mats$IncLevel2, ","), function(z) mean(as.numeric(z)))
      mats$status <-ifelse(mats$IncLevelDifference > 0, "CTRL", "TREAT")
    }
  }
  mats$chr = gsub("chr", "", mats$chr)
  mats$geneSymbol = mats$GeneID
  mats = mats[mats$chr %in% GRCm39_table$RefSeq.seq.accession,]
  mats = droplevels(mats)
  return(mats)
}


fix_gene_id <- function(rmats_output_path) {
  
  GRCm39_table = read.delim("files/sequence_report.tsv", sep = "\t")
  mats = NULL
  
  for (count_type in c("JC", "JCEC")) {
    for (event_type in c("SE", "A5SS", "A3SS", "MXE", "RI")) {
      mats = read.delim(paste0(rmats_output_path, paste0(event_type, ".MATS.", count_type,".txt")))
      
      mats$chr = gsub("chr", "", mats$chr)
      mats$geneSymbol = mats$GeneID
      #mats = mats[mats$chr %in% GRCm39_table$RefSeq.seq.accession,]
      #mats$chr <- as.factor(mats$chr)
      #levels(mats$chr) <- GRCm39_table$UCSC.style.name[match(levels(mats$chr), GRCm39_table$RefSeq.seq.accession)]
      #mats = droplevels(mats)
      save_file(mats, rmats_output_path, paste0(event_type, ".MATS.", count_type,"_processed.txt"))
    }
  }
}


read_file <- function(dir, h=TRUE, s="\t"){
  file = read.delim(dir, header = h, sep = s)
  return(file)
}

save_file <- function(obj, dir, file_name, row=F, sep="\t") {
  write.table(obj, 
              file = paste0(dir, file_name), 
              sep = sep, 
              row.names = row, 
              quote = FALSE)
}


save_GR_bed_file <- function(obj, dir, file_name, extra=TRUE) {
  obj_df <- as.data.frame(obj)
  if (extra) {
    obj_df = obj_df[, c('seqnames', 'start', 'end',  'name', 'status', 'strand')]
    colnames(obj_df) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand')
  }
  write.table(obj_df, file = paste0(dir, file_name), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}

export_svg <- function(obj, dir, file_name, h=8, w=11){
  print("Exporting svg file")
  filename = paste0(dir, file_name, '.svg')
  print(filename)
  svglite(filename, width = w, height = h)
  if (!is.null(obj))
    print(obj)
  dev.off()
}

get_name <- function(x){
  aux = unlist(strsplit(x, "/"))
  aux = aux[length(aux)]
  aux = gsub(".bam", "", aux)
  return(aux)
}

get_counts <- function(raw, col, names_txt) {
  splitdat = do.call("rbind", strsplit( raw[, colnames(raw) == col ], ","))
  splitdat = data.frame(apply(splitdat, 2, as.numeric))
  names(splitdat) = names_txt
  return(splitdat)
}

### get_sample names order from ctrl.txt and treat.txt 

get_rMATS_RAW_input <- function(rmats_output_path, count_type="JC", event_type="SE", meta_table, group = "SAMPLE_1") {
  raw = NULL
  t_samples = read.delim(paste0(rmats_output_path, "/treat.txt"), header = F, sep = ",")
  t_samples = as.character(apply(t_samples, 2, get_name))
  c_samples = read.delim(paste0(rmats_output_path, "/ctrl.txt"), header = F, sep = ",")  
  c_samples = as.character(apply(c_samples, 2, get_name))
  
  if (count_type == "JC" | count_type == "JCEC") {
    if (event_type == "SE" | event_type == "A5SS" | event_type == "A3SS" | event_type == "MXE" | event_type == "RI") {
      raw = read.delim(paste0(rmats_output_path, paste0(count_type, ".raw.input.", event_type,".txt")))
    }
  }
  print(head(raw))
  if (group == "SAMPLE_1") {
    counts  <- get_counts(raw, "IJC_SAMPLE_1", c_samples)
    #rownames(counts) <- raw$ID
    aux = get_counts(raw, "SJC_SAMPLE_1", c_samples)
    #rownames(counts) <- aux$ID
    counts  <- rbind(counts, aux)  
    rownames(counts) <- 1:nrow(counts)
    
  } else if (group == "SAMPLE_2") {
    counts  <- get_counts(raw, "IJC_SAMPLE_2", t_samples)
    #rownames(counts) <- raw$ID
    aux = get_counts(raw, "SJC_SAMPLE_2", t_samples)
    #rownames(counts) <- aux$ID
    counts  <- rbind(counts, aux)  
    rownames(counts) <- 1:nrow(counts)
  }
  
  return(counts)
}



write_rMATS_RAW_input <- function(counts, rmats_output_path, count_type, event_type) {
  
  t_samples = read.delim(paste0(rmats_output_path, "treat.txt"), header = F, sep = ",")
  t_samples = as.character(apply(t_samples, 2, get_name))
  c_samples = read.delim(paste0(rmats_output_path, "ctrl.txt"), header = F, sep = ",")  
  c_samples = as.character(apply(c_samples, 2, get_name))
  
  
  if (count_type == "JC" | count_type == "JCEC") {
    if (event_type == "SE" | event_type == "A5SS" | event_type == "A3SS" | event_type == "MXE" | event_type == "RI") {
      
      raw = read.delim(paste0(rmats_output_path, paste0(count_type, ".raw.input.", event_type,".txt")))
      file_name = paste0(count_type, ".raw.input.", event_type,"_copy.txt")
      write.table(raw, file = paste0(rmats_output_path, file_name), sep = '\t', col.names = T, row.names = FALSE, quote = FALSE)
      
      IJC_SAMPLE_1 = apply(counts[1:(nrow(counts)/2),c_samples], 1, paste, collapse=",")
      SJC_SAMPLE_1 = apply(counts[(nrow(counts)/2 + 1):nrow(counts), c_samples], 1, paste, collapse=",")
      IJC_SAMPLE_2 = apply(counts[1:(nrow(counts)/2),t_samples], 1, paste, collapse=",")
      SJC_SAMPLE_2 = apply(counts[(nrow(counts)/2 + 1):nrow(counts), t_samples], 1, paste, collapse=",")
      new_raw = cbind(ID=raw$ID, data.frame(IJC_SAMPLE_1), data.frame(SJC_SAMPLE_1), data.frame(IJC_SAMPLE_2), data.frame(SJC_SAMPLE_2), IncFormLen=raw$IncFormLen, SkipFormLen=raw$SkipFormLen)
      
      file_name = paste0(count_type, ".raw.input.", event_type,".txt")
      write.table(new_raw, file = paste0(rmats_output_path, file_name), sep = '\t', col.names = T, row.names = FALSE, quote = FALSE)
    }
  }
  
}

get_rMATS_filtered_results <- function(dir, file_name) {
  file = read.delim(paste0(dir, file_name), sep="\t")
  return(file)
}


export_event_files <- function(loc_res, dir, file_name) {
  e_type = c("SE", "A5SS", "A3SS", "MXE", "RI")
  
  for (i in seq(1, nrow(loc_res), 5)){
    comb_events = NULL
    n = 0
    celltype = loc_res$CellType[i]
    timep = loc_res$Timepoint[i]
    treat = loc_res$Treatment[i]
    while (n < 5) {
      index_f = i + n
      events = read.delim(file = loc_res_$Filtered[index_f], sep = "\t")
      ## Keep non-corrected??? 
      events = events[events$psi_passed == "yes",]
      if (celltype == c("Rbp4") & timep %in% c("1week", "1month") || celltype == "Cux2ERT" & timep == "1month"){
        events = events[events$psi_passed == "yes" & events$non_corrected_overlap == "no",]
      }
      ## add event type column
      events$EventType = e_type[n + 1]
      ## combine events
      comb_events = rbind(comb_events, events)
      n = n + 1
    }
    file_name = paste(celltype, timep, trat, "AS_events_combined", sep="_")
    save_file(obj = comb_events, dir = "results/rmats/Episode_V/rMATS_events_filtered/")
  }
}

