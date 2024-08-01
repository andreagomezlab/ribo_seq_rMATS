get_rMATS_results <- function(rmats_output_path, count_type="JC", event_type="SE", fdr=0.05) {
  
  mats = NULL
  
  if (count_type == "JC" | count_type == "JCEC") {
    if (event_type == "SE" | event_type == "A5SS" | event_type == "A3SS" | event_type == "MXE" | event_type == "RI") {
      mats = read.delim(paste0(rmats_output_path, paste0(event_type, ".MATS.", count_type,".txt")))
      mats = mats[mats$FDR <= fdr,]
      mats$avgPSI_1 <- sapply(strsplit(mats$IncLevel1, ","), function(z) mean(as.numeric(z)))
      mats$avgPSI_2 <- sapply(strsplit(mats$IncLevel2, ","), function(z) mean(as.numeric(z)))
      mats$status <-ifelse(mats$IncLevelDifference > 0, "CTRL", "TREAT")
    }
  }
  
  mats$chr = gsub("chr", "", mats$chr)
  mats = mats[mats$chr %in% GRCm39_table$RefSeq.seq.accession,]
  mats = droplevels(mats)
  return(mats)
}

save_file <- function(obj, dir, file_name, row=F, sep="\t") {
  write.table(obj, file = paste0(dir, file_name), sep = sep, row.names = row, quote = FALSE)
}


save_GR_bed_file <- function(obj, dir, file_name, extra=TRUE) {
  obj_df <- as.data.frame(obj)
  if (extra) {
    obj_df = obj_df[, c('seqnames', 'start', 'end',  'name', 'status', 'strand')]
    colnames(obj_df) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand')
  }
  write.table(obj_df, file = paste0(dir, file_name), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}
