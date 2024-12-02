library(tidyverse)
library(hrbrthemes)
library(ggVennDiagram)
library("VennDiagram")
library(ggbio)
library(foreach)
library(doMC)

registerDoMC(cores = 4)


check_UCSC_seqnames <- function(obj) {
  
  GRCm39_table = read.delim("files/sequence_report.tsv", sep = "\t")
  
  if (grepl('NC_', levels(seqnames(obj))[1])) {
    obj <- obj[as.character(seqnames(obj)) %in% GRCm39_table$RefSeq.seq.accession, ]
    obj_df <- as.data.frame(obj)
    obj_df$seqnames <- plyr::mapvalues(x = obj_df$seqnames, from = GRCm39_table$RefSeq.seq.accession, to = GRCm39_table$UCSC.style.name)
    obj <- makeGRangesFromDataFrame(obj_df, keep.extra.columns = T)
  } 
  obj = keepStandardChromosomes(obj, pruning.mode="coarse")
}

build_coordinates_gr <- function(res, et, site='single') {
  
  if (et == "SE") {
    aux = res[, c("chr", "strand", "exonStart_0base", "exonEnd", "GeneID", "source", "IncLevelDifference", "status", "psi_passed")]
  } else if (et == "A5SS") {
    aux = res[, c("chr", "strand", "shortES", "shortEE", "GeneID", "source", "IncLevelDifference","status","psi_passed")]
  } else if (et == "A3SS") {
    aux = res[, c("chr", "strand", "shortES", "shortEE", "GeneID", "source", "IncLevelDifference","status","psi_passed")]
  } else if (et == "MXE") {
    aux = res[, c("chr", "strand", "X2ndExonStart_0base", "X2ndExonEnd", "GeneID", "source", "IncLevelDifference","status","psi_passed")]
  } else if (et == "RI") {
    aux = res[, c("chr", "strand", "riExonStart_0base", "riExonEnd", "GeneID", "source", "IncLevelDifference","status","psi_passed")]
  }
  colnames(aux) <- c("chr", "strand", "start", "end", "name", "source", "IncLevelDifference","status","psi_passed")
  aux$event_type = et
  aux_gr = makeGRangesFromDataFrame(aux, keep.extra.columns = T)
  aux_gr = check_UCSC_seqnames(aux_gr)
  
  return(aux_gr)
}

mark_filter_events <- function(res, psi=0.1){
  res$IncLevelDifference_passed = abs(res$IncLevelDifference) >= psi
  res$psi_passed = ifelse(res$IncLevelDifference_passed, "yes", "no")
  return(res)
}

get_AS_associated_genes <- function(res) {
  data_ = res[, c('study','source', 'GeneID')]
  return(data_[!duplicated(data_),])
}

get_overlapping_genes <- function(x) {
  t=get.venn.partitions(x, keep.elements = T, force.unique = T)
  return(t$..values..[[1]])
}

jaccard_index <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

overlap_coef <- function(a, b){
  intersection = length(intersect(a, b))
  den = min(length(a), length(b))
  return (intersection/den)
}

count_overlapping_events <- function(gr) {
  ## baseline 
  gr_baseline = reduce(gr)
  #save_GR_bed_file(gr_baseline, dir = "files/", file_name = "baseline_events.bed", extra=FALSE)
  gr_list_source <- split(x = gr, f = gr$source)
  
  overlaps <- foreach(i=1:length(gr_list_source), .combine = cbind)%dopar%{
    gr_source <- gr_list_source[[i]]
    overlap <- countOverlaps(query = gr_baseline, subject = gr_source, type = 'within') > 0
    return(overlap)
  }
  colnames(overlaps) <- names(gr_list_source)
  gr_baseline_df = as.data.frame(gr[findOverlaps(gr_baseline, gr, select = "first"),])
  count = rowSums(overlaps)
  aux = cbind(gr_baseline_df[, -which(colnames(gr_baseline_df) == 'source')], overlaps)
  aux$counts = count
  gr_overlapps = makeGRangesFromDataFrame(aux, keep.extra.columns = T)
  return(gr_overlapps)
}


export_sequence <- function() {
  
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome)
  
  Hsapiens <- BSgenome.
  
  ## random tx subset
  tx <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
  gr <- tx[sample(seq(length(tx)),5)]
  
  ## extract sequence
  seq <- getSeq(Hsapiens, gr)
  
  ## add names
  names(seq) <- gr$tx_name
  
  writeXStringSet(seq,file="result")
  
}

numeric_from_comma_string <- function(comma_string) {
  #print(comma_string)
  aux = as.numeric(unlist(strsplit(comma_string, split=',')))
  #print(aux)
  return(aux)
}

calculate_SD_from_row <- function(row) {
  avg_values <- numeric_from_comma_string(row)
  return(sd(avg_values))
}

calculate_avg_from_row <- function(row) {
  avg_values <- numeric_from_comma_string(row)
  return(mean(avg_values))
}

calculate_median_from_row <- function(row) {
  avg_values <- numeric_from_comma_string(row)
  return(median(avg_values))
}

calculate_median_SJC <- function(x) {
  x$median_SJC1 = apply(x, 1, function(s) calculate_median_from_row(s["SJC_SAMPLE_1"]))
  x$median_SJC2 = apply(x, 1, function(s) calculate_median_from_row(s["SJC_SAMPLE_2"]))
  return(x)
}

calculate_median_IJC <- function(x) {
  x$median_IJC1 = apply(x, 1, function(s) calculate_median_from_row(s["IJC_SAMPLE_1"]))
  x$median_IJC2 = apply(x, 1, function(s) calculate_median_from_row(s["IJC_SAMPLE_2"]))
  return(x)
}


calculate_sd <- function(x) {
  x$IncSD_1 = apply(x, 1, function(s) calculate_SD_from_row(s["IJC_SAMPLE_1"])) + apply(x, 1, function(s) calculate_SD_from_row(s["SJC_SAMPLE_1"]))
  x$IncSD_2 = apply(x, 1, function(s) calculate_SD_from_row(s["IJC_SAMPLE_2"])) + apply(x, 1, function(s) calculate_SD_from_row(s["SJC_SAMPLE_2"]))
  return(x)
}

calculate_avg <- function(x) {
  x$IncAVG_1 = apply(x, 1, function(s) calculate_avg_from_row(s["IJC_SAMPLE_1"])) + apply(x, 1, function(s) calculate_avg_from_row(s["SJC_SAMPLE_1"]))
  x$IncAVG_2 = apply(x, 1, function(s) calculate_avg_from_row(s["IJC_SAMPLE_2"])) + apply(x, 1, function(s) calculate_avg_from_row(s["SJC_SAMPLE_2"]))
  return(x)
}

min_psi_from_row <- function(row) {
  sample_1_psi_values <- numeric_from_comma_string(row["IncLevel1"])
  sample_2_psi_values <- numeric_from_comma_string(row["IncLevel2"])
  min_sample_1 <- min(sample_1_psi_values)
  min_sample_2 <- min(sample_2_psi_values)
  return(min(min_sample_1, min_sample_2))
}

max_psi_from_row <- function(row) {
  sample_1_psi_values <- numeric_from_comma_string(row["IncLevel1"])
  sample_2_psi_values <- numeric_from_comma_string(row["IncLevel2"])
  max_sample_1 <- max(sample_1_psi_values)
  max_sample_2 <- max(sample_2_psi_values)
  return(max(max_sample_1, max_sample_2))
}

calculate_mim_max_diff <- function(x) {
  x$min_psi <- apply(x, 1, min_psi_from_row)
  x$max_psi <- apply(x, 1, max_psi_from_row)
  x$psi_diff <- x$max_psi - x$min_psi
  return(x)
}

mark_filter_max_mim_diff <- function(x, p=0.05) {
  x$max_min_diff = ifelse(x$psi_diff > p, "yes", "no")
  return(x)
}

mark_filter_avg <- function(x, p=10) {
  res$avg_passed = ifelse(x$IncAVG_1 >= p & x$IncAVG_2 >= p, "yes", "no")
  return(res)
}

mark_filter_IncLevel <- function(x, bt=0.05, up=0.95) {
  x$IncLevel_passed = ifelse(x$avgPSI_1>bt & x$avgPSI_1 < up & x$avgPSI_2>bt & x$avgPSI_2 < up, "yes", "no")
  return(x)
}

mark_filter_median_SJC <- function(x, p=0) {
  x$Median_passed = ifelse( (x$median_SJC1 > p & x$median_IJC1 > p) | (x$median_SJC2 > p & x$median_IJC2 > p), "yes", "no")
  return(x) 
}

adjust_pvalue <- function(x, m="BH") {
  x$pval.adjust = p.adjust(x$PValue, method=m)
  return(x)
}


jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

overlap_coef <- function(a, b){
  intersection = length(intersect(a, b))
  den = min(length(a), length(b))
  return (intersection/den)
}


# Perform quantile normalization
qn <- function(x){
  rnk <- apply(x, 2, rank, ties.method = "min")
  sorted_x <- apply(x, 2, sort)
  ranked_means <- rowMeans(sorted_x)
  cbind(sorted_x, ranked_means)
  
  x_norm <- matrix(ranked_means[rnk], ncol = ncol(x))
  dimnames(x_norm) <- dimnames(x)
  return(x_norm)
}


mark_overlaps_non_corrected <- function(res, non_corrected) {
  print(table(non_corrected$psi_passed))
  ids_sig_nc = non_corrected$ID[which(non_corrected$psi_passed == "yes")]
  res$non_corrected_overlap = "no"
  o = match(ids_sig_nc, res$ID)
  res$non_corrected_overlap[o] <- "yes"
  return(res)
}


