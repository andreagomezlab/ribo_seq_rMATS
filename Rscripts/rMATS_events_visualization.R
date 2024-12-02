library(hrbrthemes)
library(ggVennDiagram)
library("VennDiagram")
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)

plot_table <- function(table) {
  grid.newpage()
  return(gridExtra::grid.table(table))
}

plot_events_filtering <- function(res, title_txt) {
  p1 = ggplot(res, aes(avgPSI_1, avgPSI_2, color=psi_passed))  +
    geom_point(aes(shape = status ),
               stroke = 2
    ) +
    theme_ipsum(base_size = 14) +
    labs(colour = paste0("InLevellDiff > |0.1|"), shape="IncLevelDifference", title = paste(title_txt, table(res$psi_passed)[2], "of", nrow(res), "events"))
  return(p1)
}


plot_venn <- function(x) {
  p1 = ggVennDiagram(x) +
    scale_color_brewer(palette = "Paired")
  return(p1)
}

plot_PSI_hist <- function(res, title) {
  
  p1 <- ggplot(res, aes(IncLevelDifference, fill=psi_passed)) + 
    geom_histogram(color="black") + 
    xlab("IncLevelDifference") + ylab("count") +
    ggtitle(paste0(title)) +
    theme_classic2(base_size = 13)
  return(p1)
}

plot_ICounts_SAMPLE <- function(rmats, col1='IJC_SAMPLE_1', col2='IJC_SAMPLE_2', filter=TRUE, log=TRUE) {
  
  rmats_filt = rmats
  
  if (filter)
    rmats_filt = rmats[rmats$psi_passed == 'yes',]
  
  rownames(rmats_filt) <- rmats_filt$ID
  rmats_filt = rmats_filt[order(rmats_filt$IncLevelDifference),]
  
  splitdat = do.call("rbind", strsplit( rmats_filt[, colnames(rmats_filt) == col1 ], ","))
  splitdat = data.frame(apply(splitdat, 2, as.numeric))
  names(splitdat) = paste("ctrl", 1:5, sep = "")
  head(splitdat)
  matrix_IJC = splitdat
  
  splitdat = do.call("rbind", strsplit(rmats_filt[, colnames(rmats_filt) == col2 ], ","))
  splitdat = data.frame(apply(splitdat, 2, as.numeric))
  names(splitdat) = paste("treat", 1:5, sep = "")
  head(splitdat)
  matrix_IJC = cbind(matrix_IJC, splitdat)
  
  meta_table = rmats_filt[, c('ID', 'PValue', 'FDR', 'IncLevelDifference','avgPSI_1', 'avgPSI_2', 'psi_passed', 'min_psi', 'max_psi', 'psi_diff', 'IncSD_1', 'IncSD_2', 'IncAVG_1', 'IncAVG_2', "IncLevel_passed", "max_min_diff", "avg_passed", "Median_passed")]
  
  column_ha = HeatmapAnnotation(Batch = c(rep("B1", 5), rep("B2", 5)))
  
  if (log) 
    matrix_IJC = log2(matrix_IJC + 1)
  
  ## add column annotation
  ha = rowAnnotation(df = meta_table[, c(2:ncol(meta_table))])
  
  p1 <- Heatmap(as.matrix(matrix_IJC), 
                cluster_rows = F,
                name = "mat", 
                cluster_columns = T,
                col = colorRampPalette(c("white", "blue"))(10),
                #right_annotation = ha,
                top_annotation = column_ha)
  draw(p1)
  return(p1)
}

make_ratio_numeric = function(ratio) {
  ratio = as.numeric(unlist(str_split(ratio, "/"))[1])/as.numeric(unlist(str_split(ratio, "/"))[2])
  ratio
}

plot_enrich_dotplot <- function(enrich.df, ncat=10) {
  enrich.df = enrich.df[order(enrich.df$Count),]
  enrich.df <- cbind(enrich.df,order=1:nrow(enrich.df))
  enrich.df$GeneRatio_num =  sapply(enrich.df$GeneRatio, make_ratio_numeric)
  print(head(enrich.df))
  p1 <- ggplot(head(enrich.df, ncat), aes(x = GeneRatio_num, y = order, fill = p.adjust, size = Count)) +
    geom_point(shape = 21, stroke = 0.5) + # change the thickness of the boarder with stroke
    scale_fill_gradient(low = "#317ebb", 
                        high = "#de6763",
                        trans = 'reverse') +
    ylab("") +
    xlab("") + 
    scale_y_continuous(breaks=1:nrow(enrich.df),labels=enrich.df$Description) +
    guides(color = guide_colorbar(reverse = TRUE)) +
    theme_bw(base_size = 15) +
    theme(axis.text.y=element_text(colour="black"),
          axis.text.x=element_text(colour="black")) +
    labs(
      x = "Gene Ratio",
      y = "",
      fill = "P-value adjusted",
      size = "Enrichment"
    )

  
  return(p1)
}

plot_enrich_barplot <- function(enrich.df, ncat=20) {
  
  top = head(enrich.df, ncat)
  top$Description <- factor(top$Description, levels = top$Description)
  
  p1 <- ggplot(top, aes(x = Count, y = reorder(Description), fill = p.adjust)) +
    geom_bar(stat="identity", width=0.9) +
    scale_fill_gradient(low = "#317ebb", 
                        high = "#de6763",
                        trans = 'reverse') +
    guides(color = guide_colorbar(reverse = TRUE)) +
    theme_bw(base_size = 15) +
    theme(axis.text.y=element_text(colour="black"),
          axis.text.x=element_text(colour="black")) +
    labs(
      x = "Enrichment",
      y = "",
      fill = "P-value adjusted",
      #size = "Enrichment"
    )
  
  return(p1)
}


plot_pie <- function(pie_data) {
  # Compute the position of labels
  pie_data <- pie_data %>% 
    arrange(desc(group)) %>%
    mutate(prop = values / sum(pie_data$values) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  # Basic piechart
  p1 <- ggplot(pie_data, aes(x="", y=prop, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="none") +
    
    geom_text(aes(y = ypos, label = paste0(group, " ", values)), color = "black", size=6) +
    scale_fill_manual(values=c("#fc8d59", "#fee08b", "#e6f598","#3288bd"))
  return(p1)
}

plot_euler <- function(pie_data) {
  fit = euler(c('total' = pie_data$value[4],
                'total&events_passed_lc' = pie_data$value[2],
                'total&events_passed_lc&Pval.adjust' = pie_data$value[3],
                'total&events_passed_lc&Pval.adjust&PSI.filter' = pie_data$value[4]))
  plot(fit,
       fills=list(fill = c("#f7f7f7", "#fde0ef", "#f1b6da","#de77ae")),
       labels=list(col="black", font=4), quantities=T)
}

