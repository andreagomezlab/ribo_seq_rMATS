# https://eliocamp.github.io/codigo-r/en/2021/08/docker-renv/
# https://medium.com/@tea_legs/deploying-an-r-environment-in-docker-part-1-1475210ece7b
rm(list = ls())

#set seed for reproductive results
set.seed(123) 

source("src/rMATS_events_visualization.R")
source("src/rMATS_events_IO.R")
source("src/rMATS_processing_events.R")
source("src/rMTAS_events_annotation.R")


##############################


# fix_gene_id('/apps/rmats-turbo/PV_48_DOI_output_novel_var_read_len/')
# fix_gene_id('/apps/rmats-turbo/Cux2_48_DOI_output_novel_var_read_len/')
# fix_gene_id('/apps/rmats-turbo/Rbp4_48_DOI_output_novel_var_read_len/')
# 
# 
# fix_gene_id('/apps/rmats-turbo/PV_48_DOI_output_var_read_len/')
# fix_gene_id('/apps/rmats-turbo/Cux2_48_DOI_output_var_read_len/')
# fix_gene_id('/apps/rmats-turbo/Rbp4_48_DOI_output_var_read_len/')


fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_48hrs_DOI/processed/rMATs/PV_48hrs_DOI_rMATS/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_48hrs_DOI/processed/rMATs/Cux2_48hrs_DOI_rMATS/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_48hrs_DOI/processed/rMATs/Rbp4_48hrs_DOI_rMATS/')

fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_1week_DOI/processed/rMATs/PV_1week_DOI_rMATS/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_1week_DOI/processed/rMATs/Cux2_1week_DOI_rMATS/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1week_DOI/processed/rMATs/Rbp4_1week_DOI_rMATS/')

fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_1month_DOI/processed/rMATs/PV_1month_DOI_rMATS/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_1month_DOI/processed/rMATs/Cux2_1month_DOI_rMATS/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1month_DOI/processed/rMATs/Rbp4_1month_DOI_rMATS/')

fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_48hrs_Psilo/processed/rMATs/PV_48hrs_Psilo_rMATS/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_48hrs_Psilo/processed/rMATs/Cux2_48hrs_Psilo_rMATS/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_48hrs_Psilo/processed/rMATs/Rbp4_48hrs_Psilo_rMATS/')

fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1week_DOI/processed/rMATs/Rbp4_1week_DOI_rMATS/combat_combined/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1month_DOI/processed/rMATs/Rbp4_1month_DOI_rMATS/combat_combined/')
fix_gene_id('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_1month_DOI/processed/rMATs/Cux2_1month_DOI_rMATS/combat_combined/')


dirs = c('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_48hrs_DOI/processed/rMATs/PV_48hrs_DOI_rMATS/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_48hrs_DOI/processed/rMATs/Cux2_48hrs_DOI_rMATS/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_48hrs_DOI/processed/rMATs/Rbp4_48hrs_DOI_rMATS/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_1week_DOI/processed/rMATs/PV_1week_DOI_rMATS/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_1week_DOI/processed/rMATs/Cux2_1week_DOI_rMATS/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1week_DOI/processed/rMATs/Rbp4_1week_DOI_rMATS/non_corrected/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_1month_DOI/processed/rMATs/PV_1month_DOI_rMATS/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_1month_DOI/processed/rMATs/Cux2_1month_DOI_rMATS/non_corrected/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1month_DOI/processed/rMATs/Rbp4_1month_DOI_rMATS/non_corrected/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_48hrs_Psilo/processed/rMATs/PV_48hrs_Psilo_rMATS/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_48hrs_Psilo/processed/rMATs/Cux2_48hrs_Psilo_rMATS/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_48hrs_Psilo/processed/rMATs/Rbp4_48hrs_Psilo_rMATS/',
         #'/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1week_DOI/processed/rMATs/Rbp4_1week_DOI_rMATS/combat/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1week_DOI/processed/rMATs/Rbp4_1week_DOI_rMATS/combat_combined/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1month_DOI/processed/rMATs/Rbp4_1month_DOI_rMATS/combat_combined/',
         '/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_1month_DOI/processed/rMATs/Cux2_1month_DOI_rMATS/combat_combined/'
         #'/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1week_DOI/processed/rMATs/Rbp4_1week_DOI_rMATS/b1_samples/',
         #'/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1week_DOI/processed/rMATs/Rbp4_1week_DOI_rMATS/combat_qn/'
         )
         

study_table = data.frame(Study_id=rep(c("Hsiao"),15),
                         Treatment = c(rep("DOI", 9), rep("Psilo", 3), "DOI", "DOI", "DOI"),
                         CellType = c(rep(c("PV", "Cux2", "Rbp4"), 4), "Rbp4", "Rbp4", "Cux2"),
                         Timepoint = c(rep(c("48hrs"), 3), rep(c("1week"), 3), rep(c("1month"), 3), rep("48hrs", 3), "1week_combat_combined", "1month_combat_combined", "1month_combat_combined"),
                         Study_dir=dirs
                         )
study_table$Study_source = paste(study_table$CellType, study_table$Timepoint, study_table$Treatment, sep = "_")
plot_table(study_table)

et = c("SE", "A5SS", "A3SS", "MXE", "RI")
ct = c("JC", "JCEC")
index_c = 1

### 
#res_table_passed_all=NULL
index_e = 1
index_s = 13

### filter out false positive then calculate FDR

title=paste(study_table$Study_id[index_s], study_table$Study_source[index_s], et[index_e], ct[index_c], sep = "_")
title
res = get_rMATS_results(study_table$Study_dir[index_s], count_type=ct[index_c], event_type=et[index_e], fdr=1) ## Filtering FDR 5%
head(res)
res$source = study_table$Study_source[index_s]
res$study = study_table$Study_id[index_s]

res = calculate_median_SJC(res)
res = calculate_median_IJC(res)

hist(log2(res$median_SJC1 +1))
hist(log2(res$median_SJC2 +1))

hist(log2(res$median_IJC1 +1))
hist(log2(res$median_IJC2 +1))

summary(res$median_SJC1)
summary(res$median_SJC2)

summary(log2(res$median_SJC1 + 1))
summary(log2(res$median_SJC2 + 1))

summary(log2(res$median_IJC1 + 1))
summary(log2(res$median_IJC2 + 1))

res = mark_filter_median_SJC(res, 0)
table(res$Median_passed)

# ggplot(res[res$Median_passed == "yes",], aes(x = IncLevelDifference, y=PValue))+
#   geom_point() +
#   theme_bw(base_size = 17)

res_filt = res[res$Median_passed == "yes",]
res_filt = adjust_pvalue(res_filt, "BH")
head(res_filt)

table(res_filt$FDR < 0.05)
table(res_filt$pval.adjust < 0.05)

ggplot(res_filt, aes(x = FDR, y=pval.adjust))+
  geom_point() +
  theme_bw(base_size = 17)


x <- list(
  FDR = res_filt$ID[res_filt$FDR < 0.05],
  pval.adjust = res_filt$ID[res_filt$pval.adjust < 0.05]
)

plot_venn(x)

res_filt = res_filt[res_filt$pval.adjust < 0.05,]

res = mark_filter_events(res_filt, psi = 0.1)
#res_table_passed_all = rbind(res_table_passed_all,res)
#file_name = paste(title, "PSI_1.0_BH_0.05.tsv", sep = "_")
#file_name
#save_file(obj = res, dir = 'results/rmats/rMATS_events_filtered/', file_name = file_name)
pe = plot_events_filtering(res, title_txt = title)
ph = plot_PSI_hist(res, title)
cowplot::plot_grid(plotlist = list(ph, pe), cols = 1)

res = calculate_mim_max_diff(res)
res = calculate_sd(res)
res = calculate_avg(res)

res = mark_filter_avg(res, p=10)
res = mark_filter_IncLevel(res, 0.05, 0.95)
res = mark_filter_max_mim_diff(res, 0.05)

#table(res$psi_passed)
#table(res$max_min_diff, res$psi_passed)
#table(res$IncLevel_passed, res$psi_passed)
#table(res$Median_passed, res$psi_passed)

# 
#plot_ICounts_SAMPLE(rmats = res[which(res$psi_passed == "yes"),], col1 = 'SJC_SAMPLE_1', col2 = 'SJC_SAMPLE_2', filter = TRUE)
# 
# plot_ICounts_SAMPLE(rmats = res, col1 = 'IJC_SAMPLE_1', col2 = 'IJC_SAMPLE_2', filter = FALSE)
# 
# plot_ICounts_SAMPLE(rmats = res, col1 = 'IncLevel1', col2 = 'IncLevel2', log=FALSE, filter = FALSE)

file_name = paste(title, "_PSI_1.0_BH_0.05.tsv", sep = "_")
file_name
save_file(obj = res, dir = 'results/rmats/Episode_VI/rMATS_events_filtered/', file_name = file_name)

###  annotation
#res = events
# AS_genes = unique(res$geneSymbol[res$psi_passed == "yes"])
# 
# resp = get_enrichPathway(AS_genes, pcutoff=0.05, ncat=10)
# save_file(as.data.frame(resp), 
#           dir = "results/rmats/GO/",
#           file_name = paste("Reactome", study_table$Study_id[index_s], study_table$Study_source[index_s], et[index_e], ct[index_c], "genes_combined.txt", sep = "_"),
#           row = FALSE,
#           sep = "\t")
# 
# p_enrich = plot_enrich_dotplot(as.data.frame(resp), ncat = 10)
# 
# export_svg(obj = p_enrich,
#            dir = 'figures/drafts/rMATS/', 
#            file_name = paste("Reactome", study_table$Study_id[index_s], study_table$Study_source[index_s], et[index_e], ct[index_c], sep = "_"),
#            h = 4, # PV
#            w = 14  # PV
# )
# 
# # number of categories to plot (ncat), Biological Process GO terms (ont)
# ego = get_enrichGO(AS_genes, pcutoff = 0.01, qcutoff = 0.05, ncat = 20, ont = "BP")
# 
# save_file(as.data.frame(ego), 
#           dir = "results/rmats/GO/",
#           file_name = paste("GO", study_table$Study_id[index_s], study_table$Study_source[index_s], et[index_e], ct[index_c], "genes_combined.txt", sep = "_"),
#           row = FALSE,
#           sep = "\t")
# 
# as.data.frame(ego)
# p_ego = plot_enrich_barplot(ego)
# p_ego
# 
# export_svg(obj = p_ego,
#            dir = 'figures/drafts/rMATS/', 
#            file_name = paste("GO", study_table$Study_id[index_s], study_table$Study_source[index_s], et[index_e], ct[index_c], sep = "_"),
#            h = 5.5, # PV
#            w = 11.2  # PV
# )

# AS_UP = unique(res$geneSymbol[res$passed == "yes" & res$status == "TREAT"])
# get_enrichPathway(AS_UP, pcutoff=0.05, ncat=10)
# get_enrichGO(AS_UP, pcutoff=0.01, qcutoff = 0.05, ncat=20)
# 
# AS_DOWN = unique(res$geneSymbol[res$passed == "yes" & res$status == "CTRL"])
# get_enrichPathway(AS_DOWN, pcutoff=0.05, ncat=5)
# get_enrichGO(AS_DOWN, pcutoff=0.01, qcutoff = 0.05, ncat=20)

##

##### filter out non_corrected overlaps
index_nonc = 6
non_corrected = get_rMATS_filtered_results(dir = 'results/rmats/Episode_VI/rMATS_events_filtered/', 
                                           file_name = paste(paste(study_table$Study_id[index_s], study_table$Study_source[index_nonc], et[index_e], ct[index_c], sep = "_"), "non_corrected_PSI_1.0_BH_0.05.tsv", sep = "_"))
head(non_corrected)
table(non_corrected$psi_passed)

res = mark_overlaps_non_corrected(res, non_corrected)
file_name = paste(title, "marked_filt_PSI_1.0_BH_0.05.tsv", sep = "_")
file_name
save_file(obj = res, dir = 'results/rmats/Episode_VI/rMATS_events_filtered/', file_name = file_name)


x <- list(
  non_corr = non_corrected$ID[non_corrected$psi_passed == "yes"],
  combat= res$ID[res$psi_passed == "yes"]
)


pv = plot_venn(x)
pv
export_pdf(obj = pv, dir = "figures/drafts/rMATS/", file_name = "venn_batch_overlapps_Rbp4_1week", h = 4, w = 4)


ph = plot_ICounts_SAMPLE(rmats = res[which(res$psi_passed == "yes"),], col1 = 'SJC_SAMPLE_1', col2 = 'SJC_SAMPLE_2', filter = TRUE)
ph
export_pdf(obj = ph, dir = "figures/drafts/rMATS/", file_name = "heatmap_batch_corrected_all_Rbp4_1week", h = 6, w = 3)


ph = plot_ICounts_SAMPLE(rmats = res[which(res$psi_passed == "yes" & res$non_corrected_overlap == "no"),], col1 = 'SJC_SAMPLE_1', col2 = 'SJC_SAMPLE_2', filter = TRUE)
ph
export_pdf(obj = ph, dir = "figures/drafts/rMATS/", file_name = "heatmap_batch_corrected_exclusive_Rbp4_1week", h = 6, w = 3)

ph = plot_ICounts_SAMPLE(rmats = res[which(res$psi_passed == "yes" & res$non_corrected_overlap == "yes"),], col1 = 'SJC_SAMPLE_1', col2 = 'SJC_SAMPLE_2', filter = TRUE)
ph
export_pdf(obj = ph, dir = "figures/drafts/rMATS/", file_name = "heatmap_batch_corrected_over_non_corrected_Rbp4_1week", h = 6, w = 3)

ph = plot_ICounts_SAMPLE(rmats = res[which(res$ID %in% setdiff(x$non_corr, x$combat)),], col1 = 'SJC_SAMPLE_1', col2 = 'SJC_SAMPLE_2', filter = FALSE)
ph

non_corrected$ID[non_corrected$psi_passed == "yes"]


## filter out overlapping events non-corrected
res_filt = res[res$non_corrected_overlap == "no",]
res_gr = build_coordinates_gr(res_filt, et[index_e])

res_gr_filt = res_gr[res_gr$psi_passed == "yes"]
# total
table(res_gr_filt$psi_passed, res_gr_filt$source)
# status
table(res_gr_filt$source, res_gr_filt$status)
# genes total
aggregate(name~source,as.data.frame(res_gr_filt), function(x) length(unique(x)))
# genes by status
aggregate(name~source+status,as.data.frame(res_gr_filt), function(x) length(unique(x)))

#############


#### 
data_jc = get_AS_associated_genes(res)
res_gr = build_coordinates_gr(res, et[index_e])
###

res_gr_filt = res_gr[res_gr$psi_passed == "yes"]
# total
table(res_gr_filt$psi_passed, res_gr_filt$source)
# status
table(res_gr_filt$source, res_gr_filt$status)
# genes total
aggregate(name~source,as.data.frame(res_gr_filt), function(x) length(unique(x)))
# genes by status
aggregate(name~source+status,as.data.frame(res_gr_filt), function(x) length(unique(x)))

#res_gr_filt

file_name = paste(title,".bed", sep = "")
file_name
save_GR_bed_file(subset(res_gr_filt, source==study_table$Study_source[index_s]), dir = "results/rmats/Episode_VI/rMATS_events_filtered/", file_name = file_name)



###############################################################################################
###### GO union set events

loc_res = read.delim("files/AS_results_location.csv", sep = ",", comment.char = "#")
head(loc_res)

cell_type = "Rbp4" # [PV Cux2ERT Rbp4]
treat = "DOI" # [DOI Psilocybin]
timepoint = "1week" # [1month  1week  48hrs]
b = FALSE
title_txt = paste("Hsiao", cell_type, timepoint, treat, sep="_")
title_txt

AS_genes <- get_combine_gene_events(loc_res, cell_type, treat, timepoint, batch = b) # [PV Cux2ERT Rbp4] [1month  1week  48hrs] 

resp = get_enrichPathway(AS_genes, pcutoff=0.05, ncat=10)
save_file(as.data.frame(resp),
          dir = "results/rmats/Episode_VI/GO/",
          file_name = paste("Reactome", title_txt,  "AS_genes_combined.txt", sep = "_"),
          row = FALSE,
          sep = "\t")

# p_enrich = plot_enrich_dotplot(as.data.frame(resp), ncat = 10)
# p_enrich

# export_svg(obj = p_enrich,
#            dir = 'figures/drafts/rMATS/',
#            file_name = paste("Reactome", title_txt,"AS_genes_combined", sep = "_"),
#            h = 4, # PV
#            w = 14  # PV
# )

# number of categories to plot (ncat), Biological Process GO terms (ont)
ego = get_enrichGO(AS_genes, pcutoff = 0.01, qcutoff = 0.05, ncat = 20, ont = "BP")

save_file(as.data.frame(ego), 
          dir = "results/rmats/Episode_VI/GO/",
          file_name = paste("GO", title_txt, "AS_genes_combined.txt", sep = "_"),
          row = FALSE,
          sep = "\t")

# as.data.frame(ego)
# p_ego = plot_enrich_barplot(ego)
# p_ego
# 
# export_pdf(obj = p_ego,
#            dir = 'figures/drafts/rMATS/', 
#            file_name = paste("GO", title_txt, "AS_genes_combined", sep = "_"),
#            h = 5.5, # PV
#            w = 11.2  # PV)
# )
# 
# export_svg(obj = p_ego,
#            dir = 'figures/drafts/rMATS/', 
#            file_name = paste("GO", title_txt, "AS_genes_combined.txt", sep = "_"),
#            h = 5.5, # PV
#            w = 11.2  # PV
# )


####################################################### Summary Pie


loc_res = read.delim("files/AS_results_location.csv", sep = ",", comment.char = "#")
head(loc_res)

cell_type = "Rbp4" # [PV Cux2ERT Rbp4]
timep = "48hrs" # [DOI Psilocybin]
treat = "Psilocybin" # [1month  1week  48hrs]

bar_data <- get_summary_pie_values(loc_res, cell_type, timepoint = timep, treat = treat) # it handles bathc insede the function
bar_data
bar_data_ = bar_data[bar_data$group == 'PSI.filter',]
bar_data_
bar_data_$per = bar_data_$value/sum(bar_data_$value)
aux= bar_data_

bar_data <- get_summary_pie_values(loc_res, "Cux2ERT", timepoint = "48hrs", treat = "DOI")
bar_data
bar_data_ = bar_data[bar_data$group == 'PSI.filter',]
bar_data_$per = bar_data_$value/sum(bar_data_$value)
aux = rbind(aux, bar_data_)

bar_data <- get_summary_pie_values(loc_res, "Rbp4")
bar_data
bar_data_ = bar_data[bar_data$group == 'PSI.filter',]
bar_data_$per = bar_data_$value/sum(bar_data_$value)
aux = rbind(aux, bar_data_)


euler_data = aggregate(bar_data$value, by=list(group=bar_data$group), FUN=sum)
colnames(euler_data)[2] <- 'value'
euler_data = euler_data[c(4,1,3,2),]

p1 = plot_euler(euler_data)
export_svg(p1, dir = "figures/drafts/rMATS/", file_name = "PV_AS_events_euler_plot", h = 5, w = 5)

p1 = plot_euler(euler_data)
export_svg(p1, dir = "figures/drafts/rMATS/", file_name = "Cux2_AS_events_euler_plot", h = 5, w = 5)

p1 = plot_euler(euler_data)
export_svg(p1, dir = "figures/drafts/rMATS/", file_name = "Rbp4_AS_events_euler_plot", h = 5, w = 5)


p1 <- ggplot(aux, aes(x = celltype, fill = event_type, y = per)) +
  geom_bar(position="fill", stat="identity", color="black",) +
  geom_text(aes(label = paste0(signif(round(per, 4) * 100, 3), "%" )), position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values=c("#fc8d59", "#fee08b", "#e6f598","#abd9e9", "#3288bd")) +
  labs(y="Splicing events", x="") +
  theme_bw(base_size = 17)
p1

export_svg(p1, dir = "figures/drafts/rMATS/", file_name = "ASplicing_events_per_type", h = 5, w = 5)
  
#####

## baseline considering all the sources in the GR paramater
# gr_overlapping_counts = count_overlapping_events(res_gr_filt)
# gr_overlapping_counts_filt = subset(gr_overlapping_counts, counts > 1)
# 
# table(gr_overlapping_counts_filt$name, gr_overlapping_counts_filt$counts)
# 
# file_name = paste0(et[index_e], "_baseline_counts_at_least_2.bed")
# save_GR_bed_file(gr_overlapping_counts_filt, dir = "results/rmats/rMATS_events_filtered/", file_name = file_name, extra = FALSE)

## build a workflow for rmats2sashimiplot plot visualization
# 
# c("/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam1_1_len.txt",
#   "/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam1_2_len.txt",
#   "/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam1_3_len.txt",
#   "/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam2_1_len.txt",
#   "/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam2_2_len.txt",
#   "/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam2_3_len.txt")
# 
# 
# aux = read.delim("/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam1_1_len.txt", header = F)
# aux$sample = 'C1'
# read_len = aux
# 
# aux = read.delim("/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam1_2_len.txt", header = F)
# aux$sample = 'C2'
# read_len = rbind(read_len, aux)
# 
# aux = read.delim("/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam1_3_len.txt", header = F)
# aux$sample = 'C3'
# read_len = rbind(read_len, aux)
# 
# 
# aux = read.delim("/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam2_1_len.txt", header = F)
# aux$sample = 'DOI1'
# read_len = aux
# 
# aux = read.delim("/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam2_2_len.txt", header = F)
# aux$sample = 'DOI2'
# read_len = rbind(read_len, aux)
# 
# aux = read.delim("/apps/rmats-turbo/PV_48_DOI_tmp_output_novel_var_read_len/bam2_3_len.txt", header = F)
# aux$sample = 'DOI3'
# read_len = rbind(read_len, aux)
# 
# table(read_len$sample)
# 
# options(scipen = 999)
# 
# ggplot(data = read_len, aes(x = V1, fill = sample)) +
#   geom_histogram(alpha = 0.30, position = 'identity')+
#   #geom_density(alpha = 0.30) +
#   xlab("Read length") +
#   ylab("Counts") +
#   theme_gray(base_size = 16)
# 
# summary(read_len)




######################################### Compare AS genes


loc_res = read.delim("files/AS_results_location.csv", sep = ",", comment.char = "#")
head(loc_res)

data_jc = data.frame(id="PV_48hrs_DOI", Genes=get_combine_gene_events(loc_res, celltype = "PV", treat = "DOI"))
data_jc = rbind(data_jc, data.frame(id="PV_48hrs_Psilo", Genes=get_combine_gene_events(loc_res, celltype = "PV", treat = "Psilocybin")))
data_jc = rbind(data_jc, data.frame(id="Cux2_48hrs_DOI", Genes=get_combine_gene_events(loc_res, celltype = "Cux2ERT", treat = "DOI")))
data_jc = rbind(data_jc, data.frame(id="Cux2_48hrs_Psilo", Genes=get_combine_gene_events(loc_res, celltype = "Cux2ERT", treat = "Psilocybin")))
data_jc = rbind(data_jc, data.frame(id="Rbp4_48hrs_DOI", Genes=get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "DOI")))
data_jc = rbind(data_jc, data.frame(id="Rbp4_48hrs_Psilo", Genes=get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "Psilocybin")))
table(data_jc$id)

ids = unique(data_jc$id)
table_jc = NULL

for (i in 1:length(ids)) {
  for (j  in 1:length(ids)){
    x = ids[i]
    y = ids[j]
    a = data_jc$Genes[data_jc$id == x]
    b = data_jc$Genes[data_jc$id == y]
    j_index = jaccard(a, b)  
    row = data.frame(x=x, y=y, value=j_index)
    table_jc = rbind(table_jc, row)
  }
}

table_jc$x <- factor(table_jc$x, levels=rev(levels(factor(table_jc$x))))

library(hrbrthemes)
library(ggplot2)

ggplot(table_jc, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = signif(value, 1)), color = "black", size = 4) +
  coord_fixed() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  theme_ipsum() +
  ylab("") +
  labs(title="Jaccard index - DEG") +
  theme(
    axis.text.x = element_text( angle = 45,vjust = 0.5, hjust=0.5)
  )

###


ids = unique(data_jc$id)
table_oc = NULL

for (i in 1:length(ids)) {
  for (j  in 1:length(ids)){
    x = ids[i]
    y = ids[j]
    a = data_jc$Genes[data_jc$id == x]
    b = data_jc$Genes[data_jc$id == y]
    j_index = overlap_coef(a, b)  
    row = data.frame(x=x, y=y, value=j_index)
    table_oc = rbind(table_oc, row)
  }
}

table_oc$x <- factor(table_oc$x, levels=rev(levels(factor(table_oc$x))))


ggplot(table_oc, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "purple") +
  geom_text(aes(label = signif(value, 1)), color = "black", size = 3) +
  coord_fixed() +
  #theme_ipsum() +
  ylab("") +
  labs(title="Paired Overlap Coeficient - AS genes") +
  theme(axis.text.x = element_text(angle = 45,
                               vjust=0.9,
                               hjust=1,
                               lineheight=0.75)
  )



clusters = readxl::read_xlsx("data/1-s2.0-S2211124721013000-mmc8.xlsx")
clusters$id = 1:nrow(clusters)
molted=melt(clusters,id.vars=c("id"))
molted

molted = molted[,-1]
colnames(molted) = c("id", "Genes")

data_jc = rbind(data_jc, molted[, c('id', 'Genes')])
table(data_jc$id)

ids = unique(data_jc$id)
table_jc = NULL

for (i in 1:length(ids)) {
  for (j  in 1:length(ids)){
    x = ids[i]
    y = ids[j]
    a = data_jc$Genes[data_jc$id == x]
    b = data_jc$Genes[data_jc$id == y]
    j_index = jaccard(a, b)  
    row = data.frame(x=x, y=y, value=j_index)
    table_jc = rbind(table_jc, row)
  }
}

table_jc$x <- factor(table_jc$x, levels=rev(levels(factor(table_jc$x))))

library(hrbrthemes)
library(ggplot2)

ggplot(table_jc, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = signif(value, 1)), color = "black", size = 4) +
  coord_fixed() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  theme_ipsum() +
  ylab("") +
  labs(title="Jaccard index - AS genes") +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust=0.9,
                                   hjust=1,
                                   lineheight=0.75)
  )


###


ids = unique(data_jc$id)
table_oc = NULL

for (i in 1:length(ids)) {
  for (j  in 1:length(ids)){
    x = ids[i]
    y = ids[j]
    a = data_jc$Genes[data_jc$id == x]
    b = data_jc$Genes[data_jc$id == y]
    j_index = overlap_coef(a, b)  
    row = data.frame(x=x, y=y, value=j_index)
    table_oc = rbind(table_oc, row)
  }
}

table_oc$x <- factor(table_oc$x, levels=rev(levels(factor(table_oc$x))))


ggplot(table_oc, aes(x = x, y = y, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "purple") +
  geom_text(aes(label = signif(value, 1)), color = "black", size = 3) +
  coord_fixed() +
  theme_ipsum() +
  ylab("") +
  labs(title="Paired Overlap Coeficient - AS and Modules") +
  theme(axis.text.x = element_text(angle = 45,
                               vjust=0.9,
                               hjust=1,
                               lineheight=0.75)
  )


######


loc_res = read.delim("files/AS_results_location.csv", sep = ",", comment.char = "#")
head(loc_res)
table(loc_res$CellType)

x <- list(
  Rbp4_48_DOI = get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "DOI"),
  Rbp4_48_Psilo = get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "Psilocybin")
)


plot_venn(x, "DOI vs Psilo gene overlap")
plot_venn_v2(x, "DOI vs Psilo gene overlap")



##### Cluster Profile AS genes

title_txt = "AS genes Psilo DOI"
library(clusterProfiler)
library(org.Mm.eg.db)


x <- list(
  PV_48_DOI_AS = get_combine_gene_events(loc_res, celltype = "PV", treat = "DOI"),
  PV_48_Psilo_AS = get_combine_gene_events(loc_res, celltype = "PV", treat = "Psilocybin"),
  Cux2_48_DOI_AS = get_combine_gene_events(loc_res, celltype = "Cux2ERT", treat = "DOI"),
  Cux2_48_Psilo_AS = get_combine_gene_events(loc_res, celltype = "Cux2ERT", treat = "Psilocybin"),
  Rbp4_48_DOI_AS = get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "DOI"),
  Rbp4_48_Psilo_AS = get_combine_gene_events(loc_res, celltype = "Rbp4", treat = "Psilocybin")
)


gs <- lapply(x, function(symbols) {
  IDs = mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  return(IDs)
})

res <- compareCluster(gs, 
                      fun="enrichGO", 
                      pvalueCutoff=0.05, 
                      ont          = "BP",
                      pAdjustMethod="BH", 
                      OrgDb=org.Mm.eg.db) %>% 
  setReadable(., OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

p1 = clusterProfiler::dotplot(res, showCategory=20, title=title_txt) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient(low = "white", 
                      high = "#317ebb",
                      trans = 'reverse') +
  scale_y_discrete(labels=function(x) str_wrap(x, width=60)) +
  labs(
    x = "",
    y = "",
    fill = "P-value adjusted",
    size = "Enrichment"
  )
p1

export_pdf(p1, dir = "figures/drafts/", file_name = paste0("ClusterProfiler_", title_txt), h = 9, w = 8)

save_file(p1$data, 
          dir = "results/rmats/Episode V/",
          file_name = paste("ClusterProfile", title_txt, sep = "_"),
          row = FALSE,
          sep = "\t")


####################### enrichPathway

res <- compareCluster(gs, 
                      fun="enrichPathway", 
                      pvalueCutoff=0.05, 
                      pAdjustMethod="BH", 
                      readable = TRUE,
                      organism = "mouse")

p1 = clusterProfiler::dotplot(res, showCategory=20, title=title_txt) +
  theme(axis.text.x = element_text(angle = 90))
p1

export_svg(p1, dir = "figures/drafts/", file_name = paste0("ClusterProfiler_Reactome_", title_txt, ".tsv"), h = 5, w = 6)

save_file(p1$data, 
          dir = "results/rmats/Episode V/",
          file_name = paste0("ClusterProfile_Reactome_", title_txt, ".tsv"),
          row = FALSE,
          sep = "\t")


### Individual Cells treatment

cell_type = "Rbp4" # ["PV" "Cux2ERT" "Rbp4"]
timepoint = "48hrs" # [48hrs 1week 1month]
treat = "Psilocybin" # ["DOI" "Psilocybin"]
title_txt = paste(cell_type, timepoint, treat, sep = "_")
title_txt

AS_genes <- get_combine_gene_events(loc_res, cell_type, treat)

resp = get_enrichPathway(AS_genes, pcutoff=0.05, ncat=20)
save_file(as.data.frame(resp), 
          dir = "results/rmats/Episode V/",
          file_name = paste("Reactome", title_txt,  "events_genes_combined.txt", sep = "_"),
          row = FALSE,
          sep = "\t")

# number of categories to plot (ncat), Biological Process GO terms (ont)
ego = get_enrichGO(AS_genes, pcutoff = 0.01, qcutoff = 0.05, ncat = 20, ont = "BP")

save_file(as.data.frame(ego), 
          dir = "results/rmats/Episode V/",
          file_name = paste("GO", title_txt, "events_genes_combined.txt", sep = "_"),
          row = FALSE,
          sep = "\t")

#############################################################################################################################################
### batch correction using Combat


library(DESeq2)
library(Hmisc)

cell_type = "Cux2ERT" # ["PV" "Cux2ERT" "Rbp4"]
timepoint = "1month" # [48hrs 1week 1month]
treat = "DOI" # ["DOI" "Psilocybin"]
title_txt = paste(cell_type, timepoint, treat, sep = "_")
title_txt

meta_table = get_metatable(cell_type = cell_type, timepoint = timepoint, treat = treat)
plot_table(meta_table)


dir_txt = "/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1week_DOI/processed/rMATs/Rbp4_1week_DOI_rMATS/"
dir_txt = "/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Rbp4_1month_DOI/processed/rMATs/Rbp4_1month_DOI_rMATS/"
dir_txt = "/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/Cux2_1month_DOI/processed/rMATs/Cux2_1month_DOI_rMATS/"

e_type = "RI" # [SE A5SS A3SS MXE RI]

counts_ctrl = get_rMATS_RAW_input(dir_txt, count_type="JC", event_type=e_type, group = "SAMPLE_1") # [IJC SJC]
head(counts_ctrl)
dim(counts_ctrl)
counts_treat = get_rMATS_RAW_input(dir_txt, count_type="JC", event_type=e_type, group = "SAMPLE_2") # [IJC SJC])
head(counts_treat)
dim(counts_treat)

### Merge for PCA
counts = cbind(counts_ctrl, counts_treat)
head(counts, 30)
tail(counts)
dim(counts)

hist.data.frame(log2(counts))

meta_table$File.ID
colnames(counts)
all(meta_table$File.ID %in% colnames(counts))
rownames(meta_table) <- meta_table$File.ID

# Creating deseq2 object
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = meta_table, 
                              design = ~Treatment)

vst <- vst(dds, blind=T)
pca <- plotPCA(vst, "Batch")
pca + theme_classic2(base_size = 17) +
  geom_text_repel(data = pca$data, aes(x=PC1, y=PC2, label=name), box.padding = 0.8) +
  labs(title = title_txt)

################################################################################ Ctrl and Treat separated
# nonzeros_ctrl = which(!rowSums(counts_ctrl == 0) >= 4)
# counts_ctrl_filt = counts_ctrl[nonzeros_ctrl ,]
# dim(counts_ctrl)  
# dim(counts_ctrl_filt)  
# meta_table_filt = meta_table[meta_table$Treatment == "Ctrl",]
# cbat_ctrl <- ComBat(as.matrix(counts_ctrl_filt), batch=meta_table_filt$Batch)
# sapply(as.data.frame(cbat_ctrl), function(x) quantile(x, probs = seq(0, 1, 1/4)))
# 
# cbat_ctrl_qn = cbat_ctrl
# ## https://www.r-bloggers.com/2024/03/mastering-quantile-normalization-in-r-a-step-by-step-guide/#:~:text=At%20its%20core%2C%20quantile%20normalization,values%2C%20making%20meaningful%20comparisons%20possible.
# #https://www.statology.org/quantile-normalization-in-r/
# #perform quantile normalization
# #cbat_ctrl_qn = qn(cbat_ctrl)
# sapply(as.data.frame(cbat_ctrl_qn), function(x) quantile(x, probs = seq(0, 1, 1/4)))
# 
# head(counts_ctrl_filt)
# head(cbat_ctrl_qn)
# table(cbat_ctrl_qn < 0)
# cbat_ctrl_qn[cbat_ctrl_qn < 0] = 0
# cbat_ctrl_final = round(cbat_ctrl_qn, 0)
# head(cbat_ctrl_final)
# dim(cbat_ctrl_final)
# 
# #Rbp4 1 week
# plot(cbat_ctrl_final[,'AG2395_merged'], counts_ctrl_filt[,'AG2395_merged'])
# #Rbp4 1 month
# plot(cbat_ctrl_final[,'AG2068_merged'], counts_ctrl_filt[,'AG2068_merged'])
# #Cux2 1 month
# plot(cbat_ctrl_final[,'AG4469_merged'], counts_ctrl_filt[,'AG4469_merged'])
# 
# 
# nonzeros_treat = which(!rowSums(counts_treat == 0) >= 4)
# counts_treat_filt = counts_treat[nonzeros_treat ,]
# dim(counts_treat)  
# dim(counts_treat_filt)  
# meta_table_filt = meta_table[meta_table$Treatment == "DOI",]
# 
# cbat_treat <- ComBat(as.matrix(counts_treat_filt), batch=meta_table_filt$Batch)
# sapply(as.data.frame(cbat_treat), function(x) quantile(x, probs = seq(0, 1, 1/4)))
# 
# cbat_treat_qn = cbat_treat # no QN
# 
# #cbat_treat_qn = qn(cbat_treat)
# sapply(as.data.frame(cbat_treat_qn), function(x) quantile(x, probs = seq(0, 1, 1/4)))
# 
# head(counts_treat_filt)
# head(cbat_treat_qn)
# table(cbat_treat_qn < 0)
# cbat_treat_qn[cbat_treat_qn < 0] = 0
# cbat_treat_final = round(cbat_treat_qn, 0)
# head(cbat_treat_final)
# 
# dim(cbat_ctrl_final)
# dim(cbat_treat_final)
# #Rbp4 1 week
# plot(cbat_treat_final[,'AG2388_merged'], counts_treat_filt[,'AG2388_merged'])
# #Rbp4 1 month
# plot(cbat_treat_final[,'AG2064_merged'], counts_treat_filt[,'AG2064_merged'])
# #Cux2 1 month
# plot(cbat_ctrl_final[,'AG4469_merged'], counts_ctrl_filt[,'AG4469_merged'])
# 
# ## add zeros back
# head(counts_ctrl[nonzeros_ctrl,])
# counts_ctrl[nonzeros_ctrl,] = cbat_ctrl_final
# head(counts_ctrl[nonzeros_ctrl,])
# head(counts_treat[nonzeros_treat,])
# counts_treat[nonzeros_treat,] = cbat_treat_final
# head(counts_treat[nonzeros_treat,])
# 
# cbat = cbind(counts_ctrl, counts_treat)
# counts_cbat = cbat
####################################################################################
## Corrected combined groups

nonzeros = which(!rowSums(counts == 0) >= 8)
counts_filt = counts[nonzeros ,]
dim(counts)  
dim(counts_filt)  
#mod <- model.matrix(~ Treatment, meta_table)
cbat <- ComBat(as.matrix(counts_filt), batch=meta_table$Batch)
#cbat = qn(cbat)

## add zeros back
head(counts_filt)
head(cbat)

counts_cbat = counts
counts_cbat[nonzeros,] = cbat
table(counts_cbat < 0)
counts_cbat[counts_cbat < 0] = 0
counts_cbat = round(counts_cbat, 0)


dim(counts)
dim(counts_cbat)
head(counts, 20)
head(counts_cbat, 20)

tail(counts, 20)
tail(counts_cbat, 20)


pca_combat = prcomp(t(counts_cbat))

p1 = autoplot(pca_combat,
              data = meta_table,
              colour="Treatment",
              shape="Sex",
              size=5) +
  theme_classic2(base_size = 17) +
  geom_text_repel(aes(x=PC1, y=PC2, label=File.ID), box.padding = 0.8) +
  scale_color_manual(values = figure_colors$Treatment)
p1
p1 = autoplot(pca_combat,
              data = meta_table,
              colour="Batch",
              shape="Treatment",
              size=5) +
  theme_classic2(base_size = 17) +
  geom_text_repel(aes(x=PC1, y=PC2, label=File.ID), box.padding = 0.8) +
  scale_color_manual(values = figure_colors$Batch)
p1

hist.data.frame(log2(counts_cbat))

counts_cbat = as.data.frame(counts_cbat)

summary(counts)
summary(counts_cbat)

## add zeros back


write_rMATS_RAW_input(counts_cbat, dir_txt, count_type="JC", event_type=e_type)


################################################################################
## check candidates sashimi


loc_res = read.delim("files/AS_results_location.csv", sep = ",", comment.char = "#")
head(loc_res)
table(loc_res$CellType)

cell_type = "Rbp4" # ["PV" "Cux2ERT" "Rbp4"]
timepoint = "48hrs" # [48hrs 1week 1month]
treat = "Psilocybin" # ["DOI" "Psilocybin"]
title_txt = paste(cell_type, timepoint, treat, sep = "_")
title_txt

meta_table = get_metatable(cell_type, timepoint)
meta_table

loc_res_ = loc_res[loc_res$CellType == cell_type & loc_res$Timepoint == timepoint & loc_res$Treatment == treat,]
go_terms = read.delim(loc_res_$GO[1], sep = "\t")
go_genes = unique(unlist(strsplit(go_terms$geneID, "/")))

events = read.delim(loc_res_$Filtered[1])
events_ = events[events$GeneID %in% go_genes & events$psi_passed == "yes",]
dim(events_)

events_ = events_[order(events_$IncLevelDifference, decreasing = F),]
topN = 25

title_txt = gsub('ERT', "", title_txt)
title_txt = gsub('cybin', "", title_txt)
root_path = paste0('/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/',title_txt, '/processed/rMATs/', title_txt, '_rMATS/')
root_path

#G_PATH = '/apps/rmats-turbo/grouping_4.txt'
G_PATH = '/apps/rmats-turbo/grouping.txt'
#G_PATH = '/apps/rmats-turbo/grouping_5_4.txt'

for (i in 1:topN) {
  
  GENE = events_$GeneID[i]
  ID = events_$ID[i]
  
  print(GENE)
  
  b1 = paste0(' --b1 ', root_path, 'ctrl.txt ')
  b2 = paste0(' --b2 ', root_path, 'treat.txt ')
  cmd = paste0('cat ', root_path, 'SE.MATS.JC_processed.txt |  head -1 > ', root_path, GENE, '_', ID, '_SE.MATS.JC.txt')
  system(cmd, wait = TRUE)
  cmd = paste0('grep ', GENE, ' ', root_path, 'SE.MATS.JC_processed.txt | grep ', ID, ' >> ', root_path, GENE, '_', ID, '_SE.MATS.JC.txt')
  system(cmd, wait = TRUE)
  if (!dir.exists(paste0('Sashimi_plot ', root_path, 'Sashimi_plot_', GENE, '_', ID, '_group'))) {
    efile = paste0(root_path, GENE, "_", ID, '_SE.MATS.JC.txt')
  
    run_sashimi_cmd = 'conda run -n py2_env rmats2sashimiplot'
    run_sashimi_cmd = paste0(run_sashimi_cmd, b1, b2, '--event-type SE -e ', efile, ' --l1 Ctrl --l2 DOI --remove-event-chr-prefix --exon_s 1 --intron_s 5 --min-counts 1 -o ', root_path)
    run_sashimi_cmd
    system(run_sashimi_cmd, wait = TRUE)
    cmd = paste0('mv ', root_path, 'Sashimi_plot ', root_path, 'Sashimi_plot_', GENE, '_', ID, '_rep')
    system(cmd, wait = TRUE)
    cmd = paste0('rm -r ', root_path, 'Sashimi_index*')
    system(cmd, wait = TRUE)
    
    run_sashimi_cmd = paste0(run_sashimi_cmd, b1, b2, '--event-type SE -e ', efile, ' --l1 Ctrl --l2 DOI --remove-event-chr-prefix --exon_s 1 --intron_s 5 --min-counts 1 -o ', root_path, ' --group-info ', G_PATH)
    run_sashimi_cmd
    system(run_sashimi_cmd, wait = TRUE)
    cmd = paste0('mv ', root_path, 'Sashimi_plot ', root_path, 'Sashimi_plot_', GENE, '_', ID, '_group')
    system(cmd, wait = TRUE)
    cmd = paste0('rm -r ', root_path, 'Sashimi_index*')
    system(cmd, wait = TRUE)  
  } else {
    print("File exists")
  }
}

