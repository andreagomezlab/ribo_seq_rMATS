# https://eliocamp.github.io/codigo-r/en/2021/08/docker-renv/
# https://medium.com/@tea_legs/deploying-an-r-environment-in-docker-part-1-1475210ece7b
rm(list = ls())

source("src/rMATS_events_visualization.R")
source("src/rMATS_events_IO.R")
source("src/rMATS_processing_events.R")

study_table = data.frame(Study_id=c("de_la_fuente_Revenga", "Nardou", "Nardou", "Nardou", "Nardou", "Hsiao", "Hsiao", "Hsiao"),
                   Study_source = c("DOI", "LSD", "Cocaine", "Ketamine", "MDMA", "Cux2_DOI", "PV_DOI", "Rbp4_DOI"),
                   Study_dir=c("/apps/rmats-turbo/de_la_fuente_Revenga_48hrs_output/", 
                           "/apps/rmats-turbo/Nardou_48hrs_LSD_output/",
                           "/apps/rmats-turbo/Nardou_48hrs_Cocaine_output/",
                           "/apps/rmats-turbo/Nardou_48hrs_katemine_output/",
                           "/apps/rmats-turbo/Nardou_48hrs_MDMA_output/",
                           "results/rmats/cux2_48_DOI_output_novel/",
                           "results/rmats/PV_48_DOI_output_novel/",
                           "results/rmats/Rbp4_48_DOI_output_novel/"
                           ))

plot_table(study_table)

et = c("SE", "A5SS", "A3SS", "MXE", "RI")
ct = c("JC", "JCE")

### 
res_table_passed_all=NULL
index_c = 1
index_e = 1
index_s = 8

title=paste(study_table$Study_id[index_s], study_table$Study_source[index_s], et[index_e], ct[index_c], sep = "_")
title
res = get_rMATS_results(study_table$Study_dir[index_s], count_type=ct[index_c], event_type=et[index_e]) ## Filtering FDR 5%
head(res)
res$source = study_table$Study_source[index_s]
res$study = study_table$Study_id[index_s]
res = mark_filter_events(res, qt=0.8, psi = 0.1)
res_table_passed_all = rbind(res_table_passed_all,res)
file_name = paste(title, "PSI_10_FDR_005.tsv", sep = "_")
file_name
#save_file(obj = res, dir = study_table$Study_dir[index_s], file_name = file_name)
plot_events_filtering(res, title_txt = title)

#### 
data_jc = get_AS_associated_genes(res_table_passed_all)
res_gr = build_coordinates_gr(res_table_passed_all, et[index_e])
###

table(res_table_passed_all$passed, res_table_passed_all$source)
table(res_gr$passed, res_gr$source)

res_gr_filt = res_gr[res_gr$passed == "yes"]
table(res_gr_filt$passed, res_gr_filt$source)
table(res_gr_filt$source, res_gr_filt$status)
res_gr_filt

# applying aggregate function 
aggregate(name~source+status,as.data.frame(res_gr_filt), function(x) length(unique(x)))


file_name = paste(title,".bed", sep = "")
file_name
save_GR_bed_file(subset(res_gr_filt, source==study_table$Study_source[index_s]), dir = "files/", file_name = file_name)

## baseline considering all the sources in the GR paramater
gr_overlapping_counts = count_overlapping_events(res_gr_filt)
gr_overlapping_counts_filt = subset(gr_overlapping_counts, counts > 1)

table(gr_overlapping_counts_filt$name, gr_overlapping_counts_filt$counts)

file_name = paste0(et[index_e], "_baseline_counts_at_least_2.bed")
save_GR_bed_file(gr_overlapping_counts_filt, dir = "files/", file_name = file_name, extra = FALSE)


## fix chromossome name for plotting and table results
## build a workflow for rmats2sashimiplot plot visualization
## refine pipeline for BAMs input and aligment with STAR

