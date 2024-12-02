rm(list = ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)

Hsiao_cux2_rmats = read.delim("results/rmats/Cux2_48_DOI_output_novel_var_read_len/summary.txt")
Hsiao_cux2_rmats$ID = "Hisao_Cux2_DOI"
Hsiao_pv_rmats = read.delim("results/rmats/PV_48_DOI_output_novel_var_read_len//summary.txt")
Hsiao_pv_rmats$ID = "Hisao_PV_DOI"
Hsiao_rbp4_rmats = read.delim("results/rmats/Rbp4_48_DOI_output_novel_var_read_len//summary.txt")
Hsiao_rbp4_rmats$ID = "Hisao_Rbp4_DOI"


Hsiao_pv_rmats_novel = read.delim("results/rmats/PV_48_DOI_output_novel_var_read_len/summary.txt")
Hsiao_pv_rmats_novel$ID = "Hisao_PV_DOI_novel"
Hsiao_pv_rmats = read.delim("results/rmats/PV_48_DOI_output_var_read_len//summary.txt")
Hsiao_pv_rmats$ID = "Hisao_PV_DOI"
sum_rmats = rbind(Hsiao_pv_rmats_novel, Hsiao_pv_rmats)


Hsiao_cux2_rmats_novel = read.delim("results/rmats/Cux2_48_DOI_output_novel_var_read_len/summary.txt")
Hsiao_cux2_rmats_novel$ID = "Hisao_Cux2_DOI_novel"
Hsiao_cux2_rmats = read.delim("results/rmats/Cux2_48_DOI_output_var_read_len//summary.txt")
Hsiao_cux2_rmats$ID = "Hisao_Cux2_DOI_no_novel"
sum_rmats = rbind(Hsiao_cux2_rmats_novel, Hsiao_cux2_rmats)


Hsiao_rbp4_rmats_novel = read.delim("results/rmats/Rbp4_48_DOI_output_novel_var_read_len/summary.txt")
Hsiao_rbp4_rmats_novel$ID = "Hisao_Rbp4_DOI_novel"
Hsiao_rbp4_rmats = read.delim("results/rmats/Rbp4_48_DOI_output_var_read_len/summary.txt")
Hsiao_rbp4_rmats$ID = "Hisao_Rbp4_DOI_no_novel"
sum_rmats = rbind(Hsiao_rbp4_rmats_novel, Hsiao_rbp4_rmats)


Hsiao_pv_psilo_rmats = read.delim("/media/data01/gomezlab/tmp_bkp/1735-2728-09242024_140206/PV_48h/processed/rMATs/PV_48h_Psilo_rMATS/summary.txt")
Hsiao_pv_psilo_rmats$ID = "Hisao_PV_Psilo"
Hsiao_pv_doi_rmats = read.delim("results/rmats/PV_48_DOI_output_var_read_len//summary.txt")
Hsiao_pv_doi_rmats$ID = "Hisao_PV_DOI"
sum_rmats = rbind(Hsiao_pv_psilo_rmats, Hsiao_pv_doi_rmats)


# ### Nardou
# Nardou_Cocaine_rmats = read.delim("/apps/rmats-turbo/Nardou_48hrs_Cocaine_output/summary.txt")
# Nardou_Cocaine_rmats$ID = "Nardou_Cocaine"
# Nardou_LSD_rmats = read.delim("/apps/rmats-turbo/Nardou_48hrs_LSD_output/summary.txt")
# Nardou_LSD_rmats$ID = "Nardou_LSD"
# Nardou_Ketamine_rmats = read.delim("/apps/rmats-turbo/Nardou_48hrs_katemine_output/summary.txt")
# Nardou_Ketamine_rmats$ID = "Nardou_Ketamine"
# Nardou_MDMA_rmats = read.delim("/apps/rmats-turbo/Nardou_48hrs_MDMA_output/summary.txt")
# Nardou_MDMA_rmats$ID = "Nardou_MDMA"
# 
# 
# ### de la fuente revenga
# dlfr_DOI_rmats = read.delim("/apps/rmats-turbo/de_la_fuente_Revenga_48hrs_output/summary.txt")
# dlfr_DOI_rmats$ID = "de_la_Fuente_Revenga_DOI"


sum_rmats = rbind(Hsiao_cux2_rmats, Hsiao_pv_rmats)
sum_rmats = rbind(sum_rmats, Hsiao_rbp4_rmats)
# sum_rmats = rbind(sum_rmats, Nardou_Cocaine_rmats)
# sum_rmats = rbind(sum_rmats, Nardou_LSD_rmats)
# sum_rmats = rbind(sum_rmats, Nardou_Ketamine_rmats)
# sum_rmats = rbind(sum_rmats, Nardou_MDMA_rmats)
# sum_rmats = rbind(sum_rmats, dlfr_DOI_rmats)
head(sum_rmats)

table(sum_rmats$ID, sum_rmats$EventType)

sum_rmats_filt = sum_rmats[, c('EventType', 'TotalEventsJC', 'ID')]
colnames(sum_rmats_filt)[2] <- "Events"
sum_rmats_filt$Set = 'Total'
sum_rmats_filt$Class = 'JC'

sum_rmats_filt  %>% 
  ggplot( aes(x=ID, y=Events, fill=EventType)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label = Events, group = EventType), 
            position = position_dodge(width = .9),
            vjust = -1, 
            size = 3.7,
            colour="black")+
  scale_fill_viridis(discrete=TRUE, name="") +
  theme_ipsum() +
  labs(title="Total Splicing Events JC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("")

sum_rmats_filt  %>% 
  ggplot( aes(x=EventType, y=Events, fill=ID)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label = Events, group = ID), 
            position = position_dodge(width = .9),
            vjust = -1, 
            size = 3.7,
            colour="black")+            
  #scale_fill_viridis(discrete=TRUE, name="") +
  scale_fill_manual(values = c(Hisao_Rbp4_DOI= "#0f7733", Hisao_PV_DOI="#a94498", Hisao_Cux2_DOI="#45aa9a")) +
  theme_ipsum() +
  labs(title="Total Splicing Events JC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("")

aux = sum_rmats[, c('EventType', 'SignificantEventsJC', 'ID')]
colnames(aux)[2] <- "Events"
aux$Set = 'Significant'
aux$Class = 'JC'


aux  %>% 
  ggplot( aes(x=ID, y=Events, fill=EventType)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label = Events, group = EventType), 
            position = position_dodge(width = .9),
            vjust = -1, 
            size = 3.7,
            colour="black")+      
  scale_fill_viridis(discrete=TRUE, name="") +
  theme_ipsum() +
  labs(title="Significant Splicing Events JC (FDR 5%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("")



aux  %>% 
  ggplot( aes(x=EventType, y=Events, fill=ID)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label = Events, group = ID), 
            position = position_dodge(width = .9),
            vjust = -1, 
            size = 3.7,
            colour="black")+      
  scale_fill_manual(values = c(Hisao_Rbp4_DOI= "#0f7733", Hisao_PV_DOI="#a94498", Hisao_Cux2_DOI="#45aa9a")) +
  theme_ipsum() +
  labs(title="Significant Splicing Events (FDR 5%)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("")





####

sum_rmats_filt = sum_rmats[, c('EventType', 'TotalEventsJCEC', 'ID')]
colnames(sum_rmats_filt)[2] <- "Events"
sum_rmats_filt$Set = 'Total'
sum_rmats_filt$Class = 'JCEC'
aux = sum_rmats[, c('EventType', 'SignificantEventsJCEC', 'ID')]
colnames(aux)[2] <- "Events"
aux$Set = 'Significant'
aux$Class = 'JCEC'

sum_rmats_plot = rbind(sum_rmats_filt, aux)
sum_rmats_plot = sum_rmats_plot[sum_rmats_plot$Set == "Significant",]

sum_rmats_plot  %>% 
  ggplot( aes(x=ID, y=Events, fill=EventType)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label = Events, group = EventType), 
            position = position_dodge(width = .9),
            vjust = -1, size = 3.7)+
  scale_fill_viridis(discrete=TRUE, name="") +
  theme_ipsum_rc() +
  labs(title="Significant Splicing Events JCEC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("")
