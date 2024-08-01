rm(list = ls())

library(tidyverse)
library(hrbrthemes)
library(viridis)

Hsiao_cux2_rmats = read.delim("results/rmats/cux2_48_DOI_output_novel/summary.txt")
Hsiao_cux2_rmats$ID = "Hisao_Cux2_DOI"
Hsiao_pv_rmats = read.delim("results/rmats/PV_48_DOI_output_novel/summary.txt")
Hsiao_pv_rmats$ID = "Hisao_PV_DOI"
Hsiao_rbp4_rmats = read.delim("results/rmats/Rbp4_48_DOI_output_novel/summary.txt")
Hsiao_rbp4_rmats$ID = "Hisao_Rbp4_DOI"

### Nardou
Nardou_Cocaine_rmats = read.delim("/apps/rmats-turbo/Nardou_48hrs_Cocaine_output/summary.txt")
Nardou_Cocaine_rmats$ID = "Nardou_Cocaine"
Nardou_LSD_rmats = read.delim("/apps/rmats-turbo/Nardou_48hrs_LSD_output/summary.txt")
Nardou_LSD_rmats$ID = "Nardou_LSD"
Nardou_Ketamine_rmats = read.delim("/apps/rmats-turbo/Nardou_48hrs_katemine_output/summary.txt")
Nardou_Ketamine_rmats$ID = "Nardou_Ketamine"
Nardou_MDMA_rmats = read.delim("/apps/rmats-turbo/Nardou_48hrs_MDMA_output/summary.txt")
Nardou_MDMA_rmats$ID = "Nardou_MDMA"


### de la fuente revenga
dlfr_DOI_rmats = read.delim("/apps/rmats-turbo/de_la_fuente_Revenga_48hrs_output/summary.txt")
dlfr_DOI_rmats$ID = "de_la_Fuente_Revenga_DOI"


sum_rmats = rbind(Hsiao_cux2_rmats, Hsiao_pv_rmats)
sum_rmats = rbind(sum_rmats, Hsiao_rbp4_rmats)
sum_rmats = rbind(sum_rmats, Nardou_Cocaine_rmats)
sum_rmats = rbind(sum_rmats, Nardou_LSD_rmats)
sum_rmats = rbind(sum_rmats, Nardou_Ketamine_rmats)
sum_rmats = rbind(sum_rmats, Nardou_MDMA_rmats)
sum_rmats = rbind(sum_rmats, dlfr_DOI_rmats)
head(sum_rmats)

table(sum_rmats$ID, sum_rmats$EventType)

sum_rmats_filt = sum_rmats[, c('EventType', 'TotalEventsJC', 'ID')]
colnames(sum_rmats_filt)[2] <- "Events"
sum_rmats_filt$Set = 'Total'
sum_rmats_filt$Class = 'JC'
aux = sum_rmats[, c('EventType', 'SignificantEventsJC', 'ID')]
colnames(aux)[2] <- "Events"
aux$Set = 'Significant'
aux$Class = 'JC'

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
  labs(title="Significant Splicing Events JC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  xlab("")


sum_rmats_plot  %>% 
  ggplot( aes(x=EventType, y=Events, fill=ID)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label = Events, group = ID), 
            position = position_dodge(width = .9),
            vjust = -1, size = 3.7)+
  scale_fill_viridis(discrete=TRUE, name="") +
  theme_modern_rc() +
  labs(title="Significant Splicing Events JCEC") +
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
