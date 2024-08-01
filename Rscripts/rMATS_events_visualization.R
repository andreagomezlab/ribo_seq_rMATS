library(hrbrthemes)
library(ggVennDiagram)
library("VennDiagram")

plot_table <- function(table) {
  grid.newpage()
  return(gridExtra::grid.table(table))
}

plot_events_filtering <- function(res, title_txt) {
  p1 = ggplot(res, aes(avgPSI_1, avgPSI_2, color=passed))  +
    geom_point(aes(shape = status ),
               stroke = 2
    ) +
    theme_ipsum(base_size = 14) +
    labs(colour = "average PSI > 0.8 & InLevellDiff > |0.1|", shape="IncLevelDifference", title = paste(title_txt, table(res$passed)[2], "of", table(res$passed)[1], "events passed", sep = "_"))
  return(p1)
}


plot_venn <- function(x) {
  p1 = ggVennDiagram(x) +
    scale_color_brewer(palette = "Paired")
  return(p1)
}


