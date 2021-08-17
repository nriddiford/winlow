library(tidyverse)
setwd("~/GitHub/winlow/")

#' Set a clean ggplot2
#' @keywords theme
#' @import ggplot2
cleanTheme <- function(base_size = 12, hide_legend = FALSE, rotate_x = FALSE) {
  
  mod_theme <-
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5, size = 15),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5),
      axis.text = element_text(size = base_size),
      axis.title = element_text(size = base_size + base_size/10),
      strip.text = element_text(size = base_size),
      text = element_text(family = "Helvetica"),
      strip.background = element_rect(color = "#232323", fill = "#0508382C", size = .5, linetype = "solid")
    )
  if(hide_legend){
    mod_theme <- mod_theme + ggplot2::theme(legend.position="none")
  }
  if(rotate_x){
    mod_theme <- mod_theme + ggplot2::theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = base_size)
    )
  }
  return(mod_theme)
}


ec_stats <- read.delim("stats.csv", sep=',', header = F)
colnames(ec_stats) <- c('sample', 'tissue', 'total_reads', 'locus', 'n_samples', 'av_read_count', 'fc', 'log2_FC', 'p-val', 'sig')

ec_stats <- ec_stats %>%
  dplyr::filter(tissue == 'tumour') %>%
  dplyr::mutate(contamination = ifelse(abs(log2_FC) >= 1, 'high',
                                       ifelse(abs(log2_FC) >= 0.32, 'medium', 'low'))) %>%
  dplyr::mutate(group = substring(sample, 1, 1))


red <- "#FC4E07"
blue <- "#00AFBB"
yellow <- "#E7B800"

colours <- c(
  "low" = "#00AFBB",
  "medium" = "#E7B800",
  "high" = "#FC4E07"
)

ec_stats %>% 
  ggplot(.) +
  dplyr::arrange(-contamination) %>% 
  geom_bar(aes(sample, log2_FC, colour = contamination, fill = contamination), alpha = 0.5, stat = "identity") +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  scale_y_continuous("Log2(FC)") +
  cleanTheme(rotate_x = 'T', ) +
  theme(axis.title.x = element_blank()) +
  ggtitle("Log2 ratio of tumour/normal read counts in 3R:0-3000000")
