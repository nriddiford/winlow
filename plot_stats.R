library(tidyverse)
setwd("~/GitHub/winlow/")
ec_stats <- read.delim("stats.csv", sep=',', header = F)
colnames(ec_stats) <- c('sample', 'tissue', 'total_reads', 'locus', 'n_samples', 'av_read_count', 'fc', 'log2_FC', 'p-val', 'sig')

ec_stats <- ec_stats %>%
  dplyr::filter(tissue == 'tumour') %>%
  dplyr::mutate(contamination = ifelse(abs(log2_FC) >= 1, '>50%',
                                       ifelse(abs(log2_FC) >= 0.32, '>25%', '<25%'))) %>%
  dplyr::mutate(group = substring(sample, 1, 1))


red <- "#FC4E07"
blue <- "#00AFBB"
yellow <- "#E7B800"

library(ggpubr)

ggbarplot(ec_stats, "sample", "log2_FC", sort.val = 'desc', sort.by.groups = F,
          fill = 'contamination', color = 'contamination', palette = c(blue, yellow, red), alpha = 0.6,
          rotate = TRUE, x.text.angle = 90, title="Log2 ratio of tumour/normal read counts in 3R:0-3000000", ylab="Log2(FC)", xlab=F)
