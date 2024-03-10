library(ggplot2)
library(tidyverse)
library(stringr)
library(ggrepel)

rm(list = ls())

df <- read.csv("gProfiler_output.csv")

df <- subset(df, source == "KEGG")

df <- df %>% rowwise() %>%
  mutate(label = ifelse(str_detect(term_name, "infection|Infection"), term_name, NA),
         coloring = ifelse(str_detect(term_name, "infection|Infection"), "#009100", "black"))

p <- ggplot(df, aes(y = negative_log10_of_adjusted_p_value, x = "KEGG pathway")) +
  geom_vline(xintercept = "KEGG pathway") +
  geom_point(mapping = aes(size = intersection_size, color = coloring), color = df$coloring, alpha = 0.1) +
  coord_flip() +
  xlab("") +
  ylab("-Log10(p)") +
  geom_text_repel(label = df$label, max.overlaps = Inf, box.padding = 1)

p + theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 6)) +
  scale_size_continuous(range = c(0.05,10))

ggsave("Infection_enrichment.pdf", width = 140, height = 40, units = "mm")