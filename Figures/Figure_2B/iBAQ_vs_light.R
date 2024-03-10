library(ggplot2)
library(tidyverse)

retrievedf1 <- function(){
  skylineExport <- read.csv("Quant_export_CVB3.csv")
  
  skylineExport <- skylineExport %>% filter(Standard.Type != "iRT", Replicate != "20211004_TSQALTIS_UM2_4019253_SA_SNB200921863_mock_CVB3_2AM_HL_pt2") %>% 
    dplyr::select(!Standard.Type)
  
  lookUpLN <- read.csv("lookUpLN_pt2.csv")
  
  skylineExport <- skylineExport %>% left_join(lookUpLN, by = c("Protein.Name", "Peptide", "Peptide.Modified.Sequence"))
  
  print(paste("prefiltered list has ", nrow(skylineExport), " rows and ", length(unique(skylineExport$pSites)), " pSites before filtering.", collapse = ""))
  skylineExport <- skylineExport %>% filter(!str_detect(Peptide.Modified.Sequence, "\\[\\+8\\]|\\[\\+10\\]"))
  print(paste("prefiltered list has ", nrow(skylineExport), " rows and ", length(unique(skylineExport$pSites)), " pSites after filtering.", collapse = ""))
  
  skylineExport$LN_part2 <- NA
  
  return(skylineExport)
}

retrievedf2 <- function(){
  skylineExport <- read.csv("Quant_export_EMCV.csv")
  
  skylineExport <- skylineExport %>% 
    filter(Standard.Type != "iRT", 
           Replicate != "20211004_TSQALTIS_UM2_4019253_SA_SNB200921863_mock_CVB3_2AM_HL_pt2",
           Replicate != "10h_mock_r2") %>% 
    dplyr::select(!Standard.Type)

  lookUpLN <- read.csv("lookUpLN_pt2.csv")
  
  skylineExport <- skylineExport %>% left_join(lookUpLN, by = c("Protein.Name", "Peptide", "Peptide.Modified.Sequence"))
  
  print(paste("prefiltered list has ", nrow(skylineExport), " rows and ", length(unique(skylineExport$pSites)), " pSites before filtering.", collapse = ""))
  skylineExport <- skylineExport %>% filter(!str_detect(Peptide.Modified.Sequence, "\\[\\+8\\]|\\[\\+10\\]"))
  print(paste("prefiltered list has ", nrow(skylineExport), " rows and ", length(unique(skylineExport$pSites)), " pSites after filtering.", collapse = ""))
  
  skylineExport$LN_groups <- NA
  
  return(skylineExport)
}

iBAQdf <- function(){
  dfibaq <- read.delim("proteinGroups.txt")
  
  dfFig <- data.frame(Kinase = dfibaq$Gene.names, Quantification = dfibaq$iBAQ, ID = dfibaq$Protein.IDs)
  
  return(dfFig)
}

df <- rbind(retrievedf1(), retrievedf2())[c("pSites", "Total.Area", "Protein.Name")]
dfibaq <- iBAQdf()

df <- df %>% mutate(Total.Area = as.integer(Total.Area)) %>%
  arrange(desc(Total.Area)) %>%
  dplyr::distinct(pSites, .keep_all = TRUE)

df <- df %>% rowwise() %>%
  mutate(ID = str_split(Protein.Name, "\\|", simplify = TRUE)[2]) %>%
  ungroup

df <- df %>% left_join(dfibaq)

df[is.na(df$Quantification),]$Quantification <- 2
df[df$Quantification == 0,]$Quantification <- 1
df[df$Quantification == 2,]$Quantification <- 0

#Colorcoding
df <- df %>% rowwise() %>%
  mutate(coloring = case_when(Quantification == 1 ~ "Not quantifiable",
                              Quantification == 0 ~ "Below detection limit",
                              TRUE ~ "Quantified"))

df$Total.Area <- log(df$Total.Area, 10)
df$Quantification <- log(df$Quantification, 10)

p <- ggplot(df, aes(x = Total.Area, y = Quantification, color = coloring)) +
  geom_point(mapping = aes(shape = coloring), size = 0.6, fill = "white", show.legend = FALSE) +
  ylab("Protein abundance (log10)") +
  xlab("Intensity kinome assay (log10)")

p + theme_bw() +
  theme(text = element_text(size = 8),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("#BF3B00", "#000083","#009100"))

ggsave("iBAQ_vs_area_update.pdf", width = 45, height = 45, units = c("mm"))