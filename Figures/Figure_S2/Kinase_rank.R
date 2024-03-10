library(ggplot2)
library(tidyverse)
library(stringr)
library(data.table)

df <- fread("proteinGroups.txt")
dfKinases <- read.csv("allKinases.csv")
dfAll <- read.csv("All_Proteins.csv")
dfFound <- read.csv("lookUpLN_pt2.csv")

df <- subset(df, Reverse != "+")
df <- subset(df, `Only identified by site` != "+")
df <- subset(df, `Potential contaminant` != "+")

#Vectors with the kinases
includedKinaseVec <- unique(dfAll$Protein.Accession)
foundKinaseVec <- unique(dfFound$Protein.Name)
foundKinaseVec <- as.vector(unlist(sapply(foundKinaseVec, function(x) str_split(x, "\\|", simplify = TRUE)[2])))
allKinaseVec <- dfKinases$Kinases

#iBAQ plot
df <- df %>% rowwise() %>%
  mutate(Kinase = ifelse(str_detect(`Protein IDs`, paste(allKinaseVec, collapse = "|")), "Yes", "No"))

dfFig <- data.frame(Protein = df$`Gene names`, Quantification = df$iBAQ, ID = df$`Protein IDs`)

dfFig$color <- "grey"
  
dfFig <- dfFig %>% rowwise() %>%
  mutate(color = ifelse(grepl(ID, paste(allKinaseVec, collapse = "|")), "black", color))
  
dfFig <- dfFig %>% rowwise() %>%
  mutate(color = ifelse(grepl(ID, paste(includedKinaseVec, collapse = "|")), "#BF3B00", color))

dfFig <- dfFig %>% rowwise() %>%
  mutate(color = ifelse(grepl(ID, paste(foundKinaseVec, collapse = "|")), "#009100", color))

dfFig <- dfFig %>% rowwise() %>%
  mutate(size = ifelse(color %in% c("#BF3B00", "#009100"), 1.5, 
                       ifelse(color == "black" , 1, 0.5)))

dfFig <- dfFig %>% arrange(desc(Quantification))

dfFig$Rank <- seq(1, nrow(dfFig))

dfFig$Quantification = log(dfFig$Quantification, 10)

p <- ggplot(dfFig, aes(x = Rank, y = Quantification, color = color)) +
  geom_point(size = dfFig$size, color = dfFig$color) +
  ylab("Log10 Protein Abundance")

p + theme_bw() +
  theme(text = element_text(size=6),
        panel.grid = element_blank())

ggsave("proteome_and_kinases.pdf", width = 50, height = 50, units = "mm")
