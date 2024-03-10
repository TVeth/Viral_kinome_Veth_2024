library(tidyverse)
library(ggplot2)
library(data.table)
library(RColorBrewer)

df <- read.csv("SRM_traces.csv")

df <- df[df$Intensity != "",]

count <- 1
colVec <- c()
for(rows in 1:nrow(df)){
  if(val <- df[rows,1] == "Retention Time"){
    count <- count + 1
    colVec <- append(colVec, count)
  }else{
    colVec <- append(colVec, count)
  }
}

df$peakID <- colVec

df <- df[df$Retention.Time != "Retention Time",]

df <- mutate(df, Retention.Time = as.double(Retention.Time),
             Intensity = as.double(Intensity))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = rep(col_vector, 10)

p <- ggplot(df, aes(x = Retention.Time, y = Intensity, group = as.factor(peakID))) + 
  geom_area(aes(fill = as.factor(peakID)), show.legend = FALSE, alpha = 0.2) + 
  geom_line(aes(color = as.factor(peakID)), show.legend = FALSE, size = 0.1) +
  xlab("Retention Time (min)")

p + theme_bw() + 
  theme(text = element_text(size = 6),
        panel.grid = element_blank()) +
  scale_color_manual(values = col_vector[1:length(unique(df$peakID))]) +
  scale_fill_manual(values = col_vector[1:length(unique(df$peakID))]) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7e6))

#greys
p + theme_bw() + 
  theme(text = element_text(size = 6),
        panel.grid = element_blank()) +
  scale_color_manual(values = gray.colors(length(unique(df$peakID)), start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL, rev = FALSE)) +
  scale_fill_manual(values = gray.colors(length(unique(df$peakID)), start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL, rev = FALSE)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7e6))

#color
colfunc <- colorRampPalette(c("#BF3B00", "#000083","#009100", "#FFCC00"))

p + theme_bw() + 
  theme(text = element_text(size = 6),
        panel.grid = element_blank()) +
  scale_color_manual(values = colfunc(length(unique(df$peakID)))) +
  scale_fill_manual(values = colfunc(length(unique(df$peakID)))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7e6))

ggsave("MS_traces_SRM_greenredblueyellow.pdf", width = 100, height = 45, units = c("mm"))

#The following part is to identify JNK2 peaks
#JNK2 phos 1 = peakID 38
#JNK2 phos 9 = peakID 148

df <- df %>% rowwise() %>%
  mutate(label = ifelse(peakID == 38, "JNK2_1",
                        ifelse(peakID == 148, "JNK2_9", NA)))

df <- df %>% rowwise() %>%
  mutate(coloring = ifelse(peakID == 91, "red",
                        ifelse(peakID == 188, "green", "black")))

dfcolor <- df %>% distinct(peakID, .keep_all = TRUE)
         
p <- ggplot(df, aes(x = Retention.Time, y = Intensity, group = as.factor(peakID), color = as.factor(peakID)), show.legend = FALSE) + 
  #geom_area(aes(fill = as.factor(peakID)), show.legend = FALSE) + 
  geom_line(size = 0.1, show.legend = FALSE) +
  xlab("Retention Time (minutes)")

p + theme_bw() + 
  theme(text = element_text(size = 6),
        panel.grid = element_blank()) +
  scale_color_manual(values = dfcolor$coloring) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7e6))
ggsave("MS_traces_SRM_color.pdf", width = 100, height = 45, units = c("mm"))
