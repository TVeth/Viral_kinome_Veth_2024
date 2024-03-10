library(tidyverse)
library(ggplot2)

JNK1df <- read.csv("JNK2_phos_1.csv")

JNK1df <- JNK1df %>%
  pivot_longer(cols = colnames(JNK1df)[2:ncol(JNK1df)], values_to = "Intensity", names_to = "Transition")

JNK1df[JNK1df$Retention.Time < 73,]$Intensity <- 0
JNK1df[JNK1df$Retention.Time > 74.2,]$Intensity <- 0

p <- ggplot(JNK1df, aes(x = Retention.Time, y = Intensity, color = Transition)) +
  geom_line() +
  xlab("Retention Time (hours)")

p + theme_bw() +
  theme(text = element_text(size = 6),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100000)) +
  scale_x_continuous(expand = c(0,0), limits= c(72, 75))

ggsave("JNK2_phos1.pdf", width = 50, height = 45, units = c("mm"))

#JNK
JNK9df <- read.csv("JNK2_phos_9.csv")

JNK9df <- JNK9df %>%
  pivot_longer(cols = colnames(JNK9df)[2:ncol(JNK9df)], values_to = "Intensity", names_to = "Transition")

JNK9df[JNK9df$Retention.Time < 67.2,]$Intensity <- 0
JNK9df[JNK9df$Retention.Time > 68.4,]$Intensity <- 0

p <- ggplot(JNK9df, aes(x = Retention.Time, y = Intensity, color = Transition)) +
  geom_line() +
  xlab("Retention Time (hours)")

p + theme_bw() +
  theme(text = element_text(size = 6),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100000)) +
  scale_x_continuous(expand = c(0,0), limits= c(66.75, 68.75),
                     breaks = c(67, 68), labels = c(67, 68))

ggsave("JNK2_phos9.pdf", width = 50, height = 45, units = c("mm"))