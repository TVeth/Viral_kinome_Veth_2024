library(tidyverse)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(factoextra)

rm(list=ls(all=TRUE))

setwd("") #Set your work directory that contains the input files

retrieveFun1 <- function(){
  skylineExport <- read.csv("Quant_export_CVB3.csv")
  
  skylineExport <- skylineExport %>% filter(Standard.Type != "iRT", Replicate != "20211004_TSQALTIS_UM2_4019253_SA_SNB200921863_mock_CVB3_2AM_HL_pt2") %>% 
    dplyr::select(!Standard.Type)
  
  skylineExport <- skylineExport %>% rowwise() %>% mutate(Time = str_sub(str_split(Replicate, "_", simplify = TRUE)[2],1, -2),
                                                          Group = paste0(str_split(Replicate, "_", simplify = TRUE)[1:2], collapse="_"),
                                                          Rep = str_sub(Replicate, -1),
                                                          Grouping = str_split(Replicate, "_", simplify = TRUE)[1])
  
  skylineExport <- mutate(skylineExport, Time = as.integer(Time))

  lookUpLN <- read.csv("lookUpLN.csv")
  
  skylineExport <- skylineExport %>% left_join(lookUpLN, by = c("Protein.Name", "Peptide", "Peptide.Modified.Sequence"))
  
  print(paste("prefiltered list has ", nrow(skylineExport), " rows and ", length(unique(skylineExport$pSites)), " pSites before filtering.", collapse = ""))
  skylineExport <- skylineExport %>% filter(!str_detect(Modified.Sequence, "\\[\\+8\\]|\\[\\+10\\]"))
  print(paste("prefiltered list has ", nrow(skylineExport), " rows and ", length(unique(skylineExport$pSites)), " pSites after filtering.", collapse = ""))
  
  #Delete other peptides
  skylineExport <- skylineExport %>% filter(pSites != "MAPKAPK2-T222")
  
  skylineExport <- skylineExport %>% distinct(pSites, Replicate.Name, .keep_all = TRUE)
  
  skylineExport <- skylineExport %>% filter(pSites != "MAPK14-Y182" & pSites != "MAPK14-T180")
  
  skylineExport <- skylineExport %>% filter(RatioLightToHeavy != "#N/A")
  
  skylineExport <- skylineExport %>% filter(!pSites %in% c("MAPK3-T202", "NEK6-S206"))
  
  skylineExport <- mutate(skylineExport, RatioLightToHeavy = as.double(RatioLightToHeavy))
  skylineExport <- skylineExport %>% group_by(Peptide.Modified.Sequence, Group) %>% mutate(meanHL = mean(RatioLightToHeavy, na.rm = TRUE))
  
  #Center the plots at y=0 at t=2 based on the mock
  centerToZero <- function(x, y, z, q){
    tempdf <- data.frame(Time = x, Signal = y, group = z, meanSignal = q)
    
    t0 <- subset(tempdf, Time == 2 & group == "mock")[1,]$meanSignal
    
    return(t0)
  }
  
  skylineExport <- skylineExport %>% group_by(Peptide.Modified.Sequence) %>% mutate(normalizedSignalt0 = centerToZero(x = Time, y = RatioLightToHeavy, z = Grouping, q = meanHL))
  
  skylineExport <- skylineExport %>% rowwise() %>% mutate(normSignal = RatioLightToHeavy / normalizedSignalt0,
                                                          normMeanSignal = meanHL / normalizedSignalt0)

  anoveSig <- function(x, y){
    if (length(x) <= 2){
      return(1)
    }else{
      res.aov <- aov(x ~ y)
      sums <- summary(res.aov)[[1]][["Pr(>F)"]][1]
      #print(res.aov)
      return(sums)}
  }
  
  skylineExport$LN_groups <- factor(skylineExport$LN_groups, levels = unique(skylineExport$LN_groups))
  skylineExport <- skylineExport %>% group_by(Peptide.Modified.Sequence, Time) %>% mutate(anovaRes = anoveSig(x = RatioLightToHeavy, y = LN_groups))
  return(skylineExport)
}

retrieveFun2 <- function(){
  skylineExport <- read.csv("Quant_export_EMCV.csv")
    
    skylineExport <- skylineExport %>% 
      filter(Standard.Type != "iRT", 
             Replicate != "20211004_TSQALTIS_UM2_4019253_SA_SNB200921863_mock_CVB3_2AM_HL_pt2",
             Replicate != "10h_mock_r2") %>% 
      dplyr::select(!Standard.Type)
    
    skylineExport <- skylineExport %>% rowwise() %>% mutate(Time = str_sub(str_split(Replicate, "_", simplify = TRUE)[1],1, -2),
                                                            Group = paste0(str_split(Replicate, "_", simplify = TRUE)[1:2], collapse="_"),
                                                            Rep = str_sub(Replicate, -1),
                                                            Grouping = str_split(Replicate, "_", simplify = TRUE)[2])
    
    skylineExport <- mutate(skylineExport, Time = as.integer(Time))

    lookUpLN <- read.csv("lookUpLN_pt2.csv")
    
    skylineExport <- skylineExport %>% left_join(lookUpLN, by = c("Protein.Name", "Peptide", "Peptide.Modified.Sequence"))
    
    print(paste("prefiltered list has ", nrow(skylineExport), " rows and ", length(unique(skylineExport$pSites)), " pSites before filtering.", collapse = ""))
    skylineExport <- skylineExport %>% filter(!str_detect(Modified.Sequence, "\\[\\+8\\]|\\[\\+10\\]"))
    print(paste("prefiltered list has ", nrow(skylineExport), " rows and ", length(unique(skylineExport$pSites)), " pSites after filtering.", collapse = ""))
    
    #MAPKAPK2-T222 is in the noise
    #skylineExport <- skylineExport %>% filter(pSites != "MAPKAPK2-T222")
    
    #KS6A1-T573;KS6A3-T577 is listed twice. Both show identical profiles
    skylineExport <- skylineExport %>% distinct(pSites, Replicate.Name, .keep_all = TRUE)
    
    #MK14 showed RT shifting
    #skylineExport <- skylineExport %>% filter(pSites != "MAPK14???Y182" & pSites != "MAPK14???T180")
    
    skylineExport <- skylineExport %>%
      dplyr::filter(!pSites %in% c("MAPKAPK2-T222", "FAK1-Y576", "MAPK14-Y182", "MAPK14-Y182-Mox",
                       "MAPK14-T180", "MAPK14-T180-Mox", "MAPK1-Y187", "CaMKID-T180", 
                       "MAPK6-S189"))
    
    #Delete the NAs in the RatioLightToHeavy column
    skylineExport <- skylineExport %>% filter(RatioLightToHeavy != "#N/A")
  
    skylineExport <- mutate(skylineExport, RatioLightToHeavy = as.double(RatioLightToHeavy))
    skylineExport <- skylineExport %>% group_by(Peptide.Modified.Sequence, Group) %>% mutate(meanHL = mean(RatioLightToHeavy, na.rm = TRUE))
    
    centerToZero <- function(x, y, z, q){
      tempdf <- data.frame(Time = x, Signal = y, group = z, meanSignal = q)
      
      t0 <- subset(tempdf, Time == 4 & group == "mock")[1,]$meanSignal
      
      return(t0)
    }
    
    skylineExport <- skylineExport %>% group_by(Peptide.Modified.Sequence) %>% mutate(normalizedSignalt0 = centerToZero(x = Time, y = RatioLightToHeavy, z = Grouping, q = meanHL))
    
    skylineExport <- skylineExport %>% rowwise() %>% mutate(normSignal = RatioLightToHeavy / normalizedSignalt0,
                                                            normMeanSignal = meanHL / normalizedSignalt0)
  
    anoveSig <- function(x, y){
      if (length(x) <= 2){
        return(1)
      }else{
        res.aov <- aov(x ~ y)
        sums <- summary(res.aov)[[1]][["Pr(>F)"]][1]
        #print(res.aov)
        return(sums)}
    }
    
    skylineExport$LN_part2 <- factor(skylineExport$LN_part2, levels = unique(skylineExport$LN_part2))
    skylineExport <- skylineExport %>% group_by(Peptide.Modified.Sequence, Time) %>% mutate(anovaRes = anoveSig(x = RatioLightToHeavy, y = LN_part2))
    return(skylineExport)
}

exportlinePlotsCVB3 <- function(){
  dir.create("lineplots", showWarnings = FALSE)
  
  linePlotsdf <- skylineExport %>% ungroup() %>%
    filter(!is.na(LN_groups)) %>%
    dplyr::select(pSites, normSignal, Time, Grouping) %>%
    filter(Grouping %in% c("CVB3", "2AM", "mock") & Time != 2)

  linePlotsdf$Grouping <- factor(linePlotsdf$Grouping, levels = c("CVB3", "2AM", "mock"))
  
  #CDK9
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "CDK9-S175" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("CDK9-S175") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-3, 2)))
  
  ggsave("lineplots/CDK9.pdf", width = 60, height = 40, units = c("mm"))

  #CDK12
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "CDK12-S1083" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("CDK12-S1083") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-3, 2)))
  
  ggsave("lineplots/CDK12.pdf", width = 60, height = 40, units = c("mm"))
  
  #RSK1/RSK2
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "KS6A1-T573;KS6A3-T577" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("RSK1-T573;RSK2-T577") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.6, 4.5)))
  
  ggsave("lineplots/RSK1.pdf", width = 60, height = 40, units = c("mm"))
  
  #MEK1/2
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "MP2K1-S222;MP2K2-S226" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("MEK1-S222;MEK2-S226") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.6, 4.5)))
  
  ggsave("lineplots/MEK12.pdf", width = 60, height = 40, units = c("mm"))
  
  #CHK2
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "CHK2-S379" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("CHK2-S379") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.6, 4.5)))
  
  ggsave("lineplots/CHK2.pdf", width = 60, height = 40, units = c("mm"))
  
  #ERK2
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "MAPK1-T185/Y187" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("ERK2???T185/Y187") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.6, 4.5)))
  
  ggsave("lineplots/ERK2.pdf", width = 60, height = 40, units = c("mm"))
}

exportlinePlotsEMCV <- function(){
  linePlotsdf <- skylineExport %>% ungroup() %>%
    filter(!is.na(LN_part2)) %>%
    dplyr::select(pSites, normSignal, Time, Grouping) %>%
    filter(Grouping %in% c("EMCV", "L", "mock") & Time != 2)
  
  #CDK9
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "CDK9-S175" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("CDK9-S175") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-3, 2)))
  
  ggsave("lineplots/EMCV_CDK9.pdf", width = 60, height = 40, units = c("mm"))
  
  #CDK12
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "CDK12-S1083" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("CDK12-S1083") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-3, 2)))
  
  ggsave("lineplots/EMCV_CDK12.pdf", width = 60, height = 40, units = c("mm"))
  
  #RSK1/RSK2
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "KS6A1-T573;KS6A3-T577" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("RSK1-T573;RSK2-T577") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.6, 4.5)))
  
  ggsave("lineplots/EMCV_RSK1.pdf", width = 60, height = 40, units = c("mm"))
  
  #MEK1/MEK2
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "MP2K1-S222;MP2K2-S226" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("MEK1-S222;MEK2-S226") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.6, 4.5)))
  
  ggsave("lineplots/EMCV_MEK12.pdf", width = 60, height = 40, units = c("mm"))
  
  #CHK2
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "CHK2-S379" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("CHK2-S379") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.6, 4.5)))
  
  ggsave("lineplots/EMCV_CHK2.pdf", width = 60, height = 40, units = c("mm"))
  
  #ERK2
  linePlotsdfTemp <- linePlotsdf %>% filter(pSites == "MAPK1???T185/Y187" & Grouping != "E2a")
  
  p <- ggplot(linePlotsdfTemp, aes(x = Time, y = log(normSignal, 2), group = Grouping, color = Grouping)) + 
    geom_hline(yintercept = c(-log(1.5,2), log(1.5,2)), linetype = "dashed", color = "black", alpha =0.8, size = 0.1) +
    geom_point(size = 0.25) + 
    geom_smooth(span = 0.6, level = 0.9, se = TRUE, size = 0.1) +
    xlab("Time (h.p.i.)") +
    ylab("Log2 fold change(/mock)")
  
  print(p + theme_bw() + theme(text = element_text(size = 6),
                               panel.grid = element_blank()) +
          ggtitle("ERK2???T185/Y187") +
          scale_color_manual(values = c("#009100", "#BF3B00", "darkgrey")) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.6, 4.5)))
  
  ggsave("lineplots/EMCV_ERK2.pdf", width = 60, height = 40, units = c("mm"))
}

exportAllLineplots <- function(){
  linePlotsdf <- skylineExport %>% ungroup() %>%
    dplyr::select(pSites, normSignal, Time, Grouping, LN_groups) %>%
    filter(Grouping %in% c("EMCV", "CVB3", "mock") & Time != 2) 
  
  linePlotsdf <- subset(linePlotsdf, !is.na(LN_groups) | Grouping != "CVB3")
  
  pdf("lineplots/all.pdf")
  for (pep in unique(linePlotsdf$pSites)){
    p <- ggplot(subset(linePlotsdf, pSites == pep), aes(x = Time, y = normSignal, group = Grouping, color = Grouping)) + geom_point() + geom_smooth(span = 0.6, level = 0.8, se = TRUE)
    print(p + theme_bw() + theme(text = element_text(size = 20),
                                 panel.grid = element_blank()) +
            ggtitle(pep) +
            xlab("Time (hours") +
            ylab("Normalized signal (a.u.)") +
            scale_color_manual(values = c("#D95F02", "#1B9E77", "#7570B3", "grey", "black")))
  }
  dev.off()
}

#This function imports both the CVB3 and EMCV data
skylineExport <- rbind(retrieveFun1(), retrieveFun2())

#This function generates the heatmap of figure 
generateHeatmap(save = FALSE)
exportlinePlotsCVB3()
exportlinePlotsEMCV()
exportAllLineplots()