library(tidyverse)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(factoextra)

rm(list=ls(all=TRUE))

setwd("") #Set your work directory containing the input files

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

generateHeatmap <- function(save = FALSE){
  #This function creates a heatmap of the individual FC values
  sigValues <- skylineExport %>% filter(anovaRes <= 0.05)
  sigValues <- unique(sigValues$pSites)
  
  heatMapdf <- skylineExport %>% ungroup() %>% dplyr::select(pSites, normSignal, Replicate.Name) %>% filter(pSites %in% sigValues)
  
  #convert FC to log2
  heatMapdf$normSignal <- log(heatMapdf$normSignal, 2)
  
  heatMapmtrx <- heatMapdf %>% pivot_wider(names_from = Replicate.Name, values_from = normSignal)
  
  #To get things in the right order manually select the correct columns
  heatMapmtrx <- heatMapmtrx %>% dplyr::select(pSites,
                                        CVB3_4h_r1, CVB3_4h_r2, CVB3_4h_r3,
                                        CVB3_6h_r1, CVB3_6h_r2, CVB3_6h_r3,
                                        CVB3_8h_r1, CVB3_8h_r2, CVB3_8h_r3,
                                        CVB3_10h_r1, CVB3_10h_r2, CVB3_10h_r3,
                                        `2AM_4h_r1`, `2AM_4h_r2`, `2AM_4h_r3`,
                                        `2AM_6h_r1`, `2AM_6h_r2`, `2AM_6h_r3`,
                                        `2AM_8h_r1`, `2AM_8h_r2`, `2AM_8h_r3`,
                                        `2AM_10h_r1`, `2AM_10h_r2`, `2AM_10h_r3`,
                                        `4h_EMCV_r1`, `4h_EMCV_r2`, `4h_EMCV_r3`,
                                        `6h_EMCV_r1`, `6h_EMCV_r2`, `6h_EMCV_r3`,
                                        `8h_EMCV_r1`, `8h_EMCV_r2`, `8h_EMCV_r3`,
                                        `10h_EMCV_r1`, `10h_EMCV_r2`, `10h_EMCV_r3`,
                                        `4h_L_r1`, `4h_L_r2`, `4h_L_r3`,
                                        `6h_L_r1`,  `6h_L_r2`,  `6h_L_r3`,
                                        `8h_L_r1`, `8h_L_r2`, `8h_L_r3`,
                                        `10h_L_r1`, `10h_L_r2`, `10h_L_r3`)
  
  #Import lookuptable to join
  lookUpLN <- read.csv("lookUpLN_pt2.csv")
  lookUpLN <- lookUpLN[c("pSites", "Pathway")] %>% distinct(pSites, .keep_all = TRUE)
  
  decoration <- heatMapmtrx["pSites"] %>% left_join(lookUpLN, by = c("pSites"))
  
  #Change the missing values to 0 to allow clustering
  heatMapmtrx[is.na(heatMapmtrx)] <- 0
  
  #Generate the actual matrix
  mtrx <- data.matrix(heatMapmtrx[2:ncol(heatMapmtrx)])
  rownames(mtrx) <- heatMapmtrx$pSites
  
  #Create a color gradient
  col <- colorRampPalette(c("#b35806", "white", "#542788"), bias = 1.1)(256)
  col_fun = colorRamp2(c(seq(3, 0.1, length.out = 128), seq(0, -3, length.out = 128)), col)
  #col_fun = colorRamp2(c(-6, 0, 0, 0, 6), c("#b35806", "white", "white", "white", "#542788"))
  
  #Split the columns based on the separate groups
  colSplit <- c(rep("CVB3", 12), rep("CVB3-2AM", 12), rep("EMCV", 12), rep("EMCV-Lzn", 12))
  
  #Annotations
  mapAnnotation <- rowAnnotation(foo2 = decoration$Pathway,
                                 col = list("PI3K signaling" = "#1B9E77",
                                            "MAPK signaling" = "#D95F02",
                                            "AGC kinase signaling" = "#7570B3",
                                            "Cell cycle progression" = "#E7298A",
                                            "Alternative" = "#66A61E",
                                            "Cytoskeleton regulation" = "#E6AB02",
                                            "DNA damage" = "#A6761D"))
  
  set.seed(2)
  #Generate the heatmap
  p <- Heatmap(mtrx, cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, column_split = colSplit,
               cluster_row_slices = TRUE, 
               cluster_column_slices = TRUE,
               border = TRUE,
               right_annotation = mapAnnotation,
               row_km = 4)
  
  print(p)
  
  dir.create("heatmaps", showWarnings = FALSE)
  
  if(save){
    pdf("heatmaps/heatmap_all_FC_Kmeans.pdf", width = 11, height = 11) 
  }
  print(p)
  if(save){
    dev.off()
  }

}

#This function imports both the CVB3 and EMCV data
skylineExport <- rbind(retrieveFun1(), retrieveFun2())

#This function generates the heatmap of figure 
generateHeatmap(save = FALSE) #Speciy whether the Heatmap should be exported