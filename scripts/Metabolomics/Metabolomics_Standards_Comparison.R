#Read in required libraries
library("reshape")
#library(plyr)
library("dplyr")
library("tidyverse")
library("ggplot2")
library("arsenal")
library("Rmisc")
library(gridExtra)
library(ggpubr)
library(factoextra)
library(ropls)
library(mixOmics)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(car)
library(effects)
library(ggfortify)
library(cowplot)
library(vegan)
library(corrr)
library(ggcorrplot)
library(GGally)
library(broom)
library(cowplot)
library(RVAideMemoire)
library(arsenal)
library(patchwork)
library(tidyr)
library(ggrepel)
library(MetaboAnalystR)

#Import data

sample.info <- read.csv("data/Metabolomics/Porites-Kevin_SampleInfo.csv")
peak.pos <- read.csv("data/Metabolomics/Peaks_Pos.csv")
peak.neg <- read.csv("data/Metabolomics/Peaks_Neg.csv")
peak.unknown <- read.csv("data/Metabolomics/Peaks_Unknown.csv")

#Removing unncessesary columns

peak.pos2 <- peak.pos[ -c(1:9)] 
peak.pos3 <- peak.pos2[ -c(2:9)] # if we need the blanks, change the 9 to 6

peak.neg2 <- peak.neg[ -c(1:9)] 
peak.neg3 <- peak.neg2[ -c(2:9)] # if we need the blanks, change the 9 to 6

peak.unknown2 <- peak.unknown[ -c(1:9)] 
peak.unknown3 <- peak.unknown2[ -c(2:9)] # if we need the blanks, change the 9 to 6

# Selecting columns with standards for each dataset

pos.standards <- peak.pos3 %>% 
  filter(grepl("PosIS", compound, fixed=TRUE))
pos.standards$compound.ID <- c('Alanine', 'Glutamine', 'Inosine', 'Serine', 'Lysine')

neg.standards <- peak.neg3 %>% 
  filter(grepl('NegIS|SIM', compound)) 

neg.standards$compound.ID <- c('NAD+', 'NADH', 'NADP+', 'NADPH', 'Alanine', 'Glucose','Glutamine', 'Glycine', 'Inosine', 'Lysine', 'Malate', 'Serine', 'Thymine')

# Selecting corresponding metabolites

pos.metabs.list <- c('Alanine', 'Glutamine', 'Inosine', 'Serine', 'Lysine')
pos.metabs <- peak.pos3[peak.pos3$compound %in% pos.metabs.list,]
pos.metabs$compound.ID <- pos.metabs$compound

neg.metabs.list <- c('NAD+', 'Alanine', 'Glucose','Glutamine', 'Glycine', 'Inosine', 'Malate', 'Thymine','Serine', 'Lysine')
neg.metabs <- peak.neg3[peak.neg3$compound %in% neg.metabs.list,]
neg.metabs$compound.ID <- neg.metabs$compound

# Reshaping datasets and adding columns

pos.standards.2 <- melt(pos.standards, id= c("compound.ID", "compound")) #melting dataset 
pos.standards.2$standard.count <- pos.standards.2$value
pos.standards.3 <- pos.standards.2 %>%
  dplyr::select(compound.ID, variable, standard.count)

pos.metabs.2 <- melt(pos.metabs, id= c("compound.ID", "compound"))#melting dataset 
pos.metabs.2$norm.count <- pos.metabs.2$value
pos.metabs.3 <- pos.metabs.2 %>%
  dplyr::select(compound.ID, variable, norm.count)

neg.standards.2 <- melt(neg.standards, id=c("compound.ID", "compound")) #melting dataset 
neg.standards.2$standard.count <- neg.standards.2$value
neg.standards.3 <- neg.standards.2 %>%
  dplyr::select(compound.ID, variable, standard.count)

neg.metabs.2 <- melt(neg.metabs, id= c("compound.ID", "compound")) #melting dataset 
neg.metabs.2$norm.count <- neg.metabs.2$value
neg.metabs.3 <- neg.metabs.2 %>%
  dplyr::select(compound.ID, variable, norm.count)

# binding datasets

pos.all <- merge(pos.standards.3, pos.metabs.3, by = c("compound.ID", "variable"))
neg.all <- merge(neg.standards.3, neg.metabs.3, by = c("compound.ID", "variable"))

# Adding metadata

sample.info$variable <- sample.info$Sample.ID
sample.info2 <- sample.info %>%
  dplyr::select(variable, Fragment.ID, Time, Treatment)

pos.all.meta <- merge(pos.all, sample.info2, by ="variable")
neg.all.meta <- merge(neg.all, sample.info2, by ="variable")

# Plotting scatterplots

pos.scatter <- ggplot(data = pos.all.meta, aes(x=standard.count, y=norm.count))+
  ylab("Normal counts")+ xlab("Standard Counts") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ compound.ID, scales = "free")

neg.scatter <- ggplot(data = neg.all.meta, aes(x=standard.count, y=norm.count))+
  ylab("Normal counts")+ xlab("Standard Counts") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ compound.ID, scales = "free")
