#Title: Lipid Assay 
#Author: KH Wong
#Date Last Modified: 20210623
#See Readme file for details 


##### 20210622 Lipid Testing #####

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(Rmisc)

# Read in datafiles

Initial <- read.csv("data/Physiology/Lipid/20210622/20210622_Lipid_Initial.csv")
Final <- read.csv("data/Physiology/Lipid/20210622/20210622_Lipid_Final.csv")
Meta <- read.csv("data/Physiology/Lipid/20210622/20210622_Lipid_Well_Meta.csv")

# Merging Files and renaming columns

data.1 <- merge(Meta, Initial, by = "Well")
names(data.1)[11] <- "Initial.abs"

data.2 <- merge(data.1, Final, by = "Well")
names(data.2)[12] <- "Final.abs"

#### Removing the third triplicate from each samples -- this was due to errors in pipetting 

data.2 <- data.2 %>%
  filter(Replicate != "3")

# Blank correction for each run separately

Initial_Blank <- data.2 %>% 
  filter(Type == "Blank") %>%
  summarise(blk.avg = mean(Initial.abs))

Final_Blank <- data.2 %>% 
  filter(Type == "Blank") %>%
  summarise(blk.avg = mean(Final.abs))

data.2$i.abs.corr <- data.2$Initial.abs - Initial_Blank$blk.avg
data.2$f.abs.corr <- data.2$Final.abs - Final_Blank$blk.avg

# this code makes sure that any negative values on plate are converted to 0

data.2$i.abs.corr[data.2$i.abs.corr < 0] <- 0
data.2$f.abs.corr[data.2$f.abs.corr < 0] <- 0

# Subtract Final from Initial 

data.2$sub.abs <- data.2$f.abs.corr - data.2$i.abs.corr

# Standard curve 

Standard <- data.2 %>% 
  filter(Type == "Standard") 

mean.std <- summarySE(Standard, measurevar="sub.abs", groupvars="Concentration.mg.mL")

Standard.plot <- ggplot(data = mean.std, aes(x=Concentration.mg.mL, y=sub.abs))+
  ylab("Absorbance (nm)")+ xlab("Concentration (mg/mL)") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

lmstandard <- lm (Concentration.mg.mL ~ sub.abs, data = mean.std)
lmsummary <- summary(lmstandard) # R2 = 0.8397 

# Standard curve without highest concentration 

mean.std.2 <- mean.std[-c(6),]

Standard.plot2 <- ggplot(data = mean.std.2, aes(x=Concentration.mg.mL, y=sub.abs))+
  ylab("Absorbance (nm)")+ xlab("Concentration (mg/mL)") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

lmstandard2 <- lm (Concentration.mg.mL ~ sub.abs, data = mean.std.2)
lmsummary2 <- summary(lmstandard2) # R2 = 0.9565 

# Obtaining concentration values from standard curve

Samples <- data.2 %>% #subsetting Samples
  filter(Type == "Coral") 

Samples$Concentration.mg.mL <- predict(lmstandard2, newdata = Samples) #using model to get concentration

