#Title: Surface area calculations for Porites Nurition July 2019 Project (Rim)
#Author: HM Putnam and NJ Silbiger
#Edited by: Kevin Wong
#Date Last Modified: 20190711
#See Readme file for details 

library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyverse)


### 2017 colony surface area ###

### 2018 colony surface area ###
coral.sa.Rim <-read.csv("data/July_Rim/Surface_Area/rim_colony_SA_raw.csv", header=T, sep = ",")
standard <- read.csv("data/July_Rim/Surface_Area/Foil_Standard.csv", header=T, sep = ",")

#Standard calculation
ggplot(data = standard, aes(x=SA.cm2, y=weight.g))+
  ylab("Weight (g)")+ xlab("Surface Area (cm2)") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

lmstandard <- lm (SA.cm2 ~ weight.g, data = standard)
summary(lmstandard)

# Calculate SA with standard lm

coral.sa.Rim$SA.cm2 <- predict(lmstandard, newdata = coral.sa.Rim) #using model to get concentration
write.csv(coral.sa.Rim, 'data/July_Rim/Surface_Area/colony.Rim.calcsa.csv')
