#Title: Hobo Logger Tank Temperatures
#Author: KH Wong
#Date Last Modified: 20201001
#See Readme file for details

rm(list=ls()) #clears workspace 

library(lubridate)
library(tidyr)
library(plotrix)
library(tidyverse)
library(dplyr)
library(reshape2)
library(gridExtra)
library(FSA)
library(weathermetrics)
library(ggplot2)
library(Rmisc)

##### Experimental #####
Tank5 <- read.csv("data/Daily_Measurements/20190903_5_Mesocosm.csv", sep=",", skip=c(2), header=FALSE, na.strings = "NA")[ ,2:3]
Tank6 <- read.csv("data/Daily_Measurements/20190903_6_Mesocosm.csv", sep=",", skip=c(2), header=FALSE, na.strings = "NA")[ ,2:3]
Tank6 <- Tank6[1:nrow(Tank5),] #Makes data sets equal in length

data <- cbind(Tank5[,c(1,2)], Tank6$V3)
colnames(data) <- c("Date.Time", "Tank5", "Tank6")
data <- data [1:(nrow(data )-250),]
data$Date.Time <- parse_date_time(data$Date.Time, "%m/%d/%y %I:%M:%S %p")

tmp.col <- c("coral","blue")
tnks <- c("Tank 5", "Tank 6")


# Making a date column to merge lunar day 
data$date <-as.Date(data$Date.Time, "%Y-%m-%d")

exp.day <- read.csv("data/Metadata/experiment_day.csv")
exp.day$date <-as.Date(exp.day$date, "%Y-%m-%d")

data.exp <- merge(data, exp.day, by = "date")
data.exp$Experimental.Day <- as.numeric(data.exp$Experimental.Day)

data.exp_long <- gather(data.exp, Treatment, Temperature, Tank5:Tank6, factor_key=TRUE)

temp.sum <- summarySE(data.exp_long, measurevar="Temperature", groupvars="Treatment")

t.test(Temperature ~ Treatment, data = data.exp_long)





#Making actual plot 
pdf("output/HoboTempstimeline.pdf", width = 15, height = 6)
#par(mfrow=c(1,3))
par(mar=c(6,6,2,2)) #sets the bottom, left, top and right
plot(data.exp$Experimental.Day, data.exp$Tank5, cex=0.2, col="white",ylim = c(25,33), ylab="Temperature °C", xlab="Experimental Day", las=1, xaxt = "n")
#axis(1, data.lunar$lunar.day, labels = T, cex.axis = 0.7, las = 2)
axis(1, at=seq(min(data.exp$Experimental.Day), max(data.exp$Experimental.Day), 1), labels = T)
#axis.Date(1, at=seq(min(data$Date.Time), max(data$Date.Time), by="1 d"), format="%m/%d")
#legend(data$Date.Time[5], 31, legend = "tnks", col=tmp.col, cex=0.9, lty=1, box.lty=0)
dev.off()

#Making actual plot 
pdf("output/HoboTempsExperimental.pdf", width = 15, height = 6)
#par(mfrow=c(1,3))
par(mar=c(6,6,2,2)) #sets the bottom, left, top and right
plot(data.exp$Date.Time, data.exp$Tank5, cex=0.2, col="blue",ylim = c(25,33), ylab="Temperature °C", xlab="Experimental Day", las=1, xaxt = "n")
#axis(1, data.lunar$lunar.day, labels = T, cex.axis = 0.7, las = 2)
axis(1, at=seq(min(data.exp$Experimental.Day), max(data.exp$Experimental.Day), 1), labels = T)
#axis.Date(1, at=seq(min(data$Date.Time), max(data$Date.Time), by="1 d"), format="%m/%d")
points(data.exp$Date.Time, data.exp$Tank6, cex=0.2, col="coral")
#legend(data$Date.Time[5], 31, legend = "tnks", col=tmp.col, cex=0.9, lty=1, box.lty=0)
dev.off()
