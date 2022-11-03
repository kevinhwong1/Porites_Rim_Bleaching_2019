#Title: Daily Measurements for Porties Nutrition JULY experiment (RIM )
#Author: KWong (modified from HPutnam)
#Date Modified: 20190715
#See Readme file for details on data files and metadata

rm(list=ls()) # removes all prior objects

#R Version: R version 3.3.1
#RStudio Version: 1.0.44
######Read in required libraries#####
library(car) #version 2.1-4 Date: 2016-11-30 Depends: R (>= 3.2.0) Imports:MASS, mgcv, nnet, pbkrtest (>= 0.3-2), quantreg, grDevices, utils, stats, graphics, Rcpp
library(ggplot2) #version 2.2.1 Date/Publication: 2016-12-30 Depends: R (>= 3.1) Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(gridExtra) #version: 2.2.1 Date/Publication: 2016-02-29 Depends: R(>= 2.5.0) Imports: gtable, grid, grDevices, graphics, utils
library(lsmeans)  #version: 2.26-3 Date: 2017-05-09 Depends: estimability, methods, R (>= 3.0) Imports: graphics, stats, utils, nlme, coda (>= 0.17), multcomp, plyr,mvtnorm, xtable (>= 1.8-2)
library(multcomp) #version: 1.4-6 Date: 2016-07-14 Depends: stats, graphics, mvtnorm (>= 1.0-3), survival (>= 2.39-4), TH.data (>= 1.0-2)
library(nlme) #version: 3.1-131 Date: 2017-02-06 Depends: R (>= 3.0.2) Imports: graphics, stats, utils, lattice
library(plotrix) #version: 3.6-5 Date: 2017-05-09 Depends: NA Imports: grDevices, graphics, stats, utils
library(plyr) #version: 1.8.4 Date/Publication: 2016-06-08 Depends: R (>= 3.1.0) Imports: Rcpp (>= 0.11.0)
library(reshape) #version: 3.3.1 Date/Publication: 2016-06-24  Depends: R (>= 3.3.1)
library(seacarb) #version: 3.2 Date/Publication: 2017-06-19 Depends: R (>= 2.10), oce, gsw Imports: NA
library(grid) #version: 3.3.1 Date/Publication: 2016-06-24  Depends: R (>= 3.3.1)
library(xtable) #version 1.8-2 Date/Publication: 2016-01-08 Depends: R (>= 2.10.0)
library(lme4) #version: 1.1-13 Date/Publication: 2017-04-19 Depends: R (>= 3.0.2), Matrix (>= 1.1.1), methods, stats Imports: graphics, grid, splines, utils, parallel, MASS, lattice, nlme(>= 3.1-123), minqa (>= 1.1.15), nloptr (>= 1.0.4)
library(dplyr)
#library(blmeco) #version: 1.1 Date/Publication: 2015-08-22 Depends: R (>= 3.0.0), stats, MASS Imports: MuMIn, arm, lme4
#library(MuMIn) #version: 1.15.6 Date/Publication: 2016-01-07 Depends: R (>= 3.0.0) Imports: graphics, methods, Matrix, stats, stats4

#####Required Data files#####

#CRM_TA_Data.csv
#Daily_Temp_pH_Sal.csv
#~/MyProjects/Porites_Nutrition_June2019/data/pH_Calibration_Files

##### DISCRETE pH CALCULATIONS #####
path <-("~/MyProjects/Porites_Rim_Bleaching_2019/data/Daily_Measurements/pH_Calibration_Files")
file.names<-list.files(path = path, pattern = "csv$") #list all the file names in the folder to get only get the csv files
pH.cals <- data.frame(matrix(NA, nrow=length(file.names), ncol=3, dimnames=list(file.names,c("Date", "Intercept", "Slope")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Calib.Data <-read.table(file.path(path,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
  file.names[i]
  model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
  coe <- coef(model) #extracts the coeffecients
  summary(model)$r.squared
  plot(Calib.Data$mVTris, Calib.Data$TTris)
  pH.cals[i,2:3] <- coe #inserts them in the dataframe
  pH.cals[i,1] <- substr(file.names[i],1,8) #stores the file name in the Date column
}
colnames(pH.cals) <- c("Calib.Date",  "Intercept",  "Slope") #rename columns
pH.cals #view data

#constants for use in pH calculation 
R <- 8.31447215 #gas constant in J mol-1 K-1 
F <-96485.339924 #Faraday constant in coulombs mol-1

#read in probe measurements of pH, temperature, and salinity from tanks
daily <- read.csv("data/Daily_Measurements/Daily_Temp_pH_Sal_RIM_ADULT.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
daily <- daily
min(na.omit(daily$Temperature))
max(na.omit(daily$Temperature))

# Making tank number categorical
daily$Tank <- as.factor(daily$Tank)

#plot(na.omit(daily$Sample.ID[1:84]), na.omit(daily$Temperature[1:84]), las=2)

#merge with Seawater chemistry file
SW.chem <- merge(pH.cals, daily, by="Calib.Date")

SW.chem <- SW.chem[with(SW.chem, order(Date, Time, Tank)), ]

mvTris <- SW.chem$Temperature*SW.chem$Slope+SW.chem$Intercept #calculate the mV of the tris standard using the temperature mv relationships in the measured standard curves 
STris<-35 #salinity of the Tris
phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(SW.chem$Temperature+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(SW.chem$Temperature+273.15)-0.11149858*(SW.chem$Temperature+273.15) #calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem$pH.Total<-phTris+(mvTris/1000-SW.chem$pH.MV/1000)/(R*(SW.chem$Temperature+273.15)*log(10)/F) #calculate the pH on the total scale (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)

#Combining date and time
SW.chem$Date.Time <- paste(SW.chem$Date, SW.chem$Time, sep = "-")
SW.chem$Date.Time <- as.POSIXct(SW.chem$Date.Time, format="%Y%m%d-%H:%M")

##### TESTING #####

# pdf("~/MyProjects/Porites_Nutrition_June2019/output/Daily_Measures_Tank_Time.pdf")
# par(mfrow=c(1,3))
# plot(SW.chem$Date.Time, SW.chem$Temperature, xlab="Date and Time", ylab="Temperature°C", ylim=c(23,29),las=2, col=SW.chem$Tank)
# plot(SW.chem$Date.Time, SW.chem$pH.Total, xlab="Date and Time", ylab="pH Total Scale", ylim=c(7.9,8.1),las=2, col=SW.chem$Tank)
# plot(SW.chem$Date.Time, SW.chem$Salinity, xlab="Date and Time", ylab="Salinity psu", ylim=c(36,38),las=2, col=SW.chem$Tank)
# dev.off()
# 
# pdf("~/MyProjects/Porites_Nutrition_June2019/output/Daily_by_Tank.pdf")
# par(mfrow=c(1,3))
# plot(SW.chem$Tank, SW.chem$Temperature, xlab="Tank", ylab="Temperature°C", ylim=c(23,29),las=2, col=SW.chem$Tank)
# plot(SW.chem$Tank, SW.chem$pH.Total, xlab="Tank", ylab="pH Total Scale", ylim=c(7.4,8.1),las=2, col=SW.chem$Tank)
# plot(SW.chem$Tank, SW.chem$Salinity, xlab="Tank", ylab="Salinity psu", ylim=c(36,38),las=2, col=SW.chem$Tank)
# dev.off()
# 
# 
# pdf("~/MyProjects/Porites_Nutrition_June2019/output/Daily_by_Treatment.pdf")
# par(mfrow=c(1,3))
# plot(SW.chem$Treatment, SW.chem$Temperature, xlab="Treatment", ylab="Temperature°C", ylim=c(23,29),las=2, col=SW.chem$Tank)
# plot(SW.chem$Treatment, SW.chem$pH.Total, xlab="Treatment", ylab="pH Total Scale", ylim=c(7.4,8.1),las=2, col=SW.chem$Tank)
# plot(SW.chem$Treatment, SW.chem$Salinity, xlab="Treatment", ylab="Salinity psu", ylim=c(36,38),las=2, col=SW.chem$Tank)
# dev.off()
# 
# Temp <- ggplot(SW.chem, aes(x=Date.Time, y=Temperature, group = Tank, color = Tank))+
#   geom_point() +
#   geom_line() + 
#   theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
#                     panel.grid.minor = element_blank(), axis.line = element_blank()) + theme(legend.position = "none")
# pH <- ggplot(SW.chem, aes(x=Date.Time, y=pH.Total, group = Tank, color = Tank))+
#   geom_point() +
#   geom_line() + 
#   theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
#                      panel.grid.minor = element_blank(), axis.line = element_blank()) + theme(legend.position = "none")
# Salinity <- ggplot(SW.chem, aes(x=Date.Time, y=Salinity, group = Tank, color = Tank))+
#   geom_point() +
#   geom_line() + 
#   theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
#                      panel.grid.minor = element_blank(), axis.line = element_blank()) 
# sw.chem <- arrangeGrob (Temp, pH, Salinity, ncol=3)
# ggsave(file="output/Daily_Measures_Tank_Time.pdf", sw.chem, width = 11, height = 11, units = c("in"))

### Experimental 

SW.chem.Exp <- SW.chem %>%
  filter(Type == "Experimental")

pdf("~/MyProjects/Porites_Nutrition_June2019/output/Daily_Measures_AdultRIM_Tank_Time_Experimental.pdf")
par(mfrow=c(1,3))
plot(SW.chem.Exp$Date.Time, SW.chem.Exp$Temperature, xlab="Date and Time", ylab="Temperature°C", ylim=c(23,33),las=2, col=SW.chem.Exp$Tank)
plot(SW.chem.Exp$Date.Time, SW.chem.Exp$pH.Total, xlab="Date and Time", ylab="pH Total Scale", ylim=c(7.9,8.1),las=2, col=SW.chem.Exp$Tank)
plot(SW.chem.Exp$Date.Time, SW.chem.Exp$Salinity, xlab="Date and Time", ylab="Salinity psu", ylim=c(36,38),las=2, col=SW.chem.Exp$Tank)
dev.off()

pdf("~/MyProjects/Porites_Nutrition_June2019/output/Daily_by_Tank_AdultRIM_Experimental.pdf")
par(mfrow=c(1,3))
plot(SW.chem.Exp$Tank, SW.chem.Exp$Temperature, xlab="Tank", ylab="Temperature°C", ylim=c(23,33),las=2, col=SW.chem.Exp$Tank)
plot(SW.chem.Exp$Tank, SW.chem.Exp$pH.Total, xlab="Tank", ylab="pH Total Scale", ylim=c(7.8,8.1),las=2, col=SW.chem.Exp$Tank)
plot(SW.chem.Exp$Tank, SW.chem.Exp$Salinity, xlab="Tank", ylab="Salinity psu", ylim=c(36.5,37.5),las=2, col=SW.chem.Exp$Tank)
dev.off()


pdf("~/MyProjects/Porites_Nutrition_June2019/output/Daily_by_Treatment_AdultRIM_Experimental.pdf")
par(mfrow=c(1,3))
plot(SW.chem.Exp$Treatment, SW.chem.Exp$Temperature, xlab="Treatment", ylab="Temperature°C", ylim=c(23,33),las=2, col=SW.chem.Exp$Tank)
plot(SW.chem.Exp$Treatment, SW.chem.Exp$pH.Total, xlab="Treatment", ylab="pH Total Scale", ylim=c(7.8,8.1),las=2, col=SW.chem.Exp$Tank)
plot(SW.chem.Exp$Treatment, SW.chem.Exp$Salinity, xlab="Treatment", ylab="Salinity psu", ylim=c(36.5,37.5),las=2, col=SW.chem.Exp$Tank)
dev.off()

Temp <- ggplot(SW.chem.Exp, aes(x=Date.Time, y=Temperature, group = Tank, color = Tank))+
  geom_point() +
  geom_line() + 
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) + theme(legend.position = "none")
pH <- ggplot(SW.chem.Exp, aes(x=Date.Time, y=pH.Total, group = Tank, color = Tank))+
  geom_point() +
  geom_line() + 
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) + theme(legend.position = "none")
Salinity <- ggplot(SW.chem.Exp, aes(x=Date.Time, y=Salinity, group = Tank, color = Tank))+
  geom_point() +
  geom_line() + 
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) 
sw.chem <- arrangeGrob (Temp, pH, Salinity, ncol=1)
ggsave(file="output/Daily_Measures_Tank_Time_AdultRIM_Experimental.pdf", sw.chem, width = 11, height = 11, units = c("in"))


### parameters 
library(Rmisc)
library(dplyr)

SW.chem.Exp.clean <- SW.chem.Exp %>% filter(pH.Total != "NA")
SW.chem.Exp.clean$Treatment <- as.factor(SW.chem.Exp.clean$Treatment)

pH.sum <- summarySE(SW.chem.Exp.clean, measurevar="pH.Total", groupvars="Treatment")
t.test(pH.Total ~ Treatment, data = SW.chem.Exp.clean)

sal.sum <- summarySE(SW.chem.Exp.clean, measurevar="Salinity", groupvars="Treatment")
t.test(Salinity ~ Treatment, data = SW.chem.Exp.clean)



