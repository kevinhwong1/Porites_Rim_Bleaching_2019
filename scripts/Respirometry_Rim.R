#Title: Photosynthesis and Respiration Calculations for Porites Nurition July 2019 Project (Rim)
#Author: HM Putnam and NJ Silbiger
#Edited by: Kevin Wong
#Date Last Modified: 20190711
#See Readme file for details 

##Install packages
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("stringr" %in% rownames(installed.packages()) == 'FALSE') install.packages('stringr') 

#Read in required libraries

##### Include Versions of libraries
install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('stringr')
library(lsmeans)

##### PHOTOSYNTHESIS Time 0 20190708 #####
path.p<-"data/Colony_Respirometry/20190708" #the location of all your respirometry files 

# bring in the respiration files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file

#add file names that include the subdirectory name
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#generate a 3 column dataframe with specific column names
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=5))
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec","Temp.C","PR")

#Load Sample Info T1
Sample.Info <- read.csv(file="data/Colony_Respirometry/20190708_Photo_Resp_Sample_Info.csv", header=T) #read in volume and sample.info data
#Sample.Info$Vol.L <- Sample.Info$Chamber.Vol.mL/1000 #calculate volume

#PHOTOSYNTHESIS
path.p<-"data/Colony_Respirometry/20190708" #the location of all your respirometry files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=3)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec")
file.names
#Add names for photosynthesis or respiration for for loop
PR<-c('Photo','Resp')

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]), skip = 1, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  brk <- which(diff(Photo.Data1$Time) > 30) #look for breaks in time of 30 seconds or more
  Photo <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) < brk)  #subset by break in time stamp keeping everything before break
  Resp <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) > brk) #subset by break in time stamp keeping everything before break
  lt.levs <- list(Photo, Resp) #list levels of segmentation
  
  
  for(j in 1:length(lt.levs)){    
    Photo.Data <- as.data.frame(lt.levs[j])
    n<-dim(Photo.Data )[1] #identify length of data
    Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data )[1] #list length of trimmed data
    Photo.Data $sec <- 1:n #set seconds by one from start to finish of run
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("output/Colony_Photo_Resp_Output/Rim_T1",rename,"_",j,"thinning.pdf"))
    par(omi=rep(0.3, 4)) #set size of the outer margins in inches
    par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    
    #save original unthinned data
    Photo.Data.orig <- Photo.Data
    Photo.Data   <-  thinData(Photo.Data , by=20)$newData1 #thin data by every 20 points for all the O2 values
    Photo.Data$sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
    
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr') #plotting graphics using 'usr'
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) #giving specific coordinates for plot using 'usr'
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    dev.off()
    ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
    Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2, 
                         method="pc", verbose=TRUE) 
    pdf(paste0("output/Colony_Photo_Resp_Output/Rim_T1",rename,"_",j,"regression.pdf"))
    plot(Regs)
    dev.off()
    
    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
    #Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
    Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
    Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
  }
}

write.csv(Photo.R, 'data/Colony_Respirometry/Rim_T1_Photo.R.csv')

Photo.R <- read.csv('data/Colony_Respirometry/Rim_T1_Photo.R.csv')#Split up the photostynthesis and respiration data into two dataframes
PHO <- Photo.R[Photo.R$V5=='Photo', ]
RES <- Photo.R[Photo.R$V5=='Resp', ]

#Load Sample Info
Sample.Info <- read.csv(file="data/Colony_Respirometry/20190708_Photo_Resp_Sample_Info.csv", header=T) #read in volume and sample.info data
#Convert sample volume to mL
Sample.Info$Water.Height.in <- Sample.Info$Water.Height.cm*0.393701
Sample.Info$Chamber.Vol.L <- (Sample.Info$Water.Height.in - 0.3125)*15.625*8.000*0.0163871 

#Merge data with sample info
Resp <- merge(RES,Sample.Info, by="Fragment.ID" )
Photo <- merge(PHO,Sample.Info, by="Fragment.ID")

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Resp$umol.sec <- Resp$umol.L.sec*Sample.Info$Chamber.Vol.L
Photo$umol.sec <- Photo$umol.L.sec*Sample.Info$Chamber.Vol.L

## BLANKS HAVE TO BE SPECIFIC TO RESPONSE VARIABLE (I.E., PHOTO OR RESP) AND TEMP
photo.blnk <- aggregate(umol.sec ~ Sample.Type, data=Photo, mean)
photo.Blanks <- subset(photo.blnk, Sample.Type == "Blank")
Photo$Blank <- photo.Blanks[1,2]

#Calculate resp blank rate
resp.blnk <- aggregate(umol.sec ~ Sample.Type, data=Resp, mean)
resp.Blanks <- subset(resp.blnk, Sample.Type == "Blank")
Resp$Blank <- resp.Blanks[1,2]

#Account for blank rate Subtract Blank by the temperature blank
Resp$umol.sec.corr <- Resp$umol.sec-Resp$Blank
Photo$umol.sec.corr <- Photo$umol.sec-Photo$Blank

#normalize to surface area and h-1

#Merging SA file
SA <- read.csv("data/Surface_Area/colony.Rim.calcsa.csv")
Resp2 <- merge(Resp, SA, by = "Fragment.ID")
Photo2 <- merge(Photo, SA, by = "Fragment.ID")

Resp2$umol.cm2.hr <- (Resp2$umol.sec.corr*3600)/Resp2$SA.cm2
Photo2$umol.cm2.hr <- (Photo2$umol.sec.corr*3600)/Photo2$SA.cm2

##T1 Results
#remove blanks from dataset
Photo3 <- subset(Photo2, Sample.Type!= "Blank")
Resp3 <- subset(Resp2, Sample.Type!= "Blank")

write.csv(Photo3, file="output/Colony_Photo_Resp_Output/Photosynthesis.rates.T1.csv")
write.csv(Resp3, file="output/Colony_Photo_Resp_Output/Respiration.rates.T1.csv")

#calculate gross photosynthesis Pnet -- Rdark

Photo3 <- read.csv("output/Colony_Photo_Resp_Output/Photosynthesis.rates.T1.csv")
Resp3 <- read.csv("output/Colony_Photo_Resp_Output/Respiration.rates.T1.csv")

resp.data.T1 <- merge(Photo3[,c(2,18,35)],Resp3[,c(2,35)], by="Fragment.ID")

#rename the columns
names(resp.data.T1)[names(resp.data.T1) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
names(resp.data.T1)[names(resp.data.T1) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"

#Pnet plus resp (if positive) is pGross
resp.data.T1$Pgross_umol.cm2.hr <- resp.data.T1$Pnet_umol.cm2.hr-resp.data.T1$Rdark_umol.cm2.hr
resp.data.T1$Timepoint <- "T1"

write.csv(resp.data.T1, file="output/Colony_Photo_Resp_Output/Rim_Respdata.T1.csv")

# 
# #T2
# #PHOTOSYNTHESIS
# path.p<-"data/July_Rim/Colony_Respirometry/20190805" #the location of all your respirometry files
# file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
# Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=3)) #generate a 3 column dataframe with specific column names
# colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec")
# file.names
# #Add names for photosynthesis or respiration for for loop
# PR<-c('Photo','Resp')
# 
# for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
#   Photo.Data1 <-read.table(file.path(path.p,file.names[i]), skip = 1, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
#   Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
#   Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
#   brk <- which(diff(Photo.Data1$Time) > 30) #look for breaks in time of 30 seconds or more
#   Photo <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) < brk)  #subset by break in time stamp keeping everything before break
#   Resp <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) > brk) #subset by break in time stamp keeping everything before break
#   lt.levs <- list(Photo, Resp) #list levels of segmentation
#   
#   
#   for(j in 1:length(lt.levs)){    
#     Photo.Data <- as.data.frame(lt.levs[j])
#     n<-dim(Photo.Data )[1] #identify length of data
#     Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
#     n<-dim(Photo.Data )[1] #list length of trimmed data
#     Photo.Data $sec <- 1:n #set seconds by one from start to finish of run
#     
#     #Save plot prior to and after data thinning to make sure thinning is not too extreme
#     rename <- sub("_.*", "", file.names[i])
#     pdf(paste0("output/July_Rim/Colony_Photo_Resp_Output/Rim_T2",rename,"_",j,"thinning.pdf"))
#     par(omi=rep(0.3, 4)) #set size of the outer margins in inches
#     par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
#     plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
#     usr  <-  par('usr')
#     rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
#     whiteGrid()
#     box()
#     points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
#     axis(1)
#     axis(2, las=1)
#     
#     #save original unthinned data
#     Photo.Data.orig <- Photo.Data
#     Photo.Data   <-  thinData(Photo.Data , by=20)$newData1 #thin data by every 20 points for all the O2 values
#     Photo.Data$sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
#     
#     plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
#     usr  <-  par('usr') #plotting graphics using 'usr'
#     rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) #giving specific coordinates for plot using 'usr'
#     whiteGrid()
#     box()
#     points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
#     axis(1)
#     axis(2, las=1)
#     dev.off()
#     ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
#     Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2, 
#                          method="pc", verbose=TRUE) 
#     pdf(paste0("output/July_Rim/Colony_Photo_Resp_Output/Rim_T2",rename,"_",j,"regression.pdf"))
#     plot(Regs)
#     dev.off()
#     
#     s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
#     Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
#     Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
#     #Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
#     Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
#     Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
#   }
# }
# 
# write.csv(Photo.R, 'data/July_Rim/Colony_Respirometry/Rim_T2_Photo.R.csv')
# 
# Photo.R <- read.csv('data/July_Rim/Colony_Respirometry/Rim_T2_Photo.R.csv')#Split up the photostynthesis and respiration data into two dataframes
# PHO <- Photo.R[Photo.R$V5=='Photo', ]
# RES <- Photo.R[Photo.R$V5=='Resp', ]
# 
# 
# Sample.Info2 <- read.csv(file="data/July_Rim/Colony_Respirometry/20190805_Photo_Resp_Sample_Info.csv", header=T) #read in volume and sample.info data
# 
# #Convert sample volume to mL
# Sample.Info2$Water.Height.in <- Sample.Info2$Water.Height.cm*0.393701
# Sample.Info2$Chamber.Vol.L <- (Sample.Info2$Water.Height.in - 0.3125)*15.625*8.000*0.0163871 
# 
# #Merge with sample info
# Resp2 <- merge(RES,Sample.Info2, by="Fragment.ID")
# Photo2 <- merge(PHO,Sample.Info2, by="Fragment.ID")
# 
# #split sample info between photo and resp
# #Sample.Info.P2<-Sample.Info2[Sample.Info2$Light_Dark=='Light',]
# #Sample.Info.R2<-Sample.Info2[Sample.Info2$Light_Dark=='Dark',]
# 
# #Account for chamber volume umol s-1
# 
# Resp2$umol.sec <- Resp2$umol.L.sec*Resp2$Chamber.Vol.L
# Photo2$umol.sec <- Photo2$umol.L.sec*Photo2$Chamber.Vol.L
# 
# ##Getting Blank means and attaching to data
# photo.blnk <- aggregate(umol.sec ~ Sample.Type*Temp.Cat, data=Photo2, mean)
# photo.Blanks <- subset(photo.blnk, Sample.Type == "Blank") #values for blank
# Photo2$Blank <- photo.Blanks$umol.sec[match(Photo2$Temp.Cat,photo.Blanks$Temp.Cat)] #matching the right treatment blank to samples
# 
# #Calculate resp blank rate
# resp.blnk <- aggregate(umol.sec ~ Sample.Type*Temp.Cat, data=Resp2, mean)
# resp.Blanks <- subset(resp.blnk, Sample.Type == "Blank")
# Resp2$Blank <- resp.Blanks$umol.sec[match(Resp2$Temp.Cat,resp.Blanks$Temp.Cat)]
# 
# #Account for blank rate Subtract Blank by the temperature blank
# Resp2$umol.sec.corr <- Resp2$umol.sec-Resp2$Blank
# Photo2$umol.sec.corr <- Photo2$umol.sec-Photo2$Blank
# 
# #Merging SA file
# SA <- read.csv("data/July_Rim/Surface_Area/colony.Rim.calcsa.csv")
# Resp3 <- merge(Resp2, SA, by = "Fragment.ID")
# Photo3 <- merge(Photo2, SA, by = "Fragment.ID")
# 
# #normalize to surface area and h-1
# Resp3$umol.cm2.hr <- (Resp3$umol.sec.corr*3600)/Resp3$SA.cm2
# Photo3$umol.cm2.hr <- (Photo3$umol.sec.corr*3600)/Photo3$SA.cm2
# 
# write.csv(Photo2, file="output/July_Rim/Colony_Photo_Resp_Output/Photosynthesis.rates.T2.csv")
# write.csv(Resp2, file="output/July_Rim/Colony_Photo_Resp_Output/Respiration.rates.T2.csv")
# 
# #T2 Results
# Photo2 <- read.csv("output/July_Rim/Colony_Photo_Resp_Output/Photosynthesis.rates.T2.csv")
# Resp2 <- read.csv("output/July_Rim/Colony_Photo_Resp_Output/Respiration.rates.T2.csv")
# #remove blanks
# Photo3 <- subset(Photo3, Sample.Type!= "Blank")
# Resp3 <- subset(Resp3, Sample.Type!= "Blank")
# 
# #calculate gross photosynthesis Pnet -- Rdark
# 
# resp.data.T2 <- merge(Photo3[,c(1,17,34)],Resp3[,c(1,34)], by="Fragment.ID")
# 
# #rename the column
# names(resp.data.T2)[names(resp.data.T2) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
# names(resp.data.T2)[names(resp.data.T2) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"
# names(resp.data.T2)[names(resp.data.T2) == "Treatment.y"]<- "Treatment" 
# 
# #Pnet plus resp (if positive) is pGross
# resp.data.T2$Pgross_umol.cm2.hr <- resp.data.T2$Pnet_umol.cm2.hr-resp.data.T2$Rdark_umol.cm2.hr
# resp.data.T2$Timepoint <- "T2"
# 
# write.csv(resp.data.T2, file="output/July_Rim/Colony_Photo_Resp_Output/Rim_Respdata.T2.csv")


#T2 (20190815)
#PHOTOSYNTHESIS
path.p<-"data/Colony_Respirometry/20190815" #the location of all your respirometry files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=3)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec")
file.names
#Add names for photosynthesis or respiration for for loop
PR<-c('Photo','Resp')

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]), skip = 1, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  brk <- which(diff(Photo.Data1$Time) > 30) #look for breaks in time of 30 seconds or more
  Photo <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) < brk)  #subset by break in time stamp keeping everything before break
  Resp <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) > brk) #subset by break in time stamp keeping everything before break
  lt.levs <- list(Photo, Resp) #list levels of segmentation
  
  
  for(j in 1:length(lt.levs)){    
    Photo.Data <- as.data.frame(lt.levs[j])
    n<-dim(Photo.Data )[1] #identify length of data
    Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data )[1] #list length of trimmed data
    Photo.Data $sec <- 1:n #set seconds by one from start to finish of run
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("output/Colony_Photo_Resp_Output/Rim_T2",rename,"_",j,"thinning.pdf"))
    par(omi=rep(0.3, 4)) #set size of the outer margins in inches
    par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    
    #save original unthinned data
    Photo.Data.orig <- Photo.Data
    Photo.Data   <-  thinData(Photo.Data , by=20)$newData1 #thin data by every 20 points for all the O2 values
    Photo.Data$sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
    
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr') #plotting graphics using 'usr'
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) #giving specific coordinates for plot using 'usr'
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    dev.off()
    ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
    Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2, 
                         method="pc", verbose=TRUE) 
    pdf(paste0("output/Colony_Photo_Resp_Output/Rim_T2",rename,"_",j,"regression.pdf"))
    plot(Regs)
    dev.off()
    
    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
    #Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
    Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
    Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
  }
}

write.csv(Photo.R, 'data/Colony_Respirometry/Rim_T2_Photo.R.csv')

Photo.R <- read.csv('data/Colony_Respirometry/Rim_T2_Photo.R.csv')#Split up the photostynthesis and respiration data into two dataframes
PHO <- Photo.R[Photo.R$V5=='Photo', ]
RES <- Photo.R[Photo.R$V5=='Resp', ]


Sample.Info3 <- read.csv(file="data/Colony_Respirometry/20190815_Photo_Resp_Sample_Info.csv", header=T) #read in volume and sample.info data

#Convert sample volume to mL
Sample.Info3$Water.Height.in <- Sample.Info3$Water.Height.cm*0.393701
Sample.Info3$Chamber.Vol.L <- (Sample.Info3$Water.Height.in - 0.3125)*15.625*8.000*0.0163871 

#Merge with sample info
Resp2 <- merge(RES,Sample.Info3, by="Fragment.ID")
Photo2 <- merge(PHO,Sample.Info3, by="Fragment.ID")

#split sample info between photo and resp
#Sample.Info.P2<-Sample.Info2[Sample.Info2$Light_Dark=='Light',]
#Sample.Info.R2<-Sample.Info2[Sample.Info2$Light_Dark=='Dark',]

#Account for chamber volume umol s-1

Resp2$umol.sec <- Resp2$umol.L.sec*Resp2$Chamber.Vol.L
Photo2$umol.sec <- Photo2$umol.L.sec*Photo2$Chamber.Vol.L

##Getting Blank means and attaching to data
photo.blnk <- aggregate(umol.sec ~ Sample.Type*Temp.Cat, data=Photo2, mean)
photo.Blanks <- subset(photo.blnk, Sample.Type == "Blank") #values for blank
Photo2$Blank <- photo.Blanks$umol.sec[match(Photo2$Temp.Cat,photo.Blanks$Temp.Cat)] #matching the right treatment blank to samples

#Calculate resp blank rate
resp.blnk <- aggregate(umol.sec ~ Sample.Type*Temp.Cat, data=Resp2, mean)
resp.Blanks <- subset(resp.blnk, Sample.Type == "Blank")
Resp2$Blank <- resp.Blanks$umol.sec[match(Resp2$Temp.Cat,resp.Blanks$Temp.Cat)]

#Account for blank rate Subtract Blank by the temperature blank
Resp2$umol.sec.corr <- Resp2$umol.sec-Resp2$Blank
Photo2$umol.sec.corr <- Photo2$umol.sec-Photo2$Blank

#Merging SA file
SA <- read.csv("data/Surface_Area/colony.Rim.calcsa.csv")
Resp3 <- merge(Resp2, SA, by = "Fragment.ID")
Photo3 <- merge(Photo2, SA, by = "Fragment.ID")

#normalize to surface area and h-1
Resp3$umol.cm2.hr <- (Resp3$umol.sec.corr*3600)/Resp3$SA.cm2
Photo3$umol.cm2.hr <- (Photo3$umol.sec.corr*3600)/Photo3$SA.cm2

write.csv(Photo3, file="output/Colony_Photo_Resp_Output/Photosynthesis.rates.T2.csv")
write.csv(Resp3, file="output/Colony_Photo_Resp_Output/Respiration.rates.T2.csv")

#T3 Results
Photo3 <- read.csv("output/Colony_Photo_Resp_Output/Photosynthesis.rates.T2.csv")
Resp3 <- read.csv("output/Colony_Photo_Resp_Output/Respiration.rates.T2.csv")
#remove blanks
Photo3 <- subset(Photo3, Sample.Type!= "Blank")
Resp3 <- subset(Resp3, Sample.Type!= "Blank")

#calculate gross photosynthesis Pnet -- Rdark

resp.data.T3 <- merge(Photo3[,c(2,18,35)],Resp3[,c(2,35)], by="Fragment.ID")

#rename the column
names(resp.data.T3)[names(resp.data.T3) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
names(resp.data.T3)[names(resp.data.T3) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"
names(resp.data.T3)[names(resp.data.T3) == "Treatment.y"]<- "Treatment" 

#Pnet plus resp (if positive) is pGross
resp.data.T3$Pgross_umol.cm2.hr <- resp.data.T3$Pnet_umol.cm2.hr-resp.data.T3$Rdark_umol.cm2.hr
resp.data.T3$Timepoint <- "T2"

write.csv(resp.data.T3, file="output/Colony_Photo_Resp_Output/Rim_Respdata.T2.csv")

# 
# #T4
# #PHOTOSYNTHESIS
# path.p<-"data/July_Rim/Colony_Respirometry/20190823" #the location of all your respirometry files
# file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
# Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=3)) #generate a 3 column dataframe with specific column names
# colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec")
# file.names
# #Add names for photosynthesis or respiration for for loop
# PR<-c('Photo','Resp')
# 
# for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
#   Photo.Data1 <-read.table(file.path(path.p,file.names[i]), skip = 1, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
#   Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
#   Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
#   brk <- which(diff(Photo.Data1$Time) > 30) #look for breaks in time of 30 seconds or more
#   Photo <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) < brk)  #subset by break in time stamp keeping everything before break
#   Resp <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) > brk) #subset by break in time stamp keeping everything before break
#   lt.levs <- list(Photo, Resp) #list levels of segmentation
#   
#   
#   for(j in 1:length(lt.levs)){    
#     Photo.Data <- as.data.frame(lt.levs[j])
#     n<-dim(Photo.Data )[1] #identify length of data
#     Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
#     n<-dim(Photo.Data )[1] #list length of trimmed data
#     Photo.Data $sec <- 1:n #set seconds by one from start to finish of run
#     
#     #Save plot prior to and after data thinning to make sure thinning is not too extreme
#     rename <- sub("_.*", "", file.names[i])
#     pdf(paste0("output/July_Rim/Colony_Photo_Resp_Output/Rim_T4",rename,"_",j,"thinning.pdf"))
#     par(omi=rep(0.3, 4)) #set size of the outer margins in inches
#     par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
#     plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
#     usr  <-  par('usr')
#     rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
#     whiteGrid()
#     box()
#     points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
#     axis(1)
#     axis(2, las=1)
#     
#     #save original unthinned data
#     Photo.Data.orig <- Photo.Data
#     Photo.Data   <-  thinData(Photo.Data , by=20)$newData1 #thin data by every 20 points for all the O2 values
#     Photo.Data$sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
#     
#     plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
#     usr  <-  par('usr') #plotting graphics using 'usr'
#     rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) #giving specific coordinates for plot using 'usr'
#     whiteGrid()
#     box()
#     points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
#     axis(1)
#     axis(2, las=1)
#     dev.off()
#     ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
#     Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2, 
#                          method="pc", verbose=TRUE) 
#     pdf(paste0("output/July_Rim/Colony_Photo_Resp_Output/Rim_T4",rename,"_",j,"regression.pdf"))
#     plot(Regs)
#     dev.off()
#     
#     s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
#     Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
#     Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
#     #Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
#     Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
#     Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
#   }
# }
# 
# write.csv(Photo.R, 'data/July_Rim/Colony_Respirometry/Rim_T4_Photo.R.csv')
# 
# Photo.R <- read.csv('data/July_Rim/Colony_Respirometry/Rim_T4_Photo.R.csv')#Split up the photostynthesis and respiration data into two dataframes
# PHO <- Photo.R[Photo.R$V5=='Photo', ]
# RES <- Photo.R[Photo.R$V5=='Resp', ]
# 
# 
# Sample.Info3 <- read.csv(file="data/July_Rim/Colony_Respirometry/20190823_Photo_Resp_Sample_Info.csv", header=T) #read in volume and sample.info data
# 
# #Convert sample volume to mL
# Sample.Info3$Water.Height.in <- Sample.Info3$Water.Height.cm*0.393701
# Sample.Info3$Chamber.Vol.L <- (Sample.Info3$Water.Height.in - 0.3125)*15.625*8.000*0.0163871 
# 
# #Merge with sample info
# Resp2 <- merge(RES,Sample.Info3, by="Fragment.ID")
# Photo2 <- merge(PHO,Sample.Info3, by="Fragment.ID")
# 
# #split sample info between photo and resp
# #Sample.Info.P2<-Sample.Info2[Sample.Info2$Light_Dark=='Light',]
# #Sample.Info.R2<-Sample.Info2[Sample.Info2$Light_Dark=='Dark',]
# 
# #Account for chamber volume umol s-1
# 
# Resp2$umol.sec <- Resp2$umol.L.sec*Resp2$Chamber.Vol.L
# Photo2$umol.sec <- Photo2$umol.L.sec*Photo2$Chamber.Vol.L
# 
# ##Getting Blank means and attaching to data
# photo.blnk <- aggregate(umol.sec ~ Sample.Type*Temp.Cat, data=Photo2, mean)
# photo.Blanks <- subset(photo.blnk, Sample.Type == "Blank") #values for blank
# Photo2$Blank <- photo.Blanks$umol.sec[match(Photo2$Temp.Cat,photo.Blanks$Temp.Cat)] #matching the right treatment blank to samples
# 
# #Calculate resp blank rate
# resp.blnk <- aggregate(umol.sec ~ Sample.Type*Temp.Cat, data=Resp2, mean)
# resp.Blanks <- subset(resp.blnk, Sample.Type == "Blank")
# Resp2$Blank <- resp.Blanks$umol.sec[match(Resp2$Temp.Cat,resp.Blanks$Temp.Cat)]
# 
# #Account for blank rate Subtract Blank by the temperature blank
# Resp2$umol.sec.corr <- Resp2$umol.sec-Resp2$Blank
# Photo2$umol.sec.corr <- Photo2$umol.sec-Photo2$Blank
# 
# #Merging SA file
# SA <- read.csv("data/July_Rim/Surface_Area/colony.Rim.calcsa.csv")
# Resp3 <- merge(Resp2, SA, by = "Fragment.ID")
# Photo3 <- merge(Photo2, SA, by = "Fragment.ID")
# 
# #normalize to surface area and h-1
# Resp3$umol.cm2.hr <- (Resp3$umol.sec.corr*3600)/Resp3$SA.cm2
# Photo3$umol.cm2.hr <- (Photo3$umol.sec.corr*3600)/Photo3$SA.cm2
# 
# write.csv(Photo3, file="output/July_Rim/Colony_Photo_Resp_Output/Photosynthesis.rates.T4.csv")
# write.csv(Resp3, file="output/July_Rim/Colony_Photo_Resp_Output/Respiration.rates.T4.csv")
# 
# #T4 Results
# Photo3 <- read.csv("output/July_Rim/Colony_Photo_Resp_Output/Photosynthesis.rates.T4.csv")
# Resp3 <- read.csv("output/July_Rim/Colony_Photo_Resp_Output/Respiration.rates.T4.csv")
# #remove blanks
# Photo3 <- subset(Photo3, Sample.Type!= "Blank")
# Resp3 <- subset(Resp3, Sample.Type!= "Blank")
# 
# #calculate gross photosynthesis Pnet -- Rdark
# 
# resp.data.T4 <- merge(Photo3[,c(2,18,35)],Resp3[,c(2,35)], by="Fragment.ID")
# 
# #rename the column
# names(resp.data.T4)[names(resp.data.T4) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
# names(resp.data.T4)[names(resp.data.T4) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"
# names(resp.data.T4)[names(resp.data.T4) == "Treatment.y"]<- "Treatment" 
# 
# #Pnet plus resp (if positive) is pGross
# resp.data.T4$Pgross_umol.cm2.hr <- resp.data.T4$Pnet_umol.cm2.hr-resp.data.T4$Rdark_umol.cm2.hr
# resp.data.T4$Timepoint <- "T4"
# 
# write.csv(resp.data.T4, file="output/July_Rim/Colony_Photo_Resp_Output/Rim_Respdata.T4.csv")


#T5 (20190830)
#PHOTOSYNTHESIS
path.p<-"data/Colony_Respirometry/20190830" #the location of all your respirometry files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=3)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec")
file.names
#Add names for photosynthesis or respiration for for loop
PR<-c('Photo','Resp')

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]), skip = 1, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  brk <- which(diff(Photo.Data1$Time) > 30) #look for breaks in time of 30 seconds or more
  Photo <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) < brk)  #subset by break in time stamp keeping everything before break
  Resp <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) > brk) #subset by break in time stamp keeping everything before break
  lt.levs <- list(Photo, Resp) #list levels of segmentation
  
  
  for(j in 1:length(lt.levs)){    
    Photo.Data <- as.data.frame(lt.levs[j])
    n<-dim(Photo.Data )[1] #identify length of data
    Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data )[1] #list length of trimmed data
    Photo.Data $sec <- 1:n #set seconds by one from start to finish of run
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("output/Colony_Photo_Resp_Output/Rim_T3",rename,"_",j,"thinning.pdf"))
    par(omi=rep(0.3, 4)) #set size of the outer margins in inches
    par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    
    #save original unthinned data
    Photo.Data.orig <- Photo.Data
    Photo.Data   <-  thinData(Photo.Data , by=20)$newData1 #thin data by every 20 points for all the O2 values
    Photo.Data$sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
    
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr') #plotting graphics using 'usr'
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) #giving specific coordinates for plot using 'usr'
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    dev.off()
    ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
    Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2, 
                         method="pc", verbose=TRUE) 
    pdf(paste0("output/Colony_Photo_Resp_Output/Rim_T3",rename,"_",j,"regression.pdf"))
    plot(Regs)
    dev.off()
    
    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
    #Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
    Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
    Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
  }
}

write.csv(Photo.R, 'data/Colony_Respirometry/Rim_T3_Photo.R.csv')

Photo.R <- read.csv('data/Colony_Respirometry/Rim_T3_Photo.R.csv')#Split up the photostynthesis and respiration data into two dataframes
PHO <- Photo.R[Photo.R$V5=='Photo', ]
RES <- Photo.R[Photo.R$V5=='Resp', ]


Sample.Info3 <- read.csv(file="data/Colony_Respirometry/20190830_Photo_Resp_Sample_Info.csv", header=T) #read in volume and sample.info data

#Convert sample volume to mL
Sample.Info3$Water.Height.in <- Sample.Info3$Water.Height.cm*0.393701
Sample.Info3$Chamber.Vol.L <- (Sample.Info3$Water.Height.in - 0.3125)*15.625*8.000*0.0163871 

#Merge with sample info
Resp2 <- merge(RES,Sample.Info3, by="Fragment.ID")
Photo2 <- merge(PHO,Sample.Info3, by="Fragment.ID")

#split sample info between photo and resp
#Sample.Info.P2<-Sample.Info2[Sample.Info2$Light_Dark=='Light',]
#Sample.Info.R2<-Sample.Info2[Sample.Info2$Light_Dark=='Dark',]

#Account for chamber volume umol s-1

Resp2$umol.sec <- Resp2$umol.L.sec*Resp2$Chamber.Vol.L
Photo2$umol.sec <- Photo2$umol.L.sec*Photo2$Chamber.Vol.L

##Getting Blank means and attaching to data
photo.blnk <- aggregate(umol.sec ~ Sample.Type*Temp.Cat, data=Photo2, mean)
photo.Blanks <- subset(photo.blnk, Sample.Type == "Blank") #values for blank
Photo2$Blank <- photo.Blanks$umol.sec[match(Photo2$Temp.Cat,photo.Blanks$Temp.Cat)] #matching the right treatment blank to samples

#Calculate resp blank rate
resp.blnk <- aggregate(umol.sec ~ Sample.Type*Temp.Cat, data=Resp2, mean)
resp.Blanks <- subset(resp.blnk, Sample.Type == "Blank")
Resp2$Blank <- resp.Blanks$umol.sec[match(Resp2$Temp.Cat,resp.Blanks$Temp.Cat)]

#Account for blank rate Subtract Blank by the temperature blank
Resp2$umol.sec.corr <- Resp2$umol.sec-Resp2$Blank
Photo2$umol.sec.corr <- Photo2$umol.sec-Photo2$Blank

#Merging SA file
SA <- read.csv("data/July_Rim/Surface_Area/colony.Rim.calcsa.csv")
Resp3 <- merge(Resp2, SA, by = "Fragment.ID")
Photo3 <- merge(Photo2, SA, by = "Fragment.ID")

#normalize to surface area and h-1
Resp3$umol.cm2.hr <- (Resp3$umol.sec.corr*3600)/Resp3$SA.cm2
Photo3$umol.cm2.hr <- (Photo3$umol.sec.corr*3600)/Photo3$SA.cm2

write.csv(Photo3, file="output/Colony_Photo_Resp_Output/Photosynthesis.rates.T3.csv")
write.csv(Resp3, file="output/Colony_Photo_Resp_Output/Respiration.rates.T3.csv")

#T5 Results
Photo3 <- read.csv("output/Colony_Photo_Resp_Output/Photosynthesis.rates.T3.csv")
Resp3 <- read.csv("output/Colony_Photo_Resp_Output/Respiration.rates.T3.csv")
#remove blanks
Photo3 <- subset(Photo3, Sample.Type!= "Blank")
Resp3 <- subset(Resp3, Sample.Type!= "Blank")

#calculate gross photosynthesis Pnet -- Rdark

resp.data.T5 <- merge(Photo3[,c(2,18,35)],Resp3[,c(2,35)], by="Fragment.ID")

#rename the column
names(resp.data.T5)[names(resp.data.T5) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
names(resp.data.T5)[names(resp.data.T5) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"
names(resp.data.T5)[names(resp.data.T5) == "Treatment.y"]<- "Treatment" 

#Pnet plus resp (if positive) is pGross
resp.data.T5$Pgross_umol.cm2.hr <- resp.data.T5$Pnet_umol.cm2.hr-resp.data.T5$Rdark_umol.cm2.hr
resp.data.T5$Timepoint <- "T3"

write.csv(resp.data.T5, file="output/Colony_Photo_Resp_Output/Rim_Respdata.T3.csv")


### Combining datasets

resp.data.T1 <- read.csv("output/Colony_Photo_Resp_Output/Rim_Respdata.T1.csv")
resp.data.T2 <- read.csv("output/Colony_Photo_Resp_Output/Rim_Respdata.T2.csv")
resp.data.T3 <- read.csv("output/Colony_Photo_Resp_Output/Rim_Respdata.T3.csv")
sample.info.col <- read.csv("data/Metadata/Colony_Sample_Info.csv")


resp.data.all <- rbind(resp.data.T1,resp.data.T2, resp.data.T3)
resp.data.all$PR <- resp.data.all$Pgross_umol.cm2.hr / -(resp.data.all$Rdark_umol.cm2.hr)
resp.data.all.info <- merge(resp.data.all, sample.info.col, by = "Fragment.ID")
resp.data.select <- resp.data.all.info %>%
  filter(Group != "NA")

resp.data.select$Group <- factor(resp.data.select$Group, levels = c('Control', 'Bleached', 'Mortality')) #re-ordering factos
resp.data.select$abs.Rdark_umol.cm2.hr <- -(resp.data.select$Rdark_umol.cm2.hr)
#Calculate means
AllMeans <- ddply(resp.data.select, c('Group','Timepoint'), summarize,
                  N = sum(!is.na(Pnet_umol.cm2.hr)), # sample size
                  #pnet
                  Pnet.mean= mean(Pnet_umol.cm2.hr, na.rm=T), #mean pnet
                  Pnet.se = sd(Pnet_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                  #Rdark
                  Rdark.mean= mean(Rdark_umol.cm2.hr, na.rm=T), #mean rdark
                  Rdark.se = sd(Rdark_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                  #Pgross
                  Pgross.mean  = mean(Pgross_umol.cm2.hr, na.rm=T),
                  Pgross.se = sd(Pgross_umol.cm2.hr, na.rm=T)/sqrt(N),
                  #P:R
                  PR.mean  = mean(PR, na.rm=T),
                  PR.se = sd(PR, na.rm=T)/sqrt(N))
#write Results
write.csv(resp.data.select, file="data/Colony_Respirometry/SelectRates.Rim.csv") # raw data
write.csv(AllMeans, file="data/Colony_Respirometry/AllMeans.Rim.csv") # Mean data


PG_Box <- ggplot(resp.data.select, aes(x=Timepoint, y=Pgross_umol.cm2.hr, fill = Group)) +
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5)) +
  stat_summary(fun=mean, geom="line", aes(group=Group, color = Group), position = position_dodge(width = 0.5))  + 
#  stat_summary(fun=mean, geom="point")
  geom_point(pch = 21, position=position_jitterdodge(dodge.width=0.5), outlier.shape= NA) +
#  ylim(0,0.5) +
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(labels=c("T1" = "Day 0", "T2" = "Day 37", "T3" = "Day 52")) +
  xlab("Timepoint") +  ylab(expression("Gross Photosynthetic Rate " (mu*mol ~ cm^{-2} ~ h^{-1}))) + #label y axis + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, color = "black")) +
  theme(legend.position = "none")

Resp_Box <- ggplot(resp.data.select, aes(x=Timepoint, y=abs.Rdark_umol.cm2.hr, fill = Group)) +
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5)) +
  stat_summary(fun=mean, geom="line", aes(group=Group, color = Group), position = position_dodge(width = 0.5))  + 
  #  stat_summary(fun=mean, geom="point")
  geom_point(pch = 21, position=position_jitterdodge(dodge.width=0.5), outlier.shape= NA) +
  #  ylim(0,0.5) +
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(labels=c("T1" = "Day 0", "T2" = "Day 37", "T3" = "Day 52")) +
  xlab("Timepoint") +  ylab(expression("Respiration Rate " (mu*mol ~ cm^{-2} ~ h^{-1}))) + #label y axis + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, color = "black")) +
  theme(legend.position = "none")

PR_Box <- ggplot(resp.data.select, aes(x=Timepoint, y=PR, fill = Group)) +
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5)) +
  stat_summary(fun=mean, geom="line", aes(group=Group, color = Group), position = position_dodge(width = 0.5))  + 
  #  stat_summary(fun=mean, geom="point")
  geom_point(pch = 21, position=position_jitterdodge(dodge.width=0.5), outlier.shape= NA) +
  #  ylim(0,0.5) +
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(labels=c("T1" = "Day 0", "T2" = "Day 37", "T3" = "Day 52")) +
  xlab("Timepoint") +  ylab("P:R") + #label y axis + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, color = "black")) +
  theme(legend.position = c(0.78,0.88), 
        legend.title = element_text(colour="black", size=15),
        legend.text = element_text(colour="black", size = 12))


dev.off()

Resp.Figs <- arrangeGrob (PG_Box, Resp_Box, PR_Box,ncol=3)
ggsave(file="output/Rim_ColonyRespirometry.pdf", Resp.Figs, width = 11, height = 6, units = c("in"))


 ## Statsitics 
#Pgross
PG.anova <- aov(Pgross_umol.cm2.hr~Group*Timepoint, data = resp.data.select)
qqnorm(resid(PG.anova))
qqline(resid(PG.anova))

boxplot(resid(PG.anova)~resp.data.select$Group)
boxplot(resid(PG.anova)~resp.data.select$Timepoint)
anova(PG.anova)

capture.output(anova(PG.anova), file = "output/Statistics/ANOVA_PG.csv")

PG.ph <- lsmeans(PG.anova, ~Group*Timepoint, adjust="tukey")

PG.ph.pairs <- multcomp::cld(PG.ph, alpha=.05, Letters=letters) #list pairwise tests and letter display
capture.output(PG.ph.pairs, file = "output/Statistics/PH_PG.csv")

#Resp
RD.anova <- aov(abs.Rdark_umol.cm2.hr~Group*Timepoint, data = resp.data.select)
qqnorm(resid(RD.anova))
qqline(resid(RD.anova))

boxplot(resid(RD.anova)~resp.data.select$Group)
boxplot(resid(RD.anova)~resp.data.select$Timepoint)
anova(RD.anova)

capture.output(anova(RD.anova), file = "output/Statistics/ANOVA_RD.csv")

RD.ph <- lsmeans(RD.anova, ~Group*Timepoint, adjust="tukey")

RD.ph.pairs <- multcomp::cld(RD.ph, alpha=.05, Letters=letters) #list pairwise tests and letter display
capture.output(RD.ph.pairs, file = "output/Statistics/RD_PG.csv")

#PR
PR.anova <- aov(PR~Group*Timepoint, data = resp.data.select)
qqnorm(resid(PR.anova))
qqline(resid(PR.anova))

boxplot(resid(PR.anova)~resp.data.select$Group)
boxplot(resid(PR.anova)~resp.data.select$Timepoint)
anova(PR.anova)

capture.output(anova(PR.anova), file = "output/Statistics/ANOVA_PR.csv")

PR.ph <- lsmeans(PR.anova, ~Group*Timepoint, adjust="tukey")

PR.ph.pairs <- multcomp::cld(PR.ph, alpha=.05, Letters=letters) #list pairwise tests and letter display
capture.output(PR.ph.pairs, file = "output/Statistics/PR_PG.csv")
