#Title: Metabolomic Data Analysis
#Author: KH Wong
#Date Last Modified: 20210521
#See Readme file for details

rm(list=ls()) #clears workspace 

#Read in required libraries
library("reshape")
library("Rmisc")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("arsenal")
library(gridExtra)
library(ggpubr)
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/factoextra")
# 
# install.packages("factoextra")
library(factoextra)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ropls")
library(ropls)

#BiocManager::install("mixOmics")
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


# # Or for the latest dev version:
# BiocManager::install("M3C")
# BiocManager::install("pathifier")
# BiocManager::install("pathview")
# BiocManager::install("RCy3")
# 
# install.packages("rJava")
# 
# devtools::install_github("lanagarmire/lilikoi")
# library(lilikoi)
# 
# install.packages("lilikoi")
# 
# devtools::install(build_vignettes = T)
# 
# library("FELLA")
# devtools::install_github("b2slab/FELLA")

#Overall Workflow
# 1. Import Data (Sample info, negative peaks, positive peaks)
# 2. Make a polarity column for each dataset and integrate into name
# 3. Bind data sets together 
# 4. Normalize ion counts to sample input weight and make new dataframe
# 5. Select active columns
# 6. Zero inflate values by 1000
# 7. Log normalization 
# 8. Preliminary PCA and statistics 
# 9. Preliminary heatmap 
# 10. Preliminary PSL-DA
# 11. PLS-DA Validation and cutoffs 

#Import data

sample.info <- read.csv("data/Metabolomics/Porites-Kevin_SampleInfo.csv")
peak.pos <- read.csv("data/Metabolomics/Peaks_Pos.csv")
peak.neg <- read.csv("data/Metabolomics/Peaks_Neg.csv")

#Removing unncessesary columns

peak.pos2 <- peak.pos[ -c(1:9)] 
peak.pos3 <- peak.pos2[ -c(2:9)] # if we need the blanks, change the 9 to 6

peak.neg2 <- peak.neg[ -c(1:9)] 
peak.neg3 <- peak.neg2[ -c(2:9)] # if we need the blanks, change the 9 to 6

#### Selecting polarity based on higher signal intensity ###

#obtain row means 
pos.mean <- data.frame(ID=peak.pos3[,1], Means=rowMeans(peak.pos3[,-1]))
neg.mean <- data.frame(ID=peak.neg3[,1], Means=rowMeans(peak.neg3[,-1]))

# Comparing polarity means of each metabolite, x = pos, y = neg 
cmp <- comparedf(pos.mean, neg.mean, by = "ID",
                 tol.factor = "labels")        # match only factor labels

n.diffs(cmp) #58 shared metabolites

list.diffs <- as.data.frame(do.call(cbind, diffs(cmp))) #creating a dataframe with shared metabolites as a list

df.diffs <- data.frame(matrix(unlist(list.diffs), nrow=n.diffs(cmp), byrow=F),stringsAsFactors=FALSE) #converting list to dataframe

names(df.diffs)[1:7] <- c("var.x", "var.y", "ID", "values.x", "values.y", "row.x", "row.y") #changing column names

i <- c(4, 5) #specifying values column
df.diffs[ , i] <- apply(df.diffs[ , i], 2,            # changing values column to numeric
                        function(x) as.numeric(as.character(x)))

df.diffs$selected.pol <-  ifelse(df.diffs$values.x > df.diffs$values.y, 'x', 
                                 ifelse(df.diffs$values.x < df.diffs$values.y, 'y', 'tie')) #selecting polarity with highest mean

x.diff.keep <- df.diffs %>% 
  filter(selected.pol == "x") #extracting x selection

y.diff.keep <- df.diffs %>%
  filter(selected.pol == "y") #extracting y selection

x.all.keep <- peak.pos3[!(peak.pos3$compound %in% y.diff.keep$ID),] #removing rows where the metabolite is higher in neg df
y.all.keep <- peak.neg3[!(peak.neg3$compound %in% x.diff.keep$ID),] #removing rows where the metabolite is higher in pos df


# Re-shaping dataset 
peak.pos4<- melt(x.all.keep, id= "compound") #melting dataset 
peak.neg4<- melt(y.all.keep, id= "compound") #melting dataset 

# adding polarity 
peak.pos4$polarity <- "positive"
peak.neg4$polarity <- "negative"

# Binding positive and negative datasets together
peak.all <- rbind(peak.pos4, peak.neg4)

# Renaming column 
names(peak.all)[2] <- "Sample.ID"
names(peak.all)[3] <- "Raw.IonCount"

# Merging weight sample info
peak.all2 <- merge(peak.all, sample.info, by = "Sample.ID")

# Normalization by weight

peak.all2$Raw.IonCount <- as.numeric(as.character(peak.all2$Raw.IonCount))
peak.all2$Norm.IonCount <- peak.all2$Raw.IonCount / peak.all2$Weight.mg

#Selecting columns of interest
peak.all3 <- peak.all2 %>% dplyr::select(Sample.ID, compound, Fragment.ID, Time, Treatment, Norm.IonCount)

#Reformatting dataframe so compounds are listed as column headers
peak.all4 <- peak.all3 %>% spread(compound, Norm.IonCount)

#adding 1000 to all variable to account for 0 values and log normalization 
#Log normalization (https://www.intechopen.com/books/metabolomics-fundamentals-and-applications/processing-and-visualization-of-metabolomics-data-using-r)

peak.all5 <- log((peak.all4[5:184] + 1000), 2)

norm.data <- cbind(peak.all4[1:4], peak.all5)

#### PCA ####
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
# https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/integration_biological.Rmd line 2465

scaled_pca <-prcomp(norm.data[c(5:184)], scale=TRUE, center=TRUE)
fviz_eig(scaled_pca)

coral_info<-norm.data[c(3,4)]

pca_data<- scaled_pca%>%
  augment(coral_info)%>%
  group_by(Time, Treatment)%>%
  mutate(PC1.mean = mean(.fittedPC1),
         PC2.mean = mean(.fittedPC2))

detach(package:Rmisc)
detach(package:plyr)

pca.centroids<- pca_data%>% 
  select(Time, Treatment, PC1.mean, PC2.mean)%>%
  group_by(Time, Treatment)%>%
  summarise(PC1.mean = mean(PC1.mean),
            PC2.mean = mean(PC2.mean))


#Examine PERMANOVA results.  

# scale data
vegan <- scale(norm.data[c(5:184)])

# PerMANOVA 
permanova<-adonis(vegan ~ Time*Treatment, data = norm.data, method='eu')
z_pca<-permanova$aov.tab
z_pca


#Assemble plot with background points

#1. make plot with dots

#adding percentages on axis

names(pca_data)[4] <- "PCA1"
names(pca_data)[5] <- "PCA2"
percentage <- round((scaled_pca$sdev^2) / sum((scaled_pca$sdev^2)) * 100, 2)
percentage <- paste( colnames(pca_data[4:50]), "(", paste(as.character(percentage), "%", ")", sep="") )

#setting up data to add polygons
pca_data$Time.Treatment <- paste(pca_data$Time, pca_data$Treatment)
find_hull <- function(pca_data) pca_data[chull(pca_data$PCA1, pca_data$PCA2), ]
hulls <- ddply(pca_data, "Time.Treatment", find_hull)


PCA<-ggplot(pca_data, aes(PCA1, PCA2, color=Treatment)) + 
  geom_point(size = 4, alpha=0.2, aes(shape = Time))+
  scale_colour_manual(values=c("#46008B", "#8B0046", "#468B00")) +
  scale_fill_manual(values=c("#46008B", "#8B0046", "#468B00")) + 
  scale_shape_manual(values=c(15, 17, 19)) +
  theme_classic()+
  ylim(-10,10)+
  xlim(-20,20)+
  ylab(percentage[2])+
  xlab(percentage[1])+
  geom_text(x=9, y=-8.25, label=paste("p(Time)=", z_pca$`Pr(>F)`[1]), size=4, color=ifelse(z_pca$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
  geom_text(x=9, y=-9, label=paste("p(Treatment)=", z_pca$`Pr(>F)`[2]), size=4, color=ifelse(z_pca$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
  geom_text(x=9, y=-9.75, label=paste("p(Time x Treatment)=", z_pca$`Pr(>F)`[3]), size=4, color=ifelse(z_pca$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=8), 
        legend.position="none",
        plot.background = element_blank(),
        legend.title = element_text(size=10), 
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(size=18), 
        title = element_text(size=25, face="bold"), 
        axis.title = element_text(size=18));PCA

#Add centroids  

#2. add centroids 
PCAcen<-PCA +  geom_polygon(data = hulls, alpha = 0.2, aes(color = Treatment, fill = Treatment, lty = Time)) +
  geom_point(aes(x=PC1.mean, y=PC2.mean,color=Treatment, shape = Time), data=pca.centroids, size=4, show.legend=FALSE) + 
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_colour_manual(values=c("#46008B", "#8B0046", "#468B00"), breaks = c("Control_Ambient","Bleached_Hot", "Mortality_Hot"), labels = c("Control", "Bleached", "Partial Mortality")) +
  scale_fill_manual(values=c("#46008B", "#8B0046", "#468B00"), breaks = c("Control_Ambient","Bleached_Hot", "Mortality_Hot"), labels = c("Control", "Bleached", "Partial Mortality")) + 
  scale_shape_manual(values=c(15, 17, 19)) +
  theme(legend.text = element_text(size=8), 
        legend.position=c(0.95,0.8),
        plot.background = element_blank(),
        legend.title = element_text(size=10), 
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(size=18), 
        title = element_text(size=25, face="bold"), 
        axis.title = element_text(size=18));PCAcen

#Add segments

#3. add segments
segpoints<-pca.centroids%>%
  gather(variable, value, -(Time:Treatment)) %>%
  unite(temp, Time, variable) %>%
  spread(temp, value)

PCAfull<-PCAcen + 
  geom_segment(aes(x = Day0_PC1.mean, y = Day0_PC2.mean, xend = Day37_PC1.mean, yend = Day37_PC2.mean, colour = Treatment), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = Day37_PC1.mean, y = Day37_PC2.mean, xend = Day52_PC1.mean, yend = Day52_PC2.mean, colour = Treatment), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); PCAfull


ggsave(filename="output/FullPCA_metabolomics.pdf", plot=PCAfull, dpi=300, width=12, height=10, units="in")

### PLS-DA for only final timepoint 
#https://rdrr.io/cran/RVAideMemoire/man/PLSDA.VIP.html try this

#Control vs Bleached

norm.data_2 <- column_to_rownames(norm.data, 'Sample.ID')

D52_CvsB <- norm.data_2 %>%
  filter(Time == "Day52") %>%
  filter(Treatment != "Mortality_Hot")

D52_CvsB_clean <- na.omit(D52_CvsB)

#assigning datasets 
X <- D52_CvsB[4:183]
Y <- as.factor(D52_CvsB$Treatment) 
summary(Y) ## class summary
summary(X)
dim(X) ## number of samples and features
length(Y) ## length of class membership factor = number of samples


#PLSDA without variable selection
MyResult.plsda <- plsda(X, Y, ncomp = 2) # 1 Run the method
plotIndiv(MyResult.plsda)    # 2 Plot the samples

plotVar(MyResult.plsda, cutoff = 0.7)    

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,ellipse = TRUE, title="Day 52 - Bleached vs Control")

MyResult.plsda2 <- plsda(X,Y, ncomp=6) #number of components is #classes-1
selectVar(MyResult.plsda2, comp=1)$value

plotLoadings(MyResult.plsda2, comp = 1, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1
plotLoadings(MyResult.plsda2, comp = 2, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1

comp1.select.metabolites.all <- data.frame(selectVar(MyResult.plsda2, comp = 1)$value)

# component validation 
set.seed(200) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 5, #figure out why only folds =5 works
                     progressBar = FALSE, auc = TRUE, nrepeat = 10) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

#ROC Curve
auc.plsda = auroc(MyResult.plsda2, roc.comp = 2)

#VIP Extraction
D52_CvsB_VIP <- PLSDA.VIP(MyResult.plsda2)
D52_CvsB_VIP_DF <- as.data.frame(D52_CvsB_VIP[["tab"]])

######Control vs Partial Mortality

D52_CvsP <- norm.data_2 %>%
  filter(Time == "Day52") %>%
  filter(Treatment != "Bleached_Hot")

D52_CvsP_clean <- na.omit(D52_CvsP)

#assigning datasets 
X <- D52_CvsP[4:183]
Y <- as.factor(D52_CvsP$Treatment) 
summary(Y) ## class summary
summary(X)
dim(X) ## number of samples and features
length(Y) ## length of class membership factor = number of samples


#PLSDA without variable selection
MyResult.plsda <- plsda(X, Y, ncomp = 2) # 1 Run the method
plotIndiv(MyResult.plsda)    # 2 Plot the samples

plotVar(MyResult.plsda, cutoff = 0.7)    

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,ellipse = TRUE, title="Day 52 - Bleached vs Control")

MyResult.plsda2 <- plsda(X,Y, ncomp=6) #number of components is #classes-1
selectVar(MyResult.plsda2, comp=1)$value

plotLoadings(MyResult.plsda2, comp = 1, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1
plotLoadings(MyResult.plsda2, comp = 2, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1

comp1.select.metabolites.all <- data.frame(selectVar(MyResult.plsda2, comp = 1)$value)

# component validation 
set.seed(200) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 5, #figure out why only folds =5 works
                     progressBar = FALSE, auc = TRUE, nrepeat = 10) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

#ROC Curve
auc.plsda = auroc(MyResult.plsda2, roc.comp = 2)

#VIP Extraction
D52_CvsP_VIP <- PLSDA.VIP(MyResult.plsda2)
D52_CvsP_VIP_DF <- as.data.frame(D52_CvsP_VIP[["tab"]])

######Bleached vs Partial Mortality

D52_BvsP <- norm.data_2 %>%
  filter(Time == "Day52") %>%
  filter(Treatment != "Control_Ambient")

D52_BvsP_clean <- na.omit(D52_BvsP)

#assigning datasets 
X <- D52_BvsP[4:183]
Y <- as.factor(D52_BvsP$Treatment) 
summary(Y) ## class summary
summary(X)
dim(X) ## number of samples and features
length(Y) ## length of class membership factor = number of samples


#PLSDA without variable selection
MyResult.plsda <- plsda(X, Y, ncomp = 2) # 1 Run the method
plotIndiv(MyResult.plsda)    # 2 Plot the samples

plotVar(MyResult.plsda, cutoff = 0.7)    

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,ellipse = TRUE, title="Day 52 - Bleached vs Control")

MyResult.plsda2 <- plsda(X,Y, ncomp=6) #number of components is #classes-1
selectVar(MyResult.plsda2, comp=1)$value

plotLoadings(MyResult.plsda2, comp = 1, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1
plotLoadings(MyResult.plsda2, comp = 2, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1

comp1.select.metabolites.all <- data.frame(selectVar(MyResult.plsda2, comp = 1)$value)

# component validation 
set.seed(200) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 5, #figure out why only folds =5 works
                     progressBar = FALSE, auc = TRUE, nrepeat = 10) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

#ROC Curve
auc.plsda = auroc(MyResult.plsda2, roc.comp = 2)

#VIP Extraction
D52_BvsP_VIP <- PLSDA.VIP(MyResult.plsda2)
D52_BvsP_VIP_DF <- as.data.frame(D52_BvsP_VIP[["tab"]])


####### Overlaps for VIPs >1 #########

# Converting row names to column
D52_CvsB_VIP_table <- rownames_to_column(D52_CvsB_VIP_DF, var = "Metabolite")
D52_CvsP_VIP_table <- rownames_to_column(D52_CvsP_VIP_DF, var = "Metabolite")
D52_BvsP_VIP_table <- rownames_to_column(D52_BvsP_VIP_DF, var = "Metabolite")

# Filtering for VIP > 1
D52_CvsB_VIP_1 <- D52_CvsB_VIP_table %>% 
  filter(VIP >= 1)

D52_CvsP_VIP_1 <- D52_CvsP_VIP_table %>% 
  filter(VIP >= 1)

D52_BvsP_VIP_1 <- D52_BvsP_VIP_table %>% 
  filter(VIP >= 1)

D52_CvsB_VIP_1 %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_point() +
  ylab("Metabolite") +
  xlab("VIP Score") +
  ggtitle("Control vs Bleached") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

D52_CvsP_VIP_1 %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_point() +
  ylab("Metabolite") +
  xlab("VIP Score") +
  ggtitle("Control vs Partial Mortality") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




### Compare for Venn Diagram (https://github.com/gaospecial/ggVennDiagram)

nrow(D52_CvsB_VIP_1)
nrow(D52_CvsP_VIP_1)

library("VennDiagram")

venn.x <- list(
  D52_CvsB = sample(D52_CvsB_VIP_1$Metabolite), 
  D52_CvsP = sample(D52_CvsP_VIP_1$Metabolite)
)

myCol <- c("#8B0046", "#468B00")

venn.diagram(venn.x,
             filename = 'output/Metabolomics/D52_BvsP_venn.png',
             output = TRUE,
             height = 480, 
             width = 480,
             resolution = 300,
             category.names = c("Bleached", "Partial Mortality"),
             lwd = 2,
             col = c("#8B0046", "#468B00"),
             fill = c(alpha("#8B0046",0.3), alpha("#468B00", 0.3)),
             cex = 0.5,
             fontfamily = "sans",
             fontface = "bold",
             cat.cex = 0.5,
             cat.default.pos = "outer",
             cat.pos = c(12, -12),
             cat.dist = c(-0.45, -0.45),
             cat.fontfamily = "sans",
             cat.fontface = "bold"
             )



