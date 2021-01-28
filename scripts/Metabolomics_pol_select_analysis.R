#Title: Metabolomic Data Analysis
#Author: KH Wong
#Date Last Modified: 20201208
#See Readme file for details

rm(list=ls()) #clears workspace 

#Read in required libraries
library("reshape")
library("dplyr")
library("tidyverse")
library("Rmisc")
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

#subsetting for pos and neg analysis
peak.all.pos2 <- peak.all2 %>% filter(polarity == "positive")
peak.all.neg2 <- peak.all2 %>% filter(polarity == "negative")

#Selecting columns of interest
peak.all3 <- peak.all2 %>% dplyr::select(Sample.ID, compound, Fragment.ID, Time, Treatment, Norm.IonCount)

peak.all.pos3 <- peak.all.pos2 %>% dplyr::select(Sample.ID, compound, Fragment.ID, Time, Treatment, Norm.IonCount)
peak.all.neg3 <- peak.all.neg2 %>% dplyr::select(Sample.ID, compound, Fragment.ID, Time, Treatment, Norm.IonCount)

#Reformatting dataframe so compounds are listed as column headers
peak.all4 <- peak.all3 %>% spread(compound, Norm.IonCount)

peak.all.pos4 <- peak.all.pos3 %>% spread(compound, Norm.IonCount)
peak.all.neg4 <- peak.all.neg3 %>% spread(compound, Norm.IonCount)

#Making row names as sample ID
peak.all5 <- column_to_rownames(peak.all4, 'Sample.ID')

peak.all.pos5 <- column_to_rownames(peak.all.pos4, 'Sample.ID')
peak.all.neg5 <- column_to_rownames(peak.all.neg4, 'Sample.ID')

#Making treatment and time groups 
peak.all5$Treatment_Time <- paste(peak.all5$Treatment, peak.all5$Time, sep = "_")

peak.all.pos5$Treatment_Time <- paste(peak.all.pos5$Treatment, peak.all.pos5$Time, sep = "_")
peak.all.neg5$Treatment_Time <- paste(peak.all.neg5$Treatment, peak.all.neg5$Time, sep = "_")

#Selecting active rows and columns (in that order)
peak.all6 <- peak.all5[1:45, 4:184]

peak.all.pos6 <- peak.all.pos5[1:45, 4:85]
peak.all.neg6 <- peak.all.neg5[1:45, 4:102]

#adding 1000 to all variable to account for 0 values and log normalization 
peak.all7 <- peak.all6 + 1000

peak.all.pos7 <- peak.all.pos6 + 1000
peak.all.neg7 <- peak.all.neg6 + 1000

#Log normalization (https://www.intechopen.com/books/metabolomics-fundamentals-and-applications/processing-and-visualization-of-metabolomics-data-using-r)
peak.all.active <- log(peak.all7, 2)

peak.pos.active <- log(peak.all.pos7, 2)
peak.neg.active <- log(peak.all.neg7, 2)

#### PCA ####
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/



### Grouped PCA 

#All polarities
res.pca <- prcomp(peak.all.active, scale = TRUE)
fviz_eig(res.pca)

pca.all <- fviz_pca_ind(res.pca,
             habillage=peak.all5$Treatment_Time,
             palette = c("seagreen1", "seagreen3", "seagreen4", "skyblue1", "skyblue3", "skyblue4", "salmon1", "salmon3", "salmon4"),
             label = "none",
             title = "All Polarities",
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             legend = c(0.8,0.8),
             repel = FALSE
) + scale_shape_manual(values = c(17, 17, 17, 15, 15, 15, 19, 19, 19)) 

#Positive polarity
res.pca.pos <- prcomp(peak.pos.active, scale = TRUE)
fviz_eig(res.pca.pos)

pca.pos <- fviz_pca_ind(res.pca.pos,
                        habillage=peak.all5$Treatment_Time,
                        palette = c("seagreen1", "seagreen3", "seagreen4", "skyblue1", "skyblue3", "skyblue4", "salmon1", "salmon3", "salmon4"),
                        label = "none",
                        title = "Positive Polarity",
                        addEllipses = TRUE, # Concentration ellipses
                        ellipse.type = "confidence",
                        legend.title = "Groups",
                        repel = FALSE, 
                        legend = "none",
) + scale_shape_manual(values = c(17, 17, 17, 15, 15, 15, 19, 19, 19))

#Negative polarity
res.pca.neg <- prcomp(peak.neg.active, scale = TRUE)
fviz_eig(res.pca.neg)

pca.neg <- fviz_pca_ind(res.pca.neg,
                        habillage=peak.all5$Treatment_Time,
                        palette = c("seagreen1", "seagreen3", "seagreen4", "skyblue1", "skyblue3", "skyblue4", "salmon1", "salmon3", "salmon4"),
                        label = "none",
                        title = "Negative Polarity",
                        addEllipses = TRUE, # Concentration ellipses
                        ellipse.type = "confidence",
#                        legend.title = "Groups",
                        repel = FALSE,
                        legend = "none", 
) + scale_shape_manual(values = c(17, 17, 17, 15, 15, 15, 19, 19, 19))

Metab.Figs <- ggarrange(pca.pos, pca.neg, pca.all, common.legend = TRUE, legend = "right", ncol=3) 

ggsave(file="output/Metab.pcas.pdf", Metab.Figs, width = 11, height = 5, units = c("in"))

 ###Assessing PCA
# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

# -----------------------------------

### PLSDA
#https://mixomicsteam.github.io/Bookdown/plsda.html

#assigning datasets 
X <- peak.all.active
Y <- as.factor(peak.all5$Treatment_Time[1:45]) 
summary(Y) ## class summary
dim(X) ## number of samples and features
length(Y) ## length of class membership factor = number of samples

#PLSDA without variable selection
MyResult.plsda <- plsda(X, Y) # 1 Run the method
plotIndiv(MyResult.plsda)    # 2 Plot the samples

plotVar(MyResult.plsda, cutoff = 0.7)    

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,ellipse = TRUE, title="PLS-DA - ALL DATA")

MyResult.plsda2 <- plsda(X,Y, ncomp=8) #number of components is #classes-1
selectVar(MyResult.plsda2, comp=1)$value

plotLoadings(MyResult.plsda2, comp = 1, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1
plotLoadings(MyResult.plsda2, comp = 2, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1


# component validation 
MyResult.plsda2 <- plsda(X,Y, ncomp=8)
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 5, #figure out why only folds =5 works
                     progressBar = FALSE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

# Assessing optimal number of variables to keep 

list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 5, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(ncomp))

# Final model based on component and variable validation 

MyResult.splsda.final <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX)

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")

selectVar(MyResult.splsda.final, comp=2)$value

plotLoadings(MyResult.splsda.final, comp = 1, contrib = 'max', method = 'median')
plotLoadings(MyResult.splsda.final, comp = 2, contrib = 'max', method = 'median')
dev.off()

auc.splsda = auroc(MyResult.splsda.final, roc.comp = ncomp)

## Performance validation of final spls-da
perf.splsda <- perf(MyResult.splsda.final, validation = "Mfold", folds = 5, 
                    progressBar = FALSE, auc = TRUE, nrepeat = 50) 

perf.splsda$auc
perf.splsda$error.rate

plot(perf.splsda)
dev.off()

par(mfrow=c(1,4))
plot(perf.splsda$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.splsda$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)
plot(perf.splsda$features$stable[[3]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 3', las =2)
plot(perf.splsda$features$stable[[4]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 4', las =2)
dev.off()


# here we match the selected variables to the stable features
ind.match = match(selectVar(MyResult.splsda.final, comp = 1)$name, 
                  names(perf.splsda$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.splsda$features$stable[[1]][ind.match])

comp1.select.metabolites <- data.frame(selectVar(MyResult.splsda.final, comp = 1)$value, Freq)


#Heatmap of loading 1
cim(MyResult.splsda.final, comp=1, title="Component 1")
dev.off()

# 
# #Calculating Q and Y for spls
# 
# X <- ... # X data
# Y <- ... # Y data (factor or discrete data)
# Y.mat <- unmap(Y) # creates a dummy matrix
# res <- spls(X, Y.mat, …)
# val <- perf(res, criterion = c("R2", "Q2"))
# val
# 
# ncomp = 4
# result.spls <- spls(X, Y, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'regression')
# tune.spls <- perf(result.spls, validation = 'Mfold', folds = 10,
#                   criterion = 'all', progressBar = FALSE)
# 
# 
# 
# ##-----------
# 

### Try Multi-level sPLS-DA next 

############## 



