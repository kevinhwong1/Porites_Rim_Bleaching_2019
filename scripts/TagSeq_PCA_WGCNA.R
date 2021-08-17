#title: "Porites July Bleaching TagSeq WGCNA"
#author: "Erin Chille - Ariana Huffmyer - Kevin Wong"
#date: "Last updated 2021/08/016"

# Load Packages

if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse') 
if ("genefilter" %in% rownames(installed.packages()) == 'FALSE') install.packages('genefilter') 
if ("DESeq2" %in% rownames(installed.packages()) == 'FALSE') install.packages('DESeq2') 
if ("RColorBrewer" %in% rownames(installed.packages()) == 'FALSE') install.packages('RColorBrewer') 
if ("WGCNA" %in% rownames(installed.packages()) == 'FALSE') install.packages('WGCNA') 
if ("flashClust" %in% rownames(installed.packages()) == 'FALSE') install.packages('flashClust') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("ComplexHeatmap" %in% rownames(installed.packages()) == 'FALSE') install.packages('ComplexHeatmap') 
if ("goseq" %in% rownames(installed.packages()) == 'FALSE') install.packages('goseq') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("clusterProfiler" %in% rownames(installed.packages()) == 'FALSE') install.packages('clusterProfiler') 
if ("pheatmap" %in% rownames(installed.packages()) == 'FALSE') install.packages('pheatmap') 
library(BiocManager)
#if ("simplifyEnrichment" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("simplifyEnrichment") 
library("tidyverse")
library("genefilter")
library("DESeq2")
library("RColorBrewer")
library("WGCNA")
library("flashClust")
library("gridExtra")
library("ComplexHeatmap")
library("goseq")
library("dplyr")
library("clusterProfiler")
library("simplifyEnrichment") 
library("pheatmap")
library("grDevices")

#### Data input, cleaning and pre-processing

treatmentinfo <- read.csv("data/Molecular/Sample_Info.csv", header = TRUE, sep = ",")
head(treatmentinfo)
treatmentinfo$Day <- as.factor(treatmentinfo$Day)

gcount <- read.table("output/TagSeq/full_count_matrix.txt", header = TRUE, row.names="Name")
head(gcount)
dim(gcount)

#### Quality-filter gene counts

#set filter values for PoverA, P=100% percent of the samples have counts over A=10. This means that only 5 out of 45 (0.11) samples need to have counts over 10. 
filt <- filterfun(pOverA(0.1,10))

gfilt <- genefilter(gcount, filt)

#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),]) ###### look into why this is different than gfilt

#How many rows do we have before and after filtering?
nrow(gcount) #30740
nrow(gcount_filt) #5646 for 10, X for 5

### Quality-check of datasets  
#In order for the DESeq2 algorithms to work, the SampleIDs on the treatmentinfo file and count matrices have to match exactly and in the same order. The following R clump will check to make sure that these match.

#Checking that all row and column names match. Should return "TRUE"
all(rownames(treatmentinfo$sample_id_full) %in% colnames(gcount_filt))
all(rownames(treatmentinfo$sample_id_full) == colnames(gcount_filt)) 


### Read normalization
#We are now going normalize our read counts using VST-normalization in DESeq2

#### Construct the DESeq2 dataset

#Create a DESeqDataSet design from gene count matrix and labels. Here we set the design to look at time_point to test for any differences in gene expression across timepoints.



#Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = round(gcount_filt),
                               colData = treatmentinfo,
                               design = ~Group+Day)

#### Log-transform the count data
#First we are going to log-transform the data using a variance stabilizing transforamtion (VST). This is only for visualization purposes. Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.

#To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.

SF.gdds <- estimateSizeFactors(gdds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 for us to use vst
print(sizeFactors(SF.gdds)) #View size factors


#Our size factors are all less than 4, so we can use VST! VST = variance stabilizing transformation to minimize effects of small counts and normalize wrt library size

gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst), 3) #view transformed gene count data

##### Plot a heatmap of sample-to-sample distances
gsampleDists <- dist(t(assay(gvst))) #calculate distance matix
gsampleDistMatrix <- as.matrix(gsampleDists) #distance matrix
rownames(gsampleDistMatrix) <- colnames(gvst) #assign row names
colnames(gsampleDistMatrix) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pdf("output/TagSeq/Full_TagSeq_heatmap.pdf")
pht=pheatmap(gsampleDistMatrix, #plot matrix
             clustering_distance_rows=gsampleDists, #cluster rows
             clustering_distance_cols=gsampleDists, #cluster columns
             col=colors) #set colors
draw(pht)
dev.off()

##### Principal component plot of samples

gPCAdata <- plotPCA(gvst, intgroup = c("Group", "Day"), returnData=TRUE)
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
allgenesfilt_PCA <- ggplot(gPCAdata, aes(PC1, PC2, shape=Day, color = Group)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #xlim(-40,40)+ 
  #ylim(-40,40)+
  coord_fixed()+
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines 
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) # + #Set the plot background
#theme(legend.position = ("none")) #set title attributes
allgenesfilt_PCA

### PCA

pca.centroids<- gPCAdata%>% 
  dplyr::select(Day, Group, PC1, PC2)%>%
  dplyr::group_by(Day, Group)%>%
  dplyr::summarise(PC1.mean = mean(PC1),
                   PC2.mean = mean(PC2))

gPCAdata$Day.Group <- paste(gPCAdata$Day, gPCAdata$Group)
find_hull <- function(gPCAdata) gPCAdata[chull(gPCAdata$PC1, gPCAdata$PC2), ]
hulls <- plyr::ddply(gPCAdata, "Day.Group", find_hull)

PCA<-ggplot(gPCAdata, aes(PC1, PC2, color=Group)) + 
  geom_point(size = 4, alpha=0.2, aes(shape = Day))+
  scale_colour_manual(values=c("#46008B", "#8B0046", "#468B00")) +
  scale_fill_manual(values=c("#46008B", "#8B0046", "#468B00")) + 
  scale_shape_manual(values=c(15, 17, 19)) +
  theme_classic()+
#  ylim(-7,12)+
#  xlim(-6,6)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#  geom_text(x=4.5, y=-4.75, label=paste("p(Day)=", z_pca$`Pr(>F)`[1]), size=4, color=ifelse(z_pca$`Pr(>F)`[1] < 0.05, "black", "darkgray")) + 
#  geom_text(x=4.5, y=-5.5, label=paste("p(Group)=", z_pca$`Pr(>F)`[2]), size=4, color=ifelse(z_pca$`Pr(>F)`[2] < 0.05, "black", "darkgray")) + 
#  geom_text(x=4.5, y=-6.25, label=paste("p(Day x Group)=", z_pca$`Pr(>F)`[3]), size=4, color=ifelse(z_pca$`Pr(>F)`[3] < 0.05, "black", "darkgray")) + 
  theme(legend.text = element_text(size=8), 
        legend.position="none",
        plot.background = element_blank(),
        legend.title = element_text(size=10), 
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(size=18), 
        title = element_text(size=25, face="bold"), 
        axis.title = element_text(size=18))


#Add centroids  

#2. add centroids 
PCAcen<-PCA +  geom_polygon(data = hulls, alpha = 0.2, aes(color = Group, fill = Group, lty = Day)) +
  geom_point(aes(x=PC1.mean, y=PC2.mean,color=Group, shape = Day), data=pca.centroids, size=4, show.legend=FALSE) + 
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_colour_manual(values=c("#46008B", "#8B0046", "#468B00"), breaks = c("Control","Bleached", "Partial-Mortality"), labels = c("Control", "Bleached", "Partial Mortality")) +
  scale_fill_manual(values=c("#46008B", "#8B0046", "#468B00"), breaks = c("Control","Bleached", "Partial-Mortality"), labels = c("Control", "Bleached", "Partial Mortality")) + 
  scale_shape_manual(values=c(15, 17, 19)) +
  theme(legend.text = element_text(size=8), 
        legend.position=c(0.95,0.85),
        plot.background = element_blank(),
        legend.title = element_text(size=10), 
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(size=18), 
        title = element_text(size=25, face="bold"), 
        axis.title = element_text(size=18))

#Add segments

#3. add segments
segpoints<-pca.centroids%>%
  gather(variable, value, -(Day:Group)) %>%
  unite(temp, Day, variable) %>%
  spread(temp, value) 

names(segpoints)[2] <- "Day0_PC1.mean"
names(segpoints)[3] <- "Day0_PC2.mean"
names(segpoints)[4] <- "Day37_PC1.mean"
names(segpoints)[5] <- "Day37_PC2.mean"
names(segpoints)[6] <- "Day52_PC1.mean"
names(segpoints)[7] <- "Day52_PC2.mean"

PCAfull<-PCAcen + 
  geom_segment(aes(x = Day0_PC1.mean, y = Day0_PC2.mean, xend = Day37_PC1.mean, yend = Day37_PC2.mean, colour = Group), data = segpoints, size=2, show.legend=FALSE) +
  geom_segment(aes(x = Day37_PC1.mean, y = Day37_PC2.mean, xend = Day52_PC1.mean, yend = Day52_PC2.mean, colour = Group), data = segpoints, size=2, arrow = arrow(length=unit(0.5,"cm")), show.legend=FALSE); PCAfull

ggsave(filename="output/TagSeq/Full_PCA_TagSeq.pdf", plot=PCAfull, dpi=300, width=12, height=10, units="in")
