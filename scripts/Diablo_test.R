# library(mixOmics)
# data(breast.TCGA)
# # extract training data and name each data frame
# 
# str(breast.TCGA$data.train$mrna)
# X <- list(mRNA = breast.TCGA$data.train$mrna, # all rows are samples, columns are compounds/genes
#           miRNA = breast.TCGA$data.train$mirna, 
#           protein = breast.TCGA$data.train$protein)
# Y <- breast.TCGA$data.train$subtype
# summary(Y)
# 
# list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5)) # not sure what this does, maybe training set?
# 
# MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
# plotIndiv(MyResult.diablo) ## sample plot
# plotVar(MyResult.diablo) ## variable plot
# 
# MyResult.diablo2 <- block.splsda(X, Y)
# 
# ## Customization
# plotIndiv(MyResult.diablo2, 
#           ind.names = FALSE, 
#           legend=TRUE, cex=c(1,2,3),
#           title = 'BRCA with DIABLO')
# plotVar(MyResult.diablo2, var.names = c(FALSE, FALSE, TRUE),
#         legend=TRUE, pch=c(16,16,1))
# 
# #Plot diablo
# plotDiablo(MyResult.diablo2, ncomp = 1)
# 
# #circos plot
# corMat <- circosPlot(MyResult.diablo2, cutoff=0.7, line = TRUE)
# 
# #cimdiablo (heatmap)
# # minimal example with margins improved: (this does not work for me)
# #cimDiablo(MyResult.diablo, margin=c(8,20))
# # extended example:
# #cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), comp = 1, margin=c(8,20), legend.position = "right")
# 
# 
# # Plot Loadings
# #plotLoadings(MyResult.diablo, contrib = "max")
# plotLoadings(MyResult.diablo2, comp = 2, contrib = "max")
# 
# network(MyResult.diablo, blocks = c(1,2,3),
#         color.node = c('darkorchid', 'brown1', 'lightgreen'), 
#         cutoff = 0.6, save = 'jpeg', name.save = 'DIABLOnetwork')

##### MY DATA

BiocManager::install("mixOmicsTeam/mixOmics@devel")
library(mixOmics)

# Import data, make sure all rows are in order in each dataset
Metabolome <- read.csv("output/Metabolomics/metab_df.csv", check.names = FALSE)
Microbiome <- read.csv("output/16S/processed_data/ASV_df.csv", check.names = FALSE)
Gene <- read.csv("output/TagSeq/gvst_df.csv", check.names = FALSE)
metadata <- read.csv("data/Metadata/Frag_Sample_Info.csv")

metadata$GroupDay <- paste(metadata$Group, metadata$Day, sep = "-")

# Need to clean the datasets: 1) remove sample that is missing from mirco DF (R19-52)
colnames(Microbiome)[1] <- "ColDay"
colnames(Metabolome)[1] <- "ColDay"
colnames(Gene)[1] <- "ColDay"

Metabolome2 <- Metabolome %>% filter(ColDay != "R19-52")
Gene2 <- Gene %>% filter(ColDay != "R19-52")

metadata2 <- metadata %>% filter(ColDay != "R19-52")

# Filtering for just day 52
metadata3 <- metadata2 %>% filter(grepl("-52", ColDay, fixed = TRUE))

# 
# # Make ColDay rownames
# rownames(Metabolome2) <- paste0(Metabolome2$ColDay)
# rownames(Microbiome) <- paste0(Microbiome$ColDay)
# rownames(Gene2) <- paste0(Gene2$ColDay)
# 
# Metabolome_final <- Metabolome2 %>% dplyr::select(-ColDay)
# Microbiome_final <- Microbiome %>% dplyr::select(-ColDay)
# Gene_final <- Gene2 %>% dplyr::select(-ColDay)
# 
# 
# 
# # extract training data and name each data frame
# X <- list(Metabolome = Metabolome_final, # all rows are samples, columns are compounds/genes
#           Microbiome = Microbiome_final)
# Y <- as.factor(metadata2$GroupDay)
# summary(Y)
# 
# length(unique(lapply(X, function(x) rownames))) # This should output as 1
# 
# MyResult.diablo <- block.splsda(X, Y)
# 
# plotIndiv(MyResult.diablo) ## sample plot
# plotVar(MyResult.diablo) ## variable plot
# 
# ## Customization
# plotIndiv(MyResult.diablo, 
#           ind.names = FALSE, 
#           legend=TRUE, cex=c(1,2,3),
#           title = 'BRCA with DIABLO')
# plotVar(MyResult.diablo, var.names = c(FALSE, FALSE, TRUE),
#         legend=TRUE, pch=c(16,16,1))
# 
# #Plot diablo
# plotDiablo(MyResult.diablo, ncomp = 1)
# 
# #circos plot
# corMat <- circosPlot(MyResult.diablo, cutoff=0.7, showIntraLinks = TRUE, line = TRUE)
# 
# #cimdiablo (heatmap)
# # minimal example with margins improved: (this does not work for me)
# #cimDiablo(MyResult.diablo, margin=c(8,20))
# # extended example:
# #cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), comp = 1, margin=c(8,20), legend.position = "right")
# 
# 
# # Plot Loadings
# #plotLoadings(MyResult.diablo, contrib = "max")
# plotLoadings(MyResult.diablo, comp = 2, contrib = "max")
# 
# network(MyResult.diablo, blocks = c(1,2,3),
#         color.node = c('darkorchid', 'brown1', 'lightgreen'), 
#         cutoff = 0.6, save = 'jpeg', name.save = 'DIABLOnetwork')
# 


###### Filtering for just day 52

metadata3 <- metadata2 %>% filter(grepl("-52", ColDay, fixed = TRUE))
Metabolome3 <- Metabolome2 %>% filter(grepl("-52", ColDay, fixed = TRUE))
Microbiome3 <- Microbiome %>% filter(grepl("-52", ColDay, fixed = TRUE))
Gene3 <- Gene2 %>% filter(grepl("-52", ColDay, fixed = TRUE))

# Make ColDay rownames
rownames(Metabolome3) <- paste0(Metabolome3$ColDay)
rownames(Microbiome3) <- paste0(Microbiome3$ColDay)
rownames(Gene3) <- paste0(Gene3$ColDay)

Metabolome_final <- Metabolome3 %>% dplyr::select(-ColDay)
Microbiome_final <- Microbiome3 %>% dplyr::select(-ColDay)
Gene_final <- Gene3 %>% dplyr::select(-ColDay)

# extract training data and name each data frame
X <- list(Metabolome = Metabolome_final, # all rows are samples, columns are compounds/genes
          Microbiome = Microbiome_final,
          Genes = Genbe_Final)
Y <- as.factor(metadata3$GroupDay)
summary(Y)


### Metabolomics PLSDA ###

#assigning datasets 
X <- metabolomics[3:184]
Y <- as.factor(metabolomics$Lifestage) #select treatment names
Y
MyResult.plsda_metab <- plsda(X,Y, ncomp=9) #number of components is classes-1
#MyResult.splsda_metab <- splsda(X,Y, ncomp=12) #number of components is classes-1
plotIndiv(MyResult.plsda_metab) 

### 



length(unique(lapply(X, function(x) rownames))) # This should output as 1

MyResult.diablo <- block.splsda(X, Y)

plotIndiv(MyResult.diablo) ## sample plot
plotVar(MyResult.diablo) ## variable plot

## Customization
plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2,3),
          title = 'BRCA with DIABLO')
plotVar(MyResult.diablo, var.names = c(FALSE, FALSE, TRUE),
        legend=TRUE, pch=c(16,16,1))

#Plot diablo
plotDiablo(MyResult.diablo, ncomp = 1)

#circos plot
corMat <- circosPlot(MyResult.diablo, cutoff=0.7, line = TRUE)

#cimdiablo (heatmap)
# minimal example with margins improved: (this does not work for me)
#cimDiablo(MyResult.diablo, margin=c(8,20))
# extended example:
#cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), comp = 1, margin=c(8,20), legend.position = "right")


# Plot Loadings
#plotLoadings(MyResult.diablo, contrib = "max")
plotLoadings(MyResult.diablo, comp = 2, contrib = "max")

network(MyResult.diablo, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        cutoff = 0.6, save = 'jpeg', name.save = 'DIABLOnetwork')

