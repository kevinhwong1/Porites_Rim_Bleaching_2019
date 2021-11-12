

D52_uniqueB_BvsP_up_plot <- D52_uniqueB_BvsP_up_VIP %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, color = "#8B0046") +
  xlim(1.2, 2) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Bleached Unique") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17))



D52_CvsB_VIP_UP <- D52_CvsB_VIP_sig_up %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point() +
  xlim(1.2, 2) +
  ylab("Accumulated") +
  xlab("") +
  ggtitle("Bleached") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

D52_CvsB_VIP_DOWN <- D52_CvsB_VIP_sig_down %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point() +
  xlim(1.2, 2) +
  ylab("Depleted") +
  xlab("VIP Score") +
  #  ggtitle("Bleached: Depleted") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

D52_CvsP_VIP_UP <- D52_CvsP_VIP_sig_up %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point() +
  xlim(1.2, 2) +
  ylab("") +
  xlab("") +
  ggtitle("Partial Mortality") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

D52_CvsP_VIP_DOWN <- D52_CvsP_VIP_sig_down %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point() +
  xlim(1.2, 2) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Partial Mortality: Depleted") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

D52_All_VIP_plot <- plot_grid(D52_CvsB_VIP_UP, D52_CvsP_VIP_UP, D52_CvsB_VIP_DOWN, D52_CvsP_VIP_DOWN, 
                              align = "v", 
                              ncol = 2, 
                              rel_heights = c(0.95, 0.35), 
                              labels = c("A", "B", "C", "D"))

D52_All_VIP_plot
#ggsave(filename="../output/Metabolomics/Day52_VIP.pdf", plot=D52_All_VIP_plot, dpi=300, width=7, height=9, units="in")

```


### VIP plotting after t-test validation and comparing overlapping metabolites

``` {r Plotting, echo=TRUE, warning=FALSE, message=FALSE}

# Up regulated metabolites that overlap between B and P

D52_overlap_BvsP_up_1 <- as.data.frame(intersect(D52_CvsB_VIP_sig_up$Metabolite, D52_CvsP_VIP_sig_up$Metabolite))
names(D52_overlap_BvsP_up_1)[1] <- 'Metabolite'
D52_overlap_BvsP_up_2 <- merge(D52_overlap_BvsP_up_1, D52_CvsB_VIP_sig_up, by="Metabolite")
D52_overlap_BvsP_up_VIP <- merge(D52_overlap_BvsP_up_2, D52_CvsP_VIP_sig_up, by="Metabolite")
names(D52_overlap_BvsP_up_VIP)[2] <- 'Bleached'
names(D52_overlap_BvsP_up_VIP)[3] <- 'Partial Mortality'

D52_overlap_BvsP_up_VIP_comb <- gather(D52_overlap_BvsP_up_VIP, key = "Group", value = "VIP", 'Bleached', 'Partial Mortality')

# Up regulated metabolites that are unique to B

D52_uniqueB_BvsP_up <- as.data.frame(setdiff(D52_CvsB_VIP_sig_up$Metabolite, D52_CvsP_VIP_sig_up$Metabolite))
names(D52_uniqueB_BvsP_up)[1] <- 'Metabolite'
D52_uniqueB_BvsP_up_VIP <- merge(D52_uniqueB_BvsP_up, D52_CvsB_VIP_sig_up, by="Metabolite")

# Up regulated metabolites that are unique to P

D52_uniqueP_BvsP_up <- as.data.frame(setdiff(D52_CvsP_VIP_sig_up$Metabolite, D52_CvsB_VIP_sig_up$Metabolite))
names(D52_uniqueP_BvsP_up)[1] <- 'Metabolite'
D52_uniqueP_BvsP_up_VIP <- merge(D52_uniqueP_BvsP_up, D52_CvsP_VIP_sig_up, by="Metabolite")

# Down regulated metabolites that overlap between B and P

D52_overlap_BvsP_down_1 <- as.data.frame(intersect(D52_CvsB_VIP_sig_down$Metabolite, D52_CvsP_VIP_sig_down$Metabolite))
names(D52_overlap_BvsP_down_1)[1] <- 'Metabolite'
D52_overlap_BvsP_down_2 <- merge(D52_overlap_BvsP_down_1, D52_CvsB_VIP_sig_down, by="Metabolite")
D52_overlap_BvsP_down_VIP <- merge(D52_overlap_BvsP_down_2, D52_CvsP_VIP_sig_down, by="Metabolite")
names(D52_overlap_BvsP_down_VIP)[2] <- 'Bleached'
names(D52_overlap_BvsP_down_VIP)[3] <- 'Partial Mortality'

D52_overlap_BvsP_down_VIP_comb <- gather(D52_overlap_BvsP_down_VIP, key = "Group", value = "VIP", 'Bleached', 'Partial Mortality')


# Down regulated metabolites that are unique to B 
# There are none 

# Down regulated metabolites that are unique to P

D52_uniqueP_BvsP_down <- as.data.frame(setdiff(D52_CvsP_VIP_sig_down$Metabolite, D52_CvsB_VIP_sig_down$Metabolite))
names(D52_uniqueP_BvsP_down)[1] <- 'Metabolite'
D52_uniqueP_BvsP_down_VIP <- merge(D52_uniqueP_BvsP_down, D52_CvsP_VIP_sig_down, by="Metabolite")


# Plotting

D52_uniqueB_BvsP_up_plot <- D52_uniqueB_BvsP_up_VIP %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, color = "#8B0046") +
  xlim(1.2, 2) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Bleached Unique") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17))

D52_uniqueP_BvsP_up_plot <- D52_uniqueP_BvsP_up_VIP %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, color = "#468B00") +
  xlim(1.2, 2) +
  ylab("") +
  xlab("") +
  #  ggtitle("Partial Mortality Unique") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17))

D52_uniqueP_BvsP_down_plot <- D52_uniqueP_BvsP_down_VIP %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, color = "#468B00") +
  xlim(1.2, 2) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Partial Mortality: Depleted") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17))

D52_overlap_BvsP_up_plot <- D52_overlap_BvsP_up_VIP_comb %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum), fill = Group)) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, aes(color = Group))+
  scale_colour_manual(values=c("#8B0046", "#468B00")) +
  scale_fill_manual(values=c("#8B0046", "#468B00")) + 
  xlim(1.2, 2) +
  ylab("") +
  xlab("") +
  #  ggtitle("Bleached Overlap") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17),
                     legend.position = "none")

D52_overlap_BvsP_down_plot <- D52_overlap_BvsP_down_VIP_comb %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum), fill = Group)) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, aes(color = Group))+
  scale_colour_manual(values=c("#8B0046", "#468B00")) +
  scale_fill_manual(values=c("#8B0046", "#468B00")) + 
  xlim(1.2, 2) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Bleached Overlap") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17), 
                     legend.text = element_text(size = 20),
                     legend.title = element_text(size = 22)) 


D52_All_VIP_plot <- D52_uniqueB_BvsP_up_plot +
  D52_overlap_BvsP_up_plot +
  D52_uniqueP_BvsP_up_plot +
  guide_area() +
  D52_overlap_BvsP_down_plot +
  D52_uniqueP_BvsP_down_plot + 
  plot_layout(ncol=3, heights=c(0.4,0.1), guides = "collect") 

ggsave(filename="../output/Metabolomics/Day52_VIP.png", plot=D52_All_VIP_plot, dpi=300, width=20, height=9, units="in")
```

```{r, echo=FALSE, out.width="100%", fig.cap="Day_52_VIP"}
knitr::include_graphics("../output/Metabolomics/Day52_VIP.png")
```


### Pathway Analysis

#### Day 52 Control vs Bleached
``` {r ORAPathD52CvBup, echo=TRUE, warning=FALSE, message=FALSE, results=FALSE}

#library(MetaboAnalystR)

# D52 Control vs Bleached UP Overrepresentation Analysis and  via Metaboanalyst 4.0

D52_CvsB_VIP_sig_up
D52_CvsB_VIP_sig_up$Metabolite[D52_CvsB_VIP_sig_up$Metabolite == "NADPH (SIM)"] <- "NADPH"

mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-c(D52_CvsB_VIP_sig_up$Metabolite)
mSet<-Setup.MapData(mSet, cmpd.vec)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "NG-dimethyl-L-arginine")
mSet<-GetCandidateList(mSet)
mSet<-SetCandidate(mSet, "NG-dimethyl-L-arginine", "Asymmetric dimethylarginine") #Replacing metabolite
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F)

D52_CvsB_up_dataset <- mSet$dataSet #extracting input dataframe
D52_CvsB_up_KEGG <- as.data.frame(D52_CvsB_up_dataset$map.table) #extracting KEGG terms

write.csv(D52_CvsB_up_KEGG, "../output/Metabolomics/MetaboAnalyst_Output/D52_Comparisons/D52_CvsB_up_KEGG.csv")

## PATHWAY ENRICHMENT

#Model org: C elegans
#ORA analysis: Fisher's exact test
#Pathway topology analysis: relative betweenness centrality (rbc)

mSet_Path<-SetKEGG.PathLib(mSet, "cel", "current") #comparing again C.elegans, review parameters
mSet_Path<-SetMetabolomeFilter(mSet_Path, F) 
mSet_Path<-CalculateOraScore(mSet_Path, "rbc", "hyperg") #review parameters

D52_CvsB_up_analset_KEGG <- mSet_Path$api
D52_CvsB_up_PATHWAY <- as.data.frame(D52_CvsB_up_analset_KEGG$ora.results) #extracting KEGG terms

write.csv(D52_CvsB_up_PATHWAY, "../output/Metabolomics/MetaboAnalyst_Output/D52_Comparisons/D52_CvsB_up_PATHWAY.csv")

D52_CvsB_up_PATHWAY 

# PLOTTING

D52_CvsB_up_PATHWAY <- tibble::rownames_to_column(D52_CvsB_up_PATHWAY, "Pathway")
names(D52_CvsB_up_PATHWAY)[5] <- "pvalue"
names(D52_CvsB_up_PATHWAY)[6] <- "neglogp"

logsig<- -log(0.05)

pointstolabel <- D52_CvsB_up_PATHWAY %>%
  filter(pvalue < 0.05)

D52_CvsB_up_PATHWAY_plot <- D52_CvsB_up_PATHWAY %>%
  mutate(color = ifelse(neglogp>logsig, "red", "grey")) %>%
  ggplot(aes(x=Impact, y=neglogp, size=Impact, color = "black", fill = color)) +
  geom_point(shape=21) +
  scale_color_identity() +
  scale_fill_identity() +
  geom_hline(yintercept=logsig, linetype="dashed", color = "red", size=2) +
  #    geom_text(data=subset(D52_CvsP_up_PATHWAY, neglogp>logsig), aes(Impact,neglogp,label=Pathway)) +
  ylab("-log(p)") +
  xlab("Pathway Impact") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17), 
                     legend.position="none") 

D52_CvsB_up_PATHWAY_plot


D52_CvsB_up_PATHWAY_plotting <- D52_CvsB_up_PATHWAY %>%
  mutate(color_plot = ifelse(neglogp>logsig, "red", "grey")) %>%
  mutate(plotname = ifelse(Pathway %in% pointstolabel$Pathway, Pathway, ""))

set.seed(42)
D52_CvsB_up_PATHWAY_plot <- ggplot(D52_CvsB_up_PATHWAY_plotting, 
                                   aes(x=Impact, y=neglogp, size=Impact, color = "black", fill = color_plot)) +
  geom_point(shape=21) +
  geom_hline(yintercept=logsig, linetype="dashed", color = "red", size=2) +
  geom_label_repel(aes(label = plotname), 
                   min.segment.length = 0,
                   arrow = arrow(length = unit(0.015, "npc")),
                   size = 3.5, 
                   fill = "white",
                   box.padding = unit(2.5, "lines"), 
                   point.padding = unit(0.8, "lines")) +
  scale_color_identity() +
  scale_fill_identity() +
  ylab("-log(p)") +
  xlab("Pathway Impact") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17), 
                     legend.position="none")

D52_CvsB_up_PATHWAY_plot


-log(0.015242) #4.183701
```


#### Day 52 Control vs Partial Mortality

``` {r ORAPathD52CvPup, echo=TRUE, warning=FALSE, message=FALSE, results=FALSE}
# D52 Control vs Partial Mortality UP Overrepresentation Analysis and  via Metaboanalyst 4.0

D52_CvsP_VIP_sig_up
D52_CvsP_VIP_sig_up$Metabolite[D52_CvsP_VIP_sig_up$Metabolite == "NADPH (SIM)"] <- "NADPH"

mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-c(D52_CvsP_VIP_sig_up$Metabolite)
mSet<-Setup.MapData(mSet, cmpd.vec)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)
# mSet<-PerformDetailMatch(mSet, "Acetyl Proline")
# mSet<-GetCandidateList(mSet)
mSet<-PerformDetailMatch(mSet, "NG-dimethyl-L-arginine")
mSet<-GetCandidateList(mSet)
mSet<-SetCandidate(mSet, "NG-dimethyl-L-arginine", "Asymmetric dimethylarginine")
mSet<-PerformDetailMatch(mSet, "Pyrroline-5-carboxylic acid")
mSet<-GetCandidateList(mSet)
#mSet<-SetCandidate(mSet, "Pyrroline-5-carboxylic acid", "1-Pyrroline-5-carboxylic acid") #this isnt working for some reson
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F)

D52_CvsP_up_dataset <- mSet$dataSet #extracting input dataframe
D52_CvsP_up_KEGG <- as.data.frame(D52_CvsP_up_dataset$map.table) #extracting KEGG terms

write.csv(D52_CvsP_up_KEGG, "../output/Metabolomics/MetaboAnalyst_Output/D52_Comparisons/D52_CvsP_up_KEGG.csv")

## PATHWAY ENRICHMENT

#Model org: C elegans
#ORA analysis: Fisher's exact test
#Pathway topology analysis: relative betweenness centrality (rbc)

mSet_Path<-SetKEGG.PathLib(mSet, "cel", "current") #comparing again C.elegans, review parameters
mSet_Path<-SetMetabolomeFilter(mSet_Path, F) 
mSet_Path<-CalculateOraScore(mSet_Path, "rbc", "hyperg") #review parameters

D52_CvsP_up_analset_KEGG <- mSet_Path$api
D52_CvsP_up_PATHWAY <- as.data.frame(D52_CvsP_up_analset_KEGG$ora.results) #extracting KEGG terms

write.csv(D52_CvsP_up_PATHWAY, "../output/Metabolomics/MetaboAnalyst_Output/D52_Comparisons/D52_CvsP_up_PATHWAY.csv")

D52_CvsP_up_PATHWAY 

# Plotting

D52_CvsP_up_PATHWAY <- tibble::rownames_to_column(D52_CvsP_up_PATHWAY, "Pathway")
names(D52_CvsP_up_PATHWAY)[5] <- "pvalue"
names(D52_CvsP_up_PATHWAY)[6] <- "neglogp"

logsig<- -log(0.05)

pointstolabel <- D52_CvsP_up_PATHWAY %>%
  filter(pvalue < 0.05)

D52_CvsP_up_PATHWAY_plotting <- D52_CvsP_up_PATHWAY %>%
  mutate(color_plot = ifelse(neglogp>logsig, "red", "grey")) %>%
  mutate(plotname = ifelse(Pathway %in% pointstolabel$Pathway, Pathway, ""))

set.seed(42)
D52_CvsP_up_PATHWAY_plot <- ggplot(D52_CvsP_up_PATHWAY_plotting, 
                                   aes(x=Impact, y=neglogp, size=Impact, color = "black", fill = color_plot)) +
  geom_point(shape=21) +
  geom_hline(yintercept=logsig, linetype="dashed", color = "red", size=2) +
  geom_label_repel(aes(label = plotname), 
                   min.segment.length = 0,
                   arrow = arrow(length = unit(0.015, "npc")),
                   size = 3.5, 
                   fill = "white",
                   box.padding = unit(2.5, "lines"), 
                   point.padding = unit(0.8, "lines")) +
  scale_color_identity() +
  scale_fill_identity() +
  ylab("-log(p)") +
  xlab("Pathway Impact") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17), 
                     legend.position="none")

D52_CvsP_up_PATHWAY_plot

# https://ggrepel.slowkow.com/articles/examples.html

```


# D37 ANALYSIS

## PLS-DA on Timepoint 37 between each pair

### Control vs Bleached

``` {r PLSDACvsB37, echo=TRUE, warning=FALSE, message=FALSE}

#Control vs Bleached

norm.data_2 <- column_to_rownames(norm.data, 'Sample.ID')

D37_CvsB <- norm.data_2 %>%
  filter(Time == "Day37") %>%
  filter(Treatment != "Mortality_Hot")

D37_CvsB_clean <- na.omit(D37_CvsB)

#assigning datasets 
X <- D37_CvsB[4:183]
Y <- as.factor(D37_CvsB$Treatment) 
#summary(Y) ## class summary
#summary(X)
#dim(X) ## number of samples and features
#length(Y) ## length of class membership factor = number of samples


#PLSDA without variable selection
MyResult.plsda <- plsda(X, Y, ncomp = 2) # 1 Run the method
#plotIndiv(MyResult.plsda)    # 2 Plot the samples

#plotVar(MyResult.plsda, cutoff = 0.7)    

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,ellipse = TRUE, title="Day 37 - Bleached vs Control")

MyResult.plsda2 <- plsda(X,Y, ncomp=6) #number of components is #classes-1
#selectVar(MyResult.plsda2, comp=1)$value

#plotLoadings(MyResult.plsda2, comp = 1, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1
#plotLoadings(MyResult.plsda2, comp = 2, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1

comp1.select.metabolites.all <- data.frame(selectVar(MyResult.plsda2, comp = 1)$value)

# component validation 
set.seed(200) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 5, #figure out why only folds =5 works
                     progressBar = FALSE, auc = TRUE, nrepeat = 10) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

#ROC Curve
auc.plsda = auroc(MyResult.plsda2, roc.comp = 2)

#VIP Extraction
D37_CvsB_VIP <- PLSDA.VIP(MyResult.plsda2)
D37_CvsB_VIP_DF <- as.data.frame(D52_CvsB_VIP[["tab"]])

```


### Control vs Partial Mortality

``` {r PLSDACvsP37, echo=TRUE, warning=FALSE, message=FALSE}
D37_CvsP <- norm.data_2 %>%
  filter(Time == "Day37") %>%
  filter(Treatment != "Bleached_Hot")

D37_CvsP_clean <- na.omit(D37_CvsP)

#assigning datasets 
X <- D37_CvsP[4:183]
Y <- as.factor(D37_CvsP$Treatment) 
#summary(Y) ## class summary
#summary(X)
#dim(X) ## number of samples and features
#length(Y) ## length of class membership factor = number of samples

#PLSDA without variable selection
MyResult.plsda <- plsda(X, Y, ncomp = 2) # 1 Run the method
#plotIndiv(MyResult.plsda)    # 2 Plot the samples

#plotVar(MyResult.plsda, cutoff = 0.7)    

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,ellipse = TRUE, title="Day 37 - Control vs Partial Mortality")

MyResult.plsda2 <- plsda(X,Y, ncomp=6) #number of components is #classes-1
#selectVar(MyResult.plsda2, comp=1)$value

#plotLoadings(MyResult.plsda2, comp = 1, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1
#plotLoadings(MyResult.plsda2, comp = 2, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1

comp1.select.metabolites.all <- data.frame(selectVar(MyResult.plsda2, comp = 1)$value)

# component validation 
set.seed(200) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 5, #figure out why only folds =5 works
                     progressBar = FALSE, auc = TRUE, nrepeat = 10) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

#ROC Curve
auc.plsda = auroc(MyResult.plsda2, roc.comp = 2)

#VIP Extraction
D37_CvsP_VIP <- PLSDA.VIP(MyResult.plsda2)
D37_CvsP_VIP_DF <- as.data.frame(D37_CvsP_VIP[["tab"]])

```


### Bleached vs Partial Mortality

``` {r PLSDABvsP37, echo=TRUE, warning=FALSE, message=FALSE}

D37_BvsP <- norm.data_2 %>%
  filter(Time == "Day37") %>%
  filter(Treatment != "Control_Ambient")

D37_BvsP_clean <- na.omit(D37_BvsP)

#assigning datasets 
X <- D37_BvsP[4:183]
Y <- as.factor(D37_BvsP$Treatment) 
#summary(Y) ## class summary
#summary(X)
#dim(X) ## number of samples and features
#length(Y) ## length of class membership factor = number of samples


#PLSDA without variable selection
MyResult.plsda <- plsda(X, Y, ncomp = 2) # 1 Run the method
#plotIndiv(MyResult.plsda)    # 2 Plot the samples

#plotVar(MyResult.plsda, cutoff = 0.7)    

plotIndiv(MyResult.plsda, ind.names = FALSE, legend=TRUE,ellipse = TRUE, title="Day 37 - Bleached vs Partial Mortality")

MyResult.plsda2 <- plsda(X,Y, ncomp=6) #number of components is #classes-1
#selectVar(MyResult.plsda2, comp=1)$value

#plotLoadings(MyResult.plsda2, comp = 1, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1
#plotLoadings(MyResult.plsda2, comp = 2, contrib = 'max', method = 'median', ndisplay = 50) #top 50 metabolites contributing to variation on component 1

comp1.select.metabolites.all <- data.frame(selectVar(MyResult.plsda2, comp = 1)$value)

# component validation 
set.seed(200) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 5, #figure out why only folds =5 works
                     progressBar = FALSE, auc = TRUE, nrepeat = 10) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

#ROC Curve
auc.plsda = auroc(MyResult.plsda2, roc.comp = 2)

#VIP Extraction
D37_BvsP_VIP <- PLSDA.VIP(MyResult.plsda2)
D37_BvsP_VIP_DF <- as.data.frame(D37_BvsP_VIP[["tab"]])

```


### VIP Comparisons for D37 Timepoint

``` {r VIP37, echo=TRUE, warning=FALSE, message=FALSE}

####### Overlaps for VIPs >1 #########

# Converting row names to column
D37_CvsB_VIP_table <- rownames_to_column(D37_CvsB_VIP_DF, var = "Metabolite")
D37_CvsP_VIP_table <- rownames_to_column(D37_CvsP_VIP_DF, var = "Metabolite")
D37_BvsP_VIP_table <- rownames_to_column(D37_BvsP_VIP_DF, var = "Metabolite")

# Filtering for VIP > 1
D37_CvsB_VIP_1 <- D37_CvsB_VIP_table %>% 
  filter(VIP >= 1)

D37_CvsP_VIP_1 <- D37_CvsP_VIP_table %>% 
  filter(VIP >= 1)

D37_BvsP_VIP_1 <- D37_BvsP_VIP_table %>% 
  filter(VIP >= 1)

D37_CvsB_VIP_1 %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_point() +
  ylab("Metabolite") +
  xlab("VIP Score") +
  ggtitle("Control vs Bleached") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

D37_CvsP_VIP_1 %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_point() +
  ylab("Metabolite") +
  xlab("VIP Score") +
  ggtitle("Control vs Partial Mortality") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

D37_BvsP_VIP_1 %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_point() +
  ylab("Metabolite") +
  xlab("VIP Score") +
  ggtitle("Bleached vs Partial Mortality") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



### Compare for Venn Diagram (https://github.com/gaospecial/ggVennDiagram)

#Open this part in Repo

# nrow(D52_CvsB_VIP_1)
# nrow(D52_CvsP_VIP_1)
# 
# library("VennDiagram")
# 
# venn.x <- list(
#   D52_CvsB = sample(D52_CvsB_VIP_1$Metabolite), 
#   D52_CvsP = sample(D52_CvsP_VIP_1$Metabolite)
# )
# 
# myCol <- c("#8B0046", "#468B00")
# 
# venn.diagram(venn.x,
#              filename = 'output/Metabolomics/D52_BvsP_venn.png',
#              output = TRUE,
#              height = 480, 
#              width = 480,
#              resolution = 300,
#              category.names = c("Bleached", "Partial Mortality"),
#              lwd = 2,
#              col = c("#8B0046", "#468B00"),
#              fill = c(alpha("#8B0046",0.3), alpha("#468B00", 0.3)),
#              cex = 0.5,
#              fontfamily = "sans",
#              fontface = "bold",
#              cat.cex = 0.5,
#              cat.default.pos = "outer",
#              cat.pos = c(12, -12),
#              cat.dist = c(-0.45, -0.45),
#              cat.fontfamily = "sans",
#              cat.fontface = "bold"
#              )
```

### T-test validation: Control vs Bleached
``` {r t-test validationCB37, echo=TRUE, warning=FALSE, message=FALSE}

# Control vs Bleached

## Gather data frame and group by metabolite 
D37_CvB_gather <- D37_CvsB_clean %>% gather(key = "Metabolite", value = "Count", -Fragment.ID, -Time, -Treatment) %>%
  dplyr::select("Treatment", "Metabolite", "Count")

D37_CvB_gather$Treatment <- as.factor(D37_CvB_gather$Treatment)
D37_CvB_gather <- dplyr::group_by(D37_CvB_gather, Metabolite)

# Select metabolites only present from the PLS-DA with VIPs >1
D37_CvB_VIP_Select <- subset(D37_CvB_gather, D37_CvB_gather$Metabolite %in% D37_CvsB_VIP_1$Metabolite)

#Looped t-test for all metabolites with VIPs >1, estimate 1 = Bleached_Hot, estimate 2 = Control
D37_CvB_t.test <-do(D37_CvB_VIP_Select, tidy(t.test(.$Count ~ .$Treatment,
                                                    alternative = "two.sided",
                                                    mu = 0,
                                                    paired = FALSE,
                                                    var.equal = FALSE,
                                                    conf.level = 0.95
)))

#adjust p value for the number of comparisons
D37_CvB_t.test$p.adj<-p.adjust(D37_CvB_t.test$p.value, method=c("fdr"), n=length(D37_CvsB_VIP_1$Metabolite)) #length = number of metabolites tested, false discovery rate method

# Filter for significantly different metabolites (p < 0.05)
D37_CvB_t.test.fdr <- D37_CvB_t.test %>% filter(p.adj <= 0.05)
length(D37_CvB_t.test.fdr$Metabolite) #3

# # Filter for significantly different metabolites (p < 0.05)
D37_CvB_t.test.sig <- D37_CvB_t.test %>% filter(p.value <= 0.05)
length(D37_CvB_t.test.sig$Metabolite) #went from 18 to 3 with p value adjustment

# Filter for metabolites accumulated compared to control
D37_CvB_t.test.sig.up <- D37_CvB_t.test.fdr %>% filter(estimate1 > estimate2)

# Filter for metabolites depleted compared to control
D37_CvB_t.test.sig.down <- D37_CvB_t.test.fdr %>% filter(estimate1 < estimate2)

```


## T-test validation: Control vs Partial Mortality
``` {r t-test_validationCP37, echo=TRUE, warning=FALSE, message=FALSE}

# Control vs Bleached

## Gather data frame and group by metabolite 
D37_CvsP_gather <- D37_CvsP_clean %>% gather(key = "Metabolite", value = "Count", -Fragment.ID, -Time, -Treatment) %>%
  dplyr::select("Treatment", "Metabolite", "Count")

D37_CvsP_gather$Treatment <- as.factor(D37_CvsP_gather$Treatment)
D37_CvsP_gather <- dplyr::group_by(D37_CvsP_gather, Metabolite)

# Select metabolites only present from the PLS-DA with VIPs >1
D37_CvsP_VIP_Select <- subset(D37_CvsP_gather, D37_CvsP_gather$Metabolite %in% D37_CvsP_VIP_1$Metabolite)

#Looped t-test for all metabolites with VIPs >1, estimate 1 = Control_Ambient, estimate 2 = Bleached_Hot
D37_CvsP_t.test <-do(D37_CvsP_VIP_Select, tidy(t.test(.$Count ~ .$Treatment,
                                                      alternative = "two.sided",
                                                      mu = 0,
                                                      paired = FALSE,
                                                      var.equal = FALSE,
                                                      conf.level = 0.95
)))

#adjust p value for the number of comparisons
D37_CvsP_t.test$p.adj<-p.adjust(D37_CvsP_t.test$p.value, method=c("fdr"), n=length(D37_CvsP_VIP_1$Metabolite)) #length = number of metabolites tested, false discovery rate method

# Filter for significantly different metabolites (p < 0.05)
D37_CvsP_t.test.fdr <- D37_CvsP_t.test %>% filter(p.adj <= 0.05)
length(D37_CvsP_t.test.fdr$Metabolite) #12

# Filter for significantly different metabolites (p < 0.05)
D37_CvsP_t.test.sig <- D37_CvsP_t.test %>% filter(p.value <= 0.05)
length(D37_CvsP_t.test.sig$Metabolite) #went from 26 to 12 with p value adjustment

# Filter for metabolites accumulated compared to control
D37_CvsP_t.test.sig.up <- D37_CvsP_t.test.fdr %>% filter(estimate1 < estimate2)

# Filter for metabolites depleted compared to control
D37_CvsP_t.test.sig.down <- D37_CvsP_t.test.fdr %>% filter(estimate1 > estimate2)

```


## T-test validation: Bleached vs Partial Mortality
``` {r t-test_validationBP37, echo=TRUE, warning=FALSE, message=FALSE}

# Bleach vs Partial Mortality

## Gather data frame and group by metabolite 
D37_BvsP_gather <- D37_BvsP_clean %>% gather(key = "Metabolite", value = "Count", -Fragment.ID, -Time, -Treatment) %>%
  dplyr::select("Treatment", "Metabolite", "Count")

D37_BvsP_gather$Treatment <- as.factor(D37_BvsP_gather$Treatment)
D37_BvsP_gather <- dplyr::group_by(D37_BvsP_gather, Metabolite)

# Select metabolites only present from the PLS-DA with VIPs >1
D37_BvsP_VIP_Select <- subset(D37_BvsP_gather, D37_BvsP_gather$Metabolite %in% D37_BvsP_VIP_1$Metabolite)

#Looped t-test for all metabolites with VIPs >1, estimate 1 = Bleached_Hot, estimate 2 = Mortality_Hot
D37_BvsP_t.test <-do(D37_BvsP_VIP_Select, tidy(t.test(.$Count ~ .$Treatment,
                                                      alternative = "two.sided",
                                                      mu = 0,
                                                      paired = FALSE,
                                                      var.equal = FALSE,
                                                      conf.level = 0.95
)))

#adjust p value for the number of comparisons
D37_BvsP_t.test$p.adj<-p.adjust(D37_BvsP_t.test$p.value, method=c("fdr"), n=length(D37_BvsP_VIP_1$Metabolite)) #length = number of metabolites tested, false discovery rate method

# Filter for significantly different metabolites (p < 0.05)
D37_BvsP_t.test.fdr <- D37_BvsP_t.test %>% filter(p.adj <= 0.05)
length(D37_BvsP_t.test.fdr$Metabolite) #0

# Filter for significantly different metabolites (p < 0.05)
D37_CvsP_t.test.sig <- D37_BvsP_t.test %>% filter(p.value <= 0.05)
length(D37_CvsP_t.test.sig$Metabolite) #went from 10 to 0 with p value adjustment

# Filter for metabolites accumulated compared to control
D37_CvsP_t.test.sig.up <- D37_CvsP_t.test.fdr %>% filter(estimate1 < estimate2)

# Filter for metabolites depleted compared to control
D37_CvsP_t.test.sig.down <- D37_CvsP_t.test.fdr %>% filter(estimate1 > estimate2)

```


### VIP plotting after t-test validation

``` {r VIP-test37, echo=TRUE, warning=FALSE, message=FALSE}

# Selecting metabolites that were validated by the t-test for each up and down accumulated
D37_CvsB_VIP_sig_up <- subset(D37_CvsB_VIP_1, D37_CvsB_VIP_1$Metabolite %in% D37_CvB_t.test.sig.up$Metabolite)
D37_CvsB_VIP_sig_down <- subset(D37_CvsB_VIP_1, D37_CvsB_VIP_1$Metabolite %in% D37_CvB_t.test.sig.down$Metabolite)

D37_CvsP_VIP_sig_up <- subset(D37_CvsP_VIP_1, D37_CvsP_VIP_1$Metabolite %in% D37_CvsP_t.test.sig.up$Metabolite)
D37_CvsP_VIP_sig_down <- subset(D37_CvsP_VIP_1, D37_CvsP_VIP_1$Metabolite %in% D37_CvsP_t.test.sig.down$Metabolite)

write.csv(D37_CvsB_VIP_sig_up, "../output/Metabolomics/D37_CvsB_VIP_sig_up.csv")
write.csv(D37_CvsB_VIP_sig_down, "../output/Metabolomics/D37_CvsB_VIP_sig_down.csv")
write.csv(D37_CvsP_VIP_sig_up, "../output/Metabolomics/D37_CvsP_VIP_sig_up.csv")
write.csv(D37_CvsP_VIP_sig_down, "../output/Metabolomics/D37_CvsP_VIP_sig_down.csv")

# Plotting

D37_CvsB_VIP_UP <- D37_CvsB_VIP_sig_up %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point() +
  xlim(1, 2.3) +
  ylab("Accumulated") +
  xlab("") +
  ggtitle("Bleached") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

D37_CvsB_VIP_DOWN <- D37_CvsB_VIP_sig_down %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point() +
  xlim(1, 2.3) +
  ylab("Depleted") +
  xlab("VIP Score") +
  #  ggtitle("Bleached: Depleted") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

D37_CvsP_VIP_UP <- D37_CvsP_VIP_sig_up %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point() +
  xlim(1, 2.3) +
  ylab("") +
  xlab("") +
  ggtitle("Partial Mortality") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

D37_CvsP_VIP_DOWN <- D37_CvsP_VIP_sig_down %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point() +
  xlim(1, 2.3) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Partial Mortality: Depleted") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

D37_All_VIP_plot <- plot_grid(D37_CvsB_VIP_UP, D37_CvsP_VIP_UP, D37_CvsB_VIP_DOWN, D37_CvsP_VIP_DOWN, 
                              align = "v", 
                              ncol = 2, 
                              rel_heights = c(0.95, 0.35), 
                              labels = c("A", "B", "C", "D"))

D37_All_VIP_plot
#ggsave(filename="../output/Metabolomics/Day37_VIP.pdf", plot=D37_All_VIP_plot, dpi=300, width=7, height=9, units="in")

```


### VIP plotting after t-test validation and comparing overlapping metabolites

``` {r Plotting37, echo=TRUE, warning=FALSE, message=FALSE}

# Up regulated metabolites that overlap between B and P

D37_overlap_BvsP_up_1 <- as.data.frame(intersect(D37_CvsB_VIP_sig_up$Metabolite, D37_CvsP_VIP_sig_up$Metabolite))
names(D37_overlap_BvsP_up_1)[1] <- 'Metabolite'
D37_overlap_BvsP_up_2 <- merge(D37_overlap_BvsP_up_1, D37_CvsB_VIP_sig_up, by="Metabolite")
D37_overlap_BvsP_up_VIP <- merge(D37_overlap_BvsP_up_2, D37_CvsP_VIP_sig_up, by="Metabolite")
names(D37_overlap_BvsP_up_VIP)[2] <- 'Bleached'
names(D37_overlap_BvsP_up_VIP)[3] <- 'Partial Mortality'

D37_overlap_BvsP_up_VIP_comb <- gather(D37_overlap_BvsP_up_VIP, key = "Group", value = "VIP", 'Bleached', 'Partial Mortality')

# Up regulated metabolites that are unique to B

D37_uniqueB_BvsP_up <- as.data.frame(setdiff(D37_CvsB_VIP_sig_up$Metabolite, D37_CvsP_VIP_sig_up$Metabolite))
names(D37_uniqueB_BvsP_up)[1] <- 'Metabolite'
D37_uniqueB_BvsP_up_VIP <- merge(D37_uniqueB_BvsP_up, D37_CvsB_VIP_sig_up, by="Metabolite")

# Up regulated metabolites that are unique to P

D37_uniqueP_BvsP_up <- as.data.frame(setdiff(D37_CvsP_VIP_sig_up$Metabolite, D37_CvsB_VIP_sig_up$Metabolite))
names(D37_uniqueP_BvsP_up)[1] <- 'Metabolite'
D37_uniqueP_BvsP_up_VIP <- merge(D37_uniqueP_BvsP_up, D37_CvsP_VIP_sig_up, by="Metabolite")

# Down regulated metabolites that overlap between B and P

D37_overlap_BvsP_down_1 <- as.data.frame(intersect(D37_CvsB_VIP_sig_down$Metabolite, D37_CvsP_VIP_sig_down$Metabolite))
names(D37_overlap_BvsP_down_1)[1] <- 'Metabolite'
D37_overlap_BvsP_down_2 <- merge(D37_overlap_BvsP_down_1, D37_CvsB_VIP_sig_down, by="Metabolite")
D37_overlap_BvsP_down_VIP <- merge(D37_overlap_BvsP_down_2, D37_CvsP_VIP_sig_down, by="Metabolite")
names(D37_overlap_BvsP_down_VIP)[2] <- 'Bleached'
names(D37_overlap_BvsP_down_VIP)[3] <- 'Partial Mortality'

D37_overlap_BvsP_down_VIP_comb <- gather(D37_overlap_BvsP_down_VIP, key = "Group", value = "VIP", 'Bleached', 'Partial Mortality')


# Down regulated metabolites that are unique to B 
# There are none 

# Down regulated metabolites that are unique to P

D37_uniqueP_BvsP_down <- as.data.frame(setdiff(D37_CvsP_VIP_sig_down$Metabolite, D37_CvsB_VIP_sig_down$Metabolite))
names(D37_uniqueP_BvsP_down)[1] <- 'Metabolite'
D37_uniqueP_BvsP_down_VIP <- merge(D37_uniqueP_BvsP_down, D37_CvsP_VIP_sig_down, by="Metabolite")


# Plotting

D37_uniqueB_BvsP_up_plot <- D37_uniqueB_BvsP_up_VIP %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, color = "#8B0046") +
  xlim(1, 2.3) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Bleached Unique") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17))

D37_uniqueP_BvsP_up_plot <- D37_uniqueP_BvsP_up_VIP %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, color = "#468B00") +
  xlim(1, 2.3) +
  ylab("") +
  xlab("") +
  #  ggtitle("Partial Mortality Unique") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17))

D37_uniqueP_BvsP_down_plot <- D37_uniqueP_BvsP_down_VIP %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, color = "#468B00") +
  xlim(1, 2.3) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Partial Mortality: Depleted") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17))

D37_overlap_BvsP_up_plot <- D37_overlap_BvsP_up_VIP_comb %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum), fill = Group)) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, aes(color = Group))+
  scale_colour_manual(values=c("#8B0046", "#468B00")) +
  scale_fill_manual(values=c("#8B0046", "#468B00")) + 
  xlim(1, 2.3) +
  ylab("") +
  xlab("") +
  #  ggtitle("Bleached Overlap") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17),
                     legend.position = "none")

D37_overlap_BvsP_down_plot <- D37_overlap_BvsP_down_VIP_comb %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum), fill = Group)) +
  geom_hline(aes(yintercept = Metabolite), linetype = "dotted", color = "grey") +
  geom_point(size = 7, aes(color = Group))+
  scale_colour_manual(values=c("#8B0046", "#468B00")) +
  scale_fill_manual(values=c("#8B0046", "#468B00")) + 
  xlim(1, 2.3) +
  ylab("") +
  xlab("VIP Score") +
  #  ggtitle("Bleached Overlap") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(size=17, color="black"), 
                     title = element_text(size=17, face="bold"), 
                     axis.title = element_text(size=17), 
                     legend.text = element_text(size = 20),
                     legend.title = element_text(size = 22)) 


D37_All_VIP_plot <- D37_uniqueB_BvsP_up_plot +
  D37_overlap_BvsP_up_plot +
  #                     D37_uniqueP_BvsP_up_plot +
  guide_area() +
  D37_overlap_BvsP_down_plot +
  #                     D37_uniqueP_BvsP_down_plot + 
  plot_layout(ncol=2, heights=c(0.4,0.1), guides = "collect") 

ggsave(filename="../output/Metabolomics/Day37_VIP.png", plot=D37_All_VIP_plot, dpi=300, width=20, height=9, units="in")
```

