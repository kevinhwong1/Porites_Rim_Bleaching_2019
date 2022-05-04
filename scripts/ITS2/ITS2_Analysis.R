# PJB ITS2 Analysis Post Symportal
# Author: Kevin Wong
# Date modified: 20220310

# Load Packages
library("dplyr")
library("tidyverse") 
library("ggplot2")
library("janitor")
library("funrar")

# Load data
its2 <- read.delim("output/ITS2/its2_type_profiles/4_PJB_analysis_20211227T125938.profiles.absolute.abund_and_meta.txt", sep = "\t")
meta <- read.csv("data/Metadata/Frag_Sample_Info.csv")

# Cleaning data frame
its2_2 <- its2[-c(1:5, 52:53), ] #removing unncessary rows
its2_2 <- its2_2[-c(1)] # Removing first column

its2_3 <- its2_2 %>%
  row_to_names(row_number = 1) # Making a new column header

colnames(its2_3)[1] <- "Fragment_ID" # renaming column 2 

its2_4 <- its2_3[,-1] #removing the first row
rownames(its2_4) <- its2_3[,1] #making the fragment ids the row names 

its2_4 <- as.numeric(its2_4[c(1:9),])

# Making the relative abundance dataframe 
its2_5 <- data.matrix(sapply(its2_4, as.numeric))
rel_its2 <- make_relative(its2_5)
rel_its2_df <- as.data.frame(rel_its2)

# adding Fragment ID column

rel_its2_df <- cbind(rel_its2_df, its2_3$Fragment_ID) 
colnames(rel_its2_df)[10] <- "Fragment_ID" # renaming column

# Gathering data 
rel_its2_df2 <- rel_its2_df %>%
  gather(key=Clade, value=relabund, 1:9)

# Merging metadata
rel_its2_meta <- merge(rel_its2_df2, meta, by = "Fragment_ID")
rel_its2_meta2 <- rel_its2_meta %>% select(-ColDay, -Day_2, -Fragment_ID)

# Plotting

ps.melt_sum <- ps.melt %>%
  group_by(Sample,Group, Day, Phylum)

rel_abun_bar <- ggplot(rel_its2_meta2, aes(x = Colony_ID, y = relabund, fill = Clade)) + 
  geom_bar(stat = "identity", aes(fill=Clade)) + 
  labs(x="", y="Relative Abundance") +
  scale_fill_brewer(palette="Spectral") +
  facet_wrap(~Group+Day, scales= "free_x", nrow=1) +
  theme_bw() + 
  theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank(), axis.text.x.bottom = element_text(angle = -90))

ggsave(filename="output/ITS2/Relative_Abundance_Barplot.png", plot=rel_abun_bar, dpi=300, width=10, height=5, units="in")


sym <- read.csv("output/Physiology/Zoox.Calc.csv")

sym_abun_bar <- ggplot(sym, aes(x = Fragment.ID, y = Cells.cm2.x6)) + 
  geom_bar(stat = "identity" ) + 
  labs(x="", y="Cellsx10^6 cm2 ") +
  scale_fill_brewer(palette="Spectral") +
  facet_wrap(~Group+Day, scales= "free_x", nrow=1) +
  theme_bw() + 
  theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.text.x.bottom = element_text(angle = -90))


