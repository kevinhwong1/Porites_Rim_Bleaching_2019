## Creating metadata files for 16S pipeline
## Emma Strand created 1/7/2022

library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)

## 1. Sample manifest file
## filelist.csv created during pipeline in andromeda and scp to desktop to work with in R

file_names <- read.csv("output/16S/filelist.csv", header = FALSE) %>%
  select(V2) %>% ## reading in filelist as dataframe and only the second column
  dplyr::rename(`absolute-filepath` = V2) # renaming to match the verbiage of qiime2

sample_manifest <- file_names # creating a new df based on the original file_names df
sample_manifest$path <- "/data/putnamlab/kevin_wong1/PJB_16S/raw_data/" #adding the absolute file path

sample_manifest <- sample_manifest %>% unite(`absolute-filepath`, path, `absolute-filepath`, sep = "") %>% # merging the two columns to complete the file path
  mutate(direction = case_when(grepl("R1", `absolute-filepath`) ~ "forward",
                               grepl("R2", `absolute-filepath`) ~ "reverse")) # creating a new column to state whether forward or reverse based on the R value in the sequence title name

sample_manifest$`sample-id` <- substr(sample_manifest$`absolute-filepath`, 46, 51) # creating a new column based on the sample id value

sample_manifest <- sample_manifest[, c(3, 1, 2)] # reordering the columns

sample_manifest <- sample_manifest %>% spread(direction, `absolute-filepath`) %>%
  dplyr::rename(`forward-absolute-filepath` = forward) %>%
  dplyr::rename(`reverse-absolute-filepath` = reverse)

write.table(sample_manifest, "output/16S/sample_manifest.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## return to terminal to secure copy paste the sample manifest file to bluewaves/andromeda folders

## 2. Sample metadata file
metadata <- read.csv("data/Molecular/16S/16S_Seq_Metadata.csv", header = TRUE) %>%
  select(-Plate, -Well, -Vial) # removing 2 columns that are not needed for this metadata sheet

categories <- c("#q2:types", "categorical", "categorical", "categorical", "categorical") # QIIME2 needs each column to be specified

metadata <- rbind(metadata, categories)
metadata <- metadata[c(46,1:45),]

write.table(metadata, "output/16S/metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)


######## denoising statistics ######

## comparing 3 options for denoising parameters

denoise_260230 <- read.table("output/16S/processed_data/denoising-stats_260230.tsv", sep="\t", header = TRUE)
denoise_270240 <- read.table("output/16S/processed_data/denoising-stats_270240.tsv", sep="\t", header = TRUE)
denoise_280240 <- read.table("output/16S/processed_data/denoising-stats_280240.tsv", sep="\t", header = TRUE)

denoise_260230$parameter <- "260forward_230reverse"
denoise_270240$parameter <- "270forward_240reverse"
denoise_280240$parameter <- "280forward_240reverse"

denoising.stats <- union(denoise_260230, denoise_270240) %>% union(denoise_280240)

denoise.reads <- denoising.stats[, c(1,2,3,5,6,8,10)] # reordering columns to make it easier to plot
denoise.percent <- denoising.stats[, c(1,4,7,9,10)] # reordering columns to make it easier to plot

denoise.reads <- denoise.reads %>% gather(statistic, value, 2:6) # aggregates the three variables we are interested in to make it easier to plot
denoise.percent <- denoise.percent %>% gather(statistic, value, 2:4) # aggregates the three variables we are interested in to make it easier to plot

denoise.reads$statistic <- factor(denoise.reads$statistic, levels=c("input","filtered","denoised","merged","non.chimeric"))
denoise.percent$statistic <- factor(denoise.percent$statistic, levels=c("percentage.of.input.passed.filter", "percentage.of.input.merged",
                                                                        "percentage.of.input.non.chimeric"))

percent <- ggplot(data = denoise.percent, aes(x = parameter, y = value, group = parameter, color = parameter)) +
  theme_classic() + geom_boxplot() +
  facet_grid(~statistic, scales = "free") +
  theme(legend.position = "none") +
  ylab("% reads") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.2, hjust = 1.3)) #Set the text angle

reads <- ggplot(data = denoise.reads, aes(x = parameter, y = value, group = parameter, color = parameter)) +
  theme_classic() + geom_boxplot() +
  facet_grid(~statistic, scales = "free") +
  theme(legend.position = "none") +
  ylab("# reads") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.2, hjust = 1.3)) #Set the text angle

percent
reads

ggsave(file="output/16S/processed_data/denoising-percent.png", percent, width = 11, height = 6, units = c("in"))
ggsave(file="output/16S/processed_data/denoising-reads.png", reads, width = 11, height = 6, units = c("in"))
