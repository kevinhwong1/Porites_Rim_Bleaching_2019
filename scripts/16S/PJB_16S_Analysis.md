# PJB 16S Analysis

Following [E. Strad pipeline](https://emmastrand.github.io/EmmaStrand_Notebook/E5-16S-Analysis/)

Primers used:
* 515F: GTGYCAGCMGCCGCGGTAA
* 806RB: GGACTACNVGGGTWTCTAAT

Modules needed:
* FastQC/0.11.9-Java-11
* MultiQC/1.9-intel-2020a-Python-3.8.2
* QIIME2/2021.4

# Copying raw files to my directory

```bash
mkdir PJB_16S
cd PJB_16S
mkdir raw_data
cd raw_data
scp /data/putnamlab/KITT/hputnam/20211210_AmpliconSeq/AmpSeq/KW_PJB_16S/W* .

```

# Quality Control with FastQC

`nano fastqc.sh`

```bash
#!/bin/bash
#SBATCH --job-name="fastqc_16S"
#SBATCH -t 24:00:00
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/PJB_16S
#SBATCH --exclusive

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in ./raw_data/*fastq.gz
do
fastqc $file --outdir ./fastqc        
done

multiqc --interactive ./fastqc

mv multiqc_report.html ./multiqc_data/PJB_16S_raw_qc_multiqc_report.html #renames file

```

```bash
scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/multiqc_data/PJB_16S_raw_qc_multiqc_report.html  /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/.
```

# Create Metadata files

```bash
cd /data/putnamlab/kevin_wong1/PJB_16S
mkdir metadata
```

#### Sample manifest file

Print all file names:

```bash
cd /data/putnamlab/kevin_wong1/PJB_16S
find raw_data -type f -print | sed 's_/_,_g' > metadata/filelist.csv
```

Export to local computer:

```bash
scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/metadata/filelist.csv  /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/.
```

Run the Sample manifest file section in 16S_metadata.R file and then return to the following steps.

R script:
```bash
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
```

Secure copy paste the sample manifest file in a terminal window outside of andromeda.

```bash
scp /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/sample_manifest.txt kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/metadata
```

```bash
scp /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/metadata.txt kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/metadata
```

# QIIME2

#### Sample Input

```bash
mkdir scripts
cd scripts/
nano import.sh
```

```bash
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/PJB_16S

source /usr/share/Modules/init/sh # load the module function
module load QIIME2/2021.8

#### METADATA FILES ####

# Metadata path
METADATA="metadata/metadata.txt"

# Sample manifest path
MANIFEST="metadata/sample_manifest.txt"

#########################

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path PJB_16S-paired-end-sequences.qza
```

#### QIIME2 denoising

Parameters chosen:

* --p-trunc-len-f: 19 (forward is 19 bp long)
* --p-trunc-len-r: 20 (reverse is 20 bp long)
* p-trim-left-f: try 260, 270, and 280
* p-trim-left-r: try 230, 240

`mkdir processed_data`

1. forward 260, reverse 230 (most conservative)

`mkdir denoise_260230`

`nano denoise_260230.sh`

```bash
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/PJB_16S/processed_data/denoise_260230

source /usr/share/Modules/init/sh # load the module function
module load QIIME2/2021.8

#### METADATA FILES ####

# Metadata path
METADATA="../../metadata/metadata.txt"

# Sample manifest path
MANIFEST="../../metadata/sample_manifest.txt"

#########################

qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs ../../PJB_16S-paired-end-sequences.qza \
  --p-trunc-len-r 230 --p-trunc-len-f 260 \
  --p-trim-left-r 20 --p-trim-left-f 19 \
  --o-table table_260230.qza \
  --o-representative-sequences rep-seqs_260230.qza \
  --o-denoising-stats denoising-stats_260230.qza \
  --p-n-threads 20

#### CLUSTERING

# Summarize feature table and sequences
qiime metadata tabulate \
  --m-input-file denoising-stats_260230.qza \
  --o-visualization denoising-stats_260230.qzv
qiime feature-table summarize \
  --i-table table_260230.qza \
  --o-visualization table_260230.qzv \
  --m-sample-metadata-file $METADATA
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_260230.qza \
  --o-visualization rep-seqs_260230.qzv
```

2. forward 270, reverse 240 (middle)

`mkdir denoise_270240`

`nano denoise_270240.sh`

```bash
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="script_error_270240" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script_270240" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/PJB_16S/processed_data/denoise_270240

source /usr/share/Modules/init/sh # load the module function
module load QIIME2/2021.8

#### METADATA FILES ####

# Metadata path
METADATA="../../metadata/metadata.txt"

# Sample manifest path
MANIFEST="../../metadata/sample_manifest.txt"

#########################

qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs ../../PJB_16S-paired-end-sequences.qza \
  --p-trunc-len-r 240 --p-trunc-len-f 270 \
  --p-trim-left-r 20 --p-trim-left-f 19 \
  --o-table table_270240.qza \
  --o-representative-sequences rep-seqs_270240.qza \
  --o-denoising-stats denoising-stats_270240.qza \
  --p-n-threads 20

#### CLUSTERING

# Summarize feature table and sequences
qiime metadata tabulate \
  --m-input-file denoising-stats_270240.qza \
  --o-visualization denoising-stats_270240.qzv
qiime feature-table summarize \
  --i-table table_270240.qza \
  --o-visualization table_270240.qzv \
  --m-sample-metadata-file $METADATA
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_270240.qza \
  --o-visualization rep-seqs_270240.qzv
```

3. forward 280, reverse 240 (least conservative)

`mkdir denoise_280240`

`nano denoise_280240.sh`

```bash
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="script_error_280240" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script_280240" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/PJB_16S/processed_data/denoise_280240

source /usr/share/Modules/init/sh # load the module function
module load QIIME2/2021.8

#### METADATA FILES ####

# Metadata path
METADATA="../../metadata/metadata.txt"

# Sample manifest path
MANIFEST="../../metadata/sample_manifest.txt"

#########################

qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs ../../PJB_16S-paired-end-sequences.qza \
  --p-trunc-len-r 240 --p-trunc-len-f 280 \
  --p-trim-left-r 20 --p-trim-left-f 19 \
  --o-table table_280240.qza \
  --o-representative-sequences rep-seqs_280240.qza \
  --o-denoising-stats denoising-stats_280240.qza \
  --p-n-threads 20

#### CLUSTERING

# Summarize feature table and sequences
qiime metadata tabulate \
  --m-input-file denoising-stats_280240.qza \
  --o-visualization denoising-stats_280240.qzv
qiime feature-table summarize \
  --i-table table_280240.qza \
  --o-visualization table_280240.qzv \
  --m-sample-metadata-file $METADATA
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_280240.qza \
  --o-visualization rep-seqs_280240.qzv
```

##### Denoising statistics

```bash
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/denoise_260230/denoising-stats_260230.qzv /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/

scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/denoise_270240/denoising-stats_270240.qzv /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/

scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/denoise_280240/denoising-stats_280240.qzv /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/
```

Open [qiime2 view](https://view.qiime2.org/) and drop in the first file you want to view. Click 'Download metadata TSV file' and save that file to 'MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data' folder.

R script:

```bash
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
```

add denoising-percent.png
add denoising-reads.png

Based on above, I will proceed using the most conservative parameters (highest quality; forward 260; reverse 230).

#### QIIME2 taxonomic Identification

Download the database file to the metadata folder on andromeda. We chose the Silva 138 99% OTUs from 515F/806R region of sequences (MD5: e05afad0fe87542704be96ff483824d4) as the classifier because we used 515F and 806RB primers for our sequences and QIIME2 recommends the classify-sklearn classifier trainer.

`cd /data/putnamlab/kevin_wong1/PJB_16S/metadata`

`wget https://data.qiime2.org/2021.4/common/silva-138-99-515-806-nb-classifier.qza`

Create two folders for the taxonomic output:

`cd /data/putnamlab/kevin_wong1/PJB_16S/processed_data/`

`mkdir filtered_taxonomy`
`mkdir unfiltered_taxonomy`

`nano taxonomy.sh`

```bash
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="script_error_taxonomy" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script_taxonomy" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/PJB_16S/processed_data/

source /usr/share/Modules/init/sh # load the module function
module load QIIME2/2021.8

#### METADATA FILES ####

# Metadata path
METADATA="../metadata/metadata.txt"

# Sample manifest path
MANIFEST="../metadata/sample_manifest.txt"

#########################

#### TAXONOMY CLASSIFICATION

qiime feature-classifier classify-sklearn \
  --i-classifier ../metadata/silva-138-99-515-806-nb-classifier.qza \
  --i-reads denoise_260230/rep-seqs_260230.qza \
  --o-classification taxonomy_260230.qza

## UNFILTERED

qiime metadata tabulate \
    --m-input-file taxonomy_260230.qza \
    --o-visualization taxonomy_260230.qzv
qiime taxa barplot \
    --i-table denoise_260230/table_260230.qza \
    --i-taxonomy taxonomy_260230.qza \
    --m-metadata-file $METADATA \
    --o-visualization unfiltered_taxonomy/taxa-bar-plots-unfiltered.qzv
qiime metadata tabulate \
    --m-input-file denoise_260230/rep-seqs_260230.qza \
    --m-input-file taxonomy_260230.qza \
    --o-visualization unfiltered_taxonomy/tabulated-feature-metadata.qzv

## FILTERED
qiime taxa filter-table \
     --i-table denoise_260230/table_260230.qza \
     --i-taxonomy taxonomy_260230.qza \
     --p-mode contains \
     --p-exclude "Unassigned","Chloroplast","Eukaryota" \
     --o-filtered-table filtered_taxonomy/table-filtered_260230.qza

qiime metadata tabulate \
    --m-input-file taxonomy_260230.qza \
    --o-visualization taxonomy_260230.qzv
qiime taxa barplot \
    --i-table filtered_taxonomy/table-filtered_260230.qza \
    --i-taxonomy taxonomy_260230.qza \
    --m-metadata-file $METADATA \
    --o-visualization filtered_taxonomy/taxa-bar-plots-filtered.qzv
qiime metadata tabulate \
    --m-input-file denoise_260230/rep-seqs_260230.qza \
    --m-input-file taxonomy_260230.qza \
    --o-visualization filtered_taxonomy/tabulated-feature-metadata.qzv
```

```
# unfiltered
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/unfiltered_taxonomy/taxa-bar-plots-unfiltered.qzv /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/

# filtered
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/filtered_taxonomy/taxa-bar-plots-filtered.qzv /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/
```
