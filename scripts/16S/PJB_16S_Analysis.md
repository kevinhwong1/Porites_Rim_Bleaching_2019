# PJB 16S Analysis

Following [E. Strand pipeline](https://emmastrand.github.io/EmmaStrand_Notebook/E5-16S-Analysis/) and [A. Huffmyer pipeline](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Mcapitata-Development-16S-Analysis-Part-3/)

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

##### MultiQC general statistics:

|       Sample Name       | % Dups | % GC | Length | M Seqs |
|:-----------------------:|:------:|:----:|:------:|:------:|
| WSH129_S148_L001_R1_001 | 95.40% |  48% | 259 bp |   0.1  |
| WSH129_S148_L001_R2_001 | 90.40% |  48% | 262 bp |   0.1  |
| WSH130_S160_L001_R1_001 | 94.70% |  48% | 260 bp |   0.1  |
| WSH130_S160_L001_R2_001 | 88.80% |  49% | 263 bp |   0.1  |
| WSH131_S172_L001_R1_001 | 90.90% |  47% | 260 bp |    0   |
| WSH131_S172_L001_R2_001 | 78.20% |  46% | 262 bp |    0   |
| WSH132_S184_L001_R1_001 | 94.40% |  52% | 261 bp |   0.1  |
| WSH132_S184_L001_R2_001 | 89.60% |  52% | 263 bp |   0.1  |
| WSH133_S101_L001_R1_001 | 94.70% |  48% | 248 bp |    0   |
| WSH133_S101_L001_R2_001 | 90.60% |  49% | 252 bp |    0   |
| WSH134_S113_L001_R1_001 | 94.10% |  49% | 258 bp |    0   |
| WSH134_S113_L001_R2_001 | 90.80% |  50% | 262 bp |    0   |
| WSH135_S125_L001_R1_001 | 90.60% |  51% | 266 bp |    0   |
| WSH135_S125_L001_R2_001 | 83.00% |  51% | 268 bp |    0   |
| WSH136_S137_L001_R1_001 | 91.60% |  48% | 269 bp |    0   |
| WSH136_S137_L001_R2_001 | 85.20% |  49% | 272 bp |    0   |
| WSH137_S149_L001_R1_001 | 93.20% |  45% | 250 bp |    0   |
| WSH137_S149_L001_R2_001 | 84.40% |  45% | 254 bp |    0   |
| WSH138_S161_L001_R1_001 | 93.00% |  50% | 271 bp |   0.1  |
| WSH138_S161_L001_R2_001 | 86.60% |  51% | 274 bp |   0.1  |
| WSH139_S173_L001_R1_001 | 94.50% |  47% | 245 bp |    0   |
| WSH139_S173_L001_R2_001 | 88.80% |  47% | 249 bp |    0   |
| WSH140_S185_L001_R1_001 | 95.00% |  50% | 252 bp |   0.1  |
| WSH140_S185_L001_R2_001 | 90.40% |  50% | 256 bp |   0.1  |
| WSH141_S102_L001_R1_001 | 92.80% |  49% | 253 bp |    0   |
| WSH141_S102_L001_R2_001 | 85.40% |  50% | 258 bp |    0   |
| WSH142_S114_L001_R1_001 | 92.50% |  50% | 262 bp |    0   |
| WSH142_S114_L001_R2_001 | 85.50% |  51% | 265 bp |    0   |
| WSH143_S126_L001_R1_001 | 92.50% |  50% | 267 bp |    0   |
| WSH143_S126_L001_R2_001 | 84.80% |  50% | 271 bp |    0   |
| WSH144_S138_L001_R1_001 | 95.10% |  49% | 256 bp |    0   |
| WSH144_S138_L001_R2_001 | 88.40% |  49% | 261 bp |    0   |
| WSH145_S150_L001_R1_001 | 93.60% |  49% | 266 bp |    0   |
| WSH145_S150_L001_R2_001 | 86.50% |  49% | 269 bp |    0   |
| WSH146_S162_L001_R1_001 | 93.50% |  49% | 260 bp |    0   |
| WSH146_S162_L001_R2_001 | 86.00% |  49% | 265 bp |    0   |
| WSH147_S174_L001_R1_001 | 95.30% |  49% | 253 bp |    0   |
| WSH147_S174_L001_R2_001 | 90.20% |  49% | 259 bp |    0   |
| WSH148_S186_L001_R1_001 | 92.70% |  48% | 252 bp |    0   |
| WSH148_S186_L001_R2_001 | 85.50% |  48% | 256 bp |    0   |
| WSH149_S103_L001_R1_001 | 92.80% |  49% | 259 bp |    0   |
| WSH149_S103_L001_R2_001 | 86.00% |  50% | 264 bp |    0   |
| WSH150_S115_L001_R1_001 | 92.90% |  46% | 264 bp |    0   |
| WSH150_S115_L001_R2_001 | 87.70% |  47% | 267 bp |    0   |
| WSH151_S127_L001_R1_001 | 92.30% |  49% | 255 bp |    0   |
| WSH151_S127_L001_R2_001 | 84.70% |  49% | 259 bp |    0   |
| WSH152_S139_L001_R1_001 | 92.00% |  51% | 259 bp |    0   |
| WSH152_S139_L001_R2_001 | 83.20% |  51% | 263 bp |    0   |
| WSH153_S151_L001_R1_001 | 92.50% |  49% | 253 bp |    0   |
| WSH153_S151_L001_R2_001 | 82.60% |  49% | 257 bp |    0   |
| WSH154_S163_L001_R1_001 | 89.90% |  43% | 253 bp |    0   |
| WSH154_S163_L001_R2_001 | 76.90% |  43% | 257 bp |    0   |
| WSH155_S175_L001_R1_001 | 92.40% |  49% | 252 bp |    0   |
| WSH155_S175_L001_R2_001 | 84.00% |  48% | 257 bp |    0   |
| WSH156_S187_L001_R1_001 | 93.80% |  48% | 237 bp |    0   |
| WSH156_S187_L001_R2_001 | 87.70% |  48% | 242 bp |    0   |
| WSH157_S104_L001_R1_001 | 95.30% |  51% | 251 bp |    0   |
| WSH157_S104_L001_R2_001 | 91.20% |  51% | 254 bp |    0   |
| WSH158_S116_L001_R1_001 | 95.70% |  50% | 250 bp |    0   |
| WSH158_S116_L001_R2_001 | 92.10% |  50% | 254 bp |    0   |
| WSH159_S128_L001_R1_001 | 94.50% |  51% | 255 bp |    0   |
| WSH159_S128_L001_R2_001 | 89.10% |  51% | 259 bp |    0   |
| WSH160_S140_L001_R1_001 | 92.10% |  51% | 263 bp |   0.1  |
| WSH160_S140_L001_R2_001 | 82.90% |  49% | 266 bp |   0.1  |
| WSH161_S152_L001_R1_001 | 93.20% |  51% | 256 bp |   0.1  |
| WSH161_S152_L001_R2_001 | 86.00% |  51% | 259 bp |   0.1  |
| WSH162_S164_L001_R1_001 | 94.60% |  50% | 252 bp |   0.1  |
| WSH162_S164_L001_R2_001 | 88.70% |  50% | 256 bp |   0.1  |
| WSH163_S176_L001_R1_001 | 95.10% |  48% | 246 bp |    0   |
| WSH163_S176_L001_R2_001 | 89.90% |  49% | 251 bp |    0   |
| WSH164_S188_L001_R1_001 | 94.00% |  48% | 243 bp |    0   |
| WSH164_S188_L001_R2_001 | 87.60% |  48% | 248 bp |    0   |
| WSH165_S105_L001_R1_001 | 94.90% |  51% | 247 bp |    0   |
| WSH165_S105_L001_R2_001 | 89.00% |  51% | 252 bp |    0   |
| WSH166_S117_L001_R1_001 | 94.70% |  51% | 253 bp |    0   |
| WSH166_S117_L001_R2_001 | 90.00% |  51% | 256 bp |    0   |
| WSH167_S129_L001_R1_001 | 94.20% |  51% | 263 bp |    0   |
| WSH167_S129_L001_R2_001 | 88.00% |  51% | 266 bp |    0   |
| WSH168_S141_L001_R1_001 | 89.10% |  49% | 258 bp |    0   |
| WSH168_S141_L001_R2_001 | 75.30% |  49% | 264 bp |    0   |
| WSH169_S153_L001_R1_001 | 93.60% |  46% | 237 bp |    0   |
| WSH169_S153_L001_R2_001 | 83.80% |  46% | 244 bp |    0   |
| WSH170_S165_L001_R1_001 | 94.40% |  48% | 253 bp |   0.1  |
| WSH170_S165_L001_R2_001 | 88.60% |  48% | 257 bp |   0.1  |
| WSH171_S177_L001_R1_001 | 93.90% |  49% | 254 bp |   0.1  |
| WSH171_S177_L001_R2_001 | 87.70% |  49% | 258 bp |   0.1  |
| WSH172_S189_L001_R1_001 | 95.90% |  50% | 249 bp |   0.1  |
| WSH172_S189_L001_R2_001 | 91.40% |  50% | 253 bp |   0.1  |
| WSH173_S106_L001_R1_001 | 91.70% |  50% | 247 bp |    0   |
| WSH173_S106_L001_R2_001 | 84.50% |  50% | 250 bp |    0   |

##### Sequence Counts

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_sequence_counts_plot.png)

##### Sequence Quality

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_per_base_sequence_quality_plot.png)

##### Per Sequence Quality Scores

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_per_sequence_quality_scores_plot.png)

##### Per Sequence GC Content

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_per_sequence_gc_content_plot.png)

##### Per Base N Content

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_per_base_n_content_plot.png)

##### Sequence Length Distribution

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_sequence_length_distribution_plot.png)

##### Sequence Duplication Levels

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_sequence_duplication_levels_plot.png)

##### Overrepresented Sequences

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_overrepresented_sequencesi_plot.png)

##### Adapter Content

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc_adapter_content_plot.png)

##### Status Check

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/fastqc_plots/fastqc-status-check-heatmap.png)


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

**R script:**
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

##### Denoising Percentage

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/processed_data/denoising-percent.png)

##### Denoising Reads

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/processed_data/denoising-reads.png)


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

##### Unfiltered level 5 taxonomy

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/processed_data/tax_lv5_unfiltered.png)

##### Filtered level 5 taxonomy

![ ](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019/blob/master/output/16S/processed_data/tax_lv5_filtered.png)


#### Phylogenetic trees

`nano phylo_tree.sh`

``` bash
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                 
#SBATCH --error="tree_script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="tree_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/PJB_16S/processed_data/

source /usr/share/Modules/init/sh # load the module function
module load QIIME2/2021.8

# align and mask sequences
qiime alignment mafft \
  --i-sequences denoise_260230/rep-seqs_260230.qza \
  --o-alignment aligned-rep-seqs_260230.qza
qiime alignment mask \
  --i-alignment aligned-rep-seqs_260230.qza \
  --o-masked-alignment masked-aligned-rep-seqs_260230.qza

# calculate tree
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs_260230.qza \
  --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

Export files to local desktop:

```bash
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/rooted-tree.qza /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/

scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/unrooted-tree.qza /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/
```

#### Calculate diversity metrics

`nano diversity.sh`

```bash
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                 
#SBATCH --error="idfiltered_script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="idfiltered_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/PJB_16S/processed_data/

source /usr/share/Modules/init/sh # load the module function
module load QIIME2/2021.8

# Metadata path
METADATA="../metadata/metadata.txt"

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table filtered_taxonomy/table-filtered_260230.qza \
  --p-sampling-depth 95 \
  --m-metadata-file $METADATA \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --m-metadata-column Timepoint \
  --o-visualization core-metrics-results/unweighted-unifrac-station-significance.qzv \
  --p-pairwise
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA  \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/unweighted-unifrac-group-significance.qzv \
  --p-pairwise

# This script calculates the rarefaction curve for the data
  qiime diversity alpha-rarefaction \
    --i-table filtered_taxonomy/table-filtered_260230.qza \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth 800 \
    --m-metadata-file $METADATA \
    --o-visualization alpha-rarefaction.qzv
```

Export files to local desktop:

```bash
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/*qza /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/

scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/alpha-rarefaction.qzv  /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/

scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/processed_data/core-metrics-results /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/16S/processed_data/
```
