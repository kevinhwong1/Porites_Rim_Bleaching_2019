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
``
