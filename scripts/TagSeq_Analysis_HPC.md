# Porites bleaching TagSeq Bioinformatics

## Script modified from [S. Gurr](https://github.com/SamGurr/Pgenerosa_TagSeq_Metabolomics/blob/main/TagSeq/HPC_work/Geoduck_TagSeq_Bioinf.md#Initial-diagnostics-upon-sequence-upload-to-HPC)

## <span style="color:blue">**Workflow**</span>
  - [Upon upload to HPC...](#Initial-diagnostics-upon-sequence-upload-to-HPC)
	  - [Upon upload to HPC...](#Initial-diagnostics-upon-sequence-upload-to-HPC)
      - [Count raw reads](#Count-the-number-of-read-files)
      - [Digital fingerprint md5sums](#run-checksum) (HPC script; <span style="color:green">**md5_checksum.sh**<span>)
  - [1. MultiQC: Initial QC](#Quality-check-of-raw-reads) (HPC script; <span style="color:green">**mutliqc.sh**<span>)
  - [2. Trimming and QC of 'clean' reads](#Trimming-and-post-trim-quality-check-of-'clean'-reads)
  	- [Remember! polyA tail in TagSeq](#Trimming-polyA-tail)
	- [fastp - about/commands](#What-this-script-will-do...)
	- [fastp and MulitiQC: Trim and QC](#shell-script-fastp_multiqc.sh) (HPC script; <span style="color:green">**fastp_mutliqc.sh**<span>)
  - [3. Alignment of cleaned reads to reference](#HISAT2-Alignment-of-cleaned-reads-to-reference)
	- [Upload reference genome](#Reference-genome-upload-to-HPC)
	- [HISAT2 - about/commands](#HISAT2-alignment)
	- [samtools - about/commands](#samtools)
	- [HISAT2: Index reference and alignment](#HPC-Job-HISAT2-Index-Reference-and-Alignment) (HPC script; <span style="color:green">**HISAT2.sh**<span>)
  - [4. Assembly and quantification](#Assembly-and-quantification)
  	- [Upload annotation reference for assembly](#Upload-annotation-reference-gff-or-gff3-to-HPC)
	- [StringTie2 - about/commands](#StringTie)
	- [gffcomapare - about/commands](#gffcompare)
	- [prepDE.py - about/commands](#Python-step-prepDE.py) (essential prep, load open-source python script for count matrix)
		- [List HISAT2 output files (for --merge in Stringtie2)](#gtf_list.txt-run) (essential Stringtie2 ```--merge``` prep; **gtf_list.txt** file!)
		- [List HISAT2 output files (for prepDE.py)](#listGTF.txt-run) (essential prepDE.py prep; **listGTF.txt** file!)
	- [Stringtie2: Assembly step](#HPC-job-Assembly) (HPC script; <span style="color:green">**Stringtie2.sh**<span>)
	- [Stringtie2, gffcomapre, prepDE.py: Merge and build read count matrix for DEG analysis](#HPC-job-Merge-and-Build-Read-Count-Matrix-for-DEG-analysis) (HPC script; <span style="color:green">**Stringtie2_merge_prepDEpy.sh**<span>)


# Initial diagnostics upon sequence upload to HPC
--------------------------------------------

HP uploaded data to Bluewaves and we are waiting for the md5 checksum from UT.

HP ran initial multiQC

# Trimming and post-trim quality check of 'clean' reads
-------------------------------------------

# shell script: <span style="color:green">**fastp_multiqc.sh**<span>

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --output=../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/"%x_out.%j"
#SBATCH --error=../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/"%x_err.%j"
#SBATCH -D /data/putnamlab/KITT/hputnam/20210312_Pastreoides_TagSeq/Pastreoides_TagSeq_JA21015
#SBATCH --cpus-per-task=3

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.7-foss-2018b-Python-2.7.15

# Make an array of sequences to trim
array1=($(ls *.fastq.gz))

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA)
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 ../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 20 -y -Y 50
        fastqc ../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/clean.${i}
done

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

cd ../../../../kevin_wong1/20210315_Past_Tagseq/output/clean/ #The following command will be run in the /clean directory

multiqc ./ #Compile MultiQC report from FastQC files

echo "Cleaned MultiQC report generated." $(date)
```

20210318 job #1883852

### EXPORT MUTLIQC REPORT
*exit bluewaves and run from terminal*
- save to gitrepo as multiqc_clean.html
```
scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/clean/multiqc_report.html  MyProjects/Porites_Rim_Bleaching_2019/output/TagSeq
```


## HPC Job: HISAT2 Index Reference and Alignment
-----------------------------------------------------------------

- create directory output\hisat2

``` mkdir HISAT2 ```


### Mapping #1 - Porites astreodes transcriptome from (Kenkel et al. 2013)[https://matzlab.weebly.com/data--code.html]

- index reference and alignment

uploading reference transcriptome from (Kenkel et al. 2013)[https://matzlab.weebly.com/data--code.html]

```
scp Desktop/URI_PHD/pastreoides_2014_ref/pastreoides_may2014/past.fasta kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/REFS/Past/Kenkel2013_past_transcriptome.fasta
```

**input**
- Kenkel2013_past_transcriptome.fasta *= reference transcriptome*
- clean/*.fastq.gz *= all clean TagSeq reads*

**ouput**
- Past_transcriptome_ref *= indexed reference by hisat2-build; stored in the output/hisat2 folder as 1.hy2, 2.ht2... 8.ht2*
- <clean.fasta>.sam *=hisat2 output, readable text file; removed at the end of the script*
- <clean.fasta>.bam *=converted binary file complementary to the hisat sam files*

# shell script: <span style="color:green">**HISAT2.sh**<span>

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Kenkel2013_REF
#SBATCH --cpus-per-task=3


#load packages
module load HISAT2/2.1.0-foss-2018b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# symbolically link 'clean' reads to hisat2 dir
ln -s ../../../clean/clean*.fastq.gz ./

# index the reference transcriptome for Porites astreoides output index to working directory
hisat2-build -f ../../../../REFS/Past/Kenkel2013_past_transcriptome.fasta ./Past_transcriptome_ref # called the reference transcriptome (scaffolds)
echo "Referece transcriptome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --dta -x Past_transcriptome_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

```

20210319 Submitted job 1883908

Took ~2 days

This job "failed" but I think it was because I had an extra "clean*.fastq.qz" file. All of the samples had mapping percentages and sam files were converted too bam then deleted.

```
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
Extra parameter(s) specified: "clean.A10-KW10_S233_L002_R1_001.fastq.gz", "clean.A11-KW11_S236_L001_R1_001.fastq.gz", "clean.A11-KW11_S236_L002_R1_001.fastq.gz", "clean.A12-KW12_S239_L001_R1_001.fastq.gz", "clean.A12-KW12_S239_L002_R1_001.fastq.gz", "clean.A$
Error: Encountered internal HISAT2 exception (#1)
Command: /opt/software/HISAT2/2.1.0-foss-2018b/bin/hisat2-align-s --wrapper basic-0 -p 8 --dta -x Past_transcriptome_ref -S A10-KW10_S233_L001_R1_001.sam -U /tmp/31840.unp clean.A10-KW10_S233_L002_R1_001.fastq.gz clean.A11-KW11_S236_L001_R1_001.fastq.gz clea$
(ERR): hisat2-align exited with value 1
[E::hts_open_format] Failed to open file A10-KW10_S233_L001_R1_001.sam
samtools sort: can't open "A10-KW10_S233_L001_R1_001.sam": No such file or directory
rm: cannot remove `A10-KW10_S233_L001_R1_001.sam': No such file or directory
```

Exporting to report mapping percentages

```
scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Kenkel2013_REF/HISAT2.sh_err.1883908  MyProjects/Porites_Rim_Bleaching_2019/output/TagSeq
```


### Salmon

replaces HISAT2 and StringTie because I have no reference genome with a GFF3 file

https://combine-lab.github.io/salmon/getting_started/

```
#!/bin/bash
#SBATCH --job-name="Salmon"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/salmon
#SBATCH --mem=100GB

module load Salmon/1.1.0-gompi-2019b

echo "Building index..." $(date)
salmon index -t ../../../REFS/Past/Kenkel2013_past_transcriptome.fasta -i PAST_Kenkel_ref_index -k 31

echo "Running the alignment and abundance estimation." $(date)
salmon quant -i PAST_Kenkel_ref_index -l A \
        -r ../clean/clean*.fastq.gz \
        --validateMappings -p 8 -o transcripts_quant

echo "Mission complete." $(date)

```
























# Annotating Kenkel transcriptome

Following this pipeline by E. Chille:
https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/1-BLAST-GO-KO/2020-10-08-M-capitata-functional-annotation-pipeline.md

### DIAMOND
https://github.com/bbuchfink/diamond/wiki

diamond_past_kenkel.sh

```
#!/bin/bash
#SBATCH --job-name="diamond" #CHANGE_NAME
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu #CHANGE_EMAIL
#SBATCH --mem=100GB

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

cd /data/putnamlab/kevin_wong1/REFS/Past/annotations/

echo "Updating Past transcriptome annotation" $(date)

diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q ../Kenkel2013_past_transcriptome.fasta -o Past.Kenkel.annot.20210608 -f 100 -b20 --more-sensitive -e 0.00001 -k1

echo "Search complete... converting format to XML and tab"

diamond view -a Past.Kenkel.annot.20210608.daa -o Past.Kenkel.annot.20210608.xml -f 5
diamond view -a Past.Kenkel.annot.20210608.daa -o Past.Kenkel.annot.20210608.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)

```



interproscan.sh

```
#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=erin_chille@uri.edu
#SBATCH -p putnamlab

cd /data/putnamlab/erin_chille/mcap2019/annotations/

echo "START $(date)"

# Load module
module load InterProScan/5.46-81.0-foss-2019b
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i ../data/ref/Mcap.IPSprotein.fa -b ./Mcap.interpro.200824  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i ./Mcap.interpro.200824.xml -b ./Mcap.interpro.200824

# -i is the input data
# -b is the output file base
# -f is formats
# -iprlookup enables mapping
# -goterms is GO Term
# -pa is pathway mapping
# -version displays version number

echo "DONE $(date)"

```

kofam.sh

```

#!/bin/bash
#SBATCH --job-name="KofamScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=/opt/software/kofam_scan/1.3.0-foss-2019b,/opt/software/HMMER/3.3.2-gompi-2019b/bin/
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=erin_chille@uri.edu
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=500GB
#SBATCH -D /data/putnamlab/erin_chille/mcap2019/annotations/

echo "Loading modules" $(date)
module load kofam_scan/1.3.0-foss-2019b
module load libyaml/0.1.5
module unload HMMER/3.3.1-foss-2019b
module load HMMER/3.3.2-gompi-2019b
module list

#echo "Starting analysis... downloading KO database" $(date)
#wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz #download KO database
#wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
#gunzip ko_list.gz
#tar xf profiles.tar.gz

echo "Beginning mapping" $(date)
/opt/software/kofam_scan/1.3.0-foss-2019b/exec_annotation -o Mcap_KO_annot.txt -k ./ko_list -p ./profiles/eukaryote.hal -E 0.00001 -f detail-tsv --report-unannotated /data/putnamlab/erin_chille/mcap2019/data/ref/Mcap.protein.fa

echo "Analysis complete!" $(date)
```













### Mapping #2 - Porites lutea genome


**input**
- plut_final_2.1.fasta.gz *= reference genome*
- clean/*.fastq.gz *= all clean TagSeq reads*

**ouput**
- Plutea_genome_ref *= indexed reference by hisat2-build; stored in the output/hisat2 folder as 1.hy2, 2.ht2... 8.ht2*
- <clean.fasta>.sam *=hisat2 output, readable text file; removed at the end of the script*
- <clean.fasta>.bam *=converted binary file complementary to the hisat sam files*

# shell script: <span style="color:green">**HISAT2_Plutea.sh**<span>

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --output="%x_out.%j"
#SBATCH --error="%x_err.%j"
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Plutea_REF
#SBATCH --cpus-per-task=3


#load packages
module load HISAT2/2.1.0-foss-2018b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

# symbolically link 'clean' reads to hisat2 dir
ln -s ../../clean/clean*.fastq.gz ./

# index the reference genome for Porites lutea output index to working directory
hisat2-build -f ../../../../../REFS/Plutea/plut_final_2.1.fasta.gz ./Plutea_genome_ref # called the reference genome (scaffolds)
echo "Reference genome indexed. Starting alignment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the symbolically linked sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --dta -x Plutea_genome_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

```

Submitted batch job 36970 on Andromeda on March 26 at 13:22.
-- stuck - genome is not indexing properly
