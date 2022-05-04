# Porites bleaching TagSeq Bioinformatics

## Script modified from my awesome lab mates [S. Gurr](https://github.com/SamGurr/Pgenerosa_TagSeq_Metabolomics/blob/main/TagSeq/HPC_work/Geoduck_TagSeq_Bioinf.md#Initial-diagnostics-upon-sequence-upload-to-HPC), [E. Chille]() and A. Huffmyer.


# Initial diagnostics upon sequence upload to HPC
--------------------------------------------

HP uploaded data to Bluewaves and we are waiting for the md5 checksum from UT.

HP ran initial multiQC

# Trimming and post-trim quality check of 'clean' reads
-------------------------------------------

### shell script: <span style="color:green">**fastp_multiqc.sh**<span>

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

### EXPORT MULTIQC REPORT
*exit bluewaves and run from terminal*
- save to gitrepo as multiqc_clean.html
```
scp kevin_wong1@bluewaves.uri.edu:/data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/clean/multiqc_report.html  MyProjects/Porites_Rim_Bleaching_2019/output/TagSeq
```


# Merge files for each sample across lanes  
-----------------------------------------------------------------

`nano cat.clean.sh`

```bash
#!/bin/bash
#SBATCH --job-name="cat.clean"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/clean
#SBATCH --mem=100GB

module load SAMtools/1.9-foss-2018b

echo "Concatenate files" $(date)

ls *R1_001.fastq.gz | awk -F '[_]' '{print $1"_"$2}' | sort | uniq > ID

for i in `cat ./ID`;
	do cat $i\_L001_R1_001.fastq.gz $i\_L002_R1_001.fastq.gz > $i\_ALL.fastq.gz;
	done

echo "Mission complete." $(date)
```

# HISAT2
-----------------------------------------------------------------

* Genome file path: `/data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta`
* GFF genes path: `/data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gff`


`nano align.sh`

```bash
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="align_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Past_Genome_final

# load modules needed
module load HISAT2/2.2.1-foss-2019b
module load SAMtools/1.9-foss-2018b

# index the reference genome for Pocillopora acuta output index to working directory
hisat2-build -f /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta ./Pastreoides_ref # called the reference genome (scaffolds)
echo "Reference genome indexed. Starting alingment" $(date)

# symbolically link 'clean' reads to hisat2 dir
ln -s ../../clean/clean*_ALL.fastq.gz ./

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls *.fastq.gz)) # call the clean sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
	hisat2 -p 8 --dta -x Pastreoides_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
    		echo "${i} bam-ified!"
        rm ${sample_name}.sam
done
```

# StringTie2
-----------------------------------------------------------------

`nano stringtie.sh`

```bash
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="stringtie_error" #if your job fails, the error report will be put in this file
#SBATCH --output="stringtie_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Past_Genome_final

#load packages
module load StringTie/2.1.4-GCC-9.3.0

#Transcript assembly: StringTie

array=($(ls *_ALL.bam)) #Make an array of sequences to assemble

for i in ${array[@]}; do #Running with the -e option to compare output to exclude novel genes. Also output a file with the gene abundances
        sample_name=`echo $i| awk -F [_] '{print $1"_"$2"_"$3}'`
	stringtie -p 8 -e -B -G /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gff -A ${sample_name}.gene_abund.tab -o ${sample_name}.gtf ${i}
        echo "StringTie assembly for seq file ${i}" $(date)
done
echo "StringTie assembly COMPLETE, starting assembly analysis" $(date)
```

# Prep DE
-----------------------------------------------------------------

`cp /data/putnamlab/ashuffmyer/pairs-rnaseq/prepDE.py .`

`nano prepDE.sh`

```bash
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="prepDE_error" #if your job fails, the error report will be put in this file
#SBATCH --output="prepDE_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Past_Genome_final

#load packages
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.1.4-GCC-9.3.0

#Transcript assembly: StringTie
module load GffCompare/0.12.1-GCCcore-8.3.0

#Transcript assembly QC: GFFCompare

#make gtf_list.txt file
ls *.bam.gtf > gtf_list.txt

stringtie --merge -p 8 -G /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gff -o Past_merged.gtf gtf_list.txt
echo "Stringtie merge complete" $(date)

gffcompare -r /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gff -G -o merged Past_merged.gtf

echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

#make gtf list text file
for filename in *.bam.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

python prepDE.py -g PJB_gene_count_matrix.csv -i listGTF.txt #Compile the gene count matrix
echo "Gene count matrix compiled." $(date)
```

`scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/20210315_Past_Tagseq/output/hisat2/Past_Genome_final/PJB_gene_count_matrix.csv /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/output/TagSeq/
`
