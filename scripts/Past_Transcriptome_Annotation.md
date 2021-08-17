
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
