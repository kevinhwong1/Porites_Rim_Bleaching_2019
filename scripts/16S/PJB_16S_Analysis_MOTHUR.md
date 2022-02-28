# 16S analysis using MOTHUR

This script is taken from [Ariana Huffmyer's post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/16S-Analysis-in-Mothr-Part-1/).

1. Prepare directory

```
cd /data/putnamlab/kevin_wong1/PJB_16S
mkdir mothur

cp /data/putnamlab/kevin_wong1/PJB_16S/raw_data/*.gz /data/putnamlab/kevin_wong1/PJB_16S/mothur
```

2. Start Mothur

```
interactive
module load Mothur/1.46.1-foss-2020b
cd mothur/
mothur
```

3. Preparing Sequences: make.file, make.contig, and summary.seq


```
nano oligos.oligos

primer GTGCCAGCMGCCGCGGTAA GGACTACNVGGGTWTCTAAT
```

`nano contigs.sh`

```bash
#!/bin/bash
#SBATCH --job-name="contigs"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="contigs_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="contigs_output_script" #once your job is completed, any final job report comments will be put in this file

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#make.file(inputdir=., type=gz, prefix=PJB)"

mothur "#make.contigs(inputdir=., outputdir=., file=PJB.files, oligos=oligos.oligos, trimoverlap=T)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.fasta)"
```

Summary output: `less contigs_output_script`

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1       1       0       1       1
2.5%-tile:      1       252     252     0       3       31747
25%-tile:       1       253     253     0       4       317468
Median:         1       253     253     0       4       634936
75%-tile:       1       253     253     0       5       952404
97.5%-tile:     1       254     254     6       6       1238125
Maximum:        1       280     280     114     30      1269871
Mean:           1       252     252     0       4
# of Seqs:      1269871
```

Check that primers are gone:
* F GTGCCAGCMGCCGCGGTAA
* R GGACTACNVGGGTWTCTAAT

`head PJB.trim.contigs.fasta`


4. QC’ing sequences with screen.seqs

`nano screen.sh`

```bash
#!/bin/bash
#SBATCH --job-name="screen"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="screen_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="screen_output_script" #once your job is completed, any final job report comments will be put in this file

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#screen.seqs(inputdir=., outputdir=., fasta=PJB.trim.contigs.fasta, group=PJB.contigs.groups, maxambig=0, maxlength=350, minlength=250)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.good.fasta)"

```

The summary output as viewed in the screen_output_script file now reads:

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       250     250     0       3       1
2.5%-tile:      1       253     253     0       3       24682
25%-tile:       1       253     253     0       4       246819
Median:         1       253     253     0       4       493638
75%-tile:       1       253     253     0       5       740457
97.5%-tile:     1       254     254     0       6       962594
Maximum:        1       280     280     0       12      987275
Mean:           1       253     253     0       4
# of Seqs:      987275
```

5. Determining and counting unique sequences

`nano unique.sh`

```bash
#!/bin/bash
#SBATCH --job-name="unique"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="unique_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="unique_output_script" #once your job is completed, any final job report comments will be put in this file

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#unique.seqs(fasta=PJB.trim.contigs.good.fasta)"

mothur "#count.seqs(name=PJB.trim.contigs.good.names, group=PJB.contigs.good.groups)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.good.unique.fasta, count=PJB.trim.contigs.good.count_table)"

mothur "#count.groups(count= PJB.trim.contigs.good.unique.fasta)"
```

`less unique_output_script`

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       250     250     0       3       1
2.5%-tile:      1       253     253     0       3       24682
25%-tile:       1       253     253     0       4       246819
Median:         1       253     253     0       4       493638
75%-tile:       1       253     253     0       5       740457
97.5%-tile:     1       254     254     0       6       962594
Maximum:        1       280     280     0       12      987275
Mean:           1       253     253     0       4
# of unique seqs:       72883
total # of seqs:        987275
```

Warning message:

```
mothur > count.groups(count= PJB.trim.contigs.good.unique.fasta)
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAA
GTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACT
CTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAA
GTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.

**** Exceeded maximum allowed command errors, quitting ****
[ERROR]: Your count table contains a sequence named GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC with a total=0. Please correct.
fbdiffs=0(match), contains 0.
fpdiffs=0(match), contains 0.
rbdiffs=0(match) contains 0.
rpdiffs=0(match) contains 0.

Size of smallest group: 0.

Total seqs: 0.
```


6. Aligning to reference database

### Prepare the reference sequences:

```
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.bacteria.zip

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip

unzip silva.bacteria.zip
cd silva.bacteria
cp silva.bacteria.fasta ../silva.bacteria.fasta

unzip trainset9_032012.pds.zip
```

`nano silva_ref.sh`

```bash
#!/bin/bash
#SBATCH --job-name="silva_ref"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="silva_ref_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="silva_ref_output_script" #once your job is completed, any final job report comments will be put in this file

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)"

mothur "#summary.seqs(fasta=silva.bacteria.pcr.fasta)"

mothur "#rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)"
```

### Align sequences to the reference

`nano align.sh`

```bash
#!/bin/bash
#SBATCH --job-name="align"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="align_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="align_output_script" #once your job is completed, any final job report comments will be put in this file


module load Mothur/1.46.1-foss-2020b

mothur

mothur "#align.seqs(fasta=PJB.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.good.unique.align)"
```

Summary output:
`less align_output_script`

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1231    3       0       1       1
2.5%-tile:      1968    11550   252     0       3       1823
25%-tile:       1968    11550   253     0       4       18221
Median:         1968    11550   253     0       4       36442
75%-tile:       1968    11550   253     0       5       54663
97.5%-tile:     1968    11550   254     0       6       71061
Maximum:        13422   13425   275     0       12      72883
Mean:           1973    11549   252     0       4
# of Seqs:      72883
```

### QC sequences according to alignment to the reference

`nano screen2.sh`

```bash
#!/bin/bash
#SBATCH --job-name="screen2"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="screen2_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="screen2_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#screen.seqs(fasta=PJB.trim.contigs.good.unique.align, count=PJB.trim.contigs.good.count_table, start=1968, end=11550, maxhomop=8)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.good.unique.good.align, count=PJB.trim.contigs.good.good.count_table)"

mothur "#count.groups(count=PJB.trim.contigs.good.good.count_table)"

```

Summary output:
`less screen2_output_script`

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       11550   250     0       3       1
2.5%-tile:      1968    11550   253     0       3       24609
25%-tile:       1968    11550   253     0       4       246085
Median:         1968    11550   253     0       4       492170
75%-tile:       1968    11550   253     0       5       738255
97.5%-tile:     1968    11550   254     0       6       959731
Maximum:        1968    13425   275     0       8       984339
Mean:           1967    11550   253     0       4
# of unique seqs:       72472
total # of seqs:        984339
```

```
mothur > count.groups(count=PJB.trim.contigs.good.good.count_table)
WSH129 contains 61252.
WSH130 contains 32582.
WSH131 contains 7084.
WSH132 contains 83387.
WSH133 contains 14744.
WSH134 contains 11080.
WSH135 contains 25048.
WSH136 contains 12600.
WSH137 contains 8890.
WSH138 contains 66867.
WSH139 contains 16472.
WSH140 contains 29521.
WSH141 contains 16537.
WSH142 contains 18938.
WSH143 contains 24958.
WSH144 contains 20990.
WSH145 contains 25098.
WSH146 contains 19106.
WSH147 contains 17605.
WSH148 contains 16765.
WSH149 contains 19290.
WSH150 contains 12994.
WSH151 contains 19274.
WSH152 contains 13537.
WSH153 contains 11916.
WSH154 contains 1731.
WSH155 contains 7268.
WSH156 contains 14572.
WSH157 contains 21936.
WSH158 contains 16699.
WSH159 contains 13874.
WSH160 contains 22918.
WSH161 contains 31488.
WSH162 contains 37097.
WSH163 contains 16455.
WSH164 contains 21257.
WSH165 contains 14702.
WSH166 contains 12564.
WSH167 contains 11057.
WSH168 contains 7832.
WSH169 contains 10627.
WSH170 contains 33534.
WSH171 contains 35915.
WSH172 contains 34605.
WSH173 contains 11673.

Size of smallest group: 1731.

Total seqs: 984339.

Output File Names:
PJB.trim.contigs.good.good.count.summary
```


### Filter Sequences

`nano filter.sh`

```bash
#!/bin/bash
#SBATCH --job-name="filter"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="filter_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="filter_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#filter.seqs(fasta=PJB.trim.contigs.good.unique.good.align, vertical=T, trump=.)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.good.unique.good.filter.fasta, count=PJB.trim.contigs.good.good.count_table)"

mothur "#count.groups(count= PJB.trim.contigs.good.good.count_table)"

```

Summary output:
`less filter_output_script`

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       530     250     0       3       1
2.5%-tile:      1       531     253     0       3       24609
25%-tile:       1       531     253     0       4       246085
Median:         1       531     253     0       4       492170
75%-tile:       1       531     253     0       5       738255
97.5%-tile:     1       531     254     0       6       959731
Maximum:        1       531     271     0       8       984339
Mean:           1       530     253     0       4
# of unique seqs:       72472
total # of seqs:        984339
```

7. Polish the data with pre clustering

`nano precluster.sh`

```bash
#!/bin/bash
#SBATCH --job-name="precluster"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="precluster_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="precluster_output_script" #once your job is completed, any final job report comments will be put in this file


module load Mothur/1.46.1-foss-2020b

mothur

mothur "#unique.seqs(fasta=PJB.trim.contigs.good.unique.good.filter.fasta, count=PJB.trim.contigs.good.good.count_table)"

mothur "#pre.cluster(fasta=PJB.trim.contigs.good.unique.good.filter.unique.fasta, count=PJB.trim.contigs.good.unique.good.filter.count_table, diffs=1)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.count_table)"

mothur "#count.groups(count=current)"
```

Summary output:
`less precluster_output_script`

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       530     250     0       3       1
2.5%-tile:      1       531     253     0       3       24609
25%-tile:       1       531     253     0       4       246085
Median:         1       531     253     0       4       492170
75%-tile:       1       531     253     0       5       738255
97.5%-tile:     1       531     254     0       6       959731
Maximum:        1       531     271     0       8       984339
Mean:           1       530     253     0       4
# of unique seqs:       43840
total # of seqs:        984339
```

Error message:
```
mothur > count.groups(count=current)
[WARNING]: no file was saved for count parameter.
You have no current groupfile, countfile or sharedfile and one is required.
[ERROR]: did not complete count.groups.

mothur > quit()
```

8. Identify chimeras

`nano chimera.sh`

```bash
#!/bin/bash
#SBATCH --job-name="chimera"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="chimera_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="chimera_output_script" #once your job is completed, any final job report comments will be put in this file


module load Mothur/1.46.1-foss-2020b

module load VSEARCH/2.18.0-GCC-10.2.0

mothur

mothur "#chimera.vsearch(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=T)"

mothur "#remove.seqs(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"

mothur "#count.groups(count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"
```

Summary output:
`less chimera_output_script`

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       530     250     0       3       1
2.5%-tile:      1       531     253     0       3       23951
25%-tile:       1       531     253     0       4       239501
Median:         1       531     253     0       4       479002
75%-tile:       1       531     253     0       5       718503
97.5%-tile:     1       531     254     0       6       934053
Maximum:        1       531     271     0       8       958003
Mean:           1       530     253     0       4
# of unique seqs:       32766
total # of seqs:        958003


```

The count of sequences in each file are:

```
mothur > count.groups(count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
WSH129 contains 58649.
WSH130 contains 31860.
WSH131 contains 6993.
WSH132 contains 80691.
WSH133 contains 14536.
WSH134 contains 10956.
WSH135 contains 24421.
WSH136 contains 12329.
WSH137 contains 8860.
WSH138 contains 64973.
WSH139 contains 16282.
WSH140 contains 28786.
WSH141 contains 15863.
WSH142 contains 18045.
WSH143 contains 23782.
WSH144 contains 20178.
WSH145 contains 23974.
WSH146 contains 18672.
WSH147 contains 17154.
WSH148 contains 16521.
WSH149 contains 18696.
WSH150 contains 12725.
WSH151 contains 18917.
WSH152 contains 13242.
WSH153 contains 11786.
WSH154 contains 1728.
WSH155 contains 7222.
WSH156 contains 14427.
WSH157 contains 21381.
WSH158 contains 16373.
WSH159 contains 13681.
WSH160 contains 21838.
WSH161 contains 30675.
WSH162 contains 35901.
WSH163 contains 16045.
WSH164 contains 20905.
WSH165 contains 14476.
WSH166 contains 12494.
WSH167 contains 10959.
WSH168 contains 7666.
WSH169 contains 10540.
WSH170 contains 32851.
WSH171 contains 35000.
WSH172 contains 33455.
WSH173 contains 11495.

Size of smallest group: 1728.

Total seqs: 958003.

Output File Names:
PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count.summary
```


9. Classifying sequences

`nano classify.sh`

```bash
#!/bin/bash
#SBATCH --job-name="classify"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="classify_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="classify_output_script" #once your job is completed, any final job report comments will be put in this file

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#classify.seqs(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax)"

mothur "#remove.lineage(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"
```

Summary output:
`less classify_output_script`

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       530     250     0       3       1
2.5%-tile:      1       531     253     0       3       22628
25%-tile:       1       531     253     0       4       226280
Median:         1       531     253     0       4       452559
75%-tile:       1       531     253     0       5       678838
97.5%-tile:     1       531     254     0       6       882489
Maximum:        1       531     271     0       8       905116
Mean:           1       530     253     0       4
# of unique seqs:       30986
total # of seqs:        905116
```

10. Cluster for OTUs

`nano cluster.sh`

```bash
#!/bin/bash
#SBATCH --job-name="cluster"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=kevin_wong1@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="cluster_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="cluster_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#dist.seqs(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)"

mothur "#cluster(column=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03)"

mothur "#cluster.split(fasta=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, taxlevel=4, cutoff=0.03, splitmethod=classify)"

mothur "#make.shared(list=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"

mothur "#classify.otu(list=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=PJB.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)"

mothur "#rename.file(taxonomy=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=PJB.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)"

mothur "#count.groups(shared=PJB.opti_mcc.shared)"
```

The count of sequences in each file are:
`less cluster_output_script`

```
mothur > count.groups(shared=PJB.opti_mcc.shared)
WSH129 contains 57046.
WSH130 contains 30636.
WSH131 contains 6973.
WSH132 contains 74507.
WSH133 contains 14419.
WSH134 contains 10908.
WSH135 contains 20866.
WSH136 contains 11094.
WSH137 contains 8789.
WSH138 contains 58396.
WSH139 contains 15638.
WSH140 contains 27886.
WSH141 contains 15379.
WSH142 contains 17475.
WSH143 contains 23156.
WSH144 contains 20043.
WSH145 contains 21560.
WSH146 contains 15054.
WSH147 contains 17104.
WSH148 contains 13757.
WSH149 contains 17375.
WSH150 contains 12317.
WSH151 contains 14510.
WSH152 contains 12731.
WSH153 contains 11634.
WSH154 contains 1691.
WSH155 contains 7213.
WSH156 contains 13272.
WSH157 contains 21281.
WSH158 contains 16355.
WSH159 contains 13519.
WSH160 contains 20729.
WSH161 contains 28550.
WSH162 contains 33867.
WSH163 contains 15613.
WSH164 contains 19120.
WSH165 contains 14342.
WSH166 contains 12461.
WSH167 contains 10946.
WSH168 contains 7393.
WSH169 contains 10084.
WSH170 contains 32468.
WSH171 contains 32391.
WSH172 contains 33225.
WSH173 contains 11343.

Size of smallest group: 1691.

Total seqs: 905116.

Output File Names:
PJB.opti_mcc.count.summary
```

11. Subsampling for Sequencing Depth

To have no samples removed from the dataset, I can only set the maximum subsampling depth to 1691.

```
interactive
module load Mothur/1.46.1-foss-2020b
mothur

sub.sample(shared=PJB.opti_mcc.shared, size=1691) #subsampling the dataset - PJB.opti_mcc.0.03.subsample.shared

rarefaction.single(shared=PJB.opti_mcc.shared, calc=sobs, freq=100) #rarefying the dataset - PJB.opti_mcc.groups.rarefaction

summary.single(shared=PJB.opti_mcc.shared, calc=nseqs-sobs-shannon-invsimpson, subsample=1691) # summarizing dataset - PJB.opti_mcc.groups.rarefaction
```

12. Caluclate Ecological Statistics

13. Population Level Analyses

14. Export for R analyses

Files needed to export from andromeda to desktop for R analysis:

* PJB.taxonomy
* PJB.opti_mcc.0.03.subsample.shared
* PJB.opti_mcc.shared

Run outside of andromeda.

`scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/mothur/PJB.taxonomy ~/MyProjects/Porites_Rim_Bleaching_2019/output/16S/mothur
`

`scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/mothur/PJB.opti_mcc.shared ~/MyProjects/Porites_Rim_Bleaching_2019/output/16S/mothur
`

`scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_16S/mothur/PJB.opti_mcc.0.03.subsample.shared ~/MyProjects/Porites_Rim_Bleaching_2019/output/16S/mothur
`
