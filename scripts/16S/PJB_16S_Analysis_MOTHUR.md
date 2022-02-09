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

mothur "#make.contigs(inputdir=., outputdir=., file=PJB.files, oligos=oligos.oligos)"

mothur "#summary.seqs(fasta=PJB.trim.contigs.fasta)"
```

Summary output: `less contigs_output_script`

```
Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       22      22      0       2       1
2.5%-tile:      1       292     292     0       4       31747
25%-tile:       1       292     292     0       4       317470
Median:         1       293     293     0       4       634939
75%-tile:       1       301     301     0       5       952408
97.5%-tile:     1       310     310     6       6       1238130
Maximum:        1       563     563     114     270     1269876
Mean:   1       297     297     0       4
# of Seqs:      1269876
```

Check that primers are gone:
* F GTGCCAGCMGCCGCGGTAA
* R GGACTACNVGGGTWTCTAAT

`head PJB.trim.contigs.fasta`

```
>M00763_26_000000000-K4TML_1_1101_13046_1536	ee=1.77103	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match)
GGATAGGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGATACCCCGGTAGTCC
>M00763_26_000000000-K4TML_1_1101_11503_1634	ee=2.58733	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match)
TGACAGGTGCCAGCCGCCGCGGGAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGCTACCCGTGTAGTCC
>M00763_26_000000000-K4TML_1_1101_17805_1647	ee=1.04918	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match)
GTGCCAGCAGCCGCGGTAAGACGGAGGATGCAAGTGTTATCCGGAATCACTGGGCGTAAAGCGTCTGTAGGTGGTTTAATAAGTCAACTGTTAAATCTTAAGGCTCAACCTTAAAATCGCAGTCGAAACTATTAAACTAGAGTATAGTAGAGGTAAAGGGAATTTCCAGTGGAGCGGTGAAATGCGTAGAGATTGGAAAGAACACCGATGGCGAAGGCACTTTACTGGGCTATTACTAACACTCAGAGACGAAAGCTAGGGTAGCAAATGGGATTAGATACCCTAGTAGTCC
>M00763_26_000000000-K4TML_1_1101_11512_1651	ee=2.09099	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match)
GGGACCGTGGCCAGCCGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTTTGTAGGCTACTTTTTAAGTCTGCTGTTAAAGAATAAGACTTAACCTTAGAAAAGCAGTATGAAACTAAAAGGATAGAGTTCAGTAGGGGTAGAGGGAATTCTCGGTGTAGTAGTGAAATGCGTAGATATCGAGAGGAACACCAATAGCGAAAGCACTCTACTAGGCTGTAACTGACGCTAAGAAACGAAAGCTAGGGGAGCCAATGGGATTAGCTACCCGTGTAGTCC
>M00763_26_000000000-K4TML_1_1101_21801_1736	ee=2.29141	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match)
TGGACAGGTGCCAGCCGCCGCGTGAAGACGGAGGATGCTAGTGTTATCCGGAATCACTGGGCGTAAAGCGTCTGTAGNTGGTTAAATAAGTCAACTGTTAAATCTTGAAGCTCAACTTCAAAATCGCAATCGAAACTATTTGACTTGAGTATAGTAGAGGTAACGGGAATTTCCAGTGGAGTGGTGAAATGCGTAGAGATTGGAAAGAACACCGATGGCGAAGGCACTTTACTGGGCTATTACTAACACTCAGAGACGAAAGCTAGGGTAGCAAATGGGATTAGATAACCCGGTAGTCC
```

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
Minimum:        1       252     252     0       3       1
2.5%-tile:      1       292     292     0       4       24702
25%-tile:       1       292     292     0       4       247017
Median:         1       293     293     0       4       494033
75%-tile:       1       301     301     0       5       741049
97.5%-tile:     1       310     310     0       6       963364
Maximum:        1       350     350     0       27      988065
Mean:   1       295     295     0       4
# of Seqs:      988065
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
Minimum:        1       252     252     0       3       1
2.5%-tile:      1       292     292     0       4       24702
25%-tile:       1       292     292     0       4       247017
Median:         1       293     293     0       4       494033
75%-tile:       1       301     301     0       5       741049
97.5%-tile:     1       310     310     0       6       963364
Maximum:        1       350     350     0       27      988065
Mean:   1       295     295     0       4
# of unique seqs:       579759
total # of seqs:        988065
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
Minimum:        0       0       0       0       1       1
2.5%-tile:      1       13424   292     0       4       14494
25%-tile:       1       13424   292     0       4       144940
Median:         1       13424   292     0       4       289880
75%-tile:       1       13424   293     0       5       434820
97.5%-tile:     1       13425   294     0       6       565266
Maximum:        13424   13425   302     0       19      579759
Mean:   5       13419   292     0       4
# of Seqs:      579759
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
Minimum:        1       11550   269     0       3       1
2.5%-tile:      1       13424   292     0       4       24681
25%-tile:       1       13424   292     0       4       246801
Median:         1       13424   292     0       4       493602
75%-tile:       1       13424   293     0       5       740403
97.5%-tile:     1       13425   293     0       6       962523
Maximum:        1968    13425   300     0       8       987203
Mean:   1       13424   292     0       4
# of unique seqs:       579013
total # of seqs:        987203
```

```
mothur > count.groups(count=PJB.trim.contigs.good.good.count_table)
WSH129 contains 61316.
WSH130 contains 32627.
WSH131 contains 7086.
WSH132 contains 83611.
WSH133 contains 14747.
WSH134 contains 11080.
WSH135 contains 25086.
WSH136 contains 12607.
WSH137 contains 8889.
WSH138 contains 68185.
WSH139 contains 16473.
WSH140 contains 29570.
WSH141 contains 16859.
WSH142 contains 18954.
WSH143 contains 25182.
WSH144 contains 20997.
WSH145 contains 25124.
WSH146 contains 19133.
WSH147 contains 17608.
WSH148 contains 16778.
WSH149 contains 19309.
WSH150 contains 12999.
WSH151 contains 19308.
WSH152 contains 13562.
WSH153 contains 12015.
WSH154 contains 1731.
WSH155 contains 7269.
WSH156 contains 14588.
WSH157 contains 21940.
WSH158 contains 16712.
WSH159 contains 13891.
WSH160 contains 22959.
WSH161 contains 31522.
WSH162 contains 37154.
WSH163 contains 16460.
WSH164 contains 21284.
WSH165 contains 14706.
WSH166 contains 12564.
WSH167 contains 11062.
WSH168 contains 7839.
WSH169 contains 10634.
WSH170 contains 33551.
WSH171 contains 35930.
WSH172 contains 34627.
WSH173 contains 11675.

Size of smallest group: 1731.

Total seqs: 987203.

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
Minimum:        1       546     230     0       3       1
2.5%-tile:      1       547     253     0       3       24681
25%-tile:       1       547     253     0       4       246801
Median:         1       547     253     0       4       493602
75%-tile:       1       547     253     0       5       740403
97.5%-tile:     1       547     254     0       6       962523
Maximum:        4       547     271     0       8       987203
Mean:   1       546     253     0       4
# of unique seqs:       579013
total # of seqs:        987203
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
Minimum:        1       546     230     0       3       1
2.5%-tile:      1       547     253     0       3       24681
25%-tile:       1       547     253     0       4       246801
Median:         1       547     253     0       4       493602
75%-tile:       1       547     253     0       5       740403
97.5%-tile:     1       547     254     0       6       962523
Maximum:        4       547     271     0       8       987203
Mean:   1       546     253     0       4
# of unique seqs:       45194
total # of seqs:        987203
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
Minimum:        1       546     230     0       3       1
2.5%-tile:      1       547     253     0       3       24025
25%-tile:       1       547     253     0       4       240245
Median:         1       547     253     0       4       480489
75%-tile:       1       547     253     0       5       720733
97.5%-tile:     1       547     254     0       6       936953
Maximum:        4       547     271     0       8       960977
Mean:   1       546     253     0       4
# of unique seqs:       34016
total # of seqs:        960977
```

The count of sequences in each file are:

```
WSH136 contains 12335.
WSH137 contains 8858.
WSH138 contains 66295.
WSH139 contains 16281.
WSH140 contains 28837.
WSH141 contains 16184.
WSH142 contains 18058.
WSH143 contains 24014.
WSH144 contains 20191.
WSH145 contains 24021.
WSH146 contains 18695.
WSH147 contains 17150.
WSH148 contains 16534.
WSH149 contains 18714.
WSH150 contains 12731.
WSH151 contains 18946.
WSH152 contains 13264.
WSH153 contains 11885.
WSH154 contains 1728.
WSH155 contains 7223.
WSH156 contains 14442.
WSH157 contains 21383.
WSH158 contains 16393.
WSH159 contains 13701.
WSH160 contains 21881.
WSH161 contains 30710.
WSH162 contains 35950.
WSH163 contains 16051.
WSH164 contains 20936.
WSH165 contains 14481.
WSH166 contains 12496.
WSH167 contains 10965.
WSH168 contains 7668.
WSH169 contains 10542.
WSH170 contains 32868.
WSH171 contains 35009.
WSH172 contains 33486.
WSH173 contains 11497.

Size of smallest group: 1728.
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
Minimum:        1       546     230     0       3       1
2.5%-tile:      1       547     253     0       3       22686
25%-tile:       1       547     253     0       4       226854
Median:         1       547     253     0       4       453708
75%-tile:       1       547     253     0       5       680562
97.5%-tile:     1       547     254     0       6       884730
Maximum:        4       547     271     0       8       907415
Mean:   1       546     253     0       4
# of unique seqs:       32032
total # of seqs:        907415

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
WSH129 contains 57170.
WSH130 contains 30695.
WSH131 contains 6976.
WSH132 contains 74729.
WSH133 contains 14423.
WSH134 contains 10908.
WSH135 contains 20916.
WSH136 contains 11099.
WSH137 contains 8787.
WSH138 contains 59116.
WSH139 contains 15635.
WSH140 contains 27905.
WSH141 contains 15699.
WSH142 contains 17488.
WSH143 contains 23388.
WSH144 contains 20055.
WSH145 contains 21605.
WSH146 contains 15075.
WSH147 contains 17100.
WSH148 contains 13768.
WSH149 contains 17393.
WSH150 contains 12323.
WSH151 contains 14537.
WSH152 contains 12753.
WSH153 contains 11733.
WSH154 contains 1691.
WSH155 contains 7214.
WSH156 contains 13286.
WSH157 contains 21283.
WSH158 contains 16375.
WSH159 contains 13530.
WSH160 contains 20767.
WSH161 contains 28585.
WSH162 contains 33915.
WSH163 contains 15618.
WSH164 contains 19150.
WSH165 contains 14346.
WSH166 contains 12463.
WSH167 contains 10952.
WSH168 contains 7394.
WSH169 contains 10086.
WSH170 contains 32485.
WSH171 contains 32399.
WSH172 contains 33256.
WSH173 contains 11344.

Size of smallest group: 1691.

Total seqs: 907415.

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

14. Export foor R analyses

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
