# PJB ITS2

```
mkdir PJB_ITS2

cd PJB_ITS2

mkdir raw_reads
```

files live here:

`scp /data/putnamlab/KITT/hputnam/20211210_AmpliconSeq/AmpSeq/KW_PJB_ITS2/* .`

### Symportal Database set up

The database is already set up from the [Thermal Transplant ITS2](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/Symportal_ThermalTransplant.md) for KW.

#### Loading Data

```bash
scp Wong_Kevin_PJB_ITS2_Meta.csv kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/PJB_ITS2/Wong_Kevin_PJB_ITS2_Meta.csv
```

`nano sp_load.sh`

```bash
#!/bin/bash
#SBATCH --job-name="SP_load"
#SBATCH -t 500:00:00
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/SymPortal
#SBATCH --exclusive

module load Miniconda3/4.9.2

eval "$(conda shell.bash hook)"
conda activate symportal_env

module load SymPortal/0.3.21-foss-2020b
module unload SymPortal/0.3.21-foss-2020b

export PYTHONPATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/lib/python3.7/site-packages:$PYTHONPATH

export PATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/bin:$PATH

main.py --load /data/putnamlab/kevin_wong1/PJB_ITS2/raw_reads \
--name PJB_1 \
--num_proc $SLURM_CPUS_ON_NODE \
--data_sheet /data/putnamlab/kevin_wong1/PJB_ITS2/Wong_Kevin_PJB_ITS2_Meta.csv

# Checking dataset number
./main.py --display_data_sets

```

#### Running analysis


```bash
#!/bin/bash
#SBATCH --job-name="SP_analysis"
#SBATCH -t 500:00:00
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/SymPortal
#SBATCH --exclusive

module load Miniconda3/4.9.2

eval "$(conda shell.bash hook)"
conda activate symportal_env

module load SymPortal/0.3.21-foss-2020b
module unload SymPortal/0.3.21-foss-2020b

export PYTHONPATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/lib/python3.7/site-packages:$PYTHONPATH

export PATH=/data/putnamlab/kevin_wong1/SymPortal/:/data/putnamlab/kevin_wong1/SymPortal/bin:$PATH

# Running analysis
./main.py --analyse 9 --name PJB_analysis --num_proc $SLURM_CPUS_ON_NODE

# Checking data analysis instances
./main.py --display_analyses

```

```bash
scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/SymPortal/outputs/analyses/3/20211216T114913/its2_type_profiles /Users/kevinwong/MyProjects/Thermal_Transplant_Molecular/output/ITS2/.
``
