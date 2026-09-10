#!/bin/bash
#SBATCH --job-name=RM_LTF
#SBATCH --partition=regular
#SBATCH --cpus-per-task=16
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=s.a.gonzalez.mollinedo@rug.nl

module load RepeatModeler
module load RepeatMasker
module load BLAST+
module load Miniconda3/22.11.1-1

SCRATCH_DIR=/scratch/$USER
ENV_PREFIX="$SCRATCH_DIR/conda-envs/earlgrey-7.2.4"
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PREFIX"

earlGrey -g LTFr_final_renamed_hap1.fa -t 16 -s PoHecki -o ./earlGreyOutputs -m yes

