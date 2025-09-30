#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<youremailaddress>
#SBATCH --output=%x_%j.out

# ───────────────────────────────────────────────
# 1) Load modules
module load nextflow/23.10.1       # or whichever version is available
module load singularity/3.11.0     # required for nf-core containers

# ───────────────────────────────────────────────
# 2) Move to work directory
cd $SCRATCH/RNAseq_project

# ───────────────────────────────────────────────
# 3) Run nf-core/rnaseq v3.10.1
nextflow run nf-core/rnaseq -r 3.10.1 \
    -profile singularity \                
    --input samplesheet.csv \            
    --genome GRCh38 \                     
    --aligner star_rsem \                  
    --outdir results_rnaseq \              
