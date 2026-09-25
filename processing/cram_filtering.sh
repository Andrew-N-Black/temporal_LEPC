#!/bin/bash
#SBATCH --job-name=dedup_old_vs_new
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail
ml biocontainers
ml samtools

# For CRAM files (sarek default output format)
for f in ./*md.cram; do
    base=$(basename "$f" .cram)
    samtools view -@ 10 -C -F 0x400 -q 20 -T /scratch/gautschi/blackan/GROUSE/old_vs_new/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna \
        -o "./${base}.dedup_q20.cram" "$f"
    samtools index "./${base}.dedup_q20.cram"
done
