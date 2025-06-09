#!/bin/bash
#SBATCH --job-name=Kobza_Code
#SBATCH --output=dada2_output_%j.log
#SBATCH --error=dada2_error_%j.log
#SBATCH --partition=guest,batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=echandrasekara2@huskers.unl.edu
#SBATCH --licenses=common
#SBATCH --cpus-per-task=8
#SBATCH --time=05:59:00
#SBATCH --mem=32G

module load R

# Go to working directory
cd /work/samodha/sachin/Kobza

# Run R script
Rscript Kobza_Code.R	
