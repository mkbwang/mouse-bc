#!/bin/sh

#SBATCH --job-name=microarray_ensembl_annotation
#SBATCH --time=14:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=5g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x-error.out

module load Rtidyverse/4.4.0
source /home/wangmk/.bashrc
micromamba activate ML


python microarray_ensembl_annotation.py

