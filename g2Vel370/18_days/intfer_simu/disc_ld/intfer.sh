#!/bin/bash
#SBATCH --nodes=1      
#SBATCH --ntasks=28
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --output=data/sh/intferwrl_2.%J.out
#SBATCH --error=data/sh/intferwrl_2.%J.err
    
#source ~/.bashrc
module load gcc/9.2.0
conda activate work
module load pmix/2.2.2
python3 intfer.py

