#!/bin/bash
#SBATCH --nodes=1      
#SBATCH --ntasks=28
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --qos=6node_qos
#SBATCH --output=data/intfer/465/simu_2.%J.out
#SBATCH --error=data/intfer/465/simu_2.%J.err
    
#source ~/.bashrc
module swap intel gcc/9.2.0
module load pmix/2.2.2
source /etc/profile.d/conda.sh
conda activate work
python3 simu.py
