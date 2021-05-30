#!/bin/bash

#SBATCH -c 2
#SBATCH --mem=5000
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH -p shared
#SBATCH --verbose
#SBATCH -J optsymsim
#SBATCH -o %x_%j_%N.out
#SBATCH -e %x_%j_%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=debug.pie@gmail.com

module load R/3.6.3-fasrc01
module load GCC/8.2.0-2.31.1

## nind: num of individual in one condition

Rscript symsim.R  \
        --nind_per_cond ${1} \
        --brn_len ${2} \
        --bimod ${3} \
        --sigma ${4} \
        --capt_alpha ${5} \
        --ngene ${6} \
        > symsim_${1}subpop_${2}brnlen_${3}bimod_${4}sigma_${5}alpha_${6}gene.log 2>&1
        
