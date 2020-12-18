#!/bin/bash

#SBATCH -c 4
#SBATCH --mem=10000
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -p shared
#SBATCH --verbose
#SBATCH -J mssc_symsim
#SBATCH -o %x_%j_%N.out
#SBATCH -e %x_%j_%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=debug.pie@gmail.com

module load R/3.6.3-fasrc01
module load GCC/8.2.0-2.31.1

## nind: num of individual in one condition

Rscript symsim.R \
        --nind ${1} \
        --nindeff ${2} \
        --groupshift ${3} \
        --rpt ${4} \
        --ngene ${5} \
        --ratio_ind2cond ${6} \
        --scale_in_diffg ${7} \
        --scale_in_nondiffg ${8} \
        > mssc_symsim_${1}ind_${2}nindeff_${3}-gs_${5}g_ric-${6}_sid-${7}_sin-${8}.log 2>&1
        
