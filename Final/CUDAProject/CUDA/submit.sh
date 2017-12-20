#!/bin/bash

#SBATCH -p slurm_shortgpu 

#SBATCH --job-name=solution

#SBATCH -N 1 -n 40 --gres=gpu:1

#SBATCH --error=/srv/home/cguo42/test.err

#SBATCH --output=/srv/home/cguo42/test.out

#SBATCH -o solution.o%j

cd $SLURM_SUBMIT_DIR

./sudokusolver16x16 16 16 16x16puzzle/puzzle16_01.txt puzzle16.txt
