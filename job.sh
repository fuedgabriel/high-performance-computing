#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=res.txt
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --time=05:00
#SBATCH -p debug


mpirun ./a.out