#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8GB
#SBATCH -t 16:00:00
#SBATCH --gres=gpu:1
#SBATCH --dependency=singleton


host=$(env - hostname -s)
wrapper=/scratch/work/public/singularity/run-hoomd-2.9.6.bash

exe=run_dlvo_spheres_metal.py

$wrapper python -u /scratch/mc10050/work/ionic-crystal-project/nucleation2/${exe}
exit




