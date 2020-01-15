#!/bin/bash
#PBS -V
#PBS -z
#PBS -l select=1:ncpus=8
cd $PBS_O_WORKDIR
abq=/scratch/users/erik/SIMULIA/CAE/2018/linux_a64/code/bin/ABQLauncher
${abq} j=overload cpus=8 interactive user=/scratch/users/erik/python_projects/transformation_fatigue/transformation_subroutine/transformation_subroutine.o

