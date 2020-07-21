#!/bin/bash

##SBATCH --partition=maintenance

# Request an hour of runtime:
#SBATCH --time=48:30:00

# Use 2 nodes with 64 tasks each, for 8 MPI tasks:
#SBATCH --nodes=8
#SBATCH --tasks-per-node=9
##SBATCH --nodelist=node1129

#SBATCH --mem=32G

# Specify a job name:
#SBATCH -J cyl3_5 

# Specify an output file
#SBATCH -o Re11K.out
#SBATCH -e Re11K.out

# Run a command
 module load mpi/mvapich2-2.3a_intel
 module load python/3.5.2
 module load cuda/10.0.130
 module load cudnn/7.4
 module load tensorflow/1.13.1_cpu_py3

 python3 -u server.py -env CFD -fil None >& python.log &

 srun -n 64 --mpi=pmi2 ./nektarF -z128 -deal -chk -r0 -S -V cyl.rea >& nektar.log &

 wait
