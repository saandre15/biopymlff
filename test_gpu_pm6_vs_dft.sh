#!/bin/bash

#SBATCH -J test-1           # Job name
#SBATCH -o test-1.o%j       # Name of stdout output file
#SBATCH -e test-1.e%j       # Name of stderr error file
#SBATCH -p small           # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for OpenMP)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for OpenMP)
#SBATCH -t 48:00:00       # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=axs210186@utdallas.edu


export OMP_NUM_THREADS=56
module use /work2/01114/jfonner/frontera/modulefiles
module load gaussian
cd $WORK/biopymlff
pipenv run python3 -m unittest biopymlff/calculators/gebf_pm6_test.py
