#!/bin/bash
module purge
module load GCC/8.3.0
module load mpich/mpi-x-gcc8.3.0
time yhrun -N 1 -n 5 -p thcp1 ./spmv

# yhbatch -N 1 -n 5 -p thcp1 ./spmv.sh