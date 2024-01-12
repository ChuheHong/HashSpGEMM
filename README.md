# HashSpGEMM
It implements the general sparse matrix-matrix via Hash accumulator.
## Complie
To complie the code you should confirm the path of the "omp.h" and change it. For example, the path of "omp.h" is `/opt/homebrew/opt/libomp/include/omp.h` in my compter.

Run the following code to complie the code.
```
cd obj
gcc-13 -fopenmp -O3 -c ../src/*.c
gcc-13 -fopenmp -O3 ./* -o ../bin/main
```
The version of my gcc is 13, please change the version according to your own situasion.

## Run
```
../bin/main
```
The nums of thread are determined by your own CPU cores.