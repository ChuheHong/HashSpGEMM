#include <cstdio>
#include <algorithm>
#include "Matrix.hpp"
using namespace std;

int main(int argc, char *argv[])
{
    COOMatrix matrix = ReadMatrix(fopen("./494_bus.mtx", "r"));
    sortCOOMatrix(matrix);
    for (int i = 0; i < matrix.nnz; i++)
    {
        printf("%d %d %f\n",
               matrix.coo_element[i].rowIndex,
               matrix.coo_element[i].colIndex,
               matrix.coo_element[i].values);
    }
    return 0;
}
