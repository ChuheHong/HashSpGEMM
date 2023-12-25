#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstdio>
#include <algorithm>
using namespace std;

// 自定义一个CSR存储元素
class CSRFloatElement
{
public:
    float values;  // 值
    int colIndex;  // 列索引
    int rowOffset; // 行偏移量
};
// 自定义一个COO存储元素
class COOFloatElement
{
public:
    float values; // 值
    int colIndex; // 列索引
    int rowIndex; // 行偏移量
};

// CSR格式矩阵
class CSRMatrix
{
public:
    int rows; // 行数
    int cols; // 列数
    int nnz;  // 非零项
    CSRFloatElement *csr_element;
};
// 创建CSR存储格式的矩阵
CSRMatrix CreateCSRMatrix(int rows, int cols, int nnz)
{
    CSRMatrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.nnz = nnz;
    // 分配nnz个空间，并初始化为0
    matrix.csr_element = new CSRFloatElement[nnz]();
    return matrix;
}

// COO格式矩阵
class COOMatrix
{
public:
    int rows; // 行数
    int cols; // 列数
    int nnz;  // 非零项
    COOFloatElement *coo_element;
};
// 创建COO存储格式的矩阵
COOMatrix CreateCOOMatrix(int rows, int cols, int nnz)
{
    COOMatrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.nnz = nnz;
    // 分配nnz个空间，并初始化为0
    matrix.coo_element = new COOFloatElement[nnz]();
    return matrix;
}

// 以COO的方式读取mtx文件的矩阵
COOMatrix ReadMatrix(FILE *f)
{
    COOMatrix matrix;
#define MAXLEN 511
    // 读取矩阵
    for (char s[MAXLEN]; fgets(s, MAXLEN, f) && !feof(f);)
        if (s[0] != '%')
        {
            int rows, cols, nnz;
            // 读取矩阵的行数、列数、非零项个数
            sscanf(s, "%d%d%d", &rows, &cols, &nnz);
            matrix = CreateCOOMatrix(rows, cols, nnz);
            for (int line = 0; line < nnz; ++line)
            {
                int r, c;
                float v;
                // 读取非零项的行索引、列索引和值
                fscanf(f, "%d%d%f", &r, &c, &v);
                // mtx矩阵的行数和列数是从1开始，所以要减去1
                matrix.coo_element[line].rowIndex = r - 1;
                matrix.coo_element[line].colIndex = c - 1;
                matrix.coo_element[line].values = v;
            }
            break;
        }
#undef MAXLEN
    return matrix;
}

// 将COO存储矩阵转换为CSR格式，COO存储格式应该按照从左到右、从上到下整理好
CSRMatrix COOToCSR(COOMatrix coo_matrix)
{
    CSRMatrix csr_matrix;
    csr_matrix = CreateCSRMatrix(coo_matrix.rows, coo_matrix.cols, coo_matrix.nnz);
    for (int i = 0; i < csr_matrix.nnz; i++)
    {
        csr_matrix.csr_element[i].values = coo_matrix.coo_element[i].values;
        csr_matrix.csr_element[i].colIndex = coo_matrix.coo_element[i].colIndex;
        // 根据coo的行，计算出每一行的非零项，并放到下一行
        csr_matrix.csr_element[coo_matrix.coo_element[i].rowIndex + 1].rowOffset++;
    }
    for (int i = 0; i < csr_matrix.rows; i++)
    {
        // 将前面累加起来
        csr_matrix.csr_element[i + 1].rowOffset += csr_matrix.csr_element[i].rowOffset;
    }
    return csr_matrix;
}

// 读取CSR矩阵第i行第j列元素值,i和j从0开始
float getCSRMatrixElement(const CSRMatrix matrix, int i, int j)
{
    int start = matrix.csr_element[i].rowOffset;
    int end = matrix.csr_element[i + 1].rowOffset;
    // 从这一行开始寻找，注意k不等于end
    for (int k = start; k < end; k++)
    {
        // 列索引相同则找到元素
        if (matrix.csr_element[k].colIndex == j)
        {
            return matrix.csr_element[k].values;
        }
    }
    return 0.0;
}

// 针对COO自定义一个排序比较函数
bool compareCOOMatrix(COOFloatElement e1, COOFloatElement e2)
{
    if (e1.rowIndex != e2.rowIndex)
        return e1.rowIndex < e2.rowIndex;
    return e1.colIndex < e2.colIndex;
}

// 对COO矩阵按照行优先的方式排序
void sortCOOMatrix(COOMatrix coo_matrix)
{
    sort(coo_matrix.coo_element, coo_matrix.coo_element + coo_matrix.nnz, compareCOOMatrix);
}

#endif