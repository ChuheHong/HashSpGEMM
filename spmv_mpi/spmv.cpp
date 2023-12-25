#include <cstdio>
#include "Matrix.hpp"
#include "mpi.h"
using namespace std;

int myid, m_id, numprocs;
COOMatrix coo_matrix;
CSRMatrix csr_matrix;
float *vector, *result, *buffer;
MPI_Status status;
// 主进程
void master()
{
    // 接收从进程发送来的buffer
    for (int i = 0; i < numprocs - 1; i++)
    {
        MPI_Recv(buffer,
                 coo_matrix.cols,
                 MPI_FLOAT,
                 MPI_ANY_SOURCE,
                 MPI_ANY_TAG,
                 MPI_COMM_WORLD,
                 &status);
        // status.MPI_TAG表示从进程号，只对从进程处理的行累加
        for (int j = status.MPI_TAG; j < coo_matrix.cols; j += (numprocs - 1))
        {
            result[j] += buffer[j];
        }
    }
}
// 从进程
void slave()
{
    // 并行计算，主进程不计算
    for (int i = myid; i < coo_matrix.cols; i += (numprocs - 1))
    {
        // 遍历csr矩阵的第i行
        int start = csr_matrix.csr_element[i].rowOffset;
        int end = csr_matrix.csr_element[i + 1].rowOffset;
        // k代表第k个元素
        for (int k = start; k < end; k++)
        {
            // 缓冲向量第i行的值为csr矩阵对应的元素乘向量vector的元素，索引为csr矩阵所在的列
            buffer[i] = csr_matrix.csr_element[k].values * (vector[csr_matrix.csr_element[k].colIndex]);
        }
    }

    // 将buffer缓冲值发送到主进程，tag为从进程号
    MPI_Send(buffer,
             coo_matrix.cols,
             MPI_FLOAT,
             m_id, myid,
             MPI_COMM_WORLD);
}
int main(int argc, char *argv[])
{
    // MPI初始化，myid是当前进程id，m_id是主进程id，numprocs是总进程数
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    // 设置主进程
    m_id = numprocs - 1;
    // 读取COO格式存储的mtx文件
    coo_matrix = ReadMatrix(fopen("./494_bus.mtx", "r"));
    // 对COO矩阵排序
    sortCOOMatrix(coo_matrix);
    // 根据矩阵的列数，生成向量
    vector = new float[coo_matrix.cols];
    result = new float[coo_matrix.cols]();
    buffer = new float[coo_matrix.cols]();
    for (int i = 0; i < coo_matrix.cols; i++)
        vector[i] = i + 1;

    // 把COO格式转换为CSR存储格式
    csr_matrix = COOToCSR(coo_matrix);

    // 执行主进程和从进程
    if (myid == m_id)
        master();
    else
        slave();

    // 打印结果
    if (myid == m_id)
    {
        for (int i = 0; i < coo_matrix.cols; i++)
            printf("%d %f\n", i, result[i]);
    }
    MPI_Finalize();
    return 0;
}
