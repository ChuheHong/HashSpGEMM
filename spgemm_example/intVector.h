#ifndef __intVector_h__
#define __intVector_h__

#include <sys/malloc.h>
#include <string.h>
#include <stdlib.h>

typedef struct IntVector
{
    int size, cap;
    int *data;
} IntVector;

void IntVectorFree(IntVector *v) { free(v->data); }
void IntVectorInit(IntVector *v, int size)
{
    v->size = v->cap = size;
    if (v->cap < 1)
        v->cap = 1;
    v->data = (int *)malloc(sizeof(int) * v->cap);
}
void IntVectorPushBack(IntVector *v, int val)
{
    while (v->cap <= v->size)
    {
        IntVector v_new;
        IntVectorInit(&v_new, v->cap << 1);
        memcpy(v_new.data, v->data, sizeof(int) * v->size);
        IntVectorFree(v);
        *v = v_new;
        v->size >>= 1;
    }
    v->data[v->size++] = val;
}
#endif