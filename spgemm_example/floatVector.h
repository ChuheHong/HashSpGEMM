#ifndef __floatVector_h__
#define __floatVector_h__

// #include <malloc.h>
#include <stdlib.h>
#include <sys/malloc.h>
#include <string.h>

typedef struct FloatVector
{
    int size, cap;
    float *data;
} FloatVector;

void FloatVectorFree(FloatVector *v) { free(v->data); }
void FloatVectorInit(FloatVector *v, int size)
{
    v->size = v->cap = size;
    if (v->cap < 1)
        v->cap = 1;
    v->data = (float *)malloc(sizeof(float) * v->cap);
}
void FloatVectorPushBack(FloatVector *v, float val)
{
    while (v->cap <= v->size)
    {
        FloatVector v_new;
        FloatVectorInit(&v_new, v->cap << 1);
        memcpy(v_new.data, v->data, sizeof(float) * v->size);
        FloatVectorFree(v);
        *v = v_new;
        v->size >>= 1;
    }
    v->data[v->size++] = val;
}
#endif