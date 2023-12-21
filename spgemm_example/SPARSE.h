#ifndef __SPARSE_hpp__
#define __SPARSE_hpp__

#include <stdio.h>
#include "floatVector.h"
#include "intVector.h"

typedef struct Scoo
{
	int
		m,	 // 行
		n,	 // 列
		nnz; // 非零元素数量
	IntVector
		RowPtr, // 大小m+1
		ColPtr, // 大小n+1
		RowInd, // 大小nnz，行坐标
		ColInd; // 大小nnz，列坐标
	FloatVector
		Val; // 大小为nnz
} Scoo;

void ScooFree(Scoo *a)
{
	IntVectorFree(&a->RowPtr);
	IntVectorFree(&a->ColPtr);
	IntVectorFree(&a->RowInd);
	IntVectorFree(&a->ColInd);
	FloatVectorFree(&a->Val);
}
void ScooMalloc(Scoo *a)
{
	IntVectorInit(&a->RowPtr, a->m + 1);
	IntVectorInit(&a->ColPtr, a->n + 1);
	IntVectorInit(&a->RowInd, a->nnz);
	IntVectorInit(&a->ColInd, a->nnz);
	FloatVectorInit(&a->Val, a->nnz);
}
void ScooPushBack(Scoo *a, int r, int c, float v)
{
	IntVectorPushBack(&a->RowInd, r);
	IntVectorPushBack(&a->ColInd, c);
	FloatVectorPushBack(&a->Val, v);
	++a->nnz;
}
void ScooRead(Scoo *a, FILE *f)
{
#define MAXLEN 511
	for (char s[MAXLEN]; fgets(s, MAXLEN, f) && !feof(f);)
		if (s[0] != '%')
		{
			int nnz;
			sscanf(s, "%d%d%d", &a->m, &a->n, &nnz);
			a->nnz = 0;
			ScooMalloc(a);
			for (int line = 0; line < nnz; ++line)
			{
				int r, c;
				float v;
				fscanf(f, "%d%d%f", &r, &c, &v);
				if (0 <= r && r < a->m &&
					0 <= c && c < a->n)
					ScooPushBack(a, r, c, v);
			}
			break;
		}
#undef MAXLEN
}
void ScooWrite(const Scoo *a, FILE *f)
{
	fprintf(
		f,
		"%%MatrixMarket matrix coordinate real general\n%d %d %d\n",
		a->m,
		a->n,
		a->nnz);
	for (int j = 0; j < a->nnz; ++j)
		fprintf(
			f,
			"%d %d %f\n",
			a->RowInd.data[j],
			a->ColInd.data[j],
			a->Val.data[j]);
}
void ScsrInit(Scoo *a)
{
	int cur = 0;
	for (int i = a->RowPtr.data[0] = 0; i < a->nnz; ++i)
		for (; cur < a->RowInd.data[i]; ++cur)
			a->RowPtr.data[cur + 1] = i;
	for (; cur < a->m; ++cur)
		a->RowPtr.data[cur + 1] = a->nnz;
}
void ScooToScsr(const Scoo *a, Scoo *b)
{
	IntVector beg, s;
	IntVectorInit(&beg, a->nnz);
	IntVectorInit(&s, 2);
	for (int i = 0; i < a->nnz; ++i)
		beg.data[i] = i;
	s.data[0] = 0, s.data[1] = a->nnz - 1;
	while (s.size)
	{
		int right = s.data[--s.size], left = s.data[--s.size], l = left, r = right, pivot = beg.data[l];
		while (l < r)
		{
#define LESS(x, y) (a->RowInd.data[x] < a->RowInd.data[y] || a->RowInd.data[x] == a->RowInd.data[y] && a->ColInd.data[x] < a->ColInd.data[y])
			while (l < r && !LESS(beg.data[r], pivot))
				--r;
			beg.data[l] = beg.data[r];
			while (l < r && !LESS(pivot, beg.data[l]))
				++l;
			beg.data[r] = beg.data[l];
#undef LESS
		}
		beg.data[l] = pivot;
		if (l - 1 > left)
		{
			IntVectorPushBack(&s, left);
			IntVectorPushBack(&s, l - 1);
		}
		if (l + 1 < right)
		{
			IntVectorPushBack(&s, l + 1);
			IntVectorPushBack(&s, right);
		}
	}
	IntVectorFree(&s);
	*b = *a;
	ScooMalloc(b);
	for (int i = 0; i < b->nnz; ++i)
	{
		b->RowInd.data[i] = a->RowInd.data[beg.data[i]];
		b->ColInd.data[i] = a->ColInd.data[beg.data[i]];
		b->Val.data[i] = a->Val.data[beg.data[i]];
	}
	ScsrInit(b);
}
void ScscInit(Scoo *a)
{
	int cur = 0;
	for (int i = a->ColPtr.data[0] = 0; i < a->nnz; ++i)
		for (; cur < a->ColInd.data[i]; ++cur)
			a->ColPtr.data[cur + 1] = i;
	for (; cur < a->m; ++cur)
		a->ColPtr.data[cur + 1] = a->nnz;
}
void ScooToScsc(const Scoo *a, Scoo *b)
{
	IntVector beg, s;
	IntVectorInit(&beg, a->nnz);
	IntVectorInit(&s, 2);
	for (int i = 0; i < a->nnz; ++i)
		beg.data[i] = i;
	s.data[0] = 0, s.data[1] = a->nnz - 1;
	while (s.size)
	{
		int right = s.data[--s.size], left = s.data[--s.size], l = left, r = right, pivot = beg.data[l];
		while (l < r)
		{
#define LESS(x, y) (a->ColInd.data[x] < a->ColInd.data[y] || a->ColInd.data[x] == a->ColInd.data[y] && a->RowInd.data[x] < a->RowInd.data[y])
			while (l < r && !LESS(beg.data[r], pivot))
				--r;
			beg.data[l] = beg.data[r];
			while (l < r && !LESS(pivot, beg.data[l]))
				++l;
			beg.data[r] = beg.data[l];
#undef LESS
		}
		beg.data[l] = pivot;
		if (l - 1 > left)
		{
			IntVectorPushBack(&s, left);
			IntVectorPushBack(&s, l - 1);
		}
		if (l + 1 < right)
		{
			IntVectorPushBack(&s, l + 1);
			IntVectorPushBack(&s, right);
		}
	}
	IntVectorFree(&s);
	*b = *a;
	ScooMalloc(b);
	for (int i = 0; i < b->nnz; ++i)
	{
		b->RowInd.data[i] = a->RowInd.data[beg.data[i]];
		b->ColInd.data[i] = a->ColInd.data[beg.data[i]];
		b->Val.data[i] = a->Val.data[beg.data[i]];
	}
	ScscInit(b);
	IntVectorFree(&beg);
}
#endif