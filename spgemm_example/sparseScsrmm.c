#include "SPARSE.h"

void ScsrMulScsc(const Scoo *a, const Scoo *b, Scoo *c)
{
	c->m = a->m;
	c->n = b->n;
	c->nnz = 0;
	ScooMalloc(c);
	for (int i = 0; i < c->m; ++i)
		if (a->RowPtr.data[i] < a->RowPtr.data[i + 1])
			for (int j = 0; j < c->n; ++j)
			{
				float res = 0;
				for (int
						 kb = b->ColPtr.data[j],
						 ka = a->RowPtr.data[i];
					 kb < b->ColPtr.data[j + 1] &&
					 ka < a->RowPtr.data[i + 1];
					 ++kb)
					if (b->Val.data[kb])
					{
						while (a->ColInd.data[ka] < b->RowInd.data[kb])
							++ka;
						while (a->ColInd.data[ka] == b->RowInd.data[kb])
							res += a->Val.data[ka] * b->Val.data[kb], ++ka;
					}
				if (res)
					ScooPushBack(c, i, j, res);
			}
	ScsrInit(c);
}
void ScsrMulScsr(const Scoo *a_csr, const Scoo *b_csr, Scoo *c_csr)
{
	Scoo B_csc;
	ScooToScsc(b_csr, &B_csc);
	ScsrMulScsc(a_csr, &B_csc, c_csr);
	ScooFree(&B_csc);
}
int main(int argc, char **argv)
{
	Scoo a_coo, b_coo, a_csr, b_csr, c_csr;

	ScooRead(&a_coo, fopen(argv[1], "r"));
	ScooRead(&b_coo, fopen(argv[2], "r"));

	ScooToScsr(&a_coo, &a_csr);
	ScooToScsr(&b_coo, &b_csr);

	ScooFree(&a_coo);
	ScooFree(&b_coo);

	ScsrMulScsr(&a_csr, &b_csr, &c_csr);

	ScooWrite(&c_csr, stdout);

	ScooFree(&a_csr);
	ScooFree(&b_csr);
	ScooFree(&c_csr);
}