#include "hpc.h"
index jds_spmv (const jds *A, const double *x, double *y)
{
	index Andiag, *Ap, *Aj, *Aperm, k, p;
	double *Ax;
	if(!A) return (0);
	Andiag = A->ndiag; Ap = A->p; Aj = A->j; Aperm = A->perm; Ax = A->x;
	for(k = 0; k < Andiag; ++k){
		for(p = Ap[k]; p < Ap[k + 1]; ++p){
			y[Aperm[p-Ap[k]]] += Ax[p]*x[Aj[p]];
		}
	}
	return (1) ;
}
