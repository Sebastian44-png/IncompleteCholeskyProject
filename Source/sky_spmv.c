#include "hpc.h"
/* y = A * x + y */
index sky_spmv (const sky *A, const double *x, double *y)
{
	index k, n, *Ap, row;
	double *Ad, *Ax;
	if(!A) return (0);
	n = A->n; Ap = A->p; Ad = A->d; Ax = A->x;
	// Diagonal
	for(k = 0; k < n; ++k){
		y[k] += Ad[k]*x[k];
	}
	
	// other entries
	for(row = 1; row < n; ++row){
		for(k = Ap[row-1]; k < Ap[row]; ++k){
			y[row] += Ax[k]*x[row - Ap[row] + k];
			y[row - Ap[row] + k] += Ax[k]*x[row];
		}
	}

  	return (1) ;
}
