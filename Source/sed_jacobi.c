#include "hpc.h"
// One step of Jacobi Iteration. Solve: Ax = b, w as auxiliary vector
index sed_jacobi (const sed *A, const double *b, double *xk, double *w )
{
    index n, *ind;   
    double *Ax;

    n = A->n; ind = A->i; Ax = A->x;
    
    for(size_t i = 0; i < n; i++){
        w[i] = b[i];
    } 
        
    //iterate through colums: xk - aij * xj für i != j
    for(size_t j = 0; j<n; j++){
        for(size_t i = ind[j]; i< ind[j+1]; i++){
            w[ind[i]] -= Ax[i]*xk[j];
        }
    }

    // xk = xk/akk
    for(size_t k=0; k < n; k++){
        xk[k] = w[k] / Ax[k];
    }
      
    return 1;
}   

