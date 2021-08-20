#include "hpc.h"

// One step of Gauss-Seidel Iteration. Solve: Ax = b, w as auxiliary vector
index sed_gauss_seidel(const sed *A, const double *b, double *xk, double *w)
{
    index n, *ind; 
    double *Ax; 
    

    n = A->n; ind = A->i; Ax = A->x;

    //copy xk to w
    for(index i = 0; i<n; i++){
        w[i] = xk[i];
    }

    //iterate through rows and cols of matrix
    for(index k = 0; k < n; k++){
        xk[k] = b[k];
    
        for(index j = 0; j < n; j++){
            if(k == j) continue;

            // find row k in col j
            for(index i = ind[j]; i < ind[j+1]; i++){
                if(ind[i] == k){

                    xk[k] -= xk[j] * Ax[i];
                }
            }  
        }
        xk[k] = xk[k] / Ax[k];
    }
    return 1;
}
