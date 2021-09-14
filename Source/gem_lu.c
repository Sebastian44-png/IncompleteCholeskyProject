#include "hpc.h"
/* perform a LU decomposition without pivoting on matrix stored 
 * in general matrix format (row wise), the lower triangular 
 * part a A gets overwritten with L */

index gem_lu(gem *A) 
{
    index i;
    index j;
    index k;
    index n;
    double  recip;
    double *Ax;

    if ( (!A) || (A->m != A->n)) { 
        printf ("(null)\n") ; 
        return (0) ; 
    }
    n = A->n; Ax = A->x;

    for (k = 0; k < n ; k++){

        /* modify entries in A, i.e. create LU decomp*/
        for (i = k+1; i < n; i++){  
            /* compute  l_ik, store in A */
            Ax[n*i+k] = Ax[n*i+k] / Ax[k*n+k];
            
            for (j = k+1; j < n; j++){
                /* mofify (n-k-1)*(n-k-1) submatrix */
                Ax[i*n+j] -= Ax[i*n+k] * Ax[k*n+j];
            }
        }
    }
    return(1);
}