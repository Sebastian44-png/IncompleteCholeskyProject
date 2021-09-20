#include "hpc.h"

/* Gauss-Seidel Iteration with constrains x(fixed) = b(fixed) */
index sed_gs (const sed *A, const double *b, double *x, double *w, index forward) 
{
    index n ; 
    index *Ap ; 
    index *Ai ;
    double *Ax ;
    
    if ( !A || !x || !b ) return (0) ;  /* check inputs */
    n = A->n ; Ai = A->i ; Ax = A->x ;
    for (index j = 0 ; j < n; j++) w[j] = b[j] ;
    
    if (forward) 
    {
        for (index j = 0 ; j < n; j++)
        {
            for (index p = Ai[j] ; p < Ai[j+1] ; p++) 
            {
                w[j] -= Ax[p] * x [Ai[p]];
            }
        }
        for (index j = 0 ; j < n; j++)
        {
            x[j] = w[j] / Ax[j];
            for (index p = Ai[j] ; p < Ai[j+1] ; p++) 
            {
                w[Ai[p]] -= Ax[p] * x [j] ;    
            }
        }
    }
    else
    {
        for (index j = 0 ; j < n ; j++)
        {
            for (index p = Ai[j] ; p < Ai[j+1] ; p++) 
            {
                w[Ai[p]] -= Ax[p] * x [j] ;
            }
        }
        for (index j = n-1 ; j >=0; j--)
        {
            for (index p = Ai[j] ; p < Ai[j+1] ; p++)
            {
                w[j] -= Ax[p] * x [Ai[p]] ;
            }
            x[j] = w[j] / Ax[j] ;
        }
    }    
    return (1) ;
}
