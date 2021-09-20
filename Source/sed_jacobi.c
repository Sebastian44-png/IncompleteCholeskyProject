#include "hpc.h"
/* author: Benjamin Bestler */

// One step of Jacobi Iteration for lower half of symmetric matrix. Solve: Ax = b, w as auxiliary vector
index sed_jacobi (const sed *A, const double *b, double *xk, double *w )
{
    // check inputs
    if(!A || !b || !xk || !w)
    {
        return 0;
    }
    index n, *ind ;   
    double *Ax ;

    n = A->n ; 
    ind = A->i ; 
    Ax = A->x ;
    
    for(size_t i = 0; i < n; i++)
    {
        w[i] = b[i] ;
    } 
        
    //iterate through colums: xk - aij * xj for i != j
    for(size_t j = 0 ; j<n ; j++)
    {
        for(size_t i = ind[j] ; i< ind[j+1] ; i++)
        {
            w[ind[i]] -= Ax[i]*xk[j] ;
            w[j] -= Ax[i] * xk[ind[i]] ;
        }
    }

    // xk = xk/akk
    for(size_t k=0 ; k < n ; k++)
    {
        xk[k] = w[k] / Ax[k] ;
    }
      
    return 1 ;
}   

