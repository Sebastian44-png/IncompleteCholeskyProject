#include "hpc.h"

/* Jacobi Iteration with constrains x(fixed) = b(fixed) */
index sed_jacobi_constr(const sed *A, const double *b, double *x, double *w, 
                     index *fixed, index nFixed)
{
    index m;
    index n;
    index t;
    index ft;
    index *Ai;
    double *Ax;
    
    if ( !A || !x || !b ) return (0);  /* check inputs */
    n = A->n; 
    Ai = A->i;
    Ax = A->x;

    for (index j = 0; j < n; j++) 
    {
        w[j] = b[j];
    }
     
    //iterate through colums: xk - aij * xj für i != j
    for(index j = 0; j < n; j++){
        for(index i = Ai[j]; i < Ai[j+1]; i++)
        {
            w[Ai[i]] -= Ax[i] * x[j];
            w[j] -= Ax[i] * x[Ai[i]];
        }
    }

    // xk = xk/akk
    ft = fixed[0];
    t = 1;
    for(index k=0; k < n; k++)
    {
        if(k != ft)
        {
            x[k] = w[k] / Ax[k];
        }
        else
        {
            if(t < nFixed)
            {
                ft = fixed[t];
                t++;
            }
        }
    }

    return (1);
}