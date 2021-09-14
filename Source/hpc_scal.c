#include "hpc.h"

/* y = y + alpha * x */
index hpc_scal(double* y, const double *x, double alpha, index n)
{

    if ( !x || !y ) return (0) ;

    for(index i = 0; i < n; i++)
    {       
        y[i] += alpha * x[i] ;
    }
    return (1) ;
}