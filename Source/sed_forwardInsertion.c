#include "hpc.h"

/*Forward Insecertion: solves the linear System Lx=b where L is a lower triangular matrix*/

index sed_forwardInsertion(sed *L , double *x , double *b)
{
    if(!L || !x || !b)
    {
        return(0);
    }
    
    index n ;
    index *ind ;
    double *val ;

    n = L->n ;
    ind = L->i ;
    val = L->x ;
   
    x [0] = b [0] / val [0] ;
    
    for (index i = 1 ; i < n ; i++)
    {
        for (index ptr = ind[i - 1] ; ptr < ind[i] ; ptr ++)
        {
            x [ind [ptr]] += val [ptr] * x [i -1] ; 
        }

        x [i] = (b [i] - x [i]) / val [i] ;
    }
    return (1) ;
}
