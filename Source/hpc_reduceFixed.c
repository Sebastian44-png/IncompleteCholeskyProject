#include "hpc.h"

/* writes non-fixed entries from x into xShort */
index hpc_reduceFixed (double *x, double *xShort, index n, const index *fixed, const index nFixed)
{
    index kShort ;
    index ft ;
    index fptr ;

    if (!x || !xShort || !fixed)
    {
        return (0) ;
    }

    kShort = 0 ;
    fptr = 0 ;
    ft = fixed [fptr] ;
    for (index k = 0 ; k < n ; k++)
    {
        if (k == ft)
        {
            fptr++ ;
            if (fptr < nFixed)
            {
                ft = fixed [fptr] ;
            }
        }
        else
        {
            xShort [kShort] = x [k] ;
            kShort++ ;
        }
    }
    return (1) ;
}