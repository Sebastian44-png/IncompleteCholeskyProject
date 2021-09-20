#include "hpc.h"
/* author: Benjamin Bestler */

/* writes non-fixed entries from x into xShort */
index hpc_reduceFixed (double *x, double *xShort, index n, const index *fixed, const index nFixed)
{
    index kShort ;
    index ft ;
    index fptr ;

    // check inputs
    if (!x || !xShort || !fixed)
    {
        return (0) ;
    }

    // kShort is position in shortened vector
    kShort = 0 ;

    // to managed fixed/skipped entries
    fptr = 0 ;
    ft = fixed [fptr] ;

    // k is position in original vector
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