#include "hpc.h"

/* writes values from x into non-fixed entries from xLong */
index hpc_expandFixed (double *x, double *xLong, index n, const index *fixed, const index nFixed)
{
    index k ;
    index ft ;
    index fptr ;

    if (!x || !xLong || !fixed)
    {
        return (0) ;
    }

    // k is position in original vector
    k = 0 ;

    // to manage fixed entries, will be skipped for insertion
    fptr = 0 ;
    ft = fixed [fptr] ;

    // kLong is position in elongated vector
    for (index kLong = 0 ; kLong < n ; kLong++)
    {
        if (kLong == ft)
        {
            fptr++ ;
            if (fptr < nFixed)
            {
                ft = fixed [fptr] ;
            }
        }
        else
        {
            xLong [kLong] = x [k] ;
            k++ ;
        }
    }
    return (1) ;
}