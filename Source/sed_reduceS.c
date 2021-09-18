#include "hpc.h"

/* takes top level Stiffness Matrix A and reduces it to only non-fixed row/columns */
sed *sed_reduceS (sed *S, index *fixed, index nFixed)
{
    sed *SRed ;
    index nRed ;
    index fptr ;
    index kR ;
    index Rptr ;
    index ft ;
    index *isf ;

    index Sn ;
    index Snzmax ;
    index *Si ;
    double *Sx ;

    index SRn ;
    index *SRi ;
    double *SRx ;

    Sn = S->n ;
    Snzmax = S->nzmax ;
    Si = S->i ;
    Sx = S->x ;

    isf = malloc (Sn * sizeof(index));
    if (!S || !isf)
    {
        return (NULL) ;
    }

    for (index k = 0 ; k < Sn ; k++)
    {
        isf [k] = 0 ;
    }

    for (fptr = 0 ; fptr < nFixed ; fptr++)
    {
        isf [fixed [fptr]] = 1 ;
    }

    // set nRed to number of deleted diagonal entries
    nRed = nFixed ;
    fptr = 0 ;
    ft = fixed [0] ;
    for (index k = 0 ; k < Sn ; k++)
    {
        if (k == ft)
        {
            // Increas nRed by number of off-Diagonal entries of fixed columns
            nRed += Si [ft + 1] - Si [ft] ;
            fptr++ ;
            if (fptr < nFixed)
            {
                ft = fixed [fptr] ;
            }
        }
        else
        {
            for (index ptr = Si [k] ; ptr < Si [k + 1] ; k++)
            {
                if (isf [Si [ptr]])
                {
                    // if row index is fixed, increase nRed.
                    nRed++;
                }
            }
        }
    }

    SRed = sed_alloc (Sn - nFixed, Snzmax - nRed, 1) ;
    SRn = SRed->n ;
    SRi = SRed->i ;
    SRx = SRed->x ;
    if (!SRed)
    {
        return (NULL) ;
    }

    // Fill reduced matrix with entries
    kR = 0;
    Rptr = SRn + 1 ;
    SRi [0] = Rptr ;
    for(index k = 0 ; k < Sn ; k++) 
    {
        if (!isf [k])
        {
            // diagonal entry
            SRx [kR] = Sx [k] ;
            // fill column
            for (index ptr = Si [k] ; ptr < Si [k + 1] ; ptr++)
            {
                if (!isf [Si [ptr]])
                {
                    SRx [Rptr] = Sx [ptr] ;
                    SRi [Rptr] = Si [ptr] ;
                    Rptr++ ;
                }
            }
            kR++ ;
            // pointer advancment indicates start of next column
            SRi [kR] = Rptr ;
        }
    }
    return (SRed) ;
}