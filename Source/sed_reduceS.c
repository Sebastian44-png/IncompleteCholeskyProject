#include "hpc.h"

/* takes top level Stiffness Matrix A and reduces it to only non-fixed row/columns */
sed *sed_reduceS (const sed *S, const index *fixed, const index nFixed)
{
    sed *SRed ;
    index nRed ;
    index fptr ;
    index kR ;
    index Rptr ;
    index ft ;
    index *isf ;
    index *NewRow ;

    index Sn ;
    index Snzmax ;
    index *Si ;
    double *Sx ;

    index SRn ;
    index *SRi ;
    double *SRx ;

    isf = malloc (Sn * sizeof(index)) ;
    NewRow = malloc (Sn * sizeof(index)) ;

    // check inputs and if mallocs were sucessfull
    if (!S || !isf || !NewRow)
    {
        free (isf) ;
        free (NewRow) ;
        return (NULL) ;
    }

    Sn = S->n ;
    Snzmax = S->nzmax ;
    Si = S->i ;
    Sx = S->x ;

    // Store in isf at each index if it is fixed
    for (index k = 0 ; k < Sn ; k++)
    {
        isf [k] = 0 ;
    }

    for (fptr = 0 ; fptr < nFixed ; fptr++)
    {
        isf [fixed [fptr]] = 1 ;
    }

    // Set NewRow to contain row index of corresponding row in new matrix, nRed counts fixed rows
    nRed = 0;
    for (index k = 0 ; k < Sn ; k++)
    {
        if (isf [k])
        {
            nRed++ ;
        }
        else
        {
            NewRow [k] = k - nRed ;
        }
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

    // Allocate new sed Matrix and check if sucessfull
    SRed = sed_alloc (Sn - nFixed, Snzmax - nRed, 1) ;

    if (!SRed)
    {
        free (isf) ;
        free (NewRow) ;
        sed_free(SRed) ;
        return (NULL) ;
    }

    SRn = SRed->n ;
    SRi = SRed->i ;
    SRx = SRed->x ;

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
                    SRi [Rptr] = NewRow [Si [ptr]] ;
                    Rptr++ ;
                }
            }
            kR++ ;
            // pointer advancment indicates start of next column
            SRi [kR] = Rptr ;
        }
    }

    // success, free buffers and return reduced matrix
    free (isf) ;
    free (NewRow) ;
    return (SRed) ;
}