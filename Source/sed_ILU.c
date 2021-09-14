#include "hpc.h"

index sed_ILU (sed *A)
{
    index n ;
    index nzmax ;
    index *Ai ;
    double *Ax ;

    index i ;
    index j ;

    double aik ;
    double akj ;

    n = A->n ;
    nzmax = A->nzmax ;
    Ai = A->i ;
    Ax = A->x ;

    for (index k = 0 ; k < n - 1 ; k++)
    {
        /* computations for diagonal elements */
        for (index ptr = k + 1 ; ptr < n ; ptr++)
        {
            if (Ax[ptr] == 0)
            {
                return (0) ;
            }

            /* find aik */
            aik = 0 ;
            for (index iptr = Ai[k]; iptr < Ai[k + 1]; iptr++)
            {
                if (Ai[iptr] == ptr)
                {
                    aik = Ax[iptr] ;
                    break ;
                }
            }
            if (aik == 0)
            {
                continue ;
            }

            /* find akj */
            akj = 0;
            for (index iptr = Ai [ptr] ; iptr < Ai [ptr + 1] ; iptr++)
            {
                if (Ai [iptr] == k)
                {
                    akj = Ax [iptr] ;
                    break ;
                }
            }
            if (akj == 0)
            {
                continue ;
            }

            Ax [ptr] -= aik * akj / Ax [k] ;
        }

        /* other matrix entries */
        for (j = k + 1; j < n; j++)
        {
            /* find akj */
            akj = 0;
            for (index iptr = Ai[j]; iptr < Ai[j + 1]; iptr++)
            {
                if (Ai[iptr] == k)
                {
                    akj = Ax[iptr] ;
                    break ;
                }
            }
            if (akj == 0)
            {
                continue ;
            }

            for (index ptr = Ai[j]; ptr < Ai[j + 1]; ptr++)
            {
                if (Ax [ptr] == 0 || Ai [ptr] <= k)
                {
                    continue ;
                }

                i = Ai[ptr] ;

                /* find aik */
                aik = 0;
                for (index iptr = Ai [k] ; iptr < Ai [k + 1] ; iptr++)
                {
                    if (Ai [iptr] == i)
                    {
                        aik = Ax [iptr] ;
                        break ;
                    }
                }
                if (aik == 0)
                {
                    continue ;
                }

                Ax[ptr] -= aik * akj / Ax[k] ;
            }
        }
    }
    return (1) ;
}