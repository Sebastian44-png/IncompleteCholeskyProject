#include "hpc.h"
/* author: Benjamin Bestler */

/*computes incomplete LU factorization on the matrix A*/
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
        /* computations for diagonal elements i = j */
        for (i = k + 1 ; i < n ; i++)
        {
            if (Ax[i] == 0)
            {
                /* Diagonal entries must not be 0 */
                return (0) ;
            }

            /* find aik */
            aik = 0 ;
            for (index iptr = Ai[k]; iptr < Ai[k + 1]; iptr++)
            {
                if (Ai[iptr] == i)
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
            for (index iptr = Ai [i] ; iptr < Ai [i + 1] ; iptr++)
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

            /* aii = aii + fii = aii - aik * aki / akk */
            Ax [i] -= aik * akj / Ax [k] ;
        }

        /* other matrix entries, j is column index */
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

            /* Loop through j-th column */
            for (index ptr = Ai[j]; ptr < Ai[j + 1]; ptr++)
            {
                if (Ax [ptr] == 0 || Ai [ptr] <= k)
                {
                    /* Nothing to do if Entry is 0 or above or in k-th row */
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

                /* aij = aij + fij = aij - aik * akj / akk */
                Ax[ptr] -= aik * akj / Ax[k] ;
            }
        }
    }
    return (1) ;
}
