#include "hpc.h"
/* Computes relaxed modified incomplete LU-factorization with relaxation factor alpha
A is required to have only nonzero diagonal entries */

index sed_MILU (sed *A, const double alpha)
{
    if (alpha == 0)
    {
        return sed_ILU (A) ;
    }
    index n ;
    index nzmax ;
    index *Ai ;
    double *Ax ;

    index i ;
    index j ;

    double aik ;
    double akj ;

    bool found ;

    n = A->n ;
    nzmax = A->nzmax ;
    Ai = A->i ;
    Ax = A->x ;

    for (index k = 0 ; k < n - 1 ; k++)
    {
        /* computations for diagonal elements , no difference to ILU as akk != 0 */
        for (index ptr = k + 1 ; ptr < n ; ptr++)
        {
            if (Ax[ptr] == 0)
            {
                return (0) ;
            }

            /* find aik */
            aik = 0 ;
            for (index iptr = Ai[k] ; iptr < Ai[k + 1] ; iptr++)
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
            for (index iptr = Ai [ptr] ; iptr < Ai[ptr + 1] ; iptr++)
            {
                if (Ai [iptr] == k)
                {
                    akj = Ax[iptr] ;
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
        for (j = k + 1 ; j < n ; j++)
        {
            for (index akjptr = Ai [j] ; akjptr < Ai [j + 1]; akjptr++)
            {
                /* going through akj */
                if (Ai [akjptr] != k)
                {
                    continue ;
                }

                for (index aikptr = Ai [k] ; aikptr < Ai [k + j] ; aikptr++)
                {
                    /* going through aik */
                    i = Ai [aikptr] ; 

                    if (i <= k)
                    {
                        continue ;
                    }

                    /* find aij and modify if found */
                    found = false ;
                    for (index iptr = Ai[j] ; iptr < Ai[j + 1] ; iptr++)
                    {
                        if (Ai[iptr] == i)
                        {
                            found = true ;
                            Ax [iptr] -= Ax [akjptr] * Ax [aikptr] / Ax [k] ;
                            break ;
                        }
                    }

                    if (!found)
                    {
                        Ax [i] -= alpha * Ax [akjptr] * Ax [aikptr] / Ax [k] ;
                    }
                }
            }
        }
    }
    return (1) ;
}