#include "hpc.h"
/* author: Benjamin Bestler */

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
        /* iterate through k-th row */
        for (j = k + 1 ; j < n ; j++)
        {
            /* Search for k-th entry in j-th column and do computations if found */ 
            for (index akjptr = Ai [j] ; akjptr < Ai [j + 1]; akjptr++)
            {
                if (Ai [akjptr] == k)
                {
                    /* going through aik */
                    for (index aikptr = Ai [k] ; aikptr < Ai [k + j] ; aikptr++)
                    {
                        i = Ai [aikptr] ; 

                        if (i <= k)
                        {
                            continue ;
                        }

                        if (i == j)
                        {
                            /* Diagonal entry */
                            Ax [i] -= Ax [akjptr] * Ax [aikptr] / Ax [k] ;
                        }
                        else
                        {
                            /* Non-diagonal entry */

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

                            /* If not found, modify aii */
                            if (!found)
                            {
                                Ax [i] -= alpha * Ax [akjptr] * Ax [aikptr] / Ax [k] ;
                            }
                        }
                    }
                    break ;
                }
            }
        }
    }
    return (1) ;
}