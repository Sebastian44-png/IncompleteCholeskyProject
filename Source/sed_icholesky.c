#include "hpc.h"
#include <math.h>

/*author: Joachim Kröner*/

/*This method computes the incomplete cholesky factorization*/
index sed_icholesky (sed *A, sed *L)
{
    /*check input*/
    if(!A||!L) return (0) ;
    
    index nzmax ;
    index n ;
    index *ind ;
    double *x ;

    /*auxiliary variables*/
    index *anzCol ;
    index counter ;
    index alk ;
    index amk ;
    index alm ;

    nzmax = A->nzmax ;
    n = A->n ;
    ind = A->i ;
    x = A->x ;

    counter = 0 ;
    anzCol = malloc(n*sizeof(index)) ;

   /*copy the matrix A in L*/ 
    for (index j = 0 ; j < nzmax ; j++)
    {
        L->i [j] = ind [j] ;
        L->x [j] = x [j] ;
    }
    for (index j = 0 ; j < n ; j++)
    {   
        /*get the number of non zero off diagonal elements of each row*/
        anzCol[j] = ind[j+1] - ind[j] ;
    }
    
    /*apply the incomplete cholesky decomposition*/
    for (index k = 0 ; k < n ; k++)
    {
        /*calculate the diagonal element*/
        L->x[k] = sqrt(L->x[k]) ;
         
        for (index l = 0 ; l < L->i[k+1] - L->i[k] ; l++)
        {
            /*calculate the k-th column*/
            L->x[L->i[k] + l] = L->x[L->i[k] + l] / L->x[k] ;
        }

        /*apply changes on the rest of the matrix*/    
        for (index m = k + 1 ; m < n ; m++)
        {   
            /*check if the element l_mk == 0*/
            amk = sed_isoccupied (L, k , m) ;
            if (amk >= 0)
            {
                L->x [m] -=pow (L->x [L->i [k] + amk] , 2) ; 
                
                for ( index l = 0 ; l < L->i [m+1] - L->i [m] ; l++)
                {
                    /*check if the element l_lm*/
                    alk = sed_isoccupied(L, k ,L->i [L->i [m] + l]) ;

                    if (alk >= 0)
                    {
                        L->x [L->i [m] + l] -= L->x [L->i [k] + alk] * L->x [L->i [k] + amk] ; 
                        /*reset column index*/
                        alk = -1 ;
                    }
                }
            /*reset column index*/
            amk = -1 ;
            }
        }
    }

    /*memory releas*/
    free (anzCol) ;

    return (1) ;
}
