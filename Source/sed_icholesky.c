#include "hpc.h"
#include <math.h>

index sed_icholesky (sed *A, sed *L)
{

    if(!A||!L) return (0) ;
    
    index nzmax ;
    index n ;
    index *ind ;
    double *x ;
    index *anzCol ;
    index counter ;
    index alk ;
    index amk ;
    index alm ;

    nzmax = A->nzmax ;
    n = A->n ;
    ind = A->i ;
    x = A->x ;
    

    /*initialize the matrix L and get the number of elements in each column of A*/
    counter = 0 ;
    anzCol = malloc(n*sizeof(index)) ;
    /*first column starts at position n+1 in L->x*/
    L->i[0] = n+1 ;
    
    for (index j = 0 ; j < n ; j++)
    {
        anzCol[j] = ind[j+1] - ind[j] ;/*off diag nz of each col in A*/
        L->x[j] = x[j];
        for (index k = 0 ; k < anzCol[j] ; k ++)
        {
            if(ind[ind[j] + k] > j)
            {
                L->i[L->i[j] + counter] = ind[ind[j] + k] ;
                L->x[L->i[j] + counter] = x[ind[j] + k] ;
                counter += 1 ;
            }
        }
        L->i[j+1] = L->i[j] + counter ;
        counter = 0 ;
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

        /*apply on the rest of the matrix*/    
        for (index m = k + 1 ; m < n ; m++)
        {
            amk = sed_isoccupied(L, k , m) ;
            if (amk >= 0)
            {
                L->x[m] -=pow(L->x[L->i[k] + amk] , 2) ; 
                
                for ( index l = 0 ; l < L->i[m+1] - L->i[m] ; l++)
                {
                    alk = sed_isoccupied(L, k ,L->i[L->i[m] + l]) ;

                    if (alk >= 0)
                    {
                        L->x[L->i[m] + l] -= L->x[L->i[k] + alk] * L->x[L->i[k] + amk] ; 
                        alk = -1;
                    }
                }

            amk = -1;
            }
        }
    }
    return (1) ;

}
