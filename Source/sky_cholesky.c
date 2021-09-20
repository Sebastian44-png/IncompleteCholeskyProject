#include "hpc.h"
#include <math.h>

/*author: Joachim Kröner*/

/*This method calculates the Cholesky factorization of a symetric positive definite matrix.*/
index sky_cholesky(sky *A)
{
    /*check input*/
    if(!A)
    {
        return (0) ;
    }
    
    index n ;
    index *p ;
    double *d ;
    double *x ;
    double *anzrow ;
    double sum ;
    index ptrDiff ;
    index i ;
    index startItI ;

    n = A->n ;
    p = A->p ;
    d = A->d ;
    x = A->x ;
    
    /*calculate the number of non zero entries in each row*/
    anzrow = malloc (n * sizeof(index)) ;
    anzrow [0] = 0 ;
    for (index ptr = 1 ; ptr < n ; ptr++)
    {
        anzrow [ptr] = p [ptr] - p [ptr-1] ;
    }
    
    /* iteration over the rows of A*/
    for (index k = 0 ; k < n ; k++)
    {
       
        for(index ptr = 0 ; ptr < anzrow [k] ; ptr++)
        {
            
            /*get the offset of the columnindex for the first non zero element in row i
              compared to row k*/
            i = k - anzrow [k] + ptr ;
            ptrDiff = (i - anzrow [i]) - (k - anzrow [k]) ;

            if(ptrDiff >= 0)
            {
                startItI = 0 ;
            }
            else
            {
                startItI = -ptrDiff ;
            }

            /*compute the inner product of row i and k*/
            sum = 0 ;
            for(index iPtr = startItI; iPtr < anzrow[i]; iPtr++)
            {
                sum += x [p [k - 1] + iPtr + ptrDiff] * x [p [i - 1] + iPtr] ;
            }
            x [p [k - 1] + ptr] = (x [p [k - 1] + ptr] - sum) / d [i] ; 
        } 
        
        /*compute l_kk stored in d_k*/
        sum = 0;
        for (index j = 0 ; j < anzrow [k] ; j++)
        {
            sum += x [p [k - 1] + j] * x [p [k - 1] + j] ;
        }
        d [k] = sqrt (d [k] - sum) ;

    }
    /*memory release*/
    free (anzrow) ;

    return (1) ;
}
