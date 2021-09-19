#include "hpc.h"
#include <math.h>
/*This method computes the Cholesky decomp of a symetrix positive definite Matrix.
 * The original Matrix will be destroied*/

index sky_cholesky(sky *A){
    
    index n;
    index *p;
    double *d;
    double *x;
    double *anzrow;
    double sum;
    index ptrDiff;
    index i;

    n = A->n;
    p = A->p;
    d = A->d;
    x = A->x;
    
    /*check input*/
    if(!A)
    {
        return (0) ;
    }

    /*calculate the number of nz entries for each row*/
    anzrow = malloc(n*sizeof(index));
    anzrow[0] = 0;
    for (index ptr = 1; ptr < n; ptr++){
        anzrow[ptr] = p[ptr] - p[ptr-1];
    }

    for (index k = 0; k < n; k++){
       
        for(index ptr = 0; ptr < anzrow[k]; ptr++){
            i = k - anzrow[k] + ptr;
            ptrDiff = (i - anzrow[i]) - (k - anzrow[k]);

            index startItI;
            if(ptrDiff >= 0){
                startItI = 0;
            }
            else{
                startItI = -ptrDiff;
            }
            
            sum = 0;
            for(index iPtr = startItI; iPtr < anzrow[i]; iPtr++){
                sum += x[p[k - 1] + iPtr + ptrDiff] * x[p[i - 1] + iPtr];
            }
            x[p[k - 1] + ptr] = (x[p[k - 1] + ptr] - sum) / d[i]; 
        } 
        
        /*compute l_kk stored in d_k*/
        sum = 0;
        for (index j = 0; j < anzrow[k]; j++){
            sum += x[p[k - 1] + j] * x[p[k - 1] + j];
        }
        d[k] = sqrt(d[k] - sum);

    }
    return (1) ;
}
