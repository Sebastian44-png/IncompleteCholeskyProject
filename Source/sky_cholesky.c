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

    n = A->n;
    p = A->p;
    d = A->d;
    x = A->x;
    /*calculate the number of nz entries for each row*/
    anzrow = malloc(n*sizeof(index));
    anzrow[0] = 0;
    for (index i = 1; i < n; i++){
        anzrow[i] = p[i] - p[i-1];
    }

    for (index k = 0; k < n; k++){
        /*auxiliary variable*/
        double sum = 0;

        /*calculate all entries in front of a_kk in the k-th row*/
        if (anzrow[k] > 0){
            for (index j = 0; j < k - 1; j++){
                    
            }
        }
        
        
        
        /*compute l_kk stored in d_k*/
        for (index i = 0; i < anzrow[k]-1; i++){
            sum += x[p[k-1] + i]*x[p[k-1] + i];
        }
        d[k] = sqrt(d[k] - sum);

    }
}
