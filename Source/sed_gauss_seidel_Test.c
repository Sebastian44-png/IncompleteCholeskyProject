#include "hpc.h"                                                                
                                                                                
/* Do the Fauss -Seidel Algorithm on a given Matrix A and a right side b
until the the error is small enough or the maxIt is reached*/
index sed_gauss_seidel(const sed *A, const double *b, double *xk, double *w, index maxIt, double tol)    
{
    index n, *ind ;                                                              
    double *Ax ;                                                                 
    index it = 0 ;
    double error =     
                                                                                
    n = A->n; ind = A->i; Ax = A->x;                                            
   while (it < maxIt && )                                                                            
    //copy xk to w                                                              
    for(index i = 0; i<n; i++){                                                 
        w[i] = xk[i];                                                           
    }                                                                           
                                                                               
    //iterate through rows and cols of matrix                                   
    for(index k = 0; k < n; k++){                                               
       xk[k] = b[k];                                                           
                                                                               
       for(index j = 0; j < n; j++){                                           
            if(k == j) continue;                                                
                                                                               
            // find row k in col j                                              
            for(index i = ind[j]; i < ind[j+1]; i++){                           
                if(ind[i] == k){                                                
                                                                                
                    xk[k] -= xk[j] * Ax[i];                                     
                }                                                               
            }                                                                   
        }                                                                       
        xk[k] = xk[k] / Ax[k];                                                  
    }                                                                           
    return 1;
}
