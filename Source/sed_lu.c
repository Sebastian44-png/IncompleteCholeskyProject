#include "hpc.h"
/* perform the lu-decompostion on an matrix stored
   in the sed format, the original matrix gets destroyed */

index sed_lu (sed *A){

    index n;
    index nnz;
    index *ind;  
    double *Ax;
    double f;
    index akj_ind;

    n = A->n;
    ind = A->i;
    nnz = A->nnz;
    Ax = A->x;

    if (!A) return (0) ;   /* check inputs */

    for(index k = 0; k < n; k++){
        
        for(index col = 0; col < n; col++){
            for(index ind[])
            if (Ax[]){


                i = 
                j = col


                f = -Ax[] / Ax[k];
            }
        }


        /* compute a_ik = a_ik / a_kk */
        for(index i = ind[k]; i < ind[k+1]; i++){
                Ax[i] = Ax[i] / Ax[k];
        }

        for(index j = k+1; j < n; j++){

            /* find a_kj-th element*/
            double akj = 0;

            for(index i = ind[k]; i < ind[k+1]; i++){

                if(ind[i] == k){ 
                    akj = Ax[i];
                    akj_ind = i;
                }
            }
            
            /* if a_kj-th element is zero, go to next column */
            if(!akj){ continue; }
            
            /* compute a_ij = a_ij - a_ik/a_kj */
            Ax[k] = Ax[k] - akj;  
            for(index i = akj_ind + 1; i < ind[k+1]; i++){
                
                Ax[i] = Ax[i] - Ax[k] * akj_ind;
            }
        }
    }
    return(1);
}
