#include "hpc.h"
/* perform the lu-decompostion on an matrix stored
   in the sed format, the original matrix gets destroyed */

index sed_lu(sed *A, const double *b, const double *x){

    index n;
    index *ind;   
    double *Ax;
    index akj_ind;

    n = A->n;
    ind = A->i;
    Ax = A->x;

    for(index k = 0; k < n; k++){
        
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
            
            /* if a_kj-th element zero, go to next column */
            if(!akj){ continue; }

            /* compute a_ij = a_ij a_ik - a-kj */
            for(index i = akj_ind + 1; i < ind[k+1]; i++){
                
                Ax[i] = Ax[i] - Ax[k] * akj_ind;
            }
        }
    }
}
