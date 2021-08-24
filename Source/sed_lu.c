#include "hpc.h"

index sed_lu(sed *A, const double *b, const double *x){

    index n;
    index *ind;   
    double *Ax;

    n = A->n;
    ind = A->i;
    Ax = A->x;

    for(index k = 0; k < n; k++){
        for(index i = k + 1; i< n; i++){
                break;
            }
        }

}
