#include "hpc.h"
bnd *bnd_compress (const cs *T)
{
    index j, k, m, n, nz, *Ti, *Tj, w;
    double *Bx, *Tx;
    bnd *B ;
    if (!HPC_TRIPLET (T)) return (NULL) ;         /* check inputs */
    m = T->m ; n = T->n; Tx = T->x ; Ti = T->ind ; Tj = T->p ; nz = T->nz ;
    if (m != n) return (NULL) ;                   /* check square matrix */
    w = 0;                                        /* compute bandwidth */
    for ( j = 0 ; j < nz ; j++) w = HPC_MAX(w,labs(Ti[j]-Tj[j]));
    /*  printf("bandwidth = %g\n",(double) w); */
    B = bnd_alloc(n, w) ;                         /* allocate result */
    if (!B) return (NULL) ;                       /* out of memory */    
    Bx = B->x ; 
    for (k = 0 ; k < nz ; k++) Bx[(2*Ti[k]+1)*w+Tj[k]] += Tx[k];
    return (B) ;
}
