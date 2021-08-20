#include "hpc.h"
/* allocate a (n x n) band matrix, with band width m */
bnd *bnd_alloc (index n, index m)
{
    bnd *A = calloc (1, sizeof (bnd)) ;    /* allocate the sky struct */
    if (!A) return (NULL) ;                /* out of memory */
    A->n = n ;                             /* define dimension */
    A->m = m ;                             /* define bandwidth */
    A->x = calloc (n * (2*m+1), sizeof (double)) ; /* allocate matrix */
    return ( !(A->x) ? bnd_free (A) : A) ;
}

/* free a gem matrix */
bnd *bnd_free (bnd *A)
{
    if (!A) return (NULL) ;      /* do nothing if A already NULL */
    free (A->x) ;                /* free the gem struct and return NULL */
    free (A);
    return (NULL) ; 
}

