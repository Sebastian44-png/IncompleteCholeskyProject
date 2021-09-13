
#include "hpc.h"

int main(int argc, char **argv)
{
    /*
    index N;

    cs *A_cs;
    sky *A_sky;

    N = 3;
    A_cs = cs_lapmat_p1_square(N);
    gem_print(gem_compress(A_cs), 0);

    printf("\nDimension of matrix    = ( %g, %g)\n", (double)A_cs->n, (double)A_cs->m);

    A_sky = sky_compress(A_cs);

    sky_cholesky(A_sky);

    sky_print(A_sky, 0);
    */
    
    index n = 5 ;
    double *x = malloc(n*sizeof(double));
    double *b = malloc(n*sizeof(double));
 
    b[0] = 1;
    b[1] = 0;
    b[2] = 1;
    b[3] = 0;
    b[4] = 1;
 
    sed *L;
    L = calloc(5, 15);
    for (index ptr = 0 ; ptr < n ; ptr++){
        L->x[ptr] = 1;
        L->i[ptr] = n+1+ptr;
    }
    for (index j = L->i[0] ; j < L->i[n] ; j ++){
        L->x[j] = 1;
    }
    sed_print(L,1);
}
