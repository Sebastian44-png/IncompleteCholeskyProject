#include "hpc.h"

int main(int argc, char **argv)
{
    index N ;

    cs *A_cs;
    sed *A_sed;
    sed *L_sed;

    N = 3;
    A_cs = cs_lapmat_p1_square(N);
    printf("Original Matrix\n");
    gem_print(gem_compress(A_cs), 0);

    printf("==========================\n");
    printf("Incomplete Cholesky\n");


    A_sed = sed_compress(A_cs);

    index lmax = ((A_sed->nzmax - A_sed->n)/2) + A_sed->n + 1;

    L_sed = sed_alloc(A_sed->n, lmax, 1);
    
    sed_icholesky(A_sed, L_sed);
    
    sed_print(L_sed, 0);    
}
