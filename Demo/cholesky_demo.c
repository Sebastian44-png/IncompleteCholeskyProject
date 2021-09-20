
#include "hpc.h"

int main(int argc, char **argv)
{
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
}