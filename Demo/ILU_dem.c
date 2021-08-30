
#include "hpc.h"

int main(int argc, char **argv)
{
    index N;

    cs *A_cs;
    sed *A_sed;

    N = 5;
    A_cs = cs_lapmat_p1_square(N);

    printf("\nDimension of matrix    = ( %g, %g)\n", (double)A_cs->n, (double)A_cs->m);

    A_sed = sed_compress(A_cs);

    sed_ILU(A_sed);

    sed_print(A_sed, 0);
}