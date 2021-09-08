
#include "hpc.h"

int main(int argc, char **argv)
{
    index N;

    cs *A_cs;
    sed *A_sed;

    N = 5;
    A_cs = cs_lapmat_p1_square(N);

    A_sed = sed_compress(A_cs);

    sed_ILU(A_sed);

    printf("After ILU") ;
    sed_print(A_sed, 0);

    A_sed = sed_compress(A_cs);

    sed_MILU(A_sed, 1) ;

    printf("Solution of MILU:") ;
    sed_print(A_sed, 0);
}