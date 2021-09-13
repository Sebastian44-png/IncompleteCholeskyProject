#include "hpc.h"

int main(int argc, char **argv)
{
    index N ;

    cs *A_cs;
    sed *A_sed;
    sed *L_sed;


    index n = 9 ;
    double *x = malloc(n*sizeof(double));
    double *b = malloc(n*sizeof(double));
    b[0] = 1;
    b[1] = 1;
    b[2] = 1;
    b[3] = 1;
    b[4] = 1;
    b[5] = 1;
    b[6] = 1;
    b[7] = 1;
    b[8] = 1;


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
    
   /* sed_print(L_sed, 0);*/    
    
    sed_forwardInsertion(L_sed, x, b);
   
    printf("\n ==============\n Ergebnis:\n");
    for(index k = 0; k <n ; k++){
        printf("%f\n", x[k]);   
    }
    printf("\n ==============\n b:\n");                                  
    for(index k = 0; k <n ; k++){                                               
        printf("%f\n", b[k]);                                                   
    }








}
