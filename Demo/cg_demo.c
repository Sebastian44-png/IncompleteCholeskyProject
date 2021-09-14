#include "hpc.h"

#include<time.h>
#include <sys/time.h>
#include <stdlib.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

int main (int argc, char **argv)
{
    index j;
    index n;
    index N;

    double *b;
    double *x0;

    b = malloc(n*sizeof(double));
    x0 = malloc(n*sizeof(double));
    
    cs *A_cs;
    sed *A_sed;
    gem *A_gem;

    printf("\n========================================\n");
    /* Get number N as parameter */ 
    N = 0;
    if ( argc > 1 ) {
        if ( (atoi(argv[1]) > 0) & (atoi(argv[1]) < 1000) ) { 
            N = atoi(argv[1]);
        }
    } 
   
    A_cs = cs_lapmat_p1_square(N);  
    A_sed = sed_compress(A_cs);
    A_gem = gem_compress(A_cs);

    n = A_cs->n;
    printf("gem_Matrix \n");

    index maxIt = 3;
    double tol = 0.0000001;

    for(index k = 0; k < n; k++){
        x0[k] = k;
        b[k] = 2;
    }

    sed_cg(A_sed , b, x0, maxIt, tol);

    
    gem_print(A_gem, 0);
}
    