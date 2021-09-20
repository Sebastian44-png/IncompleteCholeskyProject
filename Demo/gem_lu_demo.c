#include "hpc.h"
/*
Author: Sebastian Acerbi
*/

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
    double *b_ref;
    double *b_sol;
    double *x; 

    cs *A_cs;
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
    A_gem = gem_compress(A_cs);

    printf("A = \n");
    gem_print(A_gem, 0);
    
    printf("LU-Decomposition of A = \n");
    gem_lu(A_gem);
    gem_print(A_gem, 0);

    free(b);
    free(b_ref);
    free(b_sol);
    free(x);

    gem_free(A_gem);
    cs_free(A_cs);
}