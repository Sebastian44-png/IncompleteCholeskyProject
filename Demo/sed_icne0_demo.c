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
    double *b_ref;
    double *b_sol;
    double *x; 

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

    printf("gem_Matrix \n");
    gem_print(A_gem, 0);
    
    //printf("sed_Matrix \n");
    //sed_print(A_sed, 0); 
    //printf("Buffers :");
    
    sed *L = sed_alloc(A_cs->n, 0, 1);
    sed_icne0(A_sed, 2.0, L);

    printf("Matrix after icne(0) decomposition\n");
    sed_print(L, 0);
    //printf("x: "); print_buffer_double(A_sed->x, A_sed->nzmax);
    //printf("ind: "); print_buffer_int(A_sed->i, A_sed->nzmax);

    //printf("x: "); print_buffer_double(L->x, L->nzmax);
    //printf("ind: "); print_buffer_int(L->i, L->nzmax);
    //printf("%d \n",(int)L->n);
    
    sed * LLt = sed_alloc(L->n, 0, 1);
    sed_L_to_LLt(L, LLt);
    sed_print(LLt, 0);

    //printf("x: "); print_buffer_double(L->x, L->nzmax+1);
    //printf("ind: "); print_buffer_int(L->i, L->nzmax+1);

    cs_free (A_cs); 
    sed_free(A_sed);
    gem_free(A_gem);
    sed_free(L);
    //free(L->i);
    //free(L->x);
    //sed_free(LLt);
    free(b);
    return (0) ;
}
