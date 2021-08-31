#include "hpc.h"

#include<time.h>
#include <sys/time.h>
#include <stdlib.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

int main (int argc, char **argv)
{
    index N; /* N² is the Dimension of the MAtrix */ 

    double *b;
    double *b_ref;
    double *b_sol;
    double *x; 

    cs *A_cs;
    sed *A; /* gets overwritten with the lu decompostition of A*/
    sed *A_ref;
    gem *A_gem;
    
    /* Get number N as parameter */ 
    N = 0;
    if ( argc > 1 ) {
        if ( (atoi(argv[1]) > 0) & (atoi(argv[1]) < 1000) ) { 
            N = atoi(argv[1]);
        }
    } 
   
    A_cs = cs_lapmat_p1_square(N);  
    A = sed_compress(A_cs);
    A_ref = sed_compress(A_cs);
    A_gem = gem_compress(A_cs);

    gem_print(A_gem, 1);

    if (!A) return(1);

    printf ("\nDimension of matrix    = ( %g, %g)\n", (double) A->n, (double) A->n); 
    printf("\n========================================\n");
    printf("A before lu decomposition \n");
    sed_print(A, 0);

    if (!A->n) return(1);        /* check input */

    /* ------------------------- */
    /* Allocate and initiate x and b */  
    /* ------------------------- */

    index n = A->n;
    b = malloc(n * sizeof(double)); 
    b_ref = malloc(n * sizeof(double));
    b_sol = malloc(n * sizeof(double));
    x = malloc(n * sizeof(double));

    for (index j = 0 ; j < n ; j++ ) {
        b[j] = 1;
        b_ref[j] = 1;
        x[j] = 0;
        b_sol[j] = 0;
    }
    /* ---------------------- */
    /* apply decomposition */ 
    /* ---------------------- */
    TIME_SAVE(0);
    if (!sed_lu(A)) return(1);       /* perform Gauss decomposition  */
    TIME_SAVE(1);

    printf("\n========================================\n");
    printf("A after lu decomposition \n");
    sed_print(A, 0);

    //if (!gem_gausssol(A, b)) return(1); /* compute A^(-1) *  b */
    //TIME_SAVE(2);

    /* --------------------- */    
    /* verify results */
    /* --------------------- */
    //if(!gem_spmv(A_ref, b, b_sol));
    //TIME_SAVE(3);

    // double error_gauss = 0.0;
    //for ( j = 0 ; j < n ; j++ ) { 
      //  error_gauss = HPC_MAX(error_gauss, abs(b_sol[j] - b_ref[j]));
    //}

    printf("Time lu-decomp   = %9i ns\n", (int) TIME_ELAPSED(0,1) );
    //printf("Time Gauss solve    = %9i ns\n", (int) TIME_ELAPSED(1,2) );
    //printf("result - max( sol ) = %20.16f \n", error_gauss); 

    //printf("time gem spmv  =  = %9i ns\n", (int) TIME_ELAPSED(2,3));

    /* --------------------- */
    /* clear memory */  
    /* +++++++++++++++++++++ */
    cs_free (A_cs); 
    sed_free(A); 
    sed_free(A_ref);

    free(b);
    free(b_ref);
    free(b_sol);
    free(x);

    return (0) ;
}
