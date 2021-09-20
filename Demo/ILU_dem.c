#include "hpc.h"
/* author: Benjamin Bestler */

int main (int argc, char **argv)
{
    index N ;
    N = 5 ;
    double alpha ;
    alpha = 0.5 ;
    
    /* get grid size */
    if (argc > 1)
    {                        
        if (atoi(argv[1]) > 0) 
        {
            N = atoi(argv[1]) ;
        }
    } 

    /* get alpha */
    if (argc > 2)
    {                        
        if (atof(argv[2]) > 0) 
        {
            alpha = atof(argv[2]) ;
        }
    } 

    cs *A_cs ;
    sed *A_sed_ILU ;
    sed *A_sed_MILU ;
    sed *A_sed_MILU_rel ;

    /* Create problem */
    A_cs = cs_lapmat_p1_square(N) ;

    /* Prepare matrix and apply ILU */
    A_sed_ILU = sed_compress(A_cs) ;
    sed_ILU(A_sed_ILU) ;

    printf("ILU factorization of Matrix A:\n") ;
    sed_print(A_sed_ILU, 0) ;

    /* Prepare matrix and apply MILU */
    A_sed_MILU = sed_compress(A_cs) ;
    sed_MILU(A_sed_MILU, 1) ;

    printf("MILU factorization of Matrix A:\n") ;
    sed_print(A_sed_MILU, 0) ;

    /* Prepare matrix and apply MILU */
    A_sed_MILU_rel = sed_compress(A_cs) ;
    sed_MILU(A_sed_MILU_rel, alpha) ;

    printf("relaxed MILU factorization of Matrix A with alpha = %.2f:\n", alpha) ;
    sed_print(A_sed_MILU_rel, 0) ;
}