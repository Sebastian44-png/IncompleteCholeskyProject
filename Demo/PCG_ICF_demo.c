#include "hpc.h"
#include <math.h>

int main (int argc, char **argv)
{
    index N;
    
    cs *A_cs ;
    sed *A_sed ;
    double *b ;
    double *b_tilde ;
    
    sed *L_sed ;
    index lnzmax;
    

    N = 3 ;
    index n = pow(N,2) ;
    b = malloc(n*sizeof(double)) ;
    b_tilde = malloc (n *sizeof(double)) ;
    A_cs = cs_lapmat_p1_square(N) ;
    
    for (index i = 0 ; i < n ; i++)
    {
        b[i] = 1 ;
    }

    printf("\nProblem:\n") ;
    gem_print(gem_compress(A_cs), 0);
    
    printf("\n=============================\n") ;
    printf("Incomplete Cholesky Faktorization\n") ;
    
    A_sed = sed_compress(A_cs); 
        
    index lmax = ((A_sed->nzmax - A_sed->n)/2) + A_sed->n + 1;
    L_sed = sed_alloc(A_sed->n, lmax, 1);                                       
    sed_icholesky(A_sed, L_sed);      
    
    sed_print(L_sed,0);


    printf("\n===========================\n") ;
    printf("forward insertion L*b_tilde = b\n") ;
    
    sed_forwardInsertion(L_sed, b_tilde, b );

    for (index i = 0 ; i < n ; i++)
    {
        printf("%f\n",b_tilde[i] ) ;
    }

}
