#include "hpc.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

int main (int argc, char **argv)
{
    
    index N ;
    index n ;
    double *b ;
    double *x ;
    cs *A_cs ;
    sed *A_sed ;
    sed *L_sed_icf ;
    sed *L_sed_icne ;
    index maxIt ; 
    double alpha ;
    double tol ;
    
    /*variables for analysis*/
    index *anaAnzIt ;
    double *anaAvgTime ; 
    double *errorJacobi ;
    double *errorGauss ;
    double *errorICF ;
    double *errorICNE ;
    
    /*allocate memory*/
    N = 5 ;
    n = pow (N,2) ;
    maxIt = 10*n ;
    b = malloc(n*sizeof(double)) ;
    x = malloc(n*sizeof(double)) ;
    
    anaAnzIt = malloc(5*sizeof(index)) ;
    anaAvgTime = malloc(5*sizeof(double)) ;
    errorJacobi = malloc(maxIt*sizeof(double)) ;
    errorGauss = malloc(maxIt*sizeof(double)) ;
    errorICF = malloc(maxIt*sizeof(double)) ;
    errorICNE = malloc(maxIt*sizeof(double));

    /*create problem*/
    A_cs = cs_lapmat_p1_square (N) ;
    A_sed = sed_compress (A_cs); 
    L_sed_icf = sed_alloc(A_sed->n, A_sed->nzmax , 1);                                       
    L_sed_icne = sed_alloc(A_sed->n, A_sed->nzmax, 1);                                       
    for (index i = 0 ; i < n ; i++)
    {
        x [i] = 0 ;
        b [i] = 1 ;
    }
    tol = 1e-10 ;
    
    /*apply the incomplete cholesky factorization on the problem A*/
    sed_icholesky (A_sed , L_sed_icf) ;      
    
    /*apply the icne(0) factorization on the prolblem A*/
    //alpha = 1 ;
    //sed_icne0 (A_sed , alpha , L_sed_icne ) ;
    
    /*PCG with Jacobi as preconditioner*/
    TIME_SAVE (0) ;
    anaAnzIt [0] = sed_cg_jacobi (A_sed , b , x , maxIt , tol, errorJacobi) ;
    TIME_SAVE (1) ;
    anaAvgTime [0] = TIME_ELAPSED (0 , 1) / anaAnzIt [0] ; 
    

    /*PCG with Gauss-Seidel as preconditioner*/
    for (index i = 0 ; i < n ; i++)
    {
        x [i] = 0 ;
        b [i] = 1 ; 
    }
    TIME_SAVE (2) ;
    anaAnzIt [1] = sed_cg_gauss_seidel (A_sed , b , x , maxIt , tol, errorGauss) ;
    TIME_SAVE (3) ;
    anaAvgTime [1] = TIME_ELAPSED (2 , 3) / anaAnzIt [1] ; 

    /*PCG with ICF as preconditioner*/
    for (index i = 0 ; i < n ; i++)
    {
        x [i] = 0 ;
        b [i] = 1 ;
    }
    TIME_SAVE (4) ;
    anaAnzIt [2] = sed_ccg (A_sed , L_sed_icf , b , x , maxIt , tol, errorICF) ;
    TIME_SAVE (5) ;
    anaAvgTime [2] = TIME_ELAPSED (4 , 5) / anaAnzIt [2] ; 
    /*PCG with ICNE as preconditioner*/

    /*analysis in a txt file*/
    FILE *f ;
    f = fopen("Analyse_Problem_1.txt","w") ;
    fprintf(f,"error_Jacobi error_Gauss error_ICF error_ICNE error_Multigrid\n") ;
    for (index i = 0 ; i < maxIt ; i++)
    {
        fprintf(f,"%f %f %f %f %f\n" , errorJacobi [i] , errorGauss [i] , errorICF [i] , errorICNE [i] , errorICNE [i]) ;
    }
    fclose(f);

    /*free memory*/
    free (b) ;
    free (x) ;
    sed_free (A_sed) ;
    sed_free (L_sed_icf) ;
    sed_free (L_sed_icne) ;
    free (anaAnzIt) ;
    free (anaAvgTime) ;
    free (errorJacobi) ;
    free (errorGauss) ;
    free (errorICF) ;
    free (errorICNE) ;
}
