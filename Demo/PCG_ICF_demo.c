#include "hpc.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

int main (int argc, char **argv)
{
    /*create arrays for analysis*/
    index *anaAnzIt = malloc(5*sizeof(index)) ;
    double *anaAvgTime = malloc(5 * sizeof(double)) ; 
    
    /*create the Problem*/
    /*get local varialbes for solving the problem*/
    index N ;
    index n ;
    double *b ;
    double *x ;
    cs *A_cs ;
    sed *A_sed ;
    sed *L_sed ;
    index lnzmax;
    
    /*allocate memory*/   
    N = 3;
    n = pow(N,2) ;
    b = malloc(n*sizeof(double)) ;
    x = malloc(n*sizeof(double)) ;
    
    /*create  A in sed Format and a vector b*/
    A_cs = cs_lapmat_p1_square(N) ;
    A_sed = sed_compress(A_cs); 

    /*solve the problem with jacobi method*/
    index anzIt = 30 ;
    double *w = malloc(n*sizeof(double)) ;
    for (index i = 0 ; i< n ; i++) 
    {
        x [i] = 0 ;
        b [i] = 1 ;
    }
    double anaTimeSum = 0 ;
    index counter = 0;
    for (index i = 0 ; i < anzIt; i++) 
    {   
        TIME_SAVE (0) ;
        sed_jacobi (A_sed , b , x , w) ;
        TIME_SAVE (1);
        anaTimeSum += TIME_ELAPSED (0 , 1) ; 
        counter = i;
    }
    anaAvgTime[0] = anaTimeSum / counter ;
    anaAnzIt[0] = counter ;
    
    
    /*solve the problem with guass-seidel method*/ 
    for (index i = 0 ; i< n ; i++)
    {
        x[i] = 0 ; 
    }
    counter = 0 ; 
    anaTimeSum = 0 ; 
    for (index i = 0 ; i < anzIt; i++)
    {
        TIME_SAVE (0) ;
        sed_gauss_seidel(A_sed, b, x, w);
        TIME_SAVE (1) ;
        anaTimeSum = TIME_ELAPSED (0 , 1) ;
        counter = i;    
    }
    anaAvgTime[1] = anaTimeSum / counter ;
    anaAnzIt [1] = counter ;

    /*incomplete cholesky fkatorization*/
    printf("\n=============================\n") ;
    printf("Incomplete Cholesky Faktorization\n") ;
    
    for (index i = 0 ; i < n ; i++)
    {
        x[i] = 0 ; 
    }        
    index lmax = ((A_sed->nzmax - A_sed->n)/2) + A_sed->n + 1;
    L_sed = sed_alloc(A_sed->n, lmax, 1);                                       
    
    
    TIME_SAVE (3) ;
    
    sed_icholesky(A_sed, L_sed);      
    TIME_SAVE (4) ;
    sed_print(L_sed,0);
    /*    
    printf("b = \n");
    for (index i = 0 ; i < n ; i++){
        printf("%f  ", b[i]) ;
    }
    printf("\n"); 
    */
    printf("\n===========================\n") ;
    printf("Aplly the pcg method with icf\n") ;
    double tol = 1e-10 ;
    sed_ccg(A_sed , L_sed , b , x ,n , tol) ;
   


    /*write analysis in a txt file*/
    FILE *f;    
    
    f = fopen("Analysis_Problem_1.txt", "w");
    fprintf(f, "Verfahren: Jacobi Gauss-Seidel PCG_ICF PCG_ICNE Multigrid");
    fprintf(f, "\nnumber of  Iterationen ") ;
    for (index i = 0 ; i < 5 ; i++)
    {
        fprintf(f, "%d ",(int) anaAnzIt[i]);
    }
    fprintf(f,"\naverage time per Iteration ");
    for (index i = 0 ; i < 5 ; i++)
    {
        fprintf(f, "%f ", anaAvgTime[i]);
    }
    fclose(f) ; 
        
        
        
        
        
    free (anaAvgTime) ;    
    free (anaAnzIt) ;
        
    free (b) ; 
    free (x) ;
    free (w) ;
    sed_free (L_sed) ;
    sed_free (A_sed);
    cs_free (A_cs ) ; 
}
