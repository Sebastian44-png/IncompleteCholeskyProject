#include "hpc.h"
#include<time.h>
#include<sys/time.h>"

struct timeval probelm1_dim1 [50] ;
struct timeval problem1_dim2 [50] ;

#define TIME_SAVE(j) (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j, k) (1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))


int main (void)
{
    /*init a Problem*/
    TIME_SAVE (0) ;
    index N ;
    cs *A_cs ;
    sed *A_sed ;
    sed *L_sed ;                                                                
    double *b ;
    double *x ;
    index lnzmax ;

    /*allocate memory*/
    TIME_SAVE (1) ;
    b = malloc(n*sizeof(double)) ;
    x = malloc(n*sizeof(double)) ;
    
    /*initialize memory to 0*/
    for (index i = 0 ; i < n ; i++)
    {
        b [i] = 0 ;
        x [i] = 0 ;
    }

    /*release memory*/
        
    
    
    
    
}
