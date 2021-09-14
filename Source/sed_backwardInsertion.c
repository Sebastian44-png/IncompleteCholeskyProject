#include "hpc.h"                                                                

/*Backward Insecertion: solves the linear System Lx=b where L is a lower triangular matrix*/

index sed_forwardInsertion(sed *L , double *x , double *b)                      
{                                                                               
    if(!L || !x || !b)                                                          
    {                                                                           
        return(0);                                                              
    }                                                                           
 
    index n ;                                                                   
    index *ind ;
    double *val ;                                                               
                                                                             
    n =  L->n ;                                                                  
    ind = L->i ;                                                                
    val = L->x ;                                                                
    x [n - 1] = b [n - 1] / val [ n-1 ] ;                                                   
                                                                         
    for (index i = n-1 ; i >= 0 ; i--)                                            
    {                                                                           
        for (index ptr = ind[i] ; ptr < ind[i + 1] ; ptr ++)                      
        {                                                                       
            x [ind [ptr]] += val [ptr] * x [i] ;                            
        }                                                                       
        
        x [i -1] = (b [i - 1] - x [i - 1]) / val [i - 1] ;                                     
    }                                                                           
}            
