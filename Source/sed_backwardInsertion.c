#include "hpc.h"                                                                
                                                                                
/*Backward Insecertion: solves the linear System Lx=b where L is a lower triangular matrix*/
                                                                                
index sed_backwardInsertion(sed *L , double *x , double *b)                
{                                                                               
    if(!L || !x || !b)                                                          
    {                                                                           
        return(0);                                                              
    }                                                                           
                                                                                
    index n ;                                                                   
    index *ind ;                                                                
    double *val ;                                                               
    index temp ;                                                                
                                                                                
    n =  L->n ;                                                                 
    ind = L->i ;                                                                
    val = L->x ;                                                                
    temp = -1 ;                                                                 
                                                                                
    x [n - 1] = b [n - 1] / val [n - 1] ;                                       
                                                                                
    for (index i = n-1 ; i > 0 ; i--)                                           
    {                                                                           
        for (index j = 0 ; j < i ; j++)                                         
        {                                                                       
            temp = sed_isoccupied(L, j, i) ;                                    
            if(temp >= 0)                                                       
            {                                                                   
                x[j] += val [ind[j] +temp] * x [i];                             
            }                                                                   
        }                                                                       
                                                                                
        x [i -1] = (b [i - 1] - x [i - 1]) / val [i - 1] ;                      
    }                                                                           
    return (1) ;                                                                
}                  
