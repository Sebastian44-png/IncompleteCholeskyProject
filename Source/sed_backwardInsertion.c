#include "hpc.h"                                                                
/*author: Joachim Kröner*/

/*Backward Insecertion: solves the linear System L^Tx=b where L is given as  a lower triangular matrix*/
index sed_backwardInsertion(sed *L , double *x , double *b)                
{   
    /*check input*/    
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
    
    /*calculate last entry of x*/    
    x [n - 1] = b [n - 1] / val [n - 1] ;                                       
    
    /*Iteration over the columns of L*/    
    for (index i = n-1 ; i > 0 ; i--)                                           
    {   
        /*Iteration over the rows of L*/        
        for (index j = 0 ; j < i ; j++)                                         
        {   
            /*check if element is zero*/            
            temp = sed_isoccupied(L, j, i) ;                                    
            if(temp >= 0)                                                       
            {                                                                   
                x[j] += val [ind[j] +temp] * x [i];                             
            }                                                                   
        }                                                                       
        /*calculate the next entry of x*/                                                                        
        x [i -1] = (b [i - 1] - x [i - 1]) / val [i - 1] ;                      
    }                                                                           
    return (1) ;                                                                
}                  
