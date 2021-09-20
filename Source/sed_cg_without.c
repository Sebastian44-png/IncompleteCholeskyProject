#include "hpc.h"

index sed_cg_without (sed *A ,  double *b , double *x , index maxIt , double tol, double *error)
{
    /*check input*/
    if(!A || !b || !x)
    {
        return (0) ;
    }

    index An ;
    An = A->n ;
    
    /*local variables*/
    double *r ;
    double *p ;
    double *Ap ;
    double alpha ;
    double beta ; 
    double roh ;
    double roh_next ;
    
    r = malloc(An * sizeof(double)) ;
    p = malloc(An * sizeof(double)) ;
    Ap = malloc(An * sizeof(double)) ;
    
    for (index i = 0 ; i < An ; i++)
    {
        r [i] = 0. ;
        p [i] = 0. ;
    }


    /*calculate the first residual r_0 = b - A*x(0) and store the result in r*/
    sed_gaxpy(A, x , r);
    //sed_spmv(A , x , r) ;
    for (index i = 0 ; i < An ; i++)
    {
        r [i] = b [i] - r [i] ;
    }


    /*calculate p_0 =  r_0*/
    for (index i = 0 ; i < An ; i++)
    {
        p [i] = r [i] ;
    } 
    
    /*calcluate roh*/
    roh = hpc_dot(r , r, An) ;
    
    /*cg iteration*/
    for (index k = 0 ; k < maxIt ; k ++)
    {
        /*calculate the Matrix vektor product*/
        for (index i = 0 ; i < An ; i++) 
        {
            Ap [i] = 0 ;
        }
        
        sed_gaxpy(A , p , Ap) ;
        
        /*calculate alpha*/
        alpha = hpc_dot(Ap, p, An) ;
        if (alpha == 0)
        {
            return (0) ;
        }
        alpha = roh / alpha ;

        /*calculate the next solution x*/
        hpc_scal(x , p , alpha , An) ;
        
        /*save actuall residual*/
        error [k] = hpc_dot(r , r , An);
        
        /*calculate the next residuum*/
        hpc_scal(r , Ap , -alpha , An) ;
        
        /*check if abort criterion is reached */
        if(error[k] < tol)
        {
            free (r) ;
            free (p) ;
            free (Ap) ;
            return (k);
        }

        /*caluclate the next roh*/
        roh_next = hpc_dot(r , r , An) ;

        /*check if roh_next == 0*/
        if(roh_next == 0)
        {
            free (r) ;
            free (p) ;
            free (Ap) ;
            return (0) ;
        }
        /*calculate beta*/
        beta = roh_next / roh ;

        /*calculate the next p 
         * Ap is used as a buffer, in the rest of the loop the values of Ap are not used
         * L^T * Ap = r*/
        for (index i = 0 ; i < An ; i++)
        {
            Ap [i] = 0;
        }
        for (index i = 0 ; i < An ; i++)
        {
            p [i] = r[i] + beta * p [i] ;
        }
        /*update roh*/
        roh = roh_next ;
    }



    /*memory release*/
   
    free (r) ;
    free (p) ;
    free (Ap) ;
   
   return (maxIt) ;  
}
