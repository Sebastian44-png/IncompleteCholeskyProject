#include "hpc.h"

index sed_ccg (sed *A , sed *L ,  double *b , double *x , index maxIt , double tol)
{
    /*check input*/
    if(!A || !L || !b || !x)
    {
        return (0) ;
    }

    index An ;
    An = A->n ;
    
    /*local variables*/
    double *r ;
    double *r_next ;
    double *p ;
    double *Ap ;
    double alpha ;
    double beta ; 
    double roh ;
    double roh_next ;
    double *Lap ;
    
    r = malloc(An * sizeof(double)) ;
    r_next = malloc(An * sizeof(double)) ;
    p = malloc(An * sizeof(double)) ;
    Ap = malloc(An * sizeof(double)) ;
    Lap = malloc(An *sizeof(double)) ;
    
    for (index i = 0 ; i < An ; i++)
    {
        r [i] = 0. ;
        r_next [i] = 0. ;
        p [i] = 0. ;
    }


    /*calculate the first residual r_0 = b - A*x(0) and store the result in r*/
    sed_spmv(A , x , r) ;
    for (index i = 0 ; i < An ; i++)
    {
        r [i] = b [i] - r [i] ;
    }


    /*calculate the transformed residual r = L^-1 * r_0*/
    sed_forwardInsertion(L, r_next , r) ;
    for (index i = 0 ; i < An ; i++)
    {
        r [i] = r_next [i] ;
        r_next [i] = 0 ;
    }

    /*calculate p_0 = L-T * r_0*/
    sed_backwardInsertion(L , p , r) ;
    
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
        
        sed_spmv(A , p , Ap);
        
        
        /*calculate alpha*/
        alpha = hpc_dot(Ap, p, An) ;
        if (alpha == 0)
        {
            return (0) ;
        }
        alpha = roh / alpha ;

        /*calculate the next solution x*/
        hpc_scal(x , p , alpha , An) ;
        
        /*calculate the next residuum*/
        for (index i = 0 ; i < An ; i++)
        {
            Lap [i] = 0 ;
        }
        sed_forwardInsertion (L , Lap , Ap) ;
        hpc_scal(r , Lap , -alpha , An) ;
        if(hpc_dot(r,r,An) < tol)
        {
            return (k);
        }
        /*caluclate the next roh*/
        roh_next = hpc_dot(r , r , An) ;

        /*check if roh_next == 0*/
        if(roh_next == 0)
        {
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
        sed_backwardInsertion(L, Ap, r) ; 
        for (index i = 0 ; i < An ; i++)
        {
            p [i] = Ap [i] + beta * p [i] ;
        }
        /*update roh*/
        roh = roh_next ;
    }



    /*memory release*/
   
    free (r) ;
    free (r_next) ;
    free (p) ;
    free (Ap) ;
    free (Lap) ; 
   
   return (1) ;  
}
