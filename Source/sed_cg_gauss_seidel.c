#include "hpc.h"
/*author: Joachim Kröner*/

/*This method solves the linear system Ax = b with Gauss-Mehtod as preconditioner*/
index sed_cg_gauss_seidel (sed *A ,  double *b , double *x , index maxIt , double tol, double *error)
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
    double *z ;
    double *z_next ;
    double *p ;
    double *Ap ;
    double alpha ;
    double beta ; 
    double roh ;
    double *w ;
    
    /*allocate memory*/    
    r = malloc(An * sizeof(double)) ;
    z = malloc(An * sizeof(double)) ;
    z_next= malloc(An * sizeof(double)) ;
    p = malloc(An * sizeof(double)) ;
    Ap = malloc(An * sizeof(double)) ;
    w = malloc(An * sizeof(double)) ;
    for (index i = 0 ; i < An ; i++)
    {
        r [i] = 0. ;
        z [i] = 0. ;
        z_next [i] = 0. ;
        p [i] = 0. ;
    }


    /*calculate the first residual r_0 = b - A*x(0) and store the result in r*/
    sed_gaxpy(A , x ,r) ;
    for (index i = 0 ; i < An ; i++)
    {
        r [i] = b [i] - r [i] ;
    }


    /*calculate the preconditoned residual z = M^-1 * r_0 as an Gauss-Seidel Iteration*/
    sed_gs(A, r, z , w , 1 );
    sed_gs(A, r, z , w , 0 );
    
    /*calculate p*/
    for (index i = 0 ; i < An ; i++)
    {
        p [i] = z [i] ;
    }

    /*calcluate roh*/
    roh = hpc_dot(r , z , An) ;
    
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
            free (r) ;
            free (z_next) ;
            free (z) ;
            free (p) ;
            free (Ap) ;
            free (w) ; 
            return (0) ;
        }
        alpha = roh / alpha ;

        /*calculate the next solution x*/
        hpc_scal(x , p , alpha , An) ;
        
        
        /*calculate the next residuum*/
        hpc_scal(r , Ap, -alpha , An) ;
        
        /*save actuall residual*/
        error [k] = hpc_dot(r , r , An);
        
        /*check if abort criterion is reached */
        if(error [k] < tol)
        {
            free (r) ;
            free (z_next) ;
            free (z) ;
            free (p) ;
            free (Ap) ;
            free (w) ; 
            return (k);
        }
        
        sed_gs(A, r, z_next , w , 1 );
        sed_gs(A, r, z_next , w , 0 );
        
        /*calculate beta*/
        for (index i = 0 ;i < An ; i++)
        {
            z [i] = z_next [i] - z [i] ;
        }
        beta = hpc_dot(r, z , An) / roh ;
        
        for (index i = 0 ; i < An ; i++)
        {
            p [i] = beta * p [i] + z_next [i] ;
            z [i] = z_next [i] ;
            z_next [i] = 0. ;
        }

        /*update roh*/
        roh = hpc_dot(r , z , An) ;
    }

    /*memory release*/
    free (r) ;
    free (z_next) ;
    free (z) ;
    free (p) ;
    free (Ap) ;
    free (w) ; 
   
   return (maxIt) ; 
} 
