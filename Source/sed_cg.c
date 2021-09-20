#include "hpc.h"

/* f = A * x + alpha*y */
index sed_cg ( const sed *A, double *b, double *x, index maxIt, double tol , double * error ){
    //check input
    if (!A || !b || !x || !error)
    {
        return (0) ;
    }
    index *ind, n;
    double *Ax;

    ind = A->i;
    Ax = A->x;
    n = A->n;

    // helper variables
    double *r = malloc(n*sizeof(double));
    double *p = malloc(n*sizeof(double));
    double *Ap = malloc(n*sizeof(double));

    double rho;
    double alpha;
    double beta;
    double rho_prev;
    
   // first redisual: r(0) = b - Ax(0)
    sed_gaxpy(A, x, r);
    for ( index i = 0; i < n; i++)
    {
        r [i] = b[i] - r [i];
    }
    
    copy_buffer(r, p ,n);
    rho = hpc_dot(r, r, n);

    for(index k = 0; k < maxIt; k++){

        // actual error as the square of the 2Norm of r
        error[k] = hpc_dot(r , r , n) ;
        if(error[k] < tol)
        {
            free(r);
            free(p);
            free(Ap);
            return(k);
        }
        
        for(index i = 0; i < n; i++)
        {
            Ap [i] = 0;
        }
        sed_gaxpy(A, p, Ap);

        alpha = rho / hpc_dot(Ap, p, n); 
        
        hpc_scal(x, p,  alpha, n); // x_(k) = x_(k-1) + alpha * p
        hpc_scal(r, Ap, -alpha, n); // r_(k) = r(k-1) - alpha * p

        // store rho(k-1), as it's needed for later calculations
        rho_prev = rho;
        rho = hpc_dot(r, r, n); 

        if(rho == 0) 
        {
            free(r);
            free(p);
            free(Ap);
            return (0) ; 
        }

        beta = rho / rho_prev;

        // r_(k) + rho(k) * p(k) / rho(k-1)
        for(index i = 0; i < n; i++)
        {
            p [i] = r [i] + (beta * p [i]);
        }
    }
    /*free memory*/
    free(r);
    free(p);
    free(Ap);

    return (maxIt) ;
}
