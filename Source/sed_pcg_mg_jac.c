#include "hpc.h"

/* Preconditioned Conjugate gradient method with multigrid preconditioner
sed **A needs to be grid hierarchy */
double *sed_pcg_mg_jac (sed *A, sed **Amg, double *b, double *x, double tol, index maxIt,
             mesh **H, index nLevel, index pre, index post, index gamma)
{
    index n ;
    index nFull ;
    index nFixed ;
    index *fixed ;

    n = A->n ;
    nFull = Amg [nLevel]->n ;
    nFixed = H [nLevel]->nfixed ;
    fixed = H [nLevel]->fixed ;

    // helper variables
    double *r = malloc (n * sizeof (double)) ;
    double *rExp = malloc (nFull * sizeof(double)) ;
    double *zExp = malloc (nFull*sizeof(double)) ;
    double *z = malloc (n*sizeof(double)) ; // to store preconditioned residual
    double *z_prev = malloc (n*sizeof(double)) ;
    double *p = malloc (n*sizeof(double)) ;
    double *Ap = malloc (n*sizeof(double)) ;
    double *errStep = malloc (maxIt * sizeof (double)) ;

    double rho ;
    double alpha ;
    double beta ;
    double rho_prev ;
    double error ;
    double tmp ;
    
    // r(0) = b - Ax(0)
    for ( index i = 0; i < n; i++)
    {
        r [i] = 0 ;
    }
    sed_gaxpy (A, x, r) ;
    for (index i = 0 ; i < n ; i++)
    {
        r [i] = b[i] - r [i] ;
    }
    for (index i = 0 ; i < nFull ; i++)
    {
        rExp [i] = 0 ;
        zExp [i] = 0 ;
    }
    
    hpc_expandFixed(r, rExp, nFull, fixed, nFixed) ;
    hpc_mg_jac(Amg, rExp, zExp, 0, 1, H, nLevel, pre, post, gamma) ;
    hpc_reduceFixed(zExp, z, nFull, fixed, nFixed) ;
    copy_buffer(z, p ,n) ;
    rho_prev = hpc_dot(r, z, n) ;

    for(index k = 0 ; k < maxIt ; k++){
        for(index i = 0 ; i < n ; i++)
        {
            Ap[i] = 0 ;
        }
        for (index i = 0 ; i < nFull ; i++)
        {
            zExp [i] = 0 ;
        }
        sed_gaxpy (A, p, Ap) ;
        //sed_spmv(A [nLevel], p, Ap);
        //printf("spmv done");
        alpha = rho_prev / hpc_dot (Ap, p, n) ;
        
        hpc_scal (x, p,  alpha, n) ;
        hpc_scal (r, Ap, -alpha, n) ;

        error = 0;
        // return if ||r||_2 < tol  
        for (index i = 0 ; i < n ; i++)
        {
            error += r[i] * r[i] ;
        }
         printf("Error: %g\n", error);
        errStep [k] = error ;
        if (error < tol)
        {
            free (r) ;
            free (rExp) ;
            free (z) ;
            free (zExp) ;
            free (z_prev) ;
            free (p) ;
            free (Ap) ;
            return (errStep) ;
        }
        print_buffer_double(r, n);
        hpc_expandFixed (r, rExp, nFull, fixed, nFixed) ;
        print_buffer_double(rExp, n);
        hpc_mg_jac (Amg, rExp, zExp, 0, 1, H, nLevel, pre, post, gamma) ;
        print_buffer_double(zExp, n);
        hpc_reduceFixed (zExp, z, nFull, fixed, nFixed) ;
        print_buffer_double(z, n);

        for (index i = 0 ; i < n ; i++)
        {
            tmp = z_prev [i] ;
            z_prev[i] = z[i] ;
            z [i] -= tmp ;
        }
        rho = hpc_dot (r, z, n) ;

        if (rho == 0) 
        { 
            free (r) ;
            free (rExp) ;
            free (zExp) ;
            free (z) ;
            free (z_prev) ;
            free (p) ;
            free (Ap) ;
            free (errStep) ;
            return (NULL) ;
        }

        beta = rho / rho_prev ;

        for (index i = 0 ; i < n ; i++)
        {
            p[i] = z_prev[i] + (beta * p[i]) ;
        }
        rho_prev = hpc_dot(r, z_prev, n) ;
        
        //printf(" x ="); print_buffer_double(x, n);
    }
    free (r) ;
    free (rExp) ;
    free (z) ;
    free (zExp) ;
    free (z_prev) ;
    free (p) ;
    free (Ap) ;
    return (errStep) ;
}