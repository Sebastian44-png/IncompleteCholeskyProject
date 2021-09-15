#include "hpc.h"

/* Preconditioned Conjugate gradient method with multigrid preconditioner
sed **A needs to be grid hierarchy */
index sed_pcg_mg(sed **A, double *b, double *x, double tol, index maxIt,
             mesh **H, index nLevel, index pre, index post, index gamma)
{
    index *Aind, n;
    double *Ax;

    Aind = A [nLevel]->i;
    Ax = A [nLevel]->x;
    n = A [nLevel]->n;

    // helper variables
    double *r = malloc(n*sizeof(double));
    double *z = malloc(n*sizeof(double)); // to store preconditioned residual
    double *p = malloc(n*sizeof(double));
    double *Ap = malloc(n+sizeof(double));

    double rho;
    double alpha;
    double beta;
    double rho_prev;
    double error;
    
    sed_geaxpy(A [nLevel], x, b , r, -1.0); // r(0) = b - Ax(0)
    copy_buffer(r, p ,n);
    hpc_mg(A, r, z, 0, 1, H, nLevel, pre, post, gamma);
    rho = hpc_dot(r, z, n);

    for(index k = 0; k < maxIt; k++){

        // return if ||r||_2 < tol  
        for(index i = 0; i < n; i++)
        {
            error += r[i] * r[i];
        }
        if(error < tol){
            return(k);
        }
        
        for(index i = 0; i < n; i++)
        {
            Ap[i] = 0;
        }

        sed_spmv(A [nLevel], p, Ap);

        alpha = rho / hpc_dot(Ap, p, n);
        
        hpc_scal(x, p,  alpha, n);
        
        hpc_scal(r, Ap, -alpha, n);

        hpc_mg(A, r, z, 0, 1, H, nLevel, pre, post, gamma);

        rho_prev = rho;
        rho = hpc_dot(r, z, n);

        if(rho == 0) { continue; }

        beta = rho / rho_prev;

        for(index i = 0; i < n; i++)
        {
            p[i] = r[i] + (beta * p[i]);
        }
        
        printf(" x ="); print_buffer_double(x, n);
    }
}