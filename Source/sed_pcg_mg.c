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
    double *z_prev = malloc(n*sizeof(double));
    double *p = malloc(n*sizeof(double));
    double *Ap = malloc(n*sizeof(double));

    double rho;
    double alpha;
    double beta;
    double rho_prev;
    double error;
    double tmp;
    
    // r(0) = b - Ax(0)
    for ( index i = 0; i < n; i++)
    {
        r [i] = 0;
    }
    sed_spmv(A [nLevel], x, r);
    for ( index i = 0; i < n; i++)
    {
        z [i] = 0;
        r [i] = b[i] - r [i];
    }
    printf("Entering mg");
    hpc_mg(A, r, z, 0, 1, H, nLevel, pre, post, gamma);
        print_buffer_double(z, n);
    printf("First mg done");
    copy_buffer(z, p ,n);
    //rho = hpc_dot(r, z, n);

    for(index k = 0; k < maxIt; k++){
        rho_prev = hpc_dot(r, z, n);
        for(index i = 0; i < n; i++)
        {
            Ap[i] = 0;
            z[i] = 0;
        }

        sed_spmv(A [nLevel], p, Ap);
        //printf("spmv done");
        alpha = rho_prev / hpc_dot(Ap, p, n);
        
        hpc_scal(x, p,  alpha, n);
        hpc_scal(r, Ap, -alpha, n);

        error = 0;
        // return if ||r||_2 < tol  
        for(index i = 0; i < n; i++)
        {
            error += r[i] * r[i];
        }
        printf("Error: %g\n", error);
        if(error < tol){
            free (r);
            free (z);
            free (z_prev);
            free (p);
            free (Ap);
            return(k);
        }
        // preconditioning:
        hpc_mg(A, r, z, 0, 1, H, nLevel, pre, post, gamma);
        //print_buffer_double(z, n);

        //rho_prev = rho;
        for (index i = 0; i < n; i++)
        {
            tmp = z_prev [i];
            z_prev[i] = z[i];
            z [i] -= tmp;
        }
        rho = hpc_dot(r, z, n);

        if(rho == 0) 
        { 
            free (r);
            free (z);
            free (z_prev);
            free (p);
            free (Ap);
            return (-1);
        }

        beta = rho / rho_prev;

        for(index i = 0; i < n; i++)
        {
            p[i] = z_prev[i] + (beta * p[i]);
        }
        
        //printf(" x ="); print_buffer_double(x, n);
    }
    free (r);
    free (z);
    free (z_prev);
    free (p);
    free (Ap);
    return (maxIt);
}