#include "hpc.h"
/* author: Benjamin Bestler */

/* Preconditioned Conjugate gradient method for solving A*x = b with multigrid preconditioner
 *sed **Amg needs to be grid hierarchy 
 *A does not contain rows/columns corresponding to fixed nodes, the levels of Amg do contain these rows/columns */
index sed_pcg_mg (sed *A, sed **Amg, double *b, double *x, double tol, index maxIt,
             mesh **H, index nLevel, index pre, index post, index gamma, double *error)
{
    // check inputs
    if (!A || !Amg || !b || !x || !H || !error)
    {
        return (0) ;
    }

    index n ;
    index nFull ;
    index nFixed ;
    index *fixed ;

    // problem dimension
    n = A->n ;

    // for adapting data to Multigrid
    nFull = Amg [nLevel]->n ;
    nFixed = H [nLevel]->nfixed ;
    fixed = H [nLevel]->fixed ;

    // Buffers for residual, r for CG, rExp for MG
    double *r = malloc (n * sizeof (double)) ;
    double *rExp = malloc (nFull * sizeof(double)) ;

    // Buffers for residual after preconditioningn
    double *z = malloc (n * sizeof(double)) ;
    double *zExp = malloc (nFull * sizeof(double)) ;
    double *z_prev = malloc (n * sizeof(double)) ;

    // Buffer for p and Ap = A * p
    double *p = malloc (n * sizeof(double)) ;
    double *Ap = malloc (n * sizeof(double)) ;

    // check if mallocs are successfull
    if (!r || !rExp || !z || !zExp || !z_prev || !p || !Ap)
    {
        free (r) ;
        free (rExp) ;
        free (z) ;
        free (zExp) ;
        free (z_prev) ;
        free (p) ;
        free (Ap) ;
        return (0) ;
    }

    double alpha ;
    double beta ;
    double rho ;
    
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

    // Prepare for precdonditioning
    for (index i = 0 ; i < nFull ; i++)
    {
        rExp [i] = 0 ;
        zExp [i] = 0 ;
    }
    hpc_expandFixed(r, rExp, nFull, fixed, nFixed) ;

    // first preconditioning, then adapt data for CG    
    hpc_mg(Amg, rExp, zExp, 0, 1, H, nLevel, pre, post, gamma) ;
    hpc_reduceFixed(zExp, z_prev, nFull, fixed, nFixed) ;

    // copy to p
    copy_buffer(z_prev, p ,n) ;
    rho = hpc_dot(r, z_prev, n) ;

    for (index k = 0 ; k < maxIt ; k++)
    {
        // Update Ap = A * p
        for (index i = 0 ; i < n ; i++)
        {
            Ap[i] = 0 ;
        }

        sed_gaxpy (A, p, Ap) ;

        // update alpha
        alpha = rho / hpc_dot (Ap, p, n) ;
        
        // update solution vector and residual
        hpc_scal (x, p,  alpha, n) ;
        hpc_scal (r, Ap, -alpha, n) ;

        // check for convergence
        error [k] = hpc_dot (r, r, n) ;
        if (error [k] < tol)
        {
            free (r) ;
            free (rExp) ;
            free (z) ;
            free (zExp) ;
            free (z_prev) ;
            free (p) ;
            free (Ap) ;
            return (k) ;
        }

        // prepare for preconditioning
        for (index i = 0 ; i < nFull ; i++)
        {
            zExp [i] = 0 ;
        }
        hpc_expandFixed (r, rExp, nFull, fixed, nFixed) ;

        // preconditioning, then adapt data for CG    
        hpc_mg (Amg, rExp, zExp, 0, 1, H, nLevel, pre, post, gamma) ;
        hpc_reduceFixed (zExp, z, nFull, fixed, nFixed) ;

        for (index i = 0 ; i < n ; i++)
        {
            z_prev [i] = z[i] - z_prev [i] ;
        }

        // update beta, rho in denominator is previous rho
        beta = hpc_dot (r, z_prev , n) / rho ;

        // previous rho not needed anymore -> update rho and check
        rho = hpc_dot(r, z, n) ;

        if (rho == 0) 
        { 
            free (r) ;
            free (rExp) ;
            free (zExp) ;
            free (z) ;
            free (z_prev) ;
            free (p) ;
            free (Ap) ;
            free (error) ;
            return (0) ;
        }

        // Update p and store z in z_prev
        for (index i = 0 ; i < n ; i++)
        {
            p [i] = z [i] + (beta * p [i]) ;
            z_prev [i] = z [i] ;
        }
    }
    free (r) ;
    free (rExp) ;
    free (z) ;
    free (zExp) ;
    free (z_prev) ;
    free (p) ;
    free (Ap) ;
    return (maxIt) ;
}