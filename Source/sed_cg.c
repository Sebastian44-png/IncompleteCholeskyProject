#include "hpc.h"

double hpc_dot(const double *x, const double *y, index n);
index hpc_scal(double* y, const double *x, double alpha, index n);
index sed_geaxpy(const sed* A, double *x, double *y, double *f, double alpha);
index copy_buffer(const double* a, double *b, index n);

/* f = A * x + alpha*y */
index sed_cg ( const sed *A, double *b, double *x, index maxIt, double tol ){

    index *ind, n;
    double *Ax;

    ind = A->i;
    Ax = A->x;
    n = A->n;

    // helper variables
    double *r = malloc(n*sizeof(double));
    double *p = malloc(n*sizeof(double));
    double *Ap =malloc(n+sizeof(double));

    double rho;
    double alpha;
    double beta;
    double rho_prev;
    double error;
    
    sed_geaxpy(A, x, b , r, -1.0); // r(0) = b - Ax(0)
    copy_buffer(r, p ,n);
    rho = hpc_dot(r, r, n);

    for(index k = 0; k < maxIt; k++){

        // return if ||r||_2 < tol  
        for(index i = 0; i < n; i++)
        {
            error += r[i] * r[i];
        }
        if(error < tol){
            return(1);
        }
        
        for(index i = 0; i < n; i++)
        {
            Ap[i] = 0;
        }

        sed_spmv(A, p, Ap);

        alpha = rho / hpc_dot(Ap, p, n);
        
        hpc_scal(x, p,  alpha, n);
        
        hpc_scal(r, Ap, -alpha, n);

        rho_prev = rho;
        rho = hpc_dot(r, r, n);

        if(rho == 0) { continue; }

        beta = rho / rho_prev;

        for(index i = 0; i < n; i++)
        {
            p[i] = r[i] + (beta * p[i]);
        }
        
        printf(" x ="); print_buffer_double(x, n);
    }
}

double hpc_dot(const double *x, const double *y, index n){

    double res = 0;

    for(index i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

/* y = b + alpha A * x */
index sed_geaxpy(const sed* A, double *x, double *b, double *y, double alpha){

    index p, j, m, n, nz, *Ap, *Ai ;
    double *Ax ;
    double tmp ;

    n = A->n ; 
    Ai = A->i ; 
    Ax = A->x ;


    for (j = 0 ; j < n ; j++)
    {
        y [j] = b [j] ;
    }
    
    for (j = 0 ; j < n ; j++)
    {

        y [j] += alpha * (Ax [j] * x [j]) ;

        for (p = Ai[j] ; p < Ai[j+1] ; p++)
        {
            y [Ai [p] ] += alpha * ( Ax [p] * x [j]) ;
        }
        
    }
    return (1) ;
}

/* y = y + alpha * x */
index hpc_scal(double* y, const double *x, double alpha, index n)
{

    if (!x || !y) return (0) ;

    for(index i = 0; i < n; i++)
    {       
        y [i] += alpha * x [i];
    }
    return (1) ;
}
/* copy a to b*/
index copy_buffer(const double* a, double *b, index n){

    for (index i=0; i <n; i++){
        b [i] = a[i];
    }
    return(1);
}
