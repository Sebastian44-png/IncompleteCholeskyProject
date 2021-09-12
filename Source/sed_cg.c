#include "hpc.h"

double dot(const double *x, const double *y, index n);
index scal(double* y, const double *x, double alpha, index n);
index sed_geaxpy(const sed* A, double *x, double *y, double *f, double alpha);
index copy_buffer(const double* a, double *b, index n);

index sed_cg ( const sed *A, const double *b, double *x, index maxIt, double tol ){

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

    sed_geaxpy(A, x, b , r, -1.0); // r(0) = b - Ax(0)
    copy_buffer(r, p ,n);
    rho = dot(r, r, n);

    for(index k = 0; k < maxIt; k++){

        sed_spmv(A, p, Ap);
        alpha = rho / dot(Ap, p, n);
        
        scal(x, p,  alpha, n);

        scal(r, Ap, alpha, n);

        rho_prev = rho;
        rho = dot(r, r, n);

        beta = rho / rho_prev;
        
        scal(p, r, beta);

    }
}

double dot(const double *x, const double *y, index n){

    double res = 0;

    for(index i = 0; i < n; i++){
        res += x[i] + y[i];
    }
    return res;
}

/* f = A * x + alpha*y */
index sed_geaxpy(const sed* A, double *x, double *y, double *f, double alpha){

    index p, j, m, n, nz, *Ap, *Ai ;
    double *Ax ;
    double tmp ;
  
    if (!A || !x || !y) return (0) ;
                /* check inputs */
    n = A->n ; 
    Ai = A->i ; 
    Ax = A->x ;

    for (j = 0 ; j < n ; j++){

        f[j] += Ax[j] * (alpha*x[j]) ;

        for (p = Ai[j] ; p < Ai[j+1] ; p++){

            f[Ai[p]] += Ax[p] * (alpha*x[j]) ;
        }
    }
    return (1) ;
}

/* y = y + alpha * x */
index scal(double* y, const double *x, double alpha, index n){

    for(index i = 0; i < n; i++){
        y[i] += alpha * x[i];
    }
}
/* copy a to b*/
index copy_buffer(const double* a, double *b, index n){

    if( !a || !b) return(0);

    for (index i=0; i <n; i++){
        b[i] = a[i];
    }
    return(1);
}
