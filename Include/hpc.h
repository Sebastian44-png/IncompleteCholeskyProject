#ifndef _HPC_H
#define _HPC_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#include <errno.h>
#include <string.h>

#include <stdbool.h>

#define index ptrdiff_t

/* --- primary HPC routines and data structures ------------------------- */

typedef struct cs_sparse /* matrix in compressed-row/col or triplet form */
{
    index nzmax; /* maximum number of entries */
    index m;     /* number of rows */
    index n;     /* number of columns */
    index *p;    /* col/row pointers (size n+1) or col indices (size nzmax) */
    index *ind;  /* row/col indices, size nzmax */
    double *x;   /* numerical values, size nzmax */
    index nz;    /* # of entries in triplet matrix, 
                       * -1 for compressed-col, -2 for compressed-row */
} cs;

typedef struct gem_full /* general matrix form, entries stored row wise */
{
    index m;   /* number of rows */
    index n;   /* number of columns */
    double *x; /* numerical values */
} gem;

typedef struct sky_pack /* sym. matrix in sky storage form */
{
    index   n ;       /* number of rows/columns          */
    index  *p ;       /* col pointers (size n)         */
    double *d ;       /* diagonal entries (size n)       */
    double *x ;       /* off-diagonal entries, size p[n] */
} sky ;

typedef struct sed_sparse /* matrix in sparse matrix in compressed col. */
{                         /* with extracted diagonal storage form      */
    index nzmax;          /* maximum number of entries */
    index n;              /* number of rows/columns          */
    index *i;             /* col pointers and row indices    */
    double *x;            /* numerical values, size i[n] */
} sed;

typedef struct mesh_data  /* mesh */
{
    index ncoord ;    /* number of coordinates  */
    index nelem ;     /* number of elements   */
    index nedges ;    /* number of edges  */
    index nbdry ;     /* number of boundary elements  */
    index nfixed;     /* number of fixed nodes ????    */
    double *coord ;   /* coordinates (x1,y1,x2,y2, ... ,x_ncoord,y_ncoord) */
    index *elem ;     /* elements ([e1,e2,e3,m1,m2,m3,t1], ... ) */
    index *edge2no ;  /*  */
    index *bdry ;     /* bdry ([e1,e2,m1,t1], [e3,e4,m2,t2], ...) */
    index *fixed ;    /* bdry ([e1,e2,m1,t1], [e3,e4,m2,t2], ...) */
} mesh ;

/* utilities */
void *hpc_realloc(void *p, index n, size_t size, index *ok);
double hpc_cumsum(index *p, index *c, index n);
int emax(index size, index *values);

/* multigrid scheme */
index hpc_mg(sed **A, double *b, double *x, double tol, index maxit,
             mesh **H, index nLevel, index pre, index post, index gamma);
index hpc_mg_cycle(sed **A, mesh **H, index nLevel, 
                   double **b, double **x, double **r,
                   index pre, index post, index gamma);
index hpc_mg_jac(sed **A, double *b, double *x, double tol, index maxit,
             mesh **H, index nLevel, index pre, index post, index gamma);
index hpc_mg_cycle_jac(sed **A, mesh **H, index nLevel, 
                   double **b, double **x, double **r,
                   index pre, index post, index gamma);
void hpc_prol(double *x, index nx, index *edgeno, index nEdges, double *y);
void hpc_prol_quad(double *x, double *y, index *elem, index nC, index nT, index nE);
void hpc_rest(double *x, index *edgeno, index nEdges, double *y, index ny);
index hpc_reduceFixed(double *x, double *xShort, index n, const index *fixed, const index nFixed);
index hpc_expandFixed (double *x, double *xLong, index n, const index *fixed, const index nFixed);



/* gem format */
gem *gem_alloc(index n, index m);
gem *gem_free(gem *A);
index gem_gausssol(gem *A, double *x);
index gem_gauss(gem *A);
index gem_lu(gem *A);
index gem_gaxpy (const gem *A, const double *x, double *y);
index gem_spmv(const gem *A, const double *x, double *y);
gem *gem_compress(const cs *T);
index gem_print(const gem *A, index brief);

/* cs format */
cs *cs_alloc(index m, index n, index nzmax, index values, index typ);
cs *cs_load(FILE *f, index issym);
cs *cs_compress(const cs *T, index typ);
index cs_realloc(cs *A, index nzmax);
index cs_entry(cs *T, index i, index j, double x);
index cs_print(const cs *A, index brief);
cs *cs_free(cs *A);
cs *cs_done(cs *C, void *w, void *x, index ok);

index cs_spmv(const cs *A, const double *x, double *y);
cs *cs_lapmat_p1_square(index m);

/* skyline format */
sky *sky_alloc(index n, index nzmax);
sky *sky_free(sky *A);
sky *sky_load(FILE *f);
index sky_print(const sky *A, index brief);
sky *sky_compress(const cs *T);
index sky_spmv(const sky *A, const double *x, double *y);
index sky_cholesky(sky *A);

/* sed format */
sed *sed_alloc(index n, index nzmax, index values);
index sed_realloc(sed *A, index nzmax);
sed *sed_free(sed *A);
sed *sed_done(sed *C, void *w, void *x, index ok);
sed *sed_compress(const cs *A);
index sed_print(const sed *A, index brief);
index sed_spmv(const sed *A, const double *x, double *y);

index sed_icholesky(sed *A, sed *L);
index sed_isoccupied(sed *A, index column, index row);
index sed_forwardInsertion(sed *L, double *x, double *b);
index sed_backwardInsertion(sed *L, double *x, double *b);
index sed_ccg(sed *A, sed *L, double *b, double *x, index maxIt, double tol, double *error);
index sed_cg_jacobi(sed *A, double *b, double *x, index maxIt, double tol , double *error);
index sed_cg_gauss_seidel(sed *A, double *b, double *x, index maxIt, double tol , double *error);
index sed_ILU(sed *A) ;
index sed_MILU(sed *A, double alpha) ;

index sed_icne0 (sed *A, double alpha, sed* L);
double* sed_find_in_column(sed* A, index row, index col);
index sed_L_to_LLt(sed* L, sed* LLt);

index sed_gauss_seidel(const sed *A, const double *b, double *xk, double *w);
index sed_gs(const sed *A, const double *b, double *x, double *w, index forward);
index sed_jacobi (const sed *A, const double *b, double *xk, double *w );
index sed_richardson (const sed *A, const double *b, double *xk, double *w, const double omega);
index sed_cg (const sed *A, double *b, double *x, index maxIt, double tol);
index sed_pcg_mg(sed *A, sed **Amg, double *b, double *x, double tol, index maxIt,
             mesh **H, index nLevel, index pre, index post, index gamma, double *error);
index sed_pcg_mg_jac(sed *A, sed **Amg, double *b, double *x, double tol, index maxIt,
             mesh **H, index nLevel, index pre, index post, index gamma, double *error);

index sed_gaxpy (const sed *A, const double *x, double *y);
index sed_gs_constr (const sed *A, const double *b, double *x, double *w, 
                     index *fixed, index nFixed, index forward);
index sed_jacobi_constr(const sed *A, const double *b, double *x, double *w, 
                     index *fixed, index nFixed);
index sed_dupl (sed *A);

/* For MG examples */
void stima_laplace3(double p1[2], double p2[2], double p3[2],
                    index  typ, double dx[6], double ax[9]);
sed *sed_nz_pattern(mesh *M) ; 
index sed_buildS(mesh *M, sed *T);
sed *sed_reduceS(const sed *S, const index *fixed, const index nFixed);

/* mesh operations */
mesh *mesh_alloc (index ncoord, index nelem, index nbdry);
mesh *mesh_free (mesh *M);
mesh *mesh_load (char *fname);
index *mesh_getFixed(const index nCoord, const index *bdry, 
                     const index nBdry, index *nFixed);
index mesh_print (const mesh *M, index brief);
mesh *mesh_refine(mesh *In);
index mesh_getEdge2no(const index nElem, const index *Elem, 
                      index *nEdges, index **edge2no);
void mesh_buildRhs(const mesh *M, double *b, double (*f)(double *, index), 
                   double (*g)(double *, index));

/* hpc utils */
double hpc_dot(const double *x, const double *y, index n);
index hpc_scal(double* y, const double *x, double alpha, index n);

/* utils */
void print_buffer_int(index* buffer, int len);
void print_buffer_double(double* buffer, int len);
index copy_buffer(const double* a, double *b, index n);

/* for mg/FEM, declaration in demo files */
double kappa( double x[2], index typ );
double F_vol( double x[2], index typ );

#define HPC_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define HPC_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define HPC_CSC(A) (A && (A->nz == -1))
#define HPC_CSR(A) (A && (A->nz == -2))
#define HPC_TRIPLET(A) (A && (A->nz >= 0))
#endif
