#include "hpc.h"

#include<time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

double kappa( double x[2], index typ )
{
  return ( 1.0 );
}

double F_vol( double x[2], index typ )
{
  return ( 0.0 );
}

double g_Neu( double x[2], index typ )
{
  return ( x[0] * x[1] );
}

double u_D( double x[2])
{
//  return ( 0.0 );
return ( x[0] * x[1] );
//return (1.0);
}


int main (int argc, char **argv)
{
    index n, k, ncoord, nelem, nbdry, nfixed, nedges, total, *bdry,
          cnt = 0, N = 0, MAX_ITER = 100;
    char *Pdir = "../Problem/", fname[64];
    double *b, *x, *w, *Coord, x1[2], x2[2], m[2];
    mesh **H, *T ;
    sed **A;
    double *bR;
    sed *AR;
    double *dBuff;
    index ft ;
    index fptr ;
    index *fixed ;

    
    printf("\n========================================\n");
    if (argc < 2 ){ printf("Problem not specified\n"); return(1); } 
    sprintf(fname,"%s%s",Pdir,argv[1]); /* get problem as parameter */
    if (argc>2){                        /* get no. of refinements  */ 
      if ( (atoi(argv[2]) >0) & (atoi(argv[2]) < 13)) N = atoi(argv[2]);
    } 
    printf("Load data form %s, no. refinements = %g\n", fname, (double) N);   
        
    /* Allocate memory for hierachy */
    H = malloc ( (N+1) * sizeof(mesh));
    A = malloc ( (N+1) * sizeof(sed*));
    
    /* Load problem */
    H[0] = mesh_load (fname);              /* load geometry */
    mesh_getEdge2no(H[0]->nelem, H[0]->elem, &H[0]->nedges, &H[0]->edge2no);
    H[0]->fixed = mesh_getFixed(H[0]->ncoord, H[0]->bdry, H[0]->nbdry, &H[0]->nfixed);
    printf("\nInit mesh  # dofs =  %10g\n",(double)  H[0]->ncoord+H[0]->nedges);
    /* Build stiffness matrix, refine mesh and create hierachy  */ 
    k = 0;
    while(1)
    {  
      A[k] = sed_nz_pattern(H[k]) ;            /* get pattern of matrix */
      if (!A[k]) return(1);
      if ( !sed_buildS(H[k],A[k]) ) return(1); /* assemble coefficient matrix */    
      if (k >= N) break;
      H[k+1] = mesh_refine(H[k]);
      mesh_getEdge2no(H[k+1]->nelem, H[k+1]->elem,
                      &H[k+1]->nedges, &H[k+1]->edge2no);
      H[k+1]->fixed = mesh_getFixed(H[k+1]->ncoord, H[k+1]->bdry, 
                                   H[k+1]->nbdry, &H[k+1]->nfixed);
      k++;
    }
    printf("Final mesh # dofs =  %10g\n",(double)  H[N]->ncoord+H[N]->nedges);
    printf("# refinements     =  %10g\n",(double)  N);

    n = A[N]->n;
    x = calloc (n, sizeof(double));       /* get workspace for sol*/
    w = calloc (n, sizeof(double));       /* get temporary workspace */
    b = calloc (n, sizeof(double));       /* get workspace for rhs*/
    mesh_buildRhs(H[N], b, F_vol, g_Neu); /* build rhs (volume and Neumann data */
    
    /* incorporate Dirichlet data */
    ncoord = H[N]->ncoord ; nelem = H[N]->nelem ; nbdry = H[N]->nbdry ; 
    bdry = H[N]->bdry; H[N]->fixed = mesh_getFixed(ncoord, bdry, nbdry, &H[N]->nfixed);
    nfixed = H[N]->nfixed ; nedges = H[N]->nedges; Coord = H[N]->coord;

    for ( k = 0; k < nbdry; k++)
    {
      if (!bdry[4*k+3])
      {
        x1[0] = Coord[2 * bdry[4*k]];   x1[1] = Coord[2 * bdry[4*k]+1];
        x2[0] = Coord[2 * bdry[4*k+1]]; x2[1] = Coord[2 * bdry[4*k+1]+1];
        m[0] = ( x1[0] + x2[0] ) / 2.0 ;
        m[1] = ( x1[1] + x2[1] ) / 2.0 ;
        x[bdry[4*k]  ] = u_D(x1);
        x[bdry[4*k+1]] = u_D(x2);
        x[ncoord + bdry[4*k+2]] = u_D(m);
      }
    }

    bR = malloc (n * sizeof(double)) ;
    dBuff = malloc (n * sizeof(double)) ;
    fixed = H[N]->fixed ;
    ft = fixed [0];
    fptr = 0;
    for (k = 0 ; k < n ; k++)
    {
      bR [k] = b [k] ;
      if (k == ft)
      {
        dBuff [k] = -1 * x [k] ;
        fptr++ ;
        if (fptr < nfixed)
        {
          ft = fixed [fptr] ;
        }
      }
      else
      {
        dBuff [k] = 0 ;
      }
    }
    sed_gaxpy(A [N], dBuff, bR) ;
    k = 0 ;
    ft = fixed [0];
    fptr = 0;
    for (index kb; kb < n; kb++) 
    {
      if (kb == ft)
      {
        fptr++ ;
        if (fptr < nfixed)
        {
          ft = fixed [fptr] ;
        }
      }
      else
      {
        bR [k] = bR [kb] ;
        k++ ;
      }
    }

    /* AR*x = bR is LSE for standard CG solver */
    AR = sed_reduceS(A [N], fixed, nfixed) ;
    
    /*vorbereitung für Analyse*/
    
    index *anaAnzIt ;
    double *anaAvgTime ; 
    double *errorJacobi ;
    double *errorGauss ;
    double *errorICF ;
    double *errorICNE ;
    double * errorMultigrid ;
    
    index maxIt = 2 * AR->n ;
    double tol  = 1e-16 ;

    anaAnzIt = malloc(5*sizeof(index)) ;
    anaAvgTime = malloc(5*sizeof(double)) ;
    errorJacobi = malloc(maxIt*sizeof(double)) ;
    errorGauss = malloc(maxIt*sizeof(double)) ;
    errorICF = malloc(maxIt*sizeof(double)) ;
    errorICNE = malloc(maxIt*sizeof(double));
    errorMultigrid = malloc(maxIt*sizeof(double)) ;
    
    /*PCG with multigrid as preconditioner*/
    double *x0 = malloc (AR ->n *sizeof(double)) ;
    for (k = 0 ; k < AR->n ; k ++)
    {
        x0[k] = 0. ;
    }
    
    TIME_SAVE (0) ;
    anaAnzIt [4] = sed_pcg_mg(AR, A, bR, x0, tol , maxIt, H, N, 1, 1, 1, errorMultigrid);  
    TIME_SAVE (1) ;
    anaAvgTime [4] = TIME_ELAPSED (0 , 1) / anaAnzIt [4] ; 
    
    /*PCG with jacobi as preconditioner*/
    for (k = 0 ; k < AR->n ; k++)
    {
        x0 [k] = 0. ;
    }
    TIME_SAVE (2) ;
    anaAnzIt [0] = sed_cg_jacobi (AR , bR , x0 , maxIt , tol, errorJacobi) ;
    TIME_SAVE (3) ;
    anaAvgTime [0] = TIME_ELAPSED (2 , 3) / anaAnzIt [0] ; 

    /*PCG with gauss-seidel as preconditioner*/
    for (k = 0 ; k < AR->n ; k++)
    {
        x0 [k] = 0. ;
    }
    TIME_SAVE (4) ;
    anaAnzIt [1] = sed_cg_gauss_seidel (AR , bR , x0 , maxIt , tol, errorGauss) ;
    TIME_SAVE (5) ;
    anaAvgTime [1] = TIME_ELAPSED (4 , 5) / anaAnzIt [1] ; 


    /*PCG with icf as preconditioner*/
    sed *L_sed_icf ;
    L_sed_icf = sed_alloc(AR->n, AR->nzmax , 1) ;
    sed_icholesky (AR , L_sed_icf) ;      
    for (k = 0 ; k < AR->n ; k++)
    {
        x0 [k] = 0. ;
    }
    TIME_SAVE (6) ;
    anaAnzIt [2] = sed_ccg (AR , L_sed_icf , bR , x0 , maxIt , tol, errorICF) ;
    TIME_SAVE (7) ;
    anaAvgTime [2] = TIME_ELAPSED (6 , 7) / anaAnzIt [2] ; 
    
    
    sed_free(L_sed_icf) ;
    /*PCG with icne as preconditinoer*/
    for (k = 0 ; k < AR->n ; k++)
    {
        x0 [k] = 0. ;
    }
    double alpha = 2 ; 
    sed * L_sed_icne = sed_alloc(AR->n, 0 , 1) ;
    sed * A_LLt = sed_alloc(AR->n, 0, 1) ;
    
    printf("=========  AR   =====\n") ;
    sed_print(AR , 0) ;
    printf("=========  A_LLt======\n") ;
    //sed_L_to_LLt(AR, A_LLt) ;
    //print_buffer_double(A_LLt->x, A_LLt->nzmax) ;
    //sed_print(A_LLt , 0);

    sed_icne0(AR , alpha , L_sed_icne);
    print_buffer_double(L_sed_icne->x, L_sed_icne->nzmax+5) ;
    print_buffer_int(L_sed_icne->i, L_sed_icne->nzmax+5) ;
   
    printf("=========icne0 matrix ==========\n") ;
    sed_print(L_sed_icne , 0);

    TIME_SAVE (8) ;
    anaAnzIt [3] = sed_ccg (AR , L_sed_icne , bR , x0 , maxIt , tol, errorICNE) ;
    TIME_SAVE (9) ;
    anaAvgTime [3] = TIME_ELAPSED (8 , 9) / anaAnzIt [3] ; 
    

    /*AR sed ist in L Form
     * bR in dim von AR (ergebnisvektor)*/
    
    
    
    
    
    
    
    //errStep = sed_pcg_mg_jac(AR, A, bR, x, 1e-10, 50, H, N, 1, 1, 1);  
    // cnt = sed_cg(AR , bR, x, 50, 1e-10);
                                   
    //total = ncoord*2*sizeof(double) 
    //      + (7*nelem+4*nbdry+nedges*2)*sizeof(index);
    
    //printf ("Total :       %12.6g MByte\n", (double) total/1024./1024.);
    
    /*analysis in a txt file*/
    FILE *f ;
    f = fopen("Analyse_Problem_1.txt","w") ;
    fprintf(f,"error_Jacobi error_Gauss error_ICF error_ICNE error_Multigrid\n") ;
    for (index i = 0 ; i < maxIt ; i++)
    {
        fprintf(f,"%.3e %.3e %.3e %.3e %.3e\n" , errorJacobi [i] , errorGauss [i] , errorICF [i] , errorICNE [i] , errorMultigrid [i]) ;
    }
    fclose(f);
    /*memory release*/
    sed_free(AR);
    free(bR);
    free(dBuff);
    for (k = 0 ; k <= N ; k++)
    {
        mesh_free(H[k]) ;
        sed_free(A[k]) ;
    }
    printf("hello4 \n");
    free(H); free(A);
    printf("hello5 \n");
    free (anaAnzIt) ;
    free (anaAvgTime) ;
    free (errorJacobi) ;
    free (errorGauss) ;
    free (errorICF) ;
    free (errorICNE) ;

    return (0) ;
}
