#include "hpc.h"
/* print a banded matrix; use %g for integers to avoid differences with index */
index bnd_print (const bnd *A, index brief)
{
  index i, j, m, n, nzmax, nz, *Ap, *Ai ;
  double *Ad, *Ax ;
  
  if (!A) { printf ("(null)\n") ; return (0) ; }
  n = A->n ; m = A->m ; Ax = A->x ;
    
  printf ("%g-by-%g, bandwidth = %g\n", (double) n, (double) n, (double) m) ;
  printf ("matrix entries row wise \n");  
  for (i = 0 ; i < n ; i++)
  {
    for (j = 0 ; j < 2*m+1 ; j++)
    {
      printf (" %5.3g", Ax[i*(2*m+1)+j]) ;
      if (brief && j > 10) { printf ("  ...\n") ; return (1) ; }
    }
    printf ("\n") ;
  }
  return (1) ;
}

