#include "hpc.h"
cs *cs_lapmat_p1_square (index m){
  
    index i, j;
    cs *T;
    
    T = cs_alloc (m*m, m*m, (5*m-4)*m, 1, 0) ;
    if (!T) return (NULL) ;                             /* check inputs */

    for ( i=0 ; i<m*m ; i++) cs_entry (T, i, i,  4.0);

    for ( i=0 ; i<m*(m-1) ; i++){
      cs_entry (T, i, i+m, -1.0);
      cs_entry (T, i+m, i, -1.0);
    }

    for ( i=0 ; i<m; i++)
    {
      for ( j=0 ; j<m-1 ; j++){
        cs_entry (T, i*m+j, i*m+j+1, -1.0);
        cs_entry (T, i*m+j+1, i*m+j, -1.0);
      }
    }
    return (T);
}
