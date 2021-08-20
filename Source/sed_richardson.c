#include "hpc.h"
// One step of Richardson Iteration. Solve: Ax = b, w as auxiliary vector
index sed_richardson (const sed *A, const double *b, double *xk, double *w, const double omega )
{
   index n; index *ind;
   double* Ax;
   
   n = A->n; ind = A->i; Ax = A->x;
  
   // xk+1 = omega * bi - omega*(xi*(1-aii));
   for(size_t i=0; i<n; i++){
      w[i] = xk[i];
      xk[i] -= (omega * Ax[i]* w[i]);
      xk[i] += omega * b[i];
   }
   
   //iterate through rows
   for(size_t j = 0; j < n; j++){
      //iterate throug cols 
      for(size_t i = ind[j]; i < ind [j+1]; i++){
         
         //Row Index -> ind[i]
         xk[ind[i]] -= omega * Ax[i] * w[j];
      }
   }
   return 1;
}   

