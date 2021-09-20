#include "hpc.h"

/*
Author: Sebastian Acerbi
*/

/* convert L matrix to L + L-transposed*/
index sed_L_to_LLt(sed* L, sed* LLt)
{
    if (!L) return (0);
    
    index l_n; 
    index* l_ind; 
    index l_nzmax; 
    double* l_Ax; 

    index* llt_ind;
    index llt_nzmax;
    double* llt_Ax;
    index lt_index;

    l_n = L->n;
    l_ind = L->i;
    l_nzmax = L->nzmax;
    l_Ax = L->x;

    llt_nzmax = l_nzmax + (l_nzmax - l_n);

    sed_realloc(LLt, llt_nzmax);
    
    llt_ind = LLt->i;
    llt_nzmax = LLt->nzmax;
    llt_Ax = LLt->x;
    
    index *counts = malloc(l_n * sizeof(index));
    index *col_counts = malloc((l_n+1) * sizeof(index));
    for (index i = 0; i < l_n; i++){
        counts[i] = 0;
        col_counts[i] = 0;
    }

    for (index i = 0; i < l_n + 1; i++)
    {
        LLt->x[i] = L->x[i];
        llt_ind[i] = l_ind[i];
        llt_Ax[i] = l_Ax[i];
    }
    
    // calculate col counts
    for (index j = 0; j < l_n-1; j++) 
    {
        for (index i = l_ind[j]; i < l_ind[j+1]; i++)
        {
            counts[l_ind[i]] ++;
        }
    }
    hpc_cumsum(col_counts, counts, l_n+1);

    for(index i= 0; i < l_n+1; i++){
        llt_ind[i] += col_counts[i];
    }
    llt_ind[l_n + 1] = 0;

    for(index i=0; i<l_n; i++){
        counts[i] = 0;
    }
    
    for (index j = 0; j < l_n-1; j++)
    {
        for(index i = l_ind[j]; i < l_ind[j+1]; i++)
        {
            lt_index = llt_ind[l_ind[i]] + counts[l_ind[i]];

            llt_Ax[lt_index] = l_Ax[i];
            llt_ind[lt_index] = j;
            counts[l_ind[i]]++;

            llt_ind[llt_ind[j] + counts[j]] = l_ind[i];
            llt_Ax[llt_ind[j] + counts[j]] = l_Ax[i];
            counts[j]++;
        }
    }
    free(counts) ;
    
    return (1);
}