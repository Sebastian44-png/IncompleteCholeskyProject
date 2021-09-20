#include "hpc.h"
/*author: Joachim Kröner*/

/*check if a certain line in a column is occupied*/
index sed_isoccupied( sed *A , index column , index row)
{
    index isThere = -1 ;
    index *ind ;
    /*check input*/
    if(!A)
    {
        return (-1) ;
    }    
    ind = A->i ;
    /*interation over a column of A */
    for (index i = 0 ; i < ind[column +1] - ind[column] ; i++ )
    {
        /*check if a certain row index there*/
        if(ind[ind[column] + i] == row)
        {
            isThere = i ;
            break ;
        }
    }
    return isThere ;
}
