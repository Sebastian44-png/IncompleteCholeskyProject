#include "hpc.h"

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

    for (index i = 0 ; i < ind[column +1] - ind[column] ; i++ )
    {
        if(ind[ind[column] + i] == row)
        {
            isThere = i ;
            break ;
        }
    }
    return isThere ;
}
