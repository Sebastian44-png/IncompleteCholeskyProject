/*
Utilities for working with buffers
Author: Sebasian Acerbi
*/

#include "hpc.h"

void print_buffer_double(double* buffer, int len)
{

    for( size_t i = 0; i < len; i++){
        printf ( "%2.2f " , buffer [i] );
        
    }
    printf("\n") ;
}

void print_buffer_int(index* buffer, int len)
{

    for(size_t i = 0; i < len; i++){
        printf("%ld ", buffer [i]) ;
    }
    printf("\n");
}

/* copy a to b*/
index copy_buffer(const double* a, double *b, index n)
{

    for (index i = 0; i < n; i++){
        b [i] = a [i] ;
    }
    return(1) ;
}