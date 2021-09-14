#include "hpc.h"

void print_buffer_double(double* buffer, int len){

    for(size_t i=0; i<len; i++){
        printf("%2.2f ", buffer[i]);
        
    }
    printf("\n");
}
void print_buffer_int(index* buffer, int len){

    for(size_t i=0; i<len; i++){
        printf("%ld ", buffer[i]);
    }
    printf("\n");
}