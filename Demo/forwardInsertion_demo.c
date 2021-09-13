#include "hpc.h"

int main(int argc, char **argv)
{
    index n = 5 ; 
    double *x = malloc(n*sizeof(double));
    double *b = malloc(n*sizeof(double));

    b[0] = 1;
    b[1] = 0;
    b[2] = 1;
    b[3] = 0;
    b[4] = 1;
    
    sed *L; 
    L = calloc(5, 15,1);

    sed_print(L);
    
}
