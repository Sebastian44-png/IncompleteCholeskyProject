#include "hpc.h"

//calculate the dot product of x and y
double hpc_dot(const double *x, const double *y, index n){

    double res = 0;

    for(index i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}