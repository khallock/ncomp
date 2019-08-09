#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ncomp.h>
#include <math.h>

int main(int argc, char** argv){
    srand(time(0));

    size_t n = 0;
    size_t i, j, k;
    size_t fi_index, fo_index;
    int ier;

    if(argc > 1) {
        n = atoi(argv[1]);
        n = pow(2, (int) (log(n)/ log(2)));
    }
    if(n == 0)
        n = 512;

    size_t in_size = (n / 2) + 1;
    size_t out_size = n + 1;
    size_t dims_fi[] = {3, in_size, in_size};
    size_t dims_fo[] = {3, out_size, out_size};

    double *xi = (double *) malloc(sizeof(double) * dims_fi[2]);
    double *yi = (double *) malloc(sizeof(double) * dims_fi[1]);
    double *fi = (double *) malloc(sizeof(double) * dims_fi[0] * dims_fi[1] * dims_fi[2]);
    double *xo = (double *) malloc(sizeof(double) * dims_fo[2]);
    double *yo = (double *) malloc(sizeof(double) * dims_fo[1]);
    double *fo = (double *) malloc(sizeof(double) * dims_fo[0] * dims_fo[1] * dims_fo[2]);

    for(i = 0; i <= n; i+=2) {
        xi[i / 2] = i;
        yi[i / 2] = i;
    }

    for(i = 0; i <= n; i++) {
        xo[i] = i;
        yo[i] = i;
    }

    for(i = 0; i < dims_fi[0] * dims_fi[1] * dims_fi[2]; i++) {
        fi[i] = rand();
    }

    ncomp_array* ncomp_xi = ncomp_array_alloc((void*) xi, NCOMP_DOUBLE, 1, &in_size);
    ncomp_array* ncomp_yi = ncomp_array_alloc((void*) yi, NCOMP_DOUBLE, 1, &in_size);
    ncomp_array* ncomp_fi = ncomp_array_alloc((void*) fi, NCOMP_DOUBLE, 3, dims_fi);
    ncomp_array* ncomp_xo = ncomp_array_alloc((void*) xo, NCOMP_DOUBLE, 1, &out_size);
    ncomp_array* ncomp_yo = ncomp_array_alloc((void*) yo, NCOMP_DOUBLE, 1, &out_size);
    ncomp_array* ncomp_fo = ncomp_array_alloc((void*) fo, NCOMP_DOUBLE, 3, dims_fo);

    ier = linint2(ncomp_xi, ncomp_yi, ncomp_fi, ncomp_xo, ncomp_yo, ncomp_fo, 0, 0);

    for(i = 0; i < dims_fi[0]; i++) {
        for(j = 0; j < dims_fi[1]; j++) {
            for(k = 0; k < dims_fi[2]; k++) {
                fi_index = i * dims_fi[1] * dims_fi[2] + j * dims_fi[2] + k;
                fo_index = i * dims_fo[1] * dims_fo[2] + j * 2 * dims_fo[2] + k * 2;
                if( ((double *)ncomp_fi->addr)[fi_index] != ((double *)ncomp_fo->addr)[fo_index] ) {
                    printf("ncomp_fi[%ld] != ncomp_fo[%ld]\n", fi_index, fo_index);
                    return -1;
                }
            }
        }
    }

    return ier;
}
