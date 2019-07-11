#include <stdlib.h>
#include <stdio.h>
#include "ncomp.h"
#include "util.h"
#include <time.h>

/* this name mangling works with gfortran, haven't tested other compilers */
extern void dlinint2_(int *,double *,int *,double *,
                      double *,int *,int *,double *,int *,
                      double *,double *,double *, double *,
                      int *,double *,int *,int *);

int linint2_loop_struct(
    const ncomp_array* xi, const ncomp_array* yi, const ncomp_array* fi,
    const ncomp_array* xo, const ncomp_array* yo, ncomp_array* fo,
    double* xiw, double* fxiw, double* tmp_fo,
    int icycx, double xmsg, int iopt
) {
/* Adapted from linint2_W() in linint2W.c */

    long nxi, nyi, nxi2, nfi, nxo, nyo, nfo, size_leftmost, size_fo;

    nxi  = xi->shape[xi->ndim-1];
    nyi  = yi->shape[yi->ndim-1];
    nxo  = xo->shape[0];
    nyo  = yo->shape[0];
    nxi2 = nxi + 2;
    nfi  = nxi * nyi;
    nfo  = nxo * nyo;

    int inxi = (int) nxi;
    int inyi = (int) nyi;
    int inxo = (int) nxo;
    int inyo = (int) nyo;
    int inxi2 = (int) nxi2;
    int infi = (int) nfi;
    int info = (int) nfo;

    double *tmp_xi = NULL;
    double *tmp_yi = NULL;
    double *tmp_fi = NULL;
    double *tmp_xo, *tmp_yo;
    int ier = 0;
    int tmp_ier;

    long i, j, index_xi = 0, index_yi = 0, index_fi = 0, index_fo = 0;

    tmp_xo = (double*)xo->addr;
    tmp_yo = (double*)yo->addr;

  size_leftmost = 1;
  for( i = 0; i < fi->ndim-2; i++ ) size_leftmost *= fi->shape[i];
  size_fo = size_leftmost * nfo;

    for( i = 0; i < size_leftmost; i++ ) {
        if(xi->ndim > 1 || i == 0) {
            if(xi->type != NCOMP_DOUBLE) {
                _to_double(xi->addr, index_xi, nxi, xi->type, tmp_xi);
            }
            else {
                tmp_xi = &((double*)xi->addr)[index_xi];
            }
            tmp_xi = &((double*)xi->addr)[index_xi];
        }
        if(yi->ndim > 1 || i == 0) {
            if(yi->type != NCOMP_DOUBLE) {
                _to_double(yi->addr, index_yi, nyi, yi->type, tmp_yi);
            }
            else {
                tmp_yi = &((double*)yi->addr)[index_yi];
            }
            tmp_yi = &((double*)yi->addr)[index_yi];
        }

        if(fi->type != NCOMP_DOUBLE) {
            _to_double(fi->addr, index_fi, nfi, fi->type, tmp_fi);
        }
        else {
            tmp_fi = &((double*)fi->addr)[index_fi];
        }

        tmp_fi = &((double*)fi->addr)[index_fi];

        tmp_ier = 0;

        dlinint2_(&inxi, tmp_xi, &inyi, tmp_yi, tmp_fi, &icycx, &inxo, tmp_xo, &inyo, tmp_yo, tmp_fo, xiw, fxiw, &inxi2, &xmsg, &iopt, &tmp_ier);

        ier |= tmp_ier;
        if(tmp_ier) {
            for(j = 0; j < nfo; j++) {
                if(fi->type == NCOMP_DOUBLE) {
                    ((double*)fo->addr)[index_fo+j] = xmsg;
                }
                else {
                    /* FIX THIS ********/
                    ((float*)fo->addr)[index_fo+j] = (float)xmsg;
                }
            }
        }
        else {
            _to_double(tmp_fo, 0, nfo, fo->type, ((double*)fo->addr) + index_fo);
        }
        if(xi->ndim > 1) index_xi += nxi;
        if(yi->ndim > 1) index_yi += nyi;
        index_fi += nfi;
        index_fo += nfo;
    }

    return ier;
}

int linint2(
    const ncomp_array* xi, const ncomp_array* yi, const ncomp_array* fi,
    const ncomp_array* xo, const ncomp_array* yo, ncomp_array* fo,
    int icycx, double xmsg, int iopt
) {
    double* xiw  = (double*)calloc(xi->shape[xi->ndim-1] + 2,sizeof(double));
    double* fxiw = (double*)calloc(xi->shape[xi->ndim-1] + 2,sizeof(double));
    double* tmp_fo = (double*)calloc(xo->shape[0] * yo->shape[0], sizeof(double));
    int ier = linint2_loop_struct(xi, yi, fi, xo, yo, fo,
                                  xiw, fxiw, tmp_fo, icycx, xmsg, iopt);
    free(xiw);
    free(fxiw);
    free(tmp_fo);
    return ier;
}
