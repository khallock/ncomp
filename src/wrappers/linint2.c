#include <stdlib.h>
#include <stdio.h>
#include "ncomp/types.h"
#include "ncomp/util.h"
#include <time.h>

/* this name mangling works with gfortran, haven't tested other compilers */
extern void dlinint2_(int *,double *,int *,double *,
                      double *,int *,int *,double *,int *,
                      double *,double *,double *, double *,
                      int *,double *,int *,int *);

int linint2(
    const ncomp_array* xi, const ncomp_array* yi, const ncomp_array* fi,
    const ncomp_array* xo, const ncomp_array* yo, ncomp_array* fo,
    int icycx, int iopt
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

    /* handle missing values */
    int has_missing_fi = fi->has_missing;
    ncomp_missing* missing_fi;
    missing_fi = (ncomp_missing*) &(fi->msg);
    double missing_dfi;
    float missing_rfi;
    coerce_missing(fi->type, has_missing_fi, missing_fi, &missing_dfi, &missing_rfi);

    /* create temporary work arrays */
    double* xiw  = (double*)calloc(xi->shape[xi->ndim-1] + 2,sizeof(double));
    double* fxiw = (double*)calloc(xi->shape[xi->ndim-1] + 2,sizeof(double));
    double* tmp_fo = (double*)calloc(xo->shape[0] * yo->shape[0], sizeof(double));
    if(xiw == NULL || fxiw == NULL || tmp_fo == NULL) {
      printf("linint2: Unable to allocate memory for temporary work arrays");
      return(-1);
    }

  if(xo->type == NCOMP_DOUBLE) {
    tmp_xo = (double*)xo->addr;
  } else {
    tmp_xo = (double*)calloc(nxo,sizeof(double));
    _to_double(xo->addr, 0, nxo, xo->type, tmp_xo);
  }

  if(yo->type == NCOMP_DOUBLE) {
    tmp_yo = (double*)yo->addr;
  } else {
    tmp_yo = (double*)calloc(nyo,sizeof(double));
    _to_double(yo->addr, 0, nyo, yo->type, tmp_yo);
  }

  size_leftmost = 1;
  for( i = 0; i < fi->ndim-2; i++ ) size_leftmost *= fi->shape[i];
  size_fo = size_leftmost * nfo;

  if(xi->type != NCOMP_DOUBLE) {
    tmp_xi = (double*)calloc(nxi,sizeof(double));
    if(tmp_xi == NULL) {
      printf("linint2: Unable to allocate memory for coercing xi to double precision");
      return(-1);
    }
  }

  if(yi->type != NCOMP_DOUBLE) {
    tmp_yi = (double*)calloc(nyi,sizeof(double));
    if(tmp_yi == NULL) {
      printf("linint2: Unable to allocate memory for coercing yi to double precision");
      return(-1);
    }
  }

  if(fi->type != NCOMP_DOUBLE) {
    tmp_fi = (double*)calloc(nfi,sizeof(double));
    if(tmp_fi == NULL) {
      printf("linint2: Unable to allocate memory for coercing fi to double precision");
      return(-1);
    }
  }


    for( i = 0; i < size_leftmost; i++ ) {
        if(xi->ndim > 1 || i == 0) {
            if(xi->type != NCOMP_DOUBLE) {
                _to_double(xi->addr, index_xi, nxi, xi->type, tmp_xi);
            }
            else {
                tmp_xi = &((double*)xi->addr)[index_xi];
            }
        }
        if(yi->ndim > 1 || i == 0) {
            if(yi->type != NCOMP_DOUBLE) {
                _to_double(yi->addr, index_yi, nyi, yi->type, tmp_yi);
            }
            else {
                tmp_yi = &((double*)yi->addr)[index_yi];
            }
        }

        if(fi->type != NCOMP_DOUBLE) {
            _to_double(fi->addr, index_fi, nfi, fi->type, tmp_fi);
        }
        else {
            tmp_fi = &((double*)fi->addr)[index_fi];
        }


        tmp_ier = 0;

        dlinint2_(&inxi, tmp_xi, &inyi, tmp_yi, tmp_fi, &icycx, &inxo, tmp_xo, &inyo, tmp_yo, tmp_fo, xiw, fxiw, &inxi2, &missing_dfi, &iopt, &tmp_ier);

        ier |= tmp_ier;
        if(tmp_ier) {
            for(j = 0; j < nfo; j++) {
                if(fi->type == NCOMP_DOUBLE) {
                    ((double*)fo->addr)[index_fo+j] = missing_dfi;
                }
                else {
                    ((float*)fo->addr)[index_fo+j] = missing_rfi;
                }
            }
        }
        else {
            if(fi->type == NCOMP_DOUBLE) {
                _to_double(tmp_fo, 0, nfo, NCOMP_DOUBLE, ((double*)fo->addr) + index_fo);
            } else {
                _to_float(tmp_fo, 0, nfo, NCOMP_DOUBLE, ((float*)fo->addr) + index_fo);
            }
        }
        if(xi->ndim > 1) index_xi += nxi;
        if(yi->ndim > 1) index_yi += nyi;
        index_fi += nfi;
        index_fo += nfo;
    }
    if(fi->has_missing) {
        fo->has_missing = 1;
    } else {
        fo->has_missing = 0;
    }

    if(fi->type == NCOMP_DOUBLE) {
        fo->msg.msg_double = missing_dfi;
    } else {
        fo->msg.msg_float = missing_rfi;
    }

    free(xiw);
    free(fxiw);
    free(tmp_fo);
    return ier;
}
