#include "ncomp/types.h"
#include "ncomp/util.h"
#include "ncomp_build/util.hpp"
#include <vector>

/* this name mangling works with gfortran, haven't tested other compilers */
extern "C" void dlinint2_(int *,double *,int *,double *,
                          double *,int *,int *,double *,int *,
                          double *,double *,double *, double *,
                          int *,double *,int *,int *);

extern "C" int linint2(const ncomp_array *xi, const ncomp_array *yi, const ncomp_array *fi,
                       const ncomp_array *xo, const ncomp_array *yo, ncomp_array *fo,
                       int icycx, int iopt) {
  /* Adapted from linint2_W() in linint2W.c */
  long nxi, nyi, nxi2, nfi, nxo, nyo, nfo, size_leftmost, size_fo;

  nxi = xi->shape[xi->ndim - 1];
  nyi = yi->shape[yi->ndim - 1];
  nxo = xo->shape[0];
  nyo = yo->shape[0];
  nxi2 = nxi + 2;
  nfi = nxi * nyi;
  nfo = nxo * nyo;

  int inxi = (int)nxi;
  int inyi = (int)nyi;
  int inxo = (int)nxo;
  int inyo = (int)nyo;
  int inxi2 = (int)nxi2;
  int infi = (int)nfi;
  int info = (int)nfo;

  int ier = 0;
  int tmp_ier;

  long index_xi = 0, index_yi = 0, index_fi = 0, index_fo = 0;

  /* handle missing values */
  int has_missing_fi = fi->has_missing;
  ncomp_missing *missing_fi;
  missing_fi = (ncomp_missing *)&(fi->msg);
  double missing_dfi;
  float missing_rfi;
  coerce_missing(fi->type, has_missing_fi, missing_fi, &missing_dfi,
                 &missing_rfi);

  /* create temporary work arrays */
  std::vector<double> xiw(xi->shape[xi->ndim - 1] + 2);
  std::vector<double> fxiw(xi->shape[xi->ndim - 1] + 2);
  std::vector<double> tmp_fo(xo->shape[0] * yo->shape[0]);

  /* create temporary output arrays */
  double *tmp_xo = nullptr, *tmp_yo = nullptr;
  std::vector<double> tmp_vector_xo;
  std::vector<double> tmp_vector_yo;

  if (xo->type == NCOMP_DOUBLE) {
    tmp_xo = static_cast<double *>(xo->addr);
  } else {
    tmp_vector_xo.resize(nxo);
    tmp_xo = tmp_vector_xo.data();
    convert_to<double>(xo->addr, nxo, 0, xo->type, tmp_xo);
  }

  if (yo->type == NCOMP_DOUBLE) {
    tmp_yo = static_cast<double *>(yo->addr);
  } else {
    tmp_vector_yo.resize(nyo);
    tmp_yo = tmp_vector_yo.data();
    convert_to<double>(yo->addr, nyo, 0, yo->type, tmp_yo);
  }

  size_leftmost = 1;
  for (auto i = 0; i < fi->ndim - 2; i++)
    size_leftmost *= fi->shape[i];
  size_fo = size_leftmost * nfo;

  /*
   * Coerce input arrays to double if necessary.
   */
  double *tmp_xi = nullptr, *tmp_yi = nullptr, *tmp_fi = nullptr;
  std::vector<double> tmp_vector_xi;
  std::vector<double> tmp_vector_yi;
  std::vector<double> tmp_vector_fi;

  if (xi->type != NCOMP_DOUBLE) {
    tmp_vector_xi.resize(nxi);
    tmp_xi = tmp_vector_xi.data();
  }
  if (yi->type != NCOMP_DOUBLE) {
    tmp_vector_yi.resize(nyi);
    tmp_yi = tmp_vector_yi.data();
  }
  if (fi->type != NCOMP_DOUBLE) {
    tmp_vector_fi.resize(nfi);
    tmp_fi = tmp_vector_fi.data();
  }

  for (auto i = 0; i < size_leftmost; i++) {
    if (xi->ndim > 1 || i == 0) {
      if (xi->type != NCOMP_DOUBLE)
        convert_to<double>(xi->addr, nxi, index_xi, xi->type, tmp_xi);
      else
        tmp_xi = &((double *)xi->addr)[index_xi];
    }
    if (yi->ndim > 1 || i == 0) {
      if (yi->type != NCOMP_DOUBLE)
        convert_to<double>(yi->addr, nyi, index_yi, yi->type, tmp_yi);
      else
        tmp_yi = &((double *)yi->addr)[index_yi];
    }

    if (fi->type != NCOMP_DOUBLE)
      convert_to<double>(fi->addr, nfi, index_fi, fi->type, tmp_fi);
    else
      tmp_fi = &((double *)fi->addr)[index_fi];

    tmp_ier = 0;

    dlinint2_(&inxi, tmp_xi, &inyi, tmp_yi, tmp_fi, &icycx, &inxo, tmp_xo,
              &inyo, tmp_yo, tmp_fo.data(), xiw.data(), fxiw.data(), &inxi2,
              &missing_dfi, &iopt, &tmp_ier);

    ier |= tmp_ier;
    if (tmp_ier) {
      for (auto j = 0; j < nfo; j++) {
        if (fi->type == NCOMP_DOUBLE) {
          ((double *)fo->addr)[index_fo + j] = missing_dfi;
        } else {
          ((float *)fo->addr)[index_fo + j] = missing_rfi;
        }
      }
    } else {
      if (fi->type == NCOMP_DOUBLE) {
        convert_to<double>(tmp_fo.data(), nfo, 0, NCOMP_DOUBLE,
                           ((double *)fo->addr) + index_fo);
      } else {
        convert_to<float>(tmp_fo.data(), nfo, 0, NCOMP_DOUBLE,
                          ((float *)fo->addr) + index_fo);
      }
    }
    if (xi->ndim > 1)
      index_xi += nxi;
    if (yi->ndim > 1)
      index_yi += nyi;
    index_fi += nfi;
    index_fo += nfo;
  }
  if (fi->has_missing) {
    fo->has_missing = 1;
  } else {
    fo->has_missing = 0;
  }

  if (fi->type == NCOMP_DOUBLE) {
    fo->msg.msg_double = missing_dfi;
  } else {
    fo->msg.msg_float = missing_rfi;
  }

  return ier;
}
