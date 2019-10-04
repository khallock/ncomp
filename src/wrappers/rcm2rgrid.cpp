#include "ncomp/types.h"
#include "ncomp/util.h"
#include "ncomp_internal/util.hpp"
#include <iostream> // for cerr
#include <vector>

extern "C" void drcm2rgrid_(int *ingrid, int *inlat2d, int *inlon2d,
                            double *tmp_lat2d, double *tmp_lon2d,
                            double *tmp_fi, int *inlat1d, double *tmp_lat1d,
                            int *inlon1d, double *tmp_lon1d, double *tmp_fo,
                            double *missing_dfi, int *tmp_ncrit, int *tmp_opt,
                            int *ier);

extern "C" void drgrid2rcm_(int *ingrid, int *inlat1d, int *inlon1d,
                            double *tmp_lat1d, double *tmp_lon1d,
                            double *tmp_fi, int *inlat2d, int *inlon2d,
                            double *tmp_lat2d, double *tmp_lon2d,
                            double *tmp_fo, double *missing_dfi, int *tmp_ncrit,
                            int *tmp_opt, int *ier);

//////////////////////////////////////////////////////////////////////////////////////////

extern "C" int rcm2rgrid(const ncomp_array *lat2d, const ncomp_array *lon2d, const ncomp_array *fi,
			 const ncomp_array *lat1d, const ncomp_array *lon1d, ncomp_array *fo) {
  // Other variables
  size_t nlon2d, nlat2d, nlat1d, nlon1d, nfi, nfo;

  nlat2d = lat2d->shape[0];
  nlon2d = lat2d->shape[1];
  nlat1d = lat1d->shape[0];
  nlon1d = lon1d->shape[0];
  nfi = nlon2d * nlat2d;
  nfo = nlat1d * nlon1d;

  // Compute the total size of the input/output arrays.
  size_t ngrid = 1;
  for (auto i = 0; i < fi->ndim - 2; i++)
    ngrid *= fi->shape[i];
  size_t size_fi = ngrid * nfi;
  size_t size_fo = ngrid * nfo;

  /**************** START: BASIC SANITY CHECKS ****************/
  if (lat2d->shape[0] != lon2d->shape[0] ||
      lat2d->shape[1] != lon2d->shape[1]) {
#if DEBUG
    std::cerr << "ERROR rcm2rgrid: The input lat/lon grids must be the same "
      "size ! \n";
#endif
    return 1;
  }

  if (nlon2d <= 1 || nlat2d <= 1 || nlat1d <= 1 || nlon1d <= 1) {
#if DEBUG
    std::cerr << "ERROR rcm2rgrid: The input/output lat/lon grids must have at "
      "least 2 elements ! \n";
#endif
    return 1;
  }

  // Check dimensions of fi.
  if (fi->ndim < 2) {
#if DEBUG
    std::cerr << "ERROR rcm2rgrid: fi must be at least two dimensions !\n";
#endif
    return 1;
  }
  if (fi->shape[fi->ndim - 2] != nlat2d || fi->shape[fi->ndim - 1] != nlon2d) {
#if DEBUG
    std::cerr << "ERROR rcm2rgrid: The rightmost dimensions of fi must "
      "be (nlat2d x nlon2d), where nlat2d and nlon2d are the "
      "dimensions of the lat2d/lon2d arrays !\n";
#endif
    return 1;
  }

  // Test input dimension sizes.
  if((nlon2d > INT_MAX) || (nlat2d > INT_MAX) || (ngrid > INT_MAX) ||
     (nlon1d > INT_MAX) || (nlat1d > INT_MAX)) {
#if DEBUG
    std::cerr << "ERROR rcm2rgrid: one or more input dimension sizes is greater than INT_MAX !\n";
#endif
    return 1;
  }

  /**************** END: BASIC SANITY CHECKS ****************/

  int inlon2d, inlat2d, ingrid, inlon1d, inlat1d;

  inlon2d = static_cast<int>(nlon2d);
  inlat2d = static_cast<int>(nlat2d);
  ingrid = static_cast<int>(ngrid);
  inlon1d = static_cast<int>(nlon1d);
  inlat1d = static_cast<int>(nlat1d);

  /* handle missing values */
  int has_missing_fi = fi->has_missing;
  ncomp_missing *missing_fi;
  missing_fi = (ncomp_missing *)&(fi->msg);
  double missing_dfi;
  float missing_rfi;
  coerce_missing(fi->type, has_missing_fi, missing_fi, &missing_dfi,
                 &missing_rfi);

  // Allocate space for temporary output array.
  // NOTE: since the output array (fo) type and memory is already
  //		cython, just allocate tmp_fo accordingly.
  double *tmp_fo = nullptr;
  std::vector<double> tmp_fo_vec;
  if (fo->type == NCOMP_DOUBLE) {
    tmp_fo = &((double *)fo)[0];
  } else {
    tmp_fo_vec.resize(size_fo);
    tmp_fo = tmp_fo_vec.data();
  }

  // Coerce input arrays to double if necessary.
  double *tmp_lat2d, *tmp_lon2d, *tmp_lat1d, *tmp_lon1d, *tmp_fi;

  tmp_lat2d = coerce_input_T<double>(lat2d->addr, lat2d->type, nfi, 0, nullptr,
                                     nullptr);
  tmp_lon2d = coerce_input_T<double>(lon2d->addr, lon2d->type, nfi, 0, nullptr,
                                     nullptr);
  tmp_lat1d = coerce_input_T<double>(lat1d->addr, lat1d->type, nlat1d, 0,
                                     nullptr, nullptr);
  tmp_lon1d = coerce_input_T<double>(lon1d->addr, lon1d->type, nlon1d, 0,
                                     nullptr, nullptr);
  tmp_fi = coerce_input_T<double>(fi->addr, fi->type, size_fi, has_missing_fi,
                                  &missing_fi, &missing_dfi);

  // Force opt to 0 and ncrit to 1, since they are not used yet.
  int tmp_opt = 0;
  int tmp_ncrit = 1;
  int ier;

  /**************** START: FORTRAN CALL ****************/
  drcm2rgrid_(&ingrid, &inlat2d, &inlon2d, tmp_lat2d, tmp_lon2d, tmp_fi,
              &inlat1d, tmp_lat1d, &inlon1d, tmp_lon1d, tmp_fo, &missing_dfi,
              &tmp_ncrit, &tmp_opt, &ier);
  /**************** END: FORTRAN CALL ****************/

  if (ier) {
    if (ier == 1) {
#if DEBUG
      std::cerr
	<< "ERROR rcm2rgrid: not enough points in input/output array ! \n";
#endif
      return 1;
    }
    if (2 <= ier && ier <= 5) {
#if DEBUG
      std::cerr << "ERROR rcm2rgrid: lat2d, lon2d, lat1d, lon1d must be "
	"monotonically increasing ! \n";
#endif
      return 1;
    }

    set_subset_output_missing(fo->addr, 0, fo->type, size_fo, missing_dfi);
  } else {
    if (fo->type != NCOMP_DOUBLE) {
      convert_to<float>(tmp_fo, size_fo, 0, NCOMP_DOUBLE, ((float *)fo->addr));
    }
  }

  return ier;
}

//////////////////////////////////////////////////////////////////////////////////////////

extern "C" int rgrid2rcm(const ncomp_array* lat1d, const ncomp_array* lon1d, const ncomp_array* fi,
			 const ncomp_array* lat2d, const ncomp_array* lon2d, ncomp_array* fo)
{
  // Other variables
  size_t nlon2d, nlat2d, nlat1d, nlon1d, nfi, nfo;

  nlat2d = lat2d->shape[0];
  nlon2d = lat2d->shape[1]; // same as dsizes_lon2d[1]
  nlat1d = lat1d->shape[0];
  nlon1d = lon1d->shape[0];

  // Compute the total number of elements in our arrays.
  nfi = nlat1d * nlon1d;
  nfo = nlon2d * nlat2d;

  // Compute the total size of the input/output arrays.
  size_t ngrid = 1;
  for(auto i = 0; i < fi->ndim - 2; i++) ngrid *= fi->shape[i];
  size_t size_fi = ngrid * nfi;
  size_t size_fo = ngrid * nfo;

  /**************** START: BASIC SANITY CHECKS ****************/
  if(lat2d->shape[0] != lon2d->shape[0] ||
     lat2d->shape[1] != lon2d->shape[1]) {
#if DEBUG
    std::cerr << "ERROR rgrid2rcm: The output lat/lon grids must be the same size ! \n";
#endif
    return 1;
  }

  if(nlon2d <= 1 || nlat2d <= 1 || nlat1d <= 1 || nlon1d <= 1) {
#if DEBUG
    std::cerr << "ERROR rgrid2rcm: The input/output lat/lon grids must "
      "have at least 2 elements ! \n";
#endif
    return 1;
  }

  // Check dimensions of fi.
  if (fi->ndim < 2) {
#if DEBUG
    std::cerr
      << "ERROR rgrid2rcm: fi must be at least two dimensions !\n";
#endif
    return 1;
  }

  if (fi->shape[fi->ndim - 2] != nlat1d ||
      fi->shape[fi->ndim - 1] != nlon1d) {
#if DEBUG
    std::cerr << "ERROR rgrid2rcm: The rightmost dimensions of `fi` must "
      "be (nlat1d x nlon1d), where nlat1d and nlon1d are the "
      "dimensions of the lat1d/lon1d arrays !\n";
#endif
    return 1;
  }

  // Test input dimension sizes.
  if((nlon2d > INT_MAX) || (nlat2d > INT_MAX) || (ngrid > INT_MAX) ||
     (nlon1d > INT_MAX) || (nlat1d > INT_MAX)) {
#if DEBUG
    std::cerr << "ERROR rgridrcm: one or more input dimension sizes is greater than INT_MAX !\n";
#endif
    return 1;
  }
  /**************** END: BASIC SANITY CHECKS ****************/

  int inlon2d, inlat2d, ingrid, inlon1d, inlat1d;

  inlon2d = (int) nlon2d;
  inlat2d = (int) nlat2d;
  ingrid = (int) ngrid;
  inlon1d = (int) nlon1d;
  inlat1d = (int) nlat1d;

  // Coerce missing values.
  int has_missing_fi = fi->has_missing;
  ncomp_missing* missing_fi;
  missing_fi = (ncomp_missing*) &(fi->msg);
  double missing_dfi;
  float missing_rfi;
  coerce_missing(fi->type, has_missing_fi, missing_fi, &missing_dfi, &missing_rfi);

  // Allocate space for temporary output array.
  // NOTE: since the output array (fo) type and memory is already
  //		cython, just allocate tmp_fo accordingly.
  double *tmp_fo = nullptr;
  std::vector<double> tmp_fo_vec;
  if(fo->type == NCOMP_DOUBLE) {
    tmp_fo = &((double*)fo)[0];
  }
  else {
    tmp_fo_vec.resize(size_fo);
    tmp_fo = tmp_fo_vec.data();
  }

  // Coerce input arrays to double if necessary.
  double *tmp_lat2d, *tmp_lon2d, *tmp_lat1d, *tmp_lon1d, *tmp_fi;

  tmp_lat2d = coerce_input_T<double>(lat2d->addr, lat2d->type, nfo,	 0, nullptr, nullptr);
  tmp_lon2d = coerce_input_T<double>(lon2d->addr, lon2d->type, nfo,	 0, nullptr, nullptr);
  tmp_lat1d = coerce_input_T<double>(lat1d->addr, lat1d->type, nlat1d, 0, nullptr, nullptr);
  tmp_lon1d = coerce_input_T<double>(lon1d->addr, lon1d->type, nlon1d, 0, nullptr, nullptr);
  tmp_fi    = coerce_input_T<double>(fi->addr, fi->type, size_fi, has_missing_fi, &missing_fi, &missing_dfi);

  // Force opt to zero and ncrit to 1, since they are not used yet.
  int tmp_opt = 0;
  int tmp_ncrit = 1;
  int ier;

  /**************** START: FORTRAN CALL ****************/
  drgrid2rcm_(&ingrid, &inlat1d, &inlon1d, tmp_lat1d, tmp_lon1d,
	      tmp_fi, &inlat2d, &inlon2d, tmp_lat2d,
	      tmp_lon2d, tmp_fo, &missing_dfi,
	      &tmp_ncrit, &tmp_opt, &ier);
  /**************** END: FORTRAN CALL ****************/

  if(ier) {
    if(ier == 1) {
#if DEBUG
      std::cerr << "ERROR rgrid2rcm: not enough points in input/output array ! \n";
#endif
      return 1;
    }
    if(2 <= ier && ier <= 5) {
#if DEBUG
      std::cerr << "ERROR rgrid2rcm: lat2d, lon2d, lat1d, lon1d must be monotonically increasing ! \n";
#endif
      return 1;
    }

    set_subset_output_missing(fo->addr, 0, fo->type, size_fo,
			      missing_dfi);
  } else {
    if (fo->type != NCOMP_DOUBLE) {
      convert_to<float>(tmp_fo, size_fo, 0, NCOMP_DOUBLE,
			((float *)fo->addr));
    }
  }

  return ier;
}

//////////////////////////////////////////////////////////////////////////////////////////
