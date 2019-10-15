#include "ncomp/types.h"
#include "ncomp/util.h"
#include "ncomp_internal/util.hpp"
#include <iostream>

extern "C" void mocloops_(int *, int *, int *, int *, int *,
                                       double *, double *, int *, double *,
                                       double *, double *, double *, double *,
                                       double *, double *);

/*
* NOTE: Adapted from moc_globe_atl_W() in ncl/ni/src/lib/nfp/mocloopsW.c of original NCL code.
*
* Facilitates calculating the meridional overturning circulation for the globe and Atlantic.
*
*/

extern "C" int moc_globe_atl( const ncomp_array *lat_aux_grid, const ncomp_array *a_wvel,
                              const ncomp_array *a_bolus, const ncomp_array *a_submeso,
                              const ncomp_array *tlat, const ncomp_array *rmlak,
                              ncomp_array * tmp_out ){

/*
 * Various Variables
 */
  long nyaux, kdep, nlat, mlon, nlatmlon, kdepnlatmlon, kdepnyaux2;
  long i, size_output;
  int nrx, inlat, imlon, ikdep, inyaux;
  int ret;


 /*
  * Output Variables
  */
  void *tmp;
  double *dtmp1, *dtmp2, *dtmp3;
  std::unique_ptr<size_t[]> dsizes_tmp;
  int ndims_tmp;
  int type_tmp;


  /*
   * Sanity Checks
   */

  /* Check for lat_aux_grid dimension */
  if (lat_aux_grid->ndim < 1) {
   #if DEBUG
     std::cerr<<"moc_globe_atl: Empty array!!!"<<std::endl;
   #endif
   return 1;
  }

  nyaux = lat_aux_grid->shape[0];

  /* Check for a_wvel dimension */
  if (a_wvel->ndim < 3) {
    #if DEBUG
      std::cerr<<"moc_globe_atl: The input array must be at least three-dimensional"<<std::endl;
    #endif
    return 1;
  }

  kdep  = a_wvel->shape[0];
  nlat  = a_wvel->shape[1];
  mlon  = a_wvel->shape[2];
  kdepnyaux2   = 2 * kdep * nyaux;
  nlatmlon     = nlat * mlon;
  kdepnlatmlon = kdep * nlatmlon;

  /* Test dimension sizes. */
  if((mlon > INT_MAX) || (nlat > INT_MAX) ||
     (kdep > INT_MAX) || (nyaux > INT_MAX)) {
    #if DEBUG
      std::cerr<<"moc_globe_atl: one or more input dimension sizes are greater than INT_MAX"<<std::endl;
    #endif
    return 1;
  }

  imlon = (int) mlon;
  inlat = (int) nlat;
  ikdep = (int) kdep;
  inyaux = (int) nyaux;

  for(i = 0; i <= 2; i++) {
    if(a_bolus->shape[i] != a_wvel->shape[i] ||
       a_bolus->shape[i] != a_submeso->shape[i]) {
       #if DEBUG
         std::cerr<<"moc_globe_atl: a_wvel, a_submeso, and a_bolus must have the same dimensionality"<<std::endl;
       #endif
       return 1;
    }
  }

  if(tlat->shape[0] != nlat || tlat->shape[1] != mlon) {
    #if DEBUG
      std::cerr<<"moc_globe_atl: The dimensions of tlat must be nlat x mlon"<<std::endl;
    #endif
    return 1;
  }

  if(rmlak->shape[0] != 2 ||
     rmlak->shape[1] != nlat || rmlak->shape[2] != mlon) {
      #if DEBUG
        std::cerr<<"moc_globe_atl: The dimensions of rmlak must be 2 x nlat x mlon"<<std::endl;
      #endif
      return 1;
  }


 /*
  * Handle missing values
  */
  double missing_d_a_wvel;
  float missing_f_a_wvel;
  coerce_missing(a_wvel->type, a_wvel->has_missing, (ncomp_missing *)&(a_wvel->msg),
                &missing_d_a_wvel, &missing_f_a_wvel);


/*
 * Convert input arrays to double if necassary.
 */
 size_t lat_aux_grid_numel = prod(lat_aux_grid->shape, lat_aux_grid->ndim);
 double *tmp_lat_aux_grid = convert_to_with_copy_avoiding<double>((void *)lat_aux_grid->addr,
                              lat_aux_grid_numel, 0, lat_aux_grid->type, NCOMP_DOUBLE);

size_t a_wvel_numel = prod(a_wvel->shape, a_wvel->ndim);
double *tmp_a_wvel = convert_to_with_copy_avoiding<double>((void *)a_wvel->addr,
                             a_wvel_numel, 0, a_wvel->type, NCOMP_DOUBLE);

size_t a_bolus_numel = prod(a_bolus->shape, a_bolus->ndim);
double *tmp_a_bolus = convert_to_with_copy_avoiding<double>((void *)a_bolus->addr,
                            a_bolus_numel, 0, a_bolus->type, NCOMP_DOUBLE);

size_t a_submeso_numel = prod(a_submeso->shape, a_submeso->ndim);
double *tmp_a_submeso = convert_to_with_copy_avoiding<double>((void *)a_submeso->addr,
                           a_submeso_numel, 0, a_submeso->type, NCOMP_DOUBLE);

size_t tlat_numel = prod(tlat->shape, tlat->ndim);
double *tmp_tlat = convert_to_with_copy_avoiding<double>((void *)tlat->addr,
                          tlat_numel, 0, tlat->type, NCOMP_DOUBLE);


/*
 * The output will be float unless any of the first four arguments are double.
 */
  if(a_wvel->type == NCOMP_DOUBLE || a_submeso->type == NCOMP_DOUBLE ||
     a_bolus->type == NCOMP_DOUBLE) {
    type_tmp = NCOMP_DOUBLE;
  }
  else {
    type_tmp = NCOMP_FLOAT;
  }


/*
 * Allocate space for output array.
 */
 size_output = 3 * kdepnyaux2;    /* 3 x 2 x kdep x nyaux */

  /* create temporary output arrays */
  if (type_tmp != NCOMP_DOUBLE) {
    tmp   = (void *)(new float[size_output]);

    dtmp1 = new double[kdepnyaux2];
    dtmp2 = new double[kdepnyaux2];
    dtmp3 = new double[kdepnyaux2];
  }
  else {
    tmp   = (void *)(new double[size_output]);

    dtmp1 = &((double*)tmp)[0];
    dtmp2 = &((double*)tmp)[kdepnyaux2];
    dtmp3 = &((double*)tmp)[2 * kdepnyaux2];
  }

/*
 * Allocate space for output dimension sizes and set them.
 */
  ndims_tmp  = 4;
  dsizes_tmp.reset(new size_t[ndims_tmp]);
  dsizes_tmp[0] = 3;
  dsizes_tmp[1] = 2;
  dsizes_tmp[2] = kdep;
  dsizes_tmp[3] = nyaux;


/*
 * Call the Fortran routine.
 */
  nrx = 2;

  mocloops_(&inyaux, &imlon, &inlat, &ikdep, &nrx, tmp_tlat, tmp_lat_aux_grid,
            rmlak, tmp_a_wvel, tmp_a_bolus, tmp_a_submeso, &missing_d_a_wvel,
            dtmp1, dtmp2, dtmp3);


  /*
   * Return variables.
   */

   if(type_tmp != NCOMP_DOUBLE){

     /* Convert tmp array to floats first*/
     convert_to<float>(dtmp1, kdepnyaux2, 0, NCOMP_DOUBLE,
                         (float *)tmp + 0);

     convert_to<float>(dtmp2, kdepnyaux2, 0, NCOMP_DOUBLE,
                         (float *)tmp + kdepnyaux2);

     convert_to<float>(dtmp3, kdepnyaux2, 0, NCOMP_DOUBLE,
                         (float *)tmp + 2 * kdepnyaux2);

      /* Populate output ncomp_array from tmp array */
      *tmp_out = *ncomp_array_alloc((float *)tmp, NCOMP_FLOAT, ndims_tmp, dsizes_tmp.get());
   }
   else
      /* Populate output ncomp_array from tmp array */
      *tmp_out = *ncomp_array_alloc((double *)tmp, NCOMP_DOUBLE, ndims_tmp, dsizes_tmp.get());

  /* TO-DO: Check if has_missing and msg.msg_double needed
   tmp_out->has_missing = a_wvel->has_missing;
   tmp_out->msg.msg_double = missing_d_a_wvel;
   */


/*
 * Free unneeded memory.
 */
   if(lat_aux_grid->type != NCOMP_DOUBLE)
      delete [] tmp_lat_aux_grid;

   if(a_wvel->type != NCOMP_DOUBLE)
      delete [] tmp_a_wvel;

   if(a_bolus->type != NCOMP_DOUBLE)
      delete [] tmp_a_bolus;

   if(a_submeso->type != NCOMP_DOUBLE)
      delete [] tmp_a_submeso;

   if(tlat->type != NCOMP_DOUBLE)
      delete [] tmp_tlat;

   if(type_tmp != NCOMP_DOUBLE) {
      delete [] dtmp1;
      delete [] dtmp2;
      delete [] dtmp3;
   }


/*
 * TO-DO: Revisit to return a meaningful value, if any. Intuition for returnin 0
 * for now is that if any error did not occur thorughout the code, then there
 * should not be an erro at this point.
 */
  return 0;
}
