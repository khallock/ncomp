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
                              const ncomp_array *tlat, const ncomp_array *rmlak )
{
/*
 * Various
 */
  long nyaux, kdep, nlat, mlon, nlatmlon, kdepnlatmlon, kdepnyaux2;
  long i, size_output;
  int nrx, inlat, imlon, ikdep, inyaux;
  int ret;

 // Sanity Checking
 if (lat_aux_grid->ndim < 1) {
   #if DEBUG
     std::cerr<<"moc_globe_atl: Empty array!!!"<<std::endl;
   #endif
   return 1;
 }

 nyaux = lat_aux_grid->shape[0];

  // Sanity Checking
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

/*
 * Test dimension sizes.
 */
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
 * Coerce missing value to double if necessary.
 */
  coerce_missing(type_a_wvel,has_missing_a_wvel,&missing_a_wvel,
                 &missing_dbl_a_wvel,NULL);

/*
 * Coerce input arrays to double if necassary.
 */
  tmp_lat_aux_grid = coerce_input_double(lat_aux_grid,type_lat_aux_grid,nyaux,
					 0,NULL,NULL);
  if(tmp_lat_aux_grid == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"moc_globe_atl: Unable to allocate memory for coercing lat_aux_grid to double");
    return(NhlFATAL);
  }

  tmp_a_wvel = coerce_input_double(a_wvel,type_a_wvel,kdepnlatmlon,
                                   has_missing_a_wvel,&missing_a_wvel,
                                   &missing_dbl_a_wvel);

  if(tmp_a_wvel == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"moc_globe_atl: Unable to allocate memory for coercing a_wvel to double");
    return(NhlFATAL);
  }

  tmp_a_bolus = coerce_input_double(a_bolus,type_a_bolus,kdepnlatmlon,0,
                                    NULL,NULL);
  if(tmp_a_bolus == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"moc_globe_atl: Unable to allocate memory for coercing a_bolus to double");
    return(NhlFATAL);
  }

  tmp_a_submeso = coerce_input_double(a_submeso,type_a_submeso,kdepnlatmlon,
                                      0,NULL,NULL);
  if(tmp_a_submeso == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"moc_globe_atl: Unable to allocate memory for coercing a_submeso to double");
    return(NhlFATAL);
  }

  tmp_tlat = coerce_input_double(tlat,type_tlat,nlatmlon,0,NULL,NULL);
  if(tmp_tlat == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"moc_globe_atl: Unable to allocate memory for coercing tlat to double");
    return(NhlFATAL);
  }


/*
 * The output will be float unless any of the first four arguments
 * are double.
 */
  if(type_a_wvel == NCL_double || type_a_submeso == NCL_double ||
     type_a_bolus == NCL_double) {
    type_tmp = NCL_double;
  }
  else {
    type_tmp = NCL_float;
  }

/*
 * Allocate space for output array.
 */
  size_output = 3 * kdepnyaux2;    /* 3 x 2 x kdep x nyaux */

  if(type_tmp != NCL_double) {
    tmp   = (void *)calloc(size_output, sizeof(float));
    dtmp1 = (double *)calloc(kdepnyaux2, sizeof(double));
    dtmp2 = (double *)calloc(kdepnyaux2, sizeof(double));
    dtmp3 = (double *)calloc(kdepnyaux2, sizeof(double));
    if(tmp == NULL || dtmp1 == NULL || dtmp2 == NULL ||
       dtmp3 == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"moc_globe_atl: Unable to allocate memory for temporary output arrays");
      return(NhlFATAL);
    }
  }
  else {
    tmp = (void *)calloc(size_output, sizeof(double));
    if(tmp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"moc_globe_atl: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    dtmp1 = &((double*)tmp)[0];
    dtmp2 = &((double*)tmp)[kdepnyaux2];
    dtmp3 = &((double*)tmp)[2*kdepnyaux2];
  }

/*
 * Allocate space for output dimension sizes and set them.
 */
  ndims_tmp  = 4;
  dsizes_tmp = (ng_size_t*)calloc(ndims_tmp,sizeof(ng_size_t));
  if( dsizes_tmp == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"moc_globe_atl: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  dsizes_tmp[0] = 3;
  dsizes_tmp[1] = 2;
  dsizes_tmp[2] = kdep;
  dsizes_tmp[3] = nyaux;

/*
 * Call the Fortran routine.
 */
  nrx = 2;

  NGCALLF(mocloops,MOCLOOPS)(&inyaux, &imlon, &inlat, &ikdep, &nrx, tmp_tlat,
                             tmp_lat_aux_grid, rmlak, tmp_a_wvel, tmp_a_bolus,
                             tmp_a_submeso, &missing_dbl_a_wvel.doubleval,
                             dtmp1, dtmp2, dtmp3);

  if(type_tmp != NCL_double) {
    coerce_output_float_only(tmp,dtmp1,kdepnyaux2,0);
    coerce_output_float_only(tmp,dtmp2,kdepnyaux2,kdepnyaux2);
    coerce_output_float_only(tmp,dtmp3,kdepnyaux2,2*kdepnyaux2);
  }

/*
 * Free unneeded memory.
 */
  if(type_lat_aux_grid != NCL_double) NclFree(tmp_lat_aux_grid);
  if(type_a_wvel    != NCL_double) NclFree(tmp_a_wvel);
  if(type_a_bolus   != NCL_double) NclFree(tmp_a_bolus);
  if(type_a_submeso != NCL_double) NclFree(tmp_a_submeso);
  if(type_tlat      != NCL_double) NclFree(tmp_tlat);
  if(type_tmp       != NCL_double) {
    NclFree(dtmp1);
    NclFree(dtmp2);
    NclFree(dtmp3);
  }

/*
 * Return value back to NCL script.
 */
  ret = NclReturnValue(tmp,ndims_tmp,dsizes_tmp,NULL,type_tmp,0);
  NclFree(dsizes_tmp);
  return(ret);
}
