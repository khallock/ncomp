#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ncomp.h>

int main(void) {
  printf("Testing moc_globe_atl: ...\n");

  // Creating sample data:
  int nyaux = 3;            // nyaux = lat_aux_grid->shape[0]
  int kdep = 4;             // kdep  = a_wvel/a_bolus/a_submeso->shape[0]
  int nlat = 5;             // nlat  = a_wvel/a_bolus/a_submeso->shape[1] AND tlat->shape[0] AND rmlak->shape[1]
  int mlon = 6;             // mlon  = a_wvel/a_bolus/a_submeso->shape[2] AND tlat->shape[1] AND rmlak->shape[2]

  int kdepnyaux2   = 2 * kdep * nyaux;
  int nlatmlon     = nlat * mlon;
  int kdepnlatmlon = kdep * nlatmlon;

  size_t dim_lat_aux_grid[] = {nyaux};
  double * lat_aux_grid = (double *) malloc(sizeof(double) * nyaux);

  size_t dim_tlat[] = {nlat, mlon};
  double * tlat = (double *) malloc(sizeof(double) * nlatmlon);

  size_t dim_rmlak[] = {2, nlat, mlon};
  int * rmlak = (int *) malloc(sizeof(int) * 2 * nlatmlon);

  size_t dim_a_wvel[] = {kdep, nlat, mlon};
  double * a_wvel = (double *) malloc(sizeof(double) * kdepnlatmlon);
  double * a_bolus = (double *) malloc(sizeof(double) * kdepnlatmlon);
  double * a_submeso = (double *) malloc(sizeof(double) * kdepnlatmlon);

  // filling in with some data
  for (int i = 0; i < nyaux; i++) {
    lat_aux_grid[i] = (i+1) * 5;
  }

  for (int i = 0; i < nlatmlon; i++) {
    tlat[i] = (i+1);
  }

  for (int i = 0; i < 2 * nlatmlon; i++) {
    rmlak[i] = 1;
  }

  for (int i = 0; i < kdepnlatmlon; i++) {
    a_wvel[i] = (i+1);
    a_bolus[i] = (i+1);
    a_submeso[i] = (i+1);
  }

  // Generating ncmop_array data from created data
  ncomp_array* ncomp_a_wvel = ncomp_array_alloc((void*) a_wvel, NCOMP_DOUBLE, 3, dim_a_wvel);
  ncomp_a_wvel->has_missing = 1;
  ncomp_a_wvel->msg.msg_double = -999.0;
  ncomp_array* ncomp_a_submeso = ncomp_array_alloc((void*) a_submeso, NCOMP_DOUBLE, 3, dim_a_wvel);
  ncomp_array* ncomp_a_bolus = ncomp_array_alloc((void*) a_bolus, NCOMP_DOUBLE, 3, dim_a_wvel);
  ncomp_array* ncomp_lat_aux_grid = ncomp_array_alloc((void*) lat_aux_grid, NCOMP_DOUBLE, 1, dim_lat_aux_grid);
  ncomp_array* ncomp_tlat = ncomp_array_alloc((void*) tlat, NCOMP_DOUBLE, 2, dim_tlat);
  ncomp_array* ncomp_rmlak = ncomp_array_alloc((void*) rmlak, NCOMP_INT, 3, dim_rmlak);

  // Allocating ,e,ory for output data
  ncomp_array* ncomp_arr_out = NULL;

  // Calling moc_globe_atl function
  printf("Calling moc_globe_atl: ...\n");

  int ierr = moc_globe_atl(ncomp_lat_aux_grid, ncomp_a_wvel, ncomp_a_bolus, ncomp_a_submeso, ncomp_tlat, ncomp_rmlak, &ncomp_arr_out);

  if (ierr != 0) {
    printf("ierr: %d", ierr);
    return ierr;
  }

  print_ncomp_array("ncomp_arr_out", ncomp_arr_out);

  if (ncomp_arr_out->ndim != 4) {
    printf("problem with ncomp_arr_out->ndim\n");
    return 1;
  }
  if (  (ncomp_arr_out->shape[0] != 3) ||
        (ncomp_arr_out->shape[1] != 2) ||
        (ncomp_arr_out->shape[2] != kdep) ||
        (ncomp_arr_out->shape[3] != nyaux) ) {
    printf("problem with ncomp_arr_out->shape\n");
    return 2;
  };

  /*
  for (int i = 0; i < 16; ++i) {
    if (fabs( ((double*) ncomp_arr_out->addr)[i] - 0.25) > 0.001) {
      printf("problem with ncomp_x_out->addr\n");
      return 3;
    }
  }
  */

  printf("All done!!!\n");
  return 0;
}
