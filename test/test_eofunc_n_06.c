#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ncomp.h>

int main(void) {
  printf("Testing eofunc_n (06): ...\n");

  // Creating an empty options
  ncomp_attributes* options = ncomp_attributes_allocate(0);

  printf("options->nAttribute: %d\n", options->nAttribute);

  // Creating sample data:
  int nx = 4;
  int ny = 4;
  int nt = 4;
  int nelem = nx*ny*nt;
  size_t dim_x[] = {nx, ny, nt};
  double * x_in = (double *) malloc(sizeof(double) * nelem);

  // filling in with some data
  for (int i = 0; i < (nelem); ++i) {
    x_in[i] = i;
  }

  // creating rest of the eofunc arguments:
  ncomp_array* ncomp_x_in = ncomp_array_alloc((void*) x_in, NCOMP_DOUBLE, 3, dim_x);
  ncomp_x_in->has_missing = 0;
  ncomp_x_in->msg.msg_double = -999.0;
  ncomp_array* ncomp_x_out = (ncomp_array*) malloc(sizeof(ncomp_array));
  ncomp_attributes* attr = (ncomp_attributes*) malloc(sizeof(ncomp_attributes));
  int neval = 1;

  printf("Calling eofunc_n: ...\n");
  int t_dim = 1;
  int ierr = eofunc_n(ncomp_x_in, neval, t_dim, options, ncomp_x_out, attr);

  if (ierr != 0) {
    printf("ierr: %d", ierr);
    return ierr;
  }

  print_ncomp_array("x_out", ncomp_x_out);
  print_ncomp_attributes(attr);

  if (ncomp_x_out->ndim != 3) {
    printf("problem with ncomp_x_out->ndim\n");
    return 1;
  }
  if (  (ncomp_x_out->shape[0] != 1) ||
        (ncomp_x_out->shape[1] != 4) ||
        (ncomp_x_out->shape[2] != 4) ) {
    printf("problem with ncomp_x_out->shape\n");
    return 2;
  };

  for (int i = 0; i < 16; ++i) {
    if (fabs( ((double*) ncomp_x_out->addr)[i] - 0.25) > 0.001) {
      printf("problem with ncomp_x_out->addr\n");
      return 3;
    }
  }

  int expected_nAttribute = 5;
  if (attr->nAttribute != expected_nAttribute) {
    printf("problem with attr->nAttribute\n");
    return 4;
  }

  for (int i = 0; i < expected_nAttribute; ++i) {
    ncomp_single_attribute * s_attr = attr->attribute_array[i];

    if (strcmp("eval_transpose", s_attr->name) == 0) {
      if (  (s_attr->value->type != 12) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            ( fabs( ((double*) s_attr->value->addr)[0] - 85.3333) > 0.0001) ) {
        printf("problem with eval_transpose\n");
        return 5;
      }
    }

    if (strcmp("eval", s_attr->name) == 0) {
      if (  (s_attr->value->type != 12) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            ( fabs( ((double*) s_attr->value->addr)[0] - 426.6667) > 0.0001) ) {
        printf("problem with eval\n");
        return 6;
      }
    }

    if (strcmp("pcvar", s_attr->name) == 0) {
      if (  (s_attr->value->type != 11) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            ( fabs( ((float*) s_attr->value->addr)[0] - 100.00) > 0.0001) ) {
        printf("problem with pcvar\n");
        return 7;
      }
    }

    if (strcmp("matrix", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_CHAR) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != (strlen("covariance") + 1)) ||
            strcmp((char *) s_attr->value->addr, "covariance")!=0 ) {
        printf("problem with matrix\n");
        return 8;
      }
    }

    if (strcmp("method", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_CHAR) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != (strlen("transpose") + 1)) ||
            strcmp((char *) s_attr->value->addr, "transpose")!=0 ) {
        printf("problem with method\n");
        return 9;
      }
    }
  }

  printf("All done!!!\n");
  return 0;
}
