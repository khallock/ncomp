#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ncomp.h>



int main(void) {
  printf("Testing eofunc_n (04): ...\n");

  // Creating an empty options
  ncomp_attributes* options = ncomp_attributes_allocate(2);

  int jopt = 0;
  options->attribute_array[0] = create_ncomp_single_attribute_from_scalar((char *) "jopt", &jopt, NCOMP_INT);

  double pcrit = 32.0;
  options->attribute_array[1] = create_ncomp_single_attribute_from_scalar((char *) "pcrit", &pcrit, NCOMP_DOUBLE);

  printf("options->nAttribute: %d\n", options->nAttribute);

  // Creating sample data:
  int nx = 4;
  int ny = 4;
  int nt = 4;
  int nelem = nx*ny*nt;
  size_t dim_x[] = {nx, ny, nt};
  float x_in [] = {0,1,-99,-99,4,-99,6,-99,8,9,10,-99,12,-99,14,15,16,-99,18,-99,20,21,22,-99,24,25,26,27,28,-99,30,-99,32,33,34,35,36,-99,38,39,40,-99,42,-99,44,45,46,-99,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63};


  // creating rest of the eofunc arguments:
  ncomp_array* ncomp_x_in = ncomp_array_alloc((void*) x_in, NCOMP_FLOAT, 3, dim_x);
  ncomp_x_in->has_missing = 1;
  ncomp_x_in->msg.msg_float = -99.0;
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

  printf("attr_name: %s\n", attr->attribute_array[3]->name);
  printf("attr_value: %f\n", * ((double *) attr->attribute_array[3]->value->addr) );

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

  float expected_x_out[] = {0.3139, 0.1243, 0.1274, -99.0, 0.3139, 0.0318, 0.3139, -99.0, 0.3139, 0.2821, 0.3139, 0.0303, 0.3139, 0.3139, 0.3139, 0.3139};

  for (int i = 0; i < 16; ++i) {
    if (fabs( ((float*) ncomp_x_out->addr)[i] - expected_x_out[i]) > 0.0001) {
      printf("problem with ncomp_x_out->addr\n");
      printf("Expected: x_out[%d] = %f\n", i, expected_x_out[i]);
      printf("  Actual: x_out[%d] = %f\n", i, ((float*) ncomp_x_out->addr)[i]);
      return 3;
    }
  }

  int expected_nAttribute = 6;
  if (attr->nAttribute != expected_nAttribute) {
    printf("problem with attr->nAttribute\n");
    printf("Expected: %d\n", expected_nAttribute);
    printf("  Actual: %d\n", attr->nAttribute);
    return 4;
  }

  for (int i = 0; i<5; ++i) {
    ncomp_single_attribute * s_attr = attr->attribute_array[i];

    if (strcmp("eval_transpose", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_FLOAT) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            ( fabs( ((float*) s_attr->value->addr)[0] - 84.7542) > 0.0001) ) {
        printf("problem with eval_transpose\n");
        printf("Expected: %f\n", 2.9852);
        printf("  Actual: %f\n", ((float*) s_attr->value->addr)[0]);
        return 5;
      }
    }

    if (strcmp("eval", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_FLOAT) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            ( fabs( ((float*) s_attr->value->addr)[0] - 339.0166) > 0.0001) ) {
        printf("problem with eval\n");
        printf("Expected: %f\n", 14.9260);
        printf("  Actual: %f\n", ((float*) s_attr->value->addr)[0]);
        return 6;
      }
    }

    if (strcmp("pcvar", s_attr->name) == 0) {
      if (  (s_attr->value->type != 11) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            ( fabs( ((float*) s_attr->value->addr)[0] - 102.4951) > 0.0001) ) {
        printf("problem with pcvar\n");
        printf("Expected: %f\n", 98.7163);
        printf("Actual: %f\n", ((float*) s_attr->value->addr)[0]);
        return 7;
      }
    }

    if (strcmp("pcrit", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_FLOAT) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            ( fabs( ((float*) s_attr->value->addr)[0] - 32.0000) > 0.0001) ) {
        printf("problem with pcrit\n");
        return 6;
      }
    }

    if (strcmp("matrix", s_attr->name) == 0) {
      if (  (s_attr->value->type != 0) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            strcmp((char *) s_attr->value->addr, "covariance")!=0 ) {
        printf("problem with matrix\n");
        return 8;
      }
    }

    if (strcmp("method", s_attr->name) == 0) {
      if (  (s_attr->value->type != 0) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            strcmp((char *) s_attr->value->addr, "transpose")!=0 ) {
        printf("problem with method\n");
        return 9;
      }
    }
  }



  printf("All done!!!\n");
  return 0;
}
