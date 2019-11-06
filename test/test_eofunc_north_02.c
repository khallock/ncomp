#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ncomp.h>

int main(void) {
  printf("Testing eofunc_north (01): ...\n");

  size_t dim[] = {4};
  double eval[] = {974.8881, 490.3863, 190.7456, 172.2483};

  ncomp_array* ncomp_eval_in = ncomp_array_alloc((void*) eval, NCOMP_DOUBLE, 1, dim);
  ncomp_eval_in->has_missing = 1;
  ncomp_eval_in->msg.msg_double = -99.0;

  ncomp_array* ncomp_sig_out;
  ncomp_attributes* attr = (ncomp_attributes*) malloc(sizeof(ncomp_attributes));

  int N = 80;
  int prinfo = 1;

  printf("Calling eofunc_north ... \n");
  int ierr = eofunc_north(ncomp_eval_in, N, prinfo, &ncomp_sig_out, attr);

  print_ncomp_array("ncomp_sig_out", ncomp_sig_out);
  print_ncomp_attributes(attr);

  if (ncomp_sig_out->ndim != 1) {
    printf("problem with ncomp_sig_out->ndim\n");
    return 1;
  }
  if (  (ncomp_sig_out->shape[0] != 4) ) {
    printf("problem with ncomp_sig_out->shape\n");
    return 2;
  };

  int expected_sig_out[] = {1, 1, 0, 0};

  for (int i = 0; i < 4; ++i) {
    if (((int*) ncomp_sig_out->addr)[i] != expected_sig_out[i]) {
      printf("problem with ncomp_sig_out->addr\n");
      printf("Expected: ncomp_sig_out[%d] = %d\n", i, expected_sig_out[i]);
      printf("  Actual: ncomp_sig_out[%d] = %d\n", i, ((int*) ncomp_sig_out->addr)[i]);
      return 3;
    }
  }

  int expected_nAttribute = 2;
  if (attr->nAttribute != expected_nAttribute) {
    printf("problem with attr->nAttribute\n");
    printf("Expected: %d\n", expected_nAttribute);
    printf("  Actual: %d\n", attr->nAttribute);
    return 4;
  }

  for (int i = 0; i<expected_nAttribute; ++i) {
    ncomp_single_attribute * s_attr = attr->attribute_array[i];

    if (strcmp("long_name", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_CHAR) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != (strlen("EOF separation") + 1)) ||
            strcmp((char *) s_attr->value->addr, "EOF separation")!=0 ) {
        printf("problem with long_name\n");
        return 5;
      }
    }

    if (strcmp("N", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_INT) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            ( ((int*) s_attr->value->addr)[0] != N ) ) {
        printf("problem with N\n");
        printf("Expected: %d\n", N);
        printf("Actual: %d\n", ((int*) s_attr->value->addr)[0]);
        return 7;
      }
    }

  }
  return ierr;
}
