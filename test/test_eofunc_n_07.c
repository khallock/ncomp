#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ncomp.h>

double * readfile(char * filename, int nTime, int nLat, int nLon) {
  FILE * fp;

  fp = fopen(filename, "r");
  if (fp != NULL) {
    double * data = (double *) malloc( (nTime*nLat*nLon +2 +1) * sizeof(double));

    for (int k = 0; k < nTime; ++k) {
      for (int j = 0; j < nLat; ++j) {
        for (int i = 0; i < nLon; ++i) {
          char tmp[20];
          fscanf(fp,"%s", tmp);
          int idx = k*(nLat*nLon) + j*(nLon) + i;
          data[idx] = atof(tmp);
          // printf("%f ", data[idx]);
        }
        // printf("\n");
      }
    }
    fclose(fp);
    return data;
  } else {
    printf("Could not open (%s).\n", filename);
    return NULL;
  }
}

int main(void) {
  printf("Testing eofunc_n (07): ...\n");
  int nTime = 127;
  int nLat  =  91;
  int nLon  = 180;
  char * data_filename = "resources/SampleSST.csv";
  char * expected_output_filename = "resources/SampleSST_ev.txt";

  double * data = readfile(data_filename, nTime, nLat, nLon);
  if (data == NULL) return 1;

  size_t dim_x[] = {nTime, nLat, nLon};
  ncomp_array* ncomp_x_in = ncomp_array_alloc((void*) data, NCOMP_DOUBLE, 3, dim_x);
  ncomp_x_in->has_missing = 0;
  ncomp_array* ncomp_x_out = (ncomp_array*) malloc(sizeof(ncomp_array));
  ncomp_attributes* attr = (ncomp_attributes*) malloc(sizeof(ncomp_attributes));

  printf("Calling eofunc_n: ...\n");
  int t_dim = 0;
  int neval = 5;
  int ierr = eofunc_n(ncomp_x_in, neval, t_dim, NULL, ncomp_x_out, attr);

  if (ierr != 0) {
    printf("ierr: %d", ierr);
    return ierr;
  }

  print_ncomp_attributes(attr);

  // Proceeding with checking the output
  if (ncomp_x_out->ndim != 3) {
    printf("problem with ncomp_x_out->ndim\n");
    return 2;
  }
  if (  (ncomp_x_out->shape[0] != neval) ||
        (ncomp_x_out->shape[1] != 91) ||
        (ncomp_x_out->shape[2] != 180) ) {
    printf("problem with ncomp_x_out->shape\n");
    printf("expected shape: [5, 91, 180]\n");
    printf(
      "  actual shape: [%zd, %zd, %zd]\n",
      ncomp_x_out->shape[0],
      ncomp_x_out->shape[1],
      ncomp_x_out->shape[2]
    );
    return 3;
  };

  double * expected_data = readfile(expected_output_filename, neval, nLat, nLon);
  for (int k = 0 ; k < neval; ++k) {
    for (int j = 0; j < nLat; ++j) {
      for (int i = 0; i < nLon; ++i) {
        int idx = k*(nLat*nLon) + j*(nLon) + i;
        if ( fabs( ((double*) ncomp_x_out->addr)[idx] - expected_data[idx]) > 0.001 ) {
          printf("expected ev[%d,%d,%d]: %f\n", k, j, i, expected_data[idx]);
          printf("  actual ev[%d,%d,%d]: %f\n", k, j, i, ((double*) ncomp_x_out->addr)[idx]);
          return 4;
        }
      }
    }
  }

  int expected_nAttribute = 5;
  if (attr->nAttribute != expected_nAttribute) {
    printf("problem with attr->nAttribute\n");
    return 5;
  }

  for (int i = 0; i < expected_nAttribute; ++i) {
    ncomp_single_attribute * s_attr = attr->attribute_array[i];

    if (strcmp("pcvar", s_attr->name) == 0) {
      if (  (s_attr->value->type != 11) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != neval) ) {
        printf("problem with pcvar dimensions\n");
        return 6;
      }

      float expected_pcvar[] = {86.407, 3.540, 2.655, 1.047, 0.776};
      for (int i = 0; i < 5; ++i) {
        if ( fabs( ((float*) s_attr->value->addr)[i] - expected_pcvar[i]) > 0.001) {
          printf("expected pcvar[%d]: %f\n",i, expected_pcvar[i]);
          printf("  actual pcvar[%d]: %f\n",i, ((float*) s_attr->value->addr)[i]);
          return 7;
        }
      }
    }

    if (strcmp("eval", s_attr->name) == 0) {
      if (  (s_attr->value->type != 12) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != neval) ) {
        printf("problem with eval dimensions\n");
        return 8;
      }

      float expected_eval[] = {61037.9258, 2500.7486, 1875.6576, 739.3867, 548.1788};
      for (int i = 0; i < 5; ++i) {
        if ( fabs( ((double*) s_attr->value->addr)[i] - expected_eval[i]) > 0.0001) {
          printf("expected eval[%d]: %f\n",i, expected_eval[i]);
          printf("  actual eval[%d]: %f\n",i, ((double*) s_attr->value->addr)[i]);
          return 9;
        }
      }
    }

    if (strcmp("eval_transpose", s_attr->name) == 0) {
      if (  (s_attr->value->type != 12) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != neval) ) {
        printf("problem with eval_transpose dimensions\n");
        return 10;
      }

      float expected_eval_transpose[] = {570.4479, 23.3715, 17.5295, 6.9102, 5.1232};
      for (int i = 0; i < 5; ++i) {
        if ( fabs( ((double*) s_attr->value->addr)[i] - expected_eval_transpose[i]) > 0.0001) {
          printf("expected eval_transpose[%d]: %f\n",i, expected_eval_transpose[i]);
          printf("  actual eval_transpose[%d]: %f\n",i, ((double*) s_attr->value->addr)[i]);
          return 11;
        }
      }
    }

    if (strcmp("matrix", s_attr->name) == 0) {
      if (  (s_attr->value->type != 0) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            strcmp((char *) s_attr->value->addr, "covariance")!=0 ) {
        printf("problem with matrix\n");
        return 12;
      }
    }

    if (strcmp("method", s_attr->name) == 0) {
      if (  (s_attr->value->type != 0) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != 1) ||
            strcmp((char *) s_attr->value->addr, "transpose")!=0 ) {
        printf("problem with method\n");
        return 13;
      }
    }
  }


  free(data);
  free(expected_data);
  return 0;
}
