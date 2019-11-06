#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ncomp.h>

double * read_data_rearrange(char * filename, int nTime, int nLat, int nLon) {
  FILE * fp;

  fp = fopen(filename, "r");
  if (fp != NULL) {
    double * data = (double *) malloc( (nTime*nLat*nLon) * sizeof(double));

    for (int k = 0; k < nTime; ++k) {
      for (int j = 0; j < nLat; ++j) {
        for (int i = 0; i < nLon; ++i) {
          char tmp[40];
          fscanf(fp,"%s", tmp);
          // changing nTime x nLat x nLon to
          // nLat x nLon x nTime
          int idx = j*(nLon*nTime) + i*(nTime) + k;
          data[idx] = atof(tmp);
        }
      }
    }
    fclose(fp);
    return data;
  } else {
    printf("Could not open (%s).\n", filename);
    return NULL;
  }
}

double * read_expected_ev(char * filename, int nTime, int nLat, int nLon) {
  FILE * fp;

  fp = fopen(filename, "r");
  if (fp != NULL) {
    double * data = (double *) malloc( (nTime*nLat*nLon) * sizeof(double));

    for (int k = 0; k < nTime; ++k) {
      for (int j = 0; j < nLat; ++j) {
        for (int i = 0; i < nLon; ++i) {
          char tmp[40];
          fscanf(fp,"%s", tmp);
          int idx = k*(nLat*nLon) + j*(nLon) + i;
          data[idx] = atof(tmp);
        }
      }
    }
    fclose(fp);
    return data;
  } else {
    printf("Could not open (%s).\n", filename);
    return NULL;
  }
}

double * read_expected_ev_ts(char * filename, int nEval, int nTime) {
  FILE * fp;

  fp = fopen(filename, "r");
  if (fp != NULL) {
    double * data = (double *) malloc( (nEval * nTime) * sizeof(double));

    for (int n = 0; n < nEval; ++n) {
      for (int t = 0; t < nTime; ++t) {
        char tmp[40];
        fscanf(fp,"%s", tmp);
        int idx = n*nTime + t;
        data[idx] = atof(tmp);
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
  printf("Testing eofunc_ts (01): ... \n");
  int nTime = 127;
  int nLat  =  91;
  int nLon  = 180;
  char * data_filename = "resources/SampleSST.csv";
  char * expected_ev_filename = "resources/SampleSST_ev_rearranged.txt";
  char * expected_ev_ts_filename = "resources/SampleSST_ev_rearranged_ts.txt";

  double * data = read_data_rearrange(data_filename, nTime, nLat, nLon);
  if (data == NULL) return 1;

  size_t dim_x[] = {nLat, nLon, nTime};
  ncomp_array* ncomp_x_in = ncomp_array_alloc((void*) data, NCOMP_DOUBLE, 3, dim_x);
  ncomp_x_in->has_missing = 0;
  ncomp_array* ev;
  ncomp_attributes* ev_attr = (ncomp_attributes*) malloc(sizeof(ncomp_attributes));

  printf("Calling eofunc: ...\n");
  int neval = 5;
  int ierr = eofunc(ncomp_x_in, neval, NULL, &ev, ev_attr);

  if (ierr != 0) {
    printf("ierr: %d\n", ierr);
    return ierr;
  }

  // Proceeding with checking the output
  printf("Testing eofunc output ...\n");
  if (ev->ndim != 3) {
    printf("problem with ncomp_x_out->ndim\n");
    return 2;
  }
  if (  (ev->shape[0] != neval) ||
        (ev->shape[1] != 91) ||
        (ev->shape[2] != 180) ) {
    printf("problem with ncomp_x_out->shape\n");
    printf("expected shape: [5, 91, 180]\n");
    printf(
      "  actual shape: [%zd, %zd, %zd]\n",
      ev->shape[0],
      ev->shape[1],
      ev->shape[2]
    );
    return 3;
  };

  double * expected_ev = read_expected_ev(expected_ev_filename, neval, nLat, nLon);
  for (int k = 0 ; k < neval; ++k) {
    for (int j = 0; j < nLat; ++j) {
      for (int i = 0; i < nLon; ++i) {
        int idx = k*(nLat*nLon) + j*(nLon) + i;
        if ( fabs( ((double*) ev->addr)[idx] - expected_ev[idx]) > 0.001 ) {
          printf("expected ev[%d,%d,%d]: %f\n", k, j, i, expected_ev[idx]);
          printf("  actual ev[%d,%d,%d]: %f\n", k, j, i, ((double*) ev->addr)[idx]);
          return 4;
        }
      }
    }
  }

  int expected_nAttribute = 5;
  if (ev_attr->nAttribute != expected_nAttribute) {
    printf("problem with attr->nAttribute\n");
    return 5;
  }

  for (int i = 0; i < expected_nAttribute; ++i) {
    ncomp_single_attribute * s_attr = ev_attr->attribute_array[i];

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
      if (  (s_attr->value->type != NCOMP_CHAR) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != (strlen("covariance") + 1)) ||
            strcmp((char *) s_attr->value->addr, "covariance")!=0 ) {
        printf("problem with matrix\n");
        return 12;
      }
    }

    if (strcmp("method", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_CHAR) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != (strlen("transpose") + 1)) ||
            strcmp((char *) s_attr->value->addr, "transpose")!=0 ) {
        printf("problem with method\n");
        return 13;
      }
    }
  }


  printf("Calling eofunc_ts: ...\n");
  ncomp_array* ev_ts;
  ncomp_attributes* ev_ts_attr = (ncomp_attributes*) malloc(sizeof(ncomp_attributes));
  ierr = eofunc_ts(ncomp_x_in, ev, NULL, &ev_ts, ev_ts_attr);

  if (ierr != 0) {
    printf("ierr: %d\n", ierr);
    return ierr;
  }

  printf("Testing eofunc_ts output: ...\n");
  if (ev_ts->ndim != 2) {
    printf("problem with ev_ts->ndim\n");
    return 2;
  }
  if (  (ev_ts->shape[0] != neval) ||
        (ev_ts->shape[1] != nTime) ) {
    printf("problem with ev_ts->shape\n");
    printf("expected shape: [5, 127]\n");
    printf(
      "  actual shape: [%zd, %zd]\n",
      ev_ts->shape[0],
      ev_ts->shape[1]
    );
    return 3;
  };

  double * expected_ev_ts = read_expected_ev_ts(expected_ev_ts_filename, neval, nTime);
  for (int n = 0 ; n < neval; ++n) {
    for (int t = 0; t < nTime; ++t) {
      int idx = n*nTime + t;

      if ( fabs( ((double*) ev_ts->addr)[idx] - expected_ev_ts[idx]) > 0.001 ) {
        printf("expected ev_ts[%d,%d]: %f\n", n, t, expected_ev_ts[idx]);
        printf("  actual ev_ts[%d,%d]: %f\n", n, t, ((double*) ev_ts->addr)[idx]);
        return 4;
      }
    }
  }

  expected_nAttribute = 2;
  if (ev_ts_attr->nAttribute != expected_nAttribute) {
    printf("problem with ev_ts_attr->nAttribute\n");
    return 5;
  }

  for (int i = 0; i < expected_nAttribute; ++i) {
    ncomp_single_attribute * s_attr = ev_ts_attr->attribute_array[i];

    if (strcmp("ts_mean", s_attr->name) == 0) {
      if (  (s_attr->value->type != 12) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != neval) ) {
        printf("problem with eval dimensions\n");
        return 8;
      }

      float expected_eval[] = {258.2679, 381.9555, 257.2193, 611.9719, -193.1238};
      for (int i = 0; i < neval; ++i) {
        if ( fabs( ((double*) s_attr->value->addr)[i] - expected_eval[i]) > 0.0001) {
          printf("expected ts_mean[%d]: %f\n",i, expected_eval[i]);
          printf("  actual ts_mean[%d]: %f\n",i, ((double*) s_attr->value->addr)[i]);
          return 9;
        }
      }
    }

    if (strcmp("matrix", s_attr->name) == 0) {
      if (  (s_attr->value->type != NCOMP_CHAR) ||
            (s_attr->value->ndim != 1) ||
            (s_attr->value->shape[0] != (strlen("covariance") + 1)) ||
            strcmp((char *) s_attr->value->addr, "covariance")!=0 ) {
        printf("problem with matrix\n");
        return 12;
      }
    }
  }

  free(data);
  free(expected_ev);
  free(expected_ev_ts);
  return 0;
}
