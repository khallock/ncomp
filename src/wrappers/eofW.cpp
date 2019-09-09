#include "ncomp/types.h"
#include "ncomp/util.h"
#include "ncomp_internal/util.hpp"
#include <iostream>
#include <memory>

extern "C" void crveoft_( double *dx_strip, double *dx_strip_t,
                         int *nrow, int *ncol, int *nrobs,
                         int *mcsta, double *xmsg, int *neval,
                         double *eval, double *evec,
                         double *pcvar, double *trace,
                         double *xdvar, double *xave,
                         int *jopt, int *ier);

extern "C" void ddrveof_(double *,int *,int *,int *,int *,
                         double *,int *,double *, double *,
                         float*,double *,int *,int *,double*,
                         long long int *, double *,int *,
                         double *,int *,int *,int *,int *,int *);

extern "C" void deof11_( double *, int *, int *, int *, int *,
                         double  *, double *, double *, double *,
                         double *);

extern "C" void dstat2_( double *, int *, double *, double *,
                         double *, double *, int *, int *);

extern "C" void xrveoft_( double *dx_strip, double *dx_strip_t,
                          int *nrow, int *ncol, int *nrobs,
                          int *mcsta, double *xmsg, int *neval,
                          double *eval, double *evec,
                          double *pcvar, double *trace,
                          double *xdvar, double *xave,
                          int *jopt, int *ier);

/*
* NOTE: adapted from eofunc_w() in ncl/ni/src/lib/nfp/eofW.c of original NCL code
* This routine calls three different EOF routines.
*
* The first is the original eofcov routine, which can be
* extremely slow if  nrow < mcsta.
*
* The second routine is one Dennis wrote in 2004/2005 to speed
* up the case where nrow < mcsta. This routine had a problem with
* one particular case with a French grid. Dennis spent quite a
* bit of time proving that this routine works with several
* textbook examples, but he's not certain why it is having problems
* with this one particular grid.
*
* The third routine was taken from SCRIPPS and modified by Dennis
* to handle missing values.
*
* If use_old_transpose = use_new_transpose = False, then the old
* routine is used.  If use_old_transpose = True, then Dennis' transpose
* routine is used. If use_new_transpose = True, then SCRIPPS transpose
* routine is used. Note: use_old_tranpose should only be used for
* debugging purposes. It is not intended to be used by the user.
*
*/

typedef struct { // options used for eofunc and their default values
  int jopt {0};
  double pcrit {50.0};
  bool return_pcrit {false};
  bool return_eval {false};
  bool return_trace {false};
  bool anomalies {false};
  bool use_new_transpose {false};
  bool use_old_transpose {false};
  bool tr_setbyuser {false};
  bool debug {false};
} eofunc_options;

eofunc_options* extract_eofunc_options(const attributes & options_in) {
  eofunc_options* options_out = new eofunc_options;

  options_out->jopt = *(int*) getAttributeOrDefault(options_in, "jopt", &((*options_out).jopt));
  if ((options_out->jopt != 0) && (options_out->jopt != 1)) {
    options_out->jopt = 0; // jopt must be either 0 or 1
  }

  int tmpPos = -1;
  if (hasAttribute(options_in, "pcrit", tmpPos)==1) {
    options_out->pcrit = *(double*) options_in.attribute_array[tmpPos].value->addr;
    options_out->return_pcrit = true;
    if ((options_out->pcrit < 0.0) && (options_out->pcrit > 100.0)) {
      options_out->pcrit = 50.0; // pcrit must be between 0.0 and 100.0; default value is 50
    }
  }

  options_out->return_eval = *(bool*) getAttributeOrDefault(options_in, "return_eval", &((*options_out).return_eval));

  options_out->return_trace = *(bool*) getAttributeOrDefault(options_in, "return_trace", &((*options_out).return_trace));

  options_out->anomalies = *(bool*) getAttributeOrDefault(options_in, "anomalies", &((*options_out).anomalies));

  if (hasAttribute(options_in, "transpose", tmpPos)==1) {
    options_out->use_new_transpose = *(bool*) options_in.attribute_array[tmpPos].value->addr;
    options_out->tr_setbyuser = true;
  }

  if (hasAttribute(options_in, "oldtranspose", tmpPos)==1) {
    options_out->use_old_transpose = *(bool*) options_in.attribute_array[tmpPos].value->addr;
    options_out->tr_setbyuser = true;
  }

  options_out->debug = *(bool*) getAttributeOrDefault(options_in, "debug", &((*options_out).debug));

  if (options_out->debug) {
    std::cout<<"eofunc: pcrit = "<<options_out->pcrit<<std::endl;
    if (options_out->debug) {
      std::cout << "anomalies being removed..." << '\n';
    } else {
      std::cout << "anomalies NOT being removed..." << '\n';
    }
  }


  return options_out;
}

extern "C" int eofunc(const ncomp_array & x_in, const ncomp_array & neval_in,
                      const attributes & options_in,
                      ncomp_array * x_out, attributes * attr_out) {
  int i_error = 0;

  // Sanity Checking
  if (x_in.ndim < 2) {
    std::cerr<<"eofunc: The input array must be at least two-dimensional"<<std::endl;
    return 1;
  }

  /* handle missing values */
  double missing_d_x_in;
  float missing_f_x_in;
  coerce_missing(x_in.type, x_in.has_missing, (ncomp_missing *)&(x_in.msg),
                 &missing_d_x_in,&missing_f_x_in);

  // Getting xData as double
  size_t x_nelem = prod(x_in.shape, x_in.ndim);
  double * dx = convert_to_with_copy_avoiding<double>(&x_in.addr, x_nelem, 0, x_in.type, NCOMP_DOUBLE);

  // Get number of eigenvalues and eigen vectors to be computed.
  // This is supposed to be a scalar.
  int* neval =  (int *) neval_in.addr;

  // Check Dimension sizes
  size_t msta = prod(x_in.shape, x_in.ndim-1);
  size_t ncol = msta;
  size_t nobs = x_in.shape[x_in.ndim-1];
  size_t nrow = nobs;

  // Sanity Checking
  if ( msta<1 || nobs <1) {
    std::cerr << "eofunc: The dimensions of the input array must both be at least 1" << std::endl;
    return 2;
  }

  // ignoring the following tests based on the discussion with Abhishek
  // if((nrow > INT_MAX) || (ncol > INT_MAX) ||
  //    (msta > INT_MAX) || (nobs > INT_MAX)) {
  //   NhlPError(NhlFATAL,NhlEUNKNOWN,"eofunc: one or more dimension sizes is greater than INT_MAX");
  //   return(NhlFATAL);
  // }

  int inrow = (int) nrow;
  int incol = (int) ncol;
  int inobs = (int) nobs;

  // processing options
  eofunc_options* options = extract_eofunc_options(options_in);

  /*
  * Create arrays to store non-missing data and to remove mean from
  * data before entering Fortran routines.
  */
  double * dx_strop = new double[nrow*ncol];
  double * xave = new double[ncol];
  double * xvar = new double[ncol];
  double * xdvar = new double[ncol];

  /*
   * Strip all grid points that have less than "PCRIT" valid values.
   * Create "dx_strip". This may have fewer columns/grid-pts
   * than the original "dx" array, if not all columns
   * had the minimum number of valid values.
   */
  size_t  mcsta = 0;
  double xsd, pcx, con;
  int kntx;
  for (size_t nc = 0; nc < ncol; ++nc) {
    /*
     * Statistics for this station/grid-point
     */
    dstat2_(  &dx[nrow*nc], &inrow, (double*) &x_in.msg.msg_double,
              &xave[nc], &xvar[nc], &xsd, &kntx, &i_error);

    /*
     * Eliminate stations/grid-points with less than pcrit % of data.
     */
    pcx = ((double)kntx/(double)nrow)*100.0;
    if (  (pcx < options->pcrit) ||
          (xsd <= 0.0)  ) {
      xave[nc] = x_in.msg.msg_double;
    }

    /*
     * Create anomalies. If jopt=1, then normalize the anomalies.
     * mcsta is the number of acceptable grid/station points (mcsta <= msta).
     */
    con = 1.0;
    if (  (options->jopt == 1) &&
          (xave[nc] != x_in.msg.msg_double) &&
          (xsd > 0.0) ) {
      con = 1.0/xsd;
    }

    /*
     * Work with anomalies: xdave=0.0 [or standardized anomalies]
     */
    if (xave[nc] != x_in.msg.msg_double) {
      /*
      * The following can produce too much output, so I've commented it out.
      *
      *      if(debug) {
      *          printf("nc = %d xave = %g\n", nc, xave[nc]);
      *      }
      */
      /*
      * Increment counter for acceptable points.
      */
      for( size_t nr = 0; nr < nobs; ++nr) {
        if(dx[nc*nrow+nr] != x_in.msg.msg_double) {
          if(!options->anomalies) {
            /*
             * User hasn't removed anomalies, so do it here.
             */
            dx_strip[mcsta*nrow+nr] = (dx[nc*nrow+nr] - xave[nc]) * con;
          } else {
            if(options->debug) {
              std::cout<<"anomalies NOT being removed..."<<std::endl;
            }
            /*
             * User has already removed anomalies, so leave alone.
             */
             dx_strip[mcsta*nrow+nr] = dx[nc*nrow+nr];
           }
         } else {
           dx_strip[mcsta*nrow+nr] = x_in.msg.msg_double;
         }
       }
       if(options->jopt == 0) {
         xdvar[mcsta] = xvar[nc];
       } else {
         xdvar[mcsta] = 1.0;
       }
       mcsta++;
      }
  }

  int imcsta = (int) mcsta;

  /*
   * Depending on the size of the rightmost 2D arrays being processed, and/or
   * the value of the transpose or oldtranspose attributes, we call one of
   * three different Fortran routines. These routines basically behave the
   * same, except two of them operate on a transposed version of the 2d array.
   */
  if(options->debug) {
    printf("eofunc: msta = %ld mcsta = %ld nobs = %ld\n", msta, mcsta, nobs);
  }
  /*
   * If one of the transpose attributes has not explicitly been set by the
   * user, then based on the sizes of the input array's dimensions, determine
   * whether to call a transpose routine or not.
   */
  if(!options->tr_setbyuser) {
    /*
     * If mcsta <= nrow, don't call transpose routine.
     */
    if(mcsta <= nrow) {
      options->use_new_transpose = false;    /* already the default */
      options->use_old_transpose = false;
      if(options->debug) {
        printf("eofunc: transpose set to False\n");
      }
    } else {
      /*
       * Since mcsta > nrow, call transpose routine.
       */
      options->use_new_transpose = true;
      options->use_old_transpose = false;
      if(options->debug) {
        printf("eofunc: transpose set to True\n");
      }
    }
  } else if(options->debug) {
    /*
     * User explicitly set one of the transpose attributes, so indicate
     * which one here. Note that if both oldtranspose and transpose are
     * set to True, transpose will take precedence.
     */
    if(options->use_new_transpose) {
      printf("eofunc: user set use_new_transpose to True\n");
    } else if (options->use_old_transpose) {
      printf("eofunc: user set use_old_transpose to True\n");
    } else {
      printf("eofunc: user set neither transpose attribute to True\n");
    }
  }

  size_t * dsizes_evec = new size_t[x_in.ndim];
  dsizes_evec[0] = *neval;
  for (int i = 0; i<x_in.ndim-1; ++i) {
    dsizes_evec[i+1] = x_in.shape[i];
  }
  int total_size_evec = (*neval) * ncol;

  /*
   * Allocate memory for various arrays.  Depending on which Fortran routine
   * will be called later, different quantities need to be allocated here.
   */
  double * xdatat = nullptr;
  double * wevec = nullptr;
  double * prncmp = nullptr;
  double * eval = nullptr;
  double * pcvar = nullptr;
  float * revec = nullptr;
  double * evec = nullptr;
  double * trace = nullptr;
  float * rpcvar;
  long long int llcssm;
  size_t lwork;
  int liwork, lifail;
  double * cssm;
  double * work;
  double * weval;
  int * iwork;
  int * ifail;
  int ilwork, iliwork, ilifail, lweval, icovcor;
  int iopt = 0;
  if (options->use_new_transpose) {
    xdatat = new double[nrow*mcsta];

    wevec = allocateAndInit<double>(*neval * mcsta, x_in.msg.msg_double);

    prncmp = new double[*neval*nrow];

    eval = allocateAndInit<double>(*neval, x_in.msg.msg_double);

    pcvar = allocateAndInit<double>(*neval, x_in.msg.msg_double);

    if (x_in.type != NCOMP_DOUBLE) {
      revec = new float[total_size_evec];
    } else {
      /*
       * If mcsta = msta, then we can use wevec as is. Otherwise, later we
       * need to copy wevec to locations in which the input was not missing.
       */
      if (mcsta != msta) {
        evec = new double[total_size_evec];
      }
    }

    /*
     * Transpose the input array.
     */
    size_t l1=0;
    for(size_t i = 0; i < mcsta; i++ ) {
      size_t l2 = i;
      for(size_t j = 0; j < nrow; j++ ) {
        xdatat[l2] = dx_strip[l1];
        l1++;
        l2+=mcsta;
      }
    }
  } else if (options->use_old_transpose) {
    trace = allocateAndInit<double>(1, x_in.msg.msg_double);

    evec = allocateAndInit<double>(total_size_evec, x_in.msg.msg_double);

    eval = allocateAndInit<double>(*neval, x_in.msg.msg_double);

    pcvar = allocateAndInit<double>(*neval, x_in.msg.msg_double);

    xdatat = new double[nrow*mcsta]

    if (x_in.type != NCOMP_DOUBLE) {
      revec = new float[total_size_evec];
    }
  } else {
    /*
     * eofcov routine
     *
     * Allocate space needed for various arrays.
     */
    wevec = allocateAndInit<double>(total_size_evec, x_in.msg.msg_double);

    trace = allocateAndInit<double>(1, x_in.msg.msg_double);

    eval = allocateAndInit<double>(*neval, x_in.msg.msg_double);

    rpcvar = allocateAndInit<float>(*neval, x_in.msg.msg_double);

    if (x_in.type != NCOMP_DOUBLE) {
      revec = new float[total_size_evec];
    } else {
      /*
       * If mcsta = msta, then we can use wevec as is. Otherwise, later we
       * need to copy wevec to locations in which the input was not missing.
       */
      if (mcsta != msta) {
        std::vector<double> tmp_evec(total_size_evec);
        evec = new double[total_size_evec];
      }
    }

    /*
     * Check sizes of work arrays that need to be passed to Fortran
     * routine below.
     */
    llcssm = mcsta*(mcsta+1)/2;
    lwork  = 8*mcsta;
    liwork = 5*mcsta;
    lifail = mcsta;

    ilwork  = (int) lwork;
    iliwork = (int) liwork;
    ilifail = (int) lifail;

    /*
     * Create some work arrays.  This is necessary to avoid having
     * these arrays created dynamically in the Fortran file (which makes
     * it Fortran 90, and unportable to some systems.
     */
    lweval = lifail;

    cssm = new double[llcssm];

    work = new double[lwork];

    weval = new double[lweval];

    iwork = new int[liwork];

    ifail = new in[lifail];
  }

  /*
   * Call the Fortran 77 version of appropriate routine.
   */
  if (options->use_new_transpose) {
    icovcor = 0;
    deof11_(xdatat,&imcsta,&inrow,neval,&icovcor,
            (double *)&x_in.msg.msg_double,eval,wevec,pcvar,prncmp);
  } else if (options->use_old_transpose) {
    xrveoft_( dx_strip.data(),xdatat,&inrow,&incol,&inobs,&imcsta,
              (double *)&x_in.msg.msg_double,neval,eval,evec,
              pcvar,trace,xdvar.data(),xave.data(),&options->jopt,&i_error);
  } else {
    ddrveof_( dx_strip.data(),&inrow,&incol,&inobs,&imcsta,
              (double *)&x_in.msg.msg_double,neval,eval,wevec,rpcvar,
              trace,&iopt,&options->jopt,cssm,&llcssm,work,&ilwork,
              weval,iwork,&iliwork,ifail,&ilifail,&i_error);
  }

  /*
   * If we used the "old" transpose routine, then the returned eigenvectors
   * have already been returned to the original-sized array with all the
   * missing values in the correct locations.  All we need to do here is
   * convert to float if necessary.
   */
  size_t nc2;
  if(options->use_old_transpose) {
    if(x_in.type != NCOMP_DOUBLE) {
      for( size_t i = 0; i < total_size_evec; ++i ) {
        revec[i] = (float)evec[i];
      }
    }
  } else {
    /*
     * If we are dealing with the old eofcov routine, or the new SCRIPPS
     * routine, then we need to reshape the evec (or revec if float)
     * array.  Note  that for the old eofcov routine, wevec is actually
     * the same size as evec, whereas for the new routine, it's the same
     * size only if mcsta == msta.
     */

    if (mcsta < msta) {
      /*
       * First, make sure init to missing because not all values will be
       * filled in.
       *
       * This is the floating point (single precision) case.
       */
      if(x_in.type != NCOMP_DOUBLE) {
        for(size_t i = 0; i < total_size_evec; ++i) {
          revec[i] = x_in.msg.msg_float;
        }
        /*
         * Now copy over the appropriate values in the wevec array. Since the
         * wevec array is a different size depending on which routine you are
         * using, we have two different sections of code here.
         */
        if(options->use_new_transpose) {
          nc2 = 0;
          for( size_t nc = 0; nc < ncol; ++nc) {
            if (xave[nc] != x_in.msg.msg_double) {
              for( size_t ne = 0; ne < *neval; ++ne ) {
                revec[ne*ncol+nc] = (float)wevec[ne*mcsta+nc2];
              }
              nc2++;
            }
          }
        } else {
          nc2 = 0;
          for( size_t nc = 0; nc < ncol; ++nc) {
            if (xave[nc] != x_in.msg.msg_double) {
              for( size_t ne = 0; ne < *neval; ++ne ) {
                revec[ne*ncol+nc] = (float)wevec[ne*ncol+nc2];
              }
              nc2++;
            }
          }
        }
      } else {
        /*
         * This is the double precision case.
         */
        /*
         * First, make sure init to missing because not all values will be
         * filled in.
         */
        for(size_t i = 0; i < total_size_evec; ++i) {
          evec[i] = x_in.msg.msg_double;
        }
        /*
         * Now copy over the appropriate values in the wevec array. Since the
         * wevec array is a different size depending on which routine you are
         * using, we have two different sections of code here.
         */
        if(options->use_new_transpose) {
          nc2 = 0;
          for( size_t nc = 0; nc < ncol; ++nc) {
            if (xave[nc] != x_in.msg.msg_double) {
              for( size_t ne = 0; ne < *neval; ++ne ) {
                evec[ne*ncol+nc] = wevec[ne*mcsta+nc2];
              }
              nc2++;
            }
          }
        } else {
          nc2 = 0;
          for( size_t nc = 0; nc < ncol; ++nc) {
            if (xave[nc] != x_in.msg.msg_double) {
              for( size_t ne = 0; ne < *neval; ++ne ) {
                evec[ne*ncol+nc] = wevec[ne*ncol+nc2];
              }
              nc2++;
            }
          }
        }
      }
    } else {
      /*
       * mcsta = msta, so we just need to copy stuff over. It doesn't matter
       * whether we have called the old eofcov routine or the new eof SCRIPPS
       * routine, because if mcsta==msta, then wevec is the same size for
       * both routines.
       */
      if(x_in.type != NCOMP_DOUBLE) {
        for( size_t i = 0; i < total_size_evec; i++ ) {
          revec[i] = (float)wevec[i];
        }
      } else {
        evec = wevec;
      }
    }
  }

  /*
   * Check various possible error messages. The new transpose routine doesn't
   * have an ier.
   */
  if (!options->use_new_transpose && i_error != 0) {
    if (i_error == -1) {
      std::cerr<<"eofunc: cssm contains one or more missing values.\n(One or more series contains all missing values.)"<<std::endl;
      return i_error;
    }
    else if (i_error == -88) {
      std::cerr<<"eofunc: trace is equal to zero.\nAll data entries are missing or are equal to zero."<<std::endl;
      return i_error;
    }
    else if (i_error < 0) {
      std::cerr<<"eofunc: The "<<abs(i_error)<<"-th argument had an illegal value"<<std::endl;
      return i_error;
    }
    else {
      std::cerr<<"eofunc: "<<i_error<<" eigenvectors failed to converge"<<std::endl;
      return i_error;
    }
  }

  /*
   * This is the start of a rather large if-else statement. It is based
   * on whether you are returning floats or doubles.
   */
   std::vector<single_attribute> tmp_attr_out;
  if(x_in.type != NCOMP_DOUBLE) {
    /*
     * Set up return value.
     */
     // x_out is the return_md in NCL code.
     x_out = ncomp_array_alloc((void *) revec, NCOMP_FLOAT, dsizes_evec.size(),dsizes_evec.data());
    /*
     * Only return the eigenvalues if the appropriate option has been set.
     */
    if(options->return_eval) {
      /*
       * Coerce eval to float.
       */
      float * reval =  convert_to_with_copy_avoiding(eval, *neval, 0, NCOMP_DOUBLE, NCOMP_FLOAT);

      /*
       * If we didn't use the SCRIPPS routine, then the eigenvalues
       * returned are okay as is. Otherwise, we have to apply a scale
       * factor and return both the original values and the scaled values.
       */
      if(options->use_new_transpose) {
        float * scaled_reval;
        std::vector<float> tmp_scaled_reval(*neval);
        scaled_reval = tmp_scaled_reval.data();
        float scale_factor = (mcsta-1)/(nrow-1);
        for( size_t i = 0; i < *neval; ++i ) {
          scaled_reval[i] = scale_factor * reval[i];
        }

      /*
       * First return original eigenvalues as "eval_transpose".
       */
      tmp_attr_out.push_back(create_single_attribute("eval_transpose", reval, NCOMP_FLOAT, 1, neval));

      /*
       * Now return scaled eigenvalues as simply "eval".
       */
      tmp_attr_out.push_back(create_single_attribute("eval", scaled_reval, NCOMP_FLOAT, 1, neval));
      } else {
        /*
         * We didn't call the tranpose routine, so we only need to return
         * one set of eigenvalues.
         */
        tmp_attr_out.push_back(create_single_attribute("eval", reval, NCOMP_FLOAT, 1, neval));
      }
    }
    /*
     * Only return the trace if the appropriate option has been set.
     * The new transpose routine doesn't return trace.
     */
    if(!options->use_new_transpose) {
      if(return_trace) {
        /*
         * Coerce trace to float.
         */
        std::vector<float> rtrace(1, (float) *trace)
        float * tmp_rtrace;
        tmp_rtrace = rtrace.data();
        int dims[1] {1};
        tmp_attr_out.push_back(create_single_attribute("trace", tmp_rtrace, NCOMP_FLOAT, 1, dims));
      }
    }
  } else {
    /*
     *  Return doubles.
     */
    x_out = ncomp_array_alloc((void *) evec, NCOMP_DOUBLE, dsizes_evec.size(),dsizes_evec.data());  // wait for abhishek to make sure after evec is out of scope
                                                                                                    //its memory is not released and x_out is still valid;

    /*
     * Only return the eigenvalues if the appropriate option has been set.
     */
    if(options->return_eval) {
      /*
       * If we didn't use the SCRIPPS routine, then the eigenvalues
       * returned are okay as is. Otherwise, we have to apply a scale
       * factor and return both the original values and the scaled values.
       */
      if(options->use_new_transpose) {
        double * scaled_eval;
        std::vector<double> tmp_scaled_eval(*neval);
        scaled_eval = tmp_scaled_eval.data();
        double scale_factor = (mcsta-1)/(nrow-1);
        for( size_t i = 0; i < *neval; ++i ) {
          scaled_eval[i] = scale_factor * eval[i];
        }
        /*
         * First return original eigenvalues as "eval_transpose".
         */
        tmp_attr_out.push_back(create_single_attribute("eval_transpose", eval, NCOMP_DOUBLE, 1, neval));
        /*
         * Now return scaled eigenvalues as simply "eval".
         */
        tmp_attr_out.push_back(create_single_attribute("eval", scaled_eval, NCOMP_DOUBLE, 1, neval));
      } else {
        /*
         * We didn't call the tranpose routine, so we only need to return
         * one set of eigenvalues.
         */
        tmp_attr_out.push_back(create_single_attribute("eval", eval, NCOMP_DOUBLE, 1, neval));
      }
    }

    if(!options->use_new_transpose) {
      if(options->return_trace) {
        otrace = (float *)calloc(1,sizeof(float));
        std::vector<float> tmp_rtrace(1)
        *rtrace = (float)(*trace);
        int dims[1] {1};
        tmp_attr_out.push_back(create_single_attribute("trace", trace, NCOMP_DOUBLE, 1, dims));
      }
    }
  }

  /*
   * Return pcvar as float no matter what.
   */
   float * rpcvar;
  if(options->use_old_transpose || options->use_new_transpose) {
    std::vector(float) tmp_rpcvar(*neval);
    rpcvar = tmp_rpcvar.data();
    for( size_t i = 0; i < *neval; ++i ) {
      rpcvar[i] = (float)pcvar[i];
    }
  }
  tmp_attr_out.push_back(create_single_attribute("pcvar", rpcvar, NCOMP_FLOAT, 1, neval))

  /*
   * Only return "pcrit" if it was set by the user and we called one
   * of the transpose routines. The type returned is a float or a double,
   * depending on what pcrit was set to in the input.
   */
  if( (options->use_new_transpose || options->use_old_transpose) &&
      options->return_pcrit) {
    int dims[1] {1};
    std::vector()
    tmp_attr_out.push_back(create_single_attribute("pcvar", rpcvar, NCOMP_FLOAT, 1, neval))
    dsizes[0] = 1;
    // if(type_pcrit == NCL_float) {
    //   att_md = _NclCreateVal(
    //                          NULL,
    //                          NULL,
    //                          Ncl_MultiDValData,
    //                          0,
    //                          (void*)rpcrit,
    //                          NULL,
    //                          1,
    //                          dsizes,
    //                          TEMPORARY,
    //                          NULL,
    //                          (NclObjClass)nclTypefloatClass
    //                          );
    // }
    // else {
    //   att_md = _NclCreateVal(
    //                          NULL,
    //                          NULL,
    //                          Ncl_MultiDValData,
    //                          0,
    //                          (void*)pcrit,
    //                          NULL,
    //                          1,
    //                          dsizes,
    //                          TEMPORARY,
    //                          NULL,
    //                          (NclObjClass)nclTypedoubleClass
    //                          );
    // }
    // _NclAddAtt(
    //            att_id,
    //            "pcrit",
    //            att_md,
    //            NULL
    //            );
  }

  /*
   * "matrix" indicates whether the covariance or correlation matrix
   * was used.
   */
  // if(options->jopt == 0) {
  //   cmatrix = (char *)calloc(11,sizeof(char));
  //   strcpy(cmatrix,"covariance");
  // }
  // else {
  //   cmatrix = (char *)calloc(12,sizeof(char));
  //   strcpy(cmatrix,"correlation");
  // }
  // matrix  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  // *matrix = NrmStringToQuark(cmatrix);
  //
  // dsizes[0] = 1;
  // att_md = _NclCreateVal(
  //                        NULL,
  //                        NULL,
  //                        Ncl_MultiDValData,
  //                        0,
  //                        (void*)matrix,
  //                        NULL,
  //                        1,
  //                        dsizes,
  //                        TEMPORARY,
  //                        NULL,
  //                        (NclObjClass)nclTypestringClass
  //                        );
  // _NclAddAtt(
  //            att_id,
  //            "matrix",
  //            att_md,
  //            NULL
  //            );

// /*
//  * "method" indicates whether the transpose routine was called or not.
//  */
//   if(use_new_transpose) {
//     cmethod = (char *)calloc(10,sizeof(char));
//     strcpy(cmethod,"transpose");
//   }
//   else if(use_old_transpose) {
//     cmethod = (char *)calloc(14,sizeof(char));
//     strcpy(cmethod,"old_transpose");
//   }
//   else {
//     cmethod = (char *)calloc(13,sizeof(char));
//     strcpy(cmethod,"no transpose");
//   }
//   method  = (NclQuark*)NclMalloc(sizeof(NclQuark));
//   *method = NrmStringToQuark(cmethod);

  // dsizes[0] = 1;
  // att_md = _NclCreateVal(
  //                        NULL,
  //                        NULL,
  //                        Ncl_MultiDValData,
  //                        0,
  //                        (void*)method,
  //                        NULL,
  //                        1,
  //                        dsizes,
  //                        TEMPORARY,
  //                        NULL,
  //                        (NclObjClass)nclTypestringClass
  //                        );
  // _NclAddAtt(
  //            att_id,
  //            "method",
  //            att_md,
  //            NULL
  //            );

  // tmp_var = _NclVarCreate(
  //                         NULL,
  //                         NULL,
  //                         Ncl_Var,
  //                         0,
  //                         NULL,
  //                         return_md,
  //                         NULL,
  //                         att_id,
  //                         NULL,
  //                         RETURNVAR,
  //                         NULL,
  //                         TEMPORARY
  //                         );


  /*
   * Return output grid and attributes to NCL.
   */
  // return_data.kind = NclStk_VAR;
  // return_data.u.data_var = tmp_var;
  // _NclPlaceReturn(return_data);


  // Cleaning up
  delete options;
  return i_error;
}
