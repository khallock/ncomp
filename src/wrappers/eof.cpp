#include "ncomp/types.h"
#include "ncomp/util.h"
#include "ncomp_internal/util.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <cstring>

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
  bool return_eval {true};
  bool return_trace {false};
  bool anomalies {false};
  bool use_new_transpose {true};
  bool use_old_transpose {false};
  bool tr_setbyuser {false};
  bool debug {false};
} eofunc_options;

eofunc_options* extract_eofunc_options(const ncomp_attributes * options_in) {
  eofunc_options* options_out = new eofunc_options;
  if  ((options_in == nullptr) ||
      (options_in->nAttribute == 0)) {
    return options_out;
  }

  options_out->jopt = *(int *) getAttributeOrDefault(options_in, "jopt", &(options_out->jopt));
  if ((options_out->jopt != 0) && (options_out->jopt != 1)) {
    options_out->jopt = 0; // jopt must be either 0 or 1
  }

  int tmpPos = -1;
  if (hasAttribute(options_in, "pcrit", &tmpPos)==1) {
    options_out->pcrit = *(double*) options_in->attribute_array[tmpPos]->value->addr;
    options_out->return_pcrit = true;
    if ((options_out->pcrit < 0.0) || (options_out->pcrit > 100.0)) {
      options_out->pcrit = 50.0; // pcrit must be between 0.0 and 100.0; default value is 50
    }
  }

  options_out->return_eval = *(bool*) getAttributeOrDefault(options_in, "return_eval", &(options_out->return_eval));

  options_out->return_trace = *(bool*) getAttributeOrDefault(options_in, "return_trace", &(options_out->return_trace));

  options_out->anomalies = *(bool*) getAttributeOrDefault(options_in, "anomalies", &(options_out->anomalies));

  if (hasAttribute(options_in, "transpose", &tmpPos)==1) {
    options_out->use_new_transpose = *(bool*) options_in->attribute_array[tmpPos]->value->addr;
    options_out->tr_setbyuser = true;
  } else if (hasAttribute(options_in, "oldtranspose", &tmpPos)==1) { // we should set either transpose or old-transposed
                                                                    // in this case, if transpose is already provided,
                                                                    // we are not checking if old-tranpose is set.
                                                                    // We only check if old-transpose is set, in case if
                                                                    // transpose is not set.
    options_out->use_old_transpose = *(bool*) options_in->attribute_array[tmpPos]->value->addr;
    options_out->use_new_transpose = false; // making sure that only one of them is set.
    options_out->tr_setbyuser = true;
  }

  options_out->debug = *(bool*) getAttributeOrDefault(options_in, "debug", &(options_out->debug));

  if (options_out->debug) {
    std::cout<<"eofunc: pcrit = "<<options_out->pcrit<<std::endl;
    if (options_out->anomalies) {
      std::cout << "anomalies NOT being removed..." << '\n';
    } else {
      std::cout << "anomalies being removed..." << '\n';
    }
  }


  return options_out;
}

extern "C" int eofunc(const ncomp_array * x_in, const int neval_in,
                      const ncomp_attributes * options_in,
                      ncomp_array * x_out, ncomp_attributes * attrList_out) {
  int i_error = 0;

  // Sanity Checking
  if (x_in->ndim < 2) {
    std::cerr<<"eofunc: The input array must be at least two-dimensional"<<std::endl;
    return 1;
  }

  /* handle missing values */
  double missing_d_x_in;
  float missing_f_x_in;
  coerce_missing(x_in->type, x_in->has_missing, (ncomp_missing *)&(x_in->msg),
                 &missing_d_x_in,&missing_f_x_in);

  // Getting xData as double
  size_t x_nelem = prod(x_in->shape, x_in->ndim);
  double * dx = convert_to_with_copy_avoiding<double>((void *)x_in->addr, x_nelem, 0, x_in->type, NCOMP_DOUBLE);

  // Get number of eigenvalues and eigen vectors to be computed.
  // This is supposed to be a scalar.
  int neval =  neval_in;

  // Check Dimension sizes
  size_t msta = prod(x_in->shape, x_in->ndim-1);
  size_t ncol = msta;
  size_t nobs = x_in->shape[x_in->ndim-1];
  size_t nrow = nobs;
  size_t total_size_x_in = ncol*nrow;

  // Sanity Checking
  if ( msta<1 || nobs <1) {
    std::cerr << "eofunc: The dimensions of the input array must both be at least 1" << std::endl;
    return 2;
  }

  if((nrow > INT_MAX) || (ncol > INT_MAX) ||
     (msta > INT_MAX) || (nobs > INT_MAX)) {
    std::cerr<<"eofunc: one or more dimension sizes is greater than INT_MAX"<<std::endl;
    return(NCOMP_RETURN_FATAL);
  }

  int inrow = (int) nrow;
  int incol = (int) ncol;
  int inobs = (int) nobs;

  // processing options
  std::unique_ptr<eofunc_options> options (extract_eofunc_options(options_in));

  /*
  * Create arrays to store non-missing data and to remove mean from
  * data before entering Fortran routines.
  */
  // double * dx_strip = new double[nrow*ncol];
  // double * xave = new double[ncol];
  // double * xvar = new double[ncol];
  // double * xdvar = new double[ncol];
  std::unique_ptr<double[]> dx_strip(new double[nrow*ncol]);
  std::unique_ptr<double[]> xave(new double[ncol]);
  std::unique_ptr<double[]> xvar(new double[ncol]);
  std::unique_ptr<double[]> xdvar(new double[ncol]);

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
    dstat2_(  &dx[nrow*nc], &inrow, &missing_d_x_in,
              &xave[nc], &xvar[nc], &xsd, &kntx, &i_error);

    /*
     * Eliminate stations/grid-points with less than pcrit % of data.
     */
    pcx = ((double)kntx/(double)nrow)*100.0;
    if (  (pcx < options->pcrit) ||
          (xsd <= 0.0)  ) {
      xave[nc] = missing_d_x_in;
    }

    /*
     * Create anomalies. If jopt=1, then normalize the anomalies.
     * mcsta is the number of acceptable grid/station points (mcsta <= msta).
     */
    con = 1.0;
    if (  (options->jopt == 1) &&
          (xave[nc] != missing_d_x_in) &&
          (xsd > 0.0) ) {
      con = 1.0/xsd;
    }

    /*
     * Work with anomalies: xdave=0.0 [or standardized anomalies]
     */
    if (xave[nc] != missing_d_x_in) {
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
        if(dx[nc*nrow+nr] != missing_d_x_in) {
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
           dx_strip[mcsta*nrow+nr] = missing_d_x_in;
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

  if(mcsta > INT_MAX) {
    std::cerr<<"eofunc: one or more dimension sizes is greater than INT_MAX"<<std::endl;
    return(NCOMP_RETURN_FATAL);
  }
  int imcsta = (int) mcsta;

  /*
   * Depending on the size of the rightmost 2D arrays being processed, and/or
   * the value of the transpose or oldtranspose ncomp_attributes, we call one of
   * three different Fortran routines. These routines basically behave the
   * same, except two of them operate on a transposed version of the 2d array.
   */
  if(options->debug) {
    printf("eofunc: msta = %ld mcsta = %ld nobs = %ld\n", msta, mcsta, nobs);
  }
  /*
   * If one of the transpose ncomp_attributes has not explicitly been set by the
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
     * User explicitly set one of the transpose ncomp_attributes, so indicate
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

  std::unique_ptr<size_t[]> dsizes_evec(new size_t[x_in->ndim]);
  dsizes_evec[0] = neval;
  for (int i = 0; i<x_in->ndim-1; ++i) {
    dsizes_evec[i+1] = x_in->shape[i];
  }
  int total_size_evec = (neval) * ncol;

  /*
   * Allocate memory for various arrays.  Depending on which Fortran routine
   * will be called later, different quantities need to be allocated here.
   */
  std::unique_ptr<double[]> xdatat;
  std::unique_ptr<double[]> wevec;
  std::unique_ptr<double[]> prncmp;
  std::unique_ptr<double[]> eval;
  std::unique_ptr<double[]> pcvar;
  std::unique_ptr<float[]> revec;
  std::unique_ptr<double[]> evec;
  std::unique_ptr<double[]> trace;
  std::unique_ptr<float[]> rpcvar;
  long long int llcssm;
  size_t lwork;
  int liwork, lifail;
  std::unique_ptr<double[]> cssm;
  std::unique_ptr<double[]> work;
  std::unique_ptr<double[]> weval;
  std::unique_ptr<int[]> iwork;
  std::unique_ptr<int[]> ifail;
  int ilwork, iliwork, ilifail, lweval, icovcor;
  int iopt = 0;
  if (options->use_new_transpose) {
    xdatat.reset(new double[nrow*mcsta]);

    wevec.reset(allocateAndInit<double>(neval * mcsta, missing_d_x_in));

    prncmp.reset(new double[neval * nrow]);

    eval.reset(allocateAndInit<double>(neval, missing_d_x_in));

    pcvar.reset(allocateAndInit<double>(neval, missing_d_x_in));

    if (x_in->type != NCOMP_DOUBLE) {
      revec.reset(new float[total_size_evec]);
    } else {
      /*
       * If mcsta = msta, then we can use wevec as is. Otherwise, later we
       * need to copy wevec to locations in which the input was not missing.
       */
      if (mcsta != msta) {
        evec.reset(new double[total_size_evec]);
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
    trace.reset(allocateAndInit<double>(1, missing_d_x_in));

    evec.reset(allocateAndInit<double>(total_size_evec, missing_d_x_in));

    eval.reset(allocateAndInit<double>(neval, missing_d_x_in));

    pcvar.reset(allocateAndInit<double>(neval, missing_d_x_in));

    xdatat.reset(new double[nrow*mcsta]);

    if (x_in->type != NCOMP_DOUBLE) {
      revec.reset(new float[total_size_evec]);
    }
  } else {
    /*
     * eofcov routine
     *
     * Allocate space needed for various arrays.
     */
    wevec.reset(allocateAndInit<double>(total_size_evec, missing_d_x_in));

    trace.reset(allocateAndInit<double>(1, missing_d_x_in));

    eval.reset(allocateAndInit<double>(neval, missing_d_x_in));

    rpcvar.reset(allocateAndInit<float>(neval, missing_f_x_in));

    if (x_in->type != NCOMP_DOUBLE) {
      revec.reset(new float[total_size_evec]);
    } else {
      /*
       * If mcsta = msta, then we can use wevec as is. Otherwise, later we
       * need to copy wevec to locations in which the input was not missing.
       */
      if (mcsta != msta) {
        evec.reset(new double[total_size_evec]);
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

    if((lwork > INT_MAX) || (liwork > INT_MAX) || (lifail > INT_MAX)) {
      std::cerr<<"eofunc: one or more dimension sizes is greater than INT_MAX"<<std::endl;
      return(NCOMP_RETURN_FATAL);
    }

    ilwork  = (int) lwork;
    iliwork = (int) liwork;
    ilifail = (int) lifail;

    /*
     * Create some work arrays.  This is necessary to avoid having
     * these arrays created dynamically in the Fortran file (which makes
     * it Fortran 90, and unportable to some systems.
     */
    lweval = lifail;

    cssm.reset(new double[llcssm]);

    work.reset(new double[lwork]);

    weval.reset(new double[lweval]);

    iwork.reset(new int[liwork]);

    ifail.reset(new int[lifail]);
  }

  /*
   * Call the Fortran 77 version of appropriate routine.
   */
  if (options->use_new_transpose) {
    icovcor = 0;
    deof11_(xdatat.get(),&imcsta,&inrow,&neval,&icovcor,
            &missing_d_x_in,eval.get(),wevec.get(),pcvar.get(),prncmp.get());
  } else if (options->use_old_transpose) {
    xrveoft_( dx_strip.get(),xdatat.get(),&inrow,&incol,&inobs,&imcsta,
              &missing_d_x_in,&neval,eval.get(),evec.get(),
              pcvar.get(),trace.get(),xdvar.get(),xave.get(),&options->jopt,&i_error);
  } else {
    ddrveof_( dx_strip.get(),&inrow,&incol,&inobs,&imcsta,
              &missing_d_x_in,&neval,eval.get(),wevec.get(),rpcvar.get(),
              trace.get(),&iopt,&options->jopt,cssm.get(),&llcssm,work.get(),&ilwork,
              weval.get(),iwork.get(),&iliwork,ifail.get(),&ilifail,&i_error);
  }

  /*
   * If we used the "old" transpose routine, then the returned eigenvectors
   * have already been returned to the original-sized array with all the
   * missing values in the correct locations.  All we need to do here is
   * convert to float if necessary.
   */
  size_t nc2;
  if(options->use_old_transpose) {
    if(x_in->type != NCOMP_DOUBLE) {
      for( size_t i = 0; i < total_size_evec; ++i ) {
        revec[i] = (float)evec[i];
      }
      // delete[] evec;
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
      if(x_in->type != NCOMP_DOUBLE) {
        std::fill_n(revec.get(), total_size_evec, missing_f_x_in);
        /*
         * Now copy over the appropriate values in the wevec array. Since the
         * wevec array is a different size depending on which routine you are
         * using, we have two different sections of code here.
         */
        if(options->use_new_transpose) {
          nc2 = 0;
          for( size_t nc = 0; nc < ncol; ++nc) {
            if (xave[nc] != missing_d_x_in) {
              for( size_t ne = 0; ne < neval; ++ne ) {
                revec[ne*ncol+nc] = (float)wevec[ne*mcsta+nc2];
              }
              ++nc2;
            }
          }
        } else {
          nc2 = 0;
          for( size_t nc = 0; nc < ncol; ++nc) {
            if (xave[nc] != missing_d_x_in) {
              for( size_t ne = 0; ne < neval; ++ne ) {
                revec[ne*ncol+nc] = (float)wevec[ne*ncol+nc2];
              }
              ++nc2;
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
         std::fill_n(evec.get(), total_size_evec, missing_d_x_in);

        /*
         * Now copy over the appropriate values in the wevec array. Since the
         * wevec array is a different size depending on which routine you are
         * using, we have two different sections of code here.
         */
        if(options->use_new_transpose) {
          nc2 = 0;
          for( size_t nc = 0; nc < ncol; ++nc) {
            if (xave[nc] != missing_d_x_in) {
              for( size_t ne = 0; ne < neval; ++ne ) {
                evec[ne*ncol+nc] = wevec[ne*mcsta+nc2];
              }
              ++nc2;
            }
          }
        } else {
          nc2 = 0;
          for( size_t nc = 0; nc < ncol; ++nc) {
            if (xave[nc] != missing_d_x_in) {
              for( size_t ne = 0; ne < neval; ++ne ) {
                evec[ne*ncol+nc] = wevec[ne*ncol+nc2];
              }
              ++nc2;
            }
          }
        }
      }
      // delete[] wevec;
    } else {
      /*
       * mcsta = msta, so we just need to copy stuff over. It doesn't matter
       * whether we have called the old eofcov routine or the new eof SCRIPPS
       * routine, because if mcsta==msta, then wevec is the same size for
       * both routines.
       */
      if(x_in->type != NCOMP_DOUBLE) {
        for( size_t i = 0; i < total_size_evec; ++i ) {
          revec[i] = (float)wevec[i];
        }
        // delete[] wevec;
      } else {
        evec = std::move(wevec);
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
 * Free unneeded memory.
 */
  // if(x_in->type != NCOMP_DOUBLE) delete[] dx;
  // delete[] dx_strip;
  // delete[] xave;
  // delete[] xvar;
  // delete[] xdvar;
  // if(!options->use_new_transpose && !options->use_old_transpose) {
  //   delete[] work;
  //   delete[] cssm;
  //   delete[] weval;
  //   delete[] iwork;
  //   delete[] ifail;
  // }
  // else {
  //   delete[] xdatat;
  //   if(options->use_new_transpose) delete[] prncmp;
  // }

  /*
   * This is the start of a rather large if-else statement. It is based
   * on whether you are returning floats or doubles.
   */
  std::vector<ncomp_single_attribute *> tmp_attr_out;
  if(x_in->type != NCOMP_DOUBLE) {
    /*
     * Set up return value.
     */
     // x_out is the return_md in NCL code.
     ncomp_array_copy(
       ncomp_array_alloc((void *) revec.release(), NCOMP_FLOAT, x_in->ndim,dsizes_evec.get()),
       x_out);
     x_out->has_missing = x_in->has_missing;
     x_out->msg.msg_float = missing_f_x_in;

    /*
     * Only return the eigenvalues if the appropriate option has been set.
     */
    if(options->return_eval) {
      /*
       * Coerce eval to float.
       */
      float * reval =  new float[neval]; // this is going to be returned. No need to delet or use unique_ptr
      for (size_t i = 0; i < neval; ++i) {
        reval[i] = (float) eval[i];
      }

      /*
       * If we didn't use the SCRIPPS routine, then the eigenvalues
       * returned are okay as is. Otherwise, we have to apply a scale
       * factor and return both the original values and the scaled values.
       */
      if(options->use_new_transpose) {
        float * scaled_reval = new float[neval]; // this is going to be returned. No need to delet or use unique_ptr
        float scale_factor = (mcsta-1)/(nrow-1);
        for( size_t i = 0; i < neval; ++i ) {
          scaled_reval[i] = scale_factor * reval[i];
        }

      /*
       * First return original eigenvalues as "eval_transpose".
       */
      size_t dims[1] {(size_t)neval};
      tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "eval_transpose", reval, NCOMP_FLOAT, 1, dims));

      /*
       * Now return scaled eigenvalues as simply "eval".
       */
      tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "eval", scaled_reval, NCOMP_FLOAT, 1, dims));
      } else {
        /*
         * We didn't call the tranpose routine, so we only need to return
         * one set of eigenvalues.
         */
        size_t dims[1] {(size_t)neval};
        tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "eval", reval, NCOMP_FLOAT, 1, dims));
      }
      // delete[] eval;
    }
    /*
     * Only return the trace if the appropriate option has been set.
     * The new transpose routine doesn't return trace.
     */
    if(!options->use_new_transpose) {
      if(options->return_trace) {
        /*
         * Coerce trace to float.
         */
        float * rtrace = new float[1]{(float) trace[0]};
        size_t dims[1] {1};
        tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "trace", rtrace, NCOMP_FLOAT, 1, dims));
      }
      // delete[] trace;
    }
  } else {
    /*
     *  Return doubles.
     */

    ncomp_array_copy(
      ncomp_array_alloc((void *) evec.release(), NCOMP_DOUBLE, x_in->ndim,dsizes_evec.get()),
      x_out);
    x_out->has_missing = x_in->has_missing;
    x_out->msg.msg_double = missing_d_x_in;

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
        double * scaled_eval = new double[neval];
        double scale_factor = (mcsta-1)/(nrow-1);
        for( size_t i = 0; i < neval; ++i ) {
          scaled_eval[i] = scale_factor * eval[i];
        }
        /*
         * First return original eigenvalues as "eval_transpose".
         */
        tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "eval_transpose", eval.release(), NCOMP_DOUBLE, 1, dsizes_evec.get()));
        /*
         * Now return scaled eigenvalues as simply "eval".
         */
        tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "eval", scaled_eval, NCOMP_DOUBLE, 1, dsizes_evec.get()));
      } else {
        /*
         * We didn't call the tranpose routine, so we only need to return
         * one set of eigenvalues.
         */
        tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "eval", eval.release(), NCOMP_DOUBLE, 1, dsizes_evec.get()));
      }
    }

    if(!options->use_new_transpose) {
      if(options->return_trace) {
        size_t dims[1] {1};
        tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "trace", trace.release(), NCOMP_DOUBLE, 1, dims));
      }
      // else {
      //   delete[] trace;
      // }
    }
  }

  /*
   * Return pcvar as float no matter what.
   */
  if(options->use_old_transpose || options->use_new_transpose) {
    rpcvar.reset(new float[neval]);
    // rpcvar = new float[neval]; // Possible bug in NCL as well.
    for( size_t i = 0; i < neval; ++i ) {
      rpcvar[i] = (float)pcvar[i];
    }
  }
  tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "pcvar", rpcvar.release(), NCOMP_FLOAT, 1, dsizes_evec.get()));

  /*
   * Only return "pcrit" if it was set by the user and we called one
   * of the transpose routines. The type returned is a float or a double,
   * depending on what pcrit was set to in the input.
   */
  if( (options->use_new_transpose || options->use_old_transpose) &&
      options->return_pcrit) {
    double * tmp_pcrit = new double[1];
    *tmp_pcrit = options->pcrit;
    tmp_attr_out.push_back(create_ncomp_single_attribute_from_scalar((char *) "pcrit", tmp_pcrit, NCOMP_DOUBLE));
  }

  /*
   * "matrix" indicates whether the covariance or correlation matrix
   * was used.
   */
  char * cmatrix;
  if(options->jopt == 0) {
    cmatrix = new char[11];
    strcpy(cmatrix,"covariance");
    size_t dims[1] {1};
    tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "matrix", cmatrix, NCOMP_CHAR, 1, dims));
  }
  else {
    cmatrix = new char[12];
    strcpy(cmatrix,"correlation");
    size_t dims[1] {1};
    tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "matrix", cmatrix, NCOMP_CHAR, 1, dims));
  }


  /*
   * "method" indicates whether the transpose routine was called or not.
   */
  char * cmethod;
  if(options->use_new_transpose) {
    cmethod = new char[10];
    strcpy(cmethod,"transpose");
    size_t dims[1] {1};
    tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "method", cmethod, NCOMP_CHAR, 1, dims));
  }
  else if(options->use_old_transpose) {
    cmethod = new char[14];
    strcpy(cmethod,"old_transpose");
    size_t dims[1] {1};
    tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "method", cmethod, NCOMP_CHAR, 1, dims));
  }
  else {
    cmethod = new char[13];
    strcpy(cmethod,"no transpose");
    size_t dims[1] {1};
    tmp_attr_out.push_back(create_ncomp_single_attribute((char *) "method", cmethod, NCOMP_CHAR, 1, dims));
  }

  collectAttributeList(tmp_attr_out, attrList_out);

  // Cleaning up
  // delete[] dsizes_evec;
  // delete options;
  return i_error;
}
