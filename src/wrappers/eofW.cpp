#include "ncomp/types.h"
#include "ncomp/util.h"
#include "ncomp_internal/util.hpp"

extern "C" void ddrveof_(double *,int *,int *,int *,int *,
                         double *,int *,double *, double *,
                         float*,double *,int *,int *,double*,
                         long long int *, double *,int *,
                         double *,int *,int *,int *,int *,int *);

extern "C" void crveoft_( double *dx_strip, double *dx_strip_t,
                         int *nrow, int *ncol, int *nrobs,
                         int *mcsta, double *xmsg, int *neval,
                         double *eval, double *evec,
                         double *pcvar, double *trace,
                         double *xdvar, double *xave,
                         int *jopt, int *ier);

extern "C" void deof11_( double *, int *, int *, int *, int *,
                         double  *, double *, double *, double *,
                         double *);

extern "C" void dstat2_( double *, int *, double *, double *,
                         double *, double *, int *, int *);

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
