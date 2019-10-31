! ----------------------------------------------------------------------
      SUBROUTINE DSTAT2 (X, NPTS, XMSG, XMEAN, XVAR, XSD, NPTUSED, IER)
      DOUBLEPRECISION X1
      DOUBLEPRECISION X2
      DOUBLEPRECISION XN
      DOUBLEPRECISION X3

! this routine will calculate estimates of the first two  moments
! .   of the vector x containing missing data.

! input arguments:
! .   x        - input vector
! .   npts     - length of x
! .   xmsg     - missing code: if there are no msg values
! .                            set xmsg to some value which will
! .                            not be encountered.

! output arguments:
! .   xmean    - mean of x (first moment)
! .   xvar     - sample variance
! .   xsd      - sqrt( xvar )
! .   nptused  - no. of points used to calculate the estimates
! .   ier      - if (ier.ne.0) an error has occurred

! note :
! .   uncalculated quantities are set to xmsg (nptused set to zero)

      INTEGER NPTS, NPTUSED, IER
      DOUBLEPRECISION X (1:NPTS), XMSG, XMEAN, XVAR, XSD

      XMEAN = XMSG
      XVAR = XMSG
      XSD = XMSG
      NPTUSED = 0

      IF (NPTS.LT.1) THEN
        IER = 1
        RETURN
      ENDIF

      IER = 0
      X1 = 0.D0
      X2 = 0.D0
      XN = 0.D0
      c_missing = 42.D0
      DO N = 1, NPTS
      IF (X (N) .NE.XMSG) THEN
        XN = XN + 1.0D0
        X1 = X1 + X (N)
        X2 = X2 + X (N) * X (N)
      ELSE
        c_missing = c_missing + 1.D0
      ENDIF
      ENDDO

      NPTUSED = XN

      IF ( (XN - 1.D0) .GT.0.D0) THEN
        X3 = X1 * X1 / XN
!
! prevent possible roundoff
!
        XVAR = DMAX1 ( (X2 - X3) / (XN - 1.D0), 0.D0)
        XSD = SQRT (XVAR)
        XMEAN = X1 / XN
      ELSEIF ( (XN - 1.D0) .EQ.0.D0) THEN
        XMEAN = X1
        XVAR = 0.0D0
        XSD = 0.0D0
      ELSE
!
! error code for all msg values
!
        IER = 2
      ENDIF

      RETURN
      END SUBROUTINE DSTAT2
! ----------------------------------------------------------------------
      SUBROUTINE DSTAT2T (X, NPTS, XMSG, XMEANT, XVART, XSDT, NPTUSED,  &
      WORK, PTRIM, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XMEANT
      DOUBLEPRECISION XVART
      DOUBLEPRECISION XSDT
      DOUBLEPRECISION PTRIM

! this routine will calculate trimmed estimates of the first
! .   two moments of the vector x containing missing data.

! input arguments:
! .   x        - input vector
! .   npts     - length of x  (npts.ge.5 if pctrim.gt.0.)
! .   xmsg     - missing code: if there are no msg values
! .                            set xmsg to some value which will
! .                            not be encountered.
! .   work     - work vector of length npts
! .              upon return work contains x in ascending order
! .   ptrim    - portion of the series to be trimmed (0.<= pctrim <1.)
! .              if (pctrim.eq.0.) then
! .                  calculate whole series mean and st dev
! .              else
! .                  a minimum of two pts will be trimmed
! .                  a maximum of npts-2 will be trimmed
! .              endif
! .
! .              note: if u want to be very conservative and only
! .                    trim the min and max values then set
! .                    ptrim to some very small value (e.g. 1.e-36)
! .
! .              ex: npts=45, ptrim=0.50 then
! .                  0.5*45=22 pts will be trimmed (11 on each end)

! output arguments:
! .   xmeant   - mean of trimmed series x
! .   xvart    - sample trimmed variance
! .   xsdt     - sqrt( xvart ): st deviation of trimmed series
! .   nptused  - no. of points used to calculate the trimmed estimates
! .   ier      - if (ier.ne.0) an error has occurred

! note : uncalculated quantities are set to xmsg (nptused set to zero)

      DOUBLEPRECISION X (1:NPTS), WORK (1:NPTS)

      IER = 0
      XMEANT = XMSG
      XVART = XMSG
      XSDT = XMSG
      NPTUSED = 0

      IF (NPTS.LT.1) THEN
        IER = 1
        RETURN
      ENDIF

      IF (PTRIM.EQ.0.D0) THEN
        CALL DSTAT2 (X, NPTS, XMSG, XMEANT, XVART, XSDT, NPTUSED, IER)
        RETURN
      ELSEIF (PTRIM.EQ.1.D0) THEN
!
! cannot trim whole series
!
        IER = 2
        RETURN
      ENDIF

! sort x into numerically ascending order
! .   first transfer the input array to the work array
! .   excluding missing pts

      MPTS = 0
      DO N = 1, NPTS
      IF (X (N) .NE.XMSG) THEN
        MPTS = MPTS + 1
        WORK (MPTS) = X (N)
      ENDIF
      ENDDO

      IF (MPTS.LT.1) THEN
!
! series contains all msg values
!
        IER = 3
        RETURN
      ELSEIF (MPTS.LT.NPTS) THEN
!
! fill the rest of work with msg code
!
        DO M = MPTS + 1, NPTS
        WORK (M) = XMSG
        ENDDO
      ENDIF
!
! sort series in ascending order
!
      CALL DSORTU (WORK, MPTS)
!
! no. of values trimmed from one side
!
      MTRIM = MAX0 (INT (DBLE (MPTS) * 0.5D0 * PTRIM), 1)


      IF ( (MPTS - 2 * MTRIM) .GT.0) THEN
        CALL DSTAT2 (WORK (MTRIM + 1), MPTS - 2 * MTRIM, XMSG, XMEANT,  &
        XVART, XSDT, NPTUSED, IER)
      ELSE
!
! not enough trimmed values
!
        IER = 4
      ENDIF

      RETURN
      END SUBROUTINE DSTAT2T
! ----------------------------------------------------------------------
      SUBROUTINE DSTAT4 (X, NPTS, XMSG, XMEAN, XVAR, XSD, XSKEW, XKURT, &
      NPTUSED, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XMEAN
      DOUBLEPRECISION XVAR
      DOUBLEPRECISION XSD
      DOUBLEPRECISION XSKEW
      DOUBLEPRECISION XKURT
      DOUBLEPRECISION X1
      DOUBLEPRECISION X2
      DOUBLEPRECISION XN
      DOUBLEPRECISION X3
      DOUBLEPRECISION X4

! this routine will calculate estimates of the first four moments
! .   of the vector x containing missing data.

! input arguments:
! .   x        - input vector
! .   npts     - length of x
! .   xmsg     - missing code: if there are no msg values
! .                            set xmsg to some value which will
! .                            not be encountered.

! output arguments:
! .   xmean    - mean of x (first moment)
! .   xvar     - sample variance of x (second moment)
! .   xsd      - sqrt( xvar )
! .   xskew    - coef. of skewness (third moment)
! .              departure from symmetry. if skew>0 [skew<0] the
! .              distribution trails off to the right [left] .
! .   xkurt    - coef. of kurtosis (fourth moment)
! .              the normal distribution has a kurtosis of 3 . this
! .              value is subtracted from the calculated kurtosis.
! .              thus neg values are possible and the returned value
! .              is kurtosis relative to the normal distribution.
! .              if (kurt-3) > 0 [<0] it is usually more sharply peaked
! .              [flatter] than the normal distribution. (leptokurtic
! .              and platykurtic , respectively.) e.g. , a rectangular
!.               function has a kurtosis of 1.8 or 1.8-3 = -1.2 relative
!.               to the normal distribution.
! .   nptused  - no. of points used to calculate the estimates
! .   ier      - if (ier.ne.0) an error has occurred

! note :
! .   uncalculated quantities are set to xmsg (nptused set to zero)

      DOUBLEPRECISION X (1:NPTS)

      XMEAN = XMSG
      XVAR = XMSG
      XSD = XMSG
      XSKEW = XMSG
      XKURT = XMSG
      NPTUSED = 0


      IF (NPTS.LT.1) THEN
        IER = 1
        RETURN
      ENDIF

      IER = 0
      X1 = 0.D0
      X2 = 0.D0
      XN = 0.D0
      DO N = 1, NPTS
      IF (X (N) .NE.XMSG) THEN
        XN = XN + 1.0D0
        X1 = X1 + X (N)
        X2 = X2 + X (N) * X (N)
      ENDIF
      ENDDO

      NPTUSED = XN

      IF ( (XN - 1.D0) .GT.0.D0) THEN
        X3 = X1 * X1 / XN
!
! prevent possible roundoff
!
        XVAR = DMAX1 ( (X2 - X3) / (XN - 1.D0), 0.D0)
        XSD = SQRT (XVAR)
        XMEAN = X1 / XN

        IF (XVAR.GT.0.D0) THEN
          X1 = 0.D0
          X2 = 0.D0
          X3 = 0.D0
          X4 = 0.D0
          DO N = 1, NPTS
          IF (X (N) .NE.XMSG) THEN
            X4 = X (N) - XMEAN
            X1 = X4**3
            X2 = X2 + X1
            X3 = X3 + X1 * X4
          ENDIF
          ENDDO
          XSKEW = (X2 / SQRT (XVAR) **3) / XN
          XKURT = (X3 / (XVAR * XVAR) ) / XN - 3.D0
        ENDIF
      ELSEIF ( (XN - 1.D0) .EQ.0.D0) THEN
        XMEAN = X1
        XVAR = 0.0D0
        XSD = 0.0D0
      ELSE
!
! error code for all msg values
!
        IER = 2
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE DSTAT4
! ----------------------------------------------------------------------
      SUBROUTINE DMEDMRNG (X, WORK, NPTS, XMSG, XMEDIAN, XMRANGE,       &
      XRANGE, NPTUSED, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XMEDIAN
      DOUBLEPRECISION XMRANGE
      DOUBLEPRECISION XRANGE

! this routine will find the median, range and mid-range
! .   of the input  vector x containing missing data.

! arguments :
! .   x        - input vector (series)
! .   work     - work vector of length npts. upon output
! .              work(1) thru work(nptused) are in numerically ascending
! .              order. work(nptused+1) thru work(npts) will be filled
! .              the missing value code.
! .   npts     - length of x
! .   xmsg     - missing code: if no msg values set to some number
! .                            which will not be encountered.
! .   xmedian  - median value of vector x
! .   xmrange  - midrange
! .   xrange   - range of x
! .   nptused  - no. of points used in computing various quantities
! .   ier      - if (ier.ne.0) an error has occurred

! note :
! .   uncalculated quantities are set to xmsg (nptused set to zero)

      DOUBLEPRECISION X (1:NPTS), WORK (1:NPTS)

      XMEDIAN = XMSG
      XMRANGE = XMSG
      XRANGE = XMSG
      NPTUSED = 0

      IF (NPTS.LT.1) THEN
        IER = 1
        RETURN
      ENDIF
      IER = 0

! sort x into numerically ascending order
! .   first transfer the input array to the work array
! .   excluding missing pts

      MPTS = 0
      DO N = 1, NPTS
      IF (X (N) .NE.XMSG) THEN
        MPTS = MPTS + 1
        WORK (MPTS) = X (N)
      ENDIF
      ENDDO

      IF (MPTS.LT.1) THEN
!
! series contains all msg values
!
        IER = 2
        RETURN
      ENDIF

      CALL DSORTU (WORK, MPTS)

      IF (MPTS.LT.NPTS) THEN
!
! fill the rest of work with msg code
!
        DO N = MPTS + 1, NPTS
        WORK (N) = XMSG
        ENDDO
      ENDIF

      XMRANGE = 0.5D0 * (WORK (MPTS) + WORK (1) )
      XRANGE = WORK (MPTS) - WORK (1)

      IF (DMOD (DBLE (MPTS), 2.D0) .EQ.0.D0) THEN
        XMEDIAN = 0.5D0 * (WORK (MPTS / 2) + WORK (MPTS / 2 + 1) )
      ELSE
        XMEDIAN = WORK ( (MPTS + 1) / 2)
      ENDIF

!
! for nomenclature compatibility
!
      NPTUSED = MPTS

      RETURN
      END SUBROUTINE DMEDMRNG
! ----------------------------------------------------------------------
      SUBROUTINE DESAUTO (X, NPTS, XMSG, XMEAN, XVAR, MXLAG, ACV, ACR,  &
      IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XMEAN
      DOUBLEPRECISION XVAR
      DOUBLEPRECISION XSD
      DOUBLEPRECISION X1
      DOUBLEPRECISION XN

! this routine will estimate the autocovariances and autocorrelations
! .   for the input vector x with missing data. the normalizing factor
! .   in the denominator (xvar) is kept constant.

! arguments :
! .   x        - input vector
! .   npts     - length of vector x
! .   xmsg     - missing code: if no msg values set to some number
! .                            which will not be encountered.
! .   xmean    - mean of the input vector x
! .   xvar     - sample variance  of the input vector x
! .              note: if xmean and/or xvar = xmsg this routine will
! .                    calculate the quantities and return them
! .   mxlag    - max lag to be estimated  [0 <= mxlag <= npts ]
! .   acv      - sample autocovariances  [vector of length .ge.
! .                                      (mxlag+1) ]
! .   acr      - sample autocorrelations [vector of length .ge.
! .                                       (mxlag+1) ]

! .              acv(0) = xvar (lag 0  )  : acr(0)   = 1.00 ( lag 0  )
! .              acv(1) =      (lag 1  )  : acr(1)   =      ( lag 1  )
! .                         .                              .
! .                         .                              .
! .              acv(mxlag) = (lag mxlag): acr(mxlag) = ( lag mxlag)

! .   ier      - if (ier.ne.0) an error has occurred

! note :  uncalculated quantities are set to xmsg

      DOUBLEPRECISION X (1:NPTS), ACV (0:MXLAG), ACR (0:MXLAG)

      IER = 0
      IF (NPTS.LT.2) IER = 1
      IF (MXLAG.LT.0.OR.MXLAG.GT.NPTS) IER = 2
      IF (XVAR.LE.0.D0.AND.XVAR.NE.XMSG) IER = 3

      DO LAG = 0, MAX0 (0, MIN0 (MXLAG, NPTS) )
      ACV (LAG) = XMSG
      ACR (LAG) = XMSG
      ENDDO
      IF (IER.NE.0) RETURN


      IF (XMEAN.EQ.XMSG.OR.XVAR.EQ.XMSG) THEN
        CALL DSTAT2 (X, NPTS, XMSG, XMEAN, XVAR, XSD, NPTUSED, JER)
        IF (JER.NE.0) THEN
          IER = - JER
        ELSEIF (XVAR.EQ.0.D0) THEN
          IER = - 5
        ENDIF
      ENDIF
      IF (IER.NE.0) RETURN

      ACV (0) = XVAR
      ACR (0) = 1.00D0

!
! i hope somebody wants more than this
!
      IF (MXLAG.EQ.0) RETURN

      DO LAG = 1, MXLAG
      X1 = 0.D0
      XN = 0.D0
      DO N = 1, NPTS - LAG
      IF (X (N) .NE.XMSG.AND.X (N + LAG) .NE.XMSG) THEN
        X1 = X1 + ( (X (N + LAG) - XMEAN) * (X (N) - XMEAN) )
        XN = XN + 1.D0
      ENDIF
      ENDDO
      IF (XN.GE.2.D0) THEN
        ACV (LAG) = X1 / (XN - 1.D0)
        ACR (LAG) = ACV (LAG) / XVAR
      ELSE
        ACV (LAG) = XMSG
        ACR (LAG) = XMSG
      ENDIF
      ENDDO

      RETURN
      END SUBROUTINE DESAUTO
! ----------------------------------------------------------------------
      SUBROUTINE DESCROS (X, Y, NPTS, XMSG, YMSG, XMEAN, YMEAN, XSD,    &
      YSD, MXLAG, CCV, CCR, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION YMSG
      DOUBLEPRECISION XMEAN
      DOUBLEPRECISION YMEAN
      DOUBLEPRECISION XSD
      DOUBLEPRECISION YSD
      DOUBLEPRECISION XVAR
      DOUBLEPRECISION YVAR
      DOUBLEPRECISION XSDYSD
      DOUBLEPRECISION XYN
      DOUBLEPRECISION XY1


! this routine will estimate the crosscovariances and crosscorrelations
! .   for the input series x and y. missing data is allowed.

! .   to get the (+) and (-) ccv and ccr two calls must be made. the
! .   first with the calling sequence above and the second with the
! .   x and y series reversed (also the other arguments) . note
! .   that  the second ccv(0) and ccr(0) are redundant.

! arguments :
! .   x,y      - input vectors with which the calculations are to take
! .              place.
! .   npts     - length of vectors x and y
! .   xmsg,ymsg- missing code: if no msg values set to some number
! .                            which will not be encountered.
! .              xmsg will be used to fill missing values
! .   xmean    - mean of the x vector (if xmean=xmsg: program will
! .              calculate)
! .   ymean    - mean of the y vector (if ymean=ymsg: program will
! .              calculate)
! .   xsd      - st. dev. of x vector (if xsd  =xmsg: program will
! .              calculate)
! .   ysd      - st. dev. of y vector (if ysd  =ymsg: program will
! .              calculate)
! .   mxlag    - max lag crosscovariance/correlation to be estimated
! .   ccv      - vector of length .ge. (mxlag+1) containing sample
! .              cross cov
! .   ccr      - vector of length .ge. (mxlag+1) containing sample
! .              cross cor

! .                ccv(0) = xyvar (lag 0  )  : ccr(0) =  ( lag 0   )
! .                ccv(1) =       (lag 1  )  : ccr(1) =  ( lag 1   )
! .                           .                               .
! .                           .                               .
! .                ccv(mxlag) =  (lag mxlag): ccr(mxlag) = ( lag mxlag)

! .   ier      - if (ier.ne.0) an error has occurred

! note : uncalculated quantities are set to xmsg

      DOUBLEPRECISION X (1:NPTS), Y (1:NPTS), CCV (0:MXLAG), CCR (0:    &
      MXLAG)

      IER = 0
      IF (NPTS.LT.2) IER = 1
      IF (MXLAG.LT.0.OR.MXLAG.GT.NPTS) IER = 2
      IF ( (XSD.NE.XMSG.AND.XSD.LE.0.D0) .OR. (                         &
      YSD.NE.YMSG.AND.YSD.LE.0.D0) ) IER = 3
      DO LAG = 0, MAX0 (0, MIN0 (MXLAG, NPTS) )
      CCV (LAG) = XMSG
      CCR (LAG) = XMSG
      ENDDO
      IF (IER.NE.0) RETURN

      IF (XMEAN.EQ.XMSG.OR.XSD.EQ.XMSG) THEN
        CALL DSTAT2 (X, NPTS, XMSG, XMEAN, XVAR, XSD, NPTUSED, JER)
        IF (JER.NE.0) THEN
          IER = - JER
        ELSEIF (XSD.EQ.0.D0) THEN
!
! x must be a series of constant values
!
          IER = - 5
        ELSEIF (NPTUSED.EQ.0) THEN
!
! x must be a series of missing values
!
          IER = - 6
        ENDIF
        IF (IER.NE.0) RETURN
      ENDIF

      IF (YMEAN.EQ.YMSG.OR.YSD.EQ.YMSG) THEN
        CALL DSTAT2 (Y, NPTS, YMSG, YMEAN, YVAR, YSD, NPTUSED, JER)
        IF (JER.NE.0) THEN
          IER = - (JER + 100)
        ELSEIF (YSD.EQ.0.D0) THEN
!
! y must be a series of constant values
!
          IER = - 105
        ELSEIF (NPTUSED.EQ.0) THEN
!
! y must be a series of missing values
!
          IER = - 106
        ENDIF
        IF (IER.NE.0) RETURN
      ENDIF

      XSDYSD = 1.D0 / (XSD * YSD)

      DO LAG = 0, MXLAG
      XYN = 0.D0
      XY1 = 0.D0
      DO N = 1, NPTS - LAG
      IF (X (N) .NE.XMSG.AND.Y (N + LAG) .NE.YMSG) THEN
        XY1 = XY1 + ( (Y (N + LAG) - YMEAN) * (X (N) - XMEAN) )
        XYN = XYN + 1.D0
      ENDIF
      ENDDO
      IF (XYN.GE.2.D0) THEN
        CCV (LAG) = XY1 / (XYN - 1.D0)
        CCR (LAG) = CCV (LAG) * XSDYSD
      ELSE
        CCV (LAG) = XMSG
        CCR (LAG) = XMSG
      ENDIF
      ENDDO

      RETURN
      END SUBROUTINE DESCROS
! ----------------------------------------------------------------------
      SUBROUTINE DRMVAVE (X, NPTS, XMSG, XMEAN, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XMEAN

! this routine will remove the mean of vector x from each pt.
! .   missing pts are excepted.

! arguments :
! .   x        - input vector
! .   npts     - length of x
! .   xmsg     - missing code: if there are no msg values
! .                            set xmsg to some value which will
! .                            not be encountered.
! .   xmean    - mean of the input vector (x)
! .   ier      - if (ier.ne.0) an error has occurred

      DOUBLEPRECISION X (1:NPTS)

      IF (NPTS.LT.1) THEN
        IER = 1
        RETURN
      ELSE
        IER = 0
      ENDIF

      DO N = 1, NPTS
      IF (X (N) .NE.XMSG) X (N) = X (N) - XMEAN
      ENDDO

      RETURN
      END SUBROUTINE DRMVAVE
! ----------------------------------------------------------
      SUBROUTINE DSTNDX (X, NPTS, XMSG, XMEAN, XSD, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XMEAN
      DOUBLEPRECISION XSD

! standardize the vector x

! arguments :
! .   x        - vector to be standardized
! .   npts     - length of vector x
! .   xmsg     - missing code: if there are no msg values
! .                            set xmsg to some value which will
! .                            not be encountered.
! .   xmean    - mean of the input vector
! .   xsd      - standard deviation of the series
! .   ier      - if (ier.ne.0) an error has ocurred

      DOUBLEPRECISION X (1:NPTS)

      IER = 0
      IF (NPTS.LT.1) THEN
        IER = 1
      ELSEIF (XSD.LE.0.D0) THEN
        IER = 2
      ENDIF
      IF (IER.NE.0) RETURN

      DO N = 1, NPTS
      IF (X (N) .NE.XMSG) X (N) = (X (N) - XMEAN) / XSD
      ENDDO

      RETURN
      END SUBROUTINE DSTNDX
! ----------------------------------------------------------------------
      SUBROUTINE DNMDATA (X, NPTS, XMSG, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XMAX

! normalize the vector x by the absolute value of its largest element

! arguments :
! .   x        - series (vector) to be normalized
! .   npts     - length of x
! .   xmsg     - missing code: if there are no msg values
! .                            set xmsg to some value which will
! .                            not be encountered.
! .   ier      - if (ier.ne.0) an error has occurred

      DOUBLEPRECISION X (1:NPTS)

      IER = 0
      IF (NPTS.LT.1) THEN
        IER = 1
        RETURN
      ENDIF

! find the first non-missing pt

      DO N = 1, NPTS
      IF (X (N) .NE.XMSG) THEN
        NSTRT = N
        GOTO 20
      ENDIF
      ENDDO

!
! must be all msg values
!
      IER = 1
      RETURN

   20 XMAX = ABS (X (NSTRT) )
      DO N = NSTRT, NPTS
      IF (X (N) .NE.XMSG) XMAX = DMAX1 (ABS (X (N) ), XMAX)
      ENDDO

      IF (XMAX.NE.0.D0) THEN
        DO N = NSTRT, NPTS
        IF (X (N) .NE.XMSG) X (N) = X (N) / XMAX
        ENDDO
      ELSE
!
! must be all zeros
!
        IER = 2
      ENDIF

      RETURN
      END SUBROUTINE DNMDATA
! ----------------------------------------------------------------------
      SUBROUTINE DERRCHK (X, NPTS, XMEAN, XSD, XMSG, XBAND, LOCERR,     &
      NPTER, IER)
      DOUBLEPRECISION XMEAN
      DOUBLEPRECISION XSD
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XBAND
      DOUBLEPRECISION XWIDTH

! this routine will check a series of data for potential gross errors
! .   (i.e., wild points). the series will be scanned and the locations
! .   of possibly bad points are stored. this routine
! .   assumes that the series x is normally distributed.

! input arguments:
! .   x        - series (vector) to be checked
! .   npts     - length of x
! .   xmean    - series mean (... or median is the user chooses)
! .   xsd      - series standard deviation
! .   xmsg     - missing code: if there are no msg values
! .                            set xmsg to some value which will
! .                            not be encountered.
! .   xband    - standard deviation width to be allowed
! .              (e.g., xband=3.0; check values beyond 99%)

! output arguments:
! .   locerr   - integer vector which will contain the subscript
! .              location of potential gross errors. for example: if
! .              series x has two pts more than xband from xmean and
! .              these pts are located at x(17) and x(204) then
! .              locerr(1) = 17 and locerr(2) = 204.
! .   npter    - no. of pts in locerr. in the above example it would be
! .              two. if a value of zero is returned no gross errors
! .              have been detected.
! .   ier      - if (ier.ne.0) an error has occurred

      DOUBLEPRECISION X (1:NPTS)
      INTEGER LOCERR (1:NPTS)

      IER = 0
      IF (NPTS.LT.1) THEN
        IER = 1
      ELSEIF (XSD.EQ.0.D0) THEN
        IER = 2
      ENDIF
      IF (IER.NE.0) RETURN

      NN = 0
      NPTER = 0
      XWIDTH = XBAND * XSD
      DO N = 1, NPTS
      IF (X (N) .NE.XMSG.AND.ABS (X (N) - XMEAN) .GT.XWIDTH) THEN
!
! must be a potential gross error
!
        NPTER = NPTER + 1
        LOCERR (NPTER) = N
      ENDIF
      ENDDO

      RETURN
      END SUBROUTINE DERRCHK
! ----------------------------------------------------------------------
      SUBROUTINE DDOTPR (X, Y, NPTS, XMSG, YMSG, ANS, NPTUSED, IER)
      DOUBLEPRECISION X
      DOUBLEPRECISION Y
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION YMSG
      DOUBLEPRECISION ANS

! compute the dot product for series containing msg values

! arguments :
! .   x,y       - input series
! .   npts      - length of x and y
! .   xmsg,ymsg - missing code: if there are no msg values
! .                             set xmsg,ymsg to some value which will
! .                             not be encountered.
! .   ans       - answer (dot product)
! .   nptused   - number of points used to compute the dot product
! .   ier       - if (ier.ne.0) an error has occurred

      DIMENSION X (1:NPTS), Y (1:NPTS)

      IER = 0
      IF (NPTS.LT.1) THEN
        IER = 1
        RETURN
      ENDIF

      ANS = 0.D0
      NPTUSED = 0
      DO N = 1, NPTS
      IF (X (N) .NE.XMSG.AND.Y (N) .NE.YMSG) THEN
        ANS = ANS + X (N) * Y (N)
        NPTUSED = NPTUSED+1
      ENDIF
      ENDDO

      RETURN
      END SUBROUTINE DDOTPR
! ----------------------------------------------------------------------
      SUBROUTINE DVCVMNS (X, NROW, NCOL, NRT, NCS, XMSG, VCM, LVCM, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XN
      DOUBLEPRECISION XCA
      DOUBLEPRECISION XCB
      DOUBLEPRECISION XCACB

! this routine will calculate the variance-covariance matrix (vcm)
! .   of the array x containing missing data. obviously if x does
! .   contain missing data then vcm is only an approximation.

! note : conventional storage is utilized for vcm

!     x        - input data array  ( unchanged on output)
!     nrow,ncol- exact dimensions of x in calling routine
!     nrt,ncs  - dimension of sub-matrix which contains the data
!                (nrt@nrow : ncs@ncol)
!     xmsg     - missing data code (if none set to some no. not
!                encountered)
!     vcm      - var-cov matrix
!     lvcm     - not used in this routine
!     ier      - error code (if ier=-1 then vcm contains missing entry)

      DOUBLEPRECISION X (1:NROW, 1:NCOL), VCM (1:NCOL, 1:NCOL)

      IER = 0

! calculate the var-cov between columns (stations)

      DO NCA = 1, NCS
      DO NCB = NCA, NCS
      XN = 0.D0
      XCA = 0.D0
      XCB = 0.D0
      XCACB = 0.D0
      DO I = 1, NRT
      IF (X (I, NCA) .NE.XMSG.AND.X (I, NCB) .NE.XMSG) THEN
        XN = XN + 1.D0
        XCA = XCA + X (I, NCA)
        XCB = XCB + X (I, NCB)
        XCACB = XCACB + X (I, NCA) * X (I, NCB)
      ENDIF
!
! end "nrt"
!
      ENDDO
      IF (XN.GE.2.D0) THEN
        VCM (NCA, NCB) = (XCACB - (XCA * XCB) / (XN) ) / (XN - 1.D0)
        VCM (NCB, NCA) = VCM (NCA, NCB)
      ELSE
        IER = - 1
        VCM (NCA, NCB) = XMSG
        VCM (NCB, NCA) = XMSG
      ENDIF
!
! end "nrt"
!
      ENDDO
!
! end "ncs"
!
      ENDDO

      RETURN
      END SUBROUTINE DVCVMNS
! ----------------------------------------------------------------------
      SUBROUTINE DVCMSSM (X, NROW, NCOL, NRT, NCS, XMSG, VCM, LVCM, IER)
      DOUBLEPRECISION X
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION VCM
      DOUBLEPRECISION XN
      DOUBLEPRECISION XCA
      DOUBLEPRECISION XCB
      DOUBLEPRECISION XCACB

! this routine will calculate the variance-couariance matrix (vcm)
! .   of the array x containing missing data. obviously if x does
! .   contain missing data then vcm is only an approximation.

! note : symmetric storage mode is utilized for vcm to save space.

! input:
!     x        - input data array  ( unchanged on output)
!     nrow,ncol- exact dimensions of x in calling routine
!     nrt,ncs  - dimension of sub-matrix which contains the data
!                (nrt <= nrow : ncs <= ncol)
!     xmsg     - missing data code (if none set to some no. not
!                encountered)
! output:
!     vcm      - var-cov matrix
!     lvcm     - length of vcm
!     ier      - error code (if ier=-1 then vcm contains missing entry)

      DIMENSION X (1:NROW, 1:NCOL), VCM ( * )

      IER = 0
      IF (NROW.LT.1.OR.NCOL.LT.1) IER = IER + 1
      IF (NRT.LT.1.OR.NCS.LT.1) IER = IER + 10
      IF (IER.NE.0) RETURN

      LVCM = NCS * (NCS + 1) / 2

! calculate the var-cov between columns (stations)

      NN = 0
      DO NCA = 1, NCS
      DO NCB = 1, NCA
      XN = 0.D0
      XCA = 0.D0
      XCB = 0.D0
      XCACB = 0.D0
      NN = NN + 1
      DO I = 1, NRT
      IF (X (I, NCA) .NE.XMSG.AND.X (I, NCB) .NE.XMSG) THEN
        XN = XN + 1.D0
        XCA = XCA + X (I, NCA)
        XCB = XCB + X (I, NCB)
        XCACB = XCACB + (X (I, NCA) ) * (X (I, NCB) )
      ENDIF
      ENDDO
      IF (XN.GE.2.D0) THEN
        VCM (NN) = (XCACB - (XCA * XCB) / XN) / (XN - 1.D0)
      ELSEIF (XN.EQ.1.) THEN
        VCM (NN) = (XCACB - (XCA * XCB) / XN)
      ELSE
        IER = - 1
        VCM (NN) = XMSG
      ENDIF
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DVCMSSM
!
! ----------------------------------------------------------------------
      SUBROUTINE DCORMNS (X, NROW, NCOL, NRT, NCS, XMSG, CRM, LCRM, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XN
      DOUBLEPRECISION XCA
      DOUBLEPRECISION XCB
      DOUBLEPRECISION XCA2
      DOUBLEPRECISION XCB2
      DOUBLEPRECISION XCACB
      DOUBLEPRECISION XVARA
      DOUBLEPRECISION XVARB

! this routine will calculate the correlation matrix   (crm)
! .   of the array x containing missing data. obviously if x does
! .   contain missing data then crm is only an approximation.

! note : conventional storage

!     x        - input data array  ( unchanged on output)
!     nrow,ncol- exact dimensions of x in calling routine
!     nrt,ncs  - dimension of sub-matrix which contains the data
!                (nrt@nrow : ncs@ncol)
!     xmsg     - missing data code (if none set to some no. not
!                encountered)
!     crm      - correlation matrix  [full]
!     lcrm     - dummy [not used in this routine]
!     ier      - error code (if ier=-1 then crm contains missing entry)

      DOUBLEPRECISION X (1:NROW, 1:NCOL), CRM (1:NCOL, 1:NCOL)

      IER = 0
      LCRM = NCS * (NCS + 1) / 2

! calculate the var-cov between columns (stations)
! .   then standardize

      DO NCA = 1, NCS
      DO NCB = 1, NCA
      XN = 0.D0
      XCA = 0.D0
      XCB = 0.D0
      XCA2 = 0.D0
      XCB2 = 0.D0
      XCACB = 0.D0
      DO I = 1, NRT
      IF (X (I, NCA) .NE.XMSG.AND.X (I, NCB) .NE.XMSG) THEN
        XN = XN + 1.D0
        XCA = XCA + X (I, NCA)
        XCB = XCB + X (I, NCB)
        XCA2 = XCA2 + X (I, NCA) * X (I, NCA)
        XCB2 = XCB2 + X (I, NCB) * X (I, NCB)
        XCACB = XCACB + X (I, NCA) * X (I, NCB)
      ENDIF
      ENDDO
      IF (XN.GE.2.D0) THEN
        XVARA = (XCA2 - ( (XCA * XCA) / (XN) ) ) / (XN - 1.D0)
        XVARB = (XCB2 - ( (XCB * XCB) / (XN) ) ) / (XN - 1.D0)
        IF (XVARA.GT.0.D0.AND.XVARB.GT.0.D0) THEN
          CRM (NCB, NCA) = (XCACB - ( (XCA * XCB) / (XN) ) ) / (XN -    &
          1.D0)
          CRM (NCB, NCA) = CRM (NCB, NCA) / (SQRT (XVARA) * SQRT (XVARB)&
          )
          CRM (NCA, NCB) = CRM (NCB, NCA)
        ELSE
          IER = - 1
          CRM (NCB, NCA) = XMSG
          CRM (NCA, NCB) = XMSG
        ENDIF
      ELSE
        IER = - 1
        CRM (NCB, NCA) = XMSG
        CRM (NCA, NCB) = XMSG
      ENDIF
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DCORMNS
! ----------------------------------------------------------------------
      SUBROUTINE DCRMSSM (X, NROW, NCOL, NRT, NCS, XMSG, CRM, LCRM, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION XN
      DOUBLEPRECISION XCA
      DOUBLEPRECISION XCB
      DOUBLEPRECISION XCA2
      DOUBLEPRECISION XCB2
      DOUBLEPRECISION XCACB
      DOUBLEPRECISION XVARA
      DOUBLEPRECISION XVARB

! this routine will calculate the correlation matrix   (crm)
! .   of the array x containing missing data. obviously if x does
! .   contain missing data then crm is only an approximation.

! note : symmetric storage mode is utilized for crm to save space.

! input:
!     x        - input data array  ( unchanged on output)
!     nrow,ncol- exact dimensions of x in calling routine
!     nrt,ncs  - dimension of sub-matrix which contains the data
!                (nrt <= nrow : ncs <= ncol)
!     xmsg     - missing data code (if none set to some no. not
!                encountered)
! output:
!     crm      - correlation matrix
!     lcrm     - length of crm
!     ier      - error code (if ier=-1 then crm contains missing entry)

      DOUBLEPRECISION X (1:NROW, 1:NCOL), CRM ( * )

      IER = 0
      IF (NROW.LT.1.OR.NCOL.LT.1) IER = IER + 1
      IF (NRT.LT.1.OR.NCS.LT.1) IER = IER + 10
      IF (IER.NE.0) RETURN

      LCRM = NCS * (NCS + 1) / 2

! calculate the var-cov between columns (stations)

      NN = 0
      DO NCA = 1, NCS
      DO NCB = 1, NCA
      XN = 0.D0
      XCA = 0.D0
      XCB = 0.D0
      XCA2 = 0.D0
      XCB2 = 0.D0
      XCACB = 0.D0
      NN = NN + 1
      DO I = 1, NRT
      IF (X (I, NCA) .NE.XMSG.AND.X (I, NCB) .NE.XMSG) THEN
        XN = XN + 1.D0
        XCA = XCA + X (I, NCA)
        XCB = XCB + X (I, NCB)
        XCA2 = XCA2 + X (I, NCA) * X (I, NCA)
        XCB2 = XCB2 + X (I, NCB) * X (I, NCB)
        XCACB = XCACB + X (I, NCA) * X (I, NCB)
      ENDIF
      ENDDO
      IF (XN.GE.2.D0) THEN
        XVARA = (XCA2 - ( (XCA * XCA) / (XN) ) ) / (XN - 1.D0)
        XVARB = (XCB2 - ( (XCB * XCB) / (XN) ) ) / (XN - 1.D0)
        IF (XVARA.GT.0.D0.AND.XVARB.GT.0.D0) THEN
          CRM (NN) = (XCACB - ( (XCA * XCB) / (XN) ) ) / (XN - 1.D0)
          CRM (NN) = CRM (NN) / (SQRT (XVARA) * SQRT (XVARB) )
        ELSE
          IER = - 1
          CRM (NN) = XMSG
        ENDIF
      ELSEIF (XN.EQ.1.D0) THEN
        XVARA = (XCA2 - ( (XCA * XCA) / (XN) ) )
        XVARB = (XCB2 - ( (XCB * XCB) / (XN) ) )
        IF (XVARA.GT.0.D0.AND.XVARB.GT.0.D0) THEN
          CRM (NN) = (XCACB - ( (XCA * XCB) / (XN) ) )
          CRM (NN) = CRM (NN) / (SQRT (XVARA) * SQRT (XVARB) )
        ELSE
          IER = - 1
          CRM (NN) = XMSG
        ENDIF
      ELSE
        IER = - 1
        CRM (NN) = XMSG
      ENDIF
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DCRMSSM

! ----------------------------------------------------------------------
      SUBROUTINE DCCORC (Y, X, NROW, NCOL, NRT, NCS, XMSG, CORC, IER)
      DOUBLEPRECISION XMSG
      DOUBLEPRECISION YMEAN
      DOUBLEPRECISION YVAR
      DOUBLEPRECISION YSD
      DOUBLEPRECISION XMEAN
      DOUBLEPRECISION XVAR
      DOUBLEPRECISION XSD
      DOUBLEPRECISION XN
      DOUBLEPRECISION XY
      DOUBLEPRECISION COVAR

! this routine will calculate the correlation coef between a
! .   vector (station) y and a matrix x containing observations from
! .   many different stations

!     y        - vector (station data) with which data matrix x is to
!                be correlated
!     x        - input data array  ( unchanged on output)
!     nrow,ncol- exact dimensions of x in calling routine
!     nrt,ncs  - dimension of sub-matrix which contains the data
!                (nrt@nrow : ncs@ncol)
!     xmsg     - missing data code (if none set to some no. not
!                encountered)
!     corc     - correlation coef. vector
!                corc(1) contain correlation coef between y and x(nrt,1)
!                corc(2) contain correlation coef between y and
!                    x(nrt,2) etc.
!     ier      - error code (if ier=-1 then corc contains missing entry)

      DOUBLEPRECISION X (1:NROW, 1:NCOL), Y ( * ), CORC ( * )

! calculate the mean and var of y

      CALL DSTAT2 (Y, NRT, XMSG, YMEAN, YVAR, YSD, NPTUSY, IER)
      IF (IER.NE.0) RETURN
      IF (YSD.EQ.0.D0) THEN
        DO NC = 1, NCS
        CORC (NC) = XMSG
        ENDDO
        IER = 201
        RETURN
      ENDIF

! calculate the cross correlation coef
! .   calculate the mean and variance of column nc in the data array

      DO NC = 1, NCS
      CORC (NC) = XMSG
      CALL DSTAT2 (X (1, NC), NRT, XMSG, XMEAN, XVAR, XSD, NPTUSX, JER)
      IF (JER.NE.0) THEN
        IER = 205
        GOTO 40
      ENDIF
      IF (XSD.NE.0.D0) THEN
        XN = 0.D0
        XY = 0.D0
        DO N = 1, NRT
        IF (Y (N) .NE.XMSG.AND.X (N, NC) .NE.XMSG) THEN
          XN = XN + 1.D0
          XY = XY + (Y (N) - YMEAN) * (X (N, NC) - XMEAN)
        ENDIF
        ENDDO
        IF (XN.GT.2.D0) THEN
          COVAR = XY / (XN - 1.D0)
          CORC (NC) = COVAR / (XSD * YSD)
        ENDIF
      ENDIF
   40 CONTINUE
      ENDDO

      RETURN
      END SUBROUTINE DCCORC
!------------------------------------------------------------------
      SUBROUTINE DSSMIOX (X, NCS)

! output a symmetric storage mode matrix [x]

! note: format statement have to be changed

! .   x(1)
! .   x(2)  x(3)
! .   x(4)  x(5)  x(6)
! .   etc.

      DOUBLEPRECISION X (1: * )

      DATA IPR / 6 /

      NCEND = 0
      DO NC = 1, NCS
      NCSTRT = NCEND+1
      NCEND = NCSTRT + NC - 1
      WRITE (IPR, FMT = '(2x,i5,10(1x,f12.3),/, (7x,(10(1x,f12.3))))')  &
      NC, (X (N) , N = NCSTRT, NCEND)
      ENDDO

      RETURN
      END SUBROUTINE DSSMIOX
! ------------------------------------------------------------
      DOUBLEPRECISION FUNCTION DEPSMACH (IPR)

!
! =1 print else no print
!
      INTEGER IPR
      DOUBLEPRECISION EPS

      EPS = 1.0D0
    1 EPS = 0.5D0 * EPS
      IF ( (1.0D0 + EPS) .GT.1.0D0) GOTO 1
      DEPSMACH = 2.0D0 * EPS
      IF (IPR.EQ.1) WRITE ( * , FMT = '(''DEPSMACH: eps='',e15.7)')     &
      DEPSMACH

      RETURN
      END FUNCTION DEPSMACH

! ----------------------------------------------------------------------
      SUBROUTINE DSORTU (A, LA)

! sort type "real"

      INTEGER LA
      DOUBLEPRECISION A (LA)

      INTEGER IU (21), IL (21), I, M, J, K, IJ, L
      DOUBLEPRECISION T, TT, R

      M = 1
      I = 1
      J = LA
      R = .375D0
      IF (LA.LE.0) RETURN
   10 IF (I.EQ.J) GOTO 55
   15 IF (R.GT..5898437D0) GOTO 20
      R = R + 3.90625D-2
      GOTO 25
   20 R = R - .21875D0
   25 K = I
!                          select a central element of the
!                          array and save it in location t
      IJ = I + (J - I) * R
      T = A (IJ)
!                          if first element of array is greater
!                          than t, interchange with t
      IF (A (I) .LE.T) GOTO 30
      A (IJ) = A (I)
      A (I) = T
      T = A (IJ)
   30 L = J
!                           if last element of array is less than
!                           t, interchange with t
      IF (A (J) .GE.T) GOTO 40
      A (IJ) = A (J)
      A (J) = T
      T = A (IJ)
!                          if first element of array is greater
!                          than t, interchange with t
      IF (A (I) .LE.T) GOTO 40
      A (IJ) = A (I)
      A (I) = T
      T = A (IJ)
      GOTO 40
   35 IF (A (L) .EQ.A (K) ) GOTO 40
      TT = A (L)
      A (L) = A (K)
      A (K) = TT
!                           find an element in the second half of
!                           the array which is smaller than t
   40 L = L - 1
      IF (A (L) .GT.T) GOTO 40
!                           find an element in the first half of
!                           the array which is greater than t
   45 K = K + 1
      IF (A (K) .LT.T) GOTO 45
!                           interchange these elements
      IF (K.LE.L) GOTO 35
!                           save upper and lower subscripts of
!                           the array yet to be sorted
      IF (L - I.LE.J - K) GOTO 50
      IL (M) = I
      IU (M) = L
      I = K
      M = M + 1
      GOTO 60
   50 IL (M) = K
      IU (M) = J
      J = L
      M = M + 1
      GOTO 60
!                           begin again on another portion of
!                           the unsorted array
   55 M = M - 1
      IF (M.EQ.0) RETURN
      I = IL (M)
      J = IU (M)
   60 IF (J - I.GE.11) GOTO 25
      IF (I.EQ.1) GOTO 10
      I = I - 1
   65 I = I + 1
      IF (I.EQ.J) GOTO 55
      T = A (I + 1)
      IF (A (I) .LE.T) GOTO 65
      K = I
   70 A (K + 1) = A (K)
      K = K - 1
      IF (T.LT.A (K) ) GOTO 70
      A (K + 1) = T
      GOTO 65
      END SUBROUTINE DSORTU
! ----------------------------------------------------------------------
      SUBROUTINE DISORTU (A, LA)

! sort type "integer"

      INTEGER LA
      INTEGER A (LA)

      INTEGER IU (21), IL (21), I, M, J, K, IJ, L
      INTEGER T, TT
      DOUBLEPRECISION R

      M = 1
      I = 1
      J = LA
      R = .375D0
      IF (LA.LE.0) RETURN
   10 IF (I.EQ.J) GOTO 55
   15 IF (R.GT..5898437D0) GOTO 20
      R = R + 3.90625D-2
      GOTO 25
   20 R = R - .21875D0
   25 K = I
!                                  select a central element of the
!                                  array and save it in location t
      IJ = I + (J - I) * R
      T = A (IJ)
!                                  if first element of array is greater
!                                  than t, interchange with t
      IF (A (I) .LE.T) GOTO 30
      A (IJ) = A (I)
      A (I) = T
      T = A (IJ)
   30 L = J
!                                  if last element of array is less than
!                                  t, interchange with t
      IF (A (J) .GE.T) GOTO 40
      A (IJ) = A (J)
      A (J) = T
      T = A (IJ)
!                                  if first element of array is greater
!                                  than t, interchange with t
      IF (A (I) .LE.T) GOTO 40
      A (IJ) = A (I)
      A (I) = T
      T = A (IJ)
      GOTO 40
   35 IF (A (L) .EQ.A (K) ) GOTO 40
      TT = A (L)
      A (L) = A (K)
      A (K) = TT
!                                  find an element in the second half of
!                                  the array which is smaller than t
   40 L = L - 1
      IF (A (L) .GT.T) GOTO 40
!                                  find an element in the first half of
!                                  the array which is greater than t
   45 K = K + 1
      IF (A (K) .LT.T) GOTO 45
!                                  interchange these elements
      IF (K.LE.L) GOTO 35
!                                  save upper and lower subscripts of
!                                  the array yet to be sorted
      IF (L - I.LE.J - K) GOTO 50
      IL (M) = I
      IU (M) = L
      I = K
      M = M + 1
      GOTO 60
   50 IL (M) = K
      IU (M) = J
      J = L
      M = M + 1
      GOTO 60
!                                  begin again on another portion of
!                                  the unsorted array
   55 M = M - 1
      IF (M.EQ.0) RETURN
      I = IL (M)
      J = IU (M)
   60 IF (J - I.GE.11) GOTO 25
      IF (I.EQ.1) GOTO 10
      I = I - 1
   65 I = I + 1
      IF (I.EQ.J) GOTO 55
      T = A (I + 1)
      IF (A (I) .LE.T) GOTO 65
      K = I
   70 A (K + 1) = A (K)
      K = K - 1
      IF (T.LT.A (K) ) GOTO 70
      A (K + 1) = T
      GOTO 65
      END SUBROUTINE DISORTU

      SUBROUTINE COLLAPSEXY (X, Y, NPTS, XMSG, YMSG, XX, YY, MPTS)
      IMPLICIT NONE
      INTEGER NPTS, MPTS
      DOUBLEPRECISION X (NPTS), Y (NPTS)
      DOUBLEPRECISION XX (NPTS), YY (NPTS)
      DOUBLEPRECISION XMSG, YMSG

      INTEGER N

      MPTS = 0
      DO N = 1, NPTS
      IF (X (N) .NE.XMSG.AND.Y (N) .NE.YMSG) THEN
        MPTS = MPTS + 1
        XX (MPTS) = X (N)
        YY (MPTS) = Y (N)
!              PRINT *,'n, mpts, xx(mpts), yy(mpts)=',N,MPTS,XX(MPTS),
!     +          YY(MPTS)
      ENDIF
      ENDDO

      RETURN
      END SUBROUTINE COLLAPSEXY
