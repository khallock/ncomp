! ----------------------------------------------------------------------
      SUBROUTINE stat2 (x, npts, xmsg, xmean, xvar, xsd, nptused, ier) 
                                                                        
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
                                                                        
      INTEGER npts, nptused, ier 
      REAL x (1:npts), xmsg, xmean, xvar, xsd 
      DOUBLEPRECISION x1, x2, x3, temp 
                                                                        
      xmean = xmsg 
      xvar = xmsg 
      xsd = xmsg 
      nptused = 0 
                                                                        
      IF (npts.lt.1) then 
        ier = 1 
        RETURN 
      ENDIF 
                                                                        
      ier = 0 
      x1 = 0.d0 
      x2 = 0.d0 
      xn = 0. 
      DO n = 1, npts 
      IF (x (n) .ne.xmsg) then 
        xn = xn + 1.0 
        x1 = x1 + dble (x (n) ) 
        x2 = x2 + dprod (x (n), x (n) ) 
      ENDIF 
      enddo 
                                                                        
      nptused = xn 
                                                                        
      IF ( (xn - 1.) .gt.0.) then 
        x3 = x1 * x1 / dble (xn) 
        temp = (x2 - x3) / dble (xn - 1.) 
                                             ! prevent possible roundoff
        xvar = sngl (dmax1 (temp, 0.d0) ) 
        xsd = sqrt (xvar) 
        xmean = sngl (x1 / dble (xn) ) 
      ELSEIF ( (xn - 1.) .eq.0.) then 
        xmean = sngl (x1) 
        xvar = 0.0 
        xsd = 0.0 
      ELSE 
                    ! error code for all msg values                     
        ier = 2 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE stat2                          
! ----------------------------------------------------------------------
      SUBROUTINE stat2t (x, npts, xmsg, xmeant, xvart, xsdt, nptused,   &
      work, ptrim, ier)                                                 
                                                                        
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
                                                                        
      REAL x (1:npts), work (1:npts) 
                                                                        
      ier = 0 
      xmeant = xmsg 
      xvart = xmsg 
      xsdt = xmsg 
      nptused = 0 
                                                                        
      IF (npts.lt.1) then 
        ier = 1 
        RETURN 
      ENDIF 
                                                                        
      IF (ptrim.eq.0.) then 
        CALL stat2 (x, npts, xmsg, xmeant, xvart, xsdt, nptused, ier) 
        RETURN 
      ELSEIF (ptrim.eq.1.) then 
                      ! cannot trim whole series                        
        ier = 2 
        RETURN 
      ENDIF 
                                                                        
! sort x into numerically ascending order                               
! .   first transfer the input array to the work array                  
! .   excluding missing pts                                             
                                                                        
      mpts = 0 
      DO n = 1, npts 
      IF (x (n) .ne.xmsg) then 
        mpts = mpts + 1 
        work (mpts) = x (n) 
      ENDIF 
      enddo 
                                                                        
      IF (mpts.lt.1) then 
                                    ! series contains all msg values    
        ier = 3 
        RETURN 
      ELSEIF (mpts.lt.npts) then 
                                    ! fill the rest of work with msg cod
        DO m = mpts + 1, npts 
        work (m) = xmsg 
        enddo 
      ENDIF 
                                                                        
                                    ! sort series in ascending order    
      CALL sortu (work, mpts) 
                                                                        
                                                       ! no. of values t
      mtrim = max0 (int (float (mpts) * 0.5 * ptrim), 1) 
                                                       ! from one side  
                                                                        
      IF ( (mpts - 2 * mtrim) .gt.0) then 
        CALL stat2 (work (mtrim + 1), mpts - 2 * mtrim, xmsg, xmeant,   &
        xvart, xsdt, nptused, ier)                                      
      ELSE 
                                    ! not enough trimmed values         
        ier = 4 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE stat2t                         
! ----------------------------------------------------------------------
      SUBROUTINE stat4 (x, npts, xmsg, xmean, xvar, xsd, xskew, xkurt,  &
      nptused, ier)                                                     
                                                                        
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
! .              [flatter] than the normal distribution. (leptokurtic an
! .              platykurtic , respectively.) e.g. , a rectangular funct
! .              a kurtosis of 1.8 or 1.8-3 = -1.2 relative to the norma
! .              distribution.                                          
! .   nptused  - no. of points used to calculate the estimates          
! .   ier      - if (ier.ne.0) an error has occurred                    
                                                                        
! note :                                                                
! .   uncalculated quantities are set to xmsg (nptused set to zero)     
                                                                        
      REAL x (1:npts) 
      DOUBLEPRECISION x1, x2, x3, x4, dxvar 
                                                                        
      xmean = xmsg 
      xvar = xmsg 
      xsd = xmsg 
      xskew = xmsg 
      xkurt = xmsg 
      nptused = 0 
                                                                        
                                                                        
      IF (npts.lt.1) then 
        ier = 1 
        RETURN 
      ENDIF 
                                                                        
      ier = 0 
      x1 = 0.d0 
      x2 = 0.d0 
      xn = 0. 
      DO n = 1, npts 
      IF (x (n) .ne.xmsg) then 
        xn = xn + 1.0 
        x1 = x1 + dble (x (n) ) 
        x2 = x2 + dprod (x (n), x (n) ) 
      ENDIF 
      enddo 
                                                                        
      nptused = xn 
                                                                        
      IF ( (xn - 1.) .gt.0.) then 
        x3 = x1 * x1 / dble (xn) 
                                                     !prevent possible r
        dxvar = dmax1 ( (x2 - x3) / dble (xn - 1.), 0.d0) 
        xvar = sngl (dxvar) 
        xsd = sqrt (xvar) 
        xmean = sngl (x1 / dble (xn) ) 
                                                                        
        IF (dxvar.gt.0.d0) then 
          x1 = 0.d0 
          x2 = 0.d0 
          x3 = 0.d0 
          x4 = 0.d0 
          DO n = 1, npts 
          IF (x (n) .ne.xmsg) then 
            x4 = dble (x (n) - xmean) 
            x1 = x4**3 
            x2 = x2 + x1 
            x3 = x3 + x1 * x4 
          ENDIF 
          enddo 
          xskew = sngl ( (x2 / dsqrt (dxvar) **3) / dble (xn) ) 
          xkurt = sngl ( (x3 / (dxvar * dxvar) ) / dble (xn) ) - 3. 
        ENDIF 
      ELSEIF ( (xn - 1.) .eq.0.) then 
        xmean = sngl (x1) 
        xvar = 0.0 
        xsd = 0.0 
      ELSE 
                    ! error code for all msg values                     
        ier = 2 
        RETURN 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE stat4                          
! ----------------------------------------------------------------------
      SUBROUTINE medmrng (x, work, npts, xmsg, xmedian, xmrange, xrange,&
      nptused, ier)                                                     
                                                                        
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
                                                                        
      REAL x (1:npts), work (1:npts) 
                                                                        
      xmedian = xmsg 
      xmrange = xmsg 
      xrange = xmsg 
      nptused = 0 
                                                                        
      IF (npts.lt.1) then 
        ier = 1 
        RETURN 
      ENDIF 
      ier = 0 
                                                                        
! sort x into numerically ascending order                               
! .   first transfer the input array to the work array                  
! .   excluding missing pts                                             
                                                                        
      mpts = 0 
      DO n = 1, npts 
      IF (x (n) .ne.xmsg) then 
        mpts = mpts + 1 
        work (mpts) = x (n) 
      ENDIF 
      enddo 
                                                                        
      IF (mpts.lt.1) then 
                                      ! series contains all msg values  
        ier = 2 
        RETURN 
      ENDIF 
                                                                        
      CALL sortu (work, mpts) 
                                                                        
      IF (mpts.lt.npts) then 
                                      ! fill the rest of work with msg c
        DO n = mpts + 1, npts 
        work (n) = xmsg 
        enddo 
      ENDIF 
                                                                        
      xmrange = 0.5 * (work (mpts) + work (1) ) 
      xrange = work (mpts) - work (1) 
                                                                        
      IF (amod (float (mpts), 2.) .eq.0.) then 
        xmedian = 0.5 * (work (mpts / 2) + work (mpts / 2 + 1) ) 
      ELSE 
        xmedian = work ( (mpts + 1) / 2) 
      ENDIF 
                                                                        
                              ! for nomenclature compatibility          
      nptused = mpts 
                                                                        
      RETURN 
      END SUBROUTINE medmrng                        
! ----------------------------------------------------------------------
      SUBROUTINE esauto (x, npts, xmsg, xmean, xvar, mxlag, acv, acr,   &
      ier)                                                              
                                                                        
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
! .   acv      - sample autocovariances  [vector of length .ge. (mxlag+1
! .   acr      - sample autocorrelations [vector of length .ge. (mxlag+1
                                                                        
! .              acv(0)     = xvar (lag 0  )  : acr(0)       = 1.00 ( la
! .              acv(1)     =      (lag 1  )  : acr(1)       =      ( la
! .                         .                              .            
! .                         .                              .            
! .              acv(mxlag) =      (lag mxlag): acr(mxlag) =      ( lag 
                                                                        
! .   ier      - if (ier.ne.0) an error has occurred                    
                                                                        
! note :  uncalculated quantities are set to xmsg                       
                                                                        
      REAL x (1:npts), acv (0:mxlag), acr (0:mxlag) 
      DOUBLEPRECISION x1 
                                                                        
      ier = 0 
      IF (npts.lt.2) ier = 1 
      IF (mxlag.lt.0.or.mxlag.gt.npts) ier = 2 
      IF (xvar.le.0..and.xvar.ne.xmsg) ier = 3 
      IF (ier.ne.0) then 
        DO lag = 0, max0 (0, min0 (mxlag, npts) ) 
        acv (lag) = xmsg 
        acr (lag) = xmsg 
        enddo 
        RETURN 
      ENDIF 
                                                                        
      IF (xmean.eq.xmsg.or.xvar.eq.xmsg) then 
        CALL stat2 (x, npts, xmsg, xmean, xvar, xsd, nptused, jer) 
        IF (jer.ne.0) then 
          ier = - jer 
          RETURN 
        ELSEIF (xvar.eq.0.) then 
          ier = - 5 
          RETURN 
        ENDIF 
      ENDIF 
                                                                        
      acv (0) = xvar 
      acr (0) = 1.00 
                                                                        
                                   ! i hope somebody wants more than thi
      IF (mxlag.eq.0) return 
                                                                        
      DO lag = 1, mxlag 
      x1 = 0.d0 
      xn = 0. 
      DO n = 1, npts - lag 
      IF (x (n) .ne.xmsg.and.x (n + lag) .ne.xmsg) then 
        x1 = x1 + dprod ( (x (n + lag) - xmean), (x (n) - xmean) ) 
        xn = xn + 1. 
      ENDIF 
      enddo 
      IF (xn.ge.2.) then 
        acv (lag) = sngl (x1 / dble (xn - 1.) ) 
        acr (lag) = acv (lag) / xvar 
      ELSE 
        acv (lag) = xmsg 
        acr (lag) = xmsg 
      ENDIF 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE esauto                         
! ----------------------------------------------------------------------
      SUBROUTINE escros (x, y, npts, xmsg, ymsg, xmean, ymean, xsd, ysd,&
      mxlag, ccv, ccr, ier)                                             
                                                                        
                                                                        
! this routine will estimate the crosscovariances and crosscorrelations 
! .   for the input series x and y. missing data is allowed.            
                                                                        
! .   to get the (+) and (-) ccv and ccr two calls must be made. the    
! .   first with the calling sequence above and the second with the     
! .   x and y series reversed (also the other arguments) . note         
! .   that  the second ccv(0) and ccr(0) are redundant.                 
                                                                        
! arguments :                                                           
! .   x,y      - input vectors with which the calculations are to take p
! .   npts     - length of vectors x and y                              
! .   xmsg,ymsg- missing code: if no msg values set to some number      
! .                            which will not be encountered.           
! .              xmsg will be used to fill missing values               
! .   xmean    - mean of the x vector (if xmean=xmsg: program will calcu
! .   ymean    - mean of the y vector (if ymean=ymsg: program will calcu
! .   xsd      - st. dev. of x vector (if xsd  =xmsg: program will calcu
! .   ysd      - st. dev. of y vector (if ysd  =ymsg: program will calcu
! .   mxlag    - max lag crosscovariance/correlation to be estimated    
! .   ccv      - vector of length .ge. (mxlag+1) containing sample cross
! .   ccr      - vector of length .ge. (mxlag+1) containing sample cross
                                                                        
! .                ccv(0)     = xyvar (lag 0  )  : ccr(0)     =  ( lag 0
! .                ccv(1)     =       (lag 1  )  : ccr(1)     =  ( lag 1
! .                           .                               .         
! .                           .                               .         
! .                ccv(mxlag) =       (lag mxlag): ccr(mxlag) =  ( lag m
                                                                        
! .   ier      - if (ier.ne.0) an error has occurred                    
                                                                        
! note : uncalculated quantities are set to xmsg                        
                                                                        
      REAL x (1:npts), y (1:npts), ccv (0:mxlag), ccr (0:mxlag) 
      DOUBLEPRECISION xy1 
                                                                        
      ier = 0 
      IF (npts.lt.2) ier = 1 
      IF (mxlag.lt.0.or.mxlag.gt.npts) ier = 2 
      IF ( (xsd.ne.xmsg.and.xsd.le.0.) .or. (ysd.ne.ymsg.and.ysd.le.0.) &
      ) ier = 3                                                         
      DO lag = 0, max0 (0, min0 (mxlag, npts) ) 
      ccv (lag) = xmsg 
      ccr (lag) = xmsg 
      enddo 
      IF (ier.ne.0) return 
                                                                        
      IF (xmean.eq.xmsg.or.xsd.eq.xmsg) then 
        CALL stat2 (x, npts, xmsg, xmean, xvar, xsd, nptused, jer) 
        IF (jer.ne.0) then 
          ier = - jer 
        ELSEIF (xsd.eq.0.) then 
                                   ! x must be a series of constant valu
          ier = - 5 
        ELSEIF (nptused.eq.0) then 
                                   ! x must be a series of missing value
          ier = - 6 
        ENDIF 
        IF (ier.ne.0) return 
      ENDIF 
                                                                        
      IF (ymean.eq.ymsg.or.ysd.eq.ymsg) then 
        CALL stat2 (y, npts, ymsg, ymean, yvar, ysd, nptused, jer) 
        IF (jer.ne.0) then 
          ier = - (jer + 100) 
        ELSEIF (ysd.eq.0.) then 
                                   ! y must be a series of constant valu
          ier = - 105 
        ELSEIF (nptused.eq.0) then 
                                   ! y must be a series of missing value
          ier = - 106 
        ENDIF 
        IF (ier.ne.0) return 
      ENDIF 
                                                                        
      xsdysd = sngl (1.d0 / dprod (xsd, ysd) ) 
                                                                        
      DO lag = 0, mxlag 
      xyn = 0. 
      xy1 = 0.d0 
      DO n = 1, npts - lag 
      IF (x (n) .ne.xmsg.and.y (n + lag) .ne.ymsg) then 
        xy1 = xy1 + dprod ( (y (n + lag) - ymean), (x (n) - xmean) ) 
        xyn = xyn + 1. 
      ENDIF 
      enddo 
      IF (xyn.ge.2.) then 
        ccv (lag) = sngl (xy1 / dble (xyn - 1.) ) 
        ccr (lag) = ccv (lag) * xsdysd 
      ELSE 
        ccv (lag) = xmsg 
        ccr (lag) = xmsg 
      ENDIF 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE escros                         
! ----------------------------------------------------------------------
      SUBROUTINE rmvave (x, npts, xmsg, xmean, ier) 
                                                                        
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
                                                                        
      REAL x (1:npts) 
                                                                        
      IF (npts.lt.1) then 
        ier = 1 
        RETURN 
      ELSE 
        ier = 0 
      ENDIF 
                                                                        
      DO n = 1, npts 
      IF (x (n) .ne.xmsg) x (n) = x (n) - xmean 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE rmvave                         
! ----------------------------------------------------------            
      SUBROUTINE stndx (x, npts, xmsg, xmean, xsd, ier) 
                                                                        
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
                                                                        
      REAL x (1:npts) 
                                                                        
      ier = 0 
      IF (npts.lt.1) then 
        ier = 1 
      ELSEIF (xsd.le.0.) then 
        ier = 2 
      ENDIF 
      IF (ier.ne.0) return 
                                                                        
      DO n = 1, npts 
      IF (x (n) .ne.xmsg) x (n) = (x (n) - xmean) / xsd 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE stndx                          
! ----------------------------------------------------------------------
      SUBROUTINE nmdata (x, npts, xmsg, ier) 
                                                                        
! normalize the vector x by the absolute value of its largest element   
                                                                        
! arguments :                                                           
! .   x        - series (vector) to be normalized                       
! .   npts     - length of x                                            
! .   xmsg     - missing code: if there are no msg values               
! .                            set xmsg to some value which will        
! .                            not be encountered.                      
! .   ier      - if (ier.ne.0) an error has occurred                    
                                                                        
      REAL x (1:npts) 
                                                                        
      ier = 0 
      IF (npts.lt.1) then 
        ier = 1 
        RETURN 
      ENDIF 
                                                                        
! find the first non-missing pt                                         
                                                                        
      DO n = 1, npts 
      IF (x (n) .ne.xmsg) then 
        nstrt = n 
        GOTO 20 
      ENDIF 
      enddo 
                                                                        
                   ! must be all msg values                             
      ier = 1 
      RETURN 
                                                                        
   20 xmax = abs (x (nstrt) ) 
      DO n = nstrt, npts 
      IF (x (n) .ne.xmsg) xmax = amax1 (abs (x (n) ), xmax) 
      enddo 
                                                                        
      IF (xmax.ne.0.) then 
        DO n = nstrt, npts 
        IF (x (n) .ne.xmsg) x (n) = x (n) / xmax 
        enddo 
      ELSE 
                  ! must be all zeros                                   
        ier = 2 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE nmdata                         
! ----------------------------------------------------------------------
      SUBROUTINE errchk (x, npts, xmean, xsd, xmsg, xband, locerr,      &
      npter, ier)                                                       
                                                                        
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
! .   locerr   - integer vector which will contain the subscript locatio
! .              potential gross errors. for example: if series x has   
! .              two pts more than xband from xmean and these pts are   
! .              located at x(17) and x(204) then locerr(1) = 17 and    
! .              locerr(2) = 204.                                       
! .   npter    - no. of pts in locerr. in the above example it would be 
! .              two. if a value of zero is returned no gross errors    
! .              have been detected.                                    
! .   ier      - if (ier.ne.0) an error has occurred                    
                                                                        
      REAL x (1:npts) 
      INTEGER locerr (1:npts) 
                                                                        
      ier = 0 
      IF (npts.lt.1) then 
        ier = 1 
      ELSEIF (xsd.eq.0.) then 
        ier = 2 
      ENDIF 
      IF (ier.ne.0) return 
                                                                        
      nn = 0 
      npter = 0 
      xwidth = xband * xsd 
      DO n = 1, npts 
      IF (x (n) .ne.xmsg.and.abs (x (n) - xmean) .gt.xwidth) then 
                                       ! must be a potential gross error
        npter = npter + 1 
        locerr (npter) = n 
      ENDIF 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE errchk                         
! ----------------------------------------------------------------------
      SUBROUTINE dotpr (x, y, npts, xmsg, ymsg, ans, nptused, ier) 
                                                                        
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
                                                                        
      DIMENSION x (1:npts), y (1:npts) 
      DOUBLEPRECISION temp 
                                                                        
      ier = 0 
      IF (npts.lt.1) then 
        ier = 1 
        RETURN 
      ENDIF 
                                                                        
      temp = 0.d0 
      nptused = 0 
      DO n = 1, npts 
      IF (x (n) .ne.xmsg.and.y (n) .ne.ymsg) then 
        temp = temp + dprod (x (n), y (n) ) 
        nptused = nptused+1 
      ENDIF 
      enddo 
                                                                        
      ans = sngl (temp) 
                                                                        
      RETURN 
      END SUBROUTINE dotpr                          
! ----------------------------------------------------------------------
      SUBROUTINE vcvmns (x, nrow, ncol, nrt, ncs, xmsg, vcm, lvcm, ier) 
                                                                        
! this routine will calculate the variance-covariance matrix (vcm)      
! .   of the array x containing missing data. obviously if x does contai
! .   missing data then vcm is only an approximation.                   
                                                                        
! note : conventional storage is utilized for vcm                       
                                                                        
!     x        - input data array  ( unchanged on output)               
!     nrow,ncol- exact dimensions of x in calling routine               
!     nrt,ncs  - dimension of sub-matrix which contains the data        
!                (nrt@nrow : ncs@ncol)                                  
!     xmsg     - missing data code (if none set to some no. not encounte
!     vcm      - var-cov matrix                                         
!     lvcm     - not used in this routine                               
!     ier      - error code (if ier=-1 then vcm contains missing entry) 
                                                                        
      REAL x (1:nrow, 1:ncol), vcm (1:ncol, 1:ncol) 
      DOUBLEPRECISION xca, xcb, xcacb 
                                                                        
      ier = 0 
                                                                        
! calculate the var-cov between columns (stations)                      
                                                                        
      DO nca = 1, ncs 
      DO ncb = nca, ncs 
      xn = 0. 
      xca = 0. 
      xcb = 0. 
      xcacb = 0. 
      DO i = 1, nrt 
      IF (x (i, nca) .ne.xmsg.and.x (i, ncb) .ne.xmsg) then 
        xn = xn + 1. 
        xca = xca + dble (x (i, nca) ) 
        xcb = xcb + dble (x (i, ncb) ) 
        xcacb = xcacb + dprod (x (i, nca), x (i, ncb) ) 
      ENDIF 
                 ! end "nrt"                                            
      enddo 
      IF (xn.ge.2.) then 
        vcm (nca, ncb) = sngl ( (xcacb - (xca * xcb) / dble (xn) )      &
        / dble (xn - 1.) )                                              
        vcm (ncb, nca) = vcm (nca, ncb) 
      ELSE 
        ier = - 1 
        vcm (nca, ncb) = xmsg 
        vcm (ncb, nca) = xmsg 
      ENDIF 
                 ! end "ncs"                                            
      enddo 
                 ! end "ncs"                                            
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE vcvmns                         
! ----------------------------------------------------------------------
      SUBROUTINE vcmssm (x, nrow, ncol, nrt, ncs, xmsg, vcm, lvcm, ier) 
                                                                        
! this routine will calculate the variance-couariance matrix (vcm)      
! .   of the array x containing missing data. obviously if x does contai
! .   missing data then vcm is only an approximation.                   
                                                                        
! note : symmetric storage mode is utilized for vcm to save space.      
                                                                        
! input:                                                                
!     x        - input data array  ( unchanged on output)               
!     nrow,ncol- exact dimensions of x in calling routine               
!     nrt,ncs  - dimension of sub-matrix which contains the data        
!                (nrt <= nrow : ncs <= ncol)                            
!     xmsg     - missing data code (if none set to some no. not encounte
! output:                                                               
!     vcm      - var-cov matrix                                         
!     lvcm     - length of vcm                                          
!     ier      - error code (if ier=-1 then vcm contains missing entry) 
                                                                        
      DIMENSION x (1:nrow, 1:ncol), vcm ( * ) 
      DOUBLEPRECISION xca, xcb, xcacb 
                                                                        
      ier = 0 
      IF (nrow.lt.1.or.ncol.lt.1) ier = ier + 1 
      IF (nrt.lt.1.or.ncs.lt.1) ier = ier + 10 
      IF (ier.ne.0) return 
                                                                        
      lvcm = ncs * (ncs + 1) / 2 
                                                                        
! calculate the var-cov between columns (stations)                      
                                                                        
      nn = 0 
      DO nca = 1, ncs 
      DO ncb = 1, nca 
      xn = 0. 
      xca = 0.d0 
      xcb = 0.d0 
      xcacb = 0.d0 
      nn = nn + 1 
      DO i = 1, nrt 
      IF (x (i, nca) .ne.xmsg.and.x (i, ncb) .ne.xmsg) then 
        xn = xn + 1. 
        xca = xca + dble (x (i, nca) ) 
        xcb = xcb + dble (x (i, ncb) ) 
        xcacb = xcacb + dprod (x (i, nca), x (i, ncb) ) 
      ENDIF 
      enddo 
      IF (xn.ge.2.) then 
        vcm (nn) = sngl ( (xcacb - (xca * xcb) / dble (xn) ) / dble (xn &
        - 1.) )                                                         
      ELSEIF (xn.eq.1.) THEN 
        vcm (nn) = (xcacb - (xca * xcb) / xn) 
      ELSE 
        ier = - 1 
        vcm (nn) = xmsg 
      ENDIF 
      enddo 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE vcmssm                         
!                                                                       
! ----------------------------------------------------------------------
      SUBROUTINE cormns (x, nrow, ncol, nrt, ncs, xmsg, crm, lcrm, ier) 
                                                                        
! this routine will calculate the correlation matrix   (crm)            
! .   of the array x containing missing data. obviously if x does contai
! .   missing data then crm is only an approximation.                   
                                                                        
! note : conventional storage                                           
                                                                        
!     x        - input data array  ( unchanged on output)               
!     nrow,ncol- exact dimensions of x in calling routine               
!     nrt,ncs  - dimension of sub-matrix which contains the data        
!                (nrt@nrow : ncs@ncol)                                  
!     xmsg     - missing data code (if none set to some no. not encounte
!     crm      - correlation matrix  [full]                             
!     lcrm     - dummy [not used in this routine]                       
!     ier      - error code (if ier=-1 then crm contains missing entry) 
                                                                        
      REAL x (1:nrow, 1:ncol), crm (1:ncol, 1:ncol) 
      DOUBLEPRECISION xca, xcb, xca2, xcb2, xcacb, xvara, xvarb, temp 
                                                                        
      ier = 0 
      lcrm = ncs * (ncs + 1) / 2 
                                                                        
! calculate the var-cov between columns (stations)                      
! .   then standardize                                                  
                                                                        
      DO nca = 1, ncs 
      DO ncb = 1, nca 
      xn = 0. 
      xca = 0.d0 
      xcb = 0.d0 
      xca2 = 0.d0 
      xcb2 = 0.d0 
      xcacb = 0.d0 
      DO i = 1, nrt 
      IF (x (i, nca) .ne.xmsg.and.x (i, ncb) .ne.xmsg) then 
        xn = xn + 1. 
        xca = xca + dble (x (i, nca) ) 
        xcb = xcb + dble (x (i, ncb) ) 
        xca2 = xca2 + dprod (x (i, nca), x (i, nca) ) 
        xcb2 = xcb2 + dprod (x (i, ncb), x (i, ncb) ) 
        xcacb = xcacb + dprod (x (i, nca), x (i, ncb) ) 
      ENDIF 
      enddo 
      IF (xn.ge.2.) then 
        xvara = (xca2 - ( (xca * xca) / dble (xn) ) ) / dble (xn - 1.) 
        xvarb = (xcb2 - ( (xcb * xcb) / dble (xn) ) ) / dble (xn - 1.) 
        IF (xvara.gt.0..and.xvarb.gt.0.) then 
          temp = (xcacb - ( (xca * xcb) / dble (xn) ) ) / dble (xn - 1.) 
          crm (ncb, nca) = sngl (temp / (dsqrt (xvara) * dsqrt (xvarb) )&
          )                                                             
                                              ! symmetric               
          crm (nca, ncb) = crm (ncb, nca) 
        ELSE 
          ier = - 1 
          crm (ncb, nca) = xmsg 
          crm (nca, ncb) = xmsg 
        ENDIF 
      ELSE 
        ier = - 1 
        crm (ncb, nca) = xmsg 
        crm (nca, ncb) = xmsg 
      ENDIF 
      enddo 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE cormns                         
! ----------------------------------------------------------------------
      SUBROUTINE crmssm (x, nrow, ncol, nrt, ncs, xmsg, crm, lcrm, ier) 
                                                                        
! this routine will calculate the correlation matrix   (crm)            
! .   of the array x containing missing data. obviously if x does contai
! .   missing data then crm is only an approximation.                   
                                                                        
! note : symmetric storage mode is utilized for crm to save space.      
                                                                        
! input:                                                                
!     x        - input data array  ( unchanged on output)               
!     nrow,ncol- exact dimensions of x in calling routine               
!     nrt,ncs  - dimension of sub-matrix which contains the data        
!                (nrt <= nrow : ncs <= ncol)                            
!     xmsg     - missing data code (if none set to some no. not encounte
! output:                                                               
!     crm      - correlation matrix                                     
!     lcrm     - length of crm                                          
!     ier      - error code (if ier=-1 then crm contains missing entry) 
                                                                        
      REAL x (1:nrow, 1:ncol), crm ( * ) 
      DOUBLEPRECISION xca, xcb, xca2, xcb2, xcacb, xvara, xvarb, temp 
                                                                        
      ier = 0 
      IF (nrow.lt.1.or.ncol.lt.1) ier = ier + 1 
      IF (nrt.lt.1.or.ncs.lt.1) ier = ier + 10 
      IF (ier.ne.0) return 
                                                                        
      lcrm = ncs * (ncs + 1) / 2 
                                                                        
! calculate the var-cov between columns (stations)                      
                                                                        
      nn = 0 
      DO nca = 1, ncs 
      DO ncb = 1, nca 
      xn = 0. 
      xca = 0.d0 
      xcb = 0.d0 
      xca2 = 0.d0 
      xcb2 = 0.d0 
      xcacb = 0.d0 
      nn = nn + 1 
      DO i = 1, nrt 
      IF (x (i, nca) .ne.xmsg.and.x (i, ncb) .ne.xmsg) then 
        xn = xn + 1. 
        xca = xca + dble (x (i, nca) ) 
        xcb = xcb + dble (x (i, ncb) ) 
        xca2 = xca2 + dprod (x (i, nca), x (i, nca) ) 
        xcb2 = xcb2 + dprod (x (i, ncb), x (i, ncb) ) 
        xcacb = xcacb + dprod (x (i, nca), x (i, ncb) ) 
      ENDIF 
      enddo 
      IF (xn.ge.2.) then 
        xvara = (xca2 - ( (xca * xca) / dble (xn) ) ) / dble (xn - 1.) 
        xvarb = (xcb2 - ( (xcb * xcb) / dble (xn) ) ) / dble (xn - 1.) 
        IF (xvara.gt.0.d0.and.xvarb.gt.0.d0) then 
          temp = (xcacb - ( (xca * xcb) / dble (xn) ) ) / dble (xn - 1.) 
          crm (nn) = sngl (temp / (dsqrt (xvara) * dsqrt (xvarb) ) ) 
        ELSE 
          ier = - 1 
          crm (nn) = xmsg 
        ENDIF 
      ELSEIF (xn.eq.1.) then 
        xvara = (xca2 - ( (xca * xca) / (xn) ) ) 
        xvarb = (xcb2 - ( (xcb * xcb) / (xn) ) ) 
        IF (xvara.gt.0..and.xvarb.gt.0.) then 
          crm (nn) = (xcacb - ( (xca * xcb) / (xn) ) ) 
          crm (nn) = crm (nn) / (sqrt (xvara) * sqrt (xvarb) ) 
        ELSE 
          ier = - 1 
          crm (nn) = xmsg 
        ENDIF 
      ELSE 
        ier = - 1 
        crm (nn) = xmsg 
      ENDIF 
      enddo 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE crmssm                         
                                                                        
! ----------------------------------------------------------------------
      SUBROUTINE ccorc (y, x, nrow, ncol, nrt, ncs, xmsg, corc, ier) 
                                                                        
! this routine will calculate the correlation coef between a            
! .   vector (station) y and a matrix x containing observations from    
! .   many different stations                                           
                                                                        
!     y        - vector (station data) with which data matrix x is to   
!                be correlated                                          
!     x        - input data array  ( unchanged on output)               
!     nrow,ncol- exact dimensions of x in calling routine               
!     nrt,ncs  - dimension of sub-matrix which contains the data        
!                (nrt@nrow : ncs@ncol)                                  
!     xmsg     - missing data code (if none set to some no. not encounte
!     corc     - correlation coef. vector                               
!                corc(1) contain correlation coef between y and x(nrt,1)
!                corc(2) contain correlation coef between y and x(nrt,2)
!     ier      - error code (if ier=-1 then corc contains missing entry)
                                                                        
      REAL x (1:nrow, 1:ncol), y ( * ), corc ( * ) 
      DOUBLEPRECISION xy, temp 
                                                                        
! calculate the mean and var of y                                       
                                                                        
      CALL stat2 (y, nrt, xmsg, ymean, yvar, ysd, nptusy, ier) 
      IF (ier.ne.0) return 
      IF (ysd.eq.0.) then 
        DO nc = 1, ncs 
        corc (nc) = xmsg 
        enddo 
        ier = 201 
        RETURN 
      ENDIF 
                                                                        
! calculate the cross correlation coef                                  
! .   calculate the mean and variance of column nc in the data array    
                                                                        
      DO nc = 1, ncs 
      corc (nc) = xmsg 
      CALL stat2 (x (1, nc), nrt, xmsg, xmean, xvar, xsd, nptusx, jer) 
      IF (jer.ne.0) then 
        ier = 205 
        GOTO 40 
      ENDIF 
      IF (xsd.ne.0.) then 
        xn = 0. 
        xy = 0.d0 
        DO n = 1, nrt 
        IF (y (n) .ne.xmsg.and.x (n, nc) .ne.xmsg) then 
          xn = xn + 1. 
          xy = xy + dprod ( (y (n) - ymean), (x (n, nc) - xmean) ) 
        ENDIF 
        enddo 
        IF (xn.gt.2.) then 
          temp = xy / dble (xn - 1.) 
          corc (nc) = sngl (temp / dprod (xsd, ysd) ) 
        ENDIF 
      ENDIF 
   40 CONTINUE 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE ccorc                          
! ----------------------------------------------------------------------
      SUBROUTINE ssmiox (x, ncs) 
                                                                        
! output a symmetric storage mode matrix [x]                            
                                                                        
! note: format statement have to be changed                             
                                                                        
! .   x(1)                                                              
! .   x(2)  x(3)                                                        
! .   x(4)  x(5)  x(6)                                                  
! .   etc.                                                              
                                                                        
      REAL x (1: * ) 
                                                                        
      DATA ipr / 6 / 
                                                                        
      ncend = 0 
      DO nc = 1, ncs 
      ncstrt = ncend+1 
      ncend = ncstrt + nc - 1 
      WRITE (ipr, '(2x,i5,10(1x,f12.3),/, (7x,(10(1x,f12.3))))') nc,    &
      (x (n) , n = ncstrt, ncend)                                       
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE ssmiox                         
! ------------------------------------------------------------          
      REAL function epsmach (ipr) 
                                                                        
                     ! =1 print else no print                           
      INTEGER ipr 
      REAL eps 
                                                                        
      eps = 1.0 
    1 eps = 0.5 * eps 
      IF ( (1.0 + eps) .gt.1.0) goto 1 
      epsmach = 2.0 * eps 
      IF (ipr.eq.1) write ( * , " ('EPSMACH: eps=', e15.7) ") epsmach 
                                                                        
      RETURN 
      END FUNCTION epsmach                          
                                                                        
! ----------------------------------------------------------------------
      SUBROUTINE sortu (a, la) 
                                                                        
! sort type "real"                                                      
                                                                        
      INTEGER la 
      REAL a (la) 
                                                                        
      INTEGER iu (21), il (21), i, m, j, k, ij, l 
      REAL t, tt, r 
                                                                        
      m = 1 
      i = 1 
      j = la 
      r = .375 
      IF (la.le.0) return 
   10 IF (i.eq.j) goto 55 
   15 IF (r.gt..5898437) goto 20 
      r = r + 3.90625e-2 
      GOTO 25 
   20 r = r - .21875 
   25 k = i 
!                                  select a central element of the      
!                                  array and save it in location t      
      ij = i + (j - i) * r 
      t = a (ij) 
!                                  if first element of array is greater 
!                                  than t, interchange with t           
      IF (a (i) .le.t) goto 30 
      a (ij) = a (i) 
      a (i) = t 
      t = a (ij) 
   30 l = j 
!                                  if last element of array is less than
!                                  t, interchange with t                
      IF (a (j) .ge.t) goto 40 
      a (ij) = a (j) 
      a (j) = t 
      t = a (ij) 
!                                  if first element of array is greater 
!                                  than t, interchange with t           
      IF (a (i) .le.t) goto 40 
      a (ij) = a (i) 
      a (i) = t 
      t = a (ij) 
      GOTO 40 
   35 IF (a (l) .eq.a (k) ) goto 40 
      tt = a (l) 
      a (l) = a (k) 
      a (k) = tt 
!                                  find an element in the second half of
!                                  the array which is smaller than t    
   40 l = l - 1 
      IF (a (l) .gt.t) goto 40 
!                                  find an element in the first half of 
!                                  the array which is greater than t    
   45 k = k + 1 
      IF (a (k) .lt.t) goto 45 
!                                  interchange these elements           
      IF (k.le.l) goto 35 
!                                  save upper and lower subscripts of   
!                                  the array yet to be sorted           
      IF (l - i.le.j - k) goto 50 
      il (m) = i 
      iu (m) = l 
      i = k 
      m = m + 1 
      GOTO 60 
   50 il (m) = k 
      iu (m) = j 
      j = l 
      m = m + 1 
      GOTO 60 
!                                  begin again on another portion of    
!                                  the unsorted array                   
   55 m = m - 1 
      IF (m.eq.0) return 
      i = il (m) 
      j = iu (m) 
   60 IF (j - i.ge.11) goto 25 
      IF (i.eq.1) goto 10 
      i = i - 1 
   65 i = i + 1 
      IF (i.eq.j) goto 55 
      t = a (i + 1) 
      IF (a (i) .le.t) goto 65 
      k = i 
   70 a (k + 1) = a (k) 
      k = k - 1 
      IF (t.lt.a (k) ) goto 70 
      a (k + 1) = t 
      GOTO 65 
      END SUBROUTINE sortu                          
! ----------------------------------------------------------------------
      SUBROUTINE isortu (a, la) 
                                                                        
! sort type "integer"                                                   
                                                                        
      INTEGER la 
      INTEGER a (la) 
                                                                        
      INTEGER iu (21), il (21), i, m, j, k, ij, l 
      INTEGER t, tt 
      REAL r 
                                                                        
      m = 1 
      i = 1 
      j = la 
      r = .375 
      IF (la.le.0) return 
   10 IF (i.eq.j) goto 55 
   15 IF (r.gt..5898437) goto 20 
      r = r + 3.90625e-2 
      GOTO 25 
   20 r = r - .21875 
   25 k = i 
!                                  select a central element of the      
!                                  array and save it in location t      
      ij = i + (j - i) * r 
      t = a (ij) 
!                                  if first element of array is greater 
!                                  than t, interchange with t           
      IF (a (i) .le.t) goto 30 
      a (ij) = a (i) 
      a (i) = t 
      t = a (ij) 
   30 l = j 
!                                  if last element of array is less than
!                                  t, interchange with t                
      IF (a (j) .ge.t) goto 40 
      a (ij) = a (j) 
      a (j) = t 
      t = a (ij) 
!                                  if first element of array is greater 
!                                  than t, interchange with t           
      IF (a (i) .le.t) goto 40 
      a (ij) = a (i) 
      a (i) = t 
      t = a (ij) 
      GOTO 40 
   35 IF (a (l) .eq.a (k) ) goto 40 
      tt = a (l) 
      a (l) = a (k) 
      a (k) = tt 
!                                  find an element in the second half of
!                                  the array which is smaller than t    
   40 l = l - 1 
      IF (a (l) .gt.t) goto 40 
!                                  find an element in the first half of 
!                                  the array which is greater than t    
   45 k = k + 1 
      IF (a (k) .lt.t) goto 45 
!                                  interchange these elements           
      IF (k.le.l) goto 35 
!                                  save upper and lower subscripts of   
!                                  the array yet to be sorted           
      IF (l - i.le.j - k) goto 50 
      il (m) = i 
      iu (m) = l 
      i = k 
      m = m + 1 
      GOTO 60 
   50 il (m) = k 
      iu (m) = j 
      j = l 
      m = m + 1 
      GOTO 60 
!                                  begin again on another portion of    
!                                  the unsorted array                   
   55 m = m - 1 
      IF (m.eq.0) return 
      i = il (m) 
      j = iu (m) 
   60 IF (j - i.ge.11) goto 25 
      IF (i.eq.1) goto 10 
      i = i - 1 
   65 i = i + 1 
      IF (i.eq.j) goto 55 
      t = a (i + 1) 
      IF (a (i) .le.t) goto 65 
      k = i 
   70 a (k + 1) = a (k) 
      k = k - 1 
      IF (t.lt.a (k) ) goto 70 
      a (k + 1) = t 
      GOTO 65 
      END SUBROUTINE isortu                         
