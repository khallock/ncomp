! ---------------------------------------------------------             
      SUBROUTINE DDRVEOF (X, NROW, NCOL, NOBS, MSTA, XMSG, NEVAL, EVAL, &
      EVEC, PCVAR, TRACE, IOPT, JOPT, CSSM, LCSSM, WORK, LWORK, WEVAL,  &
      IWORK, LIWORK, IFAIL, LIFAIL, IER)                                
      DOUBLEPRECISION XMSG 
      DOUBLEPRECISION TRACE 
      DOUBLEPRECISION TOL 
      DOUBLEPRECISION DEPSMACH 
      DOUBLEPRECISION VLOW 
      DOUBLEPRECISION VUP 
                                                                        
! driver to calculate :                                                 
!*PL*ERROR* Comment line too long                                       
! .   the principal components (eignvalues and corresponding eigenvector
! .       of the data array x. x may contain missing observations.      
! .       if it has msg data, this will calculate a var-cov matrix      
! .       but it may not be positive definite.                          
                                                                        
! .       USES LAPACK/BLAS ROUTINES                                     
                                                                        
! nomenclature :                                                        
!*PL*ERROR* Comment line too long                                       
! .   x         - matrix containing the data.  it contains n observation
! .               for each of m stations or grid pts.                   
!*PL*ERROR* Comment line too long                                       
! .   nrow,ncol - exact row (observation) and column (station) dimension
! .               of x in the calling routine.                          
! .   nobs      - actual number of observations (nobs <= nrow)          
! .   msta      - actual number of stations     (msta <= ncol)          
! .   xmsg      - missing code (if no obs are missing set to some       
! .               number which will not be encountered)                 
! .   neval     - no. of eigenvalues and eigen vectors to be computed.  
!*PL*ERROR* Comment line too long                                       
! .               neval <= msta , if not msta eigenvalues and eigenvecto
! .               will be computed and ier = -2.                        
! .   eval      - vector containing the eigenvalues in DESCENDING order.
! .               eval must be at least neval in length.                
! .   evec      - an array which will contains the eigenvector info.    
! .               this must be dimensioned at least (ncol,neval) in the 
! .               calling routine. There is some normalization done     
! .               but I am not sure waht it is.                         
!*PL*ERROR* Comment line too long                                       
! .   pcvar     - contains the % variance associated with each eigenvalu
! .               this must be dimensioned at least neval in length.    
! .   trace     - trace of the variance - covariance matrix.in the      
! .               case of a var-covar matrix , the trace should         
! .               equal the sum of all possible eigenvalues.            
! .   iopt      - not used; kept for compatibility with old routine     
! .   jopt      - =  0 : use var-covar matrix in prncmp                 
! .                  1 : use correlation matrix in prncmp               
! .   cssm      - real vector to hold the covariance/correlation matrix 
! .   lcssm     - length of cssm =msta*(msta+1)/2                       
! .   work      - real vector                                           
! .   weval     - = msta                                                
! .   lwork     - = 8*msta                                              
! .   iwork     - integer vector                                        
! .   liwork    - = 5*msta                                              
! .   ifail     - integer vector                                        
! .   lifail    - = msta                                                
! .   ier       - error code                                            
                                                                        
! NOTE: upon return if one want to get loadings,                        
!       something like this should be tried.                            
!                                                                       
!           to calculate the coef or amplitude vector                   
!           associated with the n th eigenvector                        
!           mult the anomaly matrix (datam) times the                   
!           transpose of the eigenvector                                
!               meval = neval                                           
!               do n=1,meval                                            
!                  evwrk(n) = 0.                                        
!                 do nyr=1,nyrs                                         
!                    evyran(nyr,n) = 0.0                                
!                    evyrpv(nyr,n) = 0.0                                
!                   do m=1,msta                                         
!                      if (datam(nyr,m).ne.xmsg) then                   
!                          temp = evec(m,n)*datam(nyr,m)                
!                          evyran(nyr,n) = evyran(nyr,n) + temp         
!                          evyrpv(nyr,n) = evyrpv(nyr,n) + temp*temp    
!                          evwrk(n) = evwrk(n) + temp*temp              
!                      endif                                            
!                   enddo                                               
!                 enddo                                                 
!               .                                                       
!               enddo                                                   
                                                                        
                                                                        
      DOUBLEPRECISION X (1:NROW, 1:NCOL), EVAL ( * ), EVEC (1:NCOL,     &
      * )                                                               
      REAL PCVAR ( * ) 
      INTEGER(8) LCSSM 
      DOUBLEPRECISION CSSM (LCSSM), WORK (LWORK), WEVAL (LIFAIL) 
      INTEGER IWORK (LIWORK), IFAIL (LIFAIL) 
                                                                        
      CHARACTER(16) LABEL 
!*PT*WARNING* Already double-precision                                  
      DOUBLEPRECISION TEMP 
                                                                        
      DATA IPR / 6 / 
      DATA IPRFLG / 1 / 
! calculate covariance or correlation matrix in symmetric storage mode  
                                                                        
      IER = 0 
      IF (JOPT.EQ.0) THEN 
        CALL DVCMSSM (X, NROW, NCOL, NOBS, MSTA, XMSG, CSSM, LSSM, IER) 
      ELSE 
        CALL DCRMSSM (X, NROW, NCOL, NOBS, MSTA, XMSG, CSSM, LSSM, IER) 
      ENDIF 
      IF (IER.NE.0) THEN 
      WRITE (IPR, FMT = '(//'' sub drveof: ier,jopt= '',2i3             &
     &                     ,'' returned from vcmssm/crmssm'')') IER, JOP&
     &T                                                                 
        IF (IER.GT.0) RETURN 
      ENDIF 
                                                                        
! activate if print of cov/cor matrix [work] is desired                 
                                                                        
      IF (IPRFLG.EQ.2) THEN 
        IF (JOPT.EQ.0) THEN 
      LABEL (1:15)  = 'covar matrix:  ' 
        ELSE 
          LABEL (1:15) = 'correl matrix: ' 
        ENDIF 
        WRITE (IPR, FMT = '(//,a15,''sym storage mode'')') LABEL (1:15) 
        CALL DSSMIOX (CSSM, MSTA) 
      ENDIF 
                                                                        
! calculate the trace  before it is destroyed by sspevx                 
                                                                        
      NA = 0 
!*PT*WARNING* Already double-precision (DBLE)                           
      TEMP = DBLE (0.D0) 
      DO NN = 1, MSTA 
      NA = NA + NN 
!*PT*WARNING* Already double-precision (DBLE)                           
      TEMP = TEMP + DBLE (CSSM (NA) ) 
      ENDDO 
!*PT*WARNING* Already double-precision (DBLE)                           
      IF (TEMP.EQ.DBLE (0.D0) ) THEN 
        IER = - 88 
      WRITE (IPR, FMT = '(//'' sub drveof: ier,jopt= '',2i3             &
     &                     ,'' trace=0.0'')') IER, JOPT                 
        RETURN 
      ENDIF 
      TRACE = TEMP 
                                                                        
! calculate the specified number of eigenvalues and the corresponding   
! .   eigenvectors.                                                     
                                                                        
      MEVAL = MIN (NEVAL, MSTA) 
                                                                        
      TOL = 10.D0 * DEPSMACH (IPR) 
      VLOW = 0.0D0 
      VUP = 0.0D0 
      ILOW = MAX (MSTA - MEVAL + 1, 1) 
      IUP = MSTA 
      CALL DSPEVX ('V', 'I', 'U', MSTA, CSSM, VLOW, VUP, ILOW, IUP, TOL,&
      MEVOUT, WEVAL, EVEC, NCOL, WORK, IWORK, IFAIL, INFO)              
                                                                        
      IF (INFO.NE.0) THEN 
        IER = IER + INFO 
      WRITE (IPR, FMT = '(//,'' sub drveof: sspevx error: info=''       &
     &                        ,   i9)') INFO                            
      WRITE (IPR, FMT = '(   '' sub drveof: sspevx error: ifail=''      &
     &                        ,   20i5)')  (IFAIL (I) , I = 1, MIN (20, &
     &MSTA) )                                                           
      ENDIF 
                                                                        
! reverse the order of things from "sspevx" so                          
! .   largest eigenvalues/vectors are first.                            
                                                                        
      DO N = 1, MEVAL 
      WORK (N) = WEVAL (N) 
      ENDDO 
      DO N = 1, MEVAL 
      EVAL (N) = WORK (NEVAL + 1 - N) 
      ENDDO 
                                                                        
      DO N = 1, MSTA 
      DO NN = 1, NEVAL 
      WORK (NN) = EVEC (N, NN) 
      ENDDO 
      DO NN = 1, NEVAL 
      EVEC (N, NN) = WORK (NEVAL + 1 - NN) 
      ENDDO 
      ENDDO 
!                                                                       
! Convert EVEC and EVAL to single precision real                        
!                                                                       
      DO N = 1, MEVAL 
      PCVAR (N) = REAL (EVAL (N) / TRACE) * 100. 
      ENDDO 
                                                                        
      RETURN 
      END SUBROUTINE DDRVEOF                        
                                                                        
      SUBROUTINE DNCLDRV (XX, X, NROW, NCOL, NOBS, NSTA, XMSG, NEVAL,   &
      EVAL, EVEC, PCVAR, TRACE, IOPT, JOPT, PCRIT, EVECX, CSSM, LCSSM,  &
      WORK, LWORK, WEVAL, IWORK, LIWORK, IFAIL, LIFAIL, IER)            
      DOUBLEPRECISION PCRITX 
                                                                        
      INTEGER NROW, NCOL, NOBS, NSTA, NEVAL, IOPT, JOPT, IER 
      INTEGER(8) LCSSM 
      INTEGER LWORK, LIWORK, LIFAIL 
      DOUBLEPRECISION XX (NROW, NCOL), EVAL (NEVAL), EVEC (NCOL, NEVAL),&
      PCRIT, XMSG, TRACE, X (NROW, NCOL), EVECX (NCOL, NEVAL)           
      DOUBLEPRECISION CSSM (LCSSM), WORK (LWORK), WEVAL (LIFAIL) 
      REAL PCVAR (NEVAL) 
      INTEGER IWORK (LIWORK), IFAIL (LIFAIL) 
                                                                        
      DATA IPR / 6 / 
      DATA IPRFLG / 1 / 
                                                                        
      PCRITX = PCRIT * 0.01D0 
!                                   ! counts the total number of        
!                                   ! locations with > pcrit non=msg    
      MSTA = 0 
      DO NC = 1, NSTA 
!                                   ! counter for this location         
      KNT = 0 
      DO NR = 1, NOBS 
      IF (XX (NR, NC) .NE.XMSG) THEN 
        KNT = KNT + 1 
      ENDIF 
      ENDDO 
      IF (DBLE (KNT) / DBLE (NOBS) .GE.PCRITX) THEN 
        MSTA = MSTA + 1 
        DO IROW = 1, NROW 
        X (IROW, MSTA) = XX (IROW, NC) 
        ENDDO 
      ENDIF 
      ENDDO 
                                                                        
!      write (*,'(//'' sub dncldrv: nrow,ncol,nobs,msta= ''             
!     1              ,4i3)') nrow,ncol,nobs,msta                        
      CALL DDRVEOF (X, NROW, NCOL, NOBS, MSTA, XMSG, NEVAL, EVAL, EVECX,&
      PCVAR, TRACE, IOPT, JOPT, CSSM, LCSSM, WORK, LWORK, WEVAL, IWORK, &
      LIWORK, IFAIL, LIFAIL, IER)                                       
                                                                        
! before returning put the evecs in the correct location                
                                                                        
! preset the return array to msg                                        
      DO NE = 1, NEVAL 
      DO NC = 1, NCOL 
      EVEC (NC, NE) = XMSG 
      ENDDO 
      ENDDO 
!                                   ! counts the total number of        
!                                   ! locations with > pcrit non=msg    
      MSTA = 0 
      DO NC = 1, NSTA 
! counter for this location                                             
      KNT = 0 
      DO NR = 1, NOBS 
      IF (XX (NR, NC) .NE.XMSG) THEN 
        KNT = KNT + 1 
      ENDIF 
      ENDDO 
      IF (DBLE (KNT) / DBLE (NOBS) .GE.PCRITX) THEN 
        MSTA = MSTA + 1 
        DO IEVAL = 1, NEVAL 
        EVEC (NC, IEVAL) = EVECX (MSTA, IEVAL) 
        ENDDO 
      ENDIF 
      ENDDO 
                                                                        
      RETURN 
      END SUBROUTINE DNCLDRV                        
