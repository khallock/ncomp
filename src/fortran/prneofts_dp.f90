      SUBROUTINE DEOFTS7(X,NROW,NCOL,NOBS,MSTA,XMSG,NEVAL,EVEC,JOPT,    &
     &                   IFLAG,XX,WRK,EVECTS,EVTSAV,IER)                
      DOUBLE PRECISION XBAR 
      DOUBLE PRECISION XVAR 
      DOUBLE PRECISION XSTD 
      DOUBLE PRECISION CON 
                                                                        
! f77                                                                   
                                                                        
! takes output of prneof.f and creates a time series of the             
! .   eof amplitudes                                                    
                                                                        
! nomenclature :                                                        
! .   x         - matrix containing the data.  it contains n            
! .               observations for each of m stations or grid pts.      
! .   nrow,ncol - exact row (observation) and column (station)          
! .               dimensions of x in the calling routine.               
! .   nobs      - actual number of observations (nobs <= nrow)          
! .   msta      - actual number of stations     (msta <= ncol)          
! .   xmsg      - missing code (if no obs are missing set to some       
! .               number which will not be encountered)                 
! .   neval     - no. of eigenvalues and eigen vectors computed by      
! .               prneof.                                               
! .   evec      - array created by prneof                               
! .               this must be dimensioned at least (ncol,neval) in the 
! .               calling routine.                                      
! .   jopt      - =0 means covariance  was used to generate evec        
! .               =1 means correlation was used to generate evec        
! .   xx        - work vector of the same dimensions as x               
! .   wrk       - work vector of length nobs                            
! .   evects    - time series of eof amplitudes for each eigenvalue     
! .   ier       - error code                                            
                                                                        
! INPUT                                                                 
                                                                        
! dimensions of x                                                       
      INTEGER NROW,NCOL,NOBS,MSTA 
                                                                        
! monthly data from station/grid pt                                     
! work array                                                            
! missing code (if any)                                                 
      DOUBLE PRECISION X(NROW,NCOL),EVEC(NCOL,NEVAL),XX(NROW,NCOL),     &
     &                 WRK(NOBS),XMSG                                   
      INTEGER NEVAL,IFLAG 
                                                                        
! OUTPUT                                                                
                                                                        
! monthly anomalies from long term monthly mean                         
      DOUBLE PRECISION EVECTS(NROW,NEVAL) 
      DOUBLE PRECISION EVTSAV(NEVAL) 
! error code                                                            
      INTEGER IER 
                                                                        
                                                                        
      IER = 0 
      IF (NROW.LE.0 .OR. NCOL.LE.0) IER = IER + 1 
      IF (NOBS.LE.0 .OR. MSTA.LE.0) IER = IER + 10 
      IF (IER.NE.0) THEN 
          WRITE (*,FMT='(/,'' sub eofts7: ier='',5i5)') IER,NROW,NCOL,  &
     &      NOBS,MSTA                                                   
!          stop                                                         
          RETURN 
      END IF 
                                                                        
! set to msg as default                                                 
                                                                        
      DO M = 1,NCOL 
          DO N = 1,NROW 
              XX(N,M) = X(N,M) 
          END DO 
      END DO 
                                                                        
! the following was added in Sept 2003 so that the                      
! 'normalization' that follows the "if" would                           
! not be done.                                                          
                                                                        
! c c if (iflag.eq.1 .and. jopt.eq.0) go to 10                          
      IF (JOPT.EQ.1) THEN 
                                                                        
! calculate mean/stddev at each station/grid-pt                         
                                                                        
                                                                        
          DO M = 1,MSTA 
              DO N = 1,NOBS 
                  WRK(N) = X(N,M) 
              END DO 
                                                                        
              CALL DSTAT2(WRK,NOBS,XMSG,XBAR,XVAR,XSTD,KNTX,IER) 
              CON = 1.0D0 
              IF (JOPT.EQ.1 .AND. XSTD.GT.0.D0) THEN 
                  CON = 1.D0/XSTD 
              END IF 
                                                                        
!         write (*,"(' eofts7: m,nobs,kntx,xbar,xstd=',3i5,2f8.2)")     
!     *                        m,nobs,kntx,xbar,xstd                    
                                                                        
              DO N = 1,NOBS 
                  IF (X(N,M).NE.XMSG .AND. XBAR.NE.XMSG) THEN 
                      XX(N,M) = (X(N,M)-XBAR)*CON 
                  ELSE 
                      XX(N,M) = XMSG 
                  END IF 
              END DO 
                                                                        
! end msta                                                              
          END DO 
                                                                        
! c c write (*,"(//,'ANOMALIES',/)")                                    
! c c do n=1,nobs                                                       
! c c     write (*,"(i5 , 10(1x,f9.3) )") n, (xx(n,m),m=1,msta)         
! c c enddo                                                             
                                                                        
      END IF 
                                                                        
! calculate the amplitude time series                                   
                                                                        
   10 DO K = 1,NEVAL 
          DO N = 1,NOBS 
              KNTX = 0 
              EVECTS(N,K) = 0.0D0 
              DO M = 1,MSTA 
                  IF (XX(N,M).NE.XMSG .AND. EVEC(M,K).NE.XMSG) THEN 
                      KNTX = KNTX + 1 
                      EVECTS(N,K) = EVECTS(N,K) + EVEC(M,K)*XX(N,M) 
                  END IF 
              END DO 
              IF (KNTX.EQ.0) EVECTS(N,K) = XMSG 
          END DO 
          CALL DSTAT2(EVECTS(1,K),NOBS,XMSG,EVTSAV(K),XVAR,XSTD,KNTX,   &
     &                IER)                                              
          IF (IFLAG.EQ.0 .AND. EVTSAV(K).NE.XMSG) THEN 
              DO N = 1,NOBS 
                  IF (EVECTS(N,K).NE.XMSG) THEN 
                      EVECTS(N,K) = EVECTS(N,K) - EVTSAV(K) 
                  END IF 
              END DO 
          END IF 
      END DO 
                                                                        
      RETURN 
      END                                           
                                                                        
      SUBROUTINE DEOFTSCA(XX,NROW,NCOL,NOBS,MSTA,XMSG                   &
     &                   ,NEVAL,EVEC,JOPT,IFLAG,EVECTS,PCRIT,IER,X,WRK) 
      IMPLICIT NONE 
                                                                        
! takes output of prneof_ca.f and creates a time series of the          
! .   eof amplitudes. It is essential that xx and pcrit not change      
! .   from the call to prneof_ca.f                                      
                                                                        
! use f90 to dynamically allocate/deallocate memory (work space)        
!     also some array syntax                                            
                                                                        
! nomenclature :                                                        
! .   xx        - matrix containing the data.  it contains n observation
! .               for each of m stations or grid pts.                   
! .   nrow,ncol - exact row (observation) and column (station)          
! .               dimensions of x in the calling routine.               
! .   nobs      - actual number of observations (nobs <= nrow)          
! .   msta      - actual number of stations     (msta <= ncol)          
! .   xmsg      - missing code (if no obs are missing set to some       
! .               number which will not be encountered)                 
! .   neval     - no. of eigenvalues and eigen vectors computed by      
! .               prneof.                                               
! .   evec      - array created by prneof_ca                            
! .               this must be dimensioned at least (ncol,neval) in the 
! .               calling routine.                                      
! .   evects    - time series of eof amplitudes for each eigenvalue     
! .   pcrit     - minimum % of non-missing values required for the      
! .               station/grid-pt be used in the covariance/correlation 
! .               matrix, must be between 0 and 100 inclusive           
! .   ier       - error code                                            
! INPUT                                                                 
                                                                        
! dimensions of xx                                                      
      INTEGER NROW,NCOL,NOBS,MSTA 
      INTEGER NEVAL, JOPT, IFLAG 
                                                                        
! monthly data from station/grid pt                                     
! missing code (if any)                                                 
      DOUBLE PRECISION XX(NROW,NCOL),EVEC(NCOL,NEVAL),XMSG 
! OUTPUT                                                                
                                                                        
! monthly anomalies from long term monthly mean                         
      DOUBLE PRECISION EVECTS(NROW,NEVAL) 
! error code                                                            
      INTEGER IER 
! temporaray arrays                                                     
      DOUBLE PRECISION X(NROW,NCOL),WRK(NOBS) 
                                                                        
      DOUBLE PRECISION PCRIT 
      DOUBLE PRECISION PCRITX 
      DOUBLE PRECISION XBAR 
      DOUBLE PRECISION XVAR 
      DOUBLE PRECISION XSTD 
      DOUBLE PRECISION CON 
      INTEGER NSTA, NC, NR, KNTX, M, N, K 
                                                                        
      IER = 0 
                                                                        
      PCRITX = PCRIT*0.01D0 
                                                                        
! counts the total number of                                            
      NSTA = 0 
      DO NC = 1,MSTA 
! counter for this location                                             
          KNTX = 0 
          DO NR = 1,NOBS 
              IF (XX(NR,NC).NE.XMSG) THEN 
                  KNTX = KNTX + 1 
              END IF 
          END DO 
          IF (DBLE(KNTX)/DBLE(NOBS).GE.PCRITX) THEN 
              NSTA = NSTA + 1 
              DO N = 1,NOBS 
                  X(N,NSTA) = XX(N,NC) 
              END DO 
          END IF 
      END DO 
                                                                        
! this if added Sept 2003 to allow computation without                  
! removing the mean.                                                    
                                                                        
      IF (JOPT.EQ.0 .AND. IFLAG.EQ.1) GO TO 10 
                                                                        
! calculate each station/grid-pt  long-term mean and standard deviation 
! use the "nsta" limit; this pertains to array x                        
                                                                        
      DO M = 1,NSTA 
          DO N = 1,NOBS 
              WRK(N) = X(N,M) 
          END DO 
                                                                        
          CALL DSTAT2(WRK,NOBS,XMSG,XBAR,XVAR,XSTD,KNTX,IER) 
          CON = 1.0D0 
          IF (JOPT.EQ.1 .AND. XSTD.GT.0.D0) THEN 
              CON = 1.D0/XSTD 
          END IF 
                                                                        
          DO N = 1,NOBS 
              IF (X(N,M).NE.XMSG .AND. XBAR.NE.XMSG) THEN 
                  X(N,M) = (X(N,M)-XBAR)*CON 
              ELSE 
                  X(N,M) = XMSG 
              END IF 
          END DO 
! end nsta                                                              
!                                                                       
      END DO 
   10 CONTINUE 
                                                                        
! calculate the amplitude time series                                   
                                                                        
      DO K = 1,NEVAL 
          DO N = 1,NOBS 
              KNTX = 0 
              EVECTS(N,K) = 0.0D0 
              DO M = 1,NSTA 
                  IF (X(N,M).NE.XMSG) THEN 
                      KNTX = KNTX + 1 
                      EVECTS(N,K) = EVECTS(N,K) + EVEC(M,K)*X(N,M) 
                  END IF 
              END DO 
              if (KNTX.EQ.0) EVECTS(N,K) = XMSG 
          END DO 
      END DO 
                                                                        
      RETURN 
      END                                           
