! NCLFORTSTART                                                          
      SUBROUTINE TDRVPRC (X, NROW, NCOL, NROBS, NCSTA, XMSG, NEVAL,     &
      EVAL, EVEC, PCVAR, TRACE, IOPT, JOPT, PCRIT, IER, XDATA, XDATAT)  
      IMPLICIT NONE 
                                                                        
! this operates on the TRANSPOSE of the array "x"                       
! .   It results in [sometimes, MUCH] faster execution                  
! Note for NCL: in *this* routine NROBS=NROW, NCSTA=NCOL                
                                                                        
      INTEGER NROW, NCOL, NROBS, NCSTA, NEVAL, IOPT, JOPT, IER 
      DOUBLEPRECISION X (NROW, NCOL), EVAL (NEVAL), EVEC (NCOL, NEVAL), &
      PCVAR (NEVAL), TRACE, XMSG, PCRIT, XDATA (NROW, NCOL), XDATAT (   &
      NCSTA, NROW)                                                      
! NCLEND                                                                
                                                                        
!======== .so info=====================                                 
! SGI/dataproc: WRAPIT -L /usr/lib64 -l complib.sgimath_mp prneofTranspo
! Sun/CGD: WRAPIT -L /opt/SUNWspro/lib -l sunmath -l fsu -l fui -lsunper
!======================================                                 
                                                                        
                                                                        
! AUTOMATIC FORTRAN ARRAYS TO HOLD DATA AND DATA STATS                  
      DOUBLEPRECISION XAVE (NCOL), XVAR (NCOL), XDVAR (NCOL) 
                                                                        
                                                                        
! LOCAL VARIABLES                                                       
      DOUBLEPRECISION PCX, CON, XSD 
      INTEGER KNTX, KER, NR, NC, MCSTA 
                                                                        
                                                                        
! strip all grid points that have less that "PCRIT" valid values        
! .   create the "XDATA". This may have fewer columns/grid-pts          
! .   than the original "X" array if not all columns                    
! .   had the minimun number of valid values.                           
                                                                        
      MCSTA = 0 
      DO NC = 1, NCOL 
                                                                        
! statistics for this station/grid-point                                
                                                                        
      CALL DSTAT2 (X (1, NC), NROW, XMSG, XAVE (NC), XVAR (NC), XSD,    &
      KNTX, KER)                                                        
                                                                        
! quick way of eliminating stations/grid-points                         
                                                                        
      PCX = (DBLE (KNTX) / DBLE (NROW) ) * 100.D0 
      IF (PCX.LT.PCRIT.OR.XSD.LE.0.0D0) THEN 
        XAVE (NC) = XMSG 
      ENDIF 
                                                                        
! possible normalizing for jopt=1                                       
                                                                        
      CON = 1.0D0 
      IF (JOPT.EQ.1.AND.XAVE (NC) .NE.XMSG.AND.XSD.GT.0.D0) CON = 1.D0 /&
      XSD                                                               
                                                                        
! WORK WITH ANOMALIES: XDAVE=0.0 [or standardized anomalies]            
                                                                        
      IF (XAVE (NC) .NE.XMSG) THEN 
        MCSTA = MCSTA + 1 
                                                                        
        DO NR = 1, NROW 
        IF (X (NR, NC) .NE.XMSG) THEN 
          XDATA (NR, MCSTA) = (X (NR, NC) - XAVE (NC) ) * CON 
! c c                 XDATA(NR,MCSTA) =  X(NR,NC)                       
        ELSE 
          XDATA (NR, MCSTA) = XMSG 
        ENDIF 
        ENDDO 
                                                                        
        IF (JOPT.EQ.0) THEN 
          XDVAR (MCSTA) = XVAR (NC) 
        ELSE 
          XDVAR (MCSTA) = 1.0D0 
        ENDIF 
      ENDIF 
                                                                        
      ENDDO 
                                                                        
! pass the selected data (XDATA) to the EOF driver                      
                                                                        
      CALL XRVEOFT (XDATA, XDATAT, NROW, NCOL, NROBS, MCSTA, XMSG,      &
      NEVAL, EVAL, EVEC, PCVAR, TRACE, XDVAR, XAVE, JOPT, IER)          
                                                                        
      RETURN 
      END SUBROUTINE TDRVPRC                        
                                                                        
! ---------------------------------------------------------             
      SUBROUTINE XRVEOFT (XDATA, XDATAT, NROW, NCOL, NROBS, NCSTA, XMSG,&
      NEVAL, EVAL, EVEC, PCVAR, TRACE, XDVAR, XAVE, JOPT, IER)          
      IMPLICIT NONE 
                                                                        
! **** USES AUTOMATIC ARRAYS SO USE WITH g77 or f90                     
                                                                        
! operate on the *TRANSPOSE* OF XDATA: then use matrix stuff to         
! .   get the desired eof information.                                  
                                                                        
! Note: NROW=NROBS but, it is possible for  MCSTA<=NCOL                 
! .                     when some stations/grid_pts do not have enough d
                                                                        
! driver to calculate :                                                 
! .   the principal components (eignvalues and eigenvectors)            
! .       of the data array XDATA. XDATA may contain                    
! .       missing observations. If it has msg data, this will           
! .       calculate a var-cov matrix but it may not be                  
! .       positive definite.                                            
                                                                        
! . The eigenvectors and eigenvalues of the *TRANSPOSED* matrix         
! .      are used to derive the corresponding eigenvectors              
! .      associated with the original matrix.                           
                                                                        
! .       USES LAPACK/BLAS ROUTINES                                     
                                                                        
! nomenclature :                                                        
! .   xdata     - matrix containing the data. it contains n observations
! .               for each of m stations or grid pts.                   
! .   nrow,ncol - exact row (observation) and column (station)          
! .               dimensions of xdata in the calling routine.           
! .   nrobs     - actual number of observations (nrobs <= nrow)         
! .   ncsta     - actual number of stations     (ncsta <= ncol)         
! .   xmsg      - missing code (if no obs are missing set to some       
! .               number which will not be encountered)                 
! .   neval     - no. of eigenvalues and eigenvectors to be computed.   
! .               neval <= min(nrobs,ncsta).                            
! .               If not, ncsta eigenvalues and eigenvectors            
! .               will be computed and ier = -2.                        
! .   eval      - vector containing the eigenvalues in DESCENDING order.
! .               eval must be at least neval in length.                
! .   evec      - an array which will contains the eigenvector info.    
! .               this must be dimensioned at least (ncol,neval) in the 
! .               calling routine. There is some normalization done     
! .               but I am not sure waht it is.                         
! .   pcvar     - contains % variance associated with each eigenvalue.  
! .               This must be dimensioned at least neval in length.    
! .   trace     - trace of the variance - covariance matrix.in the      
! .               case of a var-covar matrix , the trace should         
! .               equal the sum of all possible eigenvalues.            
! .   iopt      - not used; kept for compatibility with old routine     
! .   jopt      - =  0 : use var-covar matrix in prncmp                 
! .                  1 : use correlation matrix in prncmp               
! .   ier       - error code                                            
                                                                        
      INTEGER NROW, NCOL, NROBS, NCSTA, NEVAL, LSSM, LWORK, LIWORK,     &
      JOPT, IER                                                         
                                                                        
      DOUBLEPRECISION XDATA (NROW, NCOL), PCVAR (NEVAL), XMSG, TRACE,   &
      EVAL (NEVAL), EVEC (NCOL, NEVAL), XDATAT (NCSTA, NROW)            
                                                                        
! temporary arrays (automatic or passed in via interface)               
                                                                        
      DOUBLEPRECISION CSSM (NROW * (NROW + 1) / 2), WORK (8 * NROW) 
      INTEGER IWORK (5 * NROW) 
                                                                        
      DOUBLEPRECISION TEOFPC (NROW, NEVAL), WEVAL (NROW), WEVEC (NCSTA, &
      NEVAL), W2D (NROW, NEVAL), XAVE (NCOL), XDVAR (NCOL)              
      INTEGER IFAIL (NROBS) 
                                                                        
! local                                                                 
      CHARACTER(16) LABEL 
      DOUBLEPRECISION TEMP, TOL, EPSMACH, VLOW, VUP 
      INTEGER N, NE, NR, NC, MEVAL, IPR, IPRFLG, MEVOUT, ILOW, IUP,     &
      INFO, MCSTA                                                       
                                                                        
      DATA IPR / 6 / 
      DATA IPRFLG / 0 / 
                                                                        
! length of various temporary arrays                                    
                                                                        
      LSSM = NROW * (NROW + 1) / 2 
      LWORK = 8 * NROW 
      LIWORK = 5 * NROW 
                                                                        
! EXPLICITLY INITALIZE                                                  
                                                                        
      DO NR = 1, NROW 
      WEVAL (NR) = XMSG 
      ENDDO 
                                                                        
      DO NE = 1, NEVAL 
      EVAL (NE) = XMSG 
                                                                        
      DO NC = 1, NCOL 
      EVEC (NC, NE) = XMSG 
      ENDDO 
                                                                        
      DO NC = 1, NCSTA 
      WEVEC (NC, NE) = XMSG 
      ENDDO 
                                                                        
      DO NR = 1, NROW 
      W2D (NR, NE) = XMSG 
      TEOFPC (NR, NE) = XMSG 
      ENDDO 
      ENDDO 
                                                                        
! create transposed data array                                          
                                                                        
      DO NC = 1, NCSTA 
      DO NR = 1, NROW 
      XDATAT (NC, NR) = XDATA (NR, NC) 
      ENDDO 
      ENDDO 
                                                                        
! compute the covariance matrix using the transposed data               
                                                                        
      CALL DVCMSSM (XDATAT, NCSTA, NROW, NCSTA, NROW, XMSG, CSSM, LSSM, &
      IER)                                                              
                                                                        
      IF (IER.NE.0) THEN 
      WRITE (IPR, FMT = '(//'' sub drveoft: ier= '',i3                  &
     &  ,'' returned from vcmssm/crmssm'')') IER                        
        IF (IER.GT.0) RETURN 
      ENDIF 
                                                                        
! activate if print of cov/cor matrix [work] is desired                 
                                                                        
      IF (IPRFLG.EQ.2) THEN 
        IF (JOPT.EQ.0) THEN 
      LABEL (1:15)  = 'DRVEOFT: covar matrix:  ' 
        ELSE 
          LABEL (1:15) = 'DRVEOFT: correl matrix: ' 
        ENDIF 
        WRITE (IPR, " (//, a15, 'sym storage mode') ") LABEL (1:15) 
        CALL DSSMIOX (CSSM, NROW) 
      ENDIF 
                                                                        
! calculate the trace of transposed matrix (TRACE)                      
! before matrix is destroyed by sspevx                                  
                                                                        
      N = 0 
      TRACE = 0.d0 
      DO NR = 1, NROW 
      N = N + NR 
      TRACE = TRACE+CSSM (N) 
      ENDDO 
      IF (TRACE.LE.0.d0) THEN 
        IER = - 88 
      WRITE (IPR, FMT = '(//'' SUB DRVEOFT: ier,jopt= '',2i3            &
     &          ,'' trace='',f15.5)') IER, JOPT, TRACE                  
        RETURN 
      ENDIF 
                                                                        
                                                                        
! calculate specified number of eigenvalues and eigenvectors.           
! .   make sure that  neval <= nrow (=nrobs)                            
!                                                                       
! .   Remember the TEOFPC are the eofs of the *transposed* matrix.      
! .   This means they are the principal components of the original data.
                                                                        
      MEVAL = MIN (NEVAL, NROW) 
                                                                        
      TOL = 10.D0 * EPSMACH (IPR) 
      VLOW = 0.0D0 
      VUP = 0.0D0 
      ILOW = MAX (NROBS - MEVAL + 1, 1) 
      IUP = NROBS 
      MEVOUT = 0 
                                                                        
      CALL DSPEVX ('V', 'I', 'U', NROW, CSSM, VLOW, VUP, ILOW, IUP, TOL,&
      MEVOUT, WEVAL, TEOFPC, NROW, WORK, IWORK, IFAIL, INFO)            
                                                                        
      IF (INFO.NE.0) THEN 
        IER = IER + INFO 
      WRITE (IPR, FMT = '(//,'' SUB DRVEOFT: sspevx error: info=''      &
     &                      ,   i9)') INFO                              
      WRITE (IPR, FMT = '(   '' SUB DRVEOFT: sspevx error: ifail=''     &
     &     ,   20i5)')  (IFAIL (N) , N = 1, MIN (20, NROBS) )           
      ENDIF 
                                                                        
! reorder so eigenvalues are in descending order                        
! .   make sure the asocciated eigenvectors are also reordered          
                                                                        
      DO N = 1, MEVAL 
      WORK (N) = WEVAL (N) 
      ENDDO 
                                                                        
      DO N = 1, MEVAL 
      WEVAL (N) = WORK (MEVAL + 1 - N) 
      ENDDO 
                                                                        
      DO N = 1, MEVAL 
      DO NR = 1, NROW 
      W2D (NR, N) = TEOFPC (NR, N) 
      ENDDO 
      ENDDO 
                                                                        
      DO N = 1, MEVAL 
      DO NR = 1, NROW 
      TEOFPC (NR, N) = W2D (NR, MEVAL + 1 - N) 
      ENDDO 
      ENDDO 
                                                                        
! percent variance explained (PCVAR) of transposed matrix.              
! .   Both the eigenvalues (hence, PCVAR) are equivalent                
! .   due to "duality" od space/time.                                   
                                                                        
      DO N = 1, MEVAL 
      PCVAR (N) = (WEVAL (N) / TRACE) * 100.D0 
      ENDDO 
                                                                        
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
! section to invert from transposed data space to original data space.  
! ie: compute the eigenvectors of the original data                     
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
                                                                        
! calculate the SPATIAL eofs from the TEOFPC and XDATAT                 
                                                                        
      DO N = 1, MEVAL 
      DO NC = 1, NCSTA 
      WEVEC (NC, N) = 0.0d0 
                                                                        
      DO NR = 1, NROW 
      IF (XDATAT (NC, NR) .NE.XMSG) THEN 
        WEVEC (NC, N) = WEVEC (NC, N) + TEOFPC (NR, N) * XDATAT (NC, NR) 
      ENDIF 
      ENDDO 
                                                                        
      ENDDO 
      ENDDO 
                                                                        
! normalize spatial eofs to variance=1.0                                
! normalize the time series [principal components]                      
                                                                        
      DO N = 1, MEVAL 
                                                                        
      TEMP = 0.0d0 
      DO NC = 1, NCSTA 
      TEMP = TEMP + WEVEC (NC, N) * WEVEC (NC, N) 
      ENDDO 
      TEMP = DSQRT (TEMP) 
                                                                        
      DO NC = 1, NCSTA 
      WEVEC (NC, N) = WEVEC (NC, N) / TEMP 
      ENDDO 
                                                                        
      ENDDO 
                                                                        
! return only the requested number of eigenvalues                       
                                                                        
      DO N = 1, NEVAL 
      EVAL (N) = WEVAL (N) 
      DO NC = 1, NCOL 
      EVEC (NC, N) = XMSG 
      ENDDO 
      ENDDO 
                                                                        
! return the requested eigenvector                                      
! .   the purpose of the "if" is to reassign to correct locations       
                                                                        
      MCSTA = 0 
      DO NC = 1, NCOL 
      IF (XAVE (NC) .NE.XMSG) THEN 
        MCSTA = MCSTA + 1 
        DO N = 1, NEVAL 
        EVEC (NC, N) = WEVEC (MCSTA, N) 
        ENDDO 
      ENDIF 
      ENDDO 
                                                                        
      RETURN 
      END SUBROUTINE XRVEOFT                        
