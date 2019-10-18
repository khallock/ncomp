! -------------------------------------------------------------------
SUBROUTINE DLINMSG(X,NPTS,XMSG,MFLAG,MPTCRT)
  !     implicit none

  ! NCL: xnew = linmsg(x,mflag)

  ! given a series x of length npts : this routine will linearly
  ! .   interpolate to fill in the missing pts. missing values at the
  ! .   beginning and end of the series will be be determined
  ! .   by the sign of mprcrt.
  !
  ! nomenclature :
  ! .   x         - input series which may or may not contain msg values
  ! .   npts      - length of the series
  ! .   xmsg      - missing code
  ! .   mflag     - note: IF mflag.lt.0 THEN the missing values at the
  ! .               beginning and end of the series will be set to the
  ! .               value of the nearest non-msg value. IF mflag.ge.0
  ! .               THEN set these values to missing.
  ! .   mptcrt    - IF more than "mptcrt" consecutive values are
  ! .               encountered, the routine will not interpolate across
  ! .               that segment. IF mptcrt=npts [most common option],
  ! .               THEN the routine will interpolate as many values as
  ! .               it can.
  ! .
  ! OTHER variables
  ! .   ncode     - code number
  ! .               ncode = -1  : whole series is missing
  ! .               ncode =  0  : series has no missing points upon RETURN
  ! .                             to the calling routine. either the serie
  ! .                             had no missing points or this routine
  ! .                             has filled them with interpolated values
  ! .               ncode = nn  : series still has missing values. this
  ! .                             occurs when iabs(mptcrt) is exceeded.
  ! .                             nn is the number of missing values
  ! .                             still present.
  ! .   nitp      - No. of InTerpolated Points : user shouldcheck
  ! .               the ratio (nitp/npts)

  INTEGER NPTS,MPTCRT,NCODE,NITP,N,NEND,NSTRT
  DOUBLE PRECISION X(1:NPTS),XMSG
  INTEGER NPTCRT,NN,NBASE
  DOUBLE PRECISION SLOPE
  !
  ! This DO loop was added later to check for the special
  ! case were all values in X are missing.
  !

  nloop: DO  N=1,NPTS
     IF (X(N) /= XMSG) GO TO 10
  END DO nloop
  RETURN

  ! c c MPTCRT = NPTS   ! updated version
10 NSTRT = 0
  NEND = 0
  NCODE = 0
  NITP = 0
  NPTCRT = IABS(MPTCRT)
  DO  N = 1,NPTS
     IF (X(N) == XMSG) THEN
        ! must be a msg pt : set indices
        IF (NSTRT == 0) NSTRT = N
        NEND = N
     ELSE
        ! must be a valid pt : check for prior missing values
        !        (1) IF nstrt=0 THEN there are no msg prior values : skip out
        !        (2) IF (nend-nstrt+1).gt.nptcrt the set ncode : skip out
        !        (3) IF nstrt=1 THEN initial series values are msg : set to
        !            first non-msg value
        !        ... ELSE
        !            perform the linear interpolation

        IF (NSTRT /= 0) THEN
           IF ((NEND-NSTRT+1) > NPTCRT) THEN
              GO TO 30
           END IF

           IF (NSTRT == 1) THEN
              NITP = NITP + (NEND-NSTRT+1)
              IF (MFLAG < 0) THEN
                 X(NSTRT:NEND) = X(N)
              ELSE
                 X(NSTRT:NEND) = XMSG
              END IF
           ELSE
              NBASE = NSTRT - 1
              SLOPE = (X(N)-X(NBASE))/DBLE(N-NBASE)
              NITP = NITP + (NEND-NSTRT+1)
              DO  NN = NSTRT,NEND
                 X(NN) = X(NBASE) + SLOPE*DBLE(NN-NBASE)
              END DO
           END IF

30         NSTRT = 0
           NEND = 0
        END IF
     END IF
  END DO

  ! check the end points

  IF (NEND == NPTS) THEN
     NITP = NITP + (NEND-NSTRT+1)
     IF (MFLAG < 0) THEN
        X(NSTRT:NPTS) = X(NSTRT-1)
     ELSE
        X(NSTRT:NPTS) = XMSG
     END IF
  END IF

  RETURN
END SUBROUTINE DLINMSG
