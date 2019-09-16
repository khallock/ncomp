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
  ! .   mflag     - note: if mflag.lt.0 then the missing values at the
  ! .               beginning and end of the series will be set to the
  ! .               value of the nearest non-msg value. if mflag.ge.0
  ! .               then set these values to missing.
  ! .   mptcrt    - if more than "mptcrt" consecutive values are
  ! .               encountered, the routine will not interpolate across
  ! .               that segment. If mptcrt=npts [most common option],
  ! .               then the routine will interpolate as many values as
  ! .               it can.
  ! .
  ! OTHER variables
  ! .   ncode     - code number
  ! .               ncode = -1  : whole series is missing
  ! .               ncode =  0  : series has no missing points upon return
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
  ! This do loop was added later to check for the special
  ! case were all values in X are missing.
  !
  nloop: do  N=1,NPTS
     if (X(N) /= XMSG) GO TO 10
  end do nloop
  return

  ! c c MPTCRT = NPTS   ! updated version
10 NSTRT = 0
  NEND = 0
  NCODE = 0
  NITP = 0
  NPTCRT = IABS(MPTCRT)
  do  N = 1,NPTS
     if (X(N) == XMSG) then
        ! must be a msg pt : set indices
        if (NSTRT == 0) NSTRT = N
        NEND = N
     else

        ! must be a valid pt : check for prior missing values
        !        (1) if nstrt=0 then there are no msg prior values : skip out
        !        (2) if (nend-nstrt+1).gt.nptcrt the set ncode : skip out
        !        (3) if nstrt=1 then initial series values are msg : set to
        !            first non-msg value
        !        ... else
        !            perform the linear interpolation

        if (NSTRT /= 0) then
           if ((NEND-NSTRT+1) > NPTCRT) then
              GO TO 30
           end if

           if (NSTRT == 1) then
              NITP = NITP + (NEND-NSTRT+1)
              if (MFLAG < 0) then
                 do NN = NSTRT,NEND
                    X(NN) = X(N)
                 end do
              else
                 do NN = NSTRT,NEND
                    X(NN) = XMSG
                 end do
              end if
           else
              NBASE = NSTRT - 1
              SLOPE = (X(N)-X(NBASE))/DBLE(N-NBASE)
              NITP = NITP + (NEND-NSTRT+1)
              do  NN = NSTRT,NEND
                 X(NN) = X(NBASE) + SLOPE*DBLE(NN-NBASE)
              end do
           end if

30         NSTRT = 0
           NEND = 0
        end if
     end if
  end do

  ! check the end points

  if (NEND == NPTS) then
     NITP = NITP + (NEND-NSTRT+1)
     if (MFLAG < 0) then
        do  NN = NSTRT,NPTS
           X(NN) = X(NSTRT-1)
        end do
     else
        do NN = NSTRT,NPTS
           X(NN) = XMSG
        end do
     end if
  end if

  !     nn = 0
  !     do 60 n=1,npts
  !     if (x(n).eq.xmsg) then
  !         nn = nn+1
  !     endif
  !  60 continue

  !     ncode = nn
  !     if (nn.eq.npts) then
  ! for historical reasons
  !         ncode = -1
  !     endif

  return
END SUBROUTINE DLINMSG
