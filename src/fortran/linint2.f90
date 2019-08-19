! -----------------------------------------------------------
! NCLFORTSTART
SUBROUTINE DLININT1(NXI,XI,FI,ICYCX,NXO,XO,FO,XIW,FXIW,NXI2,XMSG,&
     &                    IOPT,IER)
  IMPLICIT NONE
  INTEGER NXI,NXO,NXI2,IOPT,ICYCX,IER
  DOUBLE PRECISION XI(NXI),FI(NXI)
  DOUBLE PRECISION XO(NXO),FO(NXO),XMSG
  DOUBLE PRECISION XIW(NXI2),FXIW(NXI2)
  ! NCLEND

  ! NCL:  fo = linint1 (xi,fi, wrapX, xo, iopt)
  !
  !            {nxi,nxo} = dimsizes( {xi,xo} )
  !            fo is the same size xo and same type as "fi"
  !            xmsg = fi@_FillValue
  !            wrapX is an NCL logical; a False ==> 0
  !                                     a True  ==> 1
  !            iopt = 0 try to preserve msg areas
  !                 = 1 try to fill in entire grid
  !
  !            The NCL wrapper should allow for multiple datasets
  !            so the user need only make one call to the function.

  ! perform piecewise linear interpolation allowing for missing data
  ! .  nothing fancy

  ! nomenclature:
  ! .   nxi     - lengths of xi and dimensions of fi (must be >= 2)
  ! .   xi      - coordinates of fi (eg, lon)
  ! .             must be monotonically (in/de)creasing [same as xo]
  ! .   fi      - functional input values [2D]
  ! .   icycx   - 0 if fi is cyclic in x (cyclic pt should NOT included)
  ! .             .ne.0 if cyclic
  ! .   nxo     - lengths of xo and dimensions of fo (must be >= 1)
  ! .   xo      - coordinates of fo (eg, lon)
  ! .             must be monotonically (in/de)creasing [same as xi]
  ! .   fo      - functional output values [interpolated]
  ! .   xmsg    - missing code
  ! .   iopt    - options
  ! .           - iopt=0  means to try to preserve missing areas
  ! .           - iopt=1  means to try to try to fill-in missing areas
  ! .   ier     - error code
  ! .             =0;   no error
  ! .             =1;   not enough points in input/output array
  ! .             =2;   xi are not monotonically increasing
  ! .             =4;   xo yo are not monotonically increasing
  !
  !                              local
  INTEGER NX,NPTS,IFLAG,NXSTRT,NXLAST
  !                              gross error checking
  IER = 0
  if (NXO < 1) then
     IER = 1
     return
  end if
  !                              initialize to msg
  FO(:) = XMSG

  if (NXI <= 1) then
     IER = 1
     return
  end if
  !                              mono (in/de)creasing ?
  CALL DMONOID2(NXI,XI,NXO,XO,IFLAG,IER)
  if (IFLAG == 0 .or. IER /= 0) return
  !                              are data to be treated cyclic?
  if (ICYCX == 0) then
     !                              data are not cyclic
     if (IOPT == 0) then

        CALL DLIN2INT1(NXI,XI,FI,NXO,XO,FO,XMSG,IFLAG)
     else if (IOPT == 1) then
        !                              collapse data array (eliminate msg)
        NPTS = 0
        do NX = 1,NXI
           if (FI(NX) /= XMSG) then
              NPTS = NPTS + 1
              XIW(NPTS+1)  = XI(NX)
              FXIW(NPTS+1) = FI(NX)
           end if
        end do

        CALL DLIN2INT1(NPTS,XIW(2),FXIW(2),NXO,XO,FO,XMSG,IFLAG)
     end if
  else
     !                              data are cyclic
     if (IOPT == 0) then
        !                              preserve msg region
        do NX = 1,NXI
           XIW(NX+1)  = XI(NX)
           FXIW(NX+1) = FI(NX)
        end do

        CALL DLINCYC(NXI,XI,FI,1,NXI,IFLAG,XIW,FXIW,NXI)
        CALL DLIN2INT1(NXI+2,XIW,FXIW,NXO,XO,FO,XMSG,IFLAG)

     else if (IOPT == 1) then
        !                              collapse data array
        NPTS = 0
        NXSTRT = 0
        NXLAST = 0
        do NX = 1,NXI
           if (FI(NX) /= XMSG) then
              NPTS = NPTS + 1
              XIW(NPTS+1)  = XI(NX)
              FXIW(NPTS+1) = FI(NX)
              if (NPTS == 1) NXSTRT = NX
              NXLAST = NX
           end if
        end do

        if (NPTS == 0) then
           IER = 1
           return
        end if

        CALL DLINCYC(NXI,XI,FI,NXSTRT,NXLAST,IFLAG,XIW,FXIW,NPTS)
        CALL DLIN2INT1(NPTS+2,XIW,FXIW,NXO,XO,FO,XMSG,IFLAG)
     end if

  end if

  return
END SUBROUTINE DLININT1
! -----------------------------------------------------------
! NCLFORTSTART
SUBROUTINE DLININT2(NXI,XI,NYI,YI,FI,ICYCX,NXO,XO,NYO,YO,FO,XIW,&
     &                    FXIW,NXI2,XMSG,IOPT,IER)
  IMPLICIT NONE
  INTEGER NXI,NYI,NXO,NYO,NXI2,ICYCX,IOPT,IER
  DOUBLE PRECISION XI(NXI),YI(NYI),FI(NXI,NYI)
  DOUBLE PRECISION XO(NXO),YO(NYO),FO(NXO,NYO),XMSG
  DOUBLE PRECISION XIW(NXI2),FXIW(NXI2)
  ! NCLEND

  ! NCL:  fo = linint2 (xi,yi,fi, wrapX, xo,yo, iopt)
  !
  !            {nxi,nyi,nxo,nyo} = dimsizes( {xi,yi,xo,yo} )
  !            fo is the same size xo, yo and same type as "fi"
  !            xmsg = fi@_FillValue
  !            wrapX is an NCL logical; a False ==> 0
  !                                     a True  ==> 1
  !
  !            The NCL wrapper should allow for multiple datasets
  !            so the user need only make one call to the function.

  ! perform bilinear interpolation allowing for missing data
  ! .  nothing fancy

  ! nomenclature:
  ! .   nxi,nyi - lengths of xi,yi and dimensions of fi (must be >= 2)
  ! .   xi      - coordinates of fi (eg, lon)
  ! .             must be monotonically increasing
  ! .   yi      - coordinates of fi (eg, lat)
  ! .             must be monotonically increasing
  ! .   fi      - functional input values [2D]
  ! .   icycx   - 0 if fi is cyclic in x (cyclic pt should NOT included)
  ! .             .ne.0 if cyclic
  ! .   nxo,nyo - lengths of xo,yo and dimensions of fo (must be >= 2)
  ! .   xo      - coordinates of fo (eg, lon)
  ! .             must be monotonically increasing
  ! .   yo      - coordinates of fo (eg, lat)
  ! .             must be monotonically increasing
  ! .   fo      - functional output values [interpolated]
  ! .   xmsg    - missing code
  ! .   iopt    - =0 try to preserve msg areas
  ! .             =1 try to fill in entire grid
  ! .   ier     - error code
  ! .             =0;   no error
  ! .             =1;   not enough points in input/output array
  ! .             =2;   xi or xo are not monotonically (in/de)creasing
  ! .             =3;   yi or yo are not monotonically (in/de)creasing
  !
  !                              local and temporary/work arrays
  INTEGER IFLAG,NPTS,NXSTRT,NXLAST
  DOUBLE PRECISION YIW(NYI),FYIW(NYI),FOYW(NYO)
  DOUBLE PRECISION FTMP(NXO,NYI)

  !                              local
  INTEGER NX,NY

  !                              error checking
  IER = 0
  if (NXO < 1 .or. NYO < 1) then
     IER = 1
     return
  end if

  !                              set to msg
  FO = XMSG

  if (NXI < 2 .or. NYI < 2) then
     IER = 1
     return
  end if

  !                              error mono increasing ?
  CALL DMONOID2(NXI,XI,NXO,XO,IFLAG,IER)
  if (IFLAG == 0 .or. IER /= 0) then
     IER = 2
     return
  end if
  CALL DMONOID2(NYI,YI,NYO,YO,IFLAG,IER)
  if (IFLAG == 0 .or. IER /= 0) then
     IER = 3
     return
  end if

  !                               is the input array cyclic in x
  if (ICYCX == 0) then
     !                               preserve missing areas
     if (IOPT == 0) then
        !                               interpolate in the x direction
        do NY = 1,NYI
           CALL DLIN2INT1(NXI,XI,FI(1,NY),NXO,XO,FTMP(1,NY),&
                &                          XMSG,IFLAG)
        end do
        !                               interpolate in the y direction
        do NX = 1,NXO
           FYIW(:) = FTMP(NX,:) ! : = NYI
           CALL DLIN2INT1(NYI,YI,FYIW,NYO,YO,FOYW,XMSG,IFLAG)
           FO(NX,:) = FOYW(:) ! : = NYO
        end do
     else if (IOPT == 1) then
        !                               interpolate in the x direction
        !                               collapse data array
        do NY = 1,NYI
           NPTS = 0
           do NX = 1,NXI
              if (FI(NX,NY) /= XMSG) then
                 NPTS = NPTS + 1
                 XIW(NPTS+1)  = XI(NX)
                 FXIW(NPTS+1) = FI(NX,NY)
              end if
           end do

           CALL DLIN2INT1(NPTS,XIW(2),FXIW(2),NXO,XO,FTMP(1,NY),&
                &                          XMSG,IFLAG)
        end do
        !                               interpolate in the y direction
        do NX = 1,NXO
           NPTS = 0
           do NY = 1,NYI
              if (FTMP(NX,NY) /= XMSG) then
                 NPTS = NPTS + 1
                 YIW(NPTS) = YI(NY)
                 FYIW(NPTS) = FTMP(NX,NY)
              end if
           end do

           CALL DLIN2INT1(NPTS,YIW,FYIW,NYO,YO,FOYW,XMSG,IFLAG)

           do NY = 1,NYO
              FO(NX,NY) = FOYW(NY)
           end do
        end do

     end if
  else
     !                               must be cyclic in x
     !                               create cyclic "x" coordinates
     if (IOPT == 0) then
        do NY = 1,NYI
           do NX = 1,NXI
              XIW(NX+1)  = XI(NX)
              FXIW(NX+1) = FI(NX,NY)
           end do

           NPTS = NXI
           CALL DLINCYC(NXI,XI,FI(1,NY),1,NXI,IFLAG,XIW,FXIW,&
                &                         NPTS)
           CALL DLIN2INT1(NXI+2,XIW,FXIW,NXO,XO,FTMP(1,NY),XMSG,&
                &                           IFLAG)
        end do
        !                               interpolate in the y direction
        do NX = 1,NXO
           FYIW(:) = FTMP(NX,:) ! : = NYI
           CALL DLIN2INT1(NYI,YI,FYIW,NYO,YO,FOYW,XMSG,IFLAG)
           FO(NX,:) = FOYW(:) ! : = NYO
        end do

     else if (IOPT == 1) then
        do NY = 1,NYI
           !                              collapse data array
           NPTS = 0
           NXSTRT = 0
           NXLAST = 0
           do NX = 1,NXI
              if (FI(NX,NY) /= XMSG) then
                 NPTS = NPTS + 1
                 XIW(NPTS+1)  = XI(NX)
                 FXIW(NPTS+1) = FI(NX,NY)
                 if (NPTS == 1) NXSTRT = NX
                 NXLAST = NX
              end if
           end do

           CALL DLINCYC(NXI,XI,FI(1,NY),NXSTRT,NXLAST,IFLAG,XIW,&
                &                         FXIW,NPTS)
           CALL DLIN2INT1(NPTS+2,XIW,FXIW,NXO,XO,FTMP(1,NY),&
                &                           XMSG,IFLAG)
        end do
        !                               interpolate in the y direction
        do NX = 1,NXO
           NPTS = 0
           do NY = 1,NYI
              if (FTMP(NX,NY) /= XMSG) then
                 NPTS = NPTS + 1
                 YIW(NPTS) = YI(NY)
                 FYIW(NPTS) = FTMP(NX,NY)
              end if
           end do

           CALL DLIN2INT1(NPTS,YIW,FYIW,NYO,YO,FOYW,XMSG,IFLAG)

           FO(NX, :) = FOYW(:) ! : = NYO
        end do
     end if
  end if

  return
END SUBROUTINE DLININT2
! -----------------------------------------------------------
SUBROUTINE DLIN2INT1(NIN,XI,FI,NOUT,XO,FO,XMSG,IFLAG)
  IMPLICIT NONE
  INTEGER NIN,NOUT,IFLAG
  DOUBLE PRECISION XI(NIN),FI(NIN),XO(NOUT),FO(NOUT),XMSG

  ! perform 1D piecewise linear interpolation allowing for missing data
  ! .  this works with lin2int [does no error chk]
  ! .  nothing fancy
  !                              local
  INTEGER NI,NO,NISTRT,NIS,NEXACT
  DOUBLE PRECISION SLOPE

  ! PRINT *, "Hello World!"
  ! write ( *, '(a)' ) '  Hello, world!'

  FO(:) = XMSG ! : = NOUT

  ! main loop [exact matches]
  ! nistrt minimizes extra checks
  NEXACT = 0
  NISTRT = 1
  NIS = NISTRT
  noloop: do NO = 1,NOUT
     niloop: do NI = NISTRT,NIN
        if (XO(NO) == XI(NI)) then
           FO(NO) = FI(NI)
           NIS = NI + 1
           NEXACT = NEXACT + 1
           exit niloop
        end if
     end do niloop
     NISTRT = NIS
  end do noloop

  !                              main loop [interpolation]
  if (IFLAG == 1 .or. IFLAG == -1) then
     do NO = 1,NOUT
        do NI = 1,NIN - 1
           if ((XO(NO) > XI(NI) .and. XO(NO) < XI(NI+1) .and. IFLAG == 1)&
                &.or. (XO(NO) < XI(NI) .and. XO(NO) > XI(NI+1) .and. IFLAG == -1)) then
              if (FI(NI) /= XMSG .and. FI(NI+1) /= XMSG) then
                 SLOPE = (FI(NI+1)-FI(NI))/ (XI(NI+1)-XI(NI))
                 FO(NO) = FI(NI) + SLOPE* (XO(NO)-XI(NI))
              end if
           end if
        end do
     end do
  end if

  return
END SUBROUTINE DLIN2INT1
! ---------------------------------------------
SUBROUTINE DLINCYC(NXI,XI,FI,NXSTRT,NXLAST,IFLAG,XIW,FXIW,NPTS)
  !
  ! handle the "x" cyclic point
  !
  IMPLICIT NONE
  INTEGER NXI,NXSTRT,NXLAST,IFLAG,NPTS
  DOUBLE PRECISION XI(NXI),FI(NXI),XIW(0:NXI+1),FXIW(0:NXI+1)
  DOUBLE PRECISION DX


  if (NXSTRT == 1 .and. NXLAST == NXI) then
     DX = ABS(XI(2)-XI(1))
     XIW(0) = XI(1) - IFLAG*DX
     FXIW(0) = FI(NXI)
     DX = ABS(XI(NXI)-XI(NXI-1))
     XIW(NPTS+1) = XI(NXI) + IFLAG*DX
     FXIW(NPTS+1) = FI(1)
  else
     !                              arbitrary
     DX = (NXSTRT+ (NXI-NXLAST))*ABS(XI(2)-XI(1))
     XIW(0) = XIW(1) - IFLAG*DX
     FXIW(0) = FXIW(NPTS)
     DX = ((NXI-NXLAST)+NXSTRT)*ABS(XI(NXI)-XI(NXI-1))
     XIW(NPTS+1) = XIW(NPTS) + IFLAG*DX
     FXIW(NPTS+1) = FXIW(1)
  end if

  return
END SUBROUTINE DLINCYC
! ---------------------------------------------
SUBROUTINE DMONOINC(X,NX,NER,IER)
  IMPLICIT NONE

  ! chk to make sure that x is monotonically increasing

  INTEGER NX,NER,IER
  DOUBLE PRECISION X(NX)
  !                          local
  INTEGER N

  IER = 0
  if (NX <= 1) return

  do N = 1,NX - 1
     if (X(N+1) <= X(N)) then
        IER = NER
        return
     end if
  end do

  return
END SUBROUTINE DMONOINC
!     ---------------------------------------------------
SUBROUTINE DMONOID1(NIN,XI,IFLAG,IER)
  !
  ! chk to see if a series is mono (in/de)creasing
  !
  IMPLICIT NONE
  INTEGER NIN,IFLAG,IER
  DOUBLE PRECISION XI(NIN)
  !                              local
  INTEGER NI

  IER = 0
  IFLAG = 0
  !
  ! if only one element, then treat it as monotonically increasing
  !
  if (NIN < 2) then
     IFLAG = +1
     return
  end if

  if (XI(2) > XI(1)) then
     !                              ? mono INcreasing
     do NI = 1,NIN - 1
        if (XI(NI+1) <= XI(NI)) then
           IER = 2
           return
        end if
     end do

     IFLAG = +1
  else
     !                              ? mono DEcreasing
     do NI = 1,NIN - 1
        if (XI(NI+1) >= XI(NI)) then
           IER = 2
           return
        end if
     end do

     IFLAG = -1
  end if

  return
END SUBROUTINE DMONOID1
!     ---------------------------------------------------
SUBROUTINE DMONOID2(NIN,XI,NOUT,XO,IFLAG,IER)
  !
  ! make sure the two series are mono (in/de)creasing
  !
  IMPLICIT NONE
  INTEGER NIN,NOUT,IFLAG,IER
  DOUBLE PRECISION XI(NIN),XO(NOUT)
  !                              local
  INTEGER NFLAG,NER

  IER = 0
  IFLAG = 0
  NFLAG = 0

  CALL DMONOID1(NIN,XI,IFLAG,IER)
  if (IFLAG == 0) then
     IER = 2
  else
     CALL DMONOID1(NOUT,XO,NFLAG,NER)
     if (NFLAG == 0) then
        IER = 3
     else if (IFLAG /= NFLAG) then
        IER = 4
     end if
  end if

  return
END SUBROUTINE DMONOID2

! -----------------------------------------------------------
! NCLFORTSTART
SUBROUTINE DLININT2PTS(NXI,XI,NYI,YI,FI,ICYCX,NXYO,XO,YO,FO,&
     &                       XIW,FIXW,NXI2,XMSG,IER)
  IMPLICIT NONE
  INTEGER NXI,NYI,NXYO,ICYCX,NXI2,IER
  DOUBLE PRECISION XI(NXI),YI(NYI),FI(NXI,NYI)
  DOUBLE PRECISION XO(NXYO),YO(NXYO),FO(NXYO)
  DOUBLE PRECISION XIW(NXI2),FIXW(NXI2,NYI),XMSG
  ! NCLEND

  ! NCL:  fo = linint2_points (xi,yi,fi, wrapX, xo,yo, foOpt)
  !
  !            {nxi,nyi} = dimsizes( {xi,yi} )
  !            {nxyo} = dimsizes( {xo} )  = dimsizes( {yo} )
  !            fo is the same size xo, yo and same type as "fi"
  !            xmsg = fi@_FillValue
  !            wrapX is an NCL logical; a False ==> 0
  !                                     a True  ==> 1
  !            foOpt unused option
  !
  !            The NCL wrapper should allow for multiple datasets
  !            so the user need only make one call to the function.

  ! perform 2D piecewise linear interpolation allowing for missing data
  ! .  nothing fancy

  ! nomenclature:
  ! .   nxi,nyi - lengths of xi,yi and dimensions of fi (must be >= 2)
  ! .   xi      - coordinates of fi (eg, lon)
  ! .             must be monotonically increasing
  ! .   yi      - coordinates of fi (eg, lat)
  ! .             must be monotonically increasing
  ! .   fi      - functional input values [2D]
  ! .   icycx   - 0 if fi is cyclic in x (cyclic pt should NOT included)
  ! .             .ne.0 if cyclic
  ! .   nxyo    - number of xo,yo coordinate pairs
  ! .   xo,yo   - coordinate pairs
  ! .   fo      - functional output values [interpolated]
  ! .   xmsg    - missing code
  ! .   nopt    - option for how to handle case when not all values
  ! .              present
  ! .              >= 0 means set the interpolated point to xmsg
  ! .              =-1 try some sort of weighted average
  ! .   ier     - error code
  ! .             =0;   no error
  ! .             =1;   not enough points in input/output array
  ! .             =2/3; xi or yi are not monotonically increasing
  !
  !                              automatic temporary/work arrays
  DOUBLE PRECISION DX

  !                              local
  INTEGER NX,NY,NOPT
  !                              this could be made an argument
  NOPT = -1
  !                              error checking
  IER = 0
  if (NXI < 2 .or. NYI < 2) then
     IER = 1
     return
  end if
  !                              error mono increasing ?
  CALL DMONOINC(XI,NXI,2,IER)
  if (IER /= 0) return
  CALL DMONOINC(YI,NYI,3,IER)
  if (IER /= 0) return
  !                               is the input array cyclic in x
  if (ICYCX == 0) then
     CALL DLINT2XY(NXI,XI,NYI,YI,FI,NXYO,XO,YO,FO,XMSG,NOPT,IER)
  else
     !                               must be cyclic in x
     !                               create cyclic "x" coordinates
     do NX = 1,NXI
        XIW(NX+1) = XI(NX)
     end do
     DX = XI(2) - XI(1)
     XIW(1)    = XI(1) - DX
     XIW(NXI2) = XI(NXI) + DX

     do NY = 1,NYI
        do NX = 1,NXI
           FIXW(NX+1,NY) = FI(NX,NY)
        end do
        FIXW(1,NY)    = FI(NXI,NY)
        FIXW(NXI2,NY) = FI(1,NY)
     end do
     CALL DLINT2XY(NXI2,XIW,NYI,YI,FIXW,NXYO,XO,YO,FO,XMSG,NOPT,&
          &                  IER)
  end if

  return
END SUBROUTINE DLININT2PTS

! -----------------------------------------------------------
SUBROUTINE DLINT2XY(NXI,XI,NYI,YI,FI,NXYO,XO,YO,FO,XMSG,NOPT,IER)
  IMPLICIT NONE
  INTEGER NXI,NYI,NXYO,NOPT,IER
  DOUBLE PRECISION XI(NXI),YI(NYI),FI(NXI,NYI),XMSG
  DOUBLE PRECISION XO(NXYO),YO(NXYO),FO(NXYO)
  !                               local
  INTEGER N,M,NXY,NN,MM
  DOUBLE PRECISION TMP1,TMP2,SLPX,SLPY

  do NXY = 1,NXYO
     FO(NXY) = XMSG

     NN = 0
     nloop: do N = 1,NXI - 1
	if (XO(NXY) >= XI(N) .and. XO(NXY) < XI(N+1)) then
	   NN = N
	   exit nloop
	end if
     end do nloop

     MM = 0
     mloop: do M = 1,NYI - 1
	if (YO(NXY) >= YI(M) .and. YO(NXY) < YI(M+1)) then
	   MM = M
	   exit mloop
	end if
     end do mloop

     if (NN /= 0 .and. MM /= 0) then
        if (XO(NXY) == XI(NN) .and. YO(NXY) == YI(MM)) then
           !                               exact location [no interpolation]
           FO(NXY) = FI(NN,MM)
        else if (FI(NN,MM) /= XMSG .and. FI(NN+1,MM) /= XMSG .and.&
             &                 FI(NN,MM+1) /= XMSG .and.&
             &                 FI(NN+1,MM+1) /= XMSG) then
           !                               must interpolate: calculate slopes
           SLPX = (XO(NXY)-XI(NN))/ (XI(NN+1)-XI(NN))
           SLPY = (YO(NXY)-YI(MM))/ (YI(MM+1)-YI(MM))
           !                               interpolate in "x" first
           TMP1 = FI(NN,MM) + SLPX* (FI(NN+1,MM)-FI(NN,MM))
           TMP2 = FI(NN,MM+1) + SLPX* (FI(NN+1,MM+1)-FI(NN,MM+1))
           !                               interpolate in "y"
           FO(NXY) = TMP1 + SLPY* (TMP2-TMP1)
        else if (NOPT == -1) then
           CALL ESTFOW(FI(NN,MM),FI(NN+1,MM),FI(NN,MM+1),&
                &                        FI(NN+1,MM+1),XI(NN),XI(NN+1),YI(MM),&
                &                        YI(MM+1),FO(NXY),XO(NXY),YO(NXY),XMSG)
        end if
     end if

  end do

  return
END SUBROUTINE DLINT2XY
! ---------------------------------------------
SUBROUTINE ESTFOW(F1,F2,F3,F4,X1,X2,Y1,Y2,F0,X0,Y0,XMSG)
  IMPLICIT NONE
  DOUBLE PRECISION F1,F2,F3,F4,X1,X2,Y1,Y2,F0,X0,Y0,XMSG
  ! local
  DOUBLE PRECISION FI(2,2),W(2,2),SUM,SWT
  INTEGER N,M

  !     f3(x1,y2)+++++++++++f4(x2,y2)
  !        f0(x0,y0)         someplace within the surrounding 4 pts
  !     f3(x1,y1)+++++++++++f4(x2,y1)

  ! if all msg return

  F0 = XMSG
  if (F1 == XMSG .and. F2 == XMSG .and. F3 == XMSG .and.&
       &    F4 == XMSG) return

  ! compute 'distance wgted average
  ! .   make assumption that dx/dy are small and
  ! .   pythag theorem is good enough

  ! inverse distance wgts
  W(1,1) = 1.D0/SQRT((X1-X0)**2+ (Y1-Y0)**2)
  W(2,1) = 1.D0/SQRT((X2-X0)**2+ (Y1-Y0)**2)
  W(1,2) = 1.D0/SQRT((X1-X0)**2+ (Y2-Y0)**2)
  W(2,2) = 1.D0/SQRT((X2-X0)**2+ (Y2-Y0)**2)

  FI(1,1) = F1
  FI(2,1) = F2
  FI(1,2) = F3
  FI(2,2) = F4

  SUM = 0.D0
  SWT = 0.D0
  do M = 1,2
     do N = 1,2
        if (FI(N,M) /= XMSG) then
           SUM = SUM + FI(N,M)*W(N,M)
           SWT = SWT + W(N,M)
        end if
     end do
  end do
  !                               wgted average
  if (SWT > 0.D0) then
     F0 = SUM/SWT
  end if

  return
END SUBROUTINE ESTFOW
