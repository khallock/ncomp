SUBROUTINE DRCM2RGRID(NGRD,NYI,NXI,YI,XI,FI,NYO,YO,NXO,XO,FO&
     ,XMSG,NCRIT,OPT,IER)
  IMPLICIT NONE
  INTEGER          NGRD,NXI,NYI,NXO,NYO,NCRIT,OPT,IER
  DOUBLE PRECISION XI(NXI,NYI),YI(NXI,NYI),FI(NXI,NYI,NGRD)
  DOUBLE PRECISION XO(NXO),YO(NYO),FO(NXO,NYO,NGRD),XMSG
  ! real, dimension(nxo, nyo, ngrd)::FO
  ! real, dimension(nxo)::XO
  ! real, dimension(nyo)::YO

  ! NCLEND

  ! NCL:  fo = rcm2rgrid (lat2d,lon2d,fi, lat, lon iopt)
  !                        yi    xi   fi  yo   xo
  !
  !            fo is the same size xo, yo and same type as "fi"
  !            xmsg = fi@_FillValue
  !            opt unused option
  !
  !            The NCL wrapper should allow for multiple datasets
  !            so the user need only make one call to the function.

  ! perform 2D interpolation allowing for missing data:  nothing fancy

  ! nomenclature:
  ! .   nxi,nyi - lengths of xi,yi and dimensions of fi (must be >= 2)
  ! .   xi      - coordinates of fi (eg, lon [2D] )
  ! .   yi      - coordinates of fi (eg, lat [2D] )
  ! .   fi      - functional input values [2D]
  ! .   nxo,nyo - lengths of xo,yo and dimensions of fo (must be >= 1)
  ! .   xo      - coordinates of fo (eg, lon [1D])
  ! .             must be monotonically increasing
  ! .   yo      - coordinates of fo (eg, lat [1D])
  ! .             must be monotonically increasing
  ! .   fo      - functional output values [interpolated]
  ! .   xmsg    - missing code
  ! .   opt     - unused
  ! .   ier     - error code
  ! .             =0;   no error
  ! .             =1;   not enough points in input/output array
  ! .             =2/3; xi or yi are not monotonically increasing
  ! .             =4/5; xo or yo are not monotonically increasing
  !
  !                              local
  INTEGER          NG, NX,NY,NEXACT,IX,IY,M,N,NW,NER,K,NCRT
  INTEGER          MFLAG, MPTCRT, MKNT
  DOUBLE PRECISION FW(2,2),W(2,2),SUMF,SUMW,CHKLAT(NYI),CHKLON(NXI)
  DOUBLE PRECISION EPS
  DOUBLE PRECISION DGCDIST
  !                              error checking
  IER = 0
  if (NXI <= 1 .OR. NYI <= 1 .OR. NXO <= 1 .OR. NYO <= 1) then
     IER = 1
     return
  end if
  if (IER /= 0) return

  CALL DMONOINC(YO,NYO,IER,NER)
  if (IER /= 0) return
  CALL DMONOINC(XO,NXO,IER,NER)
  if (IER /= 0) return

  CHKLAT(:) = YI(1,:) ! : = NYI

  CALL DMONOINC(CHKLAT,NYI,IER,NER)
  if (IER /= 0) return

  CHKLON(:) = XI(:, 1) ! : = NXI

  CALL DMONOINC(CHKLAT,NYI,IER,NER)
  if (IER /= 0) return

  K = 2
  ! c c k = opt

  if (NCRIT <= 1) then
     NCRT = 1
  else
     NCRT = MIN(4,NCRIT)
  end if
  !                              initialize to xmsg
  FO = XMSG
  !                              main loop [exact matches]
  !                              people want bit-for-bit match
  EPS    = 1.D-04
  ! NEXACT = 0

 do NY = 1,NYO
    do NX = 1,NXO
        iyloop1: do IY = 1,NYI
           do IX = 1,NXI

              if (XO(NX) >= (XI(IX,IY)-EPS) .AND.&
                   XO(NX) <= (XI(IX,IY)+EPS) .AND.&
                   YO(NY) >= (YI(IX,IY)-EPS) .AND.&
                   YO(NY) <= (YI(IX,IY)+EPS) ) then

                 do NG=1,NGRD
                    FO(NX,NY,NG) = FI(IX,IY,NG)
                    ! NEXACT = NEXACT + 1
                 end do
                 exit iyloop1
              end if

           end do
        end do iyloop1
     end do
  end do

  ! c c print *, "nexact=",nexact
  !                              main loop [interpolation]
  do NY = 1,NYO
     do NX = 1,NXO

        iyloop2: do IY = 1,NYI-K
           do IX = 1,NXI-K
              if (XO(NX) >= XI(IX,IY) .AND.&
                   XO(NX) <= XI(IX+K,IY) .AND.&
                   YO(NY) >= YI(IX,IY) .AND.&
                   YO(NY) <= YI(IX,IY+K)) then


                 W(1,1) = (1.D0/DGCDIST(YO(NY),XO(NX),&
                      YI(IX,IY),XI(IX,IY),2))**2
                 W(2,1) = (1.D0/DGCDIST(YO(NY),XO(NX),&
                      YI(IX+K,IY),XI(IX+K,IY),2))**2
                 W(1,2) = (1.D0/DGCDIST(YO(NY),XO(NX),&
                      YI(IX,IY+K),XI(IX,IY+K),2))**2
                 W(2,2) = (1.D0/DGCDIST(YO(NY),XO(NX),&
                      YI(IX+K,IY+K),XI(IX+K,IY+K),2))**2

                 do NG=1,NGRD
                    if (FO(NX,NY,NG) == XMSG) then
                       FW(1,1) = FI(IX,IY,NG)
                       FW(2,1) = FI(IX+K,IY,NG)
                       FW(1,2) = FI(IX,IY+K,NG)
                       FW(2,2) = FI(IX+K,IY+K,NG)

                       NW   = 0
                       SUMF = 0.0D0
                       SUMW = 0.0D0
                       do N = 1,2
                          do M = 1,2
                             if (FW(M,N) /= XMSG) then
                                SUMF = SUMF + FW(M,N)*W(M,N)
                                SUMW = SUMW + W(M,N)
                                NW   = NW + 1
                             end if
                          end do
                       end do
                       !                                             nw >=3 arbitrary
                       ! c c                       if (NW >= 3 .AND. SUMW > 0.D0) then
                       !                                             nw =1 nearest neighbor
                       if (NW >= NCRT .AND. SUMW > 0.D0) then
                          FO(NX,NY,NG) = SUMF/SUMW
                       end if
                    end if
                 end do

                 exit iyloop2
              end if
           end do
        end do iyloop2

     end do
  end do

  ! Since the RCM grid is curvilinear the above algorithm may not work
  ! .   for all of the locations on regular grid. Fill via linear interp.

  MKNT   =  0
  MFLAG  =  0
  MPTCRT =  2
  do NG=1,NGRD
     do NY=1,NYO
        do NX=1,NXO
           if (FO(NX,NY,NG) == XMSG) then
              CALL DLINMSG(FO(1,NY,NG),NXO,XMSG,MFLAG,MPTCRT)
              MKNT = MKNT + 1
           end if
        end do
     end do
  end do

  ! C C PRINT *,"MKNT=",MKNT
  return
END SUBROUTINE DRCM2RGRID

! -----------------------------------------------------------

SUBROUTINE DRGRID2RCM(NGRD,NYI,NXI,YI,XI,FI,NYO,NXO,YO,XO,FO&
     ,XMSG,NCRIT,OPT,IER)
  IMPLICIT NONE
  INTEGER          NGRD,NXI,NYI,NXO,NYO,OPT,NCRIT,IER
  DOUBLE PRECISION XI(NXI),YI(NYI),FI(NXI,NYI,NGRD)
  DOUBLE PRECISION XO(NXO,NYO),YO(NXO,NYO),FO(NXO,NYO,NGRD),XMSG
  ! NCLEND

  !            fo is the same size xo, yo and same type as "fi"
  !            xmsg = fi@_FillValue
  !            opt unused option
  !
  !            The NCL wrapper should allow for multiple datasets
  !            so the user need only make one call to the function.

  ! perform 2D interpolation allowing for missing data:  nothing fancy

  ! nomenclature:
  ! .   nxi,nyi - lengths of xi,yi and dimensions of fi (must be >= 2)
  ! .   xi      - coordinates of fi (eg, lon [1D])
  ! .   yi      - coordinates of fi (eg, lat [1D])
  ! .   fi      - functional input values [2D]
  ! .   nxo,nyo - lengths of xo,yo and dimensions of fo (must be >= 1)
  ! .   xo      - coordinates of fo (eg, lon [2D])
  ! .             must be monotonically increasing
  ! .   yo      - coordinates of fo (eg, lat [2D])
  ! .             must be monotonically increasing
  ! .   fo      - functional output values [interpolated]
  ! .   xmsg    - missing code
  ! .   opt     - unused
  ! .   ier     - error code
  ! .             =0;   no error
  ! .             =1;   not enough points in input/output array
  ! .             =2/3; xi or yi are not monotonically increasing
  ! .             =4/5; xo or yo are not monotonically increasing
  !
  !                              local
  INTEGER          NG,NX,NY,NEXACT,IX,IY,M,N,NW,NER,K
  DOUBLE PRECISION FW(2,2),W(2,2),SUMF,SUMW,EPS
  DOUBLE PRECISION DGCDIST

  !                              in-line functions (bilinear interp)
  DOUBLE PRECISION Z1,Z2,Z3,Z4,SLOPE,SLPX,SLPY,FLI,FBLI

  FLI(Z1,Z2,SLOPE) = Z1 + SLOPE*(Z2-Z1)
  FBLI(Z1,Z2,Z3,Z4,SLPX,SLPY) = FLI(Z1,Z2,SLPX) +&
       SLPY* (FLI(Z3,Z4,SLPX)-&
       FLI(Z1,Z2,SLPX))

  !                              error checking
  IER = 0
  if (NXI <= 1 .OR. NYI <= 1 .OR. NXO <= 1 .OR. NYO <= 1) then
     IER = 1
     return
  end if
  if (IER /= 0) return

  CALL DMONOINC(YI,NYI,IER,NER)
  if (IER /= 0) return
  CALL DMONOINC(XI,NXI,IER,NER)
  if (IER /= 0) return
  !                              Init to missing
  FO = XMSG ! (:, :, :) = (NXO, NYO, NGRD)
  !                              main loop [exact matches]
  EPS    = 1.D-03
  NEXACT = 0

  do NY = 1,NYO
     do NX = 1,NXO

        iyloop1: do IY = 1,NYI
           do IX = 1,NXI
              if (XO(NX,NY) >= (XI(IX)-EPS) .AND.&
                   XO(NX,NY) <= (XI(IX)+EPS) .AND.&
                   YO(NX,NY) >= (YI(IY)-EPS) .AND.&
                   YO(NX,NY) <= (YI(IY)+EPS) ) then

                 do NG=1,NGRD
                    FO(NX,NY,NG) = FI(IX,IY,NG)
                    NEXACT = NEXACT + 1
                 end do
                 exit iyloop1
              end if
           end do
        end do iyloop1

     end do
  end do

  ! print *, "nexact=",nexact

  K = 1
  ! c c k = opt

  !                              main loop [interpolation]
  do NY = 1,NYO
     do NX = 1,NXO

        iyloop2: do IY = 1,NYI - K
           do IX = 1,NXI - K
              if (XO(NX,NY) >= XI(IX) .AND.&
                   XO(NX,NY) < XI(IX+K) .AND.&
                   YO(NX,NY) >= YI(IY) .AND.&
                   YO(NX,NY) < YI(IY+K)) then

                 do NG = 1,NGRD
                    if (FO(NX,NY,NG) == XMSG) then
                       if (FI(IX,IY,NG) /= XMSG .AND.&
                            FI(IX+K,IY,NG) /= XMSG .AND.&
                            FI(IX,IY+K,NG) /= XMSG .AND.&
                            FI(IX+K,IY+K,NG) /= XMSG) then

                          FO(NX,NY,NG)=FBLI(FI(IX,IY,NG),FI(IX+K,IY,NG),&
                               FI(IX,IY+K,NG),FI(IX+K,IY+K,NG),&
                               (XO(NX,NY)-XI(IX))/&
                               (XI(IX+K)-XI(IX)),&
                               (YO(NX,NY)-YI(IY))/&
                               (YI(IY+K)-YI(IY)))

                       else
                          !                                            OVERKILL
                          FW(1,1) = FI(IX,IY,NG)
                          FW(2,1) = FI(IX+K,IY,NG)
                          FW(1,2) = FI(IX,IY+K,NG)
                          FW(2,2) = FI(IX+K,IY+K,NG)

                          W(1,1) = (1.D0/DGCDIST(YO(NX,NY),XO(NX,NY)&
                               ,YI(IY),XI(IX),2))**2
                          W(2,1) = (1.D0/DGCDIST(YO(NX,NY),XO(NX,NY)&
                               ,YI(IY),XI(IX+K),2))**2
                          W(1,2) = (1.D0/DGCDIST(YO(NX,NY),XO(NX,NY)&
                               ,YI(IY+K),XI(IX),2))**2
                          W(2,2) = (1.D0/DGCDIST(YO(NX,NY),XO(NX,NY)&
                               ,YI(IY+K),XI(IX+K),2))**2

                          NW = 0
                          SUMF = 0.0D0
                          SUMW = 0.0D0
                          do N = 1,2
                             do M = 1,2
                                if (FW(M,N) /= XMSG) then
                                   SUMF = SUMF + FW(M,N)*W(M,N)
                                   SUMW = SUMW + W(M,N)
                                   NW = NW + 1
                                end if
                             end do
                          end do
                          !                                             nw >=3 arbitrary
                          ! c c                  if (NCRIT >= 3 .AND. SUMW > 0.D0) then
                          !                                             nw  =1 nearest neighbor
                          if (NCRIT >= 1 .AND. SUMW > 0.D0) then
                             FO(NX,NY,NG) = SUMF/SUMW
                          end if
                       end if
                    end if
                 end do
                 exit iyloop2
              end if
           end do
        end do iyloop2

     end do
  end do

  return
END SUBROUTINE DRGRID2RCM

! -----------------------------------------------------------

DOUBLE PRECISION FUNCTION DGCDIST(RLAT1,RLON1,RLAT2,RLON2,IU)
  IMPLICIT NONE
  !
  ! calculate the great circle distance between two points
  !
  ! usage: dist = gcdist (rlat1,rlon1,rlat2,rlon2,iu)
  !
  ! nomenclature :
  ! .   rlat1,rlon1 - latitude and longtitude of the first point
  ! .   rlat2,rlon2 - latitude and longtitude of the second point
  ! .   iu          - code for the type units gcdist is to return
  ! .               = 1 : gcdist returned in radians
  ! .               = 2 : gcdist returned in degrees
  ! .               = 3 : gcdist returned in meters
  ! .               = 4 : gcdist returned in kilometers
  ! .               = 5 : gcdist returned in *not used*
  !
  ! input
  INTEGER IU
  ! input types
  DOUBLE PRECISION RLAT1,RLON1,RLAT2,RLON2

  ! local stuff
  DOUBLE PRECISION UNITS(5),RAD,DLONR,RLAT1R,RLAT2R
  DATA UNITS/1.0D0,57.29577995691645D0,6371220.D0,6371.2200D0,0.D0/
  ! change as required
  DATA RAD/0.01745329238474369D0/

  ! special test if RLAT1=RLAT2 and RLON1=RLON2
  IF(RLAT1 == RLAT2 .AND. RLON1==RLON2) THEN
     DGCDIST = 0.D0
     RETURN
  END IF
  RLAT1R = RLAT1*RAD
  RLAT2R = RLAT2*RAD
  DLONR = DMIN1(ABS(RLON1-RLON2),ABS(360.D0-RLON1+RLON2),&
       ABS(360.D0-RLON2+RLON1))*RAD

  DGCDIST = ATAN2(SQRT((COS(RLAT2R) * SIN(DLONR)) ** 2 +&
       (COS(RLAT1R) * SIN(RLAT2R) -&
       SIN(RLAT1R) * COS(RLAT2R) * COS(DLONR)) ** 2&
       ),&
       SIN(RLAT1R)*SIN(RLAT2R)+&
       COS(RLAT1R)*COS(RLAT2R)*COS(DLONR)&
       ) * UNITS(IU)

  RETURN
END FUNCTION DGCDIST
